module elmxxHistMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Minimal monthly (h0-equivalent) history output for ELMxx.
  !
  ! Thinner than a port of histFileMod.F90: per-step accumulation lives in
  ! C++ (History.cpp, driven through elmxx_mod's ELMxxHistory* bindings).
  ! This module's job is exactly histFileMod's hist_htapes_wrapup minus its
  ! per-step half: parse the namelist additions, drive one accumulate call
  ! per step, and once a month pull the normalized means and write them
  ! through a real PIO decomposed write.
  !
  ! See plans/active/2026_08_23_elmxx_history_output_plan.md for the design.
  !-----------------------------------------------------------------------

  use shr_kind_mod , only : r8 => shr_kind_r8, r4 => shr_kind_r4
  use shr_sys_mod  , only : shr_sys_abort, shr_sys_flush
  use shr_const_mod, only : SHR_CONST_REARTH
  use shr_cal_mod  , only : shr_cal_ymd2julian, shr_cal_numDaysInYear, shr_cal_noleap
  use elmxxSpmdMod , only : masterproc, iam, npes
  use elmxxIO      , only : pio_subsystem, io_type
  use elmxxSubgridMod, only : num_landunits, lun_itype, lun_wtgcell, istsoil
  use elmxx_mod    , only : ELMxxType, ELMXX_SUCCESS, &
                             ELMxxHistoryActivateDefaults, ELMxxHistoryActivate, &
                             ELMxxHistoryActiveCount, ELMxxHistoryActiveName, &
                             ELMxxHistoryActiveUnits, ELMxxHistoryActiveLongName, &
                             ELMxxHistoryActiveIsAvg, ELMxxHistoryAccumulate, &
                             ELMxxHistoryNormalize, ELMxxHistoryGetField, &
                             ELMxxHistoryReset
  use pio

  implicit none
  save
  private

  public :: elmxx_hist_init
  public :: elmxx_hist_step
  public :: elmxx_hist_write_if_month_end
  public :: elmxx_hist_final

  !--------------------------------------------------------------------------
  ! Module state. One instance only -- ELMxx runs one land model per
  ! process, same assumption elmxxMod itself makes for elmxx_state.
  !--------------------------------------------------------------------------
  logical            :: hist_init_done = .false.
  character(len=256) :: hist_caseid    = ' '
  integer            :: n_local        = 0   ! local (this-rank) natural-column count
  integer            :: num_global     = 0   ! lndgrid dimension size
  integer            :: nfields        = 0

  character(len=32) , allocatable :: field_name(:)
  character(len=128), allocatable :: field_units(:)
  character(len=256), allocatable :: field_long_name(:)
  logical            , allocatable :: field_is_avg(:)
  real(r8)           , allocatable :: hist_buf(:,:)   ! (n_local, nfields), pulled at month-end

  real(r8), allocatable :: hist_lon(:), hist_lat(:), hist_area(:)  ! (n_local), time-constant
  integer , allocatable :: hist_dof(:)                             ! (n_local), 1-based lndgrid position

  type(io_desc_t) :: hist_iodesc

  integer :: hist_start_year  = 0
  integer :: hist_start_month = 0
  integer :: hist_start_day   = 0
  integer :: hist_nacs        = 0   ! this module's own step counter, mirrors C++ hist.nacs

  real(r4), parameter :: HIST_FILL = 1.0e36_r4

contains

  !-----------------------------------------------------------------------
  subroutine elmxx_hist_init(elmxx_state, logunit, caseid, elmxx_hist_fincl, &
                             hist_year, month, day, n_local_cells, num_cells_global_in, &
                             natural_id_cells_owned, lonc_g, latc_g, areac_g)
    !
    ! Activate the C++-side default fields plus any namelist additions,
    ! build the local staging buffers and the PIO decomposed-write map, and
    ! record the interval start.
    !
    implicit none
    type(ELMxxType)    , intent(in) :: elmxx_state
    integer             , intent(in) :: logunit
    character(len=*)    , intent(in) :: caseid
    character(len=*)    , intent(in) :: elmxx_hist_fincl
    integer             , intent(in) :: hist_year, month, day
    integer             , intent(in) :: n_local_cells
    integer             , intent(in) :: num_cells_global_in
    integer             , intent(in) :: natural_id_cells_owned(:)  ! (n_local_cells,), full-grid indices
    real(r8)            , intent(in) :: lonc_g(:), latc_g(:), areac_g(:)  ! (ni*nj,)

    integer :: ierr, i, is_avg
    real(r8) :: re_km
    character(len=*), parameter :: subname = '(elmxx_hist_init) '

    if (hist_init_done) call shr_sys_abort(subname//'ERROR: elmxx_hist_init called twice')

    hist_caseid = caseid
    n_local     = n_local_cells
    num_global  = num_cells_global_in

    call ELMxxHistoryActivateDefaults(elmxx_state, ierr)
    if (ierr /= ELMXX_SUCCESS) &
         call shr_sys_abort(subname//'ERROR: ELMxxHistoryActivateDefaults failed')

    call elmxx_hist_fincl_parse(elmxx_state, elmxx_hist_fincl, logunit)

    nfields = ELMxxHistoryActiveCount(elmxx_state, ierr)
    if (ierr /= ELMXX_SUCCESS) &
         call shr_sys_abort(subname//'ERROR: ELMxxHistoryActiveCount failed')

    if (nfields == 0) then
       if (masterproc) write(logunit,*) subname,'no history fields active; history output disabled'
       return
    end if

    allocate(field_name(nfields), field_units(nfields), field_long_name(nfields), &
             field_is_avg(nfields))
    do i = 1, nfields
       field_name(i) = ' '; field_units(i) = ' '; field_long_name(i) = ' '
       call ELMxxHistoryActiveName(elmxx_state, i-1, field_name(i), ierr)
       if (ierr /= ELMXX_SUCCESS) call shr_sys_abort(subname//'ERROR: ELMxxHistoryActiveName failed')
       call ELMxxHistoryActiveUnits(elmxx_state, i-1, field_units(i), ierr)
       if (ierr /= ELMXX_SUCCESS) call shr_sys_abort(subname//'ERROR: ELMxxHistoryActiveUnits failed')
       call ELMxxHistoryActiveLongName(elmxx_state, i-1, field_long_name(i), ierr)
       if (ierr /= ELMXX_SUCCESS) call shr_sys_abort(subname//'ERROR: ELMxxHistoryActiveLongName failed')
       is_avg = ELMxxHistoryActiveIsAvg(elmxx_state, i-1, ierr)
       if (ierr /= ELMXX_SUCCESS) call shr_sys_abort(subname//'ERROR: ELMxxHistoryActiveIsAvg failed')
       field_is_avg(i) = (is_avg /= 0)
    end do

    allocate(hist_buf(n_local, nfields))
    hist_buf = 0.0_r8

    ! Natural-vegetation-landunit-weight-on-gridcell == 1 guard: PATCH_TO_COL
    ! fields are natural-column values, which coincide with ELM's gridcell h0
    ! values only when the natural-veg landunit covers the whole cell (true
    ! on 1x1_brazil/1x1_glc, not in general -- see the plan's "This is
    ! natural-landunit output, not gridcell output" section).
    call elmxx_hist_check_natveg_weight(logunit)

    ! Time-constant per-cell fields, looked up once via the full-grid arrays.
    ! area converts the coupler's radians^2 to h0's conventional km^2, the
    ! same re*re conversion CIME's domain/history tooling uses.
    re_km = SHR_CONST_REARTH * 1.0e-3_r8
    allocate(hist_lon(n_local), hist_lat(n_local), hist_area(n_local))
    do i = 1, n_local
       hist_lon(i)  = lonc_g(natural_id_cells_owned(i))
       hist_lat(i)  = latc_g(natural_id_cells_owned(i))
       hist_area(i) = areac_g(natural_id_cells_owned(i)) * re_km * re_km
    end do

    ! Decomposed-write map: rank iam's i-th owned cell sits at the
    ! iam+1+(i-1)*npes ascending position among the num_global active land
    ! cells -- the compact lndgrid ordinal, NOT natural_id_cells_owned(i)
    ! (which indexes the sparse, ocean-inclusive full grid). See the plan's
    ! "PIO schema" section for why the two are easy to conflate.
    allocate(hist_dof(n_local))
    do i = 1, n_local
       hist_dof(i) = iam + 1 + (i-1)*npes
    end do
    call PIO_initdecomp(pio_subsystem, PIO_real, (/num_global/), hist_dof, hist_iodesc)

    hist_start_year  = hist_year
    hist_start_month = month
    hist_start_day   = day
    hist_nacs        = 0

    hist_init_done = .true.

    if (masterproc) then
       write(logunit,*) subname,nfields,' history fields active:'
       do i = 1, nfields
          write(logunit,*) subname,'   ',trim(field_name(i)),' [',trim(field_units(i)),'] ', &
               trim(field_long_name(i))
       end do
       call shr_sys_flush(logunit)
    end if

  end subroutine elmxx_hist_init

  !-----------------------------------------------------------------------
  subroutine elmxx_hist_fincl_parse(elmxx_state, spec, logunit)
    !
    ! elmxx_hist_fincl names fields to activate BEYOND the C++-side
    ! defaults -- additions/overrides, not a full restatement. Same
    ! comma-separated-string shape as elmxx_kernels_parse
    ! (elmxxKernelMod.F90), reused deliberately.
    !
    implicit none
    type(ELMxxType) , intent(in) :: elmxx_state
    character(len=*), intent(in) :: spec
    integer          , intent(in) :: logunit

    integer :: ib, ie, n, ierr
    character(len=32) :: token
    character(len=*), parameter :: subname = '(elmxx_hist_fincl_parse) '

    n = len_trim(spec)
    if (n == 0) return

    ib = 1
    do while (ib <= n)
       ie = index(spec(ib:n), ',')
       if (ie == 0) then
          ie = n
       else
          ie = ib + ie - 2
       end if

       token = adjustl(spec(ib:ie))
       if (len_trim(token) > 0) then
          call ELMxxHistoryActivate(elmxx_state, trim(token), ierr)
          if (ierr /= ELMXX_SUCCESS) then
             write(logunit,*) subname,'ERROR: could not activate "',trim(token), &
                  '" from elmxx_hist_fincl -- unknown field name or already active'
             call shr_sys_flush(logunit)
             call shr_sys_abort(subname//'ERROR: elmxx_hist_fincl names an unknown '// &
                  'or duplicate field: "'//trim(token)//'"')
          end if
       end if

       ib = ie + 2
    end do

  end subroutine elmxx_hist_fincl_parse

  !-----------------------------------------------------------------------
  subroutine elmxx_hist_check_natveg_weight(logunit)
    !
    ! See elmxx_hist_init's comment at the call site.
    !
    implicit none
    integer, intent(in) :: logunit
    integer :: l
    real(r8), parameter :: tol = 1.0e-6_r8
    character(len=*), parameter :: subname = '(elmxx_hist_check_natveg_weight) '

    do l = 1, num_landunits
       if (lun_itype(l) == istsoil) then
          if (abs(lun_wtgcell(l) - 1.0_r8) > tol) then
             write(logunit,*) subname,'ERROR: natural-vegetation landunit weight-on-', &
                  'gridcell = ',lun_wtgcell(l),' != 1.0 at local landunit ',l,' -- ', &
                  'PATCH_TO_COL history fields are natural-column values, and only ', &
                  'coincide with ELM''s gridcell h0 values when this weight is 1.0'
             call shr_sys_flush(logunit)
             call shr_sys_abort(subname//'ERROR: natural-veg landunit weight-on-gridcell '// &
                  '/= 1.0; history output would silently mislabel natural-column values '// &
                  'as gridcell values')
          end if
       end if
    end do

  end subroutine elmxx_hist_check_natveg_weight

  !-----------------------------------------------------------------------
  subroutine elmxx_hist_step(elmxx_state, logunit)
    !
    ! One call per internal step. Pure device work on the C++ side -- see
    ! the plan for why this must run after phenology/SurfaceAlbedo, not
    ! before.
    !
    implicit none
    type(ELMxxType), intent(in) :: elmxx_state
    integer         , intent(in) :: logunit
    integer :: ierr
    character(len=*), parameter :: subname = '(elmxx_hist_step) '

    if (.not. hist_init_done) return

    call ELMxxHistoryAccumulate(elmxx_state, ierr)
    if (ierr /= ELMXX_SUCCESS) call shr_sys_abort(subname//'ERROR: ELMxxHistoryAccumulate failed')
    hist_nacs = hist_nacs + 1

  end subroutine elmxx_hist_step

  !-----------------------------------------------------------------------
  subroutine elmxx_hist_write_if_month_end(elmxx_state, year, month, day, tod, logunit)
    !
    ! (year, month, day, tod) is the END time of the step just taken. A
    ! month boundary is day==1 .and. tod==0 of that end time; the file
    ! written is dated for the month that just closed (the interval's
    ! START, recorded at elmxx_hist_init or the previous call here), not
    ! the month whose day-1 boundary triggered the write.
    !
    implicit none
    type(ELMxxType), intent(in) :: elmxx_state
    integer         , intent(in) :: year, month, day, tod
    integer         , intent(in) :: logunit

    integer :: ierr, ifld
    character(len=*), parameter :: subname = '(elmxx_hist_write_if_month_end) '

    if (.not. hist_init_done) return
    if (.not. (day == 1 .and. tod == 0)) return

    if (hist_nacs == 0) then
       write(logunit,*) subname,'WARNING: month boundary reached with no accumulated ', &
            'samples; skipping write for ',hist_start_year,'-',hist_start_month
       call shr_sys_flush(logunit)
       hist_start_year = year; hist_start_month = month; hist_start_day = day
       return
    end if

    call ELMxxHistoryNormalize(elmxx_state, ierr)
    if (ierr /= ELMXX_SUCCESS) call shr_sys_abort(subname//'ERROR: ELMxxHistoryNormalize failed')

    do ifld = 1, nfields
       call ELMxxHistoryGetField(elmxx_state, trim(field_name(ifld)), hist_buf(:,ifld), &
            n_local, ierr)
       if (ierr /= ELMXX_SUCCESS) &
            call shr_sys_abort(subname//'ERROR: ELMxxHistoryGetField failed for '// &
            trim(field_name(ifld)))
    end do

    call elmxx_hist_write_file(hist_start_year, hist_start_month, hist_start_day, &
                               year, month, day, tod, logunit)

    call ELMxxHistoryReset(elmxx_state, ierr)
    if (ierr /= ELMXX_SUCCESS) call shr_sys_abort(subname//'ERROR: ELMxxHistoryReset failed')
    hist_nacs = 0

    ! The interval that just closed started at hist_start_*; the next one
    ! starts exactly here, at this month boundary.
    hist_start_year = year; hist_start_month = month; hist_start_day = day

  end subroutine elmxx_hist_write_if_month_end

  !-----------------------------------------------------------------------
  subroutine elmxx_hist_write_file(start_year, start_month, start_day, &
                                   end_year, end_month, end_day, end_tod, logunit)
    !
    ! One self-contained monthly NetCDF file: dims, time-constant fields and
    ! the accumulating fields are all (re)defined here, since -- unlike
    ! ELM's single always-open tape -- a new file is created every month.
    !
    implicit none
    integer, intent(in) :: start_year, start_month, start_day
    integer, intent(in) :: end_year, end_month, end_day, end_tod
    integer, intent(in) :: logunit

    type(file_desc_t) :: file
    type(var_desc_t)  :: vid_lon, vid_lat, vid_area, vid_time, vid_time_bounds
    type(var_desc_t), allocatable :: vid_field(:)
    integer :: dimid_lndgrid, dimid_time, dimid_histint
    integer :: ierr, ifld
    character(len=256) :: fname
    character(len=16)  :: cdate
    real(r8) :: start_days, end_days, tb(2,1), tval(1)
    real(r4), allocatable :: wbuf(:)
    character(len=*), parameter :: subname = '(elmxx_hist_write_file) '

    write(cdate,'(i4.4,"-",i2.2)') start_year, start_month
    fname = trim(hist_caseid)//'.elmxx.h0.'//trim(cdate)//'.nc'

    ierr = PIO_createfile(pio_subsystem, file, io_type, trim(fname), PIO_CLOBBER)
    if (ierr /= PIO_NOERR) call shr_sys_abort(subname//'ERROR: PIO_createfile failed for '//trim(fname))

    ierr = PIO_def_dim(file, 'lndgrid', num_global, dimid_lndgrid)
    if (ierr /= PIO_NOERR) call shr_sys_abort(subname//'ERROR: PIO_def_dim lndgrid failed')
    ierr = PIO_def_dim(file, 'time', PIO_UNLIMITED, dimid_time)
    if (ierr /= PIO_NOERR) call shr_sys_abort(subname//'ERROR: PIO_def_dim time failed')
    ierr = PIO_def_dim(file, 'hist_interval', 2, dimid_histint)
    if (ierr /= PIO_NOERR) call shr_sys_abort(subname//'ERROR: PIO_def_dim hist_interval failed')

    ! Time-constant per-cell fields.
    ierr = PIO_def_var(file, 'lon', PIO_real, (/dimid_lndgrid/), vid_lon)
    if (ierr /= PIO_NOERR) call shr_sys_abort(subname//'ERROR: PIO_def_var lon failed')
    ierr = PIO_put_att(file, vid_lon, 'units', 'degrees_east')
    ierr = PIO_def_var(file, 'lat', PIO_real, (/dimid_lndgrid/), vid_lat)
    if (ierr /= PIO_NOERR) call shr_sys_abort(subname//'ERROR: PIO_def_var lat failed')
    ierr = PIO_put_att(file, vid_lat, 'units', 'degrees_north')
    ierr = PIO_def_var(file, 'area', PIO_real, (/dimid_lndgrid/), vid_area)
    if (ierr /= PIO_NOERR) call shr_sys_abort(subname//'ERROR: PIO_def_var area failed')
    ierr = PIO_put_att(file, vid_area, 'units', 'km^2')

    ! Time / time_bounds. NOT part of the lndgrid decomposition -- one
    ! value (this file's single record) written collectively below.
    ierr = PIO_def_var(file, 'time', PIO_double, (/dimid_time/), vid_time)
    if (ierr /= PIO_NOERR) call shr_sys_abort(subname//'ERROR: PIO_def_var time failed')
    ierr = PIO_put_att(file, vid_time, 'units', 'days since 0001-01-01 00:00:00')
    ierr = PIO_put_att(file, vid_time, 'calendar', trim(shr_cal_noleap))
    ierr = PIO_put_att(file, vid_time, 'bounds', 'time_bounds')
    ierr = PIO_def_var(file, 'time_bounds', PIO_double, (/dimid_histint, dimid_time/), vid_time_bounds)
    if (ierr /= PIO_NOERR) call shr_sys_abort(subname//'ERROR: PIO_def_var time_bounds failed')

    ! Accumulating fields.
    allocate(vid_field(nfields))
    do ifld = 1, nfields
       ierr = PIO_def_var(file, trim(field_name(ifld)), PIO_real, (/dimid_lndgrid, dimid_time/), &
            vid_field(ifld))
       if (ierr /= PIO_NOERR) &
            call shr_sys_abort(subname//'ERROR: PIO_def_var failed for '//trim(field_name(ifld)))
       ierr = PIO_put_att(file, vid_field(ifld), 'units', trim(field_units(ifld)))
       ierr = PIO_put_att(file, vid_field(ifld), 'long_name', trim(field_long_name(ifld)))
       ierr = PIO_put_att(file, vid_field(ifld), '_FillValue', HIST_FILL)
       ierr = PIO_put_att(file, vid_field(ifld), 'missing_value', HIST_FILL)
       if (field_is_avg(ifld)) then
          ierr = PIO_put_att(file, vid_field(ifld), 'cell_methods', 'time: mean')
       end if
    end do

    ierr = PIO_enddef(file)
    if (ierr /= PIO_NOERR) call shr_sys_abort(subname//'ERROR: PIO_enddef failed for '//trim(fname))

    ! Time / time_bounds: absolute days since the twin case's own epoch
    ! (0001-01-01 00:00:00, noleap) -- not run-relative, so an ELMxx h0
    ! file overlays its ELM twin without a conversion step.
    start_days = elmxx_hist_days_since_epoch(start_year, start_month, start_day, 0)
    end_days   = elmxx_hist_days_since_epoch(end_year, end_month, end_day, end_tod)
    tval(1)   = end_days
    tb(1,1)   = start_days
    tb(2,1)   = end_days
    ierr = PIO_put_var(file, vid_time, (/1/), tval(1))
    if (ierr /= PIO_NOERR) call shr_sys_abort(subname//'ERROR: PIO_put_var time failed')
    ierr = PIO_put_var(file, vid_time_bounds, (/1,1/), (/2,1/), tb)
    if (ierr /= PIO_NOERR) call shr_sys_abort(subname//'ERROR: PIO_put_var time_bounds failed')

    ! Time-constant, decomposed fields -- written into every file (each
    ! monthly file is self-contained), not once for the whole run.
    allocate(wbuf(n_local))
    wbuf = real(hist_lon, r4)
    call PIO_write_darray(file, vid_lon, hist_iodesc, wbuf, ierr)
    if (ierr /= PIO_NOERR) call shr_sys_abort(subname//'ERROR: PIO_write_darray lon failed')
    wbuf = real(hist_lat, r4)
    call PIO_write_darray(file, vid_lat, hist_iodesc, wbuf, ierr)
    if (ierr /= PIO_NOERR) call shr_sys_abort(subname//'ERROR: PIO_write_darray lat failed')
    wbuf = real(hist_area, r4)
    call PIO_write_darray(file, vid_area, hist_iodesc, wbuf, ierr)
    if (ierr /= PIO_NOERR) call shr_sys_abort(subname//'ERROR: PIO_write_darray area failed')

    do ifld = 1, nfields
       call PIO_setframe(file, vid_field(ifld), 1_PIO_OFFSET_KIND, ierr)
       if (ierr /= PIO_NOERR) &
            call shr_sys_abort(subname//'ERROR: PIO_setframe failed for '//trim(field_name(ifld)))
       wbuf = real(hist_buf(:,ifld), r4)
       call PIO_write_darray(file, vid_field(ifld), hist_iodesc, wbuf, ierr)
       if (ierr /= PIO_NOERR) &
            call shr_sys_abort(subname//'ERROR: PIO_write_darray failed for '//trim(field_name(ifld)))
    end do
    deallocate(wbuf)
    deallocate(vid_field)

    call PIO_closefile(file)

    if (masterproc) then
       write(logunit,*) subname,'wrote ',trim(fname),' (',nfields,' fields, ',hist_nacs,' samples)'
       call shr_sys_flush(logunit)
    end if

  end subroutine elmxx_hist_write_file

  !-----------------------------------------------------------------------
  function elmxx_hist_days_since_epoch(year, month, day, tod) result(days)
    !
    ! Days since 0001-01-01 00:00:00, calendar-aware (not hardcoded to
    ! noleap's 365, even though that is the only calendar the target twins
    ! use today -- see the plan's time-metadata note).
    !
    implicit none
    integer, intent(in) :: year, month, day, tod
    real(r8) :: days
    real(r8) :: jday
    integer  :: y

    days = 0.0_r8
    do y = 1, year - 1
       days = days + real(shr_cal_numDaysInYear(y, shr_cal_noleap), r8)
    end do
    call shr_cal_ymd2julian(year, month, day, tod, jday, shr_cal_noleap)
    days = days + (jday - 1.0_r8)

  end function elmxx_hist_days_since_epoch

  !-----------------------------------------------------------------------
  subroutine elmxx_hist_final(elmxx_state, year, month, day, tod, logunit)
    !
    ! Flush a partial trailing month (nacs > 0) so a run that stops
    ! mid-month still gets a written record. Must run BEFORE ELMxxDestroy
    ! -- elmxx_final calls this first, ahead of tearing the object down.
    !
    implicit none
    type(ELMxxType), intent(in) :: elmxx_state
    integer         , intent(in) :: year, month, day, tod
    integer         , intent(in) :: logunit
    integer :: ierr

    if (.not. hist_init_done) return

    if (hist_nacs > 0) then
       call ELMxxHistoryNormalize(elmxx_state, ierr)
       if (ierr /= ELMXX_SUCCESS) &
            call shr_sys_abort('(elmxx_hist_final) ERROR: ELMxxHistoryNormalize failed')
       block
         integer :: ifld
         do ifld = 1, nfields
            call ELMxxHistoryGetField(elmxx_state, trim(field_name(ifld)), hist_buf(:,ifld), &
                 n_local, ierr)
            if (ierr /= ELMXX_SUCCESS) &
                 call shr_sys_abort('(elmxx_hist_final) ERROR: ELMxxHistoryGetField failed for '// &
                 trim(field_name(ifld)))
         end do
       end block
       call elmxx_hist_write_file(hist_start_year, hist_start_month, hist_start_day, &
                                  year, month, day, tod, logunit)
       call ELMxxHistoryReset(elmxx_state, ierr)
       hist_nacs = 0
    end if

    call PIO_freedecomp(pio_subsystem, hist_iodesc)
    hist_init_done = .false.

  end subroutine elmxx_hist_final

end module elmxxHistMod
