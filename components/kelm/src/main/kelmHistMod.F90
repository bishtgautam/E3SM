module kelmHistMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Minimal h0-equivalent history output for kELM: monthly by default, or
  ! every N hours (kelm_hist_nhtfrq < 0), kelm_hist_mfilt records per file.
  !
  ! Thinner than a port of histFileMod.F90: per-step accumulation lives in
  ! C++ (History.cpp, driven through kelm_mod's kELMHistory* bindings).
  ! This module's job is exactly histFileMod's hist_htapes_wrapup minus its
  ! per-step half: parse the namelist additions, drive one accumulate call
  ! per step, and at each interval's end pull the normalized means and
  ! write them through a real PIO decomposed write.
  !
  ! See plans/active/2026_08_23_kelm_history_output_plan.md for the design.
  !-----------------------------------------------------------------------

  use shr_kind_mod , only : r8 => shr_kind_r8, r4 => shr_kind_r4
  use shr_sys_mod  , only : shr_sys_abort, shr_sys_flush
  use shr_const_mod, only : SHR_CONST_REARTH
  use shr_cal_mod  , only : shr_cal_ymd2julian, shr_cal_numDaysInYear, shr_cal_noleap
  use kelmSpmdMod , only : masterproc, iam, npes
  use kelmIO      , only : pio_subsystem, io_type, io_format
  use kelmSubgridMod, only : num_landunits, lun_itype, lun_wtgcell, lun_gridcell, &
                             istsoil, istdlak, istice, isturb_tbd, isturb_md, &
                             num_patches, patch_column, patch_itype, patch_wtcol, &
                             col_landunit
  use kelmSurfdataMod, only : natpft
  use kelm_mod    , only : kELMType, KELM_SUCCESS, &
                             kELMHistoryActivateDefaults, kELMHistoryActivate, &
                             kELMHistoryActiveCount, kELMHistoryActiveName, &
                             kELMHistoryActiveUnits, kELMHistoryActiveLongName, &
                             kELMHistoryActiveIsAvg, kELMHistoryAccumulate, &
                             kELMHistoryNormalize, kELMHistoryGetField, &
                             kELMHistoryReset, kELMHistoryActiveNlev, &
                             kELMHistoryActiveLevDim
  use pio

  implicit none
  save
  private

  public :: kelm_hist_init
  public :: kelm_hist_step
  public :: kelm_hist_write_if_due
  public :: kelm_hist_final
  public :: kelm_hist_get_restart_state
  public :: kelm_hist_set_restart_state
  public :: kelm_hist_set_base_date
  public :: kelm_hist_check_natveg_weight

  !--------------------------------------------------------------------------
  ! Module state. One instance only -- kELM runs one land model per
  ! process, same assumption kelmMod itself makes for kelm_state.
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
  integer            , allocatable :: field_nlev(:)   ! 0, or the field's level extent
  integer            , allocatable :: field_levdim(:) ! 0, or its index into levdim_*

  ! A multi-level field's month-end values, (n_local, nlev). Single-level
  ! fields stay in hist_buf. Each level dimension the file carries (levgrnd,
  ! levlak, and the subgrid ltype/natpft) has one size and one decomposition.
  type :: hist_lev_buf_t
     real(r8), allocatable :: v(:,:)
  end type hist_lev_buf_t
  type(hist_lev_buf_t), allocatable :: hist_buf2(:)    ! (nfields)
  integer, parameter :: MAX_LEVDIM = 8
  integer            :: nlevdim = 0
  character(len=16)  :: levdim_name(MAX_LEVDIM) = ' '
  integer            :: levdim_n(MAX_LEVDIM) = 0
  type(io_desc_t)    :: levdim_iodesc(MAX_LEVDIM)

  ! ELM's time-constant subgrid fractions, written into every file:
  ! PCT_LANDUNIT on ltype (ELM's max_lunit landunit types) and PCT_NAT_PFT
  ! on natpft, both (n_local, n) in percent.
  integer, parameter :: MAX_LUNIT = 9
  integer :: ilev_ltype = 0, ilev_natpft = 0
  real(r8), allocatable :: pct_landunit(:,:), pct_nat_pft(:,:)

  ! The run's start date (the driver clock's startTime) as days since
  ! 0001-01-01, ELM's base day for mdcur/mscur.
  real(r8) :: hist_base_days = 0.0_r8

  real(r8), allocatable :: hist_lon(:), hist_lat(:), hist_area(:)  ! (n_local), time-constant
  integer , allocatable :: hist_dof(:)                             ! (n_local), 1-based lndgrid position

  type(io_desc_t) :: hist_iodesc

  integer :: hist_start_year  = 0
  integer :: hist_start_month = 0
  integer :: hist_start_day   = 0
  integer :: hist_nacs        = 0   ! this module's own step counter, mirrors C++ hist.nacs
  ! The open interval's start time of day [s]. Not in the restart file: a
  ! restart resumes it as 0, which is exact for monthly and daily output.
  integer :: hist_start_tod   = 0

  ! Output frequency and file size, ELM's hist_nhtfrq/hist_mfilt for tape 1:
  ! 0 = monthly; < 0 = every -hist_nhtfrq hours (a divisor of 24).
  integer :: hist_nhtfrq = 0
  integer :: hist_mfilt  = 1

  ! The file being filled, open from its first record until it holds
  ! hist_mfilt records (or the run ends).
  logical            :: hist_file_open = .false.
  integer            :: hist_nrec      = 0
  character(len=256) :: hist_fname     = ' '
  type(file_desc_t)  :: hist_file
  type(var_desc_t)   :: hv_time, hv_time_bounds, hv_mcdate, hv_mcsec, hv_mdcur, hv_mscur
  type(var_desc_t), allocatable :: hv_field(:)

  real(r4), parameter :: HIST_FILL = 1.0e36_r4

contains

  !-----------------------------------------------------------------------
  ! Restart (Stage 6): the Fortran half of the history state -- the open
  ! interval's start date and this module's own sample count. The device
  ! accumulators and the C++ nacs cross through the restart field registry.
  !-----------------------------------------------------------------------
  subroutine kelm_hist_get_restart_state(nacs, start_year, start_month, start_day)
    integer, intent(out) :: nacs, start_year, start_month, start_day
    nacs = hist_nacs
    start_year = hist_start_year; start_month = hist_start_month
    start_day = hist_start_day
  end subroutine kelm_hist_get_restart_state

  subroutine kelm_hist_set_restart_state(nacs, start_year, start_month, start_day)
    integer, intent(in) :: nacs, start_year, start_month, start_day
    hist_nacs = nacs
    hist_start_year = start_year; hist_start_month = start_month
    hist_start_day = start_day
  end subroutine kelm_hist_set_restart_state

  !-----------------------------------------------------------------------
  ! The run's start date, ELM's base day for mdcur/mscur. The driver clock's
  ! startTime is the case's RUN_STARTDATE on a continue run too, so a
  ! restarted run writes the same mdcur as a straight one.
  !-----------------------------------------------------------------------
  subroutine kelm_hist_set_base_date(year, month, day, tod)
    integer, intent(in) :: year, month, day, tod
    hist_base_days = kelm_hist_days_since_epoch(year, month, day, tod)
  end subroutine kelm_hist_set_base_date

  !-----------------------------------------------------------------------
  ! Index of level dimension `name` with `n` levels, added on first use.
  !-----------------------------------------------------------------------
  function kelm_hist_levdim(name, n) result(k)
    character(len=*), intent(in) :: name
    integer, intent(in) :: n
    integer :: k
    do k = 1, nlevdim
       if (trim(levdim_name(k)) == trim(name)) then
          if (levdim_n(k) /= n) call shr_sys_abort('(kelm_hist_levdim) ERROR: level '// &
               'dimension '//trim(name)//' used with two different sizes')
          return
       end if
    end do
    if (nlevdim == MAX_LEVDIM) call shr_sys_abort('(kelm_hist_levdim) ERROR: too many '// &
         'level dimensions')
    nlevdim = nlevdim + 1
    k = nlevdim
    levdim_name(k) = name
    levdim_n(k) = n
  end function kelm_hist_levdim

  !-----------------------------------------------------------------------
  subroutine kelm_hist_init(kelm_state, logunit, caseid, kelm_hist_fincl, &
                             hist_year, month, day, n_local_cells, num_cells_global_in, &
                             natural_id_cells_owned, lonc_g, latc_g, areac_g, &
                             nhtfrq, mfilt, empty_htapes)
    !
    ! Activate the C++-side default fields plus any namelist additions,
    ! build the local staging buffers and the PIO decomposed-write map, and
    ! record the interval start.
    !
    implicit none
    type(kELMType)    , intent(in) :: kelm_state
    integer             , intent(in) :: logunit
    character(len=*)    , intent(in) :: caseid
    character(len=*)    , intent(in) :: kelm_hist_fincl
    integer             , intent(in) :: hist_year, month, day
    integer             , intent(in) :: n_local_cells
    integer             , intent(in) :: num_cells_global_in
    integer             , intent(in) :: natural_id_cells_owned(:)  ! (n_local_cells,), full-grid indices
    real(r8)            , intent(in) :: lonc_g(:), latc_g(:), areac_g(:)  ! (ni*nj,)
    integer             , intent(in) :: nhtfrq        ! 0 monthly, < 0 every -nhtfrq hours
    integer             , intent(in) :: mfilt         ! records per file
    logical             , intent(in) :: empty_htapes  ! only kelm_hist_fincl's fields

    integer :: ierr, i, is_avg, l, c, p, g
    real(r8) :: re_km
    character(len=16) :: levname
    character(len=*), parameter :: subname = '(kelm_hist_init) '

    if (hist_init_done) call shr_sys_abort(subname//'ERROR: kelm_hist_init called twice')

    hist_caseid = caseid
    n_local     = n_local_cells
    num_global  = num_cells_global_in

    if (nhtfrq > 0) call shr_sys_abort(subname//'ERROR: kelm_hist_nhtfrq > 0 (a number '// &
         'of time steps) is not supported; use 0 (monthly) or -N (every N hours)')
    if (nhtfrq < 0) then
       if (mod(24, -nhtfrq) /= 0) call shr_sys_abort(subname//'ERROR: kelm_hist_nhtfrq = -N '// &
            'needs N to divide 24, so every interval ends at the same times each day')
    end if
    if (mfilt < 1) call shr_sys_abort(subname//'ERROR: kelm_hist_mfilt must be >= 1')
    hist_nhtfrq = nhtfrq
    hist_mfilt  = mfilt

    if (.not. empty_htapes) then
       call kELMHistoryActivateDefaults(kelm_state, ierr)
       if (ierr /= KELM_SUCCESS) &
            call shr_sys_abort(subname//'ERROR: kELMHistoryActivateDefaults failed')
    end if

    call kelm_hist_fincl_parse(kelm_state, kelm_hist_fincl, logunit)

    nfields = kELMHistoryActiveCount(kelm_state, ierr)
    if (ierr /= KELM_SUCCESS) &
         call shr_sys_abort(subname//'ERROR: kELMHistoryActiveCount failed')

    if (nfields == 0) then
       if (masterproc) write(logunit,*) subname,'no history fields active; history output disabled'
       return
    end if

    allocate(field_name(nfields), field_units(nfields), field_long_name(nfields), &
             field_is_avg(nfields), field_nlev(nfields), field_levdim(nfields), &
             hist_buf2(nfields))
    do i = 1, nfields
       field_name(i) = ' '; field_units(i) = ' '; field_long_name(i) = ' '
       call kELMHistoryActiveName(kelm_state, i-1, field_name(i), ierr)
       if (ierr /= KELM_SUCCESS) call shr_sys_abort(subname//'ERROR: kELMHistoryActiveName failed')
       call kELMHistoryActiveUnits(kelm_state, i-1, field_units(i), ierr)
       if (ierr /= KELM_SUCCESS) call shr_sys_abort(subname//'ERROR: kELMHistoryActiveUnits failed')
       call kELMHistoryActiveLongName(kelm_state, i-1, field_long_name(i), ierr)
       if (ierr /= KELM_SUCCESS) call shr_sys_abort(subname//'ERROR: kELMHistoryActiveLongName failed')
       is_avg = kELMHistoryActiveIsAvg(kelm_state, i-1, ierr)
       if (ierr /= KELM_SUCCESS) call shr_sys_abort(subname//'ERROR: kELMHistoryActiveIsAvg failed')
       field_is_avg(i) = (is_avg /= 0)
       field_nlev(i) = kELMHistoryActiveNlev(kelm_state, i-1, ierr)
       if (ierr /= KELM_SUCCESS) call shr_sys_abort(subname//'ERROR: kELMHistoryActiveNlev failed')
       field_levdim(i) = 0
       if (field_nlev(i) > 0) then
          levname = ' '
          call kELMHistoryActiveLevDim(kelm_state, i-1, levname, ierr)
          if (ierr /= KELM_SUCCESS .or. len_trim(levname) == 0) &
               call shr_sys_abort(subname//'ERROR: kELMHistoryActiveLevDim failed for '// &
               trim(field_name(i)))
          field_levdim(i) = kelm_hist_levdim(levname, field_nlev(i))
          allocate(hist_buf2(i)%v(n_local, field_nlev(i)))
          hist_buf2(i)%v = 0.0_r8
       end if
    end do

    ! PCT_LANDUNIT / PCT_NAT_PFT, from the subgrid kELM built (ELM writes
    ! them every month; they are constant without dynamic land use). Local
    ! cell g of the subgrid is history cell g.
    ilev_ltype  = kelm_hist_levdim('ltype', MAX_LUNIT)
    ilev_natpft = kelm_hist_levdim('natpft', natpft)
    allocate(pct_landunit(n_local, MAX_LUNIT), pct_nat_pft(n_local, natpft))
    pct_landunit = 0.0_r8
    pct_nat_pft  = 0.0_r8
    do l = 1, num_landunits
       g = lun_gridcell(l)
       if (g < 1 .or. g > n_local) cycle
       if (lun_itype(l) >= 1 .and. lun_itype(l) <= MAX_LUNIT) &
            pct_landunit(g, lun_itype(l)) = pct_landunit(g, lun_itype(l)) + 100.0_r8*lun_wtgcell(l)
    end do
    do p = 1, num_patches
       c = patch_column(p)
       l = col_landunit(c)
       if (lun_itype(l) /= istsoil) cycle
       g = lun_gridcell(l)
       if (g < 1 .or. g > n_local) cycle
       if (patch_itype(p) >= 0 .and. patch_itype(p) < natpft) &
            pct_nat_pft(g, patch_itype(p)+1) = pct_nat_pft(g, patch_itype(p)+1) + &
            100.0_r8*patch_wtcol(p)
    end do

    allocate(hist_buf(n_local, nfields))
    hist_buf = 0.0_r8

    ! Natural-vegetation-landunit-weight-on-gridcell == 1 guard: PATCH_TO_COL
    ! fields are natural-column values, which coincide with ELM's gridcell h0
    ! values only when the natural-veg landunit covers the whole cell (true
    ! on 1x1_brazil/1x1_glc, not in general -- see the plan's "This is
    ! natural-landunit output, not gridcell output" section).
    call kelm_hist_check_natveg_weight(logunit)

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

    ! Each (lndgrid, level) decomposition: local element (i, l), stored at
    ! (l-1)*n_local + i, is global (hist_dof(i), l).
    do g = 1, nlevdim
       block
         integer, allocatable :: dof2(:)
         integer :: nl
         nl = levdim_n(g)
         allocate(dof2(n_local*nl))
         do l = 1, nl
            do i = 1, n_local
               dof2((l-1)*n_local + i) = (l-1)*num_global + hist_dof(i)
            end do
         end do
         call PIO_initdecomp(pio_subsystem, PIO_real, (/num_global, nl/), dof2, levdim_iodesc(g))
         deallocate(dof2)
       end block
    end do

    hist_start_year  = hist_year
    hist_start_month = month
    hist_start_day   = day
    hist_start_tod   = 0
    hist_nacs        = 0

    hist_init_done = .true.

    if (masterproc) then
       write(logunit,*) subname,'kelm_hist_nhtfrq = ',hist_nhtfrq,', kelm_hist_mfilt = ', &
            hist_mfilt,', kelm_hist_empty_htapes = ',empty_htapes
       write(logunit,*) subname,nfields,' history fields active:'
       do i = 1, nfields
          write(logunit,*) subname,'   ',trim(field_name(i)),' [',trim(field_units(i)),'] ', &
               trim(field_long_name(i))
       end do
       call shr_sys_flush(logunit)
    end if

  end subroutine kelm_hist_init

  !-----------------------------------------------------------------------
  subroutine kelm_hist_fincl_parse(kelm_state, spec, logunit)
    !
    ! kelm_hist_fincl names fields to activate BEYOND the C++-side
    ! defaults -- additions/overrides, not a full restatement. Same
    ! comma-separated-string shape as kelm_kernels_parse
    ! (kelmKernelMod.F90), reused deliberately.
    !
    implicit none
    type(kELMType) , intent(in) :: kelm_state
    character(len=*), intent(in) :: spec
    integer          , intent(in) :: logunit

    integer :: ib, ie, n, ierr
    character(len=32) :: token
    character(len=*), parameter :: subname = '(kelm_hist_fincl_parse) '

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
          call kELMHistoryActivate(kelm_state, trim(token), ierr)
          if (ierr /= KELM_SUCCESS) then
             write(logunit,*) subname,'ERROR: could not activate "',trim(token), &
                  '" from kelm_hist_fincl -- unknown field name or already active'
             call shr_sys_flush(logunit)
             call shr_sys_abort(subname//'ERROR: kelm_hist_fincl names an unknown '// &
                  'or duplicate field: "'//trim(token)//'"')
          end if
       end if

       ib = ie + 2
    end do

  end subroutine kelm_hist_fincl_parse

  !-----------------------------------------------------------------------
  subroutine kelm_hist_check_natveg_weight(logunit)
    !
    ! See kelm_hist_init's comment at the call site. Since Stage 6.5 a cell's
    ! h0 value is the landunit-weighted mean over its natural and LAKE
    ! landunits (kELMHistorySetLandunitWeights, pushed by kelmLakeMod), and
    ! since Stage 6.7 its URBAN ones too (kELMHistorySetUrbanWeights, pushed
    ! by kelmUrbanMod), and its GLACIER one (kELMHistorySetGlacierWeights,
    ! kelmGlacierMod), so the guard is that those cover the cell: weight on
    ! a landunit kELM does not model (wetland, crop) would make the mean a
    ! mislabelled partial one.
    !
    implicit none
    integer, intent(in) :: logunit
    integer :: l, g, ncell
    real(r8), allocatable :: wsum(:)
    real(r8), parameter :: tol = 1.0e-6_r8
    character(len=*), parameter :: subname = '(kelm_hist_check_natveg_weight) '

    ncell = maxval(lun_gridcell(1:num_landunits))
    allocate(wsum(ncell))
    wsum = 0.0_r8
    do l = 1, num_landunits
       if (lun_itype(l) == istsoil .or. lun_itype(l) == istdlak .or. &
           lun_itype(l) == istice .or. &
           (lun_itype(l) >= isturb_tbd .and. lun_itype(l) <= isturb_md)) then
          wsum(lun_gridcell(l)) = wsum(lun_gridcell(l)) + lun_wtgcell(l)
       end if
    end do
    do g = 1, ncell
       if (abs(wsum(g) - 1.0_r8) > tol) then
          write(logunit,*) subname,'ERROR: natural + lake + urban + glacier landunit weight on local cell ', &
               g,' = ',wsum(g),' != 1.0 -- h0 is the mean over those landunits only'
          call shr_sys_flush(logunit)
          call shr_sys_abort(subname//'ERROR: natural + lake + urban + glacier landunit weight-on-gridcell '// &
               '/= 1.0; history output would silently mislabel a partial mean '// &
               'as a gridcell value')
       end if
    end do
    deallocate(wsum)

  end subroutine kelm_hist_check_natveg_weight

  !-----------------------------------------------------------------------
  subroutine kelm_hist_step(kelm_state, logunit)
    !
    ! One call per internal step. Pure device work on the C++ side -- see
    ! the plan for why this must run after phenology/SurfaceAlbedo, not
    ! before.
    !
    implicit none
    type(kELMType), intent(in) :: kelm_state
    integer         , intent(in) :: logunit
    integer :: ierr
    character(len=*), parameter :: subname = '(kelm_hist_step) '

    if (.not. hist_init_done) return

    call kELMHistoryAccumulate(kelm_state, ierr)
    if (ierr /= KELM_SUCCESS) call shr_sys_abort(subname//'ERROR: kELMHistoryAccumulate failed')
    hist_nacs = hist_nacs + 1

  end subroutine kelm_hist_step

  !-----------------------------------------------------------------------
  subroutine kelm_hist_write_if_due(kelm_state, year, month, day, tod, logunit)
    !
    ! (year, month, day, tod) is the END time of the step just taken. An
    ! interval closes at a month boundary (day==1 .and. tod==0) for
    ! hist_nhtfrq == 0, or every -hist_nhtfrq hours of the day for
    ! hist_nhtfrq < 0 (ELM's hist_nhtfrq convention). The closed interval's
    ! mean is appended to the open file as its next record; the file is
    ! closed once it holds hist_mfilt records.
    !
    implicit none
    type(kELMType), intent(in) :: kelm_state
    integer         , intent(in) :: year, month, day, tod
    integer         , intent(in) :: logunit

    integer :: ierr
    logical :: due
    character(len=*), parameter :: subname = '(kelm_hist_write_if_due) '

    if (.not. hist_init_done) return
    if (hist_nhtfrq == 0) then
       due = (day == 1 .and. tod == 0)
    else
       due = (mod(tod, -hist_nhtfrq*3600) == 0)
    end if
    if (.not. due) return

    if (hist_nacs == 0) then
       write(logunit,*) subname,'WARNING: history interval closed with no accumulated ', &
            'samples; skipping the record starting ',hist_start_year,'-',hist_start_month, &
            '-',hist_start_day
       call shr_sys_flush(logunit)
       hist_start_year = year; hist_start_month = month; hist_start_day = day
       hist_start_tod = tod
       return
    end if

    call kELMHistoryNormalize(kelm_state, ierr)
    if (ierr /= KELM_SUCCESS) call shr_sys_abort(subname//'ERROR: kELMHistoryNormalize failed')

    call kelm_hist_pull(kelm_state)

    if (.not. hist_file_open) call kelm_hist_open_file(logunit)
    hist_nrec = hist_nrec + 1
    call kelm_hist_write_record(year, month, day, tod, logunit)

    call kELMHistoryReset(kelm_state, ierr)
    if (ierr /= KELM_SUCCESS) call shr_sys_abort(subname//'ERROR: kELMHistoryReset failed')
    hist_nacs = 0

    if (hist_nrec == hist_mfilt) call kelm_hist_close_file(logunit)

    ! The interval that just closed started at hist_start_*; the next one
    ! starts exactly here.
    hist_start_year = year; hist_start_month = month; hist_start_day = day
    hist_start_tod = tod

  end subroutine kelm_hist_write_if_due

  !-----------------------------------------------------------------------
  subroutine kelm_hist_pull(kelm_state)
    !
    ! Copy every active field's normalized accumulator into the staging
    ! buffers: hist_buf for single-level fields, hist_buf2 for multi-level.
    !
    implicit none
    type(kELMType), intent(in) :: kelm_state
    integer :: ierr, ifld
    character(len=*), parameter :: subname = '(kelm_hist_pull) '

    do ifld = 1, nfields
       if (field_nlev(ifld) > 0) then
          call kELMHistoryGetField(kelm_state, trim(field_name(ifld)), hist_buf2(ifld)%v, &
               n_local*field_nlev(ifld), ierr)
       else
          call kELMHistoryGetField(kelm_state, trim(field_name(ifld)), hist_buf(:,ifld), &
               n_local, ierr)
       end if
       if (ierr /= KELM_SUCCESS) &
            call shr_sys_abort(subname//'ERROR: kELMHistoryGetField failed for '// &
            trim(field_name(ifld)))
    end do

  end subroutine kelm_hist_pull

  !-----------------------------------------------------------------------
  subroutine kelm_hist_open_file(logunit)
    !
    ! Create the next history file and write what does not change between
    ! its records: dimensions, variable definitions, lon/lat/area and the
    ! subgrid fractions. Named for its first record's interval start:
    ! YYYY-MM for monthly files, ELM's YYYY-MM-DD-SSSSS otherwise.
    !
    ! PCT_LANDUNIT/PCT_NAT_PFT are time-constant. A one-record file keeps
    ! ELM's h0 layout (lndgrid, lev, time); a multi-record file stores them
    ! once, without the time dimension, rather than repeating them in every
    ! record (on a 16M-cell grid that repeat alone is ~1.7 GB per record).
    !
    implicit none
    integer, intent(in) :: logunit

    type(var_desc_t)  :: vid_lon, vid_lat, vid_area, vid_pctl, vid_pctp
    integer :: dimid_lndgrid, dimid_time, dimid_histint, dimid_lev(MAX_LEVDIM)
    integer :: ierr, ifld, k
    logical :: pct_has_time
    character(len=16)  :: cdate
    real(r4), allocatable :: wbuf(:), wbuf2(:)
    character(len=*), parameter :: subname = '(kelm_hist_open_file) '

    if (hist_nhtfrq == 0) then
       write(cdate,'(i4.4,"-",i2.2)') hist_start_year, hist_start_month
    else
       write(cdate,'(i4.4,"-",i2.2,"-",i2.2,"-",i5.5)') hist_start_year, hist_start_month, &
            hist_start_day, hist_start_tod
    end if
    hist_fname = trim(hist_caseid)//'.kelm.h0.'//trim(cdate)//'.nc'
    pct_has_time = (hist_mfilt == 1)

    ierr = PIO_createfile(pio_subsystem, hist_file, io_type, trim(hist_fname), &
         ior(PIO_CLOBBER, io_format))
    if (ierr /= PIO_NOERR) call shr_sys_abort(subname//'ERROR: PIO_createfile failed for '// &
         trim(hist_fname))

    ierr = PIO_def_dim(hist_file, 'lndgrid', num_global, dimid_lndgrid)
    if (ierr /= PIO_NOERR) call shr_sys_abort(subname//'ERROR: PIO_def_dim lndgrid failed')
    ierr = PIO_def_dim(hist_file, 'time', PIO_UNLIMITED, dimid_time)
    if (ierr /= PIO_NOERR) call shr_sys_abort(subname//'ERROR: PIO_def_dim time failed')
    ierr = PIO_def_dim(hist_file, 'hist_interval', 2, dimid_histint)
    if (ierr /= PIO_NOERR) call shr_sys_abort(subname//'ERROR: PIO_def_dim hist_interval failed')
    do k = 1, nlevdim
       ierr = PIO_def_dim(hist_file, trim(levdim_name(k)), levdim_n(k), dimid_lev(k))
       if (ierr /= PIO_NOERR) call shr_sys_abort(subname//'ERROR: PIO_def_dim '// &
            trim(levdim_name(k))//' failed')
    end do

    ! Time-constant per-cell fields.
    ierr = PIO_def_var(hist_file, 'lon', PIO_real, (/dimid_lndgrid/), vid_lon)
    if (ierr /= PIO_NOERR) call shr_sys_abort(subname//'ERROR: PIO_def_var lon failed')
    ierr = PIO_put_att(hist_file, vid_lon, 'units', 'degrees_east')
    ierr = PIO_def_var(hist_file, 'lat', PIO_real, (/dimid_lndgrid/), vid_lat)
    if (ierr /= PIO_NOERR) call shr_sys_abort(subname//'ERROR: PIO_def_var lat failed')
    ierr = PIO_put_att(hist_file, vid_lat, 'units', 'degrees_north')
    ierr = PIO_def_var(hist_file, 'area', PIO_real, (/dimid_lndgrid/), vid_area)
    if (ierr /= PIO_NOERR) call shr_sys_abort(subname//'ERROR: PIO_def_var area failed')
    ierr = PIO_put_att(hist_file, vid_area, 'units', 'km^2')

    ! Time / time_bounds. NOT part of the lndgrid decomposition -- one value
    ! per record, written collectively in kelm_hist_write_record.
    ierr = PIO_def_var(hist_file, 'time', PIO_double, (/dimid_time/), hv_time)
    if (ierr /= PIO_NOERR) call shr_sys_abort(subname//'ERROR: PIO_def_var time failed')
    ierr = PIO_put_att(hist_file, hv_time, 'units', 'days since 0001-01-01 00:00:00')
    ierr = PIO_put_att(hist_file, hv_time, 'calendar', trim(shr_cal_noleap))
    ierr = PIO_put_att(hist_file, hv_time, 'bounds', 'time_bounds')
    ierr = PIO_def_var(hist_file, 'time_bounds', PIO_double, (/dimid_histint, dimid_time/), &
         hv_time_bounds)
    if (ierr /= PIO_NOERR) call shr_sys_abort(subname//'ERROR: PIO_def_var time_bounds failed')

    ! ELM's per-record date metadata, at the record's end time.
    ierr = PIO_def_var(hist_file, 'mcdate', PIO_int, (/dimid_time/), hv_mcdate)
    ierr = PIO_put_att(hist_file, hv_mcdate, 'long_name', 'current date (YYYYMMDD)')
    ierr = PIO_def_var(hist_file, 'mcsec', PIO_int, (/dimid_time/), hv_mcsec)
    ierr = PIO_put_att(hist_file, hv_mcsec, 'long_name', 'current seconds of current date')
    ierr = PIO_put_att(hist_file, hv_mcsec, 'units', 's')
    ierr = PIO_def_var(hist_file, 'mdcur', PIO_int, (/dimid_time/), hv_mdcur)
    ierr = PIO_put_att(hist_file, hv_mdcur, 'long_name', 'current day (from base day)')
    ierr = PIO_def_var(hist_file, 'mscur', PIO_int, (/dimid_time/), hv_mscur)
    ierr = PIO_put_att(hist_file, hv_mscur, 'long_name', 'current seconds of current day')

    ! Subgrid fractions.
    if (pct_has_time) then
       ierr = PIO_def_var(hist_file, 'PCT_LANDUNIT', PIO_real, &
            (/dimid_lndgrid, dimid_lev(ilev_ltype), dimid_time/), vid_pctl)
    else
       ierr = PIO_def_var(hist_file, 'PCT_LANDUNIT', PIO_real, &
            (/dimid_lndgrid, dimid_lev(ilev_ltype)/), vid_pctl)
    end if
    if (ierr /= PIO_NOERR) call shr_sys_abort(subname//'ERROR: PIO_def_var PCT_LANDUNIT failed')
    ierr = PIO_put_att(hist_file, vid_pctl, 'units', '%')
    ierr = PIO_put_att(hist_file, vid_pctl, 'long_name', '% of each landunit on topounit')
    ierr = PIO_put_att(hist_file, vid_pctl, '_FillValue', HIST_FILL)
    ierr = PIO_put_att(hist_file, vid_pctl, 'missing_value', HIST_FILL)
    if (pct_has_time) then
       ierr = PIO_def_var(hist_file, 'PCT_NAT_PFT', PIO_real, &
            (/dimid_lndgrid, dimid_lev(ilev_natpft), dimid_time/), vid_pctp)
    else
       ierr = PIO_def_var(hist_file, 'PCT_NAT_PFT', PIO_real, &
            (/dimid_lndgrid, dimid_lev(ilev_natpft)/), vid_pctp)
    end if
    if (ierr /= PIO_NOERR) call shr_sys_abort(subname//'ERROR: PIO_def_var PCT_NAT_PFT failed')
    ierr = PIO_put_att(hist_file, vid_pctp, 'units', '%')
    ierr = PIO_put_att(hist_file, vid_pctp, 'long_name', &
         '% of each PFT on the natural vegetation (i.e., soil) landunit')
    ierr = PIO_put_att(hist_file, vid_pctp, '_FillValue', HIST_FILL)
    ierr = PIO_put_att(hist_file, vid_pctp, 'missing_value', HIST_FILL)

    ! Accumulating fields.
    if (allocated(hv_field)) deallocate(hv_field)
    allocate(hv_field(nfields))
    do ifld = 1, nfields
       if (field_nlev(ifld) > 0) then
          ierr = PIO_def_var(hist_file, trim(field_name(ifld)), PIO_real, &
               (/dimid_lndgrid, dimid_lev(field_levdim(ifld)), dimid_time/), hv_field(ifld))
       else
          ierr = PIO_def_var(hist_file, trim(field_name(ifld)), PIO_real, &
               (/dimid_lndgrid, dimid_time/), hv_field(ifld))
       end if
       if (ierr /= PIO_NOERR) &
            call shr_sys_abort(subname//'ERROR: PIO_def_var failed for '//trim(field_name(ifld)))
       ierr = PIO_put_att(hist_file, hv_field(ifld), 'units', trim(field_units(ifld)))
       ierr = PIO_put_att(hist_file, hv_field(ifld), 'long_name', trim(field_long_name(ifld)))
       ierr = PIO_put_att(hist_file, hv_field(ifld), '_FillValue', HIST_FILL)
       ierr = PIO_put_att(hist_file, hv_field(ifld), 'missing_value', HIST_FILL)
       if (field_is_avg(ifld)) then
          ierr = PIO_put_att(hist_file, hv_field(ifld), 'cell_methods', 'time: mean')
       end if
    end do

    ierr = PIO_enddef(hist_file)
    if (ierr /= PIO_NOERR) call shr_sys_abort(subname//'ERROR: PIO_enddef failed for '// &
         trim(hist_fname))

    ! Time-constant, decomposed fields -- once per file (each file is
    ! self-contained), not once for the whole run.
    allocate(wbuf(n_local))
    wbuf = real(hist_lon, r4)
    call PIO_write_darray(hist_file, vid_lon, hist_iodesc, wbuf, ierr)
    if (ierr /= PIO_NOERR) call shr_sys_abort(subname//'ERROR: PIO_write_darray lon failed')
    wbuf = real(hist_lat, r4)
    call PIO_write_darray(hist_file, vid_lat, hist_iodesc, wbuf, ierr)
    if (ierr /= PIO_NOERR) call shr_sys_abort(subname//'ERROR: PIO_write_darray lat failed')
    wbuf = real(hist_area, r4)
    call PIO_write_darray(hist_file, vid_area, hist_iodesc, wbuf, ierr)
    if (ierr /= PIO_NOERR) call shr_sys_abort(subname//'ERROR: PIO_write_darray area failed')
    deallocate(wbuf)

    if (pct_has_time) call PIO_setframe(hist_file, vid_pctl, 1_PIO_OFFSET_KIND, ierr)
    allocate(wbuf2(n_local*MAX_LUNIT))
    wbuf2 = real(reshape(pct_landunit, (/n_local*MAX_LUNIT/)), r4)
    call PIO_write_darray(hist_file, vid_pctl, levdim_iodesc(ilev_ltype), wbuf2, ierr)
    if (ierr /= PIO_NOERR) call shr_sys_abort(subname//'ERROR: PIO_write_darray PCT_LANDUNIT failed')
    deallocate(wbuf2)
    if (pct_has_time) call PIO_setframe(hist_file, vid_pctp, 1_PIO_OFFSET_KIND, ierr)
    allocate(wbuf2(n_local*natpft))
    wbuf2 = real(reshape(pct_nat_pft, (/n_local*natpft/)), r4)
    call PIO_write_darray(hist_file, vid_pctp, levdim_iodesc(ilev_natpft), wbuf2, ierr)
    if (ierr /= PIO_NOERR) call shr_sys_abort(subname//'ERROR: PIO_write_darray PCT_NAT_PFT failed')
    deallocate(wbuf2)

    hist_file_open = .true.
    hist_nrec = 0

  end subroutine kelm_hist_open_file

  !-----------------------------------------------------------------------
  subroutine kelm_hist_write_record(end_year, end_month, end_day, end_tod, logunit)
    !
    ! Append record hist_nrec (already incremented) to the open file: its
    ! time metadata and every active field's interval mean.
    !
    implicit none
    integer, intent(in) :: end_year, end_month, end_day, end_tod
    integer, intent(in) :: logunit

    integer :: ierr, ifld, ival(1)
    real(r8) :: start_days, end_days, tb(2,1), tval(1)
    real(r4), allocatable :: wbuf(:), wbuf2(:)
    character(len=*), parameter :: subname = '(kelm_hist_write_record) '

    ! Time / time_bounds: absolute days since the twin case's own epoch
    ! (0001-01-01 00:00:00, noleap) -- not run-relative, so a kELM h0 file
    ! overlays its ELM twin without a conversion step.
    start_days = kelm_hist_days_since_epoch(hist_start_year, hist_start_month, &
         hist_start_day, hist_start_tod)
    end_days   = kelm_hist_days_since_epoch(end_year, end_month, end_day, end_tod)
    tval(1)   = end_days
    tb(1,1)   = start_days
    tb(2,1)   = end_days
    ierr = PIO_put_var(hist_file, hv_time, (/hist_nrec/), tval(1))
    if (ierr /= PIO_NOERR) call shr_sys_abort(subname//'ERROR: PIO_put_var time failed')
    ierr = PIO_put_var(hist_file, hv_time_bounds, (/1,hist_nrec/), (/2,1/), tb)
    if (ierr /= PIO_NOERR) call shr_sys_abort(subname//'ERROR: PIO_put_var time_bounds failed')
    ival(1) = end_year*10000 + end_month*100 + end_day
    ierr = PIO_put_var(hist_file, hv_mcdate, (/hist_nrec/), ival(1))
    ival(1) = end_tod
    ierr = PIO_put_var(hist_file, hv_mcsec, (/hist_nrec/), ival(1))
    ! Whole days and seconds since the base date, as ELM's get_curr_time.
    ival(1) = int(end_days - hist_base_days + 1.0e-9_r8)
    ierr = PIO_put_var(hist_file, hv_mdcur, (/hist_nrec/), ival(1))
    ival(1) = nint((end_days - hist_base_days - real(int(end_days - hist_base_days + &
         1.0e-9_r8), r8)) * 86400.0_r8)
    ierr = PIO_put_var(hist_file, hv_mscur, (/hist_nrec/), ival(1))
    if (ierr /= PIO_NOERR) call shr_sys_abort(subname//'ERROR: PIO_put_var date metadata failed')

    allocate(wbuf(n_local))
    do ifld = 1, nfields
       call PIO_setframe(hist_file, hv_field(ifld), int(hist_nrec, PIO_OFFSET_KIND), ierr)
       if (ierr /= PIO_NOERR) &
            call shr_sys_abort(subname//'ERROR: PIO_setframe failed for '//trim(field_name(ifld)))
       if (field_nlev(ifld) > 0) then
          allocate(wbuf2(n_local*field_nlev(ifld)))
          wbuf2 = real(reshape(hist_buf2(ifld)%v, (/n_local*field_nlev(ifld)/)), r4)
          call PIO_write_darray(hist_file, hv_field(ifld), levdim_iodesc(field_levdim(ifld)), &
               wbuf2, ierr)
          deallocate(wbuf2)
       else
          wbuf = real(hist_buf(:,ifld), r4)
          call PIO_write_darray(hist_file, hv_field(ifld), hist_iodesc, wbuf, ierr)
       end if
       if (ierr /= PIO_NOERR) &
            call shr_sys_abort(subname//'ERROR: PIO_write_darray failed for '//trim(field_name(ifld)))
    end do
    deallocate(wbuf)

    ! Each record is complete on disk before the next interval starts, so a
    ! run that dies mid-file leaves the records before it readable.
    call PIO_syncfile(hist_file)

    if (masterproc) then
       write(logunit,*) subname,'wrote record ',hist_nrec,' of ',trim(hist_fname),' (', &
            nfields,' fields, ',hist_nacs,' samples)'
       call shr_sys_flush(logunit)
    end if

  end subroutine kelm_hist_write_record

  !-----------------------------------------------------------------------
  subroutine kelm_hist_close_file(logunit)
    implicit none
    integer, intent(in) :: logunit
    if (.not. hist_file_open) return
    call PIO_closefile(hist_file)
    hist_file_open = .false.
    if (masterproc) then
       write(logunit,*) '(kelm_hist_close_file) closed ',trim(hist_fname),' with ', &
            hist_nrec,' records'
       call shr_sys_flush(logunit)
    end if
    hist_nrec = 0
  end subroutine kelm_hist_close_file

  !-----------------------------------------------------------------------
  function kelm_hist_days_since_epoch(year, month, day, tod) result(days)
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

  end function kelm_hist_days_since_epoch

  !-----------------------------------------------------------------------
  subroutine kelm_hist_final(kelm_state, year, month, day, tod, logunit)
    !
    ! Release history's resources. A partial trailing interval (nacs > 0) is
    ! dropped, not written: ELM writes a record only when its interval
    ! closes, never at the end of a run that stops mid-interval. A file
    ! holding fewer than hist_mfilt complete records is closed as it is.
    ! Must run BEFORE kELMDestroy -- kelm_final calls this first, ahead of
    ! tearing the object down.
    !
    implicit none
    type(kELMType), intent(in) :: kelm_state
    integer         , intent(in) :: year, month, day, tod
    integer         , intent(in) :: logunit
    integer :: k

    if (.not. hist_init_done) return

    if (hist_nacs > 0 .and. masterproc) then
       write(logunit,*) '(kelm_hist_final) partial interval not written (as ELM): ', &
            hist_nacs,' samples since ',hist_start_year,'-',hist_start_month,'-',hist_start_day
       call shr_sys_flush(logunit)
    end if
    hist_nacs = 0

    call kelm_hist_close_file(logunit)

    call PIO_freedecomp(pio_subsystem, hist_iodesc)
    do k = 1, nlevdim
       call PIO_freedecomp(pio_subsystem, levdim_iodesc(k))
    end do
    hist_init_done = .false.

  end subroutine kelm_hist_final

end module kelmHistMod
