module lnd_comp_mct

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! MCT coupling layer for ELMxx.
  !
  ! Structured after components/rdycore/src/cpl/rof_comp_mct.F90. ELMxx presents
  ! itself to the coupler as a present, prognostic land model: it registers a
  ! real decomposition and a real domain, but its time stepping is a no-op and
  ! it exports zeros until the Kokkos port is wired in.
  !-----------------------------------------------------------------------

  use esmf
  use mct_mod
  use seq_flds_mod
  use seq_cdata_mod   , only : seq_cdata, seq_cdata_setptrs
  use seq_infodata_mod, only : seq_infodata_type, seq_infodata_PutData, seq_infodata_GetData
  use seq_comm_mct    , only : seq_comm_inst, seq_comm_name, seq_comm_suffix
  use seq_timemgr_mod , only : seq_timemgr_EClockDateInSync
  use shr_kind_mod    , only : IN=>SHR_KIND_IN, R8=>SHR_KIND_R8, CS=>SHR_KIND_CS, CL=>SHR_KIND_CL
  use shr_sys_mod     , only : shr_sys_abort, shr_sys_flush
  use shr_file_mod    , only : shr_file_getunit, shr_file_getlogunit, shr_file_getloglevel
  use shr_file_mod    , only : shr_file_setlogunit, shr_file_setloglevel, shr_file_setio
  use shr_file_mod    , only : shr_file_freeunit
  use elmxxSpmdMod    , only : masterproc, mpicom_lnd, iam, npes, LNDID, elmxxSpmdInit
  use elmxxMod        , only : elmxx_read_namelist, elmxx_init, elmxx_run, elmxx_final
  use elmxxMod        , only : elmxx_init_albedo
  use shr_orb_mod     , only : shr_orb_decl, SHR_ORB_UNDEF_REAL
  use elmxxMod        , only : num_cells_owned, num_cells_global, natural_id_cells_owned
  use elmxxMod        , only : elmxx_caseid
  use elmxxMod        , only : nlon_g, nlat_g, lonc_g, latc_g, areac_g, maskc_g, fracc_g
  use elmxxMod        , only : inst_name, inst_index, inst_suffix, do_elmxx
  use elmxx_cpl_indices, only : elmxx_cpl_indices_set
  use elmxxForcingMod  , only : elmxx_import

  !
  ! !PUBLIC TYPES:
  implicit none
  save
  private ! except

  !--------------------------------------------------------------------------
  ! Public interfaces
  !--------------------------------------------------------------------------

  public :: lnd_init_mct
  public :: lnd_run_mct

  ! ELMxx's own model clock. ELM keeps one (elm_time_manager) and its
  ! lnd_run_mct loops until that clock catches up with the coupler's, which is
  ! why the first coupling call runs the driver twice. ELMxx had no clock and
  ! ran exactly once per call, leaving it a physics pass behind.
  type(ESMF_Time), private :: elmxx_clock_time
  logical, private :: elmxx_clock_started = .false.
  integer, private :: elmxx_nstep = 0
  public :: lnd_final_mct

  !--------------------------------------------------------------------------
  ! Private interfaces
  !--------------------------------------------------------------------------

  private :: lnd_SetgsMap_mct
  private :: lnd_domain_mct
  private :: get_step_size, get_clock_date
  private :: chkrc

  !--------------------------------------------------------------------------
  ! Private module data
  !--------------------------------------------------------------------------

  integer :: lsize                          ! number of cells owned by this rank
  integer :: logunit_lnd = 6                ! "stdout" log file unit number

  character(*), parameter :: F00 = "('(lnd_comp_mct) ',8a)"

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CONTAINS
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  !===============================================================================
  ! !IROUTINE: lnd_init_mct
  !
  ! !DESCRIPTION:
  !     Initialize ELMxx and register its decomposition and domain with the coupler
  !===============================================================================

  subroutine lnd_init_mct( EClock, cdata, x2l_l, l2x_l, NLFilename )

    ! !INPUT/OUTPUT PARAMETERS:

    type(ESMF_Clock)            , intent(inout) :: EClock
    type(seq_cdata)             , intent(inout) :: cdata
    type(mct_aVect)             , intent(inout) :: x2l_l, l2x_l
    character(len=*), optional  , intent(in)    :: NLFilename

    !--- local ---
    type(seq_infodata_type), pointer :: infodata
    type(mct_gsMap)        , pointer :: gsMap_lnd
    type(mct_gGrid)        , pointer :: dom_l
    integer :: shrlogunit                     ! original log unit
    integer :: shrloglev                      ! original log level
    integer :: mpicom_loc                     ! local mpi communicator
    integer :: month, day, year
    logical :: exists                         ! true if file exists

    character(*), parameter :: subName = "(lnd_init_mct) "
    !-------------------------------------------------------------------------------

    ! Set cdata pointers
    call seq_cdata_setptrs(cdata, &
         id=LNDID, &
         mpicom=mpicom_loc, &
         gsMap=gsMap_lnd, &
         dom=dom_l, &
         infodata=infodata)

    ! Initialize ELMxx MPI communicator
    call elmxxSpmdInit(mpicom_loc)

    ! Determine instance information
    inst_name   = seq_comm_name(LNDID)
    inst_index  = seq_comm_inst(LNDID)
    inst_suffix = seq_comm_suffix(LNDID)

    !--- open log file ---
    call shr_file_getLogUnit (shrlogunit)
    if (masterproc) then
       inquire(file='lnd_modelio.nml'//trim(inst_suffix),exist=exists)
       if (exists) then
          logunit_lnd = shr_file_getUnit()
          call shr_file_setIO('lnd_modelio.nml'//trim(inst_suffix),logunit_lnd)
       end if
       write(logunit_lnd,*) "ELMxx model initialization"
    else
       logunit_lnd = shrlogunit
    end if

    call shr_file_getLogLevel(shrloglev)
    call shr_file_setLogUnit (logunit_lnd)

    if (masterproc) then
       write(logunit_lnd,*) ' ELMxx npes = ', npes
       write(logunit_lnd,*) ' ELMxx iam  = ', iam
       write(logunit_lnd,*) ' inst_name  = ', trim(inst_name)
    endif

    !----------------------------------------------------------------------------
    ! Read the namelist and initialize ELMxx
    !----------------------------------------------------------------------------

    call elmxx_read_namelist(logunit_lnd)

    if (.not. do_elmxx) then
       call seq_infodata_PutData( infodata, lnd_present=.false., lnd_prognostic=.false.)
       call shr_file_setLogUnit (shrlogunit)
       call shr_file_setLogLevel(shrloglev)
       return
    end if

    call get_clock_date(EClock, month, day, year)

    ! History output filenames are casename-derived (casename.elmxx.h0.*.nc),
    ! same as ELM's own h0 convention -- caseid was not read anywhere in
    ! ELMxx before history needed it.
    call seq_infodata_GetData(infodata, case_name=elmxx_caseid)

    call elmxx_init(logunit_lnd, year, month, day)

    !----------------------------------------------------------------------------
    ! Register the ELMxx decomposition and domain with the coupler
    !----------------------------------------------------------------------------

    call lnd_SetgsMap_mct( gsMap_lnd )

    call lnd_domain_mct( lsize, gsMap_lnd, dom_l )

    ! Initialize cpl -> ELMxx attribute vector
    call mct_aVect_init(x2l_l, rList=seq_flds_x2l_fields, lsize=lsize)
    call mct_aVect_zero(x2l_l)

    ! Initialize ELMxx -> cpl attribute vector
    call mct_aVect_init(l2x_l, rList=seq_flds_l2x_fields, lsize=lsize)
    call mct_aVect_zero(l2x_l)

    ! Resolve the coupler field indices now that both attribute vectors exist.
    call elmxx_cpl_indices_set(x2l_l, l2x_l)

    !----------------------------------------------------------------------------
    ! Fill infodata that needs to be returned from ELMxx
    !----------------------------------------------------------------------------

    call seq_infodata_PutData( infodata, lnd_present=.true., lnd_prognostic=.true., &
         lnd_nx=nlon_g, lnd_ny=nlat_g)

    !----------------------------------------------------------------------------
    ! Reset shr logging to original values
    !----------------------------------------------------------------------------

    if (masterproc) write(logunit_lnd,F00) 'lnd_init_mct done'
    ! Initial albedo pass, as ELM does in initialize2. SurfaceRadiation reads
    ! the PREVIOUS step's albedo, and doalb is false on the first two driver
    ! passes, so without this step 0 would run on cold-start constants and
    ! photosynthesis would have no vcmaxcint until step 2.
    block
      type(seq_infodata_type), pointer :: infodata_i
      real(r8) :: nextsw_cday_i, declin_i, eccen_i, obliqr_i, lambm0_i, mvelpp_i, eccf_i
      call seq_cdata_setptrs(cdata, infodata=infodata_i)
      call seq_infodata_GetData(infodata_i, nextsw_cday=nextsw_cday_i, &
           orb_eccen=eccen_i, orb_mvelpp=mvelpp_i, &
           orb_lambm0=lambm0_i, orb_obliqr=obliqr_i)
      if (nextsw_cday_i > 0._r8) then
         call shr_orb_decl(nextsw_cday_i, eccen_i, mvelpp_i, lambm0_i, obliqr_i, &
              declin_i, eccf_i)
         call elmxx_init_albedo(logunit_lnd, nextsw_cday_i, declin_i)
      end if
    end block

    call shr_sys_flush(logunit_lnd)

    call shr_file_setLogUnit (shrlogunit)
    call shr_file_setLogLevel(shrloglev)

  end subroutine lnd_init_mct

  !===============================================================================
  ! !IROUTINE: lnd_run_mct
  !
  ! !DESCRIPTION:
  !     Advance ELMxx by one coupling interval. Atmospheric forcing is imported
  !     from x2l_l; l2x_l is still left as initialized (zero) because no physics
  !     runs yet to produce anything to send back.
  !===============================================================================

  subroutine lnd_run_mct( EClock, cdata, x2l_l, l2x_l )

    implicit none

    ! !INPUT/OUTPUT PARAMETERS:

    type(ESMF_Clock)            ,intent(inout) :: EClock
    type(seq_cdata)             ,intent(inout) :: cdata
    type(mct_aVect)             ,intent(inout) :: x2l_l, l2x_l

    !--- local ---
    integer :: coupling_dt_in_sec
    integer :: month, day
    type(seq_infodata_type), pointer :: infodata
    real(r8) :: nextsw_cday      ! calendar day of the NEXT radiation step
    real(r8) :: declinp1         ! solar declination for that step, radians
    real(r8) :: eccen, obliqr, lambm0, mvelpp, eccf
    type(ESMF_TimeInterval) :: elmxx_step
    type(ESMF_Time) :: hist_time
    logical  :: dosend, doalb_step
    integer  :: rc, cyr, cmon, cday, ctod, cymd
    integer  :: hist_year, hist_month, hist_day, hist_tod
    !-------------------------------------------------------------------------------

    if (.not. do_elmxx) return

    coupling_dt_in_sec = get_step_size(EClock)
    call get_clock_date(EClock, month, day)

    ! Orbital state for SurfaceAlbedo. Taken from the coupler rather than
    ! recomputed from the model clock: nextsw_cday is the day of the
    ! ATMOSPHERE's next radiation step, and the land albedo has to be
    ! computed for that instant or the two drift apart. This is exactly what
    ! ELM's lnd_comp_mct hands to elm_drv.
    call seq_cdata_setptrs(cdata, infodata=infodata)
    call seq_infodata_GetData(infodata, nextsw_cday=nextsw_cday, &
         orb_eccen=eccen, orb_mvelpp=mvelpp, &
         orb_lambm0=lambm0, orb_obliqr=obliqr)
    call shr_orb_decl(nextsw_cday, eccen, mvelpp, lambm0, obliqr, declinp1, eccf)

    call elmxx_import(logunit_lnd, x2l_l)

    ! Advance over model steps until ELMxx's clock catches the coupler's, the
    ! way ELM's lnd_run_mct does. On the first coupling call the clock starts
    ! at the run start while EClock is already at the end of the interval, so
    ! this runs twice -- nstep 0 then nstep 1 -- and once per call thereafter.
    ! Running once unconditionally left ELMxx a physics pass behind ELM and
    ! also broke as soon as l_ncpl stopped matching the model timestep.
    if (.not. elmxx_clock_started) then
       call ESMF_ClockGet(EClock, startTime=elmxx_clock_time, rc=rc)
       call chkrc(rc, 'lnd::lnd_run_mct: ESMF_ClockGet startTime')
       elmxx_clock_started = .true.
    end if
    call ESMF_ClockGet(EClock, timeStep=elmxx_step, rc=rc)
    call chkrc(rc, 'lnd::lnd_run_mct: ESMF_ClockGet timeStep')

    dosend = .false.
    do while (.not. dosend)

       call ESMF_TimeGet(elmxx_clock_time, yy=cyr, mm=cmon, dd=cday, s=ctod, rc=rc)
       call chkrc(rc, 'lnd::lnd_run_mct: ESMF_TimeGet model clock')
       cymd = cyr*10000 + cmon*100 + cday
       dosend = seq_timemgr_EClockDateInSync(EClock, cymd, ctod)

       ! ELM's doalb, exactly: no albedo on the nstep-0 pass; at nstep 1 only
       ! if the atmosphere's next radiation day is this step's next day; after
       ! that whenever the coupler supplies a valid nextsw_cday.
       if (elmxx_nstep == 0) then
          doalb_step = .false.
       else if (elmxx_nstep == 1) then
          doalb_step = (abs(nextsw_cday - elmxx_caldayp1(elmxx_clock_time, elmxx_step)) < 1.e-10_r8)
       else
          doalb_step = (nextsw_cday >= -0.5_r8)
       end if

       ! History's month-end trigger needs the per-substep END time -- a
       ! separate clock read from elmxx_clock_time (this substep's START),
       ! not the once-per-coupling-interval month/day above. Deliberately a
       ! new argument rather than repurposing month/day: those already drive
       ! phenology, and changing what they carry as a side effect of adding
       ! history would risk a silent phenology behavior change.
       hist_time = elmxx_clock_time + elmxx_step
       call ESMF_TimeGet(hist_time, yy=hist_year, mm=hist_month, dd=hist_day, &
                         s=hist_tod, rc=rc)
       call chkrc(rc, 'lnd::lnd_run_mct: ESMF_TimeGet hist_time')

       call elmxx_run(logunit_lnd, coupling_dt_in_sec, month, day, &
                      nextsw_cday, declinp1, doalb=doalb_step, &
                      hist_year=hist_year, hist_month=hist_month, &
                      hist_day=hist_day, hist_tod=hist_tod)

       elmxx_nstep = elmxx_nstep + 1
       elmxx_clock_time = elmxx_clock_time + elmxx_step
    end do

  end subroutine lnd_run_mct

  !===============================================================================
  ! !IROUTINE: lnd_final_mct
  !
  ! !DESCRIPTION:
  !     Finalize ELMxx
  !===============================================================================

  subroutine lnd_final_mct( EClock, cdata, x2l_l, l2x_l)

    implicit none

    ! !INPUT/OUTPUT PARAMETERS:

    type(ESMF_Clock)            ,intent(inout) :: EClock
    type(seq_cdata)             ,intent(inout) :: cdata
    type(mct_aVect)             ,intent(inout) :: x2l_l, l2x_l
    !-------------------------------------------------------------------------------
    type(ESMF_Time) :: current_time
    integer :: rc, year, month, day, tod

    if (.not. do_elmxx) return

    ! EClock is received but was otherwise unused here; history's partial-
    ! month flush needs an end time for the trailing (possibly incomplete)
    ! interval, so extract it the same way get_clock_date does.
    call ESMF_ClockGet(EClock, currTime=current_time, rc=rc)
    call chkrc(rc, 'lnd::lnd_final_mct: error return from ESMF_ClockGet')
    call ESMF_TimeGet(current_time, yy=year, mm=month, dd=day, s=tod, rc=rc)
    call chkrc(rc, 'lnd::lnd_final_mct: error return from ESMF_TimeGet')

    call elmxx_final(year, month, day, tod)

    if (masterproc .and. logunit_lnd /= 6) close (logunit_lnd)

  end subroutine lnd_final_mct

  !===============================================================================
  ! !IROUTINE: lnd_SetgsMap_mct
  !
  ! !DESCRIPTION:
  !     Build the MCT global segment map from the ELMxx decomposition.
  !
  !     The segments cover only the active land cells (mask == 1), but the gsMap
  !     global size is the full ni*nj grid. That is exactly what ELM does
  !     (lnd_setgsmap_mct passes gindex = ldecomp%gdc2glo over numg land cells
  !     with gsize = ldomain%ni * ldomain%nj), and the coupler requires the
  !     global size to match the atm grid when the atm and lnd grids are the
  !     same (seq_domain_mct.F90 aborts on gatmsize /= glndsize).
  !===============================================================================

  subroutine lnd_SetgsMap_mct( gsMap_lnd )

    implicit none
    !
    type(mct_gsMap), intent(inout) :: gsMap_lnd   ! MCT gsmap for the land model
    !
    ! LOCAL VARIABLES
    integer :: i
    integer, allocatable :: gindex(:)
    character(len=32), parameter :: sub = 'lnd_SetgsMap_mct'
    !-----------------------------------------------------

    lsize = num_cells_owned

    allocate(gindex(lsize))
    do i = 1, lsize
       gindex(i) = natural_id_cells_owned(i)
    end do

    ! gsize is the full ni*nj grid, NOT num_cells_global (which counts only the
    ! active land cells). The segments cover just the land cells, but the gsMap
    ! global size must still match the atm grid -- see the note above.
    call mct_gsMap_init( gsMap_lnd, gindex, mpicom_lnd, LNDID, lsize, nlon_g*nlat_g )

    deallocate(gindex)

  end subroutine lnd_SetgsMap_mct

  !===============================================================================
  ! !IROUTINE: lnd_domain_mct
  !
  ! !DESCRIPTION:
  !     Send the land model domain information to the coupler.
  !
  !     lat/lon in degrees, area in radians^2, mask is 1 (land), 0 (non-land).
  !     In addition land carries around landfrac for the purposes of domain
  !     checking. aream is deliberately left at its initialized special value --
  !     it is filled in by the atm-lnd mapper, exactly as in ELM.
  !===============================================================================

  subroutine lnd_domain_mct( lsz, gsMap_lnd, dom_lnd )

    implicit none
    !
    integer        , intent(in)    :: lsz
    type(mct_gsMap), intent(in)    :: gsMap_lnd
    type(mct_gGrid), intent(inout) :: dom_lnd
    !
    ! LOCAL VARIABLES
    integer :: n, ni
    integer , pointer :: idata(:) ! temporary
    real(r8), pointer :: data(:)  ! temporary
    character(len=32), parameter :: sub = 'lnd_domain_mct'
    !-----------------------------------------------------

    call mct_gGrid_init( GGrid=dom_lnd, CoordChars=trim(seq_flds_dom_coord), &
         OtherChars=trim(seq_flds_dom_other), lsize=lsz )

    ! Allocate memory
    allocate(data(lsz))

    ! Determine global gridpoint number attribute, GlobGridNum, which is set automatically by MCT
    call mct_gsMap_orderedPoints(gsMap_lnd, iam, idata)
    call mct_gGrid_importIAttr(dom_lnd,'GlobGridNum',idata,lsz)

    ! Initialize attribute vector with special value
    data(:) = -9999.0_R8
    call mct_gGrid_importRAttr(dom_lnd,"lat"  ,data,lsz)
    call mct_gGrid_importRAttr(dom_lnd,"lon"  ,data,lsz)
    call mct_gGrid_importRAttr(dom_lnd,"area" ,data,lsz)
    call mct_gGrid_importRAttr(dom_lnd,"aream",data,lsz)
    data(:) = 0.0_R8
    call mct_gGrid_importRAttr(dom_lnd,"mask" ,data,lsz)

    ! Fill in correct values for domain components from the global grid
    do n = 1, lsz
       ni = natural_id_cells_owned(n)
       data(n) = lonc_g(ni)
    end do
    call mct_gGrid_importRattr(dom_lnd,"lon",data,lsz)

    do n = 1, lsz
       ni = natural_id_cells_owned(n)
       data(n) = latc_g(ni)
    end do
    call mct_gGrid_importRattr(dom_lnd,"lat",data,lsz)

    do n = 1, lsz
       ni = natural_id_cells_owned(n)
       data(n) = areac_g(ni)
    end do
    call mct_gGrid_importRattr(dom_lnd,"area",data,lsz)

    do n = 1, lsz
       ni = natural_id_cells_owned(n)
       data(n) = maskc_g(ni)
    end do
    call mct_gGrid_importRattr(dom_lnd,"mask",data,lsz)

    do n = 1, lsz
       ni = natural_id_cells_owned(n)
       data(n) = fracc_g(ni)
    end do
    call mct_gGrid_importRattr(dom_lnd,"frac",data,lsz)

    deallocate(data)
    deallocate(idata)

  end subroutine lnd_domain_mct

  !===============================================================================

  real(r8) function elmxx_caldayp1(now, step)

    ! Calendar day of the NEXT model step, i.e. ELM's
    ! get_curr_calday(offset=dtime). Day-of-year plus the fraction of the day
    ! elapsed, 1-based, matching the convention nextsw_cday uses.

    type(ESMF_Time)        , intent(in) :: now
    type(ESMF_TimeInterval), intent(in) :: step

    type(ESMF_Time) :: nxt
    integer :: rc, doy, hh, mm, ss

    ! ESMF_TimeGet's h/m/s are COMPONENTS, not seconds-of-day: asking for s
    ! alone on 00:30:00 returns 0, not 1800.
    nxt = now + step
    call ESMF_TimeGet(nxt, dayOfYear=doy, h=hh, m=mm, s=ss, rc=rc)
    call chkrc(rc, 'lnd::elmxx_caldayp1: ESMF_TimeGet')
    elmxx_caldayp1 = real(doy, r8) + real(hh*3600 + mm*60 + ss, r8)/86400._r8

  end function elmxx_caldayp1

  !===============================================================================

  integer function get_step_size(EClock)

    ! Return the step size in seconds.

    type(ESMF_Clock) :: EClock

    type(ESMF_TimeInterval)     :: step_size       ! timestep size
    integer                     :: rc
    character(len=*), parameter :: sub = 'lnd::get_step_size'

    call ESMF_ClockGet(EClock, timeStep=step_size, rc=rc)
    call chkrc(rc, sub//': error return from ESMF_ClockGet')

    call ESMF_TimeIntervalGet(step_size, s=get_step_size, rc=rc)
    call chkrc(rc, sub//': error return from ESMF_ClockTimeIntervalGet')

  end function get_step_size

  !===============================================================================

  subroutine get_clock_date(EClock, month, day, year)

    ! Extract the component clock date used by satellite phenology.  The
    ! EClock passed into the run phase is already at the end of this coupling
    ! interval, matching ELM's get_curr_date(offset=dtime) convention.
    !
    ! year is optional and, until history needed it, was computed here and
    ! discarded -- elmxx_init's hist_year argument now uses it.

    type(ESMF_Clock), intent(inout) :: EClock
    integer, intent(out) :: month, day
    integer, intent(out), optional :: year
    type(ESMF_Time) :: current_time
    integer :: rc, yy, seconds

    call ESMF_ClockGet(EClock, currTime=current_time, rc=rc)
    call chkrc(rc, 'lnd::get_clock_date: error return from ESMF_ClockGet')
    call ESMF_TimeGet(current_time, yy=yy, mm=month, dd=day, s=seconds, rc=rc)
    call chkrc(rc, 'lnd::get_clock_date: error return from ESMF_TimeGet')
    if (present(year)) year = yy

  end subroutine get_clock_date

  !===============================================================================

  subroutine chkrc(rc, mes)

    integer, intent(in)          :: rc   ! return code from time management library
    character(len=*), intent(in) :: mes  ! error message

    if ( rc == ESMF_SUCCESS ) return

    write(logunit_lnd,*) mes

    call shr_sys_abort ('CHKRC')

  end subroutine chkrc

  !===============================================================================

end module lnd_comp_mct
