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
  use shr_kind_mod    , only : IN=>SHR_KIND_IN, R8=>SHR_KIND_R8, CS=>SHR_KIND_CS, CL=>SHR_KIND_CL
  use shr_sys_mod     , only : shr_sys_abort, shr_sys_flush
  use shr_file_mod    , only : shr_file_getunit, shr_file_getlogunit, shr_file_getloglevel
  use shr_file_mod    , only : shr_file_setlogunit, shr_file_setloglevel, shr_file_setio
  use shr_file_mod    , only : shr_file_freeunit
  use elmxxSpmdMod    , only : masterproc, mpicom_lnd, iam, npes, LNDID, elmxxSpmdInit
  use elmxxMod        , only : elmxx_read_namelist, elmxx_init, elmxx_run, elmxx_final
  use shr_orb_mod     , only : shr_orb_decl, SHR_ORB_UNDEF_REAL
  use elmxxMod        , only : num_cells_owned, num_cells_global, natural_id_cells_owned
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
    integer :: month, day
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

    call get_clock_date(EClock, month, day)
    call elmxx_init(logunit_lnd, month, day)

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

    call elmxx_run(logunit_lnd, coupling_dt_in_sec, month, day, &
                   nextsw_cday, declinp1)

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

    if (.not. do_elmxx) return

    call elmxx_final()

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

  subroutine get_clock_date(EClock, month, day)

    ! Extract the component clock date used by satellite phenology.  The
    ! EClock passed into the run phase is already at the end of this coupling
    ! interval, matching ELM's get_curr_date(offset=dtime) convention.

    type(ESMF_Clock), intent(inout) :: EClock
    integer, intent(out) :: month, day
    type(ESMF_Time) :: current_time
    integer :: rc, year, seconds

    call ESMF_ClockGet(EClock, currTime=current_time, rc=rc)
    call chkrc(rc, 'lnd::get_clock_date: error return from ESMF_ClockGet')
    call ESMF_TimeGet(current_time, yy=year, mm=month, dd=day, s=seconds, rc=rc)
    call chkrc(rc, 'lnd::get_clock_date: error return from ESMF_TimeGet')

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
