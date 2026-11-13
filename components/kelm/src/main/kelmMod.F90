module kelmMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! kELM driver-facing state and entry points.
  !
  ! Structured after components/rdycore/src/main/rdycoreMod.F90: this module owns
  ! the decomposition (num_cells_owned / num_cells_global / natural_id_cells_owned)
  ! and the init/run/final entry points, and the MCT coupling layer in
  ! src/cpl/lnd_comp_mct.F90 reads that state to build the gsMap and domain.
  !
  ! Unlike RDycore, there is no PETSc here, and the decomposition is not handed
  ! to us by a library -- for now kELM simply partitions the global land grid
  ! round-robin across ranks. Time stepping is a no-op until the Kokkos port is
  ! wired in.
  !-----------------------------------------------------------------------

  use shr_kind_mod , only : r8 => shr_kind_r8
  use shr_const_mod          , only : SHR_CONST_STEBOL
  use shr_sys_mod  , only : shr_sys_abort, shr_sys_flush
  use shr_file_mod , only : shr_file_getunit, shr_file_freeunit
  use shr_nl_mod   , only : shr_nl_find_group_name
  use kelmSpmdMod , only : masterproc, iam, npes, mpicom_lnd
  use kelmIO      , only : kelm_pio_init, kelm_read_domain
  use kelmTimerMod, only : kelm_t_start, kelm_t_stop
  use kelmSurfdataMod, only : kelm_read_surfdata, kelm_surfdata_clean, &
                               surfdata_read, numurbl, natpft, nlevsoi, &
                               pct_natveg, pct_crop, pct_lake, pct_wetland, &
                               pct_glacier, pct_urban
  use kelmSubgridMod , only : kelm_build_subgrid, kelm_subgrid_clean, &
                               subgrid_built, num_landunits, num_columns, &
                               num_patches, lun_itype, col_landunit, &
                               istsoil, isturb_tbd, isturb_hd, isturb_md, &
                               num_landunits_global, num_columns_global, &
                               num_patches_global
  use kelmSurfaceStateMod, only : kelm_surface_state_init, &
                                   kelm_surface_state_clean
  use kelmFilterMod      , only : kelm_build_filters, kelm_filters_clean
  use kelmInitCheckMod   , only : kelm_write_init_snapshot
  use kelmForcingMod , only : kelm_forcing_init, kelm_forcing_clean

  use kelm_mod              , only : kELMType, kELMCreate, kELMDestroy, KELM_SUCCESS, &
                                      kELMStepInputs, kELMRunStep, kELMRunStepFailedKernel, &
                                      kELMSetGroundHeatFluxSb, &
                                      kELMComputeSurfaceAlbedoNatural, &
                                      kELMComputeSurfaceAlbedoLake, &
                                      kELMComputeSurfaceAlbedoUrban, &
                                      kELMComputeForcingDerived, kELMAllocateCells
  use kelmSoilPropMod       , only : kelm_soil_prop_init, kelm_soil_prop_clean, &
                                      nlevtot, nlevgrnd
  use kelmPftconMod         , only : kelm_read_pftcon, kelm_pftcon_clean, pftcon_read
  use kelmSnicarMod         , only : kelm_snicar_read, kelm_snicar_clean
  use kelmFinidatMod        , only : kelm_finidat_read, kelm_finidat_clean, &
                                      kelm_finidat_apply_soilprop
  use kelmRootMod           , only : kelm_root_init, &
                                      kelm_root_clean, root_built
  use kelmPhotosynMod     , only : kelm_photosyn_init, kelm_photosyn_seed, &
                                    photosyn_built, &
                                    kelm_push_photosyn_statics, t10_period
  use kelmSurfaceAlbedoMod, only : &
                                      kelm_surface_albedo_report
  use kelmSoilKernelMod   , only : kelm_soil_kernel_init, &
                                      kelm_soil_kernel_clean, soil_kernel_built
  use kelmDiagnosticsMod , only : kelm_diag_init, kelm_diag_finalize,   &
                                   kelm_diag_new_timestep,                &
                                   kelm_diag_snapshot_state,              &
                                   kelm_diag_snapshot_fluxes,             &
                                   kelm_diag_write_maps, kelm_diag_flush
  use kelmRestMod           , only : kelm_rest_read
  use kelmHistMod           , only : kelm_hist_init, kelm_hist_step, &
                                      kelm_hist_write_if_month_end, &
                                      kelm_hist_final
  use kelmLakeMod           , only : kelm_lake_init, kelm_lake_push_forcing, &
                                      kelm_lake_clean, lake_built
  use kelmUrbanMod          , only : kelm_urban_init, kelm_urban_diag_dump, &
                                      kelm_urban_clean, urban_built, kelm_urban_diag_step, &
                                      kelm_urban_apply_finidat
  use kelmKernelMod         , only : kelm_kernels_activate_all, kelm_kernels_push, &
                                      kelm_wbal_report, &
                                      kelm_kernels_report, kelm_report_cantemp, &
                                      kelm_report_fluxes, kelm_report_surfrad, &
                                      K_SOILTEMP, K_SOILFLUX, K_SURFRUNOFF, &
                                      K_ROOTWATER, K_HYDRODRAIN, kernel_active, &
                                      any_kernel_active, NKERNEL, kernel_name
  use kelmKokkosStateMod    , only : kelm_kokkos_state_init, &
                                      kelm_kokkos_check_map_invariants, &
                                      kelm_kokkos_seed_topology, &
                                      kelm_kokkos_seed_state, &
                                      kelm_kokkos_seed_canopy_hydrology, &
                                      kelm_kokkos_push_top_solar_rad, &
                                      kelm_kokkos_seed_albedo, &
                                      kelm_kokkos_seed_pftpar, &
                                      kelm_kokkos_seed_stomata_closed, &
                                      kelm_kokkos_seed_soil_properties, &
                                      kelm_kokkos_push_forcing, &
                                      kelm_kokkos_push_latlon, &
                                      kelm_push_monthly_phenology, &
                                      kelm_kokkos_push_root_statics, &
                                      kelm_kokkos_push_snicar_tables, &
                                      kelm_kokkos_push_finidat, &
                                      kelm_kokkos_reseed_finidat_snow, &
                                      kelm_kokkos_verify_maps, &
                                      kelm_kokkos_state_clean, &
                                      kokkos_state_built, n_kokkos_col, &
                                      n_kokkos_patch, n_kokkos_urb
  use kelm_kokkos_interface , only : kELMKokkosInitialize, kELMKokkosFinalize, &
                                      kELMKokkosPrintConfiguration, kELMTimersFinalize

  implicit none
  save
  private

#include <mpif.h>

  !--------------------------------------------------------------------------
  ! Decomposition -- read by src/cpl/lnd_comp_mct.F90
  !--------------------------------------------------------------------------
  integer , public          :: num_cells_owned           ! active land cells owned by this rank
  integer , public          :: num_cells_global          ! total active land cells (mask==1)
  integer , public, pointer :: natural_id_cells_owned(:) ! 1-based global grid IDs owned by this rank

  !--------------------------------------------------------------------------
  ! Perf mode (plan Stage 6.8). Set by the standalone perf driver
  ! (src/perf/kelm_perf.F90) before kelm_init; never by the coupled model,
  ! where it stays 0 and nothing below changes.
  !
  ! With kelm_perf_ncells = N > 0, the domain's active cells are SOURCE
  ! cells (nsrc of them) and the decomposition deals N virtual cells
  ! round-robin, exactly as it deals real ones. Virtual cell g copies source
  ! cell mod(g-1, nsrc) + 1: natural_id_cells_owned holds that source's grid
  ! ID, so it repeats, and every reader that indexes its file through it
  ! (domain, surfdata) replicates the cell with no change of its own.
  ! What assumes unique IDs is switched off: history, the per-step
  ! init snapshot, the step reports, and finidat (its counts are per file).
  !--------------------------------------------------------------------------
  integer , public          :: kelm_perf_ncells = 0
  integer , public          :: kelm_perf_nsrc   = 0
  integer , public, pointer :: kelm_perf_src_owned(:) => null() ! source ordinal of each owned cell

  !--------------------------------------------------------------------------
  ! Global grid, read from the land domain file on every rank
  !--------------------------------------------------------------------------
  integer , public          :: nlon_g                    ! ni
  integer , public          :: nlat_g                    ! nj
  real(r8), public, pointer :: lonc_g(:)                 ! cell center longitude (deg)
  real(r8), public, pointer :: latc_g(:)                 ! cell center latitude  (deg)
  real(r8), public, pointer :: areac_g(:)                ! cell area (radians^2)
  real(r8), public, pointer :: maskc_g(:)                ! domain mask (0 or 1)
  real(r8), public, pointer :: fracc_g(:)                ! land fraction

  ! Per-owned-cell centres, in DEGREES, gathered from the global arrays.
  ! SurfaceAlbedo needs them for the solar zenith angle, and it indexes by
  ! local cell like the forcing does.
  real(r8), public, pointer :: cell_lat(:) => null()
  real(r8), public, pointer :: cell_lon(:) => null()

  !--------------------------------------------------------------------------
  ! Namelist (&kelm_inparm in lnd_in)
  !--------------------------------------------------------------------------
  logical           , public :: do_kelm   = .true.
  character(len=256), public :: fatmlndfrc = ' '

  ! Run the ported SurfaceAlbedo at the end of each step. Off by default, so
  ! that turning it on is a deliberate act and a run without it keeps the
  ! frozen cold-start albedos it has been graded against.
  logical, public :: kelm_do_albedo = .false.

  ! Hold stomata closed instead of leaving them wide open. Photosynthesis is
  ! not ported, so rssun/rssha are zero and transpiration is unbounded from
  ! above; this pins them at ELM's closed limit instead. Off by default: it is
  ! a placeholder, and the two settings bracket the truth rather than either
  ! being right. See kelm_kokkos_seed_stomata_closed.
  logical, public :: kelm_stomata_closed = .false.

  ! Run the ported Photosynthesis inside the CanopyFluxes Newton iteration.
  ! On by default, like kelm_do_albedo. Requires kelm_do_albedo (vcmaxcint
  ! is SurfaceAlbedo's output) and fparamfile, both checked at init.
  logical, public :: kelm_do_photosynthesis = .true.
  ! Atmospheric CO2 [ppmv]. ELM takes this from its own namelist with
  ! co2_type = 'constant'; the I1850 twin uses 284.7.
  real(r8), public :: kelm_co2_ppmv = 284.7_r8
  character(len=256), public :: fsurdat    = ' '
  ! SNICAR lookup tables. Both must be set for SNICAR to run; leave either
  ! blank and snow albedo stays at the cold-start 0.6 constant, which is
  ! exactly how kELM behaved before SNICAR was ported.
  character(len=256), public :: fsnowoptics = ' '
  character(len=256), public :: fsnowaging  = ' '
  ! ELM restart file read as an initial condition. READ-ONLY: kELM never
  ! writes restarts. This exists so both models can be started from
  ! bit-identical state and stepped once, which is the only way to separate a
  ! coupling error from an arithmetic one -- the replay tests bypass the
  ! coupling by construction and so cannot see it.
  character(len=256), public :: finidat = ' '
  ! Stage 6 restart. kelm_nrevsn is ELM's nrevsn: the kELM restart a BRANCH
  ! run starts from. A CONTINUE run reads rpointer.lnd instead. The read
  ! itself is deferred to the top of the first kelm_run (see there).
  character(len=256), public :: kelm_nrevsn = ' '
  character(len=256), private :: rest_read_file = ' '
  logical, private :: rest_read_pending = .false.
  logical, private :: rest_read_branch  = .false.   ! branch: history starts fresh
  ! Stage 3 boundary check. OFF by default since Stage 4: the probe overwrites
  ! state fields with fingerprints and restores them from the Fortran-side
  ! seed, which silently discards a step of physics. That was harmless while
  ! every kernel was a no-op; it is not harmless now. Set it in the namelist
  ! only for a run whose results you intend to throw away.
  logical           , public :: kelm_check_boundary  = .false.
  logical           , public :: kelm_stage2_snapshot = .false.
  logical           , public :: kelm_check_soft_fail = .false.
  ! C++ timers inside the Fortran-timed steps (src/utils/Timers.h in the
  ! library): 0 off, 1 the launches and phases of each step, 2 finer.
  integer           , public :: kelm_timer_level     = 0
  ! Per-kernel diagnostic tables ("(kelm_kernels_report) rank ..." and the
  ! SurfaceAlbedo report), written on step 1 and twice a model day. Off by
  ! default: they were for bringing kernels up against ELM one at a time.
  logical           , public :: kelm_kernel_reports  = .false.
  ! TOP sub-grid topographic solar radiation, ELM's switch of the same name
  ! (Hao et al. 2021): fsurdat must carry STDEV_ELEV (or STD_ELEV), SKY_VIEW,
  ! TERRAIN_CONFIG, SINSL_COSAS and SINSL_SINAS. ELM's fineTOP
  ! (use_finetop_rad) is not ported.
  logical           , public :: use_top_solar_rad    = .false.
  ! PFT parameter file (ELM's clm_params). Required: the kernels, all of which
  ! always run, read their PFT parameters from it.
  character(len=256), public :: fparamfile = ' '
  ! History output. Additions/overrides BEYOND the C++-side default field
  ! list (RegisterHistoryField's default_active) -- empty means just the
  ! defaults. See kelmHistMod.
  character(len=256), public :: kelm_hist_fincl = ' '
  ! Case name, read from infodata in lnd_init_mct (not read anywhere in
  ! kELM before history needed it for the output filename).
  character(len=256), public :: kelm_caseid = ' '

  !--------------------------------------------------------------------------
  ! Instance information
  !--------------------------------------------------------------------------
  character(len=16), public :: inst_name
  character(len=16), public :: inst_suffix   ! e.g. "_0001" or ""
  integer          , public :: inst_index

  integer, public :: iulog = 6

  ! Starts at -1 so the first driver pass is nstep 0, matching ELM: its
  ! lnd_run_mct loops until the clock syncs and so runs nstep 0 AND 1 on
  ! the first coupling call.
  integer, private :: nstep = -1
  logical, private :: root_statics_pushed = .false.
  logical, private :: ghf_sb_pushed = .false.
  logical, private :: photosyn_statics_pushed = .false.
  logical, private :: monthly_phen_pushed = .false.
  logical, private :: kernels_pushed = .false.   ! kelm_kernels_push, once
  ! History's activation is deferred to kelm_run (see the comment at its
  ! call site there): it needs col_npfts/col_pfti/wtcol, which
  ! kelm_soil_kernel_init pushes on the first step, not kelm_init.
  logical, private :: hist_state_init_done = .false.
  integer, private :: kelm_hist_start_year = 0
  integer, private :: kelm_hist_start_month = 0
  integer, private :: kelm_hist_start_day = 0

  !--------------------------------------------------------------------------
  ! kELM Kokkos/C++ model object
  !
  ! Stage 2 materializes the host-side subgrid and surface state before this
  ! handle is created. No state crosses the C API and no kernel is called yet;
  ! that one-time persistent-state transfer belongs to Stage 3. Without a
  ! surface dataset the component still retains the Stage 1 fallback counts so
  ! a domain-only coupling smoke test remains possible.
  !--------------------------------------------------------------------------
  type(kELMType), public :: kelm_state
  logical, private        :: kelm_state_created = .false.

  public :: kelm_read_namelist
  public :: kelm_init
  public :: kelm_run
  public :: kelm_final
  public :: kelm_restart_set_start
  public :: kelm_restart_write
  public :: kelm_restart_load_export

contains

  !-----------------------------------------------------------------------
  subroutine kelm_read_namelist(logunit)
    !
    ! !DESCRIPTION:
    ! Read &kelm_inparm from lnd_in on the master task and broadcast it.
    !
    implicit none
    !
    integer, intent(in) :: logunit
    !
    character(len=256) :: nlfilename
    integer            :: ier, unitn
    logical            :: lexist
    character(len=*), parameter :: subname = '(kelm_read_namelist) '

    namelist /kelm_inparm/ do_kelm, fatmlndfrc, fsurdat, &
                            kelm_check_boundary, kelm_check_soft_fail, &
                            kelm_stage2_snapshot, fparamfile, kelm_do_albedo, &
                            kelm_stomata_closed, kelm_do_photosynthesis, &
                            kelm_co2_ppmv, fsnowoptics, fsnowaging, finidat, &
                            kelm_hist_fincl, kelm_nrevsn, kelm_timer_level, &
                            kelm_kernel_reports, use_top_solar_rad

    ! defaults
    do_kelm   = .true.
    fatmlndfrc = ' '
    fsurdat    = ' '
    fsnowoptics = ' '
    fsnowaging  = ' '
    finidat     = ' '
    kelm_nrevsn = ' '
    kelm_check_boundary  = .false.
    kelm_stage2_snapshot = .false.
    kelm_check_soft_fail = .false.
    kelm_kernel_reports  = .false.
    use_top_solar_rad    = .false.
    fparamfile            = ' '
    kelm_do_albedo       = .true.
    kelm_stomata_closed  = .false.
    kelm_do_photosynthesis = .true.
    kelm_co2_ppmv        = 284.7_r8
    kelm_hist_fincl      = ' '
    kelm_timer_level     = 0

    nlfilename = "lnd_in" // trim(inst_suffix)

    inquire (file = trim(nlfilename), exist = lexist)
    if ( .not. lexist ) then
       write(logunit,*) subname//' ERROR: namelist file does NOT exist: '//trim(nlfilename)
       call shr_sys_abort(subname//' ERROR: '//trim(nlfilename)//' does not exist')
    end if

    if (masterproc) then
       unitn = shr_file_getunit()
       write(logunit,*) 'Read in kelm_inparm namelist from: ', trim(nlfilename)
       open( unitn, file=trim(nlfilename), status='old', iostat=ier )
       if (ier /= 0) call shr_sys_abort(subname//' ERROR opening '//trim(nlfilename))

       call shr_nl_find_group_name(unitn, 'kelm_inparm', status=ier)
       if (ier /= 0) then
          call shr_sys_abort(subname//' ERROR: kelm_inparm group not found in '//trim(nlfilename))
       end if

       read(unitn, kelm_inparm, iostat=ier)
       if (ier /= 0) then
          call shr_sys_abort(subname//' ERROR reading kelm_inparm from '//trim(nlfilename))
       end if

       close(unitn)
       call shr_file_freeunit(unitn)
    end if

    call mpi_bcast (do_kelm  , 1                , MPI_LOGICAL  , 0, mpicom_lnd, ier)
    call mpi_bcast (fatmlndfrc, len(fatmlndfrc)  , MPI_CHARACTER, 0, mpicom_lnd, ier)
    call mpi_bcast (fsurdat   , len(fsurdat)     , MPI_CHARACTER, 0, mpicom_lnd, ier)
    call mpi_bcast (fsnowoptics, len(fsnowoptics) , MPI_CHARACTER, 0, mpicom_lnd, ier)
    call mpi_bcast (fsnowaging , len(fsnowaging)  , MPI_CHARACTER, 0, mpicom_lnd, ier)
    call mpi_bcast (finidat    , len(finidat)     , MPI_CHARACTER, 0, mpicom_lnd, ier)
    call mpi_bcast (kelm_nrevsn, len(kelm_nrevsn), MPI_CHARACTER, 0, mpicom_lnd, ier)
    call mpi_bcast (kelm_check_boundary , 1      , MPI_LOGICAL  , 0, mpicom_lnd, ier)
    call mpi_bcast (kelm_stage2_snapshot, 1      , MPI_LOGICAL  , 0, mpicom_lnd, ier)
    call mpi_bcast (kelm_check_soft_fail, 1      , MPI_LOGICAL  , 0, mpicom_lnd, ier)
    call mpi_bcast (kelm_kernel_reports , 1      , MPI_LOGICAL  , 0, mpicom_lnd, ier)
    call mpi_bcast (fparamfile   , len(fparamfile)   , MPI_CHARACTER, 0, mpicom_lnd, ier)
    call mpi_bcast (kelm_hist_fincl, len(kelm_hist_fincl), MPI_CHARACTER, 0, mpicom_lnd, ier)
    ! These four went unbroadcast until 2026-10-04: every rank but 0 ran with
    ! albedo and photosynthesis OFF. Every namelist member needs a line here.
    call mpi_bcast (kelm_do_albedo        , 1, MPI_LOGICAL, 0, mpicom_lnd, ier)
    call mpi_bcast (kelm_stomata_closed   , 1, MPI_LOGICAL, 0, mpicom_lnd, ier)
    call mpi_bcast (kelm_do_photosynthesis, 1, MPI_LOGICAL, 0, mpicom_lnd, ier)
    call mpi_bcast (kelm_co2_ppmv         , 1, MPI_REAL8  , 0, mpicom_lnd, ier)
    call mpi_bcast (kelm_timer_level      , 1, MPI_INTEGER, 0, mpicom_lnd, ier)
    call mpi_bcast (use_top_solar_rad     , 1, MPI_LOGICAL, 0, mpicom_lnd, ier)

    if (masterproc) then
       write(logunit,*) ' '
       write(logunit,*) 'read from namelist:'
       write(logunit,*) '   do_kelm   = ', do_kelm
       write(logunit,*) '   fatmlndfrc = ', trim(fatmlndfrc)
       write(logunit,*) '   fsurdat    = ', trim(fsurdat)
       write(logunit,*) '   kelm_check_boundary  = ', kelm_check_boundary
       write(logunit,*) '   kelm_stage2_snapshot = ', kelm_stage2_snapshot
       write(logunit,*) '   kelm_check_soft_fail = ', kelm_check_soft_fail
       write(logunit,*) '   kelm_kernel_reports  = ', kelm_kernel_reports
       write(logunit,*) '   kelm_do_albedo       = ', kelm_do_albedo
       write(logunit,*) '   kelm_do_photosynthesis = ', kelm_do_photosynthesis
       write(logunit,*) '   use_top_solar_rad      = ', use_top_solar_rad
       write(logunit,*) '   fparamfile            = ', trim(fparamfile)
       write(logunit,*) '   fsnowoptics           = ', trim(fsnowoptics)
       write(logunit,*) '   fsnowaging            = ', trim(fsnowaging)
       write(logunit,*) '   kelm_hist_fincl      = ', trim(kelm_hist_fincl)
       write(logunit,*) '   kelm_timer_level     = ', kelm_timer_level
       call shr_sys_flush(logunit)
    end if

  end subroutine kelm_read_namelist

  !-----------------------------------------------------------------------
  subroutine kelm_init(logunit, hist_year, month, day)
    !
    ! !DESCRIPTION:
    ! Initialize kELM: read the land domain and build the round-robin
    ! decomposition of the global grid.
    !
    implicit none
    !
    integer, intent(in) :: logunit, hist_year, month, day
    !
    integer :: i, n, k
    integer :: num_cells_grid                  ! ni*nj, including non-land cells
    integer :: ierr_kelm                      ! kELM C API status
    integer :: n_nat_col, n_nat_patch, n_urb_lun
    integer :: nfail_maps
    integer, allocatable :: land_ids(:)        ! global grid IDs of the active land cells
    character(len=*), parameter :: subname = '(kelm_init) '

    if (masterproc) then
       ! ELM's opening line (elm_initializeMod.F90); kELM's closing pair,
       ! "Successfully initialized the land model", is in lnd_init_mct.
       write(logunit,*)
       write(logunit,*) 'Attempting to initialize the land model .....'
       write(logunit,*)
       call shr_sys_flush(logunit)
    end if

    if (len_trim(fatmlndfrc) == 0) then
       call shr_sys_abort(subname//' ERROR: fatmlndfrc is not set in lnd_in')
    end if

    ! ---- global grid, read on every rank ----
    call kelm_pio_init(inst_name)
    call kelm_read_domain(logunit, fatmlndfrc, nlon_g, nlat_g, &
                           lonc_g, latc_g, areac_g, maskc_g, fracc_g)

    num_cells_grid = nlon_g * nlat_g

    ! ---- count total land gridcells ----
    ! Only cells with mask == 1 are active land, exactly as ELM counts numg in
    ! components/elm/src/main/decompInitMod.F90. Ocean cells carry mask == 0 and
    ! frac == 0 and take no part in the land decomposition.
    num_cells_global = 0
    do n = 1, num_cells_grid
       if (nint(maskc_g(n)) == 1) num_cells_global = num_cells_global + 1
    end do

    if (num_cells_global == 0) then
       call shr_sys_abort(subname//' ERROR: no active land cells (mask==1) in '//trim(fatmlndfrc))
    end if

    ! global grid IDs of the active land cells, in ascending order
    allocate(land_ids(num_cells_global))
    k = 0
    do n = 1, num_cells_grid
       if (nint(maskc_g(n)) == 1) then
          k = k + 1
          land_ids(k) = n
       end if
    end do

    ! Perf mode: the active cells just counted are the sources; the
    ! decomposition below deals kelm_perf_ncells virtual cells instead.
    if (kelm_perf_ncells > 0) then
       if (len_trim(finidat) > 0) then
          call shr_sys_abort(subname//' ERROR: perf mode is a cold start; '// &
               'finidat must be blank (its column/pft counts are per file)')
       end if
       kelm_perf_nsrc  = num_cells_global
       num_cells_global = kelm_perf_ncells
       if (masterproc) then
          write(logunit,*) subname,'PERF MODE: ',num_cells_global, &
                           ' virtual cells replicating ',kelm_perf_nsrc,' source cells'
       end if
    end if

    if (npes > num_cells_global) then
       write(logunit,*) subname,'ERROR: number of processes exceeds number of ', &
                        'land grid cells ',npes,num_cells_global
       call shr_sys_abort(subname//' ERROR: more MPI ranks than active land cells')
    end if

    ! ---- round-robin decomposition over the active land cells ----
    ! Rank iam owns the iam+1, iam+1+npes, iam+1+2*npes, ... land cells. Note the
    ! IDs stored are global *grid* IDs (indices into the ni*nj arrays), which is
    ! what the gsMap and the domain both need: the gsMap global size stays ni*nj
    ! so it matches the atm grid, which the coupler requires when the atm and lnd
    ! grids are the same (see seq_domain_mct.F90).
    num_cells_owned = num_cells_global / npes
    if (iam < mod(num_cells_global, npes)) num_cells_owned = num_cells_owned + 1

    allocate(natural_id_cells_owned(num_cells_owned))
    if (kelm_perf_ncells > 0) then
       allocate(kelm_perf_src_owned(num_cells_owned))
       do i = 1, num_cells_owned
          kelm_perf_src_owned(i) = mod(iam + (i-1)*npes, kelm_perf_nsrc) + 1
          natural_id_cells_owned(i) = land_ids(kelm_perf_src_owned(i))
       end do
    else
       do i = 1, num_cells_owned
          natural_id_cells_owned(i) = land_ids(iam + 1 + (i-1)*npes)
       end do
    end if

    deallocate(land_ids)

    if (masterproc) then
       write(logunit,*) subname,'grid cells = ',num_cells_grid, &
                        ' active land cells = ',num_cells_global
    end if
    write(logunit,*) subname,'rank ',iam,' owns ',num_cells_owned, &
                     ' of ',num_cells_global,' active land cells'
    call shr_sys_flush(logunit)

    nstep = -1

    ! Per-timestep atmospheric forcing lives at gridcell level, so it is sized
    ! from the decomposition and does not depend on the surface dataset.
    ! Per-cell centres for the solar zenith angle. Gathered here rather than
    ! looked up per step: the decomposition does not change.
    allocate(cell_lat(num_cells_owned), cell_lon(num_cells_owned))
    do i = 1, num_cells_owned
       cell_lat(i) = latc_g(natural_id_cells_owned(i))
       cell_lon(i) = lonc_g(natural_id_cells_owned(i))
    end do

    call kelm_forcing_init(num_cells_owned)

    !-----------------------------------------------------------------------
    ! Surface dataset.
    !
    ! Optional for now: without it kELM still has a decomposition and a domain,
    ! which is all Stage 1 needed, so a case with no fsurdat stays runnable
    ! rather than aborting. It becomes mandatory once subgrid construction
    ! depends on it.
    !-----------------------------------------------------------------------
    if (len_trim(fsurdat) > 0) then
       call kelm_t_start('kelm_read_surfdata')
       call kelm_read_surfdata(logunit, fsurdat, num_cells_grid, &
                                natural_id_cells_owned, kelm_perf_ncells == 0, &
                                read_top=use_top_solar_rad)
       call kelm_t_stop('kelm_read_surfdata')
       call kelm_report_composition(logunit)
       call kelm_build_subgrid(logunit, num_cells_owned)
       ! ELM's subgrid summary (decompInitMod.F90), for the lines kELM has:
       ! no topounits, cohorts, clumps or per-level gsMaps.
       if (masterproc) then
          write(logunit,*)
          write(logunit,*)' Surface Grid Characteristics'
          write(logunit,*)'   longitude points          = ',nlon_g
          write(logunit,*)'   latitude points           = ',nlat_g
          write(logunit,*)'   total number of gridcells = ',num_cells_global
          write(logunit,*)'   total number of landunits = ',num_landunits_global
          write(logunit,*)'   total number of columns   = ',num_columns_global
          write(logunit,*)'   total number of pfts      = ',num_patches_global
          write(logunit,*)
          call shr_sys_flush(logunit)
       end if
       call kelm_surface_state_init(logunit, month, day)
       call kelm_build_filters(logunit, num_cells_owned)
       if (kelm_stage2_snapshot .and. kelm_perf_ncells == 0) then
          call kelm_write_init_snapshot(logunit, month, day, natural_id_cells_owned)
       end if
    else
       if (masterproc) then
          write(logunit,*) subname,'no fsurdat in lnd_in; no surface dataset read'
          call shr_sys_flush(logunit)
       end if
    end if

    !-----------------------------------------------------------------------
    ! Bring up the Kokkos runtime and create the kELM model object.
    !-----------------------------------------------------------------------
    call kELMKokkosInitialize()

    if (masterproc) then
       write(logunit,*) subname,'Kokkos configuration:'
       call shr_sys_flush(logunit)
       call kELMKokkosPrintConfiguration()
       call shr_sys_flush(logunit)
    end if

    ! Counts come from the subgrid once there is one. Without a surface dataset
    ! there is no subgrid, so fall back to the Stage 1 placeholder rather than
    ! failing -- a domain-only case stays runnable.
    if (subgrid_built) then
       call kelm_count_for_create(n_nat_col, n_nat_patch, n_urb_lun)
    else
       n_nat_col   = num_cells_owned
       n_nat_patch = num_cells_owned
       n_urb_lun   = 0
    end if

    call kELMCreate(n_nat_col, n_nat_patch, n_urb_lun, kelm_state, ierr_kelm)
    if (ierr_kelm /= KELM_SUCCESS) then
       write(logunit,*) subname,'kELMCreate failed with status ',ierr_kelm
       call shr_sys_abort(subname//' ERROR: kELMCreate failed')
    end if
    kelm_state_created = .true.
    ! The cell-level container (forcing, derived quantities, statics) that
    ! every surface reads, before any surface is seeded.
    call kELMAllocateCells(kelm_state, num_cells_owned, ierr_kelm)
    if (ierr_kelm /= KELM_SUCCESS) call shr_sys_abort(subname//' ERROR: kELMAllocateCells failed')

    write(logunit,*) subname,'rank ',iam,' created kELM object: natural ', &
                     n_nat_col,' columns ',n_nat_patch,' patches, urban ', &
                     n_urb_lun,' landunits'
    call shr_sys_flush(logunit)

    !-----------------------------------------------------------------------
    ! Stage 3: build the packed Fortran<->Kokkos maps and push the one piece
    ! of topology the kernels dereference.
    !
    ! The maps are built by their own pass over the subgrid, so their extents
    ! are an INDEPENDENT count from the one kELMCreate was given. Requiring
    ! the two to agree is the check that matters: a setter whose length
    ! disagrees with the allocated view is rejected and leaves that view at
    ! zero (STATUS.md E.1), which is silent, so it must be impossible by
    ! construction rather than caught downstream.
    !-----------------------------------------------------------------------
    if (subgrid_built) then
       call kelm_kokkos_state_init(logunit)

       if (n_kokkos_col /= n_nat_col .or. n_kokkos_patch /= n_nat_patch .or. &
           n_kokkos_urb /= n_urb_lun) then
          write(logunit,*) subname,'ERROR: packed map extents ',n_kokkos_col, &
               n_kokkos_patch, n_kokkos_urb,' disagree with kELMCreate ', &
               n_nat_col, n_nat_patch, n_urb_lun
          call shr_sys_abort(subname//' ERROR: packed map extents disagree with kELMCreate')
       end if

       ! Semantics before plumbing: the invariant check needs no Kokkos, so
       ! run it before anything is pushed across the boundary. A map that is
       ! wrong about which entity is which should fail here, not survive to be
       ! round-tripped consistently by kelm_kokkos_verify_maps at step 1.
       call kelm_kokkos_check_map_invariants(logunit, nfail_maps)
       if (nfail_maps /= 0) then
          call shr_sys_abort(subname//' ERROR: packed map invariants violated')
       end if

       ! History's activation (kelm_hist_init) needs col_npfts/col_pfti/wtcol
       ! for its PATCH_TO_COL sum(wtcol)==1 check -- those are NOT seeded here
       ! (checked directly: kelmSoilKernelMod.F90 is the only place
       ! kELMSetColNpfts/ColPfti/Wtcol are ever called, from
       ! kelm_soil_kernel_init, which kelm_run's first step calls). So
       ! kelm_hist_init itself is deferred to there -- see kelm_run, right
       ! after the soil-kernel bring-up block. Interval start is recorded now.
       kelm_hist_start_year  = hist_year
       kelm_hist_start_month = month
       kelm_hist_start_day   = day

       call kelm_kokkos_seed_topology(kelm_state, logunit)
       call kelm_kokkos_seed_state(kelm_state, logunit)
       call kelm_kokkos_seed_canopy_hydrology(kelm_state, logunit)
       if (use_top_solar_rad) call kelm_kokkos_push_top_solar_rad(kelm_state, logunit)
       call kelm_kokkos_seed_albedo(kelm_state, logunit)

       ! Constant for the whole run, so seeded once here rather than per step.
       if (kelm_stomata_closed) then
          call kelm_kokkos_seed_stomata_closed(kelm_state, logunit)
       end if

       ! Soil hydraulic properties: derived from surfdata texture here, then
       ! pushed. Eight kernels read these, and none of them can run until
       ! they are real rather than zero.
       call kelm_soil_prop_init(logunit)
       call kelm_kokkos_seed_soil_properties(kelm_state, logunit)
       ! The lake surface (plan Stage 6.5): after the soil grid, which it
       ! shares, and the packed natural maps, through which it reads forcing.
       call kelm_lake_init(kelm_state, logunit)
       ! The urban surface's time-constant parameters (plan Stage 6.7, U1).
       call kelm_urban_init(kelm_state, logunit)

       ! SNICAR lookup tables. Pushed once; nothing crosses per step. Both
       ! files must be given -- with either blank, snow albedo stays at the
       ! cold-start constant and SnowAge_grain is a no-op, which is how kELM
       ! behaved before SNICAR was ported.
       if (len_trim(fsnowoptics) > 0 .and. len_trim(fsnowaging) > 0) then
          call kelm_snicar_read(fsnowoptics, fsnowaging, logunit)
          call kelm_kokkos_push_snicar_tables(kelm_state, logunit)
       else if (masterproc) then
          write(logunit,*) '(kelm_init) SNICAR tables not given; snow albedo '// &
               'held at the cold-start constant and grain aging disabled'
       end if

       ! Root profile. Needs the PFT parameters and the soil grid, so it comes
       ! after both. btran itself is per-step and is computed in kelm_run.
       if (len_trim(fparamfile) > 0) then
          call kelm_read_pftcon(logunit, fparamfile)
          call kelm_kokkos_seed_pftpar(kelm_state, logunit)
          call kelm_root_init(logunit)
       end if

       ! Photosynthesis. Needs the PFT parameters, and needs SurfaceAlbedo to
       ! be running because vcmaxcint is its output -- both checked here so a
       ! misconfiguration names itself instead of surfacing as a zero
       ! stomatal resistance, which is a legal number and looks like an
       ! answer.
       if (kelm_do_photosynthesis) then
          if (len_trim(fparamfile) == 0) then
             call shr_sys_abort(subname//' ERROR: kelm_do_photosynthesis '// &
                  'requires fparamfile')
          end if
          if (.not. kelm_do_albedo) then
             call shr_sys_abort(subname//' ERROR: kelm_do_photosynthesis '// &
                  'requires kelm_do_albedo -- vcmaxcintsun/sha are '// &
                  'SurfaceAlbedo outputs and nothing else computes them')
          end if
          ! The rest waits for the first step: kelm_photosyn_init needs
          ! dtime to size the 10-day running mean, and dtime is an argument to
          ! kelm_run, not to init. Same reason the soil kernel surface is
          ! built on step one.
       end if

       ! finidat LAST among the seeding, so it overwrites every cold-start
       ! value rather than racing them. The soil GRID is deliberately left to
       ! kelm_soil_prop_init above -- it is identical by construction and
       ! pushing it again would only add a way to get it wrong.
       if (len_trim(finidat) > 0) then
          call kelm_finidat_read(finidat, num_columns, num_patches, logunit)
          ! Soil-prop arrays FIRST: the first-step soil kernel init re-pushes
          ! from them, so seeding only the device would be undone one step in.
          call kelm_finidat_apply_soilprop(logunit)
          call kelm_kokkos_push_finidat(kelm_state, logunit)
          ! Urban state too (U8 probe); lake is still a cold start (STATUS L1).
          call kelm_urban_apply_finidat(kelm_state, logunit)
       end if

       ! Every kernel runs; there is no namelist switch. The library can
       ! still run a chosen subset -- kelm_kernels_parse takes a list of
       ! names -- for kernel-by-kernel checks against ELM. Activated after
       ! seeding, so a blocked kernel's abort names a prerequisite that
       ! genuinely could not be met rather than one not yet met at this point.
       !
       ! fparamfile is required: without it the per-patch PFT parameters are
       ! never seeded, and zero is a legal value that turns into NaN several
       ! kernels later (dleaf = 0 makes CanopyFluxes' rb zero).
       if (len_trim(fparamfile) == 0) then
          call shr_sys_abort(subname//' ERROR: fparamfile is not set; the '// &
               'kELM kernels need the PFT parameters')
       end if
       call kelm_kernels_activate_all(logunit)

       ! The five soil/hydrology kernels' surface is NOT built here. It needs
       ! the coupling timestep, and kelm_init does not have it -- dt is an
       ! argument to kelm_run. Built on the first step instead; see there.
    end if



    !-----------------------------------------------------------------------
    ! Diagnostic trace. Writes ELM's own ELMDIAG1 format so a free-running
    ! kELM run can be diffed against elm_diagnostics.bin with the existing
    ! tooling. Off unless KELM_DIAG is set in the environment.
    !-----------------------------------------------------------------------
    block
      character(len=256) :: diag_path
      integer :: dlen
      call get_environment_variable('KELM_DIAG', diag_path, dlen)
      if (dlen > 0) then
         call kelm_diag_init(trim(diag_path), .true.)
         call kelm_diag_write_maps()
         call kelm_urban_diag_dump()
         call kelm_diag_flush()
         if (masterproc) write(logunit,*) 'kELM: diagnostics -> ',trim(diag_path)
      end if
    end block

  end subroutine kelm_init

  !-----------------------------------------------------------------------
  subroutine kelm_count_for_create(n_nat_col, n_nat_patch, n_urb_lun)
    use kelmSubgridMod, only : num_patches, patch_column
    !
    ! !DESCRIPTION:
    ! Translate the subgrid into the three counts kELMCreate wants.
    !
    ! kELM's C++ side is organised per surface type, not as one flat begc:endc
    ! span, so it needs the natural-vegetation columns and patches and the urban
    ! landunit count -- not the totals. Lake and glacier are carried by the
    ! subgrid but are not part of this call: lake state is allocated separately
    ! (kELMAllocateLakeState) and glacier has no kernels at all.
    !
    implicit none
    !
    integer, intent(out) :: n_nat_col, n_nat_patch, n_urb_lun
    !
    integer :: l, p

    n_nat_col = 0; n_nat_patch = 0; n_urb_lun = 0

    do l = 1, num_landunits
       select case (lun_itype(l))
       case (istsoil)
          n_nat_col = n_nat_col + 1
       case (isturb_tbd, isturb_hd, isturb_md)
          n_urb_lun = n_urb_lun + 1
       end select
    end do

    ! Patches on natural-vegetation columns. Counted from the column side so it
    ! stays correct if a landunit ever gets more than one column.
    do p = 1, num_patches
       if (lun_itype(col_landunit(patch_column(p))) == istsoil) then
          n_nat_patch = n_nat_patch + 1
       end if
    end do

    ! kELMCreate rejects zero natural columns or patches. Every gridcell gets a
    ! natural-vegetation landunit, so this can only happen with no cells at all,
    ! which kelm_init already rules out.
    if (n_nat_col <= 0 .or. n_nat_patch <= 0) then
       call shr_sys_abort('(kelm_count_for_create) ERROR: no natural columns or patches')
    end if

  end subroutine kelm_count_for_create

  !-----------------------------------------------------------------------
  subroutine kelm_report_composition(logunit)
    !
    ! !DESCRIPTION:
    ! Summarize the subgrid composition just read, and check it is self
    ! consistent.
    !
    ! This is the first piece of the Stage 2 initialization comparison: what a
    ! cell is made of has to be right before anything is built on top of it. It
    ! is a summary only -- the per-cell dump that gets compared against ELM
    ! lands with subgrid construction itself.
    !
    implicit none
    !
    integer, intent(in) :: logunit
    !
    integer  :: i, n_nat, n_urb, n_lake, n_gla, n_wet, n_crop
    real(r8) :: total, worst
    real(r8), parameter :: pct_tol = 1.0e-6_r8   ! percentage points
    character(len=*), parameter :: subname = '(kelm_report_composition) '

    n_nat = 0; n_urb = 0; n_lake = 0; n_gla = 0; n_wet = 0; n_crop = 0
    worst = 0.0_r8

    do i = 1, num_cells_owned
       if (pct_natveg(i)  > 0.0_r8) n_nat  = n_nat  + 1
       if (pct_lake(i)    > 0.0_r8) n_lake = n_lake + 1
       if (pct_glacier(i) > 0.0_r8) n_gla  = n_gla  + 1
       if (pct_wetland(i) > 0.0_r8) n_wet  = n_wet  + 1
       if (pct_crop(i)    > 0.0_r8) n_crop = n_crop + 1
       if (sum(pct_urban(i,:)) > 0.0_r8) n_urb = n_urb + 1

       ! The landunit percentages partition the gridcell, so they must sum to
       ! 100. A cell that does not is a reading or indexing error -- most likely
       ! the wrong gridcell row -- and would otherwise surface much later as
       ! nonsensical subgrid weights.
       total = pct_natveg(i) + pct_crop(i) + pct_lake(i) + pct_wetland(i) &
             + pct_glacier(i) + sum(pct_urban(i,:))
       worst = max(worst, abs(total - 100.0_r8))
    end do

    write(logunit,*) subname,'rank ',iam,' cells with: natveg ',n_nat, &
                     ' urban ',n_urb,' lake ',n_lake,' glacier ',n_gla, &
                     ' wetland ',n_wet,' crop ',n_crop
    write(logunit,*) subname,'rank ',iam,' worst |sum(pct)-100| = ',worst
    call shr_sys_flush(logunit)

    if (worst > pct_tol) then
       write(logunit,*) subname,'ERROR: landunit percentages do not sum to 100'
       call shr_sys_abort(subname//' ERROR: inconsistent subgrid composition')
    end if

    ! Glacier is out of scope (no kernels and no Fortran fallback), so say so
    ! loudly where it appears rather than silently ignoring the area.
    if (n_gla > 0 .and. masterproc) then
       write(logunit,*) subname,'WARNING: ',n_gla,' cells have glacier area, ', &
                        'which kELM does not model; it is not yet excluded ', &
                        'from the subgrid'
       call shr_sys_flush(logunit)
    end if

  end subroutine kelm_report_composition

  !-----------------------------------------------------------------------
  subroutine kelm_run(logunit, coupling_dt_in_sec, month, day, &
                       nextsw_cday, declinp1, doalb, &
                       hist_year, hist_month, hist_day, hist_tod)
    !
    ! !DESCRIPTION:
    ! Advance kELM one coupling interval.
    !
    ! This is intentionally a no-op: the Kokkos port is not wired in yet, so
    ! kELM consumes nothing from the coupler and returns zeroed fields.
    !
    implicit none
    !
    integer, intent(in) :: logunit
    integer, intent(in) :: coupling_dt_in_sec
    real(r8), intent(in) :: nextsw_cday   ! calendar day of the next radiation step
    real(r8), intent(in) :: declinp1      ! solar declination for it, radians
    integer, intent(in) :: month, day
    logical, intent(in), optional :: doalb  ! .false. on ELM's nstep-0 pass
    ! History month-end trigger needs the per-substep END time, which is a
    ! different clock read from month,day above (that pair drives phenology
    ! and must not change meaning as a side effect of adding history --
    ! see the plan's "Two small changes outside the new module").
    integer, intent(in), optional :: hist_year, hist_month, hist_day, hist_tod

    logical :: do_albedo_this_step
    logical :: doalb_in            ! the driver's doalb, independent of config
    integer :: ierr_rs
    integer :: kfail
    type(kELMStepInputs) :: step_in

    doalb_in = .true.
    if (present(doalb)) doalb_in = doalb
    do_albedo_this_step = kelm_do_albedo .and. doalb_in

    ! RESTART READ, on the first step of a continue/branch run. It has to be
    ! HERE: before the forcing push (forcing views are restart fields, and
    ! this step's forcing must win) and before root stress (btran must come
    ! from the restored soil, not the cold start). The pieces a startup
    ! builds lazily later in its first step are built now so the read can
    ! overwrite them; their own lazy checks below then see them built.
    if (rest_read_pending) then
       if (root_built .and. .not. root_statics_pushed) then
          call kelm_kokkos_push_root_statics(kelm_state, logunit)
          root_statics_pushed = .true.
       end if
       if (kokkos_state_built .and. .not. soil_kernel_built) then
          if (kernel_active(K_SOILTEMP)   .or. kernel_active(K_SOILFLUX)  .or. &
              kernel_active(K_SURFRUNOFF) .or. kernel_active(K_ROOTWATER) .or. &
              kernel_active(K_HYDRODRAIN)) then
             call kelm_soil_kernel_init(kelm_state, &
                  real(coupling_dt_in_sec, r8), logunit)
          end if
       end if
       if (kokkos_state_built .and. soil_kernel_built .and. .not. hist_state_init_done) then
          call kelm_hist_init(kelm_state, logunit, kelm_caseid, kelm_hist_fincl, &
                               kelm_hist_start_year, kelm_hist_start_month, &
                               kelm_hist_start_day, n_kokkos_col, num_cells_global, &
                               natural_id_cells_owned, lonc_g, latc_g, areac_g)
          hist_state_init_done = .true.
       end if
       if (kokkos_state_built .and. any_kernel_active .and. kelm_do_photosynthesis) then
          if (.not. photosyn_built) then
             call kelm_photosyn_init(cell_lat, real(coupling_dt_in_sec, r8), logunit)
             call kelm_photosyn_seed(kelm_state, logunit)
          end if
          if (.not. photosyn_statics_pushed) then
             call kelm_push_photosyn_statics(kelm_state, logunit)
             photosyn_statics_pushed = .true.
          end if
       end if
       call kelm_rest_read(kelm_state, rest_read_file, rest_read_branch, &
                            num_cells_owned, num_cells_global, nstep, logunit)
       rest_read_pending = .false.
    end if

    nstep = nstep + 1
    call kelm_diag_new_timestep(nstep)

    ! ELM gates SatellitePhenology on doalb (elm_driver.F90, the non-CN,
    ! non-FATES branch), and does NOT call it during initialisation for this
    ! configuration -- the initialize2 call is behind use_fates .and.
    ! use_fates_sp. So ELM carries elai = esai = frac_veg_nosno = 0 until the
    ! first doalb step, which is nstep 2 here, and treats every patch as bare
    ! ground until then. kELM updated phenology unconditionally and so entered
    ! step 0 with a full canopy, routing patches through CanopyFluxes where ELM
    ! was still running BareGroundFluxes.
    !
    ! The Stage 2 snapshot is for checking a run against its ELM twin, so it is
    ! off unless kelm_stage2_snapshot asks for it. It is replaced every step so
    ! a run of any length can be checked at its last date, which costs about
    ! 17 KB of text per cell per step and is far too much for production.
    if (kelm_stage2_snapshot .and. subgrid_built .and. kelm_perf_ncells == 0) then
       call kelm_write_init_snapshot(logunit, month, day, natural_id_cells_owned)
    end if

    !-----------------------------------------------------------------------
    ! Stage 3 closing check, once, at the end of step 1: round-trip a
    ! per-entity fingerprint through the Kokkos views and prove the packed
    ! maps put it back where it came from.
    !
    ! Deliberately at the end of step 1 rather than inside init: the plan wants
    ! the boundary proven after a full coupling interval has been driven, with
    ! every kernel still off, so nothing between the set and the get could
    ! legitimately have changed a value.
    !-----------------------------------------------------------------------
    !-----------------------------------------------------------------------
    ! Crossing one of two: atmospheric forcing in. The forcing arrays were
    ! filled from x2l_l before this call, so this is where they reach the
    ! device. Nothing else crosses per step -- the seeded state does not
    ! change while the kernels are off.
    !-----------------------------------------------------------------------
    if (kokkos_state_built) then
       call kelm_t_start('kelm_push_forcing')
       call kelm_kokkos_push_forcing(kelm_state, nextsw_cday, declinp1, &
            cell_lat, cell_lon, logunit)
       call kelm_lake_push_forcing(kelm_state)
       call kelm_t_stop('kelm_push_forcing')
    end if

    ! Root water stress's statics, pushed once; the stress itself is the
    ! first thing kELMRunStep does.
    if (root_built .and. .not. root_statics_pushed) then
       call kelm_kokkos_push_root_statics(kelm_state, logunit)
       root_statics_pushed = .true.
    end if

    !-----------------------------------------------------------------------
    ! Stage 4: run the active kernels, in driver order.
    !-----------------------------------------------------------------------
    ! The five soil/hydrology kernels' surface, built on the first step
    ! because it needs dtime.
    !
    ! kELMInitSharedMetadata IS NOT OPTIONAL AND IS NOT ONLY ABOUT FILTERS.
    ! It is the only thing in the whole C API that assigns elm->dtime, and
    ! BuildSoilTemperatureNaturalView reads exactly that -- so skipping it
    ! would leave the soil column integrating on whatever dtime the object was
    ! constructed with. The canopy kernels do not expose this, since
    ! kELMComputeCanopyHydrology takes dtime as an explicit argument.
    if (kokkos_state_built .and. .not. soil_kernel_built) then
       if (kernel_active(K_SOILTEMP)   .or. kernel_active(K_SOILFLUX)  .or. &
           kernel_active(K_SURFRUNOFF) .or. kernel_active(K_ROOTWATER) .or. &
           kernel_active(K_HYDRODRAIN)) then
          call kelm_soil_kernel_init(kelm_state, &
               real(coupling_dt_in_sec, r8), logunit)
          ! soil_kernel_init pushes _p1 with an empty snow half; the collapse
          ! that follows would write that back over the seeded snowpack.
          if (len_trim(finidat) > 0) then
             call kelm_kokkos_reseed_finidat_snow(kelm_state, logunit)
          end if
       end if
    end if

    ! History's PATCH_TO_COL activation needs col_npfts/col_pfti/wtcol, which
    ! the block just above is what actually pushes (kelm_soil_kernel_init ->
    ! kELMSetColNpfts/ColPfti/Wtcol in kelmSoilKernelMod.F90) -- checked
    ! directly after this call failed with the invariant check seeing an
    ! all-zero wtcol when history was still being activated from kelm_init.
    ! If no soil/hydrology kernel is ever active, soil_kernel_built never
    ! becomes true and history is never activated either -- no crash, just no
    ! output, which is the same prerequisite the soil kernels themselves have.
    if (kokkos_state_built .and. soil_kernel_built .and. .not. hist_state_init_done &
        .and. kelm_perf_ncells == 0) then
       call kelm_hist_init(kelm_state, logunit, kelm_caseid, kelm_hist_fincl, &
                            kelm_hist_start_year, kelm_hist_start_month, &
                            kelm_hist_start_day, n_kokkos_col, num_cells_global, &
                            natural_id_cells_owned, lonc_g, latc_g, areac_g)
       hist_state_init_done = .true.
    end if

    ! The one-time pushes the step reads, then the step itself: since plan
    ! Stage 6.9 R2 one library call, kELMRunStep (kELMDriver.cpp), runs root
    ! stress, photosynthesis forcing, the kernels in ELM's driver order with
    ! the ground heat flux between the canopy and soil blocks, then
    ! satellite phenology and SurfaceAlbedo at the end of the step, where ELM
    ! runs them (both gated on doalb). The KELM_DIAG snapshots, KELM_WBAL
    ! marks and ELM-named timers between kernels come back through the step
    ! hook (kelm_kernels_push).
    !
    ! PHOTOSYNTHESIS: vcmaxcint is SurfaceAlbedo's output at the END of the
    ! step, and ELM runs no SurfaceAlbedo at initialization (only elm_drv
    ! calls it, on doalb steps), so neither model has one before the first
    ! albedo step; the albedos until then are InitCold's (0.2).
    if (kokkos_state_built .and. any_kernel_active) then
       if (kelm_do_photosynthesis) then
          if (.not. photosyn_built) then
             call kelm_photosyn_init(cell_lat, real(coupling_dt_in_sec, r8), logunit)
             call kelm_photosyn_seed(kelm_state, logunit)
          end if
          if (.not. photosyn_statics_pushed) then
             call kelm_push_photosyn_statics(kelm_state, logunit)
             photosyn_statics_pushed = .true.
          end if
       end if
       if (soil_kernel_built .and. .not. ghf_sb_pushed) then
          call kELMSetGroundHeatFluxSb(kelm_state, SHR_CONST_STEBOL, ierr_rs)
          ghf_sb_pushed = .true.
       end if
       ! The twelve monthly LAI/SAI fields, from the surface dataset; per step
       ! only the date crosses (the library forms the weights).
       if (subgrid_built .and. doalb_in .and. .not. monthly_phen_pushed) then
          call kelm_push_monthly_phenology(kelm_state, logunit)
          monthly_phen_pushed = .true.
       end if
    end if

    if (kokkos_state_built) then
       if (.not. kernels_pushed) then
          call kelm_kernels_push(kelm_state, logunit)
          kernels_pushed = .true.
       end if
       step_in%dtime             = real(coupling_dt_in_sec, r8)
       step_in%declinp1          = declinp1
       step_in%co2_ppmv          = kelm_co2_ppmv
       step_in%nstep             = nstep
       step_in%month             = month
       step_in%day               = day
       step_in%t10_period        = t10_period
       step_in%doalb             = merge(1, 0, doalb_in)
       step_in%do_albedo         = merge(1, 0, do_albedo_this_step)
       step_in%do_photosynthesis = merge(1, 0, kelm_do_photosynthesis)
       step_in%root_stress       = merge(1, 0, root_built)
       step_in%soil_built        = merge(1, 0, soil_kernel_built)
       step_in%subgrid_built     = merge(1, 0, subgrid_built)
       step_in%lake_built        = merge(1, 0, lake_built)
       step_in%urban_built       = merge(1, 0, urban_built)
       call kELMRunStep(kelm_state, step_in, ierr_rs)
       if (ierr_rs /= KELM_SUCCESS) then
          kfail = kELMRunStepFailedKernel(kelm_state)
          if (kfail >= 1 .and. kfail <= NKERNEL) then
             call shr_sys_abort('(kelm_run) ERROR: kELMRunStep failed in kernel "'// &
                  trim(kernel_name(kfail))//'"')
          else
             call shr_sys_abort('(kelm_run) ERROR: kELMRunStep failed (root stress, '// &
                  'photosynthesis forcing, ground heat flux, phenology or albedo)')
          end if
       end if
    end if

    if (kokkos_state_built .and. any_kernel_active) then
       if (do_albedo_this_step .and. kelm_kernel_reports .and. &
           (nstep == 1 .or. mod(nstep, 24) == 0) .and. kelm_perf_ncells == 0) then
          call kelm_surface_albedo_report(logunit)
       end if
       ! Report on the first step, then TWICE daily: step 48k lands at model
       ! midnight, where every SurfaceRadiation output is trivially zero.
       if (kelm_kernel_reports .and. (nstep == 1 .or. mod(nstep, 24) == 0) &
           .and. kelm_perf_ncells == 0) then
          call kelm_report_surfrad(kelm_state, n_kokkos_patch, logunit)
          call kelm_kernels_report(kelm_state, n_kokkos_patch, logunit)
          call kelm_report_cantemp(kelm_state, n_kokkos_col, logunit)
          call kelm_report_fluxes(kelm_state, n_kokkos_patch, logunit)
       end if
    end if

    ! Gated on the FIRST pass, not nstep == 1. The probe overwrites state with
    ! fingerprints and restores from the Fortran-side seed, so it is only
    ! harmless while the state still IS that seed. Once kELM started running
    ! ELM's extra nstep-0 pass, nstep == 1 became the second pass and this
    ! silently threw away a step of physics.
    if (kokkos_state_built .and. nstep == 0 .and. kelm_check_boundary) then
       call kelm_verify_kokkos_boundary(logunit)
    end if

    !-----------------------------------------------------------------------
    ! History: accumulate every step, write at month-end. AFTER the
    ! phenology/SurfaceAlbedo block above (not at kelm_diag_snapshot_fluxes,
    ! which runs before it) -- see the plan's "Corrections from review" for
    ! why placing this earlier would accumulate the previous step's
    ! ELAI/ESAI/LAISUN/LAISHA. Not gated on any_kernel_active: history should
    ! track whatever state exists, cold-start or not.
    !-----------------------------------------------------------------------
    if (kokkos_state_built .and. subgrid_built .and. kelm_perf_ncells == 0) then
       call kelm_t_start('elm_drv_io')
       call kelm_t_start('hbuf')
       call kelm_hist_step(kelm_state, logunit)
       call kelm_t_stop('hbuf')
       ! The nstep-0 pass is accumulated (ELM's January mean includes it --
       ! measured, H2OCAN to 1.3e-7) but never closes a tape: ELM skips the
       ! write at nstep 0 (histFileMod.F90:3415). Its stamp is the run start,
       ! day 1 00:00, which would otherwise write a one-sample file.
       if (nstep > 0 .and. present(hist_year) .and. present(hist_month) .and. &
           present(hist_day) .and. present(hist_tod)) then
          call kelm_t_start('elm_drv_io_htapes')
          call kelm_hist_write_if_month_end(kelm_state, hist_year, hist_month, &
                                             hist_day, hist_tod, logunit)
          call kelm_t_stop('elm_drv_io_htapes')
       end if
       call kelm_t_stop('elm_drv_io')
    end if

    ! The per-step log lines are lnd_run_mct's, in ELM's format.

  end subroutine kelm_run



  !-----------------------------------------------------------------------
  subroutine kelm_verify_kokkos_boundary(logunit)
    !
    ! Run the packed-map round-trip and make a failure fatal.
    !
    ! Fatal is the right default while the boundary is being brought up: a
    ! wrong index map is silent corruption, and every kernel added later builds
    ! on it. Once kernels run, this becomes the URBANxx-style `_check` mode the
    ! plan describes, with a namelist soft-fail so a long run need not abort.
    !
    implicit none
    integer, intent(in) :: logunit
    integer :: nfail
    character(len=*), parameter :: subname = '(kelm_verify_kokkos_boundary) '

    call kelm_kokkos_verify_maps(kelm_state, logunit, nfail)

    if (nfail /= 0) then
       if (kelm_check_soft_fail) then
          write(logunit,*) subname,'WARNING: ',nfail,' boundary mismatches; ', &
               'continuing because kelm_check_soft_fail is set. Every kernel ', &
               'reading this state is now suspect.'
          call shr_sys_flush(logunit)
       else
          call shr_sys_abort(subname//' ERROR: packed Fortran<->Kokkos maps do not round-trip')
       end if
    end if

    ! The probe writes fingerprints into real state fields -- t_grnd, t_veg,
    ! h2osno and friends -- so the seeded values have to be put back. Harmless
    ! today, since no kernel reads them yet, but leaving garbage in persistent
    ! state to be discovered at Stage 4 is exactly the kind of thing that gets
    ! blamed on the kernel. Re-seeding is exact here because the state is
    ! static while the kernels are off.
    call kelm_kokkos_seed_state(kelm_state, logunit)

  end subroutine kelm_verify_kokkos_boundary

  !-----------------------------------------------------------------------
  subroutine kelm_final(hist_year, hist_month, hist_day, hist_tod)
    !
    ! !DESCRIPTION:
    ! Finalize kELM.
    !
    implicit none
    integer, intent(in) :: hist_year, hist_month, hist_day, hist_tod

    integer :: ierr_kelm

    ! Stop the C++ timers before the Kokkos teardown below: they fence.
    call kELMTimersFinalize()

    ! Per-kernel water budget, if KELM_WBAL was set in the environment.
    ! Before teardown: it reads state the destroy below releases.
    call kelm_wbal_report(iulog)

    !-----------------------------------------------------------------------
    ! Tear down in reverse order of kelm_init: object first, then Kokkos.
    ! kELMKokkosFinalize must come after kELMDestroy -- the object owns
    ! Kokkos views, and destroying them after the runtime is gone is undefined.
    !-----------------------------------------------------------------------
    if (kelm_state_created) then
       ! History's partial-month flush needs kelm_state, so it must run
       ! BEFORE kELMDestroy, not after.
       if (subgrid_built .and. kelm_perf_ncells == 0) then
          call kelm_hist_final(kelm_state, hist_year, hist_month, hist_day, &
                                hist_tod, iulog)
       end if
       call kELMDestroy(kelm_state, ierr_kelm)
       if (ierr_kelm /= KELM_SUCCESS) then
          write(iulog,*) 'kelm_final: kELMDestroy failed with status ',ierr_kelm
       end if
       kelm_state_created = .false.
       call kelm_soil_kernel_clean()
       call kelm_kokkos_state_clean()
       call kelm_lake_clean()
       call kelm_urban_clean()
       call kelm_soil_prop_clean()
       call kelm_root_clean()
       call kelm_pftcon_clean()
       call kELMKokkosFinalize()
    end if

    call kelm_forcing_clean()
    call kelm_surface_state_clean()
    call kelm_filters_clean()
    call kelm_subgrid_clean()
    call kelm_surfdata_clean()

    if (associated(natural_id_cells_owned)) deallocate(natural_id_cells_owned)
    if (associated(kelm_perf_src_owned))   deallocate(kelm_perf_src_owned)
    if (associated(lonc_g))  deallocate(lonc_g)
    if (associated(latc_g))  deallocate(latc_g)
    if (associated(areac_g)) deallocate(areac_g)
    if (associated(maskc_g)) deallocate(maskc_g)
    if (associated(fracc_g)) deallocate(fracc_g)

    call kelm_diag_finalize()

  end subroutine kelm_final

  !-----------------------------------------------------------------------
  subroutine kelm_restart_set_start(start_type, inst_suffix, logunit, nstep_restart)
    !
    ! Decide, from the coupler's start type, whether this run starts from an
    ! kELM restart (ELM's nsrest): 'continue' reads rpointer.lnd, 'branch'
    ! reads kelm_nrevsn, 'startup' reads nothing (finidat, if set, is the
    ! separate ELM-restart path). Returns the restart's step counter, or -1,
    ! so the driver can carry its own count on; the read itself happens at
    ! the top of the first kelm_run.
    !
    use kelmRestMod, only : kelm_rest_rpointer_read, kelm_rest_peek_nstep
    character(len=*), intent(in)  :: start_type, inst_suffix
    integer         , intent(in)  :: logunit
    integer         , intent(out) :: nstep_restart
    character(len=*), parameter :: subname = '(kelm_restart_set_start) '

    nstep_restart = -1
    select case (trim(start_type))
    case ('continue')
       rest_read_file = kelm_rest_rpointer_read(inst_suffix)
    case ('branch')
       if (len_trim(kelm_nrevsn) == 0) call shr_sys_abort(subname// &
            'ERROR: a branch run needs kelm_nrevsn (the kELM restart to branch from)')
       rest_read_file = kelm_nrevsn
       rest_read_branch = .true.
    case default
       return
    end select
    rest_read_pending = .true.
    nstep_restart = kelm_rest_peek_nstep(rest_read_file)
    if (masterproc) then
       write(logunit,*) subname, trim(start_type), ' run from ', trim(rest_read_file), &
            ' (nstep ', nstep_restart, ')'
       call shr_sys_flush(logunit)
    end if
  end subroutine kelm_restart_set_start

  !-----------------------------------------------------------------------
  subroutine kelm_restart_write(inst_suffix, yr, mon, day, tod, logunit)
    !
    ! Write a restart dated (yr, mon, day, tod) -- the END of the step just
    ! taken, which is the time the next run starts from.
    !
    use kelmRestMod, only : kelm_rest_write
    character(len=*), intent(in) :: inst_suffix
    integer         , intent(in) :: yr, mon, day, tod, logunit
    if (.not. kokkos_state_built) return
    call kelm_rest_write(kelm_state, kelm_caseid, inst_suffix, yr, mon, day, tod, &
                          nstep, num_cells_owned, num_cells_global, logunit)
  end subroutine kelm_restart_write

  !-----------------------------------------------------------------------
  logical function kelm_restart_load_export(logunit)
    !
    ! On a continue or branch run, give the export the table the restart
    ! carries, so the coupler's first l2x is the writing run's next one (ELM
    ! exports its restart state at init). False on a startup, and for a
    ! restart written before the table existed: the init export then comes
    ! from the cold start, and lnd.log says so. Collective.
    !
    use kelmRestMod     , only : kelm_rest_read_lnd2atm
    use kelmLndExportMod, only : kelm_lnd2atm_load_cell_table, NL2A => NFIELD
    integer, intent(in) :: logunit
    real(r8), allocatable :: tab(:,:)
    logical :: found

    kelm_restart_load_export = .false.
    if (.not. rest_read_pending) return
    allocate(tab(num_cells_owned, NL2A))
    call kelm_rest_read_lnd2atm(kelm_state, rest_read_file, num_cells_owned, &
                                num_cells_global, tab, found)
    if (found) then
       call kelm_lnd2atm_load_cell_table(tab)
       kelm_restart_load_export = .true.
       if (masterproc) then
          write(logunit,*) '(kelm_restart_load_export) export at init taken from ', &
               trim(rest_read_file)
          call shr_sys_flush(logunit)
       end if
    else if (masterproc) then
       write(logunit,*) '(kelm_restart_load_export) WARNING: ',trim(rest_read_file), &
            ' carries no rst:lnd2atm; the export at init is the cold start''s'
       call shr_sys_flush(logunit)
    end if
    deallocate(tab)
  end function kelm_restart_load_export

end module kelmMod
