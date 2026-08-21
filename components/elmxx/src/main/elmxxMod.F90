module elmxxMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! ELMxx driver-facing state and entry points.
  !
  ! Structured after components/rdycore/src/main/rdycoreMod.F90: this module owns
  ! the decomposition (num_cells_owned / num_cells_global / natural_id_cells_owned)
  ! and the init/run/final entry points, and the MCT coupling layer in
  ! src/cpl/lnd_comp_mct.F90 reads that state to build the gsMap and domain.
  !
  ! Unlike RDycore, there is no PETSc here, and the decomposition is not handed
  ! to us by a library -- for now ELMxx simply partitions the global land grid
  ! round-robin across ranks. Time stepping is a no-op until the Kokkos port is
  ! wired in.
  !-----------------------------------------------------------------------

  use shr_kind_mod , only : r8 => shr_kind_r8
  use shr_const_mod          , only : SHR_CONST_STEBOL
  use shr_sys_mod  , only : shr_sys_abort, shr_sys_flush
  use shr_file_mod , only : shr_file_getunit, shr_file_freeunit
  use shr_nl_mod   , only : shr_nl_find_group_name
  use elmxxSpmdMod , only : masterproc, iam, npes, mpicom_lnd
  use elmxxIO      , only : elmxx_pio_init, elmxx_read_domain
  use elmxxSurfdataMod, only : elmxx_read_surfdata, elmxx_surfdata_clean, &
                               surfdata_read, numurbl, natpft, nlevsoi, &
                               pct_natveg, pct_crop, pct_lake, pct_wetland, &
                               pct_glacier, pct_urban
  use elmxxSubgridMod , only : elmxx_build_subgrid, elmxx_subgrid_clean, &
                               subgrid_built, num_landunits, num_columns, &
                               num_patches, lun_itype, col_landunit, &
                               istsoil, isturb_tbd, isturb_hd, isturb_md
  use elmxxSurfaceStateMod, only : elmxx_surface_state_init, &
                                   elmxx_surface_state_clean, &
                                   elmxx_push_monthly_phenology, elmxx_phenology_weights
  use elmxxFilterMod      , only : elmxx_build_filters, elmxx_filters_clean
  use elmxxInitCheckMod   , only : elmxx_write_init_snapshot
  use elmxxForcingMod , only : elmxx_forcing_init, elmxx_forcing_clean

  use elmxx_mod              , only : ELMxxType, ELMxxCreate, ELMxxDestroy, ELMXX_SUCCESS, &
                                      ELMxxComputeRootStressNatural, &
                                      ELMxxComputeGroundHeatFluxNatural, &
                                      ELMxxSetGroundHeatFluxSb, &
                                      ELMxxComputeSurfaceAlbedoNatural, &
                                      ELMxxComputePhotosynForcingNatural, &
                                      ELMxxComputePhenologyNatural
  use elmxxSoilPropMod       , only : elmxx_soil_prop_init, elmxx_soil_prop_clean, &
                                      nlevtot, nlevgrnd
  use elmxxPftconMod         , only : elmxx_read_pftcon, elmxx_pftcon_clean, pftcon_read
  use elmxxRootMod           , only : elmxx_root_init, &
                                      elmxx_root_clean, root_built
  use elmxxPhotosynMod     , only : elmxx_photosyn_init, elmxx_photosyn_seed, &
                                    photosyn_built, &
                                    elmxx_push_photosyn_statics, t10_period
  use elmxxSurfaceAlbedoMod, only : elmxx_push_coszen, &
                                      elmxx_surface_albedo_report
  use elmxxSoilKernelMod   , only : elmxx_soil_kernel_init, &
                                      elmxx_soil_kernel_clean, soil_kernel_built
  use elmxxDiagnosticsMod , only : elmxx_diag_init, elmxx_diag_finalize,   &
                                   elmxx_diag_new_timestep,                &
                                   elmxx_diag_snapshot_state,              &
                                   elmxx_diag_snapshot_fluxes,             &
                                   elmxx_diag_write_maps
  use elmxxKernelMod         , only : elmxx_kernels_parse, elmxx_kernels_run, &
                                      elmxx_kernels_report, elmxx_report_cantemp, &
                                      elmxx_report_fluxes, elmxx_report_surfrad, &
                                      K_SOILTEMP, K_SOILFLUX, K_SURFRUNOFF, &
                                      K_ROOTWATER, K_HYDRODRAIN, kernel_active, &
                                      any_kernel_active
  use elmxxKokkosStateMod    , only : elmxx_kokkos_state_init, &
                                      elmxx_kokkos_check_map_invariants, &
                                      elmxx_kokkos_seed_topology, &
                                      elmxx_kokkos_seed_state, &
                                      elmxx_kokkos_seed_canopy_hydrology, &
                                      elmxx_kokkos_seed_albedo, &
                                      elmxx_kokkos_seed_pftpar, &
                                      elmxx_kokkos_seed_stomata_closed, &
                                      elmxx_kokkos_seed_soil_properties, &
                                      elmxx_kokkos_push_forcing, &
                                      elmxx_kokkos_push_root_statics, &
                                      elmxx_kokkos_verify_maps, &
                                      elmxx_kokkos_state_clean, &
                                      kokkos_state_built, n_kokkos_col, &
                                      n_kokkos_patch, n_kokkos_urb
  use elmxx_kokkos_interface , only : ELMxxKokkosInitialize, ELMxxKokkosFinalize, &
                                      ELMxxKokkosPrintConfiguration

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
  ! Namelist (&elmxx_inparm in lnd_in)
  !--------------------------------------------------------------------------
  logical           , public :: do_elmxx   = .true.
  character(len=256), public :: fatmlndfrc = ' '

  ! Run the ported SurfaceAlbedo at the end of each step. Off by default, so
  ! that turning it on is a deliberate act and a run without it keeps the
  ! frozen cold-start albedos it has been graded against.
  logical, public :: elmxx_do_albedo = .false.

  ! Hold stomata closed instead of leaving them wide open. Photosynthesis is
  ! not ported, so rssun/rssha are zero and transpiration is unbounded from
  ! above; this pins them at ELM's closed limit instead. Off by default: it is
  ! a placeholder, and the two settings bracket the truth rather than either
  ! being right. See elmxx_kokkos_seed_stomata_closed.
  logical, public :: elmxx_stomata_closed = .false.

  ! Run the ported Photosynthesis inside the CanopyFluxes Newton iteration.
  ! Off by default, like elmxx_do_albedo. Requires elmxx_do_albedo (vcmaxcint
  ! is SurfaceAlbedo's output) and fparamfile, both checked at init.
  logical, public :: elmxx_do_photosynthesis = .false.
  ! Atmospheric CO2 [ppmv]. ELM takes this from its own namelist with
  ! co2_type = 'constant'; the I1850 twin uses 284.7.
  real(r8), public :: elmxx_co2_ppmv = 284.7_r8
  character(len=256), public :: fsurdat    = ' '
  ! Stage 3 boundary check. OFF by default since Stage 4: the probe overwrites
  ! state fields with fingerprints and restores them from the Fortran-side
  ! seed, which silently discards a step of physics. That was harmless while
  ! every kernel was a no-op; it is not harmless now. Set it in the namelist
  ! only for a run whose results you intend to throw away.
  logical           , public :: elmxx_check_boundary  = .false.
  logical           , public :: elmxx_check_soft_fail = .false.
  ! Stage 4: comma-separated kernel names; empty means the timestep is a no-op.
  character(len=256), public :: elmxx_kernels = ' '
  ! PFT parameter file (ELM's clm_params). Needed once root water stress runs.
  character(len=256), public :: fparamfile = ' '

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

  !--------------------------------------------------------------------------
  ! ELMxx Kokkos/C++ model object
  !
  ! Stage 2 materializes the host-side subgrid and surface state before this
  ! handle is created. No state crosses the C API and no kernel is called yet;
  ! that one-time persistent-state transfer belongs to Stage 3. Without a
  ! surface dataset the component still retains the Stage 1 fallback counts so
  ! a domain-only coupling smoke test remains possible.
  !--------------------------------------------------------------------------
  type(ELMxxType), public :: elmxx_state
  logical, private        :: elmxx_state_created = .false.

  public :: elmxx_read_namelist
  public :: elmxx_init
  public :: elmxx_run
  public :: elmxx_init_albedo
  public :: elmxx_final

contains

  !-----------------------------------------------------------------------
  subroutine elmxx_read_namelist(logunit)
    !
    ! !DESCRIPTION:
    ! Read &elmxx_inparm from lnd_in on the master task and broadcast it.
    !
    implicit none
    !
    integer, intent(in) :: logunit
    !
    character(len=256) :: nlfilename
    integer            :: ier, unitn
    logical            :: lexist
    character(len=*), parameter :: subname = '(elmxx_read_namelist) '

    namelist /elmxx_inparm/ do_elmxx, fatmlndfrc, fsurdat, &
                            elmxx_check_boundary, elmxx_check_soft_fail, &
                            elmxx_kernels, fparamfile, elmxx_do_albedo, &
                            elmxx_stomata_closed, elmxx_do_photosynthesis, &
                            elmxx_co2_ppmv

    ! defaults
    do_elmxx   = .true.
    fatmlndfrc = ' '
    fsurdat    = ' '
    elmxx_check_boundary  = .false.
    elmxx_check_soft_fail = .false.
    elmxx_kernels         = ' '
    fparamfile            = ' '
    elmxx_do_albedo       = .false.
    elmxx_stomata_closed  = .false.
    elmxx_do_photosynthesis = .false.
    elmxx_co2_ppmv        = 284.7_r8

    nlfilename = "lnd_in" // trim(inst_suffix)

    inquire (file = trim(nlfilename), exist = lexist)
    if ( .not. lexist ) then
       write(logunit,*) subname//' ERROR: namelist file does NOT exist: '//trim(nlfilename)
       call shr_sys_abort(subname//' ERROR: '//trim(nlfilename)//' does not exist')
    end if

    if (masterproc) then
       unitn = shr_file_getunit()
       write(logunit,*) 'Read in elmxx_inparm namelist from: ', trim(nlfilename)
       open( unitn, file=trim(nlfilename), status='old', iostat=ier )
       if (ier /= 0) call shr_sys_abort(subname//' ERROR opening '//trim(nlfilename))

       call shr_nl_find_group_name(unitn, 'elmxx_inparm', status=ier)
       if (ier /= 0) then
          call shr_sys_abort(subname//' ERROR: elmxx_inparm group not found in '//trim(nlfilename))
       end if

       read(unitn, elmxx_inparm, iostat=ier)
       if (ier /= 0) then
          call shr_sys_abort(subname//' ERROR reading elmxx_inparm from '//trim(nlfilename))
       end if

       close(unitn)
       call shr_file_freeunit(unitn)
    end if

    call mpi_bcast (do_elmxx  , 1                , MPI_LOGICAL  , 0, mpicom_lnd, ier)
    call mpi_bcast (fatmlndfrc, len(fatmlndfrc)  , MPI_CHARACTER, 0, mpicom_lnd, ier)
    call mpi_bcast (fsurdat   , len(fsurdat)     , MPI_CHARACTER, 0, mpicom_lnd, ier)
    call mpi_bcast (elmxx_check_boundary , 1      , MPI_LOGICAL  , 0, mpicom_lnd, ier)
    call mpi_bcast (elmxx_check_soft_fail, 1      , MPI_LOGICAL  , 0, mpicom_lnd, ier)
    call mpi_bcast (elmxx_kernels, len(elmxx_kernels), MPI_CHARACTER, 0, mpicom_lnd, ier)
    call mpi_bcast (fparamfile   , len(fparamfile)   , MPI_CHARACTER, 0, mpicom_lnd, ier)

    if (masterproc) then
       write(logunit,*) ' '
       write(logunit,*) 'read from namelist:'
       write(logunit,*) '   do_elmxx   = ', do_elmxx
       write(logunit,*) '   fatmlndfrc = ', trim(fatmlndfrc)
       write(logunit,*) '   fsurdat    = ', trim(fsurdat)
       write(logunit,*) '   elmxx_check_boundary  = ', elmxx_check_boundary
       write(logunit,*) '   elmxx_check_soft_fail = ', elmxx_check_soft_fail
       write(logunit,*) '   elmxx_kernels         = ', trim(elmxx_kernels)
       write(logunit,*) '   fparamfile            = ', trim(fparamfile)
       call shr_sys_flush(logunit)
    end if

  end subroutine elmxx_read_namelist

  !-----------------------------------------------------------------------
  subroutine elmxx_init(logunit, month, day)
    !
    ! !DESCRIPTION:
    ! Initialize ELMxx: read the land domain and build the round-robin
    ! decomposition of the global grid.
    !
    implicit none
    !
    integer, intent(in) :: logunit, month, day
    !
    integer :: i, n, k
    integer :: num_cells_grid                  ! ni*nj, including non-land cells
    integer :: ierr_elmxx                      ! ELMxx C API status
    integer :: n_nat_col, n_nat_patch, n_urb_lun
    integer :: nfail_maps
    integer, allocatable :: land_ids(:)        ! global grid IDs of the active land cells
    character(len=*), parameter :: subname = '(elmxx_init) '

    if (masterproc) then
       write(logunit,*) 'ELMxx model initialization'
       call shr_sys_flush(logunit)
    end if

    if (len_trim(fatmlndfrc) == 0) then
       call shr_sys_abort(subname//' ERROR: fatmlndfrc is not set in lnd_in')
    end if

    ! ---- global grid, read on every rank ----
    call elmxx_pio_init(inst_name)
    call elmxx_read_domain(logunit, fatmlndfrc, nlon_g, nlat_g, &
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

    if (npes > num_cells_global) then
       write(logunit,*) subname,'ERROR: number of processes exceeds number of ', &
                        'land grid cells ',npes,num_cells_global
       call shr_sys_abort(subname//' ERROR: more MPI ranks than active land cells')
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

    ! ---- round-robin decomposition over the active land cells ----
    ! Rank iam owns the iam+1, iam+1+npes, iam+1+2*npes, ... land cells. Note the
    ! IDs stored are global *grid* IDs (indices into the ni*nj arrays), which is
    ! what the gsMap and the domain both need: the gsMap global size stays ni*nj
    ! so it matches the atm grid, which the coupler requires when the atm and lnd
    ! grids are the same (see seq_domain_mct.F90).
    num_cells_owned = num_cells_global / npes
    if (iam < mod(num_cells_global, npes)) num_cells_owned = num_cells_owned + 1

    allocate(natural_id_cells_owned(num_cells_owned))
    do i = 1, num_cells_owned
       natural_id_cells_owned(i) = land_ids(iam + 1 + (i-1)*npes)
    end do

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

    call elmxx_forcing_init(num_cells_owned)

    !-----------------------------------------------------------------------
    ! Surface dataset.
    !
    ! Optional for now: without it ELMxx still has a decomposition and a domain,
    ! which is all Stage 1 needed, so a case with no fsurdat stays runnable
    ! rather than aborting. It becomes mandatory once subgrid construction
    ! depends on it.
    !-----------------------------------------------------------------------
    if (len_trim(fsurdat) > 0) then
       call elmxx_read_surfdata(logunit, fsurdat, num_cells_grid, &
                                natural_id_cells_owned)
       call elmxx_report_composition(logunit)
       call elmxx_build_subgrid(logunit, num_cells_owned)
       call elmxx_surface_state_init(logunit, month, day)
       call elmxx_build_filters(logunit, num_cells_owned)
       call elmxx_write_init_snapshot(logunit, month, day, natural_id_cells_owned)
    else
       if (masterproc) then
          write(logunit,*) subname,'no fsurdat in lnd_in; no surface dataset read'
          call shr_sys_flush(logunit)
       end if
    end if

    !-----------------------------------------------------------------------
    ! Bring up the Kokkos runtime and create the ELMxx model object.
    !-----------------------------------------------------------------------
    call ELMxxKokkosInitialize()

    if (masterproc) then
       write(logunit,*) subname,'Kokkos configuration:'
       call shr_sys_flush(logunit)
       call ELMxxKokkosPrintConfiguration()
       call shr_sys_flush(logunit)
    end if

    ! Counts come from the subgrid once there is one. Without a surface dataset
    ! there is no subgrid, so fall back to the Stage 1 placeholder rather than
    ! failing -- a domain-only case stays runnable.
    if (subgrid_built) then
       call elmxx_count_for_create(n_nat_col, n_nat_patch, n_urb_lun)
    else
       n_nat_col   = num_cells_owned
       n_nat_patch = num_cells_owned
       n_urb_lun   = 0
    end if

    call ELMxxCreate(n_nat_col, n_nat_patch, n_urb_lun, elmxx_state, ierr_elmxx)
    if (ierr_elmxx /= ELMXX_SUCCESS) then
       write(logunit,*) subname,'ELMxxCreate failed with status ',ierr_elmxx
       call shr_sys_abort(subname//' ERROR: ELMxxCreate failed')
    end if
    elmxx_state_created = .true.

    write(logunit,*) subname,'rank ',iam,' created ELMxx object: natural ', &
                     n_nat_col,' columns ',n_nat_patch,' patches, urban ', &
                     n_urb_lun,' landunits'
    call shr_sys_flush(logunit)

    !-----------------------------------------------------------------------
    ! Stage 3: build the packed Fortran<->Kokkos maps and push the one piece
    ! of topology the kernels dereference.
    !
    ! The maps are built by their own pass over the subgrid, so their extents
    ! are an INDEPENDENT count from the one ELMxxCreate was given. Requiring
    ! the two to agree is the check that matters: a setter whose length
    ! disagrees with the allocated view is rejected and leaves that view at
    ! zero (STATUS.md E.1), which is silent, so it must be impossible by
    ! construction rather than caught downstream.
    !-----------------------------------------------------------------------
    if (subgrid_built) then
       call elmxx_kokkos_state_init(logunit)

       if (n_kokkos_col /= n_nat_col .or. n_kokkos_patch /= n_nat_patch .or. &
           n_kokkos_urb /= n_urb_lun) then
          write(logunit,*) subname,'ERROR: packed map extents ',n_kokkos_col, &
               n_kokkos_patch, n_kokkos_urb,' disagree with ELMxxCreate ', &
               n_nat_col, n_nat_patch, n_urb_lun
          call shr_sys_abort(subname//' ERROR: packed map extents disagree with ELMxxCreate')
       end if

       ! Semantics before plumbing: the invariant check needs no Kokkos, so
       ! run it before anything is pushed across the boundary. A map that is
       ! wrong about which entity is which should fail here, not survive to be
       ! round-tripped consistently by elmxx_kokkos_verify_maps at step 1.
       call elmxx_kokkos_check_map_invariants(logunit, nfail_maps)
       if (nfail_maps /= 0) then
          call shr_sys_abort(subname//' ERROR: packed map invariants violated')
       end if

       call elmxx_kokkos_seed_topology(elmxx_state, logunit)
       call elmxx_kokkos_seed_state(elmxx_state, logunit)
       call elmxx_kokkos_seed_canopy_hydrology(elmxx_state, logunit)
       call elmxx_kokkos_seed_albedo(elmxx_state, logunit)

       ! Constant for the whole run, so seeded once here rather than per step.
       if (elmxx_stomata_closed) then
          call elmxx_kokkos_seed_stomata_closed(elmxx_state, logunit)
       end if

       ! Soil hydraulic properties: derived from surfdata texture here, then
       ! pushed. Eight kernels read these, and none of them can run until
       ! they are real rather than zero.
       call elmxx_soil_prop_init(logunit)
       call elmxx_kokkos_seed_soil_properties(elmxx_state, logunit)

       ! Root profile. Needs the PFT parameters and the soil grid, so it comes
       ! after both. btran itself is per-step and is computed in elmxx_run.
       if (len_trim(fparamfile) > 0) then
          call elmxx_read_pftcon(logunit, fparamfile)
          call elmxx_kokkos_seed_pftpar(elmxx_state, logunit)
          call elmxx_root_init(logunit)
       end if

       ! Photosynthesis. Needs the PFT parameters, and needs SurfaceAlbedo to
       ! be running because vcmaxcint is its output -- both checked here so a
       ! misconfiguration names itself instead of surfacing as a zero
       ! stomatal resistance, which is a legal number and looks like an
       ! answer.
       if (elmxx_do_photosynthesis) then
          if (len_trim(fparamfile) == 0) then
             call shr_sys_abort(subname//' ERROR: elmxx_do_photosynthesis '// &
                  'requires fparamfile')
          end if
          if (.not. elmxx_do_albedo) then
             call shr_sys_abort(subname//' ERROR: elmxx_do_photosynthesis '// &
                  'requires elmxx_do_albedo -- vcmaxcintsun/sha are '// &
                  'SurfaceAlbedo outputs and nothing else computes them')
          end if
          ! The rest waits for the first step: elmxx_photosyn_init needs
          ! dtime to size the 10-day running mean, and dtime is an argument to
          ! elmxx_run, not to init. Same reason the soil kernel surface is
          ! built on step one.
       end if

       ! Parse after seeding, so a blocked kernel's abort names a
       ! prerequisite that genuinely could not be met rather than one that
       ! merely had not been met yet at this point in init.
       call elmxx_kernels_parse(elmxx_kernels, logunit)

       ! The five soil/hydrology kernels' surface is NOT built here. It needs
       ! the coupling timestep, and elmxx_init does not have it -- dt is an
       ! argument to elmxx_run. Built on the first step instead; see there.
    end if

    if (masterproc) then
       write(logunit,*) 'ELMxx model initialization completed'
       call shr_sys_flush(logunit)
    end if


    !-----------------------------------------------------------------------
    ! Diagnostic trace. Writes ELM's own ELMDIAG1 format so a free-running
    ! ELMxx run can be diffed against elm_diagnostics.bin with the existing
    ! tooling. Off unless ELMXX_DIAG is set in the environment.
    !-----------------------------------------------------------------------
    block
      character(len=256) :: diag_path
      integer :: dlen
      call get_environment_variable('ELMXX_DIAG', diag_path, dlen)
      if (dlen > 0) then
         call elmxx_diag_init(trim(diag_path), .true.)
         call elmxx_diag_write_maps()
         if (masterproc) write(logunit,*) 'ELMxx: diagnostics -> ',trim(diag_path)
      end if
    end block

  end subroutine elmxx_init

  !-----------------------------------------------------------------------
  subroutine elmxx_count_for_create(n_nat_col, n_nat_patch, n_urb_lun)
    !
    ! !DESCRIPTION:
    ! Translate the subgrid into the three counts ELMxxCreate wants.
    !
    ! ELMxx's C++ side is organised per surface type, not as one flat begc:endc
    ! span, so it needs the natural-vegetation columns and patches and the urban
    ! landunit count -- not the totals. Lake and glacier are carried by the
    ! subgrid but are not part of this call: lake state is allocated separately
    ! (ELMxxAllocateLakeState) and glacier has no kernels at all.
    !
    implicit none
    !
    integer, intent(out) :: n_nat_col, n_nat_patch, n_urb_lun
    !
    integer :: l, c

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
    do c = 1, num_columns
       if (lun_itype(col_landunit(c)) == istsoil) then
          n_nat_patch = n_nat_patch + count_patches_on(c)
       end if
    end do

    ! ELMxxCreate rejects zero natural columns or patches. Every gridcell gets a
    ! natural-vegetation landunit, so this can only happen with no cells at all,
    ! which elmxx_init already rules out.
    if (n_nat_col <= 0 .or. n_nat_patch <= 0) then
       call shr_sys_abort('(elmxx_count_for_create) ERROR: no natural columns or patches')
    end if

  end subroutine elmxx_count_for_create

  !-----------------------------------------------------------------------
  integer function count_patches_on(c)
    !
    use elmxxSubgridMod, only : num_patches, patch_column
    implicit none
    integer, intent(in) :: c
    integer :: p

    count_patches_on = 0
    do p = 1, num_patches
       if (patch_column(p) == c) count_patches_on = count_patches_on + 1
    end do

  end function count_patches_on

  !-----------------------------------------------------------------------
  subroutine elmxx_report_composition(logunit)
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
    character(len=*), parameter :: subname = '(elmxx_report_composition) '

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
                        'which ELMxx does not model; it is not yet excluded ', &
                        'from the subgrid'
       call shr_sys_flush(logunit)
    end if

  end subroutine elmxx_report_composition

  !-----------------------------------------------------------------------
  subroutine elmxx_run(logunit, coupling_dt_in_sec, month, day, &
                       nextsw_cday, declinp1, doalb)
    !
    ! !DESCRIPTION:
    ! Advance ELMxx one coupling interval.
    !
    ! This is intentionally a no-op: the Kokkos port is not wired in yet, so
    ! ELMxx consumes nothing from the coupler and returns zeroed fields.
    !
    implicit none
    !
    integer, intent(in) :: logunit
    integer, intent(in) :: coupling_dt_in_sec
    real(r8), intent(in) :: nextsw_cday   ! calendar day of the next radiation step
    real(r8), intent(in) :: declinp1      ! solar declination for it, radians
    integer, intent(in) :: month, day
    logical, intent(in), optional :: doalb  ! .false. on ELM's nstep-0 pass

    logical :: do_albedo_this_step
    logical :: doalb_in            ! the driver's doalb, independent of config
    integer :: ierr_rs
    integer :: phm1, phm2
    real(r8) :: phw1, phw2

    doalb_in = .true.
    if (present(doalb)) doalb_in = doalb
    do_albedo_this_step = elmxx_do_albedo .and. doalb_in

    nstep = nstep + 1
    call elmxx_diag_new_timestep(nstep)

    ! ELM gates SatellitePhenology on doalb (elm_driver.F90, the non-CN,
    ! non-FATES branch), and does NOT call it during initialisation for this
    ! configuration -- the initialize2 call is behind use_fates .and.
    ! use_fates_sp. So ELM carries elai = esai = frac_veg_nosno = 0 until the
    ! first doalb step, which is nstep 2 here, and treats every patch as bare
    ! ground until then. ELMxx updated phenology unconditionally and so entered
    ! step 0 with a full canopy, routing patches through CanopyFluxes where ELM
    ! was still running BareGroundFluxes.
    if (subgrid_built) then
       call elmxx_write_init_snapshot(logunit, month, day, natural_id_cells_owned)
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
       call elmxx_kokkos_push_forcing(elmxx_state, logunit)
    end if

    !-----------------------------------------------------------------------
    ! Root water stress. btran depends on soil moisture and temperature, so it
    ! is recomputed here rather than seeded once -- it would go stale the
    ! moment hydrology starts evolving the soil column.
    !-----------------------------------------------------------------------
    if (root_built) then
       ! Root water stress now runs on the device. The statics it needs are
       ! pushed once; nothing crosses per step. This deletes the push_btran
       ! crossing -- the one that silently dropped rootr and cost a month of
       ! water.
       if (.not. root_statics_pushed) then
          call elmxx_kokkos_push_root_statics(elmxx_state, logunit)
          root_statics_pushed = .true.
       end if
       call ELMxxComputeRootStressNatural(elmxx_state, ierr_rs)
       if (ierr_rs /= ELMXX_SUCCESS) &
            call shr_sys_abort('(elmxx_run) ERROR: ComputeRootStressNatural failed')
    end if

    !-----------------------------------------------------------------------
    ! Stage 4: run the active kernels, in driver order. After forcing has
    ! crossed, before the boundary probe -- the probe overwrites state fields
    ! with fingerprints, so it has to come last in the step.
    !-----------------------------------------------------------------------
    ! The five soil/hydrology kernels' surface, built on the first step
    ! because it needs dtime.
    !
    ! ELMxxInitSharedMetadata IS NOT OPTIONAL AND IS NOT ONLY ABOUT FILTERS.
    ! It is the only thing in the whole C API that assigns elm->dtime, and
    ! BuildSoilTemperatureNaturalView reads exactly that -- so skipping it
    ! would leave the soil column integrating on whatever dtime the object was
    ! constructed with. The canopy kernels do not expose this, since
    ! ELMxxComputeCanopyHydrology takes dtime as an explicit argument.
    if (kokkos_state_built .and. .not. soil_kernel_built) then
       if (kernel_active(K_SOILTEMP)   .or. kernel_active(K_SOILFLUX)  .or. &
           kernel_active(K_SURFRUNOFF) .or. kernel_active(K_ROOTWATER) .or. &
           kernel_active(K_HYDRODRAIN)) then
          call elmxx_soil_kernel_init(elmxx_state, &
               real(coupling_dt_in_sec, r8), logunit)
       end if
    end if

    if (kokkos_state_built .and. any_kernel_active) then
       ! The ground surface energy balance has to be formed AFTER the canopy
       ! and radiation kernels of this step have run and BEFORE
       ! SoilTemperature consumes it -- but elmxx_kernels_run dispatches the
       ! whole ordered set in one call. So the canopy half runs first, then
       ! this, then the soil half. Splitting the dispatch is what makes the
       ! ordering explicit rather than implicit in a comment.
       ! PHOTOSYNTHESIS INPUTS GO HERE, before the first kernel block, because
       ! canflux is in it and canflux is what consumes them. Placing this after
       ! the kernels -- next to SurfaceAlbedo, where it superficially belongs
       ! with the other end-of-step work -- would hand the canopy the previous
       ! step's daylength and a t10 that had not seen this step's temperature.
       !
       ! vcmaxcint is SurfaceAlbedo's output and SurfaceAlbedo runs at the END
       ! of the step, so on step one it does not exist yet. That is the same
       ! temporal dependency the albedos themselves have, and ELM has it too --
       ! its SurfaceAlbedo runs in initialize2 so step one has a value. ELMxx
       ! has no such call, so step one runs with vcmaxcint = 0, which is one
       ! step of no photosynthesis and is stated rather than hidden.
       if (elmxx_do_photosynthesis) then
          if (.not. photosyn_built) then
             call elmxx_photosyn_init(cell_lat, real(coupling_dt_in_sec, r8), logunit)
             call elmxx_photosyn_seed(elmxx_state, logunit)
          end if
          ! Photosynthesis forcing runs on the device. Its statics -- gridcell
          ! latitude and maximum daylength -- are pushed once. vcmaxcint used
          ! to cross here too; SurfaceAlbedo now produces it on the device,
          ! which is what let this move.
          if (.not. photosyn_statics_pushed) then
             call elmxx_push_photosyn_statics(elmxx_state, logunit)
             photosyn_statics_pushed = .true.
          end if
          call ELMxxComputePhotosynForcingNatural(elmxx_state, nstep, t10_period, &
               declinp1, elmxx_co2_ppmv, ierr_rs)
          if (ierr_rs /= ELMXX_SUCCESS) &
               call shr_sys_abort('(elmxx_run) ERROR: ComputePhotosynForcingNatural failed')
       end if

       ! State at the top of the step, before any kernel. ELM's matching
       ! anchor is canhydro_in: -- its first kernel -- so these line up.
       call elmxx_diag_snapshot_state(elmxx_state, nlevtot, nlevgrnd, 'elmxx_in')

       call elmxx_kernels_run(elmxx_state, real(coupling_dt_in_sec, r8), logunit, 1)

       if (soil_kernel_built) then
          ! Ground surface energy balance and the three quantities that used
          ! to be pushed with it now run on the device. This deletes the
          ! soil_kernel_push crossing.
          if (.not. ghf_sb_pushed) then
             call ELMxxSetGroundHeatFluxSb(elmxx_state, SHR_CONST_STEBOL, ierr_rs)
             ghf_sb_pushed = .true.
          end if
          call ELMxxComputeGroundHeatFluxNatural(elmxx_state, ierr_rs)
          if (ierr_rs /= ELMXX_SUCCESS) &
               call shr_sys_abort('(elmxx_run) ERROR: ComputeGroundHeatFluxNatural failed')
       end if

       call elmxx_kernels_run(elmxx_state, real(coupling_dt_in_sec, r8), logunit, 2)

       call elmxx_diag_snapshot_fluxes(elmxx_state, 'elmxx_out')

       ! Read the wetted soil column back, so next step's Fortran-side btran
       ! sees what the hydrology kernels just did rather than the cold start.
       if (soil_kernel_built) then
          ! soil_kernel_pull is gone. Its only per-step consumer was the
          ! Fortran btran, which now runs on the device; the host-side
          ! col_h2osoi_* arrays are used at INIT for seeding and nowhere else,
          ! so nothing needs the soil column brought back every step.
       end if

       ! SURFACE ALBEDO RUNS HERE, AT THE END OF THE STEP, BECAUSE ELM RUNS IT
       ! HERE (elm_driver.F90, gated on doalb). SurfaceRadiation therefore
       ! never consumes an albedo computed in its own step -- it reads the
       ! previous step's, and step one reads the InitCold constants
       ! elmxx_kokkos_seed_albedo supplied. Moving this to the top of the step
       ! would also change which state the two-stream sees: t_veg, fwet and
       ! h2osoi_vol have all been updated by now.
       ! Phenology updates HERE, at the end of the step just before
       ! SurfaceAlbedo, because that is where ELM does it (elm_driver.F90, the
       ! non-CN non-FATES branch, gated on doalb -- immediately ahead of its
       ! own SurfaceAlbedo call). Updating at the top of the step instead let
       ! CanopyHydrology see leaf area a step before ELM's did, which showed up
       ! as canopy water running one timestep ahead all run.
       ! Phenology interpolates ON THE DEVICE. The twelve monthly fields are
       ! pushed once from the surface dataset; per step only the two month
       ! indices and their weights cross. This replaces push_phenology, which
       ! shipped interpolated leaf area every doalb step.
       if (subgrid_built .and. doalb_in) then
          if (.not. monthly_phen_pushed) then
             call elmxx_push_monthly_phenology(elmxx_state, logunit)
             monthly_phen_pushed = .true.
          end if
          call elmxx_phenology_weights(month, day, phm1, phm2, phw1, phw2)
          call ELMxxComputePhenologyNatural(elmxx_state, phm1, phm2, phw1, phw2, ierr_rs)
          if (ierr_rs /= ELMXX_SUCCESS) &
               call shr_sys_abort('(elmxx_run) ERROR: ComputePhenologyNatural failed')
       end if

       if (do_albedo_this_step) then
          call elmxx_push_coszen(elmxx_state, nextsw_cday, declinp1, &
               cell_lat, cell_lon, logunit)
          call ELMxxComputeSurfaceAlbedoNatural(elmxx_state, ierr_rs)
          if (ierr_rs /= ELMXX_SUCCESS) &
               call shr_sys_abort('(elmxx_run) ERROR: ComputeSurfaceAlbedoNatural failed')
          if (nstep == 1 .or. mod(nstep, 24) == 0) then
             call elmxx_surface_albedo_report(logunit)
          end if
       end if
       ! Report on the first step, then TWICE daily -- not once. The extra
       ! sample is what makes the radiation readable: step 48k lands at model
       ! midnight, where incident shortwave is zero and every SurfaceRadiation
       ! output is trivially zero, so a once-daily report would show a
       ! radiation kernel that appears to do nothing. Sampling half a day out
       ! catches daylight. It also gives the surface fluxes a day/night
       ! contrast, which is worth having for free.
       if (nstep == 1 .or. mod(nstep, 24) == 0) then
          call elmxx_report_surfrad(elmxx_state, n_kokkos_patch, logunit)
          call elmxx_kernels_report(elmxx_state, n_kokkos_patch, logunit)
          call elmxx_report_cantemp(elmxx_state, n_kokkos_col, logunit)
          call elmxx_report_fluxes(elmxx_state, n_kokkos_patch, logunit)
       end if
    end if

    ! Gated on the FIRST pass, not nstep == 1. The probe overwrites state with
    ! fingerprints and restores from the Fortran-side seed, so it is only
    ! harmless while the state still IS that seed. Once ELMxx started running
    ! ELM's extra nstep-0 pass, nstep == 1 became the second pass and this
    ! silently threw away a step of physics.
    if (kokkos_state_built .and. nstep == 0 .and. elmxx_check_boundary) then
       call elmxx_verify_kokkos_boundary(logunit)
    end if

    if (masterproc) then
       write(logunit,*) 'ELMxx step ',nstep,' dt = ',coupling_dt_in_sec,' s (no-op)'
       call shr_sys_flush(logunit)
    end if

  end subroutine elmxx_run

  !-----------------------------------------------------------------------
  subroutine elmxx_init_albedo(logunit, nextsw_cday, declinp1)
    !
    ! One SurfaceAlbedo pass at the end of initialisation, mirroring ELM's
    ! initialize2. Without it ELMxx enters step 0 with no albedo and no
    ! vcmaxcint, while ELM enters with both -- and because doalb is false on
    ! ELM's own nstep 0 and 1 passes, nothing would fill them in until step 2.
    !
    implicit none
    integer , intent(in) :: logunit
    real(r8), intent(in) :: nextsw_cday, declinp1
    integer :: ierr_ia

    if (.not. kokkos_state_built) return
    if (.not. elmxx_do_albedo) return

    ! Uses the same device kernel as the run loop -- the Fortran path read a
    ! soil moisture that is never updated (H11), so it has no business
    ! setting the initial albedo either.
    call elmxx_push_coszen(elmxx_state, nextsw_cday, declinp1, &
         cell_lat, cell_lon, logunit)
    call ELMxxComputeSurfaceAlbedoNatural(elmxx_state, ierr_ia)
    if (ierr_ia /= ELMXX_SUCCESS) &
         call shr_sys_abort('(elmxx_init_albedo) ERROR: ComputeSurfaceAlbedoNatural failed')
    if (masterproc) then
       write(logunit,*) 'ELMxx: initial SurfaceAlbedo pass complete'
       call shr_sys_flush(logunit)
    end if

  end subroutine elmxx_init_albedo


  !-----------------------------------------------------------------------
  subroutine elmxx_verify_kokkos_boundary(logunit)
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
    character(len=*), parameter :: subname = '(elmxx_verify_kokkos_boundary) '

    call elmxx_kokkos_verify_maps(elmxx_state, logunit, nfail)

    if (nfail /= 0) then
       if (elmxx_check_soft_fail) then
          write(logunit,*) subname,'WARNING: ',nfail,' boundary mismatches; ', &
               'continuing because elmxx_check_soft_fail is set. Every kernel ', &
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
    call elmxx_kokkos_seed_state(elmxx_state, logunit)

  end subroutine elmxx_verify_kokkos_boundary

  !-----------------------------------------------------------------------
  subroutine elmxx_final()
    !
    ! !DESCRIPTION:
    ! Finalize ELMxx.
    !
    implicit none

    integer :: ierr_elmxx

    !-----------------------------------------------------------------------
    ! Tear down in reverse order of elmxx_init: object first, then Kokkos.
    ! ELMxxKokkosFinalize must come after ELMxxDestroy -- the object owns
    ! Kokkos views, and destroying them after the runtime is gone is undefined.
    !-----------------------------------------------------------------------
    if (elmxx_state_created) then
       call ELMxxDestroy(elmxx_state, ierr_elmxx)
       if (ierr_elmxx /= ELMXX_SUCCESS) then
          write(iulog,*) 'elmxx_final: ELMxxDestroy failed with status ',ierr_elmxx
       end if
       elmxx_state_created = .false.
       call elmxx_soil_kernel_clean()
       call elmxx_kokkos_state_clean()
       call elmxx_soil_prop_clean()
       call elmxx_root_clean()
       call elmxx_pftcon_clean()
       call ELMxxKokkosFinalize()
    end if

    call elmxx_forcing_clean()
    call elmxx_surface_state_clean()
    call elmxx_filters_clean()
    call elmxx_subgrid_clean()
    call elmxx_surfdata_clean()

    if (associated(natural_id_cells_owned)) deallocate(natural_id_cells_owned)
    if (associated(lonc_g))  deallocate(lonc_g)
    if (associated(latc_g))  deallocate(latc_g)
    if (associated(areac_g)) deallocate(areac_g)
    if (associated(maskc_g)) deallocate(maskc_g)
    if (associated(fracc_g)) deallocate(fracc_g)

    call elmxx_diag_finalize()

  end subroutine elmxx_final

end module elmxxMod
