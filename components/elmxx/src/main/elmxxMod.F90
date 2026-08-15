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
                                   elmxx_update_phenology, elmxx_surface_state_clean
  use elmxxFilterMod      , only : elmxx_build_filters, elmxx_filters_clean
  use elmxxInitCheckMod   , only : elmxx_write_init_snapshot
  use elmxxForcingMod , only : elmxx_forcing_init, elmxx_forcing_clean

  use elmxx_mod              , only : ELMxxType, ELMxxCreate, ELMxxDestroy, ELMXX_SUCCESS
  use elmxxKokkosStateMod    , only : elmxx_kokkos_state_init, &
                                      elmxx_kokkos_check_map_invariants, &
                                      elmxx_kokkos_seed_topology, &
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

  !--------------------------------------------------------------------------
  ! Namelist (&elmxx_inparm in lnd_in)
  !--------------------------------------------------------------------------
  logical           , public :: do_elmxx   = .true.
  character(len=256), public :: fatmlndfrc = ' '
  character(len=256), public :: fsurdat    = ' '

  !--------------------------------------------------------------------------
  ! Instance information
  !--------------------------------------------------------------------------
  character(len=16), public :: inst_name
  character(len=16), public :: inst_suffix   ! e.g. "_0001" or ""
  integer          , public :: inst_index

  integer, public :: iulog = 6

  integer, private :: nstep = 0

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

    namelist /elmxx_inparm/ do_elmxx, fatmlndfrc, fsurdat

    ! defaults
    do_elmxx   = .true.
    fatmlndfrc = ' '
    fsurdat    = ' '

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

    if (masterproc) then
       write(logunit,*) ' '
       write(logunit,*) 'read from namelist:'
       write(logunit,*) '   do_elmxx   = ', do_elmxx
       write(logunit,*) '   fatmlndfrc = ', trim(fatmlndfrc)
       write(logunit,*) '   fsurdat    = ', trim(fsurdat)
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

    nstep = 0

    ! Per-timestep atmospheric forcing lives at gridcell level, so it is sized
    ! from the decomposition and does not depend on the surface dataset.
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
    end if

    if (masterproc) then
       write(logunit,*) 'ELMxx model initialization completed'
       call shr_sys_flush(logunit)
    end if

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
  subroutine elmxx_run(logunit, coupling_dt_in_sec, month, day)
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
    integer, intent(in) :: month, day

    nstep = nstep + 1

    if (subgrid_built) then
       call elmxx_update_phenology(logunit, month, day)
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
    if (kokkos_state_built .and. nstep == 1) then
       call elmxx_verify_kokkos_boundary(logunit)
    end if

    if (masterproc) then
       write(logunit,*) 'ELMxx step ',nstep,' dt = ',coupling_dt_in_sec,' s (no-op)'
       call shr_sys_flush(logunit)
    end if

  end subroutine elmxx_run

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
       call shr_sys_abort(subname//' ERROR: packed Fortran<->Kokkos maps do not round-trip')
    end if

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
       call elmxx_kokkos_state_clean()
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

  end subroutine elmxx_final

end module elmxxMod
