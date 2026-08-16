module elmxxKokkosStateMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Stage 3's persistent-state handshake: the boundary between ELMxx's Fortran
  ! subgrid and its Kokkos Views.
  !
  ! THREE LAYERS, ONE MAP (plans/STATUS.md F):
  !
  !   1. Kokkos Views        packed per surface type -- naturalCol(0:nkc-1),
  !                          naturalPatch(0:nkp-1), urbanLandunit(0:nku-1).
  !                          Lake, wetland and glacier are absent entirely.
  !   2. Fortran buffers     THIS module. Dimensioned over the whole subgrid
  !                          (num_columns, num_patches, num_landunits) in
  !                          ELMxx's construction order, which follows ELM's
  !                          begc:endc convention -- clump, landunit type,
  !                          then gridcell.
  !   3. history / restart   written from layer 2, permuted to ELM's on-disk
  !      / l2x               gridcell-major order at write time, the way ELM's
  !                          gsMap does it. Not this module's job.
  !
  ! The maps below are what the plan calls "the ELM<->ELMxx index maps,
  ! persisted in the ELMxx object". They are built once, from the subgrid's own
  ! type information, and never from loop-position coincidence -- that is the
  ! URBANxx bug (urbanxx_netShortwave_check reading through the wrong filter)
  ! this whole stage exists to prevent.
  !
  ! TWO CHECKS, NEITHER SUFFICIENT ALONE:
  !   elmxx_kokkos_check_map_invariants  semantics, no Kokkos involved --
  !                                      bijection, types, patch->column
  !                                      agreement with the subgrid.
  !   elmxx_kokkos_verify_maps           plumbing, through Kokkos -- extents,
  !                                      addressability, no silent rejection.
  ! A consistently wrong map passes the second and fails the first.
  !
  ! WHY A SEPARATE BUFFER LAYER AT ALL. The Kokkos side is packed and skips
  ! landunit types it has no kernels for; the Fortran side must keep every
  ! entity, because the Stage 2 comparison against ELM depends on carrying
  ! zero-weight and out-of-scope landunits. Aliasing one onto the other is only
  ! correct when every gridcell has exactly one landunit, which is the
  ! degenerate single-cell case. Hence: map, do not alias.
  !-----------------------------------------------------------------------

  use shr_kind_mod    , only : r8 => shr_kind_r8
  use shr_sys_mod     , only : shr_sys_abort, shr_sys_flush
  use shr_const_mod   , only : SHR_CONST_RDAIR
  use elmxxSpmdMod    , only : masterproc, iam
  use elmxxSubgridMod , only : num_landunits, num_columns, num_patches, &
                               lun_gridcell, lun_itype, col_landunit, &
                               col_itype, patch_column, patch_itype, &
                               istsoil, isturb_tbd, isturb_hd, isturb_md
  use elmxxSurfaceStateMod , only : surface_state_built, patch_lai, patch_sai, &
                                    patch_height_top
  use elmxxSurfdataMod, only : topo_std, topo_slope
  use elmxx_kokkos_interface, only : ELMxxKokkosIsLayoutRight
  use elmxxRootMod    , only : rt_btran => btran
  use elmxxSoilPropMod, only : soil_prop_built, nlevgrnd, nlevtot, &
                               sp_watsat => watsat, sp_bsw => bsw, &
                               sp_sucsat => sucsat, sp_watfc => watfc, &
                               sp_dz => col_dz, sp_tsoisno => col_t_soisno, &
                               sp_liq => col_h2osoi_liq, sp_ice => col_h2osoi_ice
  use elmxxPftconMod  , only : z0mr, displar, npft_param, pftcon_read
  use elmxxForcingMod , only : forc_z, forc_u, forc_v, forc_ptem, forc_shum, forc_pbot, &
                               forc_tbot, forc_lwrad, forc_rainc, forc_rainl, &
                               forc_snowc, forc_snowl, &
                               forc_swvdr, forc_swndr, forc_swvdf, forc_swndf
  use elmxx_mod       , only : ELMxxType, ELMXX_SUCCESS, &
                               ELMxxSetPatchColumn, &
                               ELMxxSetElai, ELMxxSetEsai, ELMxxSetHtop, &
                               ELMxxSetForcTCol, ELMxxSetForcPbotCol, &
                               ELMxxSetForcQCol, ELMxxSetForcLwradCol, &
                               ELMxxSetForcUCol, ELMxxSetForcVCol, &
                               ELMxxSetForcThCol, ELMxxSetForcT, &
                               ELMxxSetForcRain, ELMxxSetForcSnow, &
                               ELMxxSetDewmx, ELMxxSetMicroSigma, ELMxxSetNMelt, &
                               ELMxxSetFracVegNosno, ELMxxSetFwet, ELMxxSetH2ocan, &
                               ELMxxSetH2osfc, ELMxxSetIntSnow, &
                               ELMxxSetFracH2osfc, ELMxxSetFracSnoEff, &
                               ELMxxSetDoCapsnow, &
                               ELMxxSetWatsat, ELMxxSetBsw, ELMxxSetSucsat, &
                               ELMxxSetWatfc, ELMxxSetDz, ELMxxSetTSoisno, &
                               ELMxxSetH2osoiLiq, ELMxxSetH2osoiIce, &
                               ELMxxSetSmpmin, ELMxxSetTH2osfc, &
                               ELMxxSetPatchItype, ELMxxSetForcHgtPatch, &
                               ELMxxSetForcSolad, ELMxxSetForcSolai, &
                               ELMxxSetForcRhoCol, ELMxxSetBtran, &
                               ELMxxSetAlbgrd, ELMxxSetAlbgri, &
                               ELMxxSetAlbsod, ELMxxSetAlbsoi, &
                               ELMxxSetAlbsndHst, ELMxxSetAlbsniHst, &
                               ELMxxSetAlbd, ELMxxSetAlbi, &
                               ELMxxSetFabd, ELMxxSetFabi, &
                               ELMxxSetFtdd, ELMxxSetFtid, ELMxxSetFtii, &
                               ELMxxSetFlxAbsdv, ELMxxSetFlxAbsdn, &
                               ELMxxSetFlxAbsiv, ELMxxSetFlxAbsin, &
                               ELMxxSetSnl        , ELMxxGetSnl, &
                               ELMxxSetSnowDepth  , ELMxxGetSnowDepth, &
                               ELMxxSetFracSno    , ELMxxGetFracSno, &
                               ELMxxSetH2osno     , ELMxxGetH2osno, &
                               ELMxxSetTGrnd      , ELMxxGetTGrnd, &
                               ELMxxSetTVeg       , ELMxxGetTVeg, &
                               ELMxxSetFsun       , ELMxxGetFsun, &
                               ELMxxSetUrbanTaf   , ELMxxGetUrbanTaf, &
                               ELMxxSetUrbanQaf   , ELMxxGetUrbanQaf, &
                               ELMxxSetZ0mrPft    , ELMxxSetDisplarPft

  implicit none
  save
  private

  !--------------------------------------------------------------------------
  ! The maps. Forward maps are 0-based, because that is what the Kokkos side
  ! indexes with (see patch.patch_column(p) in BareGroundFluxesImpl.h); a -1
  ! marks an entity ELMxx's C++ side does not carry. Reverse maps are 1-based
  ! Fortran indices.
  !--------------------------------------------------------------------------
  integer, public, pointer :: kcol_of_col(:)     => null()  ! (num_columns)   -> 0-based, or -1
  integer, public, pointer :: col_of_kcol(:)     => null()  ! (n_nat_col)     -> 1-based column
  integer, public, pointer :: kpatch_of_patch(:) => null()  ! (num_patches)   -> 0-based, or -1
  integer, public, pointer :: patch_of_kpatch(:) => null()  ! (n_nat_patch)   -> 1-based patch
  integer, public, pointer :: kurb_of_lun(:)     => null()  ! (num_landunits) -> 0-based, or -1
  integer, public, pointer :: lun_of_kurb(:)     => null()  ! (n_urb_lun)     -> 1-based landunit

  integer, public :: n_kokkos_col   = 0
  integer, public :: n_kokkos_patch = 0
  integer, public :: n_kokkos_urb   = 0

  logical, public :: kokkos_state_built = .false.

  public :: elmxx_kokkos_state_init
  public :: elmxx_kokkos_check_map_invariants
  public :: elmxx_kokkos_seed_topology
  public :: elmxx_kokkos_seed_state
  public :: elmxx_kokkos_seed_canopy_hydrology
  public :: elmxx_kokkos_seed_soil_properties
  public :: elmxx_kokkos_seed_albedo
  public :: elmxx_kokkos_seed_pftpar
  public :: elmxx_kokkos_push_forcing
  public :: elmxx_kokkos_push_btran
  public :: elmxx_kokkos_verify_maps
  public :: elmxx_kokkos_state_clean

contains

  !-----------------------------------------------------------------------
  subroutine elmxx_kokkos_state_init(logunit)
    !
    ! Build the packed maps from the subgrid's landunit types. Must agree,
    ! entity for entity, with elmxx_count_for_create in elmxxMod.F90 -- the
    ! counts it produced are the extents ELMxxCreate allocated, and a setter
    ! whose length disagrees is rejected with ELMXX_ERR_SIZE_MISMATCH and
    ! LEAVES THE VIEW UNTOUCHED AT ZERO (STATUS.md E.1). Checked below rather
    ! than assumed.
    !
    implicit none
    integer, intent(in) :: logunit
    integer :: l, c, p, kc, kp, ku
    character(len=*), parameter :: subname = '(elmxx_kokkos_state_init) '

    if (num_columns <= 0 .or. num_patches <= 0 .or. num_landunits <= 0) then
       call shr_sys_abort(subname//'ERROR: subgrid is not ready')
    end if

    call elmxx_kokkos_state_clean()

    allocate(kcol_of_col(num_columns), kpatch_of_patch(num_patches), &
             kurb_of_lun(num_landunits))
    kcol_of_col = -1; kpatch_of_patch = -1; kurb_of_lun = -1

    ! ---- pass 1: count, in the same order the packed views are laid out ----
    kc = 0; kp = 0; ku = 0
    do c = 1, num_columns
       if (lun_itype(col_landunit(c)) == istsoil) kc = kc + 1
    end do
    do p = 1, num_patches
       if (lun_itype(col_landunit(patch_column(p))) == istsoil) kp = kp + 1
    end do
    do l = 1, num_landunits
       if (is_urban(lun_itype(l))) ku = ku + 1
    end do

    n_kokkos_col = kc; n_kokkos_patch = kp; n_kokkos_urb = ku
    allocate(col_of_kcol(max(kc,1)), patch_of_kpatch(max(kp,1)), &
             lun_of_kurb(max(ku,1)))
    col_of_kcol = 0; patch_of_kpatch = 0; lun_of_kurb = 0

    ! ---- pass 2: assign, ascending in Fortran subgrid order ----
    ! Ascending order is the contract: it is what makes the packed index a
    ! stable function of the subgrid rather than of iteration accident, and it
    ! is what lets the reverse map be a plain gather.
    kc = 0
    do c = 1, num_columns
       if (lun_itype(col_landunit(c)) == istsoil) then
          kcol_of_col(c) = kc          ! 0-based for the Kokkos side
          kc = kc + 1
          col_of_kcol(kc) = c          ! 1-based Fortran column
       end if
    end do

    kp = 0
    do p = 1, num_patches
       if (lun_itype(col_landunit(patch_column(p))) == istsoil) then
          kpatch_of_patch(p) = kp
          kp = kp + 1
          patch_of_kpatch(kp) = p
       end if
    end do

    ku = 0
    do l = 1, num_landunits
       if (is_urban(lun_itype(l))) then
          kurb_of_lun(l) = ku
          ku = ku + 1
          lun_of_kurb(ku) = l
       end if
    end do

    kokkos_state_built = .true.

    write(logunit,*) subname,'rank ',iam,' packed maps: natural columns ', &
                     n_kokkos_col,' natural patches ',n_kokkos_patch, &
                     ' urban landunits ',n_kokkos_urb
    call shr_sys_flush(logunit)

  end subroutine elmxx_kokkos_state_init

  !-----------------------------------------------------------------------
  logical function is_urban(ltype)
    implicit none
    integer, intent(in) :: ltype
    is_urban = (ltype == isturb_tbd) .or. (ltype == isturb_hd) .or. &
               (ltype == isturb_md)
  end function is_urban

  !-----------------------------------------------------------------------
  subroutine elmxx_kokkos_check_map_invariants(logunit, nfail)
    !
    ! Check the maps against the subgrid, with no Kokkos in the picture.
    !
    ! This is the half that has teeth about MEANING. Every property below is
    ! checked against elmxxSubgridMod's own arrays, so a map that is internally
    ! consistent but wrong about which entity is which fails here even though
    ! it would round-trip cleanly through the Views.
    !
    implicit none
    integer, intent(in)  :: logunit
    integer, intent(out) :: nfail
    integer :: l, c, p, kc, kp, ku, kcol
    integer, allocatable :: hits(:)
    character(len=*), parameter :: subname = '(elmxx_kokkos_check_map_invariants) '

    call require_built(subname)
    nfail = 0

    ! ---- 1. forward and reverse are mutual inverses ----------------------
    do kc = 1, n_kokkos_col
       c = col_of_kcol(kc)
       if (c < 1 .or. c > num_columns) then
          call fail(logunit, nfail, subname//'col_of_kcol out of range at packed ', kc-1)
       else if (kcol_of_col(c) /= kc-1) then
          call fail(logunit, nfail, subname//'column map is not invertible at packed ', kc-1)
       end if
    end do
    do kp = 1, n_kokkos_patch
       p = patch_of_kpatch(kp)
       if (p < 1 .or. p > num_patches) then
          call fail(logunit, nfail, subname//'patch_of_kpatch out of range at packed ', kp-1)
       else if (kpatch_of_patch(p) /= kp-1) then
          call fail(logunit, nfail, subname//'patch map is not invertible at packed ', kp-1)
       end if
    end do
    do ku = 1, n_kokkos_urb
       l = lun_of_kurb(ku)
       if (l < 1 .or. l > num_landunits) then
          call fail(logunit, nfail, subname//'lun_of_kurb out of range at packed ', ku-1)
       else if (kurb_of_lun(l) /= ku-1) then
          call fail(logunit, nfail, subname//'urban map is not invertible at packed ', ku-1)
       end if
    end do

    ! ---- 2. every packed slot is claimed exactly once --------------------
    ! Catches a map that is a function but not a bijection: two Fortran
    ! entities landing on one slot silently drops one of them.
    allocate(hits(max(n_kokkos_col,1))); hits = 0
    do c = 1, num_columns
       if (kcol_of_col(c) >= 0) hits(kcol_of_col(c)+1) = hits(kcol_of_col(c)+1) + 1
    end do
    do kc = 1, n_kokkos_col
       if (hits(kc) /= 1) call fail(logunit, nfail, &
            subname//'packed column slot claimed more or less than once: ', kc-1)
    end do
    deallocate(hits)

    allocate(hits(max(n_kokkos_patch,1))); hits = 0
    do p = 1, num_patches
       if (kpatch_of_patch(p) >= 0) hits(kpatch_of_patch(p)+1) = hits(kpatch_of_patch(p)+1) + 1
    end do
    do kp = 1, n_kokkos_patch
       if (hits(kp) /= 1) call fail(logunit, nfail, &
            subname//'packed patch slot claimed more or less than once: ', kp-1)
    end do
    deallocate(hits)

    ! ---- 3. membership matches landunit type, both ways ------------------
    ! An entity is mapped IF AND ONLY IF ELMxx's C++ side carries its type.
    ! The "only if" direction is the one that matters: a lake column leaking
    ! into the natural-column pack would be handed to soil kernels.
    do c = 1, num_columns
       if ((lun_itype(col_landunit(c)) == istsoil) .neqv. (kcol_of_col(c) >= 0)) then
          call fail(logunit, nfail, subname//'column mapped against its type: ', c)
       end if
    end do
    do p = 1, num_patches
       if ((lun_itype(col_landunit(patch_column(p))) == istsoil) .neqv. &
           (kpatch_of_patch(p) >= 0)) then
          call fail(logunit, nfail, subname//'patch mapped against its type: ', p)
       end if
    end do
    do l = 1, num_landunits
       if (is_urban(lun_itype(l)) .neqv. (kurb_of_lun(l) >= 0)) then
          call fail(logunit, nfail, subname//'landunit mapped against its type: ', l)
       end if
    end do

    ! ---- 4. patch->column agrees with the subgrid ------------------------
    ! The packed patch_column the kernels dereference must name the packed
    ! slot of the SAME column elmxxSubgridMod put the patch on.
    do kp = 1, n_kokkos_patch
       p = patch_of_kpatch(kp)
       kcol = kcol_of_col(patch_column(p))
       if (kcol < 0 .or. kcol >= n_kokkos_col) then
          call fail(logunit, nfail, subname//'packed patch has no packed column: ', kp-1)
       else if (col_of_kcol(kcol+1) /= patch_column(p)) then
          call fail(logunit, nfail, subname//'packed patch->column disagrees with subgrid: ', kp-1)
       end if
    end do

    if (nfail == 0) then
       write(logunit,*) subname,'rank ',iam,' PASSED: maps are bijective, ', &
                        'type-consistent, and agree with the subgrid'
    else
       write(logunit,*) subname,'rank ',iam,' FAILED: ',nfail,' violations'
    end if
    call shr_sys_flush(logunit)

  end subroutine elmxx_kokkos_check_map_invariants

  !-----------------------------------------------------------------------
  subroutine fail(logunit, nfail, message, idx)
    implicit none
    integer, intent(in) :: logunit, idx
    integer, intent(inout) :: nfail
    character(len=*), intent(in) :: message
    nfail = nfail + 1
    if (nfail <= 10) write(logunit,*) trim(message), idx
  end subroutine fail

  !-----------------------------------------------------------------------
  subroutine elmxx_kokkos_seed_topology(elm, logunit)
    !
    ! Push the one piece of topology the Kokkos kernels dereference directly:
    ! for every packed natural patch, the packed index of the column it sits
    ! on. Every kernel that touches both levels does `const int c =
    ! patch.patch_column(p)` and indexes a column view with it, so this map
    ! being wrong is silent corruption rather than a crash.
    !
    implicit none
    type(ELMxxType), intent(in) :: elm
    integer, intent(in) :: logunit
    integer :: kp, p, ierr
    integer, allocatable :: buf(:)
    character(len=*), parameter :: subname = '(elmxx_kokkos_seed_topology) '

    call require_built(subname)

    allocate(buf(n_kokkos_patch))
    do kp = 1, n_kokkos_patch
       p = patch_of_kpatch(kp)
       buf(kp) = kcol_of_col(patch_column(p))
       if (buf(kp) < 0) then
          call shr_sys_abort(subname//'ERROR: natural patch on an unmapped column')
       end if
    end do

    call ELMxxSetPatchColumn(elm, buf, n_kokkos_patch, ierr)
    call check(ierr, subname, 'PatchColumn')
    deallocate(buf)

    write(logunit,*) subname,'rank ',iam,' seeded patch->column for ', &
                     n_kokkos_patch,' packed patches'
    call shr_sys_flush(logunit)

  end subroutine elmxx_kokkos_seed_topology

  !-----------------------------------------------------------------------
  subroutine elmxx_kokkos_seed_state(elm, logunit)
    !
    ! Seed once. The plan's "at init and at restart only" crossing: everything
    ! here is state that does not change per timestep while the kernels are
    ! off, so it is pushed once rather than every step.
    !
    ! WHAT IS SEEDED, AND FROM WHERE:
    !   Elai, Esai, Htop  Stage 2's interpolated satellite phenology. Real
    !                     data, exact against ELM's restart already.
    !   Snl, H2osno,      ELM's cold start (ColumnDataType InitCold):
    !   SnowDepth,        no snow, t_soisno = 274 K for non-lake columns and
    !   FracSno, TGrnd    t_grnd = t_soisno(snl+1), so 274 K here. Lake's
    !                     277 K does not arise: the packed natural-column view
    !                     carries no lake.
    !   TVeg              283 K (VegetationDataType InitCold). The 297.56 /
    !                     289.46 branches are use_vancouver / use_mexicocity,
    !                     both off.
    !
    ! DELIBERATELY NOT SEEDED: the soil hydraulic properties (Watsat, Watfc,
    ! Sucsat, Bsw). ELMxx has setters for them, but they are not surfdata
    ! fields -- ELM derives them from sand, clay and organic matter through the
    ! pedotransfer functions in iniTimeConst. Porting that is real physics
    ! work, it belongs with the kernels that read those arrays, and guessing it
    ! here would put plausible-looking wrong numbers underneath Stage 4. Stage
    ! 2 already stores the raw sand/clay/organic per column, so the inputs are
    ! ready when that port happens.
    !
    implicit none
    type(ELMxxType), intent(in) :: elm
    integer, intent(in) :: logunit
    integer :: kc, kp, p, ierr
    real(r8), allocatable :: rcol(:), rpatch(:)
    integer , allocatable :: icol(:)
    character(len=*), parameter :: subname = '(elmxx_kokkos_seed_state) '
    real(r8), parameter :: t_grnd_cold = 274.0_r8   ! ELM ColumnDataType InitCold
    real(r8), parameter :: t_veg_cold  = 283.0_r8   ! ELM VegetationDataType InitCold

    call require_built(subname)
    if (.not. surface_state_built) then
       call shr_sys_abort(subname//'ERROR: surface state is not ready')
    end if

    allocate(rcol(n_kokkos_col), icol(n_kokkos_col), rpatch(n_kokkos_patch))

    ! ---- natural columns: cold start, no snow ----
    icol = 0
    call ELMxxSetSnl(elm, icol, n_kokkos_col, ierr);       call check(ierr, subname, 'Snl')
    rcol = 0.0_r8
    call ELMxxSetH2osno(elm, rcol, n_kokkos_col, ierr);    call check(ierr, subname, 'H2osno')
    call ELMxxSetSnowDepth(elm, rcol, n_kokkos_col, ierr); call check(ierr, subname, 'SnowDepth')
    call ELMxxSetFracSno(elm, rcol, n_kokkos_col, ierr);   call check(ierr, subname, 'FracSno')
    rcol = t_grnd_cold
    call ELMxxSetTGrnd(elm, rcol, n_kokkos_col, ierr);     call check(ierr, subname, 'TGrnd')

    ! ---- natural patches: cold-start canopy plus real phenology ----
    rpatch = t_veg_cold
    call ELMxxSetTVeg(elm, rpatch, n_kokkos_patch, ierr);  call check(ierr, subname, 'TVeg')

    ! Exposed LAI/SAI, not raw TLAI/TSAI. With no snow the burial term is the
    ! identity, but ELM's SatellitePhenologyMod also zeroes anything below
    ! 0.05 rather than carrying a sliver of canopy, and downstream kernels
    ! branch on exactly that threshold.
    do kp = 1, n_kokkos_patch
       rpatch(kp) = elai_seeded(kp)
    end do
    call ELMxxSetElai(elm, rpatch, n_kokkos_patch, ierr);  call check(ierr, subname, 'Elai')

    do kp = 1, n_kokkos_patch
       rpatch(kp) = esai_seeded(kp)
    end do
    call ELMxxSetEsai(elm, rpatch, n_kokkos_patch, ierr);  call check(ierr, subname, 'Esai')

    do kp = 1, n_kokkos_patch
       rpatch(kp) = patch_height_top(patch_of_kpatch(kp))
    end do
    call ELMxxSetHtop(elm, rpatch, n_kokkos_patch, ierr);  call check(ierr, subname, 'Htop')

    deallocate(rcol, icol, rpatch)

    write(logunit,*) subname,'rank ',iam,' seeded cold-start state and ', &
                     'phenology for ',n_kokkos_col,' columns ',n_kokkos_patch, &
                     ' patches'
    call shr_sys_flush(logunit)

  end subroutine elmxx_kokkos_seed_state

  !-----------------------------------------------------------------------
  subroutine elmxx_kokkos_seed_soil_properties(elm, logunit)
    !
    ! Push the derived soil hydraulic properties, which eight kernels read.
    !
    ! THE 2-D LAYOUT TRAP, RESOLVED. SetView2D wraps the caller's pointer in
    ! the DEVICE VIEW'S OWN LAYOUT and deep-copies -- it does not transpose.
    ! So a Fortran (column, layer) array is correct only when the views are
    ! LayoutLeft, which is ELMxx's default (LAYOUT_RIGHT is an option, and it
    ! is OFF). Under LayoutRight the same call would silently transpose,
    ! putting layer 2's value in column 2 and reporting success -- the same
    ! class of bug as the CDL dimension-order reversal in STATUS G.
    !
    ! Queried at runtime rather than assumed, because it is a build option of
    ! a submodule and nothing else here would notice it changing.
    !
    implicit none
    type(ELMxxType), intent(in) :: elm
    integer, intent(in) :: logunit
    integer :: kc, c, j, kp, ierr, sz(2)
    real(r8), allocatable :: buf(:,:), rcol1(:), rpatch1(:)
    integer , allocatable :: ipatch1(:)
    real(r8), parameter :: smpmin_const  = -1.0e8_r8   ! ELM SoilStateType
    real(r8), parameter :: t_h2osfc_cold = 274.0_r8    ! ELM ColumnDataType InitCold
    character(len=*), parameter :: subname = '(elmxx_kokkos_seed_soil_properties) '

    call require_built(subname)
    if (.not. soil_prop_built) then
       call shr_sys_abort(subname//'ERROR: soil properties are not built')
    end if

    if (ELMxxKokkosIsLayoutRight() /= 0) then
       call shr_sys_abort(subname//'ERROR: ELMxx was built LayoutRight; the '// &
            '(column,layer) buffers here are Fortran-ordered and would be '// &
            'silently transposed. Transpose them here before removing this.')
    end if

    sz(1) = n_kokkos_col
    sz(2) = nlevgrnd
    allocate(buf(n_kokkos_col, nlevgrnd))

    do j = 1, nlevgrnd
       do kc = 1, n_kokkos_col
          buf(kc,j) = sp_watsat(col_of_kcol(kc), j)
       end do
    end do
    call ELMxxSetWatsat(elm, buf, sz, ierr); call check(ierr, subname, 'Watsat')

    do j = 1, nlevgrnd
       do kc = 1, n_kokkos_col
          buf(kc,j) = sp_bsw(col_of_kcol(kc), j)
       end do
    end do
    call ELMxxSetBsw(elm, buf, sz, ierr);    call check(ierr, subname, 'Bsw')

    do j = 1, nlevgrnd
       do kc = 1, n_kokkos_col
          buf(kc,j) = sp_sucsat(col_of_kcol(kc), j)
       end do
    end do
    call ELMxxSetSucsat(elm, buf, sz, ierr); call check(ierr, subname, 'Sucsat')

    do j = 1, nlevgrnd
       do kc = 1, n_kokkos_col
          buf(kc,j) = sp_watfc(col_of_kcol(kc), j)
       end do
    end do
    call ELMxxSetWatfc(elm, buf, sz, ierr);  call check(ierr, subname, 'Watfc')

    deallocate(buf)

    ! ---- cold-start column state, on the snow+ground index space ----
    ! These views are (nc, NLEVTOT = 20), not NLEVGRND: they span ELM's
    ! -nlevsno+1 .. nlevgrnd. Passing a 15-wide buffer would be rejected with
    ! ELMXX_ERR_SIZE_MISMATCH and leave the view at zero, silently.
    sz(2) = nlevtot
    allocate(buf(n_kokkos_col, nlevtot))

    do j = 1, nlevtot
       do kc = 1, n_kokkos_col
          buf(kc,j) = sp_dz(col_of_kcol(kc), j)
       end do
    end do
    call ELMxxSetDz(elm, buf, sz, ierr);        call check(ierr, subname, 'Dz')

    do j = 1, nlevtot
       do kc = 1, n_kokkos_col
          buf(kc,j) = sp_tsoisno(col_of_kcol(kc), j)
       end do
    end do
    call ELMxxSetTSoisno(elm, buf, sz, ierr);   call check(ierr, subname, 'TSoisno')

    do j = 1, nlevtot
       do kc = 1, n_kokkos_col
          buf(kc,j) = sp_liq(col_of_kcol(kc), j)
       end do
    end do
    call ELMxxSetH2osoiLiq(elm, buf, sz, ierr); call check(ierr, subname, 'H2osoiLiq')

    do j = 1, nlevtot
       do kc = 1, n_kokkos_col
          buf(kc,j) = sp_ice(col_of_kcol(kc), j)
       end do
    end do
    call ELMxxSetH2osoiIce(elm, buf, sz, ierr); call check(ierr, subname, 'H2osoiIce')

    deallocate(buf)

    ! ---- the scalars CanopyTemperature reads but does not write ----
    ! Only four: everything else it appears to want -- zii, the three
    ! forc_hgt_*_patch, t_ssbef, t_h2osfc_bef -- it COMPUTES, exactly as ELM
    ! does. Seeding those would have been wasted work at best and a
    ! disagreement with the kernel at worst.
    allocate(rcol1(n_kokkos_col))
    rcol1 = smpmin_const
    call ELMxxSetSmpmin(elm, rcol1, n_kokkos_col, ierr);  call check(ierr, subname, 'Smpmin')
    rcol1 = t_h2osfc_cold
    call ELMxxSetTH2osfc(elm, rcol1, n_kokkos_col, ierr); call check(ierr, subname, 'TH2osfc')
    deallocate(rcol1)

    allocate(ipatch1(n_kokkos_patch), rpatch1(n_kokkos_patch))
    do kp = 1, n_kokkos_patch
       ipatch1(kp) = patch_itype(patch_of_kpatch(kp))
    end do
    call ELMxxSetPatchItype(elm, ipatch1, n_kokkos_patch, ierr)
    call check(ierr, subname, 'PatchItype')

    deallocate(ipatch1, rpatch1)

    ! NOTE: forc_hgt_patch is NOT seeded here. It comes from forc_z, which the
    ! coupler import fills in lnd_run_mct -- after init. Seeding it here
    ! captured the zero the array is allocated with, which made
    ! zldis = z0mg in BareGroundFluxes, so log(zldis/z0mg) = 0 and ustar
    ! divided by zero. It is crossed each step instead, with the rest of the
    ! forcing.

    write(logunit,*) subname,'rank ',iam,' seeded watsat/bsw/sucsat/watfc (', &
                     nlevgrnd,' layers) and dz/t_soisno/h2osoi_liq/h2osoi_ice (', &
                     nlevtot,' layers) for ',n_kokkos_col,' columns'
    call shr_sys_flush(logunit)

  end subroutine elmxx_kokkos_seed_soil_properties

  !-----------------------------------------------------------------------
  subroutine elmxx_kokkos_seed_albedo(elm, logunit)
    !
    ! ELM's cold-start albedos (SurfaceAlbedoType InitCold).
    !
    ! WHY THIS IS ENOUGH TO START, AND WHY IT IS NOT ENOUGH TO FINISH.
    !
    ! SurfaceAlbedo runs at the END of ELM's timestep (elm_driver line ~1500,
    ! gated on doalb), so the albedos SurfaceRadiation reads at step N were
    ! computed at the end of step N-1. The dependency is TEMPORAL, not
    ! within-step -- which is why SurfaceRadiation can run before
    ! SurfaceAlbedo exists, provided something supplies step one's values.
    !
    ! ELM supplies them exactly here: lnd_run_mct wraps elm_drv in
    ! `do while (.not. dosend)` and sets doalb = .false. at nstep == 0, so the
    ! first pass runs on InitCold constants -- 0.2 for ground and canopy, 0.6
    ! for snow, ftdd and ftii unity (fully transmitting), fabd/fabi and the
    ! SNICAR flx_abs* zero.
    !
    ! THE LIMIT, STATED PLAINLY: frozen at these values, the radiation does
    ! not evolve. No solar-zenith dependence, no snow aging, no wetness
    ! effect. Kernels downstream become exercisable and gradeable against
    ! bounds, but a multi-day run is not physically right until SurfaceAlbedo
    ! is ported. Do not mistake "runs and looks sensible" for "correct".
    !
    implicit none
    type(ELMxxType), intent(in) :: elm
    integer, intent(in) :: logunit
    integer :: ierr, szc(2), szp(2), szs(2)
    real(r8), allocatable :: bufc(:,:), bufp(:,:), bufs(:,:)
    integer, parameter :: numrad = 2               ! VIS, NIR
    integer, parameter :: nsnowlyr = 6             ! NLEVSNO + 1
    real(r8), parameter :: alb_ground = 0.2_r8     ! ELM InitCold
    real(r8), parameter :: alb_snow   = 0.6_r8
    character(len=*), parameter :: subname = '(elmxx_kokkos_seed_albedo) '

    call require_built(subname)

    szc = (/ n_kokkos_col,   numrad /)
    szp = (/ n_kokkos_patch, numrad /)
    szs = (/ n_kokkos_col,   nsnowlyr /)
    allocate(bufc(n_kokkos_col, numrad), bufp(n_kokkos_patch, numrad), &
             bufs(n_kokkos_col, nsnowlyr))

    bufc = alb_ground
    call ELMxxSetAlbgrd(elm, bufc, szc, ierr); call check(ierr, subname, 'Albgrd')
    call ELMxxSetAlbgri(elm, bufc, szc, ierr); call check(ierr, subname, 'Albgri')
    call ELMxxSetAlbsod(elm, bufc, szc, ierr); call check(ierr, subname, 'Albsod')
    call ELMxxSetAlbsoi(elm, bufc, szc, ierr); call check(ierr, subname, 'Albsoi')
    bufc = alb_snow
    call ELMxxSetAlbsndHst(elm, bufc, szc, ierr); call check(ierr, subname, 'AlbsndHst')
    call ELMxxSetAlbsniHst(elm, bufc, szc, ierr); call check(ierr, subname, 'AlbsniHst')

    bufp = alb_ground
    call ELMxxSetAlbd(elm, bufp, szp, ierr); call check(ierr, subname, 'Albd')
    call ELMxxSetAlbi(elm, bufp, szp, ierr); call check(ierr, subname, 'Albi')
    bufp = 0.0_r8
    call ELMxxSetFabd(elm, bufp, szp, ierr); call check(ierr, subname, 'Fabd')
    call ELMxxSetFabi(elm, bufp, szp, ierr); call check(ierr, subname, 'Fabi')
    call ELMxxSetFtid(elm, bufp, szp, ierr); call check(ierr, subname, 'Ftid')
    bufp = 1.0_r8
    call ELMxxSetFtdd(elm, bufp, szp, ierr); call check(ierr, subname, 'Ftdd')
    call ELMxxSetFtii(elm, bufp, szp, ierr); call check(ierr, subname, 'Ftii')

    bufs = 0.0_r8
    call ELMxxSetFlxAbsdv(elm, bufs, szs, ierr); call check(ierr, subname, 'FlxAbsdv')
    call ELMxxSetFlxAbsdn(elm, bufs, szs, ierr); call check(ierr, subname, 'FlxAbsdn')
    call ELMxxSetFlxAbsiv(elm, bufs, szs, ierr); call check(ierr, subname, 'FlxAbsiv')
    call ELMxxSetFlxAbsin(elm, bufs, szs, ierr); call check(ierr, subname, 'FlxAbsin')

    deallocate(bufc, bufp, bufs)

    write(logunit,*) subname,'rank ',iam,' seeded ELM cold-start albedos ', &
                     '(ground/canopy 0.2, snow 0.6, ftdd/ftii 1) for ', &
                     n_kokkos_col,' columns ',n_kokkos_patch,' patches'
    call shr_sys_flush(logunit)

  end subroutine elmxx_kokkos_seed_albedo

  !-----------------------------------------------------------------------
  subroutine elmxx_kokkos_seed_canopy_hydrology(elm, logunit)
    !
    ! The extra inputs CanopyHydrology needs beyond elmxx_kokkos_seed_state.
    ! Split out because they belong to one kernel: keeping them here makes it
    ! obvious what activating that kernel required, and what the next kernel
    ! will have to add.
    !
    ! TWO DERIVED PARAMETERS, from ELM's initVerticalMod:
    !   n_melt      = 200 / max(10, STD_ELEV)   -- snow-cover shape function.
    !                 ELM's istice_mec branch (a flat 10) cannot arise: the
    !                 packed natural-column view carries no glacier.
    !   micro_sigma = (SLOPE + slope0)**(-slopebeta), slopebeta = 3,
    !                 slopemax = 0.4, slope0 = slopemax**(-1/slopebeta).
    !                 Microtopography, in meters.
    ! Both come from surfdata fields ELMxx did not read until now.
    !
    ! THE REST ARE COLD-START ZEROS: no canopy water, no surface water, no
    ! snow. dewmx is ELM's constant 0.1.
    !
    implicit none
    type(ELMxxType), intent(in) :: elm
    integer, intent(in) :: logunit
    integer :: kc, kp, g, ierr
    real(r8), allocatable :: rcol(:), rpatch(:)
    integer , allocatable :: icol(:), ipatch(:)
    character(len=*), parameter :: subname = '(elmxx_kokkos_seed_canopy_hydrology) '
    real(r8), parameter :: dewmx_const = 0.1_r8    ! ELM CanopyHydrologyMod
    real(r8), parameter :: slopebeta   = 3.0_r8    ! ELM initVerticalMod
    real(r8), parameter :: slopemax    = 0.4_r8
    real(r8) :: slope0

    call require_built(subname)
    slope0 = slopemax**(-1.0_r8/slopebeta)

    allocate(rcol(n_kokkos_col), icol(n_kokkos_col), &
             rpatch(n_kokkos_patch), ipatch(n_kokkos_patch))

    ! ---- derived microtopography, per natural column ----
    do kc = 1, n_kokkos_col
       g = cell_of_kcol(kc)
       rcol(kc) = 200.0_r8 / max(10.0_r8, topo_std(g))
    end do
    call ELMxxSetNMelt(elm, rcol, n_kokkos_col, ierr);      call check(ierr, subname, 'NMelt')

    do kc = 1, n_kokkos_col
       g = cell_of_kcol(kc)
       rcol(kc) = (topo_slope(g) + slope0)**(-slopebeta)
    end do
    call ELMxxSetMicroSigma(elm, rcol, n_kokkos_col, ierr); call check(ierr, subname, 'MicroSigma')

    ! ---- cold start: no surface water, no snow beyond what seed_state set ----
    rcol = 0.0_r8
    call ELMxxSetH2osfc(elm, rcol, n_kokkos_col, ierr);     call check(ierr, subname, 'H2osfc')
    call ELMxxSetIntSnow(elm, rcol, n_kokkos_col, ierr);    call check(ierr, subname, 'IntSnow')
    call ELMxxSetFracH2osfc(elm, rcol, n_kokkos_col, ierr); call check(ierr, subname, 'FracH2osfc')
    call ELMxxSetFracSnoEff(elm, rcol, n_kokkos_col, ierr); call check(ierr, subname, 'FracSnoEff')
    icol = 0
    call ELMxxSetDoCapsnow(elm, icol, n_kokkos_col, ierr);  call check(ierr, subname, 'DoCapsnow')

    ! ---- per patch ----
    rpatch = dewmx_const
    call ELMxxSetDewmx(elm, rpatch, n_kokkos_patch, ierr);  call check(ierr, subname, 'Dewmx')
    rpatch = 0.0_r8
    call ELMxxSetH2ocan(elm, rpatch, n_kokkos_patch, ierr); call check(ierr, subname, 'H2ocan')
    call ELMxxSetFwet(elm, rpatch, n_kokkos_patch, ierr);   call check(ierr, subname, 'Fwet')

    ! frac_veg_nosno follows ELM's SatellitePhenologyMod rule: vegetated for
    ! radiation purposes iff exposed LAI + SAI reaches 0.05. elm_driver copies
    ! frac_veg_nosno_alb into frac_veg_nosno when dynamic vegetation is off,
    ! which is the SP-mode case here.
    do kp = 1, n_kokkos_patch
       if (elai_seeded(kp) + esai_seeded(kp) >= 0.05_r8) then
          ipatch(kp) = 1
       else
          ipatch(kp) = 0
       end if
    end do
    call ELMxxSetFracVegNosno(elm, ipatch, n_kokkos_patch, ierr)
    call check(ierr, subname, 'FracVegNosno')
    write(logunit,*) subname,'rank ',iam,' frac_veg_nosno = 1 on ', &
         count(ipatch == 1),' of ',n_kokkos_patch,' packed patches'

    deallocate(rcol, icol, rpatch, ipatch)

    write(logunit,*) subname,'rank ',iam,' seeded CanopyHydrology inputs for ', &
                     n_kokkos_col,' columns ',n_kokkos_patch,' patches'
    call shr_sys_flush(logunit)

  end subroutine elmxx_kokkos_seed_canopy_hydrology

  !-----------------------------------------------------------------------
  ! Exposed LAI/SAI as seeded, matching what elmxx_kokkos_seed_state pushed.
  ! ELM's SatellitePhenology zeroes anything below 0.05 rather than carrying a
  ! sliver of canopy, and frac_veg_nosno is decided on the zeroed values.
  !-----------------------------------------------------------------------
  real(r8) function elai_seeded(kp)
    implicit none
    integer, intent(in) :: kp
    elai_seeded = patch_lai(patch_of_kpatch(kp))
    if (elai_seeded < 0.05_r8) elai_seeded = 0.0_r8
  end function elai_seeded

  real(r8) function esai_seeded(kp)
    implicit none
    integer, intent(in) :: kp
    esai_seeded = patch_sai(patch_of_kpatch(kp))
    if (esai_seeded < 0.05_r8) esai_seeded = 0.0_r8
  end function esai_seeded


  !-----------------------------------------------------------------------
  subroutine elmxx_kokkos_seed_pftpar(elm, logunit)
    !
    ! The two PFT-parameter views ELMxx carries: z0mr and displar.
    !
    ! THESE ARE NOT OPTIONAL TUNING KNOBS. CanopyTemperature forms
    !     z0mv = z0mr(pft) * htop      displa = displar(pft) * htop
    ! and CanopyFluxes then builds its entire resistance chain on them --
    ! ustar, ram1, uaf, rb, and every conductance in the canopy energy
    ! balance. At zero the chain does not merely become inert: ustar is zero,
    ! ram1 is Infinity, rb is NaN, and every canopy flux on every vegetated
    ! patch is NaN.
    !
    ! That is precisely what happened, and it went unseen through several
    ! runs that reported clean ranges, because MINVAL AND MAXVAL SKIP NaN.
    ! Seventeen patches with three NaN among them printed a tidy finite
    ! interval. The kernel reports now carry non-finite counts alongside every
    ! range for exactly this reason.
    !
    ! The views are (NUMPFT+1) = 17 long -- PFT 0 is bare ground -- while the
    ! parameter file carries 25. The extra entries are crops ELMxx does not
    ! index; passing all 25 would be rejected on length.
    !
    implicit none
    type(ELMxxType), intent(in) :: elm
    integer, intent(in) :: logunit
    integer :: ierr, n
    integer, parameter :: numpft_kokkos = 17     ! NUMPFT + 1
    real(r8) :: buf(numpft_kokkos)
    character(len=*), parameter :: subname = '(elmxx_kokkos_seed_pftpar) '

    call require_built(subname)
    if (.not. pftcon_read) then
       call shr_sys_abort(subname//'ERROR: PFT parameters are not read; '// &
            'fparamfile must be set before the canopy kernels can run')
    end if

    n = min(numpft_kokkos, npft_param)
    if (n < numpft_kokkos) then
       call shr_sys_abort(subname//'ERROR: parameter file carries fewer PFTs '// &
            'than ELMxx indexes')
    end if

    buf(1:numpft_kokkos) = z0mr(0:numpft_kokkos-1)
    call ELMxxSetZ0mrPft(elm, buf, numpft_kokkos, ierr)
    call check(ierr, subname, 'Z0mrPft')

    buf(1:numpft_kokkos) = displar(0:numpft_kokkos-1)
    call ELMxxSetDisplarPft(elm, buf, numpft_kokkos, ierr)
    call check(ierr, subname, 'DisplarPft')

    write(logunit,*) subname,'rank ',iam,' seeded ',numpft_kokkos, &
         ' PFT roughness/displacement ratios; z0mr ',minval(z0mr(1:numpft_kokkos-1)), &
         ' .. ',maxval(z0mr(1:numpft_kokkos-1))
    call shr_sys_flush(logunit)

  end subroutine elmxx_kokkos_seed_pftpar

  !-----------------------------------------------------------------------
  subroutine elmxx_kokkos_push_forcing(elm, logunit)
    !
    ! The per-timestep "atmospheric forcing in" crossing -- one of the two the
    ! plan allows. Everything here genuinely changes every coupling interval;
    ! anything that does not belongs in elmxx_kokkos_seed_state.
    !
    ! Forcing arrives per gridcell, so it is broadcast down to the columns and
    ! patches through the packed maps rather than by loop position. Column c
    ! belongs to gridcell lun_gridcell(col_landunit(c)); that is the only
    ! defensible route from one level to the other.
    !
    ! ALL FORCING IS NOW ACROSS. The two that lagged are both here:
    ! ForcRhoCol, derived from vapor pressure as ELM's import derives it, and
    ! ForcSolad / ForcSolai, 2-D over the radiation bands. The latter waited on
    ! the LayoutLeft question, which elmxx_kokkos_state_init settles at startup
    ! by querying ELMxxKokkosIsLayoutRight and aborting on a LayoutRight build.
    !
    implicit none
    type(ELMxxType), intent(in) :: elm
    integer, intent(in) :: logunit
    integer :: kc, kp, g, ierr
    integer :: szp(2)
    real(r8) :: vp
    real(r8), allocatable :: rcol(:), rpatch(:)
    real(r8), allocatable :: rsol(:,:)
    logical, save :: reported = .false.
    integer, parameter :: numrad = 2   ! 1 = visible, 2 = near-IR
    ! ELM's rair is elm_varcon's, which is SHR_CONST_RDAIR. Taken from the
    ! same shared constant rather than transcribed: it is derived
    ! (RGAS/MWDAIR), not a literal, so a hand-copied value would differ in the
    ! last digits and quietly bound how close a future ELM comparison can get.
    real(r8), parameter :: rair = SHR_CONST_RDAIR
    character(len=*), parameter :: subname = '(elmxx_kokkos_push_forcing) '

    call require_built(subname)
    allocate(rcol(n_kokkos_col), rpatch(n_kokkos_patch), &
             rsol(n_kokkos_patch, numrad))

    ! ---- column-level ----
    do kc = 1, n_kokkos_col
       rcol(kc) = forc_tbot(cell_of_kcol(kc))
    end do
    call ELMxxSetForcTCol(elm, rcol, n_kokkos_col, ierr);     call check(ierr, subname, 'ForcTCol')

    do kc = 1, n_kokkos_col
       rcol(kc) = forc_pbot(cell_of_kcol(kc))
    end do
    call ELMxxSetForcPbotCol(elm, rcol, n_kokkos_col, ierr);  call check(ierr, subname, 'ForcPbotCol')

    do kc = 1, n_kokkos_col
       rcol(kc) = forc_shum(cell_of_kcol(kc))
    end do
    call ELMxxSetForcQCol(elm, rcol, n_kokkos_col, ierr);     call check(ierr, subname, 'ForcQCol')

    do kc = 1, n_kokkos_col
       rcol(kc) = forc_lwrad(cell_of_kcol(kc))
    end do
    call ELMxxSetForcLwradCol(elm, rcol, n_kokkos_col, ierr); call check(ierr, subname, 'ForcLwradCol')

    do kc = 1, n_kokkos_col
       rcol(kc) = forc_u(cell_of_kcol(kc))
    end do
    call ELMxxSetForcUCol(elm, rcol, n_kokkos_col, ierr);     call check(ierr, subname, 'ForcUCol')

    do kc = 1, n_kokkos_col
       rcol(kc) = forc_v(cell_of_kcol(kc))
    end do
    call ELMxxSetForcVCol(elm, rcol, n_kokkos_col, ierr);     call check(ierr, subname, 'ForcVCol')

    do kc = 1, n_kokkos_col
       rcol(kc) = forc_ptem(cell_of_kcol(kc))
    end do
    call ELMxxSetForcThCol(elm, rcol, n_kokkos_col, ierr);    call check(ierr, subname, 'ForcThCol')

    ! Air density. ELM derives this in its coupler import rather than
    ! receiving it, so ELMxx has to derive it too, from the same three fields
    ! and the same two-step form -- vapor pressure first, then density with
    ! the 0.378 moist-air correction. Writing it as one expression would
    ! change the rounding and make a future comparison against ELM harder to
    ! read than it needs to be.
    do kc = 1, n_kokkos_col
       g  = cell_of_kcol(kc)
       vp = forc_shum(g) * forc_pbot(g) / (0.622_r8 + 0.378_r8 * forc_shum(g))
       rcol(kc) = (forc_pbot(g) - 0.378_r8 * vp) / (rair * forc_tbot(g))
    end do
    call ELMxxSetForcRhoCol(elm, rcol, n_kokkos_col, ierr); call check(ierr, subname, 'ForcRhoCol')

    ! ---- patch-level ----
    do kp = 1, n_kokkos_patch
       rpatch(kp) = forc_tbot(cell_of_kpatch(kp))
    end do
    call ELMxxSetForcT(elm, rpatch, n_kokkos_patch, ierr);    call check(ierr, subname, 'ForcT')

    ! Convective and large-scale are separate on the coupler side and summed
    ! here, which is what ELM's import does.
    do kp = 1, n_kokkos_patch
       g = cell_of_kpatch(kp)
       rpatch(kp) = forc_rainc(g) + forc_rainl(g)
    end do
    call ELMxxSetForcRain(elm, rpatch, n_kokkos_patch, ierr); call check(ierr, subname, 'ForcRain')

    do kp = 1, n_kokkos_patch
       g = cell_of_kpatch(kp)
       rpatch(kp) = forc_snowc(g) + forc_snowl(g)
    end do
    call ELMxxSetForcSnow(elm, rpatch, n_kokkos_patch, ierr); call check(ierr, subname, 'ForcSnow')

    ! The atmospheric reference height. This is forcing, not init state:
    ! CanopyTemperature adds roughness and displacement to it to get the
    ! per-patch heights, and BareGroundFluxes then takes log(zldis/z0mg), so a
    ! zero here is not a small error -- it is a division by zero.
    do kp = 1, n_kokkos_patch
       rpatch(kp) = forc_z(cell_of_kpatch(kp))
    end do
    call ELMxxSetForcHgtPatch(elm, rpatch, n_kokkos_patch, ierr)
    call check(ierr, subname, 'ForcHgtPatch')

    ! Incident shortwave, the only 2-D forcing. Direct and diffuse are separate
    ! views, and the band index is ELM's: 1 is visible, 2 is near-IR. That is
    ! lnd_import_export's assignment, not a convention chosen here --
    !     forc_solad(g,1) = swvdr    forc_solad(g,2) = swndr
    !     forc_solai(g,1) = swvdf    forc_solai(g,2) = swndf
    ! -- and getting it backwards would swap the bands against albedos that
    ! genuinely differ between them (0.2 is the cold-start constant in both,
    ! so this would hide today and surface the moment SurfaceAlbedo lands).
    !
    ! LAYOUT: rsol is (patch, band) and crosses under LayoutLeft, which is the
    ! Fortran-native order SetView2D wraps without transposing. A LayoutRight
    ! build would transpose here AND REPORT SUCCESS; that is why startup
    ! aborts on one rather than trusting this comment.
    do kp = 1, n_kokkos_patch
       g = cell_of_kpatch(kp)
       rsol(kp,1) = forc_swvdr(g)
       rsol(kp,2) = forc_swndr(g)
    end do
    szp = (/ n_kokkos_patch, numrad /)
    call ELMxxSetForcSolad(elm, rsol, szp, ierr); call check(ierr, subname, 'ForcSolad')

    do kp = 1, n_kokkos_patch
       g = cell_of_kpatch(kp)
       rsol(kp,1) = forc_swvdf(g)
       rsol(kp,2) = forc_swndf(g)
    end do
    call ELMxxSetForcSolai(elm, rsol, szp, ierr); call check(ierr, subname, 'ForcSolai')

    if (.not. reported) then
       write(logunit,*) subname,'rank ',iam,' forc_z range ', &
            minval(forc_z(1:size(forc_z))),' .. ',maxval(forc_z(1:size(forc_z))),' m'
    end if

    deallocate(rcol, rpatch, rsol)

    if (.not. reported) then
       write(logunit,*) subname,'rank ',iam,' pushing forcing each step to ', &
                        n_kokkos_col,' columns ',n_kokkos_patch,' patches'
       call shr_sys_flush(logunit)
       reported = .true.
    end if

  end subroutine elmxx_kokkos_push_forcing

  !-----------------------------------------------------------------------
  subroutine elmxx_kokkos_push_btran(elm, logunit)
    !
    ! Cross btran, which CanopyFluxes reads and does not compute.
    !
    ! Per-step rather than seeded: it is a function of soil moisture and
    ! temperature, so it goes stale as soon as hydrology runs. Cheap enough
    ! that recomputing beats reasoning about when it last changed.
    !
    implicit none
    type(ELMxxType), intent(in) :: elm
    integer, intent(in) :: logunit
    integer :: kp, ierr
    real(r8), allocatable :: buf(:)
    character(len=*), parameter :: subname = '(elmxx_kokkos_push_btran) '

    call require_built(subname)
    allocate(buf(n_kokkos_patch))
    do kp = 1, n_kokkos_patch
       buf(kp) = rt_btran(patch_of_kpatch(kp))
    end do
    call ELMxxSetBtran(elm, buf, n_kokkos_patch, ierr)
    call check(ierr, subname, 'Btran')
    deallocate(buf)

  end subroutine elmxx_kokkos_push_btran

  !-----------------------------------------------------------------------
  integer function cell_of_kcol(kc)
    implicit none
    integer, intent(in) :: kc
    cell_of_kcol = lun_gridcell(col_landunit(col_of_kcol(kc)))
  end function cell_of_kcol

  !-----------------------------------------------------------------------
  integer function cell_of_kpatch(kp)
    implicit none
    integer, intent(in) :: kp
    cell_of_kpatch = lun_gridcell(col_landunit(patch_column(patch_of_kpatch(kp))))
  end function cell_of_kpatch

  !-----------------------------------------------------------------------
  subroutine elmxx_kokkos_verify_maps(elm, logunit, nfail)
    !
    ! Round-trip a per-entity fingerprint through the Kokkos Views.
    !
    ! WHAT THIS PROVES, EXACTLY. Both the write and the read go through the
    ! SAME map, so this cannot prove the map is semantically right -- a
    ! consistently wrong map round-trips perfectly. What it does prove is the
    ! plumbing: that every packed slot is addressable, that the extents the
    ! setters were handed match the views ELMxxCreate allocated, that no setter
    ! was silently rejected and left its view at zero (STATUS.md E.1), and that
    ! the int and double paths both work.
    !
    ! The semantics are checked separately and independently, without Kokkos in
    ! the picture at all, by elmxx_kokkos_check_map_invariants below. Neither
    ! check subsumes the other and Stage 3 needs both.
    !
    ! A fingerprint rather than physical state because a cold start gives every
    ! column the same t_grnd, so a zeroed or aliased view would still compare
    ! equal. Encoding identity also makes a failure message name the entity.
    !
    ! Runs with every kernel off, so nothing between the set and the get can
    ! legitimately change a value: any difference is the boundary's fault.
    !
    implicit none
    type(ELMxxType), intent(in) :: elm
    integer, intent(in) :: logunit
    integer, intent(out) :: nfail
    character(len=*), parameter :: subname = '(elmxx_kokkos_verify_maps) '

    call require_built(subname)
    nfail = 0

    call probe_col_real(elm, logunit, nfail, 'SnowDepth', 1)
    call probe_col_real(elm, logunit, nfail, 'FracSno'  , 2)
    call probe_col_real(elm, logunit, nfail, 'H2osno'   , 3)
    call probe_col_real(elm, logunit, nfail, 'TGrnd'    , 4)
    call probe_col_int (elm, logunit, nfail)
    call probe_patch_real(elm, logunit, nfail, 'TVeg', 1)
    call probe_patch_real(elm, logunit, nfail, 'Fsun', 2)
    if (n_kokkos_urb > 0) then
       call probe_urban_real(elm, logunit, nfail, 'UrbanTaf', 1)
       call probe_urban_real(elm, logunit, nfail, 'UrbanQaf', 2)
    end if

    if (nfail == 0) then
       write(logunit,*) subname,'rank ',iam,' PASSED: packed maps round-trip ', &
                        'exactly for ',n_kokkos_col,' columns ',n_kokkos_patch, &
                        ' patches ',n_kokkos_urb,' urban landunits'
    else
       write(logunit,*) subname,'rank ',iam,' FAILED: ',nfail,' mismatches'
    end if
    call shr_sys_flush(logunit)

  end subroutine elmxx_kokkos_verify_maps

  !-----------------------------------------------------------------------
  ! Fingerprints. Distinct per entity and decodable by eye in a log message.
  !-----------------------------------------------------------------------
  real(r8) function fp_col(c, salt)
    implicit none
    integer, intent(in) :: c, salt
    fp_col = real(1000000*salt + 1000*lun_itype(col_landunit(c)) + c, r8) &
             + real(col_itype(c), r8) / 1000.0_r8
  end function fp_col

  real(r8) function fp_patch(p, salt)
    implicit none
    integer, intent(in) :: p, salt
    fp_patch = real(1000000*salt + 1000*patch_itype(p) + p, r8) &
               + real(patch_column(p), r8) / 1000.0_r8
  end function fp_patch

  real(r8) function fp_lun(l, salt)
    implicit none
    integer, intent(in) :: l, salt
    fp_lun = real(1000000*salt + 1000*lun_itype(l) + l, r8) &
             + real(lun_gridcell(l), r8) / 1000.0_r8
  end function fp_lun

  !-----------------------------------------------------------------------
  subroutine probe_col_real(elm, logunit, nfail, field, salt)
    implicit none
    type(ELMxxType), intent(in) :: elm
    integer, intent(in) :: logunit, salt
    integer, intent(inout) :: nfail
    character(len=*), intent(in) :: field
    integer :: kc, c, ierr, bad
    real(r8), allocatable :: put(:), got(:)
    character(len=*), parameter :: subname = '(elmxx_kokkos_verify_maps) '

    allocate(put(n_kokkos_col), got(n_kokkos_col))
    do kc = 1, n_kokkos_col
       put(kc) = fp_col(col_of_kcol(kc), salt)
    end do
    got = -huge(1.0_r8)

    select case (field)
    case ('SnowDepth')
       call ELMxxSetSnowDepth(elm, put, n_kokkos_col, ierr); call check(ierr, subname, field)
       call ELMxxGetSnowDepth(elm, got, n_kokkos_col, ierr); call check(ierr, subname, field)
    case ('FracSno')
       call ELMxxSetFracSno(elm, put, n_kokkos_col, ierr);   call check(ierr, subname, field)
       call ELMxxGetFracSno(elm, got, n_kokkos_col, ierr);   call check(ierr, subname, field)
    case ('H2osno')
       call ELMxxSetH2osno(elm, put, n_kokkos_col, ierr);    call check(ierr, subname, field)
       call ELMxxGetH2osno(elm, got, n_kokkos_col, ierr);    call check(ierr, subname, field)
    case ('TGrnd')
       call ELMxxSetTGrnd(elm, put, n_kokkos_col, ierr);     call check(ierr, subname, field)
       call ELMxxGetTGrnd(elm, got, n_kokkos_col, ierr);     call check(ierr, subname, field)
    case default
       call shr_sys_abort(subname//'ERROR: unknown column probe '//trim(field))
    end select

    bad = 0
    do kc = 1, n_kokkos_col
       if (got(kc) /= put(kc)) then
          bad = bad + 1
          c = col_of_kcol(kc)
          if (bad <= 5) then
             write(logunit,*) subname,'MISMATCH ',trim(field),' packed ',kc-1, &
                  ' (column ',c,' landunit type ',lun_itype(col_landunit(c)), &
                  ') expected ',put(kc),' got ',got(kc)
          end if
       end if
    end do
    nfail = nfail + bad
    deallocate(put, got)

  end subroutine probe_col_real

  !-----------------------------------------------------------------------
  subroutine probe_col_int(elm, logunit, nfail)
    !
    ! snl is the integer column field; it exercises the int path, which has its
    ! own SetView1D overload.
    !
    implicit none
    type(ELMxxType), intent(in) :: elm
    integer, intent(in) :: logunit
    integer, intent(inout) :: nfail
    integer :: kc, c, ierr, bad
    integer, allocatable :: put(:), got(:)
    character(len=*), parameter :: subname = '(elmxx_kokkos_verify_maps) '

    allocate(put(n_kokkos_col), got(n_kokkos_col))
    do kc = 1, n_kokkos_col
       c = col_of_kcol(kc)
       put(kc) = 1000*lun_itype(col_landunit(c)) + c
    end do
    got = -huge(1)

    call ELMxxSetSnl(elm, put, n_kokkos_col, ierr); call check(ierr, subname, 'Snl')
    call ELMxxGetSnl(elm, got, n_kokkos_col, ierr); call check(ierr, subname, 'Snl')

    bad = 0
    do kc = 1, n_kokkos_col
       if (got(kc) /= put(kc)) then
          bad = bad + 1
          if (bad <= 5) then
             write(logunit,*) subname,'MISMATCH Snl packed ',kc-1,' (column ', &
                  col_of_kcol(kc),') expected ',put(kc),' got ',got(kc)
          end if
       end if
    end do
    nfail = nfail + bad
    deallocate(put, got)

  end subroutine probe_col_int

  !-----------------------------------------------------------------------
  subroutine probe_patch_real(elm, logunit, nfail, field, salt)
    implicit none
    type(ELMxxType), intent(in) :: elm
    integer, intent(in) :: logunit, salt
    integer, intent(inout) :: nfail
    character(len=*), intent(in) :: field
    integer :: kp, p, ierr, bad
    real(r8), allocatable :: put(:), got(:)
    character(len=*), parameter :: subname = '(elmxx_kokkos_verify_maps) '

    allocate(put(n_kokkos_patch), got(n_kokkos_patch))
    do kp = 1, n_kokkos_patch
       put(kp) = fp_patch(patch_of_kpatch(kp), salt)
    end do
    got = -huge(1.0_r8)

    select case (field)
    case ('TVeg')
       call ELMxxSetTVeg(elm, put, n_kokkos_patch, ierr); call check(ierr, subname, field)
       call ELMxxGetTVeg(elm, got, n_kokkos_patch, ierr); call check(ierr, subname, field)
    case ('Fsun')
       call ELMxxSetFsun(elm, put, n_kokkos_patch, ierr); call check(ierr, subname, field)
       call ELMxxGetFsun(elm, got, n_kokkos_patch, ierr); call check(ierr, subname, field)
    case default
       call shr_sys_abort(subname//'ERROR: unknown patch probe '//trim(field))
    end select

    bad = 0
    do kp = 1, n_kokkos_patch
       if (got(kp) /= put(kp)) then
          bad = bad + 1
          p = patch_of_kpatch(kp)
          if (bad <= 5) then
             write(logunit,*) subname,'MISMATCH ',trim(field),' packed ',kp-1, &
                  ' (patch ',p,' pft ',patch_itype(p),') expected ',put(kp), &
                  ' got ',got(kp)
          end if
       end if
    end do
    nfail = nfail + bad
    deallocate(put, got)

  end subroutine probe_patch_real

  !-----------------------------------------------------------------------
  subroutine probe_urban_real(elm, logunit, nfail, field, salt)
    implicit none
    type(ELMxxType), intent(in) :: elm
    integer, intent(in) :: logunit, salt
    integer, intent(inout) :: nfail
    character(len=*), intent(in) :: field
    integer :: ku, l, ierr, bad
    real(r8), allocatable :: put(:), got(:)
    character(len=*), parameter :: subname = '(elmxx_kokkos_verify_maps) '

    allocate(put(n_kokkos_urb), got(n_kokkos_urb))
    do ku = 1, n_kokkos_urb
       put(ku) = fp_lun(lun_of_kurb(ku), salt)
    end do
    got = -huge(1.0_r8)

    select case (field)
    case ('UrbanTaf')
       call ELMxxSetUrbanTaf(elm, put, n_kokkos_urb, ierr); call check(ierr, subname, field)
       call ELMxxGetUrbanTaf(elm, got, n_kokkos_urb, ierr); call check(ierr, subname, field)
    case ('UrbanQaf')
       call ELMxxSetUrbanQaf(elm, put, n_kokkos_urb, ierr); call check(ierr, subname, field)
       call ELMxxGetUrbanQaf(elm, got, n_kokkos_urb, ierr); call check(ierr, subname, field)
    case default
       call shr_sys_abort(subname//'ERROR: unknown urban probe '//trim(field))
    end select

    bad = 0
    do ku = 1, n_kokkos_urb
       if (got(ku) /= put(ku)) then
          bad = bad + 1
          l = lun_of_kurb(ku)
          if (bad <= 5) then
             write(logunit,*) subname,'MISMATCH ',trim(field),' packed ',ku-1, &
                  ' (landunit ',l,' type ',lun_itype(l),') expected ',put(ku), &
                  ' got ',got(ku)
          end if
       end if
    end do
    nfail = nfail + bad
    deallocate(put, got)

  end subroutine probe_urban_real

  !-----------------------------------------------------------------------
  subroutine check(ierr, subname, field)
    !
    ! A rejected setter leaves the view at zero rather than failing loudly, so
    ! every status is checked at the call site. STATUS.md E.1 is what happens
    ! when one is not.
    !
    implicit none
    integer, intent(in) :: ierr
    character(len=*), intent(in) :: subname, field
    if (ierr /= ELMXX_SUCCESS) then
       call shr_sys_abort(subname//'ERROR: '//trim(field)//' returned a non-success status')
    end if
  end subroutine check

  !-----------------------------------------------------------------------
  subroutine require_built(subname)
    implicit none
    character(len=*), intent(in) :: subname
    if (.not. kokkos_state_built) then
       call shr_sys_abort(subname//'ERROR: packed maps are not built')
    end if
  end subroutine require_built

  !-----------------------------------------------------------------------
  subroutine elmxx_kokkos_state_clean()
    implicit none
    if (associated(kcol_of_col))     deallocate(kcol_of_col)
    if (associated(col_of_kcol))     deallocate(col_of_kcol)
    if (associated(kpatch_of_patch)) deallocate(kpatch_of_patch)
    if (associated(patch_of_kpatch)) deallocate(patch_of_kpatch)
    if (associated(kurb_of_lun))     deallocate(kurb_of_lun)
    if (associated(lun_of_kurb))     deallocate(lun_of_kurb)
    kcol_of_col     => null(); col_of_kcol     => null()
    kpatch_of_patch => null(); patch_of_kpatch => null()
    kurb_of_lun     => null(); lun_of_kurb     => null()
    n_kokkos_col = 0; n_kokkos_patch = 0; n_kokkos_urb = 0
    kokkos_state_built = .false.
  end subroutine elmxx_kokkos_state_clean

end module elmxxKokkosStateMod
