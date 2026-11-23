module kelmGlacierMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! The glacier landunit on the coupled path (glacier plan, G4): the Fortran
  ! side of kELM's glacier surface (`glacier`, GlacierColumnData, in the C++).
  ! ELM's plain istice landunit only: no glacier MEC, no firn model, no CISM.
  !
  ! WHAT THIS OWNS
  !   maps        packed glacier column/patch <-> subgrid column/patch, and the
  !               cell of each (the glacier reads its forcing from its cell,
  !               on the device: kELMGlacierGatherForcing).
  !   init        the ice-layer geometry (the soil grid, initVerticalMod) and
  !               ELM's glacier cold start (ColumnDataType InitCold, elm_instMod,
  !               SnowHydrologyMod InitSnowLayers) -- computed here, seeded once
  !               by NAME through the restart registry ("glccol:<field>").
  !   albedo      the bare land-ice albedo, albice (kELMSetAlbice).
  !   history     each cell's glacier column and glacier landunit weight, so
  !               h0 is ELM's landunit-weighted cell mean.
  ! ELM is the reference for behaviour only; nothing is copied.
  !
  ! ONE PATCH PER GLACIER COLUMN, of weight 1 (bare ice, pft 0), so ELM's
  ! glacier patch fields are glacier fields, indexed by packed column k.
  !-----------------------------------------------------------------------

  use shr_kind_mod    , only : r8 => shr_kind_r8
  use shr_sys_mod     , only : shr_sys_abort, shr_sys_flush
  use kelmSpmdMod    , only : masterproc, iam
  use kelmSubgridMod , only : num_landunits, num_columns, num_patches, &
                               lun_gridcell, lun_itype, lun_wtgcell, &
                               col_landunit, patch_column, istsoil, istice, istdlak
  use kelmSoilPropMod, only : nlevgrnd, nlevsno, zsoi, dzsoi, zisoi
  use kelmKokkosStateMod, only : col_of_kcol, n_kokkos_col
  use kelmDiagnosticsMod, only : kelm_diag_enabled, kelm_diag_1d, kelm_diag_2d
  use kelm_mod       , only : kELMType, KELM_SUCCESS, &
                               kELMAllocateGlacierSurface, kELMRestartFieldFind, &
                               kELMRestartFieldSet, kELMRestartFieldGet, kELMRestartFieldInfo, &
                               kELMHistorySetLandunitWeights, &
                               kELMHistorySetGlacierWeights, kELMSetAlbice

  implicit none
  save
  private

  integer, public :: n_glc = 0                 ! packed glacier columns = patches
  logical, public :: glacier_built = .false.
  integer, public, pointer :: col_of_kglc(:)   => null()  ! (n_glc) -> subgrid column
  integer, public, pointer :: patch_of_kglc(:) => null()  ! (n_glc) -> subgrid patch
  integer, public, pointer :: cell_of_kglc(:)  => null()  ! (n_glc) -> 1-based local cell

  public :: kelm_glacier_init
  public :: kelm_glacier_clean
  public :: kelm_glacier_diag_step
  public :: kelm_glacier_apply_finidat

  ! ELM constants the cold start uses (elm_varcon, SnowHydrologyMod).
  real(r8), parameter :: denice      = 917._r8      ! ice density [kg/m3]
  real(r8), parameter :: h2osno_max  = 1000._r8     ! ELM: elm_varcon h2osno_max [mm]
  real(r8), parameter :: bdsno       = 250._r8      ! ELM: elm_varcon bdsno, fresh-snow density [kg/m3]
  real(r8), parameter :: snw_rds_min = 54.526_r8    ! ELM: elm_varcon snw_rds_min [um]
  real(r8), parameter :: zlnd        = 0.01_r8      ! ELM: elm_varcon zlnd [m]
  real(r8), parameter :: t_cold      = 250._r8      ! ELM InitCold: glacier ice and snow [K]
  ! ELM SnowHydrologyMod: thickest top layer with no layer beneath (dzmax_l)
  ! and with layers beneath (dzmax_u), and the thinnest top layer (dzmin).
  real(r8), parameter :: dzmin1      = 0.010_r8
  real(r8), parameter :: dzmax_l1    = 0.03_r8
  real(r8), parameter :: dzmax_u(5)  = (/ 0.02_r8, 0.05_r8, 0.11_r8, 0.23_r8, 0.47_r8 /)

  interface glc_set
     module procedure glc_set_1d, glc_set_2d, glc_set_1i
  end interface glc_set

contains

  !-----------------------------------------------------------------------
  subroutine kelm_glacier_init(elm, logunit, lake_built, urban_built, albice)
    !
    ! Build the glacier maps, allocate the glacier surface and seed it, and
    ! give history the glacier weights. A no-op on a domain without glacier.
    ! After kelm_lake_init and kelm_urban_init: when neither pushed the
    ! natural landunit weights (no lake, no urban), they are pushed here.
    !
    implicit none
    type(kELMType), intent(in) :: elm
    integer , intent(in) :: logunit
    logical , intent(in) :: lake_built, urban_built
    real(r8), intent(in) :: albice(2)
    integer :: ierr
    character(len=*), parameter :: subname = '(kelm_glacier_init) '

    call kelm_glacier_clean()
    call build_maps(logunit)
    if (n_glc == 0) return

    call kELMAllocateGlacierSurface(elm, n_glc, ierr)
    if (ierr /= KELM_SUCCESS) call shr_sys_abort(subname//'ERROR: kELMAllocateGlacierSurface failed')
    call kELMSetAlbice(elm, albice(1), albice(2), ierr)
    if (ierr /= KELM_SUCCESS) call shr_sys_abort(subname//'ERROR: kELMSetAlbice failed')

    call seed_topology(elm)
    call seed_geometry(elm)
    call seed_cold_start(elm, logunit)
    call push_history_weights(elm, logunit, lake_built .or. urban_built)

    glacier_built = .true.
    if (masterproc) then
       write(logunit,*) subname,'glacier surface built: ',n_glc,' columns on rank ',iam, &
            ', albice ',albice
       call shr_sys_flush(logunit)
    end if

  end subroutine kelm_glacier_init

  !-----------------------------------------------------------------------
  subroutine build_maps(logunit)
    implicit none
    integer, intent(in) :: logunit
    integer :: c, p, k
    integer, allocatable :: np_of_col(:), patch_of_col(:)
    character(len=*), parameter :: subname = '(kelm_glacier_init) '

    n_glc = 0
    do c = 1, num_columns
       if (lun_itype(col_landunit(c)) == istice) n_glc = n_glc + 1
    end do
    if (masterproc) write(logunit,*) subname,'rank ',iam,' glacier columns ',n_glc
    if (n_glc == 0) return

    allocate(col_of_kglc(n_glc), patch_of_kglc(n_glc), cell_of_kglc(n_glc))
    allocate(np_of_col(num_columns), patch_of_col(num_columns))
    np_of_col = 0
    do p = 1, num_patches
       np_of_col(patch_column(p)) = np_of_col(patch_column(p)) + 1
       patch_of_col(patch_column(p)) = p
    end do
    k = 0
    do c = 1, num_columns
       if (lun_itype(col_landunit(c)) /= istice) cycle
       k = k + 1
       col_of_kglc(k)  = c
       cell_of_kglc(k) = lun_gridcell(col_landunit(c))
       if (np_of_col(c) /= 1) &
            call shr_sys_abort(subname//'ERROR: a glacier column must carry exactly one patch')
       patch_of_kglc(k) = patch_of_col(c)
    end do
    deallocate(np_of_col, patch_of_col)

  end subroutine build_maps

  !-----------------------------------------------------------------------
  subroutine seed_topology(elm)
    implicit none
    type(kELMType), intent(in) :: elm
    integer :: k
    integer, allocatable :: ibuf(:)

    allocate(ibuf(n_glc))
    do k = 1, n_glc
       ibuf(k) = cell_of_kglc(k) - 1
    end do
    call glc_set(elm, 'glccol:col_gridcell', ibuf)
    call glc_set(elm, 'glccol:col_topounit', ibuf)   ! one topounit per cell
    deallocate(ibuf)

  end subroutine seed_topology

  !-----------------------------------------------------------------------
  subroutine seed_geometry(elm)
    !
    ! The ice layers are the soil grid (ELM initVerticalMod gives a glacier
    ! column the same nlevgrnd layers): the split layer views and the static
    ! ice part of the SoilTemperature layout, which puts ice layer j at slot
    ! nlevsno+1+j (the slot before it is the standing-water node).
    !
    implicit none
    type(kELMType), intent(in) :: elm
    integer :: j, nlevtot_p1
    real(r8), allocatable :: v(:,:)

    allocate(v(n_glc, nlevgrnd))
    do j = 1, nlevgrnd; v(:, j) = dzsoi(j); end do
    call glc_set(elm, 'glccol:thermal.dz_soi', v)
    do j = 1, nlevgrnd; v(:, j) = zsoi(j); end do
    call glc_set(elm, 'glccol:thermal.zc_soi', v)
    deallocate(v)
    allocate(v(n_glc, nlevgrnd + 1))
    do j = 0, nlevgrnd; v(:, j + 1) = zisoi(j); end do
    call glc_set(elm, 'glccol:thermal.zi_soi', v)
    deallocate(v)

    nlevtot_p1 = nlevsno + 1 + nlevgrnd
    allocate(v(n_glc, nlevtot_p1))
    v = 0._r8
    do j = 1, nlevgrnd; v(:, nlevsno + 1 + j) = dzsoi(j); end do
    call glc_set(elm, 'glccol:dz_p1', v)
    v = 0._r8
    do j = 1, nlevgrnd; v(:, nlevsno + 1 + j) = zsoi(j); end do
    call glc_set(elm, 'glccol:z_p1', v)
    deallocate(v)
    allocate(v(n_glc, nlevtot_p1 + 1))
    v = 0._r8
    do j = 0, nlevgrnd; v(:, nlevsno + 1 + j + 1) = zisoi(j); end do
    call glc_set(elm, 'glccol:zi_p1', v)
    deallocate(v)

  end subroutine seed_geometry

  !-----------------------------------------------------------------------
  subroutine seed_cold_start(elm, logunit)
    !
    ! ELM's glacier cold start, use_extrasnowlayers off:
    !   ice      t_soisno = 250 K on every layer; h2osoi_vol = 1, so
    !            h2osoi_ice = dz*denice and h2osoi_liq = 0 (ColumnDataType)
    !   snow     h2osno = h2osno_max, snow_depth = h2osno/bdsno (elm_instMod);
    !            the layers from snow_depth (SnowHydrologyMod InitSnowLayers),
    !            each at 250 K holding dz*bdsno of ice, grain radius
    !            snw_rds_min; int_snow = h2osno; frac_sno by Niu and Yang
    !            (= 1 for a 4 m pack); frac_sno_eff 1 (non-soil column)
    !   surface  t_grnd = the top snow layer's 250 K; t_h2osfc 250 K
    ! kELM's snow index is reversed: slot 1 here is the layer on the ice.
    !
    implicit none
    type(kELMType), intent(in) :: elm
    integer, intent(in) :: logunit
    integer :: j, snl, k
    real(r8) :: snow_depth, h2osno, dzs(nlevsno), snowbd, fmelt, frac_sno
    real(r8), allocatable :: v(:,:), s(:)
    character(len=*), parameter :: subname = '(kelm_glacier_init) '

    allocate(v(n_glc, nlevgrnd))
    v = t_cold
    call glc_set(elm, 'glccol:thermal.t_soisno_soi', v)
    do j = 1, nlevgrnd; v(:, j) = dzsoi(j) * denice; end do
    call glc_set(elm, 'glccol:snowsoil.h2osoi_ice_soi', v)
    v = 0._r8
    call glc_set(elm, 'glccol:snowsoil.h2osoi_liq_soi', v)
    deallocate(v)

    h2osno     = h2osno_max
    snow_depth = h2osno / bdsno
    call init_snow_layers(snow_depth, snl, dzs)
    snowbd   = min(400._r8, h2osno / snow_depth)
    fmelt    = (snowbd / 100._r8)**1
    frac_sno = tanh(snow_depth / (2.5_r8 * zlnd * fmelt))

    allocate(v(n_glc, nlevsno))
    do k = 1, n_glc; v(k, :) = dzs; end do
    call glc_set(elm, 'glccol:snowsoil.dz_sno', v)
    do k = 1, n_glc; v(k, :) = dzs * bdsno; end do
    call glc_set(elm, 'glccol:snowsoil.h2osoi_ice_sno', v)
    v = 0._r8
    call glc_set(elm, 'glccol:snowsoil.h2osoi_liq_sno', v)
    v = 0._r8
    do j = 1, -snl; v(:, j) = t_cold; end do
    call glc_set(elm, 'glccol:snowsoil.t_soisno_sno', v)
    v = 0._r8
    do j = 1, -snl; v(:, j) = snw_rds_min; end do
    call glc_set(elm, 'glccol:snowsoil.snw_rds', v)
    deallocate(v)

    allocate(s(n_glc))
    s = real(snl, r8);  call glc_set(elm, 'glccol:snowsoil.snl', s)
    s = h2osno;         call glc_set(elm, 'glccol:snowsoil.h2osno', s)
    s = snow_depth;     call glc_set(elm, 'glccol:snowsoil.snow_depth', s)
    s = h2osno;         call glc_set(elm, 'glccol:snowsoil.int_snow', s)
    s = frac_sno;       call glc_set(elm, 'glccol:snowsoil.frac_sno', s)
    s = 1._r8;          call glc_set(elm, 'glccol:snowsoil.frac_sno_eff', s)
    s = t_cold;         call glc_set(elm, 'glccol:t_grnd', s)
    s = t_cold;         call glc_set(elm, 'glccol:t_h2osfc', s)
    deallocate(s)

    ! ELM SurfaceAlbedoType InitCold: the albedos and transmissions the first
    ! step's radiation uses, before any SurfaceAlbedo has run. Left at zero,
    ! a cold start whose first step is sunlit absorbs and reflects nothing
    ! (errsol = -solar for that step).
    allocate(v(n_glc, 2))
    v = 0.2_r8
    call glc_set(elm, 'glccol:albgrd', v); call glc_set(elm, 'glccol:albgri', v)
    call glc_set(elm, 'glccol:albsod', v); call glc_set(elm, 'glccol:albsoi', v)
    call glc_set(elm, 'glccol:albd', v);   call glc_set(elm, 'glccol:albi', v)
    v = 0.6_r8
    call glc_set(elm, 'glccol:albsnd_hst', v); call glc_set(elm, 'glccol:albsni_hst', v)
    v = 0._r8
    call glc_set(elm, 'glccol:fabd', v); call glc_set(elm, 'glccol:fabi', v)
    call glc_set(elm, 'glccol:ftid', v)
    v = 1._r8
    call glc_set(elm, 'glccol:ftdd', v); call glc_set(elm, 'glccol:ftii', v)
    deallocate(v)

    if (masterproc) then
       write(logunit,*) subname,'cold start: h2osno ',h2osno,' snow_depth ',snow_depth, &
            ' snl ',snl,' dz (ice side first) ',dzs(1:-snl)
       call shr_sys_flush(logunit)
    end if

  end subroutine seed_cold_start

  !-----------------------------------------------------------------------
  subroutine init_snow_layers(snow_depth, snl, dz)
    !
    ! ELM: SnowHydrologyMod :: InitSnowLayers, for nlevsno = 5. Returns snl
    ! and the layer thicknesses in kELM's reversed order (dz(1) on the ice).
    !
    implicit none
    real(r8), intent(in)  :: snow_depth
    integer , intent(out) :: snl
    real(r8), intent(out) :: dz(nlevsno)
    real(r8) :: dze(-nlevsno+1:0), maxbound
    integer  :: j

    dze = 0._r8
    if (snow_depth < dzmin1) then
       snl = 0
    else if (snow_depth <= dzmax_l1) then
       snl = -1
       dze(0) = snow_depth
    else
       snl = -2
       maxbound = sum(dzmax_u(1:2))
       do while (snow_depth > maxbound .and. -snl < nlevsno)
          snl = snl - 1
          maxbound = sum(dzmax_u(1:-snl))
       end do
       do j = 1, -snl - 2
          dze(j + snl) = dzmax_u(j)
       end do
       if (snow_depth <= sum(dzmax_u(1:-snl-2)) + 2._r8 * dzmax_u(-snl-1)) then
          dze(-1) = (snow_depth - sum(dzmax_u(1:-snl-2))) / 2._r8
          dze(0)  = dze(-1)
       else
          dze(-1) = dzmax_u(-snl-1)
          dze(0)  = snow_depth - sum(dzmax_u(1:-snl-1))
       end if
    end if
    ! ELM layer j (0 = on the ice) -> kELM slot 1-j.
    dz = 0._r8
    do j = 0, snl + 1, -1
       dz(1 - j) = dze(j)
    end do

  end subroutine init_snow_layers

  !-----------------------------------------------------------------------
  subroutine push_history_weights(elm, logunit, natural_pushed)
    !
    ! Per packed natural column (one per cell): the cell's glacier column and
    ! the glacier landunit's weight. Without lake and urban, nothing has
    ! pushed the natural weights the cell mean starts from; push them too.
    !
    implicit none
    type(kELMType), intent(in) :: elm
    integer, intent(in) :: logunit
    logical, intent(in) :: natural_pushed
    integer :: kc, k, l, g, ncell, ierr
    integer , allocatable :: glc_of(:), glc_of_cell(:), lake_of(:)
    real(r8), allocatable :: wg(:), wg_cell(:), wn(:), wn_cell(:), wl(:)
    character(len=*), parameter :: subname = '(kelm_glacier_init) '

    ncell = maxval(lun_gridcell(1:num_landunits))
    allocate(glc_of_cell(ncell), wg_cell(ncell), wn_cell(ncell))
    glc_of_cell = -1; wg_cell = 0._r8; wn_cell = 0._r8
    do l = 1, num_landunits
       g = lun_gridcell(l)
       if (lun_itype(l) == istice) wg_cell(g) = lun_wtgcell(l)
       if (lun_itype(l) == istsoil) wn_cell(g) = lun_wtgcell(l)
    end do
    do k = 1, n_glc
       glc_of_cell(cell_of_kglc(k)) = k - 1
    end do
    allocate(glc_of(n_kokkos_col), wg(n_kokkos_col), wn(n_kokkos_col), &
             wl(n_kokkos_col), lake_of(n_kokkos_col))
    do kc = 1, n_kokkos_col
       g = lun_gridcell(col_landunit(col_of_kcol(kc)))
       glc_of(kc) = glc_of_cell(g)
       wg(kc) = wg_cell(g)
       wn(kc) = wn_cell(g)
    end do

    if (.not. natural_pushed) then
       lake_of = -1
       wl = 0._r8
       call kELMHistorySetLandunitWeights(elm, lake_of, wn, wl, n_kokkos_col, ierr)
       if (ierr /= KELM_SUCCESS) call shr_sys_abort(subname//'ERROR: kELMHistorySetLandunitWeights failed')
    end if
    call kELMHistorySetGlacierWeights(elm, glc_of, wg, n_kokkos_col, ierr)
    if (ierr /= KELM_SUCCESS) call shr_sys_abort(subname//'ERROR: kELMHistorySetGlacierWeights failed')
    if (masterproc) then
       write(logunit,*) subname,'history weights, natural/glacier, first cell ',wn(1),wg(1)
       call shr_sys_flush(logunit)
    end if
    deallocate(glc_of_cell, wg_cell, wn_cell, glc_of, wg, wn, wl, lake_of)

  end subroutine push_history_weights

  !-----------------------------------------------------------------------
  ! Name-addressed seeding through the restart registry (kelmLakeMod's): a
  ! flat ROW-major buffer, buf((i-1)*n2 + j) = v(i,j), hence the transpose.
  !-----------------------------------------------------------------------
  integer function field_index(elm, name, n1, n2)
    implicit none
    type(kELMType), intent(in) :: elm
    character(len=*), intent(in) :: name
    integer, intent(in) :: n1, n2
    integer :: i, m1, m2, kind, ierr
    character(len=128) :: got
    call kELMRestartFieldFind(elm, name, i, ierr)
    if (ierr /= KELM_SUCCESS) call shr_sys_abort('(kelm_glacier) ERROR: no registry field '//trim(name))
    call kELMRestartFieldInfo(elm, i, got, m1, m2, kind, ierr)
    if (m1 /= n1 .or. m2 /= n2) then
       write(*,*) '(kelm_glacier) ',trim(name),' registry extents ',m1,m2,' given ',n1,n2
       call shr_sys_abort('(kelm_glacier) ERROR: extent mismatch on '//trim(name))
    end if
    field_index = i
  end function field_index

  subroutine glc_set_1d(elm, name, v)
    implicit none
    type(kELMType), intent(in) :: elm
    character(len=*), intent(in) :: name
    real(r8), intent(in) :: v(:)
    integer :: i, ierr
    i = field_index(elm, name, size(v), 0)
    call kELMRestartFieldSet(elm, i, v, size(v), ierr)
    if (ierr /= KELM_SUCCESS) call shr_sys_abort('(kelm_glacier) ERROR: set '//trim(name))
  end subroutine glc_set_1d

  subroutine glc_set_1i(elm, name, v)
    implicit none
    type(kELMType), intent(in) :: elm
    character(len=*), intent(in) :: name
    integer, intent(in) :: v(:)
    call glc_set_1d(elm, name, real(v, r8))
  end subroutine glc_set_1i

  subroutine glc_set_2d(elm, name, v)
    implicit none
    type(kELMType), intent(in) :: elm
    character(len=*), intent(in) :: name
    real(r8), intent(in) :: v(:,:)
    real(r8), allocatable :: flat(:)
    integer :: i, ierr
    i = field_index(elm, name, size(v,1), size(v,2))
    flat = reshape(transpose(v), (/ size(v) /))
    call kELMRestartFieldSet(elm, i, flat, size(flat), ierr)
    if (ierr /= KELM_SUCCESS) call shr_sys_abort('(kelm_glacier) ERROR: set '//trim(name))
  end subroutine glc_set_2d

  !-----------------------------------------------------------------------
  subroutine kelm_glacier_apply_finidat(elm, logunit)
    !
    ! Overwrite the glacier cold start with an ELM restart's glacier state
    ! (kelmFinidatMod), for identical-IC probes and runs from ELM initial
    ! conditions. Each packed column is gathered from its subgrid column, so
    ! like the natural and urban pushes this assumes one rank (the file is
    ! global, the subgrid indices local).
    !
    !   snow     layers (thickness, temperature, water, grain radius,
    !            aerosol), mass, depth, integrated snow, cover fractions,
    !            the previous step's per-layer freezing and melt
    !   ice      layer temperature and water (the split arrays; the combined
    !            ones are refreshed from them each step)
    !   surface  t_grnd, t_h2osfc; ground and patch albedos and the SNICAR
    !            absorption factors the first step's radiation reads
    ! The _p1 snow geometry is rebuilt from the snow every step; h2osoi_vol
    ! is diagnosed by HydrologyDrainage.
    !
    use kelmFinidatMod, only : finidat_read, fi_ncol, fi_npft, &
         fi_snl, fi_t_soisno, fi_h2osoi_liq, fi_h2osoi_ice, fi_dzsno, &
         fi_snw_rds, fi_qflx_snofrz_lyr, fi_snow_depth, fi_h2osno, &
         fi_int_snow, fi_frac_sno, fi_frac_sno_eff, fi_t_grnd, fi_t_h2osfc, &
         fi_qflx_snow_melt, fi_albgrd, fi_albgri, fi_albd, fi_albi, &
         fi_flx_absdv, fi_flx_absdn, fi_flx_absiv, fi_flx_absin, &
         fi_mss_bcphi, fi_mss_bcpho, fi_mss_ocphi, fi_mss_ocpho, &
         fi_mss_dst1, fi_mss_dst2, fi_mss_dst3, fi_mss_dst4
    implicit none
    type(kELMType), intent(in) :: elm
    integer, intent(in) :: logunit
    integer :: k, c, j
    real(r8), allocatable :: r1(:), rs(:,:), rg(:,:), rr(:,:), ra(:,:)
    character(len=*), parameter :: subname = '(kelm_glacier_apply_finidat) '

    if (.not. glacier_built) return
    if (.not. finidat_read) call shr_sys_abort(subname//'ERROR: finidat not read')
    if (fi_ncol /= num_columns .or. fi_npft /= num_patches) &
         call shr_sys_abort(subname//'ERROR: finidat subgrid is not this rank''s')
    allocate(r1(n_glc), rs(n_glc, nlevsno), rg(n_glc, nlevgrnd), rr(n_glc, 2), &
             ra(n_glc, nlevsno + 1))

    call col_1d('glccol:snowsoil.snl'          , real(fi_snl, r8))
    call col_1d('glccol:snowsoil.h2osno'       , fi_h2osno)
    call col_1d('glccol:snowsoil.snow_depth'   , fi_snow_depth)
    call col_1d('glccol:snowsoil.int_snow'     , fi_int_snow)
    call col_1d('glccol:snowsoil.frac_sno'     , fi_frac_sno)
    call col_1d('glccol:snowsoil.frac_sno_eff' , fi_frac_sno_eff)
    call col_1d('glccol:snowsoil.qflx_snow_melt', fi_qflx_snow_melt)
    call col_1d('glccol:t_grnd'                , fi_t_grnd)
    call col_1d('glccol:t_h2osfc'              , fi_t_h2osfc)

    ! Snow layers: ELM order (slot nlevsno on the ice) -> kELM reversed.
    call col_sno('glccol:snowsoil.dz_sno'        , fi_dzsno)
    call col_sno('glccol:snowsoil.t_soisno_sno'  , fi_t_soisno(1:nlevsno,:))
    call col_sno('glccol:snowsoil.h2osoi_liq_sno', fi_h2osoi_liq(1:nlevsno,:))
    call col_sno('glccol:snowsoil.h2osoi_ice_sno', fi_h2osoi_ice(1:nlevsno,:))
    call col_sno('glccol:snowsoil.snw_rds'       , fi_snw_rds)
    call col_sno('glccol:snowsoil.mss_bcphi', fi_mss_bcphi)
    call col_sno('glccol:snowsoil.mss_bcpho', fi_mss_bcpho)
    call col_sno('glccol:snowsoil.mss_ocphi', fi_mss_ocphi)
    call col_sno('glccol:snowsoil.mss_ocpho', fi_mss_ocpho)
    call col_sno('glccol:snowsoil.mss_dst1' , fi_mss_dst1)
    call col_sno('glccol:snowsoil.mss_dst2' , fi_mss_dst2)
    call col_sno('glccol:snowsoil.mss_dst3' , fi_mss_dst3)
    call col_sno('glccol:snowsoil.mss_dst4' , fi_mss_dst4)
    ! ELM-ordered on the kELM side too (as the natural column's).
    do k = 1, n_glc; rs(k,:) = fi_qflx_snofrz_lyr(:, col_of_kglc(k)); end do
    call glc_set(elm, 'glccol:qflx_snofrz_lyr', rs)

    ! Ice layers.
    do k = 1, n_glc; rg(k,:) = fi_t_soisno(nlevsno+1:nlevsno+nlevgrnd, col_of_kglc(k)); end do
    call glc_set(elm, 'glccol:thermal.t_soisno_soi', rg)
    do k = 1, n_glc; rg(k,:) = fi_h2osoi_liq(nlevsno+1:nlevsno+nlevgrnd, col_of_kglc(k)); end do
    call glc_set(elm, 'glccol:snowsoil.h2osoi_liq_soi', rg)
    do k = 1, n_glc; rg(k,:) = fi_h2osoi_ice(nlevsno+1:nlevsno+nlevgrnd, col_of_kglc(k)); end do
    call glc_set(elm, 'glccol:snowsoil.h2osoi_ice_soi', rg)

    ! Radiation: ground (column), surface (the one patch), SNICAR factors.
    do k = 1, n_glc; rr(k,:) = fi_albgrd(:, col_of_kglc(k)); end do
    call glc_set(elm, 'glccol:albgrd', rr)
    do k = 1, n_glc; rr(k,:) = fi_albgri(:, col_of_kglc(k)); end do
    call glc_set(elm, 'glccol:albgri', rr)
    do k = 1, n_glc; rr(k,:) = fi_albd(:, patch_of_kglc(k)); end do
    call glc_set(elm, 'glccol:albd', rr)
    do k = 1, n_glc; rr(k,:) = fi_albi(:, patch_of_kglc(k)); end do
    call glc_set(elm, 'glccol:albi', rr)
    do k = 1, n_glc; ra(k,:) = fi_flx_absdv(:, col_of_kglc(k)); end do
    call glc_set(elm, 'glccol:flx_absdv', ra)
    do k = 1, n_glc; ra(k,:) = fi_flx_absdn(:, col_of_kglc(k)); end do
    call glc_set(elm, 'glccol:flx_absdn', ra)
    do k = 1, n_glc; ra(k,:) = fi_flx_absiv(:, col_of_kglc(k)); end do
    call glc_set(elm, 'glccol:flx_absiv', ra)
    do k = 1, n_glc; ra(k,:) = fi_flx_absin(:, col_of_kglc(k)); end do
    call glc_set(elm, 'glccol:flx_absin', ra)

    if (masterproc) then
       c = col_of_kglc(1)
       write(logunit,*) subname,'seeded ',n_glc,' glacier columns from finidat; first: snl ', &
            fi_snl(c),' h2osno ',fi_h2osno(c)
       call shr_sys_flush(logunit)
    end if
    deallocate(r1, rs, rg, rr, ra)

  contains

    subroutine col_1d(name, src)
      character(len=*), intent(in) :: name
      real(r8), intent(in) :: src(:)
      integer :: kk
      do kk = 1, n_glc; r1(kk) = src(col_of_kglc(kk)); end do
      call glc_set(elm, name, r1)
    end subroutine col_1d

    subroutine col_sno(name, src)
      ! src (nlevsno, ncol), ELM order (slot nlevsno on the ice).
      character(len=*), intent(in) :: name
      real(r8), intent(in) :: src(:,:)
      integer :: kk, jj
      do kk = 1, n_glc
         do jj = 1, nlevsno
            rs(kk, jj) = src(nlevsno+1-jj, col_of_kglc(kk))
         end do
      end do
      call glc_set(elm, name, rs)
    end subroutine col_sno

  end subroutine kelm_glacier_apply_finidat

  !-----------------------------------------------------------------------
  subroutine kelm_glacier_diag_step(elm, tag)
    !
    ! Per-step glacier state into the KELM_DIAG trace, packed glacier order
    ! (glacier column k is ELM's glacier column on a one-cell domain).
    ! Snow arrays as kELM holds them: slot 1 is the layer on the ice. tag:
    ! 'kelm_glcin' at the top of the step (pairs with ELM canhydro_in,
    ! snowlayer_in, surfrad_in), 'kelm_glcout' after the fluxes.
    !
    implicit none
    type(kELMType), intent(in) :: elm
    character(len=*), intent(in) :: tag
    real(r8), allocatable :: v(:), v5(:,:), v6(:,:), v15(:,:), v2(:,:)
    integer :: j
    character(len=24), parameter :: c1(14) = (/ &
         'snowsoil.snl            ', 'snowsoil.h2osno         ', 'snowsoil.snow_depth     ', &
         'snowsoil.frac_sno       ', 'snowsoil.int_snow       ', 't_grnd                  ', &
         'snowsoil.do_capsnow     ', 'qflx_snwcp_ice_col      ', 'sabg                    ', &
         'eflx_sh_tot             ', 'eflx_lh_tot             ', 'eflx_soil_grnd          ', &
         'qflx_evap_soi           ', 'snowsoil.qflx_snow_grnd ' /)
    character(len=24), parameter :: c5(7) = (/ &
         'snowsoil.dz_sno         ', 'snowsoil.t_soisno_sno   ', 'snowsoil.h2osoi_ice_sno ', &
         'snowsoil.h2osoi_liq_sno ', 'snowsoil.snw_rds        ', 'snowsoil.mss_bcphi      ', &
         'snowsoil.mss_dst1       ' /)
    character(len=24), parameter :: c6(2) = (/ 'flx_absdv               ', 'flx_absiv               ' /)
    character(len=24), parameter :: c2(4) = (/ 'albgrd                  ', 'albgri                  ', &
         'albd                    ', 'albi                    ' /)

    if (.not. kelm_diag_enabled .or. .not. glacier_built) return
    allocate(v(n_glc), v5(n_glc, nlevsno), v6(n_glc, nlevsno + 1), v15(n_glc, nlevgrnd), v2(n_glc, 2))
    do j = 1, size(c1)
       call glc_get_1d(elm, 'glccol:'//trim(c1(j)), v)
       call kelm_diag_1d(trim(tag)//':'//trim(c1(j)), v, n_glc)
    end do
    do j = 1, size(c5)
       call glc_get_2d(elm, 'glccol:'//trim(c5(j)), v5)
       call kelm_diag_2d(trim(tag)//':'//trim(c5(j)), v5, n_glc, nlevsno)
    end do
    do j = 1, size(c6)
       call glc_get_2d(elm, 'glccol:'//trim(c6(j)), v6)
       call kelm_diag_2d(trim(tag)//':'//trim(c6(j)), v6, n_glc, nlevsno + 1)
    end do
    do j = 1, size(c2)
       call glc_get_2d(elm, 'glccol:'//trim(c2(j)), v2)
       call kelm_diag_2d(trim(tag)//':'//trim(c2(j)), v2, n_glc, 2)
    end do
    call glc_get_2d(elm, 'glccol:thermal.t_soisno_soi', v15)
    call kelm_diag_2d(trim(tag)//':t_soisno_soi', v15, n_glc, nlevgrnd)
    call glc_get_2d(elm, 'glccol:snowsoil.h2osoi_liq_soi', v15)
    call kelm_diag_2d(trim(tag)//':h2osoi_liq_soi', v15, n_glc, nlevgrnd)
    deallocate(v, v5, v6, v15, v2)

  end subroutine kelm_glacier_diag_step

  subroutine glc_get_1d(elm, name, v)
    implicit none
    type(kELMType), intent(in) :: elm
    character(len=*), intent(in) :: name
    real(r8), intent(inout) :: v(:)
    integer :: i, ierr
    i = field_index(elm, name, size(v), 0)
    call kELMRestartFieldGet(elm, i, v, size(v), ierr)
    if (ierr /= KELM_SUCCESS) call shr_sys_abort('(kelm_glacier) ERROR: get '//trim(name))
  end subroutine glc_get_1d

  subroutine glc_get_2d(elm, name, v)
    implicit none
    type(kELMType), intent(in) :: elm
    character(len=*), intent(in) :: name
    real(r8), intent(inout) :: v(:,:)
    real(r8), allocatable :: flat(:)
    integer :: i, ierr
    i = field_index(elm, name, size(v,1), size(v,2))
    allocate(flat(size(v)))
    call kELMRestartFieldGet(elm, i, flat, size(flat), ierr)
    if (ierr /= KELM_SUCCESS) call shr_sys_abort('(kelm_glacier) ERROR: get '//trim(name))
    v = transpose(reshape(flat, (/ size(v,2), size(v,1) /)))
  end subroutine glc_get_2d

  !-----------------------------------------------------------------------
  subroutine kelm_glacier_clean()
    implicit none
    if (associated(col_of_kglc))   deallocate(col_of_kglc)
    if (associated(patch_of_kglc)) deallocate(patch_of_kglc)
    if (associated(cell_of_kglc))  deallocate(cell_of_kglc)
    n_glc = 0
    glacier_built = .false.
  end subroutine kelm_glacier_clean

end module kelmGlacierMod
