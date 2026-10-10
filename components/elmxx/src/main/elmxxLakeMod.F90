module elmxxLakeMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! The lake landunit on the coupled path (plan Stage 6.5, L1 and L4): the
  ! Fortran side of ELMxx's lake surface (lakeCol/lakePatch in the C++).
  !
  ! WHAT THIS OWNS
  !   maps        packed lake column/patch <-> subgrid column/patch, the cell
  !               of each, and the natural column of that cell (the lake reads
  !               its forcing from there, on the device).
  !   init        ELM's lake geometry (initVerticalMod), the soil properties
  !               under the lake (SoilStateType, lake block), and ELM's lake
  !               cold start (ColumnDataType, LakeStateType, FrictionVelocity
  !               InitCold) -- computed here, seeded once by NAME through the
  !               restart registry ("lakecol:<field>"), no per-field setters.
  !   per step    the reference height and shortwave per lake patch, pushed
  !               the way the natural ones are.
  !   history     each cell's natural and lake landunit weights, so h0 is
  !               ELM's landunit-weighted cell mean.
  ! ELM is the reference for behaviour only; nothing is copied.
  !
  ! ONE PATCH PER LAKE COLUMN, packed in the same ascending order, so packed
  ! lake patch k sits on packed lake column k. LakeHydrology indexes p = c,
  ! and ELMxxAllocateLakeSurface refuses anything else.
  !-----------------------------------------------------------------------

  use shr_kind_mod    , only : r8 => shr_kind_r8
  use shr_sys_mod     , only : shr_sys_abort, shr_sys_flush
  use elmxxSpmdMod    , only : masterproc, iam
  use elmxxSubgridMod , only : num_landunits, num_columns, num_patches, &
                               lun_gridcell, lun_itype, lun_wtgcell, &
                               col_landunit, patch_column, istsoil, istdlak, &
                               isturb_tbd, isturb_md
  use elmxxSurfdataMod, only : lakedepth_in, etalake_in, lakefetch_in, lake_spval
  use elmxxSurfaceStateMod, only : col_pct_sand, col_pct_clay, col_organic
  use elmxxSoilPropMod, only : nlevsoi, nlevgrnd, nlevsno, nlevtot, &
                               zsoi, dzsoi, zisoi
  use elmxxKokkosStateMod, only : kcol_of_col, col_of_kcol, n_kokkos_col
  use elmxxForcingMod , only : forc_z, forc_swvdr, forc_swndr, forc_swvdf, forc_swndf
  use elmxx_mod       , only : ELMxxType, ELMXX_SUCCESS, &
                               ELMxxAllocateLakeSurface, ELMxxRestartFieldFind, &
                               ELMxxRestartFieldSet, ELMxxRestartFieldGet, &
                               ELMxxRestartFieldInfo, &
                               ELMxxHistorySetLandunitWeights

  implicit none
  save
  private

  integer, parameter, public :: nlevlak = 10   ! use_extralakelayers off

  integer, public :: n_lake = 0                 ! packed lake columns = patches
  logical, public :: lake_built = .false.
  integer, public, pointer :: col_of_klake(:)   => null()  ! (n_lake) -> subgrid column
  integer, public, pointer :: patch_of_klake(:) => null()  ! (n_lake) -> subgrid patch
  integer, public, pointer :: cell_of_klake(:)  => null()  ! (n_lake) -> 1-based local cell

  public :: elmxx_lake_init
  public :: elmxx_lake_push_forcing
  public :: elmxx_lake_get
  public :: elmxx_lake_clean

  interface lake_set
     module procedure lake_set_1d, lake_set_2d, lake_set_1i
  end interface lake_set

contains

  !-----------------------------------------------------------------------
  subroutine elmxx_lake_init(elm, logunit)
    !
    ! Build the lake maps, compute the lake's static fields and cold start,
    ! allocate the lake surface and seed it, and give history the landunit
    ! weights. A no-op on a domain without lake.
    !
    implicit none
    type(ELMxxType), intent(in) :: elm
    integer, intent(in) :: logunit
    integer :: ierr
    character(len=*), parameter :: subname = '(elmxx_lake_init) '

    call elmxx_lake_clean()
    call build_maps(logunit)
    if (n_lake == 0) return

    call ELMxxAllocateLakeSurface(elm, n_lake, n_lake, ierr)
    if (ierr /= ELMXX_SUCCESS) call shr_sys_abort(subname//'ERROR: ELMxxAllocateLakeSurface failed')

    call seed_topology(elm)
    call seed_geometry(elm, logunit)
    call seed_soil_properties(elm, logunit)
    call seed_cold_start(elm, logunit)
    call push_history_weights(elm, logunit)

    lake_built = .true.

  end subroutine elmxx_lake_init

  !-----------------------------------------------------------------------
  subroutine build_maps(logunit)
    implicit none
    integer, intent(in) :: logunit
    integer :: c, p, k, np
    character(len=*), parameter :: subname = '(elmxx_lake_init) '

    n_lake = 0
    do c = 1, num_columns
       if (lun_itype(col_landunit(c)) == istdlak) n_lake = n_lake + 1
    end do
    if (masterproc) write(logunit,*) subname,'rank ',iam,' lake columns ',n_lake
    if (n_lake == 0) return

    allocate(col_of_klake(n_lake), patch_of_klake(n_lake), cell_of_klake(n_lake))
    k = 0
    do c = 1, num_columns
       if (lun_itype(col_landunit(c)) /= istdlak) cycle
       k = k + 1
       col_of_klake(k)  = c
       cell_of_klake(k) = lun_gridcell(col_landunit(c))
       np = 0
       do p = 1, num_patches
          if (patch_column(p) == c) then
             np = np + 1
             patch_of_klake(k) = p
          end if
       end do
       if (np /= 1) call shr_sys_abort(subname//'ERROR: a lake column must carry exactly one patch')
    end do
    ! Ascending patches on ascending columns: packed patch k is on packed
    ! column k, which is what LakeHydrology's p = c assumes.
    do k = 2, n_lake
       if (patch_of_klake(k) <= patch_of_klake(k-1)) &
            call shr_sys_abort(subname//'ERROR: lake patches are not packed like lake columns')
    end do

  end subroutine build_maps


  !-----------------------------------------------------------------------
  subroutine seed_topology(elm)
    implicit none
    type(ELMxxType), intent(in) :: elm
    integer :: k
    integer, allocatable :: ibuf(:)

    allocate(ibuf(n_lake))
    do k = 1, n_lake
       ibuf(k) = cell_of_klake(k) - 1
    end do
    call lake_set(elm, 'lakecol:col_gridcell', ibuf)
    do k = 1, n_lake
       ibuf(k) = k - 1
    end do
    call lake_set(elm, 'lakepatch:patch_column', ibuf)
    deallocate(ibuf)

  end subroutine seed_topology

  !-----------------------------------------------------------------------
  subroutine seed_geometry(elm, logunit)
    !
    ! ELM initVerticalMod. Lake layers from LAKEDEPTH by its three branches,
    ! the standard 10 layers rescaled; the ground column under the lake is the
    ! same exponential grid as every other column; no snow at a cold start.
    ! Layout is the lake surface's (LakeColumnData.h): NLEVTOT slot m <-> ELM
    ! layer m - nlevsno (1-based here), zi slot i <-> ELM zi(i - 1 - nlevsno).
    !
    implicit none
    type(ELMxxType), intent(in) :: elm
    integer, intent(in) :: logunit
    integer  :: k, g, j
    real(r8) :: depth, ratio
    real(r8), allocatable :: depth_k(:), dz_lake(:,:), z_lake(:,:)
    real(r8), allocatable :: dz(:,:), z(:,:), zi(:,:), etal(:), fetch(:)
    real(r8), parameter :: dzlak(nlevlak) = (/ 0.1_r8, 1._r8, 2._r8, 3._r8, &
         4._r8, 5._r8, 7._r8, 7._r8, 10.45_r8, 10.45_r8 /)
    real(r8), parameter :: zlak(nlevlak) = (/ 0.05_r8, 0.6_r8, 2.1_r8, 4.6_r8, &
         8.1_r8, 12.6_r8, 18.6_r8, 25.6_r8, 34.325_r8, 44.775_r8 /)
    character(len=*), parameter :: subname = '(elmxx_lake_init) '

    allocate(depth_k(n_lake), dz_lake(n_lake,nlevlak), z_lake(n_lake,nlevlak), &
             dz(n_lake,nlevtot), z(n_lake,nlevtot), zi(n_lake,nlevtot+1), &
             etal(n_lake), fetch(n_lake))

    do k = 1, n_lake
       g = cell_of_klake(k)
       depth = lakedepth_in(g)
       if (depth == lake_spval) then
          ! No LAKEDEPTH on the file: ELM's standard 50 m lake.
          depth = zlak(nlevlak) + 0.5_r8*dzlak(nlevlak)
          z_lake(k,:)  = zlak
          dz_lake(k,:) = dzlak
       else if (depth > 1._r8 .and. depth < 5000._r8) then
          ratio = depth / (zlak(nlevlak) + 0.5_r8*dzlak(nlevlak))
          z_lake(k,1)  = zlak(1)
          dz_lake(k,1) = dzlak(1)
          dz_lake(k,2:nlevlak-1) = dzlak(2:nlevlak-1)*ratio
          dz_lake(k,nlevlak) = dzlak(nlevlak)*ratio - (dz_lake(k,1) - dzlak(1)*ratio)
          do j = 2, nlevlak
             z_lake(k,j) = z_lake(k,j-1) + (dz_lake(k,j-1) + dz_lake(k,j))/2._r8
          end do
       else if (depth > 0._r8 .and. depth <= 1._r8) then
          dz_lake(k,:) = depth / nlevlak
          z_lake(k,1)  = dz_lake(k,1) / 2._r8
          do j = 2, nlevlak
             z_lake(k,j) = z_lake(k,j-1) + (dz_lake(k,j-1) + dz_lake(k,j))/2._r8
          end do
       else
          write(logunit,*) subname,'ERROR: bad lake depth ',depth,' in local cell ',g
          call shr_sys_abort(subname//'ERROR: bad lake depth')
       end if
       depth_k(k) = depth
       etal(k)  = etalake_in(g)
       fetch(k) = lakefetch_in(g)
    end do

    dz = 0._r8; z = 0._r8; zi = 0._r8
    do j = 1, nlevgrnd
       dz(:, nlevsno + j) = dzsoi(j)
       z (:, nlevsno + j) = zsoi(j)
    end do
    do j = 0, nlevgrnd
       zi(:, nlevsno + 1 + j) = zisoi(j)
    end do

    call lake_set(elm, 'lakecol:lakedepth', depth_k)
    call lake_set(elm, 'lakecol:etal', etal)
    call lake_set(elm, 'lakecol:lakefetch', fetch)
    call lake_set(elm, 'lakecol:dz_lake', dz_lake)
    call lake_set(elm, 'lakecol:z_lake', z_lake)
    call lake_set(elm, 'lakecol:dz', dz)
    call lake_set(elm, 'lakecol:z', z)
    call lake_set(elm, 'lakecol:zi', zi)

    if (masterproc) then
       write(logunit,*) subname,'lake depth [m] ',depth_k(1),' etal ',etal(1),' fetch ',fetch(1)
       write(logunit,*) subname,'dz_lake ',dz_lake(1,:)
       write(logunit,*) subname,'z_lake  ',z_lake(1,:)
       call shr_sys_flush(logunit)
    end if
    deallocate(depth_k, dz_lake, z_lake, dz, z, zi, etal, fetch)

  end subroutine seed_geometry

  !-----------------------------------------------------------------------
  subroutine seed_soil_properties(elm, logunit)
    !
    ! ELM SoilStateType, "Set soil hydraulic and thermal properties: lake".
    ! Not the natural block: the organic end members are constants
    ! (om_watsat_lake 0.9), bulk density is from the mineral porosity, and
    ! om_frac is ((organic/organic_max)**2 / organic_max * organic_max)**2 --
    ! ELM squares it once storing cellorg_col and again reading it back, so
    ! the lake's organic fraction is the fourth power. ELM's, kept.
    ! Below nlevsoi the deepest layer's texture, no organic, and bedrock csol.
    ! Only what the lake kernels read is computed: watsat, tkmg, tksatu,
    ! tkdry, csol.
    !
    implicit none
    type(ELMxxType), intent(in) :: elm
    integer, intent(in) :: logunit
    integer  :: k, c, lev, jsrc
    real(r8) :: sand, clay, om_frac, cellorg, wsat, bd, tkm
    real(r8), allocatable :: watsat(:,:), tkmg(:,:), tksatu(:,:), tkdry(:,:), csol(:,:)
    real(r8), parameter :: organic_max    = 130.0_r8  ! clm_params, kg/m3
    real(r8), parameter :: om_watsat_lake = 0.9_r8
    real(r8), parameter :: om_tkm         = 0.25_r8
    real(r8), parameter :: om_tkd         = 0.05_r8
    real(r8), parameter :: om_csol        = 2.5_r8
    real(r8), parameter :: csol_bedrock   = 2.0e6_r8
    character(len=*), parameter :: subname = '(elmxx_lake_init) '

    allocate(watsat(n_lake,nlevgrnd), tkmg(n_lake,nlevgrnd), tksatu(n_lake,nlevgrnd), &
             tkdry(n_lake,nlevgrnd), csol(n_lake,nlevgrnd))

    do k = 1, n_lake
       c = col_of_klake(k)
       do lev = 1, nlevgrnd
          jsrc = min(lev, nlevsoi)
          sand = col_pct_sand(c, jsrc)
          clay = col_pct_clay(c, jsrc)
          if (lev <= nlevsoi) then
             cellorg = (col_organic(c, lev)/organic_max)**2._r8 * organic_max
             om_frac = (cellorg/organic_max)**2._r8
          else
             om_frac = 0.0_r8
          end if
          if (sand + clay <= 0.0_r8) call shr_sys_abort(subname//'ERROR: lake soil has neither sand nor clay')

          wsat = 0.489_r8 - 0.00126_r8*sand
          bd   = (1._r8 - wsat)*2.7e3_r8
          watsat(k,lev) = (1._r8 - om_frac)*wsat + om_watsat_lake*om_frac
          tkm  = (1._r8 - om_frac)*(8.80_r8*sand + 2.92_r8*clay)/(sand + clay) + om_tkm*om_frac
          tkmg(k,lev)   = tkm ** (1._r8 - watsat(k,lev))
          tksatu(k,lev) = tkmg(k,lev)*0.57_r8**watsat(k,lev)
          tkdry(k,lev)  = ((0.135_r8*bd + 64.7_r8) / (2.7e3_r8 - 0.947_r8*bd))*(1._r8 - om_frac) &
                        + om_tkd*om_frac
          csol(k,lev)   = ((1._r8 - om_frac)*(2.128_r8*sand + 2.385_r8*clay)/(sand + clay) &
                        + om_csol*om_frac)*1.e6_r8
          if (lev > nlevsoi) csol(k,lev) = csol_bedrock
       end do
    end do

    call lake_set(elm, 'lakecol:watsat', watsat)
    call lake_set(elm, 'lakecol:tkmg', tkmg)
    call lake_set(elm, 'lakecol:tksatu', tksatu)
    call lake_set(elm, 'lakecol:tkdry', tkdry)
    call lake_set(elm, 'lakecol:csol', csol)

    if (masterproc) then
       write(logunit,*) subname,'lake soil watsat ',watsat(1,:)
       write(logunit,*) subname,'lake soil tkdry  ',tkdry(1,:)
       call shr_sys_flush(logunit)
    end if
    deallocate(watsat, tkmg, tksatu, tkdry, csol)

  end subroutine seed_soil_properties

  !-----------------------------------------------------------------------
  subroutine seed_cold_start(elm, logunit)
    !
    ! ELM's lake cold start:
    !   ColumnDataType InitCold   t_soisno, t_lake, t_grnd = 277 K; soil
    !                             h2osoi_vol = watsat above bedrock, 0 below,
    !                             liquid at 277 K; no snow
    !   LakeStateType InitCold    lake_icefrac 0, savedtke1 = tkwat 0.57,
    !                             ust_lake 0.1
    !   FrictionVelocity InitCold z0mg 0.0004 m on lake columns
    ! The soil properties are recomputed here rather than read back, from the
    ! same formula, so this does not depend on seeding order.
    !
    implicit none
    type(ELMxxType), intent(in) :: elm
    integer, intent(in) :: logunit
    integer  :: k, j, m
    real(r8), allocatable :: t(:,:), liq(:,:), ice(:,:), vol(:,:), watsat(:)
    real(r8), parameter :: t_lake_cold = 277.0_r8
    real(r8), parameter :: denh2o = 1000.0_r8
    character(len=*), parameter :: subname = '(elmxx_lake_init) '

    allocate(t(n_lake,nlevtot), liq(n_lake,nlevtot), ice(n_lake,nlevtot), &
             vol(n_lake,nlevgrnd), watsat(nlevgrnd))
    t = 0._r8; liq = 0._r8; ice = 0._r8; vol = 0._r8

    call elmxx_lake_get(elm, 'lakecol:watsat', vol)   ! (n_lake, nlevgrnd)
    do k = 1, n_lake
       watsat = vol(k,:)
       do j = 1, nlevgrnd
          m = nlevsno + j
          t(k,m) = t_lake_cold
          if (j <= nlevsoi) then
             vol(k,j) = watsat(j)
          else
             vol(k,j) = 0._r8                  ! bedrock
          end if
          liq(k,m) = dzsoi(j)*denh2o*vol(k,j)  ! 277 K: all liquid
       end do
    end do

    call lake_set(elm, 'lakecol:t_soisno', t)
    call lake_set(elm, 'lakecol:h2osoi_liq', liq)
    call lake_set(elm, 'lakecol:h2osoi_ice', ice)
    call lake_set(elm, 'lakecol:h2osoi_vol', vol)
    call lake_set(elm, 'lakecol:t_lake', spread(spread(t_lake_cold, 1, n_lake), 2, nlevlak))
    call lake_set(elm, 'lakecol:lake_icefrac', spread(spread(0._r8, 1, n_lake), 2, nlevlak))
    call lake_set(elm, 'lakecol:t_grnd', spread(t_lake_cold, 1, n_lake))
    call lake_set(elm, 'lakecol:savedtke1', spread(0.57_r8, 1, n_lake))
    call lake_set(elm, 'lakecol:ust_lake', spread(0.1_r8, 1, n_lake))
    call lake_set(elm, 'lakecol:z0mg', spread(0.0004_r8, 1, n_lake))

    if (masterproc) then
       write(logunit,*) subname,'lake cold start: 277 K, soil saturated, ', &
            'h2osoi_liq(top) ',liq(1,nlevsno+1),' kg/m2'
       call shr_sys_flush(logunit)
    end if
    deallocate(t, liq, ice, vol, watsat)

  end subroutine seed_cold_start

  !-----------------------------------------------------------------------
  subroutine push_history_weights(elm, logunit)
    !
    ! Per packed natural column (one per cell): the cell's lake column and the
    ! natural and lake landunit weights. The cell mean history forms is over
    ! natural, lake and urban (elmxxUrbanMod pushes urban's), so any other
    ! landunit with weight -- wetland, glacier, crop, none of which ELMxx
    ! models -- would make it wrong; refuse that rather than write a
    ! mislabelled h0.
    !
    implicit none
    type(ELMxxType), intent(in) :: elm
    integer, intent(in) :: logunit
    integer :: kc, k, l, g, ierr
    integer , allocatable :: lake_of(:)
    real(r8), allocatable :: wn(:), wl(:), wother(:)
    real(r8), parameter :: tol = 1.0e-6_r8
    character(len=*), parameter :: subname = '(elmxx_lake_init) '

    allocate(lake_of(n_kokkos_col), wn(n_kokkos_col), wl(n_kokkos_col), &
             wother(n_kokkos_col))
    lake_of = -1; wn = 0._r8; wl = 0._r8; wother = 0._r8
    do kc = 1, n_kokkos_col
       g = lun_gridcell(col_landunit(col_of_kcol(kc)))
       do l = 1, num_landunits
          if (lun_gridcell(l) /= g) cycle
          select case (lun_itype(l))
          case (istsoil)
             wn(kc) = lun_wtgcell(l)
          case (istdlak)
             wl(kc) = lun_wtgcell(l)
          case (isturb_tbd:isturb_md)
             ! Urban joins the mean through ELMxxHistorySetUrbanWeights
             ! (elmxxUrbanMod, Stage 6.7).
          case default
             wother(kc) = wother(kc) + lun_wtgcell(l)
          end select
       end do
       do k = 1, n_lake
          if (cell_of_klake(k) == g) lake_of(kc) = k - 1
       end do
    end do
    if (any(wother > tol)) then
       call shr_sys_abort(subname//'ERROR: a lake cell has weight on a landunit '// &
            'other than natural, lake and urban, which ELMxx history cannot represent')
    end if

    call ELMxxHistorySetLandunitWeights(elm, lake_of, wn, wl, n_kokkos_col, ierr)
    if (ierr /= ELMXX_SUCCESS) call shr_sys_abort(subname//'ERROR: ELMxxHistorySetLandunitWeights failed')
    if (masterproc) then
       write(logunit,*) subname,'history weights, natural/lake, first cell ',wn(1),wl(1)
       call shr_sys_flush(logunit)
    end if
    deallocate(lake_of, wn, wl, wother)

  end subroutine push_history_weights

  !-----------------------------------------------------------------------
  subroutine elmxx_lake_push_forcing(elm)
    !
    ! Per step: the reference height and the incident shortwave, per lake
    ! patch, with ELM's band assignment (1 visible, 2 near-IR). Everything else
    ! the lake reads, ELMxxLakeGatherForcing takes on the device from the
    ! natural column of the same cell.
    !
    implicit none
    type(ELMxxType), intent(in) :: elm
    integer :: k, g
    real(r8), allocatable :: hgt(:), sold(:,:), soli(:,:)

    if (.not. lake_built) return
    allocate(hgt(n_lake), sold(n_lake,2), soli(n_lake,2))
    do k = 1, n_lake
       g = cell_of_klake(k)
       hgt(k)    = forc_z(g)
       sold(k,1) = forc_swvdr(g); sold(k,2) = forc_swndr(g)
       soli(k,1) = forc_swvdf(g); soli(k,2) = forc_swndf(g)
    end do
    call lake_set(elm, 'lakepatch:forc_hgt', hgt)
    call lake_set(elm, 'lakepatch:forc_solad', sold)
    call lake_set(elm, 'lakepatch:forc_solai', soli)
    deallocate(hgt, sold, soli)

  end subroutine elmxx_lake_push_forcing

  !-----------------------------------------------------------------------
  ! Name-addressed seeding through the restart registry. The registry crosses
  ! a flat ROW-major buffer, buf((i-1)*n2 + j) = v(i,j); a Fortran (n1,n2)
  ! array is column-major, hence the transpose.
  !-----------------------------------------------------------------------
  integer function field_index(elm, name, n1, n2)
    implicit none
    type(ELMxxType), intent(in) :: elm
    character(len=*), intent(in) :: name
    integer, intent(in) :: n1, n2
    integer :: i, m1, m2, kind, ierr
    character(len=128) :: got
    call ELMxxRestartFieldFind(elm, name, i, ierr)
    if (ierr /= ELMXX_SUCCESS) call shr_sys_abort('(elmxx_lake) ERROR: no registry field '//trim(name))
    call ELMxxRestartFieldInfo(elm, i, got, m1, m2, kind, ierr)
    if (m1 /= n1 .or. m2 /= n2) then
       write(*,*) '(elmxx_lake) ',trim(name),' registry extents ',m1,m2,' given ',n1,n2
       call shr_sys_abort('(elmxx_lake) ERROR: extent mismatch on '//trim(name))
    end if
    field_index = i
  end function field_index

  subroutine lake_set_1d(elm, name, v)
    implicit none
    type(ELMxxType), intent(in) :: elm
    character(len=*), intent(in) :: name
    real(r8), intent(in) :: v(:)
    integer :: i, ierr
    i = field_index(elm, name, size(v), 0)
    call ELMxxRestartFieldSet(elm, i, v, size(v), ierr)
    if (ierr /= ELMXX_SUCCESS) call shr_sys_abort('(elmxx_lake) ERROR: set '//trim(name))
  end subroutine lake_set_1d

  subroutine lake_set_1i(elm, name, v)
    implicit none
    type(ELMxxType), intent(in) :: elm
    character(len=*), intent(in) :: name
    integer, intent(in) :: v(:)
    call lake_set_1d(elm, name, real(v, r8))
  end subroutine lake_set_1i

  subroutine lake_set_2d(elm, name, v)
    implicit none
    type(ELMxxType), intent(in) :: elm
    character(len=*), intent(in) :: name
    real(r8), intent(in) :: v(:,:)
    real(r8), allocatable :: flat(:)
    integer :: i, ierr
    i = field_index(elm, name, size(v,1), size(v,2))
    flat = reshape(transpose(v), (/ size(v) /))
    call ELMxxRestartFieldSet(elm, i, flat, size(flat), ierr)
    if (ierr /= ELMXX_SUCCESS) call shr_sys_abort('(elmxx_lake) ERROR: set '//trim(name))
  end subroutine lake_set_2d

  !-----------------------------------------------------------------------
  subroutine elmxx_lake_get(elm, name, v)
    ! Read a 2-D lake field into a Fortran (n1,n2) array.
    implicit none
    type(ELMxxType), intent(in) :: elm
    character(len=*), intent(in) :: name
    real(r8), intent(inout) :: v(:,:)
    real(r8), allocatable :: flat(:)
    integer :: i, ierr
    i = field_index(elm, name, size(v,1), size(v,2))
    allocate(flat(size(v)))
    call ELMxxRestartFieldGet(elm, i, flat, size(flat), ierr)
    if (ierr /= ELMXX_SUCCESS) call shr_sys_abort('(elmxx_lake) ERROR: get '//trim(name))
    v = transpose(reshape(flat, (/ size(v,2), size(v,1) /)))
  end subroutine elmxx_lake_get

  !-----------------------------------------------------------------------
  subroutine elmxx_lake_clean()
    implicit none
    if (associated(col_of_klake))   deallocate(col_of_klake)
    if (associated(patch_of_klake)) deallocate(patch_of_klake)
    if (associated(cell_of_klake))  deallocate(cell_of_klake)
    n_lake = 0
    lake_built = .false.
  end subroutine elmxx_lake_clean

end module elmxxLakeMod
