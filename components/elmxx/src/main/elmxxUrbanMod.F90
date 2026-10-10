module elmxxUrbanMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! The urban landunit on the coupled path (plan Stage 6.7): the Fortran side
  ! of ELMxx's urban surface (`urban`, one record per urban landunit, with the
  ! five surfaces as members).
  !
  ! WHAT THIS OWNS (so far: U1, time-constant parameters)
  !   params   ELM's urban parameter set per landunit, from the surfdata
  !            density-type values (UrbanParamsType :: Init): geometry,
  !            emissivities, snow-free albedos, roof/wall/road thermal
  !            properties, building temperature limits; and what ELM derives
  !            from them -- the canyon view factors and the Macdonald (1998)
  !            displacement height and roughness length. Computed here,
  !            seeded once by NAME through the restart registry
  !            ("urban:<field>"), as elmxxLakeMod seeds the lake.
  !   check    elmxx_urban_diag_dump writes every seeded value, per packed
  !            landunit, to the ELMXX_DIAG trace ("elmxxurb:" labels), so
  !            tools/compare_urban_init.py can compare them with ELM's own
  !            step-1 records exactly.
  ! ELM is the reference for behaviour only; nothing is copied.
  !
  ! PACKED ORDER is elmxxKokkosStateMod's: urban landunits ascending in subgrid
  ! order, i.e. density-type-major (all TBD, then HD, then MD), every one that
  ! exists -- a 0%-weight density type is inactive, not absent.
  !
  ! ELM namelist settings assumed, with ELM's defaults (neither twin sets
  ! them): urban_hac = 'ON' (t_building_min/max from surfdata), urban_traffic
  ! = .false. (eflx_traffic_factor = 0), use_vancouver = use_mexicocity =
  ! .false.
  !-----------------------------------------------------------------------

  use shr_kind_mod    , only : r8 => shr_kind_r8
  use shr_sys_mod     , only : shr_sys_abort, shr_sys_flush
  use shr_const_mod   , only : SHR_CONST_KARMAN
  use elmxxSpmdMod    , only : masterproc, iam
  use elmxxSubgridMod , only : lun_gridcell, lun_itype, isturb_tbd, lun_wtgcell, &
                               num_landunits, num_columns, col_landunit, col_wtlunit, &
                               istsoil, col_itype, icol_roof, icol_sunwall, icol_shadewall, &
                               icol_road_imperv, icol_road_perv
  use elmxxSoilPropMod, only : nlevsoi, nlevgrnd, nlevsno, nlevtot, zsoi, dzsoi, zisoi, &
                               sp_watsat => watsat, sp_bsw => bsw, sp_sucsat => sucsat, &
                               sp_hksat => hksat, sp_watfc => watfc, sp_tkmg => tkmg, &
                               sp_tkdry => tkdry, sp_tksatu => tksatu, sp_csol => csol
  use elmxxSurfdataMod, only : urban_params_read, nlevurb, numrad_urb, &
                               canyon_hwr_in, ht_roof_in, wind_hgt_canyon_in, &
                               em_roof_in, em_wall_in, em_improad_in, em_perroad_in, &
                               thick_roof_in, thick_wall_in, &
                               t_building_min_in, t_building_max_in, nlev_improad_in, &
                               alb_roof_dir_in, alb_roof_dif_in, &
                               alb_wall_dir_in, alb_wall_dif_in, &
                               alb_improad_dir_in, alb_improad_dif_in, &
                               alb_perroad_dir_in, alb_perroad_dif_in, &
                               tk_roof_in, tk_wall_in, tk_improad_in, &
                               cv_roof_in, cv_wall_in, cv_improad_in, &
                               wtlunit_roof, wtroad_perv
  use elmxxKokkosStateMod, only : n_kokkos_urb, lun_of_kurb, n_kokkos_col, col_of_kcol
  use elmxxLakeMod    , only : lake_built
  use elmxxDiagnosticsMod, only : elmxx_diag_enabled, elmxx_diag_1d, &
                                  elmxx_diag_2d, elmxx_diag_int_1d
  use elmxx_mod       , only : ELMxxType, ELMXX_SUCCESS, &
                               ELMxxRestartFieldFind, ELMxxRestartFieldSet, &
                               ELMxxRestartFieldInfo, ELMxxRestartFieldGet, &
                               ELMxxHistorySetLandunitWeights, ELMxxHistorySetUrbanWeights, &
                               ELMxxAllocateUrbanSurface, ELMxxUrbanBuildFilters

  implicit none
  save
  private

  integer, parameter :: numrad = 2         ! ELMxx NUMRAD (VIS, NIR)
  integer, parameter :: nlevurb_x = 5      ! ELMxx NLEVURB

  logical, public :: urban_params_built = .false.
  logical, public :: urban_built = .false.       ! the urban surface is allocated and seeded
  integer, public :: n_urb_col = 0                ! packed urban columns = patches (5 per landunit)
  integer, parameter :: nlevurb_c = 5             ! ELM nlevurb
  real(r8), parameter :: spval = 1.e36_r8

  ! Per packed urban landunit (1:n_kokkos_urb), as seeded.
  real(r8), allocatable :: canyon_hwr(:), wtroad_prv(:), ht_roof(:), wtlun_roof(:)
  real(r8), allocatable :: wind_hgt_canyon(:), z_d_town(:), z_0_town(:)
  real(r8), allocatable :: vf_sr(:), vf_wr(:), vf_sw(:), vf_rw(:), vf_ww(:)
  real(r8), allocatable :: em_roof(:), em_wall(:), em_improad(:), em_perroad(:)
  real(r8), allocatable :: thick_roof(:), thick_wall(:), t_bld_min(:), t_bld_max(:)
  real(r8), allocatable :: eflx_traffic_factor(:)
  integer , allocatable :: nlev_improad(:)
  real(r8), allocatable :: alb_roof_dir(:,:), alb_roof_dif(:,:), alb_wall_dir(:,:), alb_wall_dif(:,:)
  real(r8), allocatable :: alb_improad_dir(:,:), alb_improad_dif(:,:)
  real(r8), allocatable :: alb_perroad_dir(:,:), alb_perroad_dif(:,:)
  real(r8), allocatable :: tk_roof(:,:), tk_wall(:,:), tk_improad(:,:)
  real(r8), allocatable :: cv_roof(:,:), cv_wall(:,:), cv_improad(:,:)

  public :: elmxx_urban_init
  public :: elmxx_urban_diag_dump
  public :: elmxx_urban_diag_step
  public :: elmxx_urban_clean

contains

  !-----------------------------------------------------------------------
  subroutine elmxx_urban_init(elm, logunit)
    !
    ! Compute ELM's time-constant urban parameters per packed landunit and
    ! seed them. A no-op on a domain without urban landunits.
    !
    implicit none
    type(ELMxxType), intent(in) :: elm
    integer, intent(in) :: logunit
    character(len=*), parameter :: subname = '(elmxx_urban_init) '

    if (n_kokkos_urb <= 0) return
    if (.not. urban_params_read) then
       call shr_sys_abort(subname//'ERROR: urban landunits exist but the surface '// &
            'dataset has no urban parameters (CANYON_HWR, ALB_*, TK_*, ...)')
    end if
    if (numrad_urb /= numrad .or. nlevurb /= nlevurb_x) then
       write(logunit,*) subname,'numrad ',numrad_urb,' nlevurb ',nlevurb
       call shr_sys_abort(subname//'ERROR: surfdata numrad/nlevurb differ from ELMxx NUMRAD/NLEVURB')
    end if

    call compute_params(logunit)
    call seed_params(elm)
    call push_history_weights(elm, logunit)
    urban_params_built = .true.
    call seed_surface(elm, logunit)
    urban_built = .true.

    write(logunit,*) subname,'rank ',iam,' seeded ',n_kokkos_urb,' urban landunits'
    call shr_sys_flush(logunit)

  end subroutine elmxx_urban_init

  !-----------------------------------------------------------------------
  subroutine compute_params(logunit)
    !
    ! ELM: UrbanParamsType :: Init, the urbpoi branch.
    !
    implicit none
    integer, intent(in) :: logunit
    integer  :: k, l, g, d, n
    real(r8) :: sumvf, plan_ai, frontal_ai, build_lw_ratio
    ! ELM: UrbanParamsType :: Init parameters
    real(r8), parameter :: alpha = 4.43_r8   ! coefficient for z_d_town
    real(r8), parameter :: beta  = 1.0_r8    ! coefficient for z_d_town
    real(r8), parameter :: C_d   = 1.2_r8    ! drag coefficient (Grimmond and Oke 1999)
    real(r8), parameter :: vkc   = SHR_CONST_KARMAN  ! ELM: elm_varcon vkc
    character(len=*), parameter :: subname = '(elmxx_urban_init) '

    n = n_kokkos_urb
    allocate(canyon_hwr(n), wtroad_prv(n), ht_roof(n), wtlun_roof(n), wind_hgt_canyon(n), &
             z_d_town(n), z_0_town(n), vf_sr(n), vf_wr(n), vf_sw(n), vf_rw(n), vf_ww(n), &
             em_roof(n), em_wall(n), em_improad(n), em_perroad(n), thick_roof(n), &
             thick_wall(n), t_bld_min(n), t_bld_max(n), eflx_traffic_factor(n), nlev_improad(n))
    allocate(alb_roof_dir(n,numrad), alb_roof_dif(n,numrad), alb_wall_dir(n,numrad), &
             alb_wall_dif(n,numrad), alb_improad_dir(n,numrad), alb_improad_dif(n,numrad), &
             alb_perroad_dir(n,numrad), alb_perroad_dif(n,numrad))
    allocate(tk_roof(n,nlevurb), tk_wall(n,nlevurb), tk_improad(n,nlevurb), &
             cv_roof(n,nlevurb), cv_wall(n,nlevurb), cv_improad(n,nlevurb))

    do k = 1, n
       l = lun_of_kurb(k)
       g = lun_gridcell(l)
       d = lun_itype(l) - isturb_tbd + 1    ! ELM: dindx = itype - isturb_MIN + 1

       canyon_hwr(k)      = canyon_hwr_in(g,d)
       wtroad_prv(k)      = wtroad_perv(g,d)
       ht_roof(k)         = ht_roof_in(g,d)
       wtlun_roof(k)      = wtlunit_roof(g,d)
       wind_hgt_canyon(k) = wind_hgt_canyon_in(g,d)
       em_roof(k)         = em_roof_in(g,d)
       em_wall(k)         = em_wall_in(g,d)
       em_improad(k)      = em_improad_in(g,d)
       em_perroad(k)      = em_perroad_in(g,d)
       thick_roof(k)      = thick_roof_in(g,d)
       thick_wall(k)      = thick_wall_in(g,d)
       nlev_improad(k)    = nint(nlev_improad_in(g,d))
       ! urban_hac = 'ON': the surfdata set points stand.
       t_bld_min(k)       = t_building_min_in(g,d)
       t_bld_max(k)       = t_building_max_in(g,d)
       ! urban_traffic = .false.
       eflx_traffic_factor(k) = 0.0_r8

       alb_roof_dir(k,:)    = alb_roof_dir_in(g,d,:)
       alb_roof_dif(k,:)    = alb_roof_dif_in(g,d,:)
       alb_wall_dir(k,:)    = alb_wall_dir_in(g,d,:)
       alb_wall_dif(k,:)    = alb_wall_dif_in(g,d,:)
       alb_improad_dir(k,:) = alb_improad_dir_in(g,d,:)
       alb_improad_dif(k,:) = alb_improad_dif_in(g,d,:)
       alb_perroad_dir(k,:) = alb_perroad_dir_in(g,d,:)
       alb_perroad_dif(k,:) = alb_perroad_dif_in(g,d,:)
       tk_roof(k,:)    = tk_roof_in(g,d,:)
       tk_wall(k,:)    = tk_wall_in(g,d,:)
       tk_improad(k,:) = tk_improad_in(g,d,:)
       cv_roof(k,:)    = cv_roof_in(g,d,:)
       cv_wall(k,:)    = cv_wall_in(g,d,:)
       cv_improad(k,:) = cv_improad_in(g,d,:)

       ! Canyon view factors (Masson 2000); depend on canyon_hwr only.
       vf_sr(k) = sqrt(canyon_hwr(k)**2 + 1._r8) - canyon_hwr(k)
       vf_wr(k) = 0.5_r8 * (1._r8 - vf_sr(k))
       vf_sw(k) = 0.5_r8 * (canyon_hwr(k) + 1._r8 - sqrt(canyon_hwr(k)**2+1._r8)) / canyon_hwr(k)
       vf_rw(k) = vf_sw(k)
       vf_ww(k) = 1._r8 - vf_sw(k) - vf_rw(k)
       sumvf = vf_sr(k) + 2._r8*vf_wr(k)
       if (abs(sumvf-1._r8) > 1.e-06_r8) call shr_sys_abort(subname//'ERROR: road view factors')
       sumvf = vf_sw(k) + vf_rw(k) + vf_ww(k)
       if (abs(sumvf-1._r8) > 1.e-06_r8) call shr_sys_abort(subname//'ERROR: wall view factors')

       ! Aerodynamic constants, Macdonald (1998) as in Grimmond and Oke (1999).
       plan_ai        = canyon_hwr(k)/(canyon_hwr(k) + 1._r8)
       build_lw_ratio = plan_ai
       frontal_ai     = (1._r8 - plan_ai) * canyon_hwr(k)
       frontal_ai     = frontal_ai * sqrt(1/build_lw_ratio) * sqrt(plan_ai)
       z_d_town(k)    = (1._r8 + alpha**(-plan_ai) * (plan_ai - 1._r8)) * ht_roof(k)
       z_0_town(k)    = ht_roof(k) * (1._r8 - z_d_town(k) / ht_roof(k)) * &
            exp(-1.0_r8 * (0.5_r8 * beta * C_d / vkc**2 * &
            (1 - z_d_town(k) / ht_roof(k)) * frontal_ai)**(-0.5_r8))
    end do

  end subroutine compute_params

  !-----------------------------------------------------------------------
  subroutine seed_params(elm)
    implicit none
    type(ELMxxType), intent(in) :: elm
    real(r8), allocatable :: spv(:,:)

    call urb_set_1d(elm, 'urban:canyonHWR'     , canyon_hwr)
    call urb_set_1d(elm, 'urban:wtRoadPerv'    , wtroad_prv)
    call urb_set_1d(elm, 'urban:htRoof'        , ht_roof)
    call urb_set_1d(elm, 'urban:wtRoof'        , wtlun_roof)
    call urb_set_1d(elm, 'urban:windHgtCanyon' , wind_hgt_canyon)
    call urb_set_1d(elm, 'urban:zDTown'        , z_d_town)
    call urb_set_1d(elm, 'urban:z0Town'        , z_0_town)
    call urb_set_1d(elm, 'urban:thickRoof'     , thick_roof)
    call urb_set_1d(elm, 'urban:thickWall'     , thick_wall)
    call urb_set_1d(elm, 'urban:nlevImproad'   , real(nlev_improad, r8))
    call urb_set_1d(elm, 'urban:tBuildingMin'  , t_bld_min)
    call urb_set_1d(elm, 'urban:tBuildingMax'  , t_bld_max)
    call urb_set_1d(elm, 'urban:eflxTrafficFactor', eflx_traffic_factor)
    call urb_set_1d(elm, 'urban:viewFactors.skyFromRoad'      , vf_sr)
    call urb_set_1d(elm, 'urban:viewFactors.wallFromRoad'     , vf_wr)
    call urb_set_1d(elm, 'urban:viewFactors.skyFromWall'      , vf_sw)
    call urb_set_1d(elm, 'urban:viewFactors.roadFromWall'     , vf_rw)
    call urb_set_1d(elm, 'urban:viewFactors.otherWallFromWall', vf_ww)

    ! Snow-free emissivity and albedo per surface; one wall value for both.
    call urb_set_1d(elm, 'urban:roof.emissivity'          , em_roof)
    call urb_set_1d(elm, 'urban:sunlitWall.emissivity'    , em_wall)
    call urb_set_1d(elm, 'urban:shadedWall.emissivity'    , em_wall)
    call urb_set_1d(elm, 'urban:imperviousRoad.emissivity', em_improad)
    call urb_set_1d(elm, 'urban:perviousRoad.emissivity'  , em_perroad)
    call urb_set_2d(elm, 'urban:roof.baseAlbedoDir'          , alb_roof_dir)
    call urb_set_2d(elm, 'urban:roof.baseAlbedoDif'          , alb_roof_dif)
    call urb_set_2d(elm, 'urban:sunlitWall.baseAlbedoDir'    , alb_wall_dir)
    call urb_set_2d(elm, 'urban:sunlitWall.baseAlbedoDif'    , alb_wall_dif)
    call urb_set_2d(elm, 'urban:shadedWall.baseAlbedoDir'    , alb_wall_dir)
    call urb_set_2d(elm, 'urban:shadedWall.baseAlbedoDif'    , alb_wall_dif)
    call urb_set_2d(elm, 'urban:imperviousRoad.baseAlbedoDir', alb_improad_dir)
    call urb_set_2d(elm, 'urban:imperviousRoad.baseAlbedoDif', alb_improad_dif)
    call urb_set_2d(elm, 'urban:perviousRoad.baseAlbedoDir'  , alb_perroad_dir)
    call urb_set_2d(elm, 'urban:perviousRoad.baseAlbedoDif'  , alb_perroad_dif)

    ! Layer thermal properties. The pervious road takes the soil's, so its
    ! slot holds spval: anything that reads it is reading the wrong thing.
    call urb_set_2d(elm, 'urban:roof.tkLayer'          , tk_roof)
    call urb_set_2d(elm, 'urban:roof.cvLayer'          , cv_roof)
    call urb_set_2d(elm, 'urban:sunlitWall.tkLayer'    , tk_wall)
    call urb_set_2d(elm, 'urban:sunlitWall.cvLayer'    , cv_wall)
    call urb_set_2d(elm, 'urban:shadedWall.tkLayer'    , tk_wall)
    call urb_set_2d(elm, 'urban:shadedWall.cvLayer'    , cv_wall)
    call urb_set_2d(elm, 'urban:imperviousRoad.tkLayer', tk_improad)
    call urb_set_2d(elm, 'urban:imperviousRoad.cvLayer', cv_improad)
    allocate(spv(n_kokkos_urb, nlevurb)); spv = 1.e36_r8
    call urb_set_2d(elm, 'urban:perviousRoad.tkLayer'  , spv)
    call urb_set_2d(elm, 'urban:perviousRoad.cvLayer'  , spv)

  end subroutine seed_params

  !-----------------------------------------------------------------------
  subroutine push_history_weights(elm, logunit)
    !
    ! History's urban half of the cell mean (U6): per packed natural column
    ! (one per cell) the cell's packed urban landunits; per landunit its
    ! weight on the cell; per urban column (5k+s, s in ELM's column order)
    ! its weight on the landunit and ELM's urbanf/urbans c2l scale factors
    ! (subgridAveMod create_scale_c2l). Without lake, the natural weights are
    ! pushed here too (elmxxLakeMod pushes them when lake exists).
    !
    implicit none
    type(ELMxxType), intent(in) :: elm
    integer, intent(in) :: logunit
    integer :: kc, k, l, c, s, n, nidx, ierr
    integer , allocatable :: ptr(:), idx(:), cell_of_kc(:), lake_of(:)
    real(r8), allocatable :: wlun(:), cwt(:), csf(:), css(:), wn(:), wl(:)
    real(r8) :: hwr
    character(len=*), parameter :: subname = '(elmxx_urban_init) '

    n = n_kokkos_urb
    allocate(ptr(n_kokkos_col+1), idx(max(n,1)), cell_of_kc(n_kokkos_col), &
             wlun(n), cwt(5*n), csf(5*n), css(5*n))
    do kc = 1, n_kokkos_col
       cell_of_kc(kc) = lun_gridcell(col_landunit(col_of_kcol(kc)))
    end do

    ! CSR: urban landunits per natural column, in packed order.
    nidx = 0
    ptr(1) = 0
    do kc = 1, n_kokkos_col
       do k = 1, n
          if (lun_gridcell(lun_of_kurb(k)) == cell_of_kc(kc)) then
             nidx = nidx + 1
             idx(nidx) = k - 1
          end if
       end do
       ptr(kc+1) = nidx
    end do
    if (nidx /= n) call shr_sys_abort(subname//'ERROR: an urban landunit has no natural column in its cell')

    do k = 1, n
       l = lun_of_kurb(k)
       wlun(k) = lun_wtgcell(l)
       hwr = canyon_hwr(k)
       s = 0
       do c = 1, num_columns
          if (col_landunit(c) /= l) cycle
          s = s + 1
          if (s > 5) call shr_sys_abort(subname//'ERROR: an urban landunit has more than five columns')
          cwt(5*(k-1)+s) = col_wtlunit(c)
          select case (s)
          case (1)         ! roof
             csf(5*(k-1)+s) = 1._r8
             css(5*(k-1)+s) = 1._r8
          case (2, 3)      ! sunlit, shaded wall
             csf(5*(k-1)+s) = 3.0_r8 * hwr
             css(5*(k-1)+s) = (3.0_r8 * hwr) / (2._r8*hwr + 1._r8)
          case (4, 5)      ! impervious, pervious road
             csf(5*(k-1)+s) = 3.0_r8
             css(5*(k-1)+s) = 3.0_r8 / (2._r8*hwr + 1._r8)
          end select
       end do
       if (s /= 5) call shr_sys_abort(subname//'ERROR: an urban landunit does not have five columns')
    end do

    if (.not. lake_built) then
       allocate(lake_of(n_kokkos_col), wn(n_kokkos_col), wl(n_kokkos_col))
       lake_of = -1; wn = 0._r8; wl = 0._r8
       do kc = 1, n_kokkos_col
          do l = 1, num_landunits
             if (lun_gridcell(l) == cell_of_kc(kc) .and. lun_itype(l) == istsoil) wn(kc) = lun_wtgcell(l)
          end do
       end do
       call ELMxxHistorySetLandunitWeights(elm, lake_of, wn, wl, n_kokkos_col, ierr)
       if (ierr /= ELMXX_SUCCESS) call shr_sys_abort(subname//'ERROR: ELMxxHistorySetLandunitWeights failed')
       deallocate(lake_of, wn, wl)
    end if

    call ELMxxHistorySetUrbanWeights(elm, ptr, idx, nidx, wlun, cwt, csf, css, n, &
         n_kokkos_col, ierr)
    if (ierr /= ELMXX_SUCCESS) call shr_sys_abort(subname//'ERROR: ELMxxHistorySetUrbanWeights failed')
    if (masterproc) then
       write(logunit,*) subname,'history urban landunit weights ',wlun
       call shr_sys_flush(logunit)
    end if
    deallocate(ptr, idx, cell_of_kc, wlun, cwt, csf, css)

  end subroutine push_history_weights


  !-----------------------------------------------------------------------
  subroutine seed_surface(elm, logunit)
    !
    ! The urban surface's columns and patches (U1/U2): allocate urbanCol/
    ! urbanPatch, five columns per packed landunit k (packed column
    ! 5(k-1)+s, s in ELM's column order, which is the subgrid's), one patch
    ! each; seed their topology and ELM's urban geometry, road soil
    ! properties and cold start by NAME through the registry; then let the
    ! library build its active-only filters.
    !
    !   geometry    ELM initVerticalMod, urban branch (use_vancouver and
    !               use_mexicocity off): roof and walls on nlevurb = 5 layers
    !               from THICK_ROOF/THICK_WALL, deeper slots spval; roads on
    !               the soil grid. nlevbed = nlevsoi for every column
    !               (use_var_soil_thick = .false., as ELM falls back to).
    !   properties  ELM SoilStateType: roads take the natural pedotransfer
    !               with no organic matter (elmxxSoilPropMod zeroes om_frac
    !               on urban columns); roof and walls spval. watdry/watopt
    !               from the same; the pervious road's root fraction is ELM's
    !               uniform 0.1 over nlevsoi.
    !   cold start  ELM ColumnDataType/LandunitDataType/SoilHydrologyType
    !               InitCold, urban branches: roads 274 K over nlevgrnd, roof
    !               and walls 292 K over nlevurb; pervious road 0.3 v/v above
    !               bedrock (capped at porosity), impervious road and roof/
    !               walls dry; t_grnd from layer 1; emg the snow-free surface
    !               emissivity; taf 283 K, qaf 1e-4; pervious road wa 4800,
    !               zwt = zi(nlevsoi) + 1; zwt_perched/frost_table spval.
    !
    implicit none
    type(ELMxxType), intent(in) :: elm
    integer, intent(in) :: logunit
    integer  :: k, l, c, s, kc, j, nc, ierr
    integer , allocatable :: ib(:), ib2(:)
    real(r8), allocatable :: r1(:), dz(:,:), z(:,:), zi(:,:)
    real(r8), allocatable :: rg(:,:), rg2(:,:), rt(:,:), rt1(:,:), rp1(:,:), rp2(:,:)
    real(r8), allocatable :: tsoi(:,:), liq(:,:), ice(:,:), vol(:,:)
    real(r8) :: zu(nlevurb_c), dzu(nlevurb_c), ziu(0:nlevurb_c), thick
    integer  :: subcol(5)
    real(r8), parameter :: denh2o = 1000._r8, tkfrz = 273.15_r8
    character(len=*), parameter :: subname = '(elmxx_urban_init) '

    nc = 5 * n_kokkos_urb
    n_urb_col = nc
    call ELMxxAllocateUrbanSurface(elm, nc, nc, ierr)
    if (ierr /= ELMXX_SUCCESS) call shr_sys_abort(subname//'ERROR: ELMxxAllocateUrbanSurface failed')

    allocate(ib(n_kokkos_urb), r1(n_kokkos_urb))
    ! ---- landunit topology and the explicit column maps ----
    do k = 1, n_kokkos_urb
       l = lun_of_kurb(k)
       ib(k) = merge(1, 0, lun_wtgcell(l) > 0._r8)   ! ELM: active iff weight > 0
    end do
    call urb_set_1d(elm, 'urban:active', real(ib, r8))
    call urb_set_1d(elm, 'urban:urbpoi', spread(1._r8, 1, n_kokkos_urb))
    do k = 1, n_kokkos_urb
       ib(k) = natcol_of_cell(lun_gridcell(lun_of_kurb(k)))
    end do
    call urb_set_1d(elm, 'urban:natcol', real(ib, r8))
    do k = 1, n_kokkos_urb
       ib(k) = lun_gridcell(lun_of_kurb(k)) - 1
    end do
    call urb_set_1d(elm, 'urban:gridcell', real(ib, r8))
    do s = 1, 5
       do k = 1, n_kokkos_urb
          ib(k) = 5*(k-1) + s - 1
       end do
       select case (s)
       case (1); call urb_set_1d(elm, 'urban:colRoof'     , real(ib, r8))
       case (2); call urb_set_1d(elm, 'urban:colSunwall'  , real(ib, r8))
       case (3); call urb_set_1d(elm, 'urban:colShadewall', real(ib, r8))
       case (4); call urb_set_1d(elm, 'urban:colImproad'  , real(ib, r8))
       case (5); call urb_set_1d(elm, 'urban:colPerroad'  , real(ib, r8))
       end select
    end do
    call urb_set_1d(elm, 'urban:taf', spread(283._r8, 1, n_kokkos_urb))
    call urb_set_1d(elm, 'urban:qaf', spread(1.e-4_r8, 1, n_kokkos_urb))
    call urb_set_1d(elm, 'urban:tBuilding', spread(spval, 1, n_kokkos_urb))
    deallocate(ib, r1)

    ! ---- column topology ----
    allocate(ib(nc), ib2(nc), r1(nc))
    do k = 1, n_kokkos_urb
       call landunit_columns(lun_of_kurb(k), subcol)
       do s = 1, 5
          kc = 5*(k-1) + s
          c = subcol(s)
          ib(kc) = col_itype(c)
          if (ib(kc) /= icol_roof + s - 1) call shr_sys_abort(subname// &
               'ERROR: urban columns are not in ELM''s column order')
       end do
    end do
    call urb_set_1d(elm, 'urbcol:col_itype', real(ib, r8))
    do kc = 1, nc
       k = (kc - 1)/5 + 1
       ib2(kc) = lun_itype(lun_of_kurb(k))
    end do
    call urb_set_1d(elm, 'urbcol:lun_itype', real(ib2, r8))
    call urb_set_1d(elm, 'urbcol:lun_urbpoi', spread(1._r8, 1, nc))
    call urb_set_1d(elm, 'urbcol:col_is_soil', spread(0._r8, 1, nc))
    call urb_set_1d(elm, 'urbcol:col_is_crop', spread(0._r8, 1, nc))
    do kc = 1, nc
       k = (kc - 1)/5 + 1
       ib2(kc) = merge(1, 0, lun_wtgcell(lun_of_kurb(k)) > 0._r8)
    end do
    call urb_set_1d(elm, 'urbcol:col_active', real(ib2, r8))
    call urb_set_1d(elm, 'urbpatch:patch_active', real(ib2, r8))
    do kc = 1, nc
       k = (kc - 1)/5 + 1
       ib2(kc) = lun_gridcell(lun_of_kurb(k)) - 1
    end do
    call urb_set_1d(elm, 'urbcol:col_gridcell', real(ib2, r8))
    call urb_set_1d(elm, 'urbcol:col_topounit', real(ib2, r8))
    do kc = 1, nc
       ib2(kc) = (kc - 1)/5
    end do
    call urb_set_1d(elm, 'urbcol:col_landunit', real(ib2, r8))
    call urb_set_1d(elm, 'urbpatch:patch_landunit', real(ib2, r8))
    do kc = 1, nc
       ib2(kc) = kc - 1
    end do
    call urb_set_1d(elm, 'urbcol:col_pfti', real(ib2, r8))
    call urb_set_1d(elm, 'urbpatch:patch_column', real(ib2, r8))
    call urb_set_1d(elm, 'urbcol:col_npfts', spread(1._r8, 1, nc))
    call urb_set_1d(elm, 'urbcol:nlevbed', spread(real(nlevsoi, r8), 1, nc))
    call urb_set_1d(elm, 'urbpatch:wtcol', spread(1._r8, 1, nc))
    call urb_set_1d(elm, 'urbpatch:frac_veg_nosno', spread(0._r8, 1, nc))
    call urb_set_1d(elm, 'urbpatch:is_on_soil_col', spread(0._r8, 1, nc))
    call urb_set_1d(elm, 'urbpatch:is_on_crop_col', spread(0._r8, 1, nc))

    ! ---- geometry ----
    allocate(dz(nc,nlevgrnd), z(nc,nlevgrnd), zi(nc,0:nlevgrnd))
    do k = 1, n_kokkos_urb
       do s = 1, 5
          kc = 5*(k-1) + s
          if (s <= 3) then
             ! ELM initVerticalMod: node depths evenly spaced over the
             ! thickness, layer thicknesses and interfaces from them.
             if (s == 1) then; thick = thick_roof(k); else; thick = thick_wall(k); end if
             do j = 1, nlevurb_c
                zu(j) = (j - 0.5_r8)*(thick/real(nlevurb_c, r8))
             end do
             dzu(1) = 0.5_r8*(zu(1) + zu(2))
             do j = 2, nlevurb_c - 1
                dzu(j) = 0.5_r8*(zu(j+1) - zu(j-1))
             end do
             dzu(nlevurb_c) = zu(nlevurb_c) - zu(nlevurb_c-1)
             ziu(0) = 0._r8
             do j = 1, nlevurb_c - 1
                ziu(j) = 0.5_r8*(zu(j) + zu(j+1))
             end do
             ziu(nlevurb_c) = zu(nlevurb_c) + 0.5_r8*dzu(nlevurb_c)
             dz(kc,:) = spval; z(kc,:) = spval; zi(kc,:) = spval
             dz(kc,1:nlevurb_c) = dzu
             z (kc,1:nlevurb_c) = zu
             zi(kc,0:nlevurb_c) = ziu
          else
             dz(kc,:) = dzsoi
             z (kc,:) = zsoi
             zi(kc,:) = zisoi
          end if
       end do
    end do
    call urb_set_2d(elm, 'urbcol:dz_soi', dz)
    call urb_set_2d(elm, 'urbcol:zc_soi', z)
    call urb_set_2d(elm, 'urbcol:zi_soi', zi)
    ! The combined (ELM-ordered) and SoilTemperature (_p1) layouts: snow
    ! slots empty at a cold start; _p1 carries the standing-water node.
    allocate(rt(nc,nlevtot), rt1(nc,nlevtot+1))
    rt = 0._r8; rt(:, nlevsno+1:nlevtot) = dz
    call urb_set_2d(elm, 'urbcol:dz', rt)
    rt = 0._r8; rt(:, nlevsno+1:nlevtot) = z
    call urb_set_2d(elm, 'urbcol:z', rt)
    rt1 = 0._r8; rt1(:, nlevsno+1:nlevtot+1) = zi
    call urb_set_2d(elm, 'urbcol:zi', rt1)
    allocate(rp1(nc,nlevsno+1+nlevgrnd), rp2(nc,nlevsno+2+nlevgrnd))
    rp1 = 0._r8; rp1(:, nlevsno+2:) = dz
    call urb_set_2d(elm, 'urbcol:dz_p1', rp1)
    rp1 = 0._r8; rp1(:, nlevsno+2:) = z
    call urb_set_2d(elm, 'urbcol:z_p1', rp1)
    rp2 = 0._r8; rp2(:, nlevsno+2:) = zi
    call urb_set_2d(elm, 'urbcol:zi_p1', rp2)

    ! ---- soil properties: roads natural (no organic), roof/walls spval ----
    allocate(rg(nc,nlevgrnd), rg2(nc,nlevgrnd))
    call seed_prop('urbcol:watsat', sp_watsat)
    call seed_prop('urbcol:watsat_soi', sp_watsat)
    call seed_prop('urbcol:bsw', sp_bsw)
    call seed_prop('urbcol:sucsat', sp_sucsat)
    call seed_prop('urbcol:hksat', sp_hksat)
    call seed_prop('urbcol:watfc', sp_watfc)
    call seed_prop('urbcol:tkmg', sp_tkmg)
    call seed_prop('urbcol:tkdry', sp_tkdry)
    call seed_prop('urbcol:tksatu', sp_tksatu)
    call seed_prop('urbcol:csol', sp_csol)
    ! watdry / watopt (ELM SoilStateType), roads; the impervious road's
    ! are spval in ELM, harmless since only the pervious road reads them.
    rg = spval; rg2 = spval
    do k = 1, n_kokkos_urb
       call landunit_columns(lun_of_kurb(k), subcol)
       kc = 5*(k-1) + 5
       c = subcol(5)
       do j = 1, nlevgrnd
          rg(kc,j)  = sp_watsat(c,j) * (316230._r8/sp_sucsat(c,j)) ** (-1._r8/sp_bsw(c,j))
          rg2(kc,j) = sp_watsat(c,j) * (158490._r8/sp_sucsat(c,j)) ** (-1._r8/sp_bsw(c,j))
       end do
    end do
    call urb_set_2d(elm, 'urbcol:watdry', rg)
    call urb_set_2d(elm, 'urbcol:watopt', rg2)
    rg = 0._r8
    do k = 1, n_kokkos_urb
       rg(5*(k-1)+5, 1:nlevsoi) = 0.1_r8
    end do
    call urb_set_2d(elm, 'urbcol:rootfr_road_perv', rg)

    ! ---- cold start ----
    allocate(tsoi(nc,nlevgrnd), liq(nc,nlevgrnd), ice(nc,nlevgrnd), vol(nc,nlevgrnd))
    tsoi = spval; liq = spval; ice = spval; vol = spval
    do k = 1, n_kokkos_urb
       call landunit_columns(lun_of_kurb(k), subcol)
       do s = 1, 5
          kc = 5*(k-1) + s
          c = subcol(s)
          if (s <= 3) then
             tsoi(kc,1:nlevurb_c) = 292._r8
             vol (kc,1:nlevurb_c) = 0._r8
          else
             tsoi(kc,:) = 274._r8
             vol(kc,:) = 0._r8
             if (s == 5) then
                do j = 1, nlevsoi                       ! above bedrock
                   vol(kc,j) = 0.3_r8
                end do
             end if
             do j = 1, nlevgrnd
                vol(kc,j) = min(vol(kc,j), sp_watsat(c,j))
             end do
          end if
          do j = 1, merge(nlevurb_c, nlevgrnd, s <= 3)
             if (tsoi(kc,j) <= tkfrz) then
                ice(kc,j) = dz(kc,j)*917._r8*vol(kc,j); liq(kc,j) = 0._r8
             else
                liq(kc,j) = dz(kc,j)*denh2o*vol(kc,j); ice(kc,j) = 0._r8
             end if
          end do
       end do
    end do
    call urb_set_2d(elm, 'urbcol:t_soisno_soi', tsoi)
    call urb_set_2d(elm, 'urbcol:h2osoi_liq_soi', liq)
    call urb_set_2d(elm, 'urbcol:h2osoi_ice_soi', ice)
    call urb_set_2d(elm, 'urbcol:h2osoi_vol', vol)
    rt = 0._r8; rt(:, nlevsno+1:nlevtot) = tsoi
    call urb_set_2d(elm, 'urbcol:t_soisno', rt)
    rt = 0._r8; rt(:, nlevsno+1:nlevtot) = liq
    call urb_set_2d(elm, 'urbcol:h2osoi_liq', rt)
    rt = 0._r8; rt(:, nlevsno+1:nlevtot) = ice
    call urb_set_2d(elm, 'urbcol:h2osoi_ice', rt)
    r1 = tsoi(:,1)
    call urb_set_1d(elm, 'urbcol:t_grnd', r1)
    call urb_set_1d(elm, 'urbcol:t_h2osfc', spread(274._r8, 1, nc))
    do k = 1, n_kokkos_urb
       r1(5*(k-1)+1) = em_roof(k)
       r1(5*(k-1)+2) = em_wall(k)
       r1(5*(k-1)+3) = em_wall(k)
       r1(5*(k-1)+4) = em_improad(k)
       r1(5*(k-1)+5) = em_perroad(k)
    end do
    call urb_set_1d(elm, 'urbcol:emg', r1)
    r1 = spval
    do k = 1, n_kokkos_urb
       r1(5*(k-1)+5) = 4800._r8
    end do
    call urb_set_1d(elm, 'urbcol:wa', r1)
    r1 = spval
    do k = 1, n_kokkos_urb
       r1(5*(k-1)+5) = (25._r8 + zisoi(nlevsoi)) - 4800._r8/0.2_r8/1000._r8
    end do
    call urb_set_1d(elm, 'urbcol:zwt', r1)
    call urb_set_1d(elm, 'urbcol:zwt_perched', spread(spval, 1, nc))
    call urb_set_1d(elm, 'urbcol:frost_table', spread(spval, 1, nc))
    call urb_set_1d(elm, 'urbcol:h2osfc_thresh', spread(0._r8, 1, nc))

    call ELMxxUrbanBuildFilters(elm, ierr)
    if (ierr /= ELMXX_SUCCESS) call shr_sys_abort(subname//'ERROR: ELMxxUrbanBuildFilters failed')

    if (masterproc) then
       write(logunit,*) subname,'urban surface: ',nc,' columns; roof dz ',dz(1,1:nlevurb_c)
       write(logunit,*) subname,'urban cold start: roads 274 K, roof/walls 292 K, ', &
            'pervious-road h2osoi_vol(1) ',vol(5,1)
       call shr_sys_flush(logunit)
    end if
    deallocate(ib, ib2, r1, dz, z, zi, rt, rt1, rp1, rp2, rg, rg2, tsoi, liq, ice, vol)

  contains

    subroutine seed_prop(name, src)
      character(len=*), intent(in) :: name
      real(r8), intent(in) :: src(:,:)
      integer :: kk, ss, cc
      rg = spval
      do kk = 1, n_kokkos_urb
         call landunit_columns(lun_of_kurb(kk), subcol)
         do ss = 4, 5                                 ! roads
            cc = subcol(ss)
            rg(5*(kk-1)+ss, :) = src(cc, 1:nlevgrnd)
         end do
      end do
      call urb_set_2d(elm, name, rg)
    end subroutine seed_prop

  end subroutine seed_surface

  !-----------------------------------------------------------------------
  subroutine landunit_columns(l, cols)
    ! The five subgrid columns of urban landunit l, in subgrid (ELM) order.
    implicit none
    integer, intent(in)  :: l
    integer, intent(out) :: cols(5)
    integer :: c, s
    s = 0
    do c = 1, num_columns
       if (col_landunit(c) /= l) cycle
       s = s + 1
       if (s > 5) call shr_sys_abort('(elmxx_urban_init) ERROR: urban landunit with more than five columns')
       cols(s) = c
    end do
    if (s /= 5) call shr_sys_abort('(elmxx_urban_init) ERROR: urban landunit without five columns')
  end subroutine landunit_columns

  !-----------------------------------------------------------------------
  integer function natcol_of_cell(g)
    ! The packed (0-based) natural column of local cell g: every cell has one.
    implicit none
    integer, intent(in) :: g
    integer :: kc
    natcol_of_cell = -1
    do kc = 1, n_kokkos_col
       if (lun_gridcell(col_landunit(col_of_kcol(kc))) == g) then
          natcol_of_cell = kc - 1
          return
       end if
    end do
    call shr_sys_abort('(elmxx_urban_init) ERROR: an urban cell has no natural column')
  end function natcol_of_cell

  !-----------------------------------------------------------------------
  subroutine elmxx_urban_diag_dump()
    !
    ! Every seeded parameter, packed order, plus the packed -> subgrid
    ! landunit map. Labels name ELM's variable, so a comparison needs no
    ! translation table (tools/compare_urban_init.py).
    !
    implicit none
    integer :: n
    if (.not. elmxx_diag_enabled .or. .not. urban_params_built) return
    n = n_kokkos_urb
    call elmxx_diag_int_1d('elmxxurb:lun_of_kurb', lun_of_kurb(1:n), n)
    call elmxx_diag_1d('elmxxurb:canyon_hwr'     , canyon_hwr, n)
    call elmxx_diag_1d('elmxxurb:wtroad_perv'    , wtroad_prv, n)
    call elmxx_diag_1d('elmxxurb:ht_roof'        , ht_roof, n)
    call elmxx_diag_1d('elmxxurb:wtlunit_roof'   , wtlun_roof, n)
    call elmxx_diag_1d('elmxxurb:wind_hgt_canyon', wind_hgt_canyon, n)
    call elmxx_diag_1d('elmxxurb:z_d_town'       , z_d_town, n)
    call elmxx_diag_1d('elmxxurb:z_0_town'       , z_0_town, n)
    call elmxx_diag_1d('elmxxurb:eflx_traffic_factor', eflx_traffic_factor, n)
    call elmxx_diag_1d('elmxxurb:vf_sr', vf_sr, n)
    call elmxx_diag_1d('elmxxurb:vf_wr', vf_wr, n)
    call elmxx_diag_1d('elmxxurb:vf_sw', vf_sw, n)
    call elmxx_diag_1d('elmxxurb:vf_rw', vf_rw, n)
    call elmxx_diag_1d('elmxxurb:vf_ww', vf_ww, n)
    call elmxx_diag_1d('elmxxurb:em_roof'   , em_roof, n)
    call elmxx_diag_1d('elmxxurb:em_wall'   , em_wall, n)
    call elmxx_diag_1d('elmxxurb:em_improad', em_improad, n)
    call elmxx_diag_1d('elmxxurb:em_perroad', em_perroad, n)
    call elmxx_diag_2d('elmxxurb:alb_roof_dir'   , alb_roof_dir, n, numrad)
    call elmxx_diag_2d('elmxxurb:alb_roof_dif'   , alb_roof_dif, n, numrad)
    call elmxx_diag_2d('elmxxurb:alb_wall_dir'   , alb_wall_dir, n, numrad)
    call elmxx_diag_2d('elmxxurb:alb_wall_dif'   , alb_wall_dif, n, numrad)
    call elmxx_diag_2d('elmxxurb:alb_improad_dir', alb_improad_dir, n, numrad)
    call elmxx_diag_2d('elmxxurb:alb_improad_dif', alb_improad_dif, n, numrad)
    call elmxx_diag_2d('elmxxurb:alb_perroad_dir', alb_perroad_dir, n, numrad)
    call elmxx_diag_2d('elmxxurb:alb_perroad_dif', alb_perroad_dif, n, numrad)
  end subroutine elmxx_urban_diag_dump


  !-----------------------------------------------------------------------
  subroutine elmxx_urban_diag_step(elm, tag)
    !
    ! Per-step urban column/patch state into the ELMXX_DIAG trace, packed
    ! order (urban column c is ELM column c+1 on a one-cell domain whose
    ! natural column comes first). tag: 'elmxx_urbin' at the top of the step
    ! (pairs with ELM canhydro_in), 'elmxx_urbflx' after phase 1 (pairs with
    ! urbanflux_out), 'elmxx_urbout' at the end of the step.
    !
    implicit none
    type(ELMxxType), intent(in) :: elm
    character(len=*), intent(in) :: tag
    real(r8), allocatable :: v(:), v2(:,:)
    integer :: j
    character(len=16), parameter :: cols(6) = (/ 't_grnd          ', 'h2osno          ', &
         'frac_sno        ', 'snow_depth      ', 't_soisno_soi    ', 'h2osoi_liq_soi  ' /)
    character(len=16), parameter :: pchs(7) = (/ 'eflx_sh_grnd    ', 'qflx_evap_soi   ', &
         'sabg            ', 'eflx_lwrad_net  ', 'eflx_soil_grnd  ', 't_ref2m         ', &
         'cgrnd           ' /)
    if (.not. elmxx_diag_enabled .or. .not. urban_built) return
    allocate(v(n_urb_col), v2(n_urb_col, nlevgrnd))
    do j = 1, 4
       call urb_get_1d(elm, 'urbcol:'//trim(cols(j)), v)
       call elmxx_diag_1d(trim(tag)//':'//trim(cols(j)), v, n_urb_col)
    end do
    call urb_get_1d(elm, 'urbcol:snl', v)
    call elmxx_diag_1d(trim(tag)//':snl', v, n_urb_col)
    do j = 5, 6
       call urb_get_2d(elm, 'urbcol:'//trim(cols(j)), v2)
       call elmxx_diag_2d(trim(tag)//':'//trim(cols(j)), v2, n_urb_col, nlevgrnd)
    end do
    do j = 1, 7
       call urb_get_1d(elm, 'urbpatch:'//trim(pchs(j)), v)
       call elmxx_diag_1d(trim(tag)//':'//trim(pchs(j)), v, n_urb_col)
    end do
    call urb_get_1d(elm, 'urbpatch:cgrnds', v); call elmxx_diag_1d(trim(tag)//':cgrnds', v, n_urb_col)
    call urb_get_1d(elm, 'urbpatch:cgrndl', v); call elmxx_diag_1d(trim(tag)//':cgrndl', v, n_urb_col)
    call urb_get_1d(elm, 'urbcol:htvp', v);     call elmxx_diag_1d(trim(tag)//':htvp', v, n_urb_col)
    call urb_get_1d(elm, 'urbcol:qg', v);       call elmxx_diag_1d(trim(tag)//':qg', v, n_urb_col)
    call urb_get_1d(elm, 'urbcol:dqgdT', v);    call elmxx_diag_1d(trim(tag)//':dqgdT', v, n_urb_col)
    call urb_get_1d(elm, 'urbcol:qflx_sub_snow', v); call elmxx_diag_1d(trim(tag)//':qflx_sub_snow_col', v, n_urb_col)
    call urb_get_1d(elm, 'urbpatch:qflx_sub_snow', v); call elmxx_diag_1d(trim(tag)//':qflx_sub_snow', v, n_urb_col)
    call urb_get_1d(elm, 'urbcol:qflx_evap_grnd_col', v); call elmxx_diag_1d(trim(tag)//':qflx_evap_grnd', v, n_urb_col)
    call urb_get_2d(elm, 'urbcol:h2osoi_ice_soi', v2)
    call elmxx_diag_2d(trim(tag)//':h2osoi_ice_soi', v2, n_urb_col, nlevgrnd)
    call urb_get_1d(elm, 'urbcol:hs_top_snow', v); call elmxx_diag_1d(trim(tag)//':hs_top_snow', v, n_urb_col)
    call urb_get_1d(elm, 'urbcol:dhsdT', v);       call elmxx_diag_1d(trim(tag)//':dhsdT', v, n_urb_col)
    deallocate(v2); allocate(v2(n_urb_col, nlevsno))
    call urb_get_2d(elm, 'urbcol:t_soisno_sno', v2)
    call elmxx_diag_2d(trim(tag)//':t_soisno_sno', v2, n_urb_col, nlevsno)
    call urb_get_2d(elm, 'urbcol:dz_sno', v2)
    call elmxx_diag_2d(trim(tag)//':dz_sno', v2, n_urb_col, nlevsno)
    call urb_get_2d(elm, 'urbcol:h2osoi_ice_sno', v2)
    call elmxx_diag_2d(trim(tag)//':h2osoi_ice_sno', v2, n_urb_col, nlevsno)
    deallocate(v, v2)
  end subroutine elmxx_urban_diag_step

  subroutine urb_get_1d(elm, name, v)
    implicit none
    type(ELMxxType), intent(in) :: elm
    character(len=*), intent(in) :: name
    real(r8), intent(inout) :: v(:)
    integer :: i, ierr
    i = field_index(elm, name, size(v), 0)
    call ELMxxRestartFieldGet(elm, i, v, size(v), ierr)
    if (ierr /= ELMXX_SUCCESS) call shr_sys_abort('(elmxx_urban) ERROR: get '//trim(name))
  end subroutine urb_get_1d

  subroutine urb_get_2d(elm, name, v)
    implicit none
    type(ELMxxType), intent(in) :: elm
    character(len=*), intent(in) :: name
    real(r8), intent(inout) :: v(:,:)
    real(r8), allocatable :: flat(:)
    integer :: i, ierr
    i = field_index(elm, name, size(v,1), size(v,2))
    allocate(flat(size(v)))
    call ELMxxRestartFieldGet(elm, i, flat, size(flat), ierr)
    if (ierr /= ELMXX_SUCCESS) call shr_sys_abort('(elmxx_urban) ERROR: get '//trim(name))
    v = transpose(reshape(flat, (/ size(v,2), size(v,1) /)))
  end subroutine urb_get_2d

  !-----------------------------------------------------------------------
  ! Name-addressed seeding through the restart registry (elmxxLakeMod's
  ! pattern): a flat ROW-major buffer, buf((i-1)*n2 + j) = v(i,j).
  !-----------------------------------------------------------------------
  integer function field_index(elm, name, n1, n2)
    implicit none
    type(ELMxxType), intent(in) :: elm
    character(len=*), intent(in) :: name
    integer, intent(in) :: n1, n2
    integer :: i, m1, m2, kind, ierr
    character(len=128) :: got
    call ELMxxRestartFieldFind(elm, name, i, ierr)
    if (ierr /= ELMXX_SUCCESS) call shr_sys_abort('(elmxx_urban) ERROR: no registry field '//trim(name))
    call ELMxxRestartFieldInfo(elm, i, got, m1, m2, kind, ierr)
    if (m1 /= n1 .or. m2 /= n2) then
       write(*,*) '(elmxx_urban) ',trim(name),' registry extents ',m1,m2,' given ',n1,n2
       call shr_sys_abort('(elmxx_urban) ERROR: extent mismatch on '//trim(name))
    end if
    field_index = i
  end function field_index

  subroutine urb_set_1d(elm, name, v)
    implicit none
    type(ELMxxType), intent(in) :: elm
    character(len=*), intent(in) :: name
    real(r8), intent(in) :: v(:)
    integer :: i, ierr
    i = field_index(elm, name, size(v), 0)
    call ELMxxRestartFieldSet(elm, i, v, size(v), ierr)
    if (ierr /= ELMXX_SUCCESS) call shr_sys_abort('(elmxx_urban) ERROR: set '//trim(name))
  end subroutine urb_set_1d

  subroutine urb_set_2d(elm, name, v)
    implicit none
    type(ELMxxType), intent(in) :: elm
    character(len=*), intent(in) :: name
    real(r8), intent(in) :: v(:,:)
    real(r8), allocatable :: flat(:)
    integer :: i, ierr
    i = field_index(elm, name, size(v,1), size(v,2))
    flat = reshape(transpose(v), (/ size(v) /))
    call ELMxxRestartFieldSet(elm, i, flat, size(flat), ierr)
    if (ierr /= ELMXX_SUCCESS) call shr_sys_abort('(elmxx_urban) ERROR: set '//trim(name))
  end subroutine urb_set_2d

  !-----------------------------------------------------------------------
  subroutine elmxx_urban_clean()
    implicit none
    if (allocated(canyon_hwr)) then
       deallocate(canyon_hwr, wtroad_prv, ht_roof, wtlun_roof, wind_hgt_canyon, &
            z_d_town, z_0_town, vf_sr, vf_wr, vf_sw, vf_rw, vf_ww, em_roof, em_wall, &
            em_improad, em_perroad, thick_roof, thick_wall, t_bld_min, t_bld_max, &
            eflx_traffic_factor, nlev_improad, alb_roof_dir, alb_roof_dif, &
            alb_wall_dir, alb_wall_dif, alb_improad_dir, alb_improad_dif, &
            alb_perroad_dir, alb_perroad_dif, tk_roof, tk_wall, tk_improad, &
            cv_roof, cv_wall, cv_improad)
    end if
    urban_params_built = .false.
    urban_built = .false.
    n_urb_col = 0
  end subroutine elmxx_urban_clean

end module elmxxUrbanMod
