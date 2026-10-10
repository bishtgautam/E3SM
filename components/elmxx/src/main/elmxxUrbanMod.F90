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
                               istsoil
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
                               ELMxxRestartFieldInfo, &
                               ELMxxHistorySetLandunitWeights, ELMxxHistorySetUrbanWeights

  implicit none
  save
  private

  integer, parameter :: numrad = 2         ! ELMxx NUMRAD (VIS, NIR)
  integer, parameter :: nlevurb_x = 5      ! ELMxx NLEVURB

  logical, public :: urban_params_built = .false.

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
  end subroutine elmxx_urban_clean

end module elmxxUrbanMod
