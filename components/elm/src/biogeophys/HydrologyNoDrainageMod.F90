Module HydrologyNoDrainageMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Calculate snow and soil temperatures including phase change
  !
  use shr_kind_mod      , only : r8 => shr_kind_r8
  use shr_log_mod       , only : errMsg => shr_log_errMsg
  use decompMod         , only : bounds_type
  use elm_varctl        , only : iulog, use_vichydro, use_extrasnowlayers, use_firn_percolation_and_compaction
  use elm_varcon        , only : e_ice, denh2o, denice, rpi, spval
  use atm2lndType       , only : atm2lnd_type
  use ocn2lndType       , only : ocn2lnd_type
  use lnd2atmType       , only : lnd2atm_type
  use AerosolType       , only : aerosol_type
  use EnergyFluxType    , only : energyflux_type
  use CanopyStateType   , only  : canopystate_type
  use SoilHydrologyType , only : soilhydrology_type
  use SoilStateType     , only : soilstate_type
  use LandunitType      , only : lun_pp
  use ColumnType        , only : col_pp
  use ColumnDataType    , only : col_es, col_ws, col_wf
  use VegetationType    , only : veg_pp
  use TopounitDataType  , only : top_as, top_af ! Atmospheric state and flux variables
  use elm_instMod       , only : alm_fates , ep_betr
  use ElmDiagnostics    , only : diag_snapshot_1d, diag_snapshot_2d, diag_snapshot_int_1d, diag_snapshot_scalar, diag_enabled
  use VegetationDataType, only : veg_wf


  use timeinfoMod
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public  :: HydrologyNoDrainage ! Calculates soil/snow hydrology without drainage
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine HydrologyNoDrainage(bounds, &
       num_nolakec, filter_nolakec, &
       num_hydrologyc, filter_hydrologyc, &
       num_hydrononsoic, filter_hydrononsoic, &
       num_urbanc, filter_urbanc, &
       num_snowc, filter_snowc, &
       num_nosnowc, filter_nosnowc, canopystate_vars, &
       atm2lnd_vars, ocn2lnd_vars, lnd2atm_vars, soilstate_vars, &
       energyflux_vars, soilhydrology_vars, aerosol_vars)
    ! !DESCRIPTION:
    ! This is the main subroutine to execute the calculation of soil/snow
    ! hydrology
    ! Calling sequence is:
    !    -> SnowWater:             change of snow mass and snow water onto soil
    !    -> SurfaceRunoff:         surface runoff
    !    -> Infiltration:          infiltration into surface soil layer
    !    -> SoilWater:             soil water movement between layers
    !          -> Tridiagonal      tridiagonal matrix solution
    !    -> Drainage:              subsurface runoff
    !    -> SnowCompaction:        compaction of snow layers
    !    -> CombineSnowLayers:     combine snow layers that are thinner than minimum
    !    -> DivideSnowLayers:      subdivide snow layers that are thicker than maximum
    !    -> DivideExtraSnowLayers: subdivide up to 16 snow layers that are thicker than maximum
    !
    ! !USES:
      !$acc routine seq
    use elm_varcon           , only : denh2o, denice, hfus, grav, tfrz
    use landunit_varcon      , only : istice, istwet, istsoil, istice_mec, istcrop, istdlak
    use column_varcon        , only : icol_roof, icol_road_imperv, icol_road_perv, icol_sunwall
    use column_varcon        , only : icol_shadewall
    use elm_varctl           , only : use_cn, use_betr, use_fates, use_pflotran, pf_hmode, use_fan
    use elm_varpar           , only : nlevgrnd, nlevsno, nlevsoi, nlevurb
    use SnowHydrologyMod     , only : SnowCompaction, CombineSnowLayers, DivideSnowLayers, DivideExtraSnowLayers, SnowCapping
    use SnowHydrologyMod     , only : SnowWater, BuildSnowFilter 
    use SoilHydrologyMod     , only : ELMVICMap, SurfaceRunoff, Infiltration, WaterTable
    use SoilWaterMovementMod , only : SoilWater
    use SoilWaterRetentionCurveMod, only : soil_water_retention_curve_type
    use elm_varctl           , only : use_vsfm
    use SoilHydrologyMod     , only : DrainageVSFM
    use SoilWaterMovementMod , only : Compute_EffecRootFrac_And_VertTranSink
    !
    ! !ARGUMENTS:
    type(bounds_type)        , intent(in)    :: bounds
    integer                  , intent(in)    :: num_nolakec          ! number of column non-lake points in column filter
    integer                  , intent(in)    :: filter_nolakec(:)    ! column filter for non-lake points
    integer                  , intent(in)    :: num_hydrologyc       ! number of column soil points in column filter
    integer                  , intent(in)    :: filter_hydrologyc(:) ! column filter for soil points
    integer                  , intent(in)    :: num_hydrononsoic        ! number of non-soil landunit points in hydrology filter
    integer                  , intent(in)    :: filter_hydrononsoic(:)  ! column filter for non-soil hydrology points
    integer                  , intent(in)    :: num_urbanc           ! number of column urban points in column filter
    integer                  , intent(in)    :: filter_urbanc(:)     ! column filter for urban points
    integer                  , intent(inout) :: num_snowc            ! number of column snow points
    integer                  , intent(inout) :: filter_snowc(:)      ! column filter for snow points
    integer                  , intent(inout) :: num_nosnowc          ! number of column non-snow points
    integer                  , intent(inout) :: filter_nosnowc(:)    ! column filter for non-snow points
    type(atm2lnd_type)       , intent(in)    :: atm2lnd_vars
    type(ocn2lnd_type)       , intent(in)    :: ocn2lnd_vars
    type(lnd2atm_type)       , intent(in)    :: lnd2atm_vars
    type(soilstate_type)     , intent(inout) :: soilstate_vars
    type(energyflux_type)    , intent(in)    :: energyflux_vars
    type(canopystate_type)   , intent(in)  :: canopystate_vars
    type(aerosol_type)       , intent(inout) :: aerosol_vars
    type(soilhydrology_type) , intent(inout) :: soilhydrology_vars
    real(r8) :: dtime                         ! land model time step (sec)

    !
    ! !LOCAL VARIABLES:
    integer  :: g,l,c,j,fc                    ! indices
    integer  :: nlevbed                       ! # layers to bedrock
    real(r8) :: psi,vwc,fsattmp,psifrz        ! temporary variables for soilpsi calculation
    real(r8) :: watdry                        ! temporary
    real(r8) :: rwat(bounds%begc:bounds%endc) ! soil water wgted by depth to maximum depth of 0.5 m
    real(r8) :: swat(bounds%begc:bounds%endc) ! same as rwat but at saturation
    real(r8) :: rz(bounds%begc:bounds%endc)   ! thickness of soil layers contributing to rwat (m)
    real(r8) :: h2osoi_liq_saved(bounds%begc:bounds%endc) ! h2osoi_liq_col in topmost layer before calling SoilWater
    real(r8) :: tsw                           ! volumetric soil water to 0.5 m
    real(r8) :: stsw                          ! volumetric soil water to 0.5 m at saturation
    real(r8) :: fracl                         ! fraction of soil layer contributing to 10cm total soil water
    real(r8) :: s_node                        ! soil wetness (-)
    real(r8) :: icefrac(bounds%begc:bounds%endc,1:nlevgrnd)
    !-----------------------------------------------------------------------
    
    associate(                                                          & 
         z                  => col_pp%z                                  , & ! Input:  [real(r8) (:,:) ]  layer depth  (m)                      
         dz                 => col_pp%dz                                 , & ! Input:  [real(r8) (:,:) ]  layer thickness depth (m)             
         zi                 => col_pp%zi                                 , & ! Input:  [real(r8) (:,:) ]  interface depth (m)                   
         snl                => col_pp%snl                                , & ! Input:  [integer  (:)   ]  number of snow layers                    
         nlev2bed           => col_pp%nlevbed                           , & ! Input:  [integer  (:)   ]  number of layers to bedrock                     
         ctype              => col_pp%itype                              , & ! Input:  [integer  (:)   ]  column type                              

         forc_wind          => top_as%windbot         , & ! Input:  [real(r8) (:) ]  atmospheric wind speed (m/s)

         t_h2osfc           => col_es%t_h2osfc          , & ! Input:  [real(r8) (:)   ]  surface water temperature               
         dTdz_top           => col_es%dTdz_top          , & ! Output: [real(r8) (:)   ]  temperature gradient in top layer (col) [K m-1] !
         snot_top           => col_es%snot_top          , & ! Output: [real(r8) (:)   ]  snow temperature in top layer (col) [K]
         t_soisno           => col_es%t_soisno          , & ! Output: [real(r8) (:,:) ]  soil temperature (Kelvin)
         t_grnd             => col_es%t_grnd            , & ! Output: [real(r8) (:)   ]  ground temperature (Kelvin)
         t_grnd_u           => col_es%t_grnd_u          , & ! Output: [real(r8) (:)   ]  Urban ground temperature (Kelvin)
         t_grnd_r           => col_es%t_grnd_r          , & ! Output: [real(r8) (:)   ]  Rural ground temperature (Kelvin)
         t_soi_10cm         => col_es%t_soi10cm         , & ! Output: [real(r8) (:)   ]  soil temperature in top 10cm of soil (Kelvin)
         tsoi17             => col_es%t_soi17cm         , & ! Output: [real(r8) (:)   ]  soil temperature in top 17cm of soil (Kelvin)

         snow_depth         => col_ws%snow_depth         , & ! Input:  [real(r8) (:)   ]  snow height of snow covered area (m)
         snowdp             => col_ws%snowdp             , & ! Input:  [real(r8) (:)   ]  gridcell averaged snow height (m)
         frac_sno_eff       => col_ws%frac_sno_eff       , & ! Input:  [real(r8) (:)   ]  eff.  snow cover fraction (col) [frc]
         frac_h2osfc        => col_ws%frac_h2osfc        , & ! Input:  [real(r8) (:)   ]  fraction of ground covered by surface water (0 to 1)
         begwb              => col_ws%begwb              , & ! Input:  [real(r8) (:)   ]  water mass begining of the time step
         snw_rds            => col_ws%snw_rds            , & ! Output: [real(r8) (:,:) ]  effective snow grain radius (col,lyr) [microns, m^-6]
         snw_rds_top        => col_ws%snw_rds_top        , & ! Output: [real(r8) (:)   ]  effective snow grain size, top layer(col) [microns]
         sno_liq_top        => col_ws%sno_liq_top        , & ! Output: [real(r8) (:)   ]  liquid water fraction in top snow layer (col) [frc]
         snowice            => col_ws%snowice            , & ! Output: [real(r8) (:)   ]  average snow ice lens
         snowliq            => col_ws%snowliq            , & ! Output: [real(r8) (:)   ]  average snow liquid water
         snow_persistence   => col_ws%snow_persistence   , & ! Output: [real(r8) (:)   ]  counter for length of time snow-covered
         h2osoi_liqice_10cm => col_ws%h2osoi_liqice_10cm , & ! Output: [real(r8) (:)   ]  liquid water + ice lens in top 10cm of soil (kg/m2)
         h2osoi_ice         => col_ws%h2osoi_ice         , & ! Output: [real(r8) (:,:) ]  ice lens (kg/m2)
         h2osoi_liq         => col_ws%h2osoi_liq         , & ! Output: [real(r8) (:,:) ]  liquid water (kg/m2)
         h2osoi_vol         => col_ws%h2osoi_vol         , & ! Output: [real(r8) (:,:) ]  volumetric soil water (0<=h2osoi_vol<=watsat) [m3/m3]
         h2osno_top         => col_ws%h2osno_top         , & ! Output: [real(r8) (:)   ]  mass of snow in top layer (col) [kg]
         wf                 => col_ws%wf                 , & ! Output: [real(r8) (:)   ]  soil water as frac. of whc for top 0.05 m
         wf2                => col_ws%wf2                , & ! Output: [real(r8) (:)   ]  soil water as frac. of whc for top 0.17 m
         h2osoi_liqvol      => col_ws%h2osoi_liqvol      , & ! Output: [real(r8) (:,:) ]  volumetric liquid water content
         h2osoi_icevol      => col_ws%h2osoi_icevol      , & ! Output: [real(r8) (:,:) ]  volumetric liquid water content
         air_vol            => col_ws%air_vol            , & ! Output: [real(r8) (:,:) ]  volumetric air porosity
         eff_porosity       => soilstate_vars%eff_porosity_col        , & ! Output: [real(r8) (:,:) ]  effective soil porosity

         watsat             => soilstate_vars%watsat_col              , & ! Input:  [real(r8) (:,:) ]  volumetric soil water at saturation (porosity)
         sucsat             => soilstate_vars%sucsat_col              , & ! Input:  [real(r8) (:,:) ]  minimum soil suction (mm)
         bsw                => soilstate_vars%bsw_col                 , & ! Input:  [real(r8) (:,:) ]  Clapp and Hornberger "b"
         smp_l              => soilstate_vars%smp_l_col               , & ! Input:  [real(r8) (:,:) ]  soil matrix potential [mm]
         smpmin             => soilstate_vars%smpmin_col              , & ! Input:  [real(r8) (:)   ]  restriction for min of soil potential (mm)
         soilpsi            => soilstate_vars%soilpsi_col               & ! Output: [real(r8) (:,:) ]  soil water potential in each soil layer (MPa)
         )

         dtime = dtime_mod

      ! ========== DIAGNOSTICS: HydrologyNoDrainage INPUTS ==========
      if (diag_enabled) then
         call diag_snapshot_int_1d('hydro_in:snl', snl(bounds%begc:bounds%endc), bounds%endc-bounds%begc+1)
         call diag_snapshot_2d('hydro_in:dz', dz(bounds%begc:bounds%endc, -nlevsno+1:nlevgrnd), &
              bounds%endc-bounds%begc+1, nlevsno+nlevgrnd)
         call diag_snapshot_2d('hydro_in:t_soisno', t_soisno(bounds%begc:bounds%endc, -nlevsno+1:nlevgrnd), &
              bounds%endc-bounds%begc+1, nlevsno+nlevgrnd)
         call diag_snapshot_2d('hydro_in:h2osoi_liq', h2osoi_liq(bounds%begc:bounds%endc, -nlevsno+1:nlevgrnd), &
              bounds%endc-bounds%begc+1, nlevsno+nlevgrnd)
         call diag_snapshot_2d('hydro_in:h2osoi_ice', h2osoi_ice(bounds%begc:bounds%endc, -nlevsno+1:nlevgrnd), &
              bounds%endc-bounds%begc+1, nlevsno+nlevgrnd)
         call diag_snapshot_2d('hydro_in:watsat', watsat(bounds%begc:bounds%endc, 1:nlevgrnd), &
              bounds%endc-bounds%begc+1, nlevgrnd)
         call diag_snapshot_2d('hydro_in:sucsat', sucsat(bounds%begc:bounds%endc, 1:nlevgrnd), &
              bounds%endc-bounds%begc+1, nlevgrnd)
         call diag_snapshot_2d('hydro_in:bsw', bsw(bounds%begc:bounds%endc, 1:nlevgrnd), &
              bounds%endc-bounds%begc+1, nlevgrnd)
         call diag_snapshot_1d('hydro_in:t_grnd', t_grnd(bounds%begc:bounds%endc), bounds%endc-bounds%begc+1)
         call diag_snapshot_1d('hydro_in:t_h2osfc', t_h2osfc(bounds%begc:bounds%endc), bounds%endc-bounds%begc+1)
         call diag_snapshot_1d('hydro_in:frac_sno_eff', frac_sno_eff(bounds%begc:bounds%endc), bounds%endc-bounds%begc+1)
         call diag_snapshot_1d('hydro_in:frac_h2osfc', frac_h2osfc(bounds%begc:bounds%endc), bounds%endc-bounds%begc+1)
         call diag_snapshot_1d('hydro_in:snow_depth', snow_depth(bounds%begc:bounds%endc), bounds%endc-bounds%begc+1)
         ! Additional flux inputs needed for Python validation (added for elmpy)
         call diag_snapshot_1d('hydro_in:qflx_top_soil', col_wf%qflx_top_soil(bounds%begc:bounds%endc), bounds%endc-bounds%begc+1)
         call diag_snapshot_1d('hydro_in:qflx_evap_grnd', col_wf%qflx_evap_grnd(bounds%begc:bounds%endc), bounds%endc-bounds%begc+1)
         call diag_snapshot_1d('hydro_in:qflx_tran_veg', col_wf%qflx_tran_veg(bounds%begc:bounds%endc), bounds%endc-bounds%begc+1)
         call diag_snapshot_1d('hydro_in:qflx_dew_grnd', col_wf%qflx_dew_grnd(bounds%begc:bounds%endc), bounds%endc-bounds%begc+1)
         call diag_snapshot_1d('hydro_in:qflx_dew_snow', col_wf%qflx_dew_snow(bounds%begc:bounds%endc), bounds%endc-bounds%begc+1)
         call diag_snapshot_1d('hydro_in:qflx_sub_snow', col_wf%qflx_sub_snow(bounds%begc:bounds%endc), bounds%endc-bounds%begc+1)
         call diag_snapshot_1d('hydro_in:forc_rain', top_af%rain(bounds%begt:bounds%endt), bounds%endt-bounds%begt+1)
         call diag_snapshot_1d('hydro_in:forc_snow', top_af%snow(bounds%begt:bounds%endt), bounds%endt-bounds%begt+1)
      end if
      ! ================================================================

      ! Determine initial snow/no-snow filters (will be modified possibly by
      ! routines CombineSnowLayers and DivideSnowLayers below

      call BuildSnowFilter(bounds, num_nolakec, filter_nolakec, &
           num_snowc, filter_snowc, num_nosnowc, filter_nosnowc)

      ! Direct SnowHydrology contract: records the state consumed by
      ! SnowWater before downstream hydrology changes it.
      if (diag_enabled) then
         call diag_snapshot_int_1d('snowhydro_in:snl', snl(bounds%begc:bounds%endc), &
              bounds%endc-bounds%begc+1)
         call diag_snapshot_1d('snowhydro_in:frac_sno', col_ws%frac_sno(bounds%begc:bounds%endc), &
              bounds%endc-bounds%begc+1)
         call diag_snapshot_1d('snowhydro_in:frac_sno_eff', frac_sno_eff(bounds%begc:bounds%endc), &
              bounds%endc-bounds%begc+1)
         call diag_snapshot_1d('snowhydro_in:h2osno', col_ws%h2osno(bounds%begc:bounds%endc), &
              bounds%endc-bounds%begc+1)
         call diag_snapshot_1d('snowhydro_in:int_snow', col_ws%int_snow(bounds%begc:bounds%endc), &
              bounds%endc-bounds%begc+1)
         call diag_snapshot_1d('snowhydro_in:snow_depth', snow_depth(bounds%begc:bounds%endc), &
              bounds%endc-bounds%begc+1)
         call diag_snapshot_2d('snowhydro_in:dz', dz(bounds%begc:bounds%endc, -nlevsno+1:0), &
              bounds%endc-bounds%begc+1, nlevsno)
         call diag_snapshot_2d('snowhydro_in:h2osoi_liq', h2osoi_liq(bounds%begc:bounds%endc, -nlevsno+1:0), &
              bounds%endc-bounds%begc+1, nlevsno)
         call diag_snapshot_2d('snowhydro_in:h2osoi_ice', h2osoi_ice(bounds%begc:bounds%endc, -nlevsno+1:0), &
              bounds%endc-bounds%begc+1, nlevsno)
         call diag_snapshot_1d('snowhydro_in:qflx_snomelt', col_wf%qflx_snomelt(bounds%begc:bounds%endc), &
              bounds%endc-bounds%begc+1)
         call diag_snapshot_1d('snowhydro_in:qflx_rain_grnd', col_wf%qflx_rain_grnd(bounds%begc:bounds%endc), &
              bounds%endc-bounds%begc+1)
         call diag_snapshot_1d('snowhydro_in:qflx_sub_snow', col_wf%qflx_sub_snow(bounds%begc:bounds%endc), &
              bounds%endc-bounds%begc+1)
         call diag_snapshot_1d('snowhydro_in:qflx_dew_snow', col_wf%qflx_dew_snow(bounds%begc:bounds%endc), &
              bounds%endc-bounds%begc+1)
         call diag_snapshot_1d('snowhydro_in:qflx_evap_grnd', col_wf%qflx_evap_grnd(bounds%begc:bounds%endc), &
              bounds%endc-bounds%begc+1)
         call diag_snapshot_1d('snowhydro_in:qflx_dew_grnd', col_wf%qflx_dew_grnd(bounds%begc:bounds%endc), &
              bounds%endc-bounds%begc+1)
      end if

      ! Determine the change of snow mass and the snow water onto soil

      call SnowWater(bounds, num_snowc, filter_snowc, num_nosnowc, filter_nosnowc, &
           atm2lnd_vars, aerosol_vars)

      ! Diagnostic snapshot after SnowWater
      if (diag_enabled) then
         call diag_snapshot_1d('hydro_after_snow:qflx_top_soil', &
              col_wf%qflx_top_soil(bounds%begc:bounds%endc), bounds%endc-bounds%begc+1)
         call diag_snapshot_1d('hydro_after_snow:qflx_rain_grnd', &
              col_wf%qflx_rain_grnd(bounds%begc:bounds%endc), bounds%endc-bounds%begc+1)
         call diag_snapshot_1d('hydro_after_snow:qflx_snow_melt', &
              col_wf%qflx_snomelt(bounds%begc:bounds%endc), bounds%endc-bounds%begc+1)

         ! Direct SnowWater contract, before soil hydrology or snow-layer
         ! compaction/combine/divide alter the snow state.
         call diag_snapshot_int_1d('snowwater_out:snl', snl(bounds%begc:bounds%endc), &
              bounds%endc-bounds%begc+1)
         call diag_snapshot_1d('snowwater_out:frac_sno', col_ws%frac_sno(bounds%begc:bounds%endc), &
              bounds%endc-bounds%begc+1)
         call diag_snapshot_1d('snowwater_out:frac_sno_eff', frac_sno_eff(bounds%begc:bounds%endc), &
              bounds%endc-bounds%begc+1)
         call diag_snapshot_1d('snowwater_out:h2osno', col_ws%h2osno(bounds%begc:bounds%endc), &
              bounds%endc-bounds%begc+1)
         call diag_snapshot_1d('snowwater_out:int_snow', col_ws%int_snow(bounds%begc:bounds%endc), &
              bounds%endc-bounds%begc+1)
         call diag_snapshot_1d('snowwater_out:snow_depth', snow_depth(bounds%begc:bounds%endc), &
              bounds%endc-bounds%begc+1)
         call diag_snapshot_2d('snowwater_out:dz', dz(bounds%begc:bounds%endc, -nlevsno+1:0), &
              bounds%endc-bounds%begc+1, nlevsno)
         call diag_snapshot_2d('snowwater_out:h2osoi_liq', h2osoi_liq(bounds%begc:bounds%endc, -nlevsno+1:0), &
              bounds%endc-bounds%begc+1, nlevsno)
         call diag_snapshot_2d('snowwater_out:h2osoi_ice', h2osoi_ice(bounds%begc:bounds%endc, -nlevsno+1:0), &
              bounds%endc-bounds%begc+1, nlevsno)
         call diag_snapshot_1d('snowwater_out:qflx_snow_melt', col_wf%qflx_snow_melt(bounds%begc:bounds%endc), &
              bounds%endc-bounds%begc+1)
         call diag_snapshot_1d('snowwater_out:qflx_top_soil', col_wf%qflx_top_soil(bounds%begc:bounds%endc), &
              bounds%endc-bounds%begc+1)
      end if

      ! ========== DIAGNOSTICS: SurfaceRunoff + Infiltration INPUTS ==========
      if (diag_enabled) then
         ! Soil hydraulic properties (not in hydro_in)
         call diag_snapshot_2d('surfrun_in:hksat', &
              soilstate_vars%hksat_col(bounds%begc:bounds%endc, 1:nlevgrnd), &
              bounds%endc-bounds%begc+1, nlevgrnd)
         call diag_snapshot_1d('surfrun_in:wtfact', &
              soilstate_vars%wtfact_col(bounds%begc:bounds%endc), bounds%endc-bounds%begc+1)
         ! Soilhydrology state
         call diag_snapshot_1d('surfrun_in:zwt', &
              soilhydrology_vars%zwt_col(bounds%begc:bounds%endc), bounds%endc-bounds%begc+1)
         call diag_snapshot_1d('surfrun_in:zwt_perched', &
              soilhydrology_vars%zwt_perched_col(bounds%begc:bounds%endc), bounds%endc-bounds%begc+1)
         call diag_snapshot_1d('surfrun_in:frost_table', &
              soilhydrology_vars%frost_table_col(bounds%begc:bounds%endc), bounds%endc-bounds%begc+1)
         call diag_snapshot_1d('surfrun_in:h2osfc_thresh', &
              soilhydrology_vars%h2osfc_thresh_col(bounds%begc:bounds%endc), bounds%endc-bounds%begc+1)
         ! Column properties
         call diag_snapshot_int_1d('surfrun_in:nlevbed', &
              col_pp%nlevbed(bounds%begc:bounds%endc), bounds%endc-bounds%begc+1)
         ! snl is needed to locate the top soil layer in the combined
         ! snow+soil arrays; it is captured in hydro_in too, but SurfRunInfil
         ! consumers read the surfrun_in group on its own.
         call diag_snapshot_int_1d('surfrun_in:snl', &
              snl(bounds%begc:bounds%endc), bounds%endc-bounds%begc+1)
         call diag_snapshot_1d('surfrun_in:topo_slope', &
              col_pp%topo_slope(bounds%begc:bounds%endc), bounds%endc-bounds%begc+1)
         ! Gridcell-level fover (decay factor for fsat)
         ! fover is dimensioned by gridcell: fover(bounds%begg:bounds%endg)
         ! Snapshot gridcell-level array
         call diag_snapshot_1d('surfrun_in:fover', &
              soilhydrology_vars%fover(col_pp%gridcell(bounds%begc):col_pp%gridcell(bounds%endc)), &
              col_pp%gridcell(bounds%endc)-col_pp%gridcell(bounds%begc)+1)
         ! Fluxes not in hydro_in
         call diag_snapshot_1d('surfrun_in:qflx_floodc', &
              col_wf%qflx_floodc(bounds%begc:bounds%endc), bounds%endc-bounds%begc+1)
         call diag_snapshot_1d('surfrun_in:qflx_snow_h2osfc', &
              col_wf%qflx_snow_h2osfc(bounds%begc:bounds%endc), bounds%endc-bounds%begc+1)
         call diag_snapshot_1d('surfrun_in:qflx_ev_soil', &
              col_wf%qflx_ev_soil(bounds%begc:bounds%endc), bounds%endc-bounds%begc+1)
         call diag_snapshot_1d('surfrun_in:qflx_ev_h2osfc', &
              col_wf%qflx_ev_h2osfc(bounds%begc:bounds%endc), bounds%endc-bounds%begc+1)
         ! Water state
         call diag_snapshot_1d('surfrun_in:h2osfc', &
              col_ws%h2osfc(bounds%begc:bounds%endc), bounds%endc-bounds%begc+1)
         call diag_snapshot_1d('surfrun_in:frac_h2osfc_act', &
              col_ws%frac_h2osfc_act(bounds%begc:bounds%endc), bounds%endc-bounds%begc+1)
         ! Filters
         call diag_snapshot_int_1d('surfrun_in:filter_hydrologyc', &
              filter_hydrologyc(1:num_hydrologyc), num_hydrologyc)
         call diag_snapshot_int_1d('surfrun_in:filter_urbanc', &
              filter_urbanc(1:num_urbanc), num_urbanc)
         ! Scalar config
         if (soilhydrology_vars%origflag == 1) then
            call diag_snapshot_scalar('surfrun_in:origflag', 1.0_r8)
         else
            call diag_snapshot_scalar('surfrun_in:origflag', 0.0_r8)
         end if
         call diag_snapshot_scalar('surfrun_in:e_ice', e_ice)
         call diag_snapshot_scalar('surfrun_in:dtime', dtime)
         ! State consumed by SurfaceRunoff/Infiltration. These also exist in the
         ! hydro_in group, but are snapshotted again here on purpose: hydro_in is
         ! taken at the top of HydrologyNoDrainage, while these are the values as
         ! they stand immediately before SurfaceRunoff runs, after SnowWater and
         ! the preceding updates have already modified them.
         ! Soil layers only (1:nlevgrnd), unlike the hydro_in snapshots of the
         ! same arrays which span -nlevsno+1:nlevgrnd. SurfaceRunoff and
         ! Infiltration work purely in soil-layer index space -- snow is already
         ! accounted for upstream in qflx_top_soil -- so including the snow
         ! layers here would offset every soil layer by nlevsno.
         call diag_snapshot_2d('surfrun_in:dz', dz(bounds%begc:bounds%endc, 1:nlevgrnd), &
              bounds%endc-bounds%begc+1, nlevgrnd)
         call diag_snapshot_2d('surfrun_in:h2osoi_liq', h2osoi_liq(bounds%begc:bounds%endc, 1:nlevgrnd), &
              bounds%endc-bounds%begc+1, nlevgrnd)
         call diag_snapshot_2d('surfrun_in:h2osoi_ice', h2osoi_ice(bounds%begc:bounds%endc, 1:nlevgrnd), &
              bounds%endc-bounds%begc+1, nlevgrnd)
         call diag_snapshot_2d('surfrun_in:watsat', watsat(bounds%begc:bounds%endc, 1:nlevgrnd), &
              bounds%endc-bounds%begc+1, nlevgrnd)
         call diag_snapshot_1d('surfrun_in:frac_sno_eff', &
              frac_sno_eff(bounds%begc:bounds%endc), bounds%endc-bounds%begc+1)
         call diag_snapshot_1d('surfrun_in:frac_h2osfc', &
              frac_h2osfc(bounds%begc:bounds%endc), bounds%endc-bounds%begc+1)
         call diag_snapshot_1d('surfrun_in:qflx_top_soil', &
              col_wf%qflx_top_soil(bounds%begc:bounds%endc), bounds%endc-bounds%begc+1)
         call diag_snapshot_1d('surfrun_in:qflx_evap_grnd', &
              col_wf%qflx_evap_grnd(bounds%begc:bounds%endc), bounds%endc-bounds%begc+1)
      end if
      ! ====================================================================

      ! mapping soilmoist from CLM to VIC layers for runoff calculations
      if (use_vichydro) then
         call ELMVICMap(bounds, num_hydrologyc, filter_hydrologyc, &
              soilhydrology_vars)
      end if

      call SurfaceRunoff(bounds, num_hydrologyc, filter_hydrologyc, num_urbanc, filter_urbanc, &
           soilhydrology_vars, soilstate_vars, dtime)

      !------------------------------------------------------------------------------------
      if (use_pflotran .and. pf_hmode) then

        call Infiltration(bounds, num_hydrononsoic, filter_hydrononsoic,          &
             num_urbanc, filter_urbanc, atm2lnd_vars, ocn2lnd_vars, lnd2atm_vars, &
             energyflux_vars, soilhydrology_vars, soilstate_vars, dtime)

      else
      !------------------------------------------------------------------------------------

        call Infiltration(bounds, num_hydrologyc, filter_hydrologyc,              &
             num_urbanc, filter_urbanc, atm2lnd_vars, ocn2lnd_vars, lnd2atm_vars, &
             energyflux_vars, soilhydrology_vars, soilstate_vars, dtime)

      !------------------------------------------------------------------------------------
      end if
      !------------------------------------------------------------------------------------

      ! ========== DIAGNOSTICS: SurfaceRunoff + Infiltration OUTPUTS ==========
      if (diag_enabled) then
         ! SurfaceRunoff outputs
         call diag_snapshot_1d('surfrun_out:qflx_surf', &
              col_wf%qflx_surf(bounds%begc:bounds%endc), bounds%endc-bounds%begc+1)
         call diag_snapshot_1d('surfrun_out:fsat', &
              soilhydrology_vars%fsat_col(bounds%begc:bounds%endc), bounds%endc-bounds%begc+1)
         call diag_snapshot_1d('surfrun_out:fcov', &
              soilhydrology_vars%fcov_col(bounds%begc:bounds%endc), bounds%endc-bounds%begc+1)
         call diag_snapshot_2d('surfrun_out:fracice', &
              soilhydrology_vars%fracice_col(bounds%begc:bounds%endc, 1:nlevgrnd), &
              bounds%endc-bounds%begc+1, nlevgrnd)
         call diag_snapshot_2d('surfrun_out:icefrac', &
              soilhydrology_vars%icefrac_col(bounds%begc:bounds%endc, 1:nlevgrnd), &
              bounds%endc-bounds%begc+1, nlevgrnd)
         call diag_snapshot_1d('surfrun_out:qflx_top_soil', &
              col_wf%qflx_top_soil(bounds%begc:bounds%endc), bounds%endc-bounds%begc+1)
         ! Infiltration outputs
         call diag_snapshot_1d('surfrun_out:qflx_infl', &
              col_wf%qflx_infl(bounds%begc:bounds%endc), bounds%endc-bounds%begc+1)
         call diag_snapshot_1d('surfrun_out:qflx_h2osfc_surf', &
              col_wf%qflx_h2osfc_surf(bounds%begc:bounds%endc), bounds%endc-bounds%begc+1)
         call diag_snapshot_1d('surfrun_out:qflx_gross_infl_soil', &
              col_wf%qflx_gross_infl_soil(bounds%begc:bounds%endc), bounds%endc-bounds%begc+1)
         call diag_snapshot_1d('surfrun_out:qflx_gross_evap_soil', &
              col_wf%qflx_gross_evap_soil(bounds%begc:bounds%endc), bounds%endc-bounds%begc+1)
         call diag_snapshot_1d('surfrun_out:h2osfc', &
              col_ws%h2osfc(bounds%begc:bounds%endc), bounds%endc-bounds%begc+1)
         call diag_snapshot_2d('surfrun_out:eff_porosity', &
              eff_porosity(bounds%begc:bounds%endc, 1:nlevgrnd), &
              bounds%endc-bounds%begc+1, nlevgrnd)
      end if
      ! ====================================================================

      !!TODO:  need to fix the waterstate_vars dependence here.
#ifndef _OPENACC
      if (use_betr) then
        call ep_betr%BeTRSetBiophysForcing(bounds, col_pp, veg_pp, 1, nlevsoi, waterstate_vars=col_ws)
        call ep_betr%PreDiagSoilColWaterFlux(num_hydrologyc, filter_hydrologyc)
      endif
#endif

      if (use_vsfm) then
         call DrainageVSFM(bounds, num_hydrologyc, filter_hydrologyc, &
              num_urbanc, filter_urbanc,&
              soilhydrology_vars, soilstate_vars, dtime)
      endif

      call Compute_EffecRootFrac_And_VertTranSink(bounds, num_hydrologyc, &
           filter_hydrologyc, soilstate_vars, canopystate_vars, energyflux_vars)

      ! ========== DIAGNOSTICS: RootWaterUpdate INPUTS/OUTPUTS ==========
      if (diag_enabled) then
         ! Patch-level inputs
         call diag_snapshot_2d('rootwu_in:rootr_patch', &
              soilstate_vars%rootr_patch(bounds%begp:bounds%endp, 1:nlevgrnd), &
              bounds%endp-bounds%begp+1, nlevgrnd)
         call diag_snapshot_1d('rootwu_in:qflx_tran_veg_patch', &
              veg_wf%qflx_tran_veg(bounds%begp:bounds%endp), bounds%endp-bounds%begp+1)
         call diag_snapshot_1d('rootwu_in:qflx_tran_veg_col', &
              col_wf%qflx_tran_veg(bounds%begc:bounds%endc), bounds%endc-bounds%begc+1)
         ! Patch-column mapping
         call diag_snapshot_int_1d('rootwu_in:npfts', &
              col_pp%npfts(bounds%begc:bounds%endc), bounds%endc-bounds%begc+1)
         call diag_snapshot_int_1d('rootwu_in:pfti', &
              col_pp%pfti(bounds%begc:bounds%endc), bounds%endc-bounds%begc+1)
         call diag_snapshot_int_1d('rootwu_in:patch_column', &
              veg_pp%column(bounds%begp:bounds%endp), bounds%endp-bounds%begp+1)
         call diag_snapshot_1d('rootwu_in:patch_wtcol', &
              veg_pp%wtcol(bounds%begp:bounds%endp), bounds%endp-bounds%begp+1)
         ! Column-level outputs
         call diag_snapshot_2d('rootwu_out:rootr_col', &
              soilstate_vars%rootr_col(bounds%begc:bounds%endc, 1:nlevgrnd), &
              bounds%endc-bounds%begc+1, nlevgrnd)
         call diag_snapshot_2d('rootwu_out:qflx_rootsoi', &
              col_wf%qflx_rootsoi(bounds%begc:bounds%endc, 1:nlevgrnd), &
              bounds%endc-bounds%begc+1, nlevgrnd)
      end if
      ! ====================================================================

      if ( use_fan ) then 
         ! save the h2osoi_liq in top layer before evaluating the soilwater movement
         call store_tsl_moisture(filter_hydrologyc, num_hydrologyc) 
      end if
 
#ifndef _OPENACC
      ! If FATES plant hydraulics is turned on, over-ride default transpiration sink calculation
      if( use_fates ) call alm_fates%ComputeRootSoilFlux(bounds, num_hydrologyc, filter_hydrologyc, &
                                                      soilstate_vars)
#endif
      !------------------------------------------------------------------------------------
      if (use_pflotran .and. pf_hmode) then

        call SoilWater(bounds, num_hydrononsoic, filter_hydrononsoic, &
            num_urbanc, filter_urbanc, &
            soilhydrology_vars, soilstate_vars, dtime)

      else
      !------------------------------------------------------------------------------------

        call SoilWater(bounds, num_hydrologyc, filter_hydrologyc, num_urbanc, filter_urbanc, &
            soilhydrology_vars, soilstate_vars, dtime)

      !------------------------------------------------------------------------------------
      end if
      !------------------------------------------------------------------------------------

      if ( use_fan ) then 
         ! use the saved value to calculate the tendency
         call eval_tsl_moist_tend(filter_hydrologyc, num_hydrologyc)
      end if
 
#ifndef _OPENACC
       if (use_betr) then
          call ep_betr%BeTRSetBiophysForcing(bounds, col_pp, veg_pp, 1, nlevsoi, waterstate_vars=col_ws, &
             waterflux_vars=col_wf, soilhydrology_vars = soilhydrology_vars)

          call ep_betr%DiagAdvWaterFlux(num_hydrologyc, filter_hydrologyc)

          call ep_betr%RetrieveBiogeoFlux(bounds, 1, nlevsoi, waterflux_vars=col_wf)
       endif
#endif

      if (use_vichydro) then
         ! mapping soilmoist from CLM to VIC layers for runoff calculations
         call ELMVICMap(bounds, num_hydrologyc, filter_hydrologyc, &
              soilhydrology_vars)
      end if

      !------------------------------------------------------------------------------------
      if (use_pflotran .and. pf_hmode) then

        call WaterTable(bounds, num_hydrononsoic, filter_hydrononsoic, &
           num_urbanc, filter_urbanc, &
           soilhydrology_vars, soilstate_vars, dtime)

      else
      !------------------------------------------------------------------------------------

        call WaterTable(bounds, num_hydrologyc, filter_hydrologyc, num_urbanc, filter_urbanc, &
           soilhydrology_vars, soilstate_vars, dtime)

      !------------------------------------------------------------------------------------
      end if
      !------------------------------------------------------------------------------------


#ifndef _OPENACC
      if (use_betr) then
         !apply dew and sublimation fluxes, this is a temporary work aroud for tracking water isotope
         !Jinyun Tang, Feb 4, 2015
         call ep_betr%CalcDewSubFlux(bounds, col_pp, num_hydrologyc, filter_hydrologyc)
      endif           
#endif
      
      if (use_firn_percolation_and_compaction) then
         call SnowCapping(bounds, num_nolakec, filter_nolakec, num_snowc, filter_snowc, &
                          aerosol_vars)
      end if
      
      ! Natural compaction and metamorphosis.
      call SnowCompaction(bounds, num_snowc, filter_snowc, top_as, dtime)

      ! Combine thin snow elements
      call CombineSnowLayers(bounds, num_snowc, filter_snowc, &
           aerosol_vars, dtime)

      ! Divide thick snow elements
      if (.not. use_extrasnowlayers) then
         call DivideSnowLayers(bounds, num_snowc, filter_snowc, &
              aerosol_vars,  is_lake=.false.)
      else
         call DivideExtraSnowLayers(bounds, num_snowc, filter_snowc, &
              aerosol_vars,  is_lake=.false.)
      endif
      
      ! Set empty snow layers to zero
      do j = -nlevsno+1,0
         do fc = 1, num_snowc
            c = filter_snowc(fc)
            if (j <= snl(c) .and. snl(c) > -nlevsno) then
               h2osoi_ice(c,j) = 0._r8
               h2osoi_liq(c,j) = 0._r8
               t_soisno(c,j)  = 0._r8
               dz(c,j)    = 0._r8
               z(c,j)     = 0._r8
               zi(c,j-1)  = 0._r8
            end if
         end do
      end do

      ! Build new snow filter

      call BuildSnowFilter(bounds, num_nolakec, filter_nolakec, &
           num_snowc, filter_snowc, num_nosnowc, filter_nosnowc)

      ! Complete SnowHydrology contract: SnowWater, the intervening soil-only
      ! hydrology, and compaction/layer combine/divide have all run. The
      ! intervening soil routines do not mutate this snow-layer state.
      if (diag_enabled) then
         call diag_snapshot_int_1d('snowhydro_out:snl', snl(bounds%begc:bounds%endc), &
              bounds%endc-bounds%begc+1)
         call diag_snapshot_1d('snowhydro_out:frac_sno', col_ws%frac_sno(bounds%begc:bounds%endc), &
              bounds%endc-bounds%begc+1)
         call diag_snapshot_1d('snowhydro_out:frac_sno_eff', frac_sno_eff(bounds%begc:bounds%endc), &
              bounds%endc-bounds%begc+1)
         call diag_snapshot_1d('snowhydro_out:h2osno', col_ws%h2osno(bounds%begc:bounds%endc), &
              bounds%endc-bounds%begc+1)
         call diag_snapshot_1d('snowhydro_out:snow_depth', snow_depth(bounds%begc:bounds%endc), &
              bounds%endc-bounds%begc+1)
         call diag_snapshot_1d('snowhydro_out:int_snow', col_ws%int_snow(bounds%begc:bounds%endc), &
              bounds%endc-bounds%begc+1)
         call diag_snapshot_2d('snowhydro_out:dz', dz(bounds%begc:bounds%endc, -nlevsno+1:0), &
              bounds%endc-bounds%begc+1, nlevsno)
         call diag_snapshot_2d('snowhydro_out:h2osoi_liq', h2osoi_liq(bounds%begc:bounds%endc, -nlevsno+1:0), &
              bounds%endc-bounds%begc+1, nlevsno)
         call diag_snapshot_2d('snowhydro_out:h2osoi_ice', h2osoi_ice(bounds%begc:bounds%endc, -nlevsno+1:0), &
              bounds%endc-bounds%begc+1, nlevsno)
         call diag_snapshot_2d('snowhydro_out:snw_rds', snw_rds(bounds%begc:bounds%endc, -nlevsno+1:0), &
              bounds%endc-bounds%begc+1, nlevsno)
         call diag_snapshot_1d('snowhydro_out:qflx_snow_melt', col_wf%qflx_snow_melt(bounds%begc:bounds%endc), &
              bounds%endc-bounds%begc+1)
         call diag_snapshot_1d('snowhydro_out:qflx_top_soil', col_wf%qflx_top_soil(bounds%begc:bounds%endc), &
              bounds%endc-bounds%begc+1)
      end if

      ! For columns where snow exists, accumulate 'time-covered-by-snow' counters.
      ! Otherwise, re-zero counter, since it is bareland

      do fc = 1, num_snowc
         c = filter_snowc(fc)
         snow_persistence(c) = snow_persistence(c) + dtime
      end do
      do fc = 1, num_nosnowc
         c = filter_nosnowc(fc)
         snow_persistence(c) = 0._r8
      enddo

      ! Vertically average t_soisno and sum of h2osoi_liq and h2osoi_ice
      ! over all snow layers for history output

      do fc = 1, num_nolakec
         c = filter_nolakec(fc)
         snowice(c) = 0._r8
         snowliq(c) = 0._r8
      end do

      do j = -nlevsno+1, 0
         do fc = 1, num_snowc
            c = filter_snowc(fc)
            if (j >= snl(c)+1) then
               snowice(c) = snowice(c) + h2osoi_ice(c,j)
               snowliq(c) = snowliq(c) + h2osoi_liq(c,j)
            end if
         end do
      end do

      ! Calculate column average snow depth
      do c = bounds%begc,bounds%endc
         snowdp(c) = snow_depth(c) * frac_sno_eff(c)
      end do

      ! Determine ground temperature, ending water balance and volumetric soil water
      ! Calculate soil temperature and total water (liq+ice) in top 10cm of soil
      ! Calculate soil temperature and total water (liq+ice) in top 17cm of soil
      do fc = 1, num_nolakec
         c = filter_nolakec(fc)
         l = col_pp%landunit(c)
         if (.not. lun_pp%urbpoi(l)) then
            t_soi_10cm(c) = 0._r8
            tsoi17(c) = 0._r8
            h2osoi_liqice_10cm(c) = 0._r8
         end if
      end do
      do fc = 1, num_nolakec
         c = filter_nolakec(fc)
	       nlevbed = nlev2bed(c)
         do j = 1, nlevbed
            l = col_pp%landunit(c)
            if (.not. lun_pp%urbpoi(l)) then
               ! soil T at top 17 cm added by F. Li and S. Levis
               if (zi(c,j) <= 0.17_r8) then
                  fracl = 1._r8
                  tsoi17(c) = tsoi17(c) + t_soisno(c,j)*dz(c,j)*fracl
               else
                  if (zi(c,j) > 0.17_r8 .and. zi(c,j-1) < 0.17_r8) then
                     fracl = (0.17_r8 - zi(c,j-1))/dz(c,j)
                     tsoi17(c) = tsoi17(c) + t_soisno(c,j)*dz(c,j)*fracl
                  end if
               end if

               if (zi(c,j) <= 0.1_r8) then
                  fracl = 1._r8
                  t_soi_10cm(c) = t_soi_10cm(c) + t_soisno(c,j)*dz(c,j)*fracl
                  h2osoi_liqice_10cm(c) = h2osoi_liqice_10cm(c) + &
                       (h2osoi_liq(c,j)+h2osoi_ice(c,j))* &
                       fracl
               else
                  if (zi(c,j) > 0.1_r8 .and. zi(c,j-1) < 0.1_r8) then
                     fracl = (0.1_r8 - zi(c,j-1))/dz(c,j)
                     t_soi_10cm(c) = t_soi_10cm(c) + t_soisno(c,j)*dz(c,j)*fracl
                     h2osoi_liqice_10cm(c) = h2osoi_liqice_10cm(c) + &
                          (h2osoi_liq(c,j)+h2osoi_ice(c,j))* &
                          fracl
                  end if
               end if
            end if
         end do
      end do

      ! TODO - if this block of code is moved out of here - the SoilHydrology
      ! will NOT effect t_grnd, t_grnd_u or t_grnd_r

      do fc = 1, num_nolakec

         c = filter_nolakec(fc)
         l = col_pp%landunit(c)

         ! t_grnd is weighted average of exposed soil and snow
         if (snl(c) < 0) then
            t_grnd(c) = frac_sno_eff(c) * t_soisno(c,snl(c)+1) &
                 + (1 - frac_sno_eff(c)- frac_h2osfc(c)) * t_soisno(c,1) &
                 + frac_h2osfc(c) * t_h2osfc(c)
         else
            t_grnd(c) = (1 - frac_h2osfc(c)) * t_soisno(c,1) + frac_h2osfc(c) * t_h2osfc(c)
         endif

         if (lun_pp%urbpoi(l)) then
            t_grnd_u(c) = t_soisno(c,snl(c)+1)
         else
            t_soi_10cm(c) = t_soi_10cm(c)/0.1_r8
            tsoi17(c) =  tsoi17(c)/0.17_r8         ! F. Li and S. Levis
         end if
         if (col_pp%is_soil(c) .or. col_pp%is_crop(c)) then
            t_grnd_r(c) = t_soisno(c,snl(c)+1)
         end if

      end do

      do j = 1, nlevgrnd
         do fc = 1, num_nolakec
            c = filter_nolakec(fc)
            if ((ctype(c) == icol_sunwall .or. ctype(c) == icol_shadewall &
                 .or. ctype(c) == icol_roof) .and. j > nlevurb) then
            else
               h2osoi_vol(c,j) = h2osoi_liq(c,j)/(dz(c,j)*denh2o) + h2osoi_ice(c,j)/(dz(c,j)*denice)
               h2osoi_liqvol(c,j) = h2osoi_liq(c,j)/(dz(c,j)*denh2o)
               h2osoi_icevol(c,j) = h2osoi_ice(c,j)/(dz(c,j)*denice)
               air_vol(c,j)       = max(1.e-4_r8,watsat(c,j) - h2osoi_vol(c,j))
               eff_porosity(c,j)  = max(0.01_r8,watsat(c,j) - h2osoi_ice(c,j)/(dz(c,j)*denice))

            end if
         end do
      end do

      if ( (use_cn .or. use_fates) .and. &
         .not.(use_pflotran .and. pf_hmode) ) then
         ! Update soilpsi.
         ! ZMS: Note this could be merged with the following loop updating smp_l in the future.
         do j = 1, nlevgrnd
            do fc = 1, num_hydrologyc
               c = filter_hydrologyc(fc)

               if (h2osoi_liq(c,j) > 0._r8) then

                  vwc = h2osoi_liq(c,j)/(dz(c,j)*denh2o)

                  ! the following limit set to catch very small values of
                  ! fractional saturation that can crash the calculation of psi

                  ! use the same contants used in the supercool so that psi for frozen soils is consistent
                  fsattmp = max(vwc/watsat(c,j), 0.001_r8)
                  psi = sucsat(c,j) * (-9.8e-6_r8) * (fsattmp)**(-bsw(c,j))  ! Mpa
                  soilpsi(c,j) = min(max(psi,-15.0_r8),0._r8)

               else
                  soilpsi(c,j) = -15.0_r8
               end if
            end do
         end do
      end if

      if (use_cn .or. use_fates) then
         ! Available soil water up to a depth of 0.05 m.
         ! Potentially available soil water (=whc) up to a depth of 0.05 m.
         ! Water content as fraction of whc up to a depth of 0.05 m.

         do fc = 1, num_hydrologyc
            c = filter_hydrologyc(fc)
            rwat(c) = 0._r8
            swat(c) = 0._r8
            rz(c)   = 0._r8
         end do

         do j = 1, nlevgrnd
            do fc = 1, num_hydrologyc
               c = filter_hydrologyc(fc)
               !if (z(c,j)+0.5_r8*dz(c,j) <= 0.5_r8) then
               if (z(c,j)+0.5_r8*dz(c,j) <= 0.05_r8) then
                  watdry = watsat(c,j) * (316230._r8/sucsat(c,j)) ** (-1._r8/bsw(c,j))
                  rwat(c) = rwat(c) + (h2osoi_vol(c,j)-watdry) * dz(c,j)
                  swat(c) = swat(c) + (watsat(c,j)    -watdry) * dz(c,j)
                  rz(c) = rz(c) + dz(c,j)
               end if
            end do
         end do

         do fc = 1, num_hydrologyc
            c = filter_hydrologyc(fc)
            if (rz(c) /= 0._r8) then
               tsw  = rwat(c)/rz(c)
               stsw = swat(c)/rz(c)
            else
               watdry = watsat(c,1) * (316230._r8/sucsat(c,1)) ** (-1._r8/bsw(c,1))
               tsw = h2osoi_vol(c,1) - watdry
               stsw = watsat(c,1) - watdry
            end if
            wf(c) = tsw/stsw
         end do

         do j = 1, nlevgrnd
            do fc = 1, num_hydrologyc
               c = filter_hydrologyc(fc)
               if (z(c,j)+0.5_r8*dz(c,j) <= 0.17_r8) then
                  watdry = watsat(c,j) * (316230._r8/sucsat(c,j)) ** (-1._r8/bsw(c,j))
                  rwat(c) = rwat(c) + (h2osoi_vol(c,j)-watdry) * dz(c,j)
                  swat(c) = swat(c) + (watsat(c,j)    -watdry) * dz(c,j)
                  rz(c) = rz(c) + dz(c,j)
               end if
            end do
         end do

         do fc = 1, num_hydrologyc
            c = filter_hydrologyc(fc)
            if (rz(c) /= 0._r8) then
               tsw  = rwat(c)/rz(c)
               stsw = swat(c)/rz(c)
            else
               watdry = watsat(c,1) * (316230._r8/sucsat(c,1)) ** (-1._r8/bsw(c,1))
               tsw = h2osoi_vol(c,1) - watdry
               stsw = watsat(c,1) - watdry
            end if
            wf2(c) = tsw/stsw
         end do
      end if

      ! top-layer diagnostics
      do fc = 1, num_snowc
         c = filter_snowc(fc)
         h2osno_top(c)  = h2osoi_ice(c,snl(c)+1) + h2osoi_liq(c,snl(c)+1)
      enddo

      ! Zero variables in columns without snow
      do fc = 1, num_nosnowc
         c = filter_nosnowc(fc)

         h2osno_top(c)      = 0._r8
         snw_rds(c,:)       = 0._r8
         ! top-layer diagnostics (spval is not averaged when computing history fields)
         snot_top(c)        = spval
         dTdz_top(c)        = spval
         snw_rds_top(c)     = spval
         sno_liq_top(c)     = spval
      end do

      ! ========== DIAGNOSTICS: HydrologyNoDrainage OUTPUTS ==========
      if (diag_enabled) then
         call diag_snapshot_2d('hydro_out:h2osoi_liq', h2osoi_liq(bounds%begc:bounds%endc, -nlevsno+1:nlevgrnd), &
              bounds%endc-bounds%begc+1, nlevsno+nlevgrnd)
         call diag_snapshot_2d('hydro_out:h2osoi_ice', h2osoi_ice(bounds%begc:bounds%endc, -nlevsno+1:nlevgrnd), &
              bounds%endc-bounds%begc+1, nlevsno+nlevgrnd)
         call diag_snapshot_2d('hydro_out:h2osoi_vol', h2osoi_vol(bounds%begc:bounds%endc, 1:nlevgrnd), &
              bounds%endc-bounds%begc+1, nlevgrnd)
         call diag_snapshot_2d('hydro_out:h2osoi_liqvol', h2osoi_liqvol(bounds%begc:bounds%endc, 1:nlevgrnd), &
              bounds%endc-bounds%begc+1, nlevgrnd)
         call diag_snapshot_2d('hydro_out:h2osoi_icevol', h2osoi_icevol(bounds%begc:bounds%endc, 1:nlevgrnd), &
              bounds%endc-bounds%begc+1, nlevgrnd)
         call diag_snapshot_2d('hydro_out:air_vol', air_vol(bounds%begc:bounds%endc, 1:nlevgrnd), &
              bounds%endc-bounds%begc+1, nlevgrnd)
         call diag_snapshot_2d('hydro_out:eff_porosity', eff_porosity(bounds%begc:bounds%endc, 1:nlevgrnd), &
              bounds%endc-bounds%begc+1, nlevgrnd)
         call diag_snapshot_1d('hydro_out:t_grnd', t_grnd(bounds%begc:bounds%endc), bounds%endc-bounds%begc+1)
         call diag_snapshot_1d('hydro_out:t_soi_10cm', t_soi_10cm(bounds%begc:bounds%endc), bounds%endc-bounds%begc+1)
         call diag_snapshot_1d('hydro_out:tsoi17', tsoi17(bounds%begc:bounds%endc), bounds%endc-bounds%begc+1)
         call diag_snapshot_1d('hydro_out:snowice', snowice(bounds%begc:bounds%endc), bounds%endc-bounds%begc+1)
         call diag_snapshot_1d('hydro_out:snowliq', snowliq(bounds%begc:bounds%endc), bounds%endc-bounds%begc+1)
         call diag_snapshot_1d('hydro_out:snowdp', snowdp(bounds%begc:bounds%endc), bounds%endc-bounds%begc+1)
         if (use_cn .or. use_fates) then
            call diag_snapshot_2d('hydro_out:soilpsi', soilpsi(bounds%begc:bounds%endc, 1:nlevgrnd), &
                 bounds%endc-bounds%begc+1, nlevgrnd)
            call diag_snapshot_1d('hydro_out:wf', wf(bounds%begc:bounds%endc), bounds%endc-bounds%begc+1)
            call diag_snapshot_1d('hydro_out:wf2', wf2(bounds%begc:bounds%endc), bounds%endc-bounds%begc+1)
         end if
      end if
      ! ================================================================

    end associate

contains

    ! Subroutines for storing the time derivative of top most soil layer
    ! moisture. This is used for diagnosing the downwards moisture flux within FAN.

    subroutine store_tsl_moisture(filter, num_fc)
      ! Store the soil water within topmost layer before evaluating soil
      ! moisture
      ! transport.
      integer, intent(in) :: filter(:)
      integer, intent(in) :: num_fc

      integer :: fc, c

      do fc = 1, num_fc
         c = filter(fc)
         h2osoi_liq_saved(c) = col_ws%h2osoi_liq(c,1)
      end do

    end subroutine store_tsl_moisture

    subroutine eval_tsl_moist_tend(filter, num_fc)
      ! Evaluate the time derivative of soil liquid water due to percolation as
      ! required in FAN.
      integer, intent(in) :: filter(:)
      integer, intent(in) :: num_fc

      integer :: fc, c

      associate(h2osoi_tend_tsl => col_ws%h2osoi_tend_tsl_col, &
           h2osoi_liq => col_ws%h2osoi_liq)

      do fc = 1, num_fc
         c = filter(fc)
         h2osoi_tend_tsl(c) = (h2osoi_liq(c,1) - h2osoi_liq_saved(c)) / dtime
      end do

      end associate

    end subroutine eval_tsl_moist_tend
 
  end subroutine HydrologyNoDrainage

end Module HydrologyNoDrainageMod
