module elmxxSoilKernelMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Stage 4: the integration surface for the five soil/hydrology kernels --
  ! SoilTemperature, SoilFluxes, SurfRunInfil, RootWaterUpdate and
  ! HydrologyDrainage.
  !
  ! THIS IS A MUCH SMALLER SURFACE THAN STATUS ONCE CLAIMED, and the reason is
  ! worth stating because it changed the plan. Each of these kernels has TWO
  ! entry points in ELMxx:
  !
  !   ELMxxComputeSoilTemperature         -- the STANDALONE path. Runs on a
  !     separate allocation made by ELMxxInitST, filled through 127 distinct
  !     ST_/SF_/SRI_/RWU_/HD_ setters. This is what the validation driver uses
  !     to replay ELM diagnostic snapshots kernel by kernel.
  !
  !   ELMxxComputeSoilTemperatureNatural  -- the INTEGRATED path. Builds its
  !     view directly from elm->naturalCol / elm->naturalPatch, the same state
  !     the canopy kernels already read and write.
  !
  ! An integrated driver wants the second. Decision #8 already said "using the
  ! ...Natural kernel variants"; what was not appreciated is that this deletes
  ! the entire per-kernel setter surface. There is no ELMxxInitST, no ST_
  ! setters, no separate index space to map. What is left is: seed the
  ! naturalCol views these kernels read and the canopy kernels do not, declare
  ! the shared filters, and supply the one thing ELMxx genuinely does not
  ! compute -- the ground surface energy balance.
  !
  ! THREE PARALLEL LAYER REPRESENTATIONS, AND THEY ARE NOT INTERCHANGEABLE.
  ! This is the trap that makes this module long:
  !
  !   name        width          layout
  !   ----------  -------------  --------------------------------------------
  !   plain       NLEVTOT   = 20 snow[1..5] soil[6..20]
  !   _p1         NLEVTOT_P1= 21 snow[1..5] SSW[6] soil[7..21]
  !   _soi        NLEVGRND  = 15 soil only
  !
  ! SoilTemperature reads the _p1 views. SurfRunInfil and HydrologyDrainage
  ! read the _soi views. SoilFluxes reads the plain ones. ELMxxSetDz writes the
  ! PLAIN view; ELMxxSetDzP1 and ELMxxSetDzSoi are separate calls. So the
  ! seeding that satisfied the canopy kernels leaves the other two
  ! representations AT ZERO -- which a kernel reads happily.
  !
  ! The _p1 layout inserts a standing-surface-water slot between the snow and
  ! the soil (SoilTemperatureImpl.h: "snow[0..4] | SSW[5] | soil[6..20]"), so
  ! soil layer j sits one slot further down than in the plain view. An
  ! off-by-one here is not a crash; it is a shifted soil column.
  !-----------------------------------------------------------------------

  use shr_kind_mod        , only : r8 => shr_kind_r8
  use shr_sys_mod         , only : shr_sys_abort, shr_sys_flush
  use shr_const_mod       , only : SHR_CONST_STEBOL
  use elmxxSpmdMod        , only : masterproc, iam
  use elmxxSubgridMod     , only : num_columns, num_patches, &
                                   col_landunit, lun_gridcell, patch_column, &
                                   col_itype, lun_itype, patch_wtcol, &
                                   istsoil, istcrop
  use elmxxSurfaceStateMod, only : patch_lai, patch_sai
  use elmxxForcingMod     , only : forc_lwrad
  use elmxxGroundHeatFluxKernelMod, only : elmxx_ground_heat_flux_kernel
  use elmxxSoilPropMod    , only : nlevsoi, nlevgrnd, nlevsno, nlevtot, nlevbed, &
                                   zsoi, dzsoi, zisoi, &
                                   sp_watsat  => watsat, sp_hksat => hksat, &
                                   sp_tkmg    => tkmg,   sp_tkdry => tkdry, &
                                   sp_tksatu  => tksatu, sp_csol  => csol, &
                                   sp_dz      => col_dz, &
                                   sp_tsoisno => col_t_soisno, &
                                   sp_liq     => col_h2osoi_liq, &
                                   sp_ice     => col_h2osoi_ice, &
                                   sp_vol     => col_h2osoi_vol, &
                                   soil_prop_built
  use elmxxSurfdataMod    , only : topo_slope, fmax
  use elmxxKokkosStateMod , only : n_kokkos_col, n_kokkos_patch, &
                                   col_of_kcol, kcol_of_col, &
                                   patch_of_kpatch, kpatch_of_patch, &
                                   kokkos_state_built
  use elmxx_mod           , only : ELMxxType, ELMXX_SUCCESS, &
                                   ELMxxInitSharedMetadata, &
                                   ELMxxSetSharedFilterNolakec, &
                                   ELMxxSetSharedFilterNolakep, &
                                   ELMxxSetSharedFilterHydrologyc, &
                                   ELMxxSetSharedFilterUrbanc, &
                                   ELMxxSetSharedUrbpoi, ELMxxSetSharedFover, &
                                   ELMxxSetSharedForcRain, ELMxxSetSharedForcSnow, &
                                   ELMxxSetSharedQflxFloodg, &
                                   ELMxxSetColItype, ELMxxSetColIsSoil, &
                                   ELMxxSetColIsCrop, ELMxxSetColActive, &
                                   ELMxxSetLunItype, ELMxxSetLunUrbpoi, &
                                   ELMxxSetColTopounit, ELMxxSetColGridcell, &
                                   ELMxxSetColNpfts, ELMxxSetColPfti, &
                                   ELMxxSetNlevbed, ELMxxSetZP1, ELMxxSetZiP1, &
                                   ELMxxSetDzP1, ELMxxSetTSoisnoP1, &
                                   ELMxxSetH2osoiLiqP1, ELMxxSetH2osoiIceP1, &
                                   ELMxxSetDzSoi, ELMxxSetWatsatSoi, &
                                   ELMxxSetZcSoi, ELMxxSetZiSoi, &
                                   ELMxxSetH2osoiLiqSoi, ELMxxSetH2osoiIceSoi, &
                                   ELMxxSetH2osoiVol, &
                                   ELMxxSetTkmg, ELMxxSetTkdry, ELMxxSetTksatu, &
                                   ELMxxSetCsol, ELMxxSetHksat, &
                                   ELMxxSetWtfact, ELMxxSetZwt, ELMxxSetZwtPerched, &
                                   ELMxxSetFrostTable, ELMxxSetH2osfcThresh, &
                                   ELMxxSetTopoSlope, ELMxxSetFracH2osfcAct, &
                                   ELMxxSetDzH2osfc, ELMxxSetEflxBot, &
                                   ELMxxSetHsSoil, ELMxxSetHsTopSnow, &
                                   ELMxxSetHsH2osfc, ELMxxSetDhsdT, &
                                   ELMxxSetSabgLyrCol, ELMxxSetTssbef, &
                                   ELMxxGetQflxTranVeg, ELMxxSetQflxTranVegCol, &
                                   ELMxxSetWtcol, ELMxxSetPatchActive, &
                                   ELMxxSetPatchLandunit, &
                                   ELMxxSetIsOnSoilCol, ELMxxSetIsOnCropCol, &
                                   ELMxxSetWa, ELMxxSetBegwb, &
                                   ELMxxSetTotalPlantStoredH2o, &
                                   ELMxxSetQflxIrrig, ELMxxSetQflxGlciceFrz, &
                                   ELMxxSetQflxTopSoil, ELMxxSetQflxFloodc, &
                                   ELMxxSetQflxSnowH2osfc, &
                                   ELMxxSetH2ocanCol, ELMxxSetFSurfCol, &
                                   ELMxxGetSabg, ELMxxGetSabgSoil, ELMxxGetSabgSnow, &
                                   ELMxxGetSabgLyr, ELMxxGetDlrad, ELMxxGetCgrnd, &
                                   ELMxxGetEflxShGrnd, ELMxxGetEflxShSoil, &
                                   ELMxxGetEflxShSnow, ELMxxGetEflxShH2osfc, &
                                   ELMxxGetQflxEvapSoi, ELMxxGetQflxEvSoil, &
                                   ELMxxGetQflxEvSnow, ELMxxGetQflxEvH2osfc, &
                                   ELMxxGetQflxRainGrnd, &
                                   ELMxxGetH2osoiLiqSoi, ELMxxGetH2osoiIceSoi, &
                                   ELMxxGetEmg, ELMxxGetHtvp, ELMxxGetTGrnd, &
                                   ELMxxGetTH2osfc, ELMxxGetTSoisno, ELMxxGetSnl

  implicit none
  save
  private

  ! nSnowLyr: patch-level sabg_lyr is (P, NLEVSNO+1) -- ELM's -nlevsno+1..1.
  integer, parameter :: nsnowlyr   = nlevsno + 1        ! 6
  integer, parameter :: nlevtot_p1 = nlevsno + 1 + nlevgrnd  ! 21

  ! Packed column -> first packed patch (0-based) and patch count. Built once
  ! and checked for contiguity, which the kernels' pfti/npfts idiom requires.
  integer, allocatable :: kcol_pfti(:)
  integer, allocatable :: kcol_npfts(:)

  logical, public :: soil_kernel_built = .false.

  public :: elmxx_soil_kernel_init
  public :: elmxx_soil_kernel_clean

contains

  !-----------------------------------------------------------------------
  subroutine elmxx_soil_kernel_init(elm, dtime, logunit)
    !
    ! Everything the five kernels read that does not change with time.
    !
    implicit none
    type(ELMxxType), intent(in) :: elm
    real(r8), intent(in) :: dtime
    integer , intent(in) :: logunit
    integer :: kc, kp, c, p, l, j, ierr, sz(2)
    integer , allocatable :: icol(:), ipatch(:), filt(:)
    real(r8), allocatable :: rcol(:), rpatch(:), buf(:,:)
    character(len=*), parameter :: subname = '(elmxx_soil_kernel_init) '

    if (.not. kokkos_state_built) call shr_sys_abort(subname//'ERROR: maps not built')
    if (.not. soil_prop_built)    call shr_sys_abort(subname//'ERROR: soil properties not built')

    call build_column_patch_index(logunit)

    allocate(icol(n_kokkos_col), rcol(n_kokkos_col), filt(n_kokkos_col))
    allocate(ipatch(n_kokkos_patch), rpatch(n_kokkos_patch))

    !-----------------------------------------------------------------
    ! Shared metadata and the filters.
    !
    ! NATURAL-ONLY (decision #8), so the declaration is deliberately narrow:
    ! nolakec and hydrologyc are every packed column, nolakep is every packed
    ! patch, and urbanc is EMPTY. That is not a stub -- the packed views carry
    ! no urban column at all, so an urban filter entry would index nothing.
    ! Urban arrives as its own increment with its own kernel variants
    ! (ELMxxComputeSoilFluxesUrban and friends), not by widening these.
    !-----------------------------------------------------------------
    call ELMxxInitSharedMetadata(elm, 1, num_gridcells_local(), &
         n_kokkos_col, n_kokkos_patch, n_kokkos_col, 0, dtime, ierr)
    call check(ierr, subname, 'InitSharedMetadata')

    ! Filters are 0-based indices into the packed spaces.
    do kc = 1, n_kokkos_col
       filt(kc) = kc - 1
    end do
    call ELMxxSetSharedFilterNolakec(elm, filt, n_kokkos_col, ierr)
    call check(ierr, subname, 'FilterNolakec')
    call ELMxxSetSharedFilterHydrologyc(elm, filt, n_kokkos_col, ierr)
    call check(ierr, subname, 'FilterHydrologyc')

    do kp = 1, n_kokkos_patch
       ipatch(kp) = kp - 1
    end do
    call ELMxxSetSharedFilterNolakep(elm, ipatch, n_kokkos_patch, ierr)
    call check(ierr, subname, 'FilterNolakep')

    ! urbanc is empty; the call is still made so the count is explicit rather
    ! than left to whatever InitSharedMetadata defaulted to.
    call ELMxxSetSharedFilterUrbanc(elm, filt, 0, ierr)
    call check(ierr, subname, 'FilterUrbanc')

    ! fover: ELM's runoff decay factor, 0.5 m-1 (SurfaceRunoffMod / hydrology
    ! namelist default). One value per gridcell.
    deallocate(rcol); allocate(rcol(max(1, num_gridcells_local())))
    rcol = 0.5_r8
    call ELMxxSetSharedFover(elm, rcol, num_gridcells_local(), ierr)
    call check(ierr, subname, 'Fover')
    deallocate(rcol); allocate(rcol(n_kokkos_col))

    !-----------------------------------------------------------------
    ! Column metadata.
    !-----------------------------------------------------------------
    do kc = 1, n_kokkos_col
       c = col_of_kcol(kc)
       icol(kc) = col_itype(c)
    end do
    call ELMxxSetColItype(elm, icol, n_kokkos_col, ierr); call check(ierr, subname, 'ColItype')

    ! is_soil / is_crop follow ELM's lun_itype test, not the column type.
    do kc = 1, n_kokkos_col
       l = col_landunit(col_of_kcol(kc))
       if (lun_itype(l) == istsoil) then; icol(kc) = 1; else; icol(kc) = 0; end if
    end do
    call ELMxxSetColIsSoil(elm, icol, n_kokkos_col, ierr); call check(ierr, subname, 'ColIsSoil')

    do kc = 1, n_kokkos_col
       l = col_landunit(col_of_kcol(kc))
       if (lun_itype(l) == istcrop) then; icol(kc) = 1; else; icol(kc) = 0; end if
    end do
    call ELMxxSetColIsCrop(elm, icol, n_kokkos_col, ierr); call check(ierr, subname, 'ColIsCrop')

    icol = 1
    call ELMxxSetColActive(elm, icol, n_kokkos_col, ierr); call check(ierr, subname, 'ColActive')

    do kc = 1, n_kokkos_col
       icol(kc) = lun_itype(col_landunit(col_of_kcol(kc)))
    end do
    call ELMxxSetLunItype(elm, icol, n_kokkos_col, ierr); call check(ierr, subname, 'LunItype')

    ! No urban among the packed columns, by construction.
    icol = 0
    call ELMxxSetLunUrbpoi(elm, icol, n_kokkos_col, ierr); call check(ierr, subname, 'LunUrbpoi')
    call ELMxxSetColTopounit(elm, icol, n_kokkos_col, ierr); call check(ierr, subname, 'ColTopounit')

    ! Gridcell index, 0-based, for the per-gridcell shared forcing.
    do kc = 1, n_kokkos_col
       icol(kc) = lun_gridcell(col_landunit(col_of_kcol(kc))) - 1
    end do
    call ELMxxSetColGridcell(elm, icol, n_kokkos_col, ierr); call check(ierr, subname, 'ColGridcell')

    call ELMxxSetColNpfts(elm, kcol_npfts, n_kokkos_col, ierr); call check(ierr, subname, 'ColNpfts')
    call ELMxxSetColPfti (elm, kcol_pfti , n_kokkos_col, ierr); call check(ierr, subname, 'ColPfti')

    icol = nlevbed
    call ELMxxSetNlevbed(elm, icol, n_kokkos_col, ierr); call check(ierr, subname, 'Nlevbed')

    !-----------------------------------------------------------------
    ! Patch metadata.
    !-----------------------------------------------------------------
    do kp = 1, n_kokkos_patch
       rpatch(kp) = patch_wtcol(patch_of_kpatch(kp))
    end do
    call ELMxxSetWtcol(elm, rpatch, n_kokkos_patch, ierr); call check(ierr, subname, 'Wtcol')

    ipatch = 1
    call ELMxxSetPatchActive(elm, ipatch, n_kokkos_patch, ierr); call check(ierr, subname, 'PatchActive')

    do kp = 1, n_kokkos_patch
       l = col_landunit(patch_column(patch_of_kpatch(kp)))
       if (lun_itype(l) == istsoil) then; ipatch(kp) = 1; else; ipatch(kp) = 0; end if
    end do
    call ELMxxSetIsOnSoilCol(elm, ipatch, n_kokkos_patch, ierr); call check(ierr, subname, 'IsOnSoilCol')

    do kp = 1, n_kokkos_patch
       l = col_landunit(patch_column(patch_of_kpatch(kp)))
       if (lun_itype(l) == istcrop) then; ipatch(kp) = 1; else; ipatch(kp) = 0; end if
    end do
    call ELMxxSetIsOnCropCol(elm, ipatch, n_kokkos_patch, ierr); call check(ierr, subname, 'IsOnCropCol')

    ! Landunit index is used only to look up urbpoi, which is uniformly zero
    ! here; 0 keeps every patch pointing at the same non-urban entry.
    ipatch = 0
    call ELMxxSetPatchLandunit(elm, ipatch, n_kokkos_patch, ierr); call check(ierr, subname, 'PatchLandunit')

    !-----------------------------------------------------------------
    ! The _p1 grid: 21 slots, snow[1..5] SSW[6] soil[7..21].
    !-----------------------------------------------------------------
    sz(1) = n_kokkos_col
    sz(2) = nlevtot_p1
    allocate(buf(n_kokkos_col, nlevtot_p1))

    buf = 0.0_r8
    do j = 1, nlevgrnd
       do kc = 1, n_kokkos_col
          buf(kc, nlevsno + 1 + j) = dzsoi(j)
       end do
    end do
    call ELMxxSetDzP1(elm, buf, sz, ierr); call check(ierr, subname, 'DzP1')

    buf = 0.0_r8
    do j = 1, nlevgrnd
       do kc = 1, n_kokkos_col
          buf(kc, nlevsno + 1 + j) = zsoi(j)
       end do
    end do
    call ELMxxSetZP1(elm, buf, sz, ierr); call check(ierr, subname, 'ZP1')

    buf = 0.0_r8
    do j = 1, nlevgrnd
       do kc = 1, n_kokkos_col
          buf(kc, nlevsno + 1 + j) = sp_tsoisno(col_of_kcol(kc), nlevsno + j)
       end do
    end do
    call ELMxxSetTSoisnoP1(elm, buf, sz, ierr); call check(ierr, subname, 'TSoisnoP1')

    buf = 0.0_r8
    do j = 1, nlevgrnd
       do kc = 1, n_kokkos_col
          buf(kc, nlevsno + 1 + j) = sp_liq(col_of_kcol(kc), nlevsno + j)
       end do
    end do
    call ELMxxSetH2osoiLiqP1(elm, buf, sz, ierr); call check(ierr, subname, 'H2osoiLiqP1')

    buf = 0.0_r8
    do j = 1, nlevgrnd
       do kc = 1, n_kokkos_col
          buf(kc, nlevsno + 1 + j) = sp_ice(col_of_kcol(kc), nlevsno + j)
       end do
    end do
    call ELMxxSetH2osoiIceP1(elm, buf, sz, ierr); call check(ierr, subname, 'H2osoiIceP1')
    deallocate(buf)

    ! zi_p1 has one more slot than z_p1: interfaces, not nodes.
    sz(2) = nlevtot_p1 + 1
    allocate(buf(n_kokkos_col, nlevtot_p1 + 1))
    buf = 0.0_r8
    do j = 0, nlevgrnd
       do kc = 1, n_kokkos_col
          buf(kc, nlevsno + 1 + j + 1) = zisoi(j)
       end do
    end do
    call ELMxxSetZiP1(elm, buf, sz, ierr); call check(ierr, subname, 'ZiP1')
    deallocate(buf)

    !-----------------------------------------------------------------
    ! The _soi arrays: soil layers only, 15 wide.
    !-----------------------------------------------------------------
    sz(2) = nlevgrnd
    allocate(buf(n_kokkos_col, nlevgrnd))

    do j = 1, nlevgrnd
       do kc = 1, n_kokkos_col
          buf(kc,j) = dzsoi(j)
       end do
    end do
    call ELMxxSetDzSoi(elm, buf, sz, ierr); call check(ierr, subname, 'DzSoi')

    ! Node depths on the SAME soil-only indexing. The Richards solve needs
    ! these separately from z_p1, which carries a standing-water slot between
    ! the snow and the soil and so puts layer j one slot further down.
    do j = 1, nlevgrnd
       do kc = 1, n_kokkos_col
          buf(kc,j) = zsoi(j)
       end do
    end do
    call ELMxxSetZcSoi(elm, buf, sz, ierr); call check(ierr, subname, 'ZcSoi')

    call fill_soi(buf, sp_watsat); call ELMxxSetWatsatSoi(elm, buf, sz, ierr)
    call check(ierr, subname, 'WatsatSoi')
    call fill_soi(buf, sp_hksat);  call ELMxxSetHksat(elm, buf, sz, ierr)
    call check(ierr, subname, 'Hksat')
    call fill_soi(buf, sp_tkmg);   call ELMxxSetTkmg(elm, buf, sz, ierr)
    call check(ierr, subname, 'Tkmg')
    call fill_soi(buf, sp_tkdry);  call ELMxxSetTkdry(elm, buf, sz, ierr)
    call check(ierr, subname, 'Tkdry')
    call fill_soi(buf, sp_tksatu); call ELMxxSetTksatu(elm, buf, sz, ierr)
    call check(ierr, subname, 'Tksatu')
    call fill_soi(buf, sp_csol);   call ELMxxSetCsol(elm, buf, sz, ierr)
    call check(ierr, subname, 'Csol')
    call fill_soi(buf, sp_vol);    call ELMxxSetH2osoiVol(elm, buf, sz, ierr)
    call check(ierr, subname, 'H2osoiVol')

    do j = 1, nlevgrnd
       do kc = 1, n_kokkos_col
          buf(kc,j) = sp_liq(col_of_kcol(kc), nlevsno + j)
       end do
    end do
    call ELMxxSetH2osoiLiqSoi(elm, buf, sz, ierr); call check(ierr, subname, 'H2osoiLiqSoi')

    do j = 1, nlevgrnd
       do kc = 1, n_kokkos_col
          buf(kc,j) = sp_ice(col_of_kcol(kc), nlevsno + j)
       end do
    end do
    call ELMxxSetH2osoiIceSoi(elm, buf, sz, ierr); call check(ierr, subname, 'H2osoiIceSoi')
    deallocate(buf)

    ! Soil interfaces: one more than the node count. zisoi is 0-based in
    ! elmxxSoilPropMod (0 = the surface), so slot k holds zisoi(k-1).
    sz(2) = nlevgrnd + 1
    allocate(buf(n_kokkos_col, nlevgrnd + 1))
    do j = 0, nlevgrnd
       do kc = 1, n_kokkos_col
          buf(kc, j+1) = zisoi(j)
       end do
    end do
    call ELMxxSetZiSoi(elm, buf, sz, ierr); call check(ierr, subname, 'ZiSoi')
    deallocate(buf)

    !-----------------------------------------------------------------
    ! Hydrology statics and cold-start water table.
    !
    ! ELM's cold start (SoilHydrologyType InitCold): see the block below for
    ! the water table. fmax comes from surfdata (FMAX), defaulting to 0.38
    ! where the field is absent.
    !-----------------------------------------------------------------
    do kc = 1, n_kokkos_col
       rcol(kc) = fmax(lun_gridcell(col_landunit(col_of_kcol(kc))))
    end do
    call ELMxxSetWtfact(elm, rcol, n_kokkos_col, ierr); call check(ierr, subname, 'Wtfact')

    do kc = 1, n_kokkos_col
       rcol(kc) = topo_slope(lun_gridcell(col_landunit(col_of_kcol(kc))))
    end do
    call ELMxxSetTopoSlope(elm, rcol, n_kokkos_col, ierr); call check(ierr, subname, 'TopoSlope')

    ! ELM, SoilHydrologyType InitCold, non-urban branch with
    ! use_var_soil_thick = .false. (the default):
    !     wa   = 4000
    !     zwt  = (25 + zi(nlevsoi)) - wa/0.2/1000  =  5 + zi(nlevsoi)
    !     zwt_perched = frost_table = zi(nlevsoi)
    !
    ! These were previously 4.8 / 4800 / zi(nlevbed), which mixes two branches
    ! ELM does not take here: wa = 4800 is the urban icol_road_perv value, and
    ! nlevbed belongs to the use_var_soil_thick branch. The zwt error is not
    ! cosmetic -- fsat = wtfact*exp(-0.5*fover*zwt), so a water table 4 m too
    ! shallow inflates the saturated fraction by e^(0.5*0.5*4) ~ 2.7 and sends
    ! that much extra water to surface runoff instead of into the soil.
    rcol = 4000.0_r8
    call ELMxxSetWa(elm, rcol, n_kokkos_col, ierr); call check(ierr, subname, 'Wa')
    rcol = (25.0_r8 + zisoi(nlevsoi)) - 4000.0_r8 / 0.2_r8 / 1000.0_r8
    call ELMxxSetZwt(elm, rcol, n_kokkos_col, ierr); call check(ierr, subname, 'Zwt')
    rcol = zisoi(nlevsoi)
    call ELMxxSetZwtPerched(elm, rcol, n_kokkos_col, ierr); call check(ierr, subname, 'ZwtPerched')
    call ELMxxSetFrostTable(elm, rcol, n_kokkos_col, ierr); call check(ierr, subname, 'FrostTable')

    ! Surface-water depth threshold: ELM's micro-topography relation, with the
    ! standard 1e-3 m minimum.
    rcol = 0.0_r8
    call ELMxxSetH2osfcThresh(elm, rcol, n_kokkos_col, ierr); call check(ierr, subname, 'H2osfcThresh')
    call ELMxxSetFracH2osfcAct(elm, rcol, n_kokkos_col, ierr); call check(ierr, subname, 'FracH2osfcAct')
    call ELMxxSetDzH2osfc(elm, rcol, n_kokkos_col, ierr); call check(ierr, subname, 'DzH2osfc')

    ! Zero at a cold start, and each has a real source later: eflx_bot is the
    ! geothermal flux (zero in this configuration), the irrigation and glacier
    ! terms belong to components that are out of scope, and total_plant_stored
    ! is a hydraulic-stress term that use_hydrstress = .false. turns off.
    rcol = 0.0_r8
    call ELMxxSetEflxBot(elm, rcol, n_kokkos_col, ierr); call check(ierr, subname, 'EflxBot')
    call ELMxxSetQflxIrrig(elm, rcol, n_kokkos_col, ierr); call check(ierr, subname, 'QflxIrrig')
    call ELMxxSetQflxGlciceFrz(elm, rcol, n_kokkos_col, ierr); call check(ierr, subname, 'QflxGlciceFrz')
    call ELMxxSetTotalPlantStoredH2o(elm, rcol, n_kokkos_col, ierr)
    call check(ierr, subname, 'TotalPlantStoredH2o')
    call ELMxxSetQflxFloodc(elm, rcol, n_kokkos_col, ierr); call check(ierr, subname, 'QflxFloodc')
    call ELMxxSetFSurfCol(elm, rcol, n_kokkos_col, ierr); call check(ierr, subname, 'FSurfCol')

    deallocate(icol, rcol, ipatch, rpatch, filt)
    soil_kernel_built = .true.

    if (masterproc) then
       write(logunit,*) subname,'rank ',iam,' shared-kernel surface ready: ', &
            n_kokkos_col,' columns ',n_kokkos_patch,' patches, urbanc empty'
       call shr_sys_flush(logunit)
    end if

  end subroutine elmxx_soil_kernel_init




  !-----------------------------------------------------------------------
  subroutine build_column_patch_index(logunit)
    !
    ! Packed column -> (first packed patch, count), 0-based first index.
    !
    ! CONTIGUITY IS CHECKED, NOT ASSUMED. The kernels walk a column's patches
    ! as pfti .. pfti+npfts-1, which is only meaningful if the packed patches
    ! of a column are consecutive. They are, because the subgrid fills patches
    ! per column and the packing preserves order -- but that is a property of
    ! two modules agreeing, and it would break silently if either changed.
    !
    implicit none
    integer, intent(in) :: logunit
    integer :: kc, kp, c, cprev, n
    character(len=*), parameter :: subname = '(build_column_patch_index) '

    allocate(kcol_pfti(n_kokkos_col), kcol_npfts(n_kokkos_col))
    kcol_pfti  = 0
    kcol_npfts = 0

    cprev = -1
    do kp = 1, n_kokkos_patch
       c  = patch_column(patch_of_kpatch(kp))
       kc = kcol_of_col(c) + 1
       if (kc <= 0) then
          call shr_sys_abort(subname//'ERROR: packed patch on an unpacked column')
       end if
       if (kc /= cprev) then
          if (kcol_npfts(kc) /= 0) then
             call shr_sys_abort(subname//'ERROR: packed patches of a column are '// &
                  'not contiguous; pfti/npfts cannot address them')
          end if
          kcol_pfti(kc) = kp - 1
          cprev = kc
       end if
       kcol_npfts(kc) = kcol_npfts(kc) + 1
    end do

    n = sum(kcol_npfts)
    if (n /= n_kokkos_patch) then
       call shr_sys_abort(subname//'ERROR: patch counts do not sum to the packed total')
    end if
    do kc = 1, n_kokkos_col
       if (kcol_npfts(kc) == 0) then
          call shr_sys_abort(subname//'ERROR: a packed column carries no patch')
       end if
    end do

    write(logunit,*) subname,'rank ',iam,' column->patch index built and checked; ', &
         'patches per column ',minval(kcol_npfts),' .. ',maxval(kcol_npfts)
    call shr_sys_flush(logunit)

  end subroutine build_column_patch_index

  !-----------------------------------------------------------------------
  subroutine fill_soi(buf, src)
    implicit none
    real(r8), intent(out) :: buf(:,:)
    real(r8), intent(in)  :: src(:,:)
    integer :: kc, j
    do j = 1, nlevgrnd
       do kc = 1, n_kokkos_col
          buf(kc,j) = src(col_of_kcol(kc), j)
       end do
    end do
  end subroutine fill_soi

  !-----------------------------------------------------------------------
  real(r8) function elai_of(p)
    implicit none
    integer, intent(in) :: p
    elai_of = patch_lai(p)
    if (elai_of < 0.05_r8) elai_of = 0.0_r8
  end function elai_of

  real(r8) function esai_of(p)
    implicit none
    integer, intent(in) :: p
    esai_of = patch_sai(p)
    if (esai_of < 0.05_r8) esai_of = 0.0_r8
  end function esai_of

  !-----------------------------------------------------------------------
  integer function num_gridcells_local()
    implicit none
    integer :: c, g
    num_gridcells_local = 0
    do c = 1, num_columns
       g = lun_gridcell(col_landunit(c))
       if (g > num_gridcells_local) num_gridcells_local = g
    end do
  end function num_gridcells_local

  !-----------------------------------------------------------------------
  subroutine check(ierr, subname, what)
    implicit none
    integer, intent(in) :: ierr
    character(len=*), intent(in) :: subname, what
    if (ierr /= ELMXX_SUCCESS) then
       call shr_sys_abort(subname//'ERROR: '//trim(what)//' returned a non-success status')
    end if
  end subroutine check

  !-----------------------------------------------------------------------
  subroutine elmxx_soil_kernel_clean()
    implicit none
    if (allocated(kcol_pfti))  deallocate(kcol_pfti)
    if (allocated(kcol_npfts)) deallocate(kcol_npfts)
    soil_kernel_built = .false.
  end subroutine elmxx_soil_kernel_clean

end module elmxxSoilKernelMod
