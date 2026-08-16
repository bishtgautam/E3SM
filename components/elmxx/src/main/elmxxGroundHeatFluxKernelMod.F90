module elmxxGroundHeatFluxKernelMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! The ground surface energy balance, as a pure kernel: arrays in, arrays out.
  !
  ! ELMxx has no kernel for this. SoilTemperature READS hs_soil, hs_top_snow
  ! and hs_h2osfc and never writes them; nothing else in the library writes
  ! them either. In ELM they are built by ComputeGroundHeatFluxAndDeriv inside
  ! SoilTemperatureMod, from the radiation and turbulent fluxes the earlier
  ! kernels produced. So this is Fortran-side work, exactly as btran was, and
  ! only the result crosses.
  !
  ! WHY IT IS A SEPARATE MODULE. Same reason as elmxxSurfaceAlbedoKernelMod:
  ! nothing here depends on anything but its arguments, so
  ! tools/ground_heat_flux_replay.F90 can drive it with ELM's own recorded
  ! inputs and compare against ELM's own recorded outputs. The wrapper --
  ! elmxxSoilKernelMod -- keeps the ELMxx object, the maps and the setters,
  ! none of which can be linked into a small offline program.
  !
  ! Ported from ELM SoilTemperatureMod.F90:1789-2047, ComputeGroundHeatFluxAndDeriv,
  ! NON-URBAN branch only. The urban branch uses eflx_lwrad_net and adds
  ! wasteheat, air-conditioning and traffic fluxes; it belongs with UrbanAlbedo
  ! and the urban kernels, which are blocked for other reasons.
  !
  ! DELIBERATELY NOT COMPUTED:
  !   hs, hs_top   ELM forms both alongside these. ELMxx's C++ SoilTemperature
  !                reads only hs_soil, hs_top_snow and hs_h2osfc -- there is no
  !                hs_top view in NaturalColumnData to set -- so computing them
  !                would produce numbers with no consumer. Checked, not assumed.
  !-----------------------------------------------------------------------

  use shr_kind_mod, only : r8 => shr_kind_r8

  implicit none
  save
  private

  public :: elmxx_ground_heat_flux_kernel

contains

  !-----------------------------------------------------------------------
  subroutine elmxx_ground_heat_flux_kernel(nc, np, nlevsno, nlevtot, nsnw, sb, &
       patch_col, patch_wt, frac_veg_nosno,                                    &
       emg, htvp, t_grnd, t_h2osfc, snl, t_soisno, forc_lwrad,                 &
       sabg_soil, sabg_lyr, dlrad, cgrnd,                                      &
       eflx_sh_snow, eflx_sh_soil, eflx_sh_h2osfc,                             &
       qflx_ev_snow, qflx_ev_soil, qflx_ev_h2osfc,                             &
       hs_soil, hs_top_snow, hs_h2osfc, dhsdT, sabg_lyr_col, n_nonfinite)
    !
    ! Column-mean ground heat fluxes over nc columns and np patches.
    !
    ! patch_col(p) is the 1-based column index of patch p; <= 0 means the patch
    ! is not on one of these columns and is skipped. Weights are ELM's
    ! veg_pp%wtcol and are what turns per-patch fluxes into a column mean.
    !
    ! LAYER INDEXING, and it is the trap here. ELM's t_soisno runs
    ! -nlevsno+1..nlevgrnd, so ELM layer j lives in slot nlevsno + j:
    ! soil layer 1 is slot nlevsno+1, and the top ACTIVE layer, ELM's snl+1,
    ! is slot nlevsno+snl+1. sabg_lyr is narrower -- ELM's -nlevsno+1..1, so
    ! nsnw = nlevsno+1 slots -- but shares the same origin, so the same
    ! slot arithmetic indexes it.
    !
    implicit none
    integer , intent(in)  :: nc, np, nlevsno, nlevtot, nsnw
    real(r8), intent(in)  :: sb                          ! Stefan-Boltzmann
    integer , intent(in)  :: patch_col(np)
    real(r8), intent(in)  :: patch_wt(np)
    integer , intent(in)  :: frac_veg_nosno(np)
    real(r8), intent(in)  :: emg(nc), htvp(nc), t_grnd(nc), t_h2osfc(nc)
    integer , intent(in)  :: snl(nc)
    real(r8), intent(in)  :: t_soisno(nc, nlevtot)
    real(r8), intent(in)  :: forc_lwrad(nc)
    real(r8), intent(in)  :: sabg_soil(np)
    real(r8), intent(in)  :: sabg_lyr(np, nsnw)
    real(r8), intent(in)  :: dlrad(np), cgrnd(np)
    real(r8), intent(in)  :: eflx_sh_snow(np), eflx_sh_soil(np), eflx_sh_h2osfc(np)
    real(r8), intent(in)  :: qflx_ev_snow(np), qflx_ev_soil(np), qflx_ev_h2osfc(np)
    real(r8), intent(out) :: hs_soil(nc), hs_top_snow(nc), hs_h2osfc(nc)
    real(r8), intent(out) :: dhsdT(nc)
    real(r8), intent(out) :: sabg_lyr_col(nc, nlevtot)
    integer , intent(out) :: n_nonfinite     ! patches producing a non-finite flux
    !
    integer  :: p, c, j, lyr_top
    real(r8) :: w, fvn, lw_in, gnet_soil, gnet_h2osfc, gnet_snow
    real(r8) :: dle, le_soil, le_snow, le_h2osfc
    !-----------------------------------------------------------------------

    hs_soil = 0.0_r8; hs_top_snow = 0.0_r8; hs_h2osfc = 0.0_r8
    dhsdT = 0.0_r8; sabg_lyr_col = 0.0_r8
    n_nonfinite = 0

    do p = 1, np
       c = patch_col(p)
       if (c <= 0) cycle
       w = patch_wt(p)

       lyr_top = nlevsno + snl(c) + 1
       fvn     = real(frac_veg_nosno(p), r8)

       ! Emitted longwave, fractionated. ELM balances these against
       ! CanopyFluxes and Biogeophysics2, so the three surfaces each emit at
       ! their own temperature rather than at t_grnd.
       dle       = 4.0_r8 * emg(c) * sb * t_grnd(c)**3
       le_soil   = emg(c) * sb * t_soisno(c, nlevsno + 1)**4
       le_snow   = emg(c) * sb * t_soisno(c, lyr_top)**4
       le_h2osfc = emg(c) * sb * t_h2osfc(c)**4

       ! Incident longwave reaching the ground: only the fraction the canopy
       ! does not cover. Under a full canopy this is zero and dlrad carries
       ! the downward longwave instead.
       lw_in = (1.0_r8 - fvn) * emg(c) * forc_lwrad(c)

       gnet_soil   = sabg_soil(p) + dlrad(p) + lw_in - le_soil &
                   - (eflx_sh_soil(p) + qflx_ev_soil(p) * htvp(c))
       gnet_h2osfc = sabg_soil(p) + dlrad(p) + lw_in - le_h2osfc &
                   - (eflx_sh_h2osfc(p) + qflx_ev_h2osfc(p) * htvp(c))

       ! hs_top_snow takes its shortwave from the TOP ACTIVE LAYER of the
       ! per-layer absorption, not from the bulk sabg_snow. ELM computes an
       ! eflx_gnet_snow from sabg_snow in its first loop and then throws it
       ! away; the value that reaches hs_top_snow is this one, out of its
       ! second loop (:2020).
       gnet_snow   = sabg_lyr(p, lyr_top) + dlrad(p) + lw_in - le_snow &
                   - (eflx_sh_snow(p) + qflx_ev_snow(p) * htvp(c))

       ! Non-finite guard. The canopy kernels can hand back NaN on vegetated
       ! patches when their own inputs are unset, and a NaN here silently
       ! poisons the whole column mean -- minval/maxval would not show it.
       if (.not. (abs(gnet_soil)   >= 0.0_r8) .or. &
           .not. (abs(gnet_snow)   >= 0.0_r8) .or. &
           .not. (abs(gnet_h2osfc) >= 0.0_r8)) then
          n_nonfinite = n_nonfinite + 1
       end if

       hs_soil(c)     = hs_soil(c)     + gnet_soil   * w
       hs_h2osfc(c)   = hs_h2osfc(c)   + gnet_h2osfc * w
       hs_top_snow(c) = hs_top_snow(c) + gnet_snow   * w
       dhsdT(c)       = dhsdT(c)       + (-cgrnd(p) - dle) * w

       ! Column-mean absorbed solar by layer, top active layer down to the
       ! first soil layer -- ELM's "do j = lyr_top,1,1".
       do j = lyr_top, nlevsno + 1
          sabg_lyr_col(c, j) = sabg_lyr_col(c, j) + sabg_lyr(p, j) * w
       end do
    end do

  end subroutine elmxx_ground_heat_flux_kernel

end module elmxxGroundHeatFluxKernelMod
