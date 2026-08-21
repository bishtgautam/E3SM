module elmxxRootKernelMod
  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Argument-only root water stress kernel: ELM's
  ! SoilMoistStressMod calc_root_moist_stress_clm45default.
  !
  ! Split out of elmxxRootMod so a standalone replay can drive THE SAME SOURCE
  ! the coupled run calls, rather than a copy of it. elmxxRootMod keeps the
  ! gather/call/push around this.
  !
  ! Per layer, where there is liquid water and the soil is not too cold:
  !   eff_porosity = max(0.01, watsat - h2osoi_ice/(dz*denice))
  !   liqvol       = h2osoi_liq/(dz*denh2o)            -- NOT capped, per ELM
  !   s_node       = max(liqvol/eff_porosity, 0.01)
  !   smp_node     = max(smpsc, -sucsat*s_node**(-bsw))
  !   rresis       = min( (eff_porosity/watsat)*(smp_node-smpsc)
  !                       /(smpso-smpsc), 1 )
  !   rootr        = rootfr*rresis,  btran = sum(max(rootr,0))
  ! then rootr is normalised by btran so the layers partition the uptake.
  !-----------------------------------------------------------------------

  use shr_kind_mod , only : r8 => shr_kind_r8
  use shr_const_mod, only : SHR_CONST_TKFRZ

  implicit none
  private

  public :: elmxx_root_stress_kernel

contains

  subroutine elmxx_root_stress_kernel(np, nc, nlevbed, nlevgrnd,            &
                                      patch_itype, patch_col,              &
                                      rootfr, h2osoi_liq, h2osoi_ice, dz,  &
                                      t_soisno, watsat, bsw, sucsat,       &
                                      smpsc, smpso, tc_stress, btran0,     &
                                      denice, denh2o,                      &
                                      rootr, btran, rresis_out)
    implicit none
    integer , intent(in)  :: np, nc, nlevbed, nlevgrnd
    integer , intent(in)  :: patch_itype(np), patch_col(np)
    real(r8), intent(in)  :: rootfr(np, nlevgrnd)
    ! Soil state, indexed by ELM soil layer 1..nlevgrnd (no snow slots).
    real(r8), intent(in)  :: h2osoi_liq(nc, nlevgrnd), h2osoi_ice(nc, nlevgrnd)
    real(r8), intent(in)  :: dz(nc, nlevgrnd), t_soisno(nc, nlevgrnd)
    real(r8), intent(in)  :: watsat(nc, nlevgrnd), bsw(nc, nlevgrnd)
    real(r8), intent(in)  :: sucsat(nc, nlevgrnd)
    ! Per-PFT, indexed by patch_itype.
    real(r8), intent(in)  :: smpsc(0:), smpso(0:)
    real(r8), intent(in)  :: tc_stress, btran0, denice, denh2o
    real(r8), intent(out) :: rootr(np, nlevgrnd), btran(np)
    real(r8), intent(out) :: rresis_out(np, nlevgrnd)

    integer  :: p, c, j, ivt
    real(r8) :: tcold, eff_por, liqvol, s_node, smp_node, rresis

    tcold      = SHR_CONST_TKFRZ + tc_stress
    rootr      = 0.0_r8
    btran      = 0.0_r8
    rresis_out = 0.0_r8

    do p = 1, np
       ivt = patch_itype(p)
       if (ivt == 0) cycle            ! bare ground transpires nothing
       c = patch_col(p)

       do j = 1, nlevbed

          ! ELM floors eff_porosity at 0.01 (HydrologyNoDrainageMod, where the
          ! array the stress calc reads is actually set) rather than skipping
          ! the layer. Neither binds without ice, but the floor is what ELM does.
          eff_por = max(0.01_r8, watsat(c,j) - h2osoi_ice(c,j)/(dz(c,j)*denice))

          ! ELM does NOT cap the liquid volume at the effective porosity:
          !   h2osoi_liqvol(c,j) = h2osoi_liq(c,j)/(dz(c,j)*denh2o)
          liqvol = h2osoi_liq(c,j) / (dz(c,j)*denh2o)

          if (liqvol <= 0.0_r8 .or. t_soisno(c,j) <= tcold) cycle

          s_node   = max(liqvol/eff_por, 0.01_r8)
          smp_node = max(smpsc(ivt), -sucsat(c,j) * s_node**(-bsw(c,j)))

          rresis = min( (eff_por/watsat(c,j)) * (smp_node - smpsc(ivt))     &
                        / (smpso(ivt) - smpsc(ivt)), 1.0_r8 )

          rresis_out(p,j) = rresis
          rootr(p,j)      = rootfr(p,j) * rresis
          btran(p)        = btran(p) + max(rootr(p,j), 0.0_r8)
       end do

       ! Normalise so the layers partition the uptake rather than scale it.
       if (btran(p) > btran0) then
          rootr(p,1:nlevgrnd) = rootr(p,1:nlevgrnd) / btran(p)
       else
          rootr(p,1:nlevgrnd) = 0.0_r8
       end if
    end do

  end subroutine elmxx_root_stress_kernel

end module elmxxRootKernelMod
