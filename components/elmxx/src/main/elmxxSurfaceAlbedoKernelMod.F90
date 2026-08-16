module elmxxSurfaceAlbedoKernelMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! The surface albedo physics, as a pure kernel: arrays in, arrays out.
  !
  ! WHY THIS IS A SEPARATE MODULE. elmxxSurfaceAlbedoMod knows about the
  ! ELMxx object, the subgrid maps, the PFT parameter module and the Kokkos
  ! setters. None of that can be linked into a small offline program, so as
  ! long as the physics lived there it could only be graded by running a
  ! coupled case and eyeballing the numbers -- which is how it was graded, and
  ! which is exactly the kind of grading STATUS section J says has now failed
  ! five times.
  !
  ! Everything here depends on nothing but its arguments, so
  ! tools/surface_albedo_replay.F90 can drive it with an ELM restart's own
  ! inputs and compare against that restart's own outputs. Same source, same
  ! code path as the coupled run: the test drives production code, not a copy
  ! of it.
  !
  ! Ported from ELM SurfaceAlbedoMod.F90:
  !   SurfaceAlbedo   (:240-1012) the initialization, filters and novegsol pass
  !   SoilAlbedo      (:1015)     soil colour + surface wetness
  !   TwoStream       (:1147)     Sellers/Dickinson, nlevcan = 1 branch
  ! with the albsat/albdry colour tables from SurfaceAlbedoType.F90:201.
  !
  ! DELIBERATELY NOT PORTED, each with its consequence stated:
  !   SNICAR                 snow albedo falls back to ELM's cold-start
  !                          constant. Free where frac_sno is zero; NOT free
  !                          at f19.
  !   Albedo_TOP_Adjustment  ELM gates it on use_top_solar_rad = .false.
  !   lake / glacier / urban branches of SoilAlbedo -- out of scope.
  !   vcmaxcintsun/sha       leaf-to-canopy scaling for Photosynthesis, which
  !                          is not ported. Add here when it is; ELM computes
  !                          it in this same nlevcan = 1 block.
  !-----------------------------------------------------------------------

  use shr_kind_mod, only : r8 => shr_kind_r8

  implicit none
  save
  private

  integer, parameter, public :: numrad = 2

  ! Snow two-stream parameters, ELM elm_varcon.
  real(r8), parameter :: omegas(numrad) = (/ 0.8_r8, 0.4_r8 /)
  real(r8), parameter :: betads = 0.5_r8
  real(r8), parameter :: betais = 0.5_r8
  real(r8), parameter :: tfrz   = 273.15_r8
  real(r8), parameter :: mpe    = 1.0e-6_r8

  ! Soil albedo by colour class and band, ELM SurfaceAlbedoType.F90:201.
  ! The 20-class table; mxsoil_color = 20 on modern surface datasets.
  integer, parameter, public :: mxsoil_color = 20
  real(r8), parameter :: albsat(mxsoil_color, numrad) = reshape( (/ &
       0.25_r8,0.23_r8,0.21_r8,0.20_r8,0.19_r8,0.18_r8,0.17_r8,0.16_r8, &
       0.15_r8,0.14_r8,0.13_r8,0.12_r8,0.11_r8,0.10_r8,0.09_r8,0.08_r8, &
       0.07_r8,0.06_r8,0.05_r8,0.04_r8, &
       0.50_r8,0.46_r8,0.42_r8,0.40_r8,0.38_r8,0.36_r8,0.34_r8,0.32_r8, &
       0.30_r8,0.28_r8,0.26_r8,0.24_r8,0.22_r8,0.20_r8,0.18_r8,0.16_r8, &
       0.14_r8,0.12_r8,0.10_r8,0.08_r8 /), (/ mxsoil_color, numrad /) )
  real(r8), parameter :: albdry(mxsoil_color, numrad) = reshape( (/ &
       0.36_r8,0.34_r8,0.32_r8,0.31_r8,0.30_r8,0.29_r8,0.28_r8,0.27_r8, &
       0.26_r8,0.25_r8,0.24_r8,0.23_r8,0.22_r8,0.20_r8,0.18_r8,0.16_r8, &
       0.14_r8,0.12_r8,0.10_r8,0.08_r8, &
       0.61_r8,0.57_r8,0.53_r8,0.51_r8,0.49_r8,0.48_r8,0.45_r8,0.43_r8, &
       0.41_r8,0.39_r8,0.37_r8,0.35_r8,0.33_r8,0.31_r8,0.29_r8,0.27_r8, &
       0.25_r8,0.23_r8,0.21_r8,0.16_r8 /), (/ mxsoil_color, numrad /) )

  ! Cold-start snow albedo, ELM SurfaceAlbedoType InitCold. Held constant
  ! because SNICAR is not ported; see the header.
  real(r8), parameter, public :: albsnd_const = 0.6_r8
  real(r8), parameter, public :: albsni_const = 0.6_r8

  public :: elmxx_surface_albedo_kernel

contains

  !-----------------------------------------------------------------------
  subroutine elmxx_surface_albedo_kernel(nc, np, npft, patch_col, coszen_col, &
       soil_color, h2osoi_vol_top, frac_sno, patch_ivt, elai_in, esai_in,     &
       t_veg, fwet, rhol, rhos, taul, taus, xl,                               &
       albsod, albsoi, albgrd, albgri,                                        &
       albd, albi, fabd, fabi, ftdd, ftid, ftii,                              &
       nrad, tlai_z, fsun_z, fabd_sun_z, fabi_sun_z, fabd_sha_z, fabi_sha_z, &
       vcmaxcintsun, vcmaxcintsha)
    !
    ! One full SurfaceAlbedo pass over nc non-urban columns and np patches.
    !
    ! Column arrays are indexed 1:nc, patch arrays 1:np, and patch_col(p) is
    ! the 1-based column index of patch p (<= 0 means "not on one of these
    ! columns", and the patch keeps its no-canopy defaults).
    !
    ! PFT parameter arrays are 0-based on their first index, as ELM's are:
    ! PFT 0 is bare ground, and patch_ivt is 0-based to match.
    !
    implicit none
    integer , intent(in)  :: nc, np, npft
    integer , intent(in)  :: patch_col(np)          ! patch -> column, 1-based
    real(r8), intent(in)  :: coszen_col(nc)         ! cos solar zenith angle
    integer , intent(in)  :: soil_color(nc)         ! ELM soil colour class
    real(r8), intent(in)  :: h2osoi_vol_top(nc)     ! layer-1 water [m3/m3]
    real(r8), intent(in)  :: frac_sno(nc)           ! snow cover fraction
    integer , intent(in)  :: patch_ivt(np)          ! PFT index, 0-based
    real(r8), intent(in)  :: elai_in(np), esai_in(np)
    real(r8), intent(in)  :: t_veg(np), fwet(np)
    real(r8), intent(in)  :: rhol(0:npft-1,numrad), rhos(0:npft-1,numrad)
    real(r8), intent(in)  :: taul(0:npft-1,numrad), taus(0:npft-1,numrad)
    real(r8), intent(in)  :: xl(0:npft-1)
    real(r8), intent(out) :: albsod(nc,numrad), albsoi(nc,numrad)
    real(r8), intent(out) :: albgrd(nc,numrad), albgri(nc,numrad)
    real(r8), intent(out) :: albd(np,numrad), albi(np,numrad)
    real(r8), intent(out) :: fabd(np,numrad), fabi(np,numrad)
    real(r8), intent(out) :: ftdd(np,numrad), ftid(np,numrad), ftii(np,numrad)
    integer , intent(out) :: nrad(np)
    real(r8), intent(out) :: tlai_z(np), fsun_z(np)
    real(r8), intent(out) :: fabd_sun_z(np), fabi_sun_z(np)
    real(r8), intent(out) :: fabd_sha_z(np), fabi_sha_z(np)
    ! Leaf-to-canopy scaling for vcmax. Nothing in ELMxx reads these yet --
    ! Photosynthesis is the consumer and is not ported. They are computed here
    ! because ELM computes them here, from state that exists only here, and
    ! because leaving a hole in a validated kernel is how the next port
    ! acquires an unseeded view that reads zero. See STATUS section J.
    real(r8), intent(out) :: vcmaxcintsun(np), vcmaxcintsha(np)
    !
    integer  :: c, p, ib, ivt, isc
    real(r8) :: cosz, inc, wl, ws, laisum
    real(r8) :: omegal, asu, betadl, betail, tmp0, tmp1, tmp2, tmp3, tmp4
    real(r8) :: tmp5, tmp6, tmp7, tmp8, tmp9, betad, betai
    real(r8) :: b, c1, d, f, h, sigma, p1, p2, p3, p4, t1, s1, s2
    real(r8) :: u1, u2, u3, d1, d2, h1, h2, h3, h4, h5, h6
    real(r8) :: h7, h8, h9, h10, a1, a2, om
    real(r8) :: phi1, phi2, chil, gdir, twostext, avmu, temp0, temp1, temp2v
    real(r8) :: rho_b, tau_b, elai_p, esai_p, fabd_sun, fabd_sha
    real(r8) :: fabi_sun, fabi_sha, fsunz, extkb
    real(r8), parameter :: extkn = 0.30_r8   ! leaf nitrogen decay, ELM
    !-----------------------------------------------------------------------

    !-----------------------------------------------------------------
    ! Soil and ground albedo.
    !
    ! ELM leaves every column albedo at ZERO where coszen <= 0 -- night is not
    ! a small albedo, it is "no solar calculation was done". SurfaceRadiation
    ! gates on the same test, so the zeros are never consumed.
    !-----------------------------------------------------------------
    albsod = 0.0_r8; albsoi = 0.0_r8
    albgrd = 0.0_r8; albgri = 0.0_r8

    do ib = 1, numrad
       do c = 1, nc
          if (coszen_col(c) <= 0.0_r8) cycle
          isc = soil_color(c)
          if (isc < 1 .or. isc > mxsoil_color) cycle   ! colour 0 = no soil albedo
          ! Wetter soil is darker: ELM's linear correction on layer-1 water.
          inc = max(0.11_r8 - 0.40_r8*h2osoi_vol_top(c), 0.0_r8)
          albsod(c,ib) = min(albsat(isc,ib) + inc, albdry(isc,ib))
          albsoi(c,ib) = albsod(c,ib)

          ! Weight soil against snow. With SNICAR unported the snow end is
          ! ELM's cold-start constant, which is only defensible while
          ! frac_sno is zero.
          albgrd(c,ib) = albsod(c,ib)*(1.0_r8 - frac_sno(c)) + albsnd_const*frac_sno(c)
          albgri(c,ib) = albsoi(c,ib)*(1.0_r8 - frac_sno(c)) + albsni_const*frac_sno(c)
       end do
    end do

    !-----------------------------------------------------------------
    ! Patch albedo: canopy two-stream where there is a canopy and sun,
    ! the bare ground's own albedo where there is sun but no canopy.
    !-----------------------------------------------------------------
    ! ELM's initialization, SurfaceAlbedoMod:279-317. NOT zero: an albedo of
    ! 1 with zero transmittance is ELM's marker for "no solar calculation was
    ! done here", and it is what the coupler is handed at night. Every patch
    ! either keeps it (night) or has it overwritten below (day).
    albd = 1.0_r8; albi = 1.0_r8; fabd = 0.0_r8; fabi = 0.0_r8
    ftdd = 0.0_r8; ftid = 0.0_r8; ftii = 0.0_r8
    fsun_z = 0.0_r8
    fabd_sun_z = 0.0_r8; fabi_sun_z = 0.0_r8
    fabd_sha_z = 0.0_r8; fabi_sha_z = 0.0_r8

    ! Canopy layering, ELM SurfaceAlbedoMod:814-826. Done for EVERY patch and
    ! every step, sunlit or not -- "because layering is needed for all time
    ! steps regardless of radiation", in ELM's own comment. Gating it on the
    ! sun leaves CanopySunShadeFractions with tlai_z = 0 all night, hence
    ! laisha = 0, hence no nighttime canopy conductance at all.
    nrad = 0; tlai_z = 0.0_r8
    vcmaxcintsun = 0.0_r8; vcmaxcintsha = 0.0_r8
    do p = 1, np
       if (patch_col(p) <= 0) cycle
       elai_p = elai_in(p); if (elai_p < 0.05_r8) elai_p = 0.0_r8
       nrad(p)   = 1                 ! nlevcan = 1
       tlai_z(p) = elai_p

       ! Default leaf-to-canopy scaling, ELM SurfaceAlbedoMod:946-960. The
       ! whole canopy is shaded when there is no sun, so the nitrogen profile
       ! exp(-extkn*x) integrated over the canopy goes to the shaded leaf and
       ! the sunlit coefficient is zero. TwoStream overwrites both where the
       ! sun is up; this is what a night, or a bare patch, keeps.
       if (elai_p > 0.0_r8) then
          vcmaxcintsha(p) = (1.0_r8 - exp(-extkn*elai_p)) / extkn / elai_p
       end if
    end do

    do p = 1, np
       c = patch_col(p)
       if (c <= 0) cycle
       if (coszen_col(c) <= 0.0_r8) cycle    ! night: ELM computes nothing

       elai_p = elai_in(p); if (elai_p < 0.05_r8) elai_p = 0.0_r8
       esai_p = esai_in(p); if (esai_p < 0.05_r8) esai_p = 0.0_r8

       ! ELM's filter_novegsol, SurfaceAlbedoMod:986-1004: a sunlit patch with
       ! no canopy is transparent and reflects the GROUND's albedo. Leaving it
       ! at the initialized 1, or at zero, is not a small error -- it is the
       ! whole reflected shortwave of every bare, crop-fallow and snow-buried
       ! patch, and on these twins that is fourteen patches of seventeen.
       if (elai_p + esai_p <= 0.0_r8) then
          do ib = 1, numrad
             albd(p,ib) = albgrd(c,ib)
             albi(p,ib) = albgri(c,ib)
             ftdd(p,ib) = 1.0_r8
             ftid(p,ib) = 0.0_r8
             ftii(p,ib) = 1.0_r8
          end do
          cycle
       end if

       ! ELM's filter_vegsol: vegetated AND sunlit.
       ivt  = patch_ivt(p)
       cosz = max(0.001_r8, coszen_col(c))

       ! Leaf angle distribution (Ross index), clipped as ELM clips it.
       chil = min(max(xl(ivt), -0.4_r8), 0.6_r8)
       if (abs(chil) <= 0.01_r8) chil = 0.01_r8
       phi1     = 0.5_r8 - 0.633_r8*chil - 0.330_r8*chil*chil
       phi2     = 0.877_r8 * (1.0_r8 - 2.0_r8*phi1)
       gdir     = phi1 + phi2*cosz
       twostext = gdir/cosz
       avmu     = (1.0_r8 - phi1/phi2 * log((phi1+phi2)/phi1)) / phi2
       temp0    = gdir + phi2*cosz
       temp1    = phi1*cosz
       temp2v   = 1.0_r8 - temp1/temp0 * log((temp1+temp0)/temp1)

       ! Leaf/stem weighted optics.
       wl = elai_p / max(elai_p + esai_p, mpe)
       ws = esai_p / max(elai_p + esai_p, mpe)

       do ib = 1, numrad
          rho_b = max(rhol(ivt,ib)*wl + rhos(ivt,ib)*ws, mpe)
          tau_b = max(taul(ivt,ib)*wl + taus(ivt,ib)*ws, mpe)

          omegal = rho_b + tau_b
          asu    = 0.5_r8*omegal*gdir/temp0 * temp2v
          betadl = (1.0_r8 + avmu*twostext)/(omegal*avmu*twostext) * asu
          betail = 0.5_r8 * ((rho_b+tau_b) + (rho_b-tau_b) &
                 * ((1.0_r8+chil)/2.0_r8)**2) / omegal

          ! Intercepted snow brightens and depolarises the canopy. ELM's
          ! switch is t_veg vs freezing, not a snow mass.
          if (t_veg(p) > tfrz) then
             om    = omegal
             betad = betadl
             betai = betail
          else
             om    = (1.0_r8-fwet(p))*omegal + fwet(p)*omegas(ib)
             betad = ((1.0_r8-fwet(p))*omegal*betadl + fwet(p)*omegas(ib)*betads) / om
             betai = ((1.0_r8-fwet(p))*omegal*betail + fwet(p)*omegas(ib)*betais) / om
          end if

          b     = 1.0_r8 - om + om*betai
          c1    = om*betai
          tmp0  = avmu*twostext
          d     = tmp0 * om*betad
          f     = tmp0 * om*(1.0_r8-betad)
          tmp1  = b*b - c1*c1
          h     = sqrt(tmp1) / avmu
          sigma = tmp0*tmp0 - tmp1
          p1    = b + avmu*h
          p2    = b - avmu*h
          p3    = b + tmp0
          p4    = b - tmp0

          t1 = min(h*(elai_p+esai_p), 40.0_r8)
          s1 = exp(-t1)
          t1 = min(twostext*(elai_p+esai_p), 40.0_r8)
          s2 = exp(-t1)

          ! ---- direct beam ----
          u1   = b - c1/max(albgrd(c,ib), mpe)
          u2   = b - c1*albgrd(c,ib)
          u3   = f + c1*albgrd(c,ib)
          tmp2 = u1 - avmu*h
          tmp3 = u1 + avmu*h
          d1   = p1*tmp2/s1 - p2*tmp3*s1
          tmp4 = u2 + avmu*h
          tmp5 = u2 - avmu*h
          d2   = tmp4/s1 - tmp5*s1
          h1   = -d*p4 - c1*f
          tmp6 = d - h1*p3/sigma
          tmp7 = (d - c1 - h1/sigma*(u1+tmp0)) * s2
          h2   = (tmp6*tmp2/s1 - p2*tmp7) / d1
          h3   = -(tmp6*tmp3*s1 - p1*tmp7) / d1
          h4   = -f*p3 - c1*d
          tmp8 = h4/sigma
          tmp9 = (u3 - tmp8*(u2-tmp0)) * s2
          h5   = -(tmp8*tmp4/s1 + tmp9) / d2
          h6   = (tmp8*tmp5*s1 + tmp9) / d2

          albd(p,ib) = h1/sigma + h2 + h3
          ftid(p,ib) = h4*s2/sigma + h5*s1 + h6/s1
          ftdd(p,ib) = s2
          fabd(p,ib) = 1.0_r8 - albd(p,ib) &
                     - (1.0_r8-albgrd(c,ib))*ftdd(p,ib) &
                     - (1.0_r8-albgri(c,ib))*ftid(p,ib)

          a1 = h1/sigma * (1.0_r8 - s2*s2) / (2.0_r8*twostext) &
             + h2       * (1.0_r8 - s2*s1) / (twostext + h) &
             + h3       * (1.0_r8 - s2/s1) / (twostext - h)
          a2 = h4/sigma * (1.0_r8 - s2*s2) / (2.0_r8*twostext) &
             + h5       * (1.0_r8 - s2*s1) / (twostext + h) &
             + h6       * (1.0_r8 - s2/s1) / (twostext - h)

          fabd_sun = (1.0_r8 - om) * (1.0_r8 - s2 + 1.0_r8/avmu * (a1 + a2))
          fabd_sha = fabd(p,ib) - fabd_sun

          ! ---- diffuse ----
          u1   = b - c1/max(albgri(c,ib), mpe)
          u2   = b - c1*albgri(c,ib)
          tmp2 = u1 - avmu*h
          tmp3 = u1 + avmu*h
          d1   = p1*tmp2/s1 - p2*tmp3*s1
          tmp4 = u2 + avmu*h
          tmp5 = u2 - avmu*h
          d2   = tmp4/s1 - tmp5*s1
          h7   = (c1*tmp2) / (d1*s1)
          h8   = (-c1*tmp3*s1) / d1
          h9   = tmp4 / (d2*s1)
          h10  = (-tmp5*s1) / d2

          albi(p,ib) = h7 + h8
          ftii(p,ib) = h9*s1 + h10/s1
          fabi(p,ib) = 1.0_r8 - albi(p,ib) - (1.0_r8-albgri(c,ib))*ftii(p,ib)

          a1 = h7 * (1.0_r8 - s2*s1) / (twostext + h) + h8  * (1.0_r8 - s2/s1) / (twostext - h)
          a2 = h9 * (1.0_r8 - s2*s1) / (twostext + h) + h10 * (1.0_r8 - s2/s1) / (twostext - h)

          fabi_sun = (1.0_r8 - om) / avmu * (a1 + a2)
          fabi_sha = fabi(p,ib) - fabi_sun

          ! Sunlit fraction and per-unit-LAI absorption, VIS band only --
          ! this is absorbed PAR, and only the visible band is PAR.
          if (ib == 1) then
             t1    = min(twostext*(elai_p+esai_p), 40.0_r8)
             fsunz = (1.0_r8 - s2) / max(t1, mpe)
             fsun_z(p) = fsunz
             laisum = elai_p + esai_p
             if (fsunz > 0.0_r8 .and. laisum > 0.0_r8) then
                fabd_sun_z(p) = fabd_sun / (fsunz*laisum)
                fabi_sun_z(p) = fabi_sun / (fsunz*laisum)
             end if
             if (fsunz < 1.0_r8 .and. laisum > 0.0_r8) then
                fabd_sha_z(p) = fabd_sha / ((1.0_r8 - fsunz)*laisum)
                fabi_sha_z(p) = fabi_sha / ((1.0_r8 - fsunz)*laisum)
             end if

             ! Leaf-to-canopy scaling with the sun up, ELM TwoStream:1518-1528.
             ! The sunlit leaves see both the nitrogen decay and the direct
             ! beam extinction; the shaded leaves get whatever is left.
             extkb = twostext
             vcmaxcintsun(p) = (1.0_r8 - exp(-(extkn+extkb)*elai_p)) / (extkn + extkb)
             vcmaxcintsha(p) = (1.0_r8 - exp(-extkn*elai_p)) / extkn - vcmaxcintsun(p)
             if (elai_p > 0.0_r8) then
                vcmaxcintsun(p) = vcmaxcintsun(p) / (fsunz*elai_p)
                vcmaxcintsha(p) = vcmaxcintsha(p) / ((1.0_r8 - fsunz)*elai_p)
             else
                vcmaxcintsun(p) = 0.0_r8
                vcmaxcintsha(p) = 0.0_r8
             end if
          end if
       end do
    end do

  end subroutine elmxx_surface_albedo_kernel

end module elmxxSurfaceAlbedoKernelMod
