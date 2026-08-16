module elmxxSurfaceAlbedoMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Surface albedo: soil albedo, ground albedo, and the Sellers two-stream
  ! canopy radiative transfer.
  !
  ! A FORTRAN PORT, NOT A CROSSING. ELMxx has no albedo kernel -- there is no
  ! ELMxxComputeSurfaceAlbedo anywhere in the C API. Every one of the ~20
  ! albedo fields the radiation kernels read is caller-supplied through a
  ! setter. So this is Fortran-side work, like btran and the ground heat flux.
  !
  ! CALLED AT THE END OF THE TIMESTEP, as ELM does.
  ! elm_driver.F90 calls SurfaceAlbedo near the end of the step, gated on
  ! doalb, so the albedos SurfaceRadiation reads at step N were computed at
  ! step N-1. Step one reads SurfaceAlbedoType InitCold's constants, which
  ! elmxx_kokkos_seed_albedo supplies. Keeping that placement matters: moving
  ! it to the top of the step would change which state the two-stream sees
  ! (t_veg, fwet and h2osoi_vol have all been updated by then) and quietly
  ! diverge from ELM.
  !
  ! Ported from ELM SurfaceAlbedoMod.F90:
  !   SoilAlbedo   (:1015)  soil colour + surface wetness
  !   TwoStream    (:1147)  Sellers/Dickinson two-stream, nlevcan = 1 branch
  ! with the albsat/albdry colour tables from SurfaceAlbedoType.F90:201.
  !
  ! DELIBERATELY NOT PORTED, each with its consequence stated:
  !
  !   SNICAR         the full snow radiative transfer. Snow albedo falls back
  !                  to ELM's cold-start constants, so a snow-covered column
  !                  gets the wrong albedo. Free on the brazil twins, which
  !                  never accumulate snow; NOT free at f19.
  !   Albedo_TOP_Adjustment  sub-grid topographic correction. ELM gates it on
  !                  use_top_solar_rad, which is .false. here.
  !   lake / glacier / urban branches of SoilAlbedo -- out of scope, and the
  !                  packed views carry no lake or glacier column anyway.
  !-----------------------------------------------------------------------

  use shr_kind_mod        , only : r8 => shr_kind_r8
  use shr_sys_mod         , only : shr_sys_abort, shr_sys_flush
  use shr_orb_mod         , only : shr_orb_cosz
  use shr_const_mod       , only : SHR_CONST_PI
  use elmxxSpmdMod        , only : masterproc, iam
  use elmxxSubgridMod     , only : col_landunit, lun_gridcell, patch_column, &
                                   patch_itype, num_columns, num_patches
  use elmxxSurfaceStateMod, only : col_soil_color, patch_lai, patch_sai
  use elmxxSoilPropMod    , only : nlevsno, nlevgrnd, col_h2osoi_vol
  use elmxxPftconMod      , only : rhol, rhos, taul, taus, xl, pftcon_read
  use elmxxKokkosStateMod , only : n_kokkos_col, n_kokkos_patch, &
                                   col_of_kcol, kcol_of_col, &
                                   patch_of_kpatch, kokkos_state_built
  use elmxx_mod           , only : ELMxxType, ELMXX_SUCCESS, &
                                   ELMxxSetAlbgrd, ELMxxSetAlbgri, &
                                   ELMxxSetAlbsod, ELMxxSetAlbsoi, &
                                   ELMxxSetAlbd, ELMxxSetAlbi, &
                                   ELMxxSetFabd, ELMxxSetFabi, &
                                   ELMxxSetFtdd, ELMxxSetFtid, ELMxxSetFtii, &
                                   ELMxxSetFsunZ, ELMxxSetTlaiZ, ELMxxSetNrad, &
                                   ELMxxSetFabdSunZ, ELMxxSetFabiSunZ, &
                                   ELMxxSetFabdShaZ, ELMxxSetFabiShaZ, &
                                   ELMxxSetFsun, &
                                   ELMxxGetTVeg, ELMxxGetFwet, ELMxxGetFracSno

  implicit none
  save
  private

  integer, parameter :: numrad = 2

  ! Snow two-stream parameters, ELM elm_varcon.
  real(r8), parameter :: omegas(numrad) = (/ 0.8_r8, 0.4_r8 /)
  real(r8), parameter :: betads = 0.5_r8
  real(r8), parameter :: betais = 0.5_r8
  real(r8), parameter :: tfrz   = 273.15_r8
  real(r8), parameter :: mpe    = 1.0e-6_r8

  ! Soil albedo by colour class and band, ELM SurfaceAlbedoType.F90:201.
  ! The 20-class table; mxsoil_color = 20 on modern surface datasets.
  integer, parameter :: mxsoil_color = 20
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
  real(r8), parameter :: albsnd_const = 0.6_r8
  real(r8), parameter :: albsni_const = 0.6_r8

  logical, public :: surface_albedo_built = .false.

  public :: elmxx_surface_albedo
  public :: elmxx_surface_albedo_report

  real(r8), allocatable :: last_albd(:,:), last_albgrd(:,:), last_fsun(:)
  real(r8), allocatable :: last_coszen(:)

contains

  !-----------------------------------------------------------------------
  subroutine elmxx_surface_albedo(elm, nextsw_cday, declin, lat, lon, logunit)
    !
    ! One full SurfaceAlbedo pass, end of timestep.
    !
    implicit none
    type(ELMxxType), intent(in) :: elm
    real(r8), intent(in) :: nextsw_cday          ! calendar day of next radiation step
    real(r8), intent(in) :: declin               ! solar declination, radians
    real(r8), intent(in) :: lat(:), lon(:)       ! gridcell centres, DEGREES
    integer , intent(in) :: logunit

    integer  :: kc, kp, c, p, g, ib, ierr, ivt, sz(2)
    real(r8) :: cosz, inc, wl, ws, laisum
    real(r8) :: omegal, asu, betadl, betail, tmp0, tmp1, tmp2, tmp3, tmp4
    real(r8) :: tmp5, tmp6, tmp7, tmp8, tmp9, betad, betai
    real(r8) :: b, c1, d, f, h, sigma, p1, p2, p3, p4, t1, s1, s2
    real(r8) :: u1, u2, u3, d1, d2, h1, h2, h3, h4, h5, h6
    real(r8) :: h7, h8, h9, h10, a1, a2, om
    real(r8) :: phi1, phi2, chil, gdir, twostext, avmu, temp0, temp1, temp2v
    real(r8) :: rho_b, tau_b, elai_p, esai_p, fabd_sun, fabd_sha
    real(r8) :: fabi_sun, fabi_sha, fsunz
    integer  :: isc

    real(r8), allocatable :: coszen_col(:), albsod(:,:), albsoi(:,:)
    real(r8), allocatable :: albgrd(:,:), albgri(:,:), fracsno(:)
    real(r8), allocatable :: albd(:,:), albi(:,:), fabd(:,:), fabi(:,:)
    real(r8), allocatable :: ftdd(:,:), ftid(:,:), ftii(:,:)
    real(r8), allocatable :: tveg(:), fwet(:)
    real(r8), allocatable :: fsun_z(:), tlai_z(:), fabd_sun_z(:), fabi_sun_z(:)
    real(r8), allocatable :: fabd_sha_z(:), fabi_sha_z(:)
    integer , allocatable :: nrad(:)
    character(len=*), parameter :: subname = '(elmxx_surface_albedo) '

    if (.not. kokkos_state_built) call shr_sys_abort(subname//'ERROR: maps not built')
    if (.not. pftcon_read)        call shr_sys_abort(subname//'ERROR: PFT parameters not read')

    allocate(coszen_col(n_kokkos_col), fracsno(n_kokkos_col), &
             albsod(n_kokkos_col,numrad), albsoi(n_kokkos_col,numrad), &
             albgrd(n_kokkos_col,numrad), albgri(n_kokkos_col,numrad))
    allocate(tveg(n_kokkos_patch), fwet(n_kokkos_patch), &
             albd(n_kokkos_patch,numrad), albi(n_kokkos_patch,numrad), &
             fabd(n_kokkos_patch,numrad), fabi(n_kokkos_patch,numrad), &
             ftdd(n_kokkos_patch,numrad), ftid(n_kokkos_patch,numrad), &
             ftii(n_kokkos_patch,numrad), &
             fsun_z(n_kokkos_patch), tlai_z(n_kokkos_patch), &
             fabd_sun_z(n_kokkos_patch), fabi_sun_z(n_kokkos_patch), &
             fabd_sha_z(n_kokkos_patch), fabi_sha_z(n_kokkos_patch), &
             nrad(n_kokkos_patch))

    call ELMxxGetTVeg(elm, tveg, n_kokkos_patch, ierr);      call check(ierr, subname, 'TVeg')
    call ELMxxGetFwet(elm, fwet, n_kokkos_patch, ierr);      call check(ierr, subname, 'Fwet')
    call ELMxxGetFracSno(elm, fracsno, n_kokkos_col, ierr);  call check(ierr, subname, 'FracSno')

    !-----------------------------------------------------------------
    ! Solar zenith angle, per gridcell then broadcast to columns.
    !
    ! nextsw_cday and declin come from the coupler, exactly as ELM's
    ! lnd_comp_mct hands them to elm_drv. Using the model's own clock instead
    ! would drift against the atmosphere's radiation step.
    !-----------------------------------------------------------------
    do kc = 1, n_kokkos_col
       g = lun_gridcell(col_landunit(col_of_kcol(kc)))
       coszen_col(kc) = shr_orb_cosz(nextsw_cday, lat(g)*SHR_CONST_PI/180.0_r8, &
                                     lon(g)*SHR_CONST_PI/180.0_r8, declin)
    end do

    !-----------------------------------------------------------------
    ! Soil and ground albedo.
    !
    ! ELM leaves every albedo at ZERO where coszen <= 0 -- night is not a
    ! small albedo, it is "no solar calculation was done". SurfaceRadiation
    ! gates on the same test, so the zeros are never consumed.
    !-----------------------------------------------------------------
    albsod = 0.0_r8; albsoi = 0.0_r8
    albgrd = 0.0_r8; albgri = 0.0_r8

    do ib = 1, numrad
       do kc = 1, n_kokkos_col
          if (coszen_col(kc) <= 0.0_r8) cycle
          c   = col_of_kcol(kc)
          isc = col_soil_color(c)
          if (isc < 1 .or. isc > mxsoil_color) cycle   ! colour 0 = no soil albedo
          ! Wetter soil is darker: ELM's linear correction on layer-1 water.
          inc = max(0.11_r8 - 0.40_r8*col_h2osoi_vol(c,1), 0.0_r8)
          albsod(kc,ib) = min(albsat(isc,ib) + inc, albdry(isc,ib))
          albsoi(kc,ib) = albsod(kc,ib)

          ! Weight soil against snow. With SNICAR unported the snow end is
          ! ELM's cold-start constant, which is only defensible while
          ! frac_sno is zero -- true on the brazil twins.
          albgrd(kc,ib) = albsod(kc,ib)*(1.0_r8 - fracsno(kc)) + albsnd_const*fracsno(kc)
          albgri(kc,ib) = albsoi(kc,ib)*(1.0_r8 - fracsno(kc)) + albsni_const*fracsno(kc)
       end do
    end do

    !-----------------------------------------------------------------
    ! Two-stream canopy radiative transfer, per vegetated sunlit patch.
    !-----------------------------------------------------------------
    albd = 0.0_r8; albi = 0.0_r8; fabd = 0.0_r8; fabi = 0.0_r8
    ftdd = 1.0_r8; ftid = 0.0_r8; ftii = 1.0_r8
    fsun_z = 0.0_r8; tlai_z = 0.0_r8; nrad = 0
    fabd_sun_z = 0.0_r8; fabi_sun_z = 0.0_r8
    fabd_sha_z = 0.0_r8; fabi_sha_z = 0.0_r8

    do kp = 1, n_kokkos_patch
       p  = patch_of_kpatch(kp)
       c  = patch_column(p)
       kc = kcol_of_col(c) + 1
       if (kc <= 0) cycle

       elai_p = patch_lai(p); if (elai_p < 0.05_r8) elai_p = 0.0_r8
       esai_p = patch_sai(p); if (esai_p < 0.05_r8) esai_p = 0.0_r8

       ! ELM's filter_vegsol: vegetated AND sunlit. Everything else keeps the
       ! no-canopy defaults set above (fully transmitting, nothing absorbed).
       if (coszen_col(kc) <= 0.0_r8) cycle
       if (elai_p + esai_p <= 0.0_r8) cycle

       nrad(kp)   = 1                 ! nlevcan = 1
       tlai_z(kp) = elai_p

       ivt  = patch_itype(p)
       cosz = max(0.001_r8, coszen_col(kc))

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
          if (tveg(kp) > tfrz) then
             om    = omegal
             betad = betadl
             betai = betail
          else
             om    = (1.0_r8-fwet(kp))*omegal + fwet(kp)*omegas(ib)
             betad = ((1.0_r8-fwet(kp))*omegal*betadl + fwet(kp)*omegas(ib)*betads) / om
             betai = ((1.0_r8-fwet(kp))*omegal*betail + fwet(kp)*omegas(ib)*betais) / om
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
          u1   = b - c1/max(albgrd(kc,ib), mpe)
          u2   = b - c1*albgrd(kc,ib)
          u3   = f + c1*albgrd(kc,ib)
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

          albd(kp,ib) = h1/sigma + h2 + h3
          ftid(kp,ib) = h4*s2/sigma + h5*s1 + h6/s1
          ftdd(kp,ib) = s2
          fabd(kp,ib) = 1.0_r8 - albd(kp,ib) &
                      - (1.0_r8-albgrd(kc,ib))*ftdd(kp,ib) &
                      - (1.0_r8-albgri(kc,ib))*ftid(kp,ib)

          a1 = h1/sigma * (1.0_r8 - s2*s2) / (2.0_r8*twostext) &
             + h2       * (1.0_r8 - s2*s1) / (twostext + h) &
             + h3       * (1.0_r8 - s2/s1) / (twostext - h)
          a2 = h4/sigma * (1.0_r8 - s2*s2) / (2.0_r8*twostext) &
             + h5       * (1.0_r8 - s2*s1) / (twostext + h) &
             + h6       * (1.0_r8 - s2/s1) / (twostext - h)

          fabd_sun = (1.0_r8 - om) * (1.0_r8 - s2 + 1.0_r8/avmu * (a1 + a2))
          fabd_sha = fabd(kp,ib) - fabd_sun

          ! ---- diffuse ----
          u1   = b - c1/max(albgri(kc,ib), mpe)
          u2   = b - c1*albgri(kc,ib)
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

          albi(kp,ib) = h7 + h8
          ftii(kp,ib) = h9*s1 + h10/s1
          fabi(kp,ib) = 1.0_r8 - albi(kp,ib) - (1.0_r8-albgri(kc,ib))*ftii(kp,ib)

          a1 = h7 * (1.0_r8 - s2*s1) / (twostext + h) + h8  * (1.0_r8 - s2/s1) / (twostext - h)
          a2 = h9 * (1.0_r8 - s2*s1) / (twostext + h) + h10 * (1.0_r8 - s2/s1) / (twostext - h)

          fabi_sun = (1.0_r8 - om) / avmu * (a1 + a2)
          fabi_sha = fabi(kp,ib) - fabi_sun

          ! Sunlit fraction and per-unit-LAI absorption, VIS band only --
          ! this is absorbed PAR, and only the visible band is PAR.
          if (ib == 1) then
             t1    = min(twostext*(elai_p+esai_p), 40.0_r8)
             fsunz = (1.0_r8 - s2) / max(t1, mpe)
             fsun_z(kp) = fsunz
             laisum = elai_p + esai_p
             if (fsunz > 0.0_r8 .and. laisum > 0.0_r8) then
                fabd_sun_z(kp) = fabd_sun / (fsunz*laisum)
                fabi_sun_z(kp) = fabi_sun / (fsunz*laisum)
             end if
             if (fsunz < 1.0_r8 .and. laisum > 0.0_r8) then
                fabd_sha_z(kp) = fabd_sha / ((1.0_r8 - fsunz)*laisum)
                fabi_sha_z(kp) = fabi_sha / ((1.0_r8 - fsunz)*laisum)
             end if
          end if
       end do
    end do

    !-----------------------------------------------------------------
    ! Push. Column fields first, then patch.
    !-----------------------------------------------------------------
    sz = (/ n_kokkos_col, numrad /)
    call ELMxxSetAlbsod(elm, albsod, sz, ierr); call check(ierr, subname, 'Albsod')
    call ELMxxSetAlbsoi(elm, albsoi, sz, ierr); call check(ierr, subname, 'Albsoi')
    call ELMxxSetAlbgrd(elm, albgrd, sz, ierr); call check(ierr, subname, 'Albgrd')
    call ELMxxSetAlbgri(elm, albgri, sz, ierr); call check(ierr, subname, 'Albgri')

    sz = (/ n_kokkos_patch, numrad /)
    call ELMxxSetAlbd(elm, albd, sz, ierr); call check(ierr, subname, 'Albd')
    call ELMxxSetAlbi(elm, albi, sz, ierr); call check(ierr, subname, 'Albi')
    call ELMxxSetFabd(elm, fabd, sz, ierr); call check(ierr, subname, 'Fabd')
    call ELMxxSetFabi(elm, fabi, sz, ierr); call check(ierr, subname, 'Fabi')
    call ELMxxSetFtdd(elm, ftdd, sz, ierr); call check(ierr, subname, 'Ftdd')
    call ELMxxSetFtid(elm, ftid, sz, ierr); call check(ierr, subname, 'Ftid')
    call ELMxxSetFtii(elm, ftii, sz, ierr); call check(ierr, subname, 'Ftii')

    call ELMxxSetNrad(elm, nrad, n_kokkos_patch, ierr);        call check(ierr, subname, 'Nrad')
    call ELMxxSetTlaiZ(elm, tlai_z, n_kokkos_patch, ierr);     call check(ierr, subname, 'TlaiZ')
    call ELMxxSetFsunZ(elm, fsun_z, n_kokkos_patch, ierr);     call check(ierr, subname, 'FsunZ')
    call ELMxxSetFabdSunZ(elm, fabd_sun_z, n_kokkos_patch, ierr); call check(ierr, subname, 'FabdSunZ')
    call ELMxxSetFabiSunZ(elm, fabi_sun_z, n_kokkos_patch, ierr); call check(ierr, subname, 'FabiSunZ')
    call ELMxxSetFabdShaZ(elm, fabd_sha_z, n_kokkos_patch, ierr); call check(ierr, subname, 'FabdShaZ')
    call ELMxxSetFabiShaZ(elm, fabi_sha_z, n_kokkos_patch, ierr); call check(ierr, subname, 'FabiShaZ')

    ! Keep a copy for the report; nothing downstream reads these.
    if (.not. allocated(last_albd)) then
       allocate(last_albd(n_kokkos_patch,numrad), last_albgrd(n_kokkos_col,numrad), &
                last_fsun(n_kokkos_patch), last_coszen(n_kokkos_col))
    end if
    last_albd = albd; last_albgrd = albgrd
    last_fsun = fsun_z; last_coszen = coszen_col
    surface_albedo_built = .true.

    deallocate(coszen_col, fracsno, albsod, albsoi, albgrd, albgri)
    deallocate(tveg, fwet, albd, albi, fabd, fabi, ftdd, ftid, ftii, &
               fsun_z, tlai_z, fabd_sun_z, fabi_sun_z, fabd_sha_z, fabi_sha_z, nrad)

  end subroutine elmxx_surface_albedo

  !-----------------------------------------------------------------------
  subroutine elmxx_surface_albedo_report(logunit)
    !
    ! Graded by the bounds that must hold whatever the sun is doing: an albedo
    ! is a fraction, and a sunlit fraction is a fraction. Both abort rather
    ! than warn -- a two-stream that returns an albedo outside [0,1] has a
    ! broken solve, not a marginal input.
    !
    implicit none
    integer, intent(in) :: logunit
    character(len=*), parameter :: subname = '(elmxx_surface_albedo_report) '

    if (.not. surface_albedo_built) return

    write(logunit,*) subname,'rank ',iam,' surface albedo:'
    write(logunit,*) '    coszen      [-] ',minval(last_coszen),' .. ',maxval(last_coszen)
    write(logunit,*) '    albgrd vis  [-] ',minval(last_albgrd(:,1)),' .. ',maxval(last_albgrd(:,1))
    write(logunit,*) '    albgrd nir  [-] ',minval(last_albgrd(:,2)),' .. ',maxval(last_albgrd(:,2))
    write(logunit,*) '    albd   vis  [-] ',minval(last_albd(:,1)),' .. ',maxval(last_albd(:,1))
    write(logunit,*) '    albd   nir  [-] ',minval(last_albd(:,2)),' .. ',maxval(last_albd(:,2))
    write(logunit,*) '    fsun_z      [-] ',minval(last_fsun),' .. ',maxval(last_fsun)
    call shr_sys_flush(logunit)

    if (minval(last_albd) < 0.0_r8 .or. maxval(last_albd) > 1.0_r8) then
       call shr_sys_abort(subname//'ERROR: canopy albedo outside [0,1]')
    end if
    if (minval(last_albgrd) < 0.0_r8 .or. maxval(last_albgrd) > 1.0_r8) then
       call shr_sys_abort(subname//'ERROR: ground albedo outside [0,1]')
    end if
    if (minval(last_fsun) < 0.0_r8 .or. maxval(last_fsun) > 1.0_r8) then
       call shr_sys_abort(subname//'ERROR: sunlit fraction outside [0,1]')
    end if

  end subroutine elmxx_surface_albedo_report

  !-----------------------------------------------------------------------
  subroutine check(ierr, subname, what)
    implicit none
    integer, intent(in) :: ierr
    character(len=*), intent(in) :: subname, what
    if (ierr /= ELMXX_SUCCESS) then
       call shr_sys_abort(subname//'ERROR: '//trim(what)//' returned a non-success status')
    end if
  end subroutine check

end module elmxxSurfaceAlbedoMod
