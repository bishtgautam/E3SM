module elmxxRootMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Root distribution and the transpiration wetness factor.
  !
  ! TWO PIECES, DIFFERENT LIFETIMES:
  !   rootfr  vertical root profile per patch. Depends only on PFT and the
  !           soil grid, so it is built once at init.
  !   btran   transpiration wetness factor, 0 (fully stressed) to 1 (none).
  !           Depends on soil moisture and temperature, so it is recomputed
  !           every step. With the kernels off it happens not to change, but
  !           structurally it is per-step and is treated that way rather than
  !           seeded once and left to go stale when hydrology starts running.
  !
  ! WHY THIS IS FORTRAN-SIDE. ELMxx has no rootfr view and no setter for one;
  ! CanopyFluxes reads btran and explicitly does not compute it, and
  ! RootWaterUpdate takes rootr the same way. So the root-water-stress
  ! calculation belongs here, and only its results cross the boundary.
  !
  ! Provenance: rootfr is ELM RootBiophysMod init_vegrootfr, the Zeng (1998)
  ! two-exponential profile; btran and rootr are ELM SoilMoistStressMod
  ! calc_root_moist_stress_clm45default.
  !-----------------------------------------------------------------------

  use shr_kind_mod    , only : r8 => shr_kind_r8
  use shr_sys_mod     , only : shr_sys_abort, shr_sys_flush
  use shr_const_mod   , only : SHR_CONST_TKFRZ
  use elmxxSpmdMod    , only : masterproc, iam
  use elmxxSubgridMod , only : num_patches, patch_column, patch_itype
  use elmxxPftconMod  , only : pftcon_read, npft_param, roota_par, rootb_par, &
                               smpso, smpsc, tc_stress
  use elmxxSoilPropMod, only : soil_prop_built, nlevsoi, nlevgrnd, nlevsno, &
                               nlevbed, zisoi, watsat, bsw, sucsat, &
                               col_dz, col_t_soisno, col_h2osoi_liq, col_h2osoi_ice

  implicit none
  save
  private

  real(r8), public, pointer :: rootfr(:,:) => null()  ! (num_patches, nlevgrnd)
  real(r8), public, pointer :: rootr (:,:) => null()  ! (num_patches, nlevgrnd)
  real(r8), public, pointer :: btran (:)   => null()  ! (num_patches), 0..1

  logical, public :: root_built = .false.

  ! ELM CanopyFluxesMod: btran below this counts as no uptake at all.
  real(r8), parameter, public :: btran0 = 0.0_r8

  public :: elmxx_root_init
  public :: elmxx_compute_btran
  public :: elmxx_root_clean

contains

  !-----------------------------------------------------------------------
  subroutine elmxx_root_init(logunit)
    !
    ! Build rootfr once. ELM RootBiophysMod:
    !
    !   rootfr(p,j) = 0.5*( exp(-a*zi(j-1)) + exp(-b*zi(j-1))
    !                     - exp(-a*zi(j))   - exp(-b*zi(j)) )   j < nlevsoi
    !   rootfr(p,nlevsoi) = 0.5*( exp(-a*zi(nlevsoi-1)) + exp(-b*zi(nlevsoi-1)) )
    !
    ! The last layer takes the whole remaining tail rather than a difference,
    ! which is what makes the profile sum to one.
    !
    ! ELM's use_var_soil_thick renormalization is deliberately absent: varying
    ! depth to bedrock is out of scope, so nlevbed is nlevsoi and the branch
    ! cannot fire. Layers below nlevsoi are zero, as ELM sets them.
    !
    implicit none
    integer, intent(in) :: logunit
    integer  :: p, j, ivt
    real(r8) :: a, b, total, worst
    character(len=*), parameter :: subname = '(elmxx_root_init) '

    if (.not. pftcon_read)    call shr_sys_abort(subname//'ERROR: PFT parameters not read')
    if (.not. soil_prop_built) call shr_sys_abort(subname//'ERROR: soil grid not built')

    call elmxx_root_clean()
    allocate(rootfr(num_patches, nlevgrnd), rootr(num_patches, nlevgrnd), &
             btran(num_patches))
    rootfr = 0.0_r8; rootr = 0.0_r8; btran = 0.0_r8

    worst = 0.0_r8
    do p = 1, num_patches
       ivt = patch_itype(p)
       if (ivt < 0 .or. ivt > npft_param-1) then
          call shr_sys_abort(subname//'ERROR: PFT index outside the parameter file')
       end if

       ! Bare ground has no roots. ELM skips noveg explicitly; the parameter
       ! file also carries zeros there, so this is belt and braces -- but an
       ! exp(0) profile would otherwise put spurious roots on bare ground.
       if (ivt == 0) cycle

       a = roota_par(ivt)
       b = rootb_par(ivt)

       do j = 1, nlevsoi-1
          rootfr(p,j) = 0.5_r8 * ( exp(-a*zisoi(j-1)) + exp(-b*zisoi(j-1))   &
                                 - exp(-a*zisoi(j))   - exp(-b*zisoi(j)) )
       end do
       rootfr(p,nlevsoi) = 0.5_r8 * ( exp(-a*zisoi(nlevsoi-1)) &
                                    + exp(-b*zisoi(nlevsoi-1)) )

       total = sum(rootfr(p,1:nlevsoi))
       worst = max(worst, abs(total - 1.0_r8))
    end do

    root_built = .true.

    write(logunit,*) subname,'rank ',iam,' built rootfr for ',num_patches, &
                     ' patches over ',nlevsoi,' soil layers'
    write(logunit,*) '    worst |sum(rootfr)-1| over vegetated patches = ',worst
    call shr_sys_flush(logunit)

    ! The profile is a partition of unity by construction -- the last layer
    ! closes the tail. If it does not sum to one, either zisoi is wrong or the
    ! last-layer term was dropped, and every btran downstream is scaled wrong.
    if (worst > 1.0e-12_r8) then
       call shr_sys_abort(subname//'ERROR: root fractions do not sum to one')
    end if

  end subroutine elmxx_root_init

  !-----------------------------------------------------------------------
  subroutine elmxx_compute_btran(logunit, report)
    !
    ! ELM SoilMoistStressMod calc_root_moist_stress_clm45default.
    !
    ! Per layer, where there is liquid water and the soil is not too cold:
    !   eff_porosity = watsat - h2osoi_ice/(dz*denice)
    !   liqvol       = h2osoi_liq/(dz*denh2o), capped at eff_porosity
    !   s_node       = max(liqvol/eff_porosity, 0.01)
    !   smp_node     = max(smpsc, -sucsat*s_node**(-bsw))
    !   rresis       = min( (eff_porosity/watsat)*(smp_node-smpsc)/(smpso-smpsc), 1 )
    !   rootr        = rootfr*rresis,  btran = sum(max(rootr,0))
    ! then rootr is normalized by btran, so the layers partition the uptake.
    !
    ! The cold cutoff is tfrz + tc_stress with tc_stress negative, i.e. a
    ! couple of degrees BELOW freezing, not at it.
    !
    implicit none
    integer, intent(in) :: logunit
    logical, intent(in) :: report
    integer  :: p, c, j, ivt, m
    real(r8) :: eff_por, liqvol, s_node, smp_node, rresis, tcold
    real(r8) :: diag_s, diag_smp, diag_smpsc, diag_rresis
    logical  :: diag_taken
    real(r8), parameter :: denh2o = 1000.0_r8
    real(r8), parameter :: denice =  917.0_r8
    character(len=*), parameter :: subname = '(elmxx_compute_btran) '

    if (.not. root_built) call shr_sys_abort(subname//'ERROR: rootfr not built')

    tcold = SHR_CONST_TKFRZ + tc_stress
    rootr = 0.0_r8
    btran = 0.0_r8
    diag_s = 0.0_r8; diag_smp = 0.0_r8; diag_smpsc = 0.0_r8; diag_rresis = 0.0_r8
    diag_taken = .false.

    do p = 1, num_patches
       ivt = patch_itype(p)
       if (ivt == 0) cycle            ! bare ground transpires nothing
       c = patch_column(p)

       do j = 1, nlevbed
          m = j + nlevsno             ! packed slot for ELM layer j

          ! ELM floors eff_porosity at 0.01 (HydrologyNoDrainageMod, where the
          ! array the stress calc reads is actually set) rather than skipping
          ! the layer. Neither binds without ice, but the floor is what ELM does.
          eff_por = max(0.01_r8, watsat(c,j) - col_h2osoi_ice(c,m) / (col_dz(c,m)*denice))

          ! ELM does NOT cap the liquid volume at the effective porosity:
          !   h2osoi_liqvol(c,j) = h2osoi_liq(c,j)/(dz(c,j)*denh2o)
          ! Capping it here silently limits s_node to 1 and so understates the
          ! matric potential of a near-saturated layer.
          liqvol = col_h2osoi_liq(c,m) / (col_dz(c,m)*denh2o)

          if (liqvol <= 0.0_r8 .or. col_t_soisno(c,m) <= tcold) cycle

          s_node   = max(liqvol/eff_por, 0.01_r8)
          smp_node = max(smpsc(ivt), -sucsat(c,j) * s_node**(-bsw(c,j)))

          rresis = min( (eff_por/watsat(c,j)) * (smp_node - smpsc(ivt)) &
                        / (smpso(ivt) - smpsc(ivt)), 1.0_r8 )

          rootr(p,j) = rootfr(p,j) * rresis
          btran(p)   = btran(p) + max(rootr(p,j), 0.0_r8)

          if (.not. diag_taken .and. j == 1) then
             diag_s = s_node; diag_smp = smp_node
             diag_smpsc = smpsc(ivt); diag_rresis = rresis
             diag_taken = .true.
          end if
       end do

       ! Normalize so the layers partition the uptake rather than scale it.
       if (btran(p) > btran0) then
          rootr(p,1:nlevgrnd) = rootr(p,1:nlevgrnd) / btran(p)
       else
          rootr(p,1:nlevgrnd) = 0.0_r8
       end if
    end do

    if (report) then
       write(logunit,*) subname,'rank ',iam,' btran over ',num_patches,' patches: ', &
                        minval(btran),' .. ',maxval(btran)
       ! When btran is zero everywhere the useful question is WHY -- a dry
       ! soil and a broken formula look identical from btran alone. These are
       ! the two quantities that decide it: the top-layer matric potential and
       ! the PFT's closure threshold. If smp has been clamped up to smpsc, the
       ! soil is simply drier than the plant can extract from, and zero is the
       ! right answer rather than a bug.
       write(logunit,*) '    diag: top-layer s_node ',diag_s,' smp_node ',diag_smp, &
                        ' smpsc ',diag_smpsc,' rresis ',diag_rresis
       call shr_sys_flush(logunit)
       if (minval(btran) < 0.0_r8 .or. maxval(btran) > 1.0_r8) then
          write(logunit,*) subname,'SUSPECT: btran outside [0,1]'
       end if
    end if

  end subroutine elmxx_compute_btran

  !-----------------------------------------------------------------------
  subroutine elmxx_root_clean()
    implicit none
    if (associated(rootfr)) deallocate(rootfr)
    if (associated(rootr))  deallocate(rootr)
    if (associated(btran))  deallocate(btran)
    rootfr => null(); rootr => null(); btran => null()
    root_built = .false.
  end subroutine elmxx_root_clean

end module elmxxRootMod
