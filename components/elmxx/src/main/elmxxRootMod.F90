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
  use elmxxRootKernelMod, only : elmxx_root_stress_kernel
  use elmxxSubgridMod , only : num_patches, num_columns, patch_column, patch_itype
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
    integer  :: p, c, j
    real(r8) :: diag_s, diag_smp, diag_smpsc, diag_rresis
    real(r8), allocatable :: k_liq(:,:), k_ice(:,:), k_dz(:,:), k_tsoi(:,:)
    real(r8), allocatable :: k_rresis(:,:)
    integer , allocatable :: k_pcol(:)
    real(r8), parameter :: denh2o = 1000.0_r8
    real(r8), parameter :: denice =  917.0_r8
    character(len=*), parameter :: subname = '(elmxx_compute_btran) '

    if (.not. root_built) call shr_sys_abort(subname//'ERROR: rootfr not built')

    ! Gather into the shapes the kernel takes: ELM soil layers 1..nlevgrnd,
    ! no snow slots. The packed column arrays carry snow in slots 1..nlevsno.
    allocate(k_liq(num_columns, nlevgrnd), k_ice(num_columns, nlevgrnd), &
             k_dz(num_columns, nlevgrnd),  k_tsoi(num_columns, nlevgrnd), &
             k_rresis(num_patches, nlevgrnd), k_pcol(num_patches))
    do c = 1, num_columns
       do j = 1, nlevgrnd
          k_liq (c,j) = col_h2osoi_liq(c, j + nlevsno)
          k_ice (c,j) = col_h2osoi_ice(c, j + nlevsno)
          k_dz  (c,j) = col_dz        (c, j + nlevsno)
          k_tsoi(c,j) = col_t_soisno  (c, j + nlevsno)
       end do
    end do
    do p = 1, num_patches
       k_pcol(p) = patch_column(p)
    end do

    call elmxx_root_stress_kernel(num_patches, num_columns, nlevbed, nlevgrnd, &
         patch_itype, k_pcol, rootfr, k_liq, k_ice, k_dz, k_tsoi,              &
         watsat, bsw, sucsat, smpsc, smpso, tc_stress, btran0,                 &
         denice, denh2o, rootr, btran, k_rresis)

    diag_s = 0.0_r8; diag_smp = 0.0_r8; diag_smpsc = 0.0_r8
    diag_rresis = 0.0_r8
    do p = 1, num_patches
       if (patch_itype(p) /= 0) then
          diag_rresis = k_rresis(p,1)
          diag_smpsc  = smpsc(patch_itype(p))
          exit
       end if
    end do

    deallocate(k_liq, k_ice, k_dz, k_tsoi, k_rresis, k_pcol)

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
