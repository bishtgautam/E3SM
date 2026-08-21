module elmxxSurfaceAlbedoMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Surface albedo: gather ELMxx state, run the albedo kernel, push results.
  !
  ! A FORTRAN PORT, NOT A CROSSING. ELMxx has no albedo kernel -- there is no
  ! ELMxxComputeSurfaceAlbedo anywhere in the C API. Every one of the ~20
  ! albedo fields the radiation kernels read is caller-supplied through a
  ! setter. So this is Fortran-side work, like btran and the ground heat flux.
  !
  ! THE PHYSICS IS NOT HERE. It is in elmxxSurfaceAlbedoKernelMod, which
  ! depends on nothing but its arguments so that it can be replayed offline
  ! against an ELM restart -- see tools/validate_surface_albedo.py. This
  ! module is the part that cannot be replayed: the subgrid maps, the ELMxx
  ! object, the setters.
  !
  ! CALLED AT THE END OF THE TIMESTEP, as ELM does.
  ! elm_driver.F90 calls SurfaceAlbedo near the end of the step, gated on
  ! doalb, so the albedos SurfaceRadiation reads at step N were computed at
  ! step N-1. Step one reads SurfaceAlbedoType InitCold's constants, which
  ! elmxx_kokkos_seed_albedo supplies. Keeping that placement matters: moving
  ! it to the top of the step would change which state the two-stream sees
  ! (t_veg, fwet and h2osoi_vol have all been updated by then) and quietly
  ! diverge from ELM.
  !-----------------------------------------------------------------------

  use shr_kind_mod        , only : r8 => shr_kind_r8
  use shr_sys_mod         , only : shr_sys_abort, shr_sys_flush
  use shr_const_mod       , only : SHR_CONST_PI
  use elmxxSpmdMod        , only : masterproc, iam
  use elmxxSubgridMod     , only : col_landunit, lun_gridcell, patch_column, &
                                   patch_itype, num_columns, num_patches
  use elmxxSurfaceStateMod, only : col_soil_color, patch_lai, patch_sai
  use elmxxSoilPropMod    , only : nlevsno, nlevgrnd, col_h2osoi_vol
  use elmxxPftconMod      , only : rhol, rhos, taul, taus, xl, pftcon_read, &
                                   npft_param
  use elmxxKokkosStateMod , only : n_kokkos_col, n_kokkos_patch, &
                                   col_of_kcol, kcol_of_col, &
                                   patch_of_kpatch, kokkos_state_built
  use elmxxSurfaceAlbedoKernelMod, only : numrad, elmxx_surface_albedo_kernel
  use elmxx_mod           , only : ELMxxType, ELMXX_SUCCESS, &
                                   ELMxxSetAlbgrd, ELMxxSetAlbgri, &
                                   ELMxxSetAlbsod, ELMxxSetAlbsoi, &
                                   ELMxxSetAlbd, ELMxxSetAlbi, &
                                   ELMxxSetFabd, ELMxxSetFabi, &
                                   ELMxxSetFtdd, ELMxxSetFtid, ELMxxSetFtii, &
                                   ELMxxSetFsunZ, ELMxxSetTlaiZ, ELMxxSetNrad, &
                                   ELMxxSetFabdSunZ, ELMxxSetFabiSunZ, &
                                   ELMxxSetFabdShaZ, ELMxxSetFabiShaZ, &
                                   ELMxxGetTVeg, ELMxxGetFwet, ELMxxGetFracSno

  implicit none
  save
  private

  logical, public :: surface_albedo_built = .false.

  ! Leaf-to-canopy scaling coefficients for vcmax, packed on the Kokkos patch
  ! index. ELM's SurfaceAlbedo is where these are computed and Photosynthesis
  ! is what reads them, so they wait here until that port lands -- there is no
  ! ELMxx setter for them yet, because the C API has no Photosynthesis.
  real(r8), public, allocatable :: patch_vcmaxcintsun(:)
  real(r8), public, allocatable :: patch_vcmaxcintsha(:)

  public :: elmxx_surface_albedo_report

  real(r8), allocatable :: last_albd(:,:), last_albgrd(:,:), last_fsun(:)
  real(r8), allocatable :: last_coszen(:)

contains

  !-----------------------------------------------------------------------


  !-----------------------------------------------------------------------
  subroutine elmxx_surface_albedo_report(logunit)
    !
    ! Graded by the bounds that must hold whatever the sun is doing: an albedo
    ! is a fraction, and a sunlit fraction is a fraction. Both abort rather
    ! than warn -- a two-stream that returns an albedo outside [0,1] has a
    ! broken solve, not a marginal input.
    !
    ! Bounds are all this can say. What the port is actually graded by is
    ! tools/validate_surface_albedo.py, which replays ELM's own inputs
    ! through the same kernel and compares against ELM's own outputs.
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
