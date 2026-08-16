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
  use shr_orb_mod         , only : shr_orb_cosz
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

    integer  :: kc, kp, c, p, g, ierr, sz(2)

    real(r8), allocatable :: coszen_col(:), albsod(:,:), albsoi(:,:)
    real(r8), allocatable :: albgrd(:,:), albgri(:,:), fracsno(:)
    real(r8), allocatable :: albd(:,:), albi(:,:), fabd(:,:), fabi(:,:)
    real(r8), allocatable :: ftdd(:,:), ftid(:,:), ftii(:,:)
    real(r8), allocatable :: tveg(:), fwet(:)
    real(r8), allocatable :: fsun_z(:), tlai_z(:), fabd_sun_z(:), fabi_sun_z(:)
    real(r8), allocatable :: fabd_sha_z(:), fabi_sha_z(:)
    real(r8), allocatable :: h2osoi_top(:), elai(:), esai(:)
    integer , allocatable :: nrad(:), soil_color(:), patch_col(:), patch_ivt(:)
    character(len=*), parameter :: subname = '(elmxx_surface_albedo) '

    if (.not. kokkos_state_built) call shr_sys_abort(subname//'ERROR: maps not built')
    if (.not. pftcon_read)        call shr_sys_abort(subname//'ERROR: PFT parameters not read')

    allocate(coszen_col(n_kokkos_col), fracsno(n_kokkos_col), &
             soil_color(n_kokkos_col), h2osoi_top(n_kokkos_col), &
             albsod(n_kokkos_col,numrad), albsoi(n_kokkos_col,numrad), &
             albgrd(n_kokkos_col,numrad), albgri(n_kokkos_col,numrad))
    allocate(tveg(n_kokkos_patch), fwet(n_kokkos_patch), &
             elai(n_kokkos_patch), esai(n_kokkos_patch), &
             patch_col(n_kokkos_patch), patch_ivt(n_kokkos_patch), &
             albd(n_kokkos_patch,numrad), albi(n_kokkos_patch,numrad), &
             fabd(n_kokkos_patch,numrad), fabi(n_kokkos_patch,numrad), &
             ftdd(n_kokkos_patch,numrad), ftid(n_kokkos_patch,numrad), &
             ftii(n_kokkos_patch,numrad), &
             fsun_z(n_kokkos_patch), tlai_z(n_kokkos_patch), &
             fabd_sun_z(n_kokkos_patch), fabi_sun_z(n_kokkos_patch), &
             fabd_sha_z(n_kokkos_patch), fabi_sha_z(n_kokkos_patch), &
             nrad(n_kokkos_patch))

    if (.not. allocated(patch_vcmaxcintsun)) then
       allocate(patch_vcmaxcintsun(n_kokkos_patch), patch_vcmaxcintsha(n_kokkos_patch))
    end if

    call ELMxxGetTVeg(elm, tveg, n_kokkos_patch, ierr);      call check(ierr, subname, 'TVeg')
    call ELMxxGetFwet(elm, fwet, n_kokkos_patch, ierr);      call check(ierr, subname, 'Fwet')
    call ELMxxGetFracSno(elm, fracsno, n_kokkos_col, ierr);  call check(ierr, subname, 'FracSno')

    !-----------------------------------------------------------------
    ! Gather. Column state first, then patch.
    !
    ! Solar zenith angle: nextsw_cday and declin come from the coupler,
    ! exactly as ELM's lnd_comp_mct hands them to elm_drv. Using the model's
    ! own clock instead would drift against the atmosphere's radiation step.
    !-----------------------------------------------------------------
    do kc = 1, n_kokkos_col
       c = col_of_kcol(kc)
       g = lun_gridcell(col_landunit(c))
       coszen_col(kc) = shr_orb_cosz(nextsw_cday, lat(g)*SHR_CONST_PI/180.0_r8, &
                                     lon(g)*SHR_CONST_PI/180.0_r8, declin)
       soil_color(kc) = col_soil_color(c)
       h2osoi_top(kc) = col_h2osoi_vol(c,1)
    end do

    do kp = 1, n_kokkos_patch
       p = patch_of_kpatch(kp)
       c = patch_column(p)
       patch_col(kp) = kcol_of_col(c) + 1     ! kcol_of_col is 0-based
       patch_ivt(kp) = patch_itype(p)
       elai(kp)      = patch_lai(p)
       esai(kp)      = patch_sai(p)
    end do

    call elmxx_surface_albedo_kernel(n_kokkos_col, n_kokkos_patch, npft_param, &
         patch_col, coszen_col, soil_color, h2osoi_top, fracsno, patch_ivt,    &
         elai, esai, tveg, fwet, rhol, rhos, taul, taus, xl,                   &
         albsod, albsoi, albgrd, albgri,                                       &
         albd, albi, fabd, fabi, ftdd, ftid, ftii,                             &
         nrad, tlai_z, fsun_z, fabd_sun_z, fabi_sun_z, fabd_sha_z, fabi_sha_z, &
         patch_vcmaxcintsun, patch_vcmaxcintsha)

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

    deallocate(coszen_col, fracsno, soil_color, h2osoi_top, &
               albsod, albsoi, albgrd, albgri)
    deallocate(tveg, fwet, elai, esai, patch_col, patch_ivt, &
               albd, albi, fabd, fabi, ftdd, ftid, ftii, &
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
