module elmxxPhotosynMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! The inputs the ported Photosynthesis kernel needs that no other ELMxx
  ! kernel produced, and the push that hands them to the Kokkos side.
  !
  ! The leaf model itself is C++ (PhotosynthesisImpl.h) and runs inside the
  ! CanopyFluxes Newton iteration, because rssun must be recomputed at every
  ! iterate. What is left for the Fortran side is four small ports that ELMxx
  ! simply did not have:
  !
  !   dayl_factor   daylength scaling, ELM CanopyFluxesMod:549 over
  !                 DaylengthMod's daylength() and elm_initializeMod's
  !                 max_dayl
  !   t10           10-day running mean of 2 m temperature. ELM accumulates it
  !                 with accumulMod's 'runmean'; the recurrence is reproduced
  !                 exactly here -- see the note on it below
  !   oair / cair   O2 and CO2 partial pressures,
  !                 lnd_import_export.F90:1283 and :1360
  !   the PFT parameter push, once at init
  !
  ! vcmaxcintsun/sha come from elmxxSurfaceAlbedoMod, which already computes
  ! and validates them.
  !-----------------------------------------------------------------------

  use shr_kind_mod        , only : r8 => shr_kind_r8
  use shr_sys_mod         , only : shr_sys_abort, shr_sys_flush
  use shr_const_mod       , only : SHR_CONST_PI, SHR_CONST_CDAY, SHR_CONST_TKFRZ
  use elmxxSpmdMod        , only : masterproc, iam
  use elmxxSubgridMod     , only : col_landunit, lun_gridcell, patch_column, &
                                   patch_itype
  use elmxxForcingMod     , only : forc_pbot
  use elmxxPftconMod      , only : pftcon_read, npft_param, &
                                   c3psn, leafcn, flnr, fnitr, slatop, &
                                   qe_ps, theta_cj, bbbopt, mbbopt, &
                                   photo_uniform, n_photo_uniform
  use elmxxSurfaceAlbedoMod, only : patch_vcmaxcintsun, patch_vcmaxcintsha, &
                                    surface_albedo_built
  use elmxxKokkosStateMod , only : n_kokkos_patch, patch_of_kpatch, &
                                   kokkos_state_built
  use elmxx_mod           , only : ELMxxType, ELMXX_SUCCESS, &
                                   ELMxxSetLatRad, ELMxxSetMaxDayl, &
                                   ELMxxSetPhotoUniform, ELMxxSetPhotoPftParams, &
                                   ELMxxSetUsePhotosynthesis, &
                                   ELMxxSetDaylFactor, ELMxxSetT10, &
                                   ELMxxSetOair, ELMxxSetCair, &
                                   ELMxxSetVcmaxcintSun, ELMxxSetVcmaxcintSha, &
                                   ELMxxGetTRef2m

  implicit none
  save
  private

  ! ELM DaylengthMod: seconds per radian of Earth rotation.
  real(r8), parameter :: secs_per_radian = 13750.9871_r8
  ! ELM elm_initializeMod:674 -- maximum declination for present-day orbital
  ! parameters, +/- 23.4667 degrees, negative in the southern hemisphere.
  real(r8), parameter :: max_decl_mag = 0.409571_r8
  ! ELM elm_varcon: constant atmospheric O2 molar ratio.
  real(r8), parameter :: o2_molar_const = 0.209_r8

  ! Latitude in RADIANS, saved at init. Deliberately not fetched from
  ! elmxxMod at use time: elmxxMod uses this module, so reaching back into it
  ! is a circular dependency that will not compile.
  real(r8), allocatable :: lat_rad(:)    ! per gridcell [radians]
  real(r8), allocatable :: max_dayl(:)   ! per gridcell [s]
  real(r8), allocatable :: t10(:)        ! per packed patch [K]
  integer, public       :: t10_period = 0 ! running-mean period [steps]

  logical, public :: photosyn_built = .false.

  public :: elmxx_push_photosyn_statics
  public :: elmxx_photosyn_init
  public :: elmxx_photosyn_seed

contains

  !-----------------------------------------------------------------------
  subroutine elmxx_push_photosyn_statics(elm, logunit)
    !
    ! One-time push of the gridcell statics the device kernel needs: latitude
    ! in radians and the maximum daylength, resolved per patch.
    !
    implicit none
    type(ELMxxType), intent(in) :: elm
    integer, intent(in) :: logunit
    real(r8), allocatable :: b(:)
    integer :: kp, p, c, g, ierr
    character(len=*), parameter :: subname = '(elmxx_push_photosyn_statics) '

    if (.not. photosyn_built) return
    if (n_kokkos_patch <= 0) return
    allocate(b(n_kokkos_patch))
    do kp = 1, n_kokkos_patch
       p = patch_of_kpatch(kp); c = patch_column(p)
       g = lun_gridcell(col_landunit(c))
       b(kp) = lat_rad(g)
    end do
    call ELMxxSetLatRad(elm, b, n_kokkos_patch, ierr)
    if (ierr /= ELMXX_SUCCESS) call shr_sys_abort(subname//'ERROR: SetLatRad')
    do kp = 1, n_kokkos_patch
       p = patch_of_kpatch(kp); c = patch_column(p)
       g = lun_gridcell(col_landunit(c))
       b(kp) = max_dayl(g)
    end do
    call ELMxxSetMaxDayl(elm, b, n_kokkos_patch, ierr)
    if (ierr /= ELMXX_SUCCESS) call shr_sys_abort(subname//'ERROR: SetMaxDayl')
    deallocate(b)

  end subroutine elmxx_push_photosyn_statics


  !-----------------------------------------------------------------------
  pure function daylength(lat, decl) result(dayl)
    !
    ! ELM DaylengthMod's daylength(), radians in, seconds out.
    !
    ! The latitude clamp is not cosmetic: without it cos(lat) can go negative
    ! at the pole and the acos argument leaves [-1,1].
    !
    implicit none
    real(r8), intent(in) :: lat, decl
    real(r8) :: dayl
    real(r8) :: my_lat, temp
    real(r8), parameter :: lat_epsilon = 10.0_r8 * epsilon(1.0_r8)
    real(r8), parameter :: pole        = SHR_CONST_PI / 2.0_r8
    real(r8), parameter :: offset_pole = pole - lat_epsilon

    my_lat = min(offset_pole, max(-1.0_r8 * offset_pole, lat))
    temp   = -(sin(my_lat)*sin(decl)) / (cos(my_lat)*cos(decl))
    temp   = min(1.0_r8, max(-1.0_r8, temp))
    dayl   = 2.0_r8 * secs_per_radian * acos(temp)
  end function daylength

  !-----------------------------------------------------------------------
  subroutine elmxx_photosyn_init(lat, dtime, logunit)
    !
    ! Maximum daylength per gridcell, and the running-mean state.
    !
    implicit none
    real(r8), intent(in) :: lat(:)          ! gridcell latitudes, DEGREES
    real(r8), intent(in) :: dtime           ! model timestep [s]
    integer , intent(in) :: logunit
    integer  :: g, ng
    real(r8) :: latr, max_decl
    character(len=*), parameter :: subname = '(elmxx_photosyn_init) '

    if (.not. kokkos_state_built) call shr_sys_abort(subname//'ERROR: maps not built')

    ng = size(lat)
    if (allocated(max_dayl)) deallocate(max_dayl)
    if (allocated(lat_rad))  deallocate(lat_rad)
    allocate(max_dayl(ng), lat_rad(ng))
    do g = 1, ng
       lat_rad(g) = lat(g) * SHR_CONST_PI / 180.0_r8
       max_decl = max_decl_mag
       if (lat(g) < 0.0_r8) max_decl = -max_decl
       max_dayl(g) = daylength(lat_rad(g), max_decl)
    end do

    ! ELM's accumulator: period is given as -10 (days) and converted with
    ! accumulMod:149, period = -accum_period * seconds_per_day / dtime.
    t10_period = nint(10.0_r8 * SHR_CONST_CDAY / dtime)
    if (allocated(t10)) deallocate(t10)
    allocate(t10(n_kokkos_patch))
    ! VegetationDataType:1425 init_value = TKFRZ + 20.
    t10 = SHR_CONST_TKFRZ + 20.0_r8

    photosyn_built = .true.
    if (masterproc) then
       write(logunit,*) subname,'max_dayl [s] ',minval(max_dayl),' .. ',maxval(max_dayl)
       write(logunit,*) subname,'t10 running-mean period ',t10_period,' steps'
       call shr_sys_flush(logunit)
    end if
  end subroutine elmxx_photosyn_init

  !-----------------------------------------------------------------------
  subroutine elmxx_photosyn_seed(elm, logunit)
    !
    ! The PFT parameters, pushed once. Nine vary by PFT and go per patch; the
    ! other fourteen are uniform and go as scalars -- elmxxPftconMod aborts if
    ! the parameter file ever breaks that, so this is a checked simplification.
    !
    implicit none
    type(ELMxxType), intent(in) :: elm
    integer, intent(in) :: logunit
    integer :: kp, p, ierr, sz(2)
    real(r8), allocatable :: packed(:,:)
    character(len=*), parameter :: subname = '(elmxx_photosyn_seed) '

    if (.not. pftcon_read) then
       call shr_sys_abort(subname//'ERROR: PFT parameters are not read; '// &
            'fparamfile must be set before photosynthesis can run')
    end if

    allocate(packed(n_kokkos_patch, 9))
    do kp = 1, n_kokkos_patch
       p = patch_of_kpatch(kp)
       packed(kp,1) = c3psn(patch_itype(p))
       packed(kp,2) = leafcn(patch_itype(p))
       packed(kp,3) = flnr(patch_itype(p))
       packed(kp,4) = fnitr(patch_itype(p))
       packed(kp,5) = slatop(patch_itype(p))
       packed(kp,6) = qe_ps(patch_itype(p))
       packed(kp,7) = theta_cj(patch_itype(p))
       packed(kp,8) = bbbopt(patch_itype(p))
       packed(kp,9) = mbbopt(patch_itype(p))
    end do
    sz = (/ n_kokkos_patch, 9 /)
    call ELMxxSetPhotoPftParams(elm, packed, sz, ierr)
    call check(ierr, subname, 'PhotoPftParams')
    deallocate(packed)

    call ELMxxSetPhotoUniform(elm, photo_uniform, n_photo_uniform, ierr)
    call check(ierr, subname, 'PhotoUniform')

    call ELMxxSetUsePhotosynthesis(elm, 1, ierr)
    call check(ierr, subname, 'UsePhotosynthesis')

    if (masterproc) then
       write(logunit,*) subname,'rank ',iam,' photosynthesis ENABLED; ', &
            'CanopyFluxes now computes rssun/rssha each Newton iteration'
       call shr_sys_flush(logunit)
    end if
  end subroutine elmxx_photosyn_seed


  !-----------------------------------------------------------------------
  subroutine check(ierr, subname, what)
    implicit none
    integer, intent(in) :: ierr
    character(len=*), intent(in) :: subname, what
    if (ierr /= ELMXX_SUCCESS) then
       call shr_sys_abort(subname//'ERROR: '//trim(what)//' returned a non-success status')
    end if
  end subroutine check

end module elmxxPhotosynMod
