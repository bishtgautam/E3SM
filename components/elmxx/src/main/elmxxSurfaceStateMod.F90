module elmxxSurfaceStateMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Materializes the surface-dataset fields on ELMxx's subgrid.
  !
  ! The reader owns values in gridcell order.  This module is deliberately
  ! separate from the reader and the topology builder: it owns the two
  ! mappings that the persistent-state handshake will later pass to C++:
  !
  !   * soil texture, organic matter, and soil color on every column; and
  !   * time-interpolated satellite phenology on natural-vegetation patches.
  !
  ! Nothing crosses the C API here.  Stage 3 owns that one-time transfer.
  ! Keeping these host arrays explicit now makes the Stage 2 initialization
  ! comparison possible before physics can obscure an indexing error.
  !-----------------------------------------------------------------------

  use shr_kind_mod    , only : r8 => shr_kind_r8
  use shr_sys_mod     , only : shr_sys_abort, shr_sys_flush
  use elmxxSpmdMod    , only : masterproc, iam
  use elmxx_mod           , only : ELMxxType, ELMXX_SUCCESS, &
                                   ELMxxSetPhenActive, ELMxxSetMonthlyLai, &
                                   ELMxxSetMonthlySai, ELMxxSetMonthlyHtop, &
                                   ELMxxSetMonthlyHbot
  use elmxxKokkosStateMod , only : n_kokkos_patch, patch_of_kpatch
  use elmxxSurfdataMod, only : nlevsoi, lsmpft, nmonths, pct_sand, pct_clay, &
                               organic, soil_color, monthly_lai, monthly_sai, &
                               monthly_height_top, monthly_height_bot
  use elmxxSubgridMod , only : num_columns, num_patches, lun_gridcell, lun_itype, &
                               col_landunit, patch_column, patch_itype, patch_wtcol, istsoil

  implicit none
  save
  private

  ! Per-column copies of the source fields.  These retain the surfdata's
  ! nlevsoi layers; the 15-layer physics layout is a Stage 3 concern.
  real(r8), public, pointer :: col_pct_sand(:,:) => null() ! (column, nlevsoi)
  real(r8), public, pointer :: col_pct_clay(:,:) => null() ! (column, nlevsoi)
  real(r8), public, pointer :: col_organic(:,:)  => null() ! (column, nlevsoi)
  integer , public, pointer :: col_soil_color(:) => null() ! (column)

  ! Satellite phenology after ELM's monthly interpolation.  Non-natural
  ! patches remain zero: ELM only applies these streams to vegetated patches.
  real(r8), public, pointer :: patch_lai(:)        => null()
  real(r8), public, pointer :: patch_sai(:)        => null()
  real(r8), public, pointer :: patch_height_top(:) => null()
  real(r8), public, pointer :: patch_height_bot(:) => null()

  logical, public :: surface_state_built = .false.
  integer, private :: phenology_month = -1

  public :: elmxx_surface_state_init
  public :: elmxx_push_monthly_phenology
  public :: elmxx_phenology_weights
  public :: elmxx_surface_state_clean

contains

  !-----------------------------------------------------------------------
  subroutine elmxx_phenology_weights(month, day, m1, m2, w1, w2)
    !
    ! ELM's monthly interpolation weights. Split out so the host computes only
    ! these four scalars per step and the device does the interpolation.
    !
    implicit none
    integer , intent(in)  :: month, day
    integer , intent(out) :: m1, m2       ! 0-BASED, for the device
    real(r8), intent(out) :: w1, w2

    real(r8) :: t
    integer  :: it1
    integer, parameter :: ndaypm(12) = (/ 31, 28, 31, 30, 31, 30, &
                                          31, 31, 30, 31, 30, 31 /)

    t   = (real(day, r8) - 0.5_r8) / real(ndaypm(month), r8)
    it1 = int(t + 0.5_r8)
    m1  = month + it1 - 1
    m2  = m1 + 1
    if (m1 < 1)  m1 = 12
    if (m2 > 12) m2 = 1
    w1  = (real(it1, r8) + 0.5_r8) - t
    w2  = 1.0_r8 - w1
    m1  = m1 - 1
    m2  = m2 - 1

  end subroutine elmxx_phenology_weights

  !-----------------------------------------------------------------------
  subroutine elmxx_push_monthly_phenology(elm, logunit)
    !
    ! One-time push of the monthly LAI/SAI/height fields, resolved per patch
    ! from (gridcell, PFT). After this the surface dataset never crosses again:
    ! the device holds all twelve months and interpolates in time itself.
    !
    implicit none
    type(ELMxxType), intent(in) :: elm
    integer, intent(in) :: logunit

    real(r8), allocatable :: b(:,:)
    integer , allocatable :: act(:)
    integer :: kp, p, c, g, pft, mm, ierr, sz(2)
    character(len=*), parameter :: subname = '(elmxx_push_monthly_phenology) '

    if (n_kokkos_patch <= 0) return
    allocate(b(n_kokkos_patch, 12), act(n_kokkos_patch))
    sz(1) = n_kokkos_patch; sz(2) = 12

    ! Which patches phenology touches at all -- the host's skip conditions,
    ! evaluated once and shipped as a mask rather than re-tested every step.
    do kp = 1, n_kokkos_patch
       p = patch_of_kpatch(kp); c = patch_column(p)
       act(kp) = 1
       if (lun_itype(col_landunit(c)) /= istsoil) act(kp) = 0
       if (patch_wtcol(p) <= 0.0_r8)              act(kp) = 0
       if (patch_itype(p) == 0)                   act(kp) = 0
    end do
    call ELMxxSetPhenActive(elm, act, n_kokkos_patch, ierr)
    if (ierr /= ELMXX_SUCCESS) call shr_sys_abort(subname//'ERROR: SetPhenActive')

    call fill(monthly_lai);  call ELMxxSetMonthlyLai (elm, b, sz, ierr)
    call fill(monthly_sai);  call ELMxxSetMonthlySai (elm, b, sz, ierr)
    call fill(monthly_height_top); call ELMxxSetMonthlyHtop(elm, b, sz, ierr)
    call fill(monthly_height_bot); call ELMxxSetMonthlyHbot(elm, b, sz, ierr)

    deallocate(b, act)

  contains
    subroutine fill(src)
      real(r8), intent(in) :: src(:,:,:)
      b = 0.0_r8
      do kp = 1, n_kokkos_patch
         if (act(kp) == 0) cycle
         p = patch_of_kpatch(kp); c = patch_column(p)
         g = lun_gridcell(col_landunit(c))
         pft = patch_itype(p) + 1
         do mm = 1, 12
            b(kp,mm) = src(g,pft,mm)
         end do
      end do
    end subroutine fill

  end subroutine elmxx_push_monthly_phenology


  !-----------------------------------------------------------------------
  subroutine elmxx_surface_state_init(logunit, month, day)
    !
    ! Copy gridcell soil values to their columns and initialize patch
    ! phenology at the component clock's current date.
    !
    implicit none
    integer, intent(in) :: logunit, month, day
    integer :: c, g

    if (nlevsoi <= 0 .or. num_columns <= 0 .or. num_patches <= 0) then
       call shr_sys_abort('(elmxx_surface_state_init) ERROR: surface subgrid is not ready')
    end if

    allocate(col_pct_sand(num_columns, nlevsoi), col_pct_clay(num_columns, nlevsoi), &
             col_organic(num_columns, nlevsoi), col_soil_color(num_columns))
    allocate(patch_lai(num_patches), patch_sai(num_patches), &
             patch_height_top(num_patches), patch_height_bot(num_patches))

    do c = 1, num_columns
       g = lun_gridcell(col_landunit(c))
       col_pct_sand(c,:) = pct_sand(g,:)
       col_pct_clay(c,:) = pct_clay(g,:)
       col_organic(c,:)  = organic(g,:)
       col_soil_color(c) = soil_color(g)
    end do

    surface_state_built = .true.

    ! Leaf area starts at zero, NOT at the interpolated phenology. ELM gates
    ! SatellitePhenology on doalb (elm_driver.F90, non-CN non-FATES branch) and
    ! does not call it during initialisation for this configuration -- the
    ! initialize2 call is behind use_fates .and. use_fates_sp. So ELM carries
    ! elai = esai = 0, and hence frac_veg_nosno = 0, until the first doalb
    ! step, treating every patch as bare ground until then. Computing phenology
    ! here gave ELMxx a full canopy from step 0 and routed patches through
    ! CanopyFluxes while ELM was still running BareGroundFluxes.
    patch_lai        = 0.0_r8
    patch_sai        = 0.0_r8
    patch_height_top = 0.0_r8
    patch_height_bot = 0.0_r8

    if (masterproc) then
       write(logunit,*) '(elmxx_surface_state_init) mapped soil properties to ', &
                        num_columns,' columns on rank ',iam
       call shr_sys_flush(logunit)
    end if

  end subroutine elmxx_surface_state_init


  !-----------------------------------------------------------------------
  subroutine elmxx_surface_state_clean()
    implicit none

    if (associated(col_pct_sand))    deallocate(col_pct_sand)
    if (associated(col_pct_clay))    deallocate(col_pct_clay)
    if (associated(col_organic))     deallocate(col_organic)
    if (associated(col_soil_color))  deallocate(col_soil_color)
    if (associated(patch_lai))        deallocate(patch_lai)
    if (associated(patch_sai))        deallocate(patch_sai)
    if (associated(patch_height_top)) deallocate(patch_height_top)
    if (associated(patch_height_bot)) deallocate(patch_height_bot)

    surface_state_built = .false.
    phenology_month = -1

  end subroutine elmxx_surface_state_clean

end module elmxxSurfaceStateMod
