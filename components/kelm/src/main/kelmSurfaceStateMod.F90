module kelmSurfaceStateMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Materializes the surface-dataset fields on kELM's subgrid.
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
  use kelmSpmdMod    , only : masterproc, iam
  use kelmSurfdataMod, only : nlevsoi, lsmpft, nmonths, pct_sand, pct_clay, &
                               organic, soil_color, monthly_lai, monthly_sai, &
                               monthly_height_top, monthly_height_bot
  use kelmSubgridMod , only : num_columns, num_patches, lun_gridcell, lun_itype, &
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

  public :: kelm_surface_state_init
  public :: kelm_phenology_weights
  public :: kelm_surface_state_clean

contains

  !-----------------------------------------------------------------------
  subroutine kelm_phenology_weights(month, day, m1, m2, w1, w2)
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

  end subroutine kelm_phenology_weights


  !-----------------------------------------------------------------------
  subroutine kelm_surface_state_init(logunit, month, day)
    !
    ! Copy gridcell soil values to their columns and initialize patch
    ! phenology at the component clock's current date.
    !
    implicit none
    integer, intent(in) :: logunit, month, day
    integer :: c, g, nfill

    if (nlevsoi <= 0 .or. num_columns <= 0 .or. num_patches <= 0) then
       call shr_sys_abort('(kelm_surface_state_init) ERROR: surface subgrid is not ready')
    end if

    allocate(col_pct_sand(num_columns, nlevsoi), col_pct_clay(num_columns, nlevsoi), &
             col_organic(num_columns, nlevsoi), col_soil_color(num_columns))
    allocate(patch_lai(num_patches), patch_sai(num_patches), &
             patch_height_top(num_patches), patch_height_bot(num_patches))

    nfill = 0
    do c = 1, num_columns
       g = lun_gridcell(col_landunit(c))
       col_pct_sand(c,:) = pct_sand(g,:)
       col_pct_clay(c,:) = pct_clay(g,:)
       ! ELM's rule (SoilStateType.F90, InitCold): a cell whose top layer has
       ! neither sand nor clay must have none in any layer, and then gets 1%
       ! sand and 1% clay throughout. Some datasets carry such cells (the
       ! global 1 km urban set does); the thermal properties divide by
       ! sand + clay.
       if (col_pct_sand(c,1) + col_pct_clay(c,1) == 0.0_r8) then
          if (any(col_pct_sand(c,:) + col_pct_clay(c,:) /= 0.0_r8)) then
             call shr_sys_abort('(kelm_surface_state_init) ERROR: sand + clay is '// &
                  'zero at the surface but not below it')
          end if
          col_pct_sand(c,:) = 1.0_r8
          col_pct_clay(c,:) = 1.0_r8
          nfill = nfill + 1
       end if
       col_organic(c,:)  = organic(g,:)
       col_soil_color(c) = soil_color(g)
    end do
    if (nfill > 0) then
       write(logunit,*) '(kelm_surface_state_init) rank ', iam, ': ', nfill, &
                        ' columns with no sand or clay set to 1% each, as ELM does'
    end if

    surface_state_built = .true.

    ! Leaf area starts at zero, NOT at the interpolated phenology. ELM gates
    ! SatellitePhenology on doalb (elm_driver.F90, non-CN non-FATES branch) and
    ! does not call it during initialisation for this configuration -- the
    ! initialize2 call is behind use_fates .and. use_fates_sp. So ELM carries
    ! elai = esai = 0, and hence frac_veg_nosno = 0, until the first doalb
    ! step, treating every patch as bare ground until then. Computing phenology
    ! here gave kELM a full canopy from step 0 and routed patches through
    ! CanopyFluxes while ELM was still running BareGroundFluxes.
    patch_lai        = 0.0_r8
    patch_sai        = 0.0_r8
    patch_height_top = 0.0_r8
    patch_height_bot = 0.0_r8

    if (masterproc) then
       write(logunit,*) '(kelm_surface_state_init) mapped soil properties to ', &
                        num_columns,' columns on rank ',iam
       call shr_sys_flush(logunit)
    end if

  end subroutine kelm_surface_state_init


  !-----------------------------------------------------------------------
  subroutine kelm_surface_state_clean()
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

  end subroutine kelm_surface_state_clean

end module kelmSurfaceStateMod
