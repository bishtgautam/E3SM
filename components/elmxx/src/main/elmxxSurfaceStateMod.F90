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
  public :: elmxx_update_phenology
  public :: elmxx_surface_state_clean

contains

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
  subroutine elmxx_update_phenology(logunit, month, day)
    !
    ! Apply ELM's SatellitePhenologyMod monthly interpolation convention.
    ! `day` is the date at the end of the coupling step, matching ELM's
    ! get_curr_date(offset=dtime) call.  ELM deliberately uses a fixed
    ! no-leap month-length table here, so this routine does too.
    !
    implicit none
    integer, intent(in) :: logunit, month, day
    integer :: p, c, g, pft, first_month, second_month, it1
    integer, parameter :: ndaypm(12) = (/ 31, 28, 31, 30, 31, 30, &
                                         31, 31, 30, 31, 30, 31 /)
    real(r8) :: t, wt_first, wt_second

    if (.not. surface_state_built) then
       call shr_sys_abort('(elmxx_update_phenology) ERROR: surface state is not initialized')
    end if
    if (month < 1 .or. month > nmonths .or. day < 1 .or. day > ndaypm(month)) then
       call shr_sys_abort('(elmxx_update_phenology) ERROR: invalid calendar date')
    end if
    if (nmonths /= 12) then
       call shr_sys_abort('(elmxx_update_phenology) ERROR: satellite phenology requires 12 months')
    end if

    ! These statements mirror SatellitePhenologyMod::interpMonthlyVeg.
    t = (real(day, r8) - 0.5_r8) / real(ndaypm(month), r8)
    it1 = int(t + 0.5_r8)
    first_month  = month + it1 - 1
    second_month = first_month + 1
    if (first_month < 1) first_month = 12
    if (second_month > 12) second_month = 1
    wt_first  = (real(it1, r8) + 0.5_r8) - t
    wt_second = 1.0_r8 - wt_first

    patch_lai        = 0.0_r8
    patch_sai        = 0.0_r8
    patch_height_top = 0.0_r8
    patch_height_bot = 0.0_r8

    do p = 1, num_patches
       c = patch_column(p)
       if (lun_itype(col_landunit(c)) /= istsoil) cycle

       ! ELM's vegetation filter excludes zero-area natural PFTs.  Retain
       ! their topology for structural parity, but leave their dynamic state
       ! zero so the materialized state agrees with ELM's active patch set.
       if (patch_wtcol(p) <= 0.0_r8) cycle

       ! ELM's `noveg` PFT has index zero and receives zero values rather
       ! than values from the monthly stream.  Its 1-based counterpart here
       ! is therefore one.
       if (patch_itype(p) == 0) cycle
       pft = patch_itype(p) + 1
       if (pft < 1 .or. pft > lsmpft) then
          call shr_sys_abort('(elmxx_update_phenology) ERROR: natural PFT is outside MONTHLY_*')
       end if

       g = lun_gridcell(col_landunit(c))
       patch_lai(p) = wt_first * monthly_lai(g,pft,first_month) + &
                      wt_second * monthly_lai(g,pft,second_month)
       patch_sai(p) = wt_first * monthly_sai(g,pft,first_month) + &
                      wt_second * monthly_sai(g,pft,second_month)
       patch_height_top(p) = wt_first * monthly_height_top(g,pft,first_month) + &
                             wt_second * monthly_height_top(g,pft,second_month)
       patch_height_bot(p) = wt_first * monthly_height_bot(g,pft,first_month) + &
                             wt_second * monthly_height_bot(g,pft,second_month)
    end do

    ! Report only when a new monthly pair is selected.  The values are updated
    ! every coupling step, but logging every 30-minute interpolation obscures
    ! useful initialization diagnostics in multi-year runs.
    if (masterproc .and. first_month /= phenology_month) then
       write(logunit,*) '(elmxx_phenology) date month/day ',month,day, &
                        ' interpolates months ',first_month,second_month, &
                        ' weights ',wt_first,wt_second
       call shr_sys_flush(logunit)
    end if
    phenology_month = first_month

  end subroutine elmxx_update_phenology

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
