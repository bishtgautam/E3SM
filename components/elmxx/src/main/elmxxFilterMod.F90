module elmxxFilterMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Static topology filters over ELMxx's local subgrid.
  !
  ! A filter contains 1-based indices into the ELMxx arrays, not a second
  ! layout.  It is deliberately smaller than ELM's filterMod: snow filters
  ! depend on evolving state and belong with the kernels that maintain that
  ! state.  These are the invariant filters needed to pack and dispatch the
  ! Stage 3+ kernels.
  !
  ! `all_*` filters retain zero-weight entities.  ELM creates those entities
  ! for dynamic landunits, so they are required for an exact topology map.
  ! The remaining filters use ELM's normal active convention: every weight in
  ! the path from the gridcell to the entity must be positive.
  !-----------------------------------------------------------------------

  use shr_kind_mod    , only : r8 => shr_kind_r8
  use shr_sys_mod     , only : shr_sys_abort, shr_sys_flush
  use elmxxSpmdMod    , only : masterproc, iam
  use elmxxSubgridMod , only : num_landunits, num_columns, num_patches, &
                               lun_itype, lun_wtgcell, col_landunit, &
                               col_itype, col_wtlunit, patch_column, &
                               patch_wtcol, istsoil, istdlak, isturb_tbd, &
                               isturb_hd, isturb_md, icol_road_perv

  implicit none
  save
  private

  type, public :: elmxx_filter_type
     integer, pointer :: index(:) => null()
     integer :: count = 0
  end type elmxx_filter_type

  ! Gridcell filters. All locally-owned domain cells are active.
  type(elmxx_filter_type), public :: filter_allg

  ! Landunit filters.
  type(elmxx_filter_type), public :: filter_alll, filter_activel
  type(elmxx_filter_type), public :: filter_urbanl, filter_nourbanl

  ! Column filters.
  type(elmxx_filter_type), public :: filter_allc, filter_activec
  type(elmxx_filter_type), public :: filter_lakec, filter_nolakec
  type(elmxx_filter_type), public :: filter_soilc, filter_hydrologyc
  type(elmxx_filter_type), public :: filter_urbanc, filter_nourbanc

  ! Patch filters.
  type(elmxx_filter_type), public :: filter_allp, filter_activep
  type(elmxx_filter_type), public :: filter_lakep, filter_nolakep
  type(elmxx_filter_type), public :: filter_soilp, filter_natvegp
  type(elmxx_filter_type), public :: filter_urbanp, filter_nourbanp
  type(elmxx_filter_type), public :: filter_nolakeurbanp

  logical, public :: filters_built = .false.

  public :: elmxx_build_filters
  public :: elmxx_filters_clean

contains

  !-----------------------------------------------------------------------
  subroutine elmxx_build_filters(logunit, ncells)
    !
    ! Build every static filter after the topology and its weights exist.
    !
    implicit none
    integer, intent(in) :: logunit, ncells
    integer :: l, c, p
    logical, allocatable :: l_active(:), c_active(:), p_active(:)
    logical, allocatable :: mask(:)

    if (ncells <= 0 .or. num_landunits <= 0 .or. num_columns <= 0 .or. &
        num_patches <= 0) then
       call shr_sys_abort('(elmxx_build_filters) ERROR: subgrid is not ready')
    end if

    call elmxx_filters_clean()

    allocate(l_active(num_landunits), c_active(num_columns), p_active(num_patches))
    l_active = lun_wtgcell > 0.0_r8
    do c = 1, num_columns
       c_active(c) = l_active(col_landunit(c)) .and. col_wtlunit(c) > 0.0_r8
    end do
    do p = 1, num_patches
       p_active(p) = c_active(patch_column(p)) .and. patch_wtcol(p) > 0.0_r8
    end do

    allocate(mask(ncells)); mask = .true.; call set_filter(filter_allg, mask); deallocate(mask)

    allocate(mask(num_landunits))
    mask = .true.; call set_filter(filter_alll, mask)
    call set_filter(filter_activel, l_active)
    mask = is_urban_landunit(lun_itype)
    call set_filter(filter_urbanl, mask)
    call set_filter(filter_nourbanl, .not. mask)
    deallocate(mask)

    allocate(mask(num_columns))
    mask = .true.; call set_filter(filter_allc, mask)
    call set_filter(filter_activec, c_active)
    mask = lun_itype(col_landunit) == istdlak
    call set_filter(filter_lakec, mask)
    call set_filter(filter_nolakec, .not. mask)
    mask = lun_itype(col_landunit) == istsoil
    call set_filter(filter_soilc, mask)
    mask = (lun_itype(col_landunit) == istsoil) .or. col_itype == icol_road_perv
    call set_filter(filter_hydrologyc, mask)
    mask = is_urban_landunit(lun_itype(col_landunit))
    call set_filter(filter_urbanc, mask)
    call set_filter(filter_nourbanc, .not. mask)
    deallocate(mask)

    allocate(mask(num_patches))
    mask = .true.; call set_filter(filter_allp, mask)
    call set_filter(filter_activep, p_active)
    mask = lun_itype(col_landunit(patch_column)) == istdlak
    call set_filter(filter_lakep, mask)
    call set_filter(filter_nolakep, .not. mask)
    mask = lun_itype(col_landunit(patch_column)) == istsoil
    call set_filter(filter_soilp, mask)
    call set_filter(filter_natvegp, mask)
    mask = is_urban_landunit(lun_itype(col_landunit(patch_column)))
    call set_filter(filter_urbanp, mask)
    call set_filter(filter_nourbanp, .not. mask)
    mask = (lun_itype(col_landunit(patch_column)) /= istdlak) .and. .not. &
           is_urban_landunit(lun_itype(col_landunit(patch_column)))
    call set_filter(filter_nolakeurbanp, mask)
    deallocate(mask, l_active, c_active, p_active)

    filters_built = .true.
    call check_filters(logunit, ncells)

    if (masterproc) then
       write(logunit,*) '(elmxx_filters) rank ',iam,' active landunits/columns/patches ', &
                        filter_activel%count, filter_activec%count, filter_activep%count
       call shr_sys_flush(logunit)
    end if

  end subroutine elmxx_build_filters

  !-----------------------------------------------------------------------
  subroutine set_filter(this_filter, mask)
    implicit none
    type(elmxx_filter_type), intent(inout) :: this_filter
    logical, intent(in) :: mask(:)
    integer :: i, n

    if (associated(this_filter%index)) deallocate(this_filter%index)
    n = count(mask)
    allocate(this_filter%index(n))
    this_filter%count = n
    n = 0
    do i = 1, size(mask)
       if (mask(i)) then
          n = n + 1
          this_filter%index(n) = i
       end if
    end do
  end subroutine set_filter

  !-----------------------------------------------------------------------
  elemental logical function is_urban_landunit(itype)
    implicit none
    integer, intent(in) :: itype
    is_urban_landunit = itype == isturb_tbd .or. itype == isturb_hd .or. &
                         itype == isturb_md
  end function is_urban_landunit

  !-----------------------------------------------------------------------
  subroutine check_filters(logunit, ncells)
    ! Filters are ordered increasing by source index, and their complementary
    ! pairs must partition the complete static topology.
    implicit none
    integer, intent(in) :: logunit, ncells

    call check_partition(logunit, 'landunit urban/non-urban', num_landunits, &
                         filter_urbanl, filter_nourbanl)
    call check_partition(logunit, 'column lake/non-lake', num_columns, &
                         filter_lakec, filter_nolakec)
    call check_partition(logunit, 'column urban/non-urban', num_columns, &
                         filter_urbanc, filter_nourbanc)
    call check_partition(logunit, 'patch lake/non-lake', num_patches, &
                         filter_lakep, filter_nolakep)
    call check_partition(logunit, 'patch urban/non-urban', num_patches, &
                         filter_urbanp, filter_nourbanp)

    if (filter_allg%count /= ncells .or. filter_alll%count /= num_landunits .or. &
        filter_allc%count /= num_columns .or. filter_allp%count /= num_patches) then
       call shr_sys_abort('(elmxx_filters) ERROR: incomplete all-entity filter')
    end if
  end subroutine check_filters

  !-----------------------------------------------------------------------
  subroutine check_partition(logunit, name, nentity, first, second)
    implicit none
    integer, intent(in) :: logunit, nentity
    character(len=*), intent(in) :: name
    type(elmxx_filter_type), intent(in) :: first, second
    logical, allocatable :: seen(:)
    integer :: i

    allocate(seen(nentity)); seen = .false.
    do i = 1, first%count
       if (seen(first%index(i))) call shr_sys_abort('(elmxx_filters) ERROR: duplicate '//trim(name))
       seen(first%index(i)) = .true.
    end do
    do i = 1, second%count
       if (seen(second%index(i))) call shr_sys_abort('(elmxx_filters) ERROR: overlapping '//trim(name))
       seen(second%index(i)) = .true.
    end do
    if (.not. all(seen)) call shr_sys_abort('(elmxx_filters) ERROR: incomplete '//trim(name))
    deallocate(seen)
  end subroutine check_partition

  !-----------------------------------------------------------------------
  subroutine elmxx_filters_clean()
    implicit none
    call clean_filter(filter_allg)
    call clean_filter(filter_alll); call clean_filter(filter_activel)
    call clean_filter(filter_urbanl); call clean_filter(filter_nourbanl)
    call clean_filter(filter_allc); call clean_filter(filter_activec)
    call clean_filter(filter_lakec); call clean_filter(filter_nolakec)
    call clean_filter(filter_soilc); call clean_filter(filter_hydrologyc)
    call clean_filter(filter_urbanc); call clean_filter(filter_nourbanc)
    call clean_filter(filter_allp); call clean_filter(filter_activep)
    call clean_filter(filter_lakep); call clean_filter(filter_nolakep)
    call clean_filter(filter_soilp); call clean_filter(filter_natvegp)
    call clean_filter(filter_urbanp); call clean_filter(filter_nourbanp)
    call clean_filter(filter_nolakeurbanp)
    filters_built = .false.
  end subroutine elmxx_filters_clean

  !-----------------------------------------------------------------------
  subroutine clean_filter(this_filter)
    implicit none
    type(elmxx_filter_type), intent(inout) :: this_filter
    if (associated(this_filter%index)) deallocate(this_filter%index)
    this_filter%count = 0
  end subroutine clean_filter

end module elmxxFilterMod
