module elmxxSurfdataMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Reads the ELM surface dataset for the cells this rank owns.
  !
  ! Scope is SP mode: the subgrid composition and the soil properties. The
  ! satellite-phenology streams (MONTHLY_LAI/SAI/HEIGHT_TOP/HEIGHT_BOT) are not
  ! read yet -- they are needed by the phenology update, not by subgrid
  ! construction, and land in a later increment. CN/BGC, crop and transient
  ! land-use fields are out of scope entirely.
  !
  ! Reading strategy: like elmxx_read_domain, the whole global field is read on
  ! every rank and the owned cells are then picked out of it. That is simple and
  ! adequate at the resolutions run today, but it is O(global) memory per rank --
  ! at f19 the four MONTHLY_* fields alone would be ~90 MB per rank once they are
  ! added. Move to a decomposed read (pio_read_darray with a decomp built from
  ! natural_id_cells_owned) before this reaches production resolutions.
  !
  ! Indexing: natural_id_cells_owned holds 1-based *global grid* IDs, i.e.
  ! indices into the ni*nj arrays, which is exactly how the surface dataset is
  ! laid out (its `gridcell` dimension spans the full grid, ocean included).
  ! So a cell's surfdata is a direct lookup, no land-only compaction involved.
  !-----------------------------------------------------------------------

  use shr_kind_mod  , only : r8 => shr_kind_r8
  use shr_sys_mod   , only : shr_sys_abort, shr_sys_flush
  use elmxxSpmdMod  , only : masterproc
  use elmxxIO       , only : io_type, pio_subsystem
  use pio

  implicit none
  save
  private

  !--------------------------------------------------------------------------
  ! Dimensions read from the surface dataset
  !--------------------------------------------------------------------------
  integer, public :: numurbl = 0   ! urban density types (3)
  integer, public :: natpft  = 0   ! natural PFTs (17)
  integer, public :: nlevsoi = 0   ! hydrologically active soil layers (10)

  !--------------------------------------------------------------------------
  ! Subgrid composition, per owned cell. Percentages as stored on the file.
  !
  ! PCT_NATVEG/CROP/LAKE/WETLAND/GLACIER are percentages of the gridcell.
  ! PCT_URBAN is a percentage of the gridcell, per density type.
  ! PCT_NAT_PFT is a percentage of the *natural vegetated* landunit, not of the
  ! gridcell -- it sums to 100 within that landunit.
  !--------------------------------------------------------------------------
  real(r8), public, pointer :: pct_natveg(:)    => null()  ! (ncells)
  real(r8), public, pointer :: pct_crop(:)      => null()  ! (ncells)
  real(r8), public, pointer :: pct_lake(:)      => null()  ! (ncells)
  real(r8), public, pointer :: pct_wetland(:)   => null()  ! (ncells)
  real(r8), public, pointer :: pct_glacier(:)   => null()  ! (ncells)
  real(r8), public, pointer :: pct_urban(:,:)   => null()  ! (ncells, numurbl)
  real(r8), public, pointer :: pct_nat_pft(:,:) => null()  ! (ncells, natpft)

  !--------------------------------------------------------------------------
  ! Soil properties, per owned cell
  !--------------------------------------------------------------------------
  real(r8), public, pointer :: pct_sand(:,:)    => null()  ! (ncells, nlevsoi)
  real(r8), public, pointer :: pct_clay(:,:)    => null()  ! (ncells, nlevsoi)
  real(r8), public, pointer :: organic(:,:)     => null()  ! (ncells, nlevsoi)
  real(r8), public, pointer :: fmax(:)          => null()  ! (ncells)
  integer , public, pointer :: soil_color(:)    => null()  ! (ncells)

  logical, public :: surfdata_read = .false.

  public :: elmxx_read_surfdata
  public :: elmxx_surfdata_clean

contains

  !-----------------------------------------------------------------------
  subroutine elmxx_read_surfdata(iulog, fname, gsize, cell_ids)
    !
    ! !DESCRIPTION:
    ! Read the surface dataset and keep the rows for the cells this rank owns.
    !
    implicit none
    !
    integer         , intent(in) :: iulog        ! log unit
    character(len=*), intent(in) :: fname        ! surface dataset
    integer         , intent(in) :: gsize        ! ni*nj, for cross-checking
    integer         , intent(in) :: cell_ids(:)  ! 1-based global grid IDs owned here
    !
    type(file_desc_t) :: ncid
    integer           :: dimid, status
    integer           :: ngrid                    ! gridcell dimension on the file
    integer           :: ncells                   ! cells owned by this rank
    character(len=*), parameter :: subname = 'elmxx_read_surfdata'

    ncells = size(cell_ids)

    if (masterproc) then
       write(iulog,*) trim(subname),': reading surface dataset ',trim(fname)
       call shr_sys_flush(iulog)
    end if

    status = pio_openfile(pio_subsystem, ncid, io_type, trim(fname), pio_nowrite)
    if (status /= PIO_NOERR) then
       call shr_sys_abort(subname//' ERROR: cannot open surface dataset '//trim(fname))
    end if

    ! ---- dimensions ----
    ngrid   = get_dimlen(ncid, fname, 'gridcell')
    numurbl = get_dimlen(ncid, fname, 'numurbl')
    natpft  = get_dimlen(ncid, fname, 'natpft')
    nlevsoi = get_dimlen(ncid, fname, 'nlevsoi')

    ! The surface dataset spans the whole grid, ocean cells included, so its
    ! gridcell dimension must equal ni*nj from the domain file. If it does not,
    ! cell_ids -- which are indices into that grid -- do not address this file
    ! and every lookup below would be silently wrong.
    if (ngrid /= gsize) then
       write(iulog,*) trim(subname),': gridcell = ',ngrid,' but ni*nj = ',gsize
       call shr_sys_abort(subname//' ERROR: surface dataset grid does not match the domain')
    end if

    allocate(pct_natveg(ncells), pct_crop(ncells), pct_lake(ncells), &
             pct_wetland(ncells), pct_glacier(ncells), fmax(ncells), &
             soil_color(ncells))
    allocate(pct_urban(ncells, numurbl))
    allocate(pct_nat_pft(ncells, natpft))
    allocate(pct_sand(ncells, nlevsoi), pct_clay(ncells, nlevsoi), &
             organic(ncells, nlevsoi))

    ! ---- subgrid composition ----
    call read_gc_real1d(ncid, fname, 'PCT_NATVEG' , ngrid, cell_ids, pct_natveg)
    call read_gc_real1d(ncid, fname, 'PCT_CROP'   , ngrid, cell_ids, pct_crop)
    call read_gc_real1d(ncid, fname, 'PCT_LAKE'   , ngrid, cell_ids, pct_lake)
    call read_gc_real1d(ncid, fname, 'PCT_WETLAND', ngrid, cell_ids, pct_wetland)
    call read_gc_real1d(ncid, fname, 'PCT_GLACIER', ngrid, cell_ids, pct_glacier)
    call read_gc_real2d(ncid, fname, 'PCT_URBAN'  , ngrid, numurbl, cell_ids, pct_urban)
    call read_gc_real2d(ncid, fname, 'PCT_NAT_PFT', ngrid, natpft , cell_ids, pct_nat_pft)

    ! ---- soil properties ----
    call read_gc_real2d(ncid, fname, 'PCT_SAND', ngrid, nlevsoi, cell_ids, pct_sand)
    call read_gc_real2d(ncid, fname, 'PCT_CLAY', ngrid, nlevsoi, cell_ids, pct_clay)
    call read_gc_real2d(ncid, fname, 'ORGANIC' , ngrid, nlevsoi, cell_ids, organic)
    call read_gc_real1d(ncid, fname, 'FMAX'    , ngrid, cell_ids, fmax)
    call read_gc_int1d (ncid, fname, 'SOIL_COLOR', ngrid, cell_ids, soil_color)

    call pio_closefile(ncid)

    surfdata_read = .true.

    if (masterproc) then
       write(iulog,*) trim(subname),': gridcell = ',ngrid,' numurbl = ',numurbl, &
                      ' natpft = ',natpft,' nlevsoi = ',nlevsoi
       call shr_sys_flush(iulog)
    end if

  end subroutine elmxx_read_surfdata

  !-----------------------------------------------------------------------
  integer function get_dimlen(ncid, fname, dimname)
    !
    implicit none
    type(file_desc_t), intent(inout) :: ncid
    character(len=*) , intent(in)    :: fname, dimname
    integer :: dimid, status
    character(len=*), parameter :: subname = 'elmxx_read_surfdata::get_dimlen'

    status = pio_inq_dimid(ncid, trim(dimname), dimid)
    if (status /= PIO_NOERR) then
       call shr_sys_abort(subname//' ERROR: no '//trim(dimname)//' dimension on '//trim(fname))
    end if
    status = pio_inq_dimlen(ncid, dimid, get_dimlen)
    if (status /= PIO_NOERR) then
       call shr_sys_abort(subname//' ERROR: cannot read '//trim(dimname)//' from '//trim(fname))
    end if

  end function get_dimlen

  !-----------------------------------------------------------------------
  subroutine read_gc_real1d(ncid, fname, varname, ngrid, cell_ids, out)
    !
    ! !DESCRIPTION:
    ! Read a (gridcell) real field globally and keep the owned cells.
    !
    implicit none
    type(file_desc_t), intent(inout) :: ncid
    character(len=*) , intent(in)    :: fname, varname
    integer          , intent(in)    :: ngrid
    integer          , intent(in)    :: cell_ids(:)
    real(r8)         , intent(inout) :: out(:)
    !
    real(r8), allocatable :: glob(:)
    integer :: varid, status, i
    character(len=*), parameter :: subname = 'elmxx_read_surfdata::read_gc_real1d'

    allocate(glob(ngrid))
    status = pio_inq_varid(ncid, trim(varname), varid)
    if (status /= PIO_NOERR) then
       call shr_sys_abort(subname//' ERROR: no '//trim(varname)//' on '//trim(fname))
    end if
    status = pio_get_var(ncid, varid, glob)
    if (status /= PIO_NOERR) then
       call shr_sys_abort(subname//' ERROR: cannot read '//trim(varname)//' from '//trim(fname))
    end if

    do i = 1, size(cell_ids)
       out(i) = glob(cell_ids(i))
    end do

    deallocate(glob)

  end subroutine read_gc_real1d

  !-----------------------------------------------------------------------
  subroutine read_gc_int1d(ncid, fname, varname, ngrid, cell_ids, out)
    !
    implicit none
    type(file_desc_t), intent(inout) :: ncid
    character(len=*) , intent(in)    :: fname, varname
    integer          , intent(in)    :: ngrid
    integer          , intent(in)    :: cell_ids(:)
    integer          , intent(inout) :: out(:)
    !
    integer, allocatable :: glob(:)
    integer :: varid, status, i
    character(len=*), parameter :: subname = 'elmxx_read_surfdata::read_gc_int1d'

    allocate(glob(ngrid))
    status = pio_inq_varid(ncid, trim(varname), varid)
    if (status /= PIO_NOERR) then
       call shr_sys_abort(subname//' ERROR: no '//trim(varname)//' on '//trim(fname))
    end if
    status = pio_get_var(ncid, varid, glob)
    if (status /= PIO_NOERR) then
       call shr_sys_abort(subname//' ERROR: cannot read '//trim(varname)//' from '//trim(fname))
    end if

    do i = 1, size(cell_ids)
       out(i) = glob(cell_ids(i))
    end do

    deallocate(glob)

  end subroutine read_gc_int1d

  !-----------------------------------------------------------------------
  subroutine read_gc_real2d(ncid, fname, varname, ngrid, nsecond, cell_ids, out)
    !
    ! !DESCRIPTION:
    ! Read a real field declared VAR(nsecond, gridcell) on the file and keep the
    ! owned cells, as (cell, nsecond).
    !
    ! Dimension order reverses between CDL and Fortran. In CDL the LAST
    ! dimension varies fastest, so VAR(nsecond, gridcell) has gridcell varying
    ! fastest and the matching Fortran array is glob(ngrid, nsecond) -- gridcell
    ! FIRST. Getting this backwards does not fail, it silently transposes: with
    ! 2 cells x 3 urban types it handed each cell a mix of the other's values
    ! and the landunit percentages summed to 125 and 75 instead of 100.
    !
    implicit none
    type(file_desc_t), intent(inout) :: ncid
    character(len=*) , intent(in)    :: fname, varname
    integer          , intent(in)    :: ngrid, nsecond
    integer          , intent(in)    :: cell_ids(:)
    real(r8)         , intent(inout) :: out(:,:)
    !
    real(r8), allocatable :: glob(:,:)
    integer :: varid, status, i, k
    character(len=*), parameter :: subname = 'elmxx_read_surfdata::read_gc_real2d'

    allocate(glob(ngrid, nsecond))
    status = pio_inq_varid(ncid, trim(varname), varid)
    if (status /= PIO_NOERR) then
       call shr_sys_abort(subname//' ERROR: no '//trim(varname)//' on '//trim(fname))
    end if
    status = pio_get_var(ncid, varid, glob)
    if (status /= PIO_NOERR) then
       call shr_sys_abort(subname//' ERROR: cannot read '//trim(varname)//' from '//trim(fname))
    end if

    do i = 1, size(cell_ids)
       do k = 1, nsecond
          out(i,k) = glob(cell_ids(i), k)
       end do
    end do

    deallocate(glob)

  end subroutine read_gc_real2d

  !-----------------------------------------------------------------------
  subroutine elmxx_surfdata_clean()
    !
    implicit none

    if (associated(pct_natveg))  deallocate(pct_natveg)
    if (associated(pct_crop))    deallocate(pct_crop)
    if (associated(pct_lake))    deallocate(pct_lake)
    if (associated(pct_wetland)) deallocate(pct_wetland)
    if (associated(pct_glacier)) deallocate(pct_glacier)
    if (associated(pct_urban))   deallocate(pct_urban)
    if (associated(pct_nat_pft)) deallocate(pct_nat_pft)
    if (associated(pct_sand))    deallocate(pct_sand)
    if (associated(pct_clay))    deallocate(pct_clay)
    if (associated(organic))     deallocate(organic)
    if (associated(fmax))        deallocate(fmax)
    if (associated(soil_color))  deallocate(soil_color)

    surfdata_read = .false.

  end subroutine elmxx_surfdata_clean

end module elmxxSurfdataMod
