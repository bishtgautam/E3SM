module elmxxSubgridMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Builds the ELMxx subgrid hierarchy: gridcell -> landunit -> column -> patch.
  !
  ! The rules and the ordering below are ELM's, deliberately. Stage 2's exit
  ! criterion is that this construction matches ELM exactly from the same
  ! surface dataset, so ELM is the specification here even where a different
  ! choice would look tidier. Rules were read out of
  ! components/elm/src/main/{subgridMod,initGridCellsMod,surfrdMod}.F90 -- no
  ! code was copied.
  !
  ! ORDERING (initGridCellsMod.F90's own comment): clump index varies most
  ! slowly, then LANDUNIT TYPE, then gridcell, then column, then patch. So it is
  ! landunit-type-major, NOT gridcell-major: every gridcell's natural-vegetation
  ! landunit comes before any gridcell's urban landunit. ELMxx has one clump per
  ! rank, so the outer loop drops out.
  !
  ! WHICH LANDUNITS EXIST (subgridMod.F90):
  !   istsoil     natural veg  ALWAYS, in every gridcell regardless of weight,
  !                            with all natpft patches allocated -- ELM does this
  !                            so the landunit can come into existence later
  !                            under dynamic landunits.
  !   isturb_tbd/hd/md         when URBAN_REGION_ID /= 0, again REGARDLESS OF
  !                            WEIGHT, 5 columns each. A 0%-weight urban landunit
  !                            is inactive, not absent. This is the trap URBANxx
  !                            hit: 93% of urban gridcells in f19 have at least
  !                            one density type at exactly 0%.
  !   istdlak     lake         only where weight > 0 (no dynamic expansion)
  !   istwet      wetland      only where weight > 0
  !   istice      glacier      only where weight > 0
  !   istcrop                  only if create_crop_landunit; not in SP mode, and
  !                            ELM then requires PCT_CROP == 0 everywhere.
  !
  ! GLACIER: built into the subgrid so the comparison against ELM stays valid,
  ! but no physics will run on it -- glacier is out of scope (STATUS.md §F). Do
  ! NOT drop it and renormalize the remaining weights: that would silently give
  ! ELMxx a different land surface from ELM and destroy the one exact check
  ! Stage 2 has.
  !-----------------------------------------------------------------------

  use shr_kind_mod    , only : r8 => shr_kind_r8
  use shr_sys_mod     , only : shr_sys_abort, shr_sys_flush
  use elmxxSpmdMod    , only : masterproc, iam, mpicom_lnd
  use elmxxSurfdataMod, only : numurbl, natpft, &
                               pct_natveg, pct_crop, pct_lake, pct_wetland, &
                               pct_glacier, pct_urban, pct_nat_pft, urban_region_id

  implicit none
  save
  private

#include <mpif.h>

  !--------------------------------------------------------------------------
  ! Landunit type codes -- ELM's values (landunit_varcon.F90), so that dumps can
  ! be compared against ELM's without a translation table.
  !--------------------------------------------------------------------------
  integer, parameter, public :: istsoil    = 1
  integer, parameter, public :: istcrop    = 2
  integer, parameter, public :: istice     = 3
  integer, parameter, public :: istdlak    = 5
  integer, parameter, public :: istwet     = 6
  integer, parameter, public :: isturb_tbd = 7
  integer, parameter, public :: isturb_hd  = 8
  integer, parameter, public :: isturb_md  = 9

  !--------------------------------------------------------------------------
  ! Urban column type codes (column_varcon.F90): isturb_MIN*10 + n
  !--------------------------------------------------------------------------
  integer, parameter, public :: icol_roof        = 71
  integer, parameter, public :: icol_sunwall     = 72
  integer, parameter, public :: icol_shadewall   = 73
  integer, parameter, public :: icol_road_imperv = 74
  integer, parameter, public :: icol_road_perv   = 75

  integer, parameter, public :: maxpatch_urb = 5   ! columns per urban landunit
  integer, parameter :: urban_invalid_region = 0

  !--------------------------------------------------------------------------
  ! Subgrid, for the cells this rank owns
  !--------------------------------------------------------------------------
  integer, public :: num_landunits = 0
  integer, public :: num_columns   = 0
  integer, public :: num_patches   = 0

  integer , public, pointer :: lun_gridcell(:) => null()  ! 1-based local cell index
  integer , public, pointer :: lun_itype(:)    => null()  ! istsoil, isturb_tbd, ...
  real(r8), public, pointer :: lun_wtgcell(:)  => null()  ! weight on the gridcell (0-1)

  integer , public, pointer :: col_landunit(:) => null()
  integer , public, pointer :: col_itype(:)    => null()
  real(r8), public, pointer :: col_wtlunit(:)  => null()  ! weight on the landunit

  integer , public, pointer :: patch_column(:) => null()
  integer , public, pointer :: patch_itype(:)  => null()  ! PFT index, 0-based
  real(r8), public, pointer :: patch_wtcol(:)  => null()  ! weight on the column

  logical, public :: subgrid_built = .false.

  public :: elmxx_build_subgrid
  public :: elmxx_subgrid_clean

contains

  !-----------------------------------------------------------------------
  subroutine elmxx_build_subgrid(logunit, ncells)
    !
    ! !DESCRIPTION:
    ! Build the subgrid for this rank's cells. Two passes: count, then fill,
    ! so the arrays are sized exactly and the two loops can be checked against
    ! each other.
    !
    implicit none
    !
    integer, intent(in) :: logunit
    integer, intent(in) :: ncells    ! cells owned by this rank
    !
    integer :: nl, nc, np
    character(len=*), parameter :: subname = '(elmxx_build_subgrid) '

    ! ELM requires PCT_CROP == 0 when there is no crop landunit (cft_size == 0),
    ! and SP mode has none. Anything else means this dataset needs crop support
    ! that does not exist here, and silently folding the area elsewhere would
    ! change the land surface.
    if (any(pct_crop(1:ncells) > 0.0_r8)) then
       call shr_sys_abort(subname//' ERROR: PCT_CROP > 0 but ELMxx has no crop landunit')
    end if

    call count_subgrid(ncells, num_landunits, num_columns, num_patches)

    allocate(lun_gridcell(num_landunits), lun_itype(num_landunits), &
             lun_wtgcell(num_landunits))
    allocate(col_landunit(num_columns), col_itype(num_columns), &
             col_wtlunit(num_columns))
    allocate(patch_column(num_patches), patch_itype(num_patches), &
             patch_wtcol(num_patches))

    call fill_subgrid(ncells, nl, nc, np)

    ! The two passes must agree exactly; if they do not, one of them has a rule
    ! the other does not and every index below is suspect.
    if (nl /= num_landunits .or. nc /= num_columns .or. np /= num_patches) then
       write(logunit,*) subname,'counted ',num_landunits,num_columns,num_patches
       write(logunit,*) subname,'filled  ',nl,nc,np
       call shr_sys_abort(subname//' ERROR: subgrid count and fill disagree')
    end if

    subgrid_built = .true.

    call report_by_type(logunit)

    write(logunit,*) subname,'rank ',iam,' landunits ',num_landunits, &
                     ' columns ',num_columns,' patches ',num_patches
    call shr_sys_flush(logunit)

    call report_global_totals(logunit)

  end subroutine elmxx_build_subgrid

  !-----------------------------------------------------------------------
  subroutine report_global_totals(logunit)
    !
    ! !DESCRIPTION:
    ! Sum the subgrid across ranks and report the totals on the master.
    !
    ! This is the number that gets compared against ELM, which reports its
    ! subgrid globally ("total number of landunits = ..."). Per-rank counts
    ! cannot be compared directly: the decomposition differs, and only the
    ! master's line reaches lnd.log anyway.
    !
    implicit none
    !
    integer, intent(in) :: logunit
    !
    integer :: mine(3), total(3), ier
    character(len=*), parameter :: subname = '(elmxx_subgrid_global) '

    mine = (/ num_landunits, num_columns, num_patches /)
    call mpi_reduce(mine, total, 3, MPI_INTEGER, MPI_SUM, 0, mpicom_lnd, ier)

    if (masterproc) then
       write(logunit,*) subname,'GLOBAL total landunits = ',total(1)
       write(logunit,*) subname,'GLOBAL total columns   = ',total(2)
       write(logunit,*) subname,'GLOBAL total patches   = ',total(3)
       call shr_sys_flush(logunit)
    end if

  end subroutine report_global_totals

  !-----------------------------------------------------------------------
  subroutine report_by_type(logunit)
    !
    ! !DESCRIPTION:
    ! Break the subgrid down by landunit type.
    !
    ! This is what gets compared against ELM. Totals alone are too weak: two
    ! different compositions can add up the same, and the failure mode this
    ! guards against -- a landunit type included or omitted under the wrong
    ! condition -- moves counts between types while leaving the total intact.
    !
    implicit none
    !
    integer, intent(in) :: logunit
    !
    integer :: l, c, p, t
    integer :: nlun(9), ncol(9), npat(9)
    character(len=12) :: tname(9)
    character(len=*), parameter :: subname = '(elmxx_subgrid_by_type) '

    nlun = 0; ncol = 0; npat = 0
    tname = '            '
    tname(istsoil)    = 'natveg'
    tname(istcrop)    = 'crop'
    tname(istice)     = 'glacier'
    tname(istdlak)    = 'lake'
    tname(istwet)     = 'wetland'
    tname(isturb_tbd) = 'urban_tbd'
    tname(isturb_hd)  = 'urban_hd'
    tname(isturb_md)  = 'urban_md'

    do l = 1, num_landunits
       t = lun_itype(l)
       nlun(t) = nlun(t) + 1
    end do
    do c = 1, num_columns
       t = lun_itype(col_landunit(c))
       ncol(t) = ncol(t) + 1
    end do
    do p = 1, num_patches
       t = lun_itype(col_landunit(patch_column(p)))
       npat(t) = npat(t) + 1
    end do

    do t = 1, 9
       if (nlun(t) > 0) then
          write(logunit,*) subname,'rank ',iam,' ',trim(tname(t)), &
                           ': landunits ',nlun(t),' columns ',ncol(t), &
                           ' patches ',npat(t)
       end if
    end do
    call shr_sys_flush(logunit)

  end subroutine report_by_type

  !-----------------------------------------------------------------------
  logical function urban_valid(g)
    !
    ! !DESCRIPTION:
    ! Whether this gridcell has valid urban parameters. Note this is entirely
    ! independent of PCT_URBAN: a cell with 0% urban still gets its three urban
    ! landunits if the region ID is valid.
    !
    implicit none
    integer, intent(in) :: g

    urban_valid = (urban_region_id(g) /= urban_invalid_region)

  end function urban_valid

  !-----------------------------------------------------------------------
  subroutine count_subgrid(ncells, nlun, ncol, npat)
    !
    implicit none
    integer, intent(in)  :: ncells
    integer, intent(out) :: nlun, ncol, npat
    integer :: g

    nlun = 0; ncol = 0; npat = 0

    do g = 1, ncells
       ! natural vegetation: always
       nlun = nlun + 1
       ncol = ncol + 1
       npat = npat + natpft
       ! urban x3: on valid urban parameters, regardless of weight
       if (urban_valid(g)) then
          nlun = nlun + 3
          ncol = ncol + 3*maxpatch_urb
          npat = npat + 3*maxpatch_urb
       end if
       ! lake, wetland, glacier: only where present
       if (pct_lake(g)    > 0.0_r8) then; nlun = nlun+1; ncol = ncol+1; npat = npat+1; end if
       if (pct_wetland(g) > 0.0_r8) then; nlun = nlun+1; ncol = ncol+1; npat = npat+1; end if
       if (pct_glacier(g) > 0.0_r8) then; nlun = nlun+1; ncol = ncol+1; npat = npat+1; end if
    end do

  end subroutine count_subgrid

  !-----------------------------------------------------------------------
  subroutine fill_subgrid(ncells, nl, nc, np)
    !
    ! !DESCRIPTION:
    ! Fill the subgrid in ELM's order: landunit type outermost, then gridcell.
    !
    implicit none
    integer, intent(in)  :: ncells
    integer, intent(out) :: nl, nc, np
    integer :: g, m, u, ltype

    nl = 0; nc = 0; np = 0

    ! ---- 1. natural vegetation (every gridcell) ----
    do g = 1, ncells
       nl = nl + 1
       lun_gridcell(nl) = g
       lun_itype(nl)    = istsoil
       lun_wtgcell(nl)  = pct_natveg(g) / 100.0_r8

       nc = nc + 1
       col_landunit(nc) = nl
       col_itype(nc)    = 1              ! ELM uses ctype=1 for the soil column
       col_wtlunit(nc)  = 1.0_r8         ! one column, so it is the whole landunit

       ! All natpft patches are allocated whatever their weight, matching ELM.
       ! PCT_NAT_PFT is a percentage of the natural-vegetated landunit, not of
       ! the gridcell, so it becomes the patch weight on the column directly.
       do m = 1, natpft
          np = np + 1
          patch_column(np) = nc
          patch_itype(np)  = m - 1       ! 0-based PFT index, bare ground = 0
          patch_wtcol(np)  = pct_nat_pft(g,m) / 100.0_r8
       end do
    end do

    ! ---- 2. urban, one landunit type at a time ----
    do u = 1, 3
       select case (u)
       case (1); ltype = isturb_tbd
       case (2); ltype = isturb_hd
       case (3); ltype = isturb_md
       end select

       do g = 1, ncells
          if (.not. urban_valid(g)) cycle

          nl = nl + 1
          lun_gridcell(nl) = g
          lun_itype(nl)    = ltype
          lun_wtgcell(nl)  = pct_urban(g,u) / 100.0_r8

          call add_urban_columns(nl, nc, np)
       end do
    end do

    ! ---- 3. lake, then wetland, then glacier ----
    do g = 1, ncells
       if (pct_lake(g) > 0.0_r8) &
            call add_simple_landunit(g, istdlak, pct_lake(g), nl, nc, np)
    end do
    do g = 1, ncells
       if (pct_wetland(g) > 0.0_r8) &
            call add_simple_landunit(g, istwet, pct_wetland(g), nl, nc, np)
    end do
    do g = 1, ncells
       if (pct_glacier(g) > 0.0_r8) &
            call add_simple_landunit(g, istice, pct_glacier(g), nl, nc, np)
    end do

  end subroutine fill_subgrid

  !-----------------------------------------------------------------------
  subroutine add_urban_columns(nl, nc, np)
    !
    ! !DESCRIPTION:
    ! Add the five columns of an urban landunit, one patch each, in ELM's order:
    ! roof, sunwall, shadewall, impervious road, pervious road.
    !
    ! Column weights within the landunit come from WTLUNIT_ROOF and WTROAD_PERV,
    ! which are not read yet -- they are urban physics parameters, not subgrid
    ! composition. Until they are, the weights are left at zero rather than
    ! guessed: a wrong weight here would look plausible and compare badly for a
    ! reason that is hard to see.
    !
    implicit none
    integer, intent(inout) :: nl, nc, np
    integer :: k
    integer, parameter :: ctypes(maxpatch_urb) = &
         (/ icol_roof, icol_sunwall, icol_shadewall, icol_road_imperv, icol_road_perv /)

    do k = 1, maxpatch_urb
       nc = nc + 1
       col_landunit(nc) = nl
       col_itype(nc)    = ctypes(k)
       col_wtlunit(nc)  = 0.0_r8        ! TODO: WTLUNIT_ROOF / WTROAD_PERV

       np = np + 1
       patch_column(np) = nc
       patch_itype(np)  = 0
       patch_wtcol(np)  = 1.0_r8        ! one patch per urban column
    end do

  end subroutine add_urban_columns

  !-----------------------------------------------------------------------
  subroutine add_simple_landunit(g, ltype, pct, nl, nc, np)
    !
    ! !DESCRIPTION:
    ! Add a landunit with a single column and a single patch: lake, wetland or
    ! glacier.
    !
    implicit none
    integer , intent(in)    :: g, ltype
    real(r8), intent(in)    :: pct
    integer , intent(inout) :: nl, nc, np

    nl = nl + 1
    lun_gridcell(nl) = g
    lun_itype(nl)    = ltype
    lun_wtgcell(nl)  = pct / 100.0_r8

    nc = nc + 1
    col_landunit(nc) = nl
    col_itype(nc)    = ltype    ! ELM uses the landunit type as the column type here
    col_wtlunit(nc)  = 1.0_r8

    np = np + 1
    patch_column(np) = nc
    patch_itype(np)  = 0
    patch_wtcol(np)  = 1.0_r8

  end subroutine add_simple_landunit

  !-----------------------------------------------------------------------
  subroutine elmxx_subgrid_clean()
    !
    implicit none

    if (associated(lun_gridcell)) deallocate(lun_gridcell)
    if (associated(lun_itype))    deallocate(lun_itype)
    if (associated(lun_wtgcell))  deallocate(lun_wtgcell)
    if (associated(col_landunit)) deallocate(col_landunit)
    if (associated(col_itype))    deallocate(col_itype)
    if (associated(col_wtlunit))  deallocate(col_wtlunit)
    if (associated(patch_column)) deallocate(patch_column)
    if (associated(patch_itype))  deallocate(patch_itype)
    if (associated(patch_wtcol))  deallocate(patch_wtcol)

    num_landunits = 0; num_columns = 0; num_patches = 0
    subgrid_built = .false.

  end subroutine elmxx_subgrid_clean

end module elmxxSubgridMod
