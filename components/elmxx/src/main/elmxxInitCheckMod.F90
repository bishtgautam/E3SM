module elmxxInitCheckMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Write a compact, rank-local Stage 2 initialization snapshot.
  !
  ! ELM restart files already expose ELM's topology, weights, and initialized
  ! canopy state.  This text form exposes the matching ELMxx state without
  ! adding a one-off model I/O format.  tools/compare_stage2_init.py joins the
  ! per-rank snapshots and compares them with an ELM restart from a twin case.
  ! The file is replaced after every coupling date so a twin run of any length
  ! can be checked at its final common date.
  !-----------------------------------------------------------------------

  use shr_kind_mod         , only : r8 => shr_kind_r8
  use shr_sys_mod          , only : shr_sys_abort, shr_sys_flush
  use elmxxSpmdMod         , only : iam
  use elmxxSubgridMod      , only : num_landunits, num_columns, num_patches, &
                                    lun_gridcell, lun_itype, lun_wtgcell, &
                                    col_landunit, col_itype, col_wtlunit, &
                                    patch_column, patch_itype, patch_wtcol
  use elmxxSurfaceStateMod , only : surface_state_built, col_pct_sand, &
                                    col_pct_clay, col_organic, col_soil_color, &
                                    patch_lai, patch_sai, patch_height_top, &
                                    patch_height_bot

  implicit none
  save
  private

  logical :: snapshot_reported = .false.

  public :: elmxx_write_init_snapshot

contains

  !-----------------------------------------------------------------------
  subroutine elmxx_write_init_snapshot(logunit, month, day, global_cell_ids)
    implicit none
    integer, intent(in) :: logunit, month, day
    integer, intent(in) :: global_cell_ids(:)
    integer :: unitn, ios, l, c, p, lev, g
    character(len=64) :: filename

    if (.not. surface_state_built .or. size(global_cell_ids) <= 0) then
       call shr_sys_abort('(elmxx_init_check) ERROR: surface state is not ready')
    end if

    write(filename,'("elmxx_stage2_rank",I6.6,".dat")') iam
    open(newunit=unitn, file=trim(filename), status='replace', action='write', &
         form='formatted', iostat=ios)
    if (ios /= 0) call shr_sys_abort('(elmxx_init_check) ERROR: cannot open '//trim(filename))

    write(unitn,'(A,1X,I0,1X,I0,1X,I0)') 'ELMXX_STAGE2', iam, month, day
    write(unitn,'(A,1X,I0,1X,I0,1X,I0)') 'COUNTS', num_landunits, num_columns, num_patches

    do l = 1, num_landunits
       g = global_cell_ids(lun_gridcell(l))
       write(unitn,'(A,1X,I0,1X,I0,1X,ES24.16)') 'L', g, lun_itype(l), lun_wtgcell(l)
    end do

    do c = 1, num_columns
       g = global_cell_ids(lun_gridcell(col_landunit(c)))
       write(unitn,'(A,1X,I0,1X,I0,1X,I0,1X,ES24.16)', advance='no') &
            'C', g, lun_itype(col_landunit(c)), col_itype(c), col_wtlunit(c)
       do lev = 1, size(col_pct_sand, 2)
          write(unitn,'(1X,ES24.16)', advance='no') col_pct_sand(c,lev)
       end do
       do lev = 1, size(col_pct_clay, 2)
          write(unitn,'(1X,ES24.16)', advance='no') col_pct_clay(c,lev)
       end do
       do lev = 1, size(col_organic, 2)
          write(unitn,'(1X,ES24.16)', advance='no') col_organic(c,lev)
       end do
       write(unitn,'(1X,I0)') col_soil_color(c)
    end do

    do p = 1, num_patches
       c = patch_column(p)
       g = global_cell_ids(lun_gridcell(col_landunit(c)))
       write(unitn,'(A,1X,I0,1X,I0,1X,I0,1X,I0,5(1X,ES24.16))') 'P', g, &
            lun_itype(col_landunit(c)), col_itype(c), patch_itype(p), &
            patch_wtcol(p), patch_lai(p), patch_sai(p), patch_height_top(p), &
            patch_height_bot(p)
    end do

    close(unitn)
    if (.not. snapshot_reported) then
       write(logunit,*) '(elmxx_init_check) writing latest state to ',trim(filename), &
                        ' (month/day ',month,day,')'
       call shr_sys_flush(logunit)
       snapshot_reported = .true.
    end if
  end subroutine elmxx_write_init_snapshot

end module elmxxInitCheckMod
