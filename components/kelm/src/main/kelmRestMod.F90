module kelmRestMod
  !
  ! kELM restart files (Stage 6, step 2).
  !
  ! WHAT IS WRITTEN. All persistent state, by the design decided on
  ! 2026-10-04: every field of the restart registry (the C++ state
  ! containers' device views, generated list, plus one "hist:<NAME>"
  ! accumulator per active history field) and the few pieces of state that
  ! live in Fortran -- the step counter, which sets the t10 running-mean
  ! window, and the history interval's start date and sample count. Nothing
  ! is curated: the 10-day vs 5+5-day bit-for-bit test proves the set is
  ! complete, and trimming comes only after that passes.
  !
  ! FORMAT. One NetCDF variable per registry field, named exactly as the
  ! registry names it ("natcol:t_grnd"), double precision. Its leading
  ! dimension is its ENTITY KIND ("natcol", "natpatch", "cell", "urblun" --
  ! the registry's kELMRestartKind) over the whole domain; a rank-2 field
  ! (n1, n2) adds a level dim named by extent ("n15"), and is (n1, n2) on
  ! file in C order. Integers are stored as doubles -- exact for I4.
  !
  ! LAYOUT-INDEPENDENT, as ELM's restarts are. An entity's position on file
  ! is (global cell ordinal, then its occurrence within that cell), never
  ! its rank-local packed index, so a restart written on N ranks reads on M.
  ! Two things make that hold: within one cell the packed order is the
  ! subgrid's, whatever rank owns the cell; and per-cell counts come from an
  ! allreduce, so no rank assumes how many entities another cell has. The
  ! packed urban order is (density class, cell), not cell-major -- the
  ! reordering is exactly what this handles. TOPO fields hold rank-local
  ! indices (col_gridcell, filters): they are rebuilt at init from the
  ! subgrid and are neither written nor read.
  !
  ! KINDS NOT MAPPED. shared:urbpoi is never allocated on the coupled path;
  ! if it ever is, write/read abort naming it rather than guess its
  ! entity-to-cell map. Lake columns and patches are mapped (Stage 6.5, L6):
  ! one of each per lake cell (kelmLakeMod). Urban columns and patches are
  ! mapped (Stage 6.7, U7): five of each per urban landunit (kelmUrbanMod).
  !
  ! A BRANCH STARTS HISTORY FRESH, as ELM's does (hist_restart_ncd zeroes
  ! every tape's ntimes and reads no history buffers on nsrBranch): the
  ! "hist:" accumulators and the interval start/sample count are not read,
  ! so the first h0 holds only the branch's own steps and the branch may
  ! change kelm_hist_fincl.
  !
  ! READING IS STRICT. Every allocated, non-TOPO registry field must be on
  ! the file with the same global shape, or the run aborts: under "write all
  ! state" a field that silently fails to restore is a bug, not a tolerance.
  !
  use shr_kind_mod , only : r8 => shr_kind_r8, CL => shr_kind_cl
  use shr_sys_mod  , only : shr_sys_abort, shr_sys_flush
  use shr_file_mod , only : shr_file_getunit, shr_file_freeunit
  use kelmSpmdMod , only : masterproc, iam, npes, mpicom_lnd
  use kelmIO      , only : pio_subsystem, io_type
  use kelmHistMod , only : kelm_hist_get_restart_state, kelm_hist_set_restart_state
  use kelmSubgridMod     , only : lun_gridcell, col_landunit, patch_column
  use kelmLakeMod        , only : n_lake, cell_of_klake
  use kelmKokkosStateMod , only : n_kokkos_col, n_kokkos_patch, n_kokkos_urb, &
                                   col_of_kcol, patch_of_kpatch, lun_of_kurb
  use kelmLndExportMod  , only : kelm_lnd2atm_ready, kelm_lnd2atm_cell_table, &
                                   NL2A => NFIELD
  use kelm_mod    , only : kELMType, KELM_SUCCESS, &
                            kELMRestartFieldCount, kELMRestartFieldInfo, &
                            kELMRestartFieldGet, kELMRestartFieldSet, &
                            kELMHistoryGetNacs, kELMHistorySetNacs
  use pio
  use iso_c_binding, only : c_int, c_double
  implicit none
  private
#include <mpif.h>

  public :: kelm_rest_write
  public :: kelm_rest_read
  public :: kelm_rest_peek_nstep
  public :: kelm_rest_read_lnd2atm
  public :: kelm_rest_rpointer_read

  integer, parameter :: NAMELEN = 128

  ! kELMRestartKind (kELM.h). URBCOL/URBPATCH are retired: every urban
  ! surface record is indexed by landunit (URBLUN) since plan B.
  integer, parameter :: KIND_TOPO = 0, KIND_CELL = 1, KIND_NATCOL = 2, &
                        KIND_NATPATCH = 3, KIND_URBLUN = 4, KIND_URBCOL = 5, &
                        KIND_URBPATCH = 6, KIND_LAKECOL = 7, &
                        KIND_LAKEPATCH = 8, NKIND = 9
  character(len=10), parameter :: kind_name(0:NKIND) = (/ &
       'topo      ', 'cell      ', 'natcol    ', 'natpatch  ', 'urblun    ', &
       'urbcol    ', 'urbpatch  ', 'lakecol   ', 'lakepatch ', 'landunit  ' /)

  ! One field's layout, agreed across ranks.
  type :: field_t
     character(len=NAMELEN) :: name
     integer :: n1 = 0, n2 = 0     ! this rank's extents (0, 0 = not allocated here)
     integer :: kind = KIND_TOPO
     integer :: n2g = 0            ! global second extent (max over ranks)
     logical :: active = .false.   ! allocated somewhere, and restarted
  end type field_t

  ! Global positions of this rank's entities, per kind, and the decompositions
  ! built from them (one per kind and second extent), for one write or read.
  type :: kindmap_t
     logical :: built = .false.
     integer :: n = 0                        ! this rank's entities
     integer :: nglob = 0
     integer, allocatable :: pos(:)          ! (n local), 1-based global position
  end type kindmap_t
  type(kindmap_t) :: kmap(0:NKIND)

  integer, parameter :: MAXDEC = 64
  integer :: ndec = 0
  integer :: dec_kind(MAXDEC), dec_n2(MAXDEC)
  type(io_desc_t) :: dec_iodesc(MAXDEC)

contains

  !-----------------------------------------------------------------------
  subroutine kelm_rest_write(elm, caseid, inst_suffix, yr, mon, day, tod, &
                              nstep, ncell, ncell_global, logunit)
    !
    ! Write <caseid>.kelm.r.YYYY-MM-DD-SSSSS.nc, dated with the END of the
    ! step just taken (the time the next run starts from, as ELM names its
    ! restarts), and point rpointer.lnd<inst_suffix> at it. Collective.
    !
    type(kELMType) , intent(in) :: elm
    character(len=*), intent(in) :: caseid, inst_suffix
    integer         , intent(in) :: yr, mon, day, tod, nstep
    integer         , intent(in) :: ncell, ncell_global   ! cells owned here / in the domain
    integer         , intent(in) :: logunit

    type(file_desc_t) :: file
    type(var_desc_t)  :: vid_scalar(5), vid_l2a
    logical :: with_l2a
    real(r8), allocatable :: tab(:,:)
    integer :: e, k
    type(var_desc_t), allocatable :: vid(:)
    type(field_t), allocatable :: f(:)
    integer :: nfld, i, ierr, d1, d2, nwritten
    integer(c_int) :: st
    integer :: nacs_f, hy, hm, hd, nacs_c
    integer :: dim_ext(64), dim_id(64), ndims
    integer :: kdim_id(0:NKIND)
    real(c_double), allocatable :: buf(:)
    character(len=CL) :: fname
    character(len=16) :: rdate   ! YYYY-MM-DD-SSSSS, as ELM logs it
    character(len=*), parameter :: subname = '(kelm_rest_write) '

    write(fname,'(a,".kelm.r.",i4.4,"-",i2.2,"-",i2.2,"-",i5.5,".nc")') &
         trim(caseid), yr, mon, day, tod

    ! ELM's restart lines (restFileMod.F90) bracket kELM's own summary.
    if (masterproc) then
       write(rdate,'(i4.4,"-",i2.2,"-",i2.2,"-",i5.5)') yr, mon, day, tod
       write(logunit,*)'writing restart file ',trim(fname),' for model date = ',rdate
       write(logunit,*)
       call shr_sys_flush(logunit)
    end if

    call layout_begin(elm, f, nfld, ncell, ncell_global, .false., subname)
    allocate(vid(nfld))

    ierr = PIO_createfile(pio_subsystem, file, io_type, trim(fname), PIO_CLOBBER)
    if (ierr /= PIO_NOERR) call shr_sys_abort(subname//'ERROR: cannot create '//trim(fname))

    ndims = 0
    kdim_id = -1
    do i = 1, nfld
       if (.not. f(i)%active) cycle
       if (kdim_id(f(i)%kind) < 0) then
          ierr = PIO_def_dim(file, trim(kind_name(f(i)%kind)), kmap(f(i)%kind)%nglob, &
                             kdim_id(f(i)%kind))
          if (ierr /= PIO_NOERR) call shr_sys_abort(subname//'ERROR: def_dim '// &
               trim(kind_name(f(i)%kind)))
       end if
       d1 = kdim_id(f(i)%kind)
       if (f(i)%n2g > 0) then
          d2 = dimid_for(f(i)%n2g)
          ! Fortran dim order is the reverse of C: (n2, n1) here is (n1, n2) on file.
          ierr = PIO_def_var(file, trim(f(i)%name), PIO_double, (/d2, d1/), vid(i))
       else
          ierr = PIO_def_var(file, trim(f(i)%name), PIO_double, (/d1/), vid(i))
       end if
       if (ierr /= PIO_NOERR) call shr_sys_abort(subname//'ERROR: def_var '//trim(f(i)%name))
    end do

    ! The coupler export of the state being written (kELM-only): a continue
    ! run sends it at init, as ELM sends the export of its restart state.
    ! Not a registry field, and read only for the init export, so restarts
    ! written before it existed still read.
    with_l2a = kelm_lnd2atm_ready() .and. kdim_id(KIND_CELL) >= 0
    if (with_l2a) then
       ierr = PIO_def_var(file, 'rst:lnd2atm', PIO_double, &
                          (/dimid_for(NL2A), kdim_id(KIND_CELL)/), vid_l2a)
       if (ierr /= PIO_NOERR) call shr_sys_abort(subname//'ERROR: def_var rst:lnd2atm')
    end if

    ierr = PIO_def_var(file, 'rst:nstep'           , PIO_int, vid_scalar(1))
    ierr = PIO_def_var(file, 'rst:hist_nacs'       , PIO_int, vid_scalar(2))
    ierr = PIO_def_var(file, 'rst:hist_start_year' , PIO_int, vid_scalar(3))
    ierr = PIO_def_var(file, 'rst:hist_start_month', PIO_int, vid_scalar(4))
    ierr = PIO_def_var(file, 'rst:hist_start_day'  , PIO_int, vid_scalar(5))
    ierr = PIO_put_att(file, PIO_GLOBAL, 'caseid', trim(caseid))
    ierr = PIO_put_att(file, PIO_GLOBAL, 'title', 'kELM restart: all persistent state')
    ierr = PIO_put_att(file, PIO_GLOBAL, 'npes_written', npes)
    ierr = PIO_enddef(file)
    if (ierr /= PIO_NOERR) call shr_sys_abort(subname//'ERROR: enddef '//trim(fname))

    nwritten = 0
    do i = 1, nfld
       if (.not. f(i)%active) cycle
       allocate(buf(max(f(i)%n1 * max(f(i)%n2g, 1), 1)))
       buf = 0.0_c_double
       if (f(i)%n1 > 0) then
          call kELMRestartFieldGet(elm, int(i-1, c_int), buf, &
               int(f(i)%n1 * max(f(i)%n2, 1), c_int), st)
          if (st /= KELM_SUCCESS) call shr_sys_abort(subname//'ERROR: get '//trim(f(i)%name))
       end if
       call PIO_write_darray(file, vid(i), dec_iodesc(decomp_for(f(i)%kind, f(i)%n2g)), buf, ierr)
       if (ierr /= PIO_NOERR) call shr_sys_abort(subname//'ERROR: write '//trim(f(i)%name))
       deallocate(buf)
       nwritten = nwritten + 1
    end do

    if (with_l2a) then
       allocate(tab(ncell, NL2A), buf(max(ncell*NL2A, 1)))
       call kelm_lnd2atm_cell_table(elm, tab)
       buf = 0.0_c_double
       do e = 1, ncell
          do k = 1, NL2A
             buf((e-1)*NL2A + k) = tab(e, k)
          end do
       end do
       call PIO_write_darray(file, vid_l2a, dec_iodesc(decomp_for(KIND_CELL, NL2A)), buf, ierr)
       if (ierr /= PIO_NOERR) call shr_sys_abort(subname//'ERROR: write rst:lnd2atm')
       deallocate(tab, buf)
    end if

    call kelm_hist_get_restart_state(nacs_f, hy, hm, hd)
    nacs_c = kELMHistoryGetNacs(elm, st)
    if (nacs_c /= nacs_f) call shr_sys_abort(subname//'ERROR: C++ and Fortran '// &
         'history sample counts disagree')
    ierr = PIO_put_var(file, vid_scalar(1), nstep)
    ierr = PIO_put_var(file, vid_scalar(2), nacs_f)
    ierr = PIO_put_var(file, vid_scalar(3), hy)
    ierr = PIO_put_var(file, vid_scalar(4), hm)
    ierr = PIO_put_var(file, vid_scalar(5), hd)
    call PIO_closefile(file)
    call layout_end()

    call write_rpointer(inst_suffix, fname)

    if (masterproc) then
       write(logunit,'(a,a,a,i0,a,i0,a,i0)') subname, trim(fname), ': ', &
            nwritten, ' fields, nstep ', nstep, ', npes ', npes
       write(logunit,*) 'Successfully wrote local restart file ',trim(fname)
       write(logunit,'(72a1)') ("-",i=1,60)
       write(logunit,*)
       write(logunit,*) 'Successfully wrote local restart pointer file'
       write(logunit,*) 'Successfully wrote out restart data at nstep = ',nstep
       write(logunit,'(72a1)') ("-",i=1,60)
       call shr_sys_flush(logunit)
    end if

  contains

    integer function dimid_for(ext)
      integer, intent(in) :: ext
      integer :: k
      character(len=16) :: dname
      do k = 1, ndims
         if (dim_ext(k) == ext) then
            dimid_for = dim_id(k); return
         end if
      end do
      if (ndims == size(dim_ext)) call shr_sys_abort(subname//'ERROR: too many distinct extents')
      ndims = ndims + 1
      write(dname,'(a,i0)') 'n', ext
      ierr = PIO_def_dim(file, trim(dname), ext, dim_id(ndims))
      if (ierr /= PIO_NOERR) call shr_sys_abort(subname//'ERROR: def_dim '//trim(dname))
      dim_ext(ndims) = ext
      dimid_for = dim_id(ndims)
    end function dimid_for

  end subroutine kelm_rest_write

  !-----------------------------------------------------------------------
  subroutine kelm_rest_read(elm, fname, branch, ncell, ncell_global, nstep, logunit)
    !
    ! Restore every allocated, non-TOPO registry field and the Fortran-side
    ! state from fname. Returns the step counter of the last step the
    ! writing run took. A branch leaves history as kelm_hist_init set it
    ! (see the header). Collective.
    !
    type(kELMType) , intent(in)  :: elm
    character(len=*), intent(in)  :: fname
    logical         , intent(in)  :: branch
    integer         , intent(in)  :: ncell, ncell_global
    integer         , intent(out) :: nstep
    integer         , intent(in)  :: logunit

    type(file_desc_t) :: file
    type(var_desc_t)  :: vdesc
    type(field_t), allocatable :: f(:)
    integer :: nfld, i, ierr, varid, ndims_v, dimids(2), len1, len2, nread
    integer :: nacs, hy, hm, hd, nglob
    integer(c_int) :: st
    real(c_double), allocatable :: buf(:)
    character(len=*), parameter :: subname = '(kelm_rest_read) '

    call layout_begin(elm, f, nfld, ncell, ncell_global, branch, subname)

    ierr = PIO_openfile(pio_subsystem, file, io_type, trim(fname), PIO_NOWRITE)
    if (ierr /= PIO_NOERR) call shr_sys_abort(subname//'ERROR: cannot open '//trim(fname))

    nread = 0
    do i = 1, nfld
       if (.not. f(i)%active) cycle
       ierr = PIO_inq_varid(file, trim(f(i)%name), varid)
       if (ierr /= PIO_NOERR) call shr_sys_abort(subname//'ERROR: '//trim(f(i)%name)// &
            ' is not on '//trim(fname)//' -- a restart must carry all state')
       ierr = PIO_inq_varndims(file, varid, ndims_v)
       ierr = PIO_inq_vardimid(file, varid, dimids(1:ndims_v))
       ierr = PIO_inq_dimlen(file, dimids(1), len1)
       len2 = 0
       if (ndims_v == 2) ierr = PIO_inq_dimlen(file, dimids(2), len2)
       ! On file (Fortran order): rank 1 is (nglob); rank 2 is (n2, nglob).
       nglob = kmap(f(i)%kind)%nglob
       if ((f(i)%n2g == 0 .and. (ndims_v /= 1 .or. len1 /= nglob)) .or. &
           (f(i)%n2g >  0 .and. (ndims_v /= 2 .or. len1 /= f(i)%n2g .or. len2 /= nglob))) then
          call shr_sys_abort(subname//'ERROR: '//trim(f(i)%name)//' has a different '// &
               'shape on '//trim(fname)//' than in this configuration')
       end if
       ierr = PIO_inq_varid(file, trim(f(i)%name), vdesc)
       allocate(buf(max(f(i)%n1 * max(f(i)%n2g, 1), 1)))
       call PIO_read_darray(file, vdesc, dec_iodesc(decomp_for(f(i)%kind, f(i)%n2g)), buf, ierr)
       if (ierr /= PIO_NOERR) call shr_sys_abort(subname//'ERROR: read '//trim(f(i)%name))
       if (f(i)%n1 > 0) then
          call kELMRestartFieldSet(elm, int(i-1, c_int), buf, &
               int(f(i)%n1 * max(f(i)%n2, 1), c_int), st)
          if (st /= KELM_SUCCESS) call shr_sys_abort(subname//'ERROR: set '//trim(f(i)%name))
       end if
       deallocate(buf)
       nread = nread + 1
    end do

    nstep = read_int(file, 'rst:nstep')
    nacs  = read_int(file, 'rst:hist_nacs')
    hy    = read_int(file, 'rst:hist_start_year')
    hm    = read_int(file, 'rst:hist_start_month')
    hd    = read_int(file, 'rst:hist_start_day')
    call PIO_closefile(file)
    call layout_end()

    if (branch) then
       nacs = 0
    else
       call kelm_hist_set_restart_state(nacs, hy, hm, hd)
       call kELMHistorySetNacs(elm, int(nacs, c_int), st)
    end if

    if (masterproc) then
       write(logunit,'(a,a,a,i0,a,i0,a,i0,a,i0)') subname, trim(fname), ': restored ', &
            nread, ' fields, nstep ', nstep, ', history samples ', nacs, ', npes ', npes
       call shr_sys_flush(logunit)
    end if

  contains

    integer function read_int(fh, vname)
      type(file_desc_t), intent(inout) :: fh
      character(len=*) , intent(in)    :: vname
      integer :: vid_i, ival
      ierr = PIO_inq_varid(fh, vname, vid_i)
      if (ierr /= PIO_NOERR) call shr_sys_abort(subname//'ERROR: no '//vname//' on '//trim(fname))
      ierr = PIO_get_var(fh, vid_i, ival)
      read_int = ival
    end function read_int

  end subroutine kelm_rest_read

  !-----------------------------------------------------------------------
  subroutine kelm_rest_read_lnd2atm(elm, fname, ncell, ncell_global, tab, found)
    !
    ! The coupler export table a restart carries (rst:lnd2atm), by local
    ! cell, for a continue run's export at init: the state itself is read
    ! later, at the top of the first kelm_run. found is false for a restart
    ! written before the table existed. Collective.
    !
    type(kELMType) , intent(in)  :: elm
    character(len=*), intent(in)  :: fname
    integer         , intent(in)  :: ncell, ncell_global
    real(r8)        , intent(out) :: tab(:,:)   ! (ncell, NL2A)
    logical         , intent(out) :: found

    type(file_desc_t) :: file
    type(var_desc_t)  :: vdesc
    type(field_t), allocatable :: f(:)
    integer :: nfld, ierr, e, k
    real(c_double), allocatable :: buf(:)
    character(len=*), parameter :: subname = '(kelm_rest_read_lnd2atm) '

    tab = 0.0_r8
    call layout_begin(elm, f, nfld, ncell, ncell_global, .false., subname)
    ierr = PIO_openfile(pio_subsystem, file, io_type, trim(fname), PIO_NOWRITE)
    if (ierr /= PIO_NOERR) call shr_sys_abort(subname//'ERROR: cannot open '//trim(fname))
    call PIO_seterrorhandling(file, PIO_BCAST_ERROR)
    ierr = PIO_inq_varid(file, 'rst:lnd2atm', vdesc)
    call PIO_seterrorhandling(file, PIO_INTERNAL_ERROR)
    found = (ierr == PIO_NOERR)
    if (found) then
       if (.not. kmap(KIND_CELL)%built) call build_kindmap(KIND_CELL, ncell, ncell, ncell_global)
       allocate(buf(max(ncell*NL2A, 1)))
       call PIO_read_darray(file, vdesc, dec_iodesc(decomp_for(KIND_CELL, NL2A)), buf, ierr)
       if (ierr /= PIO_NOERR) call shr_sys_abort(subname//'ERROR: read rst:lnd2atm')
       do e = 1, ncell
          do k = 1, NL2A
             tab(e, k) = buf((e-1)*NL2A + k)
          end do
       end do
       deallocate(buf)
    end if
    call PIO_closefile(file)
    call layout_end()

  end subroutine kelm_rest_read_lnd2atm

  !-----------------------------------------------------------------------
  subroutine layout_begin(elm, f, nfld, ncell, ncell_global, skip_hist, subname)
    !
    ! Query every registry field, agree its layout across ranks, and build
    ! the global positions of each kind in use. skip_hist leaves the "hist:"
    ! accumulators inactive (a branch).
    !
    type(kELMType) , intent(in)  :: elm
    type(field_t), allocatable, intent(out) :: f(:)
    integer         , intent(out) :: nfld
    integer         , intent(in)  :: ncell, ncell_global
    logical         , intent(in)  :: skip_hist
    character(len=*), intent(in)  :: subname

    integer(c_int) :: cnt, st, a, b, k
    integer :: i, ierr, nloc
    integer, allocatable :: n1loc(:), n1sum(:), n2loc(:), n2max(:)

    call kELMRestartFieldCount(elm, cnt, st)
    if (st /= KELM_SUCCESS) call shr_sys_abort(subname//'ERROR: RestartFieldCount')
    nfld = cnt
    allocate(f(nfld), n1loc(nfld), n1sum(nfld), n2loc(nfld), n2max(nfld))
    do i = 1, nfld
       call kELMRestartFieldInfo(elm, int(i-1, c_int), f(i)%name, a, b, k, st)
       if (st /= KELM_SUCCESS) call shr_sys_abort(subname//'ERROR: RestartFieldInfo')
       f(i)%n1 = a; f(i)%n2 = b; f(i)%kind = k
       n1loc(i) = a; n2loc(i) = b
    end do

    ! The field list itself is identical on every rank (same registry, same
    ! history config); only the extents differ.
    call mpi_allreduce(n1loc, n1sum, nfld, MPI_INTEGER, MPI_SUM, mpicom_lnd, ierr)
    call mpi_allreduce(n2loc, n2max, nfld, MPI_INTEGER, MPI_MAX, mpicom_lnd, ierr)

    do i = 1, nfld
       f(i)%n2g = n2max(i)
       f(i)%active = n1sum(i) > 0 .and. f(i)%kind /= KIND_TOPO
       if (skip_hist .and. f(i)%name(1:5) == 'hist:') f(i)%active = .false.
       if (.not. f(i)%active) cycle
       ! A rank either holds the field for all its entities of that kind or
       ! has none of them; anything else would scramble the positions.
       nloc = local_count(f(i)%kind, ncell, f(i)%name, subname)
       if (f(i)%n1 /= nloc) then
          write(*,*) subname, trim(f(i)%name), ' n1 ', f(i)%n1, ' entities ', nloc, ' rank ', iam
          call shr_sys_abort(subname//'ERROR: '//trim(f(i)%name)//' extent disagrees '// &
               'with its kind ('//trim(kind_name(f(i)%kind))//')')
       end if
       if (f(i)%n1 > 0 .and. f(i)%n2 /= f(i)%n2g) call shr_sys_abort(subname// &
            'ERROR: '//trim(f(i)%name)//' second extent differs across ranks')
       if (.not. kmap(f(i)%kind)%built) &
            call build_kindmap(f(i)%kind, nloc, ncell, ncell_global)
    end do

  end subroutine layout_begin

  !-----------------------------------------------------------------------
  subroutine layout_end()
    integer :: k
    do k = 1, ndec
       call PIO_freedecomp(pio_subsystem, dec_iodesc(k))
    end do
    ndec = 0
    do k = 0, NKIND
       kmap(k)%built = .false.
       kmap(k)%n = 0
       kmap(k)%nglob = 0
       if (allocated(kmap(k)%pos)) deallocate(kmap(k)%pos)
    end do
  end subroutine layout_end

  !-----------------------------------------------------------------------
  integer function local_count(kind, ncell, name, subname)
    integer, intent(in) :: kind, ncell
    character(len=*), intent(in) :: name, subname
    select case (kind)
    case (KIND_CELL);     local_count = ncell
    case (KIND_NATCOL);   local_count = n_kokkos_col
    case (KIND_NATPATCH); local_count = n_kokkos_patch
    case (KIND_URBLUN);   local_count = n_kokkos_urb
    case (KIND_LAKECOL);  local_count = n_lake
    case (KIND_LAKEPATCH); local_count = n_lake   ! one patch per lake column
    case default
       local_count = 0
       call shr_sys_abort(subname//'ERROR: '//trim(name)//' is allocated, but kind '// &
            trim(kind_name(kind))//' has no entity-to-cell map for restarts yet')
    end select
  end function local_count

  !-----------------------------------------------------------------------
  subroutine build_kindmap(kind, n, ncell, ncell_global)
    !
    ! pos(e) for this rank's entities of one kind: entities ordered by global
    ! cell ordinal, then by occurrence within the cell. Rank iam's local
    ! cell g is global ordinal iam+1+(g-1)*npes (the round-robin
    ! decomposition in kelm_init; hist_dof uses the same map).
    !
    integer, intent(in) :: kind, n, ncell, ncell_global
    integer, allocatable :: cell(:), cnt(:), cntg(:), off(:), occ(:)
    integer :: e, g, gg, ierr
    character(len=*), parameter :: subname = '(kelm_rest build_kindmap) '

    allocate(cell(max(n,1)), cnt(ncell_global), cntg(ncell_global), &
             off(ncell_global), occ(max(ncell,1)))
    do e = 1, n
       select case (kind)
       case (KIND_CELL);     cell(e) = e
       case (KIND_NATCOL);   cell(e) = lun_gridcell(col_landunit(col_of_kcol(e)))
       case (KIND_NATPATCH); cell(e) = lun_gridcell(col_landunit(patch_column(patch_of_kpatch(e))))
       case (KIND_URBLUN);   cell(e) = lun_gridcell(lun_of_kurb(e))
       case (KIND_LAKECOL, KIND_LAKEPATCH); cell(e) = cell_of_klake(e)
       end select
       if (cell(e) < 1 .or. cell(e) > ncell) call shr_sys_abort(subname// &
            'ERROR: entity maps outside this rank''s cells')
    end do

    cnt = 0
    do e = 1, n
       gg = iam + 1 + (cell(e) - 1) * npes
       cnt(gg) = cnt(gg) + 1
    end do
    call mpi_allreduce(cnt, cntg, ncell_global, MPI_INTEGER, MPI_SUM, mpicom_lnd, ierr)
    off(1) = 0
    do gg = 2, ncell_global
       off(gg) = off(gg-1) + cntg(gg-1)
    end do

    allocate(kmap(kind)%pos(max(n,1)))
    kmap(kind)%pos = 0
    occ = 0
    do e = 1, n
       g  = cell(e)
       gg = iam + 1 + (g - 1) * npes
       occ(g) = occ(g) + 1
       kmap(kind)%pos(e) = off(gg) + occ(g)
    end do
    kmap(kind)%n = n
    kmap(kind)%nglob = sum(cntg)
    kmap(kind)%built = .true.
  end subroutine build_kindmap

  !-----------------------------------------------------------------------
  integer function decomp_for(kind, n2g)
    !
    ! Index in dec_iodesc of the PIO decomposition of a (kind, n2g) field,
    ! built on first use. The
    ! local buffer is the registry's row-major one: element (e, k) of a
    ! rank-2 field at (e-1)*m + k, m = max(n2g, 1); its global offset is
    ! (pos(e)-1)*m + k. A rank with no entities passes one dof of 0 (PIO's
    ! "no data").
    !
    integer, intent(in) :: kind, n2g
    integer, allocatable :: dof(:)
    integer :: d, e, k, m, n
    character(len=*), parameter :: subname = '(kelm_rest decomp_for) '

    do d = 1, ndec
       if (dec_kind(d) == kind .and. dec_n2(d) == n2g) then
          decomp_for = d; return
       end if
    end do
    if (ndec == MAXDEC) call shr_sys_abort(subname//'ERROR: too many decompositions')

    m = max(n2g, 1)
    n = kmap(kind)%n
    allocate(dof(max(n*m, 1)))
    dof = 0
    do e = 1, n
       do k = 1, m
          dof((e-1)*m + k) = (kmap(kind)%pos(e) - 1) * m + k
       end do
    end do

    ndec = ndec + 1
    dec_kind(ndec) = kind
    dec_n2(ndec) = n2g
    if (n2g > 0) then
       call PIO_initdecomp(pio_subsystem, PIO_double, (/n2g, kmap(kind)%nglob/), dof, &
                           dec_iodesc(ndec))
    else
       call PIO_initdecomp(pio_subsystem, PIO_double, (/kmap(kind)%nglob/), dof, &
                           dec_iodesc(ndec))
    end if
    decomp_for = ndec
  end function decomp_for

  !-----------------------------------------------------------------------
  integer function kelm_rest_peek_nstep(fname)
    !
    ! The step counter alone, for the driver: doalb's nstep-0/1 special cases
    ! must already see the restored count on the first step of a continued
    ! run, before the full read happens inside that step.
    !
    character(len=*), intent(in) :: fname
    type(file_desc_t) :: file
    integer :: ierr, vid_i, ival
    ierr = PIO_openfile(pio_subsystem, file, io_type, trim(fname), PIO_NOWRITE)
    if (ierr /= PIO_NOERR) call shr_sys_abort('(kelm_rest_peek_nstep) ERROR: cannot open '//trim(fname))
    ierr = PIO_inq_varid(file, 'rst:nstep', vid_i)
    if (ierr /= PIO_NOERR) call shr_sys_abort('(kelm_rest_peek_nstep) ERROR: no rst:nstep on '//trim(fname))
    ierr = PIO_get_var(file, vid_i, ival)
    call PIO_closefile(file)
    kelm_rest_peek_nstep = ival
  end function kelm_rest_peek_nstep

  !-----------------------------------------------------------------------
  function kelm_rest_rpointer_read(inst_suffix) result(fname)
    character(len=*), intent(in) :: inst_suffix
    character(len=CL) :: fname
    integer :: u, ios
    u = shr_file_getunit()
    open(u, file='rpointer.lnd'//trim(inst_suffix), status='old', action='read', iostat=ios)
    if (ios /= 0) call shr_sys_abort('(kelm_rest_rpointer_read) ERROR: cannot open rpointer.lnd'// &
         trim(inst_suffix)//' -- a continue run needs the previous run''s restart pointer')
    read(u,'(a)') fname
    close(u)
    call shr_file_freeunit(u)
    fname = adjustl(fname)
  end function kelm_rest_rpointer_read

  !-----------------------------------------------------------------------
  subroutine write_rpointer(inst_suffix, fname)
    character(len=*), intent(in) :: inst_suffix, fname
    integer :: u
    if (.not. masterproc) return
    u = shr_file_getunit()
    open(u, file='rpointer.lnd'//trim(inst_suffix), status='replace', action='write')
    write(u,'(a)') trim(fname)
    close(u)
    call shr_file_freeunit(u)
  end subroutine write_rpointer

end module kelmRestMod
