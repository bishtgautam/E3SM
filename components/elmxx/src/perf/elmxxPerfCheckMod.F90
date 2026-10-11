module elmxxPerfCheckMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Answer checks for the perf driver (plan Stage 6.8, gates B2-B4).
  ! ELMxx-only: no ELM counterpart.
  !
  ! In perf mode every virtual cell is a copy of a source cell and gets the
  ! same forcing, so every copy must carry bit-identical state. The checks
  ! read ALL persistent state through the restart registry (the same field
  ! list elmxxRestMod writes), so nothing has to be chosen by hand.
  !
  ! An entity (column, patch, urban landunit, ...) is identified by its cell
  ! and its occurrence within that cell, exactly as elmxxRestMod orders a
  ! restart: within one cell the packed order is the subgrid's, whatever
  ! rank owns the cell and however many cells there are.
  !
  ! Three checks:
  !  - elmxx_perf_check_replicas: on each rank, every copy equals the first
  !    local copy of its source, bit for bit; across ranks, a per-source hash
  !    of that first copy agrees (MPI min == max). Needs nothing external.
  !  - elmxx_perf_write_reference: the state of virtual cells 1..nsrc (one
  !    copy of each source) to a file.
  !  - elmxx_perf_compare_reference: every copy against such a file, bit for
  !    bit (tol = 0: same build, any N, any rank count) or within a relative
  !    tolerance (tol > 0: CPU against GPU, which cannot be BFB).
  !-----------------------------------------------------------------------

  use, intrinsic :: iso_fortran_env, only : int64
  use shr_kind_mod        , only : r8 => shr_kind_r8
  use shr_sys_mod         , only : shr_sys_abort, shr_sys_flush
  use elmxxSpmdMod        , only : masterproc, iam, npes, mpicom_lnd
  use elmxxSubgridMod     , only : lun_gridcell, col_landunit, patch_column
  use elmxxLakeMod        , only : n_lake, cell_of_klake
  use elmxxKokkosStateMod , only : n_kokkos_col, n_kokkos_patch, n_kokkos_urb, &
                                   col_of_kcol, patch_of_kpatch, lun_of_kurb
  use elmxx_mod           , only : ELMxxType, ELMXX_SUCCESS, &
                                   ELMxxRestartFieldCount, ELMxxRestartFieldInfo, &
                                   ELMxxRestartFieldGet
  use iso_c_binding       , only : c_int, c_double

  implicit none
  private
#include <mpif.h>

  public :: elmxx_perf_check_replicas
  public :: elmxx_perf_write_reference
  public :: elmxx_perf_compare_reference

  ! ELMxxRestartKind (ELMxx.h); the same constants elmxxRestMod uses.
  integer, parameter :: KIND_TOPO = 0, KIND_CELL = 1, KIND_NATCOL = 2, &
                        KIND_NATPATCH = 3, KIND_URBLUN = 4, &
                        KIND_LAKECOL = 7, KIND_LAKEPATCH = 8

  integer, parameter :: NAMELEN = 128
  character(len=8), parameter :: REF_MAGIC = 'ELMXXPRF'
  integer, parameter :: REF_VERSION = 1

  ! Report at most this many offending fields by name.
  integer, parameter :: MAX_NAMED = 10

  ! One registry field, laid out by entity.
  type :: field_layout
     character(len=NAMELEN) :: name = ' '
     integer :: kind = KIND_TOPO
     integer :: n = 0                 ! local entities
     integer :: m = 1                 ! values per entity (max(n2,1))
     integer, allocatable :: cell(:)  ! (n) local cell of each entity
     integer, allocatable :: occ(:)   ! (n) occurrence within its cell, 1-based
     integer, allocatable :: off(:)   ! (ncell+1) entities of cell c are list(off(c)+1:off(c+1))
     integer, allocatable :: list(:)  ! (n) entity indices grouped by cell, in occurrence order
  end type field_layout

contains

  !-----------------------------------------------------------------------
  logical function field_wanted(elm, i, f, ncell)
    !
    ! Info for registry field i, and its entity layout if it is checked.
    ! Skipped: unallocated views, TOPO (rank-local indices, rebuilt at init),
    ! history accumulators (history is off in perf mode), and any kind with
    ! no entity-to-cell map (none is allocated on the coupled path today).
    !
    type(ELMxxType), intent(in) :: elm
    integer, intent(in) :: i, ncell
    type(field_layout), intent(inout) :: f
    integer(c_int) :: n1, n2, kind, st
    integer :: e, c, ierr, mloc(2), mglob(2)

    field_wanted = .false.
    call ELMxxRestartFieldInfo(elm, int(i, c_int), f%name, n1, n2, kind, st)
    if (st /= ELMXX_SUCCESS) call shr_sys_abort('(elmxx_perf_check) ERROR: RestartFieldInfo failed')
    if (f%name(1:5) == 'hist:') return
    select case (kind)
    case (KIND_CELL, KIND_NATCOL, KIND_NATPATCH, KIND_URBLUN, KIND_LAKECOL, KIND_LAKEPATCH)
    case default
       return
    end select
    ! Collective from here: a rank whose cells have none of this kind (a
    ! lake-free source on its own rank) has n1 = 0 but must still join every
    ! reduction the others make. Checked if allocated anywhere; the level
    ! extent comes from the ranks that have it.
    mloc(1) = max(int(n1), 0)
    mloc(2) = max(int(n2), 1)
    if (n1 <= 0) mloc(2) = 0
    call mpi_allreduce(mloc, mglob, 2, MPI_INTEGER, MPI_MAX, mpicom_lnd, ierr)
    if (mglob(1) <= 0) return

    f%kind = kind
    f%n = mloc(1)
    f%m = mglob(2)
    if (allocated(f%cell)) deallocate(f%cell, f%occ, f%off, f%list)
    allocate(f%cell(f%n), f%occ(f%n), f%off(ncell+1), f%list(f%n))

    do e = 1, f%n
       select case (kind)
       case (KIND_CELL);      c = e
       case (KIND_NATCOL);    c = lun_gridcell(col_landunit(col_of_kcol(e)))
       case (KIND_NATPATCH);  c = lun_gridcell(col_landunit(patch_column(patch_of_kpatch(e))))
       case (KIND_URBLUN);    c = lun_gridcell(lun_of_kurb(e))
       case (KIND_LAKECOL, KIND_LAKEPATCH); c = cell_of_klake(e)
       end select
       if (c < 1 .or. c > ncell) call shr_sys_abort('(elmxx_perf_check) ERROR: '// &
            trim(f%name)//' maps an entity outside this rank''s cells')
       f%cell(e) = c
    end do

    ! Counting sort by cell; occurrence = position within the cell.
    f%off = 0
    do e = 1, f%n
       f%off(f%cell(e)+1) = f%off(f%cell(e)+1) + 1
    end do
    do c = 1, ncell
       f%off(c+1) = f%off(c+1) + f%off(c)
    end do
    block
      integer, allocatable :: placed(:)
      allocate(placed(ncell))
      placed = 0
      do e = 1, f%n
         c = f%cell(e)
         placed(c) = placed(c) + 1
         f%occ(e) = placed(c)
         f%list(f%off(c) + placed(c)) = e
      end do
    end block

    field_wanted = .true.
  end function field_wanted

  !-----------------------------------------------------------------------
  subroutine get_values(elm, i, f, buf)
    type(ELMxxType), intent(in) :: elm
    integer, intent(in) :: i
    type(field_layout), intent(in) :: f
    real(c_double), intent(inout) :: buf(:)
    integer(c_int) :: st
    call ELMxxRestartFieldGet(elm, int(i, c_int), buf, int(f%n * f%m, c_int), st)
    if (st /= ELMXX_SUCCESS) call shr_sys_abort('(elmxx_perf_check) ERROR: RestartFieldGet failed for '// &
         trim(f%name))
  end subroutine get_values

  !-----------------------------------------------------------------------
  pure integer(int64) function bits(x)
    real(r8), intent(in) :: x
    bits = transfer(x, 0_int64)
  end function bits

  !-----------------------------------------------------------------------
  subroutine canonical_cells(ncell, src_owned, nsrc, canon)
    !
    ! canon(s): the first local cell that copies source s, or 0 if none here.
    !
    integer, intent(in)  :: ncell, nsrc, src_owned(:)
    integer, intent(out) :: canon(:)
    integer :: c
    canon(1:nsrc) = 0
    do c = 1, ncell
       if (canon(src_owned(c)) == 0) canon(src_owned(c)) = c
    end do
  end subroutine canonical_cells

  !-----------------------------------------------------------------------
  subroutine elmxx_perf_check_replicas(elm, ncell, src_owned, nsrc, logunit, nbad_fields)
    !
    ! Gate B2 inside one run: every copy is bit-identical to the first local
    ! copy of its source, and that first copy hashes the same on every rank.
    !
    type(ELMxxType), intent(in) :: elm
    integer, intent(in)  :: ncell, nsrc, logunit
    integer, intent(in)  :: src_owned(:)
    integer, intent(out) :: nbad_fields   ! fields failing anywhere (all ranks)

    type(field_layout) :: f
    integer(c_int) :: nfld, st
    integer :: i, c, s, e, ce, j, k, ierr, nchecked, named
    integer :: canon(nsrc)
    integer(int64) :: nmis, nmis_g
    integer(int64), allocatable :: h(:), hmin(:), hmax(:)
    real(c_double), allocatable :: buf(:)
    logical :: bad

    call canonical_cells(ncell, src_owned, nsrc, canon)
    call ELMxxRestartFieldCount(elm, nfld, st)
    if (st /= ELMXX_SUCCESS) call shr_sys_abort('(elmxx_perf_check_replicas) ERROR: RestartFieldCount failed')

    allocate(h(nsrc), hmin(nsrc), hmax(nsrc))
    nbad_fields = 0
    nchecked = 0
    named = 0

    do i = 0, nfld - 1
       if (.not. field_wanted(elm, i, f, ncell)) cycle
       nchecked = nchecked + 1
       if (allocated(buf)) deallocate(buf)
       allocate(buf(max(f%n * f%m, 1)))
       if (f%n > 0) call get_values(elm, i, f, buf)

       ! Every local copy against its canonical copy.
       nmis = 0
       do e = 1, f%n
          c = f%cell(e)
          s = src_owned(c)
          if (c == canon(s)) cycle
          if (f%off(canon(s)+1) - f%off(canon(s)) /= f%off(c+1) - f%off(c)) then
             nmis = nmis + 1          ! different entity count: a layout bug
             cycle
          end if
          ce = f%list(f%off(canon(s)) + f%occ(e))
          do k = 1, f%m
             if (bits(buf((e-1)*f%m + k)) /= bits(buf((ce-1)*f%m + k))) nmis = nmis + 1
          end do
       end do
       call mpi_allreduce(nmis, nmis_g, 1, MPI_INTEGER8, MPI_SUM, mpicom_lnd, ierr)

       ! Canonical copies across ranks: hash each source's state (its entity
       ! count included), then min == max over the ranks that hold it.
       do s = 1, nsrc
          if (canon(s) == 0) then
             hmin(s) = huge(0_int64); hmax(s) = -huge(0_int64)
             cycle
          end if
          h(s) = int(f%off(canon(s)+1) - f%off(canon(s)), int64)
          do j = f%off(canon(s)) + 1, f%off(canon(s)+1)
             e = f%list(j)
             do k = 1, f%m
                h(s) = ieor(ishftc(h(s), 7), bits(buf((e-1)*f%m + k)))
             end do
          end do
          hmin(s) = h(s); hmax(s) = h(s)
       end do
       call mpi_allreduce(MPI_IN_PLACE, hmin, nsrc, MPI_INTEGER8, MPI_MIN, mpicom_lnd, ierr)
       call mpi_allreduce(MPI_IN_PLACE, hmax, nsrc, MPI_INTEGER8, MPI_MAX, mpicom_lnd, ierr)

       bad = nmis_g /= 0 .or. any(hmin /= hmax)
       if (bad) then
          nbad_fields = nbad_fields + 1
          if (masterproc .and. named < MAX_NAMED) then
             named = named + 1
             write(logunit,'(a,a,a,i0,a,l1)') ' (elmxx_perf_check_replicas) MISMATCH ', &
                  trim(f%name), ': copies differing in ', nmis_g, &
                  ' values; ranks disagree = ', any(hmin /= hmax)
          end if
       end if
    end do

    if (masterproc) then
       write(logunit,'(a,i0,a,i0,a,i0,a,i0,a)') ' (elmxx_perf_check_replicas) ', nchecked, &
            ' fields checked over ', ncell, ' local cells x ', npes, ' ranks: ', &
            nbad_fields, ' fail'
       call shr_sys_flush(logunit)
    end if
  end subroutine elmxx_perf_check_replicas

  !-----------------------------------------------------------------------
  subroutine elmxx_perf_write_reference(elm, ncell, src_owned, nsrc, fname, logunit)
    !
    ! Write one copy of each source (virtual cells 1..nsrc) to fname. Rank
    ! r owns virtual cell s when r = mod(s-1, npes), as its local cell
    ! (s-1)/npes + 1; that cell is also the first local copy of s. Values
    ! travel to rank 0 as bit patterns OR-reduced over zeros, which is exact
    ! (a sum would turn -0.0 into +0.0).
    !
    type(ELMxxType), intent(in) :: elm
    integer, intent(in) :: ncell, nsrc, logunit
    integer, intent(in) :: src_owned(:)
    character(len=*), intent(in) :: fname

    type(field_layout) :: f
    integer(c_int) :: nfld, st
    integer :: i, s, c, j, e, k, ierr, unitn, nwritten, nvals
    integer :: nent(nsrc), first(nsrc)
    integer(int64), allocatable :: pk(:)
    real(c_double), allocatable :: buf(:)

    if (nsrc > ncell * npes) call shr_sys_abort('(elmxx_perf_write_reference) ERROR: fewer cells than sources')

    call ELMxxRestartFieldCount(elm, nfld, st)
    if (masterproc) then
       open(newunit=unitn, file=trim(fname), access='stream', form='unformatted', &
            status='replace', action='write')
       write(unitn) REF_MAGIC, REF_VERSION, nsrc
    end if

    nwritten = 0
    do i = 0, nfld - 1
       if (.not. field_wanted(elm, i, f, ncell)) cycle
       if (allocated(buf)) deallocate(buf)
       allocate(buf(max(f%n * f%m, 1)))
       if (f%n > 0) call get_values(elm, i, f, buf)

       ! Entity count of each source, from its owner.
       nent = 0
       do s = 1, nsrc
          if (mod(s-1, npes) /= iam) cycle
          c = (s-1)/npes + 1
          nent(s) = f%off(c+1) - f%off(c)
       end do
       call mpi_allreduce(MPI_IN_PLACE, nent, nsrc, MPI_INTEGER, MPI_MAX, mpicom_lnd, ierr)
       first(1) = 0
       do s = 2, nsrc
          first(s) = first(s-1) + nent(s-1) * f%m
       end do
       nvals = sum(nent) * f%m

       allocate(pk(max(nvals,1)))
       pk = 0_int64
       do s = 1, nsrc
          if (mod(s-1, npes) /= iam) cycle
          c = (s-1)/npes + 1
          do j = 1, nent(s)
             e = f%list(f%off(c) + j)
             do k = 1, f%m
                pk(first(s) + (j-1)*f%m + k) = bits(buf((e-1)*f%m + k))
             end do
          end do
       end do
       if (nvals > 0) then
          if (masterproc) then
             call mpi_reduce(MPI_IN_PLACE, pk, nvals, MPI_INTEGER8, MPI_BOR, 0, mpicom_lnd, ierr)
          else
             call mpi_reduce(pk, pk, nvals, MPI_INTEGER8, MPI_BOR, 0, mpicom_lnd, ierr)
          end if
       end if

       if (masterproc) then
          write(unitn) f%name, f%kind, f%m, nent
          if (nvals > 0) write(unitn) pk(1:nvals)
       end if
       deallocate(pk)
       nwritten = nwritten + 1
    end do

    if (masterproc) then
       close(unitn)
       write(logunit,'(a,i0,a,i0,a,a)') ' (elmxx_perf_write_reference) ', nwritten, &
            ' fields x ', nsrc, ' source cells -> ', trim(fname)
       call shr_sys_flush(logunit)
    end if
  end subroutine elmxx_perf_write_reference

  !-----------------------------------------------------------------------
  subroutine elmxx_perf_compare_reference(elm, ncell, src_owned, nsrc, fname, tol, &
                                          logunit, nbad_fields)
    !
    ! Every local copy against the reference file. tol = 0: bit for bit.
    ! tol > 0: |a - b| <= tol * max(|a|, |b|), or both below tiny(1.0) in
    ! magnitude; the worst relative difference is reported per failing field.
    ! Fields present on only one side are counted as failures.
    !
    type(ELMxxType), intent(in) :: elm
    integer, intent(in)  :: ncell, nsrc, logunit
    integer, intent(in)  :: src_owned(:)
    character(len=*), intent(in) :: fname
    real(r8), intent(in) :: tol
    integer, intent(out) :: nbad_fields

    type ref_field
       character(len=NAMELEN) :: name
       integer :: kind, m
       integer, allocatable :: nent(:), first(:)
       integer(int64), allocatable :: pk(:)
       logical :: seen = .false.
    end type ref_field

    type(ref_field), allocatable :: rf(:), tmp(:)
    type(field_layout) :: f
    character(len=8) :: magic
    integer :: version, rnsrc, nref, ios, unitn, i, r, s, c, e, j, k, ierr, named, nchecked
    integer(c_int) :: nfld, st
    integer(int64) :: nmis, nmis_g
    real(r8) :: a, b, rel, relmax, relmax_g
    real(c_double), allocatable :: buf(:)
    logical :: bad

    ! ---- read the whole reference (small: nsrc cells) on every rank ----
    open(newunit=unitn, file=trim(fname), access='stream', form='unformatted', &
         status='old', action='read', iostat=ios)
    if (ios /= 0) call shr_sys_abort('(elmxx_perf_compare_reference) ERROR: cannot open '//trim(fname))
    read(unitn) magic, version, rnsrc
    if (magic /= REF_MAGIC .or. version /= REF_VERSION) &
         call shr_sys_abort('(elmxx_perf_compare_reference) ERROR: not a perf reference file: '//trim(fname))
    if (rnsrc /= nsrc) call shr_sys_abort('(elmxx_perf_compare_reference) ERROR: reference has a '// &
         'different number of source cells')
    allocate(rf(64))
    nref = 0
    do
       nref = nref + 1
       if (nref > size(rf)) then
          allocate(tmp(2*size(rf)))
          tmp(1:size(rf)) = rf
          call move_alloc(tmp, rf)
       end if
       allocate(rf(nref)%nent(nsrc), rf(nref)%first(nsrc))
       read(unitn, iostat=ios) rf(nref)%name, rf(nref)%kind, rf(nref)%m, rf(nref)%nent
       if (ios /= 0) then
          nref = nref - 1
          exit
       end if
       rf(nref)%first(1) = 0
       do s = 2, nsrc
          rf(nref)%first(s) = rf(nref)%first(s-1) + rf(nref)%nent(s-1) * rf(nref)%m
       end do
       allocate(rf(nref)%pk(max(sum(rf(nref)%nent) * rf(nref)%m, 1)))
       if (sum(rf(nref)%nent) > 0) read(unitn) rf(nref)%pk(1:sum(rf(nref)%nent) * rf(nref)%m)
    end do
    close(unitn)

    ! ---- compare ----
    call ELMxxRestartFieldCount(elm, nfld, st)
    nbad_fields = 0
    named = 0
    nchecked = 0
    do i = 0, nfld - 1
       if (.not. field_wanted(elm, i, f, ncell)) cycle
       nchecked = nchecked + 1
       r = 0
       do j = 1, nref
          if (rf(j)%name == f%name) then
             r = j; exit
          end if
       end do
       if (r == 0) then
          nbad_fields = nbad_fields + 1
          if (masterproc) write(logunit,*) '(elmxx_perf_compare_reference) MISSING in reference: ', trim(f%name)
          cycle
       end if
       rf(r)%seen = .true.
       if (rf(r)%m /= f%m .or. rf(r)%kind /= f%kind) then
          nbad_fields = nbad_fields + 1
          if (masterproc) write(logunit,*) '(elmxx_perf_compare_reference) SHAPE differs: ', trim(f%name)
          cycle
       end if

       if (allocated(buf)) deallocate(buf)
       allocate(buf(max(f%n * f%m, 1)))
       if (f%n > 0) call get_values(elm, i, f, buf)

       nmis = 0
       relmax = 0._r8
       do c = 1, ncell
          s = src_owned(c)
          if (f%off(c+1) - f%off(c) /= rf(r)%nent(s)) then
             nmis = nmis + 1
             cycle
          end if
          do j = 1, rf(r)%nent(s)
             e = f%list(f%off(c) + j)
             do k = 1, f%m
                a = buf((e-1)*f%m + k)
                b = transfer(rf(r)%pk(rf(r)%first(s) + (j-1)*f%m + k), 0._r8)
                if (tol <= 0._r8) then
                   if (bits(a) /= bits(b)) nmis = nmis + 1
                else
                   if (a /= a .or. b /= b) then
                      if ((a /= a) .neqv. (b /= b)) nmis = nmis + 1
                      cycle
                   end if
                   if (max(abs(a), abs(b)) < tiny(1._r8)) cycle
                   rel = abs(a - b) / max(abs(a), abs(b))
                   relmax = max(relmax, rel)
                   if (rel > tol) nmis = nmis + 1
                end if
             end do
          end do
       end do
       call mpi_allreduce(nmis, nmis_g, 1, MPI_INTEGER8, MPI_SUM, mpicom_lnd, ierr)
       call mpi_allreduce(relmax, relmax_g, 1, MPI_REAL8, MPI_MAX, mpicom_lnd, ierr)

       bad = nmis_g /= 0
       if (bad) then
          nbad_fields = nbad_fields + 1
          if (masterproc .and. named < MAX_NAMED) then
             named = named + 1
             write(logunit,'(a,a,a,i0,a,es10.3)') ' (elmxx_perf_compare_reference) MISMATCH ', &
                  trim(f%name), ': values differing ', nmis_g, '; worst rel diff ', relmax_g
          end if
       end if
    end do

    do j = 1, nref
       if (.not. rf(j)%seen) then
          nbad_fields = nbad_fields + 1
          if (masterproc) write(logunit,*) '(elmxx_perf_compare_reference) MISSING in this run: ', &
               trim(rf(j)%name)
       end if
    end do

    if (masterproc) then
       write(logunit,'(a,i0,a,a,a,es9.2,a,i0,a)') ' (elmxx_perf_compare_reference) ', nchecked, &
            ' fields against ', trim(fname), ' (tol ', tol, '): ', nbad_fields, ' fail'
       call shr_sys_flush(logunit)
    end if
  end subroutine elmxx_perf_compare_reference

end module elmxxPerfCheckMod
