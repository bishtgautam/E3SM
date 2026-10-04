module elmxxRestMod
  !
  ! ELMxx restart files (Stage 6, step 2).
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
  ! registry names it ("natcol:t_grnd"), double precision, dims named by
  ! extent ("n15"). A rank-2 field (n1, n2) is stored with NetCDF dims
  ! (n1, n2) in C order, which is the registry's row-major buffer as-is.
  ! Integers are stored as doubles -- exact for I4.
  !
  ! READING IS STRICT. Every allocated registry field must be on the file
  ! with the same shape, or the run aborts: under "write all state" a field
  ! that silently fails to restore is a bug, not a tolerance.
  !
  ! SINGLE RANK for now. Fields are written whole with PIO_put_var from the
  ! one rank; a multi-rank decomposition per subgrid entity comes with the
  ! hist_dof work. elmxx_rest_write/read abort on npes > 1 rather than write
  ! a file that only holds rank 0.
  !
  use shr_kind_mod , only : r8 => shr_kind_r8, CL => shr_kind_cl
  use shr_sys_mod  , only : shr_sys_abort, shr_sys_flush
  use shr_file_mod , only : shr_file_getunit, shr_file_freeunit
  use elmxxSpmdMod , only : masterproc, npes
  use elmxxIO      , only : pio_subsystem, io_type
  use elmxxHistMod , only : elmxx_hist_get_restart_state, elmxx_hist_set_restart_state
  use elmxx_mod    , only : ELMxxType, ELMXX_SUCCESS, &
                            ELMxxRestartFieldCount, ELMxxRestartFieldInfo, &
                            ELMxxRestartFieldGet, ELMxxRestartFieldSet, &
                            ELMxxHistoryGetNacs, ELMxxHistorySetNacs
  use pio
  use iso_c_binding, only : c_int, c_double
  implicit none
  private

  public :: elmxx_rest_write
  public :: elmxx_rest_read
  public :: elmxx_rest_peek_nstep
  public :: elmxx_rest_rpointer_read

  integer, parameter :: NAMELEN = 128

contains

  !-----------------------------------------------------------------------
  subroutine elmxx_rest_write(elm, caseid, inst_suffix, yr, mon, day, tod, &
                              nstep, logunit)
    !
    ! Write <caseid>.elmxx.r.YYYY-MM-DD-SSSSS.nc, dated with the END of the
    ! step just taken (the time the next run starts from, as ELM names its
    ! restarts), and point rpointer.lnd<inst_suffix> at it.
    !
    type(ELMxxType) , intent(in) :: elm
    character(len=*), intent(in) :: caseid, inst_suffix
    integer         , intent(in) :: yr, mon, day, tod, nstep, logunit

    type(file_desc_t) :: file
    type(var_desc_t)  :: vid_scalar(5)
    type(var_desc_t), allocatable :: vid(:)
    integer, allocatable :: n1(:), n2(:)
    character(len=NAMELEN), allocatable :: names(:)
    integer :: nfld, i, ierr, d1, d2
    integer(c_int) :: st, cnt
    integer :: nacs_f, hy, hm, hd, nacs_c
    integer :: dim_ext(64), dim_id(64), ndims
    real(c_double), allocatable :: buf(:), a2(:,:)
    character(len=CL) :: fname
    character(len=*), parameter :: subname = '(elmxx_rest_write) '

    if (npes > 1) call shr_sys_abort(subname//'ERROR: restart is single-rank only so far')

    write(fname,'(a,".elmxx.r.",i4.4,"-",i2.2,"-",i2.2,"-",i5.5,".nc")') &
         trim(caseid), yr, mon, day, tod

    call ELMxxRestartFieldCount(elm, cnt, st)
    if (st /= ELMXX_SUCCESS) call shr_sys_abort(subname//'ERROR: RestartFieldCount')
    nfld = cnt
    allocate(names(nfld), n1(nfld), n2(nfld), vid(nfld))
    do i = 1, nfld
       call field_info(elm, i, names(i), n1(i), n2(i))
    end do

    ierr = PIO_createfile(pio_subsystem, file, io_type, trim(fname), PIO_CLOBBER)
    if (ierr /= PIO_NOERR) call shr_sys_abort(subname//'ERROR: cannot create '//trim(fname))

    ndims = 0
    do i = 1, nfld
       if (n1(i) == 0) cycle
       d1 = dimid_for(n1(i))
       if (n2(i) > 0) then
          d2 = dimid_for(n2(i))
          ! Fortran dim order is the reverse of C: (n2, n1) here is (n1, n2) on file.
          ierr = PIO_def_var(file, trim(names(i)), PIO_double, (/d2, d1/), vid(i))
       else
          ierr = PIO_def_var(file, trim(names(i)), PIO_double, (/d1/), vid(i))
       end if
       if (ierr /= PIO_NOERR) call shr_sys_abort(subname//'ERROR: def_var '//trim(names(i)))
    end do

    ierr = PIO_def_var(file, 'rst:nstep'           , PIO_int, vid_scalar(1))
    ierr = PIO_def_var(file, 'rst:hist_nacs'       , PIO_int, vid_scalar(2))
    ierr = PIO_def_var(file, 'rst:hist_start_year' , PIO_int, vid_scalar(3))
    ierr = PIO_def_var(file, 'rst:hist_start_month', PIO_int, vid_scalar(4))
    ierr = PIO_def_var(file, 'rst:hist_start_day'  , PIO_int, vid_scalar(5))
    ierr = PIO_put_att(file, PIO_GLOBAL, 'caseid', trim(caseid))
    ierr = PIO_put_att(file, PIO_GLOBAL, 'title', 'ELMxx restart: all persistent state')
    ierr = PIO_enddef(file)
    if (ierr /= PIO_NOERR) call shr_sys_abort(subname//'ERROR: enddef '//trim(fname))

    do i = 1, nfld
       if (n1(i) == 0) cycle
       allocate(buf(n1(i) * max(n2(i), 1)))
       call ELMxxRestartFieldGet(elm, int(i-1, c_int), buf, int(size(buf), c_int), st)
       if (st /= ELMXX_SUCCESS) call shr_sys_abort(subname//'ERROR: get '//trim(names(i)))
       if (n2(i) > 0) then
          allocate(a2(n2(i), n1(i)))
          a2 = reshape(buf, (/n2(i), n1(i)/))
          ierr = PIO_put_var(file, vid(i), a2)
          deallocate(a2)
       else
          ierr = PIO_put_var(file, vid(i), buf)
       end if
       if (ierr /= PIO_NOERR) call shr_sys_abort(subname//'ERROR: put_var '//trim(names(i)))
       deallocate(buf)
    end do

    call elmxx_hist_get_restart_state(nacs_f, hy, hm, hd)
    nacs_c = ELMxxHistoryGetNacs(elm, st)
    if (nacs_c /= nacs_f) call shr_sys_abort(subname//'ERROR: C++ and Fortran '// &
         'history sample counts disagree')
    ierr = PIO_put_var(file, vid_scalar(1), nstep)
    ierr = PIO_put_var(file, vid_scalar(2), nacs_f)
    ierr = PIO_put_var(file, vid_scalar(3), hy)
    ierr = PIO_put_var(file, vid_scalar(4), hm)
    ierr = PIO_put_var(file, vid_scalar(5), hd)
    call PIO_closefile(file)

    call write_rpointer(inst_suffix, fname)

    if (masterproc) then
       write(logunit,'(a,a,a,i0,a,i0,a)') subname, trim(fname), ': ', &
            count(n1 > 0), ' fields, nstep ', nstep, ''
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

  end subroutine elmxx_rest_write

  !-----------------------------------------------------------------------
  subroutine elmxx_rest_read(elm, fname, nstep, logunit)
    !
    ! Restore every allocated registry field and the Fortran-side state from
    ! fname. Returns the step counter of the last step the writing run took.
    !
    type(ELMxxType) , intent(in)  :: elm
    character(len=*), intent(in)  :: fname
    integer         , intent(out) :: nstep
    integer         , intent(in)  :: logunit

    type(file_desc_t) :: file
    integer :: nfld, i, ierr, varid, ndims_v, dimids(2), len1, len2
    integer :: n1, n2, nacs, hy, hm, hd
    integer(c_int) :: st, cnt
    character(len=NAMELEN) :: name
    real(c_double), allocatable :: buf(:), a2(:,:)
    character(len=*), parameter :: subname = '(elmxx_rest_read) '

    if (npes > 1) call shr_sys_abort(subname//'ERROR: restart is single-rank only so far')

    ierr = PIO_openfile(pio_subsystem, file, io_type, trim(fname), PIO_NOWRITE)
    if (ierr /= PIO_NOERR) call shr_sys_abort(subname//'ERROR: cannot open '//trim(fname))

    call ELMxxRestartFieldCount(elm, cnt, st)
    nfld = cnt
    do i = 1, nfld
       call field_info(elm, i, name, n1, n2)
       if (n1 == 0) cycle
       ierr = PIO_inq_varid(file, trim(name), varid)
       if (ierr /= PIO_NOERR) call shr_sys_abort(subname//'ERROR: '//trim(name)// &
            ' is not on '//trim(fname)//' -- a restart must carry all state')
       ierr = PIO_inq_varndims(file, varid, ndims_v)
       ierr = PIO_inq_vardimid(file, varid, dimids(1:ndims_v))
       ierr = PIO_inq_dimlen(file, dimids(1), len1)
       len2 = 0
       if (ndims_v == 2) ierr = PIO_inq_dimlen(file, dimids(2), len2)
       ! On file (Fortran order): rank 1 is (n1); rank 2 is (n2, n1).
       if ((n2 == 0 .and. (ndims_v /= 1 .or. len1 /= n1)) .or. &
           (n2 >  0 .and. (ndims_v /= 2 .or. len1 /= n2 .or. len2 /= n1))) then
          call shr_sys_abort(subname//'ERROR: '//trim(name)//' has a different '// &
               'shape on '//trim(fname)//' than in this configuration')
       end if
       allocate(buf(n1 * max(n2, 1)))
       if (n2 > 0) then
          allocate(a2(n2, n1))
          ierr = PIO_get_var(file, varid, a2)
          buf = reshape(a2, (/n1 * n2/))
          deallocate(a2)
       else
          ierr = PIO_get_var(file, varid, buf)
       end if
       if (ierr /= PIO_NOERR) call shr_sys_abort(subname//'ERROR: get_var '//trim(name))
       call ELMxxRestartFieldSet(elm, int(i-1, c_int), buf, int(size(buf), c_int), st)
       if (st /= ELMXX_SUCCESS) call shr_sys_abort(subname//'ERROR: set '//trim(name))
       deallocate(buf)
    end do

    nstep = read_int(file, 'rst:nstep')
    nacs  = read_int(file, 'rst:hist_nacs')
    hy    = read_int(file, 'rst:hist_start_year')
    hm    = read_int(file, 'rst:hist_start_month')
    hd    = read_int(file, 'rst:hist_start_day')
    call PIO_closefile(file)

    call elmxx_hist_set_restart_state(nacs, hy, hm, hd)
    call ELMxxHistorySetNacs(elm, int(nacs, c_int), st)

    if (masterproc) then
       write(logunit,'(a,a,a,i0,a,i0)') subname, trim(fname), ': restored, nstep ', &
            nstep, ', history samples ', nacs
       call shr_sys_flush(logunit)
    end if

  contains

    integer function read_int(f, vname)
      type(file_desc_t), intent(inout) :: f
      character(len=*) , intent(in)    :: vname
      integer :: vid_i, ival
      ierr = PIO_inq_varid(f, vname, vid_i)
      if (ierr /= PIO_NOERR) call shr_sys_abort(subname//'ERROR: no '//vname//' on '//trim(fname))
      ierr = PIO_get_var(f, vid_i, ival)
      read_int = ival
    end function read_int

  end subroutine elmxx_rest_read

  !-----------------------------------------------------------------------
  integer function elmxx_rest_peek_nstep(fname)
    !
    ! The step counter alone, for the driver: doalb's nstep-0/1 special cases
    ! must already see the restored count on the first step of a continued
    ! run, before the full read happens inside that step.
    !
    character(len=*), intent(in) :: fname
    type(file_desc_t) :: file
    integer :: ierr, vid_i, ival
    ierr = PIO_openfile(pio_subsystem, file, io_type, trim(fname), PIO_NOWRITE)
    if (ierr /= PIO_NOERR) call shr_sys_abort('(elmxx_rest_peek_nstep) ERROR: cannot open '//trim(fname))
    ierr = PIO_inq_varid(file, 'rst:nstep', vid_i)
    if (ierr /= PIO_NOERR) call shr_sys_abort('(elmxx_rest_peek_nstep) ERROR: no rst:nstep on '//trim(fname))
    ierr = PIO_get_var(file, vid_i, ival)
    call PIO_closefile(file)
    elmxx_rest_peek_nstep = ival
  end function elmxx_rest_peek_nstep

  !-----------------------------------------------------------------------
  function elmxx_rest_rpointer_read(inst_suffix) result(fname)
    character(len=*), intent(in) :: inst_suffix
    character(len=CL) :: fname
    integer :: u, ios
    u = shr_file_getunit()
    open(u, file='rpointer.lnd'//trim(inst_suffix), status='old', action='read', iostat=ios)
    if (ios /= 0) call shr_sys_abort('(elmxx_rest_rpointer_read) ERROR: cannot open rpointer.lnd'// &
         trim(inst_suffix)//' -- a continue run needs the previous run''s restart pointer')
    read(u,'(a)') fname
    close(u)
    call shr_file_freeunit(u)
    fname = adjustl(fname)
  end function elmxx_rest_rpointer_read

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

  !-----------------------------------------------------------------------
  subroutine field_info(elm, i, name, n1, n2)
    ! 1-based i over the registry's 0-based index.
    type(ELMxxType) , intent(in)  :: elm
    integer         , intent(in)  :: i
    character(len=*), intent(out) :: name
    integer         , intent(out) :: n1, n2
    integer(c_int) :: a, b, st
    call ELMxxRestartFieldInfo(elm, int(i-1, c_int), name, a, b, st)
    if (st /= ELMXX_SUCCESS) call shr_sys_abort('(elmxxRestMod) ERROR: RestartFieldInfo')
    n1 = a; n2 = b
  end subroutine field_info

end module elmxxRestMod
