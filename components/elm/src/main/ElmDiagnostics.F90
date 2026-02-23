module ElmDiagnostics
  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Diagnostics module for ELMPy development.
  ! Captures state before and after each major physics routine as
  ! unformatted binary output for Python validation.
  !
  ! Usage:
  !   - Enable via namelist: elm_diagnostics = .true.
  !   - Or call diag_init() directly with enabled=.true.
  !   - Call diag_new_timestep() at start of each timestep
  !   - Call diag_snapshot_*() before/after each physics routine
  !   - Call diag_finalize() at end of run
  !
  ! Binary format per record:
  !   [timestep (int)] [label (64 chars)] [ndims (int)] [dim1 (int)] [dim2 (int, if 2d)] [data (r8 array)]
  !
  ! Labels follow convention: "routinename_in:varname" or "routinename_out:varname"
  !
  ! Python reader: see elm_python/harness/diagnostic_reader.py
  !-----------------------------------------------------------------------

  use shr_kind_mod, only : r8 => shr_kind_r8
  implicit none
  private

  ! ---- Public control variables ----
  logical, public :: diag_enabled = .false.
  integer, public :: diag_timestep = 0

  ! ---- Private module state ----
  integer, parameter :: LABEL_LEN = 64
  integer            :: diag_unit = 199    ! Fortran unit number for output
  logical            :: diag_is_open = .false.

  ! ---- Public interface ----
  public :: diag_init
  public :: diag_finalize
  public :: diag_new_timestep
  public :: diag_snapshot_scalar
  public :: diag_snapshot_1d
  public :: diag_snapshot_2d
  public :: diag_snapshot_int_1d

contains

  !-----------------------------------------------------------------------
  subroutine diag_init(output_path, enabled)
    !
    ! !DESCRIPTION:
    ! Initialize the diagnostics system. Opens the binary output file.
    !
    ! !ARGUMENTS:
    character(len=*), intent(in) :: output_path
    logical,          intent(in) :: enabled
    !-----------------------------------------------------------------------

    diag_enabled = enabled
    diag_timestep = 0

    if (.not. diag_enabled) return

    open(unit=diag_unit, file=trim(output_path), &
         form='unformatted', access='stream', action='write', status='replace')
    diag_is_open = .true.

    ! Write a header so the Python reader can verify the file format
    ! Header: magic string "ELMDIAG1" (8 bytes) + LABEL_LEN (int)
    write(diag_unit) 'ELMDIAG1'
    write(diag_unit) LABEL_LEN

  end subroutine diag_init

  !-----------------------------------------------------------------------
  subroutine diag_finalize()
    !
    ! !DESCRIPTION:
    ! Close the diagnostics output file.
    !-----------------------------------------------------------------------

    if (diag_is_open) then
       close(diag_unit)
       diag_is_open = .false.
    end if
    diag_enabled = .false.

  end subroutine diag_finalize

  !-----------------------------------------------------------------------
  subroutine diag_new_timestep()
    !
    ! !DESCRIPTION:
    ! Increment the timestep counter. Call at the beginning of each driver timestep.
    !-----------------------------------------------------------------------

    if (.not. diag_enabled) return
    diag_timestep = diag_timestep + 1

  end subroutine diag_new_timestep

  !-----------------------------------------------------------------------
  subroutine diag_snapshot_scalar(label, value)
    !
    ! !DESCRIPTION:
    ! Write a single scalar value.
    !
    ! !ARGUMENTS:
    character(len=*), intent(in) :: label
    real(r8),         intent(in) :: value
    !
    ! !LOCAL VARIABLES:
    character(len=LABEL_LEN) :: padded_label
    integer :: ndims
    !-----------------------------------------------------------------------

    if (.not. diag_enabled) return

    padded_label = label
    ndims = 0

    write(diag_unit) diag_timestep
    write(diag_unit) padded_label
    write(diag_unit) ndims
    write(diag_unit) value

  end subroutine diag_snapshot_scalar

  !-----------------------------------------------------------------------
  subroutine diag_snapshot_1d(label, array, n)
    !
    ! !DESCRIPTION:
    ! Write a 1D real(r8) array snapshot.
    !
    ! !ARGUMENTS:
    character(len=*), intent(in) :: label
    real(r8),         intent(in) :: array(:)
    integer,          intent(in) :: n      ! number of elements to write
    !
    ! !LOCAL VARIABLES:
    character(len=LABEL_LEN) :: padded_label
    integer :: ndims
    !-----------------------------------------------------------------------

    if (.not. diag_enabled) return

    padded_label = label
    ndims = 1

    write(diag_unit) diag_timestep
    write(diag_unit) padded_label
    write(diag_unit) ndims
    write(diag_unit) n
    write(diag_unit) array(1:n)

  end subroutine diag_snapshot_1d

  !-----------------------------------------------------------------------
  subroutine diag_snapshot_2d(label, array, n1, n2)
    !
    ! !DESCRIPTION:
    ! Write a 2D real(r8) array snapshot.
    ! Fortran stores column-major, so array(1:n1, 1:n2) is written as-is.
    ! Python reader must account for Fortran ordering (use order='F' or transpose).
    !
    ! !ARGUMENTS:
    character(len=*), intent(in) :: label
    real(r8),         intent(in) :: array(:,:)
    integer,          intent(in) :: n1, n2   ! dimensions to write
    !
    ! !LOCAL VARIABLES:
    character(len=LABEL_LEN) :: padded_label
    integer :: ndims
    !-----------------------------------------------------------------------

    if (.not. diag_enabled) return

    padded_label = label
    ndims = 2

    write(diag_unit) diag_timestep
    write(diag_unit) padded_label
    write(diag_unit) ndims
    write(diag_unit) n1, n2
    write(diag_unit) array(1:n1, 1:n2)

  end subroutine diag_snapshot_2d

  !-----------------------------------------------------------------------
  subroutine diag_snapshot_int_1d(label, array, n)
    !
    ! !DESCRIPTION:
    ! Write a 1D integer array snapshot.
    !
    ! !ARGUMENTS:
    character(len=*), intent(in) :: label
    integer,          intent(in) :: array(:)
    integer,          intent(in) :: n      ! number of elements to write
    !
    ! !LOCAL VARIABLES:
    character(len=LABEL_LEN) :: padded_label
    integer :: ndims_tag
    !-----------------------------------------------------------------------

    if (.not. diag_enabled) return

    padded_label = label
    ndims_tag = -1   ! negative ndims signals integer data

    write(diag_unit) diag_timestep
    write(diag_unit) padded_label
    write(diag_unit) ndims_tag
    write(diag_unit) n
    write(diag_unit) array(1:n)

  end subroutine diag_snapshot_int_1d

end module ElmDiagnostics
