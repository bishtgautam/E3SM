module elmxxSpmdMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! MPI information for ELMxx.
  !
  ! Modeled on components/rdycore/src/main/rdycoreSpmdMod.F90.
  !-----------------------------------------------------------------------

  implicit none
  save
  private

  logical, public :: masterproc  ! proc 0 logical for printing msgs
  integer, public :: iam         ! processor rank
  integer, public :: npes        ! number of processors
  integer, public :: mpicom_lnd  ! communicator group
  integer, public :: LNDID       ! mct compid

  ! Public methods
  public :: elmxxSpmdInit

#include <mpif.h>

contains

  !-----------------------------------------------------------------------
  subroutine elmxxSpmdInit(mpicom_local)
    !
    ! !DESCRIPTION:
    ! Set the ELMxx MPI communicator and determine rank and size.
    !
    implicit none
    !
    integer, intent(in) :: mpicom_local
    !
    integer :: ierr

    mpicom_lnd = mpicom_local

    ! Determine rank
    call mpi_comm_rank(mpicom_lnd, iam, ierr)
    masterproc = (iam == 0)

    ! Determine number of processors
    call mpi_comm_size(mpicom_lnd, npes, ierr)

  end subroutine elmxxSpmdInit

end module elmxxSpmdMod
