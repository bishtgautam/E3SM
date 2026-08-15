module elmxxMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! ELMxx driver-facing state and entry points.
  !
  ! Structured after components/rdycore/src/main/rdycoreMod.F90: this module owns
  ! the decomposition (num_cells_owned / num_cells_global / natural_id_cells_owned)
  ! and the init/run/final entry points, and the MCT coupling layer in
  ! src/cpl/lnd_comp_mct.F90 reads that state to build the gsMap and domain.
  !
  ! Unlike RDycore, there is no PETSc here, and the decomposition is not handed
  ! to us by a library -- for now ELMxx simply partitions the global land grid
  ! round-robin across ranks. Time stepping is a no-op until the Kokkos port is
  ! wired in.
  !-----------------------------------------------------------------------

  use shr_kind_mod , only : r8 => shr_kind_r8
  use shr_sys_mod  , only : shr_sys_abort, shr_sys_flush
  use shr_file_mod , only : shr_file_getunit, shr_file_freeunit
  use shr_nl_mod   , only : shr_nl_find_group_name
  use elmxxSpmdMod , only : masterproc, iam, npes, mpicom_lnd
  use elmxxIO      , only : elmxx_pio_init, elmxx_read_domain

  implicit none
  save
  private

#include <mpif.h>

  !--------------------------------------------------------------------------
  ! Decomposition -- read by src/cpl/lnd_comp_mct.F90
  !--------------------------------------------------------------------------
  integer , public          :: num_cells_owned           ! active land cells owned by this rank
  integer , public          :: num_cells_global          ! total active land cells (mask==1)
  integer , public, pointer :: natural_id_cells_owned(:) ! 1-based global grid IDs owned by this rank

  !--------------------------------------------------------------------------
  ! Global grid, read from the land domain file on every rank
  !--------------------------------------------------------------------------
  integer , public          :: nlon_g                    ! ni
  integer , public          :: nlat_g                    ! nj
  real(r8), public, pointer :: lonc_g(:)                 ! cell center longitude (deg)
  real(r8), public, pointer :: latc_g(:)                 ! cell center latitude  (deg)
  real(r8), public, pointer :: areac_g(:)                ! cell area (radians^2)
  real(r8), public, pointer :: maskc_g(:)                ! domain mask (0 or 1)
  real(r8), public, pointer :: fracc_g(:)                ! land fraction

  !--------------------------------------------------------------------------
  ! Namelist (&elmxx_inparm in lnd_in)
  !--------------------------------------------------------------------------
  logical           , public :: do_elmxx   = .true.
  character(len=256), public :: fatmlndfrc = ' '

  !--------------------------------------------------------------------------
  ! Instance information
  !--------------------------------------------------------------------------
  character(len=16), public :: inst_name
  character(len=16), public :: inst_suffix   ! e.g. "_0001" or ""
  integer          , public :: inst_index

  integer, public :: iulog = 6

  integer, private :: nstep = 0

  public :: elmxx_read_namelist
  public :: elmxx_init
  public :: elmxx_run
  public :: elmxx_final

contains

  !-----------------------------------------------------------------------
  subroutine elmxx_read_namelist(logunit)
    !
    ! !DESCRIPTION:
    ! Read &elmxx_inparm from lnd_in on the master task and broadcast it.
    !
    implicit none
    !
    integer, intent(in) :: logunit
    !
    character(len=256) :: nlfilename
    integer            :: ier, unitn
    logical            :: lexist
    character(len=*), parameter :: subname = '(elmxx_read_namelist) '

    namelist /elmxx_inparm/ do_elmxx, fatmlndfrc

    ! defaults
    do_elmxx   = .true.
    fatmlndfrc = ' '

    nlfilename = "lnd_in" // trim(inst_suffix)

    inquire (file = trim(nlfilename), exist = lexist)
    if ( .not. lexist ) then
       write(logunit,*) subname//' ERROR: namelist file does NOT exist: '//trim(nlfilename)
       call shr_sys_abort(subname//' ERROR: '//trim(nlfilename)//' does not exist')
    end if

    if (masterproc) then
       unitn = shr_file_getunit()
       write(logunit,*) 'Read in elmxx_inparm namelist from: ', trim(nlfilename)
       open( unitn, file=trim(nlfilename), status='old', iostat=ier )
       if (ier /= 0) call shr_sys_abort(subname//' ERROR opening '//trim(nlfilename))

       call shr_nl_find_group_name(unitn, 'elmxx_inparm', status=ier)
       if (ier /= 0) then
          call shr_sys_abort(subname//' ERROR: elmxx_inparm group not found in '//trim(nlfilename))
       end if

       read(unitn, elmxx_inparm, iostat=ier)
       if (ier /= 0) then
          call shr_sys_abort(subname//' ERROR reading elmxx_inparm from '//trim(nlfilename))
       end if

       close(unitn)
       call shr_file_freeunit(unitn)
    end if

    call mpi_bcast (do_elmxx  , 1                , MPI_LOGICAL  , 0, mpicom_lnd, ier)
    call mpi_bcast (fatmlndfrc, len(fatmlndfrc)  , MPI_CHARACTER, 0, mpicom_lnd, ier)

    if (masterproc) then
       write(logunit,*) ' '
       write(logunit,*) 'read from namelist:'
       write(logunit,*) '   do_elmxx   = ', do_elmxx
       write(logunit,*) '   fatmlndfrc = ', trim(fatmlndfrc)
       call shr_sys_flush(logunit)
    end if

  end subroutine elmxx_read_namelist

  !-----------------------------------------------------------------------
  subroutine elmxx_init(logunit)
    !
    ! !DESCRIPTION:
    ! Initialize ELMxx: read the land domain and build the round-robin
    ! decomposition of the global grid.
    !
    implicit none
    !
    integer, intent(in) :: logunit
    !
    integer :: i, n, k
    integer :: num_cells_grid                  ! ni*nj, including non-land cells
    integer, allocatable :: land_ids(:)        ! global grid IDs of the active land cells
    character(len=*), parameter :: subname = '(elmxx_init) '

    if (masterproc) then
       write(logunit,*) 'ELMxx model initialization'
       call shr_sys_flush(logunit)
    end if

    if (len_trim(fatmlndfrc) == 0) then
       call shr_sys_abort(subname//' ERROR: fatmlndfrc is not set in lnd_in')
    end if

    ! ---- global grid, read on every rank ----
    call elmxx_pio_init(inst_name)
    call elmxx_read_domain(logunit, fatmlndfrc, nlon_g, nlat_g, &
                           lonc_g, latc_g, areac_g, maskc_g, fracc_g)

    num_cells_grid = nlon_g * nlat_g

    ! ---- count total land gridcells ----
    ! Only cells with mask == 1 are active land, exactly as ELM counts numg in
    ! components/elm/src/main/decompInitMod.F90. Ocean cells carry mask == 0 and
    ! frac == 0 and take no part in the land decomposition.
    num_cells_global = 0
    do n = 1, num_cells_grid
       if (nint(maskc_g(n)) == 1) num_cells_global = num_cells_global + 1
    end do

    if (num_cells_global == 0) then
       call shr_sys_abort(subname//' ERROR: no active land cells (mask==1) in '//trim(fatmlndfrc))
    end if

    if (npes > num_cells_global) then
       write(logunit,*) subname,'ERROR: number of processes exceeds number of ', &
                        'land grid cells ',npes,num_cells_global
       call shr_sys_abort(subname//' ERROR: more MPI ranks than active land cells')
    end if

    ! global grid IDs of the active land cells, in ascending order
    allocate(land_ids(num_cells_global))
    k = 0
    do n = 1, num_cells_grid
       if (nint(maskc_g(n)) == 1) then
          k = k + 1
          land_ids(k) = n
       end if
    end do

    ! ---- round-robin decomposition over the active land cells ----
    ! Rank iam owns the iam+1, iam+1+npes, iam+1+2*npes, ... land cells. Note the
    ! IDs stored are global *grid* IDs (indices into the ni*nj arrays), which is
    ! what the gsMap and the domain both need: the gsMap global size stays ni*nj
    ! so it matches the atm grid, which the coupler requires when the atm and lnd
    ! grids are the same (see seq_domain_mct.F90).
    num_cells_owned = num_cells_global / npes
    if (iam < mod(num_cells_global, npes)) num_cells_owned = num_cells_owned + 1

    allocate(natural_id_cells_owned(num_cells_owned))
    do i = 1, num_cells_owned
       natural_id_cells_owned(i) = land_ids(iam + 1 + (i-1)*npes)
    end do

    deallocate(land_ids)

    if (masterproc) then
       write(logunit,*) subname,'grid cells = ',num_cells_grid, &
                        ' active land cells = ',num_cells_global
    end if
    write(logunit,*) subname,'rank ',iam,' owns ',num_cells_owned, &
                     ' of ',num_cells_global,' active land cells'
    call shr_sys_flush(logunit)

    nstep = 0

    if (masterproc) then
       write(logunit,*) 'ELMxx model initialization completed'
       call shr_sys_flush(logunit)
    end if

  end subroutine elmxx_init

  !-----------------------------------------------------------------------
  subroutine elmxx_run(logunit, coupling_dt_in_sec)
    !
    ! !DESCRIPTION:
    ! Advance ELMxx one coupling interval.
    !
    ! This is intentionally a no-op: the Kokkos port is not wired in yet, so
    ! ELMxx consumes nothing from the coupler and returns zeroed fields.
    !
    implicit none
    !
    integer, intent(in) :: logunit
    integer, intent(in) :: coupling_dt_in_sec

    nstep = nstep + 1

    if (masterproc) then
       write(logunit,*) 'ELMxx step ',nstep,' dt = ',coupling_dt_in_sec,' s (no-op)'
       call shr_sys_flush(logunit)
    end if

  end subroutine elmxx_run

  !-----------------------------------------------------------------------
  subroutine elmxx_final()
    !
    ! !DESCRIPTION:
    ! Finalize ELMxx.
    !
    implicit none

    if (associated(natural_id_cells_owned)) deallocate(natural_id_cells_owned)
    if (associated(lonc_g))  deallocate(lonc_g)
    if (associated(latc_g))  deallocate(latc_g)
    if (associated(areac_g)) deallocate(areac_g)
    if (associated(maskc_g)) deallocate(maskc_g)
    if (associated(fracc_g)) deallocate(fracc_g)

  end subroutine elmxx_final

end module elmxxMod
