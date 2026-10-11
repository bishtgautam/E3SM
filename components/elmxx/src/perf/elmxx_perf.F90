program elmxx_perf

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Standalone perf driver for ELMxx (plan Stage 6.8). ELMxx-only.
  !
  ! It runs the coupled model's own code path -- elmxx_init, elmxx_init_albedo,
  ! elmxx_run, elmxx_final, everything lnd_comp_mct calls -- without the
  ! coupler, DATM or ESMF. The case's lnd_in (read from the working
  ! directory) picks the source cell(s), the input files and the kernels.
  ! With --ncells N the decomposition deals N copies of the source cells
  ! round-robin over the MPI ranks (elmxxMod, perf mode), and every copy
  ! gets the same synthetic forcing (elmxxPerfForcingMod), so every copy must
  ! end with bit-identical state -- checked at the end (elmxxPerfCheckMod).
  !
  ! The step sequence mirrors lnd_run_mct on a noleap calendar: an nstep-0
  ! pass without albedo, then nsteps steps. Pass n runs with the clock at
  ! start + n*dt, and nextsw_cday is the calendar day one step later.
  !
  ! PIO is set up as the coupler does it (shr_pio_init1/2, one component,
  ! LND), from the case's drv_in and lnd_modelio.nml, so the I/O layout is
  ! the coupled run's.
  !
  ! Usage (from a directory holding lnd_in, drv_in and lnd_modelio.nml):
  !   mpirun -np P elmxx_perf --ncells N [--nsteps 48] [--warmup 2] [--dt 1800]
  !          [--start-doy 1] [--orb-year 1850] [--no-check]
  !          [--write-reference FILE] [--reference FILE] [--tolerance TOL]
  !          [--write-restart] [--perturb-cell G]
  !
  ! --perturb-cell G is the negative control: virtual cell G's air
  ! temperature gets +1e-9 K every step, so the checks must FAIL.
  !-----------------------------------------------------------------------

  use shr_kind_mod        , only : r8 => shr_kind_r8
  use shr_sys_mod         , only : shr_sys_abort, shr_sys_flush
  use shr_orb_mod         , only : shr_orb_params, shr_orb_decl
  use shr_pio_mod         , only : shr_pio_init1, shr_pio_init2, shr_pio_finalize
  use elmxxSpmdMod        , only : elmxxSpmdInit, masterproc, iam, npes, mpicom_lnd
  use elmxxMod            , only : elmxx_read_namelist, elmxx_init, elmxx_init_albedo, &
                                   elmxx_run, elmxx_final, elmxx_restart_write, &
                                   elmxx_state, elmxx_caseid, inst_name, inst_suffix, &
                                   inst_index, num_cells_owned, num_cells_global, &
                                   cell_lat, cell_lon, elmxx_perf_ncells, &
                                   elmxx_perf_nsrc, elmxx_perf_src_owned
  use elmxx_kokkos_interface, only : ELMxxKokkosFence
  use elmxxPerfForcingMod , only : elmxx_perf_forcing_set
  use elmxxPerfCheckMod   , only : elmxx_perf_check_replicas, &
                                   elmxx_perf_write_reference, &
                                   elmxx_perf_compare_reference

  implicit none
#include <mpif.h>

  integer, parameter :: logunit = 6
  integer, parameter :: days_per_month(12) = [31,28,31,30,31,30,31,31,30,31,30,31]

  ! ---- options ----
  integer  :: ncells = 0, nsteps = 48, nwarmup = 2, dt = 1800
  integer  :: start_doy = 1, start_year = 1, orb_year = 1850
  integer  :: perturb_cell = 0
  logical  :: do_check = .true., do_restart = .false.
  character(len=512) :: ref_out = ' ', ref_in = ' '
  real(r8) :: tol = 0._r8

  ! ---- orbit ----
  real(r8) :: eccen, obliq, mvelp, obliqr, lambm0, mvelpp, eccf

  ! ---- run ----
  integer  :: ierr, n, ntimed, nbad, nbad_total, comm
  integer  :: yr, mon, day, tod
  real(r8) :: calday, nextsw_cday, declin, declinp1
  real(r8) :: t0, t1, t2, t_init, t_force, t_run, t_total
  real(r8) :: tmax_run, tmin_run, tsum_run, tmax_force, tmax_init, tmax_total

  call mpi_init(ierr)
  call elmxxSpmdInit(MPI_COMM_WORLD)
  call parse_args()

  inst_name   = 'LND'
  inst_suffix = ''
  inst_index  = 1
  elmxx_caseid = 'elmxx_perf'

  ! PIO, as the coupler sets it up; elmxx_init then takes LND's iosystem
  ! through elmxx_pio_init like the coupled model does.
  comm = MPI_COMM_WORLD
  call shr_pio_init1(1, 'drv_in', comm)
  call shr_pio_init2([1], [inst_name], [.true.], [mpicom_lnd], [iam])

  call shr_orb_params(orb_year, eccen, obliq, mvelp, obliqr, lambm0, mvelpp, masterproc)

  call elmxx_read_namelist(logunit)
  elmxx_perf_ncells = ncells

  ! ---- initialization, timed as one region ----
  call mpi_barrier(mpicom_lnd, ierr)
  t0 = mpi_wtime()
  call clock_at(0, yr, mon, day, tod, calday)
  call elmxx_init(logunit, yr, mon, day)
  call next_radiation_day(dt, nextsw_cday)
  call shr_orb_decl(nextsw_cday, eccen, mvelpp, lambm0, obliqr, declinp1, eccf)
  call elmxx_init_albedo(logunit, nextsw_cday, declinp1)
  call ELMxxKokkosFence()
  t_init = mpi_wtime() - t0

  if (masterproc) then
     write(logunit,'(a,i0,a,i0,a,i0,a)') ' elmxx_perf: ', num_cells_global, ' cells (', &
          elmxx_perf_nsrc, ' sources) on ', npes, ' ranks'
     write(logunit,'(a,i0,a,i0,a,i0,a)') ' elmxx_perf: ', nsteps, ' steps of ', dt, ' s; first ', &
          nwarmup, ' passes (nstep 0 included) untimed'
     call shr_sys_flush(logunit)
  end if

  ! ---- time loop: nstep 0, then 1..nsteps ----
  t_force = 0._r8
  t_run   = 0._r8
  ntimed  = 0
  call mpi_barrier(mpicom_lnd, ierr)
  t_total = mpi_wtime()
  do n = 0, nsteps
     call next_radiation_day((n+1)*dt, nextsw_cday)
     call shr_orb_decl(nextsw_cday, eccen, mvelpp, lambm0, obliqr, declinp1, eccf)
     call clock_at(n*dt, yr, mon, day, tod, calday)
     call shr_orb_decl(calday, eccen, mvelpp, lambm0, obliqr, declin, eccf)

     t0 = mpi_wtime()
     call elmxx_perf_forcing_set(num_cells_owned, cell_lat, cell_lon, calday, tod, declin)
     if (perturb_cell > 0) call perturb_forcing(perturb_cell)
     call ELMxxKokkosFence()
     t1 = mpi_wtime()
     call elmxx_run(logunit, dt, mon, day, nextsw_cday, declinp1, doalb=(n > 0), &
                    hist_year=yr, hist_month=mon, hist_day=day, hist_tod=tod)
     call ELMxxKokkosFence()
     t2 = mpi_wtime()

     if (n >= nwarmup) then
        t_force = t_force + (t1 - t0)
        t_run   = t_run   + (t2 - t1)
        ntimed  = ntimed + 1
     end if
  end do
  t_total = mpi_wtime() - t_total

  ! ---- timing summary: max over ranks is the wall time ----
  call mpi_reduce(t_run  , tmax_run  , 1, MPI_REAL8, MPI_MAX, 0, mpicom_lnd, ierr)
  call mpi_reduce(t_run  , tmin_run  , 1, MPI_REAL8, MPI_MIN, 0, mpicom_lnd, ierr)
  call mpi_reduce(t_run  , tsum_run  , 1, MPI_REAL8, MPI_SUM, 0, mpicom_lnd, ierr)
  call mpi_reduce(t_force, tmax_force, 1, MPI_REAL8, MPI_MAX, 0, mpicom_lnd, ierr)
  call mpi_reduce(t_init , tmax_init , 1, MPI_REAL8, MPI_MAX, 0, mpicom_lnd, ierr)
  call mpi_reduce(t_total, tmax_total, 1, MPI_REAL8, MPI_MAX, 0, mpicom_lnd, ierr)
  if (masterproc) then
     write(logunit,'(a)') ' ---------------------------------------------------------------'
     write(logunit,'(a,f12.4,a)') ' init            ', tmax_init, ' s (max over ranks)'
     write(logunit,'(a,f12.4,a)') ' all passes      ', tmax_total, ' s (max over ranks)'
     if (ntimed > 0) then
        write(logunit,'(a,es12.4,a,es12.4,a,es12.4,a)') ' elmxx_run/step  ', tmax_run/ntimed, &
             ' s max, ', tsum_run/npes/ntimed, ' s mean, ', tmin_run/ntimed, ' s min'
        write(logunit,'(a,es12.4,a)') ' forcing/step    ', tmax_force/ntimed, ' s max (driver, not model)'
        write(logunit,'(a,es12.4)')   ' cell-steps/s    ', &
             real(num_cells_global, r8) * ntimed / tmax_run
        ! One line for scripts.
        write(logunit,'(a,i0,a,i0,a,i0,a,i0,a,es12.5,a,es12.5,a,es12.5)') ' PERF ncells=', &
             num_cells_global, ' nsrc=', elmxx_perf_nsrc, ' npes=', npes, ' ntimed=', ntimed, &
             ' t_step_max=', tmax_run/ntimed, ' t_step_mean=', tsum_run/npes/ntimed, &
             ' t_init=', tmax_init
     end if
     write(logunit,'(a)') ' ---------------------------------------------------------------'
     call shr_sys_flush(logunit)
  end if

  ! ---- answer checks ----
  nbad_total = 0
  if (do_check) then
     call elmxx_perf_check_replicas(elmxx_state, num_cells_owned, elmxx_perf_src_owned, &
                                    elmxx_perf_nsrc, logunit, nbad)
     nbad_total = nbad_total + nbad
  end if
  if (len_trim(ref_out) > 0) then
     call elmxx_perf_write_reference(elmxx_state, num_cells_owned, elmxx_perf_src_owned, &
                                     elmxx_perf_nsrc, trim(ref_out), logunit)
  end if
  if (len_trim(ref_in) > 0) then
     call elmxx_perf_compare_reference(elmxx_state, num_cells_owned, elmxx_perf_src_owned, &
                                       elmxx_perf_nsrc, trim(ref_in), tol, logunit, nbad)
     nbad_total = nbad_total + nbad
  end if
  if (do_restart) then
     ! Layout-independent (elmxxRestMod), so restarts from different rank
     ! counts compare directly with cprnc.
     call clock_at(nsteps*dt, yr, mon, day, tod, calday)
     call elmxx_restart_write(inst_suffix, yr, mon, day, tod, logunit)
  end if
  if (masterproc .and. (do_check .or. len_trim(ref_in) > 0)) then
     if (nbad_total == 0) then
        write(logunit,'(a)') ' PERF CHECK PASS'
     else
        write(logunit,'(a,i0,a)') ' PERF CHECK FAIL (', nbad_total, ' fields)'
     end if
     call shr_sys_flush(logunit)
  end if

  call clock_at(nsteps*dt, yr, mon, day, tod, calday)
  call elmxx_final(yr, mon, day, tod)
  call shr_pio_finalize()
  call mpi_finalize(ierr)
  if (nbad_total /= 0) error stop 1

contains

  !-----------------------------------------------------------------------
  subroutine clock_at(elapsed, y, m, d, s, cday)
    !
    ! Date at start + elapsed seconds on a noleap calendar. The run starts
    ! at 00:00 on day start_doy of year start_year. cday is the 1-based
    ! calendar day with the fraction of the day, as nextsw_cday is.
    !
    integer , intent(in)  :: elapsed
    integer , intent(out) :: y, m, d, s
    real(r8), intent(out) :: cday
    integer :: days, doy

    days = (start_doy - 1) + elapsed / 86400
    s    = mod(elapsed, 86400)
    y    = start_year + days / 365
    doy  = mod(days, 365) + 1
    cday = real(doy, r8) + real(s, r8) / 86400._r8
    m = 1
    d = doy
    do while (d > days_per_month(m))
       d = d - days_per_month(m)
       m = m + 1
    end do
  end subroutine clock_at

  !-----------------------------------------------------------------------
  subroutine perturb_forcing(g)
    ! Negative control: +1e-9 K on virtual cell g, wherever it lives.
    ! Local cell i is virtual cell iam + 1 + (i-1)*npes.
    use elmxxForcingMod, only : forc_tbot
    integer, intent(in) :: g
    integer :: i
    if (mod(g - 1, npes) /= iam) return
    i = (g - 1) / npes + 1
    if (i <= num_cells_owned) forc_tbot(i) = forc_tbot(i) + 1.e-9_r8
  end subroutine perturb_forcing

  !-----------------------------------------------------------------------
  subroutine next_radiation_day(elapsed, cday)
    ! Calendar day at start + elapsed, the date fields discarded.
    integer , intent(in)  :: elapsed
    real(r8), intent(out) :: cday
    integer :: y, m, d, s
    call clock_at(elapsed, y, m, d, s, cday)
  end subroutine next_radiation_day

  !-----------------------------------------------------------------------
  subroutine parse_args()
    integer :: i, nargs
    character(len=512) :: arg, val

    nargs = command_argument_count()
    i = 1
    do while (i <= nargs)
       call get_command_argument(i, arg)
       select case (trim(arg))
       case ('--no-check')
          do_check = .false.
       case ('--write-restart')
          do_restart = .true.
       case ('--ncells', '--nsteps', '--warmup', '--dt', '--start-doy', '--orb-year', &
             '--write-reference', '--reference', '--tolerance', '--perturb-cell')
          if (i == nargs) call usage('missing value for '//trim(arg))
          call get_command_argument(i+1, val)
          i = i + 1
          select case (trim(arg))
          case ('--ncells');          read(val, *) ncells
          case ('--nsteps');          read(val, *) nsteps
          case ('--warmup');          read(val, *) nwarmup
          case ('--dt');              read(val, *) dt
          case ('--start-doy');       read(val, *) start_doy
          case ('--orb-year');        read(val, *) orb_year
          case ('--write-reference'); ref_out = val
          case ('--reference');       ref_in  = val
          case ('--tolerance');       read(val, *) tol
          case ('--perturb-cell');    read(val, *) perturb_cell
          end select
       case ('--help', '-h')
          call usage(' ')
       case default
          call usage('unknown option '//trim(arg))
       end select
       i = i + 1
    end do
    if (ncells <= 0) call usage('--ncells N is required')
    if (nsteps < 1 .or. dt <= 0 .or. nwarmup < 0) call usage('bad --nsteps/--dt/--warmup')
    if (start_doy < 1 .or. start_doy > 365) call usage('--start-doy must be 1..365')
  end subroutine parse_args

  !-----------------------------------------------------------------------
  subroutine usage(msg)
    character(len=*), intent(in) :: msg
    integer :: ierr2
    if (masterproc) then
       if (len_trim(msg) > 0) write(logunit,'(a,a)') ' elmxx_perf: ', trim(msg)
       write(logunit,'(a)') ' usage: elmxx_perf --ncells N [--nsteps 48] [--warmup 2] [--dt 1800]'
       write(logunit,'(a)') '        [--start-doy 1] [--orb-year 1850] [--no-check]'
       write(logunit,'(a)') '        [--write-reference FILE] [--reference FILE] [--tolerance TOL]'
       write(logunit,'(a)') '        [--write-restart] [--perturb-cell G]'
       write(logunit,'(a)') '        (run where lnd_in, drv_in and lnd_modelio.nml are)'
       call shr_sys_flush(logunit)
    end if
    call mpi_finalize(ierr2)
    if (len_trim(msg) > 0) error stop 2
    stop
  end subroutine usage

end program elmxx_perf
