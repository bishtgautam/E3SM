module elmxxTimerMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Timers for ELMxx, named after ELM's where the work matches, so the two
  ! models' timing files can be read side by side.
  !
  ! Kernels launch asynchronously on a GPU, so a timer stopped right after a
  ! launch would measure the launch, not the kernel. Each stop fences first.
  !
  ! Timers stay off until lnd_init_mct turns them on, so elmxx_perf, which has
  ! its own timing and no GPTL setup, runs with no extra fences.
  !-----------------------------------------------------------------------

  use perf_mod              , only : t_startf, t_stopf
  use elmxx_kokkos_interface, only : ELMxxKokkosFence

  implicit none
  private

  logical, public :: elmxx_timers_on = .false.

  public :: elmxx_t_start
  public :: elmxx_t_stop

contains

  !-----------------------------------------------------------------------
  subroutine elmxx_t_start(name)
    character(len=*), intent(in) :: name

    if (elmxx_timers_on) call t_startf(name)

  end subroutine elmxx_t_start

  !-----------------------------------------------------------------------
  subroutine elmxx_t_stop(name)
    character(len=*), intent(in) :: name

    if (elmxx_timers_on) then
       call ELMxxKokkosFence()
       call t_stopf(name)
    end if

  end subroutine elmxx_t_stop

end module elmxxTimerMod
