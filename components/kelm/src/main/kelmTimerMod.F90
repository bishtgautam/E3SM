module kelmTimerMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Timers for kELM, named after ELM's where the work matches, so the two
  ! models' timing files can be read side by side.
  !
  ! Kernels launch asynchronously on a GPU, so a timer stopped right after a
  ! launch would measure the launch, not the kernel. Each stop fences first.
  !
  ! Timers stay off until lnd_init_mct turns them on, so kelm_perf, which has
  ! its own timing and no GPTL setup, runs with no extra fences.
  !-----------------------------------------------------------------------

  use perf_mod              , only : t_startf, t_stopf
  use kelm_kokkos_interface, only : kELMKokkosFence

  implicit none
  private

  logical, public :: kelm_timers_on = .false.

  public :: kelm_t_start
  public :: kelm_t_stop

contains

  !-----------------------------------------------------------------------
  subroutine kelm_t_start(name)
    character(len=*), intent(in) :: name

    if (kelm_timers_on) call t_startf(name)

  end subroutine kelm_t_start

  !-----------------------------------------------------------------------
  subroutine kelm_t_stop(name)
    character(len=*), intent(in) :: name

    if (kelm_timers_on) then
       call kELMKokkosFence()
       call t_stopf(name)
    end if

  end subroutine kelm_t_stop

end module kelmTimerMod
