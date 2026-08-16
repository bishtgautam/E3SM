program ground_heat_flux_replay

  !-----------------------------------------------------------------------
  ! Offline replay of ELMxx's ground surface energy balance.
  !
  ! Reads a deck of ELM's own ComputeGroundHeatFluxAndDeriv inputs, runs
  ! elmxx_ground_heat_flux_kernel -- the same source the coupled run calls,
  ! not a copy of it -- and writes every output as text for comparison against
  ! ELM's own outputs.
  !
  !   ./ground_heat_flux_replay <deck> <results>
  !
  ! The deck is written by tools/validate_ground_heat_flux.py, which extracts
  ! it from ELM's diagnostic binary. Deck format, list-directed, in this order:
  !
  !   nc np nlevsno nlevtot nsnw
  !   sb
  !   patch_col(np)   patch_wt(np)   frac_veg_nosno(np)
  !   emg(nc) htvp(nc) t_grnd(nc) t_h2osfc(nc) snl(nc) forc_lwrad(nc)
  !   t_soisno(nc,nlevtot)          column-major: one line per level
  !   sabg_soil(np) dlrad(np) cgrnd(np)
  !   eflx_sh_snow(np) eflx_sh_soil(np) eflx_sh_h2osfc(np)
  !   qflx_ev_snow(np) qflx_ev_soil(np) qflx_ev_h2osfc(np)
  !   sabg_lyr(np,nsnw)             one line per snow slot
  !-----------------------------------------------------------------------

  use shr_kind_mod, only : r8 => shr_kind_r8
  use elmxxGroundHeatFluxKernelMod, only : elmxx_ground_heat_flux_kernel

  implicit none

  integer :: nc, np, nlevsno, nlevtot, nsnw, u, j, ios, nbad
  real(r8) :: sb
  character(len=512) :: deck, results

  integer , allocatable :: patch_col(:), fvn(:), snl(:)
  real(r8), allocatable :: patch_wt(:), emg(:), htvp(:), t_grnd(:), t_h2osfc(:)
  real(r8), allocatable :: forc_lwrad(:), t_soisno(:,:)
  real(r8), allocatable :: sabg_soil(:), dlrad(:), cgrnd(:), sabg_lyr(:,:)
  real(r8), allocatable :: sh_snow(:), sh_soil(:), sh_h2osfc(:)
  real(r8), allocatable :: ev_snow(:), ev_soil(:), ev_h2osfc(:)
  real(r8), allocatable :: hs_soil(:), hs_top_snow(:), hs_h2osfc(:), dhsdT(:)
  real(r8), allocatable :: sabg_lyr_col(:,:)

  if (command_argument_count() /= 2) then
     write(*,*) 'usage: ground_heat_flux_replay <deck> <results>'
     stop 2
  end if
  call get_command_argument(1, deck)
  call get_command_argument(2, results)

  open(newunit=u, file=trim(deck), status='old', action='read', iostat=ios)
  if (ios /= 0) then
     write(*,*) 'cannot open deck ', trim(deck)
     stop 2
  end if

  read(u,*) nc, np, nlevsno, nlevtot, nsnw
  read(u,*) sb

  allocate(patch_col(np), fvn(np), patch_wt(np), snl(nc))
  allocate(emg(nc), htvp(nc), t_grnd(nc), t_h2osfc(nc), forc_lwrad(nc))
  allocate(t_soisno(nc, nlevtot), sabg_lyr(np, nsnw))
  allocate(sabg_soil(np), dlrad(np), cgrnd(np))
  allocate(sh_snow(np), sh_soil(np), sh_h2osfc(np))
  allocate(ev_snow(np), ev_soil(np), ev_h2osfc(np))
  allocate(hs_soil(nc), hs_top_snow(nc), hs_h2osfc(nc), dhsdT(nc))
  allocate(sabg_lyr_col(nc, nlevtot))

  read(u,*) patch_col
  read(u,*) patch_wt
  read(u,*) fvn
  read(u,*) emg
  read(u,*) htvp
  read(u,*) t_grnd
  read(u,*) t_h2osfc
  read(u,*) snl
  read(u,*) forc_lwrad
  do j = 1, nlevtot
     read(u,*) t_soisno(:,j)
  end do
  read(u,*) sabg_soil
  read(u,*) dlrad
  read(u,*) cgrnd
  read(u,*) sh_snow
  read(u,*) sh_soil
  read(u,*) sh_h2osfc
  read(u,*) ev_snow
  read(u,*) ev_soil
  read(u,*) ev_h2osfc
  do j = 1, nsnw
     read(u,*) sabg_lyr(:,j)
  end do
  close(u)

  call elmxx_ground_heat_flux_kernel(nc, np, nlevsno, nlevtot, nsnw, sb, &
       patch_col, patch_wt, fvn,                                         &
       emg, htvp, t_grnd, t_h2osfc, snl, t_soisno, forc_lwrad,           &
       sabg_soil, sabg_lyr, dlrad, cgrnd, sh_snow, sh_soil, sh_h2osfc,   &
       ev_snow, ev_soil, ev_h2osfc,                                      &
       hs_soil, hs_top_snow, hs_h2osfc, dhsdT, sabg_lyr_col, nbad)

  open(newunit=u, file=trim(results), status='replace', action='write')
  call put_1d(u, 'hs_soil', hs_soil)
  call put_1d(u, 'hs_top_snow', hs_top_snow)
  call put_1d(u, 'hs_h2osfc', hs_h2osfc)
  call put_1d(u, 'dhsdT', dhsdT)
  call put_1d(u, 'n_nonfinite', (/ real(nbad, r8) /))
  close(u)

contains

  subroutine put_1d(u, name, a)
    integer, intent(in) :: u
    character(len=*), intent(in) :: name
    real(r8), intent(in) :: a(:)
    integer :: i
    write(u,'(a)',advance='no') trim(name)
    do i = 1, size(a)
       write(u,'(1x,es24.16)',advance='no') a(i)
    end do
    write(u,'(a)') ''
  end subroutine put_1d

end program ground_heat_flux_replay
