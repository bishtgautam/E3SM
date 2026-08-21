program root_stress_replay

  !-----------------------------------------------------------------------
  ! Offline replay of ELMxx's root water stress calculation.
  !
  ! Reads a deck of ELM's own calc_root_moist_stress inputs, runs
  ! elmxx_root_stress_kernel -- the same source the coupled run calls, not a
  ! copy of it -- and writes rootr, btran and rresis for comparison against
  ! ELM's own outputs.
  !
  !   ./root_stress_replay <deck> <results>
  !
  ! The deck is written by tools/validate_root_stress.py from ELM's diagnostic
  ! binary. Deck format, list-directed, in this order:
  !
  !   np nc nlevbed nlevgrnd npft
  !   tc_stress btran0 denice denh2o
  !   patch_itype(np)   patch_col(np)          one line each
  !   smpsc(0:npft-1)   smpso(0:npft-1)        one line each
  !   rootfr(np,nlevgrnd)                      one line per patch
  !   h2osoi_liq(nc,nlevgrnd)                  one line per column
  !   h2osoi_ice(nc,nlevgrnd)
  !   dz(nc,nlevgrnd)
  !   t_soisno(nc,nlevgrnd)
  !   watsat(nc,nlevgrnd)
  !   bsw(nc,nlevgrnd)
  !   sucsat(nc,nlevgrnd)
  !-----------------------------------------------------------------------

  use shr_kind_mod      , only : r8 => shr_kind_r8
  use elmxxRootKernelMod, only : elmxx_root_stress_kernel

  implicit none

  integer :: np, nc, nlevbed, nlevgrnd, npft, u, p, c, ios
  real(r8) :: tc_stress, btran0, denice, denh2o
  character(len=512) :: deck, results

  integer , allocatable :: patch_itype(:), patch_col(:)
  real(r8), allocatable :: smpsc(:), smpso(:), rootfr(:,:)
  real(r8), allocatable :: h2osoi_liq(:,:), h2osoi_ice(:,:), dz(:,:)
  real(r8), allocatable :: t_soisno(:,:), watsat(:,:), bsw(:,:), sucsat(:,:)
  real(r8), allocatable :: rootr(:,:), btran(:), rresis(:,:)

  if (command_argument_count() /= 2) then
     write(*,*) 'usage: root_stress_replay <deck> <results>'
     stop 1
  end if
  call get_command_argument(1, deck)
  call get_command_argument(2, results)

  open(newunit=u, file=trim(deck), status='old', action='read', iostat=ios)
  if (ios /= 0) then
     write(*,*) 'cannot open deck ', trim(deck); stop 1
  end if

  read(u,*) np, nc, nlevbed, nlevgrnd, npft
  read(u,*) tc_stress, btran0, denice, denh2o

  allocate(patch_itype(np), patch_col(np))
  allocate(smpsc(0:npft-1), smpso(0:npft-1))
  allocate(rootfr(np,nlevgrnd))
  allocate(h2osoi_liq(nc,nlevgrnd), h2osoi_ice(nc,nlevgrnd), dz(nc,nlevgrnd))
  allocate(t_soisno(nc,nlevgrnd), watsat(nc,nlevgrnd), bsw(nc,nlevgrnd))
  allocate(sucsat(nc,nlevgrnd))
  allocate(rootr(np,nlevgrnd), btran(np), rresis(np,nlevgrnd))

  read(u,*) patch_itype
  read(u,*) patch_col
  read(u,*) smpsc
  read(u,*) smpso
  do p = 1, np
     read(u,*) rootfr(p,:)
  end do
  do c = 1, nc
     read(u,*) h2osoi_liq(c,:)
  end do
  do c = 1, nc
     read(u,*) h2osoi_ice(c,:)
  end do
  do c = 1, nc
     read(u,*) dz(c,:)
  end do
  do c = 1, nc
     read(u,*) t_soisno(c,:)
  end do
  do c = 1, nc
     read(u,*) watsat(c,:)
  end do
  do c = 1, nc
     read(u,*) bsw(c,:)
  end do
  do c = 1, nc
     read(u,*) sucsat(c,:)
  end do
  close(u)

  call elmxx_root_stress_kernel(np, nc, nlevbed, nlevgrnd,        &
       patch_itype, patch_col, rootfr,                            &
       h2osoi_liq, h2osoi_ice, dz, t_soisno,                      &
       watsat, bsw, sucsat, smpsc, smpso,                         &
       tc_stress, btran0, denice, denh2o,                         &
       rootr, btran, rresis)

  open(newunit=u, file=trim(results), status='replace', action='write')
  write(u,*) np, nlevgrnd
  write(u,*) btran
  do p = 1, np
     write(u,*) rootr(p,:)
  end do
  do p = 1, np
     write(u,*) rresis(p,:)
  end do
  close(u)

end program root_stress_replay
