program surface_albedo_replay

  !-----------------------------------------------------------------------
  ! Offline replay of ELMxx's surface albedo kernel.
  !
  ! Reads a deck of ELM's own SurfaceAlbedo inputs, runs
  ! elmxx_surface_albedo_kernel -- the same source the coupled run calls, not
  ! a copy of it -- and writes every output as text for comparison against
  ! ELM's own outputs.
  !
  !   ./surface_albedo_replay <deck> <results>
  !
  ! The deck is written by tools/validate_surface_albedo.py, which extracts
  ! it from an ELM restart. Deck format, list-directed, in this exact order:
  !
  !   nc np npft
  !   patch_col(np)          1-based column index per patch, <= 0 to skip
  !   coszen(nc)
  !   soil_color(nc)
  !   h2osoi_vol_top(nc)
  !   frac_sno(nc)
  !   patch_ivt(np)          0-based PFT index
  !   elai(np) esai(np) t_veg(np) fwet(np)
  !   rhol(0:npft-1,1) rhol(:,2) rhos(:,1) rhos(:,2)
  !   taul(:,1) taul(:,2) taus(:,1) taus(:,2) xl(:)
  !-----------------------------------------------------------------------

  use shr_kind_mod, only : r8 => shr_kind_r8
  use elmxxSurfaceAlbedoKernelMod, only : numrad, elmxx_surface_albedo_kernel

  implicit none

  integer :: nc, np, npft, u, ib, ios
  character(len=512) :: deck, results

  integer , allocatable :: patch_col(:), soil_color(:), patch_ivt(:), nrad(:)
  real(r8), allocatable :: coszen(:), h2osoi(:), fracsno(:)
  real(r8), allocatable :: elai(:), esai(:), tveg(:), fwet(:)
  real(r8), allocatable :: rhol(:,:), rhos(:,:), taul(:,:), taus(:,:), xl(:)
  real(r8), allocatable :: albsod(:,:), albsoi(:,:), albgrd(:,:), albgri(:,:)
  real(r8), allocatable :: albd(:,:), albi(:,:), fabd(:,:), fabi(:,:)
  real(r8), allocatable :: ftdd(:,:), ftid(:,:), ftii(:,:)
  real(r8), allocatable :: tlai_z(:), fsun_z(:)
  real(r8), allocatable :: fabd_sun_z(:), fabi_sun_z(:), fabd_sha_z(:), fabi_sha_z(:)
  real(r8), allocatable :: vcmaxcintsun(:), vcmaxcintsha(:)

  if (command_argument_count() /= 2) then
     write(*,*) 'usage: surface_albedo_replay <deck> <results>'
     stop 2
  end if
  call get_command_argument(1, deck)
  call get_command_argument(2, results)

  open(newunit=u, file=trim(deck), status='old', action='read', iostat=ios)
  if (ios /= 0) then
     write(*,*) 'cannot open deck ', trim(deck)
     stop 2
  end if

  read(u,*) nc, np, npft

  allocate(patch_col(np), soil_color(nc), patch_ivt(np), nrad(np))
  allocate(coszen(nc), h2osoi(nc), fracsno(nc))
  allocate(elai(np), esai(np), tveg(np), fwet(np))
  allocate(rhol(0:npft-1,numrad), rhos(0:npft-1,numrad), &
           taul(0:npft-1,numrad), taus(0:npft-1,numrad), xl(0:npft-1))
  allocate(albsod(nc,numrad), albsoi(nc,numrad), albgrd(nc,numrad), albgri(nc,numrad))
  allocate(albd(np,numrad), albi(np,numrad), fabd(np,numrad), fabi(np,numrad))
  allocate(ftdd(np,numrad), ftid(np,numrad), ftii(np,numrad))
  allocate(tlai_z(np), fsun_z(np))
  allocate(fabd_sun_z(np), fabi_sun_z(np), fabd_sha_z(np), fabi_sha_z(np))
  allocate(vcmaxcintsun(np), vcmaxcintsha(np))

  read(u,*) patch_col
  read(u,*) coszen
  read(u,*) soil_color
  read(u,*) h2osoi
  read(u,*) fracsno
  read(u,*) patch_ivt
  read(u,*) elai
  read(u,*) esai
  read(u,*) tveg
  read(u,*) fwet
  do ib = 1, numrad
     read(u,*) rhol(:,ib)
  end do
  do ib = 1, numrad
     read(u,*) rhos(:,ib)
  end do
  do ib = 1, numrad
     read(u,*) taul(:,ib)
  end do
  do ib = 1, numrad
     read(u,*) taus(:,ib)
  end do
  read(u,*) xl
  close(u)

  call elmxx_surface_albedo_kernel(nc, np, npft, patch_col, coszen,   &
       soil_color, h2osoi, fracsno, patch_ivt, elai, esai, tveg, fwet, &
       rhol, rhos, taul, taus, xl,                                     &
       albsod, albsoi, albgrd, albgri,                                 &
       albd, albi, fabd, fabi, ftdd, ftid, ftii,                       &
       nrad, tlai_z, fsun_z, fabd_sun_z, fabi_sun_z, fabd_sha_z, fabi_sha_z,  &
       vcmaxcintsun, vcmaxcintsha)

  open(newunit=u, file=trim(results), status='replace', action='write')
  call put_col(u, 'albsod', albsod)
  call put_col(u, 'albsoi', albsoi)
  call put_col(u, 'albgrd', albgrd)
  call put_col(u, 'albgri', albgri)
  call put_pft(u, 'albd', albd)
  call put_pft(u, 'albi', albi)
  call put_pft(u, 'fabd', fabd)
  call put_pft(u, 'fabi', fabi)
  call put_pft(u, 'ftdd', ftdd)
  call put_pft(u, 'ftid', ftid)
  call put_pft(u, 'ftii', ftii)
  call put_1d(u, 'tlai_z', tlai_z)
  call put_1d(u, 'fsun_z', fsun_z)
  call put_1d(u, 'fabd_sun_z', fabd_sun_z)
  call put_1d(u, 'fabi_sun_z', fabi_sun_z)
  call put_1d(u, 'fabd_sha_z', fabd_sha_z)
  call put_1d(u, 'fabi_sha_z', fabi_sha_z)
  call put_1d(u, 'vcmaxcintsun', vcmaxcintsun)
  call put_1d(u, 'vcmaxcintsha', vcmaxcintsha)
  call put_1d(u, 'nrad', real(nrad, r8))
  close(u)

contains

  subroutine put_col(u, name, a)
    integer, intent(in) :: u
    character(len=*), intent(in) :: name
    real(r8), intent(in) :: a(:,:)
    integer :: ib
    character(len=3), parameter :: band(2) = (/ 'vis', 'nir' /)
    do ib = 1, numrad
       call put_1d(u, name//'_'//band(ib), a(:,ib))
    end do
  end subroutine put_col

  subroutine put_pft(u, name, a)
    integer, intent(in) :: u
    character(len=*), intent(in) :: name
    real(r8), intent(in) :: a(:,:)
    call put_col(u, name, a)
  end subroutine put_pft

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

end program surface_albedo_replay
