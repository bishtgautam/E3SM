module elmxxForcingMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Atmospheric forcing imported from the coupler, held per owned gridcell.
  !
  ! This is the per-timestep half of the boundary described in the development
  ! plan (§A.4): forcing in, fluxes out, and nothing else crossing between
  ! Fortran and the device each step. The fields land here first, at gridcell
  ! level, and are pushed down to the subgrid and across the C API by the
  ! physics wiring in a later stage.
  !
  ! No unit conversion happens here. Values are stored exactly as the coupler
  ! delivers them, so that what ELMxx consumes can be compared against what the
  ! coupler sent without a transformation in between. Rain and snow are kept as
  ! their convective and large-scale parts rather than being summed, for the
  ! same reason.
  !-----------------------------------------------------------------------

  use shr_kind_mod, only : r8 => shr_kind_r8
  use shr_sys_mod , only : shr_sys_abort, shr_sys_flush
  use elmxxSpmdMod, only : masterproc, iam
  use mct_mod     , only : mct_aVect
  use elmxx_cpl_indices

  implicit none
  save
  private

  !--------------------------------------------------------------------------
  ! State of the atmosphere at the bottom model level
  !--------------------------------------------------------------------------
  real(r8), public, pointer :: forc_z(:)     => null()  ! reference height (m)
  real(r8), public, pointer :: forc_topo(:)  => null()  ! surface height (m)
  real(r8), public, pointer :: forc_u(:)     => null()  ! zonal wind (m/s)
  real(r8), public, pointer :: forc_v(:)     => null()  ! meridional wind (m/s)
  real(r8), public, pointer :: forc_ptem(:)  => null()  ! potential temperature (K)
  real(r8), public, pointer :: forc_shum(:)  => null()  ! specific humidity (kg/kg)
  real(r8), public, pointer :: forc_pbot(:)  => null()  ! pressure (Pa)
  real(r8), public, pointer :: forc_tbot(:)  => null()  ! temperature (K)

  !--------------------------------------------------------------------------
  ! Fluxes onto the surface
  !--------------------------------------------------------------------------
  real(r8), public, pointer :: forc_lwrad(:) => null()  ! downward longwave (W/m2)
  real(r8), public, pointer :: forc_rainc(:) => null()  ! convective rain (kg/m2/s)
  real(r8), public, pointer :: forc_rainl(:) => null()  ! large-scale rain (kg/m2/s)
  real(r8), public, pointer :: forc_snowc(:) => null()  ! convective snow (kg/m2/s)
  real(r8), public, pointer :: forc_snowl(:) => null()  ! large-scale snow (kg/m2/s)
  real(r8), public, pointer :: forc_swndr(:) => null()  ! direct near-IR shortwave (W/m2)
  real(r8), public, pointer :: forc_swvdr(:) => null()  ! direct visible shortwave (W/m2)
  real(r8), public, pointer :: forc_swndf(:) => null()  ! diffuse near-IR shortwave (W/m2)
  real(r8), public, pointer :: forc_swvdf(:) => null()  ! diffuse visible shortwave (W/m2)

  integer, private :: ncells_f = 0
  integer, private :: nimports = 0

  public :: elmxx_forcing_init
  public :: elmxx_import
  public :: elmxx_forcing_clean

contains

  !-----------------------------------------------------------------------
  subroutine elmxx_forcing_init(ncells)
    !
    implicit none
    integer, intent(in) :: ncells

    ncells_f = ncells

    allocate(forc_z(ncells), forc_topo(ncells), forc_u(ncells), forc_v(ncells), &
             forc_ptem(ncells), forc_shum(ncells), forc_pbot(ncells), &
             forc_tbot(ncells))
    allocate(forc_lwrad(ncells), forc_rainc(ncells), forc_rainl(ncells), &
             forc_snowc(ncells), forc_snowl(ncells), forc_swndr(ncells), &
             forc_swvdr(ncells), forc_swndf(ncells), forc_swvdf(ncells))

    forc_z = 0.0_r8; forc_topo = 0.0_r8; forc_u = 0.0_r8; forc_v = 0.0_r8
    forc_ptem = 0.0_r8; forc_shum = 0.0_r8; forc_pbot = 0.0_r8; forc_tbot = 0.0_r8
    forc_lwrad = 0.0_r8
    forc_rainc = 0.0_r8; forc_rainl = 0.0_r8
    forc_snowc = 0.0_r8; forc_snowl = 0.0_r8
    forc_swndr = 0.0_r8; forc_swvdr = 0.0_r8
    forc_swndf = 0.0_r8; forc_swvdf = 0.0_r8

    nimports = 0

  end subroutine elmxx_forcing_init

  !-----------------------------------------------------------------------
  subroutine elmxx_import(logunit, x2l)
    !
    ! !DESCRIPTION:
    ! Copy this timestep's forcing out of the coupler's attribute vector.
    !
    implicit none
    !
    integer        , intent(in) :: logunit
    type(mct_aVect), intent(in) :: x2l
    !
    integer :: g

    if (.not. cpl_indices_set) then
       call shr_sys_abort('(elmxx_import) ERROR: coupler indices not set')
    end if

    do g = 1, ncells_f
       forc_z(g)     = x2l%rAttr(index_x2l_Sa_z      , g)
       forc_topo(g)  = x2l%rAttr(index_x2l_Sa_topo   , g)
       forc_u(g)     = x2l%rAttr(index_x2l_Sa_u      , g)
       forc_v(g)     = x2l%rAttr(index_x2l_Sa_v      , g)
       forc_ptem(g)  = x2l%rAttr(index_x2l_Sa_ptem   , g)
       forc_shum(g)  = x2l%rAttr(index_x2l_Sa_shum   , g)
       forc_pbot(g)  = x2l%rAttr(index_x2l_Sa_pbot   , g)
       forc_tbot(g)  = x2l%rAttr(index_x2l_Sa_tbot   , g)
       forc_lwrad(g) = x2l%rAttr(index_x2l_Faxa_lwdn , g)
       forc_rainc(g) = x2l%rAttr(index_x2l_Faxa_rainc, g)
       forc_rainl(g) = x2l%rAttr(index_x2l_Faxa_rainl, g)
       forc_snowc(g) = x2l%rAttr(index_x2l_Faxa_snowc, g)
       forc_snowl(g) = x2l%rAttr(index_x2l_Faxa_snowl, g)
       forc_swndr(g) = x2l%rAttr(index_x2l_Faxa_swndr, g)
       forc_swvdr(g) = x2l%rAttr(index_x2l_Faxa_swvdr, g)
       forc_swndf(g) = x2l%rAttr(index_x2l_Faxa_swndf, g)
       forc_swvdf(g) = x2l%rAttr(index_x2l_Faxa_swvdf, g)
    end do

    nimports = nimports + 1

    ! Report the first import only. Its purpose is to show that the coupler is
    ! actually delivering data -- an all-zero first import means the forcing is
    ! not arriving, which otherwise stays invisible until physics produces
    ! nonsense. Repeating it every timestep would bury the signal.
    if (nimports == 1 .and. masterproc) then
       call report_first_import(logunit)
    end if

  end subroutine elmxx_import

  !-----------------------------------------------------------------------
  subroutine report_first_import(logunit)
    !
    implicit none
    integer, intent(in) :: logunit
    character(len=*), parameter :: subname = '(elmxx_import) '

    write(logunit,*) subname,'first forcing import, ',ncells_f,' cells'
    write(logunit,*) subname,'  tbot  range ',minval(forc_tbot) ,maxval(forc_tbot)
    write(logunit,*) subname,'  pbot  range ',minval(forc_pbot) ,maxval(forc_pbot)
    write(logunit,*) subname,'  shum  range ',minval(forc_shum) ,maxval(forc_shum)
    write(logunit,*) subname,'  lwdn  range ',minval(forc_lwrad),maxval(forc_lwrad)
    write(logunit,*) subname,'  swvdr range ',minval(forc_swvdr),maxval(forc_swvdr)
    write(logunit,*) subname,'  rain  range ', &
         minval(forc_rainc+forc_rainl),maxval(forc_rainc+forc_rainl)

    if (all(forc_tbot == 0.0_r8)) then
       write(logunit,*) subname,'WARNING: forc_tbot is all zero; ', &
                        'the coupler does not appear to be delivering forcing'
    end if
    call shr_sys_flush(logunit)

  end subroutine report_first_import

  !-----------------------------------------------------------------------
  subroutine elmxx_forcing_clean()
    !
    implicit none

    if (associated(forc_z))     deallocate(forc_z)
    if (associated(forc_topo))  deallocate(forc_topo)
    if (associated(forc_u))     deallocate(forc_u)
    if (associated(forc_v))     deallocate(forc_v)
    if (associated(forc_ptem))  deallocate(forc_ptem)
    if (associated(forc_shum))  deallocate(forc_shum)
    if (associated(forc_pbot))  deallocate(forc_pbot)
    if (associated(forc_tbot))  deallocate(forc_tbot)
    if (associated(forc_lwrad)) deallocate(forc_lwrad)
    if (associated(forc_rainc)) deallocate(forc_rainc)
    if (associated(forc_rainl)) deallocate(forc_rainl)
    if (associated(forc_snowc)) deallocate(forc_snowc)
    if (associated(forc_snowl)) deallocate(forc_snowl)
    if (associated(forc_swndr)) deallocate(forc_swndr)
    if (associated(forc_swvdr)) deallocate(forc_swvdr)
    if (associated(forc_swndf)) deallocate(forc_swndf)
    if (associated(forc_swvdf)) deallocate(forc_swvdf)

    ncells_f = 0
    nimports = 0

  end subroutine elmxx_forcing_clean

end module elmxxForcingMod
