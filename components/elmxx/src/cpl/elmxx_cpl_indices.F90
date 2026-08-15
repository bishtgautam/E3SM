module elmxx_cpl_indices

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Indices into the coupler's x2l and l2x attribute vectors.
  !
  ! Scope is the SP-mode subset: the atmospheric forcing the land physics
  ! actually consumes, and the state and fluxes it will send back. ELM's
  ! elm_cpl_indices.F90 carries 54 x2l fields and as many l2x, including
  ! aerosol deposition, dust, CO2 diagnostics, river feedbacks and glacier
  ! coupling. Those are for capability ELMxx does not have; carrying indices
  ! for fields nothing reads would be dead weight that looks like coverage.
  ! ELM's list is the reference for the NAMES -- they must match the coupler's
  ! seq_flds strings exactly -- not for the extent.
  !
  ! mct_avect_indexra returns 0 for a field the coupler does not carry, so
  ! every lookup is checked. A silent 0 would index element 0 of the attribute
  ! vector at run time.
  !-----------------------------------------------------------------------

  use shr_sys_mod, only : shr_sys_abort
  use mct_mod    , only : mct_aVect, mct_avect_indexra

  implicit none
  save
  private

  !--------------------------------------------------------------------------
  ! atm -> lnd
  !--------------------------------------------------------------------------
  integer, public :: index_x2l_Sa_z        = 0  ! bottom atm level height (m)
  integer, public :: index_x2l_Sa_topo     = 0  ! atm surface height (m)
  integer, public :: index_x2l_Sa_u        = 0  ! bottom atm level zonal wind (m/s)
  integer, public :: index_x2l_Sa_v        = 0  ! bottom atm level meridional wind (m/s)
  integer, public :: index_x2l_Sa_ptem     = 0  ! bottom atm level potential temp (K)
  integer, public :: index_x2l_Sa_shum     = 0  ! bottom atm level specific humidity (kg/kg)
  integer, public :: index_x2l_Sa_pbot     = 0  ! bottom atm level pressure (Pa)
  integer, public :: index_x2l_Sa_tbot     = 0  ! bottom atm level temperature (K)
  integer, public :: index_x2l_Faxa_lwdn   = 0  ! downward longwave (W/m2)
  integer, public :: index_x2l_Faxa_rainc  = 0  ! convective rain (kg/m2/s)
  integer, public :: index_x2l_Faxa_rainl  = 0  ! large-scale rain (kg/m2/s)
  integer, public :: index_x2l_Faxa_snowc  = 0  ! convective snow (kg/m2/s)
  integer, public :: index_x2l_Faxa_snowl  = 0  ! large-scale snow (kg/m2/s)
  integer, public :: index_x2l_Faxa_swndr  = 0  ! direct near-infrared shortwave (W/m2)
  integer, public :: index_x2l_Faxa_swvdr  = 0  ! direct visible shortwave (W/m2)
  integer, public :: index_x2l_Faxa_swndf  = 0  ! diffuse near-infrared shortwave (W/m2)
  integer, public :: index_x2l_Faxa_swvdf  = 0  ! diffuse visible shortwave (W/m2)

  !--------------------------------------------------------------------------
  ! lnd -> atm
  !--------------------------------------------------------------------------
  integer, public :: index_l2x_Sl_t        = 0  ! surface temperature (K)
  integer, public :: index_l2x_Sl_tref     = 0  ! 2m reference temperature (K)
  integer, public :: index_l2x_Sl_qref     = 0  ! 2m reference specific humidity (kg/kg)
  integer, public :: index_l2x_Sl_avsdr    = 0  ! albedo, direct visible
  integer, public :: index_l2x_Sl_anidr    = 0  ! albedo, direct near-infrared
  integer, public :: index_l2x_Sl_avsdf    = 0  ! albedo, diffuse visible
  integer, public :: index_l2x_Sl_anidf    = 0  ! albedo, diffuse near-infrared
  integer, public :: index_l2x_Sl_snowh    = 0  ! snow height (m)
  integer, public :: index_l2x_Sl_u10      = 0  ! 10m wind (m/s)
  integer, public :: index_l2x_Sl_fv       = 0  ! friction velocity (m/s)
  integer, public :: index_l2x_Sl_ram1     = 0  ! aerodynamical resistance (s/m)
  integer, public :: index_l2x_Fall_taux   = 0  ! zonal surface stress (N/m2)
  integer, public :: index_l2x_Fall_tauy   = 0  ! meridional surface stress (N/m2)
  integer, public :: index_l2x_Fall_lat    = 0  ! latent heat flux (W/m2)
  integer, public :: index_l2x_Fall_sen    = 0  ! sensible heat flux (W/m2)
  integer, public :: index_l2x_Fall_lwup   = 0  ! upward longwave (W/m2)
  integer, public :: index_l2x_Fall_evap   = 0  ! evaporation (kg/m2/s)
  integer, public :: index_l2x_Fall_swnet  = 0  ! net shortwave absorbed (W/m2)

  logical, public :: cpl_indices_set = .false.

  public :: elmxx_cpl_indices_set

contains

  !-----------------------------------------------------------------------
  subroutine elmxx_cpl_indices_set(x2l, l2x)
    !
    ! !DESCRIPTION:
    ! Look up every index once, after the attribute vectors are initialized.
    !
    implicit none
    !
    type(mct_aVect), intent(in) :: x2l, l2x

    ! ---- atm -> lnd ----
    call setx(x2l, 'Sa_z'       , index_x2l_Sa_z)
    call setx(x2l, 'Sa_topo'    , index_x2l_Sa_topo)
    call setx(x2l, 'Sa_u'       , index_x2l_Sa_u)
    call setx(x2l, 'Sa_v'       , index_x2l_Sa_v)
    call setx(x2l, 'Sa_ptem'    , index_x2l_Sa_ptem)
    call setx(x2l, 'Sa_shum'    , index_x2l_Sa_shum)
    call setx(x2l, 'Sa_pbot'    , index_x2l_Sa_pbot)
    call setx(x2l, 'Sa_tbot'    , index_x2l_Sa_tbot)
    call setx(x2l, 'Faxa_lwdn'  , index_x2l_Faxa_lwdn)
    call setx(x2l, 'Faxa_rainc' , index_x2l_Faxa_rainc)
    call setx(x2l, 'Faxa_rainl' , index_x2l_Faxa_rainl)
    call setx(x2l, 'Faxa_snowc' , index_x2l_Faxa_snowc)
    call setx(x2l, 'Faxa_snowl' , index_x2l_Faxa_snowl)
    call setx(x2l, 'Faxa_swndr' , index_x2l_Faxa_swndr)
    call setx(x2l, 'Faxa_swvdr' , index_x2l_Faxa_swvdr)
    call setx(x2l, 'Faxa_swndf' , index_x2l_Faxa_swndf)
    call setx(x2l, 'Faxa_swvdf' , index_x2l_Faxa_swvdf)

    ! ---- lnd -> atm ----
    call setx(l2x, 'Sl_t'       , index_l2x_Sl_t)
    call setx(l2x, 'Sl_tref'    , index_l2x_Sl_tref)
    call setx(l2x, 'Sl_qref'    , index_l2x_Sl_qref)
    call setx(l2x, 'Sl_avsdr'   , index_l2x_Sl_avsdr)
    call setx(l2x, 'Sl_anidr'   , index_l2x_Sl_anidr)
    call setx(l2x, 'Sl_avsdf'   , index_l2x_Sl_avsdf)
    call setx(l2x, 'Sl_anidf'   , index_l2x_Sl_anidf)
    call setx(l2x, 'Sl_snowh'   , index_l2x_Sl_snowh)
    call setx(l2x, 'Sl_u10'     , index_l2x_Sl_u10)
    call setx(l2x, 'Sl_fv'      , index_l2x_Sl_fv)
    call setx(l2x, 'Sl_ram1'    , index_l2x_Sl_ram1)
    call setx(l2x, 'Fall_taux'  , index_l2x_Fall_taux)
    call setx(l2x, 'Fall_tauy'  , index_l2x_Fall_tauy)
    call setx(l2x, 'Fall_lat'   , index_l2x_Fall_lat)
    call setx(l2x, 'Fall_sen'   , index_l2x_Fall_sen)
    call setx(l2x, 'Fall_lwup'  , index_l2x_Fall_lwup)
    call setx(l2x, 'Fall_evap'  , index_l2x_Fall_evap)
    call setx(l2x, 'Fall_swnet' , index_l2x_Fall_swnet)

    cpl_indices_set = .true.

  end subroutine elmxx_cpl_indices_set

  !-----------------------------------------------------------------------
  subroutine setx(av, fldname, idx)
    !
    ! !DESCRIPTION:
    ! Look up one field and refuse a missing one.
    !
    ! mct_avect_indexra returns 0 when the field is absent, which would index
    ! element 0 later. Fail here, naming the field, rather than at the first
    ! read.
    !
    implicit none
    type(mct_aVect) , intent(in)  :: av
    character(len=*), intent(in)  :: fldname
    integer         , intent(out) :: idx

    idx = mct_avect_indexra(av, trim(fldname), perrWith='quiet')
    if (idx <= 0) then
       call shr_sys_abort('(elmxx_cpl_indices_set) ERROR: coupler field not found: '// &
                          trim(fldname))
    end if

  end subroutine setx

end module elmxx_cpl_indices
