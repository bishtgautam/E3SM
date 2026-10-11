module elmxxPerfForcingMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Built-in synthetic atmospheric forcing for the perf driver (plan Stage
  ! 6.8). ELMxx-only: no ELM counterpart.
  !
  ! It fills the same per-cell arrays elmxx_import fills from the coupler
  ! (elmxxForcingMod), so elmxx_run cannot tell the two apart. Every value is
  ! a function of the model time and the cell's lat/lon only. Replicas of one
  ! source cell share its lat/lon, so they get bit-identical forcing -- which
  ! is what makes the replica checks in elmxxPerfCheckMod meaningful.
  !
  ! The climate is idealized, not realistic: a diurnal air-temperature cycle,
  ! sunlight from the true solar zenith angle, and an afternoon shower that
  ! falls as snow below freezing, so the hydrology and (on cold cells) snow
  ! paths are exercised. A recorded-DATM reader is the planned alternative
  ! behind the same entry point (elmxx_perf_forcing_set).
  !-----------------------------------------------------------------------

  use shr_kind_mod   , only : r8 => shr_kind_r8
  use shr_const_mod  , only : SHR_CONST_PI, SHR_CONST_TKFRZ, SHR_CONST_RDAIR, &
                              SHR_CONST_CPDAIR
  use shr_orb_mod    , only : shr_orb_cosz
  use elmxxForcingMod, only : forc_z, forc_topo, forc_u, forc_v, forc_ptem, &
                              forc_shum, forc_pbot, forc_tbot, forc_lwrad, &
                              forc_rainc, forc_rainl, forc_snowc, forc_snowl, &
                              forc_swndr, forc_swvdr, forc_swndf, forc_swvdf, &
                              forc_aer

  implicit none
  private

  public :: elmxx_perf_forcing_set

  real(r8), parameter :: z_ref    = 30._r8       ! reference height [m]
  real(r8), parameter :: pbot     = 98260._r8    ! surface pressure [Pa]
  real(r8), parameter :: t_mean   = 297._r8      ! daily-mean air temperature [K]
  real(r8), parameter :: t_amp    = 5._r8        ! diurnal half-range [K]
  real(r8), parameter :: shum     = 0.015_r8     ! specific humidity [kg/kg]
  real(r8), parameter :: lwdown   = 420._r8      ! downward longwave [W/m2]
  real(r8), parameter :: wind_u   = 2.0_r8       ! [m/s]
  real(r8), parameter :: wind_v   = 0.5_r8       ! [m/s]
  real(r8), parameter :: s0_eff   = 1361._r8 * 0.75_r8  ! TOA solar times a fixed transmissivity [W/m2]
  real(r8), parameter :: frac_dir = 0.7_r8       ! direct fraction of shortwave
  real(r8), parameter :: precip   = 2.e-4_r8     ! shower rate [kg/m2/s] (0.72 mm/h)
  real(r8), parameter :: shower_start = 14._r8   ! local solar hour
  real(r8), parameter :: shower_end   = 17._r8

contains

  !-----------------------------------------------------------------------
  subroutine elmxx_perf_forcing_set(ncell, lat_deg, lon_deg, calday, tod, declin)
    !
    ! Fill the forcing for every owned cell at one instant.
    !
    integer , intent(in) :: ncell
    real(r8), intent(in) :: lat_deg(:), lon_deg(:) ! per owned cell [degrees]
    real(r8), intent(in) :: calday                 ! 1-based calendar day of the instant
    integer , intent(in) :: tod                    ! seconds of the day, UTC
    real(r8), intent(in) :: declin                 ! solar declination [radians]

    real(r8), parameter :: deg2rad = SHR_CONST_PI / 180._r8
    integer  :: g
    real(r8) :: hour_local, coszen, sw, tbot, rain_or_snow

    do g = 1, ncell
       hour_local = modulo(real(tod, r8) / 3600._r8 + lon_deg(g) / 15._r8, 24._r8)

       ! Coldest near dawn, warmest at 15:00 local solar time.
       tbot = t_mean + t_amp * sin(2._r8 * SHR_CONST_PI * (hour_local - 9._r8) / 24._r8)

       coszen = shr_orb_cosz(calday, lat_deg(g) * deg2rad, lon_deg(g) * deg2rad, declin)
       sw = s0_eff * max(coszen, 0._r8)

       rain_or_snow = 0._r8
       if (hour_local >= shower_start .and. hour_local < shower_end) rain_or_snow = precip

       forc_z(g)     = z_ref
       forc_topo(g)  = 0._r8
       forc_u(g)     = wind_u
       forc_v(g)     = wind_v
       forc_tbot(g)  = tbot
       forc_pbot(g)  = pbot
       forc_ptem(g)  = tbot * (1.e5_r8 / pbot) ** (SHR_CONST_RDAIR / SHR_CONST_CPDAIR)
       forc_shum(g)  = shum
       forc_lwrad(g) = lwdown

       forc_swvdr(g) = 0.5_r8 * sw * frac_dir
       forc_swndr(g) = 0.5_r8 * sw * frac_dir
       forc_swvdf(g) = 0.5_r8 * sw * (1._r8 - frac_dir)
       forc_swndf(g) = 0.5_r8 * sw * (1._r8 - frac_dir)

       forc_rainc(g) = 0._r8
       forc_snowc(g) = 0._r8
       if (tbot > SHR_CONST_TKFRZ) then
          forc_rainl(g) = rain_or_snow
          forc_snowl(g) = 0._r8
       else
          forc_rainl(g) = 0._r8
          forc_snowl(g) = rain_or_snow
       end if

       forc_aer(g, :) = 0._r8   ! clean snow: no aerosol deposition
    end do

  end subroutine elmxx_perf_forcing_set

end module elmxxPerfForcingMod
