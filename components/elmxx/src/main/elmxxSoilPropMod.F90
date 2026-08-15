module elmxxSoilPropMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! The vertical soil grid, and the hydraulic properties derived from soil
  ! texture.
  !
  ! WHY THIS EXISTS. Eight of Stage 4's sixteen kernels read `watsat`,
  ! `watfc`, `sucsat` and `bsw`. Those are NOT surfdata fields -- surfdata
  ! carries sand, clay and organic matter, and ELM turns them into hydraulic
  ! properties through a pedotransfer function at initialization. Without this
  ! port those kernels would read zeros and produce plausible-looking garbage,
  ! so Stage 4 refuses to run them (elmxxKernelMod). This is the single
  ! biggest unblock available.
  !
  ! PROVENANCE, because every number here has to be traceable:
  !   vertical grid       ELM initVerticalMod
  !                       zsoi(j) = scalez*(exp(zecoeff*(j-0.5))-1)
  !                       scalez = 0.025, zecoeff = 0.50
  !   pedotransfer        ELM FuncPedotransferMod, **Cosby 1984 Table 5**.
  !                       Table 5 is what get_ipedof(0) selects via ipedof0 --
  !                       NOT Table 4, which the same module also implements
  !                       with different coefficients.
  !   organic blending    ELM SoilStateType, with om_frac = organic/organic_max
  !                       and organic_max = 130 kg/m3 from clm_params.
  !   field capacity      ELM SoilStateType: watfc defined where hk = 0.1 mm/day
  !
  ! Layers are `nlevgrnd` (15), not `nlevsoi` (10). That is deliberate and the
  ! subject of a real bug already fixed once: SoilTemperature's supercooling
  ! term indexes these over the full ground column, so truncating to nlevsoi
  ! divides by zero in layers 11-15 (STATUS.md E.1). Surfdata only supplies
  ! texture for nlevsoi layers; below that ELM carries the deepest soil
  ! layer's texture down, which is what `lev_src` does here.
  !-----------------------------------------------------------------------

  use shr_kind_mod    , only : r8 => shr_kind_r8
  use shr_sys_mod     , only : shr_sys_abort, shr_sys_flush
  use elmxxSpmdMod    , only : masterproc, iam
  use elmxxSubgridMod , only : num_columns
  use elmxxSurfaceStateMod, only : surface_state_built, col_pct_sand, &
                                   col_pct_clay, col_organic

  implicit none
  save
  private

  integer, parameter, public :: nlevsoi  = 10   ! layers surfdata supplies texture for
  integer, parameter, public :: nlevgrnd = 15   ! ground layers the kernels index
  integer, parameter, public :: nlevsno  =  5   ! maximum snow layers
  integer, parameter, public :: nlevtot  = nlevsno + nlevgrnd   ! 20, snow + ground

  ! Depth to bedrock. ELM reads it from surfdata's aveDTB when present and
  ! otherwise falls back to nlevsoi for every column; aveDTB is absent from the
  ! surfdata in use, so the ELM twins are already running that fallback.
  ! Varying depth to bedrock is out of scope (decision, 2026-08-15), so this is
  ! a parameter rather than a per-column array. Layers below it are bedrock and
  ! hold no water.
  integer, parameter, public :: nlevbed  = nlevsoi

  ! Vertical grid, shared by every column (ELM's is global too).
  real(r8), public :: zsoi (nlevgrnd)      ! node depth, m
  real(r8), public :: dzsoi(nlevgrnd)      ! layer thickness, m
  real(r8), public :: zisoi(0:nlevgrnd)    ! interface depth, m

  ! Derived per column and layer.
  real(r8), public, pointer :: watsat(:,:) => null()  ! (nc, nlevgrnd) porosity, v/v
  real(r8), public, pointer :: bsw   (:,:) => null()  ! (nc, nlevgrnd) Clapp-Hornberger b
  real(r8), public, pointer :: sucsat(:,:) => null()  ! (nc, nlevgrnd) sat. matric potential, mm
  real(r8), public, pointer :: hksat (:,:) => null()  ! (nc, nlevgrnd) sat. conductivity, mm/s
  real(r8), public, pointer :: watfc (:,:) => null()  ! (nc, nlevgrnd) field capacity, v/v

  !--------------------------------------------------------------------------
  ! Thermal properties, from the same ELM SoilStateType block. Split out here
  ! rather than folded above because only SoilTemperature reads them, and
  ! keeping them separate makes it obvious which kernel required them.
  !--------------------------------------------------------------------------
  real(r8), public, pointer :: tkmg  (:,:) => null()  ! (nc, nlevgrnd) dry-solid conductivity, W/m/K
  real(r8), public, pointer :: tkdry (:,:) => null()  ! (nc, nlevgrnd) dry-soil conductivity, W/m/K
  real(r8), public, pointer :: tksatu(:,:) => null()  ! (nc, nlevgrnd) saturated conductivity, W/m/K
  real(r8), public, pointer :: csol  (:,:) => null()  ! (nc, nlevgrnd) heat capacity, J/m3/K

  !--------------------------------------------------------------------------
  ! Cold-start column state, on the snow+ground index space the kernels use.
  ! Slot m = 1..nlevtot maps to ELM layer j = m - nlevsno, so m = 6 is ELM's
  ! first soil layer and m = 1..5 are the snow slots.
  !--------------------------------------------------------------------------
  real(r8), public, pointer :: col_dz        (:,:) => null()  ! (nc, nlevtot) m
  real(r8), public, pointer :: col_t_soisno  (:,:) => null()  ! (nc, nlevtot) K
  real(r8), public, pointer :: col_h2osoi_liq(:,:) => null()  ! (nc, nlevtot) kg/m2
  real(r8), public, pointer :: col_h2osoi_ice(:,:) => null()  ! (nc, nlevtot) kg/m2
  real(r8), public, pointer :: col_h2osoi_vol(:,:) => null()  ! (nc, nlevgrnd) v/v

  logical, public :: soil_prop_built = .false.

  public :: elmxx_soil_prop_init
  public :: elmxx_soil_prop_clean

contains

  !-----------------------------------------------------------------------
  subroutine elmxx_soil_prop_init(logunit)
    implicit none
    integer, intent(in) :: logunit
    character(len=*), parameter :: subname = '(elmxx_soil_prop_init) '

    if (.not. surface_state_built) then
       call shr_sys_abort(subname//'ERROR: surface state is not ready')
    end if

    call elmxx_soil_prop_clean()
    call build_vertical_grid()
    call build_hydraulic_properties(logunit)
    call build_cold_start_state(logunit)

    soil_prop_built = .true.

  end subroutine elmxx_soil_prop_init

  !-----------------------------------------------------------------------
  subroutine build_vertical_grid()
    !
    ! ELM initVerticalMod's default exponential grid. Note dzsoi's first and
    ! last layers are special cases, not the interior formula.
    !
    implicit none
    real(r8), parameter :: scalez  = 0.025_r8
    real(r8), parameter :: zecoeff = 0.50_r8
    integer :: j

    do j = 1, nlevgrnd
       zsoi(j) = scalez * (exp(zecoeff*(real(j,r8) - 0.5_r8)) - 1.0_r8)
    end do

    dzsoi(1) = 0.5_r8 * (zsoi(1) + zsoi(2))
    do j = 2, nlevgrnd-1
       dzsoi(j) = 0.5_r8 * (zsoi(j+1) - zsoi(j-1))
    end do
    dzsoi(nlevgrnd) = zsoi(nlevgrnd) - zsoi(nlevgrnd-1)

    zisoi(0) = 0.0_r8
    do j = 1, nlevgrnd-1
       zisoi(j) = 0.5_r8 * (zsoi(j) + zsoi(j+1))
    end do
    zisoi(nlevgrnd) = zsoi(nlevgrnd) + 0.5_r8 * dzsoi(nlevgrnd)

  end subroutine build_vertical_grid

  !-----------------------------------------------------------------------
  subroutine build_hydraulic_properties(logunit)
    implicit none
    integer, intent(in) :: logunit
    integer  :: c, j, jsrc
    real(r8) :: sand, clay, om_frac
    real(r8) :: wsat_min, b_min, suc_min, xksat
    real(r8) :: om_watsat, om_b, om_sucsat, om_hksat
    real(r8) :: perc_norm, perc_frac, uncon_frac, uncon_hksat
    real(r8) :: bd, tkm
    ! ELM SoilStateType
    real(r8), parameter :: organic_max = 130.0_r8   ! clm_params, kg/m3
    real(r8), parameter :: zsapric     = 0.5_r8     ! m
    real(r8), parameter :: pcalpha     = 0.5_r8     ! percolation threshold
    real(r8), parameter :: pcbeta      = 0.139_r8   ! percolation exponent
    real(r8), parameter :: secspday    = 86400.0_r8
    ! Thermal, ELM SoilStateType
    real(r8), parameter :: om_tkm       = 0.25_r8    ! organic conductivity, W/m/K
    real(r8), parameter :: om_tkd       = 0.05_r8    ! dry organic conductivity
    real(r8), parameter :: om_csol      = 2.5_r8     ! peat heat capacity, *1e6 J/K/m3
    real(r8), parameter :: csol_bedrock = 2.0e6_r8   ! granite/sandstone, J/m3/K
    character(len=*), parameter :: subname = '(elmxx_soil_prop_init) '

    allocate(watsat(num_columns, nlevgrnd), bsw(num_columns, nlevgrnd), &
             sucsat(num_columns, nlevgrnd), hksat(num_columns, nlevgrnd), &
             watfc(num_columns, nlevgrnd), tkmg(num_columns, nlevgrnd), &
             tkdry(num_columns, nlevgrnd), tksatu(num_columns, nlevgrnd), &
             csol(num_columns, nlevgrnd))

    do c = 1, num_columns
       do j = 1, nlevgrnd

          ! Surfdata carries texture for nlevsoi layers only; below that ELM
          ! carries the deepest soil layer's texture down the ground column.
          jsrc = min(j, nlevsoi)
          sand = col_pct_sand(c, jsrc)
          clay = col_pct_clay(c, jsrc)
          om_frac = min(col_organic(c, jsrc) / organic_max, 1.0_r8)

          ! --- mineral soil, Cosby 1984 Table 5 (ipedof0) ---
          wsat_min = 0.489_r8 - 0.00126_r8*sand
          b_min    = 2.91_r8  + 0.159_r8 *clay
          suc_min  = 10.0_r8 * (10.0_r8**(1.88_r8  - 0.0131_r8*sand))
          xksat    = 0.0070556_r8 * (10.0_r8**(-0.884_r8 + 0.0153_r8*sand))

          ! --- organic end members, depth-dependent toward sapric peat ---
          om_watsat = max(0.93_r8 - 0.1_r8   *(zsoi(j)/zsapric), 0.83_r8)
          om_b      = min(2.7_r8  + 9.3_r8   *(zsoi(j)/zsapric), 12.0_r8)
          om_sucsat = min(10.3_r8 - 0.2_r8   *(zsoi(j)/zsapric), 10.1_r8)
          om_hksat  = max(0.28_r8 - 0.2799_r8*(zsoi(j)/zsapric), 0.0001_r8)

          ! --- blend mineral and organic ---
          watsat(c,j) = (1.0_r8 - om_frac)*wsat_min + om_watsat*om_frac
          bsw   (c,j) = (1.0_r8 - om_frac)*b_min    + om_b     *om_frac
          sucsat(c,j) = (1.0_r8 - om_frac)*suc_min  + om_sucsat*om_frac

          ! --- conductivity: mineral and organic in series, with a
          !     percolating organic fraction shorting past them ---
          if (om_frac > pcalpha) then
             perc_norm = (1.0_r8 - pcalpha)**(-pcbeta)
             perc_frac = perc_norm * (om_frac - pcalpha)**pcbeta
          else
             perc_frac = 0.0_r8
          end if
          uncon_frac = (1.0_r8 - om_frac) + (1.0_r8 - perc_frac)*om_frac

          if (om_frac < 1.0_r8) then
             uncon_hksat = uncon_frac / ((1.0_r8 - om_frac)/xksat &
                         + ((1.0_r8 - perc_frac)*om_frac)/om_hksat)
          else
             uncon_hksat = 0.0_r8
          end if
          hksat(c,j) = uncon_frac*uncon_hksat + (perc_frac*om_frac)*om_hksat

          ! --- field capacity: the water content at which hk = 0.1 mm/day ---
          watfc(c,j) = watsat(c,j) * &
               (0.1_r8 / (hksat(c,j)*secspday))**(1.0_r8/(2.0_r8*bsw(c,j) + 3.0_r8))

          ! --- thermal properties (ELM SoilStateType) ---
          ! bd is bulk density, and the (sand+clay) denominator is why a
          ! column with neither would divide by zero -- surfdata always has
          ! one or the other, but guard rather than assume.
          bd  = (1.0_r8 - watsat(c,j)) * 2.7e3_r8
          if (sand + clay <= 0.0_r8) then
             call shr_sys_abort(subname//'ERROR: column has neither sand nor clay')
          end if
          tkm = (1.0_r8 - om_frac)*(8.80_r8*sand + 2.92_r8*clay)/(sand + clay) &
              + om_tkm*om_frac
          tkmg  (c,j) = tkm ** (1.0_r8 - watsat(c,j))
          tksatu(c,j) = tkmg(c,j) * 0.57_r8**watsat(c,j)
          tkdry (c,j) = ((0.135_r8*bd + 64.7_r8) / (2.7e3_r8 - 0.947_r8*bd)) &
                      * (1.0_r8 - om_frac) + om_tkd*om_frac
          csol  (c,j) = ((1.0_r8 - om_frac)*(2.128_r8*sand + 2.385_r8*clay)/(sand + clay) &
                      + om_csol*om_frac) * 1.0e6_r8
          if (j > nlevbed) csol(c,j) = csol_bedrock

       end do
    end do

    call report(logunit)

  end subroutine build_hydraulic_properties

  !-----------------------------------------------------------------------
  subroutine build_cold_start_state(logunit)
    !
    ! ELM's ColumnDataType InitCold, for natural soil columns.
    !
    ! WHICH BRANCH, AND WHY IT IS THIS ONE. ELM's h2osoi_vol cold start forks
    ! on FATES/hydrstress, arctic init, landunit type and bedrock depth. In
    ! this configuration:
    !   use_fates, use_hydrstress, use_arctic_init   all .false. (checked in
    !                                                the ELM twin's lnd_in)
    !   varying depth to bedrock                     out of scope; nlevbed is
    !                                                nlevsoi, which is also
    !                                                ELM's own fallback when
    !                                                aveDTB is absent -- and it
    !                                                is absent from the
    !                                                surfdata in use
    !   wetland, glacier, lake                       out of scope or unpacked
    ! What is left is the plain branch: 0.15 by volume above bedrock, zero
    ! below, capped at porosity.
    !
    ! Urban columns are NOT handled here (decision, 2026-08-15: natural first).
    ! Their cold start differs per column type -- pervious road 0.3, impervious
    ! road 0, roof and walls zero over nlevurb -- and belongs with the urban
    ! increment.
    !
    ! Liquid versus ice is decided by temperature, and at 274 K every layer is
    ! liquid. The branch is kept anyway: it costs nothing and the alternative
    ! is code that silently assumes a warm start.
    !
    implicit none
    integer, intent(in) :: logunit
    integer  :: c, m, j
    real(r8) :: vol
    real(r8), parameter :: t_soil_cold = 274.0_r8      ! ELM InitCold, non-lake
    real(r8), parameter :: h2osoi_vol_cold = 0.15_r8   ! ELM InitCold, plain branch
    real(r8), parameter :: tkfrz  = 273.15_r8          ! SHR_CONST_TKFRZ
    real(r8), parameter :: denh2o = 1000.0_r8          ! kg/m3
    real(r8), parameter :: denice =  917.0_r8          ! kg/m3

    allocate(col_dz        (num_columns, nlevtot), &
             col_t_soisno  (num_columns, nlevtot), &
             col_h2osoi_liq(num_columns, nlevtot), &
             col_h2osoi_ice(num_columns, nlevtot), &
             col_h2osoi_vol(num_columns, nlevgrnd))

    ! Snow slots stay zero: snl = 0 at a cold start, so no snow layer exists
    ! and every kernel gates on snl before reading them.
    col_dz = 0.0_r8; col_t_soisno = 0.0_r8
    col_h2osoi_liq = 0.0_r8; col_h2osoi_ice = 0.0_r8; col_h2osoi_vol = 0.0_r8

    do c = 1, num_columns
       do j = 1, nlevgrnd
          m = j + nlevsno                       ! ELM layer j -> packed slot m

          col_dz(c,m)       = dzsoi(j)
          col_t_soisno(c,m) = t_soil_cold

          if (j > nlevbed) then
             vol = 0.0_r8                       ! bedrock holds no water
          else
             vol = min(h2osoi_vol_cold, watsat(c,j))
          end if
          col_h2osoi_vol(c,j) = vol

          if (col_t_soisno(c,m) <= tkfrz) then
             col_h2osoi_ice(c,m) = col_dz(c,m) * denice * vol
             col_h2osoi_liq(c,m) = 0.0_r8
          else
             col_h2osoi_ice(c,m) = 0.0_r8
             col_h2osoi_liq(c,m) = col_dz(c,m) * denh2o * vol
          end if
       end do
    end do

    write(logunit,*) '(elmxx_soil_prop_init) rank ',iam,' cold-start column state:'
    write(logunit,*) '    h2osoi_vol [v/v]   ',minval(col_h2osoi_vol),' .. ',maxval(col_h2osoi_vol)
    write(logunit,*) '    h2osoi_liq [kg/m2] ',minval(col_h2osoi_liq),' .. ',maxval(col_h2osoi_liq)
    write(logunit,*) '    h2osoi_ice [kg/m2] ',minval(col_h2osoi_ice),' .. ',maxval(col_h2osoi_ice)
    write(logunit,*) '    t_soisno   [K]     ',t_soil_cold,' over ',nlevgrnd,' ground layers'
    write(logunit,*) '    nlevbed            ',nlevbed,' (layers below hold no water)'
    call shr_sys_flush(logunit)

    if (any(col_h2osoi_vol > watsat)) then
       call shr_sys_abort('(elmxx_soil_prop_init) ERROR: soil water exceeds porosity')
    end if

  end subroutine build_cold_start_state

  !-----------------------------------------------------------------------
  subroutine report(logunit)
    !
    ! Ranges, plus the bounds that must hold whatever the texture is. A
    ! porosity outside (0,1) or a field capacity above saturation is a broken
    ! pedotransfer, and saying so here is cheaper than tracing it out of a
    ! kernel later.
    !
    implicit none
    integer, intent(in) :: logunit
    character(len=*), parameter :: subname = '(elmxx_soil_prop_init) '

    write(logunit,*) subname,'rank ',iam,' derived soil properties over ', &
                     num_columns,' columns x ',nlevgrnd,' layers'
    ! Full arrays, in ELM's initVerticalMod log format, so the two can be
    ! diffed directly. They matched digit for digit on 2x1_brazil, which is
    ! the only independent check available for the grid -- it is time-constant
    ! and so appears in no restart file.
    write(logunit,*) '    zsoi:  ',zsoi(:)
    write(logunit,*) '    zisoi: ',zisoi(:)
    write(logunit,*) '    dzsoi: ',dzsoi(:)
    write(logunit,*) '    watsat [v/v]  ',minval(watsat),' .. ',maxval(watsat)
    write(logunit,*) '    bsw    [-]    ',minval(bsw)   ,' .. ',maxval(bsw)
    write(logunit,*) '    sucsat [mm]   ',minval(sucsat),' .. ',maxval(sucsat)
    write(logunit,*) '    hksat  [mm/s] ',minval(hksat) ,' .. ',maxval(hksat)
    write(logunit,*) '    watfc  [v/v]  ',minval(watfc) ,' .. ',maxval(watfc)
    write(logunit,*) '    tkmg   [W/m/K]',minval(tkmg)  ,' .. ',maxval(tkmg)
    write(logunit,*) '    tkdry  [W/m/K]',minval(tkdry) ,' .. ',maxval(tkdry)
    write(logunit,*) '    tksatu [W/m/K]',minval(tksatu),' .. ',maxval(tksatu)
    write(logunit,*) '    csol   [J/m3/K]',minval(csol) ,' .. ',maxval(csol)
    call shr_sys_flush(logunit)

    if (minval(watsat) <= 0.0_r8 .or. maxval(watsat) >= 1.0_r8) then
       call shr_sys_abort(subname//'ERROR: porosity outside (0,1)')
    end if
    if (minval(bsw) <= 0.0_r8) then
       call shr_sys_abort(subname//'ERROR: non-positive Clapp-Hornberger b')
    end if
    if (minval(sucsat) <= 0.0_r8) then
       call shr_sys_abort(subname//'ERROR: non-positive saturated matric potential')
    end if
    if (any(watfc > watsat)) then
       call shr_sys_abort(subname//'ERROR: field capacity exceeds saturation')
    end if
    ! Saturated soil conducts heat better than dry soil, always. If this
    ! inverts, tkmg's exponent or watsat is wrong.
    if (any(tksatu < tkdry)) then
       call shr_sys_abort(subname//'ERROR: saturated conductivity below dry conductivity')
    end if
    if (minval(csol) <= 0.0_r8) then
       call shr_sys_abort(subname//'ERROR: non-positive volumetric heat capacity')
    end if

  end subroutine report

  !-----------------------------------------------------------------------
  subroutine elmxx_soil_prop_clean()
    implicit none
    if (associated(watsat)) deallocate(watsat)
    if (associated(bsw))    deallocate(bsw)
    if (associated(sucsat)) deallocate(sucsat)
    if (associated(hksat))  deallocate(hksat)
    if (associated(watfc))  deallocate(watfc)
    if (associated(tkmg))   deallocate(tkmg)
    if (associated(tkdry))  deallocate(tkdry)
    if (associated(tksatu)) deallocate(tksatu)
    if (associated(csol))   deallocate(csol)
    if (associated(col_dz))         deallocate(col_dz)
    if (associated(col_t_soisno))   deallocate(col_t_soisno)
    if (associated(col_h2osoi_liq)) deallocate(col_h2osoi_liq)
    if (associated(col_h2osoi_ice)) deallocate(col_h2osoi_ice)
    if (associated(col_h2osoi_vol)) deallocate(col_h2osoi_vol)
    watsat => null(); bsw => null(); sucsat => null()
    hksat  => null(); watfc => null()
    tkmg => null(); tkdry => null(); tksatu => null(); csol => null()
    col_dz => null(); col_t_soisno => null()
    col_h2osoi_liq => null(); col_h2osoi_ice => null(); col_h2osoi_vol => null()
    soil_prop_built = .false.
  end subroutine elmxx_soil_prop_clean

end module elmxxSoilPropMod
