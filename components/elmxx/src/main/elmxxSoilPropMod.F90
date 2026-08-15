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
    ! ELM SoilStateType
    real(r8), parameter :: organic_max = 130.0_r8   ! clm_params, kg/m3
    real(r8), parameter :: zsapric     = 0.5_r8     ! m
    real(r8), parameter :: pcalpha     = 0.5_r8     ! percolation threshold
    real(r8), parameter :: pcbeta      = 0.139_r8   ! percolation exponent
    real(r8), parameter :: secspday    = 86400.0_r8
    character(len=*), parameter :: subname = '(elmxx_soil_prop_init) '

    allocate(watsat(num_columns, nlevgrnd), bsw(num_columns, nlevgrnd), &
             sucsat(num_columns, nlevgrnd), hksat(num_columns, nlevgrnd), &
             watfc(num_columns, nlevgrnd))

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

       end do
    end do

    call report(logunit)

  end subroutine build_hydraulic_properties

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

  end subroutine report

  !-----------------------------------------------------------------------
  subroutine elmxx_soil_prop_clean()
    implicit none
    if (associated(watsat)) deallocate(watsat)
    if (associated(bsw))    deallocate(bsw)
    if (associated(sucsat)) deallocate(sucsat)
    if (associated(hksat))  deallocate(hksat)
    if (associated(watfc))  deallocate(watfc)
    watsat => null(); bsw => null(); sucsat => null()
    hksat  => null(); watfc => null()
    soil_prop_built = .false.
  end subroutine elmxx_soil_prop_clean

end module elmxxSoilPropMod
