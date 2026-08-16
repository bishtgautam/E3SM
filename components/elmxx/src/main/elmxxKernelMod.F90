module elmxxKernelMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Stage 4: activate ELMxx's Kokkos kernels one at a time.
  !
  ! `elmxx_kernels` is a comma-separated list of the short names below; a
  ! kernel not named simply does not run. THERE IS NO FORTRAN ALTERNATIVE TO
  ! FALL BACK ON -- ELMxx is written from scratch -- so the timestep is
  ! incomplete until every kernel is on, and intermediate runs are diagnostic
  ! only. That is by design (plan, Stage 4), not a defect.
  !
  ! PREREQUISITES ARE ENFORCED, NOT ASSUMED. Every kernel reads state that
  ! something else has to have produced: another kernel, Stage 3's seeding, or
  ! a port that does not exist yet. A kernel run on unset inputs does not
  ! crash -- it reads zeros and produces plausible-looking garbage, which is
  ! the worst possible failure mode and exactly what Stage 2 and 3's checks
  ! were built to avoid. So each kernel carries the reason it is or is not
  ! runnable today, and asking for a blocked one aborts at init with that
  ! reason rather than at some later point with wrong physics.
  !
  ! Update `blocked_reason` as ports land. It is the honest record of what
  ! Stage 4 is waiting on.
  !-----------------------------------------------------------------------

  use shr_kind_mod    , only : r8 => shr_kind_r8
  use shr_sys_mod     , only : shr_sys_abort, shr_sys_flush
  use elmxxSpmdMod    , only : masterproc, iam
  use elmxx_mod       , only : ELMxxType, ELMXX_SUCCESS, &
                               ELMxxComputeSurfaceRadiation, &
                               ELMxxComputeCanopyHydrology, &
                               ELMxxComputeCanopySunShadeFractions, &
                               ELMxxComputeCanopyTemperature, &
                               ELMxxComputeBareGroundFluxes, &
                               ELMxxComputeCanopyFluxes, &
                               ELMxxComputeUrbanRadiation, &
                               ELMxxComputeUrbanFluxes, &
                               ELMxxComputeLakeFluxes, &
                               ELMxxComputeLakeTemperature, &
                               ELMxxComputeSoilTemperature, &
                               ELMxxComputeSoilFluxes, &
                               ELMxxComputeSurfRunInfil, &
                               ELMxxComputeRootWaterUpdate, &
                               ELMxxComputeHydrologyDrainage, &
                               ELMxxComputeLakeHydrology, &
                               ELMxxGetQflxPrecIntr, ELMxxGetQflxPrecGrnd, &
                               ELMxxGetH2ocan, ELMxxGetFwet, ELMxxGetFdry, &
                               ELMxxGetTGrnd, ELMxxGetQg, ELMxxGetThv, &
                               ELMxxGetHtvp, ELMxxGetSoilbeta, ELMxxGetZ0mg, &
                               ELMxxGetEflxShGrnd, ELMxxGetEflxShVeg, &
                               ELMxxGetQflxEvapSoi, ELMxxGetQflxTranVeg, &
                               ELMxxGetTVeg, ELMxxGetBtran, &
                               ELMxxGetFsa, ELMxxGetFsr, ELMxxGetSabv, &
                               ELMxxGetSabg, ELMxxGetSabgSoil

  implicit none
  save
  private

  ! Driver order. This is the order the kernels run in, and the order they
  ! should be activated in. It is ELM's driver order, not alphabetical.
  integer, parameter, public :: NKERNEL = 16

  integer, parameter, public :: K_SURFRAD     =  1
  integer, parameter, public :: K_CANHYDRO    =  2
  integer, parameter, public :: K_CANSUNSHADE =  3
  integer, parameter, public :: K_CANTEMP     =  4
  integer, parameter, public :: K_BAREGRND    =  5
  integer, parameter, public :: K_CANFLUX     =  6
  integer, parameter, public :: K_URBANRAD    =  7
  integer, parameter, public :: K_URBANFLUX   =  8
  integer, parameter, public :: K_LAKEFLUX    =  9
  integer, parameter, public :: K_LAKETEMP    = 10
  integer, parameter, public :: K_SOILTEMP    = 11
  integer, parameter, public :: K_SOILFLUX    = 12
  integer, parameter, public :: K_SURFRUNOFF  = 13
  integer, parameter, public :: K_ROOTWATER   = 14
  integer, parameter, public :: K_HYDRODRAIN  = 15
  integer, parameter, public :: K_LAKEHYDRO   = 16

  character(len=16), parameter, public :: kernel_name(NKERNEL) = [ &
       'surfrad         ', 'canhydro        ', 'cansunshade     ', &
       'cantemp         ', 'baregrnd        ', 'canflux         ', &
       'urbanrad        ', 'urbanflux       ', 'lakeflux        ', &
       'laketemp        ', 'soiltemp        ', 'soilflux        ', &
       'surfrunoff      ', 'rootwater       ', 'hydrodrain      ', &
       'lakehydro       ' ]

  logical, public :: kernel_active(NKERNEL) = .false.
  logical, public :: any_kernel_active      = .false.

  public :: elmxx_kernels_parse
  public :: elmxx_kernels_run
  public :: elmxx_kernels_report
  public :: elmxx_report_cantemp
  public :: elmxx_report_fluxes
  public :: elmxx_report_surfrad

contains

  !-----------------------------------------------------------------------
  function blocked_reason(k) result(why)
    !
    ! Why kernel k cannot run yet, or ' ' if it can.
    !
    ! Derived by reading each kernel's Impl header for the views it reads and
    ! comparing against what Stage 2 seeds, Stage 3 crosses, and earlier
    ! kernels produce. Recorded here rather than in a comment so that asking
    ! for a blocked kernel fails loudly with the reason.
    !
    implicit none
    integer, intent(in) :: k
    character(len=256)  :: why

    why = ' '

    select case (k)

    case (K_SURFRAD)
       ! Runnable on ELM's cold-start albedos.
       !
       ! SurfaceAlbedo runs at the END of ELM's timestep, so the albedos
       ! SurfaceRadiation reads at step N came from step N-1; the dependency
       ! is temporal, not within-step. ELM bootstraps step one from
       ! SurfaceAlbedoType InitCold, and elmxx_kokkos_seed_albedo supplies the
       ! same constants. Also needs the 2-D forc_solad/forc_solai, crossed
       ! with the rest of the forcing.
       !
       ! LIMIT: with albedos frozen there is no solar-zenith dependence, no
       ! snow aging, no wetness effect. Bounded and gradeable, but a multi-day
       ! run is not physically right until SurfaceAlbedo is ported.
       why = ' '

    case (K_CANSUNSHADE)
       ! Both read the full canopy/ground albedo set -- albd, albi, fabd,
       ! fabi, ftdd, ftid, ftii, albgrd, albgri, albso*, albsn*_hst, and the
       ! SNICAR flx_abs* factors. Every one of those is SurfaceAlbedo's
       ! output, and SurfaceAlbedo is a Stage 5 port that does not exist.
       ! They also read forc_solad/forc_solai, which Stage 3 does not cross.
       ! Needs the sunlit/shaded canopy decomposition -- fabd_sun_z,
       ! fabd_sha_z, fabi_sun_z, fabi_sha_z, tlai_z, nrad -- which is
       ! SurfaceAlbedo's per-canopy-layer output. Unlike the bulk albedos
       ! there is no meaningful cold-start constant for a vertical profile,
       ! so this one genuinely waits for the port.
       why = 'needs SurfaceAlbedo for the per-canopy-layer sunlit/shaded ' // &
             'decomposition (fabd_sun_z, fabi_sun_z, tlai_z, nrad)'

    case (K_CANTEMP)
       ! Runnable. Everything it reads is now seeded: the hydraulic properties
       ! and cold-start column state from elmxxSoilPropMod, and the four
       ! scalars it reads without writing -- smpmin, t_h2osfc, patch_itype and
       ! forc_hgt_patch. What it appears to want beyond that -- zii, the three
       ! forc_hgt_*_patch, t_ssbef, t_h2osfc_bef -- it computes itself, as ELM
       ! does; checking that rather than assuming saved seeding four fields
       ! the kernel would have overwritten.
       why = ' '

    case (K_CANFLUX)
       ! btran now has a source: elmxxRootMod computes it from rootfr and soil
       ! matric potential, and it crosses each step.
       !
       ! IT WILL BE ZERO, AND THAT IS CORRECT. The cold-start soil is 0.15 by
       ! volume, about 28% of saturation, which at bsw ~ 8 gives a matric
       ! potential near -4.7e6 mm -- roughly eighteen times drier than these
       ! PFTs' closure threshold of -2.55e5. The smp_node clamp pins it to
       ! smpsc, rresis goes to zero, and the plant transpires nothing. ELM
       ! would do the same from the same cold start.
       !
       ! So this kernel runs but cannot be meaningfully GRADED yet: nothing
       ! wets the soil until the hydrology kernels are active, so btran stays
       ! zero and every canopy flux with it. Structurally exercised,
       ! numerically trivial -- worth knowing before reading its output as
       ! evidence of anything.
       why = ' '

    case (K_BAREGRND)
       ! Their inputs are CanopyTemperature's outputs -- qg, thv, htvp, the
       ! roughness lengths, soilbeta, displa, thm -- plus forc_rho_col, which
       ! is now derived from vapor pressure and crossed each step. Runnable
       ! provided cantemp runs first, which driver order guarantees.
       !
       ! ugust is NOT provided: gustiness is out of scope for now (decision,
       ! 2026-08-15). BareGroundFluxes reads it, so it will see zero, which is
       ! ELM's no-gustiness case rather than an unset value -- but it is worth
       ! confirming against the kernel before trusting its fluxes.
       why = ' '

    case (K_SOILTEMP, K_SOILFLUX, K_SURFRUNOFF, K_ROOTWATER, K_HYDRODRAIN)
       ! CORRECTED 2026-08-15. An earlier note here claimed these "cannot
       ! usefully precede SurfaceAlbedo". That conflated two different things
       ! and is withdrawn.
       !
       ! What is true: SoilTemperature is driven by hs_soil / hs_top_snow /
       ! hs_h2osfc, the ground surface energy balance, which ELM builds in
       ! SoilTemperatureMod as
       !     eflx_gnet_soil = sabg_soil + dlrad
       !                    + (1-frac_veg_nosno)*emg*forc_lwrad
       !                    - lwrad_emit_soil
       !                    - (eflx_sh_soil + qflx_ev_soil*htvp)
       ! so they do need sabg, and sabg comes from SurfaceRadiation. That is a
       ! dependency on the SURFRAD KERNEL, which runs earlier in this same
       ! timestep and is now active.
       !
       ! What was false: that it is a dependency on the SurfaceAlbedo PORT.
       ! SurfaceAlbedo runs at the END of ELM's timestep, so SurfaceRadiation
       ! never consumes an albedo computed in its own step -- it reads the
       ! previous step's, and step one reads InitCold's constants. Frozen
       ! albedos cost realism, not runnability.
       !
       ! So the only thing left blocking these five is the integration
       ! surface, and it is a different one entirely. These are SHARED kernels:
       ! they do not read naturalCol, they read their own per-kernel state
       ! seeded through ST_/SF_/SRI_/RWU_/HD_ setters -- 163 of them -- on top
       ! of a topology declared by ELMxxInitSharedMetadata with the nolakec,
       ! nolakep, hydrologyc and urbanc filters and an urbpoi flag.
       ! elmxxFilterMod already builds all four filters, so the Fortran side
       ! fits; the seeding does not exist yet.
       why = 'needs ELMxxInitSharedMetadata, the shared filters, and ' // &
             'per-kernel ST_/SF_/SRI_/RWU_/HD_ seeding'

    case (K_URBANRAD, K_URBANFLUX)
       why = 'needs UrbanAlbedo for sabs_dir/sabs_dif, which is part of the ' // &
             'Stage 5 albedo port'

    case (K_LAKEFLUX, K_LAKETEMP, K_LAKEHYDRO)
       ! The packed Kokkos views carry no lake at all (see
       ! elmxxKokkosStateMod): ELMxxCreate takes natural columns, natural
       ! patches and urban landunits only. Lake state is a separate
       ! allocation (ELMxxAllocateLakeState) that Stage 3 does not make.
       why = 'needs the lake state allocation and its own packed maps; ' // &
             'ELMxxCreate carries natural and urban only'

    case (K_CANHYDRO)
       ! The closest to runnable, and the real first kernel of Stage 4 -- not
       ! surfrad, whatever the plan's driver order says. Its forcing (forc_rain,
       ! forc_snow, forc_t_col) is already crossed and its phenology (elai,
       ! esai) already seeded; most of the rest is zero at a cold start.
       ! What is missing is small and specific: micro_sigma and n_melt, which
       ! ELM derives in initVerticalMod from SLOPE and STD_ELEV -- both present
       ! on surfdata, neither yet read by elmxxSurfdataMod.
       why = ' '

    case default
       why = 'unknown kernel'

    end select

  end function blocked_reason

  !-----------------------------------------------------------------------
  subroutine elmxx_kernels_parse(spec, logunit)
    !
    ! Turn the namelist string into the active flags.
    !
    ! An unrecognised or blocked name aborts. Silently ignoring a name the
    ! user asked for would mean a run that looks like it exercised a kernel
    ! and did not.
    !
    implicit none
    character(len=*), intent(in) :: spec
    integer, intent(in) :: logunit
    integer :: i, k, n, ib, ie
    character(len=16) :: token
    character(len=256) :: why
    logical :: matched
    character(len=*), parameter :: subname = '(elmxx_kernels_parse) '

    kernel_active = .false.
    any_kernel_active = .false.

    n = len_trim(spec)
    if (n == 0) then
       if (masterproc) then
          write(logunit,*) subname,'no kernels active; the timestep is a no-op'
          call shr_sys_flush(logunit)
       end if
       return
    end if

    ib = 1
    do while (ib <= n)
       ie = index(spec(ib:n), ',')
       if (ie == 0) then
          ie = n
       else
          ie = ib + ie - 2
       end if

       token = adjustl(spec(ib:ie))
       if (len_trim(token) > 0) then
          matched = .false.
          do k = 1, NKERNEL
             if (trim(token) == trim(kernel_name(k))) then
                why = blocked_reason(k)
                if (len_trim(why) > 0) then
                   write(logunit,*) subname,'ERROR: kernel "',trim(token), &
                        '" cannot run yet: ',trim(why)
                   call shr_sys_flush(logunit)
                   call shr_sys_abort(subname//'ERROR: kernel "'//trim(token)// &
                        '" has unmet prerequisites')
                end if
                kernel_active(k) = .true.
                matched = .true.
                exit
             end if
          end do

          if (.not. matched) then
             write(logunit,*) subname,'ERROR: unknown kernel "',trim(token),'"'
             write(logunit,*) subname,'valid names, in driver order:'
             do i = 1, NKERNEL
                write(logunit,*) '    ',trim(kernel_name(i))
             end do
             call shr_sys_flush(logunit)
             call shr_sys_abort(subname//'ERROR: unknown kernel "'//trim(token)//'"')
          end if
       end if

       ib = ie + 2
    end do

    any_kernel_active = any(kernel_active)

    if (masterproc) then
       write(logunit,*) subname,'active kernels, in driver order:'
       do k = 1, NKERNEL
          if (kernel_active(k)) write(logunit,*) '    ',trim(kernel_name(k))
       end do
       call shr_sys_flush(logunit)
    end if

  end subroutine elmxx_kernels_parse

  !-----------------------------------------------------------------------
  subroutine elmxx_kernels_run(elm, dtime, logunit)
    !
    ! Dispatch the active kernels in driver order.
    !
    ! Order is the array order, never the order they were named on the
    ! namelist -- a kernel reads what the ones before it wrote, so letting the
    ! namelist reorder them would silently change the physics.
    !
    implicit none
    type(ELMxxType), intent(in) :: elm
    real(r8), intent(in) :: dtime
    integer, intent(in) :: logunit
    integer :: ierr
    character(len=*), parameter :: subname = '(elmxx_kernels_run) '

    if (.not. any_kernel_active) return

    if (kernel_active(K_SURFRAD)) then
       call ELMxxComputeSurfaceRadiation(elm, ierr)
       call check(ierr, logunit, K_SURFRAD)
    end if

    if (kernel_active(K_CANHYDRO)) then
       call ELMxxComputeCanopyHydrology(elm, dtime, ierr)
       call check(ierr, logunit, K_CANHYDRO)
    end if

    if (kernel_active(K_CANSUNSHADE)) then
       call ELMxxComputeCanopySunShadeFractions(elm, ierr)
       call check(ierr, logunit, K_CANSUNSHADE)
    end if

    if (kernel_active(K_CANTEMP)) then
       call ELMxxComputeCanopyTemperature(elm, ierr)
       call check(ierr, logunit, K_CANTEMP)
    end if

    if (kernel_active(K_BAREGRND)) then
       call ELMxxComputeBareGroundFluxes(elm, ierr)
       call check(ierr, logunit, K_BAREGRND)
    end if

    if (kernel_active(K_CANFLUX)) then
       call ELMxxComputeCanopyFluxes(elm, dtime, ierr)
       call check(ierr, logunit, K_CANFLUX)
    end if

    if (kernel_active(K_URBANRAD)) then
       call ELMxxComputeUrbanRadiation(elm, ierr)
       call check(ierr, logunit, K_URBANRAD)
    end if

    if (kernel_active(K_URBANFLUX)) then
       call ELMxxComputeUrbanFluxes(elm, ierr)
       call check(ierr, logunit, K_URBANFLUX)
    end if

    if (kernel_active(K_LAKEFLUX)) then
       call ELMxxComputeLakeFluxes(elm, ierr)
       call check(ierr, logunit, K_LAKEFLUX)
    end if

    if (kernel_active(K_LAKETEMP)) then
       call ELMxxComputeLakeTemperature(elm, ierr)
       call check(ierr, logunit, K_LAKETEMP)
    end if

    if (kernel_active(K_SOILTEMP)) then
       call ELMxxComputeSoilTemperature(elm, ierr)
       call check(ierr, logunit, K_SOILTEMP)
    end if

    if (kernel_active(K_SOILFLUX)) then
       call ELMxxComputeSoilFluxes(elm, ierr)
       call check(ierr, logunit, K_SOILFLUX)
    end if

    if (kernel_active(K_SURFRUNOFF)) then
       call ELMxxComputeSurfRunInfil(elm, ierr)
       call check(ierr, logunit, K_SURFRUNOFF)
    end if

    if (kernel_active(K_ROOTWATER)) then
       call ELMxxComputeRootWaterUpdate(elm, ierr)
       call check(ierr, logunit, K_ROOTWATER)
    end if

    if (kernel_active(K_HYDRODRAIN)) then
       call ELMxxComputeHydrologyDrainage(elm, ierr)
       call check(ierr, logunit, K_HYDRODRAIN)
    end if

    if (kernel_active(K_LAKEHYDRO)) then
       call ELMxxComputeLakeHydrology(elm, ierr)
       call check(ierr, logunit, K_LAKEHYDRO)
    end if

  end subroutine elmxx_kernels_run

  !-----------------------------------------------------------------------
  subroutine elmxx_kernels_report(elm, npatch, logunit)
    !
    ! Report the range of each active kernel's outputs.
    !
    ! A kernel that runs without aborting has proved almost nothing -- it
    ! reads zeros happily and writes zeros back. What makes a range useful is
    ! that it can be argued with: canopy interception must be positive when it
    ! rains on a canopy and exactly zero when it does not, h2ocan must not
    ! exceed dewmx*(elai+esai), and fwet+fdry must not leave [0,1]. That is
    ! the same standard Stage 2 held the forcing import to -- the range is the
    ! evidence, not the absence of a crash.
    !
    implicit none
    type(ELMxxType), intent(in) :: elm
    integer, intent(in) :: npatch, logunit
    integer :: ierr
    real(r8), allocatable :: intr(:), grnd(:), can(:), fwet(:), fdry(:)
    character(len=*), parameter :: subname = '(elmxx_kernels_report) '

    if (.not. kernel_active(K_CANHYDRO)) return
    if (npatch <= 0) return

    allocate(intr(npatch), grnd(npatch), can(npatch), fwet(npatch), fdry(npatch))

    call ELMxxGetQflxPrecIntr(elm, intr, npatch, ierr); call check(ierr, logunit, K_CANHYDRO)
    call ELMxxGetQflxPrecGrnd(elm, grnd, npatch, ierr); call check(ierr, logunit, K_CANHYDRO)
    call ELMxxGetH2ocan(elm, can, npatch, ierr);        call check(ierr, logunit, K_CANHYDRO)
    call ELMxxGetFwet(elm, fwet, npatch, ierr);         call check(ierr, logunit, K_CANHYDRO)
    call ELMxxGetFdry(elm, fdry, npatch, ierr);         call check(ierr, logunit, K_CANHYDRO)

    write(logunit,*) subname,'rank ',iam,' canhydro over ',npatch,' patches:'
    write(logunit,*) '    qflx_prec_intr [kg/m2/s] ',minval(intr),' .. ',maxval(intr)
    write(logunit,*) '    qflx_prec_grnd [kg/m2/s] ',minval(grnd),' .. ',maxval(grnd)
    write(logunit,*) '    h2ocan         [kg/m2]   ',minval(can) ,' .. ',maxval(can)
    write(logunit,*) '    fwet           [-]       ',minval(fwet),' .. ',maxval(fwet)
    write(logunit,*) '    fdry           [-]       ',minval(fdry),' .. ',maxval(fdry)

    ! Bounds that must hold whatever the forcing is. Reported, not asserted:
    ! this is a diagnostic stage and an abort here would stop a run that is
    ! still useful to inspect. They become assertions when the kernel set is
    ! complete enough for a run to be graded rather than watched.
    if (minval(intr) < 0.0_r8) &
         write(logunit,*) subname,'SUSPECT: negative canopy interception'
    if (minval(can) < 0.0_r8) &
         write(logunit,*) subname,'SUSPECT: negative canopy water'
    if (minval(fwet) < 0.0_r8 .or. maxval(fwet) > 1.0_r8) &
         write(logunit,*) subname,'SUSPECT: fwet outside [0,1]'
    if (minval(fdry) < 0.0_r8 .or. maxval(fdry) > 1.0_r8) &
         write(logunit,*) subname,'SUSPECT: fdry outside [0,1]'

    call shr_sys_flush(logunit)
    deallocate(intr, grnd, can, fwet, fdry)

  end subroutine elmxx_kernels_report

  !-----------------------------------------------------------------------
  subroutine elmxx_report_cantemp(elm, ncol, logunit)
    !
    ! CanopyTemperature's outputs, graded the same way: by whether the ranges
    ! can be argued with. Ground temperature should sit near the 274 K it was
    ! cold-started to, ground specific humidity must be positive and small,
    ! the latent heat of vaporization is ~2.5e6 J/kg at these temperatures and
    ! ~2.83e6 if it sublimates, and soilbeta is a fraction.
    !
    implicit none
    type(ELMxxType), intent(in) :: elm
    integer, intent(in) :: ncol, logunit
    integer :: ierr
    real(r8), allocatable :: tg(:), qg(:), thv(:), htvp(:), sbeta(:), z0mg(:)
    character(len=*), parameter :: subname = '(elmxx_kernels_report) '

    if (.not. kernel_active(K_CANTEMP)) return
    if (ncol <= 0) return

    allocate(tg(ncol), qg(ncol), thv(ncol), htvp(ncol), sbeta(ncol), z0mg(ncol))
    call ELMxxGetTGrnd(elm, tg, ncol, ierr);       call check(ierr, logunit, K_CANTEMP)
    call ELMxxGetQg(elm, qg, ncol, ierr);          call check(ierr, logunit, K_CANTEMP)
    call ELMxxGetThv(elm, thv, ncol, ierr);        call check(ierr, logunit, K_CANTEMP)
    call ELMxxGetHtvp(elm, htvp, ncol, ierr);      call check(ierr, logunit, K_CANTEMP)
    call ELMxxGetSoilbeta(elm, sbeta, ncol, ierr); call check(ierr, logunit, K_CANTEMP)
    call ELMxxGetZ0mg(elm, z0mg, ncol, ierr);      call check(ierr, logunit, K_CANTEMP)

    write(logunit,*) subname,'rank ',iam,' cantemp over ',ncol,' columns:'
    write(logunit,*) '    t_grnd   [K]     ',minval(tg)   ,' .. ',maxval(tg)
    write(logunit,*) '    qg       [kg/kg] ',minval(qg)   ,' .. ',maxval(qg)
    write(logunit,*) '    thv      [K]     ',minval(thv)  ,' .. ',maxval(thv)
    write(logunit,*) '    htvp     [J/kg]  ',minval(htvp) ,' .. ',maxval(htvp)
    write(logunit,*) '    soilbeta [-]     ',minval(sbeta),' .. ',maxval(sbeta)
    write(logunit,*) '    z0mg     [m]     ',minval(z0mg) ,' .. ',maxval(z0mg)

    if (minval(tg) < 200.0_r8 .or. maxval(tg) > 350.0_r8) &
         write(logunit,*) subname,'SUSPECT: ground temperature outside 200-350 K'
    if (minval(qg) < 0.0_r8) &
         write(logunit,*) subname,'SUSPECT: negative ground specific humidity'
    if (minval(sbeta) < 0.0_r8 .or. maxval(sbeta) > 1.0_r8) &
         write(logunit,*) subname,'SUSPECT: soilbeta outside [0,1]'
    if (minval(htvp) <= 0.0_r8) &
         write(logunit,*) subname,'SUSPECT: non-positive latent heat'

    call shr_sys_flush(logunit)
    deallocate(tg, qg, thv, htvp, sbeta, z0mg)

  end subroutine elmxx_report_cantemp

  !-----------------------------------------------------------------------
  subroutine elmxx_report_fluxes(elm, npatch, logunit)
    !
    ! BareGroundFluxes and CanopyFluxes outputs.
    !
    ! These two partition the patches between them -- BareGroundFluxes takes
    ! frac_veg_nosno == 0, CanopyFluxes takes == 1 -- so a patch appears in
    ! exactly one. That is why they are reported together: the union is the
    ! surface energy balance, and a gap or an overlap shows up here as a
    ! patch with no flux or two.
    !
    implicit none
    type(ELMxxType), intent(in) :: elm
    integer, intent(in) :: npatch, logunit
    integer :: ierr
    real(r8), allocatable :: shg(:), shv(:), evs(:), trv(:), tv(:), btr(:)
    character(len=*), parameter :: subname = '(elmxx_kernels_report) '

    if (.not. (kernel_active(K_BAREGRND) .or. kernel_active(K_CANFLUX))) return
    if (npatch <= 0) return

    allocate(shg(npatch), shv(npatch), evs(npatch), trv(npatch), &
             tv(npatch), btr(npatch))
    call ELMxxGetEflxShGrnd(elm, shg, npatch, ierr);  call check(ierr, logunit, K_BAREGRND)
    call ELMxxGetEflxShVeg(elm, shv, npatch, ierr);   call check(ierr, logunit, K_CANFLUX)
    call ELMxxGetQflxEvapSoi(elm, evs, npatch, ierr); call check(ierr, logunit, K_BAREGRND)
    call ELMxxGetQflxTranVeg(elm, trv, npatch, ierr); call check(ierr, logunit, K_CANFLUX)
    call ELMxxGetTVeg(elm, tv, npatch, ierr);         call check(ierr, logunit, K_CANFLUX)
    call ELMxxGetBtran(elm, btr, npatch, ierr);       call check(ierr, logunit, K_CANFLUX)

    write(logunit,*) subname,'rank ',iam,' surface fluxes over ',npatch,' patches:'
    write(logunit,*) '    eflx_sh_grnd  [W/m2]   ',minval(shg),' .. ',maxval(shg)
    write(logunit,*) '    eflx_sh_veg   [W/m2]   ',minval(shv),' .. ',maxval(shv)
    write(logunit,*) '    qflx_evap_soi [kg/m2/s]',minval(evs),' .. ',maxval(evs)
    write(logunit,*) '    qflx_tran_veg [kg/m2/s]',minval(trv),' .. ',maxval(trv)
    write(logunit,*) '    t_veg         [K]      ',minval(tv) ,' .. ',maxval(tv)
    write(logunit,*) '    btran         [-]      ',minval(btr),' .. ',maxval(btr)

    ! Bounds that hold whatever the forcing is.
    if (minval(tv) < 200.0_r8 .or. maxval(tv) > 350.0_r8) &
         write(logunit,*) subname,'SUSPECT: vegetation temperature outside 200-350 K'
    if (minval(btr) < 0.0_r8 .or. maxval(btr) > 1.0_r8) &
         write(logunit,*) subname,'SUSPECT: btran outside [0,1]'
    if (minval(trv) < 0.0_r8) &
         write(logunit,*) subname,'SUSPECT: negative transpiration'

    call shr_sys_flush(logunit)
    deallocate(shg, shv, evs, trv, tv, btr)

  end subroutine elmxx_report_fluxes


  !-----------------------------------------------------------------------
  subroutine elmxx_report_surfrad(elm, npatch, logunit)
    !
    ! SurfaceRadiation's outputs, graded by the one thing that makes a
    ! radiation kernel checkable: CONSERVATION. Incident shortwave is either
    ! absorbed or reflected, so fsa + fsr must equal what came in, and with the
    ! cold-start albedos frozen at 0.2 the split is not merely conserved but
    ! PREDICTABLE -- fsr/(fsa+fsr) should sit at 0.2 to rounding.
    !
    ! That is a much stronger check than a range, and it is available precisely
    ! because the albedos are constants right now. It stops being available the
    ! moment SurfaceAlbedo lands and starts varying them, so it is worth
    ! spending here: it pins the band mapping, the 2-D layout and the
    ! patch-level broadcast all at once. A transposed forc_solad would still
    ! conserve, but it would not hold 0.2 unless both bands carried equal
    ! flux -- and vis/NIR direct differ by day.
    !
    ! sabg is also split soil/snow. With no snow at a cold start, sabg_soil
    ! should account for all of sabg.
    !
    implicit none
    type(ELMxxType), intent(in) :: elm
    integer, intent(in) :: npatch, logunit
    integer :: ierr
    real(r8) :: inc, refl_frac
    real(r8), allocatable :: fsa(:), fsr(:), sabv(:), sabg(:), sabgs(:)
    character(len=*), parameter :: subname = '(elmxx_kernels_report) '

    if (.not. kernel_active(K_SURFRAD)) return
    if (npatch <= 0) return

    allocate(fsa(npatch), fsr(npatch), sabv(npatch), sabg(npatch), sabgs(npatch))
    call ELMxxGetFsa(elm, fsa, npatch, ierr);        call check(ierr, logunit, K_SURFRAD)
    call ELMxxGetFsr(elm, fsr, npatch, ierr);        call check(ierr, logunit, K_SURFRAD)
    call ELMxxGetSabv(elm, sabv, npatch, ierr);      call check(ierr, logunit, K_SURFRAD)
    call ELMxxGetSabg(elm, sabg, npatch, ierr);      call check(ierr, logunit, K_SURFRAD)
    call ELMxxGetSabgSoil(elm, sabgs, npatch, ierr); call check(ierr, logunit, K_SURFRAD)

    write(logunit,*) subname,'rank ',iam,' surfrad over ',npatch,' patches:'
    write(logunit,*) '    fsa       [W/m2] ',minval(fsa)  ,' .. ',maxval(fsa)
    write(logunit,*) '    fsr       [W/m2] ',minval(fsr)  ,' .. ',maxval(fsr)
    write(logunit,*) '    sabv      [W/m2] ',minval(sabv) ,' .. ',maxval(sabv)
    write(logunit,*) '    sabg      [W/m2] ',minval(sabg) ,' .. ',maxval(sabg)
    write(logunit,*) '    sabg_soil [W/m2] ',minval(sabgs),' .. ',maxval(sabgs)

    ! The reflected fraction, which the frozen albedos make a known constant.
    inc = maxval(fsa) + maxval(fsr)
    if (inc > 1.0_r8) then
       refl_frac = maxval(fsr) / inc
       write(logunit,*) '    incident  [W/m2] ',inc
       write(logunit,*) '    fsr/incident [-] ',refl_frac,'  (expect 0.2 while albedos are frozen)'
       if (abs(refl_frac - 0.2_r8) > 1.0e-3_r8) &
            write(logunit,*) subname,'SUSPECT: reflected fraction is not the ', &
                 'cold-start albedo -- check band order or 2-D layout'
    end if

    if (minval(fsa) < 0.0_r8 .or. minval(fsr) < 0.0_r8) &
         write(logunit,*) subname,'SUSPECT: negative absorbed or reflected shortwave'
    if (minval(sabg) < 0.0_r8) &
         write(logunit,*) subname,'SUSPECT: negative ground absorption'
    if (maxval(abs(sabg - sabgs)) > 1.0e-10_r8) &
         write(logunit,*) subname,'SUSPECT: sabg /= sabg_soil with no snow present'

    call shr_sys_flush(logunit)
    deallocate(fsa, fsr, sabv, sabg, sabgs)

  end subroutine elmxx_report_surfrad
  !-----------------------------------------------------------------------
  subroutine check(ierr, logunit, k)
    implicit none
    integer, intent(in) :: ierr, logunit, k
    if (ierr /= ELMXX_SUCCESS) then
       write(logunit,*) '(elmxx_kernels_run) ERROR: kernel ', &
            trim(kernel_name(k)),' returned status ',ierr
       call shr_sys_flush(logunit)
       call shr_sys_abort('(elmxx_kernels_run) ERROR: kernel '// &
            trim(kernel_name(k))//' failed')
    end if
  end subroutine check

end module elmxxKernelMod
