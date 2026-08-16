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
                               ELMxxComputeSoilTemperatureNatural, &
                               ELMxxComputeSoilFluxesNatural, &
                               ELMxxComputeSurfRunInfilHydroActive, &
                               ELMxxComputeRootWaterUpdateNatural, &
                               ELMxxComputeHydrologyDrainageNatural, &
                               ELMxxComputeLakeHydrology, &
                               ELMxxGetQflxPrecIntr, ELMxxGetQflxPrecGrnd, &
                               ELMxxGetH2ocan, ELMxxGetFwet, ELMxxGetFdry, &
                               ELMxxGetTGrnd, ELMxxGetQg, ELMxxGetThv, &
                               ELMxxGetHtvp, ELMxxGetSoilbeta, ELMxxGetZ0mg, &
                               ELMxxGetEflxShGrnd, ELMxxGetEflxShVeg, &
                               ELMxxGetQflxEvapSoi, ELMxxGetQflxTranVeg, &
                               ELMxxGetTVeg, ELMxxGetBtran, &
                               ELMxxGetFsa, ELMxxGetFsr, ELMxxGetSabv, &
                               ELMxxGetUstar, ELMxxGetRb1, ELMxxGetRam1, &
                               ELMxxGetDispla, ELMxxGetZ0mv, ELMxxGetFsun, &
                               ELMxxGetRssun, ELMxxGetRssha, ELMxxGetLaisun, &
                               ELMxxGetSabg, ELMxxGetSabgSoil

  implicit none
  save
  private

  ! Driver order. This is the order the kernels run in, and the order they
  ! should be activated in. It is ELM's driver order, not alphabetical.
  integer, parameter, public :: NKERNEL = 16

  !--------------------------------------------------------------------------
  ! ELM'S DRIVER ORDER, TAKEN FROM elm_driver.F90 RATHER THAN REASONED OUT.
  !
  ! The sequence of processes here must match ELM's, because each kernel reads
  ! what the ones before it wrote and ELM's ordering encodes those
  ! dependencies. Line numbers are elm_driver.F90:
  !
  !    743  CanopyHydrology
  !    780  CanopySunShadeFractions        <- BEFORE SurfaceRadiation
  !    795  SurfaceRadiation
  !    815  UrbanRadiation
  !    835  CanopyTemperature
  !    852  BareGroundFluxes
  !    863  CanopyFluxes
  !    880  UrbanFluxes
  !    899  LakeFluxes
  !    951  LakeTemperature
  !    966  SoilTemperature
  !    983  SoilFluxes
  !   1008  HydrologyNoDrainage, which contains
  !           SurfaceRunoff + Infiltration  -> surfrunoff
  !           SoilWater's root extraction   -> rootwater
  !   1045  LakeHydrology                  <- BEFORE HydrologyDrainage
  !   1345  HydrologyDrainage
  !   1500  SurfaceAlbedo                  <- END of the step; see STATUS
  !
  ! CORRECTED 2026-08-15. This list previously began with surfrad, put
  ! cansunshade after it, and put hydrodrain before lakehydro. All three were
  ! wrong: ELM runs CanopyHydrology first, CanopySunShadeFractions before
  ! SurfaceRadiation, and LakeHydrology before HydrologyDrainage.
  !--------------------------------------------------------------------------

  integer, parameter, public :: K_CANHYDRO    =  1
  integer, parameter, public :: K_CANSUNSHADE =  2
  integer, parameter, public :: K_SURFRAD     =  3
  integer, parameter, public :: K_URBANRAD    =  4
  integer, parameter, public :: K_CANTEMP     =  5
  integer, parameter, public :: K_BAREGRND    =  6
  integer, parameter, public :: K_CANFLUX     =  7
  integer, parameter, public :: K_URBANFLUX   =  8
  integer, parameter, public :: K_LAKEFLUX    =  9
  integer, parameter, public :: K_LAKETEMP    = 10
  integer, parameter, public :: K_SOILTEMP    = 11
  integer, parameter, public :: K_SOILFLUX    = 12
  integer, parameter, public :: K_SURFRUNOFF  = 13
  integer, parameter, public :: K_ROOTWATER   = 14
  integer, parameter, public :: K_LAKEHYDRO   = 15
  integer, parameter, public :: K_HYDRODRAIN  = 16

  character(len=16), parameter, public :: kernel_name(NKERNEL) = [ &
       'canhydro        ', 'cansunshade     ', 'surfrad         ', &
       'urbanrad        ', 'cantemp         ', 'baregrnd        ', &
       'canflux         ', 'urbanflux       ', 'lakeflux        ', &
       'laketemp        ', 'soiltemp        ', 'soilflux        ', &
       'surfrunoff      ', 'rootwater       ', 'lakehydro       ', &
       'hydrodrain      ' ]

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
       ! Runnable, and the earlier reason recorded here was wrong.
       !
       ! It said this kernel "needs SurfaceAlbedo for the per-canopy-layer
       ! sunlit/shaded decomposition". But CanopySunShadeFractionsPatch
       ! COMPUTES that decomposition -- it is the kernel's whole job:
       !     laisun = tlai_z * fsun_z
       !     laisha = tlai_z * (1 - fsun_z)
       ! What it reads is nrad, tlai_z and fsun_z, and the albedo set feeds
       ! only parsun_z/parsha_z, the absorbed PAR that goes on to
       ! photosynthesis. Blocking the whole kernel for the PAR half also
       ! withheld the LAI half, which nothing required.
       !
       ! nlevcan = 1 in both models, and on that path ELM does not need the
       ! two-stream for the first two either -- SurfaceAlbedoMod sets
       ! nrad = 1 and tlai_z = elai outright. Only fsun_z is genuinely a
       ! two-stream output, and ELM's cold start leaves it at zero.
       !
       ! LIMIT: fsun_z frozen at zero means the canopy is entirely shaded, so
       ! parsun/laisun stay zero and sunlit photosynthesis cannot switch on.
       ! ELM's own step one looks the same; a multi-day run does not, and that
       ! resolves when SurfaceAlbedo lands.
       why = ' '

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
       ! Runnable. Two earlier explanations for its NaN were both wrong, and
       ! both are recorded because each cost a cycle of guessing.
       !
       ! WRONG #1: "btran is zero, so canflux cannot be graded." btran only
       ! gates whether transpiration is drawn from rppdry; rppdry was already
       ! NaN by then.
       !
       ! WRONG #2: "it needs laisun/laisha from cansunshade and rssun/rssha
       ! from a Photosynthesis port." Running cansunshade changed nothing, and
       ! zero stomatal resistance is finite arithmetic, not a NaN.
       !
       ! ACTUAL CAUSE: dleaf, the leaf characteristic dimension, was never
       ! seeded. CanopyFluxes forms
       !     cf = 0.01 / (sqrt(uaf) * sqrt(dleaf))    rb = 1/(cf*uaf)
       ! so dleaf = 0 makes cf infinite and rb EXACTLY ZERO -- and rb is the
       ! denominator of
       !     rppdry = fdry*rb*(laisun/(rb+rssun) + laisha/(rb+rssha))/elai
       ! giving elai/0 and 0/0 together. It is a pftcon parameter that ELMxx
       ! stores per patch rather than per PFT, which is why it was missed
       ! alongside z0mr and displar. elmxxPftconMod now reads it and aborts on
       ! a non-positive value for a vegetated PFT.
       !
       ! THE PATTERN, THIRD OCCURRENCE: an unseeded parameter view reads zero,
       ! and zero is a legal number that produces NaN several kernels later.
       ! z0mr/displar zeroed ustar; dleaf zeroed rb. Ask of every parameter a
       ! kernel divides by whether anything actually sets it.
       !
       ! REMAINING LIMIT, real but not a blocker: rssun_iter/rssha_iter are
       ! zero because ELMxx has no Photosynthesis kernel, so stomata offer no
       ! resistance. With btran zero as well, transpiration is off entirely.
       ! Bounded and finite; not yet the right physics.
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
       ! Runnable, on the ...Natural kernel variants.
       !
       ! TWO CORRECTIONS ARE FOLDED IN HERE, both recorded because each cost
       ! real time.
       !
       ! (1) An earlier note claimed these "cannot usefully precede
       ! SurfaceAlbedo". Withdrawn. They do need sabg, but sabg comes from the
       ! SURFRAD KERNEL, which runs earlier in this same timestep and is now
       ! active. SurfaceAlbedo runs at the END of ELM's step, so it was never
       ! a within-step dependency.
       !
       ! (2) An earlier note counted "163 calls, 127 distinct" of
       ! ST_/SF_/SRI_/RWU_/HD_ seeding on top of ELMxxInitSharedMetadata. That
       ! describes the STANDALONE entry points, which run on a separate
       ! allocation and exist so the validation driver can replay ELM
       ! diagnostic snapshots. The ...Natural variants build their view
       ! straight from naturalCol/naturalPatch -- the state the canopy kernels
       ! already use -- so none of those setters is needed. Decision #8 chose
       ! the Natural variants; what was missed is how much that deletes.
       !
       ! What IS needed is in elmxxSoilKernelMod: the shared filters, the
       ! naturalCol views the canopy kernels never touched (the _p1 and _soi
       ! layer representations, thermal properties, column metadata), and the
       ! ground surface energy balance, which ELMxx genuinely does not compute.
       why = ' '

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
  subroutine elmxx_kernels_run(elm, dtime, logunit, phase)
    !
    ! Dispatch the active kernels in driver order.
    !
    ! Order is the array order, never the order they were named on the
    ! namelist -- a kernel reads what the ones before it wrote, so letting the
    ! namelist reorder them would silently change the physics.
    !
    ! TWO PHASES, because something has to happen between them. The ground
    ! surface energy balance (elmxxSoilKernelMod) is built from the radiation
    ! and turbulent fluxes phase 1 produces, and SoilTemperature in phase 2
    ! consumes it. ELMxx does not compute it, so the driver must, and it must
    ! do so at exactly this point. A single dispatch would have left that
    ! ordering as a comment; two phases make it a signature.
    !
    !   phase 1: surfrad .. laketemp   -- radiation, canopy, surface fluxes
    !   phase 2: soiltemp .. lakehydro -- the soil column and hydrology
    !
    implicit none
    type(ELMxxType), intent(in) :: elm
    real(r8), intent(in) :: dtime
    integer, intent(in) :: logunit
    integer, intent(in) :: phase
    integer :: ierr
    character(len=*), parameter :: subname = '(elmxx_kernels_run) '

    if (.not. any_kernel_active) return
    if (phase /= 1 .and. phase /= 2) then
       call shr_sys_abort(subname//'ERROR: phase must be 1 or 2')
    end if

    if (phase == 1) then

    if (kernel_active(K_CANHYDRO)) then
       call ELMxxComputeCanopyHydrology(elm, dtime, ierr)
       call check(ierr, logunit, K_CANHYDRO)
    end if

    if (kernel_active(K_CANSUNSHADE)) then
       call ELMxxComputeCanopySunShadeFractions(elm, ierr)
       call check(ierr, logunit, K_CANSUNSHADE)
    end if

    if (kernel_active(K_SURFRAD)) then
       call ELMxxComputeSurfaceRadiation(elm, ierr)
       call check(ierr, logunit, K_SURFRAD)
    end if

    if (kernel_active(K_URBANRAD)) then
       call ELMxxComputeUrbanRadiation(elm, ierr)
       call check(ierr, logunit, K_URBANRAD)
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

    end if

    if (phase == 2) then

    if (kernel_active(K_SOILTEMP)) then
       call ELMxxComputeSoilTemperatureNatural(elm, ierr)
       call check(ierr, logunit, K_SOILTEMP)
    end if

    if (kernel_active(K_SOILFLUX)) then
       call ELMxxComputeSoilFluxesNatural(elm, ierr)
       call check(ierr, logunit, K_SOILFLUX)
    end if

    if (kernel_active(K_SURFRUNOFF)) then
       call ELMxxComputeSurfRunInfilHydroActive(elm, ierr)
       call check(ierr, logunit, K_SURFRUNOFF)
    end if

    if (kernel_active(K_ROOTWATER)) then
       call ELMxxComputeRootWaterUpdateNatural(elm, ierr)
       call check(ierr, logunit, K_ROOTWATER)
    end if

    if (kernel_active(K_LAKEHYDRO)) then
       call ELMxxComputeLakeHydrology(elm, ierr)
       call check(ierr, logunit, K_LAKEHYDRO)
    end if

    if (kernel_active(K_HYDRODRAIN)) then
       call ELMxxComputeHydrologyDrainageNatural(elm, ierr)
       call check(ierr, logunit, K_HYDRODRAIN)
    end if

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
    write(logunit,*) '    eflx_sh_grnd  [W/m2]   ',minval(shg),' .. ',maxval(shg),' nan ',nonfinite(shg)
    write(logunit,*) '    eflx_sh_veg   [W/m2]   ',minval(shv),' .. ',maxval(shv),' nan ',nonfinite(shv)
    write(logunit,*) '    qflx_evap_soi [kg/m2/s]',minval(evs),' .. ',maxval(evs),' nan ',nonfinite(evs)
    write(logunit,*) '    qflx_tran_veg [kg/m2/s]',minval(trv),' .. ',maxval(trv),' nan ',nonfinite(trv)
    write(logunit,*) '    t_veg         [K]      ',minval(tv) ,' .. ',maxval(tv) ,' nan ',nonfinite(tv)
    write(logunit,*) '    btran         [-]      ',minval(btr),' .. ',maxval(btr),' nan ',nonfinite(btr)

    if (nonfinite(shg) + nonfinite(shv) + nonfinite(evs) + nonfinite(trv) + &
        nonfinite(tv) + nonfinite(btr) > 0) then
       write(logunit,*) subname,'SUSPECT: non-finite surface flux'
       call diagnose_canflux(elm, npatch, logunit)
    end if

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
  subroutine diagnose_canflux(elm, npatch, logunit)
    !
    ! Walk CanopyFluxes' resistance chain when its outputs go non-finite.
    !
    ! The chain is ustar -> ram1 -> uaf -> rb -> the conductances, and every
    ! link divides by the one before it, so a zero anywhere upstream becomes a
    ! NaN downstream. Printing the whole chain says WHICH link failed instead
    ! of leaving it to be guessed -- which is the point, given that guessing
    ! has already lost twice here.
    !
    implicit none
    type(ELMxxType), intent(in) :: elm
    integer, intent(in) :: npatch, logunit
    integer :: ierr
    real(r8), allocatable :: us(:), rb(:), ram(:), dsp(:), z0(:), rsun(:), rsha(:), lsun(:)
    character(len=*), parameter :: subname = '(diagnose_canflux) '

    allocate(us(npatch), rb(npatch), ram(npatch), dsp(npatch), z0(npatch), &
             rsun(npatch), rsha(npatch), lsun(npatch))
    call ELMxxGetUstar(elm, us, npatch, ierr)
    call ELMxxGetRb1(elm, rb, npatch, ierr)
    call ELMxxGetRam1(elm, ram, npatch, ierr)
    call ELMxxGetDispla(elm, dsp, npatch, ierr)
    call ELMxxGetZ0mv(elm, z0, npatch, ierr)
    call ELMxxGetRssun(elm, rsun, npatch, ierr)
    call ELMxxGetRssha(elm, rsha, npatch, ierr)
    call ELMxxGetLaisun(elm, lsun, npatch, ierr)

    write(logunit,*) subname,'resistance chain (range, then nan count):'
    write(logunit,*) '    ustar  ',minval(us)  ,maxval(us)  ,nonfinite(us)
    write(logunit,*) '    ram1   ',minval(ram) ,maxval(ram) ,nonfinite(ram)
    write(logunit,*) '    rb1    ',minval(rb)  ,maxval(rb)  ,nonfinite(rb)
    write(logunit,*) '    displa ',minval(dsp) ,maxval(dsp) ,nonfinite(dsp)
    write(logunit,*) '    z0mv   ',minval(z0)  ,maxval(z0)  ,nonfinite(z0)
    write(logunit,*) '    rssun  ',minval(rsun),maxval(rsun),nonfinite(rsun)
    write(logunit,*) '    rssha  ',minval(rsha),maxval(rsha),nonfinite(rsha)
    write(logunit,*) '    laisun ',minval(lsun),maxval(lsun),nonfinite(lsun)
    call shr_sys_flush(logunit)
    deallocate(us, rb, ram, dsp, z0, rsun, rsha, lsun)

  end subroutine diagnose_canflux

  !-----------------------------------------------------------------------
  integer function nonfinite(a)
    !
    ! Count NaN/Inf. This exists because MINVAL AND MAXVAL DO NOT SHOW THEM:
    ! gfortran's reductions skip NaN, so an array with three NaN among
    ! seventeen values reports a clean finite range. That blind spot hid a
    ! NaN in CanopyFluxes' component fluxes through several "graded" runs --
    ! the same failure mode as a dict-keyed comparison that cannot see
    ! ordering (STATUS G). Every range printed here is now accompanied by
    ! this count.
    !
    implicit none
    real(r8), intent(in) :: a(:)
    integer :: i
    nonfinite = 0
    do i = 1, size(a)
       if (.not. (abs(a(i)) >= 0.0_r8)) nonfinite = nonfinite + 1
    end do
  end function nonfinite

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
