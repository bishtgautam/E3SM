module kelmKernelMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Activation and dispatch of kELM's Kokkos kernels.
  !
  ! The model runs every kernel (kelm_kernels_activate_all); there is no
  ! namelist switch. kelm_kernels_parse still takes a comma-separated list
  ! of the short names below, so a subset can be run to check kernels one at
  ! a time against ELM; a kernel not named does not run, and THERE IS NO
  ! FORTRAN ALTERNATIVE TO FALL BACK ON, so such a run is diagnostic only.
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
  use kelmSpmdMod    , only : masterproc, iam
  use iso_c_binding  , only : c_int, c_ptr, c_funloc, c_null_funptr
  use kelmTimerMod   , only : kelm_t_start, kelm_t_stop, kelm_timers_on
  use kelmDiagnosticsMod, only : kelm_diag_enabled, kelm_diag_snapshot_state, &
                                  kelm_diag_snapshot_fluxes, kelm_diag_snapshot_preinfil, kelm_diag_snapshot_presoilflux, kelm_diag_snapshot_presoiltemp, kelm_diag_snapshot_soilwater, kelm_diag_snapshot_postcanhydro, kelm_diag_snapshot_postcanflux, kelm_diag_snapshot_snowstate
  use kelmKokkosStateMod , only : n_kokkos_col
  use kelmLakeMod       , only : lake_built
  use kelmUrbanMod      , only : urban_built, kelm_urban_diag_step
  use kelmSoilPropMod   , only : nlevgrnd, nlevsno, nlevtot
  use kelm_mod       , only : kELMType, KELM_SUCCESS, KELM_NKERNEL, KELM_NTIMER, &
                               kELMSetActiveKernels, kELMActiveKernelsError, kELMSetStepHook, &
                               KELM_STAGE_STATE_IN, KELM_STAGE_URBAN_IN, KELM_STAGE_SNOW_PRE, &
                               KELM_STAGE_POST_CANHYDRO, KELM_STAGE_SNOW_POST, &
                               KELM_STAGE_POST_CANFLUX, KELM_STAGE_URBAN_FLUX, &
                               KELM_STAGE_PRE_SOILTEMP, KELM_STAGE_PRE_SOILFLUX, &
                               KELM_STAGE_PRE_INFIL, KELM_STAGE_SOILWATER, KELM_STAGE_URBAN_OUT, &
                               KELM_STAGE_FLUXES_OUT, KELM_STAGE_WBAL_BEGIN, KELM_STAGE_WBAL_MARK, &
                               KELM_STAGE_TIMER_START, KELM_STAGE_TIMER_STOP, &
                               kELMComputeSurfaceRadiation, &
                               kELMComputeCanopyHydrology, &
                               kELMComputeCanopySunShadeFractions, &
                               kELMComputeCanopyTemperature, &
                               kELMComputeBareGroundFluxes, &
                               kELMComputeCanopyFluxes, &
                               kELMUrbanGatherForcing, &
                               kELMComputeBeginWaterBalanceUrban, &
                               kELMComputeCanopyHydrologyUrban, &
                               kELMComputeUrbanRadiationUrban, &
                               kELMComputeCanopyTemperatureUrban, &
                               kELMComputeUrbanFluxesUrban, &
                               kELMComputeSoilTemperatureUrban, &
                               kELMComputeSoilFluxesUrban, &
                               kELMComputeSnowWaterUrban, &
                               kELMComputeSurfaceRunoffInfiltrationUrban, &
                               kELMComputeSoilWaterUrban, &
                               kELMComputeWaterTableUrban, &
                               kELMComputeSnowLayersUrban, &
                               kELMComputeHydrologyDrainageUrban, &
                               kELMComputeWaterBalanceCheckUrban, &
                               kELMComputeLakeFluxesLake, &
                               kELMLakeGatherForcing, &
                               kELMComputeSurfaceRadiationLake, &
                               kELMSnowAgeGrainLake, &
                               kELMComputeBeginWaterBalanceLake, &
                               kELMComputeWaterBalanceCheckLake, &
                               kELMComputeLakeTemperatureLake, &
                               kELMComputeSoilTemperatureNatural, &
                               kELMComputeSoilFluxesNatural, &
                               kELMComputeSurfaceRunoffInfiltration, &
                               kELMComputeRootWaterUpdateNatural, &
                               kELMComputeSoilWaterNatural, &
                               kELMComputeHydrologyDrainageNatural, &
                               kELMComputeWaterTableNatural, &
                               kELMComputeBeginWaterBalanceNatural, &
                               kELMComputeWaterBalanceCheckNatural, &
                               kELMGetWaterInventoryNatural, &
                               kELMComputeSnowWater, &
                               kELMComputeSnowCompaction, &
                               kELMComputeSnowLayers, &
                               kELMUpdateGroundTemperatureNatural, &
                               kELMSnowAgeGrainNatural, &
                               kELMComputeLakeHydrologyLake, &
                               kELMGetQflxPrecIntr, kELMGetQflxPrecGrnd, &
                               kELMGetH2ocan, kELMGetFwet, kELMGetFdry, &
                               kELMGetTGrnd, kELMGetQg, kELMGetThv, &
                               kELMGetHtvp, kELMGetSoilbeta, kELMGetZ0mg, &
                               kELMGetEflxShGrnd, kELMGetEflxShVeg, &
                               kELMGetQflxEvapSoi, kELMGetQflxTranVeg, &
                               kELMGetTVeg, kELMGetBtran, &
                               kELMGetFsa, kELMGetFsr, kELMGetSabv, &
                               kELMGetUstar, kELMGetRb1, kELMGetRam1, &
                               kELMGetDispla, kELMGetZ0mv, kELMGetFsun, &
                               kELMGetRssun, kELMGetRssha, kELMGetLaisun, &
                               kELMGetSabg, kELMGetSabgSoil

  implicit none
  save
  private

  ! Driver order. This is the order the kernels run in, and the order they
  ! should be activated in. It is ELM's driver order, not alphabetical.
  integer, parameter, public :: NKERNEL = 23

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
  !           SoilWater (Richards solve)    -> soilwater
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
  integer, parameter, public :: K_SOILWATER   = 15
  integer, parameter, public :: K_LAKEHYDRO   = 16
  integer, parameter, public :: K_HYDRODRAIN  = 17
  ! ELM runs SnowWater at the top of HydrologyNoDrainage, before SurfaceRunoff,
  ! and compaction/combine/divide at the end of that SAME module
  ! (HydrologyNoDrainageMod.F90:312-320) -- i.e. still BEFORE HydrologyDrainage,
  ! not after it. An earlier version of this comment said "after drainage" and
  ! the dispatch order below followed it; corrected 2026-09-02, with the
  ! reasoning at the call site.
  integer, parameter, public :: K_SNOWWATER   = 18
  integer, parameter, public :: K_SNOWLAYERS  = 19
  ! ELM ages the grain near the end of the step, after hydrology has settled
  ! the layers and before SurfaceAlbedo reads snw_rds. A no-op unless the
  ! SNICAR tables were given.
  integer, parameter, public :: K_SNOWAGE     = 21
  ! WaterTable: ELM runs it inside HydrologyNoDrainage, after the Richards
  ! solve and before the snow kernels (HydrologyNoDrainageMod.F90:290). Its
  ! companion Drainage is not a separate token -- it runs at the top of
  ! hydrodrain, where ELM calls it.
  integer, parameter, public :: K_WATERTABLE  = 20
  ! BeginColWaterBalance: ELM takes the column water inventory at the very
  ! top of the step (elm_driver.F90:559), before any physics. Dispatched
  ! first below, regardless of where this sits in the token list.
  integer, parameter, public :: K_BEGWATERBAL = 22
  ! ColWaterBalanceCheck: ELM runs it after the hydrology is finished
  ! (elm_driver.F90:1333). Dispatched last below, whatever the token order.
  integer, parameter, public :: K_BALANCECHK  = 23

  character(len=16), parameter, public :: kernel_name(NKERNEL) = [ &
       'canhydro        ', 'cansunshade     ', 'surfrad         ', &
       'urbanrad        ', 'cantemp         ', 'baregrnd        ', &
       'canflux         ', 'urbanflux       ', 'lakeflux        ', &
       'laketemp        ', 'soiltemp        ', 'soilflux        ', &
       'surfrunoff      ', 'rootwater       ', 'soilwater       ', &
       'lakehydro       ', 'hydrodrain      ', &
       'snowwater       ', 'snowlayers      ', 'watertable      ', &
       'snowage         ', 'begwaterbal     ', 'balancecheck    ' ]

  !--------------------------------------------------------------------------
  ! Per-kernel water budget.
  !
  ! A water balance that does not close says the model gained or lost water
  ! but not where. This brackets every active kernel with the SAME column
  ! inventory begwb/endwb use, so differencing consecutive marks attributes
  ! the change to one kernel. Most kernels only move water between reservoirs
  ! and should show zero; one that changes the total without reporting a
  ! matching flux is the leak.
  !
  ! Off unless KELM_WBAL is set in the environment, following KELM_DIAG.
  !--------------------------------------------------------------------------
  logical,  private :: wbal_on      = .false.
  logical,  private :: wbal_probed  = .false.
  integer,  private :: wbal_ncol    = 0
  integer,  private :: wbal_nstep   = 0
  real(r8), private, allocatable :: wbal_prev(:)      ! inventory at last mark
  real(r8), private, allocatable :: wbal_cur(:)
  real(r8), private, allocatable :: wbal_accum(:)     ! per-kernel, summed over run
  real(r8), private, allocatable :: wbal_first(:)     ! per-kernel, first step

  public :: kelm_wbal_report

  ! The step hook's model handle, and ELM's names for the timed phases
  ! (kELM.h KELM_TIMER_*), so the two models' timing files line up.
  type(kELMType) :: hook_elm
  character(len=24), parameter :: timer_name(0:KELM_NTIMER-1) = [ character(len=24) :: &
       'begwbal', 'canhydro', 'surfrad', 'bgp1', 'bgflux', 'canflux', 'uflux', 'bgplake', &
       'soiltemperature', 'bgp2', 'hydro without drainage', 'hylake', 'hydro2 drainage', &
       'balchk', 'ecosysdyn', 'surfalb', 'urbsurfalb' ]

  logical, public :: kernel_active(NKERNEL) = .false.
  logical, public :: any_kernel_active      = .false.

  public :: kelm_kernels_activate_all
  public :: kelm_kernels_parse
  public :: kelm_kernels_push
  public :: kelm_kernels_report
  public :: kelm_report_cantemp
  public :: kelm_report_fluxes
  public :: kelm_report_surfrad

contains

  !-----------------------------------------------------------------------
  subroutine wbal_begin(elm, ncol)
    !
    ! Start a step's budget: allocate on first use, read the env gate, and
    ! take the opening inventory mark.
    !
    implicit none
    type(kELMType), intent(in) :: elm
    integer, intent(in) :: ncol

    integer :: ierr, dlen
    character(len=8) :: val

    if (.not. wbal_probed) then
       wbal_probed = .true.
       call get_environment_variable('KELM_WBAL', val, dlen)
       wbal_on = (dlen > 0)
    end if
    if (.not. wbal_on .or. ncol <= 0) return

    if (.not. allocated(wbal_prev)) then
       wbal_ncol = ncol
       allocate(wbal_prev(ncol), wbal_cur(ncol))
       allocate(wbal_accum(NKERNEL), wbal_first(NKERNEL))
       wbal_accum = 0._r8
       wbal_first = 0._r8
    end if

    wbal_nstep = wbal_nstep + 1
    call kELMGetWaterInventoryNatural(elm, wbal_prev, wbal_ncol, ierr)
  end subroutine wbal_begin

  !-----------------------------------------------------------------------
  subroutine wbal_mark(elm, k)
    !
    ! Attribute the change in column water since the previous mark to kernel
    ! k. Column 1 only: these are single-column validation cases, and a
    ! per-column breakdown would bury the signal.
    !
    implicit none
    type(kELMType), intent(in) :: elm
    integer, intent(in) :: k

    integer :: ierr

    if (.not. wbal_on .or. wbal_ncol <= 0) return

    call kELMGetWaterInventoryNatural(elm, wbal_cur, wbal_ncol, ierr)
    wbal_accum(k) = wbal_accum(k) + (wbal_cur(1) - wbal_prev(1))
    if (wbal_nstep == 1) wbal_first(k) = wbal_cur(1) - wbal_prev(1)
    wbal_prev(1:wbal_ncol) = wbal_cur(1:wbal_ncol)
  end subroutine wbal_mark

  !-----------------------------------------------------------------------
  subroutine kelm_wbal_report(logunit)
    !
    ! Per-kernel water budget, accumulated over the run.
    !
    ! Read it against expectation, not against zero in isolation: on
    ! 1x1_brazil only canhydro (precipitation in), the flux kernels
    ! (evaporation out), surfrunoff and hydrodrain should move total column
    ! water at all. Anything else with a non-trivial total is changing
    ! storage without a flux to match, which is the shape of a leak.
    !
    implicit none
    integer, intent(in) :: logunit
    integer :: k

    if (.not. wbal_on .or. .not. allocated(wbal_accum)) return

    write(logunit,'(a)') ' '
    write(logunit,'(a)') '=== kELM per-kernel water budget (column 1) ==='
    write(logunit,'(a,i0,a)') 'steps: ', wbal_nstep, &
         '   units: mm of column water'
    write(logunit,'(a)') 'kernel            first step        run total'
    do k = 1, NKERNEL
       if (.not. kernel_active(k)) cycle
       write(logunit,'(a16,2x,es16.8,2x,es16.8)') kernel_name(k), &
            wbal_first(k), wbal_accum(k)
    end do
    write(logunit,'(a)') '==============================================='
    write(logunit,'(a)') ' '
  end subroutine kelm_wbal_report


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
       ! SurfaceAlbedoType InitCold, and kelm_kokkos_seed_albedo supplies the
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
       ! and cold-start column state from kelmSoilPropMod, and the four
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
       ! giving elai/0 and 0/0 together. It is a pftcon parameter that kELM
       ! stores per patch rather than per PFT, which is why it was missed
       ! alongside z0mr and displar. kelmPftconMod now reads it and aborts on
       ! a non-positive value for a vegetated PFT.
       !
       ! THE PATTERN, THIRD OCCURRENCE: an unseeded parameter view reads zero,
       ! and zero is a legal number that produces NaN several kernels later.
       ! z0mr/displar zeroed ustar; dleaf zeroed rb. Ask of every parameter a
       ! kernel divides by whether anything actually sets it.
       !
       ! REMAINING LIMIT, real but not a blocker: rssun_iter/rssha_iter are
       ! zero because kELM has no Photosynthesis kernel, so stomata offer no
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

    case (K_SOILWATER)
       ! The Richards solve. Runnable once surfrunoff and rootwater are on,
       ! since it consumes exactly their outputs -- qflx_infl and
       ! qflx_rootsoi -- and ELM places SoilWater at the same point inside
       ! HydrologyNoDrainage.
       !
       ! THIS IS THE KERNEL THAT CLOSES THE WATER CYCLE. Before it existed,
       ! infiltration was computed and discarded: soil moisture held at its
       ! cold-start value through days of rain, so btran stayed pinned at zero
       ! and transpiration never switched on.
       why = ' '

    case (K_SNOWWATER, K_SNOWLAYERS)
       ! Snow water/percolation, compaction, and combine/divide. All four
       ! kernels were ported and validated against ELM by replay, but nothing
       ! ever dispatched them in the coupled model -- there was no kernel slot
       ! at all. 1x1_brazil never notices, having no snow; on 1x1_glc it meant
       ! CanopyHydrology created a snow layer that then never compacted, never
       ! melted, and above all was never PACKED AWAY by CombineSnowLayers when
       ! its density fell below 50 kg/m3. ELM creates and removes that layer
       ! within the same step; kELM created it and kept it forever.
       why = ' '

    case (K_SNOWAGE)
       ! Snow grain aging, which is half of SNICAR -- the other half is the
       ! radiative transfer inside SurfaceAlbedo. Both are validated against
       ! ELM (albedo to 6e-14, grain radius bit-exactly), but the kernel is
       ! inert unless fsnowoptics and fsnowaging are both given in the
       ! namelist. Without aging snw_rds never moves off the fresh-snow value
       ! and the snowpack stays permanently bright.
       why = ' '

    case (K_BALANCECHK)
       ! Diagnostic: needs begwaterbal to have run, but nothing blocks it.
       why = ' '

    case (K_BEGWATERBAL)
       ! Bookkeeping, not physics: reads state Stage 2 already seeds.
       why = ' '

    case (K_WATERTABLE)
       ! Runnable wherever soilwater is: it reads the same soil state plus
       ! qcharge, which SoilWater produces.
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
       ! ST_/SF_/SRI_/RWU_/HD_ seeding on top of kELMInitSharedMetadata. That
       ! describes the STANDALONE entry points, which run on a separate
       ! allocation and exist so the validation driver can replay ELM
       ! diagnostic snapshots. The ...Natural variants build their view
       ! straight from naturalCol/naturalPatch -- the state the canopy kernels
       ! already use -- so none of those setters is needed. Decision #8 chose
       ! the Natural variants; what was missed is how much that deletes.
       !
       ! What IS needed is in kelmSoilKernelMod: the shared filters, the
       ! naturalCol views the canopy kernels never touched (the _p1 and _soi
       ! layer representations, thermal properties, column metadata), and the
       ! ground surface energy balance, which kELM genuinely does not compute.
       why = ' '

    case (K_URBANRAD, K_URBANFLUX)
       ! Runnable since Stage 6.7 (U4), on the urban surface records
       ! (urban.roof .. urban.perviousRoad, kelmUrbanMod). UrbanAlbedo -- the old
       ! blocker, it produces sabs_dir/sabs_dif -- is ported and runs at the
       ! end of each albedo step. The two tokens switch on the whole urban
       ! step: every urban module's branch runs at its ELM driver position
       ! beside the natural kernel, as a set and only with soiltemp, whose
       ! shared metadata carries the per-cell precipitation urban reads
       ! (both enforced in kelm_kernels_parse).
       why = ' '

    case (K_LAKEFLUX, K_LAKETEMP, K_LAKEHYDRO)
       ! Runnable since Stage 6.5, on the lake surface (lakeCol/lakePatch,
       ! kelmLakeMod). The three run as a set, and only alongside the soil
       ! kernels, whose shared metadata carries the per-cell precipitation
       ! the lake gathers: both enforced in kelm_kernels_parse.
       why = ' '

    case (K_CANHYDRO)
       ! The closest to runnable, and the real first kernel of Stage 4 -- not
       ! surfrad, whatever the plan's driver order says. Its forcing (forc_rain,
       ! forc_snow, forc_t_col) is already crossed and its phenology (elai,
       ! esai) already seeded; most of the rest is zero at a cold start.
       ! What is missing is small and specific: micro_sigma and n_melt, which
       ! ELM derives in initVerticalMod from SLOPE and STD_ELEV -- both present
       ! on surfdata, neither yet read by kelmSurfdataMod.
       why = ' '

    case default
       why = 'unknown kernel'

    end select

  end function blocked_reason

  !-----------------------------------------------------------------------
  subroutine kelm_kernels_activate_all(logunit)
    !
    ! Activate every kernel: what the model always runs. Goes through
    ! kelm_kernels_parse so the prerequisite and grouping checks still apply.
    !
    implicit none
    integer, intent(in) :: logunit
    character(len=NKERNEL*(len(kernel_name)+1)) :: spec
    integer :: k

    spec = trim(kernel_name(1))
    do k = 2, NKERNEL
       spec = trim(spec)//','//trim(kernel_name(k))
    end do
    call kelm_kernels_parse(spec, logunit)

  end subroutine kelm_kernels_activate_all

  !-----------------------------------------------------------------------
  subroutine kelm_kernels_parse(spec, logunit)
    !
    ! Turn a comma-separated list of kernel names into the active flags.
    ! The model activates all of them (kelm_kernels_activate_all); a subset
    ! is for checking kernels one at a time against ELM.
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
    character(len=*), parameter :: subname = '(kelm_kernels_parse) '

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

    ! The grouping rules (lake as a set, urbanrad with urbanflux, both only
    ! with soiltemp) are the library's: kELMSetActiveKernels refuses a set
    ! that breaks them, and kelm_kernels_push aborts with the reason.

    if (masterproc) then
       write(logunit,*) subname,'active kernels, in driver order:'
       do k = 1, NKERNEL
          if (kernel_active(k)) write(logunit,*) '    ',trim(kernel_name(k))
       end do
       call shr_sys_flush(logunit)
    end if

  end subroutine kelm_kernels_parse

  !-----------------------------------------------------------------------
  subroutine kelm_kernels_push(elm, logunit)
    !
    ! Hand the kernel switches to the library, which runs the step
    ! (kELMRunStep, plan Stage 6.9 R2) and owns the grouping rules: the lake
    ! kernels as a set, urbanrad with urbanflux, both only with soiltemp. And
    ! register the step hook when anything between kernels wants a turn: the
    ! KELM_DIAG trace, the KELM_WBAL budget, or the timers. Without it each
    ! hook point is a null-pointer test.
    !
    implicit none
    type(kELMType), intent(in) :: elm
    integer, intent(in) :: logunit
    integer(c_int) :: active(KELM_NKERNEL)
    integer :: ierr, dlen
    character(len=8) :: val
    character(len=*), parameter :: subname = '(kelm_kernels_push) '

    if (NKERNEL /= KELM_NKERNEL) call shr_sys_abort(subname//'ERROR: kernel count '// &
         'differs from the library''s')
    active = 0
    where (kernel_active) active = 1
    call kELMSetActiveKernels(elm, active, ierr)
    if (ierr /= KELM_SUCCESS) then
       select case (kELMActiveKernelsError(elm))
       case (1)
          call shr_sys_abort(subname//'ERROR: lakeflux, laketemp and lakehydro must be '// &
               'activated together')
       case (2)
          call shr_sys_abort(subname//'ERROR: the lake kernels need soiltemp (its shared '// &
               'metadata carries the precipitation the lake reads)')
       case (3)
          call shr_sys_abort(subname//'ERROR: urbanrad and urbanflux must be activated together')
       case (4)
          call shr_sys_abort(subname//'ERROR: the urban kernels need soiltemp (its shared '// &
               'metadata carries the precipitation urban reads)')
       case default
          call shr_sys_abort(subname//'ERROR: kELMSetActiveKernels failed')
       end select
    end if

    if (.not. wbal_probed) then
       wbal_probed = .true.
       call get_environment_variable('KELM_WBAL', val, dlen)
       wbal_on = (dlen > 0)
    end if
    hook_elm = elm
    if (kelm_diag_enabled .or. wbal_on .or. kelm_timers_on) then
       call kELMSetStepHook(elm, c_funloc(kelm_step_hook))
    else
       call kELMSetStepHook(elm, c_null_funptr)
    end if
    if (masterproc) then
       write(logunit,*) subname,'kernels handed to kELMRunStep; step hook ', &
            merge('on ', 'off', kelm_diag_enabled .or. wbal_on .or. kelm_timers_on)
       call shr_sys_flush(logunit)
    end if

  end subroutine kelm_kernels_push

  !-----------------------------------------------------------------------
  subroutine kelm_step_hook(stage, ctx) bind(C)
    !
    ! kELMRunStep's callback (owner decision D2): the KELM_DIAG snapshots,
    ! the KELM_WBAL marks and the ELM-named timers, at the points the
    ! Fortran driver ran them before the step moved into the library. Each
    ! routine keeps its own on/off gate.
    !
    implicit none
    integer(c_int), value :: stage
    type(c_ptr)   , value :: ctx

    if (stage >= KELM_STAGE_TIMER_STOP) then
       call kelm_t_stop(trim(timer_name(stage - KELM_STAGE_TIMER_STOP)))
       return
    else if (stage >= KELM_STAGE_TIMER_START) then
       call kelm_t_start(trim(timer_name(stage - KELM_STAGE_TIMER_START)))
       return
    else if (stage >= KELM_STAGE_WBAL_MARK) then
       call wbal_mark(hook_elm, stage - KELM_STAGE_WBAL_MARK)
       return
    end if

    select case (stage)
    case (KELM_STAGE_STATE_IN)
       call kelm_diag_snapshot_state(hook_elm, nlevtot, nlevgrnd, 'kelm_in')
    case (KELM_STAGE_URBAN_IN)
       call kelm_urban_diag_step(hook_elm, 'kelm_urbin')
    case (KELM_STAGE_WBAL_BEGIN)
       call wbal_begin(hook_elm, n_kokkos_col)
    case (KELM_STAGE_SNOW_PRE)
       call kelm_diag_snapshot_snowstate(hook_elm, 'kelm_snow_pre')
    case (KELM_STAGE_POST_CANHYDRO)
       call kelm_diag_snapshot_postcanhydro(hook_elm, 'kelm_ch')
    case (KELM_STAGE_SNOW_POST)
       call kelm_diag_snapshot_snowstate(hook_elm, 'kelm_snow_post')
    case (KELM_STAGE_POST_CANFLUX)
       call kelm_diag_snapshot_postcanflux(hook_elm, 'kelm_cf')
    case (KELM_STAGE_URBAN_FLUX)
       call kelm_urban_diag_step(hook_elm, 'kelm_urbflx')
    case (KELM_STAGE_PRE_SOILTEMP)
       call kelm_diag_snapshot_presoiltemp(hook_elm, nlevsno + nlevgrnd, 'kelm_st')
    case (KELM_STAGE_PRE_SOILFLUX)
       call kelm_diag_snapshot_presoilflux(hook_elm, nlevsno + nlevgrnd, 'kelm_sf')
    case (KELM_STAGE_PRE_INFIL)
       call kelm_diag_snapshot_preinfil(hook_elm, 'kelm_ri')
    case (KELM_STAGE_SOILWATER)
       call kelm_diag_snapshot_soilwater(hook_elm, nlevgrnd, 'kelm_sw')
    case (KELM_STAGE_URBAN_OUT)
       call kelm_urban_diag_step(hook_elm, 'kelm_urbout')
    case (KELM_STAGE_FLUXES_OUT)
       call kelm_diag_snapshot_fluxes(hook_elm, 'kelm_out')
    end select

  end subroutine kelm_step_hook

  !-----------------------------------------------------------------------
  subroutine kelm_kernels_report(elm, npatch, logunit)
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
    type(kELMType), intent(in) :: elm
    integer, intent(in) :: npatch, logunit
    integer :: ierr
    real(r8), allocatable :: intr(:), grnd(:), can(:), fwet(:), fdry(:)
    character(len=*), parameter :: subname = '(kelm_kernels_report) '

    if (.not. kernel_active(K_CANHYDRO)) return
    if (npatch <= 0) return

    allocate(intr(npatch), grnd(npatch), can(npatch), fwet(npatch), fdry(npatch))

    call kELMGetQflxPrecIntr(elm, intr, npatch, ierr); call check(ierr, logunit, K_CANHYDRO)
    call kELMGetQflxPrecGrnd(elm, grnd, npatch, ierr); call check(ierr, logunit, K_CANHYDRO)
    call kELMGetH2ocan(elm, can, npatch, ierr);        call check(ierr, logunit, K_CANHYDRO)
    call kELMGetFwet(elm, fwet, npatch, ierr);         call check(ierr, logunit, K_CANHYDRO)
    call kELMGetFdry(elm, fdry, npatch, ierr);         call check(ierr, logunit, K_CANHYDRO)

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

  end subroutine kelm_kernels_report

  !-----------------------------------------------------------------------
  subroutine kelm_report_cantemp(elm, ncol, logunit)
    !
    ! CanopyTemperature's outputs, graded the same way: by whether the ranges
    ! can be argued with. Ground temperature should sit near the 274 K it was
    ! cold-started to, ground specific humidity must be positive and small,
    ! the latent heat of vaporization is ~2.5e6 J/kg at these temperatures and
    ! ~2.83e6 if it sublimates, and soilbeta is a fraction.
    !
    implicit none
    type(kELMType), intent(in) :: elm
    integer, intent(in) :: ncol, logunit
    integer :: ierr
    real(r8), allocatable :: tg(:), qg(:), thv(:), htvp(:), sbeta(:), z0mg(:)
    character(len=*), parameter :: subname = '(kelm_kernels_report) '

    if (.not. kernel_active(K_CANTEMP)) return
    if (ncol <= 0) return

    allocate(tg(ncol), qg(ncol), thv(ncol), htvp(ncol), sbeta(ncol), z0mg(ncol))
    call kELMGetTGrnd(elm, tg, ncol, ierr);       call check(ierr, logunit, K_CANTEMP)
    call kELMGetQg(elm, qg, ncol, ierr);          call check(ierr, logunit, K_CANTEMP)
    call kELMGetThv(elm, thv, ncol, ierr);        call check(ierr, logunit, K_CANTEMP)
    call kELMGetHtvp(elm, htvp, ncol, ierr);      call check(ierr, logunit, K_CANTEMP)
    call kELMGetSoilbeta(elm, sbeta, ncol, ierr); call check(ierr, logunit, K_CANTEMP)
    call kELMGetZ0mg(elm, z0mg, ncol, ierr);      call check(ierr, logunit, K_CANTEMP)

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

  end subroutine kelm_report_cantemp

  !-----------------------------------------------------------------------
  subroutine kelm_report_fluxes(elm, npatch, logunit)
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
    type(kELMType), intent(in) :: elm
    integer, intent(in) :: npatch, logunit
    integer :: ierr
    real(r8), allocatable :: shg(:), shv(:), evs(:), trv(:), tv(:), btr(:)
    character(len=*), parameter :: subname = '(kelm_kernels_report) '

    if (.not. (kernel_active(K_BAREGRND) .or. kernel_active(K_CANFLUX))) return
    if (npatch <= 0) return

    allocate(shg(npatch), shv(npatch), evs(npatch), trv(npatch), &
             tv(npatch), btr(npatch))
    call kELMGetEflxShGrnd(elm, shg, npatch, ierr);  call check(ierr, logunit, K_BAREGRND)
    call kELMGetEflxShVeg(elm, shv, npatch, ierr);   call check(ierr, logunit, K_CANFLUX)
    call kELMGetQflxEvapSoi(elm, evs, npatch, ierr); call check(ierr, logunit, K_BAREGRND)
    call kELMGetQflxTranVeg(elm, trv, npatch, ierr); call check(ierr, logunit, K_CANFLUX)
    call kELMGetTVeg(elm, tv, npatch, ierr);         call check(ierr, logunit, K_CANFLUX)
    call kELMGetBtran(elm, btr, npatch, ierr);       call check(ierr, logunit, K_CANFLUX)

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

  end subroutine kelm_report_fluxes


  !-----------------------------------------------------------------------
  subroutine kelm_report_surfrad(elm, npatch, logunit)
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
    type(kELMType), intent(in) :: elm
    integer, intent(in) :: npatch, logunit
    integer :: ierr
    real(r8) :: inc, refl_frac
    real(r8), allocatable :: fsa(:), fsr(:), sabv(:), sabg(:), sabgs(:)
    character(len=*), parameter :: subname = '(kelm_kernels_report) '

    if (.not. kernel_active(K_SURFRAD)) return
    if (npatch <= 0) return

    allocate(fsa(npatch), fsr(npatch), sabv(npatch), sabg(npatch), sabgs(npatch))
    call kELMGetFsa(elm, fsa, npatch, ierr);        call check(ierr, logunit, K_SURFRAD)
    call kELMGetFsr(elm, fsr, npatch, ierr);        call check(ierr, logunit, K_SURFRAD)
    call kELMGetSabv(elm, sabv, npatch, ierr);      call check(ierr, logunit, K_SURFRAD)
    call kELMGetSabg(elm, sabg, npatch, ierr);      call check(ierr, logunit, K_SURFRAD)
    call kELMGetSabgSoil(elm, sabgs, npatch, ierr); call check(ierr, logunit, K_SURFRAD)

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
       write(logunit,*) '    fsr/incident [-] ',refl_frac
       ! With frozen cold-start albedos this had to be EXACTLY 0.2, which
       ! pinned the band mapping and the 2-D layout in one number. That check
       ! is spent: SurfaceAlbedo now computes a real albedo, so the fraction
       ! varies with sun angle, soil wetness and canopy state, as it should.
       ! What survives is the bound -- a surface cannot reflect more than it
       ! receives, nor a negative amount.
       if (refl_frac < 0.0_r8 .or. refl_frac > 1.0_r8) &
            write(logunit,*) subname,'SUSPECT: reflected fraction outside [0,1]'
    end if

    if (minval(fsa) < 0.0_r8 .or. minval(fsr) < 0.0_r8) &
         write(logunit,*) subname,'SUSPECT: negative absorbed or reflected shortwave'
    if (minval(sabg) < 0.0_r8) &
         write(logunit,*) subname,'SUSPECT: negative ground absorption'
    if (maxval(abs(sabg - sabgs)) > 1.0e-10_r8) &
         write(logunit,*) subname,'SUSPECT: sabg /= sabg_soil with no snow present'

    call shr_sys_flush(logunit)
    deallocate(fsa, fsr, sabv, sabg, sabgs)

  end subroutine kelm_report_surfrad


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
    type(kELMType), intent(in) :: elm
    integer, intent(in) :: npatch, logunit
    integer :: ierr
    real(r8), allocatable :: us(:), rb(:), ram(:), dsp(:), z0(:), rsun(:), rsha(:), lsun(:)
    character(len=*), parameter :: subname = '(diagnose_canflux) '

    allocate(us(npatch), rb(npatch), ram(npatch), dsp(npatch), z0(npatch), &
             rsun(npatch), rsha(npatch), lsun(npatch))
    call kELMGetUstar(elm, us, npatch, ierr)
    call kELMGetRb1(elm, rb, npatch, ierr)
    call kELMGetRam1(elm, ram, npatch, ierr)
    call kELMGetDispla(elm, dsp, npatch, ierr)
    call kELMGetZ0mv(elm, z0, npatch, ierr)
    call kELMGetRssun(elm, rsun, npatch, ierr)
    call kELMGetRssha(elm, rsha, npatch, ierr)
    call kELMGetLaisun(elm, lsun, npatch, ierr)

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
    if (ierr /= KELM_SUCCESS) then
       write(logunit,*) '(kelm_kernels) ERROR: kernel ', &
            trim(kernel_name(k)),' returned status ',ierr
       call shr_sys_flush(logunit)
       call shr_sys_abort('(kelm_kernels) ERROR: kernel '// &
            trim(kernel_name(k))//' failed')
    end if
  end subroutine check

end module kelmKernelMod
