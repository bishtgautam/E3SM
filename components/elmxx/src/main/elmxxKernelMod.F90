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
                               ELMxxComputeLakeHydrology

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

    case (K_SURFRAD, K_CANSUNSHADE)
       ! Both read the full canopy/ground albedo set -- albd, albi, fabd,
       ! fabi, ftdd, ftid, ftii, albgrd, albgri, albso*, albsn*_hst, and the
       ! SNICAR flx_abs* factors. Every one of those is SurfaceAlbedo's
       ! output, and SurfaceAlbedo is a Stage 5 port that does not exist.
       ! They also read forc_solad/forc_solai, which Stage 3 does not cross.
       why = 'needs SurfaceAlbedo (Stage 5) for albd/albi/fabd/fabi/ftdd/' // &
             'ftid/ftii and albgrd/albgri, plus 2-D forc_solad/forc_solai'

    case (K_CANTEMP, K_BAREGRND, K_CANFLUX, K_SOILTEMP, K_SOILFLUX, &
          K_SURFRUNOFF, K_ROOTWATER, K_HYDRODRAIN)
       ! All read the soil hydraulic properties (watsat, watfc, sucsat, bsw,
       ! smpmin) and the soil column state (h2osoi_liq, h2osoi_ice, dz,
       ! t_soisno). The hydraulic properties are NOT surfdata fields: ELM
       ! derives them from sand, clay and organic matter through
       ! iniTimeConst's pedotransfer functions, which has not been ported.
       ! Stage 2 stores the raw inputs per column, so the port has what it
       ! needs -- it just has not happened.
       why = 'needs the iniTimeConst pedotransfer port for watsat/watfc/' // &
             'sucsat/bsw and an initialized soil column (h2osoi_*, dz, t_soisno)'

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
       why = 'needs SLOPE/STD_ELEV read from surfdata and its cold-start ' // &
             'inputs seeded (micro_sigma, n_melt, frac_veg_nosno, dewmx, ' // &
             'h2ocan, fwet, h2osfc, int_snow, frac_h2osfc, frac_sno_eff)'

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
