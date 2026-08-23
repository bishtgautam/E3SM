module elmxxDiagnosticsMod
  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Writes ELMxx state in ELM's own ELMDIAG1 binary format, so a free-running
  ! ELMxx run can be diffed against ELM's elm_diagnostics.bin with the tooling
  ! that already exists for the replay tests.
  !
  ! This is deliberately a byte-for-byte mirror of ElmDiagnostics.F90 rather
  ! than a new format: everything that reads ELM's binary then reads this one
  ! unchanged.
  !
  ! Records are written in ELMxx's PACKED index space (0..n_kokkos_col-1), not
  ! ELM's global column space.  The maps needed to align the two are written
  ! once per run under 'elmxxmap:' so the comparison tool can do it.
  !
  ! Snapshot point: the top of the timestep, before any kernel has run.  ELM's
  ! equivalent anchor is 'canhydro_in:' — CanopyHydrology is its first kernel —
  ! so 'elmxx_in:<var>' and 'canhydro_in:<var>' describe the same instant and
  ! are directly comparable.  That is what makes this an error-growth trace
  ! rather than another single-step replay.
  !-----------------------------------------------------------------------

  use shr_kind_mod, only : r8 => shr_kind_r8
  use shr_sys_mod , only : shr_sys_abort
  use elmxx_mod   , only : ELMxxType, ELMXX_SUCCESS,                        &
                           ELMxxGetTSoisno, ELMxxGetH2osoiLiqSoi,           &
                           ELMxxGetTSoisnoSno, ELMxxGetH2osoiLiqSno,        &
                           ELMxxGetH2osoiIceSno, ELMxxGetDzSno,             &
                           ELMxxGetSnwRds,                                  &
                           ELMxxGetDz, ELMxxGetH2osoiIce,                   &
                           ELMxxGetH2osoiLiq, ELMxxGetSnl,                  &
                           ELMxxGetFracSnoEff,                              &
                           ELMxxGetHsSoil, ELMxxGetHsTopSnow,               &
                           ELMxxGetHsH2osfc, ELMxxGetDhsdT,                 &
                           ELMxxGetEflxShSnow, ELMxxGetQflxEvSnow,          &
                           ELMxxGetDlrad, ELMxxGetEflxShSoil,                &
                           ELMxxGetQflxEvSoil, ELMxxGetTVeg,                 &
                           ELMxxGetQgSnow, ELMxxGetQgSoil, ELMxxGetThm,      &
                           ELMxxGetQg, ELMxxGetHtvp, ELMxxGetSoilbeta,       &
                           ELMxxGetZ0mg, ELMxxGetThv,                        &
                           ELMxxGetQflxEvapGrndCol, ELMxxGetQflxEvapGrnd,    &
                           ELMxxGetTSsbef, ELMxxGetCgrndl, ELMxxGetCgrnds,   &
                           ELMxxGetQflxTopSoilCol, ELMxxGetFracH2osfc,       &
                           ELMxxGetForcRhoCol, ELMxxGetDqgdT, ELMxxGetZii,   &
                           ELMxxGetSweOld, ELMxxGetImeltReal,                &
                           ELMxxGetQflxDewSnowCol,                          &
                           ELMxxGetSabgSoil, ELMxxGetSabgSnow,               &
                           ELMxxGetSabgLyr,                                  &
                           ELMxxGetH2osoiIceSoi, ELMxxGetTGrnd,             &
                           ELMxxGetTH2osfc, ELMxxGetH2osfc,                 &
                           ELMxxGetH2osno, ELMxxGetSnowDepth,               &
                           ELMxxGetFracSno, ELMxxGetFracH2osfc,             &
                           ELMxxGetIntSnow, ELMxxGetSnl,                    &
                           ELMxxGetTVeg, ELMxxGetBtran, ELMxxGetH2ocan,     &
                           ELMxxGetEflxShTot, ELMxxGetEflxLhTot,            &
                           ELMxxGetQflxEvapTot, ELMxxGetFsa, ELMxxGetFsr,   &
                           ELMxxGetTRef2m, ELMxxGetQflxTranVeg,             &
                           ELMxxGetQflxInflCol, ELMxxGetQflxSurfCol,        &
                           ELMxxGetQflxDrainCol, ELMxxGetQflxEvapSoi,       &
                           ELMxxGetQflxTopSoilCol, ELMxxGetFsat,            &
                           ELMxxGetFcov, ELMxxGetZwt, ELMxxGetWtfact,       &
                           ELMxxGetEffPorosity, ELMxxGetAlbgrd,             &
                           ELMxxGetAlbgri, ELMxxGetAlbsod, ELMxxGetAlbd,     &
                           ELMxxGetElai
  use elmxxKokkosStateMod, only : n_kokkos_col, n_kokkos_patch,             &
                                  col_of_kcol, patch_of_kpatch,            &
                                  kcol_of_col
  use elmxxSoilPropMod   , only : watsat, bsw, sucsat, hksat, nlevgrnd, nlevsno
  use elmxxRootMod       , only : rootr, btran_root => btran
  use elmxxSurfaceStateMod, only : patch_lai, patch_sai
  implicit none
  private

  logical, public :: elmxx_diag_enabled  = .false.
  integer, public :: elmxx_diag_timestep = 0

  integer, parameter :: LABEL_LEN = 64
  integer            :: dunit     = 211
  logical            :: is_open   = .false.

  public :: elmxx_diag_init
  public :: elmxx_diag_finalize
  public :: elmxx_diag_new_timestep
  public :: elmxx_diag_1d
  public :: elmxx_diag_2d
  public :: elmxx_diag_int_1d
  public :: elmxx_diag_snapshot_state
  public :: elmxx_diag_snapshot_fluxes
  public :: elmxx_diag_snapshot_preinfil
  public :: elmxx_diag_snapshot_presoilflux
  public :: elmxx_diag_snapshot_presoiltemp
  public :: elmxx_diag_snapshot_soilwater
  public :: elmxx_diag_write_maps

contains

  subroutine elmxx_diag_init(output_path, enabled)
    character(len=*), intent(in) :: output_path
    logical,          intent(in) :: enabled

    elmxx_diag_enabled  = enabled
    elmxx_diag_timestep = 0
    if (.not. elmxx_diag_enabled) return

    open(unit=dunit, file=trim(output_path), &
         form='unformatted', access='stream', action='write', status='replace')
    is_open = .true.

    write(dunit) 'ELMDIAG1'
    write(dunit) LABEL_LEN
  end subroutine elmxx_diag_init

  subroutine elmxx_diag_finalize()
    if (is_open) then
       close(dunit)
       is_open = .false.
    end if
    elmxx_diag_enabled = .false.
  end subroutine elmxx_diag_finalize

  subroutine elmxx_diag_new_timestep(nstep)
    ! Label records by the driver's nstep, matching what ElmDiagnostics now
    ! does, so ELMxx step N and ELM step N are the same instant.
    integer, intent(in), optional :: nstep
    if (.not. elmxx_diag_enabled) return
    if (present(nstep)) then
       elmxx_diag_timestep = nstep
    else
       elmxx_diag_timestep = elmxx_diag_timestep + 1
    end if
  end subroutine elmxx_diag_new_timestep

  subroutine elmxx_diag_1d(label, array, n)
    character(len=*), intent(in) :: label
    real(r8),         intent(in) :: array(:)
    integer,          intent(in) :: n
    character(len=LABEL_LEN) :: padded
    if (.not. elmxx_diag_enabled) return
    if (n <= 0) return
    padded = label
    write(dunit) elmxx_diag_timestep
    write(dunit) padded
    write(dunit) 1
    write(dunit) n
    write(dunit) array(1:n)
  end subroutine elmxx_diag_1d

  subroutine elmxx_diag_2d(label, array, n1, n2)
    character(len=*), intent(in) :: label
    real(r8),         intent(in) :: array(:,:)
    integer,          intent(in) :: n1, n2
    character(len=LABEL_LEN) :: padded
    if (.not. elmxx_diag_enabled) return
    if (n1 <= 0 .or. n2 <= 0) return
    padded = label
    write(dunit) elmxx_diag_timestep
    write(dunit) padded
    write(dunit) 2
    write(dunit) n1, n2
    write(dunit) array(1:n1, 1:n2)
  end subroutine elmxx_diag_2d

  subroutine elmxx_diag_int_1d(label, array, n)
    character(len=*), intent(in) :: label
    integer,          intent(in) :: array(:)
    integer,          intent(in) :: n
    character(len=LABEL_LEN) :: padded
    if (.not. elmxx_diag_enabled) return
    if (n <= 0) return
    padded = label
    write(dunit) elmxx_diag_timestep
    write(dunit) padded
    write(dunit) -1
    write(dunit) n
    write(dunit) array(1:n)
  end subroutine elmxx_diag_int_1d

  !-----------------------------------------------------------------------
  ! Write the packed -> ELM index maps once, so the comparison tool can put
  ! ELMxx's columns and patches back where ELM has them.
  !-----------------------------------------------------------------------
  subroutine elmxx_diag_write_maps()
    integer :: kc, c, j
    real(r8), allocatable :: tmp(:,:)
    if (.not. elmxx_diag_enabled) return
    if (associated(col_of_kcol)) &
         call elmxx_diag_int_1d('elmxxmap:col_of_kcol', col_of_kcol, n_kokkos_col)
    if (associated(patch_of_kpatch)) &
         call elmxx_diag_int_1d('elmxxmap:patch_of_kpatch', patch_of_kpatch, n_kokkos_patch)

    ! Soil hydraulic properties, once. Static, but they set the matric
    ! potential that drives btran and root extraction, so a wrong value here
    ! shows up as a moisture drift and nowhere else.
    if (associated(watsat) .and. n_kokkos_col > 0) then
       allocate(tmp(n_kokkos_col, nlevgrnd))
       do kc = 1, n_kokkos_col
          c = col_of_kcol(kc)
          do j = 1, nlevgrnd
             tmp(kc,j) = watsat(c,j)
          end do
       end do
       call elmxx_diag_2d('elmxxsoil:watsat', tmp, n_kokkos_col, nlevgrnd)
       do kc = 1, n_kokkos_col
          c = col_of_kcol(kc)
          tmp(kc,1:nlevgrnd) = bsw(c,1:nlevgrnd)
       end do
       call elmxx_diag_2d('elmxxsoil:bsw', tmp, n_kokkos_col, nlevgrnd)
       do kc = 1, n_kokkos_col
          c = col_of_kcol(kc)
          tmp(kc,1:nlevgrnd) = sucsat(c,1:nlevgrnd)
       end do
       call elmxx_diag_2d('elmxxsoil:sucsat', tmp, n_kokkos_col, nlevgrnd)
       do kc = 1, n_kokkos_col
          c = col_of_kcol(kc)
          tmp(kc,1:nlevgrnd) = hksat(c,1:nlevgrnd)
       end do
       call elmxx_diag_2d('elmxxsoil:hksat', tmp, n_kokkos_col, nlevgrnd)
       deallocate(tmp)
    end if
  end subroutine elmxx_diag_write_maps

  !-----------------------------------------------------------------------
  ! Snapshot the carried state at the top of a timestep, before any kernel
  ! has run.  ELM's matching anchor is 'canhydro_in:' — its first kernel — so
  ! these records line up instant-for-instant with ELM's own.
  !-----------------------------------------------------------------------
  subroutine elmxx_diag_snapshot_state(elm, nlevtot, nlevgrnd, tag)
    type(ELMxxType) , intent(in) :: elm
    integer         , intent(in) :: nlevtot, nlevgrnd
    character(len=*), intent(in) :: tag

    integer  :: ierr, sz(2), szg(2), szs(2)
    real(r8), allocatable :: c1(:), p1(:), c2(:,:), cg(:,:), cs(:,:)
    integer , allocatable :: ci(:)

    if (.not. elmxx_diag_enabled) return
    if (n_kokkos_col <= 0) return

    allocate(c1(n_kokkos_col), ci(n_kokkos_col))
    allocate(c2(n_kokkos_col, nlevtot))
    allocate(cg(n_kokkos_col, nlevgrnd))
    allocate(cs(n_kokkos_col, nlevsno))
    allocate(p1(max(n_kokkos_patch,1)))

    ! ---- column scalars ----
    call get_c('t_grnd',      ELMxxGetTGrnd)
    call get_c('t_h2osfc',    ELMxxGetTH2osfc)
    call get_c('h2osfc',      ELMxxGetH2osfc)
    call get_c('h2osno',      ELMxxGetH2osno)
    call get_c('snow_depth',  ELMxxGetSnowDepth)
    call get_c('frac_sno',    ELMxxGetFracSno)
    call get_c('frac_h2osfc', ELMxxGetFracH2osfc)
    call get_c('int_snow',    ELMxxGetIntSnow)

    call ELMxxGetSnl(elm, ci, n_kokkos_col, ierr)
    if (ierr == ELMXX_SUCCESS) call elmxx_diag_int_1d(tag//':snl', ci, n_kokkos_col)

    ! ---- column profiles ----
    sz(1) = n_kokkos_col; sz(2) = nlevtot
    call get_c2('t_soisno',   ELMxxGetTSoisno)
    szg(1) = n_kokkos_col; szg(2) = nlevgrnd
    call get_cg('h2osoi_liq_soi', ELMxxGetH2osoiLiqSoi)
    call get_cg('h2osoi_ice_soi', ELMxxGetH2osoiIceSoi)

    ! The snow half. This used to be omitted, on the reasoning that
    ! "1x1_brazil never carries snow, so nothing is lost" -- which stopped
    ! being true the moment 1x1_glc was used, and left the snowpack
    ! unobservable exactly where the divergence lives. Reported in ELM order
    ! so it lines up with the restart file and elm_diagnostics.bin.
    szs(1) = n_kokkos_col; szs(2) = nlevsno
    call get_cs('t_soisno_sno',   ELMxxGetTSoisnoSno)
    call get_cs('h2osoi_liq_sno', ELMxxGetH2osoiLiqSno)
    call get_cs('h2osoi_ice_sno', ELMxxGetH2osoiIceSno)
    call get_cs('dz_sno',         ELMxxGetDzSno)
    call get_cs('snw_rds',        ELMxxGetSnwRds)

    ! ---- patch scalars ----
    if (n_kokkos_patch > 0) then
       call get_p('t_veg',         ELMxxGetTVeg)
       call get_p('btran',         ELMxxGetBtran)
       call get_p('h2ocan',        ELMxxGetH2ocan)
       call get_p('t_ref2m',       ELMxxGetTRef2m)
       call get_p('eflx_sh_tot',   ELMxxGetEflxShTot)
       call get_p('eflx_lh_tot',   ELMxxGetEflxLhTot)
       call get_p('qflx_evap_tot', ELMxxGetQflxEvapTot)
       ! Water budget terms. Recorded so a moisture drift can be attributed to
       ! a flux rather than inferred from the storage profile.
       call get_p('qflx_tran_veg', ELMxxGetQflxTranVeg)
       call get_p('qflx_evap_soi', ELMxxGetQflxEvapSoi)
    end if

    deallocate(c1, ci, c2, cg, p1)

  contains

    subroutine get_c(name, getter)
      character(len=*), intent(in) :: name
      external :: getter
      call getter(elm, c1, n_kokkos_col, ierr)
      if (ierr == ELMXX_SUCCESS) call elmxx_diag_1d(tag//':'//name, c1, n_kokkos_col)
    end subroutine get_c

    subroutine get_p(name, getter)
      character(len=*), intent(in) :: name
      external :: getter
      call getter(elm, p1, n_kokkos_patch, ierr)
      if (ierr == ELMXX_SUCCESS) call elmxx_diag_1d(tag//':'//name, p1, n_kokkos_patch)
    end subroutine get_p

    subroutine get_cg(name, getter)
      character(len=*), intent(in) :: name
      external :: getter
      call getter(elm, cg, szg, ierr)
      if (ierr == ELMXX_SUCCESS) &
           call elmxx_diag_2d(tag//':'//name, cg, n_kokkos_col, nlevgrnd)
    end subroutine get_cg

    subroutine get_cs(name, getter)
      character(len=*), intent(in) :: name
      external :: getter
      call getter(elm, cs, szs, ierr)
      if (ierr == ELMXX_SUCCESS) &
           call elmxx_diag_2d(tag//':'//name, cs, n_kokkos_col, nlevsno)
    end subroutine get_cs

    subroutine get_c2(name, getter)
      character(len=*), intent(in) :: name
      external :: getter
      call getter(elm, c2, sz, ierr)
      if (ierr == ELMXX_SUCCESS) &
           call elmxx_diag_2d(tag//':'//name, c2, n_kokkos_col, nlevtot)
    end subroutine get_c2

  end subroutine elmxx_diag_snapshot_state

  !-----------------------------------------------------------------------
  ! Water fluxes, sampled at the END of the step. The HydrologyDrainage
  ! getters read elm->hydrologyDrainage, which is a separate structure from
  ! the natural column and holds nothing at the top of a step -- sampling
  ! there returns zeros and reads as "no infiltration ever happened".
  !-----------------------------------------------------------------------
  subroutine elmxx_diag_snapshot_fluxes(elm, tag)
    type(ELMxxType) , intent(in) :: elm
    character(len=*), intent(in) :: tag
    integer :: ierr
    real(r8), allocatable :: c1(:)

    if (.not. elmxx_diag_enabled) return
    if (n_kokkos_col <= 0) return
    allocate(c1(n_kokkos_col))

    call ELMxxGetQflxInflCol(elm, c1, n_kokkos_col, ierr)
    if (ierr == ELMXX_SUCCESS) call elmxx_diag_1d(tag//':qflx_infl', c1, n_kokkos_col)
    call ELMxxGetQflxSurfCol(elm, c1, n_kokkos_col, ierr)
    if (ierr == ELMXX_SUCCESS) call elmxx_diag_1d(tag//':qflx_surf', c1, n_kokkos_col)
    call ELMxxGetQflxDrainCol(elm, c1, n_kokkos_col, ierr)
    if (ierr == ELMXX_SUCCESS) call elmxx_diag_1d(tag//':qflx_drain', c1, n_kokkos_col)

    ! Absorbed/reflected shortwave, sampled here rather than at the top of the
    ! step: ELM records surfrad_out: DURING the step, so a top-of-step ELMxx
    ! value is one timestep stale and reads as a whole-diurnal-cycle shift.
    block
      real(r8), allocatable :: p1b(:)
      if (n_kokkos_patch > 0) then
         allocate(p1b(n_kokkos_patch))
         call ELMxxGetFsa(elm, p1b, n_kokkos_patch, ierr)
         if (ierr == ELMXX_SUCCESS) call elmxx_diag_1d(tag//':fsa', p1b, n_kokkos_patch)
         call ELMxxGetFsr(elm, p1b, n_kokkos_patch, ierr)
         if (ierr == ELMXX_SUCCESS) call elmxx_diag_1d(tag//':fsr', p1b, n_kokkos_patch)
         deallocate(p1b)
      end if
    end block

    ! Phenology, so leaf area can be compared directly rather than inferred
    ! from whether the canopy intercepted anything.
    if (associated(patch_lai) .and. n_kokkos_patch > 0) then
       block
         real(r8), allocatable :: pl(:)
         integer :: kp
         allocate(pl(n_kokkos_patch))
         ! Read the DEVICE elai, not the Fortran patch_lai -- phenology now
         ! interpolates on the device and patch_lai is no longer maintained.
         call ELMxxGetElai(elm, pl, n_kokkos_patch, ierr)
         if (ierr == ELMXX_SUCCESS) call elmxx_diag_1d(tag//':elai', pl, n_kokkos_patch)
         do kp = 1, n_kokkos_patch
            pl(kp) = patch_sai(patch_of_kpatch(kp))
         end do
         call elmxx_diag_1d(tag//':esai', pl, n_kokkos_patch)
         deallocate(pl)
       end block
    end if

    ! Root water uptake partition. rootr is SoilWater's per-layer sink, so it
    ! needs to be gradeable on its own rather than through btran.
    if (associated(rootr) .and. n_kokkos_patch > 0) then
       block
         real(r8), allocatable :: pg(:,:)
         integer :: kp, szp(2)
         allocate(pg(n_kokkos_patch, nlevgrnd))
         do kp = 1, n_kokkos_patch
            pg(kp,1:nlevgrnd) = rootr(patch_of_kpatch(kp), 1:nlevgrnd)
         end do
         szp(1) = n_kokkos_patch; szp(2) = nlevgrnd
         call elmxx_diag_2d(tag//':rootr', pg, n_kokkos_patch, nlevgrnd)
         deallocate(pg)
       end block
    end if

    ! Albedo, so either path -- Fortran or C++ -- can be diffed against ELM's
    ! surfrad_in records rather than against each other.
    block
      real(r8), allocatable :: c2(:,:), p2(:,:)
      integer :: szc(2), szp(2)
      allocate(c2(n_kokkos_col, 2), p2(max(n_kokkos_patch,1), 2))
      szc(1) = n_kokkos_col; szc(2) = 2
      call ELMxxGetAlbgrd(elm, c2, szc, ierr)
      if (ierr == ELMXX_SUCCESS) call elmxx_diag_2d(tag//':albgrd', c2, n_kokkos_col, 2)
      call ELMxxGetAlbgri(elm, c2, szc, ierr)
      if (ierr == ELMXX_SUCCESS) call elmxx_diag_2d(tag//':albgri', c2, n_kokkos_col, 2)
      call ELMxxGetAlbsod(elm, c2, szc, ierr)
      if (ierr == ELMXX_SUCCESS) call elmxx_diag_2d(tag//':albsod', c2, n_kokkos_col, 2)
      if (n_kokkos_patch > 0) then
         szp(1) = n_kokkos_patch; szp(2) = 2
         call ELMxxGetAlbd(elm, p2, szp, ierr)
         if (ierr == ELMXX_SUCCESS) call elmxx_diag_2d(tag//':albd', p2, n_kokkos_patch, 2)
      end if
      deallocate(c2, p2)
    end block

    ! The quantities that set the infiltration/runoff split.
    call ELMxxGetQflxTopSoilCol(elm, c1, n_kokkos_col, ierr)
    if (ierr == ELMXX_SUCCESS) call elmxx_diag_1d(tag//':qflx_top_soil', c1, n_kokkos_col)
    call ELMxxGetFsat(elm, c1, n_kokkos_col, ierr)
    if (ierr == ELMXX_SUCCESS) call elmxx_diag_1d(tag//':fsat', c1, n_kokkos_col)
    call ELMxxGetFcov(elm, c1, n_kokkos_col, ierr)
    if (ierr == ELMXX_SUCCESS) call elmxx_diag_1d(tag//':fcov', c1, n_kokkos_col)
    call ELMxxGetZwt(elm, c1, n_kokkos_col, ierr)
    if (ierr == ELMXX_SUCCESS) call elmxx_diag_1d(tag//':zwt', c1, n_kokkos_col)
    call ELMxxGetWtfact(elm, c1, n_kokkos_col, ierr)
    if (ierr == ELMXX_SUCCESS) call elmxx_diag_1d(tag//':wtfact', c1, n_kokkos_col)

    block
      real(r8), allocatable :: cg(:,:)
      integer :: szg(2)
      allocate(cg(n_kokkos_col, nlevgrnd))
      szg(1) = n_kokkos_col; szg(2) = nlevgrnd
      call ELMxxGetEffPorosity(elm, cg, szg, ierr)
      if (ierr == ELMXX_SUCCESS) &
           call elmxx_diag_2d(tag//':eff_porosity', cg, n_kokkos_col, nlevgrnd)
      deallocate(cg)
    end block

    deallocate(c1)
  end subroutine elmxx_diag_snapshot_fluxes

  !-----------------------------------------------------------------------
  ! Soil liquid immediately after the Richards solve, before HydrologyDrainage
  ! touches it. ELM's matching record is soilwater_out:h2osoi_liq.
  !-----------------------------------------------------------------------
  ! Temporary (G5): the state the SoilTemperature expand reads, captured
  ! immediately before the solve, so it can be diffed against ELM's
  ! soiltemp_in:* record for the same step.
  ! Temporary (G6): the state SoilFluxes differences against, captured
  ! between soiltemp and soilflux so it lines up with ELM's soilflx_in record.
  ! Temporary (G6): the terms that build qflx_infl, captured after snowwater
  ! and before surfrunoff.
  subroutine elmxx_diag_snapshot_preinfil(elm, tag)
    type(ELMxxType) , intent(in) :: elm
    character(len=*), intent(in) :: tag
    integer :: ierr
    real(r8), allocatable :: c1(:)
    if (.not. elmxx_diag_enabled) return
    if (n_kokkos_col <= 0) return
    allocate(c1(n_kokkos_col))
    call ELMxxGetQflxEvapGrndCol(elm, c1, n_kokkos_col, ierr)
    if (ierr == ELMXX_SUCCESS) call elmxx_diag_1d(tag//':qflx_evap_grnd_col', c1, n_kokkos_col)
    call ELMxxGetForcRhoCol(elm, c1, n_kokkos_col, ierr)
    if (ierr == ELMXX_SUCCESS) call elmxx_diag_1d(tag//':forc_rho', c1, n_kokkos_col)
    call ELMxxGetDqgdT(elm, c1, n_kokkos_col, ierr)
    if (ierr == ELMXX_SUCCESS) call elmxx_diag_1d(tag//':dqgdT', c1, n_kokkos_col)
    call ELMxxGetZii(elm, c1, n_kokkos_col, ierr)
    if (ierr == ELMXX_SUCCESS) call elmxx_diag_1d(tag//':zii', c1, n_kokkos_col)
    call ELMxxGetQflxTopSoilCol(elm, c1, n_kokkos_col, ierr)
    if (ierr == ELMXX_SUCCESS) call elmxx_diag_1d(tag//':qflx_top_soil', c1, n_kokkos_col)
    call ELMxxGetFracH2osfc(elm, c1, n_kokkos_col, ierr)
    if (ierr == ELMXX_SUCCESS) call elmxx_diag_1d(tag//':frac_h2osfc', c1, n_kokkos_col)
    call ELMxxGetFracSnoEff(elm, c1, n_kokkos_col, ierr)
    if (ierr == ELMXX_SUCCESS) call elmxx_diag_1d(tag//':frac_sno_eff', c1, n_kokkos_col)
    call ELMxxGetQflxDewSnowCol(elm, c1, n_kokkos_col, ierr)
    if (ierr == ELMXX_SUCCESS) call elmxx_diag_1d(tag//':qflx_dew_snow_col', c1, n_kokkos_col)
    block
      real(r8), allocatable :: cg(:,:)
      integer :: szt(2)
      allocate(cg(n_kokkos_col, nlevsno+nlevgrnd))
      szt(1) = n_kokkos_col; szt(2) = nlevsno + nlevgrnd
      call ELMxxGetSweOld(elm, cg, szt, ierr)
      if (ierr == ELMXX_SUCCESS) call elmxx_diag_2d(tag//':swe_old', cg, n_kokkos_col, nlevsno+nlevgrnd)
      call ELMxxGetImeltReal(elm, cg, szt, ierr)
      if (ierr == ELMXX_SUCCESS) call elmxx_diag_2d(tag//':imelt', cg, n_kokkos_col, nlevsno+nlevgrnd)
      deallocate(cg)
    end block
    block
      real(r8), allocatable :: p1(:)
      allocate(p1(n_kokkos_patch))
      call ELMxxGetQflxEvapGrnd(elm, p1, n_kokkos_patch, ierr)
      if (ierr == ELMXX_SUCCESS) call elmxx_diag_1d(tag//':qflx_evap_grnd', p1, n_kokkos_patch)
      deallocate(p1)
    end block
    deallocate(c1)
  end subroutine elmxx_diag_snapshot_preinfil

  subroutine elmxx_diag_snapshot_presoilflux(elm, nlevtot_in, tag)
    type(ELMxxType) , intent(in) :: elm
    integer         , intent(in) :: nlevtot_in
    character(len=*), intent(in) :: tag
    integer :: ierr, szg(2)
    real(r8), allocatable :: cg(:,:)
    real(r8), allocatable :: c1(:), p1(:)
    if (.not. elmxx_diag_enabled) return
    if (n_kokkos_col <= 0) return
    allocate(cg(n_kokkos_col, nlevtot_in), c1(n_kokkos_col), p1(n_kokkos_patch))
    szg(1) = n_kokkos_col; szg(2) = nlevtot_in
    call ELMxxGetTSsbef(elm, cg, szg, ierr)
    if (ierr == ELMXX_SUCCESS) call elmxx_diag_2d(tag//':tssbef', cg, n_kokkos_col, nlevtot_in)
    call ELMxxGetTSoisno(elm, cg, szg, ierr)
    if (ierr == ELMXX_SUCCESS) call elmxx_diag_2d(tag//':t_soisno', cg, n_kokkos_col, nlevtot_in)
    call ELMxxGetTGrnd(elm, c1, n_kokkos_col, ierr)
    if (ierr == ELMXX_SUCCESS) call elmxx_diag_1d(tag//':t_grnd', c1, n_kokkos_col)
    call ELMxxGetCgrndl(elm, p1, n_kokkos_patch, ierr)
    if (ierr == ELMXX_SUCCESS) call elmxx_diag_1d(tag//':cgrndl', p1, n_kokkos_patch)
    call ELMxxGetCgrnds(elm, p1, n_kokkos_patch, ierr)
    if (ierr == ELMXX_SUCCESS) call elmxx_diag_1d(tag//':cgrnds', p1, n_kokkos_patch)
    deallocate(cg, c1, p1)
  end subroutine elmxx_diag_snapshot_presoilflux

  subroutine elmxx_diag_snapshot_presoiltemp(elm, nlevtot_in, tag)
    type(ELMxxType) , intent(in) :: elm
    integer         , intent(in) :: nlevtot_in
    character(len=*), intent(in) :: tag
    integer :: ierr, szg(2)
    real(r8), allocatable :: cg(:,:)
    real(r8), allocatable :: c1(:)
    integer , allocatable :: i1(:)
    if (.not. elmxx_diag_enabled) return
    if (n_kokkos_col <= 0) return
    allocate(cg(n_kokkos_col, nlevtot_in), c1(n_kokkos_col), i1(n_kokkos_col))
    szg(1) = n_kokkos_col; szg(2) = nlevtot_in
    call ELMxxGetTSoisno(elm, cg, szg, ierr)
    if (ierr == ELMXX_SUCCESS) call elmxx_diag_2d(tag//':t_soisno', cg, n_kokkos_col, nlevtot_in)
    call ELMxxGetDz(elm, cg, szg, ierr)
    if (ierr == ELMXX_SUCCESS) call elmxx_diag_2d(tag//':dz', cg, n_kokkos_col, nlevtot_in)
    call ELMxxGetH2osoiIce(elm, cg, szg, ierr)
    if (ierr == ELMXX_SUCCESS) call elmxx_diag_2d(tag//':h2osoi_ice', cg, n_kokkos_col, nlevtot_in)
    call ELMxxGetH2osoiLiq(elm, cg, szg, ierr)
    if (ierr == ELMXX_SUCCESS) call elmxx_diag_2d(tag//':h2osoi_liq', cg, n_kokkos_col, nlevtot_in)
    call ELMxxGetSnl(elm, i1, n_kokkos_col, ierr)
    if (ierr == ELMXX_SUCCESS) call elmxx_diag_int_1d(tag//':snl', i1, n_kokkos_col)
    call ELMxxGetFracSnoEff(elm, c1, n_kokkos_col, ierr)
    if (ierr == ELMXX_SUCCESS) call elmxx_diag_1d(tag//':frac_sno_eff', c1, n_kokkos_col)
    call ELMxxGetHsSoil(elm, c1, n_kokkos_col, ierr)
    if (ierr == ELMXX_SUCCESS) call elmxx_diag_1d(tag//':hs_soil', c1, n_kokkos_col)
    call ELMxxGetHsTopSnow(elm, c1, n_kokkos_col, ierr)
    if (ierr == ELMXX_SUCCESS) call elmxx_diag_1d(tag//':hs_top_snow', c1, n_kokkos_col)
    call ELMxxGetHsH2osfc(elm, c1, n_kokkos_col, ierr)
    if (ierr == ELMXX_SUCCESS) call elmxx_diag_1d(tag//':hs_h2osfc', c1, n_kokkos_col)
    call ELMxxGetDhsdT(elm, c1, n_kokkos_col, ierr)
    if (ierr == ELMXX_SUCCESS) call elmxx_diag_1d(tag//':dhsdT', c1, n_kokkos_col)
    call ELMxxGetQgSnow(elm, c1, n_kokkos_col, ierr)
    if (ierr == ELMXX_SUCCESS) call elmxx_diag_1d(tag//':qg_snow', c1, n_kokkos_col)
    call ELMxxGetQgSoil(elm, c1, n_kokkos_col, ierr)
    if (ierr == ELMXX_SUCCESS) call elmxx_diag_1d(tag//':qg_soil', c1, n_kokkos_col)
    call ELMxxGetTGrnd(elm, c1, n_kokkos_col, ierr)
    if (ierr == ELMXX_SUCCESS) call elmxx_diag_1d(tag//':t_grnd', c1, n_kokkos_col)
    call ELMxxGetQg(elm, c1, n_kokkos_col, ierr)
    if (ierr == ELMXX_SUCCESS) call elmxx_diag_1d(tag//':qg', c1, n_kokkos_col)
    call ELMxxGetHtvp(elm, c1, n_kokkos_col, ierr)
    if (ierr == ELMXX_SUCCESS) call elmxx_diag_1d(tag//':htvp', c1, n_kokkos_col)
    call ELMxxGetSoilbeta(elm, c1, n_kokkos_col, ierr)
    if (ierr == ELMXX_SUCCESS) call elmxx_diag_1d(tag//':soilbeta', c1, n_kokkos_col)
    call ELMxxGetZ0mg(elm, c1, n_kokkos_col, ierr)
    if (ierr == ELMXX_SUCCESS) call elmxx_diag_1d(tag//':z0mg', c1, n_kokkos_col)
    call ELMxxGetThv(elm, c1, n_kokkos_col, ierr)
    if (ierr == ELMXX_SUCCESS) call elmxx_diag_1d(tag//':thv', c1, n_kokkos_col)
    call ELMxxGetQflxEvapGrndCol(elm, c1, n_kokkos_col, ierr)
    if (ierr == ELMXX_SUCCESS) call elmxx_diag_1d(tag//':qflx_evap_grnd_col', c1, n_kokkos_col)
    block
      real(r8), allocatable :: p1(:)
      allocate(p1(n_kokkos_patch))
      call ELMxxGetEflxShSnow(elm, p1, n_kokkos_patch, ierr)
      if (ierr == ELMXX_SUCCESS) call elmxx_diag_1d(tag//':eflx_sh_snow', p1, n_kokkos_patch)
      call ELMxxGetQflxEvSnow(elm, p1, n_kokkos_patch, ierr)
      if (ierr == ELMXX_SUCCESS) call elmxx_diag_1d(tag//':qflx_ev_snow', p1, n_kokkos_patch)
      call ELMxxGetDlrad(elm, p1, n_kokkos_patch, ierr)
      if (ierr == ELMXX_SUCCESS) call elmxx_diag_1d(tag//':dlrad', p1, n_kokkos_patch)
      call ELMxxGetEflxShSoil(elm, p1, n_kokkos_patch, ierr)
      if (ierr == ELMXX_SUCCESS) call elmxx_diag_1d(tag//':eflx_sh_soil', p1, n_kokkos_patch)
      call ELMxxGetQflxEvSoil(elm, p1, n_kokkos_patch, ierr)
      if (ierr == ELMXX_SUCCESS) call elmxx_diag_1d(tag//':qflx_ev_soil', p1, n_kokkos_patch)
      call ELMxxGetTVeg(elm, p1, n_kokkos_patch, ierr)
      if (ierr == ELMXX_SUCCESS) call elmxx_diag_1d(tag//':t_veg', p1, n_kokkos_patch)
      call ELMxxGetThm(elm, p1, n_kokkos_patch, ierr)
      if (ierr == ELMXX_SUCCESS) call elmxx_diag_1d(tag//':thm', p1, n_kokkos_patch)
      call ELMxxGetQflxEvapGrnd(elm, p1, n_kokkos_patch, ierr)
      if (ierr == ELMXX_SUCCESS) call elmxx_diag_1d(tag//':qflx_evap_grnd', p1, n_kokkos_patch)
      call ELMxxGetSabgSoil(elm, p1, n_kokkos_patch, ierr)
      if (ierr == ELMXX_SUCCESS) call elmxx_diag_1d(tag//':sabg_soil', p1, n_kokkos_patch)
      call ELMxxGetSabgSnow(elm, p1, n_kokkos_patch, ierr)
      if (ierr == ELMXX_SUCCESS) call elmxx_diag_1d(tag//':sabg_snow', p1, n_kokkos_patch)
      deallocate(p1)
    end block
    block
      real(r8), allocatable :: plyr(:,:)
      integer :: szp(2)
      allocate(plyr(n_kokkos_patch, nlevsno+1))
      szp(1) = n_kokkos_patch; szp(2) = nlevsno + 1
      call ELMxxGetSabgLyr(elm, plyr, szp, ierr)
      if (ierr == ELMXX_SUCCESS) call elmxx_diag_2d(tag//':sabg_lyr', plyr, n_kokkos_patch, nlevsno+1)
      deallocate(plyr)
    end block
    deallocate(cg, c1, i1)
  end subroutine elmxx_diag_snapshot_presoiltemp

  subroutine elmxx_diag_snapshot_soilwater(elm, nlevgrnd_in, tag)
    type(ELMxxType) , intent(in) :: elm
    integer         , intent(in) :: nlevgrnd_in
    character(len=*), intent(in) :: tag
    integer :: ierr, szg(2)
    real(r8), allocatable :: cg(:,:)
    if (.not. elmxx_diag_enabled) return
    if (n_kokkos_col <= 0) return
    allocate(cg(n_kokkos_col, nlevgrnd_in))
    szg(1) = n_kokkos_col; szg(2) = nlevgrnd_in
    call ELMxxGetH2osoiLiqSoi(elm, cg, szg, ierr)
    if (ierr == ELMXX_SUCCESS) &
         call elmxx_diag_2d(tag//':h2osoi_liq', cg, n_kokkos_col, nlevgrnd_in)
    deallocate(cg)
  end subroutine elmxx_diag_snapshot_soilwater

end module elmxxDiagnosticsMod
