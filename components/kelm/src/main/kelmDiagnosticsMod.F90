module kelmDiagnosticsMod
  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Writes kELM state in ELM's own ELMDIAG1 binary format, so a free-running
  ! kELM run can be diffed against ELM's elm_diagnostics.bin with the tooling
  ! that already exists for the replay tests.
  !
  ! This is deliberately a byte-for-byte mirror of ElmDiagnostics.F90 rather
  ! than a new format: everything that reads ELM's binary then reads this one
  ! unchanged.
  !
  ! Records are written in kELM's PACKED index space (0..n_kokkos_col-1), not
  ! ELM's global column space.  The maps needed to align the two are written
  ! once per run under 'kelmmap:' so the comparison tool can do it.
  !
  ! Snapshot point: the top of the timestep, before any kernel has run.  ELM's
  ! equivalent anchor is 'canhydro_in:' — CanopyHydrology is its first kernel —
  ! so 'kelm_in:<var>' and 'canhydro_in:<var>' describe the same instant and
  ! are directly comparable.  That is what makes this an error-growth trace
  ! rather than another single-step replay.
  !-----------------------------------------------------------------------

  use shr_kind_mod, only : r8 => shr_kind_r8
  use shr_sys_mod , only : shr_sys_abort
  use kelm_mod   , only : kELMType, KELM_SUCCESS,                        &
                           kELMGetTSoisno, kELMGetH2osoiLiqSoi,           &
                           kELMGetTSoisnoSno, kELMGetH2osoiLiqSno,        &
                           kELMGetH2osoiIceSno, kELMGetDzSno,             &
                           kELMGetSnwRds,                                  &
                           kELMGetDz, kELMGetH2osoiIce,                   &
                           kELMGetH2osoiLiq, kELMGetSnl,                  &
                           kELMGetFracSnoEff,                              &
                           kELMGetHsSoil, kELMGetHsTopSnow,               &
                           kELMGetHsH2osfc, kELMGetDhsdT,                 &
                           kELMGetEflxShSnow, kELMGetQflxEvSnow,          &
                           kELMGetDlrad, kELMGetEflxShSoil,                &
                           kELMGetQflxEvSoil, kELMGetTVeg,                 &
                           kELMGetQgSnow, kELMGetQgSoil, kELMGetThm,      &
                           kELMGetQg, kELMGetHtvp, kELMGetSoilbeta,       &
                           kELMGetZ0mg, kELMGetThv,                        &
                           kELMGetQflxEvapGrndCol, kELMGetQflxEvapGrnd,    &
                           kELMGetTSsbef, kELMGetCgrndl, kELMGetCgrnds,   &
                           kELMGetQflxTopSoilCol, kELMGetFracH2osfc,       &
                           kELMGetForcRhoCol, kELMGetDqgdT, kELMGetZii,   &
                           kELMGetSweOld, kELMGetImeltReal,                &
                           kELMGetQflxDewSnowCol,                          &
                           kELMGetSabgSoil, kELMGetSabgSnow,               &
                           kELMGetSabg, kELMGetSabgPen,                    &
                           kELMGetSabgLyr,                                  &
                           kELMGetMssBcphi, kELMGetMssBcpho,               &
                           kELMGetMssOcphi, kELMGetMssOcpho,               &
                           kELMGetMssDst1, kELMGetMssDst2,                 &
                           kELMGetMssDst3, kELMGetMssDst4,                 &
                           kELMGetQflxSnowGrnd, kELMGetQflxRainGrnd,       &
                           kELMGetQflxPrecIntr, kELMGetQflxEvapVeg,        &
                           kELMGetH2osoiIceSoi, kELMGetTGrnd,             &
                           kELMGetTH2osfc, kELMGetH2osfc,                 &
                           kELMGetH2osno, kELMGetSnowDepth,               &
                           kELMGetFracSno, kELMGetFracH2osfc,             &
                           kELMGetIntSnow, kELMGetSnl,                    &
                           kELMGetQflxSnowMelt,                            &
                           kELMGetTVeg, kELMGetBtran, kELMGetH2ocan,     &
                           kELMGetEflxShTot, kELMGetEflxLhTot,            &
                           kELMGetQflxEvapTot, kELMGetFsa, kELMGetFsr,   &
                           kELMGetTRef2m, kELMGetQflxTranVeg,             &
                           kELMGetQflxInflCol, kELMGetQflxSurfCol,        &
                           kELMGetQflxDrainCol, kELMGetQflxEvapSoi,       &
                           kELMGetQflxTopSoilCol, kELMGetFsat,            &
                           kELMGetFcov, kELMGetZwt, kELMGetWtfact,       &
                           kELMGetEffPorosity, kELMGetAlbgrd,             &
                           kELMGetAlbgri, kELMGetAlbsod, kELMGetAlbd,     &
                           kELMGetElai
  use kelmKokkosStateMod, only : n_kokkos_col, n_kokkos_patch,             &
                                  col_of_kcol, patch_of_kpatch,            &
                                  kcol_of_col
  use kelmSoilPropMod   , only : watsat, bsw, sucsat, hksat, nlevgrnd, nlevsno
  use kelmRootMod       , only : rootr, btran_root => btran
  use kelmSurfaceStateMod, only : patch_lai, patch_sai
  implicit none
  private

  logical, public :: kelm_diag_enabled  = .false.
  integer, public :: kelm_diag_timestep = 0

  integer, parameter :: LABEL_LEN = 64
  integer            :: dunit     = 211
  logical            :: is_open   = .false.

  public :: kelm_diag_init
  public :: kelm_diag_finalize
  public :: kelm_diag_new_timestep
  public :: kelm_diag_1d
  public :: kelm_diag_2d
  public :: kelm_diag_int_1d
  public :: kelm_diag_snapshot_state
  public :: kelm_diag_snapshot_fluxes
  public :: kelm_diag_snapshot_preinfil
  public :: kelm_diag_snapshot_presoilflux
  public :: kelm_diag_snapshot_presoiltemp
  public :: kelm_diag_snapshot_soilwater
  public :: kelm_diag_snapshot_snowstate
  public :: kelm_diag_snapshot_postcanhydro
  public :: kelm_diag_snapshot_postcanflux
  public :: kelm_diag_write_maps
  public :: kelm_diag_flush

contains

  subroutine kelm_diag_init(output_path, enabled)
    character(len=*), intent(in) :: output_path
    logical,          intent(in) :: enabled

    kelm_diag_enabled  = enabled
    kelm_diag_timestep = 0
    if (.not. kelm_diag_enabled) return

    open(unit=dunit, file=trim(output_path), &
         form='unformatted', access='stream', action='write', status='replace')
    is_open = .true.

    write(dunit) 'ELMDIAG1'
    write(dunit) LABEL_LEN
  end subroutine kelm_diag_init

  subroutine kelm_diag_finalize()
    if (is_open) then
       close(dunit)
       is_open = .false.
    end if
    kelm_diag_enabled = .false.
  end subroutine kelm_diag_finalize

  subroutine kelm_diag_flush()
    ! Push buffered records to disk. Called after the init-time dumps, so a
    ! run that aborts at its first step still leaves them readable.
    if (is_open) flush(dunit)
  end subroutine kelm_diag_flush

  subroutine kelm_diag_new_timestep(nstep)
    ! Label records by the driver's nstep, matching what ElmDiagnostics now
    ! does, so kELM step N and ELM step N are the same instant.
    integer, intent(in), optional :: nstep
    if (.not. kelm_diag_enabled) return
    if (present(nstep)) then
       kelm_diag_timestep = nstep
    else
       kelm_diag_timestep = kelm_diag_timestep + 1
    end if
  end subroutine kelm_diag_new_timestep

  subroutine kelm_diag_1d(label, array, n)
    character(len=*), intent(in) :: label
    real(r8),         intent(in) :: array(:)
    integer,          intent(in) :: n
    character(len=LABEL_LEN) :: padded
    if (.not. kelm_diag_enabled) return
    if (n <= 0) return
    padded = label
    write(dunit) kelm_diag_timestep
    write(dunit) padded
    write(dunit) 1
    write(dunit) n
    write(dunit) array(1:n)
  end subroutine kelm_diag_1d

  subroutine kelm_diag_2d(label, array, n1, n2)
    character(len=*), intent(in) :: label
    real(r8),         intent(in) :: array(:,:)
    integer,          intent(in) :: n1, n2
    character(len=LABEL_LEN) :: padded
    if (.not. kelm_diag_enabled) return
    if (n1 <= 0 .or. n2 <= 0) return
    padded = label
    write(dunit) kelm_diag_timestep
    write(dunit) padded
    write(dunit) 2
    write(dunit) n1, n2
    write(dunit) array(1:n1, 1:n2)
  end subroutine kelm_diag_2d

  subroutine kelm_diag_int_1d(label, array, n)
    character(len=*), intent(in) :: label
    integer,          intent(in) :: array(:)
    integer,          intent(in) :: n
    character(len=LABEL_LEN) :: padded
    if (.not. kelm_diag_enabled) return
    if (n <= 0) return
    padded = label
    write(dunit) kelm_diag_timestep
    write(dunit) padded
    write(dunit) -1
    write(dunit) n
    write(dunit) array(1:n)
  end subroutine kelm_diag_int_1d

  !-----------------------------------------------------------------------
  ! Write the packed -> ELM index maps once, so the comparison tool can put
  ! kELM's columns and patches back where ELM has them.
  !-----------------------------------------------------------------------
  subroutine kelm_diag_write_maps()
    integer :: kc, c, j
    real(r8), allocatable :: tmp(:,:)
    if (.not. kelm_diag_enabled) return
    if (associated(col_of_kcol)) &
         call kelm_diag_int_1d('kelmmap:col_of_kcol', col_of_kcol, n_kokkos_col)
    if (associated(patch_of_kpatch)) &
         call kelm_diag_int_1d('kelmmap:patch_of_kpatch', patch_of_kpatch, n_kokkos_patch)

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
       call kelm_diag_2d('kelmsoil:watsat', tmp, n_kokkos_col, nlevgrnd)
       do kc = 1, n_kokkos_col
          c = col_of_kcol(kc)
          tmp(kc,1:nlevgrnd) = bsw(c,1:nlevgrnd)
       end do
       call kelm_diag_2d('kelmsoil:bsw', tmp, n_kokkos_col, nlevgrnd)
       do kc = 1, n_kokkos_col
          c = col_of_kcol(kc)
          tmp(kc,1:nlevgrnd) = sucsat(c,1:nlevgrnd)
       end do
       call kelm_diag_2d('kelmsoil:sucsat', tmp, n_kokkos_col, nlevgrnd)
       do kc = 1, n_kokkos_col
          c = col_of_kcol(kc)
          tmp(kc,1:nlevgrnd) = hksat(c,1:nlevgrnd)
       end do
       call kelm_diag_2d('kelmsoil:hksat', tmp, n_kokkos_col, nlevgrnd)
       deallocate(tmp)
    end if
  end subroutine kelm_diag_write_maps

  !-----------------------------------------------------------------------
  ! Snapshot the carried state at the top of a timestep, before any kernel
  ! has run.  ELM's matching anchor is 'canhydro_in:' — its first kernel — so
  ! these records line up instant-for-instant with ELM's own.
  !-----------------------------------------------------------------------
  subroutine kelm_diag_snapshot_state(elm, nlevtot, nlevgrnd, tag)
    type(kELMType) , intent(in) :: elm
    integer         , intent(in) :: nlevtot, nlevgrnd
    character(len=*), intent(in) :: tag

    integer  :: ierr, sz(2), szg(2), szs(2)
    real(r8), allocatable :: c1(:), p1(:), c2(:,:), cg(:,:), cs(:,:)
    integer , allocatable :: ci(:)

    if (.not. kelm_diag_enabled) return
    if (n_kokkos_col <= 0) return

    allocate(c1(n_kokkos_col), ci(n_kokkos_col))
    allocate(c2(n_kokkos_col, nlevtot))
    allocate(cg(n_kokkos_col, nlevgrnd))
    allocate(cs(n_kokkos_col, nlevsno))
    allocate(p1(max(n_kokkos_patch,1)))

    ! ---- column scalars ----
    call get_c('t_grnd',      kELMGetTGrnd)
    call get_c('t_h2osfc',    kELMGetTH2osfc)
    call get_c('h2osfc',      kELMGetH2osfc)
    call get_c('h2osno',      kELMGetH2osno)
    call get_c('snow_depth',  kELMGetSnowDepth)
    call get_c('frac_sno',    kELMGetFracSno)
    call get_c('frac_h2osfc', kELMGetFracH2osfc)
    call get_c('int_snow',    kELMGetIntSnow)

    call kELMGetSnl(elm, ci, n_kokkos_col, ierr)
    if (ierr == KELM_SUCCESS) call kelm_diag_int_1d(tag//':snl', ci, n_kokkos_col)

    ! ---- column profiles ----
    sz(1) = n_kokkos_col; sz(2) = nlevtot
    call get_c2('t_soisno',   kELMGetTSoisno)
    szg(1) = n_kokkos_col; szg(2) = nlevgrnd
    call get_cg('h2osoi_liq_soi', kELMGetH2osoiLiqSoi)
    call get_cg('h2osoi_ice_soi', kELMGetH2osoiIceSoi)

    ! The snow half. This used to be omitted, on the reasoning that
    ! "1x1_brazil never carries snow, so nothing is lost" -- which stopped
    ! being true the moment 1x1_glc was used, and left the snowpack
    ! unobservable exactly where the divergence lives. Reported in ELM order
    ! so it lines up with the restart file and elm_diagnostics.bin.
    szs(1) = n_kokkos_col; szs(2) = nlevsno
    call get_cs('t_soisno_sno',   kELMGetTSoisnoSno)
    call get_cs('h2osoi_liq_sno', kELMGetH2osoiLiqSno)
    call get_cs('h2osoi_ice_sno', kELMGetH2osoiIceSno)
    call get_cs('dz_sno',         kELMGetDzSno)
    call get_cs('snw_rds',        kELMGetSnwRds)
    ! Snow aerosol mass per layer [kg/m2], ELM order. ELM records the same at
    ! snowlayer_in:mss_*. Traced so melt-season scavenging (S3) can be graded
    ! against ELM rather than by construction -- the suspect for April's snow
    ! being brighter than ELM's.
    call get_cs('mss_bcphi', kELMGetMssBcphi)
    call get_cs('mss_bcpho', kELMGetMssBcpho)
    call get_cs('mss_ocphi', kELMGetMssOcphi)
    call get_cs('mss_ocpho', kELMGetMssOcpho)
    call get_cs('mss_dst1',  kELMGetMssDst1)
    call get_cs('mss_dst2',  kELMGetMssDst2)
    call get_cs('mss_dst3',  kELMGetMssDst3)
    call get_cs('mss_dst4',  kELMGetMssDst4)

    ! ---- patch scalars ----
    if (n_kokkos_patch > 0) then
       call get_p('t_veg',         kELMGetTVeg)
       call get_p('btran',         kELMGetBtran)
       call get_p('h2ocan',        kELMGetH2ocan)
       call get_p('t_ref2m',       kELMGetTRef2m)
       call get_p('eflx_sh_tot',   kELMGetEflxShTot)
       call get_p('eflx_lh_tot',   kELMGetEflxLhTot)
       call get_p('qflx_evap_tot', kELMGetQflxEvapTot)
       ! Water budget terms. Recorded so a moisture drift can be attributed to
       ! a flux rather than inferred from the storage profile.
       call get_p('qflx_tran_veg', kELMGetQflxTranVeg)
       call get_p('qflx_evap_soi', kELMGetQflxEvapSoi)
    end if

    deallocate(c1, ci, c2, cg, p1)

  contains

    subroutine get_c(name, getter)
      character(len=*), intent(in) :: name
      external :: getter
      call getter(elm, c1, n_kokkos_col, ierr)
      if (ierr == KELM_SUCCESS) call kelm_diag_1d(tag//':'//name, c1, n_kokkos_col)
    end subroutine get_c

    subroutine get_p(name, getter)
      character(len=*), intent(in) :: name
      external :: getter
      call getter(elm, p1, n_kokkos_patch, ierr)
      if (ierr == KELM_SUCCESS) call kelm_diag_1d(tag//':'//name, p1, n_kokkos_patch)
    end subroutine get_p

    subroutine get_cg(name, getter)
      character(len=*), intent(in) :: name
      external :: getter
      call getter(elm, cg, szg, ierr)
      if (ierr == KELM_SUCCESS) &
           call kelm_diag_2d(tag//':'//name, cg, n_kokkos_col, nlevgrnd)
    end subroutine get_cg

    subroutine get_cs(name, getter)
      character(len=*), intent(in) :: name
      external :: getter
      call getter(elm, cs, szs, ierr)
      if (ierr == KELM_SUCCESS) &
           call kelm_diag_2d(tag//':'//name, cs, n_kokkos_col, nlevsno)
    end subroutine get_cs

    subroutine get_c2(name, getter)
      character(len=*), intent(in) :: name
      external :: getter
      call getter(elm, c2, sz, ierr)
      if (ierr == KELM_SUCCESS) &
           call kelm_diag_2d(tag//':'//name, c2, n_kokkos_col, nlevtot)
    end subroutine get_c2

  end subroutine kelm_diag_snapshot_state

  !-----------------------------------------------------------------------
  ! Water fluxes, sampled at the END of the step. The HydrologyDrainage
  ! getters read elm->hydrologyDrainage, which is a separate structure from
  ! the natural column and holds nothing at the top of a step -- sampling
  ! there returns zeros and reads as "no infiltration ever happened".
  !-----------------------------------------------------------------------
  subroutine kelm_diag_snapshot_fluxes(elm, tag)
    type(kELMType) , intent(in) :: elm
    character(len=*), intent(in) :: tag
    integer :: ierr
    real(r8), allocatable :: c1(:)

    if (.not. kelm_diag_enabled) return
    if (n_kokkos_col <= 0) return
    allocate(c1(n_kokkos_col))

    call kELMGetQflxInflCol(elm, c1, n_kokkos_col, ierr)
    if (ierr == KELM_SUCCESS) call kelm_diag_1d(tag//':qflx_infl', c1, n_kokkos_col)
    call kELMGetQflxSurfCol(elm, c1, n_kokkos_col, ierr)
    if (ierr == KELM_SUCCESS) call kelm_diag_1d(tag//':qflx_surf', c1, n_kokkos_col)
    call kELMGetQflxDrainCol(elm, c1, n_kokkos_col, ierr)
    if (ierr == KELM_SUCCESS) call kelm_diag_1d(tag//':qflx_drain', c1, n_kokkos_col)

    ! Absorbed/reflected shortwave, sampled here rather than at the top of the
    ! step: ELM records surfrad_out: DURING the step, so a top-of-step kELM
    ! value is one timestep stale and reads as a whole-diurnal-cycle shift.
    block
      real(r8), allocatable :: p1b(:)
      if (n_kokkos_patch > 0) then
         allocate(p1b(n_kokkos_patch))
         call kELMGetFsa(elm, p1b, n_kokkos_patch, ierr)
         if (ierr == KELM_SUCCESS) call kelm_diag_1d(tag//':fsa', p1b, n_kokkos_patch)
         call kELMGetFsr(elm, p1b, n_kokkos_patch, ierr)
         if (ierr == KELM_SUCCESS) call kelm_diag_1d(tag//':fsr', p1b, n_kokkos_patch)
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
         call kELMGetElai(elm, pl, n_kokkos_patch, ierr)
         if (ierr == KELM_SUCCESS) call kelm_diag_1d(tag//':elai', pl, n_kokkos_patch)
         do kp = 1, n_kokkos_patch
            pl(kp) = patch_sai(patch_of_kpatch(kp))
         end do
         call kelm_diag_1d(tag//':esai', pl, n_kokkos_patch)
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
         call kelm_diag_2d(tag//':rootr', pg, n_kokkos_patch, nlevgrnd)
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
      call kELMGetAlbgrd(elm, c2, szc, ierr)
      if (ierr == KELM_SUCCESS) call kelm_diag_2d(tag//':albgrd', c2, n_kokkos_col, 2)
      call kELMGetAlbgri(elm, c2, szc, ierr)
      if (ierr == KELM_SUCCESS) call kelm_diag_2d(tag//':albgri', c2, n_kokkos_col, 2)
      call kELMGetAlbsod(elm, c2, szc, ierr)
      if (ierr == KELM_SUCCESS) call kelm_diag_2d(tag//':albsod', c2, n_kokkos_col, 2)
      if (n_kokkos_patch > 0) then
         szp(1) = n_kokkos_patch; szp(2) = 2
         call kELMGetAlbd(elm, p2, szp, ierr)
         if (ierr == KELM_SUCCESS) call kelm_diag_2d(tag//':albd', p2, n_kokkos_patch, 2)
      end if
      deallocate(c2, p2)
    end block

    ! The quantities that set the infiltration/runoff split.
    call kELMGetQflxTopSoilCol(elm, c1, n_kokkos_col, ierr)
    if (ierr == KELM_SUCCESS) call kelm_diag_1d(tag//':qflx_top_soil', c1, n_kokkos_col)
    call kELMGetFsat(elm, c1, n_kokkos_col, ierr)
    if (ierr == KELM_SUCCESS) call kelm_diag_1d(tag//':fsat', c1, n_kokkos_col)
    call kELMGetFcov(elm, c1, n_kokkos_col, ierr)
    if (ierr == KELM_SUCCESS) call kelm_diag_1d(tag//':fcov', c1, n_kokkos_col)
    call kELMGetZwt(elm, c1, n_kokkos_col, ierr)
    if (ierr == KELM_SUCCESS) call kelm_diag_1d(tag//':zwt', c1, n_kokkos_col)
    call kELMGetWtfact(elm, c1, n_kokkos_col, ierr)
    if (ierr == KELM_SUCCESS) call kelm_diag_1d(tag//':wtfact', c1, n_kokkos_col)

    block
      real(r8), allocatable :: cg(:,:)
      integer :: szg(2)
      allocate(cg(n_kokkos_col, nlevgrnd))
      szg(1) = n_kokkos_col; szg(2) = nlevgrnd
      call kELMGetEffPorosity(elm, cg, szg, ierr)
      if (ierr == KELM_SUCCESS) &
           call kelm_diag_2d(tag//':eff_porosity', cg, n_kokkos_col, nlevgrnd)
      deallocate(cg)
    end block

    deallocate(c1)
  end subroutine kelm_diag_snapshot_fluxes

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
  subroutine kelm_diag_snapshot_preinfil(elm, tag)
    type(kELMType) , intent(in) :: elm
    character(len=*), intent(in) :: tag
    integer :: ierr
    real(r8), allocatable :: c1(:)
    if (.not. kelm_diag_enabled) return
    if (n_kokkos_col <= 0) return
    allocate(c1(n_kokkos_col))
    call kELMGetQflxEvapGrndCol(elm, c1, n_kokkos_col, ierr)
    if (ierr == KELM_SUCCESS) call kelm_diag_1d(tag//':qflx_evap_grnd_col', c1, n_kokkos_col)
    call kELMGetForcRhoCol(elm, c1, n_kokkos_col, ierr)
    if (ierr == KELM_SUCCESS) call kelm_diag_1d(tag//':forc_rho', c1, n_kokkos_col)
    call kELMGetDqgdT(elm, c1, n_kokkos_col, ierr)
    if (ierr == KELM_SUCCESS) call kelm_diag_1d(tag//':dqgdT', c1, n_kokkos_col)
    call kELMGetZii(elm, c1, n_kokkos_col, ierr)
    if (ierr == KELM_SUCCESS) call kelm_diag_1d(tag//':zii', c1, n_kokkos_col)
    call kELMGetQflxTopSoilCol(elm, c1, n_kokkos_col, ierr)
    if (ierr == KELM_SUCCESS) call kelm_diag_1d(tag//':qflx_top_soil', c1, n_kokkos_col)
    call kELMGetFracH2osfc(elm, c1, n_kokkos_col, ierr)
    if (ierr == KELM_SUCCESS) call kelm_diag_1d(tag//':frac_h2osfc', c1, n_kokkos_col)
    call kELMGetFracSnoEff(elm, c1, n_kokkos_col, ierr)
    if (ierr == KELM_SUCCESS) call kelm_diag_1d(tag//':frac_sno_eff', c1, n_kokkos_col)
    call kELMGetQflxDewSnowCol(elm, c1, n_kokkos_col, ierr)
    if (ierr == KELM_SUCCESS) call kelm_diag_1d(tag//':qflx_dew_snow_col', c1, n_kokkos_col)
    block
      real(r8), allocatable :: cg(:,:)
      integer :: szt(2)
      allocate(cg(n_kokkos_col, nlevsno+nlevgrnd))
      szt(1) = n_kokkos_col; szt(2) = nlevsno + nlevgrnd
      call kELMGetSweOld(elm, cg, szt, ierr)
      if (ierr == KELM_SUCCESS) call kelm_diag_2d(tag//':swe_old', cg, n_kokkos_col, nlevsno+nlevgrnd)
      call kELMGetImeltReal(elm, cg, szt, ierr)
      if (ierr == KELM_SUCCESS) call kelm_diag_2d(tag//':imelt', cg, n_kokkos_col, nlevsno+nlevgrnd)
      deallocate(cg)
    end block
    block
      real(r8), allocatable :: p1(:)
      allocate(p1(n_kokkos_patch))
      call kELMGetQflxEvapGrnd(elm, p1, n_kokkos_patch, ierr)
      if (ierr == KELM_SUCCESS) call kelm_diag_1d(tag//':qflx_evap_grnd', p1, n_kokkos_patch)
      deallocate(p1)
    end block
    deallocate(c1)
  end subroutine kelm_diag_snapshot_preinfil

  subroutine kelm_diag_snapshot_presoilflux(elm, nlevtot_in, tag)
    type(kELMType) , intent(in) :: elm
    integer         , intent(in) :: nlevtot_in
    character(len=*), intent(in) :: tag
    integer :: ierr, szg(2)
    real(r8), allocatable :: cg(:,:)
    real(r8), allocatable :: c1(:), p1(:)
    if (.not. kelm_diag_enabled) return
    if (n_kokkos_col <= 0) return
    allocate(cg(n_kokkos_col, nlevtot_in), c1(n_kokkos_col), p1(n_kokkos_patch))
    szg(1) = n_kokkos_col; szg(2) = nlevtot_in
    call kELMGetTSsbef(elm, cg, szg, ierr)
    if (ierr == KELM_SUCCESS) call kelm_diag_2d(tag//':tssbef', cg, n_kokkos_col, nlevtot_in)
    call kELMGetTSoisno(elm, cg, szg, ierr)
    if (ierr == KELM_SUCCESS) call kelm_diag_2d(tag//':t_soisno', cg, n_kokkos_col, nlevtot_in)
    call kELMGetTGrnd(elm, c1, n_kokkos_col, ierr)
    if (ierr == KELM_SUCCESS) call kelm_diag_1d(tag//':t_grnd', c1, n_kokkos_col)
    call kELMGetCgrndl(elm, p1, n_kokkos_patch, ierr)
    if (ierr == KELM_SUCCESS) call kelm_diag_1d(tag//':cgrndl', p1, n_kokkos_patch)
    call kELMGetCgrnds(elm, p1, n_kokkos_patch, ierr)
    if (ierr == KELM_SUCCESS) call kelm_diag_1d(tag//':cgrnds', p1, n_kokkos_patch)
    deallocate(cg, c1, p1)
  end subroutine kelm_diag_snapshot_presoilflux

  subroutine kelm_diag_snapshot_presoiltemp(elm, nlevtot_in, tag)
    type(kELMType) , intent(in) :: elm
    integer         , intent(in) :: nlevtot_in
    character(len=*), intent(in) :: tag
    integer :: ierr, szg(2)
    real(r8), allocatable :: cg(:,:)
    real(r8), allocatable :: c1(:)
    integer , allocatable :: i1(:)
    if (.not. kelm_diag_enabled) return
    if (n_kokkos_col <= 0) return
    allocate(cg(n_kokkos_col, nlevtot_in), c1(n_kokkos_col), i1(n_kokkos_col))
    szg(1) = n_kokkos_col; szg(2) = nlevtot_in
    call kELMGetTSoisno(elm, cg, szg, ierr)
    if (ierr == KELM_SUCCESS) call kelm_diag_2d(tag//':t_soisno', cg, n_kokkos_col, nlevtot_in)
    call kELMGetDz(elm, cg, szg, ierr)
    if (ierr == KELM_SUCCESS) call kelm_diag_2d(tag//':dz', cg, n_kokkos_col, nlevtot_in)
    call kELMGetH2osoiIce(elm, cg, szg, ierr)
    if (ierr == KELM_SUCCESS) call kelm_diag_2d(tag//':h2osoi_ice', cg, n_kokkos_col, nlevtot_in)
    call kELMGetH2osoiLiq(elm, cg, szg, ierr)
    if (ierr == KELM_SUCCESS) call kelm_diag_2d(tag//':h2osoi_liq', cg, n_kokkos_col, nlevtot_in)
    call kELMGetSnl(elm, i1, n_kokkos_col, ierr)
    if (ierr == KELM_SUCCESS) call kelm_diag_int_1d(tag//':snl', i1, n_kokkos_col)
    call kELMGetFracSnoEff(elm, c1, n_kokkos_col, ierr)
    if (ierr == KELM_SUCCESS) call kelm_diag_1d(tag//':frac_sno_eff', c1, n_kokkos_col)
    call kELMGetHsSoil(elm, c1, n_kokkos_col, ierr)
    if (ierr == KELM_SUCCESS) call kelm_diag_1d(tag//':hs_soil', c1, n_kokkos_col)
    call kELMGetHsTopSnow(elm, c1, n_kokkos_col, ierr)
    if (ierr == KELM_SUCCESS) call kelm_diag_1d(tag//':hs_top_snow', c1, n_kokkos_col)
    call kELMGetHsH2osfc(elm, c1, n_kokkos_col, ierr)
    if (ierr == KELM_SUCCESS) call kelm_diag_1d(tag//':hs_h2osfc', c1, n_kokkos_col)
    call kELMGetDhsdT(elm, c1, n_kokkos_col, ierr)
    if (ierr == KELM_SUCCESS) call kelm_diag_1d(tag//':dhsdT', c1, n_kokkos_col)
    call kELMGetQgSnow(elm, c1, n_kokkos_col, ierr)
    if (ierr == KELM_SUCCESS) call kelm_diag_1d(tag//':qg_snow', c1, n_kokkos_col)
    call kELMGetQgSoil(elm, c1, n_kokkos_col, ierr)
    if (ierr == KELM_SUCCESS) call kelm_diag_1d(tag//':qg_soil', c1, n_kokkos_col)
    call kELMGetTGrnd(elm, c1, n_kokkos_col, ierr)
    if (ierr == KELM_SUCCESS) call kelm_diag_1d(tag//':t_grnd', c1, n_kokkos_col)
    call kELMGetQg(elm, c1, n_kokkos_col, ierr)
    if (ierr == KELM_SUCCESS) call kelm_diag_1d(tag//':qg', c1, n_kokkos_col)
    call kELMGetHtvp(elm, c1, n_kokkos_col, ierr)
    if (ierr == KELM_SUCCESS) call kelm_diag_1d(tag//':htvp', c1, n_kokkos_col)
    call kELMGetSoilbeta(elm, c1, n_kokkos_col, ierr)
    if (ierr == KELM_SUCCESS) call kelm_diag_1d(tag//':soilbeta', c1, n_kokkos_col)
    call kELMGetZ0mg(elm, c1, n_kokkos_col, ierr)
    if (ierr == KELM_SUCCESS) call kelm_diag_1d(tag//':z0mg', c1, n_kokkos_col)
    call kELMGetThv(elm, c1, n_kokkos_col, ierr)
    if (ierr == KELM_SUCCESS) call kelm_diag_1d(tag//':thv', c1, n_kokkos_col)
    call kELMGetQflxEvapGrndCol(elm, c1, n_kokkos_col, ierr)
    if (ierr == KELM_SUCCESS) call kelm_diag_1d(tag//':qflx_evap_grnd_col', c1, n_kokkos_col)
    block
      real(r8), allocatable :: p1(:)
      allocate(p1(n_kokkos_patch))
      call kELMGetEflxShSnow(elm, p1, n_kokkos_patch, ierr)
      if (ierr == KELM_SUCCESS) call kelm_diag_1d(tag//':eflx_sh_snow', p1, n_kokkos_patch)
      call kELMGetQflxEvSnow(elm, p1, n_kokkos_patch, ierr)
      if (ierr == KELM_SUCCESS) call kelm_diag_1d(tag//':qflx_ev_snow', p1, n_kokkos_patch)
      call kELMGetDlrad(elm, p1, n_kokkos_patch, ierr)
      if (ierr == KELM_SUCCESS) call kelm_diag_1d(tag//':dlrad', p1, n_kokkos_patch)
      call kELMGetEflxShSoil(elm, p1, n_kokkos_patch, ierr)
      if (ierr == KELM_SUCCESS) call kelm_diag_1d(tag//':eflx_sh_soil', p1, n_kokkos_patch)
      call kELMGetQflxEvSoil(elm, p1, n_kokkos_patch, ierr)
      if (ierr == KELM_SUCCESS) call kelm_diag_1d(tag//':qflx_ev_soil', p1, n_kokkos_patch)
      call kELMGetTVeg(elm, p1, n_kokkos_patch, ierr)
      if (ierr == KELM_SUCCESS) call kelm_diag_1d(tag//':t_veg', p1, n_kokkos_patch)
      call kELMGetThm(elm, p1, n_kokkos_patch, ierr)
      if (ierr == KELM_SUCCESS) call kelm_diag_1d(tag//':thm', p1, n_kokkos_patch)
      call kELMGetQflxEvapGrnd(elm, p1, n_kokkos_patch, ierr)
      if (ierr == KELM_SUCCESS) call kelm_diag_1d(tag//':qflx_evap_grnd', p1, n_kokkos_patch)
      call kELMGetSabgSoil(elm, p1, n_kokkos_patch, ierr)
      if (ierr == KELM_SUCCESS) call kelm_diag_1d(tag//':sabg_soil', p1, n_kokkos_patch)
      call kELMGetSabgSnow(elm, p1, n_kokkos_patch, ierr)
      if (ierr == KELM_SUCCESS) call kelm_diag_1d(tag//':sabg_snow', p1, n_kokkos_patch)
      ! sabg and sabg_pen: ELM records both at surfrad_out, and sabg_pen is a
      ! DIFFERENCE of two large numbers (sabg - lyr[top]), so it cannot be
      ! reconstructed from the other traces without error. Added 2026-09-11
      ! while chasing April's +3.72 W/m2.
      call kELMGetSabg(elm, p1, n_kokkos_patch, ierr)
      if (ierr == KELM_SUCCESS) call kelm_diag_1d(tag//':sabg', p1, n_kokkos_patch)
      call kELMGetSabgPen(elm, p1, n_kokkos_patch, ierr)
      if (ierr == KELM_SUCCESS) call kelm_diag_1d(tag//':sabg_pen', p1, n_kokkos_patch)
      deallocate(p1)
    end block
    block
      real(r8), allocatable :: plyr(:,:)
      integer :: szp(2)
      allocate(plyr(n_kokkos_patch, nlevsno+1))
      szp(1) = n_kokkos_patch; szp(2) = nlevsno + 1
      ! ELM order (top snow first, soil surface last), like every snow field
      ! in this trace: the getter reflects the H7 storage. Until 2026-09-27 it
      ! did not, and this one label was H7 among ELM-ordered neighbours (F3).
      call kELMGetSabgLyr(elm, plyr, szp, ierr)
      if (ierr == KELM_SUCCESS) call kelm_diag_2d(tag//':sabg_lyr', plyr, n_kokkos_patch, nlevsno+1)
      deallocate(plyr)
    end block
    ! Snow aerosol mass (mss_*) is traced with the rest of the snow state, in
    ! the per-step snapshot above.
    deallocate(cg, c1, i1)
  end subroutine kelm_diag_snapshot_presoiltemp

  subroutine kelm_diag_snapshot_soilwater(elm, nlevgrnd_in, tag)
    type(kELMType) , intent(in) :: elm
    integer         , intent(in) :: nlevgrnd_in
    character(len=*), intent(in) :: tag
    integer :: ierr, szg(2)
    real(r8), allocatable :: cg(:,:)
    if (.not. kelm_diag_enabled) return
    if (n_kokkos_col <= 0) return
    allocate(cg(n_kokkos_col, nlevgrnd_in))
    szg(1) = n_kokkos_col; szg(2) = nlevgrnd_in
    call kELMGetH2osoiLiqSoi(elm, cg, szg, ierr)
    if (ierr == KELM_SUCCESS) &
         call kelm_diag_2d(tag//':h2osoi_liq', cg, n_kokkos_col, nlevgrnd_in)
    deallocate(cg)
  end subroutine kelm_diag_snapshot_soilwater

  subroutine kelm_diag_snapshot_snowstate(elm, tag)
    ! The snow-cover state CanopyHydrology both reads and writes, dumped so it
    ! can be differenced step by step against ELM's own canhydro_in /
    ! canhydro_out records. frac_sno is the field of interest: ELM never
    ! modifies it between two CanopyHydrology calls (verified, 0 of 1488 steps
    ! of the 1x1_glc binary), so its whole trajectory is that one routine's,
    ! and a divergence has to come from one of the inputs dumped alongside.
    type(kELMType) , intent(in) :: elm
    character(len=*), intent(in) :: tag
    integer :: ierr
    real(r8), allocatable :: c1(:)
    integer,  allocatable :: i1(:)
    if (.not. kelm_diag_enabled) return
    if (n_kokkos_col <= 0) return
    allocate(c1(n_kokkos_col), i1(n_kokkos_col))
    call kELMGetFracSno(elm, c1, n_kokkos_col, ierr)
    if (ierr == KELM_SUCCESS) call kelm_diag_1d(tag//':frac_sno', c1, n_kokkos_col)
    call kELMGetFracSnoEff(elm, c1, n_kokkos_col, ierr)
    if (ierr == KELM_SUCCESS) call kelm_diag_1d(tag//':frac_sno_eff', c1, n_kokkos_col)
    call kELMGetH2osno(elm, c1, n_kokkos_col, ierr)
    if (ierr == KELM_SUCCESS) call kelm_diag_1d(tag//':h2osno', c1, n_kokkos_col)
    call kELMGetIntSnow(elm, c1, n_kokkos_col, ierr)
    if (ierr == KELM_SUCCESS) call kelm_diag_1d(tag//':int_snow', c1, n_kokkos_col)
    call kELMGetSnowDepth(elm, c1, n_kokkos_col, ierr)
    if (ierr == KELM_SUCCESS) call kelm_diag_1d(tag//':snow_depth', c1, n_kokkos_col)
    call kELMGetQflxSnowMelt(elm, c1, n_kokkos_col, ierr)
    if (ierr == KELM_SUCCESS) call kelm_diag_1d(tag//':qflx_snow_melt', c1, n_kokkos_col)
    call kELMGetSnl(elm, i1, n_kokkos_col, ierr)
    if (ierr == KELM_SUCCESS) then
       c1 = real(i1, r8)
       call kelm_diag_1d(tag//':snl', c1, n_kokkos_col)
    end if
    deallocate(c1, i1)
  end subroutine kelm_diag_snapshot_snowstate

  subroutine kelm_diag_snapshot_postcanhydro(elm, tag)
    ! CanopyHydrology's own outputs, right after it runs. No existing
    ! snapshot point captures this -- the earliest is 'kelm_st', which
    ! runs just before SoilTemperature, several kernels later. Added to
    ! trace whether canopy interception delivers the same snowfall to the
    ! ground as ELM (G11, STATUS.md): a per-event bias here, unlike the
    ! melt-diversion path checked and ruled out already, would fire on
    ! every snowfall step and so could plausibly compound over a month.
    type(kELMType) , intent(in) :: elm
    character(len=*), intent(in) :: tag
    integer :: ierr
    real(r8), allocatable :: c1(:), p1(:)
    if (.not. kelm_diag_enabled) return
    if (n_kokkos_col <= 0) return
    allocate(c1(n_kokkos_col))
    call kELMGetQflxSnowGrnd(elm, c1, n_kokkos_col, ierr)
    if (ierr == KELM_SUCCESS) call kelm_diag_1d(tag//':qflx_snow_grnd', c1, n_kokkos_col)
    call kELMGetQflxRainGrnd(elm, c1, n_kokkos_col, ierr)
    if (ierr == KELM_SUCCESS) call kelm_diag_1d(tag//':qflx_rain_grnd', c1, n_kokkos_col)
    deallocate(c1)
    if (n_kokkos_patch > 0) then
       allocate(p1(n_kokkos_patch))
       call kELMGetQflxPrecIntr(elm, p1, n_kokkos_patch, ierr)
       if (ierr == KELM_SUCCESS) call kelm_diag_1d(tag//':qflx_prec_intr', p1, n_kokkos_patch)
       deallocate(p1)
    end if
  end subroutine kelm_diag_snapshot_postcanhydro

  subroutine kelm_diag_snapshot_postcanflux(elm, tag)
    ! CanopyFluxes' own outputs, right after it runs -- the dew/evaporation
    ! term (h2ocan += (qflx_tran_veg - qflx_evap_veg)*dtime) is where G11's
    ! h2ocan divergence was traced to (STATUS.md): the formula is validated
    ! identical to ELM's, so if h2ocan still diverges here, the cause is
    ! qflx_evap_veg/qflx_tran_veg themselves differing, not this update.
    type(kELMType) , intent(in) :: elm
    character(len=*), intent(in) :: tag
    integer :: ierr
    real(r8), allocatable :: p1(:)
    if (.not. kelm_diag_enabled) return
    if (n_kokkos_patch <= 0) return
    allocate(p1(n_kokkos_patch))
    call kELMGetQflxEvapVeg(elm, p1, n_kokkos_patch, ierr)
    if (ierr == KELM_SUCCESS) call kelm_diag_1d(tag//':qflx_evap_veg', p1, n_kokkos_patch)
    call kELMGetQflxTranVeg(elm, p1, n_kokkos_patch, ierr)
    if (ierr == KELM_SUCCESS) call kelm_diag_1d(tag//':qflx_tran_veg', p1, n_kokkos_patch)
    call kELMGetH2ocan(elm, p1, n_kokkos_patch, ierr)
    if (ierr == KELM_SUCCESS) call kelm_diag_1d(tag//':h2ocan', p1, n_kokkos_patch)
    call kELMGetTVeg(elm, p1, n_kokkos_patch, ierr)
    if (ierr == KELM_SUCCESS) call kelm_diag_1d(tag//':t_veg', p1, n_kokkos_patch)
    deallocate(p1)
  end subroutine kelm_diag_snapshot_postcanflux

end module kelmDiagnosticsMod
