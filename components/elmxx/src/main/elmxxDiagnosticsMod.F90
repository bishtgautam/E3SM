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
                           ELMxxGetH2osoiIceSoi, ELMxxGetTGrnd,             &
                           ELMxxGetTH2osfc, ELMxxGetH2osfc,                 &
                           ELMxxGetH2osno, ELMxxGetSnowDepth,               &
                           ELMxxGetFracSno, ELMxxGetFracH2osfc,             &
                           ELMxxGetIntSnow, ELMxxGetSnl,                    &
                           ELMxxGetTVeg, ELMxxGetBtran, ELMxxGetH2ocan,     &
                           ELMxxGetEflxShTot, ELMxxGetEflxLhTot,            &
                           ELMxxGetQflxEvapTot, ELMxxGetFsa, ELMxxGetFsr,   &
                           ELMxxGetTRef2m
  use elmxxKokkosStateMod, only : n_kokkos_col, n_kokkos_patch,             &
                                  col_of_kcol, patch_of_kpatch,            &
                                  kcol_of_col
  use elmxxSoilPropMod   , only : watsat, bsw, sucsat, hksat, nlevgrnd
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

    integer  :: ierr, sz(2), szg(2)
    real(r8), allocatable :: c1(:), p1(:), c2(:,:), cg(:,:)
    integer , allocatable :: ci(:)

    if (.not. elmxx_diag_enabled) return
    if (n_kokkos_col <= 0) return

    allocate(c1(n_kokkos_col), ci(n_kokkos_col))
    allocate(c2(n_kokkos_col, nlevtot))
    allocate(cg(n_kokkos_col, nlevgrnd))
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
    ! Only the soil-only (NLEVGRND) water getters are exposed to Fortran.
    ! 1x1_brazil never carries snow, so nothing is lost here.
    szg(1) = n_kokkos_col; szg(2) = nlevgrnd
    call get_cg('h2osoi_liq_soi', ELMxxGetH2osoiLiqSoi)
    call get_cg('h2osoi_ice_soi', ELMxxGetH2osoiIceSoi)

    ! ---- patch scalars ----
    if (n_kokkos_patch > 0) then
       call get_p('t_veg',         ELMxxGetTVeg)
       call get_p('btran',         ELMxxGetBtran)
       call get_p('h2ocan',        ELMxxGetH2ocan)
       call get_p('t_ref2m',       ELMxxGetTRef2m)
       call get_p('eflx_sh_tot',   ELMxxGetEflxShTot)
       call get_p('eflx_lh_tot',   ELMxxGetEflxLhTot)
       call get_p('qflx_evap_tot', ELMxxGetQflxEvapTot)
       call get_p('fsa',           ELMxxGetFsa)
       call get_p('fsr',           ELMxxGetFsr)
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

    subroutine get_c2(name, getter)
      character(len=*), intent(in) :: name
      external :: getter
      call getter(elm, c2, sz, ierr)
      if (ierr == ELMXX_SUCCESS) &
           call elmxx_diag_2d(tag//':'//name, c2, n_kokkos_col, nlevtot)
    end subroutine get_c2

  end subroutine elmxx_diag_snapshot_state

end module elmxxDiagnosticsMod
