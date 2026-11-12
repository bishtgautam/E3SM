module kelmLndExportMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Land-to-coupler export: ELM's lnd_export (lnd_import_export.F90).
  !
  ! The library forms the gridcell means (kELMComputeLnd2Atm, ELM's lnd2atm)
  ! on the device, in lnd2atm's sign; one copy brings them to the host
  ! (kelm_kokkos_pull_lnd2atm) and kelm_lnd_export fills l2x with ELM's
  ! lnd_export signs and units, line for line.
  !
  ! Fields kELM does not produce are zero, and say so in lnd.log at init:
  ! dust (Fall_flxdst1..4) and Sl_u10withgusts are not supported;
  ! Flrl_Tqsur/Tqsub (MOSART heat only) are not produced yet. Fields that only exist when a
  ! process kELM does not model is configured (dry deposition, MEGAN,
  ! methane, NH3, soil water to the atmosphere, sediment) abort at init.
  !-----------------------------------------------------------------------

  use shr_kind_mod    , only : r8 => shr_kind_r8
  use shr_sys_mod     , only : shr_sys_abort, shr_sys_flush
  use shr_const_mod   , only : SHR_CONST_STEBOL, SHR_CONST_SPVAL
  use mct_mod         , only : mct_aVect, mct_avect_indexra
  use kelmSpmdMod     , only : masterproc
  use kelmSubgridMod  , only : lun_gridcell, col_landunit, patch_column, patch_wtcol
  use kelmKokkosStateMod, only : n_kokkos_col, n_kokkos_patch, col_of_kcol, &
                                 kcol_of_col, patch_of_kpatch
  use kelm_mod        , only : kELMType, KELM_SUCCESS, kELMLnd2AtmFieldCount, &
                               kELMLnd2AtmSetPatchMap, kELMComputeLnd2Atm, kELMGetLnd2Atm
  use kelm_cpl_indices
  use kelmHistMod     , only : kelm_hist_check_natveg_weight
  use seq_drydep_mod  , only : lnd_drydep
  use shr_megan_mod   , only : shr_megan_mechcomps_n
  use shr_fan_mod     , only : shr_fan_to_atm

  implicit none
  private

  public :: kelm_lnd_export_init
  public :: kelm_kokkos_pull_lnd2atm
  public :: kelm_lnd_export
  public :: kelm_lnd2atm_ready
  public :: kelm_lnd2atm_cell_table
  public :: kelm_lnd2atm_load_cell_table
  public :: NFIELD

  ! Field order of the library's table (Lnd2AtmData, HistData.h :: L2A).
  integer, parameter :: L_FIRE = 1, L_H2OSNO = 2, &
                        L_ALBD_VIS = 3, L_ALBD_NIR = 4, L_ALBI_VIS = 5, L_ALBI_NIR = 6, &
                        L_TREF = 7, L_QREF = 8, L_RAM1 = 9, L_TAUX = 10, L_TAUY = 11, &
                        L_LH = 12, L_SH = 13, L_EVAP = 14, L_FSA = 15, &
                        L_QSUR = 16, L_QSURP = 17, L_QSUB = 18, L_QSUBP = 19, &
                        L_QGWL = 20, L_QICE = 21, L_COSZEN = 22, &
                        L_FV = 23, L_U10 = 24, NFIELD = 24

  real(r8), allocatable :: lnd2atm(:,:)   ! (n_kokkos_col, NFIELD), lnd2atm's sign
  integer , allocatable :: cell_of_kcol(:) ! natural column -> local cell (l2x index)
  integer :: ncells = 0
  logical :: export_ready = .false.

contains

  !-----------------------------------------------------------------------
  subroutine kelm_lnd_export_init(elm, num_cells_owned, l2x, logunit)
    !
    ! Hand the library its own copy of the natural patch map (naturalCol's
    ! arrives on the device only at the first step, and the coupler takes an
    ! export at init), map natural columns to l2x cells, and log every l2x
    ! field as filled, zero-as-ELM, or not produced.
    !
    implicit none
    type(kELMType) , intent(in) :: elm
    integer        , intent(in) :: num_cells_owned
    type(mct_aVect), intent(in) :: l2x
    integer        , intent(in) :: logunit

    integer :: kc, kp, c, cprev, g, ierr
    integer , allocatable :: npfts(:), pfti(:), seen(:)
    real(r8), allocatable :: wt(:)
    character(len=*), parameter :: subname = '(kelm_lnd_export_init) '

    if (kELMLnd2AtmFieldCount() /= NFIELD) &
         call shr_sys_abort(subname//'ERROR: field count differs from the library''s')

    ! Processes kELM does not model: the coupler would get a silent zero.
    if (lnd_drydep) call shr_sys_abort(subname//'ERROR: dry deposition is configured; '// &
         'kELM does not compute deposition velocities')
    if (shr_megan_mechcomps_n > 0) call shr_sys_abort(subname//'ERROR: MEGAN is configured; '// &
         'kELM does not compute VOC emissions')
    if (shr_fan_to_atm) call shr_sys_abort(subname//'ERROR: FAN NH3 to the atmosphere is '// &
         'configured; kELM does not model it')
    call refuse(l2x, 'Fall_methane')
    call refuse(l2x, 'Sl_soilw')
    call refuse(l2x, 'Flrl_rofmud')

    ! Glacier, wetland and crop are not supported (owner, 2026-11-11): every
    ! cell's natural + lake + urban weight must be 1, or the means below would
    ! be partial ones. History aborts on the same condition.
    call kelm_hist_check_natveg_weight(logunit)

    ncells = num_cells_owned
    allocate(npfts(n_kokkos_col), pfti(n_kokkos_col), wt(n_kokkos_patch))
    allocate(cell_of_kcol(n_kokkos_col), seen(ncells))

    ! Patches of a natural column are contiguous (build_column_patch_index
    ! checks the same at the first step).
    npfts = 0; pfti = 0; cprev = -1
    do kp = 1, n_kokkos_patch
       c  = patch_column(patch_of_kpatch(kp))
       kc = kcol_of_col(c) + 1
       if (kc <= 0) call shr_sys_abort(subname//'ERROR: packed patch on an unpacked column')
       if (kc /= cprev) then
          if (npfts(kc) /= 0) call shr_sys_abort(subname//'ERROR: patches not contiguous')
          pfti(kc) = kp - 1
          cprev = kc
       end if
       npfts(kc) = npfts(kc) + 1
       wt(kp) = patch_wtcol(patch_of_kpatch(kp))
    end do
    call kELMLnd2AtmSetPatchMap(elm, npfts, pfti, n_kokkos_col, wt, n_kokkos_patch, ierr)
    if (ierr /= KELM_SUCCESS) call shr_sys_abort(subname//'ERROR: kELMLnd2AtmSetPatchMap failed')

    ! One natural column per cell: its cell is the l2x row.
    seen = 0
    do kc = 1, n_kokkos_col
       g = lun_gridcell(col_landunit(col_of_kcol(kc)))
       if (g < 1 .or. g > ncells) call shr_sys_abort(subname//'ERROR: column cell out of range')
       cell_of_kcol(kc) = g
       seen(g) = seen(g) + 1
    end do
    if (any(seen /= 1)) &
         call shr_sys_abort(subname//'ERROR: every cell must have exactly one natural column')

    allocate(lnd2atm(n_kokkos_col, NFIELD))
    lnd2atm = 0._r8
    export_ready = .true.

    if (masterproc) then
       write(logunit,*) subname,'l2x fields filled: Sl_t Sl_snowh Sl_avsdr Sl_anidr Sl_avsdf ', &
            'Sl_anidf Sl_tref Sl_qref Sl_ram1 Sl_fv Sl_u10 Fall_taux Fall_tauy Fall_lat Fall_sen Fall_lwup ', &
            'Fall_evap Fall_swnet Flrl_rofsur Flrl_rofsub Flrl_rofgwl Flrl_rofi coszen_str'
       write(logunit,*) subname,'l2x fields zero, as ELM sends them in SP without irrigation: ', &
            'Flrl_wslake Flrl_demand Flrl_rofdto Fall_fco2_lnd Flrl_inundinf'
       write(logunit,*) subname,'l2x fields zero, NOT PRODUCED by kELM: Fall_flxdst1..4 and ', &
            'Sl_u10withgusts (not supported); Flrl_Tqsur Flrl_Tqsub (not yet)'
       call shr_sys_flush(logunit)
    end if

    deallocate(npfts, pfti, wt, seen)

  end subroutine kelm_lnd_export_init

  !-----------------------------------------------------------------------
  subroutine kelm_kokkos_pull_lnd2atm(elm)
    !
    ! Form this step's gridcell means on the device and copy them to the
    ! host: one device-to-host copy per coupling interval.
    !
    implicit none
    type(kELMType), intent(in) :: elm
    integer :: ierr

    if (.not. export_ready) call shr_sys_abort('(kelm_kokkos_pull_lnd2atm) ERROR: not initialized')
    call kELMComputeLnd2Atm(elm, ierr)
    if (ierr /= KELM_SUCCESS) call shr_sys_abort('(kelm_kokkos_pull_lnd2atm) ERROR: kELMComputeLnd2Atm failed')
    call kELMGetLnd2Atm(elm, lnd2atm, n_kokkos_col, NFIELD, ierr)
    if (ierr /= KELM_SUCCESS) call shr_sys_abort('(kelm_kokkos_pull_lnd2atm) ERROR: kELMGetLnd2Atm failed')

  end subroutine kelm_kokkos_pull_lnd2atm

  !-----------------------------------------------------------------------
  subroutine kelm_lnd_export(l2x)
    !
    ! ELM: lnd_import_export :: lnd_export, for the fields kELM produces.
    ! Coupler convention: fluxes positive into the surface.
    !
    implicit none
    real(r8), intent(out) :: l2x(:,:)   ! (field, cell): l2x_l%rAttr
    integer :: kc, i

    l2x(:,:) = 0.0_r8

    do kc = 1, size(cell_of_kcol)
       i = cell_of_kcol(kc)
       l2x(index_l2x_Sl_t,i)      =  sqrt(sqrt(lnd2atm(kc,L_FIRE)/SHR_CONST_STEBOL))
       l2x(index_l2x_Sl_snowh,i)  =  lnd2atm(kc,L_H2OSNO)/1000._r8
       l2x(index_l2x_Sl_avsdr,i)  =  lnd2atm(kc,L_ALBD_VIS)
       l2x(index_l2x_Sl_anidr,i)  =  lnd2atm(kc,L_ALBD_NIR)
       l2x(index_l2x_Sl_avsdf,i)  =  lnd2atm(kc,L_ALBI_VIS)
       l2x(index_l2x_Sl_anidf,i)  =  lnd2atm(kc,L_ALBI_NIR)
       l2x(index_l2x_Sl_tref,i)   =  lnd2atm(kc,L_TREF)
       l2x(index_l2x_Sl_qref,i)   =  lnd2atm(kc,L_QREF)
       l2x(index_l2x_Sl_ram1,i)   =  lnd2atm(kc,L_RAM1)
       l2x(index_l2x_Sl_fv,i)     =  lnd2atm(kc,L_FV)
       l2x(index_l2x_Sl_u10,i)    =  lnd2atm(kc,L_U10)
       l2x(index_l2x_Fall_taux,i) = -lnd2atm(kc,L_TAUX)
       l2x(index_l2x_Fall_tauy,i) = -lnd2atm(kc,L_TAUY)
       l2x(index_l2x_Fall_lat,i)  = -lnd2atm(kc,L_LH)
       l2x(index_l2x_Fall_sen,i)  = -lnd2atm(kc,L_SH)
       l2x(index_l2x_Fall_lwup,i) = -lnd2atm(kc,L_FIRE)
       l2x(index_l2x_Fall_evap,i) = -lnd2atm(kc,L_EVAP)
       l2x(index_l2x_Fall_swnet,i)=  lnd2atm(kc,L_FSA)
       l2x(index_l2x_Flrl_rofi,i)   = lnd2atm(kc,L_QICE)
       l2x(index_l2x_Flrl_rofsur,i) = lnd2atm(kc,L_QSUR) + lnd2atm(kc,L_QSURP)
       l2x(index_l2x_Flrl_rofsub,i) = lnd2atm(kc,L_QSUB) + lnd2atm(kc,L_QSUBP)
       l2x(index_l2x_Flrl_rofgwl,i) = lnd2atm(kc,L_QGWL)
       l2x(index_l2x_coszen_str,i)  = lnd2atm(kc,L_COSZEN)
    end do

    ! A missing value reaching the coupler is a bug in the rules above.
    if (any(abs(l2x) >= 0.5_r8*SHR_CONST_SPVAL)) &
         call shr_sys_abort('(kelm_lnd_export) ERROR: spval in the export')

  end subroutine kelm_lnd_export

  !-----------------------------------------------------------------------
  logical function kelm_lnd2atm_ready()
    kelm_lnd2atm_ready = export_ready
  end function kelm_lnd2atm_ready

  !-----------------------------------------------------------------------
  subroutine kelm_lnd2atm_cell_table(elm, tab)
    !
    ! The table the coupler would get from the current state, by local cell:
    ! what the restart carries (kelmRestMod), so a continue run's export at
    ! init is the one the writing run would have sent next.
    !
    implicit none
    type(kELMType), intent(in)  :: elm
    real(r8)      , intent(out) :: tab(:,:)   ! (ncells, NFIELD)
    integer :: kc

    call kelm_kokkos_pull_lnd2atm(elm)
    do kc = 1, size(cell_of_kcol)
       tab(cell_of_kcol(kc), :) = lnd2atm(kc, :)
    end do

  end subroutine kelm_lnd2atm_cell_table

  !-----------------------------------------------------------------------
  subroutine kelm_lnd2atm_load_cell_table(tab)
    !
    ! Take the table from a restart (by local cell) in place of computing it.
    !
    implicit none
    real(r8), intent(in) :: tab(:,:)   ! (ncells, NFIELD)
    integer :: kc

    if (.not. export_ready) call shr_sys_abort('(kelm_lnd2atm_load_cell_table) ERROR: not initialized')
    do kc = 1, size(cell_of_kcol)
       lnd2atm(kc, :) = tab(cell_of_kcol(kc), :)
    end do

  end subroutine kelm_lnd2atm_load_cell_table

  !-----------------------------------------------------------------------
  subroutine refuse(l2x, fldname)
    type(mct_aVect) , intent(in) :: l2x
    character(len=*), intent(in) :: fldname
    if (mct_avect_indexra(l2x, trim(fldname), perrWith='quiet') > 0) &
         call shr_sys_abort('(kelm_lnd_export_init) ERROR: coupler field '//trim(fldname)// &
              ' is configured, but kELM does not model its process')
  end subroutine refuse

end module kelmLndExportMod
