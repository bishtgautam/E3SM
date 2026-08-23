module elmxxFinidatMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Read an ELM restart file as an initial condition (finidat).
  !
  ! WHY THIS EXISTS. Every ELMxx kernel matches ELM to ~1e-14 when replayed
  ! against ELM's recorded inputs, and the coupled model is still ~0.9 K off.
  ! The replay harness cannot close that gap: it seeds each kernel from ELM
  ! and so bypasses the coupling, which is exactly where the divergence is.
  ! Starting both models from bit-identical state and stepping ONCE removes
  ! accumulation from the picture -- whatever differs was caused in that step.
  !
  ! READ-ONLY, DELIBERATELY. ELMxx does not write restarts. This is a
  ! debugging instrument, not restart capability; do not mistake it for one.
  !
  ! NETCDF DIMENSION ORDER. A CDL declaration H2OSOI_LIQ(column, levtot) means
  ! levtot varies FASTEST in memory, so the matching Fortran array is
  ! (levtot, column) -- reversed. Declaring it (column, levtot) reads the file
  ! transposed, which is silent: the shapes are both 2D and PIO does not
  ! complain. It surfaced as spval appearing mid-array in the soil profile.
  !
  ! ORDERING. The restart stores snow layers ELM-ordered: index NLEVSNO-1 is
  ! the soil-adjacent layer, and for snl = -n only the last n slots are live.
  ! ELMxx's `_sno` arrays are H7 (slot 0 soil-adjacent). The setters flip on
  ! the way in, so everything here stays in the file's own ordering.
  !
  ! WHAT IS NOT SEEDED, AND WHY IT IS SAFE. frac_iceold and do_capsnow are not
  ! in the restart; ELM recomputes both before first use each step. Anything
  ! else absent is a hole in the premise, which is why elmxx_finidat_verify
  ! exists -- see the note there.
  !-----------------------------------------------------------------------

  use shr_kind_mod  , only : r8 => shr_kind_r8
  use shr_sys_mod   , only : shr_sys_abort, shr_sys_flush
  use elmxxSpmdMod  , only : masterproc
  use elmxxIO       , only : io_type, pio_subsystem
  use pio

  implicit none
  save
  private

  integer, parameter :: nlevsno_r = 5
  integer, parameter :: nlevgrnd_r = 15
  integer, parameter :: nlevtot_r  = nlevsno_r + nlevgrnd_r   ! 20

  logical, public :: finidat_read = .false.

  ! ---- column state, ELM ordering ----
  integer , allocatable, public :: fi_snl(:)               ! (ncol)
  real(r8), allocatable, public :: fi_t_soisno(:,:)        ! (nlevtot, ncol) -- see the ordering note
  real(r8), allocatable, public :: fi_h2osoi_liq(:,:)      ! (nlevtot, ncol) -- see the ordering note
  real(r8), allocatable, public :: fi_h2osoi_ice(:,:)      ! (nlevtot, ncol) -- see the ordering note
  real(r8), allocatable, public :: fi_dzsno(:,:)           ! (nlevsno, ncol) -- see the ordering note
  real(r8), allocatable, public :: fi_zsno(:,:)            ! (nlevsno, ncol) -- see the ordering note
  real(r8), allocatable, public :: fi_zisno(:,:)           ! (nlevsno, ncol) -- see the ordering note
  real(r8), allocatable, public :: fi_snw_rds(:,:)         ! (nlevsno, ncol) -- see the ordering note
  real(r8), allocatable, public :: fi_qflx_snofrz_lyr(:,:) ! (nlevsno, ncol) -- see the ordering note
  real(r8), allocatable, public :: fi_snow_depth(:), fi_h2osno(:), fi_int_snow(:)
  real(r8), allocatable, public :: fi_frac_sno(:), fi_frac_sno_eff(:)
  real(r8), allocatable, public :: fi_t_grnd(:), fi_t_h2osfc(:), fi_h2osfc(:)
  real(r8), allocatable, public :: fi_frac_h2osfc(:)
  real(r8), allocatable, public :: fi_coszen(:), fi_wa(:), fi_zwt(:)
  real(r8), allocatable, public :: fi_albgrd(:,:), fi_albgri(:,:)
  real(r8), allocatable, public :: fi_flx_absdv(:,:), fi_flx_absdn(:,:)
  real(r8), allocatable, public :: fi_flx_absiv(:,:), fi_flx_absin(:,:)

  ! ---- patch state ----
  real(r8), allocatable, public :: fi_t_veg(:), fi_h2ocan(:), fi_fwet(:)
  real(r8), allocatable, public :: fi_elai(:), fi_esai(:), fi_htop(:)
  real(r8), allocatable, public :: fi_albd(:,:), fi_albi(:,:)

  integer, public :: fi_ncol = 0, fi_npft = 0

  public :: elmxx_finidat_read
  public :: elmxx_finidat_apply_soilprop
  public :: elmxx_finidat_clean

contains

  !-----------------------------------------------------------------------
  subroutine elmxx_finidat_read(fname, ncol_expect, npft_expect, logunit)
    !
    implicit none
    character(len=*), intent(in) :: fname
    integer         , intent(in) :: ncol_expect, npft_expect, logunit
    type(file_desc_t) :: ncid
    integer :: status
    character(len=*), parameter :: subname = '(elmxx_finidat_read) '

    call elmxx_finidat_clean()

    status = pio_openfile(pio_subsystem, ncid, io_type, trim(fname), pio_nowrite)
    if (status /= PIO_NOERR) then
       call shr_sys_abort(subname//'ERROR: cannot open finidat '//trim(fname))
    end if

    fi_ncol = get_dimlen(ncid, fname, 'column')
    fi_npft = get_dimlen(ncid, fname, 'pft')

    ! A restart from a different grid would read cleanly and seed nonsense, so
    ! the shapes are checked rather than assumed.
    if (fi_ncol /= ncol_expect .or. fi_npft /= npft_expect) then
       write(logunit,*) subname,'finidat has column=',fi_ncol,' pft=',fi_npft, &
            ' but ELMxx has ',ncol_expect,' and ',npft_expect
       call shr_sys_abort(subname//'ERROR: finidat does not match this configuration')
    end if
    call check_dim(ncid, fname, 'levsno',  nlevsno_r)
    call check_dim(ncid, fname, 'levtot',  nlevtot_r)

    allocate(fi_snl(fi_ncol))
    allocate(fi_t_soisno  (nlevtot_r, fi_ncol), &
             fi_h2osoi_liq(nlevtot_r, fi_ncol), &
             fi_h2osoi_ice(nlevtot_r, fi_ncol))
    allocate(fi_dzsno(nlevsno_r, fi_ncol), fi_zsno(nlevsno_r, fi_ncol), &
             fi_zisno(nlevsno_r, fi_ncol), fi_snw_rds(nlevsno_r, fi_ncol), &
             fi_qflx_snofrz_lyr(nlevsno_r, fi_ncol))
    allocate(fi_snow_depth(fi_ncol), fi_h2osno(fi_ncol), fi_int_snow(fi_ncol), &
             fi_frac_sno(fi_ncol), fi_frac_sno_eff(fi_ncol), &
             fi_t_grnd(fi_ncol), fi_t_h2osfc(fi_ncol), fi_h2osfc(fi_ncol), &
             fi_frac_h2osfc(fi_ncol), fi_coszen(fi_ncol), &
             fi_wa(fi_ncol), fi_zwt(fi_ncol))
    allocate(fi_albgrd(2, fi_ncol), fi_albgri(2, fi_ncol))
    allocate(fi_flx_absdv(nlevsno_r+1, fi_ncol), fi_flx_absdn(nlevsno_r+1, fi_ncol), &
             fi_flx_absiv(nlevsno_r+1, fi_ncol), fi_flx_absin(nlevsno_r+1, fi_ncol))
    allocate(fi_t_veg(fi_npft), fi_h2ocan(fi_npft), fi_fwet(fi_npft), &
             fi_elai(fi_npft), fi_esai(fi_npft), fi_htop(fi_npft))
    allocate(fi_albd(2, fi_npft), fi_albi(2, fi_npft))

    call read_int1d (ncid, fname, 'SNLSNO'      , fi_snl)
    call read_real2d(ncid, fname, 'T_SOISNO'    , fi_t_soisno)
    call read_real2d(ncid, fname, 'H2OSOI_LIQ'  , fi_h2osoi_liq)
    call read_real2d(ncid, fname, 'H2OSOI_ICE'  , fi_h2osoi_ice)
    call read_real2d(ncid, fname, 'DZSNO'       , fi_dzsno)
    call read_real2d(ncid, fname, 'ZSNO'        , fi_zsno)
    call read_real2d(ncid, fname, 'ZISNO'       , fi_zisno)
    call read_real2d(ncid, fname, 'snw_rds'     , fi_snw_rds)
    call read_real2d(ncid, fname, 'qflx_snofrz_lyr', fi_qflx_snofrz_lyr)
    call read_real1d(ncid, fname, 'SNOW_DEPTH'  , fi_snow_depth)
    call read_real1d(ncid, fname, 'H2OSNO'      , fi_h2osno)
    call read_real1d(ncid, fname, 'INT_SNOW'    , fi_int_snow)
    call read_real1d(ncid, fname, 'frac_sno'    , fi_frac_sno)
    call read_real1d(ncid, fname, 'frac_sno_eff', fi_frac_sno_eff)
    call read_real1d(ncid, fname, 'T_GRND'      , fi_t_grnd)
    call read_real1d(ncid, fname, 'TH2OSFC'     , fi_t_h2osfc)
    call read_real1d(ncid, fname, 'H2OSFC'      , fi_h2osfc)
    call read_real1d(ncid, fname, 'FH2OSFC'     , fi_frac_h2osfc)
    call read_real1d(ncid, fname, 'coszen'      , fi_coszen)
    call read_real1d(ncid, fname, 'WA'          , fi_wa)
    call read_real1d(ncid, fname, 'ZWT'         , fi_zwt)
    call read_real2d(ncid, fname, 'albgrd'      , fi_albgrd)
    call read_real2d(ncid, fname, 'albgri'      , fi_albgri)
    call read_real2d(ncid, fname, 'flx_absdv'   , fi_flx_absdv)
    call read_real2d(ncid, fname, 'flx_absdn'   , fi_flx_absdn)
    call read_real2d(ncid, fname, 'flx_absiv'   , fi_flx_absiv)
    call read_real2d(ncid, fname, 'flx_absin'   , fi_flx_absin)

    call read_real1d(ncid, fname, 'T_VEG'  , fi_t_veg)
    call read_real1d(ncid, fname, 'H2OCAN' , fi_h2ocan)
    call read_real1d(ncid, fname, 'FWET'   , fi_fwet)
    call read_real1d(ncid, fname, 'elai'   , fi_elai)
    call read_real1d(ncid, fname, 'esai'   , fi_esai)
    call read_real1d(ncid, fname, 'htop'   , fi_htop)
    call read_real2d(ncid, fname, 'albd'   , fi_albd)
    call read_real2d(ncid, fname, 'albi'   , fi_albi)

    call pio_closefile(ncid)
    finidat_read = .true.

    if (masterproc) then
       write(logunit,*) subname,'read finidat ',trim(fname)
       write(logunit,*) '    columns = ',fi_ncol,'  patches = ',fi_npft
       write(logunit,*) '    snl(1)  = ',fi_snl(1),'  h2osno(1) = ',fi_h2osno(1)
       write(logunit,*) '    t_grnd(1) = ',fi_t_grnd(1)
       call shr_sys_flush(logunit)
    end if

  end subroutine elmxx_finidat_read

  !-----------------------------------------------------------------------
  subroutine elmxx_finidat_clean()
    implicit none
    if (allocated(fi_snl))               deallocate(fi_snl)
    if (allocated(fi_t_soisno))          deallocate(fi_t_soisno)
    if (allocated(fi_h2osoi_liq))        deallocate(fi_h2osoi_liq)
    if (allocated(fi_h2osoi_ice))        deallocate(fi_h2osoi_ice)
    if (allocated(fi_dzsno))             deallocate(fi_dzsno)
    if (allocated(fi_zsno))              deallocate(fi_zsno)
    if (allocated(fi_zisno))             deallocate(fi_zisno)
    if (allocated(fi_snw_rds))           deallocate(fi_snw_rds)
    if (allocated(fi_qflx_snofrz_lyr))   deallocate(fi_qflx_snofrz_lyr)
    if (allocated(fi_snow_depth))        deallocate(fi_snow_depth)
    if (allocated(fi_h2osno))            deallocate(fi_h2osno)
    if (allocated(fi_int_snow))          deallocate(fi_int_snow)
    if (allocated(fi_frac_sno))          deallocate(fi_frac_sno)
    if (allocated(fi_frac_sno_eff))      deallocate(fi_frac_sno_eff)
    if (allocated(fi_t_grnd))            deallocate(fi_t_grnd)
    if (allocated(fi_t_h2osfc))          deallocate(fi_t_h2osfc)
    if (allocated(fi_h2osfc))            deallocate(fi_h2osfc)
    if (allocated(fi_frac_h2osfc))       deallocate(fi_frac_h2osfc)
    if (allocated(fi_coszen))            deallocate(fi_coszen)
    if (allocated(fi_wa))                deallocate(fi_wa)
    if (allocated(fi_zwt))               deallocate(fi_zwt)
    if (allocated(fi_albgrd))            deallocate(fi_albgrd)
    if (allocated(fi_albgri))            deallocate(fi_albgri)
    if (allocated(fi_flx_absdv))         deallocate(fi_flx_absdv)
    if (allocated(fi_flx_absdn))         deallocate(fi_flx_absdn)
    if (allocated(fi_flx_absiv))         deallocate(fi_flx_absiv)
    if (allocated(fi_flx_absin))         deallocate(fi_flx_absin)
    if (allocated(fi_t_veg))             deallocate(fi_t_veg)
    if (allocated(fi_h2ocan))            deallocate(fi_h2ocan)
    if (allocated(fi_fwet))              deallocate(fi_fwet)
    if (allocated(fi_elai))              deallocate(fi_elai)
    if (allocated(fi_esai))              deallocate(fi_esai)
    if (allocated(fi_htop))              deallocate(fi_htop)
    if (allocated(fi_albd))              deallocate(fi_albd)
    if (allocated(fi_albi))              deallocate(fi_albi)
    finidat_read = .false.
    fi_ncol = 0; fi_npft = 0
  end subroutine elmxx_finidat_clean

  !-----------------------------------------------------------------------
  integer function get_dimlen(ncid, fname, dimname)
    implicit none
    type(file_desc_t), intent(inout) :: ncid
    character(len=*) , intent(in)    :: fname, dimname
    integer :: dimid, status
    character(len=*), parameter :: subname = '(elmxx_finidat_read::get_dimlen) '
    status = pio_inq_dimid(ncid, trim(dimname), dimid)
    if (status /= PIO_NOERR) call shr_sys_abort(subname//'ERROR: no '// &
         trim(dimname)//' on '//trim(fname))
    status = pio_inq_dimlen(ncid, dimid, get_dimlen)
    if (status /= PIO_NOERR) call shr_sys_abort(subname//'ERROR: cannot read '// &
         trim(dimname)//' from '//trim(fname))
  end function get_dimlen

  !-----------------------------------------------------------------------
  subroutine check_dim(ncid, fname, dimname, expected)
    implicit none
    type(file_desc_t), intent(inout) :: ncid
    character(len=*) , intent(in)    :: fname, dimname
    integer          , intent(in)    :: expected
    character(len=*), parameter :: subname = '(elmxx_finidat_read::check_dim) '
    if (get_dimlen(ncid, fname, dimname) /= expected) then
       call shr_sys_abort(subname//'ERROR: '//trim(dimname)//' on '// &
            trim(fname)//' is not the expected length')
    end if
  end subroutine check_dim

  !-----------------------------------------------------------------------
  subroutine read_real1d(ncid, fname, varname, out)
    implicit none
    type(file_desc_t), intent(inout) :: ncid
    character(len=*) , intent(in)    :: fname, varname
    real(r8)         , intent(out)   :: out(:)
    integer :: varid, status
    character(len=*), parameter :: subname = '(elmxx_finidat_read::read_real1d) '
    status = pio_inq_varid(ncid, trim(varname), varid)
    if (status /= PIO_NOERR) call shr_sys_abort(subname//'ERROR: no '// &
         trim(varname)//' on '//trim(fname))
    status = pio_get_var(ncid, varid, out)
    if (status /= PIO_NOERR) call shr_sys_abort(subname//'ERROR: cannot read '// &
         trim(varname)//' from '//trim(fname))
  end subroutine read_real1d

  !-----------------------------------------------------------------------
  subroutine read_real2d(ncid, fname, varname, out)
    implicit none
    type(file_desc_t), intent(inout) :: ncid
    character(len=*) , intent(in)    :: fname, varname
    real(r8)         , intent(out)   :: out(:,:)
    integer :: varid, status
    character(len=*), parameter :: subname = '(elmxx_finidat_read::read_real2d) '
    status = pio_inq_varid(ncid, trim(varname), varid)
    if (status /= PIO_NOERR) call shr_sys_abort(subname//'ERROR: no '// &
         trim(varname)//' on '//trim(fname))
    status = pio_get_var(ncid, varid, out)
    if (status /= PIO_NOERR) call shr_sys_abort(subname//'ERROR: cannot read '// &
         trim(varname)//' from '//trim(fname))
  end subroutine read_real2d

  !-----------------------------------------------------------------------
  subroutine read_int1d(ncid, fname, varname, out)
    implicit none
    type(file_desc_t), intent(inout) :: ncid
    character(len=*) , intent(in)    :: fname, varname
    integer          , intent(out)   :: out(:)
    integer :: varid, status
    character(len=*), parameter :: subname = '(elmxx_finidat_read::read_int1d) '
    status = pio_inq_varid(ncid, trim(varname), varid)
    if (status /= PIO_NOERR) call shr_sys_abort(subname//'ERROR: no '// &
         trim(varname)//' on '//trim(fname))
    status = pio_get_var(ncid, varid, out)
    if (status /= PIO_NOERR) call shr_sys_abort(subname//'ERROR: cannot read '// &
         trim(varname)//' from '//trim(fname))
  end subroutine read_int1d


  !-----------------------------------------------------------------------
  subroutine elmxx_finidat_apply_soilprop(logunit)
    !
    ! !DESCRIPTION:
    ! Overwrite elmxxSoilPropMod's column state with the restart's.
    !
    ! WHY THIS IS SEPARATE FROM THE KOKKOS PUSH, AND WHY IT IS NOT OPTIONAL.
    ! elmxx_soil_kernel_init runs on the FIRST STEP, not at init, because it
    ! needs the coupling timestep. It re-pushes t_soisno / h2osoi_liq /
    ! h2osoi_ice to the device from THESE arrays. So seeding the device at
    ! init and stopping there gets silently undone one step later: the
    ! snowpack comes from ELM and the soil column reverts to cold start.
    !
    ! That is not a hypothetical. It is what happened: soil liquid held its
    ! cold-start profile (2.63, 4.14, 6.82 ...) against ELM's (0.77, 2.07,
    ! 3.45 ...) with all the ice missing, and the resulting inconsistent
    ! column blew up in the snow-layer kernel.
    !
    ! These arrays are on the FULL column grid, the same shape as the restart,
    ! so no gather is needed here -- unlike the Kokkos push, which packs.
    !
    use elmxxSoilPropMod, only : soil_prop_built, col_t_soisno, &
                                 col_h2osoi_liq, col_h2osoi_ice
    implicit none
    integer, intent(in) :: logunit
    integer :: c, j
    character(len=*), parameter :: subname = '(elmxx_finidat_apply_soilprop) '

    if (.not. finidat_read) then
       call shr_sys_abort(subname//'ERROR: finidat not read')
    end if
    if (.not. soil_prop_built) then
       call shr_sys_abort(subname//'ERROR: soil properties not built yet')
    end if
    if (size(col_t_soisno,1) /= fi_ncol .or. &
        size(col_t_soisno,2) /= nlevtot_r) then   ! soil-prop side is (col, lev)
       call shr_sys_abort(subname//'ERROR: soil-prop arrays are not the '// &
            'shape of the restart')
    end if

    do c = 1, fi_ncol
       do j = 1, nlevtot_r
          col_t_soisno  (c,j) = fi_t_soisno  (j,c)
          col_h2osoi_liq(c,j) = fi_h2osoi_liq(j,c)
          col_h2osoi_ice(c,j) = fi_h2osoi_ice(j,c)
       end do
    end do

    if (masterproc) then
       write(logunit,*) subname,'overwrote soil-prop column state from finidat'
       call shr_sys_flush(logunit)
    end if

  end subroutine elmxx_finidat_apply_soilprop

end module elmxxFinidatMod
