module elmxxIO

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Minimal PIO-based input for ELMxx.
  !
  ! At present ELMxx only needs to read the land domain (fraction) file so it
  ! can hand a domain to the coupler. The whole global grid is read on every
  ! rank -- the same thing RDycore does with latc_g/lonc_g/areac_g in
  ! components/rdycore/src/cpl/rof_comp_mct.F90 -- which is adequate at the
  ! resolutions we run today. Move to a decomposed read when that stops being
  ! true.
  !
  ! This is deliberately NOT a port of components/rdycore/src/main/rdycoreIO.F90,
  ! which is mostly a generic ncd_io wrapper suite that ELMxx does not need.
  !-----------------------------------------------------------------------

  use shr_kind_mod  , only : r8 => shr_kind_r8
  use shr_sys_mod   , only : shr_sys_abort, shr_sys_flush
  use elmxxSpmdMod  , only : masterproc
  use pio

  implicit none
  save
  private

  integer                       , public :: io_type
  type(iosystem_desc_t), pointer, public :: pio_subsystem

  public :: elmxx_pio_init
  public :: elmxx_read_domain

contains

  !-----------------------------------------------------------------------
  subroutine elmxx_pio_init(inst_name)
    !
    ! !DESCRIPTION:
    ! Grab the PIO subsystem the driver set up for this component instance.
    !
    ! !USES:
    use shr_pio_mod, only : shr_pio_getiosys, shr_pio_getiotype
    !
    implicit none
    !
    character(len=*), intent(in) :: inst_name

    pio_subsystem => shr_pio_getiosys(inst_name)
    io_type       =  shr_pio_getiotype(inst_name)

  end subroutine elmxx_pio_init

  !-----------------------------------------------------------------------
  subroutine elmxx_read_domain(iulog, fname, ni, nj, lonc, latc, areac, maskc, fracc)
    !
    ! !DESCRIPTION:
    ! Read a CIME domain file and return the flattened (ni*nj) cell center
    ! coordinates, cell areas, mask and land fraction.
    !
    ! Units are left exactly as they are on the file: lon/lat in degrees and
    ! area in radians^2. That is what the coupler expects for the MCT domain --
    ! ELM multiplies area by re**2 when reading and divides by re*re when
    ! exporting, which is a round trip (see surfrdMod.F90 and lnd_domain_mct in
    ! components/elm/src/cpl/lnd_comp_mct.F90).
    !
    implicit none
    !
    integer           , intent(in)  :: iulog        ! log unit
    character(len=*)  , intent(in)  :: fname        ! domain file name
    integer           , intent(out) :: ni, nj       ! grid dimensions
    real(r8), pointer , intent(out) :: lonc(:)      ! cell center longitude (deg)
    real(r8), pointer , intent(out) :: latc(:)      ! cell center latitude  (deg)
    real(r8), pointer , intent(out) :: areac(:)     ! cell area (radians^2)
    real(r8), pointer , intent(out) :: maskc(:)     ! domain mask (0 or 1)
    real(r8), pointer , intent(out) :: fracc(:)     ! land fraction
    !
    ! !LOCAL VARIABLES:
    type(file_desc_t) :: ncid
    integer           :: dimid, varid
    integer           :: status
    integer           :: gsize
    real(r8), allocatable :: tmp2d(:,:)
    integer , allocatable :: itmp2d(:,:)
    character(len=*), parameter :: subname = 'elmxx_read_domain'

    if (masterproc) then
       write(iulog,*) trim(subname),': reading land domain from ',trim(fname)
       call shr_sys_flush(iulog)
    end if

    status = pio_openfile(pio_subsystem, ncid, io_type, trim(fname), pio_nowrite)
    if (status /= PIO_NOERR) then
       call shr_sys_abort(subname//' ERROR: cannot open domain file '//trim(fname))
    end if

    ! ---- grid dimensions ----
    status = pio_inq_dimid(ncid, 'ni', dimid)
    if (status /= PIO_NOERR) call shr_sys_abort(subname//' ERROR: no ni dimension on '//trim(fname))
    status = pio_inq_dimlen(ncid, dimid, ni)

    status = pio_inq_dimid(ncid, 'nj', dimid)
    if (status /= PIO_NOERR) call shr_sys_abort(subname//' ERROR: no nj dimension on '//trim(fname))
    status = pio_inq_dimlen(ncid, dimid, nj)

    gsize = ni*nj

    allocate(lonc(gsize), latc(gsize), areac(gsize), maskc(gsize), fracc(gsize))
    allocate(tmp2d(ni,nj))
    allocate(itmp2d(ni,nj))

    ! ---- real fields ----
    call read_real2d(ncid, fname, 'xc'  , tmp2d); lonc(:)  = reshape(tmp2d, (/gsize/))
    call read_real2d(ncid, fname, 'yc'  , tmp2d); latc(:)  = reshape(tmp2d, (/gsize/))
    call read_real2d(ncid, fname, 'area', tmp2d); areac(:) = reshape(tmp2d, (/gsize/))
    call read_real2d(ncid, fname, 'frac', tmp2d); fracc(:) = reshape(tmp2d, (/gsize/))

    ! ---- mask is an integer field on CIME domain files ----
    status = pio_inq_varid(ncid, 'mask', varid)
    if (status /= PIO_NOERR) call shr_sys_abort(subname//' ERROR: no mask variable on '//trim(fname))
    status = pio_get_var(ncid, varid, itmp2d)
    if (status /= PIO_NOERR) call shr_sys_abort(subname//' ERROR: cannot read mask from '//trim(fname))
    maskc(:) = real(reshape(itmp2d, (/gsize/)), r8)

    deallocate(tmp2d)
    deallocate(itmp2d)

    call pio_closefile(ncid)

    if (masterproc) then
       write(iulog,*) trim(subname),': ni = ',ni,' nj = ',nj,' ncells = ',gsize
       write(iulog,*) trim(subname),': lon range   ',minval(lonc) ,maxval(lonc)
       write(iulog,*) trim(subname),': lat range   ',minval(latc) ,maxval(latc)
       write(iulog,*) trim(subname),': area range  ',minval(areac),maxval(areac)
       write(iulog,*) trim(subname),': frac range  ',minval(fracc),maxval(fracc)
       write(iulog,*) trim(subname),': mask range  ',minval(maskc),maxval(maskc)
       call shr_sys_flush(iulog)
    end if

  end subroutine elmxx_read_domain

  !-----------------------------------------------------------------------
  subroutine read_real2d(ncid, fname, varname, data)
    !
    ! !DESCRIPTION:
    ! Read a 2D double precision variable in its entirety on every rank.
    !
    implicit none
    !
    type(file_desc_t), intent(inout) :: ncid
    character(len=*) , intent(in)    :: fname
    character(len=*) , intent(in)    :: varname
    real(r8)         , intent(inout) :: data(:,:)
    !
    integer :: varid, status
    character(len=*), parameter :: subname = 'elmxx_read_domain::read_real2d'

    status = pio_inq_varid(ncid, trim(varname), varid)
    if (status /= PIO_NOERR) then
       call shr_sys_abort(subname//' ERROR: no '//trim(varname)//' variable on '//trim(fname))
    end if

    status = pio_get_var(ncid, varid, data)
    if (status /= PIO_NOERR) then
       call shr_sys_abort(subname//' ERROR: cannot read '//trim(varname)//' from '//trim(fname))
    end if

  end subroutine read_real2d

end module elmxxIO
