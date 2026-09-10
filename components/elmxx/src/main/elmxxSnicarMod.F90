module elmxxSnicarMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Reads the two SNICAR lookup tables and pushes them to the device.
  !
  ! WHY THIS IS FORTRAN-SIDE. ELMxx's C++ has no NetCDF dependency and is not
  ! getting one for two static tables. The tables are read here with PIO, the
  ! same way elmxxSurfdataMod reads the surface dataset, and pushed across once
  ! at init. Nothing crosses the boundary per timestep.
  !
  ! THE TWO TABLES:
  !   fsnowoptics  Mie parameters for pure ice grains -- single-scatter
  !                albedo, asymmetry parameter and mass extinction
  !                cross-section, for direct and diffuse incidence, each
  !                (radius_um, wvl) = (1471, 5).
  !   fsnowaging   Best-fit grain growth parameters tau, kappa and drdsdt0,
  !                each (sno_dns, dTdz, T) = (8, 31, 11) as Fortran sees them.
  !
  ! MEMORY LAYOUT. Both C++ sides expect C-ordered arrays -- optics as
  ! (band, radius) and aging as (T, dTdz, density). A Fortran array declared
  ! (radius, band) has exactly the memory of a C (band, radius) array, and
  ! likewise (density, dTdz, T) matches C (T, dTdz, density). So the flat
  ! buffers map straight across and no transpose happens anywhere.
  !
  ! AEROSOLS ARE READ (2026-09-09). The eight species carry no grain-radius
  ! axis -- an aerosol particle's optics do not depend on the snow grain it
  ! sits in -- so each is five numbers, one per band. They are stored
  ! (band, species) here, which is the memory of a C (species, band) array,
  ! matching ELMxxSetSnicarAerosolOptics. Species order is ELM's and is the
  ! contract shared with the deposition and the mixing; see SnicarData.h.
  !-----------------------------------------------------------------------

  use shr_kind_mod  , only : r8 => shr_kind_r8
  use shr_sys_mod   , only : shr_sys_abort, shr_sys_flush
  use elmxxSpmdMod  , only : masterproc
  use elmxxIO       , only : io_type, pio_subsystem
  use pio

  implicit none
  save
  private

  ! Table dimensions, mirroring the C++ constants in SnicarData.h/SnowAgeData.h.
  integer, parameter, public :: snicar_nbnd  = 5
  integer, parameter, public :: snicar_nrds  = 1471
  integer, parameter, public :: snicar_nT    = 11
  integer, parameter, public :: snicar_nTgrd = 31
  integer, parameter, public :: snicar_nrhos = 8

  logical, public :: snicar_tables_read = .false.

  ! (radius, band) -- see the memory-layout note above.
  real(r8), allocatable, public :: ss_alb_drc(:,:), asm_prm_drc(:,:), ext_cff_drc(:,:)
  real(r8), allocatable, public :: ss_alb_dfs(:,:), asm_prm_dfs(:,:), ext_cff_dfs(:,:)
  ! (band, species) -- see the memory-layout note above.
  integer, parameter, public :: snicar_naer = 8
  real(r8), allocatable, public :: ss_alb_aer(:,:), asm_prm_aer(:,:), ext_cff_aer(:,:)
  ! (density, dTdz, T)
  real(r8), allocatable, public :: snowage_tau(:,:,:), snowage_kappa(:,:,:), &
                                   snowage_drdt0(:,:,:)

  public :: elmxx_snicar_read
  public :: elmxx_snicar_clean

contains

  !-----------------------------------------------------------------------
  subroutine elmxx_snicar_read(fsnowoptics, fsnowaging, logunit)
    !
    implicit none
    character(len=*), intent(in) :: fsnowoptics, fsnowaging
    integer         , intent(in) :: logunit
    type(file_desc_t) :: ncid
    integer :: status
    character(len=*), parameter :: subname = '(elmxx_snicar_read) '

    call elmxx_snicar_clean()

    allocate(ss_alb_drc (snicar_nrds, snicar_nbnd), &
             asm_prm_drc(snicar_nrds, snicar_nbnd), &
             ext_cff_drc(snicar_nrds, snicar_nbnd), &
             ss_alb_dfs (snicar_nrds, snicar_nbnd), &
             asm_prm_dfs(snicar_nrds, snicar_nbnd), &
             ext_cff_dfs(snicar_nrds, snicar_nbnd))
    allocate(snowage_tau  (snicar_nrhos, snicar_nTgrd, snicar_nT), &
             snowage_kappa(snicar_nrhos, snicar_nTgrd, snicar_nT), &
             snowage_drdt0(snicar_nrhos, snicar_nTgrd, snicar_nT))

    ! ---- snow optics ----
    status = pio_openfile(pio_subsystem, ncid, io_type, trim(fsnowoptics), pio_nowrite)
    if (status /= PIO_NOERR) then
       call shr_sys_abort(subname//'ERROR: cannot open '//trim(fsnowoptics))
    end if
    call check_dim(ncid, fsnowoptics, 'wvl',       snicar_nbnd)
    call check_dim(ncid, fsnowoptics, 'radius_um', snicar_nrds)
    call read_var2d(ncid, fsnowoptics, 'ss_alb_ice_drc',      ss_alb_drc)
    call read_var2d(ncid, fsnowoptics, 'asm_prm_ice_drc',     asm_prm_drc)
    call read_var2d(ncid, fsnowoptics, 'ext_cff_mss_ice_drc', ext_cff_drc)
    call read_var2d(ncid, fsnowoptics, 'ss_alb_ice_dfs',      ss_alb_dfs)
    call read_var2d(ncid, fsnowoptics, 'asm_prm_ice_dfs',     asm_prm_dfs)
    call read_var2d(ncid, fsnowoptics, 'ext_cff_mss_ice_dfs', ext_cff_dfs)

    ! ---- aerosol Mie parameters, ELM species order 1..8 ----
    allocate(ss_alb_aer (snicar_nbnd, snicar_naer), &
             asm_prm_aer(snicar_nbnd, snicar_naer), &
             ext_cff_aer(snicar_nbnd, snicar_naer))
    call read_aerosol(ncid, fsnowoptics, 'bcphil', 1)
    call read_aerosol(ncid, fsnowoptics, 'bcphob', 2)
    call read_aerosol(ncid, fsnowoptics, 'ocphil', 3)
    call read_aerosol(ncid, fsnowoptics, 'ocphob', 4)
    call read_aerosol(ncid, fsnowoptics, 'dust01', 5)
    call read_aerosol(ncid, fsnowoptics, 'dust02', 6)
    call read_aerosol(ncid, fsnowoptics, 'dust03', 7)
    call read_aerosol(ncid, fsnowoptics, 'dust04', 8)

    call pio_closefile(ncid)

    ! ---- snow aging ----
    status = pio_openfile(pio_subsystem, ncid, io_type, trim(fsnowaging), pio_nowrite)
    if (status /= PIO_NOERR) then
       call shr_sys_abort(subname//'ERROR: cannot open '//trim(fsnowaging))
    end if
    call check_dim(ncid, fsnowaging, 'T',       snicar_nT)
    call check_dim(ncid, fsnowaging, 'dTdz',    snicar_nTgrd)
    call check_dim(ncid, fsnowaging, 'sno_dns', snicar_nrhos)
    call read_var3d(ncid, fsnowaging, 'tau',     snowage_tau)
    call read_var3d(ncid, fsnowaging, 'kappa',   snowage_kappa)
    call read_var3d(ncid, fsnowaging, 'drdsdt0', snowage_drdt0)
    call pio_closefile(ncid)

    snicar_tables_read = .true.

    if (masterproc) then
       write(logunit,*) subname,'read SNICAR tables'
       write(logunit,*) '    optics: ',trim(fsnowoptics)
       write(logunit,*) '    aging : ',trim(fsnowaging)
       ! ELM prints these three at startup for T=263K, dTdz=100 K/m,
       ! rhos=150 kg/m3. Printing the same numbers makes a transposed or
       ! mis-indexed read visible immediately instead of as a wrong albedo
       ! three thousand timesteps later.
       write(logunit,*) '    tau  (3,11,9) = ',snowage_tau  (3,11,9)
       write(logunit,*) '    kappa(3,11,9) = ',snowage_kappa(3,11,9)
       write(logunit,*) '    drdt0(3,11,9) = ',snowage_drdt0(3,11,9)
       call shr_sys_flush(logunit)
    end if

  end subroutine elmxx_snicar_read

  !-----------------------------------------------------------------------
  subroutine elmxx_snicar_clean()
    implicit none
    if (allocated(ss_alb_aer   )) deallocate(ss_alb_aer)
    if (allocated(asm_prm_aer  )) deallocate(asm_prm_aer)
    if (allocated(ext_cff_aer  )) deallocate(ext_cff_aer)
    if (allocated(ss_alb_drc   )) deallocate(ss_alb_drc)
    if (allocated(asm_prm_drc  )) deallocate(asm_prm_drc)
    if (allocated(ext_cff_drc  )) deallocate(ext_cff_drc)
    if (allocated(ss_alb_dfs   )) deallocate(ss_alb_dfs)
    if (allocated(asm_prm_dfs  )) deallocate(asm_prm_dfs)
    if (allocated(ext_cff_dfs  )) deallocate(ext_cff_dfs)
    if (allocated(snowage_tau  )) deallocate(snowage_tau)
    if (allocated(snowage_kappa)) deallocate(snowage_kappa)
    if (allocated(snowage_drdt0)) deallocate(snowage_drdt0)
    snicar_tables_read = .false.
  end subroutine elmxx_snicar_clean

  !-----------------------------------------------------------------------
  subroutine check_dim(ncid, fname, dimname, expected)
    !
    ! A table whose shape is not what the C++ side allocates would be pushed
    ! into the wrong slots and produce a plausible-looking albedo, so this
    ! aborts rather than trusting the file.
    !
    implicit none
    type(file_desc_t), intent(inout) :: ncid
    character(len=*) , intent(in)    :: fname, dimname
    integer          , intent(in)    :: expected
    integer :: dimid, dlen, status
    character(len=*), parameter :: subname = '(elmxx_snicar_read::check_dim) '

    status = pio_inq_dimid(ncid, trim(dimname), dimid)
    if (status /= PIO_NOERR) then
       call shr_sys_abort(subname//'ERROR: no '//trim(dimname)//' on '//trim(fname))
    end if
    status = pio_inq_dimlen(ncid, dimid, dlen)
    if (status /= PIO_NOERR .or. dlen /= expected) then
       call shr_sys_abort(subname//'ERROR: '//trim(dimname)//' on '// &
            trim(fname)//' is not the expected length')
    end if
  end subroutine check_dim

  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  subroutine read_aerosol(ncid, fname, species, idx)
    !
    ! One aerosol species' three Mie parameters into column `idx`.
    !
    implicit none
    type(file_desc_t), intent(inout) :: ncid
    character(len=*) , intent(in)    :: fname, species
    integer          , intent(in)    :: idx
    real(r8) :: buf(snicar_nbnd)

    call read_var1d(ncid, fname, 'ss_alb_'//trim(species),      buf)
    ss_alb_aer(:, idx)  = buf
    call read_var1d(ncid, fname, 'asm_prm_'//trim(species),     buf)
    asm_prm_aer(:, idx) = buf
    call read_var1d(ncid, fname, 'ext_cff_mss_'//trim(species), buf)
    ext_cff_aer(:, idx) = buf

  end subroutine read_aerosol

  !-----------------------------------------------------------------------
  subroutine read_var1d(ncid, fname, varname, out)
    !
    implicit none
    type(file_desc_t), intent(inout) :: ncid
    character(len=*) , intent(in)    :: fname, varname
    real(r8)         , intent(out)   :: out(:)
    type(var_desc_t) :: vardesc
    integer :: status

    status = pio_inq_varid(ncid, trim(varname), vardesc)
    if (status /= PIO_NOERR) then
       call shr_sys_abort('(read_var1d) ERROR: '//trim(varname)// &
                          ' not found in '//trim(fname))
    end if
    status = pio_get_var(ncid, vardesc, out)
    if (status /= PIO_NOERR) then
       call shr_sys_abort('(read_var1d) ERROR: cannot read '//trim(varname))
    end if

  end subroutine read_var1d

  !-----------------------------------------------------------------------
  subroutine read_var2d(ncid, fname, varname, out)
    implicit none
    type(file_desc_t), intent(inout) :: ncid
    character(len=*) , intent(in)    :: fname, varname
    real(r8)         , intent(out)   :: out(:,:)
    integer :: varid, status
    character(len=*), parameter :: subname = '(elmxx_snicar_read::read_var2d) '

    status = pio_inq_varid(ncid, trim(varname), varid)
    if (status /= PIO_NOERR) then
       call shr_sys_abort(subname//'ERROR: no '//trim(varname)//' on '//trim(fname))
    end if
    status = pio_get_var(ncid, varid, out)
    if (status /= PIO_NOERR) then
       call shr_sys_abort(subname//'ERROR: cannot read '//trim(varname)//' from '//trim(fname))
    end if
  end subroutine read_var2d

  !-----------------------------------------------------------------------
  subroutine read_var3d(ncid, fname, varname, out)
    implicit none
    type(file_desc_t), intent(inout) :: ncid
    character(len=*) , intent(in)    :: fname, varname
    real(r8)         , intent(out)   :: out(:,:,:)
    integer :: varid, status
    character(len=*), parameter :: subname = '(elmxx_snicar_read::read_var3d) '

    status = pio_inq_varid(ncid, trim(varname), varid)
    if (status /= PIO_NOERR) then
       call shr_sys_abort(subname//'ERROR: no '//trim(varname)//' on '//trim(fname))
    end if
    status = pio_get_var(ncid, varid, out)
    if (status /= PIO_NOERR) then
       call shr_sys_abort(subname//'ERROR: cannot read '//trim(varname)//' from '//trim(fname))
    end if
  end subroutine read_var3d

end module elmxxSnicarMod
