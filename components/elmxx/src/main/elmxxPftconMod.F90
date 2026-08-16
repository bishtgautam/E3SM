module elmxxPftconMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! PFT parameters, read from ELM's parameter file.
  !
  ! This is the first thing in ELMxx to read `clm_params`. Everything before
  ! it came from surfdata or the domain; root distribution and stomatal
  ! response are properties of a plant functional type, not of a gridcell, so
  ! they live here instead.
  !
  ! DELIBERATELY MINIMAL. The parameter file carries hundreds of fields, and
  ! reading them all would be a large surface with nothing to check it
  ! against. Only what a ported kernel actually needs is read, so that every
  ! field here has a consumer. Add one when a kernel needs it, not before.
  !
  ! Indexing is ELM's: PFT 0 is bare ground (`noveg`) and gets no roots and no
  ! stomatal response, which is why its roota_par and smpso are zero on file.
  ! `patch_itype` from elmxxSubgridMod is already 0-based, so it indexes these
  ! arrays directly with no translation.
  !-----------------------------------------------------------------------

  use shr_kind_mod  , only : r8 => shr_kind_r8
  use shr_sys_mod   , only : shr_sys_abort, shr_sys_flush
  use elmxxSpmdMod  , only : masterproc, iam
  use elmxxIO       , only : io_type, pio_subsystem
  use pio

  implicit none
  save
  private

  integer, public :: npft_param = 0    ! PFTs on the parameter file

  ! Root distribution, Zeng (1998) two-exponential profile [1/m].
  real(r8), public, pointer :: roota_par(:) => null()   ! (0:npft_param-1)
  real(r8), public, pointer :: rootb_par(:) => null()

  ! Soil matric potential at full stomatal opening and at closure [mm].
  real(r8), public, pointer :: smpso(:) => null()
  real(r8), public, pointer :: smpsc(:) => null()

  ! Canopy roughness and displacement, as RATIOS of canopy top height.
  ! CanopyTemperature forms z0mv = z0mr*htop and displa = displar*htop, and
  ! every aerodynamic resistance in CanopyFluxes divides by something derived
  ! from them -- so leaving these at zero does not damp the canopy, it makes
  ! ustar zero and every flux NaN. See STATUS.
  real(r8), public, pointer :: z0mr(:)   => null()
  real(r8), public, pointer :: displar(:) => null()

  ! Leaf characteristic dimension [m]. CanopyFluxes forms the leaf boundary
  ! layer conductance as cf = 0.01/(sqrt(uaf)*sqrt(dleaf)), so a zero here
  ! makes cf infinite and rb exactly zero -- and rb sits in the denominator
  ! of rppdry. Zero dleaf is therefore not a small leaf; it is a NaN.
  real(r8), public, pointer :: dleaf(:) => null()

  ! Canopy optical properties, per PFT and waveband (1 = VIS, 2 = NIR).
  ! rhol/rhos are leaf and stem reflectance, taul/taus leaf and stem
  ! transmittance, xl the leaf/stem orientation index. These are the inputs to
  ! the Sellers two-stream in elmxxSurfaceAlbedoMod. On file each band is its
  ! own variable (rholvis/rholnir and so on) rather than a banded array.
  real(r8), public, pointer :: rhol(:,:) => null()   ! (0:npft-1, numrad)
  real(r8), public, pointer :: rhos(:,:) => null()
  real(r8), public, pointer :: taul(:,:) => null()
  real(r8), public, pointer :: taus(:,:) => null()
  real(r8), public, pointer :: xl(:)     => null()

  ! Photosynthesis PFT parameters. Nine that genuinely vary between PFTs...
  real(r8), public, pointer :: c3psn(:)    => null()   ! 1 = C3, 0 = C4
  real(r8), public, pointer :: leafcn(:)   => null()   ! leaf C:N [gC/gN]
  real(r8), public, pointer :: flnr(:)     => null()   ! leaf N in Rubisco
  real(r8), public, pointer :: fnitr(:)    => null()   ! foliage N limitation
  real(r8), public, pointer :: slatop(:)   => null()   ! specific leaf area [m2/gC]
  real(r8), public, pointer :: qe_ps(:)    => null()   ! quantum efficiency, C4
  real(r8), public, pointer :: theta_cj(:) => null()   ! ac/aj co-limitation
  real(r8), public, pointer :: bbbopt(:)   => null()   ! Ball-Berry intercept
  real(r8), public, pointer :: mbbopt(:)   => null()   ! Ball-Berry slope

  ! ...and fourteen that are UNIFORM across every PFT on clm_params. They are
  ! read per PFT anyway and the uniformity is ASSERTED, because the C++ side
  ! carries them as scalars and a parameter file that broke the assumption
  ! would otherwise be silently wrong for every PFT but the first.
  ! Order matches ELMxxSetPhotoUniform.
  integer, parameter, public :: n_photo_uniform = 14
  real(r8), public :: photo_uniform(n_photo_uniform) = 0.0_r8
  character(len=8), parameter, public :: photo_uniform_names(n_photo_uniform) = &
       (/ 'fnr     ', 'act25   ', 'kcha    ', 'koha    ', 'cpha    ', &
          'vcmaxha ', 'jmaxha  ', 'tpuha   ', 'lmrha   ', 'vcmaxhd ', &
          'jmaxhd  ', 'tpuhd   ', 'lmrhd   ', 'lmrse   ' /)

  ! Critical soil temperature for soil water stress [C]. Scalar on file
  ! (dimension allpfts = 1), not per PFT, despite living with the PFT params.
  real(r8), public :: tc_stress = 0.0_r8

  logical, public :: pftcon_read = .false.

  public :: elmxx_read_pftcon
  public :: elmxx_pftcon_clean

contains

  !-----------------------------------------------------------------------
  subroutine elmxx_read_pftcon(logunit, fname)
    implicit none
    integer         , intent(in) :: logunit
    character(len=*), intent(in) :: fname
    !
    type(file_desc_t) :: ncid
    integer :: status, dimid, varid
    real(r8) :: scalar(1)
    character(len=*), parameter :: subname = '(elmxx_read_pftcon) '

    if (len_trim(fname) == 0) then
       call shr_sys_abort(subname//'ERROR: no parameter file given (fparamfile)')
    end if

    status = pio_openfile(pio_subsystem, ncid, io_type, trim(fname), pio_nowrite)
    if (status /= PIO_NOERR) then
       call shr_sys_abort(subname//'ERROR: cannot open '//trim(fname))
    end if

    status = pio_inq_dimid(ncid, 'pft', dimid)
    if (status /= PIO_NOERR) then
       call shr_sys_abort(subname//'ERROR: no pft dimension on '//trim(fname))
    end if
    status = pio_inq_dimlen(ncid, dimid, npft_param)

    call elmxx_pftcon_clean()
    allocate(roota_par(0:npft_param-1), rootb_par(0:npft_param-1), &
             smpso(0:npft_param-1), smpsc(0:npft_param-1), &
             z0mr(0:npft_param-1), displar(0:npft_param-1), &
             dleaf(0:npft_param-1), xl(0:npft_param-1), &
             rhol(0:npft_param-1,2), rhos(0:npft_param-1,2), &
             taul(0:npft_param-1,2), taus(0:npft_param-1,2), &
             c3psn(0:npft_param-1), leafcn(0:npft_param-1), &
             flnr(0:npft_param-1), fnitr(0:npft_param-1), &
             slatop(0:npft_param-1), qe_ps(0:npft_param-1), &
             theta_cj(0:npft_param-1), bbbopt(0:npft_param-1), &
             mbbopt(0:npft_param-1))

    call read_pft_real(ncid, fname, 'roota_par', roota_par)
    call read_pft_real(ncid, fname, 'rootb_par', rootb_par)
    call read_pft_real(ncid, fname, 'smpso'    , smpso)
    call read_pft_real(ncid, fname, 'smpsc'    , smpsc)
    call read_pft_real(ncid, fname, 'z0mr'     , z0mr)
    call read_pft_real(ncid, fname, 'displar'  , displar)
    call read_pft_real(ncid, fname, 'dleaf'    , dleaf)
    call read_pft_real(ncid, fname, 'xl'       , xl)
    call read_pft_real(ncid, fname, 'rholvis'  , rhol(:,1))
    call read_pft_real(ncid, fname, 'rholnir'  , rhol(:,2))
    call read_pft_real(ncid, fname, 'rhosvis'  , rhos(:,1))
    call read_pft_real(ncid, fname, 'rhosnir'  , rhos(:,2))
    call read_pft_real(ncid, fname, 'taulvis'  , taul(:,1))
    call read_pft_real(ncid, fname, 'taulnir'  , taul(:,2))
    call read_pft_real(ncid, fname, 'tausvis'  , taus(:,1))
    call read_pft_real(ncid, fname, 'tausnir'  , taus(:,2))

    ! ---- Photosynthesis: the nine that vary by PFT ----
    call read_pft_real(ncid, fname, 'c3psn'    , c3psn)
    call read_pft_real(ncid, fname, 'leafcn'   , leafcn)
    call read_pft_real(ncid, fname, 'flnr'     , flnr)
    call read_pft_real(ncid, fname, 'fnitr'    , fnitr)
    call read_pft_real(ncid, fname, 'slatop'   , slatop)
    call read_pft_real(ncid, fname, 'qe'       , qe_ps)
    call read_pft_real(ncid, fname, 'theta_cj' , theta_cj)
    call read_pft_real(ncid, fname, 'bbbopt'   , bbbopt)
    call read_pft_real(ncid, fname, 'mbbopt'   , mbbopt)

    ! ---- Photosynthesis: the fourteen claimed uniform ----
    ! Read per PFT and CHECKED, not assumed. The C++ kernel takes these as
    ! scalars; if a parameter file ever varies one by PFT, that treatment
    ! becomes wrong for every PFT but the first, and silently so.
    block
      real(r8), allocatable :: tmp(:)
      integer :: k, ipft
      real(r8) :: lo, hi
      allocate(tmp(0:npft_param-1))
      do k = 1, n_photo_uniform
         call read_pft_real(ncid, fname, trim(photo_uniform_names(k)), tmp)
         lo = tmp(0); hi = tmp(0)
         do ipft = 0, npft_param-1
            lo = min(lo, tmp(ipft))
            hi = max(hi, tmp(ipft))
         end do
         if (hi /= lo) then
            call shr_sys_abort(subname//'ERROR: '// &
                 trim(photo_uniform_names(k))//' varies by PFT on '// &
                 trim(fname)//'; ELMxx carries it as a scalar and that is '// &
                 'no longer valid -- make it a per-patch view')
         end if
         photo_uniform(k) = tmp(0)
      end do
      deallocate(tmp)
    end block

    ! tc_stress is dimensioned allpfts = 1, so it reads as a length-1 array.
    status = pio_inq_varid(ncid, 'tc_stress', varid)
    if (status /= PIO_NOERR) then
       call shr_sys_abort(subname//'ERROR: no tc_stress on '//trim(fname))
    end if
    status = pio_get_var(ncid, varid, scalar)
    if (status /= PIO_NOERR) then
       call shr_sys_abort(subname//'ERROR: cannot read tc_stress from '//trim(fname))
    end if
    tc_stress = scalar(1)

    call pio_closefile(ncid)
    pftcon_read = .true.

    if (masterproc) then
       write(logunit,*) subname,'read ',npft_param,' PFTs from ',trim(fname)
       write(logunit,*) '    roota_par [1/m] ',minval(roota_par),' .. ',maxval(roota_par)
       write(logunit,*) '    rootb_par [1/m] ',minval(rootb_par),' .. ',maxval(rootb_par)
       write(logunit,*) '    smpso     [mm]  ',minval(smpso),' .. ',maxval(smpso)
       write(logunit,*) '    smpsc     [mm]  ',minval(smpsc),' .. ',maxval(smpsc)
       write(logunit,*) '    z0mr      [-]   ',minval(z0mr),' .. ',maxval(z0mr)
       write(logunit,*) '    displar   [-]   ',minval(displar),' .. ',maxval(displar)
       write(logunit,*) '    dleaf     [m]   ',minval(dleaf),' .. ',maxval(dleaf)
       write(logunit,*) '    xl        [-]   ',minval(xl),' .. ',maxval(xl)
       write(logunit,*) '    rhol vis  [-]   ',minval(rhol(:,1)),' .. ',maxval(rhol(:,1))
       write(logunit,*) '    taul vis  [-]   ',minval(taul(:,1)),' .. ',maxval(taul(:,1))
       write(logunit,*) '    tc_stress [C]   ',tc_stress
       call shr_sys_flush(logunit)
    end if

    ! Bounds that must hold for the btran formula to mean anything: the
    ! stress range smpso - smpsc is its denominator, so equal values would
    ! divide by zero, and closure must be drier (more negative) than opening.
    ! Checked over vegetated PFTs only -- bare ground is zero by design.
    if (any(smpso(1:npft_param-1) <= smpsc(1:npft_param-1))) then
       call shr_sys_abort(subname//'ERROR: smpso must exceed smpsc for vegetated PFTs')
    end if

    ! Vegetated PFTs must have a positive roughness ratio. Zero here is not a
    ! smooth canopy -- it drives ustar to zero and NaNs every canopy flux, and
    ! it is exactly what an unseeded view looks like. Aborting makes the
    ! failure name itself instead of surfacing as NaN four kernels later.
    if (any(z0mr(1:npft_param-1) <= 0.0_r8)) then
       call shr_sys_abort(subname//'ERROR: z0mr must be positive for vegetated PFTs')
    end if
    if (any(dleaf(1:npft_param-1) <= 0.0_r8)) then
       call shr_sys_abort(subname//'ERROR: dleaf must be positive for vegetated PFTs')
    end if

  end subroutine elmxx_read_pftcon

  !-----------------------------------------------------------------------
  subroutine read_pft_real(ncid, fname, varname, out)
    !
    ! Read one (pft) field whole. Unlike surfdata there is no decomposition
    ! here: every rank needs every PFT.
    !
    implicit none
    type(file_desc_t), intent(inout) :: ncid
    character(len=*) , intent(in)    :: fname, varname
    real(r8)         , intent(inout) :: out(0:)
    integer :: varid, status
    character(len=*), parameter :: subname = '(elmxx_read_pftcon) '

    status = pio_inq_varid(ncid, trim(varname), varid)
    if (status /= PIO_NOERR) then
       call shr_sys_abort(subname//'ERROR: no '//trim(varname)//' on '//trim(fname))
    end if
    status = pio_get_var(ncid, varid, out)
    if (status /= PIO_NOERR) then
       call shr_sys_abort(subname//'ERROR: cannot read '//trim(varname)//' from '//trim(fname))
    end if

  end subroutine read_pft_real

  !-----------------------------------------------------------------------
  subroutine elmxx_pftcon_clean()
    implicit none
    if (associated(roota_par)) deallocate(roota_par)
    if (associated(rootb_par)) deallocate(rootb_par)
    if (associated(smpso))     deallocate(smpso)
    if (associated(smpsc))     deallocate(smpsc)
    if (associated(z0mr))      deallocate(z0mr)
    if (associated(displar))   deallocate(displar)
    if (associated(dleaf))     deallocate(dleaf)
    if (associated(xl))        deallocate(xl)
    if (associated(rhol))      deallocate(rhol)
    if (associated(rhos))      deallocate(rhos)
    if (associated(taul))      deallocate(taul)
    if (associated(taus))      deallocate(taus)
    roota_par => null(); rootb_par => null()
    smpso => null(); smpsc => null()
    z0mr => null(); displar => null(); dleaf => null(); xl => null()
    rhol => null(); rhos => null(); taul => null(); taus => null()
    pftcon_read = .false.
  end subroutine elmxx_pftcon_clean

end module elmxxPftconMod
