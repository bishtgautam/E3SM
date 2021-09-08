module ScalarVarianceMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! calculate scalar variance via HOM and HET methods
  !
  ! !USES:
  use shr_kind_mod         , only : r8 => shr_kind_r8
  use shr_log_mod          , only : errMsg => shr_log_errMsg
  use abortutils           , only : endrun
  use elm_varcon           , only : spval, nameg
  use elm_varctl           , only : iulog
  use decompMod            , only : bounds_type
  use subgridAveMod        , only : build_scale_l2g_gpu
  use domainMod            , only : ldomain
  use GridcellType         , only : grc_pp
  use LandunitType         , only : lun_pp
  use ColumnDataType       , only : col_pp
  use VegetationDataType   , only : veg_pp
  use column_varcon        , only : icol_roof, icol_sunwall, icol_shadewall
  use column_varcon        , only : icol_road_perv , icol_road_imperv
  !
  ! !PUBLIC TYPES:
  implicit none
  private
  save
  !
  ! !PUBLIC MEMBER FUNCTIONS:

  public :: calculate_scalar_covaiance_het
  public :: calculate_scalar_covaiance_hom
  public :: calculate_scalar_covaiance_pft_hom

  integer, parameter :: unity = 0, urbanf = 1, urbans = 2
  integer, parameter :: natveg = 3, veg =4, ice=5, nonurb=6, lake=7
  !------------------------------------------------------------------------

contains

!!! HET method
subroutine calculate_scalar_covaiance_het(bounds, &
         eflx_sh_tot, &
         eflx_lh_tot, &
         q_ref2m, &
         t_ref2m, &
         taux_patch, &
         tauy_patch, &
         forc_pbot_not_downscaled_grc, &
         forc_rho_not_downscaled_grc, &
         q_ref2m_grc, &
         t_ref2m_grc, &
         thlp2_het_grc, &
         rtp2_het_grc, &
         rtpthlp_het_grc, &
         p2c_scale_type, c2l_scale_type, l2g_scale_type)
    ! !DESCRIPTION:
    ! Perfrom subgrid-average from pfts to gridcells.
    ! Averaging is only done for points that are not equal to "spval".
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds
	!!!
    real(r8), intent(in)  :: eflx_sh_tot( bounds%begp: )  ! input pft array
	real(r8), intent(in)  :: eflx_lh_tot( bounds%begp: )  ! input pft array
	real(r8), intent(in)  :: q_ref2m( bounds%begp: )  ! input pft array
	real(r8), intent(in)  :: t_ref2m( bounds%begp: )  ! input pft array
	!real(r8), intent(in)  :: fv_patch( bounds%begp: )  ! input pft array
    real(r8), intent(in)  :: taux_patch( bounds%begp: )  ! input pft array
    real(r8), intent(in)  :: tauy_patch( bounds%begp: )  ! input pft array
    
	real(r8), intent(in)  :: forc_pbot_not_downscaled_grc( bounds%begg: )  ! input pft array
	real(r8), intent(in)  :: forc_rho_not_downscaled_grc( bounds%begg: )  ! input pft array
	real(r8), intent(in)  :: q_ref2m_grc( bounds%begg: )  ! input pft array
	real(r8), intent(in)  :: t_ref2m_grc( bounds%begg: )  ! input pft array
    
    real(r8), intent(out) :: thlp2_het_grc( bounds%begg: )  ! output gridcell array
    real(r8), intent(out) :: rtp2_het_grc( bounds%begg: )  ! output gridcell array
    real(r8), intent(out) :: rtpthlp_het_grc( bounds%begg: )  ! output gridcell array
    
    !character(len=*), intent(in) :: p2c_scale_type ! scale factor type for averaging
    !character(len=*), intent(in) :: c2l_scale_type ! scale factor type for averaging
    !character(len=*), intent(in) :: l2g_scale_type ! scale factor type for averaging
    integer , intent(in) :: p2c_scale_type ! scale factor type for averaging
    integer , intent(in) :: c2l_scale_type ! scale factor type for averaging |||  unity = 0, urbanf = 1, urbans = 2
    integer , intent(in) :: l2g_scale_type ! unity =0, natveg = 3, veg =4, ice=5, nonurb=6, lake=7
    !
    !  !LOCAL VARIABLES:
    integer  :: p,c,l,g,index                   ! indices
    logical  :: found                              ! temporary for error check
    real(r8) :: scale_p2c(bounds%begp:bounds%endp) ! scale factor
    real(r8) :: scale_c2l(bounds%begc:bounds%endc) ! scale factor
    real(r8) :: scale_l2g(bounds%begl:bounds%endl) ! scale factor
    real(r8) :: sumwt(bounds%begg:bounds%endg)     ! sum of weights
    real(r8) :: tsa_pft,pressure,lh_pft,sh_pft,taux_pft,tauy_pft,rho_pft,t_pft,q_pft,twa,qwa,thlp2_pft, rtp2_pft, rtpthlp_pft,p_weight
    !------------------------------------------------------------------------

    ! Enforce expected array sizes

    !call build_scale_l2g(bounds, l2g_scale_type, &
    !     scale_l2g)

    !call create_scale_c2l(bounds,c2l_scale_type,scale_c2l)

    !if (p2c_scale_type == 'unity') then
     !  do p = bounds%begp,bounds%endp
     !     scale_p2c(p) = 1.0_r8
    !   end do
    !else
     !  write(iulog,*)'p2g_1d error: scale type ',c2l_scale_type,' not supported'
     !  call endrun(msg=errMsg(__FILE__, __LINE__))
    !end if


    call build_scale_l2g_gpu(bounds, l2g_scale_type, &
         scale_l2g)

    if (c2l_scale_type == unity) then
       do c = bounds%begc,bounds%endc
          scale_c2l(c) = 1.0_r8
       end do
    else if (c2l_scale_type == urbanf) then
       do c = bounds%begc,bounds%endc
          l = col_pp%landunit(c)
          if (lun_pp%urbpoi(l)) then
             if (col_pp%itype(c) == icol_sunwall) then
                scale_c2l(c) = 3.0 * lun_pp%canyon_hwr(l)
             else if (col_pp%itype(c) == icol_shadewall) then
                scale_c2l(c) = 3.0 * lun_pp%canyon_hwr(l)
             else if (col_pp%itype(c) == icol_road_perv .or. col_pp%itype(c) == icol_road_imperv) then
                scale_c2l(c) = 3.0_r8
             else if (col_pp%itype(c) == icol_roof) then
                scale_c2l(c) = 1.0_r8
             end if
          else
             scale_c2l(c) = 1.0_r8
          end if
       end do
    else if (c2l_scale_type == urbans) then
       do c = bounds%begc,bounds%endc
          l = col_pp%landunit(c)
          if (lun_pp%urbpoi(l)) then
             if (col_pp%itype(c) == icol_sunwall) then
                scale_c2l(c) = (3.0 * lun_pp%canyon_hwr(l)) / (2.*lun_pp%canyon_hwr(l) + 1.)
             else if (col_pp%itype(c) == icol_shadewall) then
                scale_c2l(c) = (3.0 * lun_pp%canyon_hwr(l)) / (2.*lun_pp%canyon_hwr(l) + 1.)
             else if (col_pp%itype(c) == icol_road_perv .or. col_pp%itype(c) == icol_road_imperv) then
                scale_c2l(c) = 3.0 / (2.*lun_pp%canyon_hwr(l) + 1.)
             else if (col_pp%itype(c) == icol_roof) then
                scale_c2l(c) = 1.0_r8
             end if
          else
             scale_c2l(c) = 1.0_r8
          end if
       end do
    else
       print *, 'p2g_1d error: scale type ',c2l_scale_type,' not supported'
    end if

    if (p2c_scale_type == unity) then
       do p = bounds%begp,bounds%endp
          scale_p2c(p) = 1.0_r8
       end do
    else
       print *, 'p2g_1d error: scale type ',c2l_scale_type,' not supported'
    end if
    
    !!! start
    thlp2_het_grc(bounds%begg : bounds%endg) = spval
    rtp2_het_grc(bounds%begg : bounds%endg) = spval
    rtpthlp_het_grc(bounds%begg : bounds%endg) = spval
    sumwt(bounds%begg : bounds%endg) = 0._r8
    
    do p = bounds%begp,bounds%endp
       if (veg_pp%active(p) .and. veg_pp%wtgcell(p) /= 0._r8) then
          c = veg_pp%column(p)
          l = veg_pp%landunit(p)
          g = veg_pp%gridcell(p)
          if (ldomain%frac(g) > 0.999999_r8)  then ! only calculate pure land grid
          if (eflx_sh_tot(p) /= spval .and. eflx_sh_tot(p) /= spval .and. q_ref2m(p) /= spval .and. t_ref2m(p) /= spval &
          .and. taux_patch(p) /= spval .and. tauy_patch(p) /= spval .and. forc_pbot_not_downscaled_grc(g) /= spval .and. forc_rho_not_downscaled_grc(g) /= spval) then
             if (sumwt(g) == 0._r8) then
               thlp2_het_grc(g) = 0._r8
               rtp2_het_grc(g) = 0._r8
               rtpthlp_het_grc(g) = 0._r8
             end if
                lh_pft = eflx_lh_tot(p);
                sh_pft = eflx_sh_tot(p);
                rho_pft = forc_rho_not_downscaled_grc(g);
                !uf_pft = fv_patch(p);  !! need to be updated use parameterizaed uf
                taux_pft = taux_patch(p); 
                tauy_pft = tauy_patch(p); 
                tsa_pft = t_ref2m(p);
                pressure =forc_pbot_not_downscaled_grc(g);
                t_pft = tsa_pft*(100000/pressure)**0.286;
                q_pft = q_ref2m(p);
                
                twa = t_ref2m_grc(g)*(100000/pressure)**0.286;
                qwa = q_ref2m_grc(g);
                
                p_weight = scale_p2c(p) * scale_c2l(c) * scale_l2g(l) * veg_pp%wtgcell(p);
                
                call calculate_scalar_covaiance_pft_hom(lh_pft, sh_pft, rho_pft, taux_pft,tauy_pft,thlp2_pft, rtp2_pft, rtpthlp_pft)

                !wp2_het = nansum(p_weight*wp2_pft, 1);
                thlp2_het_grc(g) = thlp2_het_grc(g) + p_weight*thlp2_pft + p_weight*((t_pft-twa)**2);
                rtp2_het_grc(g) = rtp2_het_grc(g) + p_weight*rtp2_pft + p_weight*(q_pft-qwa)**2;
                rtpthlp_het_grc(g) = rtpthlp_het_grc(g) + p_weight*rtpthlp_pft + p_weight*(t_pft-twa)*(q_pft-qwa);

                sumwt(g) = sumwt(g) + veg_pp%wtgcell(p)
          end if
          end if
       end if
    end do
    found = .false.
    do g = bounds%begg, bounds%endg
       if (sumwt(g) > 1.0_r8 + 1.e-6_r8) then
          found = .true.
          index = g
       else if (sumwt(g) /= 0._r8) then
          thlp2_het_grc(g) = thlp2_het_grc(g)/sumwt(g)
          rtp2_het_grc(g) = rtp2_het_grc(g)/sumwt(g)
          rtpthlp_het_grc(g) = rtpthlp_het_grc(g)/sumwt(g)
       end if
    end do
    
    if (found) then
       write(iulog,*)'p2g_1d error: sumwt is greater than 1.0 at g= ',index
       call endrun(decomp_index=index, elmlevel=nameg, msg=errMsg(__FILE__, __LINE__))
    end if

  end subroutine calculate_scalar_covaiance_het
  
  
  
  
!!! HOM method  
  subroutine calculate_scalar_covaiance_hom(bounds, &
         eflx_sh_tot_grc , &
         eflx_lh_tot_grc, &
         !fv_grc, &
         taux_grc, &
         tauy_grc, &
         forc_rho_not_downscaled_grc, &
         thlp2_hom_grc, &
         rtp2_hom_grc, &
         rtpthlp_hom_grc)
    ! !DESCRIPTION:
    ! Perfrom subgrid-average from pfts to gridcells.
    ! Averaging is only done for points that are not equal to "spval".
    !
    ! !ARGUMENTS:
      
    type(bounds_type), intent(in) :: bounds
	!!!
    real(r8), intent(in)  :: eflx_sh_tot_grc( bounds%begg: )  ! input gridcell array
	real(r8), intent(in)  :: eflx_lh_tot_grc( bounds%begg: )  ! input gridcell array
	!real(r8), intent(in)  :: fv_grc( bounds%begg: )  ! input gridcell array
    real(r8), intent(in)  :: taux_grc( bounds%begg: )  ! input gridcell array
    real(r8), intent(in)  :: tauy_grc( bounds%begg: )  ! input gridcell array
	real(r8), intent(in)  :: forc_rho_not_downscaled_grc( bounds%begg: )  ! input gridcell array
	
    real(r8), intent(out) :: thlp2_hom_grc( bounds%begg: )  ! output gridcell array
    real(r8), intent(out) :: rtp2_hom_grc( bounds%begg: )  ! output gridcell array
    real(r8), intent(out) :: rtpthlp_hom_grc( bounds%begg: )  ! output gridcell array
    
    !
    !  !LOCAL VARIABLES:
    integer  :: g,index                   ! indices
    real(r8), parameter :: cpair = 1004.64_r8; !heat capacity dry air at const pres (J/kg/k)
    real(r8), parameter :: lv = 2.501e6 !latent heat of evaporation (J/Kg)
    real(r8), parameter :: thl_tol=1.e-2 !(K)
    real(r8), parameter :: rt_tol=1.e-8 !(kg/kg)
    real(r8), parameter :: w_tol=2.e-2 !(m/s)
    real(r8), parameter :: a_const=1.8_r8
    real(r8), parameter :: t0=300._r8
    real(r8), parameter :: grav=9.80616_r8
    real(r8), parameter :: z_const=1._r8
    real(r8), parameter :: ufmin=0.01_r8

    real(r8) :: lh, sh, uf, rho, ustar2, wstar
    real(r8) :: upwp, vpwp, wprtp, wpthlp, wp2_sfc


    !------------------------------------------------------------------------

    ! Enforce expected array sizes
    thlp2_hom_grc(bounds%begg : bounds%endg) = spval
    rtp2_hom_grc(bounds%begg : bounds%endg) = spval
    rtpthlp_hom_grc(bounds%begg : bounds%endg) = spval
    
    do g = bounds%begg,bounds%endg
             
             
        if (ldomain%frac(g) > 0.999999_r8)  then ! only calculate pure land grid
          if (eflx_sh_tot_grc(g) /= spval .and. eflx_lh_tot_grc(g) /= spval .and. taux_grc(g) /= spval  .and. tauy_grc(g) /= spval .and.  forc_rho_not_downscaled_grc(g) /= spval) then
            lh = eflx_lh_tot_grc(g);
            sh = eflx_sh_tot_grc(g);
            rho = forc_rho_not_downscaled_grc(g);
            !uf = fv_grc(g); !! need to be updated use parameterizaed uf
            upwp = taux_grc(g)/rho;
            vpwp = tauy_grc(g)/rho;
            
            wprtp = lh/(lv*rho);
            wpthlp= sh/(cpair*rho);
            
            ustar2 = sqrt(upwp**2 + vpwp**2);
            
            if (wpthlp > 0) then 
                 wstar = (1/t0 * grav * wpthlp * z_const) ** (1._r8/3._r8)
            else
                wstar = 0;
            endif
            
            uf = sqrt(ustar2 + 0.3 * wstar**2)
            
            if (uf <= ufmin) then
                uf = ufmin;
            endif
            
            wp2_sfc = a_const * uf**2;

            thlp2_hom_grc(g)=0.4_r8 * a_const * (wpthlp/uf)**2;
            if (thlp2_hom_grc(g) < thl_tol**2) then
                thlp2_hom_grc(g) = thl_tol**2;
            end if
             
            rtp2_hom_grc(g)=0.4_r8 * a_const * (wprtp/uf)**2;
            if (rtp2_hom_grc(g) < rt_tol**2) then
                rtp2_hom_grc(g) = rt_tol**2;
            end if
            rtpthlp_hom_grc(g)=0.2_r8*a_const*(wpthlp/uf)*(wprtp/uf);
        end if
       end if
    end do
   

  end subroutine calculate_scalar_covaiance_hom
  
  
!!! HOM method at pft level
  subroutine calculate_scalar_covaiance_pft_hom(lh_pft, sh_pft, rho_pft, taux_pft, tauy_pft, &
         thlp2_pft, &
         rtp2_pft, &
         rtpthlp_pft)
    ! !DESCRIPTION:
    ! Perfrom subgrid-average from pfts to gridcells.
    ! Averaging is only done for points that are not equal to "spval".
    !
    ! !ARGUMENTS:
   
	!!!
    real(r8), intent(in)  :: lh_pft ! input pft array
	real(r8), intent(in)  :: sh_pft  ! input pft array
	real(r8), intent(in)  :: rho_pft  ! input pft array
	!real(r8), intent(in)  :: uf_pft  ! input pft array // !! need to be updated use parameterizaed uf
	real(r8), intent(in)  :: taux_pft 
    real(r8), intent(in)  :: tauy_pft 
    
    real(r8), intent(out) :: thlp2_pft  ! output gridcell array
    real(r8), intent(out) :: rtp2_pft  ! output gridcell array
    real(r8), intent(out) :: rtpthlp_pft ! output gridcell array
    
    !
    !  !LOCAL VARIABLES:
    real(r8), parameter :: cpair = 1004.64_r8 !heat capacity dry air at const pres (J/kg/k)
    real(r8), parameter :: lv = 2.501e6 !latent heat of evaporation (J/Kg)
    real(r8), parameter :: thl_tol=1.e-2 !(K)
    real(r8), parameter :: rt_tol=1.e-8 !(kg/kg)
    real(r8), parameter :: w_tol=2.e-2 !(m/s)
    real(r8), parameter :: a_const=1.8_r8
    real(r8), parameter :: t0=300._r8
    real(r8), parameter :: grav=9.80616_r8
    real(r8), parameter :: z_const=1._r8
    real(r8), parameter :: ufmin=0.01_r8
    
    real(r8) :: upwp, vpwp, wprtp, wpthlp, wp2_sfc, uf_pft, ustar2_pft, wstar_pft, ustar2, wstar, rho

    !------------------------------------------------------------------------

    ! Enforce expected array sizes
    thlp2_pft= spval
    rtp2_pft = spval
    rtpthlp_pft = spval
    
   

      if (lh_pft /= spval .and. sh_pft /= spval .and. rho_pft /= spval  .and. taux_pft /= spval  .and. tauy_pft /= spval) then
             
             upwp = taux_pft/rho_pft;
             vpwp = tauy_pft/rho_pft;
            
             wprtp = lh_pft/(lv*rho_pft);
             wpthlp= sh_pft/(cpair*rho_pft);
             
             ustar2 = sqrt(upwp**2 + vpwp**2);
            
            if (wpthlp > 0) then 
                 wstar = (1/t0 * grav * wpthlp * z_const) ** (1._r8/3._r8)
            else
                wstar = 0;
            endif
            
            uf_pft = sqrt(ustar2 + 0.3 * wstar**2)
            
            if (uf_pft <= ufmin) then
                uf_pft = ufmin;
            endif
            
             wp2_sfc = a_const * uf_pft**2;

             thlp2_pft=0.4_r8 * a_const * (wpthlp/uf_pft)**2;
             if (thlp2_pft < thl_tol**2) then
                thlp2_pft = thl_tol**2;
             end if
             
             rtp2_pft=0.4_r8 * a_const * (wprtp/uf_pft)**2;
             if (rtp2_pft < rt_tol**2) then
                rtp2_pft = rt_tol**2;
             end if
             rtpthlp_pft=0.2_r8 * a_const*(wpthlp/uf_pft)*(wprtp/uf_pft);

       end if
 
   

  end subroutine calculate_scalar_covaiance_pft_hom
  
end module ScalarVarianceMod
