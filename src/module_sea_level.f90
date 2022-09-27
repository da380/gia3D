module module_sea_level

  use module_constants
  use module_physical_constants
  use module_spherical_model
  use module_spherical_harmonics
  use module_linear_solver
  implicit none

contains


  !=================================================================!
  !                      elastic fingerprint code                   !
  !=================================================================!


  subroutine sl_fingerprint(grid,lln,tln,ice1,ice2,sl1,sl2,model_parms,fix_shorelines,verb)
    type(gauss_legendre_grid), intent(in) :: grid
    type(love_number), dimension(0:grid%lmax), intent(in) :: lln
    type(love_number), dimension(2:2), intent(in) :: tln
    type(real_scalar_gauss_legendre_field), intent(in)  :: ice1
    type(real_scalar_gauss_legendre_field), intent(in)  ::  sl1
    type(real_scalar_gauss_legendre_field), intent(in)  :: ice2
    type(real_scalar_gauss_legendre_field), intent(inout) ::  sl2
    type(model_parameters), intent(in), optional :: model_parms
    logical, intent(in), optional :: fix_shorelines
    logical, intent(in), optional :: verb

    integer(i4b), parameter :: maxit = 50
    real(dp), parameter :: eps_tol = 1.e-6_dp

    logical :: fixed,write_error
    integer(i4b) :: it
    real(dp) :: area,int,eps,g,b,Om,psi_20,psi_21,phi_20,phi_21
    type(model_parameters) :: mp
    type(real_scalar_gauss_legendre_field) :: sls,ofun,sigma
    type(real_scalar_spherical_harmonic_expansion) :: sigma_lm

    ! deal with optional arguments
    if(present(model_parms)) mp = model_parms

    if(present(fix_shorelines)) then
       fixed = fix_shorelines
    else
       fixed = .true.
    end if

    if(present(verb)) then
       write_error = verb
    else
       write_error = .false.
    end if

    ! extract required model parameters
    g = mp%g()
    
    ! allocate the local variables
    call sls%allocate(grid)
    call sigma%allocate(grid)
    call sigma_lm%allocate(grid)
    
    ! precompute ocean function if needed
    if(fixed) then
       call ofun%allocate(grid)
       call ocean_function(sl1,ice1,ofun)
       area = ofun%integrate(grid)
    end if

    ! set the initial values
    sl2 = sl1
    psi_20 = 0.0_dp
    psi_21 = 0.0_dp

    ! store the initial guess
    sls = sl2
    
    ! start the iterations
    do it = 1,maxit


       ! compute the water_and_ice_load
       if(fixed) then
          call water_and_ice_load(sl1,ice1,sl2,ice2,sigma,ofun)
       else
          call water_and_ice_load(sl1,ice1,sl2,ice2,sigma)
       end if

       ! transform the load
       call grid%SH_trans(sigma,sigma_lm)       

       ! impose mass conservation
       if(.not.fixed) area = ocean_area(grid,sl2,ice2)
       int = sigma%integrate(grid)
       sl2%rdata = sl2%rdata - int/(rho_water*area)

       ! estimate the error
       eps = maxval(abs(sl2%rdata-sls%rdata))/maxval(abs(sl2%rdata))

       if(write_error) then
          write(6,'(a,i2,a,e12.6)') "iteration = ",it,', relative error = ',eps
       end if
       if(eps < eps_tol) exit

       ! store the old sea level
       sls%rdata = sl2%rdata

       ! compute the rotational perturbation
       psi_20 = 0.0_dp
       psi_21 = 0.0_dp
       
       ! solve the loading problem spectrally
       call sigma_lm%filter(-lln(:)%ku - lln%kp/g)

       ! add in rotationall contributions
       call sigma_lm%set(2,0,sigma_lm%get(2,0)-psi_20/g)
       call sigma_lm%set(2,1,sigma_lm%get(2,1)-psi_20/g)	

       ! update sea level spatially
       call grid%SH_itrans(sigma_lm,sigma)
       sl2%rdata =  sl1%rdata + sigma%rdata
     
       if(it == maxit) stop 'fingerprint: no convergence'
       
    end do

    call sls%delete()
    call sigma%delete()
    call sigma_lm%delete()
    if(fixed) call ofun%delete()
    
    return
  end subroutine sl_fingerprint




  

  subroutine ocean_function(sl,ice,ofun)
    type(real_scalar_gauss_legendre_field), intent(in)  :: sl
    type(real_scalar_gauss_legendre_field), intent(in)  :: ice
    type(real_scalar_gauss_legendre_field), intent(inout)  :: ofun
    where(rho_water*sl%rdata - rho_ice*ice%rdata > 0.0_dp)
       ofun%rdata = 1.0_dp
    elsewhere
       ofun%rdata = 0.0_dp
    end where    
    return
  end subroutine ocean_function

  subroutine ocean_function_mask(sl,ice,fun,comp)
    type(real_scalar_gauss_legendre_field), intent(in)  :: sl
    type(real_scalar_gauss_legendre_field), intent(in)  :: ice
    type(real_scalar_gauss_legendre_field), intent(inout)  :: fun
    logical, intent(in), optional :: comp
    logical :: complement,test
    if(present(comp)) then
       complement = comp
    else
       complement  = .false.
    end if
    where(rho_water*sl%rdata - rho_ice*ice%rdata <= 0.0_dp .eqv. complement)
       fun%rdata = 0.0_dp
    end where    
    return
  end subroutine ocean_function_mask


  real(dp) function ocean_area(grid,sl,ice) result(area)
    type(gauss_legendre_grid), intent(in) :: grid
    type(real_scalar_gauss_legendre_field), intent(in)  :: sl
    type(real_scalar_gauss_legendre_field), intent(in)  :: ice
    integer(i4b) :: ith,iph,nth,nph
    real(dp) :: atmp,fac
    nth = grid%nth
    nph = grid%nph
    area = 0.0_dp
    do ith = 1,nth
       atmp = 0.0_dp
       do iph = 1,nph          
          if(rho_water*sl%get(iph,ith) - rho_ice*ice%get(iph,ith) > 0.0_dp) then
             atmp = atmp + 1.0_dp             
          end if          
       end do
       fac = grid%w(ith)*twopi/nph
       area = area + atmp*fac
    end do
    return
  end function ocean_area
  
  subroutine water_and_ice_load(sl1,ice1,sl2,ice2,sigma,ofun)
    type(real_scalar_gauss_legendre_field), intent(in)  :: sl1
    type(real_scalar_gauss_legendre_field), intent(in)  :: ice1
    type(real_scalar_gauss_legendre_field), intent(in)  :: sl2
    type(real_scalar_gauss_legendre_field), intent(in)  :: ice2
    type(real_scalar_gauss_legendre_field), intent(inout)  :: sigma
    type(real_scalar_gauss_legendre_field), intent(in), optional  :: ofun

    logical :: ofp
    integer(i4b) :: ith,iph
    real(dp) :: sig,f1,f2,c

    ofp = present(ofun)    
    do ith = 1,sl1%nth
       do iph = 1,sl1%nph          
          f1  = rho_water*sl1%get(iph,ith)
          f2  = rho_ice*ice1%get(iph,ith)
          if(ofp) then
             c = ofun%get(iph,ith)
             sig = c*f1 + (1.0_dp-c)*f2
          else
             if(f1 > f2) then
                sig = f1
             else
                sig = f2
             end if
          end if
          f1 = rho_water*sl2%get(iph,ith)
          f2 = rho_ice*ice2%get(iph,ith)
          if(ofp) then
             sig = c*f1 + (1.0_dp-c)*f2 - sig
          else
             if(f1 > f2) then
                sig = f1-sig
             else
                sig = f2-sig
             end if
          end if
          call sigma%set(iph,ith,sig)          
       end do
    end do
    

       
    return
  end subroutine water_and_ice_load
     
  
  
end module module_sea_level
