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


  subroutine sl_fingerprint(grid,lln,tln,ice1,ice2,sl1,sl2,model_parms,shorelines,rotation,verb)
    type(gauss_legendre_grid), intent(in) :: grid
    type(love_number), dimension(0:grid%lmax), intent(in) :: lln,tln
    real(dp), dimension(grid%nph,grid%nth), intent(in)  :: ice1
    real(dp), dimension(grid%nph,grid%nth), intent(in)  :: ice2
    real(dp), dimension(grid%nph,grid%nth), intent(in)  ::  sl1
    real(dp), dimension(grid%nph,grid%nth), intent(out) ::  sl2
    type(model_parameters), intent(in), optional :: model_parms
    logical, intent(in), optional :: shorelines,rotation
    logical, intent(in), optional :: verb

    integer(i4b), parameter :: maxit = 50
    real(dp), parameter :: eps_tol = 1.e-6_dp

    logical :: shorelines_move,write_error,rot
    integer(i4b) :: it,ilm
    real(dp) :: area,int,eps,g
    type(model_parameters) :: mp
    real(dp), dimension(grid%nph,grid%nph) :: sls,sigma
    real(dp), dimension(:,:), allocatable :: ofun
    complex(dpc) :: psi_20,psi_21
    complex(dpc), dimension(grid%ncoef_r) :: sigma_lm,u_lm,phi_lm
    
    ! deal with optional arguments
    if(present(model_parms)) mp = model_parms

    if(present(shorelines)) then
       shorelines_move = shorelines
    else
       shorelines_move = .true.
    end if

    if(present(verb)) then
       write_error = verb
    else
       write_error = .false.
    end if

    if(present(rotation)) then
       rot = rotation
    else
       rot = .true.
    end if

    ! extract required model parameters
    g = mp%g()
    
    
    ! precompute ocean function if needed
    if(.not.shorelines_move) then
       allocate(ofun(grid%nph,grid%nth))
       call ocean_function(sl1,ice1,ofun)
       area = grid%integrate(ofun)       
    end if

    ! set the initial values
    sl2 = sl1
    if(rot) then
       psi_20 = 0.0_dp
       psi_21 = 0.0_dp
    end if

    ! store the initial guess
    sls = sl2

    
    ! start the iterations
    do it = 1,maxit


       ! compute the water_and_ice_load
       if(shorelines_move) then
          call water_and_ice_load(sl1,ice1,sl2,ice2,sigma)
       else
          call water_and_ice_load(sl1,ice1,sl2,ice2,sigma,ofun)
       end if

       ! transform the load
       call grid%SH_trans(sigma,sigma_lm)

       ! impose mass conservation
       if(shorelines_move) area = ocean_area(grid,sl2,ice2)
       int = grid%integrate(sigma)
       sl2 = sl2 - int/(rho_water*area)
       ilm = grid%rindex(0,0)
       sigma_lm(ilm) = 0.0_dp

       ! estimate the error
       eps = maxval(abs(sl2-sls))/maxval(abs(sl2))

       if(write_error) then
          write(6,'(a,i2,a,e12.6)') "iteration = ",it,', relative error = ',eps
       end if
       if(eps < eps_tol) exit

       ! store the old sea level
       sls = sl2

       ! solve the loading problem spectrally
       call grid%filter(lln(:)%ku,sigma_lm,u_lm)
       call grid%filter(lln(:)%kp,sigma_lm,phi_lm)

       ! add in rotational contribution to deformation
       if(rot) then
          ilm = grid%rindex(2,0)
          u_lm(ilm)   = u_lm(ilm)   + tln(2)%ku*psi_20
          phi_lm(ilm) = phi_lm(ilm) + tln(2)%kp*psi_20
          ilm = grid%rindex(2,1)
          u_lm(ilm)   = u_lm(ilm)   + tln(2)%ku*psi_21
          phi_lm(ilm) = phi_lm(ilm) + tln(2)%kp*psi_21
          call centrifugal_potential_perturbation(mp,phi_lm,psi_20,psi_21)
       end if
       
       ! update sea level
       sigma_lm = -(u_lm + phi_lm/g)
       if(rot) then
          ilm = grid%rindex(2,0)
          sigma_lm(ilm) = sigma_lm(ilm) - psi_20/g
          ilm = grid%rindex(2,1)
          sigma_lm(ilm) = sigma_lm(ilm) - psi_21/g
       end if
       call grid%SH_itrans(sigma_lm,sigma)
       sl2 = sl1 + sigma
     
       if(it == maxit) print *, 'fingerprint: no convergence'
       
    end do

    
    return
  end subroutine sl_fingerprint


  !=================================================!
  !            sea level utility routines           !
  !=================================================!
  

  subroutine ocean_function(sl,ice,ofun)
    real(dp), dimension(:,:), intent(in) :: sl
    real(dp), dimension(:,:), intent(in) :: ice
    real(dp), dimension(:,:), intent(out) :: ofun
    where(rho_water*sl - rho_ice*ice > 0.0_dp)
       ofun = 1.0_dp
    elsewhere
       ofun = 0.0_dp
    end where    
    return
  end subroutine ocean_function

  subroutine ocean_function_mask(sl,ice,fun,comp)
    real(dp), dimension(:,:), intent(in)  :: sl
    real(dp), dimension(:,:), intent(in)  :: ice
    real(dp), dimension(:,:), intent(inout)  :: fun
    logical, intent(in), optional :: comp
    logical :: complement
    if(present(comp)) then
       complement = comp
    else
       complement  = .false.
    end if
    where(rho_water*sl - rho_ice*ice <= 0.0_dp .eqv. complement)
       fun = 0.0_dp
    end where    
    return
  end subroutine ocean_function_mask


  real(dp) function ocean_area(grid,sl,ice) result(area)
    type(gauss_legendre_grid), intent(in) :: grid
    real(dp), dimension(:,:), intent(in) :: sl
    real(dp), dimension(:,:), intent(in) :: ice
    integer(i4b) :: ith,iph,nth,nph
    real(dp) :: atmp,fac
    nth = grid%nth
    nph = grid%nph
    area = 0.0_dp
    do ith = 1,nth
       atmp = 0.0_dp
       do iph = 1,nph          
          if(rho_water*sl(iph,ith) - rho_ice*ice(iph,ith) > 0.0_dp) then
             atmp = atmp + 1.0_dp             
          end if          
       end do
       fac = grid%w(ith)*twopi/nph
       area = area + atmp*fac
    end do
    return
  end function ocean_area
  

  subroutine water_and_ice_load(sl1,ice1,sl2,ice2,sigma,ofun)
    real(dp), dimension(:,:), intent(in) :: sl1
    real(dp), dimension(:,:), intent(in) :: ice1
    real(dp), dimension(:,:), intent(in) :: sl2
    real(dp), dimension(:,:), intent(in) :: ice2
    real(dp), dimension(:,:), intent(out) :: sigma
    real(dp), dimension(:,:), intent(in), optional :: ofun
    
    if(present(ofun)) then
       where(ofun == 1.0_dp) 
          sigma = rho_water*(sl2-sl1)
       elsewhere
          sigma = rho_ice*(ice2-ice1)
       end where
    else
       where(rho_water*sl1 + rho_ice*ice1 > 0.0_dp)
          sigma = rho_water*sl1
       elsewhere
          sigma = rho_ice*ice1
       end where
       where(rho_water*sl2 + rho_ice*ice2 > 0.0_dp)
          sigma = rho_water*sl2-sigma
       elsewhere
          sigma = rho_ice*ice2-sigma
       end where
    end if
       
    return
  end subroutine water_and_ice_load
     

  subroutine centrifugal_potential_perturbation(mp,phi_lm,psi_20,psi_21)
    type(model_parameters), intent(in) :: mp
    complex(dpc), dimension(:), intent(in) :: phi_lm
    complex(dpc), intent(out) :: psi_20,psi_21

    integer(i4b) :: i20,i21
    real(dp) :: DI13,DI23,DI33,I1,I2,I3,b,Om,Dom1,Dom2,Dom3

    ! get necessary model parameters
    I1 = mp%I1()
    I2 = mp%I2()
    I3 = mp%I3()
    b  = mp%b()
    Om = mp%Om()

    ! indices for spherical harmonic arrays
    i20 = index_real_spherical_harmonic_expansion(2,0)
    i21 = index_real_spherical_harmonic_expansion(2,1)
    
    ! compute the inertia tensor perturbations
    DI33 =  sqrt(5.0_dp/(9.0_dp*pi))*(b**3/bigg)*real(phi_lm(i20))
    DI13 = -sqrt(5.0_dp/(6.0_dp*pi))*(b**3/bigg)*real(phi_lm(i21))
    DI23 =  sqrt(5.0_dp/(6.0_dp*pi))*(b**3/bigg)*imag(phi_lm(i21))    

    ! compute perturbed rotation vector components
    Dom1 =  Om*DI13/(I3-I1)
    Dom2 =  Om*DI23/(I3-I2)
    Dom3 = -Om*DI33/I3

    ! compute the centrifugal potential perturbation
    psi_20 = b*b*Om*sqrt((16.0_dp*pi)/45.0_dp)*Dom3
    psi_21 = b*b*Om*sqrt((2.0_dp*pi)/15.0_dp)*(-Dom1 + ii*Dom2)

    
    return
  end subroutine centrifugal_potential_perturbation
  
end module module_sea_level
