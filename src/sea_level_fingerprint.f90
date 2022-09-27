program sea_level_fingerprint
  
  use module_constants
  use module_physical_constants
  use module_util
  use module_spherical_harmonics
  use module_spherical_model
  use module_PREM
  use module_linear_solver
  use module_sea_level
  implicit none
  
  logical :: found
  integer(i4b) :: lmax,ith,iph,io,l,m,it
  real(dp) :: th,ph,f,th1,th2,g,int,area,fac,eps,start,finish
  complex(dpc) :: ctmp

  type(spherical_model), allocatable :: model
  type(gauss_legendre_grid) :: grid
  type(love_number), dimension(:), allocatable :: lln,tln
  type(real_scalar_gauss_legendre_field) :: ice1,sl1,ice2,sl2,sigma,sls,ofun
  type(real_scalar_spherical_harmonic_expansion) :: sigma_lm
  
  ! set up the GL grid
  call check_arguments(1,'-lmax [maximum degree]')
  found = found_command_argument('-lmax',lmax)
  call grid%allocate(lmax,0)

  ! compute the love numbers
  model = elastic_PREM(.false.)
  
  allocate(lln(0:lmax))
  call make_love_numbers(model,0,lmax,lln = lln)
  allocate(tln(2:2))
  call make_love_numbers(model,2,2,tln = tln)


  ! allocate the GL fields
  call ice1%allocate(grid)
  call  sl1%allocate(grid)
  call ice2%allocate(grid)
  call  sl2%allocate(grid)

  ! set values for the input fields
  do ith = 1,grid%nth
     th = grid%th(ith)
     do iph = 1,grid%nph
        ph = grid%ph(iph)
        f = initial_sea_level(ph,th)
        call sl1%set(iph,ith,f)
        f = initial_ice(ph,th)
        call ice1%set(iph,ith,f)
        call ice2%set(iph,ith,0.0_dp)
     end do
  end do

  call cpu_time(start)
  call sl_fingerprint(grid,lln,tln,ice1,ice2,sl1,sl2)
  call cpu_time(finish)
  print *, finish-start
  
  ! write out the new sea level
  open(newunit = io,file='sl.out')
  write(io,*) grid%nth,grid%nph,0
  do ith = 1,grid%nth
     do iph = 1,grid%nph
        write(io,*) grid%ph(iph),grid%th(ith),(sl1%get(iph,ith)-0*sl1%get(iph,ith))*length_norm
     end do
  end do
  
  close(io)
  
  
contains


  real(dp) function initial_sea_level(ph,th) result(sl)
    real(dp), intent(in) :: ph
    real(dp), intent(in) :: th

    real(dp), parameter :: amp1  = 0.3_dp
    real(dp), parameter :: amp2  = 0.1_dp
    real(dp), parameter :: th1  = 20.0_dp*deg2rad
    real(dp), parameter :: th2  = 30.0_dp*deg2rad
    real(dp), parameter :: th3  = 140.0_dp*deg2rad
    real(dp), parameter :: th4  = 160.0_dp*deg2rad
    real(dp), parameter :: sl1  = -100.0_dp/length_norm
    real(dp), parameter :: sl2  =  3000.0_dp/length_norm
    real(dp), parameter :: sl3  = -5.0_dp/length_norm

    real(dp) :: th11,th22,th33,th44

    th11 = th1*(1.0_dp+amp1*sin(4.0_dp*ph))
    th22 = th2*(1.0_dp+amp1*sin(4.0_dp*ph))
    th33 = th3*(1.0_dp+amp2*cos(5.0_dp*ph))
    th44 = th4*(1.0_dp+amp2*cos(5.0_dp*ph))
    
    if(th <= th11) then
       sl = sl1
    else if(th > th11 .and. th < th22) then
       sl = sl1 + (th-th11)*(sl2-sl1)/(th22-th11)
    else if(th >= th22 .and. th < th33) then
       sl = sl2
    else if(th >= th33 .and. th < th44) then
       sl = sl2 + (th-th33)*(sl3-sl2)/(th44-th33)
    else
       sl = sl3
    end if
    
    return
  end function initial_sea_level

  
  real(dp) function initial_ice(ph,th) result(ice)
    real(dp), intent(in) :: ph
    real(dp), intent(in) :: th

    real(dp), parameter :: th1   = 18.0_dp*deg2rad
    real(dp), parameter :: th2   = 20.0_dp*deg2rad
    real(dp), parameter :: ice1  = 1000.0_dp/length_norm
    real(dp), parameter :: ice2  =    0.0_dp/length_norm 


    if(th <= th1) then
       ice = ice1
    else if(th > th1 .and. th < th2) then
       ice = ice1 + (ice2-ice1)*(th-th1)/(th2-th1)
    else
       ice = ice2
    end if
    
    return
  end function initial_ice
  
end program sea_level_fingerprint
