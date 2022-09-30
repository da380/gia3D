program sea_level_fingerprint
  
  use module_constants
  use module_physical_constants
  use module_util
  use module_spherical_harmonics
  use module_spherical_model
  use module_PREM
  use module_DECK
  use module_linear_solver
  use module_sea_level
  implicit none

  integer(i4b), parameter :: nreq = 2,nopt = 4
  logical, dimension(nreq) :: dreq
  logical, dimension(nopt) :: dopt
  character(len=256), dimension(nreq) :: ireq
  character(len=256), dimension(nopt) :: iopt
  
  logical :: verb,rotation,shorelines
  character(len=:), allocatable :: model_file,output_file
  integer(i4b) :: lmax,ith,iph,io,l,m,it,nph,nth
  real(dp) :: th,ph,start,finish

  type(spherical_model), allocatable :: model
  type(gauss_legendre_grid) :: grid
  type(love_number), dimension(:), allocatable :: lln,tln
  real(dp), dimension(:,:), allocatable :: ice1,sl1,ice2,sl2, &
                                           sigma,sls,ofun
  

  ! deal with inputs
  ireq(1) = '-l [int] -- maximum degree'
  dreq(1) = .true.
  ireq(2) = '-f [str] -- output file name'
  dreq(2) = .true.
  iopt(1) = '-m [str] -- model file (optional, default = PREM)'
  dopt(1) = .true.
  iopt(2) = '-r       -- rotation on (optional)'
  dopt(2) = .false.
  iopt(3) = '-s       -- shoreline migration on (optional)'
  dopt(3) = .false.
  iopt(4) = '-v       -- verbose (optional)'
  dopt(4) = .false.
  call print_argument_info(nreq,ireq,dreq,nopt,iopt,dopt)   
  if (.not. found_command_argument('-lmax',lmax)) stop 'lmax not set'
  if (.not.found_command_argument('-f',output_file)) stop 'output file not set'    
  if(found_command_argument('-m',model_file)) then
     model = elastic_DECK(model_file)     
  else
     model = elastic_PREM(.false.)     
  end if  
  rotation = found_command_argument('-r')
  shorelines = found_command_argument('-s')
  verb = found_command_argument('-v')
  
  
  ! set up the GLL mesh
  call grid%build(lmax,0)
  nph = grid%nph
  nth = grid%nth

  ! compute the Love numbers
  allocate(lln(0:lmax),tln(0:lmax))
  call make_love_numbers(model,0,lmax,lln = lln,tln = tln)

  
  ! allocate the arrays
  allocate(sl1(nph,nth), ice1(nph,nth), &
           sl2(nph,nth), ice2(nph,nth),ofun(nph,nth))

  ! set values for the input fields
  do ith = 1,nth
     th = grid%th(ith)
     do iph = 1,nph
        ph = grid%ph(iph)
        sl1(iph,ith)  = initial_sea_level(ph,th)
        ice1(iph,ith) = initial_ice(ph,th)
        ice2(iph,ith) = 0.95_dp*ice1(iph,ith)
     end do
  end do

  call sl_fingerprint(grid,lln,tln,ice1,ice2,sl1,sl2,verb = verb,            &
                                                     rotation = rotation,    &
                                                     shorelines = shorelines)

  ! write out the new sea level
  sl2 = sl2-sl1
  open(newunit = io,file=output_file)
  write(io,*) grid%nth,grid%nph,0
  do ith = 1,grid%nth
     do iph = 1,grid%nph
        write(io,*) grid%ph(iph),grid%th(ith),sl2(iph,ith)*length_norm

     end do
  end do
  
  close(io)
  
  
contains


  real(dp) function initial_sea_level(ph,th) result(sl)
    real(dp), intent(in) :: ph
    real(dp), intent(in) :: th

    real(dp), parameter :: th1 = pio4
    real(dp), parameter :: ph1 = pio2
    real(dp), parameter :: w11 = 0.15_dp*pi
    real(dp), parameter :: w12 = 0.25_dp*pi
    real(dp), parameter :: t1 = 1000.0_dp

    real(dp), parameter :: th2 = pio2
    real(dp), parameter :: ph2 = 3.0_dp*pi/2.0_dp
    real(dp), parameter :: w21 = 0.2_dp*pi
    real(dp), parameter :: w22 = 0.3_dp*pi
    real(dp), parameter :: t2 = 500.0_dp

    real(dp), parameter :: sl0 = 3000.0_dp
    
    real(dp) :: del1,del2
    
    del1 = acos(cos(th)*cos(th1) + sin(th)*sin(th1)*cos(ph-ph1))
    del2 = acos(cos(th)*cos(th2) + sin(th)*sin(th2)*cos(ph-ph2))

    if(del1 <= w11) then
       sl = -t1
    else if(del1 > w11 .and. del1 <= w12) then
       sl = sl0 - (t1+sl0)*((w12-del1)/(w12-w11))**2
    else if(del2 <= w21) then
       sl = -t2
    else if(del2 > w21 .and. del2 <= w22) then
       sl = sl0 - (t2+sl0)*((w22-del2)/(w22-w21))**2
    else
       sl = sl0
    end if
       
    sl = sl/length_norm

    
    return
  end function initial_sea_level



  real(dp) function initial_ice(ph,th) result(ice)
    real(dp), intent(in) :: ph
    real(dp), intent(in) :: th

    real(dp), parameter :: th1 = pio4
    real(dp), parameter :: ph1 = pio2
    real(dp), parameter :: w11 = 0.1_dp*pi
    real(dp), parameter :: w12 = 0.2_dp*pi
    real(dp), parameter :: ice1 = 1000.0_dp
    
    real(dp) :: del1
    
    del1 = acos(cos(th)*cos(th1) + sin(th)*sin(th1)*cos(ph-ph1))


    if(del1 <= w11) then
       ice = ice1
    else if(del1 > w11 .and. del1 <= w12) then
       ice =  ice1*((w12-del1)/(w12-w11))**2
    else
       ice = 0.0_dp
    end if
       
    ice = ice/length_norm

    
    return
  end function initial_ice


  
  
end program sea_level_fingerprint
