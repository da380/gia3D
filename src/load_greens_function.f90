program load_greens_function


  use module_constants
  use module_util
  use module_special_functions
  use module_spherical_model
  use module_PREM
  use module_DECK
!  use module_mesh
!  use module_matrix
!  use module_force
  use module_linear_solver
  implicit none
  


  character(len=:), allocatable :: model_file,output_file
  integer(i4b) :: lmax,l,io,nth,ith
  real(dp) :: th,th1,th2,dth,ku,kp,pl,fac
  type(spherical_model), allocatable :: model
  type(love_number), dimension(:), allocatable :: ln
  type(legendre_value) :: p

  if (.not. found_command_argument('-f',output_file)) stop 'output file not set'  
  if (.not. found_command_argument('-lmax',lmax)) stop 'lmax not set'
  allocate(ln(1:lmax)) 
  if(found_command_argument('-m',model_file)) then
     model = elastic_DECK(model_file)     
  else
     model = elastic_PREM(.false.)     
  end if
  call make_love_numbers(model,1,lmax,lln = ln)


  th1 = 0.0_dp
  th2 = pi
  nth = 180
  dth = (th2-th1)/(nth-1)

  open(newunit = io,file =trim(output_file))
  do ith = 1,nth
     th = th1 + (ith-1)*dth
     ku = 0.0_dp
     kp = 0.0_dp
     call p%init(th,0)
     call p%next()
     do l = 1,lmax
        fac = exp(-10.0_dp*l**2/lmax**2)
        print *, fac
        call p%next()
        pl = p%get(0)*fac
        ku  = ku + ln(l)%ku*pl
        kp  = kp + ln(l)%kp*pl
     end do
     write(io,*) th,ku,kp
  end do
  close(io)
  
  


     
end program load_greens_function
