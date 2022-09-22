program load_love_numbers
 
  use module_constants
  use module_util
  use module_spherical_model
  use module_PREM
  use module_DECK
  use module_mesh
  use module_matrix
  use module_force
  use module_linear_solver
  implicit none

  character(len=:), allocatable :: model_file,output_file
  integer(i4b) :: lmax,l,io
  type(spherical_model), allocatable :: model
  type(love_number), dimension(:), allocatable :: lln
 
  ! get the inputs
  call check_arguments(2,'-lmax [lmax], -f [output file]',1,'-m [model file]')
  if (.not. found_command_argument('-lmax',lmax)) stop 'lmax not set'
  if (.not. found_command_argument('-f',output_file)) stop 'output file not set'
  if(found_command_argument('-m',model_file)) then
     model = elastic_DECK(model_file)     
  else
     model = elastic_PREM(.false.)     
  end if

  allocate(lln(0:lmax))

  call make_love_numbers(model,lmax,lln = lln)
  
  ! loop over the degrees
  open(newunit = io,file=trim(output_file))
  do l = 1,lmax
            
     write(io,'(i6,4e20.8)') l,lln(l)%ku*length_norm, &
                               lln(l)%kv*length_norm, &
                               lln(l)%kw*length_norm, &          
                               lln(l)%kp*gravitational_potential_norm
            
  end do
  close(io)
  
end program load_love_numbers
