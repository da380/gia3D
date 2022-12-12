program love_numbers
  
  use module_constants
  use module_util
  use module_spherical_model
  use module_PREM
  use module_DECK
  use module_linear_solver
  implicit none

  character(len=:), allocatable :: model_file,output_file
  integer(i4b) :: lmax,l,io
  type(spherical_model), allocatable :: model
  type(love_number), dimension(:), allocatable :: ln

  
  if (.not. found_command_argument('-lmax',lmax)) stop 'lmax not set'
  allocate(ln(1:lmax))
  
  if(found_command_argument('-m',model_file)) then
     model = elastic_DECK(model_file)     
  else
     model = elastic_PREM(.false.)     
  end if

  if (found_command_argument('-lln',output_file)) then
     call make_love_numbers(model,1,lmax,lln = ln)
     open(newunit = io,file=trim(output_file))
     do l = 1,lmax
        write(io,'(i6,3e20.8)') l,ln(l)%ku*length_norm/load_norm, &
                                  ln(l)%kv*length_norm/load_norm, &
                                  ln(l)%kp*gravitational_potential_norm/load_norm
     end do
     close(io)
  end if


  if (found_command_argument('-tln',output_file)) then
     call make_love_numbers(model,1,lmax,tln = ln)
     open(newunit = io,file=trim(output_file))
     do l = 1,lmax
        write(io,'(i6,3e20.8)') l,ln(l)%ku*length_norm/load_norm, &
                                  ln(l)%kv*length_norm/load_norm, &
                                  ln(l)%kp*gravitational_potential_norm/load_norm
     end do
     close(io)
  end if

  
     
end program love_numbers
