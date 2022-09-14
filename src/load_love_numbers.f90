program load_love_numbers
 
  use module_constants
  use module_util
  use module_spherical_model
  use module_PREM
  use module_DECK
  use module_mesh
  use module_matrix
  implicit none

  character(len=:), allocatable :: model_file,output_file
  integer(i4b) :: lmax,ngll,l,ndim,kd,ldab,isection, &
                  ilayer,ispec,inode,i,info,io
  real(dp) :: drmax,u,v,p,g
  real(dp), dimension(:,:), allocatable :: b
  type(spherical_model), allocatable :: model
  type(spherical_model_mesh) :: mesh
  type(radial_matrix) :: tormat,sphmat

  ! get the inputs
  call check_arguments(2,'-lmax [lmax], -f [output file]',1,'-m [model file]')
  if (.not. found_command_argument('-lmax',lmax)) stop 'lmax not set'
  if (.not. found_command_argument('-f',output_file)) stop 'output file not set'
  if(found_command_argument('-m',model_file)) then
     model = elastic_DECK(model_file)     
  else
     model = elastic_PREM(.false.)     
  end if
  
  ! build the mesh
  drmax = 0.1_dp*model%r2/(lmax+1)
  ngll = 5
  mesh = spherical_mesh(ngll,model,drmax)

  ! loop over the degrees
  open(newunit = io,file=trim(output_file))
  do l = 1,lmax

     ! build the matrix
     sphmat = build_spheroidal_matrix(mesh,l)

     ! set the RHS
     ndim = sphmat%ndim
     kd = sphmat%kd
     ldab = sphmat%ldab
     if(allocated(b)) deallocate(b); allocate(b(ndim,1))
     b = 0.0_dp

     isection = mesh%nsections
     ilayer = mesh%section(isection)%nlayers
     associate(layer => mesh%section(isection)%layer(ilayer),         &
               ibool => sphmat%ibool%section(isection)%layer(ilayer))
       ispec = layer%nspec
       inode = layer%ngll
       i = ibool%get(1,inode,ispec)
       g = layer%g(inode,ispec)
       b(i,1) = g
       if(l > 1) then
          i = ibool%get(3,inode,ispec)
          b(i,1) = 1.0_dp
       end if
       b = b/load_norm
       

       ! solve the linear system
       call dpbtrs('U',ndim,kd,1,sphmat%a,ldab,b,ndim,info)
       call error(info /= 0,'test_matrix','problem with spheroidal substitution')  

       i = ibool%get(1,inode,ispec)
       u = b(i,1)
       i = ibool%get(2,inode,ispec)
       v = b(i,1)
       if(l > 1) then
          i = ibool%get(3,inode,ispec)
          p = b(i,1)/g
       else
          p = 0.0_dp
       end if

       write(io,'(i6,3e20.8)') l,u*length_norm, &
                                 v*length_norm, &
                                 p*length_norm
       
     end associate

     
  end do
  close(io)
  
end program load_love_numbers
