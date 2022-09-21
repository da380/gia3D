program test_matrix

  use module_constants
  use module_util
  use module_spherical_model
  use module_PREM
  use module_DECK
  use module_mesh
  use module_matrix
  use module_force
  implicit none

  logical :: spectop,ibnd,found
  integer(i4b) :: ngll,nspec,inode,ispec,ilayer,nlayers, &
                  isection,i,count,ndim,l,isection1,ilayer1, &
                  nsections,ispec1,ldab,kd,info,i1,i2,i3,k,j
  real(dp) :: drmax,rstart,tmp
  type(spherical_model), allocatable :: model
  type(spherical_model_mesh) :: mesh
  type(boolean_array) :: ibool
  type(radial_matrix) :: tormat,sphmat,sphmat2
  real(dp), dimension(:,:), allocatable :: b
  type(boolean_array) :: bool
  type(radial_matrices) :: mat

  call check_arguments(1,'-l [int > 0]')
  if(.not. found_command_argument('-l',l) .or. l < 0) stop 'l not set'
  
  model = elastic_PREM(.false.)
  drmax = 0.1_dp*model%r2/(l+1)
  ngll = 5
  mesh = spherical_mesh(ngll,model,drmax)



  mat = build_radial_matrices(mesh,l,spheroidals = .true. ,toroidals = .true.)

  stop
  
  !============================================================!
  !                       toroidal section                     !
  !============================================================!
  
  ! build the matrix
  tormat = build_toroidal_matrix(mesh,l)

  ! set the RHS
  ndim = tormat%ndim
  kd = tormat%kd
  ldab = tormat%ldab
  allocate(b(ndim,1))

  b = 0.0_dp
  if(l /= 1) then
     b(ndim,1) = 1.0_dp
  else
     b(ndim/2,1) = 1.0_dp
  end if

  call dpbtrs('U',ndim,kd,1,tormat%a,ldab,b,ndim,info)
  call error(info /= 0,'test_matrix','problem with toroidal substitution')  
  
  isection1 = tormat%ibool%isection1
  nsections = mesh%nsections

  open(99,file='tor.out')
  do isection = isection1,nsections
     ilayer1 = tormat%ibool%section(isection)%ilayer1
     nlayers = mesh%section(isection)%nlayers
     do ilayer = ilayer1,nlayers
        associate(layer => mesh%section(isection)%layer(ilayer), &
                  ibool => tormat%ibool%section(isection)%layer(ilayer))
          ispec1 = ibool%ispec1
          nspec  = layer%nspec
          ngll   = layer%ngll
          do ispec = ispec1,nspec
             do inode = 1,ngll

                i = ibool%get(1,inode,ispec)
                if(i == ndim .and. l == 1 ) then
                   write(99,*)  layer%r(inode,ispec)*length_norm,0.0_dp
                else
                   write(99,*)  layer%r(inode,ispec)*length_norm,b(i,1)
                end if
                
             end do             
          end do
        end associate
     end do
  end do
  close(99)

  
  !============================================================!
  !                      spheroidal section                    !
  !============================================================!
  
  sphmat = build_spheroidal_matrix(mesh,l)

  ! set the RHS
  ndim = sphmat%ndim  
  kd = sphmat%kd
  ldab = sphmat%ldab
  deallocate(b); allocate(b(ndim,1))
  b = 0.0_dp

!  isection = mesh%nsections
!  ilayer = mesh%section(isection)%nlayers
!  associate(layer => mesh%section(isection)%layer(ilayer),         &
!            ibool => sphmat%ibool%section(isection)%layer(ilayer))
!    ispec = layer%nspec
!    inode = layer%ngll
!    i = ibool%get(1,inode,ispec)
!    b(i,1) = layer%g(inode,ispec)
!    if(l > 1) then
!       i = ibool%get(3,inode,ispec)
!       b(i,1) = 1.0_dp
!    end if
  !  end associate

  call force_for_unit_harmonic_tide(mesh,sphmat%ibool,l,b)

  
  call dpbtrs('U',ndim,kd,1,sphmat%a,ldab,b,ndim,info)
  call error(info /= 0,'test_matrix','problem with spheroidal substitution')  


  
  isection1 = sphmat%ibool%isection1
  nsections = mesh%nsections  

  
  open(99,file='sph.out')
  do isection = isection1,nsections
     ilayer1 = sphmat%ibool%section(isection)%ilayer1
     nlayers = mesh%section(isection)%nlayers
     do ilayer = ilayer1,nlayers
        associate(layer => mesh%section(isection)%layer(ilayer), &
                  ibool => sphmat%ibool%section(isection)%layer(ilayer))
          ispec1 = ibool%ispec1
          nspec  = layer%nspec
          ngll   = layer%ngll
          do ispec = ispec1,nspec
             spectop = (l == 1) .and. &
                       (isection == mesh%nsections) .and.  & 
                       (ilayer == mesh%section(isection)%nlayers) .and. &
                       (ispec == nspec)
             do inode = 1,ngll

                
                ibnd = spectop .and. (inode == ngll)                              

                select type(layer)

                class is (spherical_solid_elastic_layer_mesh)

                   i1 = ibool%get(1,inode,ispec)
                   i2 = ibool%get(2,inode,ispec)
                   i3 = ibool%get(3,inode,ispec)
                   
                   if(ibnd) then
                      write(99,*)  layer%r(inode,ispec)*length_norm,b(i1,1),b(i2,1),0.0_dp
                   else
                      write(99,*)  layer%r(inode,ispec)*length_norm,b(i1,1),b(i2,1),b(i3,1)
                   end if

                class is (spherical_fluid_elastic_layer_mesh)

                   i = ibool%get(1,inode,ispec)

                   if(ibnd) then
                      write(99,*)  layer%r(inode,ispec)*length_norm,0.0_dp,0.0_dp,0.0_dp
                   else
                      write(99,*)  layer%r(inode,ispec)*length_norm,0.0_dp,0.0_dp,b(i,1)
                   end if
                   
                end select

                
                
             end do             
          end do
        end associate
     end do
  end do
  close(99)

  
  
end program test_matrix
