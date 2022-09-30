program test_mesh

  use module_constants
  use module_spherical_model
  use module_PREM
  use module_DECK
  use module_mesh
  use module_matrix
  use module_special_functions
  implicit none

  integer(i4b) :: ngll,nspec,inode,ispec,ilayer,nlayers,isection,i,count,ndim,l,io
  real(dp) :: drmax,rstart,Om
  type(spherical_model), allocatable :: model
  type(spherical_model_mesh) :: mesh
  type(boolean_array) :: ibool
  type(model_parameters) :: mp
  
  model = elastic_PREM(.false.)
!  model = anelastic_DECK('prem.200')
  l = 2
  drmax = 0.01_dp*model%r2/(l+1)
  ngll = 5

  Om = mp%Om()
  
  mesh = spherical_mesh(ngll,model,drmax,Om)

  l = 2
  rstart = spheroidal_start(mesh,l)
  ibool = build_boolean_spheroidal(mesh,rstart)
!  ibool = build_boolean_toroidal(mesh,rstart)  
  ndim = ibool%ndim

  open(newunit=io,file='mesh.out')
  count = 0  
  do isection = 1,mesh%nsections
     
       do ilayer = 1,mesh%section(isection)%nlayers

          associate(r   => mesh%section(isection)%layer(ilayer)%r, &
                    rho => mesh%section(isection)%layer(ilayer)%rho,    &
                   drho => mesh%section(isection)%layer(ilayer)%drho,  &
                      g => mesh%section(isection)%layer(ilayer)%g,     &
                     ep => mesh%section(isection)%layer(ilayer)%ep,    &
                  nspec => mesh%section(isection)%layer(ilayer)%nspec, &
                  ngll  => mesh%section(isection)%layer(ilayer)%ngll)

            do ispec = 1,nspec
               do inode = 1,ngll
            

                  write(io,*) r(inode,ispec)*length_norm,                  &
                              rho(inode,ispec)*density_norm,               &
                              drho(inode,ispec)*density_norm/length_norm,  & 
                              g(inode,ispec)*acceleration_norm,            &
                              ep(inode,ispec)
               end do

            end do
                        
          end associate

       end do
          
  end do
  close(io)

  
print *, ndim
  
  
end program test_mesh
