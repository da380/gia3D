program test_mesh

  use module_constants
  use module_spherical_model
  use module_PREM
  use module_DECK
  use module_meshing
  use module_special_functions
  implicit none

  integer(i4b) :: ngll,nspec,inode,ispec,ilayer,nlayers,isection,i,count
  real(dp) :: drmax,rstart
  type(spherical_model), allocatable :: model
  type(spherical_model_mesh) :: mesh
  type(boolean_array) :: ibool
  
  
  model = elastic_PREM(.false.)
  drmax = model%r2/100.0_dp
  ngll = 5 
  mesh = spherical_mesh(ngll,model,drmax)

  rstart = 0.6*mesh%r2
  ibool = build_boolean_toroidal(mesh,rstart)

  
  count = 0
  do isection = ibool%isection1,ibool%nsections

     
     
     associate(ibool => ibool%section(isection))

       do ilayer = ibool%ilayer1,ibool%nlayers

          associate(ibool => ibool%layer(ilayer), &
                    mesh => mesh%section(isection)%layer(ilayer))

            ngll = mesh%ngll
            nspec = mesh%nspec

            do ispec = ibool%ispec1,nspec

               do inode = 1,ngll
                  count = count + 1
                  print *, isection,ilayer,ispec,inode,ibool%get(1,inode,ispec)-count
               end do
               count = count -1
               
            end do
            
            
            
          end associate
          
       end do

       
     end associate

     
  end do
  
  
end program test_mesh
