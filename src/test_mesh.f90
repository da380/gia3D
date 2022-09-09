program test_mesh

  use module_constants
  use module_spherical_model
  use module_PREM
  use module_DECK
  use module_meshing
  use module_special_functions
  implicit none

  integer(i4b) :: ngll,nspec,inode,ispec,ilayer,nlayers,isection,i,count,ndim
  real(dp) :: drmax,rstart
  type(spherical_model), allocatable :: model
  type(spherical_model_mesh) :: mesh
  type(boolean_array) :: ibool
  
  
  model = elastic_PREM(.false.)
  drmax = model%r2/100.0_dp
  ngll = 5 
  mesh = spherical_mesh(ngll,model,drmax)

  rstart = 0.6*mesh%r2
  !  ibool = build_boolean_spheroidal(mesh,rstart)
  ibool = build_boolean_toroidal(mesh,rstart)  
  ndim = ibool%ndim

  
  count = 0
  do isection = ibool%isection1,mesh%nsections

     
     
     associate(ibool => ibool%section(isection), &
               section => mesh%section(isection))

       do ilayer = ibool%ilayer1,section%nlayers

          associate(ibool => ibool%layer(ilayer), &
                    layer => section%layer(ilayer))

            ngll = layer%ngll
            nspec = layer%nspec

            do ispec = ibool%ispec1,nspec

               do inode = 1,ngll

                  select type(layer)

                  class is(spherical_solid_elastic_layer_mesh)
                     
!                     print *, isection,ilayer,ispec,inode,ibool%get(1,inode,ispec), &
!                                                          ibool%get(2,inode,ispec), &
!                                                          ibool%get(3,inode,ispec)

                     print *, isection,ilayer,ispec,inode,ibool%get(1,inode,ispec)
                     
                  class is(spherical_fluid_elastic_layer_mesh)

                     print *, isection,ilayer,ispec,inode,ibool%get(1,inode,ispec)

                     
                     
                  end select
                  

               end do
               count = count -1
               
            end do
            
            
            
          end associate

       print *, '--------------------------------------------------------------------'
          
       end do


       
     end associate

     print *, '===================================================================='
     
  end do

  print *, ndim
  
  
end program test_mesh
