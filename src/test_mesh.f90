program test_mesh

  use module_constants
  use module_spherical_model
  use module_PREM
  use module_meshing
  use module_special_functions
  implicit none

  integer(i4b) :: ngll,nspec,inode,ispec,ilayer,nlayers,isection,i
  real(dp) :: drmax
  type(spherical_model), allocatable :: model
  type(spherical_model_mesh) :: mesh
  type(boolean_array) :: ibool
  
  
!  model = anelastic_PREM()
  model = maxwell_PREM(1000.0_dp,500.0_dp)
  drmax = model%r2/100.0_dp
  ngll = 5  
  mesh = spherical_mesh(ngll,model,drmax)
  ibool = build_boolean_array_scalar_field(mesh)


  
  do isection = 1,mesh%nsections

     do ilayer = 1,mesh%section(isection)%nlayers

        do ispec = 1,mesh%section(isection)%layer(ilayer)%nspec

           do inode = 1,mesh%section(isection)%layer(ilayer)%ngll

              i = ibool%get(1,inode,ispec,ilayer,isection)

              print *, i,ibool%nglob

              
           end do
           
        end do
        
        
     end do
     
  end do
  
end program test_mesh
