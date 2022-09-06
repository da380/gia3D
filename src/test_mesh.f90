program test_mesh

  use module_constants
  use module_spherical_model
  use module_PREM
  use module_meshing
  use module_special_functions
  implicit none

  integer(i4b) :: ngll,nspec,inode,ispec,ilayer,nlayers,isection
  real(dp) :: drmax
  type(spherical_model), allocatable :: model
  type(spherical_model_mesh) :: mesh
  class(spherical_layer_mesh), allocatable :: mesh_tmp

  !  model = maxwell_PREM(1000.0_dp,500.0_dp)
  model = anelastic_PREM()
  drmax = model%r2/100.0_dp
  ngll = 5  
  mesh = make_spherical_mesh(ngll,model,drmax)

  associate(mesh => mesh%section(3)%layer(2))

    select type(mesh)

    class is(spherical_solid_anelastic_layer_mesh)
       print *, mesh%qA
       
    end select
    
  end associate
  
end program test_mesh
