program test_mesh

  use module_constants
  use module_spherical_model
  use module_PREM
  use module_DECK
  use module_meshing
  use module_special_functions
  implicit none

  integer(i4b) :: ngll,nspec,inode,ispec,ilayer,nlayers,isection,i
  real(dp) :: drmax
  type(spherical_model), allocatable :: model
  type(spherical_model_mesh) :: mesh
  type(boolean_array) :: ibool
  
  
  model = elastic_PREM()
  drmax = model%r2/100.0_dp
  ngll = 5 
  mesh = spherical_mesh(ngll,model,drmax)

  
  
end program test_mesh
