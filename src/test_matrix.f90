program test_matrix

  use module_constants
  use module_spherical_model
  use module_PREM
  use module_DECK
  use module_mesh
  use module_SEM_matrix
  use module_special_functions
  implicit none
 
  integer(i4b) :: ngll,nspec,inode,ispec,ilayer,nlayers, &
                  isection,i,count,ndim,l
  real(dp) :: drmax,rstart
  type(spherical_model), allocatable :: model
  type(spherical_model_mesh) :: mesh
  type(boolean_array) :: ibool
  type(radial_matrix) :: tormat
  
  model = elastic_PREM(.false.)
  drmax = model%r2/100.0_dp
  ngll = 5 
  mesh = spherical_mesh(ngll,model,drmax)

  l = 1
  tormat = build_toroidal_matrix(mesh,l)
  
  
  
end program test_matrix
