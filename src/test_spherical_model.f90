program test_spherical_model
 
  use module_constants
  use module_spherical_model
  use module_PREM
  use module_DECK
  use module_util
  implicit none

  logical :: ok
  character(len=:), allocatable :: line
  integer(i4b) :: ilayer,ir,isection,nr,io
  real(dp) :: r,r1,r2,dr,tau_LM,tau_UM,tau_C

  type(spherical_section) :: in,upper,lower
  type(spherical_model) :: model

  
  model = elastic_PREM(.false.)
  call model%write_elastic('prem.out',200,isotropic=.true.)
  model = elastic_DECK('prem.out')
  call model%write_elastic('prem2.out',200)


  

  
end program test_spherical_model


