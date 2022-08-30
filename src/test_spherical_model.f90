program test_spherical_model
 
  use module_constants
  use module_spherical_model
  implicit none

  integer(i4b) :: ilayer,ir,isection,nr,io
  real(dp) :: r,r1,r2,dr,eta_LM,eta_UM

  type(spherical_model) :: model


  model = elastic_PREM()





  
  
  
end program test_spherical_model


