program test_spherical_model
 
  use module_constants
  use module_spherical_model
  use module_util
  implicit none

  logical :: ok
  character(len=:), allocatable :: line
  integer(i4b) :: ilayer,ir,isection,nr,io
  real(dp) :: r,r1,r2,dr,tau_LM,tau_UM,tau_C

  type(spherical_model) :: model


  
!  model = elastic_PREM(.false.)
!  call model%write_elastic('prem.out',200)

  

!  model = elastic_DECK('prem.out')
!  call model%write_elastic('prem2.out',200)

  tau_LM =  10000.0_dp*yr2sec/time_norm
  tau_UM =   1000.0_dp*yr2sec/time_norm
  tau_C  = 5000000.0_dp*yr2sec/time_norm
  model = maxwell_PREM(tau_LM,tau_UM)
  call model%write_maxwell('prem_maxwell.out',200)

  model = maxwell_DECK('prem_maxwell.out')
  
end program test_spherical_model


