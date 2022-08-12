program test_wigner

  use module_constants
  use module_special_functions

  integer(i4b) :: l,lmax,nmax,mmax
  real(dp) :: beta
  type(legendre_value) :: plm
  type(wigner_value) :: dlm

  beta = 0.2_dp
  lmax = 1000
  nmax = 0
  mmax = 5
  call dlm%init(beta,nmax,mmax)
  do l = 0,lmax
     call dlm%next()
     print *, l,dlm%get(0,0,0,mmax)
  end do

  
end program test_wigner
