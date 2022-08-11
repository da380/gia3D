program test_special_funtions

  use module_constants
  use module_special_functions
  implicit none

  integer(i4b) :: n,ix,nx
  real(dp) :: alpha,beta,x,x1,x2,dx,y1,y2
  class(orthogonal_polynomial), allocatable :: poly1,poly2

  ! set the polynomial type
  poly1 = legendre()

  alpha = 1.0_dp
  beta  = 0._dp
  poly2 = jacobi(alpha,beta)

  nx = 1000
  x1 = -1.0_dp
  x2 =  1.0_dp
  dx = (x2-x1)/(nx-1)

  n = 5

  open(99,file='test.out')
  do ix = 1,nx
     x = x1+(ix-1)*dx
     y1 = poly1%eval(n,x)
     y2 = poly2%eval(n,x)
     write(99,*) x,y1,y2,y1-y2
  end do
  close(99)
  
  
end program test_special_funtions
