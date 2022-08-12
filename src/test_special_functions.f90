program test_special_funtions

  use module_constants
  use module_special_functions
  implicit none

  integer(i4b) :: n,ix,nx,l
  real(dp) :: alpha,beta,x,x1,x2,dx,y1,y2
  class(orthogonal_polynomial), allocatable :: poly1,poly2

  type(legendre_value) :: p
  type(wigner_value) :: d
  
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


  call p%init(0.2_dp,4)
  call d%init(0.2_dp,0,4)
  do l = 0,4

     call p%next()
     print *, p%get(0,l)
     call d%next()
     print *, d%get(0,0,0,l,norm=.true.)
  end do
  
end program test_special_funtions

