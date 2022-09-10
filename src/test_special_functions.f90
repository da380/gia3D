program test_special_funtions

  use module_constants
  use module_special_functions
  implicit none

  integer(i4b) :: n,ix,nx,l
  real(dp) :: alpha,beta,x,x1,x2,dx,y1,y2

  type(legendre_value) :: p
  type(wigner_value) :: d
  


  call p%init(0.2_dp,4)
  call d%init(0.2_dp,0,4)
  do l = 0,4

     call p%next()
     print *, p%get(0,l)
     call d%next()
     print *, d%get(0,0,0,l)
  end do
  
end program test_special_funtions

