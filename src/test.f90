program test

  use module_constants
  use module_util
  use module_special_functions
  implicit none

  integer(i4b) :: i,n = 10
  real(dp) :: x
  real(dp), dimension(:), allocatable :: xx 

  allocate(xx(0))

  do i = 1,n
     x = 2*i
     xx = [xx,x]
  end do

  xx = [xx,2*xx]
  
  print *, xx
  

  
end program test
