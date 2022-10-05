program test

  use module_constants
  use module_util
  use module_special_functions
  implicit none

  integer(i4b), parameter :: n = 1000001
  integer(i4b) :: i
  real(dp) :: r
  real(dp), dimension(n) :: rv

  call random_seed()
  call normal_random_variable(rv)
  
  open(99,file='test.out')
  do i = 1,n
!     call normal_random_variable(r)
     write(99,*) rv(i)
  end do
  close(99)
  
  

  
end program test
