program test_model


  use module_constants
  use module_spherical_harmonics
  use module_quadrature
  implicit none

  character(len=256) :: string
  integer(i4b) :: l,m,mmax,lmax,nth,ith,i
  real(dp) :: th,th1,th2,dth,start,finish,f
  real(dp), dimension(:,:), allocatable :: pa,fa
  

  type(legendre_value) :: p

  lmax = 512
  mmax = lmax
  th = 0.3_dp

  call cpu_time(start)
  call p%init(th,mmax)
  do l = 1,lmax
     call p%next()     
  end do
  call cpu_time(finish)
  print *, finish-start


  
end program test_model

