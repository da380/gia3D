program test



  use module_constants
  use module_spherical_harmonics
  use module_quadrature
  implicit none

  character(len=256) :: string
  integer(i4b) :: l,m,mmax,lmax,nth,ith,i
  real(dp) :: th,th1,th2,dth,start,finish,f
  real(dp), dimension(:,:), allocatable :: pa,fa
  

  type(legendre_value) :: p

  lmax = 20
  mmax = 5
  th = 0.3_dp

  call p%init(th,mmax)  
  print *, p%deg(),p%get()
  do l = 1,lmax
     call p%next()     
     print *, p%deg(),p%get()
  end do




  
  
  
    
end program test
