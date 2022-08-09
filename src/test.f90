program test



  use module_constants
  use module_spherical_harmonics
  use module_quadrature
  implicit none

  
  integer(i4b) :: l,m,n,lmax,mmax,nmax
  real(dp) :: beta,start,finish

  type(wigner_value) :: p

  type(legendre_value) :: q
  
  lmax = 1000
  nmax = 0
  mmax = lmax
  beta = 0.3_dp
  call cpu_time(start)
  call p%init(beta,nmax,mmax)
!  print *, p%v(p%ind(0,0))
  do l = 1,lmax
     call p%next()
!     print *, p%v(p%ind(0,l))
  end do
  call cpu_time(finish)
  print *, finish-start

  call cpu_time(start)
  call q%init(beta,mmax)
!  print *, p%v(p%ind(0,0))
  do l = 1,lmax
     call q%next()
!     print *, p%v(p%ind(0,l))
  end do
  call cpu_time(finish)
  print *, finish-start


  




  
  
  
    
end program test
