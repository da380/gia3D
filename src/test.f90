program test



  use module_constants
  use module_spherical_harmonics
  use module_quadrature
  implicit none

  
  integer(i4b) :: l,m,n,lmax,mmax,nmax
  real(dp) :: beta

  type(wigner_value) :: p

  type(legendre_value) :: q
  
  lmax = 1000
  nmax = 0
  mmax = 5
  beta = 0.2_dp
  call p%init(beta,nmax,mmax)
!  call p%init(beta,mmax)
  print *, p%v(:)!*sifourpi
  do l = 1,lmax
     call p%next()
     print *, p%v(:)!*sqrt((2*l+1.0_dp)/fourpi)
  end do

  print *, '====================================='

  call q%init(beta,mmax)
  print *, q%v(:)
  do l = 1,lmax
     call q%next()
     print *, q%v(:)
  end do
  
!  
!  call p%init(beta,nmax,mmax)

!  do m = 0,mmax
!     do n = -min(m,nmax),min(m,nmax)
!        print *, m,n,p%ind(n,m)
!     end do
!  end do
  




  
  
  
    
end program test
