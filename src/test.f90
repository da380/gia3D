program test



  use module_constants
  use module_spherical_harmonics
  use module_quadrature
  implicit none

  
  integer(i4b) :: l,m,n,lmax,mmax,nmax
  real(dp) :: th
  type(legendre_value) :: p

  lmax = 5
  mmax = 5  
  th = 0.3_dp
  call p%init(th,mmax)
  print *, p%deg(),p%get()
  do l = 1,lmax
     call p%next()
     print *, p%deg(),p%get()     
  end do
  
  
!  type(wigner_value) :: p
  


!  l = 5
!  mmax = 3
!  nmax = l
!  call p%allocate(mmax,nmax)
!  do n = 0,nmax
!     do m = -min(mmax,n),min(n,mmax)
!        print *, n,m, p%ind(m,n)
!     end do
!  end do




  
  
  
    
end program test
