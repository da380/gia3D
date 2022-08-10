program test

  use module_constants
  use module_spherical_harmonics
  use module_quadrature
  implicit none

  
  integer(i4b) :: l,lp,m,n,lmax,mmax,nmax,i
  real(dp) :: beta,start,finish,time1,time2,v
  real(dp), dimension(:,:), allocatable :: d
  
  type(wigner_value) :: p
  type(legendre_value) :: q
  
  lmax = 3
  nmax = lmax
  mmax = lmax
  beta = pio4
  call p%init(beta,nmax,mmax)  
  do l = 0,lmax-1
     call p%next()
  end do



  allocate(d(2*nmax+1,2*mmax+1))
  d = p%get(-nmax,nmax,-mmax,mmax)

  
  do n = 1,2*nmax+1
     print *, d(n,:)
  end do

  

  
 
!  call q%init(beta,mmax)
!  print *, q%v(:)
!  do l = 1,lmax
!     call q%next()
!     print *, q%v(:)
!  end do




  




  
  
  
    
end program test
