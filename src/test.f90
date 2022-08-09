program test



  use module_constants
  use module_spherical_harmonics
  use module_quadrature
  implicit none

  
  integer(i4b) :: l,lp,m,n,lmax,mmax,nmax,i
  real(dp) :: beta,start,finish,time1,time2
  real(dp), dimension(:,:), allocatable :: d
  
  type(wigner_value) :: p

  type(legendre_value) :: q
  
  l = 2
  nmax = 1
  mmax = l
  beta = 0.200_dp
  call p%init(beta,nmax,mmax)
  do lp = 0,l-1
     call p%next()
  end do


  allocate(d(-nmax:nmax,-l:l))
  d = 0.0_dp
  do m = 0,l
     do n = -min(nmax,m),min(nmax,m)
        i = p%ind(n,m)
        d(n,m) = p%v(i)
     end do
  end do

  do n = -nmax,nmax
     print *, d(n,:)
  end do

  

  
 
!  call q%init(beta,mmax)
!  print *, q%v(:)
!  do l = 1,lmax
!     call q%next()
!     print *, q%v(:)
!  end do




  




  
  
  
    
end program test
