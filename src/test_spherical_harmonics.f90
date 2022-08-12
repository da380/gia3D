program test_spherical_harmonics

  use module_constants
  use module_spherical_harmonics
  use module_fftw3
  implicit none

  integer(i4b) :: lmax,ith,iph,ncoef,ilm,l,m
  real(dp) :: th,ph,start,finish
  complex(dpc), dimension(:), allocatable :: ulm
  type(gauss_legendre_grid) :: grid
  type(scalar_gauss_legendre_grid) :: u
  type(C_PTR) :: plan_forward, plan_backward
  real(dp), dimension(:,:), allocatable :: xlm
  
  ! build the grid
  lmax = 512
  call u%build(lmax)

  ! make the FFT plans
  plan_forward  = u%plan(.true.)
  plan_backward = u%plan(.false.)
    

  ! set the values for the field
  do ith = 1,u%nth
     th = u%th(ith)
     do iph = 1,u%nph
        ph = u%ph(iph)
        u%val(iph,ith) = 0.25_dp*sqrt(5.0_dp/pi)*(3.0_dp*cos(th)**2-1.0_dp)
     end do
  end do

  ! precompute the Legendre polynomial values
  xlm = u%build_legendre()

  
  ! allocate the coefficient array
  ncoef = (lmax+1)**2
  allocate(ulm(ncoef))

  call cpu_time(start)
  call trans_scalar_gauss_legendre_grid(u,plan_forward,xlm,ulm)
  call cpu_time(finish)
  print *, finish-start
  
  

  
!  ilm = 0
!  do l = 0,lmax
!     ilm = ilm+1
!     print *, l,0,ulm(ilm)
!     do m = 1,l
!        ilm =ilm+1
!        print *, l,m,ulm(ilm)
!        ilm =ilm+1
!        print *, l,-m,ulm(ilm)
!     end do
!  end do
  
end program test_spherical_harmonics
