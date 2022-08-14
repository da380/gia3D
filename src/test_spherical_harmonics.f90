program test_spherical_harmonics

  use module_constants
  use module_spherical_harmonics
  use module_fftw3
  implicit none

  integer(i4b) :: lmax,nmax,n,l,m
  real(dp), dimension(:,:), allocatable :: y
  type(gauss_legendre_grid) :: grid
  type(scalar_spherical_harmonic_expansion) :: u,v,w

  real(dp) :: start,finish
  
  lmax = 512
  call u%allocate(lmax)
  u%data = 1.0_dp
  v = ii*u
  w = 2.0_dp*u+v*ii
  call u%saxpy(2.0_dp,v)
  
  
  
end program test_spherical_harmonics




