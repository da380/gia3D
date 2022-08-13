program test_spherical_harmonics

  use module_constants
  use module_spherical_harmonics
  use module_fftw3
  implicit none

  integer(i4b) :: lmax,nmax,n,l,m
  real(dp), dimension(:,:), allocatable :: y
  type(gauss_legendre_grid) :: grid
  type(spherical_harmonic_coefficient) :: u
  
  lmax = 128
  nmax = 2
  call grid%build(lmax,nmax,precomp=.false.)

  call u%allocate(lmax,nmax)



  
end program test_spherical_harmonics
