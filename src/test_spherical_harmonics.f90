program test_spherical_harmonics

  use module_constants
  use module_spherical_harmonics
  use module_fftw3
  implicit none

  integer(i4b) :: lmax,nmax,n,l
  type(gauss_legendre_grid) :: grid

  
  lmax = 256
  nmax = 2
  call grid%build(lmax,nmax)



end program test_spherical_harmonics
