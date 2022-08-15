program test_spherical_harmonics

  use module_constants
  use module_special_functions
  use module_spherical_harmonics
  use module_fftw3
  implicit none

  integer(i4b) :: lmax,nmax,n,l,m,ith,iph,ilm
  real(dp) :: start,finish,th,ph
  type(gauss_legendre_grid) :: grid
  type(scalar_gauss_legendre_field) :: u
  type(scalar_spherical_harmonic_expansion) :: ulm
  type(wigner_array) :: d

  ! make the GL-grid
  lmax = 512
  nmax = 0
  
  call grid%allocate(lmax,nmax)

  ! allocate a scalar field
  call u%allocate(grid)


 
  ! set the values for the field
  do ith = 1,grid%nth()
     th = grid%th(ith)
     do iph = 1,grid%nph()
        ph = grid%ph(iph)
        u%data(u%index(iph,ith)) = 0.25_dp*sqrt(5.0_dp/pi)*(3.0_dp*cos(th)**2-1.0_dp) + sifourpi
     end do
  end do
  
  
  ! allocate the expansion
  call ulm%allocate(lmax)
  
  call cpu_time(start)
  call grid%SH_trans(u,ulm)
  call cpu_time(finish)
  print *, finish-start

!  ilm = 0
!  do l = 0,lmax
!     ilm = ilm+1
!     print *, l,0,ulm%data(ilm)
!     do m = 1,l
!        ilm = ilm+1
!        print *, l,m,ulm%data(ilm)
!        ilm = ilm+1
!        print *, l,-m,ulm%data(ilm)
!     end do
!  end do
  
  
end program test_spherical_harmonics




