program test_spherical_harmonics

  use module_constants
  use module_special_functions
  use module_spherical_harmonics
  use module_fftw3
  implicit none
  
  integer(i4b) :: lmax,nmax,n,l,m,ith,iph,ilm,io
  real(dp) :: start,finish,th,ph
  real(dp), dimension(:,:), allocatable :: u
  complex(dpc), dimension(:,:), allocatable :: w
  complex(dpc), dimension(:), allocatable :: ulm,vlm,wlm,xlm
  type(gauss_legendre_grid) :: grid
  
  
  ! make the grid  
  lmax = 256
  nmax = 2
  print *, 'building the grid for degree ',lmax
  call cpu_time(start)
  call grid%allocate(lmax,nmax)  
  call cpu_time(finish)
  print *, 'assembly time = ',finish-start


  !==================================!
  !     test real transformations    !
  !==================================!  
  
  ! allocate the spatial array
  allocate(u(grid%nph,grid%nth))

  ! allocate the coefficient array
  allocate(ulm(grid%ncoef_r))

  ! set the coefficients
  ulm = 0.0_dp
  ilm = grid%rindex(2,2)
  ulm(ilm) = 1.0_dp

  ! form the spatial field
  call grid%SH_itrans(ulm,u)

  ! allocate a new coefficient array
  allocate(vlm(grid%ncoef_r))

  ! transform back
  call grid%SH_trans(u,vlm)
  print *, ulm(ilm),vlm(ilm)
  
  open(newunit = io,file='test_spherical_harmonics.out')
  write(io,*) grid%nth,grid%nph,0.0_dp
  do ith = 1,grid%nth
     th = grid%th(ith)
     do iph = 1,grid%nph
        ph = grid%ph(iph)
        write(io,*) ph,th,u(iph,ith)
     end do
  end do
  close(io)


  !=====================================!
  !     test complex transformations    !
  !=====================================!  
  
  ! allocate the spatial array
  allocate(w(grid%nph,grid%nth))

  ! allocate the coefficient array
  allocate(wlm(grid%ncoef_c))

  ! set the coefficients
  wlm = 0.0_dp
  ilm = grid%cindex(2,2)
  wlm(ilm) = 1.0_dp

  ! form the spatial field
  call grid%SH_itrans(wlm,w)

  ! allocate a new coefficient array
  allocate(xlm(grid%ncoef_c))

  ! transform back
  call grid%SH_trans(w,xlm)
  print *, wlm(ilm),xlm(ilm)
  
  open(newunit = io,file='test_spherical_harmonics_complex.out')
  write(io,*) grid%nth,grid%nph,0.0_dp
  do ith = 1,grid%nth
     th = grid%th(ith)
     do iph = 1,grid%nph
        ph = grid%ph(iph)
        write(io,*) ph,th,real(w(iph,ith))
     end do
  end do
  close(io)


  
  
end program test_spherical_harmonics




