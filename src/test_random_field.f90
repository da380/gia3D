program test_random_field

  use module_constants
  use module_util
  use module_spherical_harmonics
  use module_random_fields
  implicit none

  integer(i4b) :: lmax,l,ith,iph,io,inode,ispec,i,nx,n
  real(dp) :: lambda,s,th,ph,x1,x2,x,y,sigma,dx,x0
  real(dp), dimension(:), allocatable :: xx,ff
  real(dp), dimension(:,:), allocatable :: v
  type(gauss_legendre_grid) :: grid
  type(gaussain_random_scalar_field_sphere) :: u
  type(interp_1D_cubic) :: fun,cfun
  class(GRF_1D), allocatable :: rfun


  !---------------------------------------------!
  !     test random functions on an interval    !
  !---------------------------------------------!
  
  x1 = 0.0_dp
  x2 = 1.0_dp
  if(.not.found_command_argument('-s',s)) stop 's missing'
  if(.not.found_command_argument('-lambda',lambda)) stop 'lambda missing'
  if(.not.found_command_argument('-sigma',sigma)) stop 'sigma missing'


  ! set the mean function
  n  = 100
  dx = (x2-x1)/(n-1)
  allocate(xx(n),ff(n))
  do i = 1,n
     x = x1 + (i-1)*dx
     xx(i) = x
     ff(i) = x**2*sin(3*x)
  end do

!  rfun = GRF_1D_SEM(x1,x2,lambda,s,sigma)
  rfun = GRF_1D_Fourier(x1,x2,lambda,s,sigma)

  x0 = 0.5_dp
  call rfun%corr(x0,cfun)
  call rfun%realise(fun)
  
  nx = 20*(x2-x1)/lambda
  dx = (x2-x1)/(nx-1)
  x0 = 0.3_dp
  open(newunit = io,file='random.out')
  do i = 1,nx
     x = x1 + (i-1)*dx
     write(io,*) x,fun%f(x),cfun%f(x)
  end do
  close(io)
  


  
 
  
  
  !----------------------------------------!
  !   test random functions on a sphere    !
  !----------------------------------------!

  ! get arguments
!  if(.not.found_command_argument('-lmax',lmax)) stop 'lmax missing'
!  if(.not.found_command_argument('-s',s)) stop 's missing'
!  if(.not.found_command_argument('-lambda',lambda)) stop 'lambda missing'
  
  ! build the grid
!  call grid%build(lmax)  

  ! allocate the random field
!  call u%build(grid,lambda,s)

  ! allocate coefficient array and spatial array
!  allocate(v(grid%nph,grid%nth))
  
  ! form realisation of the random field
!  call u%realise()

!  call grid%SH_itrans(u%ulm,v)
  
  ! write out the field
!  open(newunit = io,file='random.out')
!  write(io,*) grid%nth,grid%nph,0.0_dp
!  do ith = 1,grid%nth
!     th = grid%th(ith)
!     do iph = 1,grid%nph
!        ph = grid%ph(iph)
!        write(io,*) ph,th,v(iph,ith)
!     end do
!  end do
!  close(io)  

  

  
  
  
end program test_random_field
