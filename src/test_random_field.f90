program test_random_field

  use module_constants
  use module_util
  use module_spherical_harmonics
  use module_random_fields
  implicit none

  integer(i4b) :: lmax,l,ith,iph,io,inode,ispec,i
  real(dp) :: lambda,s,th,ph,x1,x2,x,y,sigma
  real(dp), dimension(:,:), allocatable :: v
  type(gauss_legendre_grid) :: grid
  type(gaussain_random_scalar_field_sphere) :: u
  type(gaussian_random_scalar_field_interval) :: w


  !---------------------------------------------!
  !     test random functions on an interval    !
  !---------------------------------------------!
  
  x1 = 0.0_dp
  x2 = 1.0_dp
  if(.not.found_command_argument('-s',s)) stop 's missing'
  if(.not.found_command_argument('-lambda',lambda)) stop 'lambda missing'
  if(.not.found_command_argument('-sigma',sigma)) stop 'sigma missing'
  
  call w%build(x1,x2,lambda,s,sigma)
!  call w%set_mean(f)  
  call w%realise()
  
  
  open(newunit = io,file='random.out')
  do ispec = w%ispec1,w%ispec2
     do inode = 1,w%ngll

        i = w%ibool(inode,ispec)
        write(io,*) w%x(inode,ispec),w%u(i)
        
     end do
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

  

contains


  real(dp) function f(x) result(y)
    real(dp), intent(in) :: x
    y = 4.0_dp*x*(1-x)
    return
  end function f
  
  
  
end program test_random_field
