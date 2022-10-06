program test_random_field

  use module_constants
  use module_util
  use module_spherical_harmonics
  use module_random_fields
  implicit none

  integer(i4b) :: lmax,l,ith,iph,io,inode,ispec
  real(dp) :: lambda,s,th,ph,x1,x2
  real(dp), dimension(:,:), allocatable :: v
  type(gauss_legendre_grid) :: grid
  type(gaussain_random_scalar_field_sphere) :: u
  type(gaussian_random_scalar_field_interval) :: w


  !---------------------------------------------!
  !     test random functions on an interval    !
  !---------------------------------------------!
  
  x1 = 0.0_dp
  x2 = 1.0_dp
  lambda = 0.1_dp
  s = 2.0_dp
 
  call w%set(x1,x2,lambda,s)

  do ispec = 1,w%nspec

     do inode = 1,w%ngll

!        print *, w%x(inode,ispec)
        
     end do
     
  end do


  
  !----------------------------------------!
  !   test random functions on a sphere    !
  !----------------------------------------!

  ! get arguments
  if(.not.found_command_argument('-lmax',lmax)) stop 'lmax missing'
  if(.not.found_command_argument('-s',s)) stop 's missing'
  if(.not.found_command_argument('-lambda',lambda)) stop 'lambda missing'
  
  ! build the grid
  call grid%build(lmax)  

  ! allocate the random field
  call u%set(grid,lambda,s,real=.true.)

  ! allocate coefficient array and spatial array
  allocate(v(grid%nph,grid%nth))
  
  ! form realisation of the random field
  call u%realise()

  call grid%SH_itrans(u%ulm,v)
  
  ! write out the field
  open(newunit = io,file='random.out')
  write(io,*) grid%nth,grid%nph,0.0_dp
  do ith = 1,grid%nth
     th = grid%th(ith)
     do iph = 1,grid%nph
        ph = grid%ph(iph)
        write(io,*) ph,th,v(iph,ith)
     end do
  end do
  close(io)  

  

  
  
end program test_random_field
