program test_interp

  use module_constants
  use module_interp
  implicit none

  logical :: store
  integer(i4b) :: nx,ny,ix,iy,io,jx,jy
  real(dp) :: x,x1,x2,dx,y,y1,y2,dy,f,start,finish
  real(dp), dimension(:), allocatable :: xx,yy
  real(dp), dimension(:,:), allocatable :: ff

  ! make the function to be interpolated

  nx = 20
  x1 = 0.0_dp
  x2 = 1.0_dp
  dx = (x2-x1)/(nx-1)
  
  ny = 40
  y1 = 0.0_dp
  y2 = 2.0_dp
  dy = (y2-y1)/(ny-1)

  allocate(xx(nx),yy(ny),ff(nx,ny))


  do ix = 1,nx
     do iy = 1,ny

        x = x1 + (ix-1)*dx
        y = y1 + (iy-1)*dy
        f = x - y

        xx(ix)    = x
        yy(iy)    = y
        ff(ix,iy) = f
        
     end do
  end do



  ! set the values for plotting
  nx = nx*100
  dx = (x2-x1)/(nx-1)  
  ny = ny*100
  dy = (y2-y1)/(ny-1)

  jx = 0
  jy = 0
  store = .true.

  call cpu_time(start)
  do ix = 1,nx
     do iy = 1,ny

        x = x1 + (ix-1)*dx
        y = y1 + (iy-1)*dy

        if(store) then
           f = bilinear_interp(xx,yy,ff,x,y,jx,jy)
        else
           f = bilinear_interp(xx,yy,ff,x,y)
        end if
        
     end do
  end do
  call cpu_time(finish)
  print *, finish-start
  
  
end program test_interp
