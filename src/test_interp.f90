program test_interp

  use module_constants
  use module_interp
  implicit none


  integer(i4b) :: nx,ny,ix,iy,io,jo
  real(dp) :: x,x1,x2,dx,y,y1,y2,dy,f,start,finish,x0 = 0.2_dp
  real(dp), dimension(:), allocatable :: xx,yy,gg
  real(dp), dimension(:,:), allocatable :: ff

  type(interp_1D_cubic) :: k
  type(interp_2D) :: g
  type(interp_2D_bicubic_spline) :: h

  
  ! make the function to be interpolated

  nx = 100
  x1 = 0.0_dp
  x2 = 1.0_dp
  dx = (x2-x1)/(nx-1)
  
  ny = 100
  y1 = 0.0_dp
  y2 = 1.0_dp
  dy = (y2-y1)/(ny-1)

  allocate(xx(nx),yy(ny),ff(nx,ny),gg(ny))

  do iy = 1,ny
     do ix = 1,nx

        x = x1 + (ix-1)*dx
        y = y1 + (iy-1)*dy

        xx(ix)    = x
        yy(iy)    = y
        ff(ix,iy) = func(x,y)
        
     end do
     gg(iy) = func(x0,y)
  end do

  call k%parms(lnat=.false.,rnat = .false.)
  call g%set(xx,yy,ff)
  call h%set(xx,yy,ff)
  call k%set(yy,gg)
  
  ! set the values for plotting
  nx = nx*3
  dx = (x2-x1)/(nx-1)  
  ny = ny*6
  dy = (y2-y1)/(ny-1)

  open(newunit = io,file='interp.out')
  open(newunit = jo,file='interp_1D.out')
  write(io,*) ny,nx,0
  do iy = 1,ny
     do ix = 1,nx

        x = x1 + (ix-1)*dx
        y = y1 + (iy-1)*dy
        
        write(io,*) x,y,func(x,y)-g%f(x,y)
        
     end do
     write(jo,*) y,func(x0,y)- k%f(y)
  end do
  close(io)
  close(jo)


contains

  real(dp) function func(x,y) result(f)
    real(dp), intent(in) :: x,y
    f = 2.0_dp*x + 5.0_dp*y + x**2 + y**2
    return
  end function func
  
end program test_interp
