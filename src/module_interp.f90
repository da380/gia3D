module module_interp

  use module_constants
  use module_error
  implicit none

  type :: interp_1D
     private
     logical :: built = .false.
     integer(i4b) :: n
     integer(i4b) :: i_save = 0
     real(dp), dimension(:), allocatable :: xx
     real(dp), dimension(:), allocatable :: ff
     real(dp) :: x_save,dx_save
     real(dp) :: fp_save
   contains
     procedure :: set    => set_interp_1D
     procedure :: delete => delete_interp_1D
     procedure :: find   => find_interp_1D
     procedure :: f    => value_interp_1D
     procedure :: fp => derivative_interp_1D
  end type interp_1D
  private :: set_interp_1D,find_interp_1D, &
             value_interp_1D,delete_interp_1D
  
  
  type, extends(interp_1D) :: interp_1D_cubic
     private
     real(dp), dimension(:), allocatable :: ff2
     real(dp) :: fpp_save
   contains
     procedure :: set => set_interp_1D_cubic
     procedure :: delete => delete_interp_1D_cubic
     procedure :: f => value_interp_1D_cubic
     procedure :: fp => derivative_interp_1D_cubic
     procedure :: fpp => derivative2_interp_1D_cubic
  end type interp_1D_cubic
 private :: set_interp_1D_cubic,value_interp_1D_cubic, &
             delete_interp_1D_cubic

 

 type :: interp_2D
    private
    logical :: built = .false.
    integer(i4b) :: nx
    integer(i4b) :: ny
    real(dp), dimension(:), allocatable :: xx
    real(dp), dimension(:), allocatable :: yy
    real(dp), dimension(:,:), allocatable :: ff
    integer(i4b) :: ix_save = 0
    integer(i4b) :: iy_save = 0
    real(dp) :: x_save,dx_save
    real(dp) :: y_save,dy_save
  contains
    procedure :: set => set_interp_2D
    procedure :: delete => delete_interp_2D
!    procedure :: find   => find_interp_2D
!    procedure :: f    => value_interp_2D    
 end type interp_2D

 
contains


  !======================================================================!
  !             general procedures for interpolation in 1D               !
  !======================================================================!
  
  
  subroutine set_interp_1D(self,xx,ff)
    implicit none
    class(interp_1D) :: self
    real(dp), dimension(:), intent(in) :: xx,ff
    
    integer(i4b) :: n,j,i
    
    ! check the size of the arrays
    n = size(xx)
    call error(n /= size(ff),'set_linear_interp_1D', &
                             'dimensions of arrays do not match')

    
    ! check the arrays are consecutive
    j = 0
    do i = 1,n-1
       if(xx(i+1) >  xx(i)) then
          j = j+1
       else if(xx(i+1) > xx(i)) then
          j = j-1
       end if
    end do
    call error(j /= n-1 .and. j /= -n+1,'set_linear_interp_1D','array not consecutive')

    
    ! store the dimensions
    self%n = n

    ! store the data
    self%xx = xx
    self%ff = ff

    ! set dx_save
    self%dx_save = abs(xx(n)-xx(1))/5.0_dp

    ! build done
    self%built = .true.
    
    return
  end subroutine set_interp_1D

  integer(i4b) function find_interp_1D(self,x)
    use module_util, only: bisect_list,hunt_list
    implicit none
    class(interp_1D) :: self
    real(dp), intent(in) :: x

    integer(i4b) :: il,iu,im,n

    ! is it worth hunting? 
    if(self%i_save /= 0 .and. abs(x-self%x_save) > self%dx_save) then
       self%i_save = 0
    end if

    ! find the point
    find_interp_1D = hunt_list(self%xx,x,self%i_save)

    ! _save the previous values
    self%i_save = find_interp_1D
    self%x_save = x
    
    return
  end function find_interp_1D

  function value_interp_1D(self,x) result(f)
    implicit none
    class(interp_1D) :: self
    real(dp), intent(in) :: x

    integer(i4b) :: i
    real(dp) :: x1,x2,f1,f2,f
    
    i  = self%find(x)
    x1 = self%xx(i)
    x2 = self%xx(i+1)

    f1 = self%ff(i)
    f2 = self%ff(i+1)

    self%fp_save = (f2-f1)/(x2-x1)
    f = f1 + self%fp_save*(x-x1)
    
    return
  end function value_interp_1D


  function derivative_interp_1D(self,x) result(fp)
    implicit none
    class(interp_1D) :: self
    real(dp), intent(in) :: x

    integer(i4b) :: i
    real(dp) :: f,fp
    
    if(x /= self%x_save) f = self%f(x)
    fp = self%fp_save
    
    return
  end function derivative_interp_1D


  subroutine delete_interp_1D(self)
    implicit none
    class(interp_1D) :: self
    deallocate(self%xx,self%ff)
    self%built = .false.    
    return
  end subroutine delete_interp_1D
  
    
  !==================================================================!
  !              routines for cubic spline interpolation             !
  !==================================================================!


  subroutine set_interp_1D_cubic(self,xx,ff)

    ! sets up a cubic spline from the data using natural
    ! boundary conditions.
    
    use module_constants
    implicit none
    class(interp_1D_cubic) :: self
    real(dp), dimension(:), intent(in) :: xx,ff

    integer(i4b) :: i,n,info
    real(dp), dimension(:), allocatable :: dl,d,du
    real(dp), dimension(:,:), allocatable :: b
     
    ! set up the base type
    call self%interp_1D%set(xx,ff)
    
    ! dimension of the arrays
    n = self%n

    ! build the matrix
    allocate(dl(n-1),d(n),du(n-1),b(n,1))

    ! upper diagonal
    du(1) = 0.0_dp
    do i = 2,n-1
       du(i) = (xx(i+1)-xx(i))/6.0_dp
    end do

    ! lower diagonl
    do i = 2,n-1
       dl(i-1) = (xx(i)-xx(i-1))/6.0_dp       
    end do
    dl(n-1) = 0.0_dp

    ! diagonal and right hand side
    do i = 2,n-1
       d(i)   =  (xx(i+1)-xx(i-1))/3.0_dp
       b(i,1) =  (ff(i+1)-ff(i))/(xx(i+1)-xx(i)) &
                -(ff(i)-ff(i-1))/(xx(i)-xx(i-1))
    end do

    ! fix the boundary conditions
    b(1,1) = 0.0_dp
    b(n,1) = 0.0_dp
    d(1) = 1.0_dp
    d(n) = 1.0_dp
    
    ! solve the tridiagonal system
    call dgtsv (n, 1, dl, d, du, b, n, info)

    ! store the computed second derivatives
    self%ff2 = b(:,1)

    ! deallocate temporary arrays
    deallocate(dl,d,du)

    
    return
  end subroutine set_interp_1D_cubic
    

    
  function value_interp_1D_cubic(self,x) result(f)

    ! evaluates the cubic spline at a given location.
    ! Also computes and stores first and second derivatives, 
    ! with the results stored as type data  
    implicit none
    class(interp_1D_cubic) :: self
    real(dp), intent(in) :: x

    integer(i4b) :: i1,i2
    real(dp) :: a,b,h,x1,x2,f1,f2,f21,f22,f

    
    ! get the indices and points
    i1 = self%find(x)
    i2 = i1+1

    ! compute some parameters
    x1  = self%xx(i1)
    x2  = self%xx(i2)
    f1  = self%ff(i1)
    f2  = self%ff(i2)
    f21 = self%ff2(i1)
    f22 = self%ff2(i2)
    h   = x2-x1
    a   = (x2-x)/h
    b   = (x-x1)/h

    ! compute the function  
    f = a*f1 + b*f2 + ((a*a*a-a)*f21 +  (b*b*b-b)*f22)*(h*h)/6.0_dp

    ! compute the first and second derivatives and store values   
    self%fp_save  = (f2-f1)/h  - (3.0_dp*a*a-1)*h*f21/6.0_dp  &
                          + (3.0_dp*b*b-1)*h*f22/6.0_dp
    self%fpp_save = a*f21 + b*f22
     
    
    return
  endfunction value_interp_1D_cubic


  function derivative_interp_1D_cubic(self,x) result(fp)

    ! evaluates the cubic spline at a given location.
    ! Also computes and stores first and second derivatives, 
    ! with the results stored as type data  
    implicit none
    class(interp_1D_cubic) :: self
    real(dp), intent(in) :: x
    real(dp) :: f,fp

    if(x /= self%x_save) f = self%f(x)
    fp = self%fp_save

    return
  end function derivative_interp_1D_cubic


  function derivative2_interp_1D_cubic(self,x) result(fpp)

    ! evaluates the cubic spline at a given location.
    ! Also computes and stores first and second derivatives, 
    ! with the results stored as type data  
    implicit none
    class(interp_1D_cubic) :: self
    real(dp), intent(in) :: x
    real(dp) :: f,fpp

    if(x /= self%x_save) f = self%f(x)
    fpp = self%fpp_save

    return
  end function derivative2_interp_1D_cubic
  
  
  subroutine delete_interp_1D_cubic(self)
    implicit none
    class(interp_1D_cubic) :: self
    deallocate(self%xx,self%ff,self%ff2)
    self%built = .false.    
    return
  end subroutine delete_interp_1D_cubic



  !==================================================================!
  !               routines for 2D bilinear interpolation             !
  !==================================================================!

  subroutine set_interp_2D(self,xx,yy,ff)
    implicit none
    class(interp_2D) :: self
    real(dp), dimension(:), intent(in) :: xx,yy
    real(dp), dimension(:,:), intent(in) :: ff
    
    integer(i4b) :: nx,ny,i,j
    
    ! check the size of the arrays
    nx = size(xx)
    ny = size(yy)
    call error(nx /= size(ff,1) .or. ny /= size(ff,2), &
               'set_linear_interp_1D','dimensions of arrays do not match')

    
    ! check the arrays are consecutive
    j = 0
    do i = 1,nX-1
       if(xx(i+1) >  xx(i)) then
          j = j+1
       else if(xx(i+1) > xx(i)) then
          j = j-1
       end if
    end do
    call error(j /= nx-1 .and. j /= -nx+1,'set_linear_interp_1D','x-array not consecutive')

    
    ! store the dimensions
    self%nx = nx
    self%ny = ny

    ! store the data
    self%xx = xx
    self%yy = yy
    self%ff = ff

    ! set dx_save
    self%dx_save = abs(xx(nx)-xx(1))/5.0_dp
    self%dy_save = abs(yy(ny)-yy(1))/5.0_dp

    ! build done
    self%built = .true.
    
    return
  end subroutine set_interp_2D


  subroutine delete_interp_2D(self)
    implicit none
    class(interp_2D) :: self
    deallocate(self%xx,self%yy,self%ff)
    self%built = .false.    
    return
  end subroutine delete_interp_2D
  
end module module_interp



