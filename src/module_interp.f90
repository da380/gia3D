module module_interp

  use module_constants
  use module_error
  implicit none

  type :: interp_1D
     private
     logical :: built = .false.
     logical :: ascnd
     integer(i4b) :: n
     integer(i4b) :: isave = 0
     real(dp), dimension(:), allocatable :: xx
     real(dp), dimension(:), allocatable :: yy
     real(dp) :: xsave,dxsave
     real(dp) :: yp
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
     real(dp), dimension(:), allocatable :: yy2
     real(dp) :: ypp
   contains
     procedure :: set => set_interp_1D_cubic
     procedure :: delete => delete_interp_1D_cubic
     procedure :: f => value_interp_1D_cubic
     procedure :: fp => derivative_interp_1D_cubic
     procedure :: fpp => derivative2_interp_1D_cubic
  end type interp_1D_cubic
 private :: set_interp_1D_cubic,value_interp_1D_cubic, &
             delete_interp_1D_cubic
  
contains


  !======================================================================!
  !             general procedures for interpolation in 1D               !
  !======================================================================!
  
  
  subroutine set_interp_1D(self,xx,yy)
    implicit none
    class(interp_1D) :: self
    real(dp), dimension(:), intent(in) :: xx,yy
    
    integer(i4b) :: n,j,i
    
    ! check the size of the arrays
    n = size(xx)
    call error(n /= size(yy),'set_linear_interp_1D', &
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
    if(j == n-1) then
       self%ascnd = .true.
    else if(j == -n+1) then
       self%ascnd = .false.
    else
       call error(n /= size(yy),'set_linear_interp_1D', &
                                'array not consecutive')
    end if

    
    ! store the dimensions
    self%n = n

    ! store the data
    self%xx = xx
    self%yy = yy

    ! set dxsave
    self%dxsave = abs(xx(n)-xx(1))/5.0_dp

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
    if(self%isave /= 0 .and. abs(x-self%xsave) > self%dxsave) then
       self%isave = 0
    end if

    ! find the point
    find_interp_1D = hunt_list(self%xx,x,self%isave)

    ! save the previous values
    self%isave = find_interp_1D
    self%xsave = x
    
    return
  end function find_interp_1D

  function value_interp_1D(self,x) result(y)
    implicit none
    class(interp_1D) :: self
    real(dp), intent(in) :: x

    integer(i4b) :: i
    real(dp) :: x1,x2,y1,y2,y
    
    i  = self%find(x)
    x1 = self%xx(i)
    x2 = self%xx(i+1)

    y1 = self%yy(i)
    y2 = self%yy(i+1)

    self%yp = (y2-y1)/(x2-x1)
    y = y1 + self%yp*(x-x1)
    
    return
  end function value_interp_1D


  function derivative_interp_1D(self,x) result(yp)
    implicit none
    class(interp_1D) :: self
    real(dp), intent(in) :: x

    integer(i4b) :: i
    real(dp) :: y,yp
    
    if(x /= self%xsave) y = self%f(x)
    yp = self%yp
    
    return
  end function derivative_interp_1D


  subroutine delete_interp_1D(self)
    implicit none
    class(interp_1D) :: self
    deallocate(self%xx,self%yy)
    self%built = .false.    
    return
  end subroutine delete_interp_1D
  
    
  !==================================================================!
  !              routines for cubic spline interpolation             !
  !==================================================================!


  subroutine set_interp_1D_cubic(self,xx,yy)

    ! sets up a cubic spline from the data using natural
    ! boundary conditions.
    
    use module_constants
    implicit none
    class(interp_1D_cubic) :: self
    real(dp), dimension(:), intent(in) :: xx,yy

    integer(i4b) :: i,n,info
    real(dp), dimension(:), allocatable :: dl,d,du
    real(dp), dimension(:,:), allocatable :: b
     
    ! set up the base type
    call self%interp_1D%set(xx,yy)
    
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
       b(i,1) =  (yy(i+1)-yy(i))/(xx(i+1)-xx(i)) &
                -(yy(i)-yy(i-1))/(xx(i)-xx(i-1))
    end do

    ! fix the boundary conditions
    b(1,1) = 0.0_dp
    b(n,1) = 0.0_dp
    d(1) = 1.0_dp
    d(n) = 1.0_dp
    
    ! solve the tridiagonal system
    call dgtsv (n, 1, dl, d, du, b, n, info)

    ! store the computed second derivatives
    self%yy2 = b(:,1)

    ! deallocate temporary arrays
    deallocate(dl,d,du)

    
    return
  end subroutine set_interp_1D_cubic
    

    
  function value_interp_1D_cubic(self,x) result(y)

    ! evaluates the cubic spline at a given location.
    ! Also computes and stores first and second derivatives, 
    ! with the results stored as type data  
    implicit none
    class(interp_1D_cubic) :: self
    real(dp), intent(in) :: x

    integer(i4b) :: i1,i2
    real(dp) :: a,b,h,x1,x2,y1,y2,y21,y22,y

    ! get the indices and points
    i1 = self%find(x)
    i2 = i1+1

    ! compute some parameters
    x1  = self%xx(i1)
    x2  = self%xx(i2)
    y1  = self%yy(i1)
    y2  = self%yy(i2)
    y21 = self%yy2(i1)
    y22 = self%yy2(i2)
    h   = x2-x1
    a   = (x2-x)/h
    b   = (x-x1)/h

    ! compute the function  
    y = a*y1 + b*y2 + ((a*a*a-a)*y21 +  (b*b*b-b)*y22)*(h*h)/6.0_dp

    ! compute the first and second derivatives and store values   
    self%yp  = (y2-y1)/h  - (3.0_dp*a*a-1)*h*y21/6.0_dp  &
                          + (3.0_dp*b*b-1)*h*y22/6.0_dp
    self%ypp = a*y21 + b*y22
     
    
    return
  endfunction value_interp_1D_cubic


  function derivative_interp_1D_cubic(self,x) result(yp)

    ! evaluates the cubic spline at a given location.
    ! Also computes and stores first and second derivatives, 
    ! with the results stored as type data  
    implicit none
    class(interp_1D_cubic) :: self
    real(dp), intent(in) :: x
    real(dp) :: y,yp

    if(x /= self%xsave) y = self%f(x)
    yp = self%yp

    return
  end function derivative_interp_1D_cubic


  function derivative2_interp_1D_cubic(self,x) result(ypp)

    ! evaluates the cubic spline at a given location.
    ! Also computes and stores first and second derivatives, 
    ! with the results stored as type data  
    implicit none
    class(interp_1D_cubic) :: self
    real(dp), intent(in) :: x
    real(dp) :: y,ypp

    if(x /= self%xsave) y = self%f(x)
    ypp = self%ypp

    return
  end function derivative2_interp_1D_cubic
  
  
  subroutine delete_interp_1D_cubic(self)
    implicit none
    class(interp_1D_cubic) :: self
    deallocate(self%xx,self%yy,self%yy2)
    self%built = .false.    
    return
  end subroutine delete_interp_1D_cubic

  
end module module_interp



