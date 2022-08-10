module module_quadrature

  
  use module_constants
  use module_error
  
  type, abstract :: quadrature
     private
     integer(i4b) :: n
     real(dp) :: a,b
     real(dp), dimension(:), allocatable :: x
     real(dp), dimension(:), allocatable :: w
   contains
     procedure :: order => order_quadrature     
     procedure :: points => points_quadrature
     procedure :: weights => weights_quadrature
     procedure :: trans => trans_quadrature
     procedure :: intr => int_quadrature_r
     procedure :: intc => int_quadrature_c
     generic :: int => intr,intc
  end type quadrature
  private :: int_quadrature_r,int_quadrature_c, &
             points_quadrature, weights_quadrature, &
             trans_quadrature,order_quadrature

  
  type, extends(quadrature) :: trapezoid_quadrature
   contains     
     procedure :: set => set_trapezoid_quadrature
  end type trapezoid_quadrature
  private :: set_trapezoid_quadrature

  
  type, extends(quadrature) :: gaussian_quadrature
   contains     
     procedure :: set => set_gaussian_quadrature
  end type gaussian_quadrature
  private :: set_gaussian_quadrature

  
  abstract interface
     
     function f_r_r(x) result(f)
       use module_constants
       implicit none
       real(dp), intent(in) :: x
       real(dp) :: f
     end function f_r_r
     
     function f_r_c(x) result(f)
       use module_constants
       implicit none
       real(dp), intent(in) :: x
       complex(dpc) :: f
     end function f_r_c

  end interface
  
contains


  !=========================================================================!
  !                      general quadrature procedures                      !
  !=========================================================================!
  
  function order_quadrature(self) result(n)
    implicit none
    class(quadrature), intent(in) :: self
    integer(i4b) :: n
    n = self%n
    return
  end function order_quadrature

  function points_quadrature(self) result(x)
    implicit none
    class(quadrature), intent(in) :: self
    real(dp), dimension(:), allocatable :: x
    x = self%x
    return
  end function points_quadrature


  function weights_quadrature(self) result(w)
    implicit none
    class(quadrature), intent(in) :: self
    real(dp), dimension(:), allocatable :: w
    w = self%w
    return
  end function weights_quadrature
  
  function int_quadrature_r(self,f) result(int)
    implicit none
    class(quadrature), intent(inout) :: self
    procedure(f_r_r) :: f
    real(dp) :: int
    integer(i4b) :: i
    int = 0.0_dp
    do i = 1,self%n
       int = int + f(self%x(i))*self%w(i)
    end do    
    return
  end function int_quadrature_r

  function int_quadrature_c(self,f) result(int)
    implicit none
    class(quadrature), intent(inout) :: self
    procedure(f_r_c) :: f
    complex(dpc) :: int
    integer(i4b) :: i
    int = 0.0_dp
    do i = 1,self%n
       int = int + f(self%x(i))*self%w(i)
    end do    
    return
  end function int_quadrature_c
  
  subroutine trans_quadrature(self,c,d)
    implicit none
    class(quadrature), intent(inout) :: self
    real(dp), intent(in) :: c,d
    real(dp) :: jac,a,b
    a = self%a
    b = self%b
    jac = (d-c)/(b-a)
    self%x = (b-self%x)*c + (self%x-a)*d
    self%x = self%x/(b-a)
    self%w = self%w*jac
    self%a = c
    self%b = d        
    return
  end subroutine trans_quadrature


  !=========================================================================!
  !                      specific quadrature procedures                     !
  !=========================================================================!
  

  subroutine set_trapezoid_quadrature(self,n,a,b)
    implicit none
    class(trapezoid_quadrature), intent(inout) :: self
    integer(i4b), intent(in) :: n
    real(dp), intent(in) :: a,b

    integer(i4b) :: i
    real(dp) :: dx
    real(dp), dimension(n) :: x,w

    dx = (b-a)/(n-1)
    do i = 1,n
       x(i) = a+(i-1)*dx
    end do

    w = dx
    w(1) = 0.5_dp*dx
    w(n) = 0.5_dp*dx

    self%n = n
    self%a = a
    self%b = b
    self%x = x
    self%w = w
    
    return
  end subroutine set_trapezoid_quadrature
  
  
  subroutine set_gaussian_quadrature(self,n,poly)
    implicit none
    class(gaussian_quadrature), intent(inout) :: self
    integer(i4b), intent(in) :: n
    character(len=*), intent(in) :: poly

    integer(i4b) :: il,iu,abstol,m,ldz,lwork,liwork,info,i
    integer(i4b), dimension(10*n) :: iwork
    integer(i4b), dimension(2*n) :: isuppz
    real(dp) :: vl,vu,mu,a,b
    real(dp), dimension(n) :: d,e,w
    real(dp), dimension(18*n) :: work
    real(dp), dimension(n,n) :: z

    self%n = n    
    if(poly == 'legendre') then
       self%a = -1.0_dp
       self%b =  1.0_dp
       mu = 2.0_dp
    else 
       call error(.true.,' set_gaussian_quadrature','polynomial not defined')       
    end if
    
    ! build the matrix    
    do i = 1,n
       if(poly == 'legendre') then          
          d(i) = 0.0_dp
          e(i) = sqrt(i*i/(4.0_dp*i*i-1.0_dp))
       end if       
    end do

    ! solve the eigenvalue problem
    lwork  = 18*n
    liwork = 10*n
    call dstegr('V','A',n,d,e,vl,vu,il,iu,abstol, &
                 m,w,z,n,isuppz,work,lwork,iwork,liwork,info)		

    ! store the points and weights
    self%x = w
    self%w = mu*z(1,:)**2    

    return
  end subroutine set_gaussian_quadrature


  
end module module_quadrature
