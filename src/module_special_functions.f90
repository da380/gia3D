module module_special_functions
  
  use module_constants
  use module_error
  implicit none
  
  type, abstract :: orthogonal_polynomial
     !
     ! defines an orthogonal polynomial by its three
     ! term recursion relation.
     !
     ! This relation has the form used in Goulb and Welsch:
     !
     ! p_{j} = (a_{j}x + b_{j})p_{j-1} - c_{j}p_{j-1}
     !
     ! with p_{-1} = 0, p_{0} = 1
     !
     private
     real(dp) :: x1
     real(dp) :: x2
     real(dp) :: mu
   contains
     procedure :: xl => xl_orthogonal_polynomial
     procedure :: xr => xr_orthogonal_polynomial
     procedure :: mu0 => mu0_orthogonal_polynomial
     procedure(op_f_sr_r), deferred :: weight
     procedure(op_f_si_r), deferred :: a
     procedure(op_f_si_r), deferred :: b
     procedure(op_f_si_r), deferred :: c
     procedure :: eval => evaluate_orthogonal_polynomial
  end type orthogonal_polynomial

  type, extends(orthogonal_polynomial) :: legendre_polynomial
   contains     
     procedure :: weight => legendre_weight
     procedure :: a      => legendre_a
     procedure :: b      => legendre_b
     procedure :: c      => legendre_c
  end type legendre_polynomial
  
  
  abstract interface
     function op_f_si_r(p,i) result(f)
       use module_constants
       import :: orthogonal_polynomial
       implicit none
       class(orthogonal_polynomial), intent(in) :: p
       integer(i4b), intent(in) :: i
       real(dp) :: f
     end function op_f_si_r
     function op_f_sr_r(p,x) result(f)
       use module_constants
       import :: orthogonal_polynomial
       implicit none
       class(orthogonal_polynomial), intent(in) :: p
       real(dp), intent(in) :: x
       real(dp) :: f
     end function op_f_sr_r     
  end interface

contains

  !=========================================================================!
  !                     orthogonal polynomial procedures                    !
  !=========================================================================!  


  function xl_orthogonal_polynomial(poly) result(xl)
    implicit none
    class(orthogonal_polynomial), intent(in) :: poly
    real(dp) :: xl
    xl = poly%x1
    return
  end function xl_orthogonal_polynomial

  function xr_orthogonal_polynomial(poly) result(xr)
    implicit none
    class(orthogonal_polynomial), intent(in) :: poly
    real(dp) :: xr
    xr = poly%x2
    return
  end function xr_orthogonal_polynomial

  function mu0_orthogonal_polynomial(poly) result(mu0)
    implicit none
    class(orthogonal_polynomial), intent(in) :: poly
    real(dp) :: mu0
    mu0 = poly%mu
    return
  end function mu0_orthogonal_polynomial

  function evaluate_orthogonal_polynomial(poly,n,x) result(p)    
    implicit none
    class(orthogonal_polynomial), intent(in) :: poly
    integer(i4b), intent(in) :: n
    real(dp), intent(in) :: x
    real(dp) :: p

    integer(i4b) :: i
    real(dp) :: pm1,pp1

    pm1 = 0.0_dp
    p   = 1.0_dp
    do i = 1,n
       pp1 = (poly%a(i)*x+poly%b(i))*p - poly%c(i)*pm1
       pm1 = p
       p = pp1       
    end do
    return
  end function evaluate_orthogonal_polynomial


  !=========================================================================!
  !                        legendre polynomial procedures                   !
  !=========================================================================!  

  function legendre() result(p)
    implicit none
    type(legendre_polynomial) :: p
    p%x1 = -1.0_dp
    p%x2 =  1.0_dp
    p%mu =  2.0_dp
    return
  end function legendre
  
  
  function legendre_weight(p,x) result(w)
    implicit none
    class(legendre_polynomial), intent(in) :: p
    real(dp), intent(in) :: x
    real(dp) :: w
    w = 1.0_dp
    return
  end function legendre_weight

  function legendre_a(p,i) result(a)
    implicit none
    class(legendre_polynomial), intent(in) :: p
    integer(i4b), intent(in) :: i
    real(dp) :: a
    if(i < 1) then
       a = 0
    else       
       a = (2*i-1)
       a = a/i
    end if
    return
  end function legendre_a

  function legendre_b(p,i) result(b)
    implicit none
    class(legendre_polynomial), intent(in) :: p
    integer(i4b), intent(in) :: i
    real(dp) :: b
    b = 0.0_dp
    return
  end function legendre_b

  function legendre_c(p,i) result(c)
    implicit none
    class(legendre_polynomial), intent(in) :: p
    integer(i4b), intent(in) :: i
    real(dp) :: c
    if(i < 1) then
       c = 0
    else       
       c = (i-1)
       c = c/i
    end if
    return
  end function legendre_c
  


  !=========================================================================!
  !                         jacobi polynomial procedures                    !
  !=========================================================================!
  


  
end module module_special_functions
