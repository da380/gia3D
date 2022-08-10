module module_special_functions
  
  use module_constants
  use module_error
  implicit none
  
  type :: orthogonal_polynomial
     !
     ! defines an orthogonal polynomial by its three
     ! term recursion relation.
     !
     ! This relation has the form used in Goulb and Welsch:
     !
     ! p_{j} = (a_{j} + b_{j})p_{j-1} - c_{j}p_{j-1}
     !
     ! with p_{-1} = 0, p_{0} = 1
     !
     real(dp) :: x1
     real(dp) :: x2
     real(dp) :: mu
     procedure(f_sr_r), pointer :: weight
     procedure(f_si_r), pointer :: a
     procedure(f_si_r), pointer :: b
     procedure(f_si_r), pointer :: c
   contains
     procedure :: eval => evaluate_orthogonal_polynomial
  end type orthogonal_polynomial

  
  abstract interface
     function f_si_r(p,i) result(f)
       use module_constants
       import :: orthogonal_polynomial
       implicit none
       class(orthogonal_polynomial), intent(in) :: p
       integer(i4b), intent(in) :: i
       real(dp) :: f
     end function f_si_r
     function f_sr_r(p,x) result(f)
       use module_constants
       import :: orthogonal_polynomial
       implicit none
       class(orthogonal_polynomial), intent(in) :: p
       real(dp), intent(in) :: x
       real(dp) :: f
     end function f_sr_r     
  end interface

contains

  !=========================================================================!
  !                     orthogonal polynomial procedures                    !
  !=========================================================================!  


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
       pp1 = (poly%a(i)+poly%b(i))*x*p - poly%c(i)*pm1
       pm1 = p
       p = pp1       
    end do
    return
  end function evaluate_orthogonal_polynomial


  !=========================================================================!
  !                      legendre polynomial procedures                     !
  !=========================================================================!  
  
  subroutine set_legendre_polynomial(p)
    implicit none
    type(orthogonal_polynomial), intent(inout) :: p
    p%x1 = -1.0_dp
    p%x2 =  1.0_dp
    p%mu =  2.0_dp
    p%weight => legendre_weight
    p%a      => legendre_a
    p%b      => legendre_b
    p%c      => legendre_c
    return
  end subroutine set_legendre_polynomial
  
  function legendre_weight(p,x) result(w)
    implicit none
    class(orthogonal_polynomial), intent(in) :: p
    real(dp), intent(in) :: x
    real(dp) :: w
    w = 1.0_dp
    return
  end function legendre_weight

  function legendre_a(p,i) result(a)
    implicit none
    class(orthogonal_polynomial), intent(in) :: p
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
    class(orthogonal_polynomial), intent(in) :: p
    integer(i4b), intent(in) :: i
    real(dp) :: b
    b = 0.0_dp
    return
  end function legendre_b

  function legendre_c(p,i) result(c)
    implicit none
    class(orthogonal_polynomial), intent(in) :: p
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
  



  


  
end module module_special_functions
