module module_util

  use module_constants
  use module_error

contains

  function poly_eval(n,coef,x) result(y)
    real(dp), dimension(n), intent(in) :: coef
    real(dp), intent(in) :: X
    real(dp) :: y
    integer(i4b) :: i
    y = coef(n)
    do i = n-1,1,-1
       y = coef(i) + x*y
    end do
  end function poly_eval


  function bisect_list(xx,x) result(il)
    implicit none
    real(dp), dimension(:), intent(in) :: xx
    real(dp), intent(in) :: x
    integer(i4b) :: il
    
    logical :: ascnd
    integer(i4b) :: n,im,iu
    
    n = size(xx)
    ascnd = xx(n) > xx(1)      
    il = 1
    iu = n
    do while(iu-il > 1)
       im = (il+iu)/2       
       if((x >= xx(im)) .eqv.  ascnd) then
          il = im
       else
          iu = im
       endif
    end do
    
    return
  end function bisect_list


  function hunt_list(xx,x,it) result(il)
    implicit none
    real(dp), dimension(:), intent(in) :: xx
    real(dp), intent(in) :: x
    integer(i4b), intent(in) :: it
    integer(i4b) :: il

    logical :: ascnd
    integer(i4b) :: n,im,iu,inc
    
    n = size(xx)
    ascnd = xx(n) > xx(1)      
    if(it < 1 .or. it >= n) then
       il =  bisect_list(xx,x)
       return
    end if

    il = it
    inc = 1

    if(x >= xx(il) .eqv. ascnd) then
       
       do
          iu = min(il + inc,n)
          if(iu == n) exit
          if(x < xx(iu) .eqv. ascnd) exit
          il = iu
          inc = inc+inc
       end do
       il = il + bisect_list(xx(il:iu),x)-1
       
    else
       
       iu = il
       do
          il = max(iu-inc,1)
          if(il == 1) exit
          if(x > xx(il) .eqv. ascnd) exit
          iu = il
          inc = inc+inc          
       end do
       il = il + bisect_list(xx(il:iu),x)-1
       
    end if
    
    return
  end function hunt_list

    
  
end module module_util
