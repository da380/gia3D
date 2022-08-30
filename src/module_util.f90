module module_util

  use module_constants
  use module_error


  type list_data
     integer(i4b) :: a
   contains
     procedure :: print_data
  end type list_data

  type, extends(list_data) :: extended_list_data
     real(dp) :: b
  end type extended_list_data


  
  type linked_list
     class(list_data), allocatable :: data
     type(linked_list), pointer :: next
  end type linked_list

  
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


  
  subroutine print_data(data)
    class(list_data), intent(in) :: data
    select type(data)
    type is(list_data)
       print *, data%a
    type is(extended_list_data)
       print *, data%a,data%b       
    end select

    return
  end subroutine print_data

  subroutine list_create( list, data )
    type(LINKED_LIST), pointer  :: list
    type(LIST_DATA), intent(in) :: data
    allocate( list )
    list%next => null()
    list%data =  data
  end subroutine list_create


  subroutine list_destroy( list )
    type(LINKED_LIST), pointer  :: list
    type(LINKED_LIST), pointer  :: current
    type(LINKED_LIST), pointer  :: next

    current => list
    do while ( associated(current%next) )
       next => current%next
        deallocate( current )
        current => next
    end do
  end subroutine list_destroy


  integer function list_count( list )
    type(LINKED_LIST), pointer  :: list
    type(LINKED_LIST), pointer  :: current
    type(LINKED_LIST), pointer  :: next
    
    if ( associated(list) ) then
       list_count = 1
       current => list
       do while ( associated(current%next) )
          current => current%next
          list_count = list_count + 1
       end do
    else
       list_count = 0
    end if
  end function list_count

  function list_next( elem ) result(next)
    type(LINKED_LIST), pointer :: elem
    type(LINKED_LIST), pointer :: next
    next => elem%next
  end function list_next

  subroutine list_insert( elem, data )
    type(LINKED_LIST), pointer  :: elem
    type(LIST_DATA), intent(in) :: data
    type(LINKED_LIST), pointer :: next
    allocate(next)
    next%next => elem%next
    elem%next => next
    next%data =  data
  end subroutine list_insert

  subroutine list_insert_head( list, data )
    type(LINKED_LIST), pointer  :: list
    type(LIST_DATA), intent(in) :: data
    type(LINKED_LIST), pointer :: elem
    allocate(elem)
    elem%data =  data
    elem%next => list
    list      => elem
  end subroutine list_insert_head


  
  subroutine list_insert_end( list, data )
    type(LINKED_LIST), pointer  :: list
    class(LIST_DATA), intent(in) :: data
    type(LINKED_LIST), pointer :: elem,current
    allocate(elem)
    elem%data =  data
    elem%next => null()
    current => list
    do while(associated(current%next))
       current => current%next
    end do
    current%next => elem
  end subroutine list_insert_end


  subroutine list_delete_element( list, elem )
    type(LINKED_LIST), pointer  :: list
    type(LINKED_LIST), pointer  :: elem

    type(LINKED_LIST), pointer  :: current
    type(LINKED_LIST), pointer  :: prev

    if ( associated(list,elem) ) then
        list => elem%next
        deallocate( elem )
    else
        current => list
        prev    => list
        do while ( associated(current) )
            if ( associated(current,elem) ) then
                prev%next => current%next
                deallocate( current ) ! Is also "elem"
                exit
            endif
            prev    => current
            current => current%next
        enddo
     endif
     !    allocate(next)
!
     !    next%next => elem%next
     !    elem%next => next
     !    next%data =  data
   end subroutine list_delete_element
  
  
end module module_util
