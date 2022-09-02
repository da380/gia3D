module module_util

  use module_constants
  use module_error

  interface count_columns
     procedure :: count_columns,count_columns_opened
  end interface count_columns

  
contains

  function poly(n,coef,x) result(y)
    real(dp), dimension(n), intent(in) :: coef
    real(dp), intent(in) :: X
    real(dp) :: y
    integer(i4b) :: i
    y = coef(n)
    do i = n-1,1,-1
       y = coef(i) + x*y
    end do
  end function poly


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


  function count_columns(file) result(ncol)
    implicit none
    character(len = *), intent(in) :: file
    integer(i4b) :: ncol

    character(len=:), allocatable :: line
    integer(i4b) :: io,ios
    real(dp), dimension(:), allocatable :: tmp
    open(newunit = io,file = trim(file))
    line = readline(io)
    allocate(tmp(len(line)))
    close(io)
    ncol = 0
    do
       read(line,*,iostat = ios) tmp(1:ncol)
       if(ios /= 0) exit
       ncol = ncol + 1
    end do 
    ncol = ncol-1    
    return
  end function count_columns

  function count_columns_opened(io) result(ncol)
    implicit none
    integer(i4b), intent(in) :: io
    integer(i4b) :: ncol

    logical :: ok
    character(len=:), allocatable :: line
    integer(i4b) :: ios
    real(dp), dimension(:), allocatable :: tmp
    line = readline(io)
    allocate(tmp(len(line)))
    backspace(io)
    ncol = 0
    do
       read(line,*,iostat = ios) tmp(1:ncol)
       if(ios /= 0) exit
       ncol = ncol + 1
    end do 
    ncol = ncol-1
    
    return
  end function count_columns_opened

  
  function readline(io) result(line)
    implicit none
    integer, intent(IN) :: io
    character(LEN=:), allocatable :: line
    integer, parameter :: line_buf_len= 1024*4
    character(LEN=line_buf_len) :: InS
    logical :: OK, set
    integer status, size    
    OK = .false.
    set = .true.
    do
       read (io,'(a)',advance='NO',iostat=status, size=size) InS
       OK = .not. IS_IOSTAT_END(status)
       if (.not. OK) return
       if (set) then
          line = InS(1:size)
          set=.false.
       else
          line = line // InS(1:size)
       end if
       if (IS_IOSTAT_EOR(status)) exit
    end do
  end function readline

end module module_util
