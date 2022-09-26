module module_util

  use module_constants
  use module_error
  implicit none

  interface count_columns
     procedure :: count_columns,count_columns_opened
  end interface count_columns


  interface found_command_argument
     procedure :: found_command_argument_character, &
                  found_command_argument_integer,   &
                  found_command_argument_real 
  end interface found_command_argument
  
contains

  !==============================================!
  !             polynomial evaluation            !
  !==============================================!
  
  function poly(n,coef,x) result(y)
    integer(i4b), intent(in) :: n
    real(dp), dimension(n), intent(in) :: coef
    real(dp), intent(in) :: X
    real(dp) :: y
    integer(i4b) :: i
    y = coef(n)
    do i = n-1,1,-1
       y = coef(i) + x*y
    end do
  end function poly


  !==============================================!
  !              search ordered list             !
  !==============================================!
  
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



  !==============================================!
  !                  IO routines                 !
  !==============================================!
  
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
    character(len=:), allocatable :: line
    integer, parameter :: buf_len= max_string_length
    character(len=buf_len) :: buf
    logical :: okay, set
    integer status, size    
    okay = .false.
    set = .true.
    do
       read (io,'(a)',advance='NO',iostat=status, size=size) buf
       okay = .not. IS_IOSTAT_END(status)
       if (.not. okay) return
       if (set) then
          line = buf(1:size)
          set=.false.
       else
          line = line // buf(1:size)
       end if
       if (IS_IOSTAT_EOR(status)) exit
    end do
  end function readline


  !==============================================!
  !             command line routines            !
  !==============================================!

  subroutine check_arguments(narg,larg,nopt,lopt) 
    integer(i4b), intent(in) :: narg    
    character(len=*) :: larg
    integer(i4b), intent(in), optional :: nopt
    character(len=*), intent(in), optional :: lopt
    character(len=:), allocatable :: loptl
    integer(i4b) :: n,noptl
    logical :: okay
    if(present(nopt) .and. present(lopt)) then
       noptl = nopt
       loptl = lopt
    else
       noptl = 0
       loptl = ''
    end if
    n = command_argument_count()
    okay = (n  >= 2*narg) .and. (n <= 2*(narg+noptl))
    if(.not.okay) then
       print *, 'Required: ',larg
       if(noptl > 0) then
          print *, 'Optional: ', lopt
       end if
       stop
    end if
    return
  end subroutine check_arguments
  
  logical function found_command_argument_character(tag,val) result(found)
    character(len=*), intent(in) :: tag
    character(len=:), allocatable, intent(out) :: val
    character(len=max_string_length) :: rtag,string
    integer(i4b) :: narg,iarg,jarg
    narg = command_argument_count()
    if(mod(narg,2) /= 0) then
       stop 'found_command_argument: number of arguments needs to be divisible by two'
    end if
    found = .false.
    do iarg = 1,narg/2
       jarg = 2*iarg-1
       call get_command_argument(jarg,rtag)
       call get_command_argument(jarg+1,string)
       if(trim(rtag) == tag) then
          found = .true.
          val = trim(string)
          return
       end if
    end do
    return
  end function found_command_argument_character

  logical function found_command_argument_integer(tag,val) result(found)
    character(len=*), intent(in) :: tag
    integer(i4b), intent(out) :: val
    character(len=max_string_length) :: rtag,string
    integer(i4b) :: narg,iarg,jarg,ios
    narg = command_argument_count()
    if(mod(narg,2) /= 0) then
       stop 'found_command_argument: number of arguments needs to be divisible by two'
    end if
    found = .false.
    do iarg = 1,narg/2
       jarg = 2*iarg-1
       call get_command_argument(jarg,rtag)
       call get_command_argument(jarg+1,string)
       if(trim(rtag) == tag) then
          read(string,*,iostat = ios) val
          if(ios == 0)  then
             found = .true.
             return
          else
             found = .false.
             return
          end if
       end if
    end do
    return
  end function found_command_argument_integer


  logical function found_command_argument_real(tag,val) result(found)
    character(len=*), intent(in) :: tag
    real(dp), intent(out) :: val
    character(len=max_string_length) :: rtag,string
    integer(i4b) :: narg,iarg,jarg,ios
    narg = command_argument_count()
    if(mod(narg,2) /= 0) then
       stop 'found_command_argument: number of arguments needs to be divisible by two'
    end if
    found = .false.
    do iarg = 1,narg/2
       jarg = 2*iarg-1
       call get_command_argument(jarg,rtag)
       call get_command_argument(jarg+1,string)
       if(trim(rtag) == tag) then
          read(string,*,iostat = ios) val
          if(ios == 0)  then
             found = .true.
             return
          else
             found = .false.
             return
          end if
       end if
    end do
    return
  end function found_command_argument_real


  
end module module_util
