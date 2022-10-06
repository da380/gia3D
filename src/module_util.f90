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
                  found_command_argument_real,      &
                  found_command_argument_logical,   &
                  found_command_argument_flag 
  end interface found_command_argument


  interface normal_random_variable
     procedure :: single_normal_random_variable, &
                  pair_normal_random_variable,    &
                  vector_normal_random_variable
  end interface normal_random_variable
  
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

  subroutine print_argument_info(nreq,ireq,dreq,nopt,iopt,dopt) 
    integer(i4b), intent(in) :: nreq
    character(len=*), dimension(nreq), intent(in) :: ireq
    logical, dimension(nreq), intent(in) :: dreq
    integer(i4b), intent(in), optional :: nopt
    character(len=*), dimension(:), intent(in), optional :: iopt
    logical, dimension(:), intent(in), optional :: dopt

    integer(i4b) :: i,nc
    
    ! work out the minimum number of command line arguments
    nc = 0
    do i = 1,nreq
       if(dreq(i)) then
          nc = nc+2
       else
          nc = nc+1
       end if
    end do

    if(command_argument_count() < nc) then
       do i = 1,nreq
          print *, trim(ireq(i))
       end do
       if(present(nopt) .and. present(iopt) .and. present(dopt)) then
          do i = 1,nopt
             print *, trim(iopt(i))
          end do
       end if
       stop
    end if
  end subroutine print_argument_info

  logical function found_command_argument_flag(tag) result(found)
    character(len=*), intent(in) :: tag
    character(len=max_string_length) :: rtag,string
    integer(i4b) :: narg,iarg,jarg,ios
    narg = command_argument_count()
    found = .false.
    do iarg = 1,narg
       call get_command_argument(iarg,rtag)
       if(trim(rtag) == tag) then
          found = .true.          
          return
       end if
    end do
    return
  end function found_command_argument_flag

  
  logical function found_command_argument_character(tag,val) result(found)
    character(len=*), intent(in) :: tag
    character(len=:), allocatable, intent(out) :: val
    character(len=max_string_length) :: rtag,string
    integer(i4b) :: narg,iarg,ios
    narg = command_argument_count()
    found = .false.
    do iarg = 1,narg-1
       call get_command_argument(iarg,rtag)
       if(trim(rtag) == tag) then
          call get_command_argument(iarg+1,string,status = ios)
          if(ios /= 0) return
          found = .true.
          val = trim(string)         
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
    found = .false.
    do iarg = 1,narg-1
       call get_command_argument(iarg,rtag)
       if(trim(rtag) == tag) then
          call get_command_argument(iarg+1,string,status = ios)
          if(ios /= 0) return
          read(string,*,iostat = ios) val
          if(ios /= 0)  return
          found = .true.
          return
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
    found = .false.
    do iarg = 1,narg-1
       call get_command_argument(iarg,rtag)
       if(trim(rtag) == tag) then
          call get_command_argument(iarg+1,string,status = ios)
          if(ios /= 0) return
          read(string,*,iostat = ios) val
          if(ios /= 0)  return
          found = .true.
          return
       end if
    end do
    return
  end function found_command_argument_real


  logical function found_command_argument_logical(tag,val) result(found)
    character(len=*), intent(in) :: tag
    logical, intent(out) :: val
    character(len=max_string_length) :: rtag,string
    integer(i4b) :: narg,iarg,jarg,ios
    narg = command_argument_count()
    found = .false.
    do iarg = 1,narg-1
       call get_command_argument(iarg,rtag)
       if(trim(rtag) == tag) then
          call get_command_argument(iarg+1,string,status = ios)
          if(ios /= 0) return
          read(string,*,iostat = ios) val
          if(ios /= 0)  return
          found = .true.
          return
       end if
    end do
    return
  end function found_command_argument_logical
  


  !==============================================!
  !             command line routines            !
  !==============================================!
  

  subroutine pair_normal_random_variable(ran1,ran2)
    real(dp), intent(out) :: ran1,ran2
    real(dp) :: u1,u2,r
    call random_number(u1)
    call random_number(u2)
    r = -2.0_dp*log(u1)
    if(r > 0.0_dp) then
       r = sqrt(r)
    else
       print *, 'hello'
       r = 0.0_dp
    end if
    ran1 = r*cos(twopi*u2)
    ran2 = r*sin(twopi*u2)    
    return
  end subroutine pair_normal_random_variable
  
  subroutine single_normal_random_variable(ran)
    real(dp), intent(out) :: ran
    real(dp) :: u1,u2,r
    call random_number(u1)
    call random_number(u2)
    r = -2.0_dp*log(u1)
    if(r > 0.0_dp) then
       r = sqrt(r)
    else
       print *, 'hello'
       r = 0.0_dp
    end if
    ran = r*cos(twopi*u2)    
    return
  end subroutine single_normal_random_variable

  subroutine vector_normal_random_variable(ran)
    real(dp), dimension(:), intent(out) :: ran
    integer(i4b) :: n,m,i
    real(dp) :: u1,u2,z1,z2
    n = size(ran)
    m = n/2
    do i = 1,m
       call pair_normal_random_variable(z1,z2)
       ran(i) = z1
       ran(m+i) = z2
    end do
    if(2*m < n) then
       call single_normal_random_variable(z1)
       ran(n) = z1
    end if
    return
  end subroutine vector_normal_random_variable



  
end module module_util
