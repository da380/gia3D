module module_error


contains

  subroutine check(test,routine,string)
    implicit none
    logical :: test
    character(len=*) :: routine,string
    if(.not.test) then
       print *, 'Error in routine "',trim(routine),'": ', trim(string)       
       stop
    end if
    return
  end subroutine check


end module module_error
