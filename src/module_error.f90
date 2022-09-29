module module_error


contains

  subroutine error(check,routine,string)
    implicit none
    logical :: check
    character(len=*) :: routine,string
    if(check) then
       print *, 'Error in routine "',trim(routine),'": ', trim(string)       
       stop
    end if
    return
  end subroutine error


end module module_error
