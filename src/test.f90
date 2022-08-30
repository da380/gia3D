program test

  use module_constants
  use module_util

  
  integer(i4b) :: i
  type(linked_list), pointer :: list,current,next
  type(list_data) :: data
  type(extended_list_data) :: extended_data

  
  data%a = 0
  call list_create(list,data)
  current => list
  do while(associated(current))
     call current%data%print_data()
     current => current%next
  end do
  print *, '======================'
  
  do i = 1,5
     data%a = data%a + 1
     call list_insert_end( list, data )
  end do

  extended_data%a = data%a
  extended_data%b = 0.0_dp
  do i = 1,5
     extended_data%a = extended_data%a + 1
     extended_data%b = extended_data%b + 1
     call list_insert_end( list, extended_data )
  end do

  

  
  current => list
  do while(associated(current))
     call current%data%print_data()
     current => current%next
  end do
  print *, '======================'

  call list_destroy( list )

  
end program test
