program fftw3_test

  use module_constants
  use module_fftw3
  implicit none

  integer(i4b), parameter :: n = 256, m = 3
  integer(i4b) :: i,j
  integer(i4b) :: rank = 1
  integer(i4b), dimension(1) :: nn 
  integer(i4b) :: howmany = m
  integer(i4b) :: idist = 1,odist = 1
  integer(i4b) :: istride = m,ostride = m
  integer(i4b), dimension(1) :: inembed,onembed

  complex(C_DOUBLE_COMPLEX), pointer :: in(:,:),out(:,:),rec(:,:)
  type(C_PTR) :: pin,pout,prec
  type(C_PTR) :: plan
  
  nn(1) = n
  inembed = nn
  onembed = nn
  
  ! allocate the arrays
  pin  = fftw_alloc_complex(int(n*m, C_SIZE_T))
  pout = fftw_alloc_complex(int(n*m, C_SIZE_T))
  prec = fftw_alloc_complex(int(n*m, C_SIZE_T))
  call c_f_pointer(pin, in, [n,m])
  call c_f_pointer(pout, out, [n,m])
  call c_f_pointer(prec, rec, [n,m])


  ! set the forward array
  do i = 1,n
     do j = 1,m
        in(i,j) = cos(3*2*pi*(i+j)/n)        
     end do     
  end do


  ! make the forward plan
  plan =  fftw_plan_many_dft(rank, nn,  howmany,   &
                             in, inembed, istride, & 
			     idist, out, onembed,  & 
			     ostride, odist, FFTW_FORWARD, FFTW_ESTIMATE)
  

  ! exucute the plan
  call fftw_execute_dft(plan,in,out)



  ! make the backwards plan
  plan =  fftw_plan_many_dft(rank, nn,  howmany,   &
                             out, inembed, istride, & 
			     idist, rec, onembed,  & 
			     ostride, odist, FFTW_BACKWARD, FFTW_ESTIMATE)

  ! exucute the plan
  call fftw_execute_dft(plan,out,rec)



  open(99,file='fftw_test.out')

  do i = 1,n
     write(99,*) i,real(in(i,:)),real(in(i,:))-real(rec(i,:))/n
  end do
  
  close(99)
  
!  integer, parameter :: n = 4
!  integer :: i
!  type(C_PTR) :: plan
!  complex(C_DOUBLE_COMPLEX), pointer :: in(:),out(:)
!  type(C_PTR) :: pin,pout


!  pin  = fftw_alloc_complex(int(n, C_SIZE_T))
!  pout = fftw_alloc_complex(int(n, C_SIZE_T))
!  call c_f_pointer(pin, in, [n])
!  call c_f_pointer(pout, out, [n])
!  plan = fftw_plan_dft_1d(n, in, out, FFTW_FORWARD, FFTW_ESTIMATE)


!  do i = 1,n
!     in(i) = sin(pi*i)
!  end do

!  print *, real(in)

!  call fftw_execute_dft(plan, in, out)


!  plan = fftw_plan_dft_1d(n, out, in, FFTW_BACKWARD, FFTW_ESTIMATE)

!  call fftw_execute_dft(plan, out, in)

!  print *, real(in/n)
  
  
  
end program fftw3_test
