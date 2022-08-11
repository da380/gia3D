module module_fftw3

  use module_constants
  use module_error
  use, intrinsic :: iso_c_binding
  include 'fftw3.f03'

  
  interface dft_1D
     module procedure :: dft_complex_scalar_1D,       &
                         dft_complex_vector_1D,       &
                         dft_real_scalar_1D_forward,  &
                         dft_real_scalar_1D_backward, &
                         dft_real_vector_1D_forward,  &
                         dft_real_vector_1D_backward     
  end interface dft_1D

contains

  subroutine dft_real_scalar_1D_forward(n,in,out)
    ! simple driver routine for performing an FFT of
    ! a  real scalar using FFTW3.
    !
    ! Variables:
    !
    ! n       -- integer. Number of sample points. This does not need
    !                     to be a power of 2, but the routines are
    !                     fastest if it is.
    !
    ! in      -- real. Array to be transformed. Dimension n
    !
    ! out     -- complex. Array containing the result. Dimension n/2+1
    !
    ! Note that the backward transformation does NOT include the
    ! 1/n normalisation.    
    implicit none
    integer(i4b), intent(in) :: n
    real(dp), dimension(n), intent(in) :: in
    complex(dpc), dimension(n/2+1), intent(out) :: out

    type(C_PTR) :: plan
    real(C_DOUBLE), pointer :: lin(:)
    complex(C_DOUBLE_COMPLEX), pointer :: lout(:)
    type(C_PTR) :: pin,pout

    ! set up the C pointers
    pin  = fftw_alloc_real(int(n, C_SIZE_T))
    pout = fftw_alloc_complex(int(n/2+1, C_SIZE_T))
    call c_f_pointer(pin,   lin, [n])
    call c_f_pointer(pout, lout, [n/2+1])

    ! copy the input data
    lin = in

    ! make the plan
    plan = fftw_plan_dft_r2c_1d(n,lin,lout,FFTW_ESTIMATE);

    ! do the transform
    call fftw_execute_dft_r2c(plan, lin, lout)

    ! copy the result
    out = lout

    ! delete the pointers
    call fftw_free(pin)
    call fftw_free(pout)
    lin  => null()
    lout => null()
    
    return
  end subroutine dft_real_scalar_1D_forward


  subroutine dft_real_vector_1D_forward(n,m,in,out)
    ! simple driver routine for performing an FFT of
    ! a  real scalar using FFTW3.
    !
    ! Variables:
    !
    ! n       -- integer. Number of sample points. This does not need
    !                     to be a power of 2, but the routines are
    !                     fastest if it is.
    !
    ! m       -- integer. Number of vector components. Note that the
    !
    ! in      -- real. Array to be transformed. Dimension n
    !
    ! out     -- complex. Array containing the result. Dimension n/2+1
    !
    ! Note that the backward transformation does NOT include the
    ! 1/n normalisation.    
    implicit none
    integer(i4b), intent(in) :: n
    integer(i4b), intent(in) :: m
    real(dp), dimension(n,m), intent(in) :: in
    complex(dpc), dimension(n/2+1,m), intent(out) :: out

    integer(i4b) :: j
    type(C_PTR) :: plan
    real(C_DOUBLE), pointer :: lin(:)
    complex(C_DOUBLE_COMPLEX), pointer :: lout(:)
    type(C_PTR) :: pin,pout

    ! set up the C pointers
    pin  = fftw_alloc_real(int(n, C_SIZE_T))
    pout = fftw_alloc_complex(int(n/2+1, C_SIZE_T))
    call c_f_pointer(pin,   lin, [n])
    call c_f_pointer(pout, lout, [n/2+1])

   
    ! make the plan
    plan = fftw_plan_dft_r2c_1d(n,lin,lout,FFTW_ESTIMATE);

    do j = 1,m

       ! copy the input data
       lin = in(:,j)

       ! do the transform
       call fftw_execute_dft_r2c(plan, lin, lout)

       ! copy the result
       out(:,j) = lout
    end do

    ! delete the pointers
    call fftw_free(pin)
    call fftw_free(pout)
    lin  => null()
    lout => null()
    
    return
  end subroutine dft_real_vector_1D_forward


  subroutine dft_real_scalar_1D_backward(n,in,out)
    ! simple driver routine for performing an inverse FFT of
    ! a  real vector using FFTW3.
    !
    ! Variables:
    !
    ! n       -- integer. Number of sample points. This does not need
    !                     to be a power of 2, but the routines are
    !                     fastest if it is.
    !
    ! in      -- complex. Array to be transformed. Dimension n/2+1
    !
    ! out     -- real. Array containing the result. Dimension n
    !
    ! Note that the backward transformation does NOT include the
    ! 1/n normalisation.    
    implicit none
    integer(i4b), intent(in) :: n
    complex(dpc), dimension(n/2+1), intent(in) :: in
    real(dp), dimension(n), intent(out) :: out

    integer(i4b) :: j
    type(C_PTR) :: plan
    complex(C_DOUBLE_COMPLEX), pointer :: lin(:)
    real(C_DOUBLE), pointer :: lout(:)
    type(C_PTR) :: pin,pout

    ! set up the C pointers
    pin = fftw_alloc_complex(int(n/2+1, C_SIZE_T))
    pout  = fftw_alloc_real(int(n, C_SIZE_T))
    call c_f_pointer(pin,   lin, [n/2+1])
    call c_f_pointer(pout, lout, [n])

    ! make the plan
    plan = fftw_plan_dft_c2r_1d(n,lin,lout,FFTW_ESTIMATE);

    ! copy the input data
    lin = in

    ! do the transform
    call fftw_execute_dft_c2r(plan, lin, lout)

    out = lout

    ! delete the pointers
    call fftw_free(pin)
    call fftw_free(pout)
    lin  => null()
    lout => null()
    
    return
  end subroutine dft_real_scalar_1D_backward
  

  subroutine dft_real_vector_1D_backward(n,m,in,out)
    ! simple driver routine for performing an inverse FFT of
    ! a  real vector using FFTW3.
    !
    ! Variables:
    !
    ! n       -- integer. Number of sample points. This does not need
    !                     to be a power of 2, but the routines are
    !                     fastest if it is.
    !
    ! m       -- integer. Number of vector components. Note that the
    !
    ! in      -- complex. Array to be transformed. Dimension n/2+1
    !
    ! out     -- real. Array containing the result. Dimension n
    !
    ! Note that the backward transformation does NOT include the
    ! 1/n normalisation.    
    implicit none
    integer(i4b), intent(in) :: n
    integer(i4b), intent(in) :: m
    complex(dpc), dimension(n/2+1,m), intent(in) :: in
    real(dp), dimension(n,m), intent(out) :: out

    integer(i4b) :: j
    type(C_PTR) :: plan
    complex(C_DOUBLE_COMPLEX), pointer :: lin(:)
    real(C_DOUBLE), pointer :: lout(:)
    type(C_PTR) :: pin,pout

    ! set up the C pointers
    pin = fftw_alloc_complex(int(n/2+1, C_SIZE_T))
    pout  = fftw_alloc_real(int(n, C_SIZE_T))
    call c_f_pointer(pin,   lin, [n/2+1])
    call c_f_pointer(pout, lout, [n])


    ! make the plan
    plan = fftw_plan_dft_c2r_1d(n,lin,lout,FFTW_ESTIMATE);

    do j = 1,m

       ! copy the input data
       lin = in(:,j)

       ! do the transform
       call fftw_execute_dft_c2r(plan, lin, lout)

       ! copy the result
       out(:,j) = lout
       
    end do

    ! delete the pointers
    call fftw_free(pin)
    call fftw_free(pout)
    lin  => null()
    lout => null()
    
    return
  end subroutine dft_real_vector_1D_backward
  

  subroutine dft_complex_scalar_1D(forward,n,in,out)
    ! simple driver routine for performing an FFT of
    ! a complex scalar using FFTW3.
    !
    ! Variables:
    !
    ! forward -- logical. Equals .true. for a forward transformation
    !                     (meaning + sign in exponent) or .false.
    !                     for a backward transformation.
    !
    ! n       -- integer. Number of sample points. This does not need
    !                     to be a power of 2, but the routines are
    !                     fastest if it is.
    !
    ! in      -- complex. Array to be transformed.
    !
    ! out     -- complex. Array containing the result
    !
    ! Note that the backward transformation does NOT include the
    ! 1/n normalisation.    
    implicit none
    logical, intent(in) :: forward
    integer(i4b), intent(in) :: n
    complex(dpc), dimension(n), intent(in) :: in
    complex(dpc), dimension(n), intent(out) :: out

    type(C_PTR) :: plan
    complex(C_DOUBLE_COMPLEX), pointer :: lin(:),lout(:)
    type(C_PTR) :: pin,pout

    ! set up the C pointers
    pin  = fftw_alloc_complex(int(n, C_SIZE_T))
    pout = fftw_alloc_complex(int(n, C_SIZE_T))
    call c_f_pointer(pin,   lin, [n])
    call c_f_pointer(pout, lout, [n])

    ! copy over the input data
    lin  = in

    ! make the plan
    if(forward) then
       plan = fftw_plan_dft_1d(n, lin, lout, FFTW_FORWARD, FFTW_ESTIMATE)
    else
       plan = fftw_plan_dft_1d(n, lin, lout, FFTW_BACKWARD, FFTW_ESTIMATE)
    end if

    ! do the transform
    call fftw_execute_dft(plan, lin, lout)

    ! copy the result
    out = lout

    ! delete the pointers
    call fftw_free(pin)
    call fftw_free(pout)
    lin  => null()
    lout => null()

    
    return
  end subroutine dft_complex_scalar_1D



  subroutine dft_complex_vector_1D(forward,n,m,in,out)
    ! simple driver routine for performing an FFT of
    ! a complex vector using FFTW3.
    !
    ! Variables:
    !
    ! forward -- logical. Equals .true. for a forward transformation
    !                     (meaning + sign in exponent) or .false.
    !                     for a backward transformation.
    !
    ! n       -- integer. Number of sample points. This does not need
    !                     to be a power of 2, but the routines are
    !                     fastest if it is.
    !
    ! m       -- integer. Number of vector components. Note that the
    !
    ! in      -- complex. Array to be transformed. Note that each
    !                     column denotes a component of the vector.
    !
    ! out     -- complex. Array containing the result
    !
    ! Note that the backward transformation does NOT include the
    ! 1/n normalisation.
    implicit none
    logical, intent(in) :: forward
    integer(i4b), intent(in) :: n,m
    complex(dpc), dimension(n,m), intent(in) :: in
    complex(dpc), dimension(n,m), intent(out) :: out       

    integer(i4b) :: j 
    type(C_PTR) :: plan
    complex(C_DOUBLE_COMPLEX), pointer :: lin(:),lout(:)
    type(C_PTR) :: pin,pout
    
    ! set up the C pointers
    pin  = fftw_alloc_complex(int(n, C_SIZE_T))
    pout = fftw_alloc_complex(int(n, C_SIZE_T))
    call c_f_pointer(pin,   lin, [n])
    call c_f_pointer(pout, lout, [n])

    ! make the plan
    if(forward) then
       plan = fftw_plan_dft_1d(n, lin, lout, FFTW_FORWARD, FFTW_ESTIMATE)
    else
       plan = fftw_plan_dft_1d(n, lin, lout, FFTW_BACKWARD, FFTW_ESTIMATE)
    end if

    do j = 1,m

       lin = in(:,j)

       ! do the transform
       call fftw_execute_dft(plan,lin,lout)
       
       ! copy over the result
       out(:,j) = lout
       
    end do

    ! delete the pointers
    call fftw_free(pin)
    call fftw_free(pout)
    lin  => null()
    lout => null()
    
    return
  end subroutine dft_complex_vector_1D




  
end module module_fftw3
