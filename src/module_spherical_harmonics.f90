module module_spherical_harmonics
  
  use module_constants
  use module_special_functions
  use, intrinsic :: iso_c_binding


  !===============================================================!
  !               type declaration for the GL-grid                !
  !===============================================================!
  
  type gauss_legendre_grid
     logical :: allocated = .false.
     integer(i4b) :: lmax
     integer(i4b) :: nmax
     integer(i4b) :: nth
     integer(i4b) :: nph
     real(dp), dimension(:), allocatable :: th
     real(dp), dimension(:), allocatable :: w
     real(dp), dimension(:,:,:), allocatable :: dlm
     type(C_PTR) :: plan_forward
     type(C_PTR) :: plan_backward
     type(C_PTR) :: plan_r2c
     type(C_PTR) :: plan_c2r     
   contains
     procedure :: delete => delete_gauss_legendre_grid
     procedure :: allocate =>  allocate_gauss_legendre_grid
     procedure :: ph =>  ph_gauss_legednre_grid
     procedure :: SH_trans_gauss_legendre_grid
     procedure :: wrapper_SH_trans_gauss_legendre_grid
     generic :: SH_trans => SH_trans_gauss_legendre_grid,        &
                            wrapper_SH_trans_gauss_legendre_grid
     procedure :: SH_itrans_gauss_legendre_grid
     procedure :: wrapper_SH_itrans_gauss_legendre_grid
     generic :: SH_itrans => SH_itrans_gauss_legendre_grid,        &
                            wrapper_SH_itrans_gauss_legendre_grid
     
  end type gauss_legendre_grid


  !===============================================================!
  !              type declarations for the GL-fields              !
  !===============================================================!

  type, abstract :: gauss_legendre_field
     logical :: allocated = .false.
     integer(i4b) :: lmax
     integer(i4b) :: nth
     integer(i4b) :: nph
     integer(i4b) :: nc = 0
     integer(i4b) :: cdim
     integer(i4b), dimension(:), allocatable :: cnval
     complex(dpc), dimension(:), allocatable :: cdata
     integer(i4b) :: nr = 0
     integer(i4b) :: rdim
     real(dp), dimension(:), allocatable :: rdata
     
   contains
     procedure :: delete => delete_gauss_legendre_field
     procedure :: index =>    index_gauss_legendre_field
     procedure :: scale_gauss_legendre_field
     procedure :: real_scale_gauss_legendre_field
     generic   :: scale => scale_gauss_legendre_field,     &
                           real_scale_gauss_legendre_field
     procedure :: saxpy_gauss_legendre_field
     procedure :: real_saxpy_gauss_legendre_field          
     generic   :: saxpy => saxpy_gauss_legendre_field,     &
                           real_saxpy_gauss_legendre_field
  end type gauss_legendre_field

  type, extends(gauss_legendre_field) :: scalar_gauss_legendre_field
     logical :: real = .false.
   contains
     procedure :: allocate => allocate_scalar_gauss_legendre_field
     procedure :: set_complex_scalar_gauss_legendre_field
     procedure :: set_real_scalar_gauss_legendre_field
     generic   :: set => set_complex_scalar_gauss_legendre_field, &
                         set_real_scalar_gauss_legendre_field
     procedure :: getc => get_complex_scalar_gauss_legendre_field
     procedure :: getr => get_real_scalar_gauss_legendre_field
  end type scalar_gauss_legendre_field


  type, extends(gauss_legendre_field) :: vector_gauss_legendre_field
     logical :: real = .false.
   contains
     procedure :: allocate => allocate_vector_gauss_legendre_field
     procedure :: set => set_value_vector_gauss_legendre_field
     procedure :: get => get_value_vector_gauss_legendre_field
  end type vector_gauss_legendre_field

  
  type, extends(gauss_legendre_field) :: second_tensor_gauss_legendre_field
     logical :: symmetric  = .false.
     logical :: trace_free = .false.
     logical :: real       = .false.
   contains
     procedure :: allocate => allocate_second_tensor_gauss_legendre_field
     procedure :: set => set_value_second_tensor_gauss_legendre_field
     procedure :: get => get_value_second_tensor_gauss_legendre_field
  end type second_tensor_gauss_legendre_field


  
  !===============================================================!
  !                 type declarations SH expansions               !
  !===============================================================!
  
  type, abstract :: spherical_harmonic_expansion
     logical :: allocated
     integer(i4b) :: lmax
     real(dp) :: s = 0.0_dp
     real(dp) :: mu = 1.0_dp
     integer(i4b) :: nc = 0
     integer(i4b) :: cncoef
     integer(i4b) :: cdim
     integer(i4b), dimension(:), allocatable :: cnval     
     complex(dpc), dimension(:), allocatable :: cdata
     integer(i4b) :: nr = 0
     integer(i4b) :: rncoef
     integer(i4b) :: rdim
     complex(dpc), dimension(:), allocatable :: rdata
   contains
     procedure :: delete    => delete_spherical_harmonic_expansion
     procedure :: cindex     => cindex_spherical_harmonic_expansion
     procedure :: rindex     => rindex_spherical_harmonic_expansion
     procedure :: set_sobolev => set_sobolev_spherical_harmonic_expansion
     procedure :: scale_spherical_harmonic_expansion
     procedure :: real_scale_spherical_harmonic_expansion
     generic   :: scale => scale_spherical_harmonic_expansion,     &
                           real_scale_spherical_harmonic_expansion
     procedure :: saxpy_spherical_harmonic_expansion
     procedure :: real_saxpy_spherical_harmonic_expansion          
     generic   :: saxpy => saxpy_spherical_harmonic_expansion,     &
                           real_saxpy_spherical_harmonic_expansion
  end type spherical_harmonic_expansion

  type, extends(spherical_harmonic_expansion) :: scalar_spherical_harmonic_expansion
     logical :: real = .false.
   contains
     procedure :: allocate => allocate_scalar_spherical_harmonic_expansion
     procedure :: set => set_complex_scalar_spherical_harmonic_expansion
     procedure :: get => get_complex_scalar_spherical_harmonic_expansion
  end type scalar_spherical_harmonic_expansion

  
  type, extends(spherical_harmonic_expansion) :: vector_spherical_harmonic_expansion
     logical :: real = .false.
   contains
     procedure :: allocate => allocate_vector_spherical_harmonic_expansion
     procedure :: set => set_value_vector_spherical_harmonic_expansion
     procedure :: get => get_value_vector_spherical_harmonic_expansion
  end type vector_spherical_harmonic_expansion

  
  type, extends(spherical_harmonic_expansion) :: second_tensor_spherical_harmonic_expansion
     logical :: symmetric  = .false.
     logical :: trace_free = .false.
     logical :: real       = .false.
   contains
     procedure :: allocate => allocate_second_tensor_spherical_harmonic_expansion
     procedure :: set => set_value_second_tensor_spherical_harmonic_expansion
     procedure :: get => get_value_second_tensor_spherical_harmonic_expansion
  end type second_tensor_spherical_harmonic_expansion
  
     
contains

  !=======================================================================!
  !                     procedures for the grid type                      !
  !=======================================================================!


  !--------------------------------------!
  !             basic routines           !
  !--------------------------------------!
  
  subroutine delete_gauss_legendre_grid(grid)
    use module_fftw3
    implicit none
    class(gauss_legendre_grid), intent(inout) :: grid
    if(.not.grid%allocated) return
    deallocate(grid%th,grid%w,grid%dlm)
    call fftw_destroy_plan(grid%plan_forward)
    call fftw_destroy_plan(grid%plan_backward)
    grid%allocated = .false.
    return
  end subroutine delete_gauss_legendre_grid

  
  subroutine allocate_gauss_legendre_grid(grid,lmax,nmax,fftw_flag)
    use module_special_functions
    use module_quadrature
    use module_fftw3    
    implicit none    
    class(gauss_legendre_grid), intent(inout) :: grid
    integer(i4b), intent(in) :: lmax
    integer(i4b), intent(in) :: nmax
    integer(C_INT), intent(in), optional :: fftw_flag

    
    integer(i4b) :: l,ith,n,nth,nph,cdim
    real(dp) :: fac
    class(orthogonal_polynomial), allocatable :: poly
    type(gauss_quadrature) :: quad
    real(C_DOUBLE), pointer :: rin(:)
    complex(C_DOUBLE_COMPLEX), pointer ::  in(:)
    complex(C_DOUBLE_COMPLEX), pointer :: out(:),rout(:)
    type(C_PTR) :: pin,pout,plan1,plan2,plan3,plan4
    integer(C_INT) :: plan_flag

    ! deal with optional arguments
    if(present(fftw_flag)) then
       plan_flag = fftw_flag
    else
       plan_flag = FFTW_MEASURE
    end if

    ! check it has not already been allocated
    call grid%delete()

    ! store the basic parameters
    grid%lmax = lmax
    grid%nmax = nmax
    nth = lmax+1
    nph = 2*lmax
    grid%nth = nth
    grid%nph = nph

    
    ! make the quadrature points and weights
    allocate(grid%th(lmax+1))
    allocate(grid%w(lmax+1))
    poly = legendre()
    call quad%set(lmax+1,poly)
    grid%th = acos(quad%points())
    grid%w = quad%weights()


    ! make the wigner d-functions
    allocate(grid%dlm((lmax+1)*(lmax+2)/2,-nmax:nmax,nth))
    do ith = 1,nth
       call set_wigner_array(grid%th(ith),lmax,nmax,grid%dlm(:,:,ith),norm=.true.)
    end do

    !-----------------------------!
    !    make the FFTW3 plans     !
    !-----------------------------!
    n = nph

    ! complex transforms
    pin  = fftw_alloc_complex(int(n, C_SIZE_T))
    pout = fftw_alloc_complex(int(n, C_SIZE_T))
    call c_f_pointer(pin,   in, [n])
    call c_f_pointer(pout, out, [n])
    plan1 = fftw_plan_dft_1d(n, in, out, FFTW_FORWARD, plan_flag)
    plan2 = fftw_plan_dft_1d(n, out, in, FFTW_BACKWARD, plan_flag)
    call fftw_free(pin)
    call fftw_free(pout)

    ! real transforms
    pin  = fftw_alloc_real(int(n, C_SIZE_T))
    pout = fftw_alloc_complex(int(n/2+1, C_SIZE_T))
    call c_f_pointer(pin,   rin, [n])
    call c_f_pointer(pout, rout, [n/2+1])
    plan3 = fftw_plan_dft_r2c_1d(n,rin,rout,plan_flag);
    plan4 = fftw_plan_dft_c2r_1d(n,rout,rin,plan_flag);
    call fftw_free(pin)
    call fftw_free(pout)
    
    ! store the plans
    grid%plan_forward  = plan1
    grid%plan_backward = plan2
    grid%plan_r2c      = plan3
    grid%plan_c2r      = plan4
    
    grid%allocated =.true.
    return
  end subroutine allocate_gauss_legendre_grid

  
  function ph_gauss_legednre_grid(grid,iph) result(ph)
    implicit none
    class(gauss_legendre_grid), intent(in) :: grid
    integer(i4b), intent(in) :: iph
    real(dp) :: ph
    ph = (iph-1)*twopi/grid%nph
    return
  end function ph_gauss_legednre_grid
  




  !-----------------------------------------------------------------!
  !            spherical harmonic transformation routines           !
  !-----------------------------------------------------------------!

  subroutine SH_trans_gauss_legendre_grid(grid,n,u,ulm)
    use module_fftw3
    implicit none
    class(gauss_legendre_grid), intent(in) :: grid
    integer(i4b), intent(in) :: n
    complex(dpc), dimension(grid%nph*grid%nth), intent(in) :: u
    complex(dpc), dimension((grid%lmax+1)**2), intent(out) :: ulm

    logical :: even
    integer(i4b) :: l,m,nth,nph,ith,ilm,i1,i2,klm,na,lmax
    real(dp) :: fac,w,sign
    complex(dpc) :: tmp
    complex(C_DOUBLE_COMPLEX), pointer :: in(:),out(:)
    type(C_PTR) :: pin,pout
    
    ! initialise the coefficients
    ulm = 0.0_dp
    
    ! get some parameters
    lmax = grid%lmax
    nth = grid%nth
    nph = grid%nph
    
    ! set up the C pointers
    pin  = fftw_alloc_complex(int(nph, C_SIZE_T))
    pout = fftw_alloc_complex(int(nph, C_SIZE_T))
    call c_f_pointer(pin,   in, [nph])
    call c_f_pointer(pout, out, [nph])

    ! do the transformation
    na = abs(n)
    even = modulo(n,2) == 0

    if(even) then

       i2 = 0
       do ith = 1,nth
          i1 = i2+1
          i2 = i1+nph-1
          in = u(i1:i2)
          call fftw_execute_dft(grid%plan_forward,in,out)
          ilm = na*na
          klm = na*(na+1)/2
          w = grid%w(ith)*twopi/grid%nph
          do l = na,lmax-1
             sign = 1.0_dp
             klm = klm+1
             tmp =  grid%dlm(klm,n,ith)*w
             tmp = tmp*out(1)
             ilm = ilm+1
             ulm(ilm) = ulm(ilm) + tmp
             do m = 1,l
                klm = klm+1             
                tmp = grid%dlm(klm,n,ith)*w
                tmp = tmp*out(m+1)
                ilm = ilm+1
                ulm(ilm) = ulm(ilm)+tmp             
                sign = -sign
                tmp = sign*grid%dlm(klm,-n,ith)*w
                tmp = tmp*out(nph-m+1)
                ilm = ilm+1
                ulm(ilm) = ulm(ilm)+tmp
             end do
          end do
          l = lmax
          sign = 1.0_dp
          klm = klm+1
          tmp =  grid%dlm(klm,n,ith)*w
          tmp = tmp*out(1)
          ilm = ilm+1
          ulm(ilm) = ulm(ilm) + tmp
          do m = 1,l-1
             klm = klm+1             
             tmp = grid%dlm(klm,n,ith)*w
             tmp = tmp*out(m+1)
             ilm = ilm+1
             ulm(ilm) = ulm(ilm)+tmp             
             sign = -sign
             tmp = sign*grid%dlm(klm,-n,ith)*w
             tmp = tmp*out(nph-m+1)
             ilm = ilm+1
             ulm(ilm) = ulm(ilm)+tmp             
          end do
          m = l
          klm = klm+1             
          tmp = grid%dlm(klm,n,ith)*w
          tmp = tmp*out(m+1)
          ilm = ilm+1
          ulm(ilm) = ulm(ilm)+tmp
          sign = -sign
       end do

       
    else

       i2 = 0
       do ith = 1,nth
          i1 = i2+1
          i2 = i1+nph-1
          in = u(i1:i2)
          call fftw_execute_dft(grid%plan_forward,in,out)
          ilm = na*na
          klm = na*(na+1)/2
          w = grid%w(ith)*twopi/grid%nph
          do l = na,lmax-1
             sign = -1.0_dp
             klm = klm+1
             tmp =  grid%dlm(klm,n,ith)*w
             tmp = tmp*out(1)
             ilm = ilm+1
             ulm(ilm) = ulm(ilm) + tmp
             do m = 1,l
                klm = klm+1             
                tmp = grid%dlm(klm,n,ith)*w
                tmp = tmp*out(m+1)
                ilm = ilm+1
                ulm(ilm) = ulm(ilm)+tmp             
                sign = -sign
                tmp = sign*grid%dlm(klm,-n,ith)*w
                tmp = tmp*out(nph-m+1)
                ilm = ilm+1
                ulm(ilm) = ulm(ilm)+tmp
             end do
          end do
          l = lmax
          sign = -1.0_dp
          klm = klm+1
          tmp =  grid%dlm(klm,n,ith)*w
          tmp = tmp*out(1)
          ilm = ilm+1
          ulm(ilm) = ulm(ilm) + tmp
          do m = 1,l-1
             klm = klm+1             
             tmp = grid%dlm(klm,n,ith)*w
             tmp = tmp*out(m+1)
             ilm = ilm+1
             ulm(ilm) = ulm(ilm)+tmp             
             sign = -sign
             tmp = sign*grid%dlm(klm,-n,ith)*w
             tmp = tmp*out(nph-m+1)
             ilm = ilm+1
             ulm(ilm) = ulm(ilm)+tmp             
          end do
          m = l
          klm = klm+1             
          tmp = grid%dlm(klm,n,ith)*w
          tmp = tmp*out(m+1)
          ilm = ilm+1
          ulm(ilm) = ulm(ilm)+tmp
          sign = -sign
       end do
       
    end if

    

    ! nullify the pointers
    call fftw_free(pin)
    call fftw_free(pout)
    in  => null()
    out => null()
    
    return
  end subroutine SH_trans_gauss_legendre_grid


  subroutine SH_itrans_gauss_legendre_grid(grid,n,ulm,u)
    use module_fftw3
    implicit none
    class(gauss_legendre_grid), intent(in) :: grid
    integer(i4b), intent(in) :: n
    complex(dpc), dimension((grid%lmax+1)**2), intent(in) :: ulm
    complex(dpc), dimension(grid%nph*grid%nth), intent(out) :: u

    logical :: even
    integer(i4b) :: l,m,nth,nph,ith,ilm,i1,i2,klm,na,lmax,i,j
    real(dp) :: fac,w,sign
    complex(dpc) :: tmp
    complex(dpc), dimension(grid%nph) :: lout
    complex(C_DOUBLE_COMPLEX), pointer :: in(:),out(:)
    type(C_PTR) :: pin,pout
    
    ! initialise the function
    u = 0.0_dp
    
    ! get some parameters
    lmax = grid%lmax
    nth = grid%nth
    nph = grid%nph
    
    ! set up the C pointers
    pin  = fftw_alloc_complex(int(nph, C_SIZE_T))
    pout = fftw_alloc_complex(int(nph, C_SIZE_T))
    call c_f_pointer(pin,   in, [nph])
    call c_f_pointer(pout, out, [nph])

    ! do the transformation
    na = abs(n)
    even = modulo(n,2) == 0


    if(even) then

       i2 = 0
       do ith = 1,nth
          in = 0.0_dp
          ilm = na*na
          klm = na*(na+1)/2
          do l = na,lmax-1
             sign = 1.0_dp
             klm = klm+1
             tmp = grid%dlm(klm,n,ith)
             ilm = ilm+1
             tmp = tmp*ulm(ilm)
             in(1) = in(1) + tmp
             do m = 1,l
                klm = klm+1
                tmp = grid%dlm(klm,n,ith)
                ilm = ilm+1
                tmp = tmp*ulm(ilm)
                in(m+1) = in(m+1) + tmp
                sign = -sign
                tmp = sign*grid%dlm(klm,-n,ith)
                ilm = ilm+1             
                tmp = tmp*ulm(ilm)
                in(nph-m+1) = in(nph-m+1) + tmp
             end do
          end do
          l = lmax
          sign = 1.0_dp
          klm = klm+1
          tmp = grid%dlm(klm,n,ith)
          ilm = ilm+1
          tmp = tmp*ulm(ilm)
          in(1) = in(1) + tmp
          do m = 1,l-1
             klm = klm+1
             tmp = grid%dlm(klm,n,ith)
             ilm = ilm+1
             tmp = tmp*ulm(ilm)
             in(m+1) = in(m+1) + tmp
             sign = -sign
             tmp = sign*grid%dlm(klm,-n,ith)
             ilm = ilm+1             
             tmp = tmp*ulm(ilm)
             in(nph-m+1) = in(nph-m+1) + tmp
          end do
          m = l
          klm = klm+1
          tmp = grid%dlm(klm,n,ith)
          ilm = ilm+1
          tmp = tmp*ulm(ilm)
          in(m+1) = in(m+1) + tmp       
          call fftw_execute_dft(grid%plan_backward,in,out)
          i1 = i2+1
          i2 = i1+nph-1
          u(i1:i2) = out          
       end do

    else

       i2 = 0
       do ith = 1,nth
          in = 0.0_dp
          ilm = na*na
          klm = na*(na+1)/2
          do l = na,lmax-1
             sign = -1.0_dp
             klm = klm+1
             tmp = grid%dlm(klm,n,ith)
             ilm = ilm+1
             tmp = tmp*ulm(ilm)
             in(1) = in(1) + tmp
             do m = 1,l
                klm = klm+1
                tmp = grid%dlm(klm,n,ith)
                ilm = ilm+1
                tmp = tmp*ulm(ilm)
                in(m+1) = in(m+1) + tmp
                sign = -sign
                tmp = sign*grid%dlm(klm,-n,ith)
                ilm = ilm+1             
                tmp = tmp*ulm(ilm)
                in(nph-m+1) = in(nph-m+1) + tmp
             end do
          end do
          l = lmax
          sign = -1.0_dp
          klm = klm+1
          tmp = grid%dlm(klm,n,ith)
          ilm = ilm+1
          tmp = tmp*ulm(ilm)
          in(1) = in(1) + tmp
          do m = 1,l-1
             klm = klm+1
             tmp = grid%dlm(klm,n,ith)
             ilm = ilm+1
             tmp = tmp*ulm(ilm)
             in(m+1) = in(m+1) + tmp
             sign = -sign
             tmp = sign*grid%dlm(klm,-n,ith)
             ilm = ilm+1             
             tmp = tmp*ulm(ilm)
             in(nph-m+1) = in(nph-m+1) + tmp
          end do
          m = l
          klm = klm+1
          tmp = grid%dlm(klm,n,ith)
          ilm = ilm+1
          tmp = tmp*ulm(ilm)
          in(m+1) = in(m+1) + tmp       
          call fftw_execute_dft(grid%plan_backward,in,out)
          i1 = i2+1
          i2 = i1+nph-1
          u(i1:i2) = out
       end do

    end if

    ! nullify pointers
    call fftw_free(pin)
    call fftw_free(pout)
    in  => null()
    out => null()
    
    return
  end subroutine SH_itrans_gauss_legendre_grid



  subroutine real_SH_trans_gauss_legendre_grid(grid,u,ulm)
    use module_fftw3
    implicit none
    class(gauss_legendre_grid), intent(in) :: grid
    real(dp), dimension(grid%nph*grid%nth), intent(in) :: u
    complex(dpc), dimension(((grid%lmax+1)*(grid%lmax+2))/2), intent(out) :: ulm

    integer(i4b) :: l,m,nth,nph,ith,ilm,i1,i2,klm,lmax
    real(dp) :: fac,w
    complex(dpc) :: tmp
    real(C_DOUBLE), pointer :: in(:)
    complex(C_DOUBLE_COMPLEX), pointer :: out(:)
    type(C_PTR) :: pin,pout
    
    ! initialise the coefficients
    ulm = 0.0_dp
    
    ! get some parameters
    lmax = grid%lmax
    nth = grid%nth
    nph = grid%nph
    
    ! set up the C pointers
    pin  = fftw_alloc_real(int(nph, C_SIZE_T))
    pout = fftw_alloc_complex(int(nph/2+1, C_SIZE_T))
    call c_f_pointer(pin,   in, [nph])
    call c_f_pointer(pout, out, [nph/2+1])
    
    i2 = 0
    do ith = 1,nth
       i1 = i2+1
       i2 = i1+nph-1
       in = u(i1:i2)
       call fftw_execute_dft_r2c(grid%plan_r2c,in,out)
       ilm = 0
       klm = 0
       w = grid%w(ith)*twopi/grid%nph
       do l = 0,lmax
          klm = klm+1
          tmp =  grid%dlm(klm,0,ith)*w
          tmp = tmp*out(1)
          ilm = ilm+1
          ulm(ilm) = ulm(ilm) + tmp
          do m = 1,l
             klm = klm+1             
             tmp = grid%dlm(klm,0,ith)*w
             tmp = tmp*out(m+1)
             ilm = ilm+1
             ulm(ilm) = ulm(ilm)+tmp             
          end do
       end do
    end do

    ! nullify the pointers
    call fftw_free(pin)
    call fftw_free(pout)
    in  => null()
    out => null()
    
    return
  end subroutine real_SH_trans_gauss_legendre_grid


  subroutine real_SH_itrans_gauss_legendre_grid(grid,ulm,u)
    use module_fftw3
    implicit none
    class(gauss_legendre_grid), intent(in) :: grid
    complex(dpc), dimension(((grid%lmax+1)*(grid%lmax+2))/2), intent(in) :: ulm
    real(dp), dimension(grid%nph*grid%nth), intent(out) :: u

    integer(i4b) :: l,m,nth,nph,ith,ilm,i1,i2,klm,lmax
    real(dp) :: fac,w
    complex(dpc) :: tmp
    real(C_DOUBLE), pointer :: out(:)
    complex(C_DOUBLE_COMPLEX), pointer :: in(:)
    type(C_PTR) :: pin,pout
    
    ! initialise the function
    u = 0.0_dp
    
    ! get some parameters
    lmax = grid%lmax
    nth = grid%nth
    nph = grid%nph
    
    ! set up the C pointers
    pin = fftw_alloc_complex(int(nph/2+1, C_SIZE_T))
    pout  = fftw_alloc_real(int(nph, C_SIZE_T))
    call c_f_pointer(pin, in, [nph/2+1])
    call c_f_pointer(pout,   out, [nph])

    
    i2 = 0
    do ith = 1,nth
       in = 0.0_dp
       ilm = 0
       klm = 0
       do l = 0,lmax
          klm = klm+1
          tmp = grid%dlm(klm,0,ith)
          ilm = ilm+1
          tmp = tmp*ulm(ilm)
          in(1) = in(1) + tmp
          do m = 1,l
             klm = klm+1
             tmp = grid%dlm(klm,0,ith)
             ilm = ilm+1
             tmp = tmp*ulm(ilm)
             in(m+1) = in(m+1) + tmp
          end do
       end do
       call fftw_execute_dft_c2r(grid%plan_c2r,in,out)
       i1 = i2+1
       i2 = i1+nph-1
       u(i1:i2) = out
    end do
    
    ! nullify the pointers
    call fftw_free(pin)
    call fftw_free(pout)
    in  => null()
    out => null()
    
    return
  end subroutine real_SH_itrans_gauss_legendre_grid


  
  subroutine wrapper_SH_trans_gauss_legendre_grid(grid,u,ulm)
    implicit none
    class(gauss_legendre_grid), intent(in) :: grid
    class(gauss_legendre_field), intent(in) :: u
    class(spherical_harmonic_expansion), intent(inout) :: ulm
    integer(i4b) :: icomp,i1,i2,j1,j2
    do icomp = 1,u%nr
       i1 = u%index(1,1,icomp)
       i2 = u%index(u%nph,u%nth,icomp)
       j1 = ulm%rindex(0,0,icomp)
       j2 = ulm%rindex(grid%lmax,grid%lmax,icomp)
       call real_SH_trans_gauss_legendre_grid(grid,u%rdata(i1:i2),ulm%rdata(j1:j2))       
    end do
    do icomp = 1,u%nc
       i1 = u%index(1,1,icomp)
       i2 = u%index(u%nph,u%nth,icomp)
       j1 = ulm%cindex(0,0,icomp)
       j2 = ulm%cindex(-grid%lmax,grid%lmax,icomp)
       call SH_trans_gauss_legendre_grid(grid,u%cnval(icomp),u%cdata(i1:i2),ulm%cdata(j1:j2))       
    end do
    return
  end subroutine wrapper_SH_trans_gauss_legendre_grid

  
  subroutine wrapper_SH_itrans_gauss_legendre_grid(grid,ulm,u)
    implicit none
    class(gauss_legendre_grid), intent(in) :: grid
    class(spherical_harmonic_expansion), intent(in) :: ulm
    class(gauss_legendre_field), intent(inout) :: u
    integer(i4b) :: icomp,i1,i2,j1,j2    
    do icomp = 1,u%nr
       i1 = u%index(1,1,icomp)
       i2 = u%index(u%nph,u%nth,icomp)
       j1 = ulm%rindex(0,0,icomp)
       j2 = ulm%rindex(grid%lmax,grid%lmax,icomp)
       call real_SH_itrans_gauss_legendre_grid(grid,ulm%rdata(j1:j2),u%rdata(i1:i2))       
    end do
    do icomp = 1,u%nc
       i1 = u%index(1,1,icomp)
       i2 = u%index(u%nph,u%nth,icomp)
       j1 = ulm%cindex(0,0,icomp)
       j2 = ulm%cindex(-ulm%lmax,ulm%lmax,icomp)
       call SH_itrans_gauss_legendre_grid(grid,u%cnval(icomp),ulm%cdata(j1:j2),u%cdata(i1:i2))
    end do
    return
  end subroutine wrapper_SH_itrans_gauss_legendre_grid

    
  !=======================================================================!
  !                    procedures for the field types                     !
  !=======================================================================!

  !-----------------------------------------------!
  !                 the basic type                !
  !-----------------------------------------------!

  subroutine delete_gauss_legendre_field(self)
    implicit none
    class(gauss_legendre_field), intent(inout) :: self
    if(.not.self%allocated) return
    deallocate(self%cdata)
    deallocate(self%cnval)
    deallocate(self%rdata)
    self%allocated = .false.
    return
  end subroutine delete_gauss_legendre_field

  
  subroutine allocate_gauss_legendre_field(self,lmax,nc,cnval,nr)
    implicit none
    class(gauss_legendre_field), intent(inout) :: self
    integer(i4b), intent(in) :: lmax
    integer(i4b), intent(in) :: nc
    integer(i4b), intent(in), dimension(nc) :: cnval
    integer(i4b), intent(in) :: nr
    call self%delete()    
    self%lmax = lmax
    self%nc = nc
    self%nth = lmax+1
    self%nph = 2*lmax
    self%cdim = self%nc*self%nth*self%nph
    allocate(self%cnval(self%nc))
    allocate(self%cdata(self%cdim))
    if(nc > 0) then
       self%cdata = 0.0_dp
       self%cnval = cnval
    end if
    self%nr = nr
    self%rdim = self%nr*self%nth*self%nph
    allocate(self%rdata(self%rdim))
    if(nr > 0) self%rdata = 0.0_dp
    self%allocated = .true.    
    return
  end subroutine allocate_gauss_legendre_field
  
  
  function index_gauss_legendre_field(self,iph,ith,icomp) result(i)
    implicit none
    class(gauss_legendre_field), intent(in) :: self
    integer(i4b), intent(in) :: iph,ith,icomp
    integer(i4b) :: i
    i =  self%nph*self%nth*(icomp-1)+self%nph*(ith-1)+iph
    return
  end function index_gauss_legendre_field

  subroutine scale_gauss_legendre_field(self,a)
    implicit none
    class(gauss_legendre_field), intent(inout) :: self
    complex(dpc), intent(in) :: a
    call zscal(self%cdim,a,self%cdata,1)
    return
  end subroutine scale_gauss_legendre_field

  subroutine real_scale_gauss_legendre_field(self,a)
    implicit none
    class(gauss_legendre_field), intent(inout) :: self
    real(dp), intent(in) :: a
    call zdscal(self%cdim,a,self%cdata,1)
    return
  end subroutine real_scale_gauss_legendre_field

  subroutine saxpy_gauss_legendre_field(self,a,other)
    implicit none
    class(gauss_legendre_field), intent(inout) :: self
    complex(dpc), intent(in) :: a
    class(gauss_legendre_field), intent(in) :: other
    call zaxpy(self%cdim,a,other%cdata,1,self%cdata,1)
    return
  end subroutine saxpy_gauss_legendre_field

  subroutine real_saxpy_gauss_legendre_field(self,a,other)
    implicit none
    class(gauss_legendre_field), intent(inout) :: self
    real(dp), intent(in) :: a
    class(gauss_legendre_field), intent(in) :: other
    complex(dpc) :: aloc
    aloc = a
    call zaxpy(self%cdim,aloc,other%cdata,1,self%cdata,1)
    return
  end subroutine real_saxpy_gauss_legendre_field


  !------------------------------------------------------!
  !            procedures for derived types              !
  !------------------------------------------------------!

  ! scalar fields
  
  subroutine allocate_scalar_gauss_legendre_field(self,grid,real)
    implicit none
    class(scalar_gauss_legendre_field), intent(inout) :: self
    type(gauss_legendre_grid), intent(in) :: grid
    logical, intent(in), optional :: real
    integer(i4b), dimension(0) :: null
    if(present(real)) then
       self%real = real
    end if
    if(self%real) then
       call allocate_gauss_legendre_field(self,grid%lmax,0,null,1)
    else
       call allocate_gauss_legendre_field(self,grid%lmax,1,(/0/),0)
    end if
    return
  end subroutine allocate_scalar_gauss_legendre_field

  subroutine set_complex_scalar_gauss_legendre_field(self,iph,ith,val)
    implicit none
    class(scalar_gauss_legendre_field), intent(inout) :: self
    integer(i4b), intent(in) :: ith,iph
    complex(dpc), intent(in) :: val
    self%cdata(self%index(iph,ith,1)) = val    
    return
  end subroutine set_complex_scalar_gauss_legendre_field


  subroutine set_real_scalar_gauss_legendre_field(self,iph,ith,val)
    implicit none
    class(scalar_gauss_legendre_field), intent(inout) :: self
    integer(i4b), intent(in) :: ith,iph
    real(dp), intent(in) :: val
    self%rdata(self%index(iph,ith,1)) = val    
    return
  end subroutine set_real_scalar_gauss_legendre_field

  function get_complex_scalar_gauss_legendre_field(self,iph,ith) result(val)
    implicit none
    class(scalar_gauss_legendre_field), intent(in) :: self
    integer(i4b), intent(in) :: ith,iph
    complex(dpc) :: val
    val = self%cdata(self%index(iph,ith,1)) 
    return
  end function get_complex_scalar_gauss_legendre_field


  function get_real_scalar_gauss_legendre_field(self,iph,ith) result(val)
    implicit none
    class(scalar_gauss_legendre_field), intent(in) :: self
    integer(i4b), intent(in) :: ith,iph
    real(dp) :: val
    val = self%rdata(self%index(iph,ith,1)) 
    return
  end function get_real_scalar_gauss_legendre_field

  ! vector fields
  
  subroutine allocate_vector_gauss_legendre_field(self,grid,real)
    implicit none
    class(vector_gauss_legendre_field), intent(inout) :: self
    type(gauss_legendre_grid), intent(in) :: grid
    logical, intent(in), optional :: real
    if(present(real)) then
       self%real = real
    end if
    if(self%real) then
       call allocate_gauss_legendre_field(self,grid%lmax,1,(/-1/),1)
    else
       call allocate_gauss_legendre_field(self,grid%lmax,3,(/-1,0,1/),0)
    end if
    return
  end subroutine allocate_vector_gauss_legendre_field

  subroutine set_value_vector_gauss_legendre_field(self,iph,ith,alpha,val)
    implicit none
    class(vector_gauss_legendre_field), intent(inout) :: self
    integer(i4b), intent(in) :: ith,iph,alpha
    complex(dpc), intent(in) :: val
    integer(i4b) :: icomp
    if(self%real) then
       select case(alpha)
       case(-1)
          self%cdata(self%index(iph,ith,1)) = val    
       case(0)
          self%rdata(self%index(iph,ith,1)) = real(val)
       case(1)
          self%cdata(self%index(iph,ith,1)) = -conjg(val)
       end select
    else
       icomp = alpha+2
       self%cdata(self%index(iph,ith,icomp)) = val
    end if
    return
  end subroutine set_value_vector_gauss_legendre_field

  function get_value_vector_gauss_legendre_field(self,iph,ith,alpha) result(val)
    implicit none
    class(vector_gauss_legendre_field), intent(in) :: self
    integer(i4b), intent(in) :: ith,iph,alpha
    complex(dpc) :: val
    integer(i4b) :: icomp
    if(self%real) then
       select case(alpha)
       case(-1)
          val = self%cdata(self%index(iph,ith,1)) 
       case(0)
          val = self%rdata(self%index(iph,ith,1))
       case(1)
          val = -conjg(self%cdata(self%index(iph,ith,1)))
       end select
    else
       icomp = alpha+2
       val = self%cdata(self%index(iph,ith,icomp))
    end if
    return
  end function get_value_vector_gauss_legendre_field


  ! second order tensor fields
  
  subroutine allocate_second_tensor_gauss_legendre_field(self,grid,symmetric,trace_free,real)
    implicit none
    class(second_tensor_gauss_legendre_field), intent(inout) :: self
    type(gauss_legendre_grid), intent(in) :: grid
    logical, intent(in), optional :: symmetric,trace_free,real
    if(present(symmetric)) then
       self%symmetric = symmetric
    end if
    if(present(trace_free)) then
       self%trace_free = trace_free
    end if
    if(present(real)) then
       self%real = real
    end if
    if(.not.(self%symmetric .or. self%trace_free .or. self%real)) then       
       call allocate_gauss_legendre_field(self,grid%lmax,9,(/-2,-1,0,-1,0,1,0,1,2/),0)       
    else
       stop 'allocate_second_tensor_gauss_legendre_field: invalid allocation choice'
    end if
    
    return
  end subroutine allocate_second_tensor_gauss_legendre_field

  subroutine set_value_second_tensor_gauss_legendre_field(self,iph,ith,alpha,beta,val)
    implicit none
    class(second_tensor_gauss_legendre_field), intent(inout) :: self
    integer(i4b), intent(in) :: ith,iph,alpha,beta
    complex(dpc), intent(in) :: val
    integer(i4b) :: icomp
    icomp = 3*(beta+1) + alpha + 2
    self%cdata(self%index(iph,ith,icomp)) = val    
    return
  end subroutine set_value_second_tensor_gauss_legendre_field

  function get_value_second_tensor_gauss_legendre_field(self,iph,ith,alpha,beta) result(val)
    implicit none
    class(second_tensor_gauss_legendre_field), intent(in) :: self
    integer(i4b), intent(in) :: ith,iph,alpha,beta
    complex(dpc) :: val
    integer(i4b) :: icomp
    icomp = 3*(beta+1) + alpha + 2
    val = self%cdata(self%index(iph,ith,icomp)) 
    return
  end function get_value_second_tensor_gauss_legendre_field




  
  
  
  !=======================================================================!
  !          procedures for the spherical_harmonic_expansion type         !
  !=======================================================================!


  !-----------------------------------------------!
  !                 the basic type                !
  !-----------------------------------------------!
  
  subroutine delete_spherical_harmonic_expansion(self)
    implicit none
    class(spherical_harmonic_expansion), intent(inout) :: self
    if(.not.self%allocated) return
    deallocate(self%cdata)
    deallocate(self%cnval)
    deallocate(self%rdata)
    self%allocated = .false.
    return
  end subroutine delete_spherical_harmonic_expansion


  subroutine allocate_spherical_harmonic_expansion(self,lmax,nc,cnval,nr)
    implicit none
    class(spherical_harmonic_expansion), intent(inout) :: self
    integer(i4b), intent(in) :: lmax
    integer(i4b), intent(in) :: nc
    integer(i4b), dimension(nc), intent(in) :: cnval
    integer(i4b), intent(in) :: nr
    call self%delete()
    self%lmax= lmax
    self%nc = nc
    self%cncoef = (lmax+1)**2
    self%cdim = nc*self%cncoef
    allocate(self%cnval(nc))
    allocate(self%cdata(self%cdim))
    if(nc > 0) then
       self%cdata = 0.0_dp
       self%cnval = cnval
    end if
    self%nr = nr
    self%rncoef = ((lmax+1)*(lmax+2))/2
    self%rdim = nr*self%rncoef
    allocate(self%rdata(self%rdim))
    if(nr > 0) self%rdata = 0.0_dp
    self%allocated = .true.
    return
  end subroutine allocate_spherical_harmonic_expansion

  function cindex_spherical_harmonic_expansion(self,m,l,icomp) result(i)
    implicit none
    class(spherical_harmonic_expansion), intent(in) :: self
    integer(i4b), intent(in) :: m,l,icomp
    integer(i4b) :: i
    i = self%cncoef*(icomp-1) + l**2
    if(m == 0) then
       i = i+1
    else if(m > 0) then
       i = i + 2*m
    else
       i = i -2*m+1
    end if
    return
  end function cindex_spherical_harmonic_expansion


  function rindex_spherical_harmonic_expansion(self,m,l,icomp) result(i)
    implicit none
    class(spherical_harmonic_expansion), intent(in) :: self
    integer(i4b), intent(in) :: m,l,icomp
    integer(i4b) :: i
    i = self%rncoef*(icomp-1) + l*(l+1)/2 + m + 1
    return
  end function rindex_spherical_harmonic_expansion
  
  
  subroutine set_sobolev_spherical_harmonic_expansion(u,s,mu)
    implicit none
    class(spherical_harmonic_expansion), intent(inout) :: u
    real(dp), intent(in) :: s,mu
    u%s  = s
    u%mu = mu
    return
  end subroutine set_sobolev_spherical_harmonic_expansion

  subroutine scale_spherical_harmonic_expansion(self,a)
    implicit none
    class(spherical_harmonic_expansion), intent(inout) :: self
    complex(dpc), intent(in) :: a
    call zscal(self%cdim,a,self%cdata,1)
    return
  end subroutine scale_spherical_harmonic_expansion

  subroutine real_scale_spherical_harmonic_expansion(self,a)
    implicit none
    class(spherical_harmonic_expansion), intent(inout) :: self
    real(dp), intent(in) :: a
    call zdscal(self%cdim,a,self%cdata,1)
    return
  end subroutine real_scale_spherical_harmonic_expansion

  subroutine saxpy_spherical_harmonic_expansion(self,a,other)
    implicit none
    class(spherical_harmonic_expansion), intent(inout) :: self
    complex(dpc), intent(in) :: a
    class(spherical_harmonic_expansion), intent(in) :: other
    call zaxpy(self%cdim,a,other%cdata,1,self%cdata,1)
    return
  end subroutine saxpy_spherical_harmonic_expansion

  subroutine real_saxpy_spherical_harmonic_expansion(self,a,other)
    implicit none
    class(spherical_harmonic_expansion), intent(inout) :: self
    real(dp), intent(in) :: a
    class(spherical_harmonic_expansion), intent(in) :: other
    complex(dpc) :: aloc
    aloc = a
    call zaxpy(self%cdim,aloc,other%cdata,1,self%cdata,1)
    return
  end subroutine real_saxpy_spherical_harmonic_expansion


  !------------------------------------------------------!
  !            procedures for derived types              !
  !------------------------------------------------------!


  ! scalar fields
  
  subroutine allocate_scalar_spherical_harmonic_expansion(self,grid,real)
    implicit none
    class(scalar_spherical_harmonic_expansion), intent(inout) :: self
    type(gauss_legendre_grid), intent(in) :: grid
    logical, intent(in), optional :: real
    integer(i4b), dimension(0) :: null
    if(present(real)) then
       self%real = real
    end if
    if(self%real) then
       call allocate_spherical_harmonic_expansion(self,grid%lmax,0,null,1)
    else
       call allocate_spherical_harmonic_expansion(self,grid%lmax,1,(/0/),0)
    end if
    return
  end subroutine allocate_scalar_spherical_harmonic_expansion

  subroutine set_complex_scalar_spherical_harmonic_expansion(self,l,m,val)
    implicit none
    class(scalar_spherical_harmonic_expansion), intent(inout) :: self
    integer(i4b), intent(in) :: l,m
    complex(dpc), intent(in) :: val
    if(self%real) then
       self%rdata(self%rindex(m,l,1)) = val
    else       
       self%cdata(self%cindex(m,l,1)) = val
    end if
    return
  end subroutine set_complex_scalar_spherical_harmonic_expansion

  function get_complex_scalar_spherical_harmonic_expansion(self,l,m) result(val)
    implicit none
    class(scalar_spherical_harmonic_expansion), intent(in) :: self
    integer(i4b), intent(in) :: l,m
    complex(dpc) :: val
    if(self%real) then
       val = self%rdata(self%rindex(m,l,1))
    else
       val = self%cdata(self%cindex(m,l,1))
    end if
    return
  end function get_complex_scalar_spherical_harmonic_expansion


  ! vector fields
  
  subroutine allocate_vector_spherical_harmonic_expansion(self,grid,real)
    implicit none
    class(vector_spherical_harmonic_expansion), intent(inout) :: self
    type(gauss_legendre_grid), intent(in) :: grid
    logical, intent(in), optional :: real
    if(present(real)) then
       self%real = real
    end if
    if(self%real) then
       call allocate_spherical_harmonic_expansion(self,grid%lmax,1,(/-1/),1)
    else
       call allocate_spherical_harmonic_expansion(self,grid%lmax,3,(/-1,0,1/),0)
    end if
    return
  end subroutine allocate_vector_spherical_harmonic_expansion

  subroutine set_value_vector_spherical_harmonic_expansion(self,l,m,alpha,val)
    implicit none
    class(vector_spherical_harmonic_expansion), intent(inout) :: self
    integer(i4b), intent(in) :: l,m,alpha
    complex(dpc), intent(in) :: val
    integer(i4b) :: icomp,sign
    if(self%real) then
       select case(alpha)
       case(-1)
          self%cdata(self%cindex(m,l,1)) = val    
       case(0)
          self%rdata(self%rindex(m,l,1)) = val
       case(1)
          if(modulo(m,2) == 0) then
             sign = 1
          else
             sign = -1
          end if
          self%cdata(self%cindex(-m,l,1)) = sign*conjg(val)
       end select
    else
       icomp = alpha+2
       self%cdata(self%cindex(m,l,icomp)) = val
    end if
    return
  end subroutine set_value_vector_spherical_harmonic_expansion

  function get_value_vector_spherical_harmonic_expansion(self,l,m,alpha) result(val)
    implicit none
    class(vector_spherical_harmonic_expansion), intent(in) :: self
    integer(i4b), intent(in) :: l,m,alpha
    complex(dpc) :: val
    integer(i4b) :: icomp,sign
    if(self%real) then
       select case(alpha)
       case(-1)
          val = self%cdata(self%cindex(m,l,1))
       case(0)
          val = self%rdata(self%rindex(m,l,1))
       case(1)
          if(modulo(m,2) == 0) then
             sign = 1
          else
             sign = -1
          end if
          val = sign*conjg(self%cdata(self%cindex(-m,l,1)))
       end select
    else
       icomp = alpha+2
       val = self%cdata(self%cindex(m,l,icomp))
    end if
    return
  end function get_value_vector_spherical_harmonic_expansion
  
  ! second order tensor fields
  
  subroutine allocate_second_tensor_spherical_harmonic_expansion(self,grid)
    implicit none
    class(second_tensor_spherical_harmonic_expansion), intent(inout) :: self
    type(gauss_legendre_grid), intent(in) :: grid
    call allocate_spherical_harmonic_expansion(self,grid%lmax,9,(/-2,-1,0,-1,0,1,0,1,2/),0)
    return
  end subroutine allocate_second_tensor_spherical_harmonic_expansion

  subroutine set_value_second_tensor_spherical_harmonic_expansion(self,l,m,alpha,beta,val)
    implicit none
    class(second_tensor_spherical_harmonic_expansion), intent(inout) :: self
    integer(i4b), intent(in) :: l,m,alpha,beta
    complex(dpc), intent(in) :: val
    integer(i4b) :: icomp
    icomp = 3*(beta+1) + alpha + 2
    self%cdata(self%cindex(m,l,icomp)) = val    
    return
  end subroutine set_value_second_tensor_spherical_harmonic_expansion

  function get_value_second_tensor_spherical_harmonic_expansion(self,l,m,alpha,beta) result(val)
    implicit none
    class(second_tensor_spherical_harmonic_expansion), intent(in) :: self
    integer(i4b), intent(in) :: l,m,alpha,beta
    complex(dpc) :: val
    integer(i4b) :: icomp
    icomp = 3*(beta+1) + alpha + 2
    val = self%cdata(self%cindex(m,l,icomp)) 
    return
  end function get_value_second_tensor_spherical_harmonic_expansion



  
  
end module module_spherical_harmonics




