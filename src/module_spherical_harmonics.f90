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
     type(wigner_array), dimension(:), allocatable :: d
     type(C_PTR) :: plan_forward
     type(C_PTR) :: plan_backward     
   contains
     procedure :: delete  => delete_gauss_legendre_grid
     procedure :: allocate   =>  allocate_gauss_legendre_grid
     procedure :: ph  =>  ph_gauss_legednre_grid
     procedure :: check_field => check_field_gauss_legendre_grid
     procedure :: SH_trans_gauss_legendre_grid
     procedure :: SH_trans_scalar_gauss_legendre_grid
     procedure :: SH_trans_vector_gauss_legendre_grid
     generic :: SH_trans => SH_trans_gauss_legendre_grid,        &
                            SH_trans_scalar_gauss_legendre_grid, &
                            SH_trans_vector_gauss_legendre_grid
     procedure :: real_SH_trans_gauss_legendre_grid
     procedure :: real_SH_trans_scalar_gauss_legendre_grid
     generic :: real_SH_trans => real_SH_trans_gauss_legendre_grid, & 
                                 real_SH_trans_scalar_gauss_legendre_grid
  end type gauss_legendre_grid


  !===============================================================!
  !              type declarations for the GL-fields              !
  !===============================================================!

  type, abstract :: gauss_legendre_field
     logical :: allocated = .false.
     integer(i4b) :: lmax
     integer(i4b) :: ndim
     integer(i4b) :: nth
     integer(i4b) :: nph
     complex(dpc), dimension(:), allocatable :: data
   contains
     procedure :: delete => delete_gauss_legendre_field
     procedure :: check => check_gauss_legendre_field
     procedure :: assign_gauss_legendre_field
     generic, public :: assignment(=) => assign_gauss_legendre_field
     procedure, pass(self) :: add_gauss_legendre_field     
     generic, public :: operator(+)  => add_gauss_legendre_field     
     procedure, pass(self) :: subtract_gauss_legendre_field     
     generic, public :: operator(-)  => subtract_gauss_legendre_field     
     procedure, pass(self) :: left_multiply_gauss_legendre_field
     procedure, pass(self) :: right_multiply_gauss_legendre_field
     procedure, pass(self) :: real_left_multiply_gauss_legendre_field
     procedure, pass(self) :: real_right_multiply_gauss_legendre_field
     generic, public :: operator(*)  => left_multiply_gauss_legendre_field,       &
                                        right_multiply_gauss_legendre_field,      &
                                        real_left_multiply_gauss_legendre_field,  &
                                        real_right_multiply_gauss_legendre_field
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
   contains
     procedure :: allocate => allocate_scalar_gauss_legendre_field
     procedure :: index => index_scalar_gauss_legendre_field
  end type scalar_gauss_legendre_field

  
  type, extends(gauss_legendre_field) :: vector_gauss_legendre_field
   contains
     procedure :: allocate => allocate_vector_gauss_legendre_field
     procedure :: index => index_vector_gauss_legendre_field
  end type vector_gauss_legendre_field

  
  type, extends(gauss_legendre_field) :: real_vector_gauss_legendre_field
   contains
     procedure :: allocate => allocate_real_vector_gauss_legendre_field
     procedure :: index => index_real_vector_gauss_legendre_field
  end type real_vector_gauss_legendre_field


  !===============================================================!
  !                 type declarations SH expansions               !
  !===============================================================!
  
  type, abstract :: spherical_harmonic_expansion
     logical :: allocated
     integer(i4b) :: lmax
     integer(i4b) :: nmax
     integer(i4b) :: ndim
     complex(dpc), dimension(:), allocatable :: data
     real(dp) :: s = 0.0_dp
     real(dp) :: mu = 1.0_dp
   contains
     procedure :: delete    => delete_spherical_harmonic_expansion
     procedure :: set_sobolev => set_sobolev_spherical_harmonic_expansion
     procedure :: check => check_spherical_harmonic_expansion 
     procedure :: assign_spherical_harmonic_expansion     
     generic, public :: assignment(=) => assign_spherical_harmonic_expansion     
     procedure, pass(self) :: add_spherical_harmonic_expansion     
     generic, public :: operator(+)  => add_spherical_harmonic_expansion     
     procedure, pass(self) :: subtract_spherical_harmonic_expansion     
     generic, public :: operator(-)  => subtract_spherical_harmonic_expansion     
     procedure, pass(self) :: left_multiply_spherical_harmonic_expansion
     procedure, pass(self) :: right_multiply_spherical_harmonic_expansion
     procedure, pass(self) :: real_left_multiply_spherical_harmonic_expansion
     procedure, pass(self) :: real_right_multiply_spherical_harmonic_expansion
     generic, public :: operator(*)  => left_multiply_spherical_harmonic_expansion,       &
                                        right_multiply_spherical_harmonic_expansion,      &
                                        real_left_multiply_spherical_harmonic_expansion,  &
                                        real_right_multiply_spherical_harmonic_expansion
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
   contains
     procedure :: allocate => allocate_scalar_spherical_harmonic_expansion
     procedure :: index    => index_scalar_spherical_harmonic_expansion
  end type scalar_spherical_harmonic_expansion


  type, extends(spherical_harmonic_expansion) :: real_scalar_spherical_harmonic_expansion
   contains
     procedure :: allocate => allocate_real_scalar_spherical_harmonic_expansion
     procedure :: index    => index_real_scalar_spherical_harmonic_expansion
  end type real_scalar_spherical_harmonic_expansion


  type, extends(spherical_harmonic_expansion) :: vector_spherical_harmonic_expansion
   contains
     procedure :: allocate => allocate_vector_spherical_harmonic_expansion
     procedure :: index    => index_vector_spherical_harmonic_expansion
  end type vector_spherical_harmonic_expansion


  type, extends(spherical_harmonic_expansion) :: real_vector_spherical_harmonic_expansion
   contains
     procedure :: allocate => allocate_real_vector_spherical_harmonic_expansion
     procedure :: index    => index_real_vector_spherical_harmonic_expansion
  end type real_vector_spherical_harmonic_expansion
  
     
contains

  !=======================================================================!
  !                     procedures for the grid type                      !
  !=======================================================================!


  !--------------------------------------!
  !             basic routines           !
  !--------------------------------------!
  
  subroutine delete_gauss_legendre_grid(grid)
    implicit none
    class(gauss_legendre_grid), intent(inout) :: grid
    if(.not.grid%allocated) return
    deallocate(grid%th,grid%w)
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

    
    integer(i4b) :: l,ith,n,nth,nph,ndim
    real(dp) :: fac
    class(orthogonal_polynomial), allocatable :: poly
    type(gauss_quadrature) :: quad
    complex(C_DOUBLE_COMPLEX), pointer ::  in(:)
    complex(C_DOUBLE_COMPLEX), pointer :: out(:)
    type(C_PTR) :: pin,pout,plan1,plan2
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
    allocate(grid%d(nth))
    do ith = 1,nth
       call grid%d(ith)%set(grid%th(ith),lmax,nmax)
    end do
    
    ! make the FFTW3 plans
    n = 2*lmax
    pin  = fftw_alloc_complex(int(n, C_SIZE_T))
    pout = fftw_alloc_complex(int(n, C_SIZE_T))
    call c_f_pointer(pin,   in, [n])
    call c_f_pointer(pout, out, [n])
    plan1 = fftw_plan_dft_1d(n, in, out, FFTW_FORWARD, plan_flag)
    plan2 = fftw_plan_dft_1d(n, out, in, FFTW_BACKWARD, plan_flag)
    grid%plan_forward  = plan1
    grid%plan_backward = plan2
    call fftw_free(pin)
    call fftw_free(pout)
    
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
  

  function check_field_gauss_legendre_grid(grid,u) result(check)
    implicit none
    class(gauss_legendre_grid), intent(in) :: grid
    type(scalar_gauss_legendre_field), intent(in) :: u
    logical :: check
    check = grid%allocated .and. u%allocated
    check = check .and. (grid%lmax == u%lmax)    
    return
  end function check_field_gauss_legendre_grid


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
    integer(i4b) :: l,m,nth,nph,ith,ilm,im,jlm,i1,i2,klm,na,lmax
    real(dp) :: fac,w,sign
    complex(dpc) :: tmp1,tmp2
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
    
    if(n == 0) then
       
       ! perform the transformation
       i2 = 0
       do ith = 1,nth

          i1 = i2+1
          i2 = i1+nph-1
          in = u(i1:i2)
          call fftw_execute_dft(grid%plan_forward,in,out)
          sign = 1.0_dp
          ilm = 0
          klm = 0
          w = grid%w(ith)*twopi/grid%nph
          do l = 0,lmax
             klm = klm+1
             tmp1 =  grid%d(ith)%data(klm,0)*w
             tmp1 = tmp1*out(1)
             ilm = ilm+1
             ulm(ilm) = ulm(ilm) + tmp1
             do m = 1,l
                klm = klm+1             
                tmp1 = grid%d(ith)%data(klm,0)*w
                tmp2 = sign*tmp1
                tmp1 = tmp1*out(m+1)
                ilm = ilm+1
                ulm(ilm) = ulm(ilm)+tmp1
                tmp2 = tmp2*out(nph-m+1)
                ilm = ilm+1
                ulm(ilm) = ulm(ilm)+tmp2
                sign = -sign
             end do
          end do
       end do

    else

       na = abs(n)
       even = modulo(n,2) == 0       
       i2 = 0
       do ith = 1,nth
          i1 = i2+1
          i2 = i1+nph-1
          in = u(i1:i2)
          call fftw_execute_dft(grid%plan_forward,in,out)
          if(even) then
             sign = 1.0_dp
          else
             sign = -1.0_dp
          end if
          ilm = na*na
          klm = na*(na+1)/2
          w = grid%w(ith)*twopi/grid%nph
          do l = na,lmax
             klm = klm+1
             tmp1 =  grid%d(ith)%data(klm,n)*w
             tmp1 = tmp1*out(1)
             ilm = ilm+1
             ulm(ilm) = ulm(ilm) + tmp1
             do m = 1,l
                klm = klm+1             
                tmp1 = grid%d(ith)%data(klm,n)*w
                tmp1 = tmp1*out(m+1)
                ilm = ilm+1
                ulm(ilm) = ulm(ilm)+tmp1
                tmp2 = sign*grid%d(ith)%data(klm,-n)*w
                tmp2 = tmp2*out(nph-m+1)
                ilm = ilm+1
                ulm(ilm) = ulm(ilm)+tmp2
                sign = -sign
             end do
          end do
       end do
       
    end if
       
    ! add in normalising factors
    ilm = 0
    do l = 0,lmax
       fac = sqrt((2*l+1)/fourpi)
       do m = -l,l
          ilm = ilm+1
          ulm(ilm) = fac*ulm(ilm)
       end do
    end do
       
    return
  end subroutine SH_trans_gauss_legendre_grid


  subroutine real_SH_trans_gauss_legendre_grid(grid,n,u,ulm)
    use module_fftw3
    implicit none
    class(gauss_legendre_grid), intent(in) :: grid
    integer(i4b), intent(in) :: n
    complex(dpc), dimension(grid%nph*grid%nth), intent(in) :: u
    complex(dpc), dimension(((grid%lmax+1)*(grid%lmax+2))/2), intent(out) :: ulm

    logical :: even
    integer(i4b) :: l,m,nth,nph,ith,ilm,im,jlm,i1,i2,klm,na,lmax
    real(dp) :: fac,w,sign
    complex(dpc) :: tmp1,tmp2
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

    na = abs(n)
    i2 = 0
    do ith = 1,nth
       i1 = i2+1
       i2 = i1+nph-1
       in = u(i1:i2)
       call fftw_execute_dft(grid%plan_forward,in,out)
       ilm = na*(na+1)/2
       klm = na*(na+1)/2
       w = grid%w(ith)*twopi/grid%nph
       do l = na,lmax
          do m = 0,l
             klm = klm+1             
             tmp1 = grid%d(ith)%data(klm,n)*w
             tmp1 = tmp1*out(m+1)
             ilm = ilm+1
             ulm(ilm) = ulm(ilm)+tmp1
          end do
       end do
    end do
    
    ! add in normalising factors
    ilm = 0
    do l = 0,lmax
       fac = sqrt((2*l+1)/fourpi)
       do m = -l,l
          ilm = ilm+1
          ulm(ilm) = fac*ulm(ilm)
       end do
    end do
       
    return
  end subroutine real_SH_trans_gauss_legendre_grid


  !-----------------------------------------------------------------!
  !            spherical harmonic transformation wrappers           !
  !-----------------------------------------------------------------!
  
  subroutine SH_trans_scalar_gauss_legendre_grid(grid,u,ulm)
    implicit none
    class(gauss_legendre_grid), intent(in) :: grid
    type(scalar_gauss_legendre_field), intent(in) :: u
    type(scalar_spherical_harmonic_expansion), intent(inout) :: ulm
    call grid%SH_trans(0,u%data,ulm%data)
    return
  end subroutine SH_trans_scalar_gauss_legendre_grid

  subroutine real_SH_trans_scalar_gauss_legendre_grid(grid,u,ulm)
    implicit none
    class(gauss_legendre_grid), intent(in) :: grid
    type(scalar_gauss_legendre_field), intent(in) :: u
    type(real_scalar_spherical_harmonic_expansion), intent(inout) :: ulm
    call grid%real_SH_trans(0,u%data,ulm%data)
    return
  end subroutine real_SH_trans_scalar_gauss_legendre_grid
  
  subroutine SH_trans_vector_gauss_legendre_grid(grid,u,ulm)
    implicit none
    class(gauss_legendre_grid), intent(in) :: grid
    type(vector_gauss_legendre_field), intent(in) :: u
    type(vector_spherical_harmonic_expansion), intent(inout) :: ulm
    integer(i4b) :: i1,i2,j1,j2,alpha
    do alpha = -1,1
       i1 = u%index(alpha,1,1)
       i2 = u%index(alpha,grid%nph,grid%nth)
       j1 = ulm%index(alpha,0,0)
       j2 = ulm%index(alpha,grid%lmax,-grid%lmax)
       call grid%SH_trans(alpha,u%data(i1:i2),ulm%data(j1:j2))
    end do
    return
  end subroutine SH_trans_vector_gauss_legendre_grid

  
  subroutine SH_trans_real_vector_gauss_legendre_grid(grid,u,ulm)
    implicit none
    class(gauss_legendre_grid), intent(in) :: grid
    type(real_vector_gauss_legendre_field), intent(in) :: u
    type(real_vector_spherical_harmonic_expansion), intent(inout) :: ulm
    integer(i4b) :: i1,i2,j1,j2,alpha
    do alpha = -1,0
       i1 = u%index(alpha,1,1)
       i2 = u%index(alpha,grid%nph,grid%nth)
       j1 = ulm%index(alpha,0,0)
       j2 = ulm%index(alpha,grid%lmax,-grid%lmax)
       call grid%real_SH_trans(alpha,u%data(i1:i2),ulm%data(j1:j2))
    end do
    return
  end subroutine SH_trans_real_vector_gauss_legendre_grid





  !=======================================================================!
  !                    procedures for the field types                     !
  !=======================================================================!

  !-----------------------------------------------!
  !               the abstract type               !
  !-----------------------------------------------!

  subroutine delete_gauss_legendre_field(self)
    implicit none
    class(gauss_legendre_field), intent(inout) :: self
    if(.not.self%allocated) return
    deallocate(self%data)
    self%allocated = .false.
    return
  end subroutine delete_gauss_legendre_field
  
  function nth_gauss_legendre_field(self) result(nth)
    implicit none
    class(gauss_legendre_field), intent(in) :: self
    integer(i4b) :: nth
    nth = self%lmax+1
    return
  end function nth_gauss_legendre_field

  function nph_gauss_legendre_field(self) result(nph)
    implicit none
    class(gauss_legendre_field), intent(in) :: self
    integer(i4b) :: nph
    nph = 2*self%lmax
    return
  end function nph_gauss_legendre_field
  
  function check_gauss_legendre_field(self,other) result(check)
    implicit none
    class(gauss_legendre_field), intent(in) :: self
    class(gauss_legendre_field), intent(in) :: other
    logical :: check
    check = self%allocated .and. other%allocated
    if(.not.check) return
    check = check .and. (self%lmax == other%lmax)
    check = check .and. (self%ndim == other%ndim)
    return
  end function check_gauss_legendre_field

  subroutine assign_gauss_legendre_field(self,other)
    implicit none
    class(gauss_legendre_field), intent(inout) :: self
    class(gauss_legendre_field), intent(in) :: other
    if(.not.other%allocated) then
       call self%delete()
       return
    end if
    if(self%check(other)) then
       self%lmax = other%lmax
       self%ndim = other%ndim
       self%data = other%data               
    else
       self%allocated = .true.
       self%lmax = other%lmax
       self%ndim = other%ndim
       allocate(self%data(self%ndim))
       self%data = other%data
    end if
    return
  end subroutine assign_gauss_legendre_field


  function add_gauss_legendre_field(self,other) result(new)
    implicit none
    class(gauss_legendre_field), intent(in) :: self
    class(gauss_legendre_field), intent(in) :: other
    class(gauss_legendre_field), allocatable :: new
    allocate(new,source = self)
    new%allocated = .true.
    new%lmax = self%lmax
    new%ndim = self%ndim
    new%data = self%data + other%data
    return
  end function add_gauss_legendre_field

  function subtract_gauss_legendre_field(self,other) result(new)
    implicit none
    class(gauss_legendre_field), intent(in) :: self
    class(gauss_legendre_field), intent(in) :: other
    class(gauss_legendre_field), allocatable :: new
    allocate(new,source = self)
    new%allocated = .true.
    new%lmax = self%lmax
    new%ndim = self%ndim
    new%data = self%data - other%data
    return
  end function subtract_gauss_legendre_field

  function left_multiply_gauss_legendre_field(a,self) result(new)
    implicit none
    complex(dpc), intent(in) :: a
    class(gauss_legendre_field), intent(in) :: self
    class(gauss_legendre_field), allocatable :: new
    allocate(new,source = self)
    new%allocated = .true.
    new%lmax = self%lmax
    new%ndim = self%ndim
    call zcopy(new%ndim,self%data,1,new%data,1)
    call zscal(new%ndim,a,new%data,1)
    return
  end function left_multiply_gauss_legendre_field

  function right_multiply_gauss_legendre_field(self,a) result(new)
    implicit none
    complex(dpc), intent(in) :: a
    class(gauss_legendre_field), intent(in) :: self
    class(gauss_legendre_field), allocatable :: new
    allocate(new,source = self)
    new%allocated = .true.
    new%lmax = self%lmax
    new%ndim = self%ndim
    call zcopy(new%ndim,self%data,1,new%data,1)
    call zscal(new%ndim,a,new%data,1)
    return
  end function right_multiply_gauss_legendre_field

  function real_left_multiply_gauss_legendre_field(a,self) result(new)
    implicit none
    real(dp), intent(in) :: a
    class(gauss_legendre_field), intent(in) :: self
    class(gauss_legendre_field), allocatable :: new
    allocate(new,source = self)
    new%allocated = .true.
    new%lmax = self%lmax
    new%ndim = self%ndim
    call zcopy(new%ndim,self%data,1,new%data,1)
    call zdscal(new%ndim,a,new%data,1)
    return
  end function real_left_multiply_gauss_legendre_field

  function real_right_multiply_gauss_legendre_field(self,a) result(new)
    implicit none
    real(dp), intent(in) :: a
    class(gauss_legendre_field), intent(in) :: self
    class(gauss_legendre_field), allocatable :: new
    allocate(new,source = self)
    new%allocated = .true.
    new%lmax = self%lmax
    new%ndim = self%ndim
    call zcopy(new%ndim,self%data,1,new%data,1)
    call zdscal(new%ndim,a,new%data,1)
    return
  end function real_right_multiply_gauss_legendre_field

  subroutine scale_gauss_legendre_field(self,a)
    implicit none
    class(gauss_legendre_field), intent(inout) :: self
    complex(dpc), intent(in) :: a
    call zscal(self%ndim,a,self%data,1)
    return
  end subroutine scale_gauss_legendre_field

  subroutine real_scale_gauss_legendre_field(self,a)
    implicit none
    class(gauss_legendre_field), intent(inout) :: self
    real(dp), intent(in) :: a
    call zdscal(self%ndim,a,self%data,1)
    return
  end subroutine real_scale_gauss_legendre_field

  subroutine saxpy_gauss_legendre_field(self,a,other)
    implicit none
    class(gauss_legendre_field), intent(inout) :: self
    complex(dpc), intent(in) :: a
    class(gauss_legendre_field), intent(in) :: other
    call zaxpy(self%ndim,a,other%data,1,self%data,1)
    return
  end subroutine saxpy_gauss_legendre_field

  subroutine real_saxpy_gauss_legendre_field(self,a,other)
    implicit none
    class(gauss_legendre_field), intent(inout) :: self
    real(dp), intent(in) :: a
    class(gauss_legendre_field), intent(in) :: other
    complex(dpc) :: aloc
    aloc = a
    call zaxpy(self%ndim,aloc,other%data,1,self%data,1)
    return
  end subroutine real_saxpy_gauss_legendre_field
  
  !-----------------------------------------------!
  !                 scalar fields                 !
  !-----------------------------------------------!


  subroutine allocate_scalar_gauss_legendre_field(self,grid)
    implicit none
    class(scalar_gauss_legendre_field), intent(inout) :: self
    type(gauss_legendre_grid), intent(in) :: grid
    call self%delete()
    self%lmax = grid%lmax
    self%nth = grid%nth
    self%nph = grid%nph
    self%ndim = self%nph*self%nth 
    allocate(self%data(self%ndim))
    self%allocated = .true.
    return
  end subroutine allocate_scalar_gauss_legendre_field

  function index_scalar_gauss_legendre_field(self,iph,ith) result(i)
    implicit none
    class(scalar_gauss_legendre_field), intent(in) :: self
    integer(i4b), intent(in) :: iph,ith
    integer(i4b) :: i
    i = self%nph*(ith-1)+iph
    return
  end function index_scalar_gauss_legendre_field



  !-----------------------------------------------!
  !                 vector fields                 !
  !-----------------------------------------------!


  subroutine allocate_vector_gauss_legendre_field(self,grid)
    implicit none
    class(vector_gauss_legendre_field), intent(inout) :: self
    type(gauss_legendre_grid), intent(in) :: grid
    call self%delete()
    self%lmax = grid%lmax
    self%nth = grid%nth
    self%nph = grid%nph
    self%ndim = self%nph*self%nth 
    self%ndim = 3*self%ndim
    allocate(self%data(self%ndim))
    self%allocated = .true.
    return
  end subroutine allocate_vector_gauss_legendre_field

  function index_vector_gauss_legendre_field(self,alpha,iph,ith) result(i)
    implicit none
    class(vector_gauss_legendre_field), intent(in) :: self
    integer(i4b), intent(in) :: alpha,iph,ith
    integer(i4b) :: i
    i = self%nph*(ith-1)+iph+(alpha+1)*self%nph*self%nth
    return
  end function index_vector_gauss_legendre_field


  !-----------------------------------------------!
  !                real vector fields             !
  !-----------------------------------------------!


  subroutine allocate_real_vector_gauss_legendre_field(self,grid)
    implicit none
    class(real_vector_gauss_legendre_field), intent(inout) :: self
    type(gauss_legendre_grid), intent(in) :: grid
    call self%delete()
    self%lmax = grid%lmax
    self%nth = grid%nth
    self%nph = grid%nph
    self%ndim = self%nph*self%nth 
    self%ndim = 2*self%ndim
    allocate(self%data(self%ndim))
    self%allocated = .true.
    return
  end subroutine allocate_real_vector_gauss_legendre_field

  function index_real_vector_gauss_legendre_field(self,alpha,iph,ith) result(i)
    implicit none
    class(real_vector_gauss_legendre_field), intent(in) :: self
    integer(i4b), intent(in) :: alpha,iph,ith
    integer(i4b) :: i
    i = self%nph*(ith-1)+iph+(alpha+1)*self%nph*self%nth
    return
  end function index_real_vector_gauss_legendre_field

  
  !=======================================================================!
  !          procedures for the spherical_harmonic_expansion type         !
  !=======================================================================!


  !-----------------------------------------------!
  !               the abstract type               !
  !-----------------------------------------------!
  
  subroutine delete_spherical_harmonic_expansion(u)
    implicit none
    class(spherical_harmonic_expansion), intent(inout) :: u
    if(.not.u%allocated) return
    deallocate(u%data)
    u%allocated = .false.
    return
  end subroutine delete_spherical_harmonic_expansion

  subroutine set_sobolev_spherical_harmonic_expansion(u,s,mu)
    implicit none
    class(spherical_harmonic_expansion), intent(inout) :: u
    real(dp), intent(in) :: s,mu
    u%s  = s
    u%mu = mu
    return
  end subroutine set_sobolev_spherical_harmonic_expansion

  function check_spherical_harmonic_expansion(self,other) result(check)
    implicit none
    class(spherical_harmonic_expansion), intent(in) :: self
    class(spherical_harmonic_expansion), intent(in) :: other
    logical :: check
    check = self%allocated .and. other%allocated
    if(.not.check) return
    check = check .and. (self%lmax == other%lmax)
    check = check .and. (self%nmax == other%nmax)
    check = check .and. (self%ndim == other%ndim)
    return
  end function check_spherical_harmonic_expansion

  subroutine assign_spherical_harmonic_expansion(self,other)
    implicit none
    class(spherical_harmonic_expansion), intent(inout) :: self
    class(spherical_harmonic_expansion), intent(in) :: other
    if(.not.other%allocated) then
       call self%delete()
       return
    end if
    if(self%check(other)) then
       self%lmax = other%lmax
       self%nmax = other%nmax
       self%ndim = other%ndim
       self%s    = other%s
       self%mu   = other%mu       
       self%data = other%data               
    else
       self%allocated = .true.
       self%lmax = other%lmax
       self%nmax = other%nmax
       self%ndim = other%ndim
       self%s    = other%s
       self%mu   = other%mu
       allocate(self%data(self%ndim))
       self%data = other%data
    end if
    return
  end subroutine assign_spherical_harmonic_expansion

  function add_spherical_harmonic_expansion(self,other) result(new)
    implicit none
    class(spherical_harmonic_expansion), intent(in) :: self
    class(spherical_harmonic_expansion), intent(in) :: other
    class(spherical_harmonic_expansion), allocatable :: new
    allocate(new,source = self)
    new%allocated = .true.
    new%lmax = self%lmax
    new%nmax = self%nmax
    new%ndim = self%ndim
    new%s    = self%s
    new%mu   = self%mu
    new%data = self%data + other%data
    return
  end function add_spherical_harmonic_expansion

  function subtract_spherical_harmonic_expansion(self,other) result(new)
    implicit none
    class(spherical_harmonic_expansion), intent(in) :: self
    class(spherical_harmonic_expansion), intent(in) :: other
    class(spherical_harmonic_expansion), allocatable :: new
    allocate(new,source = self)
    new%allocated = .true.
    new%lmax = self%lmax
    new%nmax = self%nmax
    new%ndim = self%ndim
    new%s    = self%s
    new%mu   = self%mu
    new%data = self%data - other%data
    return
  end function subtract_spherical_harmonic_expansion

  function left_multiply_spherical_harmonic_expansion(a,self) result(new)
    implicit none
    complex(dpc), intent(in) :: a
    class(spherical_harmonic_expansion), intent(in) :: self
    class(spherical_harmonic_expansion), allocatable :: new
    allocate(new,source = self)
    new%allocated = .true.
    new%lmax = self%lmax
    new%nmax = self%nmax
    new%ndim = self%ndim
    new%s    = self%s
    new%mu   = self%mu
    call zcopy(self%ndim,self%data,1,new%data,1)
    call zscal(new%ndim,a,new%data,1)
    return
  end function left_multiply_spherical_harmonic_expansion

  function right_multiply_spherical_harmonic_expansion(self,a) result(new)
    implicit none
    complex(dpc), intent(in) :: a
    class(spherical_harmonic_expansion), intent(in) :: self
    class(spherical_harmonic_expansion), allocatable :: new
    allocate(new,source = self)
    new%allocated = .true.
    new%lmax = self%lmax
    new%nmax = self%nmax
    new%ndim = self%ndim
    new%s    = self%s
    new%mu   = self%mu
    call zcopy(self%ndim,self%data,1,new%data,1)
    call zscal(new%ndim,a,new%data,1)
    return
  end function right_multiply_spherical_harmonic_expansion

  function real_left_multiply_spherical_harmonic_expansion(a,self) result(new)
    implicit none
    real(dp), intent(in) :: a
    class(spherical_harmonic_expansion), intent(in) :: self
    class(spherical_harmonic_expansion), allocatable :: new
    allocate(new,source = self)
    new%allocated = .true.
    new%lmax = self%lmax
    new%nmax = self%nmax
    new%ndim = self%ndim
    new%s    = self%s
    new%mu   = self%mu
    call zcopy(self%ndim,self%data,1,new%data,1)
    call zdscal(new%ndim,a,new%data,1)
    return
  end function real_left_multiply_spherical_harmonic_expansion

  function real_right_multiply_spherical_harmonic_expansion(self,a) result(new)
    implicit none
    real(dp), intent(in) :: a
    class(spherical_harmonic_expansion), intent(in) :: self
    class(spherical_harmonic_expansion), allocatable :: new
    allocate(new,source = self)
    new%allocated = .true.
    new%lmax = self%lmax
    new%nmax = self%nmax
    new%ndim = self%ndim
    new%s    = self%s
    new%mu   = self%mu
    call zcopy(self%ndim,self%data,1,new%data,1)
    call zdscal(new%ndim,a,new%data,1)
    return
  end function real_right_multiply_spherical_harmonic_expansion

  subroutine scale_spherical_harmonic_expansion(self,a)
    implicit none
    class(spherical_harmonic_expansion), intent(inout) :: self
    complex(dpc), intent(in) :: a
    call zscal(self%ndim,a,self%data,1)
    return
  end subroutine scale_spherical_harmonic_expansion

  subroutine real_scale_spherical_harmonic_expansion(self,a)
    implicit none
    class(spherical_harmonic_expansion), intent(inout) :: self
    real(dp), intent(in) :: a
    call zdscal(self%ndim,a,self%data,1)
    return
  end subroutine real_scale_spherical_harmonic_expansion

  subroutine saxpy_spherical_harmonic_expansion(self,a,other)
    implicit none
    class(spherical_harmonic_expansion), intent(inout) :: self
    complex(dpc), intent(in) :: a
    class(spherical_harmonic_expansion), intent(in) :: other
    call zaxpy(self%ndim,a,other%data,1,self%data,1)
    return
  end subroutine saxpy_spherical_harmonic_expansion

  subroutine real_saxpy_spherical_harmonic_expansion(self,a,other)
    implicit none
    class(spherical_harmonic_expansion), intent(inout) :: self
    real(dp), intent(in) :: a
    class(spherical_harmonic_expansion), intent(in) :: other
    complex(dpc) :: aloc
    aloc = a
    call zaxpy(self%ndim,aloc,other%data,1,self%data,1)
    return
  end subroutine real_saxpy_spherical_harmonic_expansion


  !-----------------------------------------------!
  !                 scalar fields                 !
  !-----------------------------------------------!

  subroutine allocate_scalar_spherical_harmonic_expansion(self,lmax)
    implicit none
    class(scalar_spherical_harmonic_expansion), intent(inout) :: self
    integer(i4b), intent(in) :: lmax
    call self%delete()
    self%lmax = lmax
    self%nmax = 0
    self%ndim = (lmax+1)**2
    allocate(self%data(self%ndim))
    self%data = 0.0_dp
    self%allocated = .true.
    return
  end subroutine allocate_scalar_spherical_harmonic_expansion


  function index_scalar_spherical_harmonic_expansion(self,l,m) result(i)
    implicit none
    class(scalar_spherical_harmonic_expansion), intent(in) :: self
    integer(i4b), intent(in) :: l,m
    integer(i4b) :: i
    i = l**2
    if(m == 0) then
       i = i+1
    else if(m > 0) then
       i = i + 2*m
    else
       i = i -2*m+1
    end if
    return
  end function index_scalar_spherical_harmonic_expansion




  !-----------------------------------------------!
  !               real scalar fields              !
  !-----------------------------------------------!


  subroutine allocate_real_scalar_spherical_harmonic_expansion(self,lmax)
    implicit none
    class(real_scalar_spherical_harmonic_expansion), intent(inout) :: self
    integer(i4b), intent(in) :: lmax
    call self%delete()
    self%lmax = lmax
    self%nmax = 0
    self%ndim = ((lmax+1)*(lmax+2))/2
    allocate(self%data(self%ndim))
    self%data = 0.0_dp
    self%allocated = .true.
    return
  end subroutine allocate_real_scalar_spherical_harmonic_expansion

  
  function index_real_scalar_spherical_harmonic_expansion(self,l,m) result(i)
    implicit none
    class(real_scalar_spherical_harmonic_expansion), intent(in) :: self
    integer(i4b), intent(in) :: l,m
    integer(i4b) :: i
    i = (l*(l+1))/2 + m + 1
    return
  end function index_real_scalar_spherical_harmonic_expansion


  !-----------------------------------------------!
  !                 vector fields                 !
  !-----------------------------------------------!

  subroutine allocate_vector_spherical_harmonic_expansion(self,lmax)
    implicit none
    class(vector_spherical_harmonic_expansion), intent(inout) :: self
    integer(i4b), intent(in) :: lmax
    call self%delete()
    self%lmax = lmax
    self%nmax = 0
    self%ndim = (lmax+1)**2
    self%ndim = 3*self%ndim
    allocate(self%data(self%ndim))
    self%data = 0.0_dp
    self%allocated = .true.
    return
  end subroutine allocate_vector_spherical_harmonic_expansion


  function index_vector_spherical_harmonic_expansion(self,alpha,l,m) result(i)
    implicit none
    class(vector_spherical_harmonic_expansion), intent(in) :: self
    integer(i4b), intent(in) :: alpha,l,m
    integer(i4b) :: i
    i = l**2
    if(m == 0) then
       i = i+1
    else if(m > 0) then
       i = i + 2*m
    else
       i = i -2*m+1
    end if
    i = i + (1+alpha)*(self%ndim/3)
    return
  end function index_vector_spherical_harmonic_expansion


  !-----------------------------------------------!
  !             real vector fields                !
  !-----------------------------------------------!

  
  subroutine allocate_real_vector_spherical_harmonic_expansion(self,lmax)
    implicit none
    class(real_vector_spherical_harmonic_expansion), intent(inout) :: self
    integer(i4b), intent(in) :: lmax
    call self%delete()
    self%lmax = lmax
    self%nmax = 0
    self%ndim = ((lmax+1)*(lmax+2))/2
    self%ndim = 2*self%ndim
    allocate(self%data(self%ndim))
    self%data = 0.0_dp
    self%allocated = .true.
    return
  end subroutine allocate_real_vector_spherical_harmonic_expansion


  function index_real_vector_spherical_harmonic_expansion(self,alpha,l,m) result(i)
    implicit none
    class(real_vector_spherical_harmonic_expansion), intent(in) :: self
    integer(i4b), intent(in) :: alpha,l,m
    integer(i4b) :: i
    i = (l*(l+1))/2 + m + 1 + (alpha+1)*(self%ndim/2)
    return
  end function index_real_vector_spherical_harmonic_expansion



  
end module module_spherical_harmonics




