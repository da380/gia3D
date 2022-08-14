module module_spherical_harmonics
  
  use module_constants
  use module_error
  use module_special_functions, only : wigner_value
  use, intrinsic :: iso_c_binding

  type gauss_legendre_grid
     private
     logical :: allocated = .false.
     logical :: precomp   = .true.
     integer(i4b) :: lmax
     integer(i4b) :: nmax
     integer(i4b) :: nth
     integer(i4b) :: nph
     real(dp), dimension(:), allocatable :: th
     real(dp), dimension(:), allocatable :: w
     type(wigner_value), dimension(:,:), allocatable :: dlm
     type(C_PTR) :: plan_forward
     type(C_PTR) :: plan_backward
     type(C_PTR) :: plan_r2c
     type(C_PTR) :: plan_c2r
   contains
     procedure :: delete => delete_gauss_legendre_grid
     procedure :: build =>  build_gauss_legendre_grid
     procedure :: ph    =>  ph_gauss_legednre_grid
  end type gauss_legendre_grid


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
  end type scalar_spherical_harmonic_expansion

  



     
contains

  !=======================================================================!
  !                     procedures for the grid type                      !
  !=======================================================================!
  
  subroutine delete_gauss_legendre_grid(grid)
    implicit none
    class(gauss_legendre_grid), intent(inout) :: grid
    if(.not.grid%allocated) return
    deallocate(grid%th,grid%w)
    grid%allocated = .false.
    return
  end subroutine delete_gauss_legendre_grid

  
  subroutine build_gauss_legendre_grid(grid,lmax,nmax,fftw_flag,precomp)
    use module_special_functions
    use module_quadrature
    use module_fftw3    
    implicit none    
    class(gauss_legendre_grid), intent(inout) :: grid
    integer(i4b), intent(in) :: lmax
    integer(i4b), intent(in) :: nmax
    integer(C_INT), intent(in), optional :: fftw_flag
    logical, intent(in), optional :: precomp
    
    integer(i4b) :: l,ith,n
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
    
    if(present(precomp)) then
       grid%precomp = precomp
    end if
       
    
    ! check it has not already been allocated
    call grid%delete()

    ! store the basic parameters
    grid%lmax = lmax
    grid%nmax = nmax
    grid%nth  = lmax+1
    grid%nph  = 2*lmax

    ! build the quadrature points and weights
    allocate(grid%th(lmax+1))
    allocate(grid%w(lmax+1))
    poly = legendre()
    call quad%set(lmax+1,poly)
    grid%th = acos(quad%points())
    grid%w = quad%weights()

    
    ! precompute the wigner d-functions
    if(grid%precomp) then
       allocate(grid%dlm(0:lmax,grid%nth/2+1))
       do ith = 1,grid%nth/2+1
          call grid%dlm(0,ith)%init(grid%th(ith),nmax,lmax)
          do l = 0,lmax-1
             call grid%dlm(l,ith)%next()
             grid%dlm(l+1,ith) = grid%dlm(l,ith)
             call grid%dlm(l,ith)%freeze()
          end do
          call grid%dlm(lmax,ith)%next()
          call grid%dlm(lmax,ith)%freeze()
       end do
    end if

    ! make the FFTW3 pland
    n = grid%nph
    pin  = fftw_alloc_complex(int(n, C_SIZE_T))
    pout = fftw_alloc_complex(int(n, C_SIZE_T))
    call c_f_pointer(pin,   in, [n])
    call c_f_pointer(pout, out, [n])
    plan1 = fftw_plan_dft_1d(n, in, out, FFTW_FORWARD, plan_flag)
    plan2 = fftw_plan_dft_1d(n, out, in, FFTW_BACKWARD, plan_flag)
    pin  = fftw_alloc_real(int(n, C_SIZE_T))
    pout = fftw_alloc_complex(int(n/2+1, C_SIZE_T))
    call c_f_pointer(pin,   rin, [n])
    call c_f_pointer(pout, rout, [n/2+1])
    plan3 = fftw_plan_dft_r2c_1d(n,rin,rout,plan_flag);
    plan4 = fftw_plan_dft_c2r_1d(n,rout,rin,plan_flag);
    grid%plan_forward = plan1
    grid%plan_backward = plan2
    grid%plan_r2c = plan3
    grid%plan_c2r = plan4
    call fftw_free(pin)
    call fftw_free(pout)
    
    grid%allocated =.true.
    return
  end subroutine build_gauss_legendre_grid

  function ph_gauss_legednre_grid(grid,iph) result(ph)
    implicit none
    class(gauss_legendre_grid), intent(in) :: grid
    integer(i4b), intent(in) :: iph
    real(dp) :: ph
    ph = (iph-1)*twopi/grid%nph
    return
  end function ph_gauss_legednre_grid

  !=======================================================================!
  !          procedures for the spherical_harmonic_expansion type         !
  !=======================================================================!
  
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
    call error(.not.self%check(other),'add_spherical_harmonic_expansion','incompatible inputs')
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
    call error(.not.self%check(other),'subtract_spherical_harmonic_expansion','incompatible inputs')
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
    call error(.not.self%allocated,'left_multiply_spherical_harmonic_expansion','bad input')
    allocate(new,source = self)
    new%allocated = .true.
    new%lmax = self%lmax
    new%nmax = self%nmax
    new%ndim = self%ndim
    new%s    = self%s
    new%mu   = self%mu
    new%data = a*self%data
    return
  end function left_multiply_spherical_harmonic_expansion

  function right_multiply_spherical_harmonic_expansion(self,a) result(new)
    implicit none
    complex(dpc), intent(in) :: a
    class(spherical_harmonic_expansion), intent(in) :: self
    class(spherical_harmonic_expansion), allocatable :: new
    call error(.not.self%allocated,'left_multiply_spherical_harmonic_expansion','bad input')
    allocate(new,source = self)
    new%allocated = .true.
    new%lmax = self%lmax
    new%nmax = self%nmax
    new%ndim = self%ndim
    new%s    = self%s
    new%mu   = self%mu
    new%data = a*self%data
    return
  end function right_multiply_spherical_harmonic_expansion

  function real_left_multiply_spherical_harmonic_expansion(a,self) result(new)
    implicit none
    real(dp), intent(in) :: a
    class(spherical_harmonic_expansion), intent(in) :: self
    class(spherical_harmonic_expansion), allocatable :: new
    call error(.not.self%allocated,'left_multiply_spherical_harmonic_expansion','bad input')
    allocate(new,source = self)
    new%allocated = .true.
    new%lmax = self%lmax
    new%nmax = self%nmax
    new%ndim = self%ndim
    new%s    = self%s
    new%mu   = self%mu
    new%data = a*self%data
    return
  end function real_left_multiply_spherical_harmonic_expansion

  function real_right_multiply_spherical_harmonic_expansion(self,a) result(new)
    implicit none
    real(dp), intent(in) :: a
    class(spherical_harmonic_expansion), intent(in) :: self
    class(spherical_harmonic_expansion), allocatable :: new
    call error(.not.self%allocated,'left_multiply_spherical_harmonic_expansion','bad input')
    allocate(new,source = self)
    new%allocated = .true.
    new%lmax = self%lmax
    new%nmax = self%nmax
    new%ndim = self%ndim
    new%s    = self%s
    new%mu   = self%mu
    new%data = a*self%data
    return
  end function real_right_multiply_spherical_harmonic_expansion

  subroutine scale_spherical_harmonic_expansion(self,a)
    implicit none
    class(spherical_harmonic_expansion), intent(inout) :: self
    complex(dpc), intent(in) :: a
    call error(.not.self%allocated,'scale_spherical_harmonic_expansion','not allocated')    
    self%data = a*self%data
    return
  end subroutine scale_spherical_harmonic_expansion

  subroutine real_scale_spherical_harmonic_expansion(self,a)
    implicit none
    class(spherical_harmonic_expansion), intent(inout) :: self
    real(dp), intent(in) :: a
    call error(.not.self%allocated,'real_scale_spherical_harmonic_expansion','not allocated')
    self%data = a*self%data
    return
  end subroutine real_scale_spherical_harmonic_expansion

  subroutine saxpy_spherical_harmonic_expansion(self,a,other)
    implicit none
    class(spherical_harmonic_expansion), intent(inout) :: self
    complex(dpc), intent(in) :: a
    class(spherical_harmonic_expansion), intent(in) :: other
    call error(.not.self%check(other),'saxpy_spherical_harmonic_expansion', &
                                      'incompatible inputs')
    self%data = a*self%data + other%data
    return
  end subroutine saxpy_spherical_harmonic_expansion

  subroutine real_saxpy_spherical_harmonic_expansion(self,a,other)
    implicit none
    class(spherical_harmonic_expansion), intent(inout) :: self
    real(dp), intent(in) :: a
    class(spherical_harmonic_expansion), intent(in) :: other
    call error(.not.self%check(other),'real_saxpy_spherical_harmonic_expansion', &
                                      'incompatible inputs')
    self%data = a*self%data + other%data
    return
  end subroutine real_saxpy_spherical_harmonic_expansion

  


  !=======================================================================!
  !                      procedures for scalar fields                     !
  !=======================================================================!


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



  
!  function index_spherical_harmonic_expansion(u,l,n,m) result(i)
!    implicit none
!    class(spherical_harmonic_expansion), intent(in) :: u
!    integer(i4b), intent(in) :: l,n,m
!    integer(i4b) :: i

!    call error( l < 0 .or. l > u%lmax, &
!         'index_spherical_harmonic_expansion','l out of range')
!    call error(abs(n) > l .or. abs(n) > u%nmax, &
!         'index_spherical_harmonic_expansion','n out of range')
!    call error(abs(m) > l, &
!         'index_spherical_harmonic_expansion','m out of range')

!    if(l > u%nmax) then
!       i = (u%nmax+1)*(2*u%nmax+1)*(2*u%nmax+3)/3 + (2*u%nmax + 1)*(l-1-u%nmax)*(l+u%nmax+1)
!    else
!       i = l*(4*l*l-1)/3
!    end if
!    if(l >= u%nmax) then
!       i = i + (n+u%nmax)*(2*l+1) + m+l + 1
!    else
!       i = i + (n+l)*(2*l+1) + m+l+1
!    end if
    
!    return
!  end function index_spherical_harmonic_expansion


!  subroutine allocate_spherical_harmonic_expansion(u,lmax,nmax,s,mu)
!    implicit none
!    class(spherical_harmonic_expansion), intent(inout) :: u
!    integer(i4b), intent(in) :: lmax,nmax
!    real(dp), intent(in),  optional :: s,mu
!    call u%delete()
!    u%lmax = lmax
!    u%nmax = nmax
!    u%ndim =   (nmax+1)*(2*nmax+1)*(2*nmax+3)/3           &
!             + (2*nmax + 1)*(lmax-nmax)*(lmax+nmax+2)
!    allocate(u%data(u%ndim))
!    u%data = 1.0_dp
!    if(present(s)) then
!       u%s = s
!    end if
!    if(present(mu)) then
!       u%mu = mu
!    end if    
!    u%allocated = .true.    
!    return
!  end subroutine allocate_spherical_harmonic_expansion

  

  !  function get_spherical_harmonic_expansion(u,l,n,m) result(v)
!    implicit none
!    class(spherical_harmonic_expansion), intent(in) :: u
!    integer(i4b), intent(in) :: l,n,m
!    complex(dp) :: v
!    integer(i4b) ::i
!    i = u%index(l,n,m)
!    v = u%data(i)    
!    return    
!  end function get_spherical_harmonic_expansion

!  function mslice_spherical_harmonic_expansion(u,l,n,m1,m2) result(v)
!    implicit none
!    class(spherical_harmonic_expansion), intent(in) :: u
!    integer(i4b), intent(in) :: l,n,m1,m2
!    complex(dp), dimension(m2-m1+1) :: v
!    integer(i4b) :: i1,i2
!    call error(m1 > m2,'get_mslice_spherical_harmonic_expansion',' m1 should be less than m2')
!    i1 = u%index(l,n,m1)
!    i2 = u%index(l,n,m2)
!    v = u%data(i1:i2)    
!    return
!  end function mslice_spherical_harmonic_expansion


!  function nslice_spherical_harmonic_expansion(u,l,n1,n2,m) result(v)
!    implicit none
!    class(spherical_harmonic_expansion), intent(in) :: u
!    integer(i4b), intent(in) :: l,n1,n2,m
!    real(dp), dimension(n2-n1+1) :: v
!    integer(i4b) :: i,n
!    call error(n1 > n2,'get_nslice_spherical_harmonic_expansion',' n1 should be less than n2')
!    do n = n1,n2
!       i = u%index(l,n,m)
!       v(n+u%nmax+1) = u%data(i)    
!    end do
!    return
!  end function nslice_spherical_harmonic_expansion

!  function nmslice_spherical_harmonic_expansion(u,l,n1,n2,m1,m2) result(v)
!    implicit none
!    class(spherical_harmonic_expansion), intent(in) :: u
!    integer(i4b), intent(in) :: l,n1,n2,m1,m2
!    complex(dp), dimension(n2-n1+1,m2-m1+1) :: v
!    integer(i4b) :: i1,i2,n
!    call error(m1 > m2,'get_mslice_spherical_harmonic_expansion',' m1 should be less than m2')
!    call error(n1 > n2,'get_nslice_spherical_harmonic_expansion',' n1 should be less than n2')
!    do n = n1,n2
!       i1 = u%index(l,n,m1)
!       i2 = u%index(l,n,m2)
!       v(n+u%nmax+1,:) = u%data(i1:i2)    
!    end do
!    return
!  end function nmslice_spherical_harmonic_expansion


  !function conjugate_spherical_harmonic_expansion(u) result(v)
  !  implicit none
  !  class(spherical_harmonic_expansion), intent(in) :: u
  !  class(spherical_harmonic_expansion), allocatable :: v
  !  v = u
  !  v%data = conjg(v%data)
  !  return
  !end function conjugate_spherical_harmonic_expansion
  
end module module_spherical_harmonics




