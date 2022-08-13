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


  type spherical_harmonic_coefficient
     private
     logical :: allocated
     integer(i4b) :: lmax
     integer(i4b) :: nmax
     integer(i4b) :: ndim
     complex(dpc), dimension(:), allocatable :: data
     real(dp) :: s = 0.0_dp
     real(dp) :: mu = 1.0_dp
   contains
     procedure :: delete    => delete_spherical_harmonic_coefficient
     procedure :: allocate  => allocate_spherical_harmonic_coefficient
     procedure :: set_sobolev => set_sobolev_spherical_harmonic_coefficient
     procedure :: index     => index_spherical_harmonic_coefficient
     procedure :: get       => get_spherical_harmonic_coefficient
     procedure :: mslice    => mslice_spherical_harmonic_coefficient
     procedure :: nslice    => nslice_spherical_harmonic_coefficient
     procedure :: nmslice   => nmslice_spherical_harmonic_coefficient
     procedure, pass(u) :: complex_left_multiply_spherical_harmonic_coefficient
     procedure, pass(u) :: complex_right_multiply_spherical_harmonic_coefficient
     procedure, pass(u) :: real_left_multiply_spherical_harmonic_coefficient
     procedure, pass(u) :: real_right_multiply_spherical_harmonic_coefficient
     procedure, pass(u) :: real_sp_left_multiply_spherical_harmonic_coefficient
     procedure, pass(u) :: real_sp_right_multiply_spherical_harmonic_coefficient
     procedure, pass(u) :: int_left_multiply_spherical_harmonic_coefficient
     procedure, pass(u) :: int_right_multiply_spherical_harmonic_coefficient
     procedure, pass(u) :: dot_spherical_harmonic_coefficient
     generic, public :: operator(*) => complex_left_multiply_spherical_harmonic_coefficient,   &
                                       complex_right_multiply_spherical_harmonic_coefficient,  &
                                       real_left_multiply_spherical_harmonic_coefficient,      &
                                       real_right_multiply_spherical_harmonic_coefficient,     &
                                       real_sp_left_multiply_spherical_harmonic_coefficient,   &
                                       real_sp_right_multiply_spherical_harmonic_coefficient,  &
                                       int_left_multiply_spherical_harmonic_coefficient,       &
                                       int_right_multiply_spherical_harmonic_coefficient,      &
                                       dot_spherical_harmonic_coefficient     
     procedure :: add_spherical_harmonic_coefficient
     generic, public :: operator(+)  => add_spherical_harmonic_coefficient
     procedure :: subtract_spherical_harmonic_coefficient
     generic, public :: operator(-)  => subtract_spherical_harmonic_coefficient
     procedure :: conjg => conjugate_spherical_harmonic_coefficient
     procedure :: riesz_map => riesz_map_spherical_harmonic_coefficient
     procedure :: inverse_riesz_map => inverse_riesz_map_spherical_harmonic_coefficient
  end type spherical_harmonic_coefficient

  type, extends(spherical_harmonic_coefficient)  :: real_spherical_harmonic_coefficient
   contains
     procedure :: allocate  => allocate_real_spherical_harmonic_coefficient
  end type real_spherical_harmonic_coefficient



     
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
  !                 procedures for the coefficient types                  !
  !=======================================================================!
  

  subroutine delete_spherical_harmonic_coefficient(u)
    implicit none
    class(spherical_harmonic_coefficient), intent(inout) :: u
    if(.not.u%allocated) return
    deallocate(u%data)
    u%allocated = .false.
    return
  end subroutine delete_spherical_harmonic_coefficient


  subroutine allocate_spherical_harmonic_coefficient(u,lmax,nmax,s,mu)
    implicit none
    class(spherical_harmonic_coefficient), intent(inout) :: u
    integer(i4b), intent(in) :: lmax,nmax
    real(dp), intent(in),  optional :: s,mu
    call u%delete()
    u%lmax = lmax
    u%nmax = nmax
    u%ndim =   (nmax+1)*(2*nmax+1)*(2*nmax+3)/3           &
             + (2*nmax + 1)*(lmax-nmax)*(lmax+nmax+2)
    allocate(u%data(u%ndim))
    u%data = 1.0_dp
    if(present(s)) then
       u%s = s
    end if
    if(present(mu)) then
       u%mu = mu
    end if    
    u%allocated = .true.    
    return
  end subroutine allocate_spherical_harmonic_coefficient


  subroutine allocate_real_spherical_harmonic_coefficient(u,lmax,nmax,s,mu)
    implicit none
    class(real_spherical_harmonic_coefficient), intent(inout) :: u
    integer(i4b), intent(in) :: lmax,nmax
    real(dp), intent(in),  optional :: s,mu
    call u%delete()
    u%lmax = lmax
    u%nmax = nmax
    u%ndim =   (nmax+1)*(2*nmax+1)*(2*nmax+3)/3           &
             + (2*nmax + 1)*(lmax-nmax)*(lmax+nmax+2)
    allocate(u%data(u%ndim))
    u%data = 1.0_dp
    if(present(s)) then
       u%s = s
    end if
    if(present(mu)) then
       u%mu = mu
    end if    
    u%allocated = .true.    
    return
  end subroutine allocate_real_spherical_harmonic_coefficient


  subroutine set_sobolev_spherical_harmonic_coefficient(u,s,mu)
    implicit none
    class(spherical_harmonic_coefficient), intent(inout) :: u
    real(dp), intent(in) :: s,mu
    u%s  = s
    u%mu = mu
    return
  end subroutine set_sobolev_spherical_harmonic_coefficient
  
  function index_spherical_harmonic_coefficient(u,l,n,m) result(i)
    implicit none
    class(spherical_harmonic_coefficient), intent(in) :: u
    integer(i4b), intent(in) :: l,n,m
    integer(i4b) :: i

    call error( l < 0 .or. l > u%lmax, &
         'index_spherical_harmonic_coefficient','l out of range')
    call error(abs(n) > l .or. abs(n) > u%nmax, &
         'index_spherical_harmonic_coefficient','n out of range')
    call error(abs(m) > l, &
         'index_spherical_harmonic_coefficient','m out of range')

    if(l > u%nmax) then
       i = (u%nmax+1)*(2*u%nmax+1)*(2*u%nmax+3)/3 + (2*u%nmax + 1)*(l-1-u%nmax)*(l+u%nmax+1)
    else
       i = l*(4*l*l-1)/3
    end if
    if(l >= u%nmax) then
       i = i + (n+u%nmax)*(2*l+1) + m+l + 1
    else
       i = i + (n+l)*(2*l+1) + m+l+1
    end if
    
    return
  end function index_spherical_harmonic_coefficient


  function get_spherical_harmonic_coefficient(u,l,n,m) result(v)
    implicit none
    class(spherical_harmonic_coefficient), intent(in) :: u
    integer(i4b), intent(in) :: l,n,m
    complex(dp) :: v
    integer(i4b) ::i
    i = u%index(l,n,m)
    v = u%data(i)    
    return    
  end function get_spherical_harmonic_coefficient

  function mslice_spherical_harmonic_coefficient(u,l,n,m1,m2) result(v)
    implicit none
    class(spherical_harmonic_coefficient), intent(in) :: u
    integer(i4b), intent(in) :: l,n,m1,m2
    complex(dp), dimension(m2-m1+1) :: v
    integer(i4b) :: i1,i2
    call error(m1 > m2,'get_mslice_spherical_harmonic_coefficient',' m1 should be less than m2')
    i1 = u%index(l,n,m1)
    i2 = u%index(l,n,m2)
    v = u%data(i1:i2)    
    return
  end function mslice_spherical_harmonic_coefficient


  function nslice_spherical_harmonic_coefficient(u,l,n1,n2,m) result(v)
    implicit none
    class(spherical_harmonic_coefficient), intent(in) :: u
    integer(i4b), intent(in) :: l,n1,n2,m
    real(dp), dimension(n2-n1+1) :: v
    integer(i4b) :: i,n
    call error(n1 > n2,'get_nslice_spherical_harmonic_coefficient',' n1 should be less than n2')
    do n = n1,n2
       i = u%index(l,n,m)
       v(n+u%nmax+1) = u%data(i)    
    end do
    return
  end function nslice_spherical_harmonic_coefficient

  function nmslice_spherical_harmonic_coefficient(u,l,n1,n2,m1,m2) result(v)
    implicit none
    class(spherical_harmonic_coefficient), intent(in) :: u
    integer(i4b), intent(in) :: l,n1,n2,m1,m2
    complex(dp), dimension(n2-n1+1,m2-m1+1) :: v
    integer(i4b) :: i1,i2,n
    call error(m1 > m2,'get_mslice_spherical_harmonic_coefficient',' m1 should be less than m2')
    call error(n1 > n2,'get_nslice_spherical_harmonic_coefficient',' n1 should be less than n2')
    do n = n1,n2
       i1 = u%index(l,n,m1)
       i2 = u%index(l,n,m2)
       v(n+u%nmax+1,:) = u%data(i1:i2)    
    end do
    return
  end function nmslice_spherical_harmonic_coefficient


  function complex_left_multiply_spherical_harmonic_coefficient(a,u) result(v)
    implicit none
    complex(dpc), intent(in) :: a
    class(spherical_harmonic_coefficient), intent(in) :: u
    type(spherical_harmonic_coefficient) :: v
    v = u
    v%data = a*v%data
    return
  end function complex_left_multiply_spherical_harmonic_coefficient

  function complex_right_multiply_spherical_harmonic_coefficient(u,a) result(v)
    implicit none
    class(spherical_harmonic_coefficient), intent(in) :: u
    complex(dpc), intent(in) :: a
    type(spherical_harmonic_coefficient) :: v
    v = u
    v%data = a*v%data
    return
  end function complex_right_multiply_spherical_harmonic_coefficient
  
  function real_left_multiply_spherical_harmonic_coefficient(a,u) result(v)
    implicit none
    real(dp), intent(in) :: a
    class(spherical_harmonic_coefficient), intent(in) :: u
    type(spherical_harmonic_coefficient) :: v
    v = u
    v%data = a*v%data
    return
  end function real_left_multiply_spherical_harmonic_coefficient

  function real_right_multiply_spherical_harmonic_coefficient(u,a) result(v)
    implicit none
    class(spherical_harmonic_coefficient), intent(in) :: u
    real(dp), intent(in) :: a
    type(spherical_harmonic_coefficient) :: v
    v = u
    v%data = a*v%data
    return
  end function real_right_multiply_spherical_harmonic_coefficient

  
  function real_sp_left_multiply_spherical_harmonic_coefficient(a,u) result(v)
    implicit none
    real(sp), intent(in) :: a
    class(spherical_harmonic_coefficient), intent(in) :: u
    type(spherical_harmonic_coefficient) :: v
    v = u
    v%data = a*v%data
    return
  end function real_sp_left_multiply_spherical_harmonic_coefficient

  function real_sp_right_multiply_spherical_harmonic_coefficient(u,a) result(v)
    implicit none
    class(spherical_harmonic_coefficient), intent(in) :: u
    real(sp), intent(in) :: a
    type(spherical_harmonic_coefficient) :: v
    v = u
    v%data = a*v%data
    return
  end function real_sp_right_multiply_spherical_harmonic_coefficient

  function int_left_multiply_spherical_harmonic_coefficient(a,u) result(v)
    implicit none
    integer(i4b), intent(in) :: a
    class(spherical_harmonic_coefficient), intent(in) :: u
    type(spherical_harmonic_coefficient) :: v
    v = u
    v%data = a*v%data
    return
  end function int_left_multiply_spherical_harmonic_coefficient

  function int_right_multiply_spherical_harmonic_coefficient(u,a) result(v)
    implicit none
    class(spherical_harmonic_coefficient), intent(in) :: u
    integer(i4b), intent(in) :: a
    type(spherical_harmonic_coefficient) :: v
    v = u
    v%data = a*v%data
    return
  end function int_right_multiply_spherical_harmonic_coefficient
  

  function add_spherical_harmonic_coefficient(u,v) result(w)
    implicit none
    class(spherical_harmonic_coefficient), intent(in) :: u
    type(spherical_harmonic_coefficient), intent(in) :: v
    type(spherical_harmonic_coefficient) :: w
    w = u
    w%data = u%data + v%data
    return
  end function add_spherical_harmonic_coefficient

  function subtract_spherical_harmonic_coefficient(u,v) result(w)
    implicit none
    class(spherical_harmonic_coefficient), intent(in) :: u
    type(spherical_harmonic_coefficient), intent(in) :: v
    type(spherical_harmonic_coefficient) :: w
    w = u
    w%data = u%data - v%data    
    return
  end function subtract_spherical_harmonic_coefficient


  function conjugate_spherical_harmonic_coefficient(u) result(v)
    implicit none
    class(spherical_harmonic_coefficient), intent(in) :: u
    type(spherical_harmonic_coefficient) :: v
    v = u
    v%data = conjg(v%data)
    return
  end function conjugate_spherical_harmonic_coefficient

  function dot_spherical_harmonic_coefficient(v,u) result(f)
    implicit none
    type(spherical_harmonic_coefficient), intent(in) :: v
    class(spherical_harmonic_coefficient), intent(in) :: u
    complex(dpc) :: f

    integer(i4b) :: l,n,m,i
    real(dp) :: fac

    i = 0
    f = 0.0_dp
    do l = 0,u%lmax
       fac = (1.0_dp + u%mu**2*l*(l+1))**(u%s)
       do n = -min(l,u%nmax),min(l,u%nmax)
          do m = -l,l
             i = i+1
             f = f + fac*conjg(v%data(i))*u%data(i)             
          end do          
       end do       
    end do
    
    return
  end function dot_spherical_harmonic_coefficient

  function riesz_map_spherical_harmonic_coefficient(u) result(v)
    implicit none
    class(spherical_harmonic_coefficient), intent(in) :: u
    type(spherical_harmonic_coefficient) :: v

    integer(i4b) :: l,n,m,i
    real(dp) :: fac

    v = u    
    i = 0
    do l = 0,u%lmax
       fac = (1.0_dp + u%mu**2*l*(l+1))**(u%s)
       do n = -min(l,u%nmax),min(l,u%nmax)
          do m = -l,l
             i = i+1
             v%data(i) = fac*u%data(i)
          end do          
       end do       
    end do
    
    
    return
  end function riesz_map_spherical_harmonic_coefficient

  function inverse_riesz_map_spherical_harmonic_coefficient(u) result(v)
    implicit none
    class(spherical_harmonic_coefficient), intent(in) :: u
    type(spherical_harmonic_coefficient) :: v

    integer(i4b) :: l,n,m,i
    real(dp) :: fac

    v = u    
    i = 0
    do l = 0,u%lmax
       fac = (1.0_dp + u%mu**2*l*(l+1))**(u%s)
       do n = -min(l,u%nmax),min(l,u%nmax)
          do m = -l,l
             i = i+1
             v%data(i) = u%data(i)/fac
          end do          
       end do       
    end do
    
    
    return
  end function inverse_riesz_map_spherical_harmonic_coefficient
  

  !=======================================================================!
  !                     procedures for transformations                    !
  !=======================================================================!

  
  
end module module_spherical_harmonics
