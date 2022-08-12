module module_spherical_harmonics
  
  use module_constants
  use module_error
  use module_special_functions, only : wigner_value
  use, intrinsic :: iso_c_binding

  type gauss_legendre_grid
     logical :: allocated = .false.
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





  
contains

  !=======================================================================!
  !                   procedures for the basic grid type                  !
  !=======================================================================!
  
  subroutine delete_gauss_legendre_grid(grid)
    implicit none
    class(gauss_legendre_grid), intent(inout) :: grid
    if(.not.grid%allocated) return
    deallocate(grid%th,grid%w)
    grid%allocated = .false.
    return
  end subroutine delete_gauss_legendre_grid

  subroutine build_gauss_legendre_grid(grid,lmax,nmax,flag)
    use module_special_functions
    use module_quadrature
    use module_fftw3    
    implicit none    
    class(gauss_legendre_grid), intent(inout) :: grid
    integer(i4b), intent(in) :: lmax
    integer(i4b), intent(in) :: nmax
    integer(C_INT), intent(in), optional :: flag
    
    integer(i4b) :: l,ith,n
    class(orthogonal_polynomial), allocatable :: poly
    type(gauss_quadrature) :: quad
    real(C_DOUBLE), pointer :: rin(:)
    complex(C_DOUBLE_COMPLEX), pointer ::  in(:)
    complex(C_DOUBLE_COMPLEX), pointer :: out(:),rout(:)
    type(C_PTR) :: pin,pout,plan1,plan2,plan3,plan4
    integer(C_INT) :: plan_flag
    
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

    ! make the FFTW3 pland
    if(present(flag)) then
       plan_flag = flag
    else
       plan_flag = FFTW_ESTIMATE
    end if
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
  








  
end module module_spherical_harmonics
