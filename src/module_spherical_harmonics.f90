module module_spherical_harmonics
  
  use module_constants
  use module_error


  type gauss_legendre_grid
     logical :: allocated = .false.
     integer(i4b) :: lmax
     integer(i4b) :: nth
     integer(i4b) :: nph
     integer(i4b) :: nplm
     real(dp) :: dph
     real(dp), dimension(:), allocatable :: th
     real(dp), dimension(:), allocatable :: w
   contains
     procedure :: delete => delete_gauss_legendre_grid
     procedure :: build =>  build_gauss_legendre_grid
     procedure :: precompute_legendre 
     procedure :: ph    =>  ph_gauss_legednre_grid
  end type gauss_legendre_grid

  type, extends(gauss_legendre_grid) :: scalar_gauss_legendre_grid
     integer(i4b) :: ncoef
     complex(dpc), dimension(:,:), allocatable :: val
   contains
     procedure :: delete => delete_scalar_gauss_legendre_grid
     procedure :: build => build_scalar_gauss_legendre_grid
     procedure :: plan => plan_scalar_gauss_legendre_grid
  end type scalar_gauss_legendre_grid 




  
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

  subroutine build_gauss_legendre_grid(grid,lmax)
    use module_special_functions
    use module_quadrature
    implicit none    
    class(gauss_legendre_grid), intent(inout) :: grid
    integer(i4b), intent(in) :: lmax
    class(orthogonal_polynomial), allocatable :: poly
    type(gauss_quadrature) :: quad
    call grid%delete()    
    grid%lmax = lmax
    grid%nth  = lmax+1
    grid%nph  = 2*lmax    
    grid%dph = twopi/grid%nph
    allocate(grid%th(lmax+1))
    allocate(grid%w(lmax+1))
    poly = legendre()
    call quad%set(lmax+1,poly)
    grid%th = acos(quad%points())
    grid%w = quad%weights()    
    grid%allocated =.true.
    return
  end subroutine build_gauss_legendre_grid

  function ph_gauss_legednre_grid(grid,iph) result(ph)
    implicit none
    class(gauss_legendre_grid), intent(in) :: grid
    integer(i4b), intent(in) :: iph
    real(dp) :: ph
    ph = (iph-1)*grid%dph
    return
  end function ph_gauss_legednre_grid
  


  function precompute_legendre(grid) result(xlm)
    use module_special_functions
    implicit none
    class(gauss_legendre_grid), intent(in) :: grid
    real(dp), dimension(:,:), allocatable :: xlm
    integer(i4b) :: lmax,nth,ith,ndim,l,i1,i2
    real(dp) :: th
    type(legendre_value) :: p    
    lmax = grid%lmax
    nth = grid%nth
    ndim = grid%nplm
    allocate(xlm(ndim,nth))
    do ith = 1,nth
       th = grid%th(ith)
       call p%init(th,lmax)
       i1 = 1       
       do l = 0,lmax
          i2 = i1+l          
          call p%next()
          xlm(i1:i2,ith) = p%get(0,l)
          i1 = i2+1
       end do
    end do
    return
  end function precompute_legendre

  

  !=======================================================================!
  !                 procedures for complex scalar fields                  !
  !=======================================================================!

  subroutine delete_scalar_gauss_legendre_grid(grid)
    implicit none    
    class(scalar_gauss_legendre_grid), intent(inout) :: grid
    if(.not.grid%allocated) return
    call grid%gauss_legendre_grid%delete()
    deallocate(grid%val)
    grid%allocated = .false.
    return
  end subroutine delete_scalar_gauss_legendre_grid


  subroutine build_scalar_gauss_legendre_grid(grid,lmax)
    implicit none    
    class(scalar_gauss_legendre_grid), intent(inout) :: grid
    integer(i4b), intent(in) :: lmax
    call grid%delete()
    call grid%gauss_legendre_grid%build(lmax)
    grid%ncoef = (lmax+1)**2
    grid%nplm = lmax*(lmax+1)/2 + lmax + 1
    allocate(grid%val(grid%nph,grid%nth))
    grid%val = 0.0_dp    
    grid%allocated = .true.
    return
  end subroutine build_scalar_gauss_legendre_grid


  
  
  function plan_scalar_gauss_legendre_grid(u,forward,flag) result(plan)
    use module_fftw3
    implicit none    
    class(scalar_gauss_legendre_grid), intent(in) :: u
    logical, intent(in) :: forward
    integer(C_INT), intent(in), optional :: flag
    integer(C_INT) :: plan_flag
    type(C_PTR) :: plan

    integer(i4b) :: n
    complex(C_DOUBLE_COMPLEX), pointer ::  in(:)
    complex(C_DOUBLE_COMPLEX), pointer :: out(:)
    type(C_PTR) :: pin,pout
    
    ! set up C pointers
    n = u%nph
    pin  = fftw_alloc_complex(int(n, C_SIZE_T))
    pout = fftw_alloc_complex(int(n, C_SIZE_T))
    call c_f_pointer(pin,   in, [n])
    call c_f_pointer(pout, out, [n])

    ! make the plan
    if(present(flag)) then
       plan_flag = flag
    else
       plan_flag = FFTW_MEASURE
    end if
    
    if(forward) then
       plan = fftw_plan_dft_1d(n, in, out, FFTW_FORWARD, plan_flag)
    else
       plan = fftw_plan_dft_1d(n, out, in, FFTW_BACKWARD, plan_flag)
    end if

    ! delete the pointers
    call fftw_free(pin)
    call fftw_free(pout)
    in  => null()
    out => null()
    
    return
  end function plan_scalar_gauss_legendre_grid

  

 subroutine trans_scalar_gauss_legendre_grid(u,plan,xlm,ulm)
    use module_fftw3
    use module_special_functions
    implicit none    
    class(scalar_gauss_legendre_grid), intent(in) :: u
    type(C_PTR), intent(in) :: plan
    real(dp), dimension(u%nplm,u%nth), intent(in) :: xlm    
    complex(dpc), dimension(:), intent(out) :: ulm

    integer(i4b) :: n,ith,l,m,ilm,im,ix,sign
    real(dp) :: fac
    complex(C_DOUBLE_COMPLEX), pointer :: in(:),out(:)
    type(C_PTR) :: pin,pout
     
    ! set up the C pointers
    pin   = fftw_alloc_complex(int(u%nph, C_SIZE_T))
    pout  = fftw_alloc_complex(int(u%nph, C_SIZE_T))
    call c_f_pointer(pin,   in, [u%nph])
    call c_f_pointer(pout,  out, [u%nph])



    
    ! initialise the coefficients
    ulm = 0.0_dp

    
    ! loop over colatitude
    do ith = 1,u%nth


       ! transform the ith column
       in = u%val(:,ith)
       call fftw_execute_dft(plan, in, out)
       
       ilm = 0
       ix = 0
       do l = 0,u%lmax

          ! deal with m = 0
          ilm = ilm+1
          im  = 1
          ix = ix+1
          ulm(ilm) = ulm(ilm) + out(im)*xlm(ix,ith)*u%w(ith)


          ! deal with m /= 0
          sign = 1
          do m = 1,l

             ! positive m
             ilm = ilm+1
             im = m+1
             ix = ix+1
             ulm(ilm) = ulm(ilm) + out(im)*xlm(ix,ith)*u%w(ith)

             ! negative m
             ilm = ilm+1
             im = -m+u%nph+1
             sign = -sign
             ulm(ilm) = ulm(ilm) + sign*out(im)*xlm(ix,ith)*u%w(ith)

                          
          end do
          
       end do
       
    end do


    ulm = (twopi/u%nph)*ulm


    call fftw_free(pin)
    call fftw_free(pout)
    in  => null() 
    out => null()


    
    return
  end subroutine trans_scalar_gauss_legendre_grid









  
end module module_spherical_harmonics
