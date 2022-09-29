module module_spherical_harmonics
  
  use module_constants
  use module_special_functions
  use module_interp
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
     integer(i4b) :: ncoef_c
     integer(i4b) :: ncoef_r
     real(dp), dimension(:), allocatable :: th
     real(dp), dimension(:), allocatable :: w
     real(dp), dimension(:,:,:), allocatable :: dlm
     type(C_PTR) :: plan_forward
     type(C_PTR) :: plan_backward
     type(C_PTR) :: plan_r2c
     type(C_PTR) :: plan_c2r     
   contains
     procedure :: delete => delete_gauss_legendre_grid
     procedure :: build =>  build_gauss_legendre_grid
     procedure :: ph =>  ph_gauss_legednre_grid
     procedure, private :: SH_trans_gauss_legendre_grid
     procedure, private :: scalar_SH_trans_gauss_legendre_grid
     procedure, private :: real_SH_trans_gauss_legendre_grid
     procedure, private :: real_scalar_SH_trans_gauss_legendre_grid
     generic :: SH_trans => SH_trans_gauss_legendre_grid,        &
                            scalar_SH_trans_gauss_legendre_grid, & 
                            real_SH_trans_gauss_legendre_grid,   &
                            real_scalar_SH_trans_gauss_legendre_grid
     procedure, private :: SH_itrans_gauss_legendre_grid
     procedure, private :: scalar_SH_itrans_gauss_legendre_grid
     procedure, private :: real_SH_itrans_gauss_legendre_grid
     procedure, private :: real_scalar_SH_itrans_gauss_legendre_grid
     generic :: SH_itrans => SH_itrans_gauss_legendre_grid,        &
                             scalar_SH_itrans_gauss_legendre_grid, & 
                             real_SH_itrans_gauss_legendre_grid,   &
                             real_scalar_SH_itrans_gauss_legendre_grid
     procedure, private :: integrate_gauss_legendre_grid_complex
     procedure, private :: integrate_gauss_legendre_grid_real
     generic   :: integrate => integrate_gauss_legendre_grid_complex, &
                               integrate_gauss_legendre_grid_real
     procedure, nopass :: cindex => index_complex_spherical_harmonic_expansion
     procedure, nopass :: rindex => index_real_spherical_harmonic_expansion
     procedure, private :: filter_spherical_harmonic_coefficients_complex
     procedure, private :: filter_spherical_harmonic_coefficients_real
     generic :: filter => filter_spherical_harmonic_coefficients_complex, &
                          filter_spherical_harmonic_coefficients_real
  end type gauss_legendre_grid

  
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
    call fftw_destroy_plan(grid%plan_r2c)
    call fftw_destroy_plan(grid%plan_c2r)
    grid%allocated = .false.
    return
  end subroutine delete_gauss_legendre_grid

  
  subroutine build_gauss_legendre_grid(self,lmax,nmax,fftw_flag)
    use module_special_functions
    use module_quadrature
    use module_fftw3    
    implicit none    
    class(gauss_legendre_grid), intent(inout) :: self
    integer(i4b), intent(in) :: lmax
    integer(i4b), intent(in) :: nmax
    integer(C_INT), intent(in), optional :: fftw_flag

    
    integer(i4b) :: l,ith,n,nth,nph,cdim
    real(dp) :: fac
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
    call self%delete()

    ! store the basic parameters
    self%lmax = lmax
    self%nmax = nmax
    nth = lmax+1
    nph = 2*lmax
    self%nth = nth
    self%nph = nph
    self%ncoef_c = (lmax+1)**2
    self%ncoef_r = (lmax+1)*(lmax+2)/2
    
    ! make the quadrature points and weights
    allocate(self%th(lmax+1))
    allocate(self%w(lmax+1))
    call quad%set(lmax+1)
    self%th = acos(quad%points())
    self%w = quad%weights()

    ! make the wigner d-functions
    allocate(self%dlm(self%ncoef_r,-nmax:nmax,nth))
    do ith = 1,nth
       call set_wigner_array(self%th(ith),lmax,nmax,self%dlm(:,:,ith),norm=.true.)
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
    self%plan_forward  = plan1
    self%plan_backward = plan2
    self%plan_r2c      = plan3
    self%plan_c2r      = plan4
    
    self%allocated =.true.
    return
  end subroutine build_gauss_legendre_grid

  
  function ph_gauss_legednre_grid(self,iph) result(ph)
    implicit none
    class(gauss_legendre_grid), intent(in) :: self
    integer(i4b), intent(in) :: iph
    real(dp) :: ph
    ph = (iph-1)*twopi/self%nph
    return
  end function ph_gauss_legednre_grid
  


  !-----------------------------------------------------------------!
  !            spherical harmonic transformation routines           !
  !-----------------------------------------------------------------!

  
  subroutine SH_trans_gauss_legendre_grid(self,n,u,ulm)
    use module_fftw3
    implicit none
    class(gauss_legendre_grid), intent(in) :: self
    integer(i4b), intent(in) :: n
    complex(dpc), dimension(self%nph,self%nth), intent(in) :: u
    complex(dpc), dimension(self%ncoef_c), intent(out) :: ulm

    logical :: even
    integer(i4b) :: l,m,nth,nph,ith,ilm,klm,na,lmax
    real(dp) :: fac,w,sign
    complex(dpc) :: tmp
    complex(C_DOUBLE_COMPLEX), pointer :: in(:),out(:)
    type(C_PTR) :: pin,pout
    
    ! initialise the coefficients
    ulm = 0.0_dp
    
    ! get some parameters
    lmax = self%lmax
    nth = self%nth
    nph = self%nph
    
    ! set up the C pointers
    pin  = fftw_alloc_complex(int(nph, C_SIZE_T))
    pout = fftw_alloc_complex(int(nph, C_SIZE_T))
    call c_f_pointer(pin,   in, [nph])
    call c_f_pointer(pout, out, [nph])

    ! do the transformation
    na = abs(n)
    even = modulo(n,2) == 0

    if(even) then

       do ith = 1,nth
          in = u(:,ith)
          call fftw_execute_dft(self%plan_forward,in,out)
          ilm = na*na
          klm = na*(na+1)/2
          w = self%w(ith)*twopi/self%nph
          do l = na,lmax-1
             sign = 1.0_dp
             klm = klm+1
             tmp =  self%dlm(klm,n,ith)*w
             tmp = tmp*out(1)
             ilm = ilm+1
             ulm(ilm) = ulm(ilm) + tmp
             do m = 1,l
                klm = klm+1             
                tmp = self%dlm(klm,n,ith)*w
                tmp = tmp*out(m+1)
                ilm = ilm+1
                ulm(ilm) = ulm(ilm)+tmp             
                sign = -sign
                tmp = sign*self%dlm(klm,-n,ith)*w
                tmp = tmp*out(nph-m+1)
                ilm = ilm+1
                ulm(ilm) = ulm(ilm)+tmp
             end do
          end do
          l = lmax
          sign = 1.0_dp
          klm = klm+1
          tmp =  self%dlm(klm,n,ith)*w
          tmp = tmp*out(1)
          ilm = ilm+1
          ulm(ilm) = ulm(ilm) + tmp
          do m = 1,l-1
             klm = klm+1             
             tmp = self%dlm(klm,n,ith)*w
             tmp = tmp*out(m+1)
             ilm = ilm+1
             ulm(ilm) = ulm(ilm)+tmp             
             sign = -sign
             tmp = sign*self%dlm(klm,-n,ith)*w
             tmp = tmp*out(nph-m+1)
             ilm = ilm+1
             ulm(ilm) = ulm(ilm)+tmp             
          end do
          m = l
          klm = klm+1             
          tmp = self%dlm(klm,n,ith)*w
          tmp = tmp*out(m+1)
          ilm = ilm+1
          ulm(ilm) = ulm(ilm)+tmp
          sign = -sign
       end do

       
    else

       do ith = 1,nth
          in = u(:,ith)
          call fftw_execute_dft(self%plan_forward,in,out)
          ilm = na*na
          klm = na*(na+1)/2
          w = self%w(ith)*twopi/self%nph
          do l = na,lmax-1
             sign = -1.0_dp
             klm = klm+1
             tmp =  self%dlm(klm,n,ith)*w
             tmp = tmp*out(1)
             ilm = ilm+1
             ulm(ilm) = ulm(ilm) + tmp
             do m = 1,l
                klm = klm+1             
                tmp = self%dlm(klm,n,ith)*w
                tmp = tmp*out(m+1)
                ilm = ilm+1
                ulm(ilm) = ulm(ilm)+tmp             
                sign = -sign
                tmp = sign*self%dlm(klm,-n,ith)*w
                tmp = tmp*out(nph-m+1)
                ilm = ilm+1
                ulm(ilm) = ulm(ilm)+tmp
             end do
          end do
          l = lmax
          sign = -1.0_dp
          klm = klm+1
          tmp =  self%dlm(klm,n,ith)*w
          tmp = tmp*out(1)
          ilm = ilm+1
          ulm(ilm) = ulm(ilm) + tmp
          do m = 1,l-1
             klm = klm+1             
             tmp = self%dlm(klm,n,ith)*w
             tmp = tmp*out(m+1)
             ilm = ilm+1
             ulm(ilm) = ulm(ilm)+tmp             
             sign = -sign
             tmp = sign*self%dlm(klm,-n,ith)*w
             tmp = tmp*out(nph-m+1)
             ilm = ilm+1
             ulm(ilm) = ulm(ilm)+tmp             
          end do
          m = l
          klm = klm+1             
          tmp = self%dlm(klm,n,ith)*w
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


  subroutine scalar_SH_trans_gauss_legendre_grid(self,u,ulm)
    implicit none
    class(gauss_legendre_grid), intent(in) :: self
    complex(dpc), dimension(self%nph,self%nth), intent(in) :: u
    complex(dpc), dimension(self%ncoef_c), intent(out) :: ulm
    call SH_trans_gauss_legendre_grid(self,0,u,ulm)
    return
  end subroutine scalar_SH_trans_gauss_legendre_grid

  
  subroutine SH_itrans_gauss_legendre_grid(self,n,ulm,u)
    use module_fftw3
    implicit none
    class(gauss_legendre_grid), intent(in) :: self
    integer(i4b), intent(in) :: n
    complex(dpc), dimension(self%ncoef_c), intent(in) :: ulm
    complex(dpc), dimension(self%nph,self%nth), intent(out) :: u

    logical :: even
    integer(i4b) :: l,m,nth,nph,ith,ilm,klm,na,lmax,i,j
    real(dp) :: fac,w,sign
    complex(dpc) :: tmp
    complex(dpc), dimension(self%nph) :: lout
    complex(C_DOUBLE_COMPLEX), pointer :: in(:),out(:)
    type(C_PTR) :: pin,pout
    
    ! initialise the function
    u = 0.0_dp
    
    ! get some parameters
    lmax = self%lmax
    nth = self%nth
    nph = self%nph
    
    ! set up the C pointers
    pin  = fftw_alloc_complex(int(nph, C_SIZE_T))
    pout = fftw_alloc_complex(int(nph, C_SIZE_T))
    call c_f_pointer(pin,   in, [nph])
    call c_f_pointer(pout, out, [nph])

    ! do the transformation
    na = abs(n)
    even = modulo(n,2) == 0

    if(even) then

       do ith = 1,nth
          in = 0.0_dp
          ilm = na*na
          klm = na*(na+1)/2
          do l = na,lmax-1
             sign = 1.0_dp
             klm = klm+1
             tmp = self%dlm(klm,n,ith)
             ilm = ilm+1
             tmp = tmp*ulm(ilm)
             in(1) = in(1) + tmp
             do m = 1,l
                klm = klm+1
                tmp = self%dlm(klm,n,ith)
                ilm = ilm+1
                tmp = tmp*ulm(ilm)
                in(m+1) = in(m+1) + tmp
                sign = -sign
                tmp = sign*self%dlm(klm,-n,ith)
                ilm = ilm+1             
                tmp = tmp*ulm(ilm)
                in(nph-m+1) = in(nph-m+1) + tmp
             end do
          end do
          l = lmax
          sign = 1.0_dp
          klm = klm+1
          tmp = self%dlm(klm,n,ith)
          ilm = ilm+1
          tmp = tmp*ulm(ilm)
          in(1) = in(1) + tmp
          do m = 1,l-1
             klm = klm+1
             tmp = self%dlm(klm,n,ith)
             ilm = ilm+1
             tmp = tmp*ulm(ilm)
             in(m+1) = in(m+1) + tmp
             sign = -sign
             tmp = sign*self%dlm(klm,-n,ith)
             ilm = ilm+1             
             tmp = tmp*ulm(ilm)
             in(nph-m+1) = in(nph-m+1) + tmp
          end do
          m = l
          klm = klm+1
          tmp = self%dlm(klm,n,ith)
          ilm = ilm+1
          tmp = tmp*ulm(ilm)
          in(m+1) = in(m+1) + tmp       
          call fftw_execute_dft(self%plan_backward,in,out)
          u(:,ith) = out          
       end do

    else

       do ith = 1,nth
          in = 0.0_dp
          ilm = na*na
          klm = na*(na+1)/2
          do l = na,lmax-1
             sign = -1.0_dp
             klm = klm+1
             tmp = self%dlm(klm,n,ith)
             ilm = ilm+1
             tmp = tmp*ulm(ilm)
             in(1) = in(1) + tmp
             do m = 1,l
                klm = klm+1
                tmp = self%dlm(klm,n,ith)
                ilm = ilm+1
                tmp = tmp*ulm(ilm)
                in(m+1) = in(m+1) + tmp
                sign = -sign
                tmp = sign*self%dlm(klm,-n,ith)
                ilm = ilm+1             
                tmp = tmp*ulm(ilm)
                in(nph-m+1) = in(nph-m+1) + tmp
             end do
          end do
          l = lmax
          sign = -1.0_dp
          klm = klm+1
          tmp = self%dlm(klm,n,ith)
          ilm = ilm+1
          tmp = tmp*ulm(ilm)
          in(1) = in(1) + tmp
          do m = 1,l-1
             klm = klm+1
             tmp = self%dlm(klm,n,ith)
             ilm = ilm+1
             tmp = tmp*ulm(ilm)
             in(m+1) = in(m+1) + tmp
             sign = -sign
             tmp = sign*self%dlm(klm,-n,ith)
             ilm = ilm+1             
             tmp = tmp*ulm(ilm)
             in(nph-m+1) = in(nph-m+1) + tmp
          end do
          m = l
          klm = klm+1
          tmp = self%dlm(klm,n,ith)
          ilm = ilm+1
          tmp = tmp*ulm(ilm)
          in(m+1) = in(m+1) + tmp       
          call fftw_execute_dft(self%plan_backward,in,out)
          u(:,ith) = out
       end do

    end if

    ! nullify pointers
    call fftw_free(pin)
    call fftw_free(pout)
    in  => null()
    out => null()
    
    return
  end subroutine SH_itrans_gauss_legendre_grid


  subroutine scalar_SH_itrans_gauss_legendre_grid(self,ulm,u)
    use module_fftw3
    implicit none
    class(gauss_legendre_grid), intent(in) :: self
    complex(dpc), dimension(self%ncoef_c), intent(in) :: ulm
    complex(dpc), dimension(self%nph,self%nth), intent(out) :: u
    call SH_itrans_gauss_legendre_grid(self,0,ulm,u)
    return
  end subroutine scalar_SH_itrans_gauss_legendre_grid

  
  subroutine real_SH_trans_gauss_legendre_grid(self,n,u,ulm)
    use module_fftw3
    implicit none
    class(gauss_legendre_grid), intent(in) :: self
    integer(i4b), intent(in) :: n
    real(dp), dimension(self%nph,self%nth), intent(in) :: u
    complex(dpc), dimension(self%ncoef_r), intent(out) :: ulm

    integer(i4b) :: l,m,nth,nph,ith,ilm,klm,lmax
    real(dp) :: fac,w
    complex(dpc) :: tmp
    real(C_DOUBLE), pointer :: in(:)
    complex(C_DOUBLE_COMPLEX), pointer :: out(:)
    type(C_PTR) :: pin,pout
    
    ! initialise the coefficients
    ulm = 0.0_dp
    
    ! get some parameters
    lmax = self%lmax
    nth = self%nth
    nph = self%nph
    
    ! set up the C pointers
    pin  = fftw_alloc_real(int(nph, C_SIZE_T))
    pout = fftw_alloc_complex(int(nph/2+1, C_SIZE_T))
    call c_f_pointer(pin,   in, [nph])
    call c_f_pointer(pout, out, [nph/2+1])
    

    do ith = 1,nth
       in = u(:,ith)
       call fftw_execute_dft_r2c(self%plan_r2c,in,out)
       ilm = 0
       klm = 0
       w = self%w(ith)*twopi/self%nph
       do l = 0,lmax
          klm = klm+1
          tmp =  self%dlm(klm,n,ith)*w
          tmp = tmp*out(1)
          ilm = ilm+1
          ulm(ilm) = ulm(ilm) + tmp
          do m = 1,l
             klm = klm+1             
             tmp = self%dlm(klm,n,ith)*w
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


  subroutine real_scalar_SH_trans_gauss_legendre_grid(self,u,ulm)
    class(gauss_legendre_grid), intent(in) :: self
    real(dp), dimension(self%nph,self%nth), intent(in) :: u
    complex(dpc), dimension(self%ncoef_r), intent(out) :: ulm
    call real_SH_trans_gauss_legendre_grid(self,0,u,ulm)
    return
  end subroutine real_scalar_SH_trans_gauss_legendre_grid

  
  subroutine real_SH_itrans_gauss_legendre_grid(self,n,ulm,u)
    use module_fftw3
    implicit none
    class(gauss_legendre_grid), intent(in) :: self
    integer(i4b), intent(in) :: n
    complex(dpc), dimension(self%ncoef_r), intent(in) :: ulm
    real(dp), dimension(self%nph,self%nth), intent(out) :: u

    integer(i4b) :: l,m,nth,nph,ith,ilm,klm,lmax
    real(dp) :: fac,w
    complex(dpc) :: tmp
    real(C_DOUBLE), pointer :: out(:)
    complex(C_DOUBLE_COMPLEX), pointer :: in(:)
    type(C_PTR) :: pin,pout
    
    ! initialise the function
    u = 0.0_dp
    
    ! get some parameters
    lmax = self%lmax
    nth = self%nth
    nph = self%nph
    
    ! set up the C pointers
    pin = fftw_alloc_complex(int(nph/2+1, C_SIZE_T))
    pout  = fftw_alloc_real(int(nph, C_SIZE_T))
    call c_f_pointer(pin, in, [nph/2+1])
    call c_f_pointer(pout,   out, [nph])
    
    do ith = 1,nth
       in = 0.0_dp
       ilm = 0
       klm = 0
       do l = 0,lmax
          klm = klm+1
          tmp = self%dlm(klm,n,ith)
          ilm = ilm+1
          tmp = tmp*ulm(ilm)
          in(1) = in(1) + tmp
          do m = 1,l
             klm = klm+1
             tmp = self%dlm(klm,n,ith)
             ilm = ilm+1
             tmp = tmp*ulm(ilm)
             in(m+1) = in(m+1) + tmp
          end do
       end do
       call fftw_execute_dft_c2r(self%plan_c2r,in,out)
       u(:,ith) = out
    end do
    
    ! nullify the pointers
    call fftw_free(pin)
    call fftw_free(pout)
    in  => null()
    out => null()
    
    return
  end subroutine real_SH_itrans_gauss_legendre_grid

  subroutine real_scalar_SH_itrans_gauss_legendre_grid(self,ulm,u)
    class(gauss_legendre_grid), intent(in) :: self
    complex(dpc), dimension(self%ncoef_r), intent(in) :: ulm
    real(dp), dimension(self%nph,self%nth), intent(out) :: u
    call real_SH_itrans_gauss_legendre_grid(self,0,ulm,u)
    return
  end subroutine real_scalar_SH_itrans_gauss_legendre_grid


  complex(dpc) function integrate_gauss_legendre_grid_complex(self,u) result(int)
    class(gauss_legendre_grid), intent(in) :: self
    complex(dpc), dimension(self%nph,self%nth), intent(in) :: u
    integer(i4b) :: ith,iph,nth,nph
    real(dp) :: fac
    complex(dpc) :: tmp
    nth = self%nth
    nph = self%nph
    fac = twopi/nph
    area = 0.0_dp
    do ith = 1,nth
       tmp = sum(u(:,ith))
       int = int + tmp*self%w(ith)*fac
    end do    
    return
  end function integrate_gauss_legendre_grid_complex



  real(dp) function integrate_gauss_legendre_grid_real(self,u) result(int)
    class(gauss_legendre_grid), intent(in) :: self
    real(dp), dimension(self%nph,self%nth), intent(in) :: u
    integer(i4b) :: ith,iph,nth,nph
    real(dp) :: fac
    real(dp) :: tmp
    nth = self%nth
    nph = self%nph
    fac = twopi/nph
    area = 0.0_dp
    do ith = 1,nth
       tmp = sum(u(:,ith))
       int = int + tmp*self%w(ith)*fac
    end do    
    return
  end function integrate_gauss_legendre_grid_real

  
  
  
  integer(i4b) function index_complex_spherical_harmonic_expansion(l,m) result(i)
    implicit none
    integer(i4b), intent(in) :: l,m
    i = l**2
    if(m == 0) then
       i = i+1
    else if(m > 0) then
       i = i + 2*m
    else
       i = i -2*m+1
    end if
    return
  end function index_complex_spherical_harmonic_expansion

  integer(i4b) function index_real_spherical_harmonic_expansion(l,m) result(i)
    implicit none
    integer(i4b), intent(in) :: l,m
    i = l*(l+1)/2 + m + 1
    return
  end function index_real_spherical_harmonic_expansion

  subroutine filter_spherical_harmonic_coefficients_complex(grid,fac,u,v)
    class(gauss_legendre_grid), intent(in) :: grid
    complex(dpc), dimension(0:grid%lmax) :: fac
    complex(dpc), dimension(grid%ncoef_c), intent(in) :: u
    complex(dpc), dimension(grid%ncoef_c), intent(out) :: v
    integer(i4b) :: l,i1,i2
    do l = 0,grid%lmax
       i1 = grid%cindex(l,0)
       i2 = grid%cindex(l,-l)
       v(i1:i2) = fac(l)*u(i1:i2)
    end do    
    return
  end subroutine filter_spherical_harmonic_coefficients_complex

  subroutine filter_spherical_harmonic_coefficients_real(grid,fac,u,v)
    class(gauss_legendre_grid), intent(in) :: grid
    real(dp), dimension(0:grid%lmax) :: fac
    complex(dpc), dimension(grid%ncoef_r), intent(in) :: u
    complex(dpc), dimension(grid%ncoef_r), intent(out) :: v
    integer(i4b) :: l,i1,i2
    do l = 0,grid%lmax
       i1 = grid%rindex(l,0)
       i2 = grid%rindex(l,l)
       v(i1:i2) = fac(l)*u(i1:i2)
    end do    
    return
  end subroutine filter_spherical_harmonic_coefficients_real



end module module_spherical_harmonics




