module module_random_fields

  use module_constants
  use module_util
  use module_error
  use module_quadrature
  use module_interp
  use module_fftw3
  use module_spherical_harmonics
  use, intrinsic :: iso_c_binding
  implicit none
 
  type, abstract :: GRF_1D
   contains
     procedure(GRF_1D_delete),  deferred :: delete
     procedure(GRF_1D_realise), deferred :: realise
     procedure(GRF_1D_corr), deferred :: corr
  end type GRF_1D

  abstract interface
     subroutine GRF_1D_delete(self)
       import :: GRF_1D
       class(GRF_1D), intent(inout) :: self
     end subroutine GRF_1D_delete
     subroutine GRF_1D_realise(self,fun)
       use module_interp
       import :: GRF_1D
       class(GRF_1D), intent(inout) :: self
       class(interp_1D), intent(inout) :: fun
     end subroutine GRF_1D_realise
     subroutine GRF_1D_corr(self,y,fun)
       use module_constants
       use module_interp
       import :: GRF_1D
       class(GRF_1D), intent(in) :: self
       real(dp), intent(in) :: y
       class(interp_1D), intent(inout) :: fun
     end subroutine GRF_1D_corr
  end interface

  type, extends(GRF_1D) :: GRF_1D_SEM
     logical :: allocated = .false.
     integer(i4b) :: ndim
     integer(i4b) :: mdim
     real(dp), dimension(:), allocatable :: Q
     real(dp), dimension(:), allocatable :: x
     real(dp), dimension(:,:), allocatable :: evec     
   contains
     procedure :: delete => delete_GRF_1D_SEM
     procedure :: realise => realise_GRF_1D_SEM
     procedure :: corr => corr_GRF_1D_SEM
  end type GRF_1D_SEM

  interface GRF_1D_SEM
     procedure :: build_GRF_1D_SEM
  end interface GRF_1D_SEM

  type, extends(GRF_1D) :: GRF_1D_Fourier
     logical :: allocated = .false.
     integer(i4b) :: n,i1,i2     
     real(dp), dimension(:), allocatable :: x
     real(dp), dimension(:), allocatable :: Q
     type(C_PTR) :: plan_c2r
   contains
     procedure :: delete => delete_GRF_1D_Fourier
     procedure :: realise => realise_GRF_1D_Fourier
     procedure :: corr => corr_GRF_1D_Fourier
  end type GRF_1D_Fourier

  interface GRF_1D_Fourier
          procedure :: build_GRF_1D_Fourier
  end interface GRF_1D_Fourier

   
  type gaussain_random_scalar_field_sphere
     logical :: allocated = .false.
     integer(i4b) :: lmax
     real(dp), dimension(:), allocatable :: Q_l
     complex(dpc), dimension(:), allocatable :: mlm
     complex(dpc), dimension(:), allocatable :: ulm
   contains
     procedure :: delete => delete_gaussain_random_scalar_field_sphere
     procedure :: build => build_gaussain_random_scalar_field_sphere
     procedure :: set_mean_coefficients_gaussain_random_scalar_field_sphere
     procedure :: set_mean_values_gaussain_random_scalar_field_sphere
     generic   :: set_mean => set_mean_coefficients_gaussain_random_scalar_field_sphere, &
                              set_mean_values_gaussain_random_scalar_field_sphere
     procedure :: realise => realise_gaussain_random_scalar_field_sphere
  end type gaussain_random_scalar_field_sphere

  
contains


  !====================================================!
  !                 1D fields using SEM                !
  !====================================================!
  
  subroutine delete_GRF_1D_SEM(self)
    class(GRF_1D_SEM), intent(inout) :: self
    if(.not.self%allocated) return
    deallocate(self%Q,   & 
               self%x,    &          
               self%evec)
    self%allocated = .false.
    return
  end subroutine delete_GRF_1D_SEM


  subroutine realise_GRF_1D_SEM(self,fun)
    class(GRF_1D_SEM), intent(inout) :: self
    class(interp_1D), intent(inout) :: fun

    integer(i4b) :: ispec,inode,i,j
    real(dp), dimension(self%ndim) :: un
    real(dp), dimension(self%mdim) :: u

    call random_seed()
    call normal_random_variable(un)
    un = un*sqrt(self%Q)
    u = 0.0_dp
    do i = 1,self%ndim
      u(:) = u(:) + un(i)*self%evec(:,i)
    end do
    call fun%set(self%x,u)    
    return
  end subroutine realise_GRF_1D_SEM

  
  subroutine corr_GRF_1D_SEM(self,y,fun)
    class(GRF_1D_SEM), intent(in) :: self
    real(dp), intent(in) :: y
    class(interp_1D), intent(inout) :: fun

    integer(i4b) :: iy,i,n
    real(dp) :: v1,v2,v,y1,y2
    real(dp), dimension(self%mdim) :: u

    iy = bisect_list(self%x,y)
    y1 = self%x(iy)
    y2 = self%x(iy+1)
    u = 0.0_dp
    do i = 1,self%ndim
       v1 = self%evec(iy,i)
       v2 = self%evec(iy+1,i)
       v = v1 + (v2-v1)*(y-y1)/(y2-y1)
       u(:) = u(:) + self%Q(i)*v*self%evec(:,i)
    end do
    call fun%set(self%x,u)    
    
    return
  end subroutine corr_GRF_1D_SEM
  

  type(GRF_1D_SEM) function build_GRF_1D_SEM(x1,x2,lambda,s,sigma,eps) &
                                             result(rfun)
    
    real(dp), intent(in) :: x1,x2,lambda,s,sigma
    real(dp), intent(in), optional :: eps

    integer(i4b), parameter :: ngll = 5
    integer(i4b), parameter :: npad = 10
    real(dp), parameter :: eps_default = 1.e-5_dp

    integer(i4b) :: inode,ispec,count,kda,ldab,kdb,ldbb, &
                    ndim,i,j,k,jnode,knode,info,nmax,    &
                    nspec,ispec1,ispec2,i1,i2
    integer(i4b), dimension(:,:), allocatable :: ibool
    real(dp) :: x11,x22,xl,xr,dx,fac,sum,x,eps_loc,x0
    real(dp), dimension(:), allocatable :: eval,work,jac
    real(dp), dimension(:,:), allocatable :: aa,bb,evec,hp,xx
    type(gauss_lobatto_quadrature) :: quad

        
    if(present(eps)) then
       eps_loc = eps
    else
       eps_loc = eps_default
    end if

    
    ! estimate the cutoff eigenvalue
    nmax = 0
    sum = 0.0_dp
    do       
       fac = 1.0_dp + (lambda*pi*nmax/(x2-x1))**2
       fac = fac**(-s)
       sum = sum + fac
       if(fac/sum < eps_loc) exit
       nmax = nmax+1
    end do
    
    ! work out the mesh size
    dx = (x2-x1)/nmax
    nspec = (x2-x1)/dx
    dx = (x2-x1)/nspec
    nspec = nspec + 2*npad
    x11 = x1 - npad*dx
    x22 = x2 + npad*dx
    ispec1 = 1 + npad
    ispec2 = nspec - npad
    
    ! allocate the mesh arrays
    allocate(hp(ngll,ngll))
    allocate(jac(nspec))
    allocate(xx(ngll,nspec))
    allocate(ibool(ngll,nspec))

    
    ! get the gll points and weights  
    call quad%set(ngll)
    do inode = 1,ngll
       call lagrange_polynomial(quad%x(inode),ngll,quad%x,xx(:,1),hp(inode,:))
    end do
    
    ! build up the mesh
    xl = x11
    count = 0
    do ispec = 1,nspec
       xr = xl + dx
       jac(ispec) = 0.5_dp*(xr-xl)
       do inode = 1,ngll
          xx(inode,ispec) = xl + 0.5_dp*(xr-xl)*(quad%x(inode)+1.0_dp)
          count = count + 1
          ibool(inode,ispec) = count
       end do       
       xl = xr
       count = count-1
    end do
    
    ! allocate the matrices
    ndim = ibool(ngll,nspec)
    kda  = ngll-1
    ldab = kda+1
    kdb  = 0
    ldbb = kdb+1
    allocate(aa(ldab,ndim),bb(ldbb,ndim))
    aa = 0.0_dp
    bb = 0.0_dp
    
    ! build the matrices
    do ispec = 1,nspec
       do inode = 1,ngll
          i = ibool(inode,ispec)
          k = kdb+1
          bb(k,i) = bb(k,i) + quad%w(inode)*jac(ispec)
          do jnode = inode,ngll
             j = ibool(jnode,ispec)
             k = kda+1+i-j
             do knode = 1,ngll
                aa(k,j) = aa(k,j) + hp(knode,inode) &
                                  * hp(knode,jnode) &
                                  * quad%w(knode)   & 
                                  / jac(ispec)
             end do
          end do
       end do
    end do
    
    ! solve the eigenvalue problem
    allocate(eval(ndim))
    allocate(evec(ndim,ndim))
    allocate(work(3*ndim-2))
    call dsbgv('V','U',ndim,kda,kdb,aa,ldab,bb,ldbb,eval,evec,ndim,work,info)
    call check(info == 0,'set_gaussian_random_scalar_field_1D','problem with eigendecomposition')
    eval(1) = 0.0_dp

    
    ! work out eigenvalue cutoff
    rfun%ndim = min(nmax,ndim)
    
    ! store the eigenvectors
    i1 = ibool(1,ispec1)
    i2 = ibool(ngll,ispec2)
    rfun%mdim = i2-i1+1

    
    allocate(rfun%evec(rfun%mdim,rfun%ndim),rfun%x(rfun%mdim))
    do i = 1,rfun%ndim
       rfun%evec(:,i) = evec(i1:i2,i)       
    end do
    do ispec = ispec1,ispec2
       do inode = 1,ngll
          i = ibool(inode,ispec)-i1+1
          rfun%x(i) = xx(inode,ispec)
       end do
    end do
    
    ! set the covariance
    allocate(rfun%Q(rfun%ndim))
    rfun%Q(:) = (1.0_dp + lambda*lambda*eval(1:rfun%ndim))**(-s)
    j = rfun%mdim/2
    sum = 0.0_dp
    do i = 1,rfun%ndim
       sum = sum + rfun%Q(i)*rfun%evec(j,i)**2
    end do
    rfun%Q = sigma*sigma*rfun%Q/sum

    ! finish up
    rfun%allocated = .true.
 
    return
  end function build_GRF_1D_SEM


  !====================================================!
  !                 1D fields using FFT                !
  !====================================================!
  
  
  subroutine delete_GRF_1D_Fourier(self)
    class(GRF_1D_Fourier), intent(inout) :: self
    if(.not.self%allocated) return
    deallocate(self%Q)
    call fftw_destroy_plan(self%plan_c2r)
    self%allocated = .false.
    return
  end subroutine delete_GRF_1D_Fourier


  subroutine realise_GRF_1D_Fourier(self,fun)
    class(GRF_1D_Fourier), intent(inout) :: self
    class(interp_1D), intent(inout) :: fun

    integer(i4b) :: i,n,i1,i2
    real(dp) :: r1,r2
    real(C_DOUBLE), pointer :: out(:)
    complex(C_DOUBLE_COMPLEX), pointer :: in(:)
    type(C_PTR) :: pin,pout

    n = self%n
    i1 = self%i1
    i2 = self%i2
    pin = fftw_alloc_complex(int(n/2+1, C_SIZE_T))
    pout  = fftw_alloc_real(int(n, C_SIZE_T))
    call c_f_pointer(pin, in, [n/2+1])
    call c_f_pointer(pout,   out, [n])
    call random_seed()
    call normal_random_variable(r1)
    in(1) = r1*sqrt(self%Q(1))
    do i = 2,n/2+1
       call normal_random_variable(r1,r2)
       in(i) = (r1+ii*r2)*sqrt(self%Q(i))
    end do
    call fftw_execute_dft_c2r(self%plan_c2r,in,out)
    call fun%set(self%x(i1:i2),out(i1:i2))
    
    return
  end subroutine realise_GRF_1D_Fourier


  subroutine corr_GRF_1D_Fourier(self,y,fun)
    class(GRF_1D_Fourier), intent(in) :: self
    real(dp), intent(in) :: y
    class(interp_1D), intent(inout) :: fun


    integer(i4b) :: m,n,i1,i2
    real(dp) :: th,y1,y2
    real(C_DOUBLE), pointer :: out(:)
    complex(C_DOUBLE_COMPLEX), pointer :: in(:)
    type(C_PTR) :: pin,pout

    n = self%n
    i1 = self%i1
    i2 = self%i2
    y1 = self%x(1)
    y2 = self%x(n)
    th = twopi*(y-y1)/(y2-y1)
    pin = fftw_alloc_complex(int(n/2+1, C_SIZE_T))
    pout  = fftw_alloc_real(int(n, C_SIZE_T))
    call c_f_pointer(pin, in, [n/2+1])
    call c_f_pointer(pout,   out, [n])
    in = 0.0_dp
    do m = 1,n/2+1
       in(m) = self%Q(m)*exp(-ii*(m-1)*th)
    end do
    call fftw_execute_dft_c2r(self%plan_c2r,in,out)
    call fun%set(self%x(i1:i2),out(i1:i2))
    
    return
  end subroutine corr_GRF_1D_Fourier
  


  
  type(GRF_1D_Fourier) function build_GRF_1D_Fourier(x1,x2,lambda,s, &
                                                    sigma,eps,fftw_flag) result(rfun)
    real(dp), intent(in) :: x1,x2,lambda,s,sigma
    real(dp), intent(in), optional :: eps
    integer(C_INT), intent(in), optional :: fftw_flag
    
    integer(i4b), parameter :: npad = 10
    real(dp), parameter :: eps_default = 1.0e-5_dp

    integer(i4b) :: i,n,i1,i2,ne
    real(dp) :: dx,x11,x22,k,eps_loc,fac,sum

    real(C_DOUBLE), pointer :: out(:)
    complex(C_DOUBLE_COMPLEX), pointer :: in(:)
    type(C_PTR) :: pin,pout,plan
    integer(C_INT) :: plan_flag

    if(present(eps)) then
       eps_loc = eps
    else
       eps_loc = eps_default
    end if

    if(present(fftw_flag)) then
       plan_flag = fftw_flag
    else
       plan_flag = FFTW_MEASURE
    end if

    ! estimate the step size
    k = 0.0_dp
    sum = 0.0_dp    
    do
       fac = (1.0_dp+(lambda*k)**2)**(-s)
       sum  = sum + fac
       if(fac/sum < eps_loc) exit
       k = k+1.0_dp
    end do
    dx = 0.5_dp/k

    ! number of points and range
    n = (x2-x1)/dx + 1
    dx = (x2-x1)/(n-1)
    n = n + 2*npad
    x11 = x1-npad*dx
    x22 = x2+npad*dx
    rfun%i1 = 1+npad
    rfun%i2 = n-npad

    ! extend to next power of 2
    ne = log(1.0_dp*n)/log(2.0_dp) + 1
    n = 2**ne
    rfun%n = n

    ! set the points
    allocate(rfun%x(n))
    do i = 1,n
       rfun%x(i) = x11 + (i-1)*dx
    end do

    ! set the covariance
    allocate(rfun%Q(n/2+1))
    sum = 0.0_dp
    do i = 1,n/2+1
       k = twopi*(i-1)/(dx*n)
       fac = (1.0_dp+(lambda*k)**2)**(-s)
       if(i == 1) then
          sum = sum + fac
       else
          sum = sum + 2.0_dp*fac
       end if
       rfun%Q(i) = fac
    end do
    rfun%Q = sigma*sigma*rfun%Q/sum

    ! set up fft plan
    pin = fftw_alloc_complex(int(n/2+1, C_SIZE_T))
    pout  = fftw_alloc_real(int(n, C_SIZE_T))
    call c_f_pointer(pin,   in, [n/2+1])
    call c_f_pointer(pout, out, [n])
    plan = fftw_plan_dft_c2r_1d(n,in,out,plan_flag)
    call fftw_free(pin)
    call fftw_free(pout)
    rfun%plan_c2r = plan
    

    return
  end function build_GRF_1D_Fourier
  
  !====================================================!
  !              random fields on a sphere             !
  !====================================================!


  
  subroutine delete_gaussain_random_scalar_field_sphere(self)
    class(gaussain_random_scalar_field_sphere), intent(inout) :: self
    if(.not.self%allocated) return
    deallocate(self%mlm,self%Q_l)
    self%allocated = .false.
    return
  end subroutine delete_gaussain_random_scalar_field_sphere
  
  subroutine build_gaussain_random_scalar_field_sphere(self,grid,lambda,s)
    class(gaussain_random_scalar_field_sphere), intent(inout) :: self
    type(gauss_legendre_grid), intent(in) :: grid
    real(dp), intent(in) :: lambda,s

    integer(i4b) :: l
    
    call self%delete()
    self%lmax = grid%lmax
    ! allocate arrays
    allocate(self%mlm(grid%ncoef_r))
    allocate(self%ulm(grid%ncoef_r))
    self%mlm = 0.0_dp

    allocate(self%Q_l(0:grid%lmax))
    do l = 0,self%lmax
       self%Q_l(l) = (1.0_dp + lambda*lambda*l*(l+1))**(-s)
    end do
    
    return
  end subroutine build_gaussain_random_scalar_field_sphere
  

  subroutine realise_gaussain_random_scalar_field_sphere(self)
    class(gaussain_random_scalar_field_sphere), intent(inout) :: self

    integer(i4b) :: l,m,lmax,ncoef,ilm
    real(dp) :: r1,r2,std

    lmax = self%lmax
    ncoef = (lmax+1)*(lmax+2)/2

    call random_seed()    
    ilm = 0
    do l = 0,lmax
       std = sqrt(self%Q_l(l))
       ilm = ilm + 1
       call normal_random_variable(r1)
       self%ulm(ilm) = self%mlm(ilm) + std*r1       
       do m = 1,l
          ilm = ilm+1
          call normal_random_variable(r1,r2)
          self%ulm(ilm) = self%mlm(ilm) + std*(r1+ii*r2)
       end do
    end do

    return
  end subroutine realise_gaussain_random_scalar_field_sphere

  subroutine set_mean_coefficients_gaussain_random_scalar_field_sphere(self,grid,mlm)
    class(gaussain_random_scalar_field_sphere), intent(inout) :: self
    type(gauss_legendre_grid), intent(in) :: grid
    complex(dpc), dimension(grid%ncoef_r), intent(in) :: mlm
    self%mlm = mlm
    return
  end subroutine set_mean_coefficients_gaussain_random_scalar_field_sphere


  subroutine set_mean_values_gaussain_random_scalar_field_sphere(self,grid,m)
    class(gaussain_random_scalar_field_sphere), intent(inout) :: self
    type(gauss_legendre_grid), intent(in) :: grid
    real(dp), dimension(grid%nph,grid%nth), intent(in) :: m
    call grid%SH_trans(m,self%mlm)
    return
  end subroutine set_mean_values_gaussain_random_scalar_field_sphere


  
end module module_random_fields
