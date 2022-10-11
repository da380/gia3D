module module_random_fields

  use module_constants
  use module_util
  use module_error
  use module_quadrature
  use module_interp
  use module_spherical_harmonics
  implicit none


  type, abstract :: GRSF_interval
     real(dp) :: x1
     real(dp) :: x2
   contains
     procedure(GRSF_I_delete),  deferred :: delete
     procedure(GRSF_I_realise), deferred :: realise
  end type GRSF_interval

  abstract interface
     subroutine GRSF_I_delete(self)
       import :: GRSF_interval
       class(GRSF_interval), intent(inout) :: self
     end subroutine GRSF_I_delete
     subroutine GRSF_I_realise(self,fun)
       use module_interp
       import :: GRSF_interval
       class(GRSF_interval), intent(inout) :: self
       class(interp_1D), intent(inout) :: fun
     end subroutine GRSF_I_realise
  end interface

  type, abstract, extends(GRSF_interval) :: GRSF_interval_SEM
!     logical :: allocated = .false.
!     logical :: pad       = .true.
!     integer(i4b) :: n
!     integer(i4b) :: ngll=5
!     integer(i4b) :: nspec
!     integer(i4b) :: ispec1
!     integer(i4b) :: ispec2
!     integer(i4b) :: npad = 2
!     integer(i4b), dimension(:,:), allocatable :: ibool
!     real(dp) :: eps = 1.e-4_dp
!     real(dp), dimension(:), allocatable :: mn
!     real(dp), dimension(:), allocatable :: u
!     real(dp), dimension(:), allocatable :: Qn
!     real(dp), dimension(:,:), allocatable :: evec
!     real(dp), dimension(:), allocatable :: w
!     real(dp), dimension(:), allocatable :: jac
!     real(dp), dimension(:,:), allocatable :: hp
!     real(dp), dimension(:,:), allocatable :: x
   contains
!     procedure :: delete => delete_gaussian_random_scalar_field_interval
!     procedure :: build => build_gaussian_random_scalar_field_interval
!     procedure :: set_mean_coefficients_gaussian_random_scalar_field_interval
!     procedure :: set_mean_values_gaussian_random_scalar_field_interval
!     generic   :: set_mean => set_mean_coefficients_gaussian_random_scalar_field_interval, &
!                              set_mean_values_gaussian_random_scalar_field_interval
!     procedure :: realise => realise_gaussian_random_scalar_field_interval
  end type GRSF_interval_SEM

  
  type :: gaussian_random_scalar_field_interval
     ! main parameters
     logical :: allocated = .false.
     logical :: pad       = .true.
     integer(i4b) :: n
     real(dp) :: eps = 1.e-4_dp
     real(dp), dimension(:), allocatable :: mn
     real(dp), dimension(:), allocatable :: u
     real(dp), dimension(:), allocatable :: Qn
     real(dp), dimension(:,:), allocatable :: evec
     ! mesh parameters
     integer(i4b) :: ngll=5
     integer(i4b) :: nspec
     integer(i4b) :: ispec1
     integer(i4b) :: ispec2
     integer(i4b) :: npad = 2
     integer(i4b), dimension(:,:), allocatable :: ibool
     real(dp) :: x1,x2
     real(dp), dimension(:), allocatable :: w
     real(dp), dimension(:), allocatable :: jac
     real(dp), dimension(:,:), allocatable :: hp
     real(dp), dimension(:,:), allocatable :: x
   contains
     procedure :: delete => delete_gaussian_random_scalar_field_interval
     procedure :: build => build_gaussian_random_scalar_field_interval
     procedure :: set_mean_coefficients_gaussian_random_scalar_field_interval
     procedure :: set_mean_values_gaussian_random_scalar_field_interval
     generic   :: set_mean => set_mean_coefficients_gaussian_random_scalar_field_interval, &
                              set_mean_values_gaussian_random_scalar_field_interval
     procedure :: realise => realise_gaussian_random_scalar_field_interval
  end type gaussian_random_scalar_field_interval


   
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
  !            random fields on an interval            !
  !====================================================!
  

  subroutine delete_gaussian_random_scalar_field_interval(self)
    class(gaussian_random_scalar_field_interval), intent(inout) :: self
    if(.not.self%allocated) return
    deallocate(self%w,self%jac,self%hp,self%x)
    self%allocated = .false.
    return
  end subroutine delete_gaussian_random_scalar_field_interval


  subroutine build_gaussian_random_scalar_field_interval(self,x1,x2,lambda,s,sigma,ngll,pad,npad,eps)
    class(gaussian_random_scalar_field_interval), intent(inout) :: self
    real(dp), intent(in) :: x1,x2,lambda,s,sigma
    integer(i4b), intent(in), optional :: ngll
    logical, intent(in), optional :: pad
    integer(i4b), intent(in), optional :: npad
    real(dp), intent(in), optional :: eps
    
    integer(i4b) :: inode,ispec,count,kda,ldab,kdb,ldbb, &
                    ndim,i,j,k,jnode,knode,info,nmax
    real(dp) :: x11,x22,xl,xr,dx,fac,sum
    real(dp), dimension(:), allocatable :: eval,work
    real(dp), dimension(:,:), allocatable :: aa,bb,evec
    type(gauss_lobatto_quadrature) :: quad
    
    ! deallocate if necessary
    call self%delete()
    
    ! deal with optional arguments
    if(present(ngll)) self%ngll = ngll
    if(present(pad)) self%pad = pad
    if(present(npad)) self%npad = npad
    if(present(npad)) self%eps = eps


    ! estimate the cutoff eigenvalue
    nmax = 0
    sum = 0.0_dp
    do       
       fac = 1.0_dp + (lambda*pi*nmax/(x2-x1))**2
       fac = fac**(-s)
       sum = sum + fac
       if(fac/sum < self%eps) exit
       nmax = nmax+1
    end do
    
    ! work out the mesh size
    dx = 0.5_dp*(x2-x1)/nmax
    self%nspec = (x2-x1)/dx
    dx = (x2-x1)/self%nspec
    
    if(self%pad) then
       self%nspec = self%nspec + 2*self%npad
       x11 = x1 - self%npad*dx
       x22 = x2 + self%npad*dx
       self%ispec1 = 1 + self%npad
       self%ispec2 = self%nspec - self%npad
    else
       x11 = x1
       x22 = x2
       self%ispec1 = 1
       self%ispec2 = self%nspec
    end if
        
    ! allocate the mesh arrays
    allocate(self%w(self%ngll))
    allocate(self%hp(self%ngll,self%ngll))
    allocate(self%jac(self%nspec))
    allocate(self%x(self%ngll,self%nspec))
    allocate(self%ibool(self%ngll,self%nspec))
    
    ! get the gll points and weights  
    call quad%set(self%ngll)
    self%w = quad%w
    do inode = 1,self%ngll
       call lagrange_polynomial(quad%x(inode),self%ngll,quad%x,self%x(:,1),self%hp(inode,:))
    end do

    ! build up the mesh
    xl = x11
    count = 0
    do ispec = 1,self%nspec
       xr = xl + dx
       self%jac(ispec) = 0.5_dp*(xr-xl)
       do inode = 1,self%ngll
          self%x(inode,ispec) = xl + 0.5_dp*(xr-xl)*(quad%x(inode)+1.0_dp)
          count = count + 1
          self%ibool(inode,ispec) = count
       end do       
       xl = xr
       count = count-1
    end do
    
    ! allocate the matrices
    ndim = self%ibool(self%ngll,self%nspec)    
    kda  = self%ngll-1
    ldab = kda+1
    kdb  = 0
    ldbb = kdb+1
    allocate(aa(ldab,ndim),bb(ldbb,ndim))
    aa = 0.0_dp
    bb = 0.0_dp

    ! build the matrices
    do ispec = 1,self%nspec
       do inode = 1,self%ngll
          i = self%ibool(inode,ispec)
          k = kdb+1
          bb(k,i) = bb(k,i) + self%w(inode)*self%jac(ispec)
          do jnode = inode,self%ngll
             j = self%ibool(jnode,ispec)
             k = kda+1+i-j
             do knode = 1,self%ngll
                aa(k,j) = aa(k,j) + self%hp(knode,inode) &
                                  * self%hp(knode,jnode) &
                                  * self%w(knode)        & 
                                  / self%jac(ispec)
             end do
          end do
       end do
    end do
    
    ! solve the eigenvalue problem
    allocate(eval(ndim))
    allocate(evec(ndim,ndim))
    allocate(work(3*ndim-2))
    call dsbgv('V','U',ndim,kda,kdb,aa,ldab,bb,ldbb,eval,evec,ndim,work,info)
    call check(info == 0,'set_gaussian_random_scalar_field_interval','problem with eigendecomposition')
    eval(1) = 0.0_dp
    
    ! work out eigenvalue cutoff
    self%n = min(nmax,ndim)

    
    ! store the eigenvectors
    allocate(self%evec(ndim,self%n))
    self%evec(:,:) = evec(:,1:self%n)

    ! set the covariance
    allocate(self%Qn(self%n))
    self%Qn(:) = (1.0_dp + lambda*lambda*eval(1:self%n))**(-s)
    sum = 0.0_dp
    j = ndim/2
    do i = 1,self%n
       sum = sum + self%Qn(i)*self%evec(j,i)**2
    end do
    self%Qn = sigma*sigma*self%Qn/sum
    
    ! allocate array for the mean coefficients and realised function
    allocate(self%mn(self%n),self%u(ndim))
    self%mn = 0.0_dp
    self%u  = 0.0_dp
    
    ! finish up
    self%allocated = .true.
    
    return
  end subroutine build_gaussian_random_scalar_field_interval


  

  
  subroutine set_mean_coefficients_gaussian_random_scalar_field_interval(self,mn)
    class(gaussian_random_scalar_field_interval), intent(inout) :: self
    real(dp), dimension(self%n), intent(in) :: mn
    self%mn = mn
    return
  end subroutine set_mean_coefficients_gaussian_random_scalar_field_interval


  subroutine set_mean_values_gaussian_random_scalar_field_interval(self,mfun)
    class(gaussian_random_scalar_field_interval), intent(inout) :: self
    interface r2r
       real(dp) function r2r(x) result(f)
         use module_constants
         real(dp), intent(in) :: x
       end function r2r
    end interface r2r
    procedure(r2r) :: mfun
    integer(i4b) :: i,inode,ispec,j
    real(dp) :: x,sum
    do i = 1,self%n
       sum = 0.0_dp
       do ispec = 1,self%nspec
          do inode = 1,self%ngll
             j = self%ibool(inode,ispec)
             x = self%x(inode,ispec)
             sum = sum + mfun(x)*self%evec(j,i)*self%w(inode)*self%jac(ispec)
          end do
       end do
       self%mn(i) = sum
    end do
    return
  end subroutine set_mean_values_gaussian_random_scalar_field_interval
  

  subroutine realise_gaussian_random_scalar_field_interval(self)
    class(gaussian_random_scalar_field_interval), intent(inout) :: self

    integer(i4b) :: i
    real(dp), dimension(self%n) :: un
    call random_seed()
    call normal_random_variable(un)
    un = self%mn + un*sqrt(self%Qn)
    self%u = 0.0_dp
    do i = 1,self%n
       self%u(:) = self%u(:) + un(i)*self%evec(:,i)
    end do
    return
  end subroutine realise_gaussian_random_scalar_field_interval



  
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
