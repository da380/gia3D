module module_random_fields

  use module_constants
  use module_util
  use module_error
  use module_quadrature
  use module_spherical_harmonics
  implicit none


  type gaussian_random_scalar_field_interval
     ! basic parameters
     logical :: allocated = .false.
     logical :: pad       = .true.
     integer(i4b) :: n
     real(dp) :: eps = 1.e-6_dp
     real(dp), dimension(:), allocatable :: mn
     real(dp), dimension(:), allocatable :: un
     real(dp), dimension(:), allocatable :: eval
     real(dp), dimension(:,:), allocatable :: evec
     ! mesh parameters
     integer(i4b) :: ngll=5
     integer(i4b) :: nspec
     integer(i4b) :: npad = 2
     integer(i4b), dimension(:,:), allocatable :: ibool
     real(dp) :: x1,x2
     real(dp), dimension(:), allocatable :: w
     real(dp), dimension(:), allocatable :: jac
     real(dp), dimension(:,:), allocatable :: hp
     real(dp), dimension(:,:), allocatable :: x
   contains
     procedure :: delete => delete_gaussian_random_scalar_field_interval
     procedure :: set => set_gaussian_random_scalar_field_interval
     procedure :: realise => realise_gaussian_random_scalar_field_interval
  end type gaussian_random_scalar_field_interval


   
  type gaussain_random_scalar_field_sphere
     logical :: allocated = .false.
     logical :: real=.false.
     integer(i4b) :: lmax
     real(dp), dimension(:), allocatable :: Q_l
     complex(dpc), dimension(:), allocatable :: mlm
     complex(dpc), dimension(:), allocatable :: ulm
   contains
     procedure :: delete => delete_gaussain_random_scalar_field_sphere
     procedure :: set => set_gaussain_random_scalar_field_sphere
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


  subroutine set_gaussian_random_scalar_field_interval(self,x1,x2,lambda,s,ngll,pad,npad)
    class(gaussian_random_scalar_field_interval), intent(inout) :: self
    real(dp), intent(in) :: x1,x2,lambda,s
    integer(i4b), intent(in), optional :: ngll
    logical, intent(in), optional :: pad
    integer(i4b), intent(in), optional :: npad
    
    integer(i4b) :: inode,ispec,count,kda,ldab,kdb,ldbb, &
                    ndim,i,j,k,jnode,knode,info
    real(dp) :: x11,x22,xl,xr,dx,fac
    real(dp), dimension(:), allocatable :: eval,work
    real(dp), dimension(:,:), allocatable :: aa,bb,evec
    type(gauss_lobatto_quadrature) :: quad
    

    ! deallocate if necessary
    call self%delete()
    
    ! deal with optional arguments
    if(present(ngll)) self%ngll = ngll
    if(present(pad)) self%pad = pad
    if(present(npad)) self%npad = npad

    ! work out the mesh size
    self%nspec = 6*(x2-x1)/lambda
    dx = (x2-x1)/self%nspec

    if(self%pad) then
       self%nspec = self%nspec + 2*self%npad
       x11 = x1 - self%npad*dx
       x22 = x2 + self%npad*dx
    else
       x11 = x1
       x22 = x2
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
          self%x(inode,ispec) = xl + 0.5_dp*(xr-xl)*quad%x(inode)
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

    ! work out eigenvalue cutoff
    self%n = 0
    do i = 1,ndim
       fac = 1.0_dp+lambda*lambda*eval(i)
       fac = fac**(-s)
       if(fac < self%eps) then
          self%n = i
          exit
       end if
    end do

    ! store the eigenvalues and vectors
    allocate(self%eval(self%n))
    allocate(self%evec(ndim,self%n))
    self%eval(:) = eval(1:self%n)
    self%evec(:,:) = evec(:,1:self%n)

    ! allocate array for the mean and realised coefficients
    allocate(self%mn(self%n),self%un(self%n))
    self%mn = 0.0_dp
    self%un = 0.0_dp
    
    ! finish up
    self%allocated = .true.
    
    return
  end subroutine set_gaussian_random_scalar_field_interval
  

  subroutine realise_gaussian_random_scalar_field_interval(self)
    class(gaussian_random_scalar_field_interval), intent(inout) :: self
    call normal_random_variable(self%un)
    self%un = self%mn + sqrt(self%eval)*self%un
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
  
  subroutine set_gaussain_random_scalar_field_sphere(self,grid,lambda,s,real)
    class(gaussain_random_scalar_field_sphere), intent(inout) :: self
    type(gauss_legendre_grid), intent(in) :: grid
    real(dp), intent(in) :: lambda,s
    logical, intent(in), optional :: real

    integer(i4b) :: l
    
    call self%delete()
    self%lmax = grid%lmax
    if(present(real)) self%real = real
    
    if(self%real) then
       allocate(self%mlm(grid%ncoef_r))
       allocate(self%ulm(grid%ncoef_r))
    else
       allocate(self%mlm(grid%ncoef_c))
       allocate(self%ulm(grid%ncoef_c))
    end if
    self%mlm = 0.0_dp

    allocate(self%Q_l(0:grid%lmax))
    do l = 0,self%lmax
       self%Q_l(l) = (1.0_dp + lambda*lambda*l*(l+1))**(-s)
    end do
    
    return
  end subroutine set_gaussain_random_scalar_field_sphere
  

  subroutine realise_gaussain_random_scalar_field_sphere(self)
    class(gaussain_random_scalar_field_sphere), intent(inout) :: self

    integer(i4b) :: l,m,lmax,ncoef,ilm
    real(dp) :: r1,r2,std

    lmax = self%lmax
    if(self%real) then
       ncoef = (lmax+1)*(lmax+2)/2
    else
       ncoef = (lmax+1)**2
    end if

    call random_seed()    
    ilm = 0
    do l = 0,lmax
       std = sqrt(self%Q_l(l))
       ilm = ilm + 1
       call normal_random_variable(r1,r2)
       if(self%real) then
          self%ulm(ilm) = self%mlm(ilm) + std*r1
       else
          self%ulm(ilm) = self%mlm(ilm) + std*(r1+ii*r2)
       end if
       do m = 1,l
          ilm = ilm+1
          call normal_random_variable(r1,r2)
          self%ulm(ilm) = self%mlm(ilm) + std*(r1+ii*r2)
          if(.not.self%real) then
             call normal_random_variable(r1,r2)
             ilm = ilm+1
             self%ulm(ilm) = self%mlm(ilm) + std*(r1+ii*r2)
          end if
       end do
    end do

    return
  end subroutine realise_gaussain_random_scalar_field_sphere


!  subroutine realise_gaussain_random_scalar_field_sphere_spatial_real(self,grid,u)
!    class(gaussain_random_scalar_field_sphere), intent(in) :: self
!    type(gauss_legendre_grid), intent(in) :: grid
!    real(dp), dimension(grid%nph,grid%nth), intent(out) :: u
!    complex(dpc), dimension(grid%ncoef_r) :: ulm
!    call check(self%real,'realise_gaussain_random_scalar_field_sphere_spatial_real', &
!                         'the random field is complex!')
!    call self%realise_gaussain_random_scalar_field_sphere(ulm)
!    call grid%SH_itrans(ulm,u)
!    return
!  end subroutine realise_gaussain_random_scalar_field_sphere_spatial_real

!  subroutine realise_gaussain_random_scalar_field_sphere_spatial_complex(self,grid,u)
!    class(gaussain_random_scalar_field_sphere), intent(in) :: self
!    type(gauss_legendre_grid), intent(in) :: grid
!    complex(dpc), dimension(grid%nph,grid%nth), intent(out) :: u
!    complex(dpc), dimension(grid%ncoef_c) :: ulm
!    call check(.not.self%real,'realise_gaussain_random_scalar_field_sphere_spatial_complex', &
!                              'the random field is real!')
!    call self%realise_gaussain_random_scalar_field_sphere(ulm)
!    call grid%SH_itrans(ulm,u)
!    return
!  end subroutine realise_gaussain_random_scalar_field_sphere_spatial_complex

  
end module module_random_fields
