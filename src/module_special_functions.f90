module module_special_functions  
  use module_constants
  implicit none
 



  type legendre_value
     logical  :: allocated = .false.
     integer(i4b) :: l
     integer(i4b) :: mmax
     real(dp) :: th
     real(dp), dimension(:), allocatable :: vm1
     real(dp), dimension(:), allocatable :: v
     real(dp), dimension(:), allocatable :: vp1
   contains
     procedure :: delete => delete_legendre_value
     procedure :: allocate => allocate_legendre_value
     procedure :: init => initialise_legendre_value
     procedure :: next => next_legendre_value
     procedure :: get_single_legendre_value
     procedure :: get_slice_legendre_value
     generic   :: get => get_single_legendre_value, &
                         get_slice_legendre_value     
  end type legendre_value
  private :: delete_legendre_value,     &
             allocate_legendre_value,   &
             initialise_legendre_value, &
             next_legendre_value,       &
             get_single_legendre_value, &
             get_slice_legendre_value



  type wigner_value
     logical :: allocated = .false.
     logical :: transposed = .false.
     integer(i4b) :: l
     integer(i4b) :: nmax
     integer(i4b) :: mmax
     real(dp) :: beta
     real(dp), dimension(:), allocatable :: vm1
     real(dp), dimension(:), allocatable :: v
     real(dp), dimension(:), allocatable :: vp1
   contains
     procedure :: delete => delete_wigner_value
     procedure :: allocate => allocate_wigner_value
     procedure :: index => index_wigner_value
     procedure :: init => initialise_wigner_value
     procedure :: next => next_wigner_value
     procedure :: get_single_wigner_value
     procedure :: get_slice_wigner_value
     generic   :: get => get_single_wigner_value, &
                         get_slice_wigner_value
  end type wigner_value
  private :: delete_wigner_value,     &
             allocate_wigner_value,   &
             initialise_wigner_value, &
             next_wigner_value,       &
             get_single_wigner_value



  

contains

  
  !============================================================================!
  !============================================================================!
  !                     procedures for legendre polynomials                    !
  !============================================================================!
  !============================================================================!
  
  subroutine delete_legendre_value(p)
    implicit none
    class(legendre_value), intent(inout) :: p
    if(.not.p%allocated) return
    deallocate(p%vm1)
    deallocate(p%v)
    deallocate(p%vp1)
    p%allocated = .false.
    return
  end subroutine delete_legendre_value

  subroutine allocate_legendre_value(p,mmax)
    implicit none
    class(legendre_value), intent(inout) :: p
    integer(i4b), intent(in) :: mmax
    if(p%allocated) call p%delete()
    allocate(p%vm1(0:mmax))
    allocate(p%v(0:mmax))
    allocate(p%vp1(0:mmax))
    p%vm1 = 0.0_dp
    p%v   = 0.0_dp
    p%vp1 = 0.0_dp
    p%allocated = .true.
    return
  end subroutine allocate_legendre_value

  subroutine initialise_legendre_value(p,th,mmax)
    implicit none
    class(legendre_value), intent(inout) :: p
    real(dp), intent(in) :: th
    integer(i4b), intent(in) :: mmax
    call p%delete()    
    p%l  = -1
    p%th = th
    p%mmax = mmax
    call p%allocate(mmax)       
    return
  end subroutine initialise_legendre_value


  subroutine next_legendre_value(p)
    implicit none
    class(legendre_value), intent(inout) :: p
    integer(i4b) :: l,lold,m,mmax,im
    real(dp) :: th,ct,st,fac1,fac2
    ! get local parameters
    lold = p%l
    l = lold+1
    if(l == 0) then
       p%v(0) = sifourpi
       p%l = l
       return
    end if
    th = p%th
    ct = cos(th)
    st = sin(th)
    mmax = p%mmax
    ! apply three-term recursion
    do m = 0,min(mmax,l-1)
       fac1 = (4.0_dp*l*l-1.0_dp)/(l*l-m*m)
       fac1 = sqrt(fac1)
       fac2 = ( (l-1.0_dp)*(l-1.0_dp)-m*m)            &
               /(4.0_dp*(l-1.0_dp)*(l-1.0_dp)-1.0_dp)
       fac2 = sqrt(fac2)
       p%vp1(m) = fac1*(ct*p%v(m) - fac2*p%vm1(m))
    end do
    ! two term recursion for m = l
    if(l <= mmax) then
       fac1 = (l+0.5_dp)/l
       fac1 = sqrt(fac1)
       p%vp1(l) = -fac1*st*p%v(l-1)
    end if
    ! update the data arrays
    p%l = p%l+1
    im = min(lold,mmax)
    p%vm1(0:im) = p%v(0:im)
    im = min(l,mmax)
    p%v(0:im) = p%vp1(0:im)
    return
  end subroutine next_legendre_value

  
  function get_single_legendre_value(p,m,flip) result(v)
    implicit none
    class(legendre_value), intent(in) :: p
    integer(i4b), intent(in) :: m
    logical, intent(in), optional :: flip
    real(dp) :: v
    v = p%v(abs(m))
    if(m < 0 .and. modulo(m,2) /= 0) v = -v
    if(present(flip)) then
       if(flip .and. modulo(p%l+m,2) /= 0) then
          v = -v          
       end if
    end if   
    return
  end function get_single_legendre_value

    
  function get_slice_legendre_value(p,m1,m2,flip) result(v)
    implicit none
    class(legendre_value), intent(in) :: p
    integer(i4b), intent(in) :: m1,m2
    real(dp), dimension(m2-m1+1) :: v
    logical, intent(in), optional :: flip
    integer(i4b) :: m,im
    im = 0
    do m = m1,m2
       im = im+1
       v(im) = p%v(abs(m))
       if(m < 0 .and. modulo(m,2) /= 0) v(im) = -v(im)
       if(present(flip)) then
          if(flip .and. modulo(p%l+m,2) /= 0) then
             v(im) = -v(im)
          end if
       end if      
    end do   
    return
  end function get_slice_legendre_value

  



  !============================================================================!
  !============================================================================!
  !              procedures for generalised legendre polynomials               !
  !============================================================================!
  !============================================================================!


  subroutine delete_wigner_value(p)
    implicit none
    class(wigner_value), intent(inout) :: p
    if(.not.p%allocated) return
    deallocate(p%vm1)
    deallocate(p%v)
    deallocate(p%vp1)
    p%allocated = .false.
    return
  end subroutine delete_wigner_value
  
  
  subroutine allocate_wigner_value(p,nmax,mmax)
    implicit none
    class(wigner_value), intent(inout) :: p
    integer(i4b), intent(in) :: nmax,mmax
    integer(i4b) :: pdim,n,m
    if(p%allocated) call p%delete()
    if(mmax >= nmax) then
       p%mmax = mmax
       p%nmax = nmax
    else
       p%transposed = .true.
       p%mmax = nmax
       p%nmax = mmax
    end if
    pdim = (p%nmax+1)**2 + (p%mmax-p%nmax)*(2*p%nmax+1)
    allocate(p%vm1(pdim))
    allocate(p%v(pdim))
    allocate(p%vp1(pdim))
    p%vm1 = 0.0_dp
    p%v   = 0.0_dp
    p%vp1 = 0.0_dp
    p%allocated = .true.
    return
  end subroutine allocate_wigner_value

  function index_wigner_value(p,n,m) result(i)
    implicit none
    class(wigner_value), intent(in) :: p
    integer(i4b), intent(in) :: m,n
    integer(i4b) :: i,mmax
    if(m <= p%nmax) then
       i = m**2 + n+m+1
    else
       i = (p%nmax+1)**2 + (m-p%nmax-1)*(2*p%nmax+1) + n+p%nmax+1
    end if     
    return
  end function index_wigner_value

  
  subroutine initialise_wigner_value(p,beta,nmax,mmax)
    implicit none
    class(wigner_value), intent(inout) :: p
    real(dp), intent(in) :: beta
    integer(i4b), intent(in) :: mmax,nmax
    integer(i4b) :: i
    p%l = -1
    p%beta = beta
    call p%allocate(nmax,mmax)
    return
  end subroutine initialise_wigner_value

  
  subroutine next_wigner_value(p)
    implicit none
    class(wigner_value), intent(inout) :: p
    integer(i4b) :: l,m,n,mmax,nmax,i,ip1,j,nm
    real(dp) :: beta,cb,chb,shb,fac1,fac2,fac3,xl,xm,xn
    ! set local paramters
    mmax = p%mmax
    nmax = p%nmax
    l = p%l +1
    if(l == 0) then
       p%v(1) = 1.0_dp
       p%l = l
       return
    end if
    xl = l
    beta = p%beta
    cb = cos(beta)    
    i = 0    
    if(l == 1) then
       ! deal with l = 1 separately
       i = i+1
       p%vp1(i) = cb       
    else
       ! apply three-term recursion
       do m = 0,min(mmax,l-1)
          xm = m
          nm = min(m,nmax)
          do n = -nm,nm
             xn = n                 
             i = i+1
             fac1 = (2*xl-1)*(-xn*xm+(xl-1)*xl*cb)
             fac2 = (xl-1-xn)*(xl-1+xn)*(xl-1-xm)*(xl-1+xm)
             fac2 = -xl*sqrt(fac2)
             fac3 = (-xn+xl)*(xn+xl)*(-xm+xl)*(xm+xl)
             fac3 = (xl-1)*sqrt(fac3)
             fac1 = fac1/fac3
             fac2 = fac2/fac3
             ! update values
             p%vp1(i) = fac1*p%v(i) + fac2*p%vm1(i)
          end do
       end do
    end if
    ! set edge values directly
    if(l <= mmax) then      
       chb = cos(0.5_dp*beta)
       shb = sin(0.5_dp*beta)
       if(shb <= 0.0_dp) then
          ! set special values at beta = 0.0
          do n = -min(l,nmax),min(l-1,nmax)
             i = i+1
             p%vp1(i) = 0.0_dp
          end do
          if(l <= nmax) then
             i = i+1
             p%vp1(i) = 1.0_dp
          end if                    
       else if(chb <= 0.0_dp) then         
          ! set special value at beta == pi
          if(l <= nmax) then
             i = i+1
             p%vp1(i) = 1.0_dp
          end if
          do n = -min(l-1,nmax),min(l,nmax)
             i = i+1
             p%vp1(i) = 0.0_dp
          end do                    
       else
          ! deal with interior arguments
          chb = log(chb)
          shb = log(shb)          
          nm = min(l,nmax)
          do n = -nm,nm
             xn = n
             i = i+1
             fac1 =  0.5_dp*(log_gamma(2*xl+1)-log_gamma(xl-xn+1)  &
                           - log_gamma(xl+xn+1)) + (l+n)*chb + (l-n)*shb
             fac1 = exp(fac1)
             if(modulo(l-n,2) /= 0) fac1 = -fac1
             p%vp1(i) = fac1
          end do
       end if
    end if
    p%l = l
    p%vm1 = p%v
    p%v = p%vp1
    return
  end subroutine next_wigner_value

  
  function get_single_wigner_value(p,nin,min,flip) result(v)
    implicit none
    class(wigner_value), intent(in) :: p
    integer(i4b), intent(in) :: nin,min
    logical, intent(in), optional  :: flip
    real(dp) :: v
    logical ::  check
    integer(i4b) :: l,n,m,i,sign,mmax,nmax,j
    real(dp) :: fac
    l = p%l    
    ! transpose the indices if needed
    if(p%transposed) then
       n = min
       m = nin
    else
       n = nin       
       m = min
    end if
    mmax = p%mmax
    nmax = p%nmax
    sign = 1
    if(present(flip)) then
       if(flip) then
          ! map beta to pi-beta
          if(modulo(l-n,2) /= 0) sign = -1
          j = n
          n = m
          m = -j
       end if
    end if
    if(m < 0 .and. abs(n) <= abs(m)) then
       ! parity transform
       if(modulo(n-m,2) /= 0) sign = -sign
       m = -m
       n = -n
    else if(n > abs(m)) then
       ! reflection in diagonal
       if(modulo(n-m,2) /= 0) sign = -sign
       j = n
       n = m
       m = j
    else if(n < -abs(m)) then
       ! parity and reflection
       j = n
       n = -m
       m = -j
    end if
    if(p%transposed) then
       ! fix sign if the transposed matrix is stored
       if(modulo(n-m,2) /= 0) sign = -sign
    end if
    i = p%index(n,m)
    v = sign*p%v(i)
    return
  end function get_single_wigner_value

  
  function get_slice_wigner_value(p,n1,n2,m1,m2,flip) result(v)
    implicit none
    class(wigner_value), intent(in) :: p
    integer(i4b), intent(in) :: n1,n2,m1,m2
    real(dp), dimension(n2-n1+1,m2-m1+1) :: v
    logical, intent(in), optional :: flip
    integer(i4b) :: n,m,i,j
    i = 0
    do n = n1,n2
       i = i+1
       j = 0
       do m = m1,m2
          j = j+1
          v(i,j) = p%get(n,m,flip)          
       end do
    end do
    return
  end function get_slice_wigner_value

  
  subroutine set_wigner_array(beta,lmax,nmax,dlm,norm)
    implicit none
    real(dp), intent(in) :: beta
    integer(i4b), intent(in) :: nmax,lmax
    real(dp), dimension((lmax+1)*(lmax+2)/2,-nmax:nmax), intent(out) :: dlm
    logical, intent(in), optional :: norm
    integer(i4b) :: l,i,n,m
    real(dp) :: fac
    type(wigner_value) :: p    
    call p%init(beta,nmax,lmax)
    do l = 0,lmax
       call p%next()
       fac = 1.0_dp
       if(present(norm)) then
          if(norm) fac = sqrt((2*l+1)/fourpi)
       end if
       do n = -min(nmax,l),min(nmax,l)
          do m = 0,l
             i = l*(l+1)/2 + m + 1
             dlm(i,n) = fac*p%get(n,m)
          end do
       end do
    end do
    call p%delete()    
    return
  end subroutine set_wigner_array


  !==================================================================!
  !                    Lagrange polynomial routines                  !
  !==================================================================!

  subroutine lagrange_polynomial(x,n,xn,h,hp)
    real(dp), intent(in) :: x
    integer(i4b), intent(in) :: n
    real(dp), dimension(n), intent(in) :: xn
    real(dp), dimension(n), intent(out) :: h
    real(dp), dimension(n), intent(out) :: hp

    integer(i4b) :: i,j,k
    real(dp) :: prod1,prod2


    do i = 1,n
       prod1 = 1.0_dp
       prod2 = 1.0_dp
       do j = 1,n
          if(j /= i) then
             prod1 = prod1*(x-xn(j))
             prod2 = prod2*(xn(i)-xn(j))
          endif
       end do
       h(i) = prod1/prod2

       hp(i) = 0.0_dp
       do j = 1,n
          if(j /= i) then
             prod1 = 1.0_dp
             do k = 1,n
                if(k /= i .and. k/= j) then
                   prod1 = prod1*(x-xn(k))
                end if
             end do
             hp(i) = hp(i) + prod1
          end if
       end do
       hp(i) = hp(i)/prod2
    end do


    
    return
  end subroutine lagrange_polynomial
  
  
end module module_special_functions
