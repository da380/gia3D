module module_spherical_harmonics
  
  use module_constants
  use module_error

  type legendre_value
     private
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
     procedure :: deg => degree_legendre_value
     procedure :: get_single_legendre_value
     procedure :: get_slice_legendre_value
     generic   :: get => get_single_legendre_value, &
                         get_slice_legendre_value     
  end type legendre_value
  private :: delete_legendre_value,     &
             allocate_legendre_value,   &
             initialise_legendre_value, &
             next_legendre_value,       &
             degree_legendre_value,     &
             get_single_legendre_value, &
             get_slice_legendre_value



  type wigner_value
     private
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
     procedure :: ind => index_wigner_value
     procedure :: init => initialise_wigner_value
     procedure :: next => next_wigner_value
     procedure :: deg => degree_wigner_value
     procedure :: get_single_wigner_value
     procedure :: get_slice_wigner_value
     generic   :: get => get_single_wigner_value, &
                         get_slice_wigner_value
  end type wigner_value
  private :: delete_wigner_value,     &
             allocate_wigner_value,   &
             initialise_wigner_value, &
             next_wigner_value,       &
             degree_wigner_value,     &
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
    p%l  = 0
    p%th = th
    p%mmax = mmax
    call p%allocate(mmax)       
    p%v(0) = sifourpi
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

  function degree_legendre_value(p) result(l)
    implicit none
    class(legendre_value), intent(in) :: p
    integer(i4b) :: l
    l = p%l
    return
  end function degree_legendre_value
  
  function get_single_legendre_value(p,m,flip) result(v)
    implicit none
    class(legendre_value), intent(in) :: p
    integer(i4b), intent(in) :: m
    logical, intent(in), optional :: flip
    real(dp) :: v
    call error(.not.p%allocated,'get_single_legendre_value','not allocated')
    call error(abs(m) > p%mmax, 'get_single_legendre_value','m out of range')  
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
    call error(.not.p%allocated,'get_single_legendre_value','not allocated')
    call error(abs(m1) > p%mmax, 'get_single_legendre_value','m out of range')
    call error(abs(m2) > p%mmax, 'get_single_legendre_value','m out of range')
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

  
  subroutine legendre(l,mmax,th,p)
    implicit none
    integer(i4b), intent(in) :: l
    integer(i4b), intent(in) :: mmax
    real(dp), intent(in)     :: th
    real(dp), dimension(0:mmax), intent(out) :: p

    integer(i4b) :: n,m
    real(dp) :: x,y,fac1,fac2
    real(dp), dimension(0:mmax) :: pm1,pm2


    call error( (mmax > l) .or. (mmax < 0),  &
               'legendre','invalid mmax value')
   
    ! set some trig functions
    x = cos(th)
    y = sqrt(1.0_dp-x*x)
    
    ! initialise
    p   = 0.0_dp    
    pm1 = 0.0_dp
    pm2 = 0.0_dp
    
    ! initilaise l = 0
    p(0) = sifourpi
    
    ! are we done? 
    if(l == 0) return
    
    ! save previous degree
    pm1(0) = p(0)
    
    ! initialise l = 1
    p(0) =  x*sqrt3*p(0)
    if(mmax > 0) then
       fac1 = sqrt3/sqrt2
       p(1) = -fac1*y*pm1(0)
    end if
    
    ! are we done? 
    if(l == 1) return
    
    ! save previous degrees
    pm2(0)   = pm1(0)
    pm1(0:min(mmax,1)) = p(0:min(mmax,1))

        
    do n = 2,l

       ! apply the three-term recursion for m < min(mmax,n-2)
       do m = 0,min(mmax,n-2)
          fac1 = (4.0_dp*n*n-1.0_dp)/(n*n-m*m)
          fac1 = sqrt(fac1)
          fac2 = ( (n-1.0_dp)*(n-1.0_dp)-m*m)            &
                  /(4.0_dp*(n-1.0_dp)*(n-1.0_dp)-1.0_dp)
          fac2 = sqrt(fac2)
          p(m) = fac1*(x*pm1(m) - fac2*pm2(m))
       end do

       ! apply two-term recursion P_{n-1,n-1} -> P_{n,n-1}
       if(n-1 <= mmax) then
          fac1 = 2.0_dp*n+1
          fac1 = sqrt(fac1)
          p(n-1) = x*fac1*pm1(n-1)
       end if

       ! apply two-term recursion P_{n-1,n-1} -> P_{n,n}
       if(n <= mmax) then
          fac1 = (n+0.5_dp)/n
          fac1 = sqrt(fac1)
          p(n) = -fac1*y*pm1(n-1)
       end if

       ! update the stored values
       pm2(0:min(mmax,n-1)) = pm1(0:min(mmax,n-1))
       pm1(0:min(mmax,n))   =   p(0:min(mmax,n))
       
    end do
    
    return
  end subroutine legendre



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
    p%l = 0
    p%beta = beta
    call p%allocate(nmax,mmax)
    p%v(1) = 1.0_dp
    return
  end subroutine initialise_wigner_value

  
  subroutine next_wigner_value(p)
    implicit none
    class(wigner_value), intent(inout) :: p
    integer(i4b) :: l,m,n,mmax,nmax,i,ip1,j,nm
    real(dp) :: beta,cb,chb,shb,fac1,fac2,fac3,xl,xm,xn,tlm1,lm1
    ! set local paramters
    mmax = p%mmax
    nmax = p%nmax
    l = p%l +1
    lm1 = l-1
    tlm1 = 2*l-1
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
             fac1 = tlm1*(-n*m+lm1*l*cb)
             fac2 = (lm1-n)*(lm1+n)*(lm1-m)*(lm1+m)
             fac2 = -l*sqrt(fac2)
             fac3 = (-n+l)*(n+l)*(-m+l)*(m+l)
             fac3 = lm1*sqrt(fac3)
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

  function degree_wigner_value(p) result(l)
    implicit none
    class(wigner_value), intent(in) :: p
    integer(i4b) :: l
    l = p%l
    return
  end function degree_wigner_value



  
  function get_single_wigner_value(p,nin,min,norm,flip) result(v)
    implicit none
    class(wigner_value), intent(in) :: p
    integer(i4b), intent(in) :: nin,min
    logical, intent(in), optional  :: norm
    logical, intent(in), optional  :: flip
    real(dp) :: v
    logical ::  check
    integer(i4b) :: l,n,m,i,sign,mmax,nmax,j
    real(dp) :: fac
    call error(.not.p%allocated,'get_single_wigner_value','not allocated')
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
    call error(abs(m) > mmax, 'get_single_wigner_value','m out of range')
    call error(abs(n) > nmax, 'get_single_wigner_value','n out of range')  
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
    i = p%ind(n,m)
    fac = 1.0_dp
    if(present(norm)) then
       if(norm) then
          fac = sqrt((2*l+1)/fourpi)
       end if       
    end if
    v = sign*fac*p%v(i)
    return
  end function get_single_wigner_value

  
  function get_slice_wigner_value(p,n1,n2,m1,m2,norm,flip) result(v)
    implicit none
    class(wigner_value), intent(in) :: p
    integer(i4b), intent(in) :: n1,n2,m1,m2
    real(dp), dimension(n2-n1+1,m2-m1+1) :: v
    logical, intent(in), optional :: norm
    logical, intent(in), optional :: flip
    integer(i4b) :: n,m,i,j
    i = 0
    do n = n1,n2
       i = i+1
       j = 0
       do m = m1,m2
          j = j+1
          v(i,j) = p%get(n,m,norm,flip)          
       end do
    end do
    return
  end function get_slice_wigner_value

  


  
end module module_spherical_harmonics
