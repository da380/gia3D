module module_spherical_harmonics

  use module_constants
  use module_error

  
  

  type legendre_value
!     private
     logical  :: allocated = .false.
     integer(i4b) :: l
     integer(i4b) :: mmax
     real(dp) :: th,x,y
     real(dp), dimension(:), allocatable :: vm1
     real(dp), dimension(:), allocatable :: v
     real(dp), dimension(:), allocatable :: vp1
   contains
     procedure :: delete => delete_legendre_value
     procedure :: allocate => allocate_legendre_value
     procedure :: initialise_zero_legendre_value
     procedure :: initialise_legendre_value
     generic   :: init => initialise_zero_legendre_value, & 
                          initialise_legendre_value
     procedure :: next => next_legendre_value
     !     procedure :: get_single_legendre_value
     !     procedure :: get_slice_legendre_value     
     !     generic   :: get => get_single_legendre_value, &
     !                         get_slice_legendre_value
     
  end type legendre_value

  
contains

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
    p%allocated = .true.
    return
  end subroutine allocate_legendre_value

  subroutine initialise_zero_legendre_value(p,th,mmax)
    implicit none
    class(legendre_value), intent(inout) :: p
    real(dp), intent(in) :: th
    integer(i4b), intent(in) :: mmax
    call p%delete()    
    p%l  = 0
    p%th = 0
    p%x = cos(th)
    p%y = sin(th)    
    p%mmax = mmax
    call p%allocate(mmax)       
    p%v(0) = sifourpi    
    return
  end subroutine initialise_zero_legendre_value

  subroutine initialise_legendre_value(p,l,th,mmax)
    implicit none
    class(legendre_value), intent(inout) :: p
    integer(i4b), intent(in) :: l
    real(dp), intent(in) :: th    
    integer(i4b), intent(in) :: mmax    
    p%l  = l
    p%th = th
    p%x = cos(th)
    p%y = sin(th)    
    p%mmax = min(l,mmax)
    call p%allocate(mmax)       
    call legendre(l,p%mmax,th,p%v)           
    return
  end subroutine initialise_legendre_value

  subroutine next_legendre_value(p)
    implicit none
    class(legendre_value), intent(inout) :: p
    integer(i4b) :: l,lold,m,mmax
    real(dp) :: x,y,fac1,fac2
    lold = p%l
    l = lold+1
    x = p%x
    y = p%y        
    mmax = p%mmax
    if(lold == 0) then       
       p%vp1(0) = x*sqrt3*p%v(0)
       if(mmax > 0) then
          p%vp1(1) = -sqrt3*y*p%v(0)/sqrt2
       end if       
    else
       do m = 0,min(mmax,l-2)
          fac1 = (4.0_dp*l*l-1.0_dp)/(l*l-m*m)
          fac1 = sqrt(fac1)
          fac2 = ( (l-1.0_dp)*(l-1.0_dp)-m*m)            &
                  /(4.0_dp*(l-1.0_dp)*(l-1.0_dp)-1.0_dp)
          fac2 = sqrt(fac2)
          p%vp1(m) = fac1*(x*p%v(m) - fac2*p%vm1(m))
       end do
       if(l-1 <= mmax) then
          fac1 = 2.0_dp*l+1
          fac1 = sqrt(fac1)
          p%vp1(l-1) = x*fac1*p%v(l-1)
       end if
       if(l <= mmax) then
          fac1 = (l+0.5_dp)/l
          fac1 = sqrt(fac1)
          p%vp1(l) = -fac1*y*p%v(l-1)
       end if
    end if
    p%l = p%l+1
    p%vm1 = p%v
    p%v = p%vp1
    return
  end subroutine next_legendre_value
  
  
!  function get_single_legendre_value(self,m,flip) result(p)
!    implicit none
!    class(legendre_value), intent(in) :: self
!    integer(i4b), intent(in) :: m
!    logical, intent(in), optional :: flip
!    real(dp) :: p

!    call error(.not.self%allocated,'get_single_legendre_value','not allocated')
!    call error(abs(m) > self%mmax, 'get_single_legendre_value','m out of range')  
!    p = self%p(abs(m))
!    if(m < 0 .and. modulo(m,2) /= 0) p = -p

!    if(present(flip)) then
!       if(flip .and. modulo(self%l+m,2) /= 0) then
!          p = -p          
!       end if
!    end if
    
!    return
!  end function get_single_legendre_value

    
!  function get_slice_legendre_value(self,m1,m2,flip) result(p)
!    implicit none
!    class(legendre_value), intent(in) :: self
!    integer(i4b), intent(in) :: m1,m2
!    real(dp), dimension(m2-m1+1) :: p
!    logical, intent(in), optional :: flip

!    integer(i4b) :: m,im
 
!    call error(.not.self%allocated,'get_single_legendre_value','not allocated')
!    call error(abs(m1) > self%mmax, 'get_single_legendre_value','m out of range')
!    call error(abs(m2) > self%mmax, 'get_single_legendre_value','m out of range')

!    im = 0
!    do m = m1,m2
!       im = im+1
!       p(im) = self%p(abs(m))
!       if(m < 0 .and. modulo(m,2) /= 0) p(im) = -p(im)

!       if(present(flip)) then
!          if(flip .and. modulo(self%l+m,2) /= 0) then
!             p(im) = -p(im)
!          end if
!       end if
       
!    end do
    
!    return
!  end function get_slice_legendre_value
  


  
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



  subroutine legendre_woodhouse(theta,l,m,x,xp,xc)
    
    ! This routine computes the associated legendre 
    ! polynomials P_{lm}(cos(theta)) along with
    ! their angular derivative and  products with
    ! cosec(theta). This routine is a f90 translation
    ! of the routine legndr by John Woodhouse      
    ! the output arrays should have dimension at least as 
    ! big as m+1
    
    ! inputs: 
    !          theta = angular argument in radians
    !          l     = angular degree of the Legendre functions
    !          m     = maximum angular order to calculate
    !
    ! outputs:
    !          x  = array of the Legendre functions; the array
    !               order is such that x(1) = P_{l0}, x(2) = P_{l1},
    !               etc upto x(m+1) = P_{lm}. 
    !          xp = theta derivatives of the Legendre functions
    !          xc = legendre functions times by cosec(theta).
    
    ! the legendre functions returned are such that the fully normalized 
    ! spherical harmonic functions Y_{lm} are given by
    ! 
    ! Y_{lm}(theta,phi) = P_{lm}(cos(theta))*exp(i*m*phi),
    !
    ! The P_{lm}(cos(theta)) calculated are equal to the
    ! X_{lm}(theta) described in appendix B of Dahlen & Tromp (1998).
    
    

    implicit none
    
    ! inputs:
    real(dp), intent(in) :: theta
    integer(i4b), intent(in) :: l,m
    ! outputs:
    real(dp), dimension(:), intent(out) :: x,xp,xc
    
    ! local variables
    integer(i4b) :: lp1,mp1,i,k
    real(dp) :: sum,th,ct,st,fct,sfl3,compar,dsfl3, &
         cot,cosec,x1,x2,x3,f1,f2,xm,small,lsign
    
    sum = 0.0_dp
    lp1 = l+1
    mp1 = m+1
    
    th = theta
    ct = cos(th)
    st = sin(th)
      
    fct    = sqrt(real(2*l+1)/fourpi)
    sfl3   = sqrt(real(l*(l+1)))
    compar = real(2*l+1)/fourpi
    dsfl3  = sfl3
    small  = 1.0e-16_dp*compar
    
    x      = 0.0_dp
    xp     = 0.0_dp
    xc = 0.0_dp
      
    if(l <= 1 .or. abs(theta) <= 1.0e-5_dp) then
       x(1) = fct
       if(l == 0) return     
       x(1)  = ct*fct
       x(2)  = -0.5_dp*st*fct*dsfl3
       xp(1) = -0.5_dp*st*fct*dsfl3*dsfl3
       xp(2)  = -0.5_dp*ct*fct*dsfl3
       if(abs(theta) <  1.0e-5_dp) xc(2) = xp(2)
       if(abs(theta) >= 1.0e-5_dp) xc(2) = x(2)/st
       return
    end if
    
    if(abs(pi-theta) <= 1.0e-4_dp) then
       
       lsign = (-1)**l
       x(1) = -lsign*fct*ct
       x(2) = lsign*0.5_dp*fct*dsfl3*st
       xp(1) = lsign*0.5_dp*fct*dsfl3**2*st
       xp(2) = lsign*0.5_dp*fct*dsfl3*ct
       xc(2) = -xp(2)
       
       return
    end if
    
    x1 = 1.0_dp
    x2 = ct
    
    do i = 2,l
       x3 = (real(2*i-1)*ct*x2-real(i-1)*x1)/real(i)
       x1 = x2
       x2 = x3
    end do
    
    cot   = ct/st;
    cosec = 1.0_dp/st
    
    x3 = x2*fct
    x2 = real(l)*(x1-ct*x2)*fct/st
    
    x(1) = x3
    x(2) = x2
    sum  = x3*x3
    
    xp(1) = -x2
    xp(2) = real(l*(l+1))*x3-cot*x2
    
    x(2)      = -x(2)/sfl3
    xc(2) = x(2)*cosec
    xp(2)     = -xp(2)/sfl3
    
    sum = sum+2.0_dp*x(2)*x(2)
    if(sum-compar > small) return
    
    x1 =  x3
    x2 = -x2/sqrt(real(l*(l+1)))
    do i = 3,mp1
       k   = i-1
       f1  = sqrt(real(l+i-1)*(l-i+2))
       f2  = sqrt(real(l+i-2)*(l-i+3))
       xm  = k
       x3  = -(2.0_dp*cot*(xm-1.0_dp)*x2+f2*x1)/f1
       sum = sum+2.0_dp*x3*x3
       if(sum-compar > small .and. i /= lp1) return
       x(i)      = x3
       xc(i) = x(i)*cosec
       x1        = x2
       xp(i)     = -(f1*x2+xm*cot*x3)
       x2        = x3
    end do
    
    return
  end subroutine legendre_woodhouse



  
  
  
end module module_spherical_harmonics
