module module_spherical_harmonics

  use module_constants
  use module_error


  type legendre_value
     private
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
     procedure :: init => initialise_legendre_value
     procedure :: next => next_legendre_value
     procedure :: deg => degree_legendre_value
     procedure :: get_single_legendre_value
     procedure :: get_slice_legendre_value
     procedure :: get_all_legendre_value     
     generic   :: get => get_single_legendre_value, &
                         get_slice_legendre_value,  &
                         get_all_legendre_value
     
  end type legendre_value

  
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
    p%th = 0
    p%x = cos(th)
    p%y = sin(th)    
    p%mmax = mmax
    call p%allocate(mmax)       
    p%v(0) = sifourpi    
    return
  end subroutine initialise_legendre_value


  subroutine next_legendre_value(p)
    implicit none
    class(legendre_value), intent(inout) :: p
    integer(i4b) :: l,lold,m,mmax,im
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


  function get_all_legendre_value(p,flip) result(v)
    implicit none
    class(legendre_value), intent(in) :: p
    real(dp), dimension(0:min(p%l,p%mmax)) :: v
    logical, intent(in), optional :: flip
    integer(i4b) :: m,im
    call error(.not.p%allocated,'get_single_legendre_value','not allocated')
    v = p%v(0:min(p%l,p%mmax))
    return
  end function get_all_legendre_value
  
  
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






  
  
  
end module module_spherical_harmonics
