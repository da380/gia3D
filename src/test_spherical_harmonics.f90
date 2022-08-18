program test_spherical_harmonics

  use module_constants
  use module_special_functions
  use module_spherical_harmonics
  use module_fftw3
  implicit none

  integer(i4b) :: lmax,nmax,n,l,m,ith,iph,ilm,k,klm,alpha,i
  real(dp) :: start,finish,th,ph,rand1,rand2,rfun
  complex(dpc) :: cfun
  type(wigner_value) :: d
  type(gauss_legendre_grid) :: grid

  ! scalar fields
  type(scalar_gauss_legendre_field) :: u,v
  type(scalar_spherical_harmonic_expansion) :: ulm,vlm

  ! vector fields
  type(vector_gauss_legendre_field) :: w,x
  type(vector_spherical_harmonic_expansion) :: wlm,xlm

  ! real vector fields
  type(vector_gauss_legendre_field) :: y,z
  type(vector_spherical_harmonic_expansion) :: ylm,zlm
 
  ! make the grid  
  lmax = 512
  nmax = 1
  print *, 'building the grid for degree ',lmax
  call cpu_time(start)
  call grid%allocate(lmax,nmax)  
  call cpu_time(finish)
  print *, 'assembly time = ',finish-start

  
  !=============================================================!
  !=============================================================!
  !                   tests for scalar fields                   !
  !=============================================================!
  !=============================================================!

  print *, '#############################################'
  print *, '#             scalar field tests            #'
  print *, '#############################################'
  
  !=============================================================!
  !            test coefficient to function and back            !
  !=============================================================!
  
  ! allocate a coefficient array
  call ulm%allocate(grid)

  ! set values for the coefficients
  do l = 0,lmax
     do m = -l,l
        if(m == -lmax) cycle
        call random_number(rand1)
        call random_number(rand2)
        rand1 = 2.0_dp*(rand1-0.5_dp)
        rand2 = 2.0_dp*(rand2-0.5_dp)
        cfun = (rand1+ii*rand2)*(1.0_dp+0.1*l*(l+1))**(-1.5)
        call ulm%set(l,m,cfun)
     end do
  end do



  ! allocate a scalar field
  call u%allocate(grid)

  ! do an inverse transformation
  call cpu_time(start)
  call grid%SH_itrans(ulm,u)
  call cpu_time(finish)
  print *, 'transformation time = ',finish-start


  ! now transform back again
  call vlm%allocate(grid)
  call cpu_time(start)
  call grid%SH_trans(u,vlm)
  call cpu_time(finish)
  print *, 'transformation time = ',finish-start
  
  print *, 'error in coefficients = ', maxval(abs(ulm%cdata-vlm%cdata))
  print *, '-----------------------------------------'



  !=============================================================!
  !            test function to coefficient and back            !
  !=============================================================!
  
  ! set a function
  l = lmax-2
  m = -lmax/4
  klm = l*(l+1)/2 + abs(m) + 1
  do ith = 1,grid%nth
     th = grid%th(ith)
     do iph = 1,grid%nph
        ph = grid%ph(iph)
        if(m >= 0) then
           cfun = grid%dlm(klm,0,ith)*exp(ii*m*ph)
        else
           cfun = (-1)**m*grid%dlm(klm,0,ith)*exp(ii*m*ph)
        end if
        call u%set(iph,ith,cfun)
     end do
  end do
 
  ! compute the transform
  call cpu_time(start)
  call grid%SH_trans(u,ulm)
  call cpu_time(finish)
  print *, 'transformation time = ',finish-start

  print *, '(l,m)th coefficient  = ', ulm%get(l,m)
  
  ! transform back
  call v%allocate(grid)
  call cpu_time(start)
  call grid%SH_itrans(ulm,v)
  call cpu_time(finish)
  print *, 'transformation time = ',finish-start
  
  print *, 'error in point values = ', maxval(abs(u%cdata-v%cdata))


  print *, '#############################################'
  print *, '#           real  scalar field tests        #'
  print *, '#############################################'
  
  !=============================================================!
  !            test coefficient to function and back            !
  !=============================================================!
  
  ! allocate a coefficient array
  call ulm%allocate(grid,real = .true.)

  ! set values for the coefficients
  do l = 0,lmax
     do m = 0,l
        call random_number(rand1)
        call random_number(rand2)
        rand1 = 2.0_dp*(rand1-0.5_dp)
        if(m == 0 .or. m == lmax) then
           rand2 = 0
        else
           call random_number(rand2)
           rand2 = 2.0_dp*(rand2-0.5_dp)
        end if
        cfun = (rand1+ii*rand2)*(1.0_dp+0.1*l*(l+1))**(-1.5)
        call ulm%set(l,m,cfun)
     end do
  end do

  ! allocate a scalar field
  call u%allocate(grid,real=.true.)

  
  ! do an inverse transformation
  call cpu_time(start)
  call grid%SH_itrans(ulm,u)
  call cpu_time(finish)
  print *, 'transformation time = ',finish-start

  
  ! now transform back again
  call vlm%allocate(grid,real=.true.)
  call cpu_time(start)
  call grid%SH_trans(u,vlm)
  call cpu_time(finish)
  print *, 'transformation time = ',finish-start
  
  print *, 'error in coefficients = ', maxval(abs(ulm%rdata-vlm%rdata))
  print *, '-----------------------------------------'



  !=============================================================!
  !            test function to coefficient and back            !
  !=============================================================!
  
  ! set a function
  l = lmax
  m = lmax/2
  klm = l*(l+1)/2 + abs(m) + 1
  do ith = 1,grid%nth
     th = grid%th(ith)
     do iph = 1,grid%nph
        ph = grid%ph(iph)
        if(m >= 0) then
           cfun = grid%dlm(klm,0,ith)*exp(ii*m*ph)
        else
           cfun = (-1)**m*grid%dlm(klm,0,ith)*exp(ii*m*ph)
        end if
        rfun = real(cfun)
        call u%set(iph,ith,rfun)
     end do
  end do

  ! compute the transform
  call cpu_time(start)
  call grid%SH_trans(u,ulm)
  call cpu_time(finish)
  print *, 'transformation time = ',finish-start

  print *, '(l,m)th coefficient  = ', ulm%get(l,m)
  
  ! transform back
  call v%allocate(grid,real=.true.)
  call cpu_time(start)
  call grid%SH_itrans(ulm,v)
  call cpu_time(finish)
  print *, 'transformation time = ',finish-start
  
  print *, 'error in point values = ', maxval(abs(u%rdata-v%rdata))

  
  print *, '#############################################'
  print *, '#             vector field tests            #'
  print *, '#############################################'

  !=============================================================!
  !=============================================================!
  !                   tests for vector fields                   !
  !=============================================================!
  !=============================================================!




  !=============================================================!
  !            test coefficient to function and back            !
  !=============================================================!
  
  ! allocate a coefficient array
  call wlm%allocate(grid)

  
  ! set values for the coefficients
  do l = 0,lmax
     do m = -l,l
        if(m == -lmax) cycle
        ! set alpha = 0 component
        call random_number(rand1)
        call random_number(rand2)
        rand1 = 2.0_dp*(rand1-0.5_dp)
        rand2 = 2.0_dp*(rand2-0.5_dp)
        cfun = (rand1+ii*rand2)*(1.0_dp+0.1*l*(l+1))**(-1.5)
        call wlm%set(l,m,0,cfun)
        if(l == 0) cycle
        ! set alpha = -1 component
        call random_number(rand1)
        call random_number(rand2)
        rand1 = 2.0_dp*(rand1-0.5_dp)
        rand2 = 2.0_dp*(rand2-0.5_dp)
        cfun = (rand1+ii*rand2)*(1.0_dp+0.1*l*(l+1))**(-1.5)
        call wlm%set(l,m,-1,cfun)
        ! set alpha = +1 component
        call random_number(rand1)
        call random_number(rand2)
        rand1 = 2.0_dp*(rand1-0.5_dp)
        rand2 = 2.0_dp*(rand2-0.5_dp)
        cfun = (rand1+ii*rand2)*(1.0_dp+0.1*l*(l+1))**(-1.5)
        call wlm%set(l,m,1,cfun)
     end do
  end do
  
  ! allocate a vector field
  call w%allocate(grid)

  ! do an inverse transformation
  call cpu_time(start)
  call grid%SH_itrans(wlm,w)
  call cpu_time(finish)
  print *, 'transformation time = ',finish-start


  ! now transform back again
  call xlm%allocate(grid)
  call cpu_time(start)
  call grid%SH_trans(w,xlm)
  call cpu_time(finish)
  print *, 'transformation time = ',finish-start

  print *, 'error in coefficients = ', maxval(abs(wlm%cdata-xlm%cdata))
  print *, '-----------------------------------------'


  !=============================================================!
  !            test function to coefficient and back            !
  !=============================================================!
  
  ! set a function
  l = lmax/4-1
  m = -lmax/8
  alpha = -1
  klm = l*(l+1)/2 + abs(m) + 1
  do ith = 1,grid%nth
     th = grid%th(ith)
     do iph = 1,grid%nph
        ph = grid%ph(iph)
        if(m >= 0) then
           cfun = grid%dlm(klm,alpha,ith)*exp(ii*m*ph)
        else
           cfun = (-1)**(m+alpha)*grid%dlm(klm,-alpha,ith)*exp(ii*m*ph)
        end if
        call w%set(iph,ith,alpha,cfun)
     end do
  end do

  ! compute the transform
  call cpu_time(start)
  call grid%SH_trans(w,wlm)
  call cpu_time(finish)
  print *, 'transformation time = ',finish-start


  print *, '(l,m,alpha)th coefficient  = ', wlm%get(l,m,alpha)
  
  ! transform back
  call x%allocate(grid)
  call cpu_time(start)
  call grid%SH_itrans(wlm,x)
  call cpu_time(finish)
  print *, 'transformation time = ',finish-start

  
  print *, 'error in point values = ', maxval(abs(w%cdata-x%cdata))

  

  print *, '#############################################'
  print *, '#          real vector field tests          #'
  print *, '#############################################'


  
  !=============================================================!
  !=============================================================!
  !                   tests for vector fields                   !
  !=============================================================!
  !=============================================================!


  !=============================================================!
  !            test coefficient to function and back            !
  !=============================================================!


  
  ! allocate a coefficient array
  call ylm%allocate(grid,real = .true.)

  
  ! set values for the coefficients
  do l = 0,lmax
     ! alpha = 0, m = 0
     call random_number(rand1)
     rand1 = 2.0_dp*(rand1-0.5_dp)
     rand2 = 0.0_dp
     cfun = (rand1+ii*rand2)*(1.0_dp+0.1*l*(l+1))**(-1.5)
     call ylm%set(l,0,0,cfun)
     ! alpha = -1, m = 0
     call random_number(rand1)
     call random_number(rand2)
     rand1 = 2.0_dp*(rand1-0.5_dp)
     rand2 = 2.0_dp*(rand2-0.5_dp)
     cfun = (rand1+ii*rand2)*(1.0_dp+0.1*l*(l+1))**(-1.5)
     if(l > 0) call ylm%set(l,0,-1,cfun)        
     do m = 1,l

        ! do alpha = 0, m positive
        call random_number(rand1)
        rand1 = 2.0_dp*(rand1-0.5_dp)
        if(m == lmax) then
           rand2 = 0.0_dp
        else
           call random_number(rand2)
           rand2 = 2.0_dp*(rand2-0.5_dp)
        end if
        cfun = (rand1+ii*rand2)*(1.0_dp+0.1*l*(l+1))**(-1.5)
        call ylm%set(l,m,0,cfun)

        ! do alpha = -1, m positive
        call random_number(rand1)
        call random_number(rand2)
        rand1 = 2.0_dp*(rand1-0.5_dp)
        rand2 = 2.0_dp*(rand2-0.5_dp)
        cfun = (rand1+ii*rand2)*(1.0_dp+0.1*l*(l+1))**(-1.5)
        call ylm%set(l,m,-1,cfun)

        ! do alpha = -1, m negative
        call random_number(rand1)
        call random_number(rand2)
        rand1 = 2.0_dp*(rand1-0.5_dp)
        rand2 = 2.0_dp*(rand2-0.5_dp)
        cfun = (rand1+ii*rand2)*(1.0_dp+0.1*l*(l+1))**(-1.5)
        if(m < lmax) call ylm%set(l,-m,-1,cfun)
        
     end do
     
  end do



  
  ! allocate a real vector field
  call y%allocate(grid,real = .true.)

  
  ! do an inverse transformation
  call cpu_time(start)
  call grid%SH_itrans(ylm,y)
  call cpu_time(finish)
  print *, 'transformation time = ',finish-start


  ! now transform back again
  call zlm%allocate(grid,real = .true.)
  call cpu_time(start)
  call grid%SH_trans(y,zlm)
  call cpu_time(finish)
  print *, 'transformation time = ',finish-start

  print *, 'error in coefficients = ', maxval(abs(ylm%cdata-zlm%cdata))
  print *, '-----------------------------------------'


  !=============================================================!
  !            test function to coefficient and back            !
  !=============================================================!
  
  ! set a function
  alpha = 0
  l = lmax-1
  m = -lmax/4
  if(alpha == 0) m = abs(m)
  klm = l*(l+1)/2 + abs(m) + 1
  do ith = 1,grid%nth
     th = grid%th(ith)
     do iph = 1,grid%nph
        ph = grid%ph(iph)
        if(m >= 0) then
           cfun = grid%dlm(klm,alpha,ith)*exp(ii*m*ph)
        else
           cfun = (-1)**(m+alpha)*grid%dlm(klm,-alpha,ith)*exp(ii*m*ph)
        end if
        call y%set(iph,ith,alpha,cfun)
     end do
  end do

  ! compute the transform
  call cpu_time(start)
  call grid%SH_trans(y,ylm)
  call cpu_time(finish)
  print *, 'transformation time = ',finish-start


  print *, '(l,m,alpha)th coefficient  = ', ylm%get(l,m,alpha)

  
  ! transform back
  call z%allocate(grid,real = .true.)
  call cpu_time(start)
  call grid%SH_itrans(ylm,z)
  call cpu_time(finish)
  print *, 'transformation time = ',finish-start


  print *, 'error in point values = ', maxval(abs(y%cdata-z%cdata))



  
end program test_spherical_harmonics




