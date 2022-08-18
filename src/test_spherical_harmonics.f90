program test_spherical_harmonics

  use module_constants
  use module_special_functions
  use module_spherical_harmonics
  use module_fftw3
  implicit none

  integer(i4b) :: lmax,nmax,n,l,m,ith,iph,ilm,k,klm,alpha
  real(dp) :: start,finish,th,ph,rand1,rand2
  complex(dpc) :: fun
  type(wigner_value) :: d
  type(gauss_legendre_grid) :: grid

  ! scalar fields
  type(scalar_gauss_legendre_field) :: u,v
  type(scalar_spherical_harmonic_expansion) :: ulm,vlm

  ! vector fields
  type(vector_gauss_legendre_field) :: w,x
  type(vector_spherical_harmonic_expansion) :: wlm,xlm

  ! real vector fields
  type(real_vector_gauss_legendre_field) :: y,z
  type(real_vector_spherical_harmonic_expansion) :: ylm,zlm


  ! set the expansion degree
  lmax = 512
  

  !=============================================================!
  !=============================================================!
  !                   tests for scalar fields                   !
  !=============================================================!
  !=============================================================!

  print *, '#############################################'
  print *, '#             scalar field tests            #'
  print *, '#############################################'
  
  ! make the GL-grid
  nmax = 0  
  call grid%allocate(lmax,nmax)


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
        fun = (rand1+ii*rand2)*(1.0_dp+0.1*l*(l+1))**(-1.5)
        call ulm%set(l,m,fun)
     end do
  end do



  ! allocate a scalar field
  call u%allocate(grid)

  ! do an inverse transformation
  call cpu_time(start)
  call grid%SH_itrans(ulm,u)
  call cpu_time(finish)
  print *, 'transformation time = ',finish-start

  
  ! write out the scalar field
  open(99,file='scalar.dat')
  write(99,*) grid%nth,grid%nph,0.0_dp
  do ith = 1,grid%nth
     th = grid%th(ith)
     do iph = 1,grid%nph
        ph = grid%ph(iph)
        write(99,*) ph,th,real(u%get(iph,ith))
     end do
  end do
  close(99)

  ! now transform back again
  call vlm%allocate(grid)
  call cpu_time(start)
  call grid%SH_trans(u,vlm)
  call cpu_time(finish)
  print *, 'transformation time = ',finish-start
  
  print *, 'error in coefficients = ', maxval(abs(ulm%data-vlm%data))
  print *, '-----------------------------------------'



  !=============================================================!
  !            test function to coefficient and back            !
  !=============================================================!
  
  ! set a function
  l = 5
  m = -3
  klm = l*(l+1)/2 + abs(m) + 1
  do ith = 1,grid%nth
     th = grid%th(ith)
     do iph = 1,grid%nph
        ph = grid%ph(iph)
        if(m >= 0) then
           fun = grid%dlm(klm,0,ith)*exp(ii*m*ph)
        else
           fun = (-1)**m*grid%dlm(klm,0,ith)*exp(ii*m*ph)
        end if
        call u%set(iph,ith,fun)
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

  ! write out the difference field
  open(99,file='scalar_dif.dat')
  write(99,*) grid%nth,grid%nph,0.0_dp
  do ith = 1,grid%nth
     th = grid%th(ith)
     do iph = 1,grid%nph
        ph = grid%ph(iph)       
        write(99,*) ph,th,abs(u%get(iph,ith)- v%get(iph,ith))
     end do
  end do
  close(99)
  
  print *, 'error in point values = ', maxval(abs(u%data-v%data))


  
  print *, '#############################################'
  print *, '#             vector field tests            #'
  print *, '#############################################'

  !=============================================================!
  !=============================================================!
  !                   tests for vector fields                   !
  !=============================================================!
  !=============================================================!


  ! make the GL-grid
  nmax = 1
  call grid%allocate(lmax,nmax)


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
        fun = (rand1+ii*rand2)*(1.0_dp+0.1*l*(l+1))**(-1.5)
        call wlm%set(l,m,0,fun)
        if(l == 0) cycle
        ! set alpha = -1 component
        call random_number(rand1)
        call random_number(rand2)
        rand1 = 2.0_dp*(rand1-0.5_dp)
        rand2 = 2.0_dp*(rand2-0.5_dp)
        fun = (rand1+ii*rand2)*(1.0_dp+0.1*l*(l+1))**(-1.5)
        call wlm%set(l,m,-1,fun)
        ! set alpha = +1 component
        call random_number(rand1)
        call random_number(rand2)
        rand1 = 2.0_dp*(rand1-0.5_dp)
        rand2 = 2.0_dp*(rand2-0.5_dp)
        fun = (rand1+ii*rand2)*(1.0_dp+0.1*l*(l+1))**(-1.5)
        call wlm%set(l,m,1,fun)
     end do
  end do
  
  ! allocate a vector field
  call w%allocate(grid)

  ! do an inverse transformation
  call cpu_time(start)
  call grid%SH_itrans(wlm,w)
  call cpu_time(finish)
  print *, 'transformation time = ',finish-start


  ! write out the vector field
  open(99,file='vector.dat')
  write(99,*) grid%nth,grid%nph,0.0_dp
  do ith = 1,grid%nth
     th = grid%th(ith)
     do iph = 1,grid%nph
        ph = grid%ph(iph)
        write(99,*) ph,th,real(w%get(iph,ith,0))
     end do
  end do
  close(99)

  ! now transform back again
  call xlm%allocate(grid)
  call cpu_time(start)
  call grid%SH_trans(w,xlm)
  call cpu_time(finish)
  print *, 'transformation time = ',finish-start

  print *, 'error in coefficients = ', maxval(abs(wlm%data-xlm%data))
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
           fun = grid%dlm(klm,alpha,ith)*exp(ii*m*ph)
        else
           fun = (-1)**(m+alpha)*grid%dlm(klm,-alpha,ith)*exp(ii*m*ph)
        end if
        call w%set(iph,ith,alpha,fun)
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

  ! write out the difference field
  open(99,file='vector_dif.dat')
  write(99,*) grid%nth,grid%nph,0.0_dp
  do ith = 1,grid%nth
     th = grid%th(ith)
     do iph = 1,grid%nph
        ph = grid%ph(iph)       
        write(99,*) ph,th,real(w%get(iph,ith,alpha)- x%get(iph,ith,alpha))
     end do
  end do
  close(99)
  
  print *, 'error in point values = ', maxval(abs(w%data-x%data))

  

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
  call ylm%allocate(grid)

  
  ! set values for the coefficients
  do l = 0,lmax
     do m = -l,l
        if(m == -lmax) cycle
        ! set alpha = 0 component
        call random_number(rand1)
        call random_number(rand2)
        rand1 = 2.0_dp*(rand1-0.5_dp)
        rand2 = 2.0_dp*(rand2-0.5_dp)
        fun = (rand1+ii*rand2)*(1.0_dp+0.1*l*(l+1))**(-1.5)
        call ylm%set(l,m,0,fun)
        if(l == 0) cycle
        ! set alpha = -1 component
        call random_number(rand1)
        call random_number(rand2)
        rand1 = 2.0_dp*(rand1-0.5_dp)
        rand2 = 2.0_dp*(rand2-0.5_dp)
        fun = (rand1+ii*rand2)*(1.0_dp+0.1*l*(l+1))**(-1.5)
        call ylm%set(l,m,-1,fun)        
     end do
  end do
  
  ! allocate a real vector field
  call y%allocate(grid)

  ! do an inverse transformation
  call cpu_time(start)
  call grid%SH_itrans(ylm,y)
  call cpu_time(finish)
  print *, 'transformation time = ',finish-start


  ! write out the vector field
  open(99,file='real_vector.dat')
  write(99,*) grid%nth,grid%nph,0.0_dp
  do ith = 1,grid%nth
     th = grid%th(ith)
     do iph = 1,grid%nph
        ph = grid%ph(iph)
        write(99,*) ph,th,real(y%get(iph,ith,0))
     end do
  end do
  close(99)

  ! now transform back again
  call zlm%allocate(grid)
  call cpu_time(start)
  call grid%SH_trans(y,zlm)
  call cpu_time(finish)
  print *, 'transformation time = ',finish-start

  print *, 'error in coefficients = ', maxval(abs(ylm%data-zlm%data))
  print *, '-----------------------------------------'


  !=============================================================!
  !            test function to coefficient and back            !
  !=============================================================!
  
  ! set a function
  l = lmax-1
  m = -lmax/4
  alpha = 0
  klm = l*(l+1)/2 + abs(m) + 1
  do ith = 1,grid%nth
     th = grid%th(ith)
     do iph = 1,grid%nph
        ph = grid%ph(iph)
        if(m >= 0) then
           fun = grid%dlm(klm,alpha,ith)*exp(ii*m*ph)
        else
           fun = (-1)**(m+alpha)*grid%dlm(klm,-alpha,ith)*exp(ii*m*ph)
        end if
        call y%set(iph,ith,alpha,fun)
     end do
  end do

  ! compute the transform
  call cpu_time(start)
  call grid%SH_trans(y,ylm)
  call cpu_time(finish)
  print *, 'transformation time = ',finish-start


  print *, '(l,m,alpha)th coefficient  = ', ylm%get(l,m,alpha)

  
  ! transform back
  call z%allocate(grid)
  call cpu_time(start)
  call grid%SH_itrans(ylm,z)
  call cpu_time(finish)
  print *, 'transformation time = ',finish-start

  ! write out the difference field
  open(99,file='real_vector_dif.dat')
  write(99,*) grid%nth,grid%nph,0.0_dp
  do ith = 1,grid%nth
     th = grid%th(ith)
     do iph = 1,grid%nph
        ph = grid%ph(iph)       
        write(99,*) ph,th,real(y%get(iph,ith,alpha)-z%get(iph,ith,alpha))
     end do
  end do
  close(99)
  
  print *, 'error in point values = ', maxval(abs(y%data-z%data))



  
end program test_spherical_harmonics




