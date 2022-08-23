program test_spherical_harmonics

  use module_constants
  use module_special_functions
  use module_spherical_harmonics
  use module_fftw3
  implicit none

  integer(i4b) :: lmax,nmax,n,l,m,ith,iph,ilm,k,klm,alpha,i,beta
  real(dp) :: start,finish,th,ph,rand1,rand2,rfun
  complex(dpc) :: cfun
  type(gauss_legendre_grid) :: grid

  ! scalar fields
  type(scalar_gauss_legendre_field) :: u,v
  type(scalar_spherical_harmonic_expansion) :: ulm,vlm

  ! real scalar fields
  type(real_scalar_gauss_legendre_field) :: w,x
  type(real_scalar_spherical_harmonic_expansion) :: wlm,xlm

  ! vector fields
  type(vector_gauss_legendre_field) :: y,z
  type(vector_spherical_harmonic_expansion) :: ylm,zlm


  ! real vector fields
  type(real_vector_gauss_legendre_field) :: a,b
  type(real_vector_spherical_harmonic_expansion) :: alm,blm

  
  ! internal variable fields
  type(internal_variable_gauss_legendre_field) :: c,d
  type(internal_variable_spherical_harmonic_expansion) :: clm,dlm
  

  
  ! make the grid  
  lmax = 512
  nmax = 2
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

  ! set random values
  call ulm%set_sobolev(2.0_dp,0.5_dp)
  call ulm%random()
  
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
  
  ! set function as single harmonic
  l = lmax
  m = -lmax/2
  call u%harmonic(grid,l,m)
  
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
  call wlm%allocate(grid)

  ! set random values
  call wlm%set_sobolev(2.0_dp,0.5_dp)
  call wlm%random()

  ! allocate a scalar field
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
  
  print *, 'error in coefficients = ', maxval(abs(wlm%rdata-xlm%rdata))
  print *, '-----------------------------------------'


  !=============================================================!
  !            test function to coefficient and back            !
  !=============================================================!


  ! set function as single harmonic
  l = lmax
  m = -lmax/2
  call w%harmonic(grid,l,m)
  
  ! compute the transform
  call cpu_time(start)
  call grid%SH_trans(w,wlm)
  call cpu_time(finish)
  print *, 'transformation time = ',finish-start

  print *, '(l,m)th coefficient  = ', wlm%get(l,m)
  
  ! transform back
  call x%allocate(grid)
  call cpu_time(start)
  call grid%SH_itrans(wlm,x)
  call cpu_time(finish)
  print *, 'transformation time = ',finish-start
  
  print *, 'error in point values = ', maxval(abs(w%rdata-x%rdata))


  print *, '#############################################'
  print *, '#             vector field tests            #'
  print *, '#############################################'

  
  !=============================================================!
  !            test coefficient to function and back            !
  !=============================================================!
  
  ! allocate a coefficient array
  call ylm%allocate(grid)

  ! set random values
  call ylm%set_sobolev(2.0_dp,0.5_dp)
  call ylm%random()
    
  ! allocate a vector field
  call y%allocate(grid)

  ! do an inverse transformation
  call cpu_time(start)
  call grid%SH_itrans(ylm,y)
  call cpu_time(finish)
  print *, 'transformation time = ',finish-start

  ! now transform back again
  call zlm%allocate(grid)
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
  l = lmax/4-1
  m = -lmax/8
  alpha = -1
  call y%harmonic(grid,l,m,alpha)

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

  
  print *, 'error in point values = ', maxval(abs(y%cdata-z%cdata))


  print *, '#############################################'
  print *, '#          real vector field tests          #'
  print *, '#############################################'



  !=============================================================!
  !            test coefficient to function and back            !
  !=============================================================!
  
  ! allocate a coefficient array
  call alm%allocate(grid)
  
  ! set random values
  call alm%set_sobolev(2.0_dp,0.5_dp)
  call alm%random()
  
  ! allocate a real vector field
  call a%allocate(grid)
  
  ! do an inverse transformation
  call cpu_time(start)
  call grid%SH_itrans(alm,a)
  call cpu_time(finish)
  print *, 'transformation time = ',finish-start

  ! now transform back again
  call blm%allocate(grid)
  call cpu_time(start)
  call grid%SH_trans(a,blm)
  call cpu_time(finish)
  print *, 'transformation time = ',finish-start

  print *, 'error in complex coefficients = ', maxval(abs(alm%cdata-blm%cdata))
  print *, 'error in real coefficients = ', maxval(abs(alm%rdata-blm%rdata))
  print *, '-----------------------------------------'


  !=============================================================!
  !            test function to coefficient and back            !
  !=============================================================!
  
  ! set a function
  alpha = 0
  l = lmax-1
  m = lmax/4
  call a%harmonic(grid,l,m,alpha)

  ! compute the transform
  call cpu_time(start)
  call grid%SH_trans(a,alm)
  call cpu_time(finish)
  print *, 'transformation time = ',finish-start

  print *, '(l,m,alpha)th coefficient  = ', alm%get(l,m,alpha)

  ! transform back
  call b%allocate(grid)
  call cpu_time(start)
  call grid%SH_itrans(alm,b)
  call cpu_time(finish)
  print *, 'transformation time = ',finish-start


  print *, 'error in complex point values = ', maxval(abs(a%cdata-b%cdata))
  print *, 'error in real point values = ', maxval(abs(a%rdata-b%rdata))



  print *, '#############################################'
  print *, '#       internal variable field tests       #'
  print *, '#############################################'


  
  !=============================================================!
  !            test coefficient to function and back            !
  !=============================================================!
  
  ! allocate a coefficient array
  call clm%allocate(grid)
  
  ! set random values
  call clm%set_sobolev(2.0_dp,0.5_dp)
  call clm%random()
  
  ! allocate a real vector field
  call c%allocate(grid)
  
  ! do an inverse transformation
  call cpu_time(start)
  call grid%SH_itrans(clm,c)
  call cpu_time(finish)
  print *, 'transformation time = ',finish-start

  ! now transform back again
  call dlm%allocate(grid)
  call cpu_time(start)
  call grid%SH_trans(c,dlm)
  call cpu_time(finish)
  print *, 'transformation time = ',finish-start

  print *, 'error in complex coefficients = ', maxval(abs(clm%cdata-dlm%cdata))
  print *, 'error in real coefficients = ', maxval(abs(clm%rdata-dlm%rdata))
  print *, '-----------------------------------------'


  !=============================================================!
  !            test function to coefficient and back            !
  !=============================================================!
  
  ! set a function
  alpha = -1
  beta  = -1
  l = lmax-1
  m = lmax/4
  call c%harmonic(grid,l,m,alpha,beta)

  ! compute the transform
  call cpu_time(start)
  call grid%SH_trans(c,clm)
  call cpu_time(finish)
  print *, 'transformation time = ',finish-start

  print *, '(l,m,alpha,beta)th coefficient  = ', clm%get(l,m,alpha,beta)

  ! transform back
  call d%allocate(grid)
  call cpu_time(start)
  call grid%SH_itrans(clm,d)
  call cpu_time(finish)
  print *, 'transformation time = ',finish-start


  print *, 'error in complex point values = ', maxval(abs(c%cdata-d%cdata))
  print *, 'error in real point values = ', maxval(abs(c%rdata-d%rdata))




  
  
end program test_spherical_harmonics




