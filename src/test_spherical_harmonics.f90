program test_spherical_harmonics

  use module_constants
  use module_special_functions
  use module_spherical_harmonics
  use module_fftw3
  implicit none

  integer(i4b) :: lmax,nmax,n,l,m,ith,iph,ilm,k
  real(dp) :: start,finish,th,ph,rand1,rand2
  complex(dpc) :: fun
  type(wigner_value) :: d
  type(gauss_legendre_grid) :: grid
  type(scalar_gauss_legendre_field) :: u,v
  type(scalar_spherical_harmonic_expansion) :: ulm,vlm

  ! make the GL-grid
  lmax = 512
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
        call random_number(rand1)
        call random_number(rand2)
        rand1 = 2.0_dp*(rand1-0.5_dp)
        rand2 = 2.0_dp*(rand2-0.5_dp)
        if(m == -lmax) cycle
        ulm%data(ulm%index(m,l,1)) = (rand1+ii*rand2)*(1.0_dp+0.4*l*(l+1))**(-0.5)
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
  open(99,file='initial.dat')
  write(99,*) grid%nth,grid%nph,0.0_dp
  do ith = 1,grid%nth
     th = grid%th(ith)
     do iph = 1,grid%nph
        ph = grid%ph(iph)
        write(99,*) ph,th,real(u%data(u%index(iph,ith,1)))        
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
  do ith = 1,grid%nth
     th = grid%th(ith)
     do iph = 1,grid%nph
        ph = grid%ph(ith)
        u%data(u%index(iph,ith,1)) =  exp(-4*(th-pio2)**2)*sin(6*ph)*cos(th)
     end do
  end do

  ! compute the transform
  call cpu_time(start)
  call grid%SH_trans(u,ulm)
  call cpu_time(finish)
  print *, 'transformation time = ',finish-start


  ! transform back
  call v%allocate(grid)
  call cpu_time(start)
  call grid%SH_itrans(ulm,v)
  call cpu_time(finish)
  print *, 'transformation time = ',finish-start


  ! write out the difference field
  open(99,file='dif.dat')
  write(99,*) grid%nth,grid%nph,0.0_dp
  do ith = 1,grid%nth
     th = grid%th(ith)
     do iph = 1,grid%nph
        ph = grid%ph(iph)
        write(99,*) ph,th,abs(u%data(u%index(iph,ith,1)) &
                             -v%data(u%index(iph,ith,1)))        
     end do
  end do
  close(99)
  
  print *, 'error in point values = ', maxval(abs(u%data-v%data))
  
end program test_spherical_harmonics




