program fftw3_test

  use module_constants
  use module_fftw3
  implicit none
  
  integer, parameter :: n = 8, m = 2
  integer :: i
  real(dp) :: x,start,finish
  real(dp), dimension(n) :: inr,recr
  real(dp), dimension(n,m) :: invr,recvr
  complex(dpc), dimension(n/2+1) :: outr
  complex(dpc), dimension(n/2+1,m) :: outvr
  complex(dpc), dimension(n) :: in,out,rec
  complex(dpc), dimension(n,m) :: inv,outv,recv


  do i = 1,n
     x = (i-1)*pi/n
     in(i) = sin(5*x)
     inv(i,:) = in(i)
  end do
  
  print *, real(in)  
  call dft_1D (.true., n,  in, out)
  call dft_1D(.false., n, out, rec)
  print *, real(rec/n)
  print *, '============================================'

  
  print *, real(inv(:,1))  
  call dft_1D( .true., n, m,  inv, outv)
  call dft_1D(.false., n, m, outv, recv)  
  print *, real(recv(:,1)/n)
  print *, '============================================'


  inr = real(in)
  print *, inr
  call dft_1D(n,inr,outr)
  call dft_1D(n,outr,recr)
  print *, recr/n
  print *, '============================================'

  do i = 1,m
     invr(:,i) = inr
  end do
  print *, invr(:,1)
  call dft_1D(n,m,invr,outvr)
  call dft_1D(n,m,outvr,recvr)
  print *, recvr(:,1)/n
  print *, '============================================'
  
end program fftw3_test
