program fftw3_test

  use module_constants
  use module_fftw3
  implicit none
  
  integer, parameter :: n = 8
  integer :: i
  real(dp) :: x,start,finish
  complex(dpc), dimension(n) :: in,out,rec
  

  do i = 1,8
     x = twopi*(i-1)/n
     in(i) = exp(-4*ii*x)
  end do

  call dft_1D(.true.,n,in,out)
  out = out/n

  do i = 1,n/2+1
     print *, i-1,out(i)
  end do
  do i = n/2+2,n
     print *, -(n-i+1),out(i)
  end do

  
end program fftw3_test
