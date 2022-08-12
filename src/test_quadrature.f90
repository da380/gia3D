program quadrature_test


  use module_constants
  use module_special_functions
  use module_quadrature
  implicit none

  integer(i4b), parameter :: n = 2**2+1
  class(orthogonal_polynomial), allocatable :: poly
  type(gauss_quadrature) :: quad1
  type(gauss_radau_quadrature) :: quad2
  type(gauss_lobatto_quadrature) :: quad3   

  ! set the polynomial type
  poly = legendre()


  print *, '===================================='
  print *, 'Gauss quadrature'
  
  call quad1%set(n,poly)   
  print *, quad1%points()
  print *, quad1%weights()
  
  print *, '===================================='

  print *, 'Gauss-Radau quadrature (left endpoint)'
  
  call quad2%set(n,poly)
  print *, quad2%points()
  print *, quad2%weights()

  print *, '===================================='

  print *, 'Gauss-Radau quadrature (right endpoint)'
  
  call quad2%set(n,poly,right=.true.)
  print *, quad2%points()
  print *, quad2%weights()

  print *, '===================================='

  print *, 'Gauss-Lobatto quadrature'
  
  call quad3%set(n,poly)
  print *, quad3%points()
  print *, quad3%weights()

  
end program quadrature_test
