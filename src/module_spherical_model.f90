module module_spherical_model

  use module_constants
  implicit none


  type :: spherical_earth_model
     integer(i4b) :: nlayer
     class(*), dimension(:), allocatable :: layer
  end type spherical_earth_model
  
  type, abstract :: spherical_layer
     real(dp) :: r1
     real(dp) :: r2
  end type spherical_layer

  type, extends(spherical_layer) :: spherical_elastic_layer
   contains
     procedure :: setup
  end type spherical_elastic_layer
  

contains

  subroutine setup(self,r1,r2)
    implicit none
    class(spherical_elastic_layer), intent(inout) :: self
    real(dp), intent(in) :: r1,r2
    self%r1 = r1
    self%r2 = r2    
    return
  end subroutine setup
  
end module module_spherical_model
