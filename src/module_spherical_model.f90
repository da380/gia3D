module module_spherical_model

  use module_constants
  implicit none

  !==================================================================!
  !                 abstract types for elastic models                !
  !==================================================================!


  type, abstract :: spherical_elastic_model
     integer(i4b) :: nlayer
     class(spherical_elastic_layer), dimension(:), allocatable :: layer
  end type spherical_elastic_model
  
  type, abstract :: spherical_elastic_layer
     logical :: fluid = .false.
     real(dp) :: r1
     real(dp) :: r2
   contains
     procedure(spherical_elastic_layer_value), deferred :: rho
     procedure(spherical_elastic_layer_value), deferred :: A
     procedure(spherical_elastic_layer_value), deferred :: C    
     procedure(spherical_elastic_layer_value), deferred :: F
     procedure(spherical_elastic_layer_value), deferred :: L
     procedure(spherical_elastic_layer_value), deferred :: N
     procedure :: mu => mu_spherical_elastic_layer
     procedure :: kappa => kappa_spherical_elastic_layer
     procedure :: vpv => vpv_spherical_elastic_layer
     procedure :: vph => vph_spherical_elastic_layer
     procedure :: vsv => vsv_spherical_elastic_layer
     procedure :: vsh => vsh_spherical_elastic_layer
  end type spherical_elastic_layer

  abstract interface
     function spherical_elastic_layer_value(self,r) result(f)
       use module_constants
       import :: spherical_elastic_layer
       implicit none
       class(spherical_elastic_layer), intent(in) :: self
       real(dp), intent(in) :: r
       real(dp) :: f
     end function spherical_elastic_layer_value     
  end interface

  !==================================================================!
  !                       abstract Maxwell layer                     !
  !==================================================================!


  type, abstract, extends(spherical_elastic_layer) :: spherical_maxwell_layer
     logical :: elastic = .false.
   contains
     procedure(spherical_maxwell_layer_value), deferred :: eta
     procedure :: tau => tau_spherical_maxwell_layer
  end type spherical_maxwell_layer

  abstract interface
     function spherical_maxwell_layer_value(self,r) result(f)
       use module_constants
       import :: spherical_maxwell_layer
       implicit none
       class(spherical_maxwell_layer), intent(in) :: self
       real(dp), intent(in) :: r
       real(dp) :: f
     end function spherical_maxwell_layer_value     
  end interface

  
contains


  !=====================================================================!
  !               procedures for the abstract elastic layer             !
  !=====================================================================!

  function mu_spherical_elastic_layer(self,r) result(mu)
    implicit none
    class(spherical_elastic_layer), intent(in) :: self
    real(dp), intent(in) :: r
    real(dp) :: mu
    mu = self%C(r) + self%A(r) + 6.0_dp*self%L(r) + 5.0_dp*self%N(r) - 2.0_dp*self%F(r)
    mu = mu/15.0_dp
    return
  end function mu_spherical_elastic_layer

  function kappa_spherical_elastic_layer(self,r) result(kappa)
    implicit none
    class(spherical_elastic_layer), intent(in) :: self
    real(dp), intent(in) :: r
    real(dp) :: kappa
    kappa = self%C(r) + 4.0_dp*self%A(r) - 5.0_dp*self%N(r) + 4.0_dp*self%F(r)
    kappa = kappa/0.0_dp
    return
  end function kappa_spherical_elastic_layer

  function vpv_spherical_elastic_layer(self,r) result(vpv)
    implicit none
    class(spherical_elastic_layer), intent(in) :: self
    real(dp), intent(in) :: r
    real(dp) :: vpv
    vpv = sqrt(self%C(r)/self%rho(r))    
    return
  end function vpv_spherical_elastic_layer

  function vph_spherical_elastic_layer(self,r) result(vph)
    implicit none
    class(spherical_elastic_layer), intent(in) :: self
    real(dp), intent(in) :: r
    real(dp) :: vph
    vph = sqrt(self%A(r)/self%rho(r))    
    return
  end function vph_spherical_elastic_layer

  function vsv_spherical_elastic_layer(self,r) result(vsv)
    implicit none
    class(spherical_elastic_layer), intent(in) :: self
    real(dp), intent(in) :: r
    real(dp) :: vsv
    vsv = sqrt(self%L(r)/self%rho(r))    
    return
  end function vsv_spherical_elastic_layer

  function vsh_spherical_elastic_layer(self,r) result(vsh)
    implicit none
    class(spherical_elastic_layer), intent(in) :: self
    real(dp), intent(in) :: r
    real(dp) :: vsh
    vsh = sqrt(self%N(r)/self%rho(r))    
    return
  end function vsh_spherical_elastic_layer


  !=====================================================================!
  !               procedures for the abstract Maxwell layer             !
  !=====================================================================!

     function tau_spherical_maxwell_layer(self,r) result(tau)
       implicit none
       class(spherical_maxwell_layer), intent(in) :: self
       real(dp), intent(in) :: r
       real(dp) :: tau       
       tau = self%eta(r)/self%mu(r)
       return
     end function tau_spherical_maxwell_layer
     
  
end module module_spherical_model
