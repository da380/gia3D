module module_spherical_model

  use module_constants  
  use module_physical_constants
  use module_util, only : poly_eval, count_columns
  use module_interp
  implicit none


  !==============================================================!
  !                          basic types                         !
  !==============================================================!

  type, abstract :: spherical_layer
     real(dp) :: r1
     real(dp) :: r2
   contains
     procedure(spherical_layer_function), deferred :: rho
  end type spherical_layer

  abstract interface
     function spherical_layer_function(self,r) result(f)
       use module_constants
       import :: spherical_layer
       class(spherical_layer), intent(in) :: self
       real(dp), intent(in) :: r
       real(dp) :: f
     end function spherical_layer_function
  end interface

  type, abstract, extends(spherical_layer) :: spherical_solid_elastic_layer
   contains
     procedure(spherical_solid_elastic_layer_function), deferred :: A
     procedure(spherical_solid_elastic_layer_function), deferred :: C
     procedure(spherical_solid_elastic_layer_function), deferred :: F
     procedure(spherical_solid_elastic_layer_function), deferred :: L
     procedure(spherical_solid_elastic_layer_function), deferred :: N
     procedure :: kappa => kappa_spherical_solid_elastic_layer
     procedure :: mu => mu_spherical_solid_elastic_layer
     procedure :: vpv => vpv_spherical_solid_elastic_layer
     procedure :: vph => vph_spherical_solid_elastic_layer
     procedure :: vsv => vsv_spherical_solid_elastic_layer
     procedure :: vsh => vsh_spherical_solid_elastic_layer
     procedure :: vp => vp_spherical_solid_elastic_layer
     procedure :: vs => vs_spherical_solid_elastic_layer
  end type spherical_solid_elastic_layer


  abstract interface
     function spherical_solid_elastic_layer_function(self,r) result(f)
       use module_constants
       import :: spherical_solid_elastic_layer
       class(spherical_solid_elastic_layer), intent(in) :: self
       real(dp), intent(in) :: r
       real(dp) :: f
     end function spherical_solid_elastic_layer_function
  end interface


  type, abstract, extends(spherical_layer) :: spherical_fluid_elastic_layer
   contains
     procedure(spherical_fluid_elastic_layer_function), deferred :: kappa
     procedure :: vp => vp_spherical_fluid_elastic_layer
  end type spherical_fluid_elastic_layer

  
  abstract interface
     function spherical_fluid_elastic_layer_function(self,r) result(f)
       use module_constants
       import :: spherical_fluid_elastic_layer
       class(spherical_fluid_elastic_layer), intent(in) :: self
       real(dp), intent(in) :: r
       real(dp) :: f
     end function spherical_fluid_elastic_layer_function
  end interface
  
  
  type, abstract, extends(spherical_solid_elastic_layer) :: spherical_maxwell_layer
   contains
     procedure(spherical_maxwell_layer_function), deferred :: eta
  end type spherical_maxwell_layer

  abstract interface
     function spherical_maxwell_layer_function(self,r) result(f)
       use module_constants
       import :: spherical_maxwell_layer
       class(spherical_maxwell_layer), intent(in) :: self
       real(dp), intent(in) :: r
       real(dp) :: f
     end function spherical_maxwell_layer_function
  end interface

  type spherical_section
     integer(i4b) :: nlayers
     real(dp) :: r1
     real(dp) :: r2
     class(spherical_layer), dimension(:), allocatable :: layer
  end type spherical_section


  type spherical_model
     real(dp) :: r1
     real(dp) :: r2
     integer(i4b) :: nsections = 0
     class(spherical_section), dimension(:), allocatable :: section
   contains
     procedure :: write_elastic => write_elastic_spherical_model
     procedure :: write_maxwell => write_maxwell_spherical_model
  end type spherical_model


  !==============================================================!
  !                     specific realisations                    !
  !==============================================================!

  
  type, extends(spherical_solid_elastic_layer) :: PREM_solid_elastic_layer
     integer(i4b) :: i
   contains
     procedure :: rho => rho_PREM_solid_elastic_layer
     procedure :: A => A_PREM_solid_elastic_layer
     procedure :: C => C_PREM_solid_elastic_layer
     procedure :: F => F_PREM_solid_elastic_layer
     procedure :: L => L_PREM_solid_elastic_layer
     procedure :: N => N_PREM_solid_elastic_layer
  end type PREM_solid_elastic_layer


  type, extends(spherical_fluid_elastic_layer) :: PREM_fluid_elastic_layer
     integer(i4b) :: i
   contains
     procedure :: rho => rho_PREM_fluid_elastic_layer
     procedure :: kappa => kappa_PREM_fluid_elastic_layer
  end type PREM_fluid_elastic_layer


  type, extends(spherical_maxwell_layer) :: PREM_maxwell_layer
     integer(i4b) :: i
     real(dp) :: tau_value
   contains
     procedure :: rho => rho_PREM_maxwell_layer
     procedure :: A => A_PREM_maxwell_layer
     procedure :: C => C_PREM_maxwell_layer
     procedure :: F => F_PREM_maxwell_layer
     procedure :: L => L_PREM_maxwell_layer
     procedure :: N => N_PREM_maxwell_layer
     procedure :: eta => eta_PREM_maxwell_layer
  end type PREM_maxwell_layer

  

  type, extends(spherical_solid_elastic_layer) :: DECK_solid_elastic_layer
     type(interp_1D_cubic) :: rho_cubic
     type(interp_1D_cubic) :: A_cubic
     type(interp_1D_cubic) :: C_cubic
     type(interp_1D_cubic) :: F_cubic
     type(interp_1D_cubic) :: L_cubic
     type(interp_1D_cubic) :: N_cubic
   contains
     procedure :: rho => rho_DECK_solid_elastic_layer
     procedure :: A => A_DECK_solid_elastic_layer
     procedure :: C => C_DECK_solid_elastic_layer
     procedure :: F => F_DECK_solid_elastic_layer
     procedure :: L => L_DECK_solid_elastic_layer
     procedure :: N => N_DECK_solid_elastic_layer
  end type DECK_solid_elastic_layer



  type, extends(spherical_fluid_elastic_layer) :: DECK_fluid_elastic_layer
     type(interp_1D_cubic) :: rho_cubic
     type(interp_1D_cubic) :: kappa_cubic
   contains
     procedure :: rho => rho_DECK_fluid_elastic_layer
     procedure :: kappa => kappa_DECK_fluid_elastic_layer
  end type DECK_fluid_elastic_layer


  type, extends(spherical_maxwell_layer) :: DECK_maxwell_layer
     type(interp_1D_cubic) :: rho_cubic
     type(interp_1D_cubic) :: A_cubic
     type(interp_1D_cubic) :: C_cubic
     type(interp_1D_cubic) :: F_cubic
     type(interp_1D_cubic) :: L_cubic
     type(interp_1D_cubic) :: N_cubic
     type(interp_1D_cubic) :: eta_cubic
   contains
     procedure :: rho => rho_DECK_maxwell_layer
     procedure :: A => A_DECK_maxwell_layer
     procedure :: C => C_DECK_maxwell_layer
     procedure :: F => F_DECK_maxwell_layer
     procedure :: L => L_DECK_maxwell_layer
     procedure :: N => N_DECK_maxwell_layer
     procedure :: eta => eta_DECK_maxwell_layer     
  end type DECK_maxwell_layer


  
  !==============================================================!
  !                          PREM data                           !
  !==============================================================!


  real(dp), dimension(2,13), parameter :: r_PREM = reshape((/   0.0e3_dp, 1221.5e3_dp, &
                                                             1221.5e3_dp, 3480.0e3_dp, & 
                                                             3480.0e3_dp, 3630.0e3_dp, &
                                                             3630.0e3_dp, 5600.0e3_dp, &
                                                             5600.0e3_dp, 5701.0e3_dp, &
                                                             5701.0e3_dp, 5771.0e3_dp, &
                                                             5771.0e3_dp, 5971.0e3_dp, & 
                                                             5971.0e3_dp, 6151.0e3_dp, & 
                                                             6151.0e3_dp, 6291.0e3_dp, &
                                                             6291.0e3_dp, 6346.6e3_dp, &
                                                             6346.6e3_dp, 6356.0e3_dp, &
                                                             6356.0e3_dp, 6368.0e3_dp, &
                                                             6368.0e3_dp, 6371.0e3_dp/)/length_norm,(/2,13/))
                                                             
  real(dp), dimension(4,13), parameter :: rho_coef_PREM = reshape((/ 13.08850e3_dp, 0.0e3_dp,-8.83810e3_dp, 0.0e3_dp,            &
                                                                12.58150e3_dp, -1.26380e3_dp, -3.64260e3_dp, -5.52810e3_dp, &
                                                                7.95650e3_dp, -6.47610e3_dp,  5.52830e3_dp, -3.08070e3_dp,  &
                                                                7.95650e3_dp, -6.47610e3_dp,  5.52830e3_dp, -3.08070e3_dp,  &
                                                                7.95650e3_dp, -6.47610e3_dp,  5.52830e3_dp, -3.08070e3_dp,  &
                                                                5.31970e3_dp, -1.48360e3_dp, 0.0e3_dp, 0.0e3_dp,            &
                                                                11.24940e3_dp, -8.02980e3_dp , 0.0e3_dp, 0.0e3_dp,          &
                                                                7.10890e3_dp, -3.80450e3_dp, 0.0e3_dp, 0.0e3_dp,            &
                                                                2.69100e3_dp,  0.69240e3_dp, 0.0e3_dp, 0.0e3_dp,            &
                                                                2.69100e3_dp,  0.69240e3_dp, 0.0e3_dp, 0.0e3_dp,            &
                                                                2.90000e3_dp, 0.0e3_dp, 0.0e3_dp, 0.0e3_dp,                 &
                                                                2.60000e3_dp, 0.0e3_dp, 0.0e3_dp, 0.0e3_dp,                 &
                                                                1.02000e3_dp, 0.0e3_dp, 0.0e3_dp, 0.0e3_dp/)/density_norm,  &
                                                                (/4,13/))

  real(dp), dimension(4,13), parameter :: vpv_coef_PREM = reshape((/ 11.26220e3_dp,  0.0e3_dp, -6.36400e3_dp, 0.0e3_dp,          &
                                                                11.04870e3_dp, -4.03620e3_dp,  4.8023e3_dp, -13.57320e3_dp, &
                                                                15.38910e3_dp, -5.31810e3_dp,  5.52420e3_dp, -2.55140e3_dp, &
                                                                24.9520e3_dp, -40.4673e3_dp,  51.4832e3_dp, -26.64190e3_dp, &
                                                                29.2766e3_dp, -23.6027e3_dp,   5.52420e3_dp, -2.55140e3_dp, &
                                                                19.09570e3_dp, -9.86720e3_dp, 0.0e3_dp, 0.0e3_dp,           &
                                                                39.7027e3_dp, -32.61660e3_dp, 0.0e3_dp, 0.0e3_dp,           &
                                                                20.3926e3_dp, -12.25690e3_dp, 0.0e3_dp, 0.0e3_dp,           &
                                                                0.83170e3_dp,  7.21800e3_dp, 0.0e3_dp, 0.0e3_dp,            &
                                                                0.83170e3_dp,  7.21800e3_dp, 0.0e3_dp, 0.0e3_dp,            &
                                                                6.80000e3_dp, 0.0e3_dp,  0.0e3_dp, 0.0e3_dp,                &
                                                                5.80000e3_dp, 0.0e3_dp,  0.0e3_dp, 0.0e3_dp,                &
                                                                1.45000e3_dp, 0.0e3_dp,  0.0e3_dp, 0.0e3_dp/)/velocity_norm,&
                                                                (/4,13/))
  
  
  real(dp), dimension(4,13), parameter :: vsv_coef_PREM = reshape((/ 3.66780e3_dp,  0.0e3_dp, -4.44750e3_dp, 0.0e3_dp,           &
                                                                0.0e3_dp,0.0e3_dp,0.0e3_dp,0.0e3_dp,                        &
                                                                6.92540e3_dp,  1.46720e3_dp, -2.08340e3_dp,  0.97830e3_dp,  &
                                                                11.1671e3_dp, -13.7818e3_dp,  17.4575e3_dp,  -9.2777e3_dp,  & 
                                                                22.3459e3_dp, -17.2473e3_dp,  -2.08340e3_dp,  0.97830e3_dp, & 
                                                                9.98390e3_dp, -4.93240e3_dp, 0.0e3_dp, 0.0e3_dp,            &
                                                                22.3512e3_dp, -18.58560e3_dp, 0.0e3_dp, 0.0e3_dp,           &
                                                                8.94960e3_dp, -4.45970e3_dp, 0.0e3_dp, 0.0e3_dp,            &
                                                                5.85820e3_dp, -1.46780e3_dp , 0.0e3_dp, 0.0e3_dp,           &
                                                                5.85820e3_dp, -1.46780e3_dp , 0.0e3_dp, 0.0e3_dp,           &
                                                                3.90000e3_dp, 0.0e3_dp , 0.0e3_dp, 0.0e3_dp,                &
                                                                3.20000e3_dp, 0.0e3_dp , 0.0e3_dp, 0.0e3_dp,                &
                                                                0.0e3_dp,0.0e3_dp,0.0e3_dp,0.0e3_dp/)/velocity_norm,(/4,13/))

  real(dp), dimension(4,13), parameter :: qk_coef_PREM  = reshape((/ 1327.7_dp, 0.0_dp, 0.0_dp, 0.0_dp,                  &
                                                                57823.0_dp, 0.0_dp, 0.0_dp, 0.0_dp,                 &
                                                                57823.0_dp, 0.0_dp, 0.0_dp, 0.0_dp,                 &
                                                                57823.0_dp, 0.0_dp, 0.0_dp, 0.0_dp,                 &
                                                                57823.0_dp, 0.0_dp, 0.0_dp, 0.0_dp,                 &
                                                                57823.0_dp, 0.0_dp, 0.0_dp, 0.0_dp,                 &
                                                                57823.0_dp, 0.0_dp, 0.0_dp, 0.0_dp,                 &
                                                                57823.0_dp, 0.0_dp, 0.0_dp, 0.0_dp,                 &
                                                                57823.0_dp, 0.0_dp, 0.0_dp, 0.0_dp,                 &
                                                                57823.0_dp, 0.0_dp, 0.0_dp, 0.0_dp,                 &
                                                                57823.0_dp, 0.0_dp, 0.0_dp, 0.0_dp,                 &
                                                                57823.0_dp, 0.0_dp, 0.0_dp, 0.0_dp,                 &
                                                                57823.0_dp, 0.0_dp, 0.0_dp, 0.0_dp/),(/4,13/))



  real(dp), dimension(4,13), parameter :: qm_coef_PREM  = reshape((/ 84.6_dp, 0.0_dp, 0.0_dp, 0.0_dp,                    &
                                                                0.0_dp,  0.0_dp, 0.0_dp, 0.0_dp,                    &
                                                                312.0_dp, 0.0_dp, 0.0_dp, 0.0_dp,                   &
                                                                312.0_dp, 0.0_dp, 0.0_dp, 0.0_dp,                   &
                                                                312.0_dp, 0.0_dp, 0.0_dp, 0.0_dp,                   &
                                                                143.0_dp, 0.0_dp, 0.0_dp, 0.0_dp,                   &
                                                                143.0_dp, 0.0_dp, 0.0_dp, 0.0_dp,                   &
                                                                143.0_dp, 0.0_dp, 0.0_dp, 0.0_dp,                   &
                                                                80.0_dp, 0.0_dp, 0.0_dp, 0.0_dp,                    &
                                                                600.0_dp, 0.0_dp, 0.0_dp, 0.0_dp,                   &
                                                                600.0_dp, 0.0_dp, 0.0_dp, 0.0_dp,                   &
                                                                600.0_dp, 0.0_dp, 0.0_dp, 0.0_dp,                   &
                                                                0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp/),(/4,13/))


  real(dp), dimension(4,13), parameter :: vph_coef_PREM = reshape((/ 11.26220e3_dp,  0.0e3_dp, -6.36400e3_dp, 0.0e3_dp,          &
                                                                11.04870e3_dp, -4.03620e3_dp,  4.8023e3_dp, -13.57320e3_dp, &
                                                                15.38910e3_dp, -5.31810e3_dp,  5.52420e3_dp, -2.55140e3_dp, &
                                                                24.9520e3_dp, -40.4673e3_dp,  51.4832e3_dp, -26.64190e3_dp, &
                                                                29.2766e3_dp, -23.6027e3_dp,   5.52420e3_dp, -2.55140e3_dp, &
                                                                19.09570e3_dp, -9.86720e3_dp, 0.0e3_dp, 0.0e3_dp,           &
                                                                39.7027e3_dp, -32.61660e3_dp, 0.0e3_dp, 0.0e3_dp,           &
                                                                20.3926e3_dp, -12.25690e3_dp, 0.0e3_dp, 0.0e3_dp,           &
                                                                3.59080e3_dp,  4.61720e3_dp, 0.0e3_dp, 0.0e3_dp,            &
                                                                3.59080e3_dp,  4.61720e3_dp, 0.0e3_dp, 0.0e3_dp,            &
                                                                6.80000e3_dp, 0.0e3_dp,  0.0e3_dp, 0.0e3_dp,                &
                                                                5.80000e3_dp, 0.0e3_dp,  0.0e3_dp, 0.0e3_dp,                &
                                                                1.45000e3_dp, 0.0e3_dp,  0.0e3_dp, 0.0e3_dp/)/velocity_norm,&
                                                                (/4,13/))

  real(dp), dimension(4,13), parameter :: vsh_coef_PREM = reshape((/ 3.66780e3_dp,  0.0e3_dp, -4.44750e3_dp, 0.0e3_dp,           &
                                                                0.0e3_dp,0.0e3_dp,0.0e3_dp,0.0e3_dp,                        &
                                                                6.92540e3_dp,  1.46720e3_dp, -2.08340e3_dp,  0.97830e3_dp,  &
                                                                11.1671e3_dp, -13.7818e3_dp,  17.4575e3_dp,  -9.2777e3_dp,  & 
                                                                22.3459e3_dp, -17.2473e3_dp,  -2.08340e3_dp,  0.97830e3_dp, & 
                                                                9.98390e3_dp, -4.93240e3_dp, 0.0e3_dp, 0.0e3_dp,            &
                                                                22.3512e3_dp, -18.58560e3_dp, 0.0e3_dp, 0.0e3_dp,           &
                                                                8.94960e3_dp, -4.45970e3_dp, 0.0e3_dp, 0.0e3_dp,            &
                                                                -1.08390e3_dp,  5.71760e3_dp , 0.0e3_dp, 0.0e3_dp,          &
                                                                -1.08390e3_dp,  5.71760e3_dp , 0.0e3_dp, 0.0e3_dp,          &
                                                                3.90000e3_dp, 0.0e3_dp , 0.0e3_dp, 0.0e3_dp,                &
                                                                3.20000e3_dp, 0.0e3_dp , 0.0e3_dp, 0.0e3_dp,                &
                                                                0.0e3_dp,0.0e3_dp,0.0e3_dp,0.0e3_dp/)/velocity_norm,(/4,13/))


  real(dp), dimension(4,13), parameter :: eta_coef_PREM = reshape((/ 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp,                     &
                                                                1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp,                     &
                                                                1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp,                     &
                                                                1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp,                     &
                                                                1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp,                     &
                                                                1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp,                     &
                                                                1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp,                     &
                                                                1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp,                     &
                                                                3.36870_dp, -2.47780_dp ,0.0_dp,0.0_dp,             &
                                                                3.36870_dp, -2.47780_dp ,0.0_dp,0.0_dp,             &
                                                                1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp,                     &
                                                                1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp,                     &
                                                                1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp/),(/4,13/))
  



contains


  
  !===========================================================================!
  !                        procedures for basic types                         !
  !===========================================================================!

  function kappa_spherical_solid_elastic_layer(self,r) result(kappa)
    class(spherical_solid_elastic_layer), intent(in) :: self
    real(dp), intent(in) :: r
    real(dp) :: kappa
    kappa = self%C(r) + 4.0_dp*self%A(r) - 4.0_dp*self%N(r) + 4.0_dp*self%F(r)
    kappa = kappa/9.0_dp
    return
  end function kappa_spherical_solid_elastic_layer


  function mu_spherical_solid_elastic_layer(self,r) result(mu)
    class(spherical_solid_elastic_layer), intent(in) :: self
    real(dp), intent(in) :: r
    real(dp) :: mu
    mu = self%C(r) + self%A(r) + 6.0_dp*self%L(r) + 5.0_dp*self%N(r) - 2.0_dp*self%F(r)
    mu = mu/15.0_dp
    return
  end function mu_spherical_solid_elastic_layer


  function vpv_spherical_solid_elastic_layer(self,r) result(vpv)
    class(spherical_solid_elastic_layer), intent(in) :: self
    real(dp), intent(in) :: r
    real(dp) :: vpv
    vpv = self%C(r)/self%rho(r)
    vpv = sqrt(vpv)
    return
  end function vpv_spherical_solid_elastic_layer


  function vph_spherical_solid_elastic_layer(self,r) result(vph)
    class(spherical_solid_elastic_layer), intent(in) :: self
    real(dp), intent(in) :: r
    real(dp) :: vph
    vph = self%A(r)/self%rho(r)
    vph = sqrt(vph)
    return
  end function vph_spherical_solid_elastic_layer

  
  function vsv_spherical_solid_elastic_layer(self,r) result(vsv)
    class(spherical_solid_elastic_layer), intent(in) :: self
    real(dp), intent(in) :: r
    real(dp) :: vsv
    vsv = self%L(r)/self%rho(r)
    vsv = sqrt(vsv)
    return
  end function vsv_spherical_solid_elastic_layer


  function vsh_spherical_solid_elastic_layer(self,r) result(vsh)
    class(spherical_solid_elastic_layer), intent(in) :: self
    real(dp), intent(in) :: r
    real(dp) :: vsh
    vsh = self%N(r)/self%rho(r)
    vsh = sqrt(vsh)
    return
  end function vsh_spherical_solid_elastic_layer


  function vp_spherical_solid_elastic_layer(self,r) result(vp)
    class(spherical_solid_elastic_layer), intent(in) :: self
    real(dp), intent(in) :: r
    real(dp) :: vp
    vp = (self%kappa(r)+4.0_dp*self%mu(r)/3.0_dp)/self%rho(r)
    vp = sqrt(vp)
    return
  end function vp_spherical_solid_elastic_layer


  function vs_spherical_solid_elastic_layer(self,r) result(vs)
    class(spherical_solid_elastic_layer), intent(in) :: self
    real(dp), intent(in) :: r
    real(dp) :: vs
    vs = self%mu(r)/self%rho(r)
    vs = sqrt(vs)
    return
  end function vs_spherical_solid_elastic_layer


  function vp_spherical_fluid_elastic_layer(self,r) result(vp)
    class(spherical_fluid_elastic_layer), intent(in) :: self
    real(dp), intent(in) :: r
    real(dp) :: vp
    vp = self%kappa(r)/self%rho(r)
    vp = sqrt(vp)
    return
  end function vp_spherical_fluid_elastic_layer


  
    
  subroutine write_elastic_spherical_model(self,model_file,nknot,isotropic)
    class(spherical_model), intent(in) :: self
    character(len=*), intent(in) :: model_file
    integer(i4b), intent(in) :: nknot
    logical, intent(in), optional :: isotropic

    logical :: iso
    integer(i4b) :: io,isection,ilayer,ir,nr
    real(dp) :: r1,r2,r,dr,eta,rho,vpv,vph,vsv,vsh

    iso = .false.
    if(present(isotropic)) then
       if(isotropic) iso = isotropic
    end if
          
    
    r1 = self%r1
    r2 = self%r2
    open(newunit = io,file=trim(model_file))
    do isection = 1,self%nsections

       do ilayer = 1,self%section(isection)%nlayers

          associate(layer => self%section(isection)%layer(ilayer))

            nr = nknot*(layer%r2-layer%r1)/(r2-r1)
            nr = max(2,nr)
            dr = (layer%r2-layer%r1)/(nr-1)


            do ir = 1,nr
               
               r = layer%r1 + (ir-1)*dr
               rho = layer%rho(r)
               
               select type(layer)
                  
               class is(spherical_solid_elastic_layer)
               
                  if(iso) then
                     vpv = layer%vp(r)
                     vph = layer%vp(r)
                     vsv = layer%vs(r)
                     vsh = layer%vs(r)
                     eta = 1.0_dp
                  else
                     vpv = layer%vpv(r)
                     vph = layer%vph(r)
                     vsv = layer%vsv(r)
                     vsh = layer%vsh(r)
                     eta = layer%F(r)/(layer%A(r)-2.0_dp*layer%L(r))
                  end if

               class is(spherical_fluid_elastic_layer)

                  vpv = layer%vp(r)
                  vph = layer%vp(r)
                  vsv = 0.0_dp
                  vsh = 0.0_dp
                  eta = 1.0_dp
                  
               class default
                  stop 'write_elastic_spherical_model: layer of an invalid type'
               end select

               
               if(iso) then
                  
                  write(io,'(4e18.8)') r*length_norm,      &
                                       rho*density_norm,   & 
                                       vpv*velocity_norm,  & 
                                       vsv*velocity_norm
                  
               else
                  
                  write(io,'(7e18.8)') r*length_norm,      &
                                       rho*density_norm,   & 
                                       vpv*velocity_norm,  & 
                                       vsv*velocity_norm,  & 
                                       vph*velocity_norm,  & 
                                       vsh*velocity_norm,  &
                                       eta
                  
               end if
            end do

          end associate

          
       end do
       
    end do
    close(io)
    
    return
  end subroutine write_elastic_spherical_model



  subroutine write_maxwell_spherical_model(self,model_file,nknot,isotropic)
    class(spherical_model), intent(in) :: self
    character(len=*), intent(in) :: model_file
    integer(i4b), intent(in) :: nknot
    logical, intent(in), optional :: isotropic

    logical :: iso
    integer(i4b) :: io,isection,ilayer,ir,nr
    real(dp) :: r1,r2,r,dr,eta,rho,vpv,vph,vsv,vsh,visco

    iso = .false.
    if(present(isotropic)) then
       if(isotropic) iso = isotropic
    end if
          
    
    r1 = self%r1
    r2 = self%r2
    open(newunit = io,file=trim(model_file))
    do isection = 1,self%nsections

       do ilayer = 1,self%section(isection)%nlayers

          associate(layer => self%section(isection)%layer(ilayer))

            nr = nknot*(layer%r2-layer%r1)/(r2-r1)
            nr = max(2,nr)
            dr = (layer%r2-layer%r1)/(nr-1)


            do ir = 1,nr
               
               r = layer%r1 + (ir-1)*dr
               rho = layer%rho(r)
               
               select type(layer)
                  
               class is(spherical_solid_elastic_layer)
               
                  if(iso) then
                     vpv = layer%vp(r)
                     vph = layer%vp(r)
                     vsv = layer%vs(r)
                     vsh = layer%vs(r)
                     eta = 1.0_dp
                     visco = 0.0_dp
                  else
                     vpv = layer%vpv(r)
                     vph = layer%vph(r)
                     vsv = layer%vsv(r)
                     vsh = layer%vsh(r)
                     eta = layer%F(r)/(layer%A(r)-2.0_dp*layer%L(r))
                     visco = 0.0_dp
                  end if
                  
               class is(spherical_fluid_elastic_layer)

                  vpv = layer%vp(r)
                  vph = layer%vp(r)
                  vsv = 0.0_dp
                  vsh = 0.0_dp
                  eta = 1.0_dp
                  visco = 0.0_dp

               class is(spherical_maxwell_layer)

                  if(iso) then
                     vpv = layer%vp(r)
                     vph = layer%vp(r)
                     vsv = layer%vs(r)
                     vsh = layer%vs(r)
                     eta = 1.0_dp
                     visco = log10(layer%eta(r)*viscosity_norm)
                  else
                     vpv = layer%vpv(r)
                     vph = layer%vph(r)
                     vsv = layer%vsv(r)
                     vsh = layer%vsh(r)
                     eta = layer%F(r)/(layer%A(r)-2.0_dp*layer%L(r))
                     visco = log10(layer%eta(r)*viscosity_norm)
                  end if                  
                  
               class default
                  stop 'write_maxwell_spherical_model: layer of an invalid type'
               end select

               
               
               if(iso) then
                  
                  write(io,'(5e18.8)') r*length_norm,      &
                                       rho*density_norm,   & 
                                       vpv*velocity_norm,  & 
                                       vsv*velocity_norm,  &
                                       visco
                  
               else
                  
                  write(io,'(8e18.8)') r*length_norm,      &
                                       rho*density_norm,   & 
                                       vpv*velocity_norm,  & 
                                       vsv*velocity_norm,  & 
                                       vph*velocity_norm,  & 
                                       vsh*velocity_norm,  &
                                       eta,                &
                                       visco
                  
               end if
            end do

          end associate

          
       end do
       
    end do
    close(io)
    
    return
  end subroutine write_maxwell_spherical_model

  
  
  !=============================================================!
  !                         PREM routines                       !
  !=============================================================!


  function elastic_PREM(ocean) result(model)
    logical, intent(in), optional :: ocean
    type(spherical_model) :: model

    logical :: ocean_local
    integer(i4b) :: ilayer

    if(present(ocean)) then
       ocean_local = ocean
    else
       ocean_local = .true.
    end if

    
    ! set the section structure
    model%r1 = r_PREM(1,1)
    if(ocean_local) then
       model%r2 = r_PREM(2,13)
       model%nsections = 4
    else
       model%r2 = r_PREM(2,12)
       model%nsections = 3
    end if
    allocate(model%section(model%nsections))

    
    ! build the inner core
    model%section(1)%nlayers = 1
    model%section(1)%r1 = r_PREM(1,1)
    model%section(1)%r2 = r_PREM(2,1)
    allocate(PREM_solid_elastic_layer :: model%section(1)%layer(1))
    associate(layer => model%section(1)%layer(1))
      select type(layer)
      type is(PREM_solid_elastic_layer)
         layer%i = 1
         layer%r1 = r_PREM(1,layer%i)
         layer%r2 = r_PREM(2,layer%i)
      end select
    end associate

    ! build the outer core
    model%section(2)%nlayers = 1
    model%section(2)%r1 = r_PREM(1,2)
    model%section(2)%r2 = r_PREM(2,2)
    allocate(PREM_fluid_elastic_layer :: model%section(2)%layer(1))
    associate(layer => model%section(2)%layer(1))
      select type(layer)
      type is(PREM_fluid_elastic_layer)
         layer%i = 2
         layer%r1 = r_PREM(1,layer%i)
         layer%r2 = r_PREM(2,layer%i)
      end select
    end associate

    ! build the mantle
    model%section(3)%nlayers = 10
    model%section(3)%r1 = r_PREM(1,3)
    model%section(3)%r2 = r_PREM(2,12)
    allocate(PREM_solid_elastic_layer :: model%section(3)%layer(model%section(3)%nlayers))
    do ilayer = 1,model%section(3)%nlayers
       associate(layer => model%section(3)%layer(ilayer))
         select type(layer)
         type is(PREM_solid_elastic_layer)
            layer%i = ilayer + 2
            layer%r1 = r_PREM(1,layer%i)
            layer%r2 = r_PREM(2,layer%i)
         end select
       end associate
    end do

    if(ocean_local) then
       ! build the ocean
       model%section(4)%nlayers = 1
       model%section(4)%r1 = r_PREM(1,13)
       model%section(4)%r2 = r_PREM(2,13)
       allocate(PREM_fluid_elastic_layer :: model%section(4)%layer(1))
       associate(layer => model%section(4)%layer(1))
         select type(layer)
         type is(PREM_fluid_elastic_layer)
            layer%i = 13
            layer%r1 = r_PREM(1,layer%i)
            layer%r2 = r_PREM(2,layer%i)
         end select
       end associate       
    end if
    
    return
  end function elastic_PREM

  function rho_PREM_solid_elastic_layer(self,r) result(rho)
    class(PREM_solid_elastic_layer), intent(in) :: self
    real(dp), intent(in) :: r
    real(dp) :: rho
    rho = poly_eval(4,rho_coef_PREM(:,self%i),r/r_PREM(2,13))
    return
  end function rho_PREM_solid_elastic_layer

  function A_PREM_solid_elastic_layer(self,r) result(A)
    class(PREM_solid_elastic_layer), intent(in) :: self
    real(dp), intent(in) :: r
    real(dp) :: A
    A = self%rho(r)*poly_eval(4,vpv_coef_PREM(:,self%i),r/r_PREM(2,13))**2
    return
  end function A_PREM_solid_elastic_layer
  
  function C_PREM_solid_elastic_layer(self,r) result(C)
    class(PREM_solid_elastic_layer), intent(in) :: self
    real(dp), intent(in) :: r
    real(dp) :: C
    C = self%rho(r)*poly_eval(4,vpv_coef_PREM(:,self%i),r/r_PREM(2,13))**2
    return
  end function C_PREM_solid_elastic_layer

  function F_PREM_solid_elastic_layer(self,r) result(F)
    class(PREM_solid_elastic_layer), intent(in) :: self
    real(dp), intent(in) :: r
    real(dp) :: F
    F = poly_eval(4,eta_coef_PREM(:,self%i),r/r_PREM(2,13))*(self%A(r) - 2.0_dp*self%L(r))
    return
  end function F_PREM_solid_elastic_layer

  function L_PREM_solid_elastic_layer(self,r) result(L)
    class(PREM_solid_elastic_layer), intent(in) :: self
    real(dp), intent(in) :: r
    real(dp) :: L
    L = self%rho(r)*poly_eval(4,vsv_coef_PREM(:,self%i),r/r_PREM(2,13))**2
    return
  end function L_PREM_solid_elastic_layer

  function N_PREM_solid_elastic_layer(self,r) result(N)
    class(PREM_solid_elastic_layer), intent(in) :: self
    real(dp), intent(in) :: r
    real(dp) :: N
    N = self%rho(r)*poly_eval(4,vsh_coef_PREM(:,self%i),r/r_PREM(2,13))**2
    return
  end function N_PREM_solid_elastic_layer


  function rho_PREM_fluid_elastic_layer(self,r) result(rho)
    class(PREM_fluid_elastic_layer), intent(in) :: self
    real(dp), intent(in) :: r
    real(dp) :: rho
    rho = poly_eval(4,rho_coef_PREM(:,self%i),r/r_PREM(2,13))
    return
  end function rho_PREM_fluid_elastic_layer

  
  function kappa_PREM_fluid_elastic_layer(self,r) result(kappa)
    class(PREM_fluid_elastic_layer), intent(in) :: self
    real(dp), intent(in) :: r
    real(dp) :: kappa
    kappa = self%rho(r)*poly_eval(4,vpv_coef_PREM(:,self%i),r/r_PREM(2,13))**2
    return
  end function kappa_PREM_fluid_elastic_layer



  function maxwell_PREM(tau_LM,tau_UM,tau_C) result(model)
    real(dp), intent(in) :: tau_LM,tau_UM
    real(dp), intent(in), optional :: tau_C
    type(spherical_model) :: model

    logical :: elastic_crust
    integer(i4b) :: ilayer

    ! check for optional argument
    elastic_crust = .not.present(tau_C)
    
    ! set the section structure
    model%r1 = r_PREM(1,1)
    model%r2 = r_PREM(2,12)
    if(elastic_crust) then
       model%nsections = 4
    else
       model%nsections = 3
    end if
    allocate(model%section(model%nsections))

    
    ! build the inner core
    model%section(1)%nlayers = 1
    model%section(1)%r1 = r_PREM(1,1)
    model%section(1)%r2 = r_PREM(2,1)
    allocate(PREM_solid_elastic_layer :: model%section(1)%layer(1))
    associate(layer => model%section(1)%layer(1))
      select type(layer)
      type is(PREM_solid_elastic_layer)
         layer%i = 1
         layer%r1 = r_PREM(1,layer%i)
         layer%r2 = r_PREM(2,layer%i)
      end select
    end associate

    
    ! build the outer core
    model%section(2)%nlayers = 1
    model%section(2)%r1 = r_PREM(1,2)
    model%section(2)%r2 = r_PREM(2,2)
    allocate(PREM_fluid_elastic_layer :: model%section(2)%layer(1))
    associate(layer => model%section(2)%layer(1))
      select type(layer)
      type is(PREM_fluid_elastic_layer)
         layer%i = 2
         layer%r1 = r_PREM(1,layer%i)
         layer%r2 = r_PREM(2,layer%i)
      end select
    end associate

    ! build the mantle
    if(elastic_crust) then
       model%section(3)%nlayers = 7
       model%section(3)%r1 = r_PREM(1,3)
       model%section(3)%r2 = r_PREM(2,9)
    else
       model%section(3)%nlayers = 10
       model%section(3)%r1 = r_PREM(1,3)
       model%section(3)%r2 = r_PREM(2,12)
    end if       
    allocate(PREM_maxwell_layer :: model%section(3)%layer(model%section(3)%nlayers))
    do ilayer = 1,model%section(3)%nlayers
       associate(layer => model%section(3)%layer(ilayer))
         select type(layer)
         type is(PREM_maxwell_layer)
            layer%i = ilayer + 2
            layer%r1 = r_PREM(1,layer%i)
            layer%r2 = r_PREM(2,layer%i)
            if(ilayer <= 3) then
               layer%tau_value = tau_LM
            else
               layer%tau_value = tau_UM
            end if
            if(.not.elastic_crust .and. ilayer > 7) then
               layer%tau_value = tau_C
            end if
         end select
       end associate
    end do



    if(elastic_crust) then
       model%section(4)%nlayers = 3
       model%section(4)%r1 = r_PREM(1,9)
       model%section(4)%r2 = r_PREM(2,12)
       allocate(PREM_solid_elastic_layer :: model%section(4)%layer(model%section(4)%nlayers))
       do ilayer = 1,model%section(4)%nlayers
          associate(layer => model%section(4)%layer(ilayer))
            select type(layer)
            type is(PREM_solid_elastic_layer)
               layer%i = ilayer + 9
               layer%r1 = r_PREM(1,layer%i)
               layer%r2 = r_PREM(2,layer%i)
            end select
          end associate
       end do
    end if


    
    return
  end function maxwell_PREM

  

  function rho_PREM_maxwell_layer(self,r) result(rho)
    class(PREM_maxwell_layer), intent(in) :: self
    real(dp), intent(in) :: r
    real(dp) :: rho
    rho = poly_eval(4,rho_coef_PREM(:,self%i),r/r_PREM(2,13))
    return
  end function rho_PREM_maxwell_layer

  function A_PREM_maxwell_layer(self,r) result(A)
    class(PREM_maxwell_layer), intent(in) :: self
    real(dp), intent(in) :: r
    real(dp) :: A
    A = self%rho(r)*poly_eval(4,vpv_coef_PREM(:,self%i),r/r_PREM(2,13))**2
    return
  end function A_PREM_maxwell_layer
  
  function C_PREM_maxwell_layer(self,r) result(C)
    class(PREM_maxwell_layer), intent(in) :: self
    real(dp), intent(in) :: r
    real(dp) :: C
    C = self%rho(r)*poly_eval(4,vpv_coef_PREM(:,self%i),r/r_PREM(2,13))**2
    return
  end function C_PREM_maxwell_layer

  function F_PREM_maxwell_layer(self,r) result(F)
    class(PREM_maxwell_layer), intent(in) :: self
    real(dp), intent(in) :: r
    real(dp) :: F
    F = poly_eval(4,eta_coef_PREM(:,self%i),r/r_PREM(2,13))*(self%A(r) - 2.0_dp*self%L(r))
    return
  end function F_PREM_maxwell_layer

  function L_PREM_maxwell_layer(self,r) result(L)
    class(PREM_maxwell_layer), intent(in) :: self
    real(dp), intent(in) :: r
    real(dp) :: L
    L = self%rho(r)*poly_eval(4,vsv_coef_PREM(:,self%i),r/r_PREM(2,13))**2
    return
  end function L_PREM_maxwell_layer

  function N_PREM_maxwell_layer(self,r) result(N)
    class(PREM_maxwell_layer), intent(in) :: self
    real(dp), intent(in) :: r
    real(dp) :: N
    N = self%rho(r)*poly_eval(4,vsh_coef_PREM(:,self%i),r/r_PREM(2,13))**2
    return
  end function N_PREM_maxwell_layer


  function eta_PREM_maxwell_layer(self,r) result(eta)
    class(PREM_maxwell_layer), intent(in) :: self
    real(dp), intent(in) :: r
    real(dp) :: eta
    eta = self%tau_value*self%mu(r)
    return
  end function eta_PREM_maxwell_layer


  
  !=============================================================!
  !                         DECK routines                       !
  !=============================================================!


  function elastic_DECK(file) result(model)
    character(len=*), intent(in) :: file
    type(spherical_model) :: model

    logical :: isotropic,change
    integer(i4b) :: io,iknot,nknot,ios,ncol,isec,i1,i2,ilay,i11,i22
    integer(i4b), dimension(:,:), allocatable :: sind
    real(dp) :: tmp
    real(dp), dimension(:), allocatable :: r,rho,vpv,vsv,vph,vsh,eta

    !------------------------------------!
    !           read in the file         !
    !------------------------------------!
    open(newunit = io,file = trim(file),iostat = ios)
    if(ios /= 0) stop 'elastic_DECK: problem opening model file'
    ncol = count_columns(io)
    if(ncol == 4) then
       isotropic = .true.
    else if(ncol == 7) then
       isotropic = .false.
    else
       stop 'elastic_DECK: unexpected format'
    end if
    nknot = 0
    do
       read(io,*,iostat = ios) tmp
       if(ios /= 0) exit
       nknot = nknot + 1
    end do
    allocate(r(nknot),rho(nknot),vpv(nknot),vsv(nknot),vph(nknot),vsh(nknot),eta(nknot))
    rewind(io)    
    do iknot = 1,nknot
       if(isotropic) then
          read(io,*) r(iknot),rho(iknot),vpv(iknot),vsv(iknot)
          vph(iknot) = vpv(iknot)
          vsh(iknot) = vsv(iknot)
          eta(iknot) = 1.0_dp
       else
          read(io,*) r(iknot),rho(iknot),vpv(iknot),vsv(iknot),vph(iknot),vsh(iknot),eta(iknot)
       end if
    end do
    close(io)

    ! non-dimensionalise
    r   = r/length_norm
    rho = rho/density_norm
    vpv = vpv/velocity_norm
    vsv = vsv/velocity_norm
    vph = vph/velocity_norm
    vsh = vsh/velocity_norm

    ! convert to moduli
    vpv = rho*vpv*vpv  ! C
    vph = rho*vph*vph  ! A
    vsv = rho*vsv*vsv  ! L
    vsh = rho*vsh*vsh  ! N
    eta = eta*(vph-2.0_dp*vsv) ! F
    
    ! set the radii limits
    model%r1 = r(1)
    model%r2 = r(nknot)
    
    ! work out the number of sections
    model%nsections = 1
    do iknot = 2,nknot
       change = (vsv(iknot-1) == 0.0_dp) .neqv. (vsv(iknot) == 0.0_dp)
       if(change) model%nsections = model%nsections + 1
    end do
    allocate(model%section(model%nsections))
    

    ! locally index the sections
    allocate(sind(3,model%nsections))
    isec = 1
    i1 = 1
    do iknot = 2,nknot
       change = (vsv(iknot-1) == 0.0_dp) .neqv. (vsv(iknot) == 0.0_dp)
       if(change) then
          i2 = iknot-1
          sind(1,isec) = i1
          sind(2,isec) = i2
          if(vsv(i1) == 0.0_dp) then
             sind(3,isec) = 1
          else
             sind(3,isec) = 0
          end if
          i1 = i2+1
          isec = isec + 1
       end if
    end do
    i2 = nknot
    sind(1,isec) = i1
    sind(2,isec) = i2
    if(vsv(i1) == 0.0_dp) then
       sind(3,isec) = 1
    else
       sind(3,isec) = 0
    end if

    ! build up the layers for each section
    do isec = 1,model%nsections

       i1 = sind(1,isec)
       i2 = sind(2,isec)
       associate(nlayers => model%section(isec)%nlayers, section => model%section(isec))

         ! store the section radii
         section%r1 = r(i1)
         section%r2 = r(i2)
         
         nlayers = 1
         do iknot = i1+1,i2
            if(r(iknot-1) == r(iknot)) nlayers = nlayers + 1
         end do


         ! allocate the layers
         if(sind(3,isec) == 0) then
            allocate(DECK_solid_elastic_layer:: section%layer(nlayers))
         else
            allocate(DECK_fluid_elastic_layer:: section%layer(nlayers))
         end if

         ilay = 1
         section%layer(ilay)%r1 = section%r1
         i11 = i1
         do iknot = i1+1,i2
            if(r(iknot-1) == r(iknot)) then
               section%layer(ilay)%r2 = r(iknot)
               i22 = iknot-1
               associate(layer => section%layer(ilay))
                 layer%r1 = r(i11)
                 layer%r2 = r(i22)
                 select type(layer)
                 type is(DECK_solid_elastic_layer)
                    call layer%rho_cubic%set(r(i11:i22),rho(i11:i22))
                    call layer%A_cubic%set(r(i11:i22),vph(i11:i22))
                    call layer%C_cubic%set(r(i11:i22),vpv(i11:i22))
                    call layer%F_cubic%set(r(i11:i22),eta(i11:i22))
                    call layer%L_cubic%set(r(i11:i22),vsv(i11:i22))
                    call layer%N_cubic%set(r(i11:i22),vsh(i11:i22))
                 type is(DECK_fluid_elastic_layer)                    
                    call layer%rho_cubic%set(r(i11:i22),rho(i11:i22))
                    call layer%kappa_cubic%set(r(i11:i22),vpv(i11:i22))
                 end select
                 i11 = i22+1
               end associate
               ilay = ilay + 1
            end if
         end do
         section%layer(ilay)%r2 = section%r2
         i22 = i2
         associate(layer => section%layer(ilay))
           layer%r1 = r(i11)
           layer%r2 = r(i22)
           select type(layer)
           type is(DECK_solid_elastic_layer)
              call layer%rho_cubic%set(r(i11:i22),rho(i11:i22))
              call layer%A_cubic%set(r(i11:i22),vph(i11:i22))
              call layer%C_cubic%set(r(i11:i22),vpv(i11:i22))
              call layer%F_cubic%set(r(i11:i22),eta(i11:i22))
              call layer%L_cubic%set(r(i11:i22),vsv(i11:i22))
              call layer%N_cubic%set(r(i11:i22),vsh(i11:i22))
           type is(DECK_fluid_elastic_layer)                    
              call layer%rho_cubic%set(r(i11:i22),rho(i11:i22))
              call layer%kappa_cubic%set(r(i11:i22),vpv(i11:i22))
           end select
         end associate
         
         
       end associate
       

    end do
    
    
    return
  end function elastic_DECK
  

  function rho_DECK_solid_elastic_layer(self,r) result(rho)
    class(DECK_solid_elastic_layer), intent(in) :: self
    real(dp), intent(in) :: r
    real(dp) :: rho
    rho = self%rho_cubic%val(r)
    return
  end function rho_DECK_solid_elastic_layer

  function A_DECK_solid_elastic_layer(self,r) result(A)
    class(DECK_solid_elastic_layer), intent(in) :: self
    real(dp), intent(in) :: r
    real(dp) :: A,vpv    
    A = self%A_cubic%val(r)
    return
  end function A_DECK_solid_elastic_layer
  
  function C_DECK_solid_elastic_layer(self,r) result(C)
    class(DECK_solid_elastic_layer), intent(in) :: self
    real(dp), intent(in) :: r
    real(dp) :: C,rho,vsv
    C = self%C_cubic%val(r)
    return
  end function C_DECK_solid_elastic_layer

  function F_DECK_solid_elastic_layer(self,r) result(F)
    class(DECK_solid_elastic_layer), intent(in) :: self
    real(dp), intent(in) :: r
    real(dp) :: F
    F = self%F_cubic%val(r)
    return
  end function F_DECK_solid_elastic_layer

  function L_DECK_solid_elastic_layer(self,r) result(L)
    class(DECK_solid_elastic_layer), intent(in) :: self
    real(dp), intent(in) :: r
    real(dp) :: L
    L = self%L_cubic%val(r)
    return
  end function L_DECK_solid_elastic_layer

  function N_DECK_solid_elastic_layer(self,r) result(N)
    class(DECK_solid_elastic_layer), intent(in) :: self
    real(dp), intent(in) :: r
    real(dp) :: N
    N = self%N_cubic%val(r)
    return
  end function N_DECK_solid_elastic_layer

  
  function rho_DECK_fluid_elastic_layer(self,r) result(rho)
    class(DECK_fluid_elastic_layer), intent(in) :: self
    real(dp), intent(in) :: r
    real(dp) :: rho
    rho = self%rho_cubic%val(r)
    return
  end function rho_DECK_fluid_elastic_layer

  function kappa_DECK_fluid_elastic_layer(self,r) result(kappa)
    class(DECK_fluid_elastic_layer), intent(in) :: self
    real(dp), intent(in) :: r
    real(dp) :: kappa,vpv    
    kappa = self%kappa_cubic%val(r)
    return
  end function kappa_DECK_fluid_elastic_layer

  
  function maxwell_DECK(file) result(model)
    character(len=*), intent(in) :: file
    type(spherical_model) :: model

    logical :: isotropic,change,viscous,fluid
    integer(i4b) :: io,iknot,nknot,ios,ncol,isec,i1,i2,ilay,i11,i22
    integer(i4b), dimension(:,:), allocatable :: sind
    real(dp) :: tmp
    real(dp), dimension(:), allocatable :: r,rho,vpv,vsv,vph,vsh,eta,visco

    !------------------------------------!
    !           read in the file         !
    !------------------------------------!
    open(newunit = io,file = trim(file),iostat = ios)
    if(ios /= 0) stop 'maxwell_DECK: problem opening model file'
    ncol = count_columns(io)
    if(ncol == 5) then
       isotropic = .true.
    else if(ncol == 8) then
       isotropic = .false.
    else
       stop 'maxwell_DECK: unexpected format'
    end if
    nknot = 0
    do
       read(io,*,iostat = ios) tmp
       if(ios /= 0) exit
       nknot = nknot + 1
    end do
    allocate(r(nknot),rho(nknot),vpv(nknot),vsv(nknot),vph(nknot),vsh(nknot),eta(nknot),visco(nknot))
    rewind(io)    
    do iknot = 1,nknot
       if(isotropic) then
          read(io,*) r(iknot),rho(iknot),vpv(iknot),vsv(iknot),visco(iknot)
          vph(iknot) = vpv(iknot)
          vsh(iknot) = vsv(iknot)
          eta(iknot) = 1.0_dp
       else
          read(io,*) r(iknot),rho(iknot),vpv(iknot),vsv(iknot),vph(iknot),vsh(iknot),eta(iknot),visco(iknot)
       end if
    end do
    close(io)

    ! non-dimensionalise
    r   = r/length_norm
    rho = rho/density_norm
    vpv = vpv/velocity_norm
    vsv = vsv/velocity_norm
    vph = vph/velocity_norm
    vsh = vsh/velocity_norm
    where(visco /= 0.0_dp) visco = (10**visco)/viscosity_norm

    
    ! convert to moduli
    vpv = rho*vpv*vpv  ! C
    vph = rho*vph*vph  ! A
    vsv = rho*vsv*vsv  ! L
    vsh = rho*vsh*vsh  ! N
    eta = eta*(vph-2.0_dp*vsv) ! F
    
    ! set the radii limits
    model%r1 = r(1)
    model%r2 = r(nknot)

    
    ! work out the number of sections
    model%nsections = 1
    fluid = (vsv(1) == 0.0_dp)
    viscous = (visco(1) /= 0.0_dp)
    do iknot = 2,nknot
       change = ((vsv(iknot) == 0.0_dp) .neqv. fluid) .or. ((visco(iknot) /= 0.0_dp) .neqv. viscous)
       if(change) then
          fluid = (vsv(iknot) == 0.0_dp)
          viscous = (visco(iknot) /= 0.0_dp)
          model%nsections = model%nsections + 1
       end if
    end do
    allocate(model%section(model%nsections))


    
    ! locally index the sections
    allocate(sind(3,model%nsections))
    isec = 1
    i1 = 1
    fluid = (vsv(1) == 0.0_dp)
    viscous = (visco(1) /= 0.0_dp)
    do iknot = 2,nknot
       change = ((vsv(iknot) == 0.0_dp) .neqv. fluid) .or. ((visco(iknot) /= 0.0_dp) .neqv. viscous)
       if(change) then
          i2 = iknot-1
          sind(1,isec) = i1
          sind(2,isec) = i2
          if(fluid) then
             sind(3,isec) = 1
          else if(viscous) then
             sind(3,isec) = 2
          else
             sind(3,isec) = 0
          end if
          i1 = i2+1
          isec = isec + 1
          fluid = (vsv(iknot) == 0.0_dp)
          viscous = (visco(iknot) /= 0.0_dp)
       end if
    end do
    i2 = nknot
    sind(1,isec) = i1
    sind(2,isec) = i2
    fluid = (vsv(nknot) == 0.0_dp)
    viscous = (visco(nknot) /= 0.0_dp)
    if(fluid) then
       sind(3,isec) = 1
    else if(viscous)  then
       sind(3,isec) = 2
    else
       sind(3,isec) = 0
    end if

    
    ! build up the layers for each section
    do isec = 1,model%nsections

       i1 = sind(1,isec)
       i2 = sind(2,isec)
       associate(nlayers => model%section(isec)%nlayers, section => model%section(isec))

         ! store the section radii
         section%r1 = r(i1)
         section%r2 = r(i2)
         
         nlayers = 1
         do iknot = i1+1,i2
            if(r(iknot-1) == r(iknot)) nlayers = nlayers + 1
         end do


         ! allocate the layers
         if(sind(3,isec) == 0) then
            allocate(DECK_solid_elastic_layer:: section%layer(nlayers))
         else if(sind(3,isec) == 1) then
            allocate(DECK_fluid_elastic_layer:: section%layer(nlayers))
         else
            allocate(DECK_maxwell_layer:: section%layer(nlayers))
         end if

         ilay = 1
         section%layer(ilay)%r1 = section%r1
         i11 = i1
         do iknot = i1+1,i2
            if(r(iknot-1) == r(iknot)) then
               section%layer(ilay)%r2 = r(iknot)
               i22 = iknot-1
               associate(layer => section%layer(ilay))
                 layer%r1 = r(i11)
                 layer%r2 = r(i22)
                 select type(layer)
                 type is(DECK_solid_elastic_layer)
                    call layer%rho_cubic%set(r(i11:i22),rho(i11:i22))
                    call layer%A_cubic%set(r(i11:i22),vph(i11:i22))
                    call layer%C_cubic%set(r(i11:i22),vpv(i11:i22))
                    call layer%F_cubic%set(r(i11:i22),eta(i11:i22))
                    call layer%L_cubic%set(r(i11:i22),vsv(i11:i22))
                    call layer%N_cubic%set(r(i11:i22),vsh(i11:i22))
                 type is(DECK_fluid_elastic_layer)                    
                    call layer%rho_cubic%set(r(i11:i22),rho(i11:i22))
                    call layer%kappa_cubic%set(r(i11:i22),vpv(i11:i22))
                 type is(DECK_maxwell_layer)
                    call layer%rho_cubic%set(r(i11:i22),rho(i11:i22))
                    call layer%A_cubic%set(r(i11:i22),vph(i11:i22))
                    call layer%C_cubic%set(r(i11:i22),vpv(i11:i22))
                    call layer%F_cubic%set(r(i11:i22),eta(i11:i22))
                    call layer%L_cubic%set(r(i11:i22),vsv(i11:i22))
                    call layer%N_cubic%set(r(i11:i22),vsh(i11:i22))
                    call layer%eta_cubic%set(r(i11:i22),visco(i11:i22))
                 end select
                 i11 = i22+1
               end associate
               ilay = ilay + 1
            end if
         end do
         section%layer(ilay)%r2 = section%r2
         i22 = i2
         associate(layer => section%layer(ilay))
           layer%r1 = r(i11)
           layer%r2 = r(i22)
           select type(layer)
           type is(DECK_solid_elastic_layer)
              call layer%rho_cubic%set(r(i11:i22),rho(i11:i22))
              call layer%A_cubic%set(r(i11:i22),vph(i11:i22))
              call layer%C_cubic%set(r(i11:i22),vpv(i11:i22))
              call layer%F_cubic%set(r(i11:i22),eta(i11:i22))
              call layer%L_cubic%set(r(i11:i22),vsv(i11:i22))
              call layer%N_cubic%set(r(i11:i22),vsh(i11:i22))
           type is(DECK_fluid_elastic_layer)                    
              call layer%rho_cubic%set(r(i11:i22),rho(i11:i22))
              call layer%kappa_cubic%set(r(i11:i22),vpv(i11:i22))
           type is(DECK_maxwell_layer)
              call layer%rho_cubic%set(r(i11:i22),rho(i11:i22))
              call layer%A_cubic%set(r(i11:i22),vph(i11:i22))
              call layer%C_cubic%set(r(i11:i22),vpv(i11:i22))
              call layer%F_cubic%set(r(i11:i22),eta(i11:i22))
              call layer%L_cubic%set(r(i11:i22),vsv(i11:i22))
              call layer%N_cubic%set(r(i11:i22),vsh(i11:i22))
              call layer%eta_cubic%set(r(i11:i22),visco(i11:i22))
           end select
         end associate
         
       end associate
       
    end do
        
    return
  end function maxwell_DECK

  
  function rho_DECK_maxwell_layer(self,r) result(rho)
    class(DECK_maxwell_layer), intent(in) :: self
    real(dp), intent(in) :: r
    real(dp) :: rho
    rho = self%rho_cubic%val(r)
    return
  end function rho_DECK_maxwell_layer

  function A_DECK_maxwell_layer(self,r) result(A)
    class(DECK_maxwell_layer), intent(in) :: self
    real(dp), intent(in) :: r
    real(dp) :: A,vpv    
    A = self%A_cubic%val(r)
    return
  end function A_DECK_maxwell_layer
  
  function C_DECK_maxwell_layer(self,r) result(C)
    class(DECK_maxwell_layer), intent(in) :: self
    real(dp), intent(in) :: r
    real(dp) :: C,rho,vsv
    C = self%C_cubic%val(r)
    return
  end function C_DECK_maxwell_layer

  function F_DECK_maxwell_layer(self,r) result(F)
    class(DECK_maxwell_layer), intent(in) :: self
    real(dp), intent(in) :: r
    real(dp) :: F
    F = self%F_cubic%val(r)
    return
  end function F_DECK_maxwell_layer

  function L_DECK_maxwell_layer(self,r) result(L)
    class(DECK_maxwell_layer), intent(in) :: self
    real(dp), intent(in) :: r
    real(dp) :: L
    L = self%L_cubic%val(r)
    return
  end function L_DECK_maxwell_layer

  function N_DECK_maxwell_layer(self,r) result(N)
    class(DECK_maxwell_layer), intent(in) :: self
    real(dp), intent(in) :: r
    real(dp) :: N
    N = self%N_cubic%val(r)
    return
  end function N_DECK_maxwell_layer

  function eta_DECK_maxwell_layer(self,r) result(eta)
    class(DECK_maxwell_layer), intent(in) :: self
    real(dp), intent(in) :: r
    real(dp) :: eta
    eta = self%eta_cubic%val(r)
    return
  end function eta_DECK_maxwell_layer

  
end module module_spherical_model
