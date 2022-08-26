module module_spherical_model

  use module_constants
  use module_interp
  implicit none


  !==============================================================!
  !            constants for non-dimensionalisation              !
  !==============================================================!

  real(dp), parameter :: length_norm   = 6371000.0_dp
  real(dp), parameter :: mass_norm     = 5.972e24_dp
  real(dp), parameter :: time_norm     = 3600.0_dp
  real(dp), parameter :: frequency_norm = 1.0_dp/time_norm
  real(dp), parameter :: velocity_norm = length_norm/time_norm
  real(dp), parameter :: acceleration_norm = length_norm/time_norm**2
  real(dp), parameter :: density_norm = mass_norm/length_norm**3
  real(dp), parameter :: force_norm = mass_norm*length_norm/time_norm**2
  real(dp), parameter :: gravitational_potential_norm = length_norm**2/time_norm**2
  real(dp), parameter :: gravitational_constant_norm = length_norm**3/(mass_norm*time_norm**2)
  real(dp), parameter :: modulus_norm = mass_norm/(time_norm**2*length_norm)
  real(dp), parameter :: viscosity_norm = mass_norm/(time_norm*length_norm)

  !==============================================================!
  !           type decalaration for elastic models               !
  !==============================================================!


  type, abstract :: spherical_elastic_model
     logical :: allocated = .false.
     logical, dimension(:), allocatable :: region_fluid
     logical, dimension(:), allocatable :: layer_fluid
     integer(i4b) :: nregions = 0
     integer(i4b) :: nlayers  = 0
     integer(i4b), dimension(:,:), allocatable :: region_index
     real(dp) :: rsurf
     real(dp), dimension(:,:), allocatable :: region_radius
     real(dp), dimension(:,:), allocatable :: layer_radius
   contains
     procedure(spherical_elastic_model_function), deferred :: rho
     procedure(spherical_elastic_model_function), deferred :: A
     procedure(spherical_elastic_model_function), deferred :: C
     procedure(spherical_elastic_model_function), deferred :: F
     procedure(spherical_elastic_model_function), deferred :: L
     procedure(spherical_elastic_model_function), deferred :: N    
     procedure :: kappa => kappa_spherical_elastic_model
     procedure :: mu => mu_spherical_elastic_model
     procedure :: vp => vp_spherical_elastic_model
     procedure :: vs => vs_spherical_elastic_model
     procedure :: vpv => vpv_spherical_elastic_model
     procedure :: vph => vph_spherical_elastic_model
     procedure :: vsv => vsv_spherical_elastic_model
     procedure :: vsh => vsh_spherical_elastic_model
     procedure :: mass => mass_spherical_elastic_model
     procedure :: gsurf => gsurf_spherical_elastic_model
     procedure :: write => write_spherical_elastic_model
  end type spherical_elastic_model

  
  abstract interface
     function spherical_elastic_model_function(self,i,r) result(f)
       use module_constants
       import :: spherical_elastic_model
       class(spherical_elastic_model), intent(in) :: self
       integer(i4b), intent(in) :: i
       real(dp), intent(in) :: r
       real(dp) :: f       
     end function spherical_elastic_model_function
  end interface

  
  type, abstract, extends(spherical_elastic_model) :: spherical_maxwell_model
     logical, dimension(:), allocatable :: region_viscous
     logical, dimension(:), allocatable :: layer_viscous
   contains
     procedure(spherical_maxwell_model_function), deferred :: eta
     procedure :: write => write_spherical_maxwell_model
  end type spherical_maxwell_model


  abstract interface
     function spherical_maxwell_model_function(self,i,r) result(f)
       use module_constants
       import :: spherical_maxwell_model
       class(spherical_maxwell_model), intent(in) :: self
       integer(i4b), intent(in) :: i
       real(dp), intent(in) :: r
       real(dp) :: f       
     end function spherical_maxwell_model_function
  end interface
  
  !==============================================================!
  !                          Deck models                         !
  !==============================================================!

  type, extends(spherical_elastic_model) :: elastic_deck_model
     type(interp_1D_cubic), dimension(:), allocatable :: rho_cubic
     type(interp_1D_cubic), dimension(:), allocatable :: A_cubic
     type(interp_1D_cubic), dimension(:), allocatable :: C_cubic
     type(interp_1D_cubic), dimension(:), allocatable :: F_cubic
     type(interp_1D_cubic), dimension(:), allocatable :: L_cubic
     type(interp_1D_cubic), dimension(:), allocatable :: N_cubic
   contains
     procedure :: rho =>  rho_elastic_deck_model
     procedure :: A =>  A_elastic_deck_model
     procedure :: C =>  C_elastic_deck_model
     procedure :: F =>  F_elastic_deck_model
     procedure :: L =>  L_elastic_deck_model
     procedure :: N =>  N_elastic_deck_model
  end type elastic_deck_model
  
  interface elastic_deck_model
     procedure :: set_elastic_deck_model
  end interface elastic_deck_model


  type, extends(spherical_maxwell_model) :: maxwell_deck_model
     type(interp_1D_cubic), dimension(:), allocatable :: rho_cubic
     type(interp_1D_cubic), dimension(:), allocatable :: A_cubic
     type(interp_1D_cubic), dimension(:), allocatable :: C_cubic
     type(interp_1D_cubic), dimension(:), allocatable :: F_cubic
     type(interp_1D_cubic), dimension(:), allocatable :: L_cubic
     type(interp_1D_cubic), dimension(:), allocatable :: N_cubic
     type(interp_1D_cubic), dimension(:), allocatable :: eta_cubic
   contains
     procedure :: rho =>  rho_maxwell_deck_model
     procedure :: A =>  A_maxwell_deck_model
     procedure :: C =>  C_maxwell_deck_model
     procedure :: F =>  F_maxwell_deck_model
     procedure :: L =>  L_maxwell_deck_model
     procedure :: N =>  N_maxwell_deck_model
     procedure :: eta => eta_maxwell_deck_model
  end type maxwell_deck_model


  interface maxwell_deck_model
     procedure :: set_maxwell_deck_model
  end interface maxwell_deck_model
  
  !==============================================================!
  !                          PREM model                          !
  !==============================================================!


  
  type, extends(spherical_elastic_model) :: elastic_PREM
   contains
     procedure :: rho =>  rho_elastic_PREM
     procedure :: A =>  A_elastic_PREM
     procedure :: C =>  C_elastic_PREM
     procedure :: F =>  F_elastic_PREM
     procedure :: L =>  L_elastic_PREM
     procedure :: N =>  N_elastic_PREM
  end type elastic_PREM

  interface elastic_PREM
     procedure :: set_elastic_PREM
  end interface elastic_PREM

  
  
  type, extends(spherical_maxwell_model) :: maxwell_PREM
     real(dp) :: eta_lithosphere
     real(dp) :: eta_upper_mantle
     real(dp) :: eta_lower_mantle
   contains
     procedure :: rho =>  rho_maxwell_PREM
     procedure :: A =>  A_maxwell_PREM
     procedure :: C =>  C_maxwell_PREM
     procedure :: F =>  F_maxwell_PREM
     procedure :: L =>  L_maxwell_PREM
     procedure :: N =>  N_maxwell_PREM
     procedure :: eta => eta_maxwell_PREM
  end type maxwell_PREM

  interface maxwell_PREM
     procedure :: set_maxwell_PREM
  end interface maxwell_PREM

  
  

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
  !                       procedures for elastic models                       !
  !===========================================================================!

  function kappa_spherical_elastic_model(self,i,r) result(kappa)
    class(spherical_elastic_model), intent(in) :: self
    integer(i4b), intent(in) :: i
    real(dp), intent(in) :: r
    real(dp) :: kappa
    kappa = self%C(i,r) + 4.0_dp*self%A(i,r) - 4.0_dp*self%N(i,r) + 4.0_dp*self%F(i,r)
    kappa = kappa/9.0_dp
    return
  end function kappa_spherical_elastic_model


  function mu_spherical_elastic_model(self,i,r) result(mu)
    class(spherical_elastic_model), intent(in) :: self
    integer(i4b), intent(in) :: i
    real(dp), intent(in) :: r
    real(dp) :: mu
    mu = self%C(i,r) + self%A(i,r) + 6.0_dp*self%L(i,r) + 5.0_dp*self%N(i,r) - 2.0_dp*self%F(i,r)
    mu = mu/15.0_dp
    return
  end function mu_spherical_elastic_model


  function vpv_spherical_elastic_model(self,i,r) result(vpv)
    class(spherical_elastic_model), intent(in) :: self
    integer(i4b), intent(in) :: i
    real(dp), intent(in) :: r
    real(dp) :: vpv
    vpv = sqrt(self%C(i,r)/self%rho(i,r))
    return
  end function vpv_spherical_elastic_model


  function vph_spherical_elastic_model(self,i,r) result(vph)
    class(spherical_elastic_model), intent(in) :: self
    integer(i4b), intent(in) :: i
    real(dp), intent(in) :: r
    real(dp) :: vph
    vph = sqrt(self%A(i,r)/self%rho(i,r))
    return
  end function vph_spherical_elastic_model

  
  function vsv_spherical_elastic_model(self,i,r) result(vsv)
    class(spherical_elastic_model), intent(in) :: self
    integer(i4b), intent(in) :: i
    real(dp), intent(in) :: r
    real(dp) :: vsv
    vsv = sqrt(self%L(i,r)/self%rho(i,r))
    return
  end function vsv_spherical_elastic_model


  function vsh_spherical_elastic_model(self,i,r) result(vsh)
    class(spherical_elastic_model), intent(in) :: self
    integer(i4b), intent(in) :: i
    real(dp), intent(in) :: r
    real(dp) :: vsh
    vsh = sqrt(self%N(i,r)/self%rho(i,r))
    return
  end function vsh_spherical_elastic_model


  function vp_spherical_elastic_model(self,i,r) result(vp)
    class(spherical_elastic_model), intent(in) :: self
    integer(i4b), intent(in) :: i
    real(dp), intent(in) :: r
    real(dp) :: vp
    vp = sqrt((self%kappa(i,r) + 4.0_dp*self%mu(i,r)/3.0_dp)/self%rho(i,r))
    return
  end function vp_spherical_elastic_model


  
  function vs_spherical_elastic_model(self,i,r) result(vs)
    class(spherical_elastic_model), intent(in) :: self
    integer(i4b), intent(in) :: i
    real(dp), intent(in) :: r
    real(dp) :: vs
    vs = sqrt(self%mu(i,r)/self%rho(i,r))
    return
  end function vs_spherical_elastic_model
  

  function mass_spherical_elastic_model(self,n) result(mass)
    use module_special_functions
    use module_quadrature
    class(spherical_elastic_model), intent(in) :: self
    integer(i4b), intent(in), optional :: n
    real(dp) :: mass
    integer(i4b) :: ilayer,i,m
    real(dp) :: r1,r2,r,w
    class(orthogonal_polynomial), allocatable :: poly
    type(gauss_quadrature) :: quad
    poly = legendre()
    if(present(n)) then
       m = n
    else
       m = 5
    end if
    call quad%set(m,poly)
    mass = 0.0_dp
    do ilayer = 1,self%nlayers
       r1 = self%layer_radius(1,ilayer)
       r2 = self%layer_radius(2,ilayer)
       call quad%trans(r1,r2)
       do i = 1,m
          r = quad%x(i)
          w = quad%w(i)
          mass = mass +  fourpi*self%rho(ilayer,r)*r*r*w
       end do
    end do
    return
  end function mass_spherical_elastic_model


  function gsurf_spherical_elastic_model(self,n) result(gsurf)
    use module_special_functions
    use module_quadrature
    use module_physical_constants, only : bigg => G
    class(spherical_elastic_model), intent(in) :: self
    integer(i4b), intent(in), optional :: n
    real(dp) :: gsurf
    gsurf = (bigg/gravitational_constant_norm)*self%mass(n)/self%rsurf**2
    return
  end function gsurf_spherical_elastic_model


  subroutine write_spherical_elastic_model(self,nknot,model_name,isotropic)
    class(spherical_elastic_model), intent(in) :: self
    integer(i4b), intent(in) :: nknot
    character(len=*), intent(in) :: model_name
    logical, intent(in), optional :: isotropic
    logical :: iso
    integer(i4b) :: ilayer,nloc,iloc,io
    real(dp) :: r1,r2,r,dr
    if(present(isotropic)) then
       iso = isotropic
    else
       iso = .false.
    end if
    open(newunit = io,file=trim(model_name))
    write(io,*) iso
    do ilayer = 1,self%nlayers
       r1 = self%layer_radius(1,ilayer)
       r2 = self%layer_radius(2,ilayer)
       nloc = max(int(nknot*(r2-r1)/self%rsurf),2)
       dr = (r2-r1)/(nloc-1)
       do iloc = 1,nloc
          r = r1 + (iloc-1)*dr
          if(iso) then             
             write(io,'(4e20.8)') r*length_norm,            &
                         self%rho(ilayer,r)*density_norm,   &
                         self%kappa(ilayer,r)*modulus_norm, &
                         self%mu(ilayer,r)*modulus_norm
          else
             write(io,'(7e20.8)') r*length_norm,            &
                         self%rho(ilayer,r)*density_norm,   &
                         self%A(ilayer,r)*modulus_norm,     &
                         self%C(ilayer,r)*modulus_norm,     &
                         self%F(ilayer,r)*modulus_norm,     &
                         self%L(ilayer,r)*modulus_norm,     &
                         self%N(ilayer,r)*modulus_norm
          end if
       end do       
    end do
    close(io)    
    return
  end subroutine write_spherical_elastic_model



  subroutine write_spherical_maxwell_model(self,nknot,model_name,isotropic)
    class(spherical_maxwell_model), intent(in) :: self
    integer(i4b), intent(in) :: nknot
    character(len=*), intent(in) :: model_name
    logical, intent(in), optional :: isotropic
    logical :: iso
    integer(i4b) :: ilayer,nloc,iloc,io
    real(dp) :: r1,r2,r,dr,eta
    if(present(isotropic)) then
       iso = isotropic
    else
       iso = .false.
    end if
    open(newunit = io,file=trim(model_name))
    write(io,*) iso
    do ilayer = 1,self%nlayers
       r1 = self%layer_radius(1,ilayer)
       r2 = self%layer_radius(2,ilayer)
       nloc = max(int(nknot*(r2-r1)/self%rsurf),2)
       dr = (r2-r1)/(nloc-1)
       do iloc = 1,nloc
          r = r1 + (iloc-1)*dr
          eta = self%eta(ilayer,r)*viscosity_norm
          if(eta > 0.0_dp) eta = log10(eta)
          if(iso) then             
             write(io,'(5e20.8)') r*length_norm,            &
                         self%rho(ilayer,r)*density_norm,   &
                         self%kappa(ilayer,r)*modulus_norm, &
                         self%mu(ilayer,r)*modulus_norm,    &
                         eta
          else
             write(io,'(8e20.8)') r*length_norm,            &
                         self%rho(ilayer,r)*density_norm,   &
                         self%A(ilayer,r)*modulus_norm,     &
                         self%C(ilayer,r)*modulus_norm,     &
                         self%F(ilayer,r)*modulus_norm,     &
                         self%L(ilayer,r)*modulus_norm,     &
                         self%N(ilayer,r)*modulus_norm,     &
                         eta
          end if
       end do       
    end do
    close(io)    
    return
  end subroutine write_spherical_maxwell_model



  !===========================================================================!
  !                     procedures for elastic deck models                    !
  !===========================================================================!

  
  function set_elastic_deck_model(model_file) result(model)
    use module_error
    character(len=*), intent(in) :: model_file
    type(elastic_deck_model) :: model
    logical :: ltmp,iso,fluid
    integer(i4b) :: io,nknot,ios,i,nlayers,nregions,i1,i2,ilayer,iregion
    real(dp), dimension(:), allocatable :: r,rho,A,C,F,L,N
    
    ! check the model file exists
    inquire(file = trim(model_file),exist = ltmp)
    call error(.not.ltmp,'set_elastic_deck_model','model file does not exist')
    
    ! open the file and work out the format
    open(newunit = io,file = trim(model_file))
    read(io,*) iso
    nknot = 0
    do
       read(io,*,iostat = ios)
       if(ios /= 0) exit
       nknot = nknot + 1
    end do
    rewind(io)
    read(io,*) iso
    
    ! read in the values and non-dimensionalise
    allocate(r(nknot),rho(nknot),A(nknot),C(nknot),F(nknot),L(nknot),N(nknot))
    do i = 1,nknot
       if(iso) then
          read(io,*) r(i),rho(i),A(i),L(i)
          A(i) = A(i) + 4.0_dp*L(i)/3.0_dp
          C(i) = A(i)
          N(i) = L(i)
          F(i) = A(i)-2.0_dp*L(i)
       else
          read(io,*) r(i),rho(i),A(i),C(i),F(i),L(i),N(i)
       end if
       r(i)   = r(i)/length_norm
       rho(i) = rho(i)/density_norm
       A(i)   = A(i)/modulus_norm
       C(i)   = C(i)/modulus_norm
       F(i)   = F(i)/modulus_norm
       L(i)   = L(i)/modulus_norm
       N(i)   = N(i)/modulus_norm
    end do
    close(io)

    
    ! work out the number of layers
    model%nlayers = 1
    do i = 1,nknot-1
       if(r(i+1) == r(i)) model%nlayers = model%nlayers + 1       
    end do

    ! build the cubic splines
    allocate(model%rho_cubic(model%nlayers))
    allocate(model%A_cubic(model%nlayers))
    allocate(model%C_cubic(model%nlayers))
    allocate(model%F_cubic(model%nlayers))
    allocate(model%L_cubic(model%nlayers))
    allocate(model%N_cubic(model%nlayers))
    allocate(model%layer_radius(2,model%nlayers))
    allocate(model%layer_fluid(model%nlayers))
    model%layer_fluid(:) = .false.
    i1 = 1
    ilayer = 0
    do i = 1,nknot-1
       if(r(i+1) == r(i)) then
          ilayer = ilayer+1
          i2 = i
          model%layer_radius(1,ilayer) = r(i1)
          model%layer_radius(2,ilayer) = r(i2)
          call model%rho_cubic(ilayer)%set(r(i1:i2),rho(i1:i2))
          call model%A_cubic(ilayer)%set(r(i1:i2),A(i1:i2))
          call model%C_cubic(ilayer)%set(r(i1:i2),C(i1:i2))
          call model%F_cubic(ilayer)%set(r(i1:i2),F(i1:i2))
          call model%L_cubic(ilayer)%set(r(i1:i2),L(i1:i2))
          call model%N_cubic(ilayer)%set(r(i1:i2),N(i1:i2))
          if(all(L(i1:i2) == 0.0_dp)) model%layer_fluid(ilayer) = .true.
          i1 = i+1
       end if
    end do
    ilayer = ilayer + 1 
    i2 = nknot
    model%layer_radius(1,ilayer) = r(i1)
    model%layer_radius(2,ilayer) = r(i2)
    call model%rho_cubic(ilayer)%set(r(i1:i2),rho(i1:i2))
    call model%A_cubic(ilayer)%set(r(i1:i2),A(i1:i2))
    call model%C_cubic(ilayer)%set(r(i1:i2),C(i1:i2))
    call model%F_cubic(ilayer)%set(r(i1:i2),F(i1:i2))
    call model%L_cubic(ilayer)%set(r(i1:i2),L(i1:i2))
    call model%N_cubic(ilayer)%set(r(i1:i2),N(i1:i2))
    if(all(L(i1:i2) == 0.0_dp)) model%layer_fluid(ilayer) = .true.


    ! work out the number of regions
    model%nregions = 1
    do ilayer = 1,model%nlayers-1
       if(model%layer_fluid(ilayer) .neqv. model%layer_fluid(ilayer+1)) model%nregions = model%nregions + 1
    end do

    ! set the region indices etc
    allocate(model%region_index(2,model%nregions))
    allocate(model%region_radius(2,model%nregions))
    allocate(model%region_fluid(model%nregions))       
    iregion = 1
    model%region_index(1,iregion) = 1
    do ilayer = 1,model%nlayers-1
       if(model%layer_fluid(ilayer) .neqv. model%layer_fluid(ilayer+1)) then
          model%region_index(2,iregion) = ilayer
          model%region_fluid(iregion) = model%layer_fluid(ilayer)
          iregion = iregion + 1
          model%region_index(1,iregion) = ilayer+1
       end if
    end do
    if(model%layer_fluid(model%nlayers)) model%region_fluid(model%nregions) = .true.
    model%region_index(2,iregion) = model%nlayers
    do iregion = 1,model%nregions
       model%region_radius(1,iregion) = model%layer_radius(1,model%region_index(1,iregion))
       model%region_radius(2,iregion) = model%layer_radius(2,model%region_index(2,iregion))
    end do

    ! finalise the set up
    model%rsurf = model%layer_radius(2,model%nlayers)
    model%allocated = .true.
    
    return
  end function set_elastic_deck_model

  function rho_elastic_deck_model(self,i,r) result(rho)
    class(elastic_deck_model), intent(in) :: self
    integer(i4b), intent(in) :: i
    real(dp), intent(in) :: r
    real(dp) :: rho
    rho = self%rho_cubic(i)%val(r)
    return
  end function rho_elastic_deck_model


  function A_elastic_deck_model(self,i,r) result(A)
    class(elastic_deck_model), intent(in) :: self
    integer(i4b), intent(in) :: i
    real(dp), intent(in) :: r
    real(dp) :: A
    A = self%A_cubic(i)%val(r)
    return
  end function A_elastic_deck_model


  function C_elastic_deck_model(self,i,r) result(C)
    class(elastic_deck_model), intent(in) :: self
    integer(i4b), intent(in) :: i
    real(dp), intent(in) :: r
    real(dp) :: C
    C = self%C_cubic(i)%val(r)
    return
  end function C_elastic_deck_model


  function F_elastic_deck_model(self,i,r) result(F)
    class(elastic_deck_model), intent(in) :: self
    integer(i4b), intent(in) :: i
    real(dp), intent(in) :: r
    real(dp) :: F
    F = self%F_cubic(i)%val(r)
    return
  end function F_elastic_deck_model


  function L_elastic_deck_model(self,i,r) result(L)
    class(elastic_deck_model), intent(in) :: self
    integer(i4b), intent(in) :: i
    real(dp), intent(in) :: r
    real(dp) :: L
    L = self%L_cubic(i)%val(r)
    return
  end function L_elastic_deck_model


  
  function N_elastic_deck_model(self,i,r) result(N)
    class(elastic_deck_model), intent(in) :: self
    integer(i4b), intent(in) :: i
    real(dp), intent(in) :: r
    real(dp) :: N
    N = self%N_cubic(i)%val(r)
    return
  end function N_elastic_deck_model



  !===========================================================================!
  !                     procedures for maxwell deck models                    !
  !===========================================================================!

  function set_maxwell_deck_model(model_file) result(model)
    use module_error
    character(len=*), intent(in) :: model_file
    type(maxwell_deck_model) :: model
    logical :: ltmp,iso,fluid
    integer(i4b) :: io,nknot,ios,i,nlayers,nregions,i1,i2,ilayer,iregion
    real(dp), dimension(:), allocatable :: r,rho,A,C,F,L,N,eta
    
    ! check the model file exists
    inquire(file = trim(model_file),exist = ltmp)
    call error(.not.ltmp,'set_maxwell_deck_model','model file does not exist')
    
    ! open the file and work out the format
    open(newunit = io,file = trim(model_file))
    read(io,*) iso
    nknot = 0
    do
       read(io,*,iostat = ios)
       if(ios /= 0) exit
       nknot = nknot + 1
    end do
    rewind(io)
    read(io,*) iso
    
    ! read in the values and non-dimensionalise
    allocate(r(nknot),rho(nknot),A(nknot),C(nknot), &
             F(nknot),L(nknot),N(nknot),eta(nknot))
    do i = 1,nknot
       if(iso) then
          read(io,*) r(i),rho(i),A(i),L(i),eta(i)
          A(i) = A(i) + 4.0_dp*L(i)/3.0_dp
          C(i) = A(i)
          N(i) = L(i)
          F(i) = A(i)-2.0_dp*L(i)
       else
          read(io,*) r(i),rho(i),A(i),C(i),F(i),L(i),N(i),eta(i)
       end if
       r(i)   = r(i)/length_norm
       rho(i) = rho(i)/density_norm
       A(i)   = A(i)/modulus_norm
       C(i)   = C(i)/modulus_norm
       F(i)   = F(i)/modulus_norm
       L(i)   = L(i)/modulus_norm
       N(i)   = N(i)/modulus_norm
       if(eta(i) > 0.0_dp) then
          eta(i) = (10**eta(i))/viscosity_norm
       end if
    end do
    close(io)

    
    ! work out the number of layers
    model%nlayers = 1
    do i = 1,nknot-1
       if(r(i+1) == r(i)) model%nlayers = model%nlayers + 1       
    end do

    ! build the cubic splines
    allocate(model%rho_cubic(model%nlayers))
    allocate(model%A_cubic(model%nlayers))
    allocate(model%C_cubic(model%nlayers))
    allocate(model%F_cubic(model%nlayers))
    allocate(model%L_cubic(model%nlayers))
    allocate(model%N_cubic(model%nlayers))
    allocate(model%eta_cubic(model%nlayers))
    allocate(model%layer_radius(2,model%nlayers))
    allocate(model%layer_fluid(model%nlayers))
    allocate(model%layer_viscous(model%nlayers))
    model%layer_fluid(:) = .false.
    model%layer_viscous(:) = .true.
    i1 = 1
    ilayer = 0
    do i = 1,nknot-1
       if(r(i+1) == r(i)) then
          ilayer = ilayer+1
          i2 = i
          model%layer_radius(1,ilayer) = r(i1)
          model%layer_radius(2,ilayer) = r(i2)
          call model%rho_cubic(ilayer)%set(r(i1:i2),rho(i1:i2))
          call model%A_cubic(ilayer)%set(r(i1:i2),A(i1:i2))
          call model%C_cubic(ilayer)%set(r(i1:i2),C(i1:i2))
          call model%F_cubic(ilayer)%set(r(i1:i2),F(i1:i2))
          call model%L_cubic(ilayer)%set(r(i1:i2),L(i1:i2))
          call model%N_cubic(ilayer)%set(r(i1:i2),N(i1:i2))
          call model%eta_cubic(ilayer)%set(r(i1:i2),eta(i1:i2))
          if(all(L(i1:i2) == 0.0_dp)) model%layer_fluid(ilayer) = .true.
          if(all(eta(i1:i2) == 0.0_dp)) model%layer_viscous(ilayer) = .false.
          i1 = i+1
       end if
    end do
    ilayer = ilayer + 1 
    i2 = nknot
    model%layer_radius(1,ilayer) = r(i1)
    model%layer_radius(2,ilayer) = r(i2)
    call model%rho_cubic(ilayer)%set(r(i1:i2),rho(i1:i2))
    call model%A_cubic(ilayer)%set(r(i1:i2),A(i1:i2))
    call model%C_cubic(ilayer)%set(r(i1:i2),C(i1:i2))
    call model%F_cubic(ilayer)%set(r(i1:i2),F(i1:i2))
    call model%L_cubic(ilayer)%set(r(i1:i2),L(i1:i2))
    call model%N_cubic(ilayer)%set(r(i1:i2),N(i1:i2))
    call model%eta_cubic(ilayer)%set(r(i1:i2),eta(i1:i2))
    if(all(L(i1:i2) == 0.0_dp)) model%layer_fluid(ilayer) = .true.
    if(all(eta(i1:i2) == 0.0_dp)) model%layer_viscous(ilayer) = .false.

    
    ! work out the number of regions
    model%nregions = 1
    do ilayer = 1,model%nlayers-1
       if(model%layer_fluid(ilayer) .neqv. model%layer_fluid(ilayer+1)) model%nregions = model%nregions + 1
    end do

    ! set the region indices etc
    allocate(model%region_index(2,model%nregions))
    allocate(model%region_radius(2,model%nregions))
    allocate(model%region_fluid(model%nregions))
    allocate(model%region_viscous(model%nregions))       
    iregion = 1
    model%region_index(1,iregion) = 1
    do ilayer = 1,model%nlayers-1
       if(model%layer_fluid(ilayer) .neqv. model%layer_fluid(ilayer+1)) then
          model%region_index(2,iregion) = ilayer
          model%region_fluid(iregion) = model%layer_fluid(ilayer)
          model%region_viscous(iregion) = model%layer_viscous(ilayer)
          iregion = iregion + 1
          model%region_index(1,iregion) = ilayer+1
       end if
    end do
    if(model%layer_fluid(model%nlayers)) model%region_fluid(model%nregions) = .true.
    model%region_index(2,iregion) = model%nlayers
    do iregion = 1,model%nregions
       model%region_radius(1,iregion) = model%layer_radius(1,model%region_index(1,iregion))
       model%region_radius(2,iregion) = model%layer_radius(2,model%region_index(2,iregion))
    end do

    
    ! finalise the set up
    model%rsurf = model%layer_radius(2,model%nlayers)
    model%allocated = .true.
    
    return
  end function set_maxwell_deck_model

  function rho_maxwell_deck_model(self,i,r) result(rho)
    class(maxwell_deck_model), intent(in) :: self
    integer(i4b), intent(in) :: i
    real(dp), intent(in) :: r
    real(dp) :: rho
    rho = self%rho_cubic(i)%val(r)
    return
  end function rho_maxwell_deck_model


  function A_maxwell_deck_model(self,i,r) result(A)
    class(maxwell_deck_model), intent(in) :: self
    integer(i4b), intent(in) :: i
    real(dp), intent(in) :: r
    real(dp) :: A
    A = self%A_cubic(i)%val(r)
    return
  end function A_maxwell_deck_model


  function C_maxwell_deck_model(self,i,r) result(C)
    class(maxwell_deck_model), intent(in) :: self
    integer(i4b), intent(in) :: i
    real(dp), intent(in) :: r
    real(dp) :: C
    C = self%C_cubic(i)%val(r)
    return
  end function C_maxwell_deck_model


  function F_maxwell_deck_model(self,i,r) result(F)
    class(maxwell_deck_model), intent(in) :: self
    integer(i4b), intent(in) :: i
    real(dp), intent(in) :: r
    real(dp) :: F
    F = self%F_cubic(i)%val(r)
    return
  end function F_maxwell_deck_model


  function L_maxwell_deck_model(self,i,r) result(L)
    class(maxwell_deck_model), intent(in) :: self
    integer(i4b), intent(in) :: i
    real(dp), intent(in) :: r
    real(dp) :: L
    L = self%L_cubic(i)%val(r)
    return
  end function L_maxwell_deck_model

  
  function N_maxwell_deck_model(self,i,r) result(N)
    class(maxwell_deck_model), intent(in) :: self
    integer(i4b), intent(in) :: i
    real(dp), intent(in) :: r
    real(dp) :: N
    N = self%N_cubic(i)%val(r)
    return
  end function N_maxwell_deck_model
  
  
  function eta_maxwell_deck_model(self,i,r) result(eta)
    class(maxwell_deck_model), intent(in) :: self
    integer(i4b), intent(in) :: i
    real(dp), intent(in) :: r
    real(dp) :: eta
    eta = self%eta_cubic(i)%val(r)
    return
  end function eta_maxwell_deck_model
  
  
  !===========================================================================!
  !                              procedures for PREM                          !
  !===========================================================================!


  
  function set_elastic_PREM(ocean) result(model)
    logical, intent(in), optional :: ocean
    type(elastic_PREM) :: model
    logical :: locean
    if(present(ocean)) then
       locean = ocean
    else
       locean = .true.
    end if
    if(locean) then
       model%nregions = 4
       model%nlayers = 13      
    else
       model%nregions = 3
       model%nlayers = 12
    end if
    allocate(model%region_fluid(model%nregions))
    allocate(model%region_index(2,model%nregions))
    allocate(model%region_radius(2,model%nregions))
    allocate(model%layer_fluid(model%nlayers))
    allocate(model%layer_radius(2,model%nlayers))
    ! inner core
    model%region_index(1,1) = 1
    model%region_index(2,1) = 1
    model%region_radius(1,1) = r_PREM(1,1)
    model%region_radius(2,1) = r_PREM(1,1)
    model%region_fluid(1) = .false.
    model%layer_fluid(1) = .false.
    ! outer core
    model%region_index(1,2) = 2
    model%region_index(2,2) = 2
    model%region_radius(1,2) = r_PREM(1,2)
    model%region_radius(2,2) = r_PREM(1,2)
    model%region_fluid(2) = .true.
    model%layer_fluid(2) = .true.
    ! mantle
    model%region_index(1,3) = 3
    model%region_index(2,3) = 12
    model%region_radius(1,3) = r_PREM(1,3)
    model%region_radius(2,3) = r_PREM(1,12)
    model%region_fluid(3) = .false.
    model%layer_fluid(3:12) = .false.
    if(locean) then
       model%region_index(1,4) = 13
       model%region_index(2,4) = 13
       model%region_radius(1,4) = r_PREM(1,13)
       model%region_radius(2,4) = r_PREM(1,13)
       model%region_fluid(4) = .true.
       model%layer_fluid(13) = .true.
    end if
    model%layer_radius(:,1:model%nlayers) = r_PREM(:,1:model%nlayers)
    model%rsurf = model%layer_radius(2,model%nlayers)
    model%allocated = .true.
    return
  end function set_elastic_PREM
    
  function rho_elastic_PREM(self,i,r) result(rho)
    class(elastic_PREM), intent(in) :: self
    integer(i4b), intent(in) :: i
    real(dp), intent(in) :: r
    real(dp) :: rho
    rho = get_rho_PREM(i,r)
    return
  end function rho_elastic_PREM


  function A_elastic_PREM(self,i,r) result(A)
    class(elastic_PREM), intent(in) :: self
    integer(i4b), intent(in) :: i
    real(dp), intent(in) :: r
    real(dp) :: A
    A = get_rho_PREM(i,r)*get_vph_PREM(i,r)**2
    return
  end function A_elastic_PREM


  function C_elastic_PREM(self,i,r) result(C)
    class(elastic_PREM), intent(in) :: self
    integer(i4b), intent(in) :: i
    real(dp), intent(in) :: r
    real(dp) :: C
    C = get_rho_PREM(i,r)*get_vpv_PREM(i,r)**2
    return
  end function C_elastic_PREM


  function F_elastic_PREM(self,i,r) result(F)
    class(elastic_PREM), intent(in) :: self
    integer(i4b), intent(in) :: i
    real(dp), intent(in) :: r
    real(dp) :: F
    F = get_eta_PREM(i,r)*(self%A(i,r) - 2.0_dp*self%L(i,r))
    return
  end function F_elastic_PREM

  
  function L_elastic_PREM(self,i,r) result(L)
    class(elastic_PREM), intent(in) :: self
    integer(i4b), intent(in) :: i
    real(dp), intent(in) :: r
    real(dp) :: L
    L = get_rho_PREM(i,r)*get_vsv_PREM(i,r)**2
    return
  end function L_elastic_PREM

  
  function N_elastic_PREM(self,i,r) result(N)
    class(elastic_PREM), intent(in) :: self
    integer(i4b), intent(in) :: i
    real(dp), intent(in) :: r
    real(dp) :: N
    N = get_rho_PREM(i,r)*get_vsh_PREM(i,r)**2
    return
  end function N_elastic_PREM


  function set_maxwell_PREM(eta_lithosphere,eta_upper_mantle,eta_lower_mantle) result(model)
    type(maxwell_PREM) :: model
    real(dp), intent(in) :: eta_lithosphere
    real(dp), intent(in) :: eta_upper_mantle
    real(dp), intent(in) :: eta_lower_mantle
    model%eta_lithosphere  = (10**eta_lithosphere)/viscosity_norm
    model%eta_upper_mantle = (10**eta_upper_mantle)/viscosity_norm
    model%eta_lower_mantle = (10**eta_lower_mantle)/viscosity_norm
    model%nregions = 3
    model%nlayers = 12
    allocate(model%region_fluid(model%nregions))
    allocate(model%region_viscous(model%nregions))
    allocate(model%region_index(2,model%nregions))
    allocate(model%region_radius(2,model%nregions))
    allocate(model%layer_fluid(model%nlayers))
    allocate(model%layer_viscous(model%nregions))
    allocate(model%layer_radius(2,model%nlayers))
    ! inner core
    model%region_index(1,1) = 1
    model%region_index(2,1) = 1
    model%region_radius(1,1) = r_PREM(1,1)
    model%region_radius(2,1) = r_PREM(1,1)
    model%region_fluid(1) = .false.
    model%layer_fluid(1) = .false.
    model%region_viscous(1) = .false.
    model%layer_viscous(1) = .false.
    ! outer core
    model%region_index(1,2) = 2
    model%region_index(2,2) = 2
    model%region_radius(1,2) = r_PREM(1,2)
    model%region_radius(2,2) = r_PREM(1,2)
    model%region_fluid(2) = .true.
    model%layer_fluid(2) = .true.
    model%region_viscous(2) = .false.
    model%layer_viscous(2) = .false.
    ! mantle
    model%region_index(1,3) = 3
    model%region_index(2,3) = 12
    model%region_radius(1,3) = r_PREM(1,3)
    model%region_radius(2,3) = r_PREM(1,12)
    model%region_fluid(3) = .false.
    model%layer_fluid(3:12) = .false.
    model%region_viscous(3) = .true.
    model%layer_viscous(3:) = .true.
    model%layer_radius(:,1:model%nlayers) = r_PREM(:,1:model%nlayers)    
    model%rsurf = model%layer_radius(2,model%nlayers)
    model%allocated = .true.
    return
  end function set_maxwell_PREM
    
  function rho_maxwell_PREM(self,i,r) result(rho)
    class(maxwell_PREM), intent(in) :: self
    integer(i4b), intent(in) :: i
    real(dp), intent(in) :: r
    real(dp) :: rho
    rho = get_rho_PREM(i,r)
    return
  end function rho_maxwell_PREM


  function A_maxwell_PREM(self,i,r) result(A)
    class(maxwell_PREM), intent(in) :: self
    integer(i4b), intent(in) :: i
    real(dp), intent(in) :: r
    real(dp) :: A
    A = get_rho_PREM(i,r)*get_vph_PREM(i,r)**2
    return
  end function A_maxwell_PREM


  function C_maxwell_PREM(self,i,r) result(C)
    class(maxwell_PREM), intent(in) :: self
    integer(i4b), intent(in) :: i
    real(dp), intent(in) :: r
    real(dp) :: C
    C = get_rho_PREM(i,r)*get_vpv_PREM(i,r)**2
    return
  end function C_maxwell_PREM


  function F_maxwell_PREM(self,i,r) result(F)
    class(maxwell_PREM), intent(in) :: self
    integer(i4b), intent(in) :: i
    real(dp), intent(in) :: r
    real(dp) :: F
    F = get_eta_PREM(i,r)*(self%A(i,r) - 2.0_dp*self%L(i,r))
    return
  end function F_maxwell_PREM

  
  function L_maxwell_PREM(self,i,r) result(L)
    class(maxwell_PREM), intent(in) :: self
    integer(i4b), intent(in) :: i
    real(dp), intent(in) :: r
    real(dp) :: L
    L = get_rho_PREM(i,r)*get_vsv_PREM(i,r)**2
    return
  end function L_maxwell_PREM

  
  function N_maxwell_PREM(self,i,r) result(N)
    class(maxwell_PREM), intent(in) :: self
    integer(i4b), intent(in) :: i
    real(dp), intent(in) :: r
    real(dp) :: N
    N = get_rho_PREM(i,r)*get_vsh_PREM(i,r)**2
    return
  end function N_maxwell_PREM

  

  function eta_maxwell_PREM(self,i,r) result(eta)
    class(maxwell_PREM), intent(in) :: self
    integer(i4b), intent(in) :: i
    real(dp), intent(in) :: r
    real(dp) :: eta
    if(i >= 3 .and. i <= 7) then
       eta = self%eta_lower_mantle
    else if(i >= 7 .and. i <= 10) then
       eta = self%eta_upper_mantle
    else if(i >= 10) then
       eta = self%eta_lithosphere
    else
       eta = 0.0_dp
    end if
    return
  end function eta_maxwell_PREM



  function get_rho_PREM(i,r) result(rho)
    use module_util, only : poly_eval
    integer(i4b), intent(in) :: i
    real(dp), intent(in) :: r
    real(dp) :: rho
    rho = poly_eval(4,rho_coef_PREM(:,i),r/r_PREM(2,13))
    return
  end function get_rho_PREM


  function get_vpv_PREM(i,r) result(vpv)
    use module_util, only : poly_eval
    integer(i4b), intent(in) :: i
    real(dp), intent(in) :: r
    real(dp) :: vpv
    vpv = poly_eval(4,vpv_coef_PREM(:,i),r/r_PREM(2,13))
    return
  end function get_vpv_PREM


  function get_vsv_PREM(i,r) result(vsv)
    use module_util, only : poly_eval
    integer(i4b), intent(in) :: i
    real(dp), intent(in) :: r
    real(dp) :: vsv
    vsv = poly_eval(4,vsv_coef_PREM(:,i),r/r_PREM(2,13))
    return
  end function get_vsv_PREM


  function get_qk_PREM(i,r) result(qk)
    use module_util, only : poly_eval
    integer(i4b), intent(in) :: i
    real(dp), intent(in) :: r
    real(dp) :: qk
    qk = poly_eval(4,qk_coef_PREM(:,i),r/r_PREM(2,13))
    return
  end function get_qk_PREM


  function get_qm_PREM(i,r) result(qm)
    use module_util, only : poly_eval
    integer(i4b), intent(in) :: i
    real(dp), intent(in) :: r
    real(dp) :: qm
    qm = poly_eval(4,qm_coef_PREM(:,i),r/r_PREM(2,13))
    return
  end function get_qm_PREM

  
  function get_vph_PREM(i,r) result(vph)
    use module_util, only : poly_eval
    integer(i4b), intent(in) :: i
    real(dp), intent(in) :: r
    real(dp) :: vph
    vph = poly_eval(4,vph_coef_PREM(:,i),r/r_PREM(2,13))
    return
  end function get_vph_PREM


  function get_vsh_PREM(i,r) result(vsh)
    use module_util, only : poly_eval
    integer(i4b), intent(in) :: i
    real(dp), intent(in) :: r
    real(dp) :: vsh
    vsh = poly_eval(4,vsh_coef_PREM(:,i),r/r_PREM(2,13))
    return
  end function get_vsh_PREM

  
  function get_eta_PREM(i,r) result(eta)
    use module_util, only : poly_eval
    integer(i4b), intent(in) :: i
    real(dp), intent(in) :: r
    real(dp) :: eta
    eta = poly_eval(4,eta_coef_PREM(:,i),r/r_PREM(2,13))
    return
  end function get_eta_PREM

  
  

  
end module module_spherical_model
