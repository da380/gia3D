module module_spherical_model

  use module_constants  
  use module_physical_constants
  use module_util, only : poly_eval
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
     character(len=256) :: layer_type
     real(dp) :: r1
     real(dp) :: r2
     class(spherical_layer), dimension(:), allocatable :: layer
  end type spherical_section


  type spherical_model
     real(dp) :: r1
     real(dp) :: r2
     integer(i4b) :: nsections = 0
     class(spherical_section), dimension(:), allocatable :: section     
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
    model%section(1)%layer_type = 'solid_elastic'
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
    model%section(2)%layer_type = 'fluid_elastic'
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
    model%section(3)%layer_type = 'solid_elastic'
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
       model%section(4)%layer_type = 'fluid_elastic'
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
    real(dp) :: A,vpv    
    A = self%rho(r)*poly_eval(4,vpv_coef_PREM(:,self%i),r/r_PREM(2,13))**2
    return
  end function A_PREM_solid_elastic_layer
  
  function C_PREM_solid_elastic_layer(self,r) result(C)
    class(PREM_solid_elastic_layer), intent(in) :: self
    real(dp), intent(in) :: r
    real(dp) :: C,rho,vsv
    C = self%rho(r)*poly_eval(4,vsv_coef_PREM(:,self%i),r/r_PREM(2,13))**2
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





  
end module module_spherical_model
