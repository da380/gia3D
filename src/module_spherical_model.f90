module module_spherical_model

  use module_constants  
  use module_physical_constants
  use module_error
  use module_quadrature
  implicit none

  !========================================!
  !                basic types             !
  !========================================!

  type, abstract :: spherical_layer
     real(dp) :: r1
     real(dp) :: r2
   contains
     procedure(spherical_layer_function), deferred :: rho
  end type spherical_layer
 
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
     procedure :: M => mass_spherical_model
     procedure :: g => surface_gravity_spherical_model
     procedure :: write_elastic => write_elastic_spherical_model
     procedure :: write_maxwell => write_maxwell_spherical_model
     procedure :: write_anelastic => write_anelastic_spherical_model
  end type spherical_model

  abstract interface
     function spherical_layer_function(self,r) result(f)
       use module_constants
       import :: spherical_layer
       class(spherical_layer), intent(in) :: self
       real(dp), intent(in) :: r
       real(dp) :: f
     end function spherical_layer_function
  end interface


  !========================================!
  !                 elastic                !
  !========================================!
  
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

  !========================================!
  !                   Maxwell              !
  !========================================!  
  
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

  !========================================!
  !                anelastic               !
  !========================================!
  
  type, abstract, extends(spherical_solid_elastic_layer) :: spherical_solid_anelastic_layer
     real(dp) :: tref = 1.0_dp/time_norm
   contains
     procedure(spherical_solid_anelastic_layer_function), deferred :: qA
     procedure(spherical_solid_anelastic_layer_function), deferred :: qC
     procedure(spherical_solid_anelastic_layer_function), deferred :: qF
     procedure(spherical_solid_anelastic_layer_function), deferred :: qL
     procedure(spherical_solid_anelastic_layer_function), deferred :: qN
     procedure :: qkappa => qkappa_spherical_solid_anelastic_layer
     procedure :: qmu => qmu_spherical_solid_anelastic_layer
  end type spherical_solid_anelastic_layer

  abstract interface
     function spherical_solid_anelastic_layer_function(self,r) result(f)
       use module_constants
       import :: spherical_solid_anelastic_layer
       class(spherical_solid_anelastic_layer), intent(in) :: self
       real(dp), intent(in) :: r
       real(dp) :: f
     end function spherical_solid_anelastic_layer_function
  end interface


  type, abstract, extends(spherical_fluid_elastic_layer) :: spherical_fluid_anelastic_layer
   contains
     procedure(spherical_fluid_anelastic_layer_function), deferred :: qkappa
  end type spherical_fluid_anelastic_layer

  abstract interface
     function spherical_fluid_anelastic_layer_function(self,r) result(f)
       use module_constants
       import :: spherical_fluid_anelastic_layer
       class(spherical_fluid_anelastic_layer), intent(in) :: self
       real(dp), intent(in) :: r
       real(dp) :: f
     end function spherical_fluid_anelastic_layer_function
  end interface
  


  !-------------------------------------------!
  !            model_parameters type          !
  !-------------------------------------------!
  
  type model_parameters
     private
     real(dp) ::  a_data = 6357.0_dp/length_norm
     real(dp) ::  b_data = 6371.0_dp/length_norm
     real(dp) ::  c_data = 6371.0_dp/length_norm
     real(dp) ::  M_data = 5.972e24_dp/mass_norm
     real(dp) ::  g_data = 9.807_dp/acceleration_norm
     real(dp) :: Om_data = twopi/(24.0*3600.0_dp)*time_norm
     real(dp) :: I1_data = 8.012e37_dp/(mass_norm*length_norm**2)
     real(dp) :: I2_data = 8.012e37_dp/(mass_norm*length_norm**2)
     real(dp) :: I3_data = 8.038e37_dp/(mass_norm*length_norm**2)
   contains
     procedure ::  a => get_a
     procedure ::  b => get_b
     procedure ::  c => get_c
     procedure ::  g => get_g
     procedure :: Om => get_Om
     procedure :: I1 => get_I1
     procedure :: I2 => get_I2
     procedure :: I3 => get_I3
  end type model_parameters
  
  

contains

  real(dp)  function mass_spherical_model(self,n) result(M)
    class(spherical_model), intent(in) :: self
    integer(i4b), intent(in), optional :: n
    integer(i4b), parameter :: ndef = 5
    integer(i4b) :: isection,ilayer,i,nuse
    real(dp) :: r1,r2,r,w
    type(gauss_quadrature) :: quad
    if(present(n)) then
       nuse = n
    else
       nuse = ndef
    end if
    call quad%set(nuse)   
    M = 0.0_dp
    do isection = 1,self%nsections
       do ilayer = 1,self%section(isection)%nlayers
          associate(layer => self%section(isection)%layer(ilayer))
            r1 = layer%r1
            r2 = layer%r2
            call quad%trans(r1,r2)
            do i = 1,nuse
               r = quad%x(i)
               w = quad%w(i)
               M = M + 4.0_dp*pi*layer%rho(r)*r*r*w
            end do
          end associate
       end do
    end do
    return
  end function mass_spherical_model

  real(dp)  function surface_gravity_spherical_model(self,n) result(g)
    class(spherical_model), intent(in) :: self
    integer(i4b), intent(in), optional :: n
    g = self%M(n)*bigg/self%r2**2
    return
  end function surface_gravity_spherical_model



  
  !--------------------------------------------!
  !              elastic procedures            !
  !--------------------------------------------!

  
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
       iso = isotropic
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


  !----------------------------------------------------!
  !                  Maxwell procedures                !
  !----------------------------------------------------!


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
       iso = isotropic
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

  
  !----------------------------------------------------!
  !                 anelastic procedures               !
  !----------------------------------------------------!


  subroutine write_anelastic_spherical_model(self,model_file,nknot,isotropic,qisotropic)
    class(spherical_model), intent(in) :: self
    character(len=*), intent(in) :: model_file
    integer(i4b), intent(in) :: nknot
    logical, intent(in), optional :: isotropic,qisotropic

    logical :: iso,qiso
    integer(i4b) :: io,isection,ilayer,ir,nr
    real(dp) :: r1,r2,r,dr,eta,rho,vpv,vph,vsv,vsh, &
                qk,qm,qA,qC,qF,qL,qN

    iso = .false.
    if(present(isotropic)) then
       iso = isotropic
    end if

    qiso = .true.
    if(present(qisotropic)) then
       qiso = qisotropic
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
                  
               class is(spherical_solid_anelastic_layer)
               
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

                  if(qiso) then
                     qk = layer%qkappa(r)
                     qm = layer%qmu(r)                     
                  else
                     qA = layer%qA(r)
                     qC = layer%qC(r)
                     qF = layer%qF(r)
                     qL = layer%qL(r)
                     qN = layer%qN(r)                     
                  end if

               class is(spherical_fluid_anelastic_layer)

                  vpv = layer%vp(r)
                  vph = layer%vp(r)
                  vsv = 0.0_dp
                  vsh = 0.0_dp
                  eta = 1.0_dp

                  if(qiso) then
                     
                     qk = layer%qkappa(r)
                     qm = 0.0_dp
                     
                  else
                     
                     qA = layer%qkappa(r)
                     qC = layer%qkappa(r)
                     qF = layer%qkappa(r)
                     qL = 0.0_dp
                     qN = 0.0_dp
                     
                  end if
                  
               class default
                  stop 'write_anelastic_spherical_model: layer of an invalid type'
               end select

               
               if(iso) then
                  
                  write(io,'(6e18.8)') r*length_norm,      &
                                       rho*density_norm,   & 
                                       vpv*velocity_norm,  & 
                                       vsv*velocity_norm,  &
                                       qk,                 &
                                       qm
                  
               else

                  if(qiso) then


                     write(io,'(9e18.8)') r*length_norm,      &
                                          rho*density_norm,   & 
                                          vpv*velocity_norm,  & 
                                          vsv*velocity_norm,  &
                                          qk,                 &
                                          qm,                 &
                                          vph*velocity_norm,  & 
                                          vsh*velocity_norm,  &
                                          eta


                  else

                     
                     write(io,'(12e18.8)') r*length_norm,     &
                                          rho*density_norm,   & 
                                          vpv*velocity_norm,  & 
                                          vsv*velocity_norm,  &
                                          qC,                 &
                                          qL,                 &
                                          vph*velocity_norm,  & 
                                          vsh*velocity_norm,  &
                                          qA,                 &
                                          qN,                 &
                                          eta,                &
                                          qF
                     

                  end if
                  
               end if
            end do

          end associate

          
       end do
       
    end do
    close(io)
    
    return
  end subroutine write_anelastic_spherical_model
     
  function qkappa_spherical_solid_anelastic_layer(self,r) result(qkappa)
    class(spherical_solid_anelastic_layer), intent(in) :: self
    real(dp), intent(in) :: r
    real(dp) :: qkappa
    qkappa =   self%C(r)*self%qC(r) + 4.0_dp*self%A(r)*self%qA(r) &
             - 4.0_dp*self%N(r)*self%qN(r) + 4.0_dp*self%F(r)*self%qF(r)
    qkappa = qkappa/(9.0_dp*self%kappa(r))
    return
  end function qkappa_spherical_solid_anelastic_layer


  function qmu_spherical_solid_anelastic_layer(self,r) result(qmu)
    class(spherical_solid_anelastic_layer), intent(in) :: self
    real(dp), intent(in) :: r
    real(dp) :: qmu
    qmu =   self%C(r)*self%qC(r) + self%A(r)*self%qA(r) + 6.0_dp*self%L(r)*self%qL(r) &
          + 5.0_dp*self%N(r)*self%qN(r) - 2.0_dp*self%F(r)*self%qF(r)
    qmu = qmu/(15.0_dp*self%mu(r))
    return
  end function qmu_spherical_solid_anelastic_layer


  !======================================================!
  !       routines for the earth_parameters type         !
  !======================================================!

  real(dp) function get_a(self) result(a)
    class(model_parameters), intent(in) :: self
    a = self%a_data
    return
  end function get_a
  
  real(dp) function get_b(self) result(b)
    class(model_parameters), intent(in) :: self
    b = self%b_data
    return
  end function get_b


  real(dp) function get_c(self) result(c)
    class(model_parameters), intent(in) :: self
    c = self%c_data
    return
  end function get_c
  
  real(dp) function get_g(self) result(g)
    class(model_parameters), intent(in) :: self
    g = self%g_data
    return
  end function get_g
  
  real(dp) function get_Om(self) result(Om)
    class(model_parameters), intent(in) :: self
    Om = self%Om_data
    return
  end function get_Om

  real(dp) function get_I1(self) result(I1)
    class(model_parameters), intent(in) :: self
    I1 = self%I1_data
    return
  end function get_I1

  real(dp) function get_I2(self) result(I2)
    class(model_parameters), intent(in) :: self
    I2 = self%I2_data
    return
  end function get_I2

  real(dp) function get_I3(self) result(I3)
    class(model_parameters), intent(in) :: self
    I3 = self%I3_data
    return
  end function get_I3
  
  real(dp) function get_M(self) result(M)
    class(model_parameters), intent(in) :: self
    M = self%M_data
    return
  end function get_M

  subroutine set_model_parameters(self,a,b,c,M,g,Om,I1,I2,I3)
    class(model_parameters), intent(inout) :: self
    real(dp), intent(in), optional :: a,b,c,g,Om,I1,I2,I3,M
    if(present(a)) self%a_data = a
    if(present(b)) self%b_data = b
    if(present(c)) self%c_data = c
    if(present(g)) self%M_data = M
    if(present(g)) self%g_data = g    
    if(present(Om)) self%Om_data = Om
    if(present(I1)) self%I1_data = I1
    if(present(I2)) self%I2_data = I2
    if(present(I3)) self%I3_data = I3    
    return
  end subroutine set_model_parameters
  

  
  
end module module_spherical_model
