module module_spherical_model

  use module_constants  
  use module_physical_constants
  use module_error
  implicit none

  !========================================!
  !                basic types             !
  !========================================!
  
  type spherical_model
     real(dp) :: r1
     real(dp) :: r2
     integer(i4b) :: nsections = 0
     class(spherical_section), dimension(:), allocatable :: section
   contains
     procedure :: write_elastic => write_elastic_spherical_model
     procedure :: write_maxwell => write_maxwell_spherical_model
     procedure :: write_anelastic => write_anelastic_spherical_model
  end type spherical_model
  
  type spherical_section
     integer(i4b) :: nlayers
     real(dp) :: r1
     real(dp) :: r2
     class(spherical_layer), dimension(:), allocatable :: layer
  end type spherical_section

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
  


contains


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
  
end module module_spherical_model
