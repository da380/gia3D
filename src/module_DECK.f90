module module_DECK

  use module_constants
  use module_spherical_model
  use module_util, only : count_columns
  use module_interp
  implicit none


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


  type, extends(spherical_solid_anelastic_layer) :: DECK_solid_anelastic_layer
     type(interp_1D_cubic) :: rho_cubic
     type(interp_1D_cubic) :: A_cubic
     type(interp_1D_cubic) :: C_cubic
     type(interp_1D_cubic) :: F_cubic
     type(interp_1D_cubic) :: L_cubic
     type(interp_1D_cubic) :: N_cubic
     type(interp_1D_cubic) :: qA_cubic
     type(interp_1D_cubic) :: qC_cubic
     type(interp_1D_cubic) :: qF_cubic
     type(interp_1D_cubic) :: qL_cubic
     type(interp_1D_cubic) :: qN_cubic
   contains
     procedure :: rho => rho_DECK_solid_anelastic_layer
     procedure :: A => A_DECK_solid_anelastic_layer
     procedure :: C => C_DECK_solid_anelastic_layer
     procedure :: F => F_DECK_solid_anelastic_layer
     procedure :: L => L_DECK_solid_anelastic_layer
     procedure :: N => N_DECK_solid_anelastic_layer
     procedure :: qA => qA_DECK_solid_anelastic_layer
     procedure :: qC => qC_DECK_solid_anelastic_layer
     procedure :: qF => qF_DECK_solid_anelastic_layer
     procedure :: qL => qL_DECK_solid_anelastic_layer
     procedure :: qN => qN_DECK_solid_anelastic_layer
  end type DECK_solid_anelastic_layer


  type, extends(spherical_fluid_anelastic_layer) :: DECK_fluid_anelastic_layer
     type(interp_1D_cubic) :: rho_cubic
     type(interp_1D_cubic) :: kappa_cubic
     type(interp_1D_cubic) :: qkappa_cubic
   contains
     procedure :: rho => rho_DECK_fluid_anelastic_layer
     procedure :: kappa => kappa_DECK_fluid_anelastic_layer
     procedure :: qkappa => qkappa_DECK_fluid_anelastic_layer
  end type DECK_fluid_anelastic_layer


  
contains


  !---------------------------------------------------!
  !                  elastic models                   !
  !---------------------------------------------------!

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
    rho = self%rho_cubic%f(r)
    return
  end function rho_DECK_solid_elastic_layer

  function A_DECK_solid_elastic_layer(self,r) result(A)
    class(DECK_solid_elastic_layer), intent(in) :: self
    real(dp), intent(in) :: r
    real(dp) :: A,vpv    
    A = self%A_cubic%f(r)
    return
  end function A_DECK_solid_elastic_layer
  
  function C_DECK_solid_elastic_layer(self,r) result(C)
    class(DECK_solid_elastic_layer), intent(in) :: self
    real(dp), intent(in) :: r
    real(dp) :: C,rho,vsv
    C = self%C_cubic%f(r)
    return
  end function C_DECK_solid_elastic_layer

  function F_DECK_solid_elastic_layer(self,r) result(F)
    class(DECK_solid_elastic_layer), intent(in) :: self
    real(dp), intent(in) :: r
    real(dp) :: F
    F = self%F_cubic%f(r)
    return
  end function F_DECK_solid_elastic_layer

  function L_DECK_solid_elastic_layer(self,r) result(L)
    class(DECK_solid_elastic_layer), intent(in) :: self
    real(dp), intent(in) :: r
    real(dp) :: L
    L = self%L_cubic%f(r)
    return
  end function L_DECK_solid_elastic_layer

  function N_DECK_solid_elastic_layer(self,r) result(N)
    class(DECK_solid_elastic_layer), intent(in) :: self
    real(dp), intent(in) :: r
    real(dp) :: N
    N = self%N_cubic%f(r)
    return
  end function N_DECK_solid_elastic_layer

  
  function rho_DECK_fluid_elastic_layer(self,r) result(rho)
    class(DECK_fluid_elastic_layer), intent(in) :: self
    real(dp), intent(in) :: r
    real(dp) :: rho
    rho = self%rho_cubic%f(r)
    return
  end function rho_DECK_fluid_elastic_layer

  function kappa_DECK_fluid_elastic_layer(self,r) result(kappa)
    class(DECK_fluid_elastic_layer), intent(in) :: self
    real(dp), intent(in) :: r
    real(dp) :: kappa,vpv    
    kappa = self%kappa_cubic%f(r)
    return
  end function kappa_DECK_fluid_elastic_layer


  !---------------------------------------------------!
  !                  Maxwell models                   !
  !---------------------------------------------------!

  
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
    rho = self%rho_cubic%f(r)
    return
  end function rho_DECK_maxwell_layer

  function A_DECK_maxwell_layer(self,r) result(A)
    class(DECK_maxwell_layer), intent(in) :: self
    real(dp), intent(in) :: r
    real(dp) :: A,vpv    
    A = self%A_cubic%f(r)
    return
  end function A_DECK_maxwell_layer
  
  function C_DECK_maxwell_layer(self,r) result(C)
    class(DECK_maxwell_layer), intent(in) :: self
    real(dp), intent(in) :: r
    real(dp) :: C,rho,vsv
    C = self%C_cubic%f(r)
    return
  end function C_DECK_maxwell_layer

  function F_DECK_maxwell_layer(self,r) result(F)
    class(DECK_maxwell_layer), intent(in) :: self
    real(dp), intent(in) :: r
    real(dp) :: F
    F = self%F_cubic%f(r)
    return
  end function F_DECK_maxwell_layer

  function L_DECK_maxwell_layer(self,r) result(L)
    class(DECK_maxwell_layer), intent(in) :: self
    real(dp), intent(in) :: r
    real(dp) :: L
    L = self%L_cubic%f(r)
    return
  end function L_DECK_maxwell_layer

  function N_DECK_maxwell_layer(self,r) result(N)
    class(DECK_maxwell_layer), intent(in) :: self
    real(dp), intent(in) :: r
    real(dp) :: N
    N = self%N_cubic%f(r)
    return
  end function N_DECK_maxwell_layer

  function eta_DECK_maxwell_layer(self,r) result(eta)
    class(DECK_maxwell_layer), intent(in) :: self
    real(dp), intent(in) :: r
    real(dp) :: eta
    eta = self%eta_cubic%f(r)
    return
  end function eta_DECK_maxwell_layer



  !---------------------------------------------------!
  !                 anelastic models                  !
  !---------------------------------------------------!

  function anelastic_DECK(file) result(model)
    character(len=*), intent(in) :: file
    type(spherical_model) :: model

    logical :: isotropic,change,qisotropic
    integer(i4b) :: io,iknot,nknot,ios,ncol,isec,i1,i2,ilay,i11,i22
    integer(i4b), dimension(:,:), allocatable :: sind
    real(dp) :: tmp,qk,qm,k,m,a,c,f,l,n
    real(dp), dimension(:), allocatable :: r,rho,vpv,vsv,vph,vsh,eta, &
                                           qA,qC,qF,qL,qN

    !------------------------------------!
    !           read in the file         !
    !------------------------------------!
    open(newunit = io,file = trim(file),iostat = ios)
    if(ios /= 0) stop 'anelastic_DECK: problem opening model file'
    ncol = count_columns(io)
    if(ncol == 6) then
       isotropic = .true.
       qisotropic = .true.
    else if(ncol == 9) then
       isotropic = .false.
       qisotropic = .true.
    else if(ncol == 12) then
       isotropic = .false.
       qisotropic = .false.
    else
       stop 'anelastic_DECK: unexpected format'
    end if
    nknot = 0
    do
       read(io,*,iostat = ios) tmp
       if(ios /= 0) exit
       nknot = nknot + 1
    end do
    allocate(r(nknot),rho(nknot),vpv(nknot),vsv(nknot),vph(nknot),vsh(nknot), &
             eta(nknot),qA(nknot),qC(nknot),qF(nknot),qL(nknot),qN(nknot))
    rewind(io)    
    do iknot = 1,nknot
       if(isotropic .and. qisotropic) then
          read(io,*) r(iknot),rho(iknot),vpv(iknot),vsv(iknot),qA(iknot),qL(iknot)
          vph(iknot) = vpv(iknot)
          vsh(iknot) = vsv(iknot)
          eta(iknot) = 1.0_dp
       else if(.not. isotropic .and. qisotropic) then
          read(io,*) r(iknot),rho(iknot),vpv(iknot),vsv(iknot),qA(iknot),qL(iknot), &
                    vph(iknot),vsh(iknot),eta(iknot)
       else
          read(io,*) r(iknot),rho(iknot),vpv(iknot),vsv(iknot),qC(iknot),qL(iknot), &
                    vph(iknot),vsh(iknot),qA(iknot),qN(iknot),eta(iknot),qF(iknot)
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

    ! set the quality factors if needed
    if(qisotropic) then

       qN = qL
       do iknot = 1,nknot

          qk = qA(iknot)
          qm = qL(iknot)

          a = vph(iknot)
          c = vpv(iknot)
          f = eta(iknot)
          l = vsv(iknot)
          n = vsh(iknot)

          k = (c + 4.0_dp*a - 4.0_dp*n + 4.0_dp*f)/9.0_dp
          m = (c + a + 6.0_dp*l + 5.0_dp*n - 2.0_dp*f)/15.0_dp

          qA(iknot) = (k*qk + 4.0_dp*m*qm/3.0_dp)/a
          qC(iknot) = (k*qk + 4.0_dp*m*qm/3.0_dp)/c
          qF(iknot) = (k*qk - 2.0_dp*m*qm/3.0_dp)/f

       end do
       
    end if
    
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
            allocate(DECK_solid_anelastic_layer:: section%layer(nlayers))
         else
            allocate(DECK_fluid_anelastic_layer:: section%layer(nlayers))
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
                 type is(DECK_solid_anelastic_layer)
                    call layer%rho_cubic%set(r(i11:i22),rho(i11:i22))
                    call layer%A_cubic%set(r(i11:i22),vph(i11:i22))
                    call layer%C_cubic%set(r(i11:i22),vpv(i11:i22))
                    call layer%F_cubic%set(r(i11:i22),eta(i11:i22))
                    call layer%L_cubic%set(r(i11:i22),vsv(i11:i22))
                    call layer%N_cubic%set(r(i11:i22),vsh(i11:i22))
                    call layer%qA_cubic%set(r(i11:i22),qA(i11:i22))
                    call layer%qC_cubic%set(r(i11:i22),qC(i11:i22))
                    call layer%qF_cubic%set(r(i11:i22),qF(i11:i22))
                    call layer%qL_cubic%set(r(i11:i22),qL(i11:i22))
                    call layer%qN_cubic%set(r(i11:i22),qN(i11:i22))
                 type is(DECK_fluid_anelastic_layer)                    
                    call layer%rho_cubic%set(r(i11:i22),rho(i11:i22))
                    call layer%kappa_cubic%set(r(i11:i22),vpv(i11:i22))
                    call layer%qkappa_cubic%set(r(i11:i22),qA(i11:i22))
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
           type is(DECK_solid_anelastic_layer)
              call layer%rho_cubic%set(r(i11:i22),rho(i11:i22))
              call layer%A_cubic%set(r(i11:i22),vph(i11:i22))
              call layer%C_cubic%set(r(i11:i22),vpv(i11:i22))
              call layer%F_cubic%set(r(i11:i22),eta(i11:i22))
              call layer%L_cubic%set(r(i11:i22),vsv(i11:i22))
              call layer%N_cubic%set(r(i11:i22),vsh(i11:i22))
              call layer%qA_cubic%set(r(i11:i22),qA(i11:i22))
              call layer%qC_cubic%set(r(i11:i22),qC(i11:i22))
              call layer%qF_cubic%set(r(i11:i22),qF(i11:i22))
              call layer%qL_cubic%set(r(i11:i22),qL(i11:i22))
              call layer%qN_cubic%set(r(i11:i22),qN(i11:i22))
           type is(DECK_fluid_anelastic_layer)                    
              call layer%rho_cubic%set(r(i11:i22),rho(i11:i22))
              call layer%kappa_cubic%set(r(i11:i22),vpv(i11:i22))
              call layer%qkappa_cubic%set(r(i11:i22),qA(i11:i22))
           end select
         end associate
         
         
       end associate
       
    end do
    
    
    return
  end function anelastic_DECK


  function rho_DECK_solid_anelastic_layer(self,r) result(rho)
    class(DECK_solid_anelastic_layer), intent(in) :: self
    real(dp), intent(in) :: r
    real(dp) :: rho
    rho = self%rho_cubic%f(r)
    return
  end function rho_DECK_solid_anelastic_layer

  function A_DECK_solid_anelastic_layer(self,r) result(A)
    class(DECK_solid_anelastic_layer), intent(in) :: self
    real(dp), intent(in) :: r
    real(dp) :: A,vpv    
    A = self%A_cubic%f(r)
    return
  end function A_DECK_solid_anelastic_layer
  
  function C_DECK_solid_anelastic_layer(self,r) result(C)
    class(DECK_solid_anelastic_layer), intent(in) :: self
    real(dp), intent(in) :: r
    real(dp) :: C,rho,vsv
    C = self%C_cubic%f(r)
    return
  end function C_DECK_solid_anelastic_layer

  function F_DECK_solid_anelastic_layer(self,r) result(F)
    class(DECK_solid_anelastic_layer), intent(in) :: self
    real(dp), intent(in) :: r
    real(dp) :: F
    F = self%F_cubic%f(r)
    return
  end function F_DECK_solid_anelastic_layer

  function L_DECK_solid_anelastic_layer(self,r) result(L)
    class(DECK_solid_anelastic_layer), intent(in) :: self
    real(dp), intent(in) :: r
    real(dp) :: L
    L = self%L_cubic%f(r)
    return
  end function L_DECK_solid_anelastic_layer

  function N_DECK_solid_anelastic_layer(self,r) result(N)
    class(DECK_solid_anelastic_layer), intent(in) :: self
    real(dp), intent(in) :: r
    real(dp) :: N
    N = self%N_cubic%f(r)
    return
  end function N_DECK_solid_anelastic_layer

  
  function rho_DECK_fluid_anelastic_layer(self,r) result(rho)
    class(DECK_fluid_anelastic_layer), intent(in) :: self
    real(dp), intent(in) :: r
    real(dp) :: rho
    rho = self%rho_cubic%f(r)
    return
  end function rho_DECK_fluid_anelastic_layer

  function kappa_DECK_fluid_anelastic_layer(self,r) result(kappa)
    class(DECK_fluid_anelastic_layer), intent(in) :: self
    real(dp), intent(in) :: r
    real(dp) :: kappa,vpv    
    kappa = self%kappa_cubic%f(r)
    return
  end function kappa_DECK_fluid_anelastic_layer


  function qA_DECK_solid_anelastic_layer(self,r) result(qA)
    class(DECK_solid_anelastic_layer), intent(in) :: self
    real(dp), intent(in) :: r
    real(dp) :: qA,vpv    
    qA = self%qA_cubic%f(r)
    return
  end function qA_DECK_solid_anelastic_layer
  
  function qC_DECK_solid_anelastic_layer(self,r) result(qC)
    class(DECK_solid_anelastic_layer), intent(in) :: self
    real(dp), intent(in) :: r
    real(dp) :: qC,rho,vsv
    qC = self%qC_cubic%f(r)
    return
  end function qC_DECK_solid_anelastic_layer

  function qF_DECK_solid_anelastic_layer(self,r) result(qF)
    class(DECK_solid_anelastic_layer), intent(in) :: self
    real(dp), intent(in) :: r
    real(dp) :: qF
    qF = self%qF_cubic%f(r)
    return
  end function qF_DECK_solid_anelastic_layer

  function qL_DECK_solid_anelastic_layer(self,r) result(qL)
    class(DECK_solid_anelastic_layer), intent(in) :: self
    real(dp), intent(in) :: r
    real(dp) :: qL
    qL = self%qL_cubic%f(r)
    return
  end function qL_DECK_solid_anelastic_layer

  function qN_DECK_solid_anelastic_layer(self,r) result(qN)
    class(DECK_solid_anelastic_layer), intent(in) :: self
    real(dp), intent(in) :: r
    real(dp) :: qN
    qN = self%qN_cubic%f(r)
    return
  end function qN_DECK_solid_anelastic_layer
  
  function qkappa_DECK_fluid_anelastic_layer(self,r) result(qkappa)
    class(DECK_fluid_anelastic_layer), intent(in) :: self
    real(dp), intent(in) :: r
    real(dp) :: qkappa,vpv    
    qkappa = self%qkappa_cubic%f(r)
    return
  end function qkappa_DECK_fluid_anelastic_layer

  
  
end module module_DECK
