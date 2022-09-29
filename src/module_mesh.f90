module module_mesh

  use module_constants
  use module_physical_constants, only : bigg
  use module_spherical_model
  use module_quadrature
  use module_special_functions
  implicit none

  type spherical_layer_mesh
     logical :: bottom = .false.
     logical :: top = .false.
     integer(i4b) :: ngll
     integer(i4b) :: nspec
     real(dp) :: r1
     real(dp) :: r2
     real(dp), dimension(:,:), allocatable :: r
     real(dp), dimension(:), allocatable :: w
     real(dp), dimension(:), allocatable :: jac
     real(dp), dimension(:,:), allocatable :: hp
     real(dp), dimension(:,:), allocatable :: rho
     real(dp), dimension(:,:), allocatable :: g
     real(dp), dimension(:,:), allocatable :: ep
  end type spherical_layer_mesh
  

  type spherical_section_mesh
     integer(i4b) :: nlayers
     real(dp) :: r1
     real(dp) :: r2
     class(spherical_layer_mesh), dimension(:), allocatable :: layer     
  end type spherical_section_mesh

  
  type spherical_model_mesh
     integer(i4b) :: nsections
     real(dp) :: r1
     real(dp) :: r2
     type(spherical_section_mesh), dimension(:), allocatable :: section
   contains
     procedure :: print_summary => print_mesh_summary
  end type spherical_model_mesh
  



  type, extends(spherical_layer_mesh) ::  spherical_solid_elastic_layer_mesh
     logical :: below_fluid = .false.
     logical :: above_fluid = .false.
     real(dp) :: below_fluid_density = 0.0_dp
     real(dp) :: above_fluid_density = 0.0_dp
     real(dp), dimension(:,:), allocatable :: A
     real(dp), dimension(:,:), allocatable :: C
     real(dp), dimension(:,:), allocatable :: F
     real(dp), dimension(:,:), allocatable :: L
     real(dp), dimension(:,:), allocatable :: N
     real(dp), dimension(:,:), allocatable :: kappa
     real(dp), dimension(:,:), allocatable :: mu
  end type spherical_solid_elastic_layer_mesh


  type, extends(spherical_layer_mesh) ::  spherical_fluid_elastic_layer_mesh
     logical :: below_solid = .false.
     logical :: above_solid = .false.
     real(dp), dimension(:,:), allocatable :: drho
     real(dp), dimension(:,:), allocatable :: kappa
  end type spherical_fluid_elastic_layer_mesh
  

  type, extends(spherical_solid_elastic_layer_mesh) ::  spherical_solid_anelastic_layer_mesh     
     real(dp), dimension(:,:), allocatable :: qA
     real(dp), dimension(:,:), allocatable :: qC
     real(dp), dimension(:,:), allocatable :: qF
     real(dp), dimension(:,:), allocatable :: qL
     real(dp), dimension(:,:), allocatable :: qN
     real(dp), dimension(:,:), allocatable :: qkappa
     real(dp), dimension(:,:), allocatable :: qmu
  end type spherical_solid_anelastic_layer_mesh

  
  type, extends(spherical_fluid_elastic_layer_mesh) ::  spherical_fluid_anelastic_layer_mesh     
     real(dp), dimension(:,:), allocatable :: qkappa
  end type spherical_fluid_anelastic_layer_mesh

  
  type, extends(spherical_solid_elastic_layer_mesh) ::  spherical_maxwell_layer_mesh     
     real(dp), dimension(:,:), allocatable :: eta
  end type spherical_maxwell_layer_mesh


  type boolean_array_layer
     integer(i4b) :: nvar
     integer(i4b) :: ispec1
     integer(i4b), dimension(:,:,:), allocatable :: data
   contains
     procedure :: get => get_boolean_array_layer
  end type boolean_array_layer


  type boolean_array_section
     integer(i4b) :: ilayer1
     type(boolean_array_layer), dimension(:), allocatable :: layer
   contains
     procedure :: get => get_boolean_array_section
  end type boolean_array_section


  type boolean_array
     integer(i4b) :: ndim     
     integer(i4b) :: ngll
     integer(i4b) :: isection1
     type(boolean_array_section), dimension(:), allocatable :: section
   contains
     procedure :: get => get_boolean_array
  end type boolean_array


  interface spherical_mesh
     procedure :: make_spherical_mesh
  end interface spherical_mesh


  
contains
  

  !===========================================================!
  !                routines for mesh building                 !
  !===========================================================!
  
  function make_spherical_mesh(ngll,model,drmax,Om) result(mesh)
    integer(i4b), intent(in) :: ngll
    type(spherical_model), intent(in) :: model
    real(dp), intent(in) :: drmax
    real(dp), intent(in), optional :: Om
    type(spherical_model_mesh) :: mesh

    logical :: fluid
    integer(i4b) :: ilayer,nlayers,isection
    class(spherical_layer_mesh), allocatable :: mesh_tmp


    mesh%nsections = model%nsections
    mesh%r1 = model%r1
    mesh%r2 = model%r2
    allocate(mesh%section(mesh%nsections))

    do isection = 1,mesh%nsections
       nlayers = model%section(isection)%nlayers
       mesh%section(isection)%nlayers = nlayers
       associate(layer => model%section(isection)%layer(1))
         select type(layer)
         class is(spherical_layer)
            allocate(spherical_layer_mesh :: mesh%section(isection)%layer(nlayers))
         class is(spherical_solid_elastic_layer)
            allocate(spherical_solid_elastic_layer_mesh :: mesh%section(isection)%layer(nlayers))
         class is(spherical_fluid_elastic_layer)
            allocate(spherical_fluid_elastic_layer_mesh :: mesh%section(isection)%layer(nlayers))
         class is(spherical_solid_anelastic_layer)
            allocate(spherical_solid_anelastic_layer_mesh :: mesh%section(isection)%layer(nlayers))
         class is(spherical_fluid_anelastic_layer)
            allocate(spherical_fluid_anelastic_layer_mesh :: mesh%section(isection)%layer(nlayers))          
         class is(spherical_maxwell_layer)
            allocate(spherical_maxwell_layer_mesh :: mesh%section(isection)%layer(nlayers))          
         end select
       end associate
    end do
      
    do isection = 1,mesh%nsections
       mesh%section(isection)%r1 = model%section(isection)%r1
       mesh%section(isection)%r2 = model%section(isection)%r2       
       do ilayer = 1,mesh%section(isection)%nlayers

          associate( layer => model%section(isection)%layer(ilayer), & 
                     layer_mesh =>  mesh%section(isection)%layer(ilayer))
            select type(layer)
             

            class is(spherical_layer)
               select type(layer_mesh)
               class is(spherical_layer_mesh)
                  call make_spherical_layer_mesh(ngll,layer,drmax,layer_mesh)
               end select
               
            class is(spherical_solid_elastic_layer)
               select type(layer_mesh)
               class is(spherical_solid_elastic_layer_mesh)
                  call make_spherical_solid_elastic_layer_mesh(ngll,layer,drmax,layer_mesh)
               end select

               
            class is(spherical_fluid_elastic_layer)
               select type(layer_mesh)
               class is(spherical_fluid_elastic_layer_mesh)
                  call make_spherical_fluid_elastic_layer_mesh(ngll,layer,drmax,layer_mesh)
               end select               


            class is(spherical_solid_anelastic_layer)
               select type(layer_mesh)
               class is(spherical_solid_anelastic_layer_mesh)
                  call make_spherical_solid_anelastic_layer_mesh(ngll,layer,drmax,layer_mesh)
               end select

               
            class is(spherical_fluid_anelastic_layer)
               select type(layer_mesh)
               class is(spherical_fluid_anelastic_layer_mesh)
                  call make_spherical_fluid_anelastic_layer_mesh(ngll,layer,drmax,layer_mesh)
               end select               
               

            class is(spherical_maxwell_layer)
               select type(layer_mesh)
               class is(spherical_maxwell_layer_mesh)
                  call make_spherical_maxwell_layer_mesh(ngll,layer,drmax,layer_mesh)
               end select

               
            end select
            
            
          end associate
          
       end do
       
    end do

    call calculate_gravity_and_ellipticity(mesh,Om)
    
    ! flag the top and bottom layers
    mesh%section(1)%layer(1)%bottom = .true.
    mesh%section(mesh%nsections)%layer(mesh%section(mesh%nsections)%nlayers)%top =  .true.


    ! label solid-fluid and fluid-solid boundaries
    do isection = 2,mesh%nsections       
       associate(layer => mesh%section(isection)%layer(1), &
                 blayer => mesh%section(isection-1)%layer(mesh%section(isection-1)%nlayers))

         select type(layer)

         class is(spherical_solid_elastic_layer_mesh)


            
            select type(blayer)
               
            class is(spherical_fluid_elastic_layer_mesh)

               layer%below_fluid = .true.
               blayer%above_solid = .true.
               layer%below_fluid_density = blayer%rho(blayer%ngll,blayer%nspec)

            end select
            
         class is(spherical_fluid_elastic_layer_mesh)

            select type(blayer)

            class is(spherical_solid_elastic_layer_mesh)

               layer%below_solid = .true.
               blayer%above_fluid = .true.
               blayer%above_fluid_density = layer%rho(1,1)
               
            end select

         end select
         
       end associate
    end do

    

    
    return
  end function make_spherical_mesh


  subroutine print_mesh_summary(mesh)
    class(spherical_model_mesh), intent(in) :: mesh

    integer(i4b) :: isection,ilayer,nsections,nlayers
    do isection = 1,mesh%nsections

       print *, '====================================='
       print *, 'Section = ',isection
       print *, '====================================='
       
       nlayers = mesh%section(isection)%nlayers
       do ilayer = 1,nlayers
          
          associate(layer => mesh%section(isection)%layer(ilayer))

            print *, ''
            print *, 'Layer = ',ilayer
            print *, 'Bottom = ',layer%bottom
            print *, 'Top = ',layer%top

            select type(layer)

            class is(spherical_solid_elastic_layer_mesh)

               print *, 'Solid layer '
               print *, 'Fluid below = ',layer%below_fluid
               print *, 'Fluid above = ',layer%above_fluid
               if(layer%below_fluid) print *, 'Fluid density below = ', layer%below_fluid_density*density_norm
               if(layer%above_fluid) print *, 'Fluid density above = ', layer%above_fluid_density*density_norm

            class is(spherical_fluid_elastic_layer_mesh)

               print *, 'Fluid layer'
               print *, 'Solid below = ',layer%below_solid
               print *, 'Solid above = ',layer%above_solid
               
            end select

            
            
          end associate
          
       end do
       print *, ''

    end do

    
  end subroutine print_mesh_summary
  

  
  subroutine make_spherical_layer_mesh(ngll,layer,drmax,mesh)
    integer(i4b), intent(in) :: ngll
    class(spherical_layer), intent(in) :: layer
    real(dp), intent(in) :: drmax
    type(spherical_layer_mesh), intent(out) :: mesh
    integer(i4b) :: ispec,nspec,inode,jspec
    real(dp) :: r1,r2,r11,r22,dr,r
    real(dp), dimension(ngll) :: h
    type(gauss_lobatto_quadrature) :: quad
    
    ! store basic data
    mesh%ngll = ngll
    mesh%r1 = layer%r1
    mesh%r2 = layer%r2

    ! set the quadrature scheme
    call quad%set(ngll)
    mesh%w = quad%w

    ! store the lagrange derivatives
    allocate(mesh%hp(ngll,ngll))
    do inode = 1,ngll
       call lagrange_polynomial(quad%x(inode),ngll,quad%x,h,mesh%hp(inode,:))
    end do

    
    ! work out number of spectral elements
    r1 = layer%r1
    r2 = layer%r2
    nspec = (r2-r1)/drmax
    nspec = max(1,nspec)
    mesh%nspec = nspec
    allocate(mesh%r(ngll,mesh%nspec))
    allocate(mesh%jac(mesh%nspec))
    allocate(mesh%rho(ngll,mesh%nspec))
    allocate(mesh%g(ngll,mesh%nspec))
    allocate(mesh%ep(ngll,mesh%nspec))

    ! initialise values for g and ep
    mesh%g = 0.0_dp
    mesh%ep = 0.0_dp
    
    dr = (r2-r1)/nspec
    r11 = r1
    do ispec = 1,nspec    
       r22 = r11+dr
       do inode = 1,ngll
          r = r11 + 0.5_dp*(quad%x(inode)+1.0_dp)*(r22-r11)
          mesh%r(inode,ispec) = r
          mesh%rho(inode,ispec) = layer%rho(r)
       end do
       mesh%jac(ispec) = 0.5_dp*(r22-r11)
       r11 = r22
    end do
    
    return
  end subroutine make_spherical_layer_mesh

    
  subroutine make_spherical_solid_elastic_layer_mesh(ngll,layer,drmax,mesh)
    integer(i4b), intent(in) :: ngll
    class(spherical_solid_elastic_layer), intent(in) :: layer
    real(dp), intent(in) :: drmax
    type(spherical_solid_elastic_layer_mesh), intent(out) :: mesh

    integer(i4b) :: ispec,nspec,inode,jspec
    real(dp) :: r
        
    call make_spherical_layer_mesh(ngll,layer,drmax,mesh%spherical_layer_mesh)        
    allocate(mesh%A(mesh%ngll,mesh%nspec))
    allocate(mesh%C(mesh%ngll,mesh%nspec))
    allocate(mesh%F(mesh%ngll,mesh%nspec))
    allocate(mesh%L(mesh%ngll,mesh%nspec))
    allocate(mesh%N(mesh%ngll,mesh%nspec))
    allocate(mesh%kappa(mesh%ngll,mesh%nspec))
    allocate(mesh%mu(mesh%ngll,mesh%nspec))
    do ispec = 1,mesh%nspec
       do inode = 1,mesh%ngll
          r = mesh%r(inode,ispec)
          mesh%A(inode,ispec) = layer%A(r)
          mesh%C(inode,ispec) = layer%C(r)
          mesh%F(inode,ispec) = layer%F(r)
          mesh%L(inode,ispec) = layer%L(r)
          mesh%N(inode,ispec) = layer%N(r)
          mesh%kappa(inode,ispec) = layer%kappa(r)
          mesh%mu(inode,ispec) = layer%mu(r)          
       end do
    end do
    
    
    return
  end subroutine make_spherical_solid_elastic_layer_mesh



  subroutine make_spherical_fluid_elastic_layer_mesh(ngll,layer,drmax,mesh)
    integer(i4b), intent(in) :: ngll
    class(spherical_fluid_elastic_layer), intent(in) :: layer
    real(dp), intent(in) :: drmax
    type(spherical_fluid_elastic_layer_mesh), intent(out) :: mesh

    integer(i4b) :: ispec,nspec,inode,jnode
    real(dp) :: r
        
    call make_spherical_layer_mesh(ngll,layer,drmax,mesh%spherical_layer_mesh)        
    allocate(mesh%kappa(mesh%ngll,mesh%nspec))
    allocate(mesh%drho(mesh%ngll,mesh%nspec))
    do ispec = 1,mesh%nspec
       do inode = 1,mesh%ngll
          r = mesh%r(inode,ispec)
          mesh%kappa(inode,ispec) = layer%kappa(r)
          mesh%drho(inode,ispec) = 0.0_dp
          do jnode = 1,mesh%ngll
             mesh%drho(inode,ispec) = mesh%drho(inode,ispec) + mesh%rho(jnode,ispec)  &
                                                             * mesh%hp(inode,jnode)   &
                                                             / mesh%jac(ispec)
          end do          
       end do
    end do
    
    return
  end subroutine make_spherical_fluid_elastic_layer_mesh


  subroutine make_spherical_solid_anelastic_layer_mesh(ngll,layer,drmax,mesh)
    integer(i4b), intent(in) :: ngll
    class(spherical_solid_anelastic_layer), intent(in) :: layer
    real(dp), intent(in) :: drmax
    type(spherical_solid_anelastic_layer_mesh), intent(out) :: mesh

    integer(i4b) :: ispec,nspec,inode,jspec
    real(dp) :: r
        
    call make_spherical_solid_elastic_layer_mesh(ngll,layer,drmax,mesh%spherical_solid_elastic_layer_mesh)
    allocate(mesh%qA(mesh%ngll,mesh%nspec))
    allocate(mesh%qC(mesh%ngll,mesh%nspec))
    allocate(mesh%qF(mesh%ngll,mesh%nspec))
    allocate(mesh%qL(mesh%ngll,mesh%nspec))
    allocate(mesh%qN(mesh%ngll,mesh%nspec))
    allocate(mesh%qkappa(mesh%ngll,mesh%nspec))
    allocate(mesh%qmu(mesh%ngll,mesh%nspec))
    do ispec = 1,mesh%nspec
       do inode = 1,mesh%ngll
          r = mesh%r(inode,ispec)
          mesh%qA(inode,ispec) = layer%qA(r)
          mesh%qC(inode,ispec) = layer%qC(r)
          mesh%qF(inode,ispec) = layer%qF(r)
          mesh%qL(inode,ispec) = layer%qL(r)
          mesh%qN(inode,ispec) = layer%qN(r)
          mesh%qkappa(inode,ispec) = layer%qkappa(r)
          mesh%qmu(inode,ispec) = layer%qmu(r)          
       end do
    end do
    
    
    return
  end subroutine make_spherical_solid_anelastic_layer_mesh


  
  subroutine make_spherical_fluid_anelastic_layer_mesh(ngll,layer,drmax,mesh)
    integer(i4b), intent(in) :: ngll
    class(spherical_fluid_anelastic_layer), intent(in) :: layer
    real(dp), intent(in) :: drmax
    type(spherical_fluid_anelastic_layer_mesh), intent(out) :: mesh

    integer(i4b) :: ispec,nspec,inode,jspec
    real(dp) :: r
        
    call make_spherical_fluid_elastic_layer_mesh(ngll,layer,drmax,mesh%spherical_fluid_elastic_layer_mesh)
    allocate(mesh%qkappa(mesh%ngll,mesh%nspec))
    do ispec = 1,mesh%nspec
       do inode = 1,mesh%ngll
          r = mesh%r(inode,ispec)
          mesh%qkappa(inode,ispec) = layer%qkappa(r)
       end do
    end do
        
    return
  end subroutine make_spherical_fluid_anelastic_layer_mesh

  

  subroutine make_spherical_maxwell_layer_mesh(ngll,layer,drmax,mesh)
    integer(i4b), intent(in) :: ngll
    class(spherical_maxwell_layer), intent(in) :: layer
    real(dp), intent(in) :: drmax
    type(spherical_maxwell_layer_mesh), intent(out) :: mesh

    integer(i4b) :: ispec,nspec,inode,jspec
    real(dp) :: r
        
    call make_spherical_solid_elastic_layer_mesh(ngll,layer,drmax,mesh%spherical_solid_elastic_layer_mesh)
    allocate(mesh%eta(mesh%ngll,mesh%nspec))
    do ispec = 1,mesh%nspec
       do inode = 1,mesh%ngll
          r = mesh%r(inode,ispec)
          mesh%eta(inode,ispec) = layer%eta(r)          
       end do
    end do
    
    
    return
  end subroutine make_spherical_maxwell_layer_mesh



  !===========================================================!
  !               routines for Boolean arrays                 !
  !===========================================================!

  
  
  integer(i4b) function get_boolean_array_layer(self,ivar,inode,ispec) result(i)
    class(boolean_array_layer), intent(in) :: self
    integer(i4b), intent(in) :: ivar,inode,ispec
    i = self%data(ivar,inode,ispec)
    return
  end function get_boolean_array_layer


  integer(i4b) function get_boolean_array_section(self,ivar,inode,ispec,ilayer) result(i)
    class(boolean_array_section), intent(in) :: self
    integer(i4b), intent(in) :: ivar,inode,ispec,ilayer
    i = self%layer(ilayer)%get(ivar,inode,ispec)
    return
  end function get_boolean_array_section


  integer(i4b) function get_boolean_array(self,ivar,inode,ispec,ilayer,isection) result(i)
    class(boolean_array), intent(in) :: self
    integer(i4b), intent(in) :: ivar,inode,ispec,ilayer,isection
    i = self%section(isection)%get(ivar,inode,ispec,ilayer)
    return
  end function get_boolean_array


  type(boolean_array) function build_boolean_simple(mesh) result(ibool)
    type(spherical_model_mesh), intent(in) :: mesh
    integer(i4b) :: isection,ilayer,ngll,nspec,count,ispec,inode,nsections,nlayers
    count = 0
    nsections = mesh%nsections
    ibool%isection1 = 1
    ibool%ngll = 0    
    allocate(ibool%section(nsections))
    do isection = 1,nsections
       nlayers = mesh%section(isection)%nlayers
       ibool%section(isection)%ilayer1 = 1
       allocate(ibool%section(isection)%layer(nlayers))
       do ilayer = 1,nlayers
          ibool%section(isection)%layer(ilayer)%nvar = 1
          ibool%section(isection)%layer(ilayer)%ispec1 = 1
          ngll = mesh%section(isection)%layer(ilayer)%ngll
          if(ngll > ibool%ngll) ibool%ngll = ngll
          nspec = mesh%section(isection)%layer(ilayer)%nspec
          allocate(ibool%section(isection)%layer(ilayer)%data(1,ngll,nspec))

          do ispec = 1,nspec
             do inode = 1,ngll
                count = count + 1
                ibool%section(isection)%layer(ilayer)%data(1,inode,ispec) = count
             end do
             count = count-1
          end do
       end do
    end do
    ibool%ndim = count+1
 
    return
  end function build_boolean_simple



  !================================================================================!
  !                     routine to calculate gravity and ellipticity               !
  !================================================================================!


  subroutine calculate_gravity_and_ellipticity(mesh,Om)
    class(spherical_model_mesh), intent(inout) :: mesh
    real(dp), optional :: Om

    integer(i4b) :: ndim,kd,ldab,isection,ilayer,inode, &
                    ispec,jnode,knode,ngll,nspec,i,j,k,info,l
    real(dp) :: tmp,rtmp,gtmp,rotfac
    real(dp), dimension(:), allocatable :: drho
    real(dp), dimension(:,:), allocatable :: a,b
    type(boolean_array) :: ibool



    !=======================================================!
    !                   gravity calculation                 !
    !=======================================================!
    
    ! make the Boolean array
    ibool =  build_boolean_simple(mesh)

    ! set the dimensions
    ndim = ibool%ndim
    kd   = ibool%ngll-1
    ldab = kd+1

    ! allocate arrays for matrix and RHS
    allocate(a(ldab,ndim),b(ndim,1))


    ! build up the system matrix and RHS
    a = 0.0_dp
    b = 0.0_dp
    do isection = 1,mesh%nsections
       do ilayer = 1,mesh%section(isection)%nlayers
          associate(ngll  => mesh%section(isection)%layer(ilayer)%ngll,  &
                    nspec => mesh%section(isection)%layer(ilayer)%nspec, &
                    w     => mesh%section(isection)%layer(ilayer)%w,     &
                    hp    => mesh%section(isection)%layer(ilayer)%hp,    &
                    jac   => mesh%section(isection)%layer(ilayer)%jac,   &
                    r     => mesh%section(isection)%layer(ilayer)%r,     &
                    ibool => ibool%section(isection)%layer(ilayer),      &
                    rho => mesh%section(isection)%layer(ilayer)%rho)
            do ispec = 1,nspec
               do inode = 1,ngll
                  i = ibool%get(1,inode,ispec)
                  do jnode = inode,ngll
                     j = ibool%get(1,jnode,ispec)
                     k = kd + 1 + i - j
                     do knode = 1,ngll
                        a(k,j) = a(k,j) + hp(knode,inode)*hp(knode,jnode) &
                                        * r(knode,ispec)**2*w(knode)/jac(ispec)                      
                     end do
                  end do
                  b(i,1) = b(i,1) + rho(inode,ispec)*r(inode,ispec)**2 &
                                  * w(inode)*jac(ispec)
               end do
            end do
          end associate
       end do
    end do


    ! add in DTN term at the surface
    k = kd+1
    a(k,ndim) = a(k,ndim) + mesh%r2

    ! scale the RHS
    b = -4.0_dp*pi*bigg*b
    
    ! compute the factorisation
    call dpbtrf('U',ndim,kd,a,ldab,info)
    call check(info == 0,'calculate_gravity_and_ellipticity','problem with factorisation')

    ! solve the linear system
    call dpbtrs	('U',ndim,kd,1,a,ldab,b,ndim,info)
    call check(info == 0,'calculate_gravity_and_ellipticity','problem with substitution')

    ! convert from potential to gravity
    do isection = 1,mesh%nsections
       do ilayer = 1,mesh%section(isection)%nlayers
          associate(ngll  => mesh%section(isection)%layer(ilayer)%ngll,  &
                    nspec => mesh%section(isection)%layer(ilayer)%nspec, &
                    hp    => mesh%section(isection)%layer(ilayer)%hp,    &
                    jac   => mesh%section(isection)%layer(ilayer)%jac,   &
                    r     => mesh%section(isection)%layer(ilayer)%r,     &
                    rho   => mesh%section(isection)%layer(ilayer)%rho,   &
                    g     => mesh%section(isection)%layer(ilayer)%g,     &
                    ibool => ibool%section(isection)%layer(ilayer))
            do ispec = 1,nspec
               do inode = 1,ngll
                  g(inode,ispec) = 0.0_dp
                  do jnode = 1,ngll
                     j = ibool%get(1,jnode,ispec)
                     g(inode,ispec) = g(inode,ispec) + b(j,1)*hp(inode,jnode)/jac(ispec)
                  end do
               end do
            end do
          end associate
       end do
    end do
    

    !=======================================================!
    !                   ellipticity calculation             !
    !=======================================================!

    if(.not.present(Om)) return
    rotfac = sqrt(fourpi/5.0_dp)*Om*Om/3.0_dp

    
    ! build up the system matrix and RHS
    l = 2
    a = 0.0_dp
    b = 0.0_dp
    do isection = 1,mesh%nsections       
       do ilayer = 1,mesh%section(isection)%nlayers
          associate(layer => mesh%section(isection)%layer(ilayer),       &
                    ngll  => mesh%section(isection)%layer(ilayer)%ngll,  &
                    nspec => mesh%section(isection)%layer(ilayer)%nspec, &
                    w     => mesh%section(isection)%layer(ilayer)%w,     &
                    hp    => mesh%section(isection)%layer(ilayer)%hp,    &
                    jac   => mesh%section(isection)%layer(ilayer)%jac,   &
                    r     => mesh%section(isection)%layer(ilayer)%r,     &
                    g     => mesh%section(isection)%layer(ilayer)%g,     &
                    ep    => mesh%section(isection)%layer(ilayer)%ep,    &
                    ibool => ibool%section(isection)%layer(ilayer),      &
                    rho => mesh%section(isection)%layer(ilayer)%rho)
            do ispec = 1,nspec
               if(allocated(drho)) then
                  if(ngll /= size(drho)) then
                     deallocate(drho)
                     allocate(drho(ngll))
                  end if
               else
                  allocate(drho(ngll))
               end if
               select type(layer)
               class is(spherical_solid_elastic_layer_mesh)                  
                  do inode = 1,ngll
                     drho(inode) = 0.0_dp
                     do jnode = 1,ngll
                        drho(inode) = drho(inode) + rho(jnode,ispec)*hp(inode,jnode)/jac(ispec)
                     end do
                  end do
               class is(spherical_fluid_elastic_layer_mesh)
                  drho(:) = layer%drho(:,ispec)                  
               end select
               do inode = 1,ngll
                  i = ibool%get(1,inode,ispec)
                  j = i
                  k = kd+1 + i-j
                  tmp = fourpi*bigg*drho(inode)*r(inode,ispec)**2
                  if(tmp /= 0.0_dp) tmp = tmp/g(inode,ispec)
                  tmp = tmp + l*(l+1)
                  a(k,j) = a(k,j) + tmp*w(inode)*jac(ispec) 
                  do jnode = inode,ngll
                     j = ibool%get(1,jnode,ispec)
                     k = kd+1+i-j
                     do knode = 1,ngll
                        a(k,j) = a(k,j) + hp(knode,inode)*hp(knode,jnode) &
                                        * r(knode,ispec)**2*w(knode)/jac(ispec)                      
                     end do
                  end do
                  tmp = rotfac*r(inode,ispec)**2
                  tmp = -fourpi*bigg*drho(inode)*tmp
                  if(tmp /= 0.0_dp) tmp = tmp/g(inode,ispec)
                  b(i,1) = b(i,1) + tmp*r(inode,ispec)**2 &
                                  * w(inode)*jac(ispec)
               end do
            end do
          end associate
       end do

       if(isection < mesh%nsections) then
          nspec = mesh%section(isection+1)%layer(1)%nspec
          ngll  = mesh%section(isection+1)%layer(1)%ngll
          tmp   = mesh%section(isection+1)%layer(1)%rho(1,1)
       else
          tmp = 0.0_dp
       end if
       ilayer = mesh%section(isection)%nlayers
       nspec  = mesh%section(isection)%layer(ilayer)%nspec
       ngll  = mesh%section(isection)%layer(ilayer)%ngll
       tmp = tmp - mesh%section(isection)%layer(ilayer)%rho(ngll,nspec)
       rtmp = mesh%section(isection)%layer(ilayer)%r(ngll,nspec)
       gtmp = mesh%section(isection)%layer(ilayer)%g(ngll,nspec)
       tmp = (fourpi*bigg*tmp*rtmp**2)/gtmp
       i = ibool%section(isection)%layer(ilayer)%get(1,ngll,nspec)
       j = i
       k = kd+1+i-j
       a(k,j) = a(k,j) + tmp
       tmp = tmp*rotfac*rtmp*rtmp
       b(i,1) = b(i,1) - tmp
       
    end do

    
    ! add in DTN term at the surface
    k = kd+1
    a(k,ndim) = a(k,ndim) + (l+1)*mesh%r2

    ! compute the factorisation
    call dpbtrf('U',ndim,kd,a,ldab,info)
    call check(info == 0,'calculate_gravity_and_ellipticity','problem with factorisation')

    ! solve the linear system
    call dpbtrs	('U',ndim,kd,1,a,ldab,b,ndim,info)
    call check(info == 0,'calculate_gravity_and_ellipticity','problem with substitution')

    ! extract the ellipticity
    do isection = 1,mesh%nsections
       do ilayer = 1,mesh%section(isection)%nlayers
          associate(ngll  => mesh%section(isection)%layer(ilayer)%ngll,  &
                    nspec => mesh%section(isection)%layer(ilayer)%nspec, &
                    hp    => mesh%section(isection)%layer(ilayer)%hp,    &
                    jac   => mesh%section(isection)%layer(ilayer)%jac,   &
                    r     => mesh%section(isection)%layer(ilayer)%r,     &
                    rho   => mesh%section(isection)%layer(ilayer)%rho,   &
                    g     => mesh%section(isection)%layer(ilayer)%g,     &
                    ep     => mesh%section(isection)%layer(ilayer)%ep,   &
                    ibool => ibool%section(isection)%layer(ilayer))
            do ispec = 1,nspec
               do inode = 1,ngll
                  i = ibool%get(1,inode,ispec)
                  tmp = rotfac*r(inode,ispec)**2
                  tmp = tmp + b(i,1)
                  if(r(inode,ispec) /= 0.0_dp) then
                     ep(inode,ispec) = tmp/(r(inode,ispec)*g(inode,ispec))
                  else
                     ep(inode,ispec) = 0.0_dp
                  end if
               end do
            end do
          end associate
       end do
    end do

    
    return
  end subroutine calculate_gravity_and_ellipticity



  
end module module_mesh
