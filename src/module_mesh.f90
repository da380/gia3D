module module_meshing

  use module_constants
  use module_physical_constants
  use module_spherical_model
  use module_quadrature
  use module_special_functions
  implicit none


  type spherical_model_mesh
     integer(i4b) :: nsections
     real(dp) :: r1
     real(dp) :: r2
     type(spherical_section_mesh), dimension(:), allocatable :: section     
  end type spherical_model_mesh
  
  type spherical_section_mesh
     integer(i4b) :: nlayers
     real(dp) :: r1
     real(dp) :: r2
     class(spherical_layer_mesh), dimension(:), allocatable :: layer     
  end type spherical_section_mesh
  
  type spherical_layer_mesh
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
  end type spherical_layer_mesh


  type, extends(spherical_layer_mesh) ::  spherical_solid_elastic_layer_mesh     
     real(dp), dimension(:,:), allocatable :: A
     real(dp), dimension(:,:), allocatable :: C
     real(dp), dimension(:,:), allocatable :: F
     real(dp), dimension(:,:), allocatable :: L
     real(dp), dimension(:,:), allocatable :: N
     real(dp), dimension(:,:), allocatable :: kappa
     real(dp), dimension(:,:), allocatable :: mu
  end type spherical_solid_elastic_layer_mesh


  type, extends(spherical_layer_mesh) ::  spherical_fluid_elastic_layer_mesh
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


  type boolean_array
     integer(i4b) :: ndim     
     integer(i4b) :: ngll
     integer(i4b) :: isection1
     type(boolean_array_section), dimension(:), allocatable :: section
   contains
     procedure :: get => get_boolean_array
  end type boolean_array

  type boolean_array_section
     integer(i4b) :: ilayer1
     type(boolean_array_layer), dimension(:), allocatable :: layer
   contains
     procedure :: get => get_boolean_array_section
  end type boolean_array_section
  
  type boolean_array_layer
     integer(i4b) :: nvar
     integer(i4b) :: ispec1
     integer(i4b), dimension(:,:,:), allocatable :: data
   contains
     procedure :: get => get_boolean_array_layer
  end type boolean_array_layer

  interface spherical_mesh
     procedure :: make_spherical_mesh
  end interface spherical_mesh


  
contains
  

  !===========================================================!
  !                routines for mesh building                 !
  !===========================================================!
  
  function make_spherical_mesh(ngll,model,drmax) result(mesh)
    integer(i4b), intent(in) :: ngll
    type(spherical_model), intent(in) :: model
    real(dp), intent(in) :: drmax
    type(spherical_model_mesh) :: mesh

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
               if(allocated(mesh_tmp)) deallocate(mesh_tmp)
               mesh_tmp = make_spherical_layer_mesh(ngll,layer,drmax)
               select type(mesh_tmp)
               class is(spherical_layer_mesh)
                  layer_mesh = mesh_tmp                
               end select
               
            class is(spherical_solid_elastic_layer)
               if(allocated(mesh_tmp)) deallocate(mesh_tmp)            
               mesh_tmp = make_spherical_solid_elastic_layer_mesh(ngll,layer,drmax)
               select type(mesh_tmp)
               class is(spherical_solid_elastic_layer_mesh)
                  layer_mesh = mesh_tmp                
               end select
               
               
            class is(spherical_fluid_elastic_layer)
               if(allocated(mesh_tmp)) deallocate(mesh_tmp)
               mesh_tmp = make_spherical_fluid_elastic_layer_mesh(ngll,layer,drmax)
               select type(mesh_tmp)
               class is(spherical_fluid_elastic_layer_mesh)
                  layer_mesh = mesh_tmp
               end select


            class is(spherical_solid_anelastic_layer)
               if(allocated(mesh_tmp)) deallocate(mesh_tmp)            
               mesh_tmp = make_spherical_solid_anelastic_layer_mesh(ngll,layer,drmax)
               select type(mesh_tmp)
               class is(spherical_solid_anelastic_layer_mesh)
                  layer_mesh = mesh_tmp                
               end select
               
               
            class is(spherical_fluid_anelastic_layer)
               if(allocated(mesh_tmp)) deallocate(mesh_tmp)
               mesh_tmp = make_spherical_fluid_anelastic_layer_mesh(ngll,layer,drmax)
               select type(mesh_tmp)
               class is(spherical_fluid_anelastic_layer_mesh)
                  layer_mesh = mesh_tmp
               end select
               
               
            class is(spherical_maxwell_layer)
               if(allocated(mesh_tmp)) deallocate(mesh_tmp)
               mesh_tmp = make_spherical_maxwell_layer_mesh(ngll,layer,drmax)
               select type(mesh_tmp)
               class is(spherical_maxwell_layer_mesh)
                  layer_mesh = mesh_tmp                
               end select
               
            end select
            
            
          end associate
          
       end do
       
    end do

    call calculate_gravity(mesh)

    return
  end function make_spherical_mesh

  
  function make_spherical_layer_mesh(ngll,layer,drmax) result(mesh)
    integer(i4b), intent(in) :: ngll
    class(spherical_layer), intent(in) :: layer
    real(dp), intent(in) :: drmax
    type(spherical_layer_mesh) :: mesh
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
  end function make_spherical_layer_mesh




  function make_spherical_solid_elastic_layer_mesh(ngll,layer,drmax) result(mesh)
    integer(i4b), intent(in) :: ngll
    class(spherical_solid_elastic_layer), intent(in) :: layer
    real(dp), intent(in) :: drmax
    type(spherical_solid_elastic_layer_mesh) :: mesh

    integer(i4b) :: ispec,nspec,inode,jspec
    real(dp) :: r
        
    mesh%spherical_layer_mesh = make_spherical_layer_mesh(ngll,layer,drmax)
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
  end function make_spherical_solid_elastic_layer_mesh



  function make_spherical_fluid_elastic_layer_mesh(ngll,layer,drmax) result(mesh)
    integer(i4b), intent(in) :: ngll
    class(spherical_fluid_elastic_layer), intent(in) :: layer
    real(dp), intent(in) :: drmax
    type(spherical_fluid_elastic_layer_mesh) :: mesh

    integer(i4b) :: ispec,nspec,inode,jnode
    real(dp) :: r
        
    mesh%spherical_layer_mesh = make_spherical_layer_mesh(ngll,layer,drmax)
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
  end function make_spherical_fluid_elastic_layer_mesh



  function make_spherical_solid_anelastic_layer_mesh(ngll,layer,drmax) result(mesh)
    integer(i4b), intent(in) :: ngll
    class(spherical_solid_anelastic_layer), intent(in) :: layer
    real(dp), intent(in) :: drmax
    type(spherical_solid_anelastic_layer_mesh) :: mesh

    integer(i4b) :: ispec,nspec,inode,jspec
    real(dp) :: r
        
    mesh%spherical_solid_elastic_layer_mesh = make_spherical_solid_elastic_layer_mesh(ngll,layer,drmax)
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
  end function make_spherical_solid_anelastic_layer_mesh


  
  function make_spherical_fluid_anelastic_layer_mesh(ngll,layer,drmax) result(mesh)
    integer(i4b), intent(in) :: ngll
    class(spherical_fluid_anelastic_layer), intent(in) :: layer
    real(dp), intent(in) :: drmax
    type(spherical_fluid_anelastic_layer_mesh) :: mesh

    integer(i4b) :: ispec,nspec,inode,jspec
    real(dp) :: r
        
    mesh%spherical_fluid_elastic_layer_mesh = make_spherical_fluid_elastic_layer_mesh(ngll,layer,drmax)
    allocate(mesh%qkappa(mesh%ngll,mesh%nspec))
    do ispec = 1,mesh%nspec
       do inode = 1,mesh%ngll
          r = mesh%r(inode,ispec)
          mesh%qkappa(inode,ispec) = layer%qkappa(r)
       end do
    end do
    
    
    return
  end function make_spherical_fluid_anelastic_layer_mesh

  

  function make_spherical_maxwell_layer_mesh(ngll,layer,drmax) result(mesh)
    integer(i4b), intent(in) :: ngll
    class(spherical_maxwell_layer), intent(in) :: layer
    real(dp), intent(in) :: drmax
    type(spherical_maxwell_layer_mesh) :: mesh

    integer(i4b) :: ispec,nspec,inode,jspec
    real(dp) :: r
        
    mesh%spherical_solid_elastic_layer_mesh = make_spherical_solid_elastic_layer_mesh(ngll,layer,drmax)
    allocate(mesh%eta(mesh%ngll,mesh%nspec))
    do ispec = 1,mesh%nspec
       do inode = 1,mesh%ngll
          r = mesh%r(inode,ispec)
          mesh%eta(inode,ispec) = layer%eta(r)          
       end do
    end do
    
    
    return
  end function make_spherical_maxwell_layer_mesh



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


  type(boolean_array) function build_boolean_toroidal(mesh,rstart) result(ibool)
    type(spherical_model_mesh), intent(in) :: mesh
    real(dp), intent(in), optional :: rstart

    integer(i4b) :: isection1,isection,nsections,ilayer,ilayer1,nlayers, &
                    ngll,nspec,ispec,ispec1,inode,count,ngllm,isection11
    real(dp) :: r1,r2,rs


    isection11 = mesh%nsections+1
    ! work out which section to start looking in
    do isection = mesh%nsections,1,-1

       associate(layer => mesh%section(isection)%layer(1))

         select type(layer)

         class is(spherical_solid_elastic_layer_mesh)

            isection11 = isection
            
         class is(spherical_fluid_elastic_layer_mesh)

            exit 
            
         class default

            stop 'build_boolean_toroidal: invalid mesh'

         end select
         
       end associate
       
    end do

    
    if(isection11 == mesh%nsections+1) then
       ibool%isection1 = isection11
       return
    end if
    
    r1 = mesh%section(isection11)%r1
    r2 = mesh%r2
    rs = r1
    if(present(rstart)) then
       rs = max(r1,rstart)
       rs = min(rs,r2)
    end if


    ! work out what section to start in
    nsections = mesh%nsections
    isection1 = nsections
    do isection = isection11,nsections
       associate(section => mesh%section(isection))
         r1 = section%r1
         r2 = section%r2
         if(rs >= r1 .and. rs < r2) then
            isection1 = isection
            exit
         end if
       end associate
    end do

    
    ibool%isection1 = isection1
    allocate(ibool%section(isection1:nsections))

    ! work out what layer to start in
    nlayers = mesh%section(isection1)%nlayers
    ilayer1 = nlayers
    do ilayer = 1,nlayers
       associate(layer => mesh%section(isection1)%layer(ilayer))
         r1 = layer%r1
         r2 = layer%r2
         if(rs >= r1 .and. rs < r2) then
            ilayer1 = ilayer
            exit
         end if
       end associate
    end do
    
    ibool%section(isection1)%ilayer1 = ilayer1
    allocate(ibool%section(isection1)%layer(ilayer1:nlayers))                
    do isection = isection1+1,nsections
       ibool%section(isection)%ilayer1 = 1
       nlayers = mesh%section(isection)%nlayers
       allocate(ibool%section(isection)%layer(nlayers))                
    end do


    ! work out what element to start in
    associate(layer => mesh%section(isection1)%layer(ilayer1))
      ngll   = layer%ngll
      nspec  = layer%nspec
      ispec1 = nspec
      do ispec = 1,nspec
         r1 = layer%r(1,ispec)
         r2 = layer%r(ngll,ispec)
         if(rs >= r1 .and. rs < r2) then
            ispec1 = ispec
            exit
         end if
      end do
    end associate
    ibool%section(isection1)%layer(ilayer1)%ispec1 = ispec1
    nlayers = mesh%section(isection1)%nlayers
    do ilayer = ilayer+1,nlayers
       ibool%section(isection1)%layer(ilayer)%ispec1 = 1       
    end do
    do isection = isection1 + 1,nsections
       ilayer1 = ibool%section(isection)%ilayer1
       nlayers = mesh%section(isection)%nlayers
       do ilayer = ilayer1,nlayers
          ibool%section(isection)%layer(ilayer)%ispec1 = 1       
       end do       
    end do

    
    ! build up the boolean array
    count = 0
    ngllm = 0
    do isection = isection1,nsections

       ilayer1 = ibool%section(isection)%ilayer1
       nlayers = mesh%section(isection)%nlayers
       
       do ilayer = ilayer1,nlayers

          ispec1 = ibool%section(isection)%layer(ilayer)%ispec1
          nspec  = mesh%section(isection)%layer(ilayer)%nspec
          ngll   = mesh%section(isection)%layer(ilayer)%ngll
          if(ngll > ngllm) ngllm = ngll
          
          
          associate(layer => mesh%section(isection)%layer(ilayer), &
                    ibool => ibool%section(isection)%layer(ilayer))

                        
            allocate(ibool%data(1,ngll,ispec1:nspec))
            do ispec = ispec1,nspec
               do inode = 1,ngll
                  count = count + 1
                  ibool%data(1,inode,ispec) = count
               end do
               count = count-1
            end do
                        
          end associate
          
       end do
       
    end do

    ibool%ndim = count+1
    ibool%ngll = ngllm

  
    return
  end function build_boolean_toroidal



  type(boolean_array) function build_boolean_spheroidal(mesh,rstart) result(ibool)
    type(spherical_model_mesh), intent(in) :: mesh
    real(dp), intent(in), optional :: rstart

    logical :: solid
    integer(i4b) :: isection1,isection,nsections,ilayer,ilayer1,nlayers, &
                    ngll,nspec,ispec,ispec1,inode,count,ngllm
    real(dp) :: r1,r2,rs
    
    r1 = mesh%r1
    r2 = mesh%r2
    rs = r1
    if(present(rstart)) then
       rs = max(r1,rstart)
       rs = min(rs,r2)
    end if

    ! work out what section to start in
    nsections = mesh%nsections
    isection1 = nsections
    do isection = 1,nsections
       associate(section => mesh%section(isection))
         r1 = section%r1
         r2 = section%r2
         if(rs >= r1 .and. rs < r2) then
            isection1 = isection
            exit
         end if
       end associate
    end do

    ibool%isection1 = isection1
    allocate(ibool%section(isection1:nsections))


    ! work out what layer to start in
    nlayers = mesh%section(isection1)%nlayers
    ilayer1 = nlayers
    do ilayer = 1,nlayers
       associate(layer => mesh%section(isection1)%layer(ilayer))
         r1 = layer%r1
         r2 = layer%r2
         if(rs >= r1 .and. rs < r2) then
            ilayer1 = ilayer
            exit
         end if
       end associate
    end do
    
    ibool%section(isection1)%ilayer1 = ilayer1
    allocate(ibool%section(isection1)%layer(ilayer1:nlayers))                
    do isection = isection1+1,nsections
       ibool%section(isection)%ilayer1 = 1
       nlayers = mesh%section(isection)%nlayers
       allocate(ibool%section(isection)%layer(nlayers))                
    end do


    ! work out what element to start in
    associate(layer => mesh%section(isection1)%layer(ilayer1))
      ngll   = layer%ngll
      nspec  = layer%nspec
      ispec1 = nspec
      do ispec = 1,nspec
         r1 = layer%r(1,ispec)
         r2 = layer%r(ngll,ispec)
         if(rs >= r1 .and. rs < r2) then
            ispec1 = ispec
            exit
         end if
      end do
    end associate
    ibool%section(isection1)%layer(ilayer1)%ispec1 = ispec1
    nlayers = mesh%section(isection1)%nlayers
    do ilayer = ilayer+1,nlayers
       ibool%section(isection1)%layer(ilayer)%ispec1 = 1       
    end do
    do isection = isection1 + 1,nsections
       ilayer1 = ibool%section(isection)%ilayer1
       nlayers = mesh%section(isection)%nlayers
       do ilayer = ilayer1,nlayers
          ibool%section(isection)%layer(ilayer)%ispec1 = 1       
       end do       
    end do

    
    ! build up the boolean array
    solid = .false.
    count = 0
    ngllm = 0
    do isection = isection1,nsections

       ilayer1 = ibool%section(isection)%ilayer1
       nlayers = mesh%section(isection)%nlayers
       
       do ilayer = ilayer1,nlayers

          ispec1 = ibool%section(isection)%layer(ilayer)%ispec1
          nspec  = mesh%section(isection)%layer(ilayer)%nspec
          ngll   = mesh%section(isection)%layer(ilayer)%ngll
          if(ngll > ngllm) ngllm = ngll
          
          
          associate(layer => mesh%section(isection)%layer(ilayer), &
                    ibool => ibool%section(isection)%layer(ilayer))


            select type(layer)
               
            class is(spherical_solid_elastic_layer_mesh)
               solid = .true.
               allocate(ibool%data(3,ngll,ispec1:nspec))
               do ispec = ispec1,nspec
                  do inode = 1,ngll
                     count = count + 1
                     ibool%data(1,inode,ispec) = count
                     count = count + 1
                     ibool%data(2,inode,ispec) = count
                     count = count + 1
                    ibool%data(3,inode,ispec) = count
                  end do
                  count = count-1
               end do
               
               
            class is(spherical_fluid_elastic_layer_mesh)
               allocate(ibool%data(1,ngll,ispec1:nspec))
               do ispec = ispec1,nspec                  
                  do inode = 1,ngll
                     if(solid) then
                        count = count-1
                        ibool%data(1,inode,ispec) = count
                        count = count + 2
                        solid = .false.
                     else
                        count = count+1
                        ibool%data(1,inode,ispec) = count
                     end if
                  end do
                  count = count-1
               end do

            class default
               stop 'build_boolean_spheroidal: invalid mesh'
            end select
            
          end associate
          
       end do
       
    end do

    ibool%ndim = count+1
    ibool%ngll = ngllm
    
    return
  end function build_boolean_spheroidal
  

  !================================================================================!
  !                            routine to calculate gravity                        !
  !================================================================================!


  subroutine calculate_gravity(mesh)
    class(spherical_model_mesh), intent(inout) :: mesh

    integer(i4b) :: ndim,kd,ldab,isection,ilayer,inode, &
                    ispec,jnode,knode,ngll,nspec,i,j,k,info
    type(boolean_array) :: ibool
    real(dp), dimension(:,:), allocatable :: a,b

    
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
    call error(info /= 0,'calculate_gravity','problem with factorisation')

    ! solve the linear system
    call dpbtrs	('U',ndim,kd,1,a,ldab,b,ndim,info)
    call error(info /= 0,'calculate_gravity','problem with substitution')

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

    return
  end subroutine calculate_gravity


  
end module module_meshing
