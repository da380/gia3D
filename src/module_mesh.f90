module module_meshing

  use module_constants
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


contains


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
    call quad%set(ngll,legendre())
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

    integer(i4b) :: ispec,nspec,inode,jspec
    real(dp) :: r
        
    mesh%spherical_layer_mesh = make_spherical_layer_mesh(ngll,layer,drmax)
    allocate(mesh%kappa(mesh%ngll,mesh%nspec))
    do ispec = 1,mesh%nspec
       do inode = 1,mesh%ngll
          r = mesh%r(inode,ispec)
          mesh%kappa(inode,ispec) = layer%kappa(r)
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




  



  
end module module_meshing
