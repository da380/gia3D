module module_matrix

  use module_constants
  use module_error
  use module_mesh
  use module_physical_constants
  implicit none

  ! small ratio for determining starting depths
  real(dp), parameter :: eps_def  = 1.0e-6_dp
  
  type radial_matrix
     logical :: factorised  = .false.
     integer(i4b) :: l
     integer(i4b) :: ndim
     integer(i4b) :: kd
     integer(i4b) :: ldab
     type(boolean_array) :: ibool
     real(dp), dimension(:,:), allocatable :: a
  end type radial_matrix

  type radial_matrices
     logical :: factorised = .false.
     logical :: spheroidals = .false.
     logical :: toroidals   = .false.
     integer(i4b) :: lmax
     integer(i4b) :: ndim_sph
     integer(i4b) :: ndim_tor
     type(radial_matrix), dimension(:), allocatable :: sph
     type(radial_matrix), dimension(:), allocatable :: tor
  end type radial_matrices
  
contains


  type(radial_matrices) function build_radial_matrices(mesh,lmax,spheroidals,toroidals,factor) result(mat)
    type(spherical_model_mesh), intent(in) :: mesh
    integer(i4b), intent(in) :: lmax
    logical, intent(in), optional :: spheroidals
    logical, intent(in), optional :: toroidals
    logical, intent(in), optional :: factor

    integer(i4b) :: l

    if(present(factor)) then
       mat%factorised = factor
    else
       mat%factorised = .true.
    end if
    
    if(lmax > 0) then
       mat%lmax = lmax
    else
       stop 'build_radial_matrices: lmax < 1'
    end if
    if(present(spheroidals)) then
       mat%spheroidals = spheroidals
       allocate(mat%sph(1:lmax))
    end if
    if(present(toroidals)) then
       mat%toroidals = toroidals
       allocate(mat%tor(1:lmax))
    end if

    mat%ndim_sph = 0
    mat%ndim_tor = 0
    do l = 1,lmax

       if(mat%spheroidals) then
          mat%sph(l) =  build_spheroidal_matrix(mesh,l,factor)
          if(mat%sph(l)%ndim > mat%ndim_sph) mat%ndim_sph = mat%sph(l)%ndim
       end if

       if(mat%toroidals) then
          mat%tor(l) =  build_toroidal_matrix(mesh,l,factor)
          if(mat%tor(l)%ndim > mat%ndim_tor) mat%ndim_tor = mat%tor(l)%ndim
       end if
       
    end do
    

    return
  end function build_radial_matrices
  
  !=================================================================!
  !                        Toroidal routines                        !
  !=================================================================!


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
    do ilayer = ilayer1+1,nlayers
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



  type(radial_matrix) function build_toroidal_matrix(mesh,l,factor) result(mat)
    type(spherical_model_mesh), intent(in) :: mesh
    integer(i4b) :: l
    logical, intent(in), optional :: factor

    logical :: factor_local
    integer(i4b) :: ndim,kd,ldab,isection,ilayer,info

    if(present(factor)) then
       factor_local = factor
    else
       factor_local = .true.
    end if       
    call check(l > 0,'build_toroidal_matrix','l<1')    
    mat%l = l
    mat%ibool = build_boolean_toroidal(mesh,toroidal_start(mesh,l))
    if(l == 1) mat%ibool%ndim = mat%ibool%ndim-1
    ndim = mat%ibool%ndim
    kd = mat%ibool%ngll-1
    ldab = kd+1    
    mat%ndim = ndim
    mat%kd = kd
    mat%ldab = ldab
    if(allocated(mat%a)) deallocate(mat%a); allocate(mat%a(ldab,ndim))
    mat%a = 0.0_dp
    do isection = mat%ibool%isection1,mesh%nsections
       do ilayer = mat%ibool%section(isection)%ilayer1,mesh%section(isection)%nlayers
          associate(layer => mesh%section(isection)%layer(ilayer),      &
                    ibool => mat%ibool%section(isection)%layer(ilayer), &
                    a => mat%a)
            select type(layer)               
            class is(spherical_solid_elastic_layer_mesh)               
               call build_toroidal_matrix_layer(layer,ibool,l,a)            
            class default
               stop 'build_toroidal_matrix: invalid mesh'
            end select
          end associate
       end do
    end do
    
    if(factor_local) then
       call dpbtrf('U',ndim,kd,mat%a,ldab,info)
       call check(info == 0,'build_toroidal_matrix','problem with factorisation')
       mat%factorised = .true.
    end if
    return
  end function build_toroidal_matrix


  subroutine build_toroidal_matrix_layer(layer,ibool,l,a)
    class(spherical_solid_elastic_layer_mesh), intent(in) :: layer
    type(boolean_array_layer), intent(in) :: ibool
    integer(i4b), intent(in) :: l
    real(dp), dimension(:,:), intent(inout) :: a

    logical :: itop,jtop,top
    integer(i4b) :: ndim,kd,ldab,isection,ilayer,inode,  &
                    ispec,jnode,knode,ngll,nspec,i,j,k,ispec1
    real(dp) :: tmp1,tmp2,zeta2,zeta2m2,rt,mut
    
    ! set some parameters
    zeta2 = l*(l+1)
    zeta2m2 = zeta2-2.0_dp
    ngll = layer%ngll
    nspec = layer%nspec
    ispec1 = ibool%ispec1
    kd = ngll-1

    ! assemble the matrix
    associate(  r => layer%r,   &
               mu => layer%mu,  &
              rho => layer%rho, &         
                w => layer%w,   &
               hp => layer%hp,  &
              jac => layer%jac)    
      do ispec = ispec1,nspec
         top = (l == 1) .and. layer%top .and. (ispec == nspec)
         do inode = 1,ngll
            itop = top .and. (inode == ngll)
            if(itop) cycle
            rt   =  r(inode,ispec)
            mut  = mu(inode,ispec)
            i = ibool%get(1,inode,ispec)
            k = kd+1
            a(k,i) = a(k,i) + zeta2*zeta2m2*mut*w(inode)*jac(ispec)
            do jnode = inode,ngll
               jtop = top .and. (jnode == ngll)
               if(jtop) cycle
               j = ibool%get(1,jnode,ispec)
               k = kd+1+i-j
               do knode = 1,ngll
                  rt  =  r(knode,ispec)
                  mut = mu(knode,ispec)
                  tmp1 = rt*hp(knode,inode)/jac(ispec)
                  if(inode == knode) tmp1 = tmp1-1.0_dp
                  tmp2 = rt*hp(knode,jnode)/jac(ispec)
                  if(jnode == knode) tmp2 = tmp2-1.0_dp
                  a(k,j) = a(k,j) + zeta2*mut*tmp1*tmp2*w(knode)*jac(ispec)
               end do
            end do
         end do
      end do
    end associate
    
    return
  end subroutine build_toroidal_matrix_layer

  
  real(dp) function toroidal_start(mesh,l,eps_in) result(rs)
    type(spherical_model_mesh), intent(in) :: mesh
    integer(i4b), intent(in) :: l
    real(dp), intent(in), optional :: eps_in
    real(dp) :: eps
    if(present(eps_in)) then
       eps = eps_in
    else
       eps = eps_def
    end if
    rs = log(eps)/(l+1)
    rs = mesh%r2*exp(rs)
    return
  end function toroidal_start

  !=================================================================!
  !                       Spheroidal routines                       !
  !=================================================================!


  type(boolean_array) function build_boolean_spheroidal(mesh,rstart) result(ibool)
    type(spherical_model_mesh), intent(in) :: mesh
    real(dp), intent(in), optional :: rstart

    logical :: fluid
    integer(i4b) :: isection1,isection,nsections,ilayer,ilayer1,nlayers, &
                    ngll,nspec,ispec,ispec1,inode,count,ngllm,ndim
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
    do ilayer = ilayer1+1,nlayers
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
    fluid = .false.
    count = 0
    ngllm = 0
    ndim = 0
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
               allocate(ibool%data(3,ngll,ispec1:nspec))
               do ispec = ispec1,nspec
                  do inode = 1,ngll
                     if(fluid) then
                        count = count + 2
                        ibool%data(1,inode,ispec) = count
                        count = count + 1
                        ibool%data(2,inode,ispec) = count
                        count = count - 2
                        ibool%data(3,inode,ispec) = count
                        count = count + 2
                        fluid = .false.
                     else                     
                        count = count + 1
                        ibool%data(1,inode,ispec) = count
                        count = count + 1
                        ibool%data(2,inode,ispec) = count
                        count = count + 1
                        ibool%data(3,inode,ispec) = count
                     end if
                  end do
                  count = count-3
               end do
               
               
            class is(spherical_fluid_elastic_layer_mesh)
               allocate(ibool%data(1,ngll,ispec1:nspec))
               if(.not.fluid .and. count /= 0) count = count + 2               
               do ispec = ispec1,nspec                  
                  do inode = 1,ngll
                     count = count+1
                     ibool%data(1,inode,ispec) = count
                  end do
                  count = count-1
               end do
               fluid = .true.
            class default
               stop 'build_boolean_spheroidal: invalid mesh'
            end select
            
          end associate
          
       end do

    end do

    ! get the dimension of the system
    associate(layer => mesh%section(mesh%nsections)%layer(1))
      select type(layer)
      class is(spherical_solid_elastic_layer_mesh)
         ibool%ndim = count+3
      class is(spherical_fluid_elastic_layer_mesh)
         ibool%ndim = count+1
      end select
    end associate
    ! store the largest spectral element order
    ibool%ngll = ngllm
    
    return
  end function build_boolean_spheroidal


  type(radial_matrix) function build_spheroidal_matrix(mesh,l,factor) result(mat)
    type(spherical_model_mesh), intent(in) :: mesh
    integer(i4b) :: l
    logical, intent(in), optional :: factor

    logical :: top,itop,jtop,factor_local,centre
    integer(i4b) :: ndim,kd,ldab,isection,ilayer,inode,  &
                    ispec,jnode,knode,ngll,nspec,i,j,k,info,ispec1, &
                    i1,i2,i3,j1,j2,j3,nlayers,ilayer1,isection1,nsections
    real(dp) :: tmp1,tmp2,zeta2,zeta2m2,zetac,ifpibigg,rt,mut,kapt,gt, &
                rhot,drhot,fac

    if(present(factor)) then
       factor_local = factor
    else
       factor_local = .true.
    end if       
    call check(l > 0,'build_spheroidal_matrix','l<1')    
    mat%l = l
    mat%ibool = build_boolean_spheroidal(mesh,spheroidal_start(mesh,l))
    if(l == 1) mat%ibool%ndim = mat%ibool%ndim-1
    ndim = mat%ibool%ndim
    ngll = mat%ibool%ngll
    kd = 3*(ngll-1)+2
    ldab = kd+1
    mat%ndim = ndim
    mat%kd = kd
    mat%ldab = ldab
    if(allocated(mat%a)) deallocate(mat%a); allocate(mat%a(ldab,ndim))
    mat%a = 0.0_dp

    
    isection1 = mat%ibool%isection1
    nsections = mesh%nsections
    do isection = isection1,nsections
       ilayer1 = mat%ibool%section(isection)%ilayer1
       nlayers = mesh%section(isection)%nlayers
       do ilayer = ilayer1,nlayers
          associate(layer => mesh%section(isection)%layer(ilayer),      &
                    ibool => mat%ibool%section(isection)%layer(ilayer), &
                    a => mat%a)
            select type(layer)               
            class is(spherical_solid_elastic_layer_mesh)
               call build_spheroidal_matrix_solid_layer(layer,ibool,l,a)
            class is(spherical_fluid_elastic_layer_mesh)
               call build_spheroidal_matrix_fluid_layer(layer,ibool,l,a)
            class default
               stop 'build_spheroidal_mesh: invalid mesh'
            end select
          end associate
       end do
    end do
    
    if(factor_local) then
       call dpbtrf('U',ndim,kd,mat%a,ldab,info)
       call check(info == 0,'build_spheroidal_matrix','problem with factorisation')
       mat%factorised = .true.
    end if
    
    return
  end function build_spheroidal_matrix
  

  subroutine build_spheroidal_matrix_solid_layer(layer,ibool,l,a)
    class(spherical_solid_elastic_layer_mesh), intent(in) :: layer
    type(boolean_array_layer), intent(in) :: ibool
    integer(i4b), intent(in) :: l
    real(dp), dimension(:,:), intent(inout) :: a

    logical :: top,itop,jtop
    integer(i4b) :: ngll,nspec,ispec1,kd,ispec,i,j,k, &
                    inode,jnode,knode
    real(dp) :: zeta2,zeta2m2,zetac,ifpibigg,rt,rhot, &
                kapt,mut,tmp1,tmp2,gt


    
    ! set some parameters
    zeta2 = l*(l+1)
    zeta2m2 = zeta2-2.0_dp
    zetac = (4.0_dp/3.0_dp)*zeta2-2.0_dp
    ifpibigg = 1.0_dp/(4.0_dp*pi*bigg)
    ngll = layer%ngll
    nspec = layer%nspec
    ispec1 = ibool%ispec1
    kd = 3*(ngll-1)+2

    
    ! do we need to deal with the DNT at the bottom of the mesh?
    if(layer%bottom .and. layer%r(1,ispec1) /= 0.0_dp) then
       rt = layer%r(1,ispec1)
       i = ibool%get(3,1,ispec1)
       k = kd+1
       a(k,i) = a(k,i) + l*ifpibigg*rt       
    end if


    
    
    ! assemble the matrix
    associate(r     => layer%r,     &
              rho   => layer%rho,   &
              mu    => layer%mu,    &
              kappa => layer%kappa, &                                        
              g     => layer%g,     &
              jac   => layer%jac,   &
              hp    => layer%hp,    &                         
              w     => layer%w)
                 
      do ispec = ispec1,nspec
         top = (l == 1) .and. layer%top .and. (ispec == nspec)
                    
         do inode = 1,ngll
            itop = top .and. (inode == ngll)

             rt  =     r(inode,ispec)
            rhot =   rho(inode,ispec)
            mut  =    mu(inode,ispec)
            kapt = kappa(inode,ispec)
            gt =       g(inode,ispec)
                       
            ! u-u'
            i = ibool%get(1,inode,ispec)
            j = i
            k = kd+1+i-j
            a(k,j) = a(k,j) + zeta2*mut*w(inode)*jac(ispec)
            a(k,j) = a(k,j) + 4.0_dp*rhot*(pi*bigg*rhot*rt-gt)  & 
                            * rt*w(inode)*jac(ispec)
            ! u-v'
            i = ibool%get(1,inode,ispec)
            j = ibool%get(2,inode,ispec)
            k = kd+1+i-j
            a(k,j) = a(k,j) + zeta2*rhot*gt*rt*w(inode)*jac(ispec)
                       
            ! v-v'
            i = ibool%get(2,inode,ispec)
            j = ibool%get(2,inode,ispec)
            k = kd+1+i-j
            a(k,j) = a(k,j) + zeta2*(zeta2*kapt+zetac*mut)*w(inode)*jac(ispec)

                       
            if(.not. itop) then

               ! v-phi'
               i = ibool%get(2,inode,ispec)
               j = ibool%get(3,inode,ispec)
               if(j >= i) then
                  k = kd+1+i-j
                  a(k,j) = a(k,j) + zeta2*rhot*rt*w(inode)*jac(ispec)
               end if

               ! phi-v'
               i = ibool%get(3,inode,ispec)
               j = ibool%get(2,inode,ispec)
               if(j >= i) then
                  k = kd+1+i-j
                  a(k,j) = a(k,j) + zeta2*rhot*rt*w(inode)*jac(ispec)
               end if
                          
               ! phi-phi'
               i = ibool%get(3,inode,ispec)
               j = i
               k = kd+1+i-j
               a(k,j) = a(k,j) + ifpibigg*zeta2*w(inode)*jac(ispec)
               
            end if

                       
                          
            do jnode = inode,ngll                          
               jtop = top .and. (jnode == ngll)


               ! u-u'
               i = ibool%get(1,inode,ispec)
               j = ibool%get(1,jnode,ispec)
               k = kd+1+i-j
               do knode = 1,ngll
                  rt = r(knode,ispec)
                  kapt = kappa(knode,ispec)
                  mut = mu(knode,ispec)
                  tmp1 = rt*hp(knode,inode)/jac(ispec)
                  if(inode == knode) tmp1 = tmp1+2.0_dp
                  tmp2 = rt*hp(knode,jnode)/jac(ispec)
                  if(jnode == knode) tmp2 = tmp2+2.0_dp
                  a(k,j) = a(k,j) + kapt*tmp1*tmp2*w(knode)*jac(ispec)
                  tmp1 = rt*hp(knode,inode)/jac(ispec)
                  if(inode == knode) tmp1 = tmp1-1.0_dp
                  tmp2 = rt*hp(knode,jnode)/jac(ispec)
                  if(jnode == knode) tmp2 = tmp2-1.0_dp
                  a(k,j) = a(k,j) + (4.0_dp/3.0_dp)*mut*tmp1*tmp2*w(knode)*jac(ispec)
               end do
               

               ! u-v'
               i = ibool%get(1,inode,ispec)
               j = ibool%get(2,jnode,ispec)
               k = kd+1+i-j
               kapt = kappa(jnode,ispec)
               rt = r(jnode,ispec)
               mut = mu(jnode,ispec)
               tmp1 = rt*hp(jnode,inode)/jac(ispec)
               tmp2 = tmp1
               if(inode == jnode) tmp1 = tmp1+2.0_dp
               a(k,j) = a(k,j) - zeta2*kapt*tmp1*w(jnode)*jac(ispec)                
               if(inode == jnode) tmp2 = tmp2-1.0_dp
               a(k,j) = a(k,j) + (2.0_dp/3.0_dp)*zeta2*mut*tmp2*w(jnode)*jac(ispec)
               mut = mu(inode,ispec)
               rt = r(inode,ispec)
               tmp1 = rt*hp(inode,jnode)/jac(ispec)
               if(inode == jnode) tmp1 = tmp1-1.0_dp
               a(k,j) = a(k,j) + zeta2*mut*tmp1*w(inode)*jac(ispec)                          


               ! u-phi'
               if(.not.jtop) then
                  i = ibool%get(1,inode,ispec)
                  j = ibool%get(3,jnode,ispec)
                  if(j >= i) then
                     k = kd+1+i-j
                     rt = r(inode,ispec)
                     rhot = rho(inode,ispec)
                     a(k,j) = a(k,j) + rhot*rt*rt*hp(inode,jnode)*w(inode)
                  end if
               end if
               
               ! v-u'
               i = ibool%get(2,inode,ispec)
               j = ibool%get(1,jnode,ispec)
               if(j >= i) then
                  k = kd+1+i-j
                  kapt = kappa(inode,ispec)
                  rt = r(inode,ispec) 
                  mut = mu(inode,ispec)
                  tmp1 = rt*hp(inode,jnode)/jac(ispec)
                  tmp2 = tmp1
                  if(inode == jnode) tmp1 = tmp1+2.0_dp
                  a(k,j) = a(k,j) - zeta2*kapt*tmp1*w(inode)*jac(ispec)
                  if(inode == jnode) tmp2 = tmp2-1.0_dp
                  a(k,j) = a(k,j) + (2.0_dp/3.0_dp)*zeta2*mut*tmp2*w(inode)*jac(ispec)
                  mut = mu(jnode,ispec)
                  rt = r(jnode,ispec)
                  tmp1 = rt*hp(jnode,inode)/jac(ispec)
                  if(inode == jnode) tmp1 = tmp1-1.0_dp
                  a(k,j) = a(k,j) + zeta2*mut*tmp1*w(jnode)*jac(ispec)
               end if
               
                          
               ! v-v'
               i = ibool%get(2,inode,ispec)
               j = ibool%get(2,jnode,ispec)
               k = kd+1+i-j
               do knode = 1,ngll
                  rt  = r(knode,ispec)
                  mut = mu(knode,ispec)                
                  tmp1 = rt*hp(knode,inode)/jac(ispec)
                  if(knode == inode) tmp1 = tmp1-1.0_dp
                  tmp2 = rt*hp(knode,jnode)/jac(ispec)
                  if(knode == jnode) tmp2 = tmp2-1.0_dp          
                  a(k,j) = a(k,j) + zeta2*mut*tmp1*tmp2*w(knode)*jac(ispec)
               end do
               
                          
               ! phi-u'
               if(.not.itop) then
                  i = ibool%get(3,inode,ispec)
                  j = ibool%get(1,jnode,ispec)
                  if(j >= i) then
                     k = kd+1+i-j
                     rt = r(jnode,ispec)
                     rhot = rho(jnode,ispec)
                     a(k,j) = a(k,j) + rhot*rt*rt*hp(jnode,inode)*w(jnode)
                  end if
               end if
               
               ! phi-phi'
               if(.not.itop .and. .not.jtop) then
                  i = ibool%get(3,inode,ispec)
                  j = ibool%get(3,jnode,ispec)
                  k = kd+1+i-j                          
                  do knode = 1,ngll
                     rt = r(knode,ispec)
                     a(k,j) = a(k,j) + ifpibigg*hp(knode,inode)*hp(knode,jnode)  & 
                                     * rt*rt*w(knode)/jac(ispec)                             
                  end do
               end if
               
            end do
            
         end do
      end do
      
    end associate


    ! deal with fluid-solid boundary
    if(layer%below_fluid .and. ispec1 == 1) then

       inode = 1
       ispec = 1
       rt = layer%r(inode,ispec)
       rhot = layer%below_fluid_density
       gt = layer%g(inode,ispec)
       ! u-u'
       i = ibool%get(1,inode,ispec)
       j = i
       k = kd+1+i-j
       a(k,j) = a(k,j) + rt*rt*gt*rhot
       ! phi-u'
       i = ibool%get(3,inode,ispec)
       j = ibool%get(1,inode,ispec)  
       k = kd+1+i-j                    
       a(k,j) = a(k,j) + rt*rt*rhot
       
    end if

    
    ! deal with solid-fluid boundary
    if(layer%above_fluid) then

       inode = ngll
       ispec = nspec
       rt = layer%r(inode,ispec)
       rhot = layer%above_fluid_density
       gt = layer%g(inode,ispec)
       ! u-u'
       i = ibool%get(1,inode,ispec)
       j = i
       k = kd+1+i-j
       a(k,j) = a(k,j) - rt*rt*gt*rhot
       ! u-phi'
       i = ibool%get(1,inode,ispec)
       j = ibool%get(3,inode,ispec)
       k = kd+1+i-j                    
       a(k,j) = a(k,j) - rt*rt*rhot
       
    end if

    ! DNT term at the surface
    if(layer%top .and. l  > 1) then
       inode = ngll
       ispec = nspec
       i = ibool%get(3,inode,ispec)
       j = i
       k = kd+1+i-j
       rt = layer%r(inode,ispec)
       a(k,j) = a(k,j)+(l+1)*ifpibigg*rt
    end if
    
    return
  end subroutine build_spheroidal_matrix_solid_layer



  subroutine build_spheroidal_matrix_fluid_layer(layer,ibool,l,a)
    class(spherical_fluid_elastic_layer_mesh), intent(in) :: layer
    type(boolean_array_layer), intent(in) :: ibool
    integer(i4b), intent(in) :: l
    real(dp), dimension(:,:), intent(inout) :: a


    logical :: top,itop,jtop
    integer(i4b) :: ngll,nspec,ispec1,kd,ispec,i,j,k, &
                    inode,jnode,knode
    real(dp) :: zeta2,zeta2m2,zetac,ifpibigg,rt, &
                gt,drhot,fac

    ! set some parameters
    zeta2 = l*(l+1)
    zeta2m2 = zeta2-2.0_dp
    zetac = (4.0_dp/3.0_dp)*zeta2-2.0_dp
    ifpibigg = 1.0_dp/(4.0_dp*pi*bigg)
    ngll = layer%ngll
    nspec = layer%nspec
    ispec1 = ibool%ispec1
    kd = 3*(ngll-1)+2

    ! do we need to deal with the DNT at the bottom of the mesh?
    if(layer%bottom .and. layer%r(1,ispec1) /= 0.0_dp) then
       rt = layer%r(1,ispec1)
       i = ibool%get(1,1,ispec1)
       k = kd+1
       a(k,i) = a(k,i) + l*ifpibigg*rt       
    end if


    
    associate(r    => layer%r,     &
              drho => layer%drho,  &
              g    => layer%g,     &
              jac  => layer%jac,   &
              hp   => layer%hp,    &                         
              w    => layer%w)
      
      do ispec = ispec1,nspec
         top = (l == 1) .and. layer%top .and. (ispec == nspec)
         
         do inode = 1,ngll
            itop = top .and. (inode == ngll)                       
            if(itop) cycle

            !phi-phi'
            i = ibool%get(1,inode,ispec)
            j = i
            k = kd+1+i-j
            a(k,j) = a(k,j) + ifpibigg*zeta2*w(inode)*jac(ispec)

            ! phi-phi'
            drhot = drho(inode,ispec)
            gt = g(inode,ispec)
            fac = rt*rt*drhot
            if(gt /= 0.0_dp) fac = fac/gt
            a(k,j) = a(k,j) + fac*w(inode)*jac(ispec)
            
            do jnode = inode,ngll                          
               jtop = top .and. (jnode == ngll)
               if(jtop) cycle

               !phi-phi'
               i = ibool%get(1,inode,ispec)                          
               j = ibool%get(1,jnode,ispec)                          
               k = kd+1+i-j                          
               do knode = 1,ngll
                  rt = r(knode,ispec)
                  a(k,j) = a(k,j) + ifpibigg*hp(knode,inode)*hp(knode,jnode)  & 
                                  * rt*rt*w(knode)/jac(ispec)                             
               end do
            end do
            
         end do
      end do
    end associate


    ! DNT term at the surface
    if(layer%top .and. l  > 1) then
       inode = ngll
       ispec = nspec
       i = ibool%get(1,inode,ispec)
       j = i
       k = kd+1+i-j
       rt = layer%r(inode,ispec)
       a(k,j) = a(k,j)+(l+1)*ifpibigg*rt
    end if
    
    return
  end subroutine build_spheroidal_matrix_fluid_layer
  
  real(dp) function spheroidal_start(mesh,l,eps_in) result(rs)
    type(spherical_model_mesh), intent(in) :: mesh
    integer(i4b), intent(in) :: l
    real(dp), intent(in), optional :: eps_in
    real(dp) :: eps
    if(present(eps_in)) then
       eps = eps_in
    else
       eps = eps_def
    end if
    rs = log(eps)/(l+1)
    rs = mesh%r2*exp(rs)
    return
  end function spheroidal_start


  
end module module_matrix

