module module_SEM_matrix

  use module_constants
  use module_error
  use module_mesh
  use module_physical_constants, only : bigg
  implicit none

  
  type radial_matrix
     logical :: factorised  = .false.
     integer(i4b) :: l
     type(boolean_array) :: ibool
     real(dp), dimension(:,:), allocatable :: a
  end type radial_matrix


  
contains


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



  type(radial_matrix) function build_toroidal_matrix(mesh,l,factor) result(mat)
    type(spherical_model_mesh), intent(in) :: mesh
    integer(i4b) :: l
    logical, intent(in), optional :: factor

    logical :: spectop,ibnd,jbnd,factor_local
    integer(i4b) :: ndim,kd,ldab,isection,ilayer,inode,  &
                    ispec,jnode,knode,ngll,nspec,i,j,k,info
    real(dp) :: tmp1,tmp2,zeta2,zeta2m2,rt,mut

    if(present(factor)) then
       factor_local = factor
    else
       factor_local = .true.
    end if       
    call error(l < 1,'build_toroidal_matrix','l<1')    
    mat%l = l
    zeta2 = l*(l+1)
    zeta2m2 = zeta2-2.0_dp
    mat%ibool = build_boolean_toroidal(mesh,toroidal_start(mesh,l))
    if(l == 1) mat%ibool%ndim = mat%ibool%ndim-1
    ndim = mat%ibool%ndim
    kd = mat%ibool%ngll-1
    ldab = kd+1
    if(allocated(mat%a)) deallocate(mat%a); allocate(mat%a(ldab,ndim))
    mat%a = 0.0_dp
    do isection = mat%ibool%isection1,mesh%nsections
       spectop = (l == 1) .and. (isection == mesh%nsections)
       do ilayer = mat%ibool%section(isection)%ilayer1,mesh%section(isection)%nlayers
          spectop = spectop .and. (ilayer == mesh%section(isection)%nlayers)
          associate(layer => mesh%section(isection)%layer(ilayer),       &
                    ibool => mat%ibool%section(isection)%layer(ilayer), &
                    a => mat%a)
            select type(layer)               
            class is(spherical_solid_elastic_layer_mesh)
               associate(ngll  => layer%ngll,   &
                    nspec => layer%nspec,  &
                    w     => layer%w,      &
                    hp    => layer%hp,     &
                    jac   => layer%jac,    &
                    r     => layer%r,      &
                    mu    => layer%mu)
                 do ispec = ibool%ispec1,nspec
                    spectop = spectop .and. (ispec == nspec)
                    do inode = 1,ngll
                       ibnd = spectop .and. (inode == ngll)
                       if(ibnd) cycle
                       rt   =  r(inode,ispec)
                       mut  = mu(inode,ispec)
                       i = ibool%get(1,inode,ispec)
                       k = kd+1
                       a(k,i) = a(k,i) + zeta2*zeta2m2*mut*w(inode)*jac(ispec)
                       do jnode = inode,ngll
                          jbnd = spectop .and. (jnode == ngll)
                          if(jbnd) cycle
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
            class default
               stop 'build_toroidal_matrix: invalid mesh'
            end select
          end associate
       end do
    end do
    if(factor_local) then
       call dpbtrf('U',ndim,kd,mat%a,ldab,info)
       call error(info /= 0,'build_toroidal_matrix','problem with factorisation')
       mat%factorised = .true.
    end if
    return
  end function build_toroidal_matrix


  real(dp) function toroidal_start(mesh,l,eps_in) result(rs)
    type(spherical_model_mesh), intent(in) :: mesh
    integer(i4b), intent(in) :: l
    real(dp), intent(in), optional :: eps_in
    real(dp), parameter :: eps_def = 1.0e-8_dp
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
    fluid = .false.
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

    ibool%ndim = count+1
    ibool%ngll = ngllm
    
    return
  end function build_boolean_spheroidal


  type(radial_matrix) function build_spheroidal_matrix(mesh,l,factor) result(mat)
    type(spherical_model_mesh), intent(in) :: mesh
    integer(i4b) :: l
    logical, intent(in), optional :: factor

    logical :: spectop,ibnd,jbnd,factor_local,centre
    integer(i4b) :: ndim,kd,ldab,isection,ilayer,inode,  &
                    ispec,jnode,knode,ngll,nspec,i,j,k,info,ispec1, &
                    i1,i2,i3,j1,j2,j3,nlayers,ilayer1
    real(dp) :: tmp1,tmp2,zeta2,zeta2m2,zetac,ifpibigg,rt,mut,kapt,gt, &
                rhot,drhot,fac

    if(present(factor)) then
       factor_local = factor
    else
       factor_local = .true.
    end if       
    call error(l < 1,'build_spheroidal_matrix','l<1')    
    mat%l = l
    zeta2 = l*(l+1)
    zeta2m2 = zeta2-2.0_dp
    zetac = (4.0_dp/3.0_dp)*zeta2-2.0_dp
    ifpibigg = 1.0_dp/(4.0_dp*pi*bigg)
    mat%ibool = build_boolean_spheroidal(mesh,spheroidal_start(mesh,l))
    if(l == 1) mat%ibool%ndim = mat%ibool%ndim-1
    ndim = mat%ibool%ndim
    kd = 3*(ngll-1)+2
    ldab = kd+1
    if(allocated(mat%a)) deallocate(mat%a); allocate(mat%a(ldab,ndim))
    mat%a = 0.0_dp    

    centre = (mesh%r1 == 0.0_dp) .and. (mat%ibool%isection1 == 1)
    if(centre) then
       centre = centre .and. (mat%ibool%section(1)%ilayer1 == 1)
       if(centre) then
          centre = centre .and. (mat%ibool%section(1)%layer(1)%ispec1 == 1)
       end if
    end if
     
    if(.not.centre) then
       isection = mat%ibool%isection1
       ilayer = mat%ibool%section(isection)%ilayer1
       ispec = mat%ibool%section(isection)%layer(ilayer)%ispec1
       rt = mesh%section(isection)%layer(ilayer)%r1
       associate(layer => mesh%section(isection)%layer(ilayer), &
                 ibool => mat%ibool%section(isection)%layer(ilayer), &
                 a => mat%a)         
         select type(layer)
         class is(spherical_solid_elastic_layer_mesh)
            i = ibool%get(3,1,ispec)
            j = i
         class is(spherical_fluid_elastic_layer_mesh)
            i = ibool%get(1,1,ispec)
            j = i
         class default
            stop 'build_spheroidal_mesh: invalid mesh'            
         end select
         k = kd+1+i-j
         a(k,j) = a(k,j) + l*ifpibigg*rt
       end associate
    end if
    
    
    do isection = mat%ibool%isection1,mesh%nsections
       spectop = (l == 1) .and. (isection == mesh%nsections)
       ilayer1 = mat%ibool%section(isection)%ilayer1
       nlayers = mesh%section(isection)%nlayers
       do ilayer = ilayer1,nlayers
          spectop = spectop .and. (ilayer == nlayers)

          associate(layer => mesh%section(isection)%layer(ilayer),      &
                    ibool => mat%ibool%section(isection)%layer(ilayer), &
                    a => mat%a)

            ngll  = layer%ngll
            nspec = layer%nspec
            ispec1 = ibool%ispec1
            
            select type(layer)
               
            class is(spherical_solid_elastic_layer_mesh)

               associate(r     => layer%r,     &
                         rho   => layer%rho,   &
                         mu    => layer%mu,    &
                         kappa => layer%kappa, &                                        
                         g     => layer%g,     &
                         jac   => layer%jac,   &
                         hp    => layer%hp,    &                         
                         w     => layer%w)
                 
                 do ispec = ispec1,nspec
                    spectop = spectop .and. (ispec == nspec)
                    do inode = 1,ngll
                       ibnd = spectop .and. (inode == ngll)

                       rhot =   rho(inode,ispec)
                       mut  =    mu(inode,ispec)
                       kapt = kappa(inode,ispec)
                       
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
                       
                       
                       if(.not. ibnd) then
                          ! phi-v'
                          i = ibool%get(3,inode,ispec)
                          j = ibool%get(2,inode,ispec)
                          if(i <= j) then
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
                          jbnd = spectop .and. (jnode == ngll)


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
                          if(.not.jbnd) then
                             i = ibool%get(1,inode,ispec)
                             j = ibool%get(3,jnode,ispec)
                             k = kd+1+i-j
                             rt = r(inode,ispec)
                             rhot = rho(inode,ispec)
                             a(k,j) = a(k,j) + rhot*rt*rt*hp(inode,jnode)*w(inode)
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


                          ! v-u'
                          i = ibool%get(2,inode,ispec)
                          j = ibool%get(1,jnode,ispec)
                          if(i <= j) then
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
                          
                          ! phi-u'
                          if(.not.ibnd) then
                             i = ibool%get(3,inode,ispec)
                             j = ibool%get(1,jnode,ispec)
                             if(i <= j) then
                                k = kd+1+i-j
                                rt = r(jnode,ispec)
                                rhot = rho(jnode,ispec)
                                a(k,j) = a(k,j) + rhot*rt*rt*hp(jnode,inode)*w(jnode)
                             end if
                          end if
                             
                          ! phi-phi'
                          if(.not.ibnd .and. .not.jbnd) then
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


               
               
            class is(spherical_fluid_elastic_layer_mesh)

               associate(r    => layer%r,     &
                         drho => layer%drho,  &
                         g    => layer%g,     &
                         jac  => layer%jac,   &
                         hp   => layer%hp,    &                         
                         w    => layer%w)
                 
                 do ispec = ispec1,nspec
                    spectop = spectop .and. (ispec == nspec)
                    
                    do inode = 1,ngll
                       ibnd = spectop .and. (inode == ngll)                       
                       if(ibnd) cycle

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
                          jbnd = spectop .and. (jnode == ngll)
                          if(jbnd) cycle

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

               
            class default
               stop 'build_spheroidal_mesh: invalid mesh'
               
            end select
            
            
          end associate
            
       end do

       ! deal with boundary terms if needed
       if(isection < mesh%nsections) then
          ngll = mesh%section(isection)%layer(nlayers)%ngll
          nspec = mesh%section(isection)%layer(nlayers)%nspec
          associate(layer1 => mesh%section(isection)%layer(nlayers),      &
                    layer2 => mesh%section(isection+1)%layer(1),          &
                    a => mat%a)

            select type(layer1)
               
            class is(spherical_solid_elastic_layer_mesh)
               select type(layer2)
               class is(spherical_solid_elastic_layer_mesh)
               class is(spherical_fluid_elastic_layer_mesh)
                  associate(ibool =>  mat%ibool%section(isection)%layer(nlayers))                  
                    ! solid-fluid boundary                  
                    rt   = layer2%r(1,1)
                    rhot = layer2%rho(1,1)
                    gt   = layer2%g(1,1)
                    ! u-u'
                    i = ibool%get(1,ngll,nspec)
                    j = i
                    k = kd+1+i-j
                    a(k,j) = a(k,j) - rt*rt*gt*rhot
                    ! u-phi'
                    i = ibool%get(1,ngll,nspec)
                    j = ibool%get(3,ngll,nspec)
                    k = kd+1+i-j
                    a(k,j) = a(k,j) - rt*rt*rhot
                  end associate
               end select
               
            class is(spherical_fluid_elastic_layer_mesh)                
               select type(layer2)
               class is(spherical_solid_elastic_layer_mesh)
                  associate(ibool =>  mat%ibool%section(isection+1)%layer(1))                  
                    ! fluid-solid boundary
                    rt   = layer1%r(ngll,nspec)
                    rhot = layer1%rho(ngll,nspec)
                    gt   = layer1%g(ngll,nspec)
                    ! u-u'
                    i = ibool%get(1,1,1)
                    j = i
                    k = kd+1+i-j
                    a(k,j) = a(k,j) + rt*rt*gt*rhot
                    ! u-phi'
                    i = ibool%get(1,1,1)
                    j = ibool%get(3,1,1)
                    k = kd+1+i-j
                    a(k,j) = a(k,j) + rt*rt*rhot
                  end associate
               class is(spherical_fluid_elastic_layer_mesh)
               end select
               
            end select
            
          end associate
          
       end if
       
    end do


    ! add in DNT terms
    if(l > 1) then
       isection = mesh%nsections
       ilayer = mesh%section(ilayer)%nlayers
       associate(layer => mesh%section(isection)%layer(ilayer),      &
                 ibool => mat%ibool%section(isection)%layer(ilayer), &
                 a => mat%a)
         select type(layer)
         class is(spherical_solid_elastic_layer_mesh)
            nspec = layer%nspec
            ngll = layer%ngll
            i = ibool%get(3,ngll,nspec)
            j = i
            k = kd+1+i-j
            rt = layer%r(ngll,nspec)
            a(k,j) = a(k,j)+(l+1)*ifpibigg*rt
         end select
       end associate
    end if
    

    if(factor_local) then
       call dpbtrf('U',ndim,kd,mat%a,ldab,info)
       call error(info /= 0,'build_spheroidal_matrix','problem with factorisation')
       mat%factorised = .true.
    end if
    
    return
  end function build_spheroidal_matrix
  

  real(dp) function spheroidal_start(mesh,l,eps_in) result(rs)
    type(spherical_model_mesh), intent(in) :: mesh
    integer(i4b), intent(in) :: l
    real(dp), intent(in), optional :: eps_in
    real(dp), parameter :: eps_def = 1.0e-8_dp
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


  
end module module_SEM_matrix

