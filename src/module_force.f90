module module_force

  use module_constants
  use module_physical_constants
  use module_mesh
  use module_matrix
  implicit none

contains


  !==========================================================!
  !              forces associated with loading              !
  !==========================================================!
  
  subroutine force_for_unit_harmonic_load(layer,ibool,l,b)
    class(spherical_layer_mesh), intent(in) :: layer
    type(boolean_array_layer), intent(in) :: ibool
    integer(i4b), intent(in) :: l
    real(dp), dimension(:,:), intent(inout) :: b
    integer(i4b) :: i,inode,ispec
    real(dp) :: g,r
    inode = layer%ngll
    ispec = layer%nspec
    g = layer%g(inode,ispec)
    r = layer%r(inode,ispec)
    i = ibool%get(1,inode,ispec)
    b(i,1) = -g*r*r
    if(l > 1) then
       i = ibool%get(3,inode,ispec)
       b(i,1) = -r*r
    end if
    return
  end subroutine force_for_unit_harmonic_load


  subroutine force_for_real_degree_l_load(layer,ibool,l,sigma,b)
    class(spherical_layer_mesh), intent(in) :: layer
    type(boolean_array_layer), intent(in) :: ibool
    integer(i4b), intent(in) :: l
    complex(dpc), dimension(l+1), intent(in) :: sigma
    real(dp), dimension(:,:), intent(inout) :: b
    integer(i4b) :: i,j,inode,ispec
    real(dp) :: g,r
    inode = layer%ngll
    ispec = layer%nspec
    g = layer%g(inode,ispec)
    r = layer%r(inode,ispec)
    i = ibool%get(1,inode,ispec)
    b(i,1)         = -g*r*r*real(sigma(1))
    b(i,2:l+1)     = -g*r*r*real(sigma(2:l+1))
    b(i,l+2:2*l+1) = -g*r*r*imag(sigma(2:l+1))
    if(l > 1) then
       j = ibool%get(3,inode,ispec)
       b(j,:) = b(i,:)/g
    end if
    return
  end subroutine force_for_real_degree_l_load
  

  !==========================================================!
  !                forces associated with tides              !
  !==========================================================!


  subroutine force_for_unit_harmonic_tide(mesh,ibool,l,b)
    class(spherical_model_mesh), intent(in) :: mesh
    type(boolean_array), intent(in) :: ibool
    integer(i4b), intent(in) :: l
    real(dp), dimension(:,:), intent(inout) :: b
    
    integer(i4b) :: ilayer,ilayer1,nlayers,isection,isection1,nsections
    real(dp) :: rn,fac


    
    ! build up the force layer-by-layer
    isection1 = ibool%isection1
    nsections = mesh%nsections
    do isection = isection1,nsections
       ilayer1 = ibool%section(isection)%ilayer1
       nlayers = mesh%section(isection)%nlayers
       do ilayer = ilayer1,nlayers
          associate(layer => mesh%section(isection)%layer(ilayer),      &
                    ibool => ibool%section(isection)%layer(ilayer))
            select type(layer)               
            class is(spherical_solid_elastic_layer_mesh)
               call force_for_unit_harmonic_tide_solid_layer(layer,ibool,l,b)
            class is(spherical_fluid_elastic_layer_mesh)
               call force_for_unit_harmonic_tide_fluid_layer(layer,ibool,l,b)
            class default
               stop 'force_for_unit_harmonic_tide: invalid mesh'
            end select
          end associate
       end do
    end do

    ! non-dimensionalise
    rn = mesh%r2
    fac = (rn**l)*gravitational_potential_norm
    b = b/fac
    
    return
  end subroutine force_for_unit_harmonic_tide

  subroutine force_for_unit_harmonic_tide_solid_layer(layer,ibool,l,b)
    class(spherical_solid_elastic_layer_mesh), intent(in) :: layer
    type(boolean_array_layer), intent(in) :: ibool
    integer(i4b), intent(in) :: l
    real(dp), dimension(:,:), intent(inout) :: b

    integer(i4b) :: inode,ngll,ispec,ispec1,nspec,i
    real(dp) :: zeta2,rt,rhot

    ! set some parameters
    zeta2 = l*(l+1)
    ngll = layer%ngll
    nspec = layer%nspec
    ispec1 = ibool%ispec1

    ! deal with volumetric terms
    associate(r   => layer%r,   &
              rho => layer%rho, &
              jac => layer%jac, &
              w  =>  layer%w)

      do ispec = ispec1,nspec

         do inode = 1,ngll
         
            rt   = r(inode,ispec)
            rt   = rt**(l+1)
            rhot = rho(inode,ispec)
            
            ! u'
            i = ibool%get(1,inode,ispec)
            b(i,1) = b(i,1) - l*rhot*rt*w(inode)*jac(ispec)
            
            ! v'
            i = ibool%get(2,inode,ispec)
            b(i,1) = b(i,1) - zeta2*rhot*rt*w(inode)*jac(ispec)
            
         end do

      end do
         
    end associate

    ! deal with fluid-solid terms if needed
    if(layer%below_fluid .and. ispec1 == 1) then

       inode = 1
       ispec = 1
       rt = layer%r(inode,ispec)
       rt = rt**(l+2)
       rhot = layer%below_fluid_density
       i = ibool%get(1,inode,ispec)
       b(i,1) = b(i,1) - rhot*rt
       
    end if

    ! deal with solid-fluid terms if needed
    if(layer%above_fluid) then

       inode = ngll
       ispec = nspec
       rt = layer%r(inode,ispec)
       rt = rt**(l+2)
       rhot = layer%above_fluid_density
       i = ibool%get(1,inode,ispec)
       b(i,1) = b(i,1) + rhot*rt
       
    end if
    
    return
  end subroutine force_for_unit_harmonic_tide_solid_layer

  
  subroutine force_for_unit_harmonic_tide_fluid_layer(layer,ibool,l,b)
    class(spherical_fluid_elastic_layer_mesh), intent(in) :: layer
    type(boolean_array_layer), intent(in) :: ibool
    integer(i4b), intent(in) :: l
    real(dp), dimension(:,:), intent(inout) :: b


    integer(i4b) :: inode,ngll,ispec,ispec1,nspec,i
    real(dp) :: rt,drhot,gt,fac

    ! set some parameters
    ngll = layer%ngll
    nspec = layer%nspec
    ispec1 = ibool%ispec1

    ! deal with volumetric terms
    associate(r    => layer%r,    &
              drho => layer%drho, &
              g    => layer%g,    &
              jac  => layer%jac,  &
              w    => layer%w)

      do ispec = ispec1,nspec

         do inode = 1,ngll
         
            rt    = r(inode,ispec)
            rt    = rt**(l+2)
            drhot = drho(inode,ispec)
            gt    = g(inode,ispec)
            
            ! phi'
            i = ibool%get(1,inode,ispec)
            b(i,1) = b(i,1) - drhot*rt*w(inode)*jac(ispec)/gt
            
            
         end do

      end do
         
    end associate

    
    return
  end subroutine force_for_unit_harmonic_tide_fluid_layer


  subroutine force_for_real_degree_l_tide(mesh,ibool,l,psi,b)
    class(spherical_model_mesh), intent(in) :: mesh
    type(boolean_array), intent(in) :: ibool
    integer(i4b), intent(in) :: l
    complex(dpc), dimension(l+1), intent(in) :: psi
    real(dp), dimension(:,:), intent(inout) :: b
    
    integer(i4b) :: ilayer,ilayer1,nlayers,isection,isection1,nsections
    
    ! build up the force layer-by-layer
    isection1 = ibool%isection1
    nsections = mesh%nsections
    do isection = isection1,nsections
       ilayer1 = ibool%section(isection)%ilayer1
       nlayers = mesh%section(isection)%nlayers
       do ilayer = ilayer1,nlayers
          associate(layer => mesh%section(isection)%layer(ilayer),      &
                    ibool => ibool%section(isection)%layer(ilayer))
            select type(layer)               
            class is(spherical_solid_elastic_layer_mesh)
               call force_for_real_degree_l_tide_solid_layer(layer,ibool,l,psi,b)
            class is(spherical_fluid_elastic_layer_mesh)
               call force_for_real_degree_l_tide_fluid_layer(layer,ibool,l,psi,b)
            class default
               stop 'force_for_real_degree_l_tide: invalid mesh'
            end select
          end associate
       end do
    end do
    
    return
  end subroutine force_for_real_degree_l_tide

  subroutine force_for_real_degree_l_tide_solid_layer(layer,ibool,l,psi,b)
    class(spherical_solid_elastic_layer_mesh), intent(in) :: layer
    type(boolean_array_layer), intent(in) :: ibool
    integer(i4b), intent(in) :: l
    complex(dpc), dimension(l+1), intent(in) :: psi
    real(dp), dimension(:,:), intent(inout) :: b

    integer(i4b) :: inode,ngll,ispec,ispec1,nspec,i
    real(dp) :: zeta2,rt,rhot,fac

    ! set some parameters
    zeta2 = l*(l+1)
    ngll = layer%ngll
    nspec = layer%nspec
    ispec1 = ibool%ispec1

    ! deal with volumetric terms
    associate(r   => layer%r,   &
              rho => layer%rho, &
              jac => layer%jac, &
              w  =>  layer%w)

      do ispec = ispec1,nspec

         do inode = 1,ngll
         
            rt   = r(inode,ispec)
            rt   = rt**(l+1)
            rhot = rho(inode,ispec)
            
            ! u'
            i = ibool%get(1,inode,ispec)
            fac = -l*rhot*rt*w(inode)*jac(ispec)
            b(i,1)         = b(i,1)         + real(psi(1))*fac
            b(i,2:l+1)     = b(i,2:l+1)     + real(psi(2:l+1))*fac
            b(i,l+2:2*l+1) = b(i,l+2:2*l+1) + imag(psi(2:l+1))*fac
            
            ! v'
            i = ibool%get(2,inode,ispec)
            fac = - zeta2*rhot*rt*w(inode)*jac(ispec)
            b(i,1)         = b(i,1)         + real(psi(1))*fac
            b(i,2:l+1)     = b(i,2:l+1)     + real(psi(2:l+1))*fac
            b(i,l+2:2*l+1) = b(i,l+2:2*l+1) + imag(psi(2:l+1))*fac
            
         end do

      end do
         
    end associate

    ! deal with fluid-solid terms if needed
    if(layer%below_fluid .and. ispec1 == 1) then

       inode = 1
       ispec = 1
       rt = layer%r(inode,ispec)
       rt = rt**(l+2)
       rhot = layer%below_fluid_density
       i = ibool%get(1,inode,ispec)
       fac = - rhot*rt
       b(i,1)         = b(i,1)         + real(psi(1))*fac
       b(i,2:l+1)     = b(i,2:l+1)     + real(psi(2:l+1))*fac
       b(i,l+2:2*l+1) = b(i,l+2:2*l+1) + imag(psi(2:l+1))*fac
       
    end if

    ! deal with solid-fluid terms if needed
    if(layer%above_fluid) then

       inode = ngll
       ispec = nspec
       rt = layer%r(inode,ispec)
       rt = rt**(l+2)
       rhot = layer%above_fluid_density
       i = ibool%get(1,inode,ispec)
       fac = rhot*rt
       b(i,1)         = b(i,1)         + real(psi(1))*fac
       b(i,2:l+1)     = b(i,2:l+1)     + real(psi(2:l+1))*fac
       b(i,l+2:2*l+1) = b(i,l+2:2*l+1) + imag(psi(2:l+1))*fac
       
    end if
    
    return
  end subroutine force_for_real_degree_l_tide_solid_layer

  
  subroutine force_for_real_degree_l_tide_fluid_layer(layer,ibool,l,psi,b)
    class(spherical_fluid_elastic_layer_mesh), intent(in) :: layer
    type(boolean_array_layer), intent(in) :: ibool
    integer(i4b), intent(in) :: l
    complex(dpc), dimension(l+1), intent(in) :: psi
    real(dp), dimension(:,:), intent(inout) :: b


    integer(i4b) :: inode,ngll,ispec,ispec1,nspec,i
    real(dp) :: rt,drhot,gt,fac

    ! set some parameters
    ngll = layer%ngll
    nspec = layer%nspec
    ispec1 = ibool%ispec1

    ! deal with volumetric terms
    associate(r    => layer%r,    &
              drho => layer%drho, &
              g    => layer%g,    &
              jac  => layer%jac,  &
              w    => layer%w)

      do ispec = ispec1,nspec

         do inode = 1,ngll
         
            rt    = r(inode,ispec)
            rt    = rt**(l+2)
            drhot = drho(inode,ispec)
            gt    = g(inode,ispec)
            
            ! phi'
            i = ibool%get(1,inode,ispec)
            fac = - drhot*rt*w(inode)*jac(ispec)/gt
            b(i,1)         = b(i,1)         + real(psi(1))*fac
            b(i,2:l+1)     = b(i,2:l+1)     + real(psi(2:l+1))*fac
            b(i,l+2:2*l+1) = b(i,l+2:2*l+1) + imag(psi(2:l+1))*fac
                        
         end do

      end do
         
    end associate

    
    return
  end subroutine force_for_real_degree_l_tide_fluid_layer


  
  
end module module_force


