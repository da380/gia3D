module module_linear_solver

  use module_constants
  use module_error
  use module_spherical_harmonics
  use module_mesh
  use module_matrix
  use module_force
  implicit none

  type love_number
     real(dp) :: ku = 0.0_dp
     real(dp) :: kv = 0.0_dp
     real(dp) :: kp = 0.0_dp
  end type love_number

contains

  subroutine make_love_numbers(model,l1,l2,lln,tln)
    type(spherical_model), intent(in) :: model
    integer(i4b), intent(in) :: l1,l2
    type(love_number), dimension(l1:l2), intent(out), optional :: lln
    type(love_number), dimension(l1:l2), intent(out), optional :: tln

    integer(i4b) :: ngll,l,ndim,kd,ldab,isection,ilayer,info,inode,ispec,i
    real(dp) :: drmax,u,v,p
    real(dp), dimension(:,:), allocatable :: b
    type(spherical_model_mesh) :: mesh
    type(radial_matrix) :: sphmat
    
    ! build the mesh
    ngll = 5
    drmax = 0.1_dp*model%r2/(l2+1)
    mesh = spherical_mesh(ngll,model,drmax)

    
    ! loop over the degrees    
    do l = l1,l2
       
       if(l == 0) cycle
       
       ! build the matrix
       sphmat = build_spheroidal_matrix(mesh,l)
       ndim = sphmat%ndim
       kd = sphmat%kd
       ldab = sphmat%ldab       
       if(allocated(b)) deallocate(b); allocate(b(ndim,1))       
       
       if(present(lln)) then
       
          b = 0.0_dp
          isection = mesh%nsections
          ilayer = mesh%section(isection)%nlayers
          associate(layer => mesh%section(isection)%layer(ilayer),         &
                    ibool => sphmat%ibool%section(isection)%layer(ilayer))            
            call force_for_unit_harmonic_load(layer,ibool,l,b)
            call dpbtrs('U',ndim,kd,1,sphmat%a,ldab,b,ndim,info)
            call check(info == 0,'test_matrix','problem with spheroidal substitution')  
            inode = layer%ngll
            ispec = layer%nspec
            i = ibool%get(1,inode,ispec)
            u = b(i,1)
            i = ibool%get(2,inode,ispec)
            v = b(i,1)
            if(l > 1) then
               i = ibool%get(3,inode,ispec)
               p = b(i,1)
            else
               p = 0.0_dp
            end if            
          end associate

          lln(l)%ku = u
          lln(l)%kv = v
          lln(l)%kp = p

       end if

       if(present(tln)) then

          b = 0.0_dp     
          call force_for_unit_harmonic_tide(mesh,sphmat%ibool,l,b)
          isection = mesh%nsections
          ilayer = mesh%section(isection)%nlayers
          associate(layer => mesh%section(isection)%layer(ilayer),         &
                    ibool => sphmat%ibool%section(isection)%layer(ilayer))            
            call dpbtrs('U',ndim,kd,1,sphmat%a,ldab,b,ndim,info)
            call check(info == 0,'test_matrix','problem with spheroidal substitution')  
            inode = layer%ngll
            ispec = layer%nspec
            i = ibool%get(1,inode,ispec)
            u = b(i,1)
            i = ibool%get(2,inode,ispec)
            v = b(i,1)
            if(l > 1) then
               i = ibool%get(3,inode,ispec)
               p = b(i,1)
            else
               p = 0.0_dp
            end if            
          end associate

          tln(l)%ku = u
          tln(l)%kv = v
          tln(l)%kp = p
          
       end if
       
    end do

    
    return
  end subroutine make_love_numbers

  



  
end module module_linear_solver
