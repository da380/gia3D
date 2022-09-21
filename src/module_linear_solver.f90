module module_linear_solver

  use module_constants
  use module_error
  use module_spherical_harmonics
  use module_mesh
  use module_matrix
  use module_force
  implicit none
  

contains


  subroutine direct_elastic_solver(mesh,mat,sigma,psi,u,phi)
    type(spherical_model_mesh), intent(in) :: mesh
    type(radial_matrices), intent(in) :: mat
    type(real_scalar_spherical_harmonic_expansion), intent(in) :: sigma
    type(real_scalar_spherical_harmonic_expansion), intent(in) :: psi
    type(real_scalar_spherical_harmonic_expansion), intent(in) :: u
    type(real_scalar_spherical_harmonic_expansion), intent(in) :: phi

    integer(i4b) :: l,lmax,ndim,isection,ilayer,i1,i2,info
    real(dp), dimension(mat%ndim_sph,mat%lmax+1) :: b

    call error(.not.mat%factorised,'direct_elastic_solver','matrices not factorised')
    call error(.not.mat%spheroidals,'direct_elastic_solver','spheroidal matrics not set up')

    lmax = mat%lmax

    do l = 1,lmax

       ! build the force term
       
       ndim = mat%ndim_sph
       b(1:ndim,1:2*l+1) = 0.0_dp
       isection = mesh%nsections
       ilayer = mesh%section(isection)%nlayers
       associate(layer => mesh%section(isection)%layer(ilayer), & 
                 ibool => mat%sph(l)%ibool%section(isection)%layer(ilayer))

         i1 = sigma%rindex(l,0,1)
         i2 = sigma%rindex(l,l,1)
         call force_for_real_degree_l_load(layer,ibool,l,sigma%rdata(i1:i2),b(1:ndim,1:2*l+1))
         
       end associate

       i1 = psi%rindex(l,0,1)
       i2 = psi%rindex(l,l,1)
       call force_for_real_degree_l_tide(mesh,mat%sph(l)%ibool,l,psi%rdata(i1:i2),b(1:ndim,1:2*l+1))

       ! solve the linear system
       call dpbtrs('U',ndim,mat%sph(l)%kd,1,mat%sph(l)%a,mat%sph(l)%ldab,b,mat%ndim_sph,info)
       call error(info /= 0,'direct_elastic_solver','problem backsubstitution')  


       ! unpack the solution vector into u and phi
       
    end do
    
    return
  end subroutine direct_elastic_solver
  
end module module_linear_solver
