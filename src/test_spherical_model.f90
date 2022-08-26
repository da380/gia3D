program test_spherical_model

  use module_constants
  use module_spherical_model
  implicit none

  integer(i4b) :: ilayer,ilayer1,ilayer2,ir,nr = 50,iregion
  real(dp) :: r,r1,r2,dr
  class(spherical_elastic_model), allocatable :: model
  class(spherical_elastic_model), allocatable :: deck_model

  class(spherical_maxwell_model), allocatable :: maxwell_model
  
  model = elastic_PREM()
  call model%write(200,'prem.200')!,isotropic = .true.)
  deck_model = elastic_deck_model('prem.200') 


  open(99,file='test_spherical_model.out')
  do ilayer = 1,model%nlayers
     r1 = model%layer_radius(1,ilayer)
     r2 = model%layer_radius(2,ilayer)
     dr = (r2-r1)/(nr-1)
     do ir = 1,nr
        r = r1+(ir-1)*dr
        write(99,*) r*length_norm,model%kappa(ilayer,r)-deck_model%kappa(ilayer,r)
     end do
  end do
  close(99)

  
  maxwell_model = maxwell_PREM(25.0_dp,21.0_dp,22.0_dp)
  call maxwell_model%write(200,'prem_maxwell.out')
  
  deallocate(maxwell_model)
  maxwell_model = maxwell_deck_model('prem_maxwell.out')
  call maxwell_model%write(200,'deck_maxwell.out')
  
  
end program test_spherical_model


