program test_spherical_model

  use module_constants
  use module_spherical_model
  implicit none

  integer(i4b) :: ilayer,ilayer1,ilayer2,ir,nr = 50,iregion
  real(dp) :: r,r1,r2,dr
  class(spherical_elastic_model), allocatable :: model
  class(spherical_elastic_model), allocatable :: deck_model

  model =  set_elastic_PREM()
  call model%write(200,'prem.200')!,isotropic = .true.)
  deck_model = set_elastic_deck_model('prem.200') 


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
  

  
end program test_spherical_model


