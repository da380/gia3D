module module_sea_level

  use module_constants
  use module_physical_constants
  use module_spherical_harmonics
  use module_linear_solver
  implicit none

contains


  !=================================================================!
  !                      elastic fingerprint code                   !
  !=================================================================!


  subroutine fingerprint(grid,lln,tln,ice1,ice2,sl1,sl2)
    type(gauss_legendre_grid), intent(in) :: grid
    type(love_number), dimension(1:grid%lmax), intent(in) :: lln
    type(love_number), intent(in) :: tln
    type(real_scalar_gauss_legendre_field), intent(in)  :: ice1
    type(real_scalar_gauss_legendre_field), intent(in)  ::  sl1
    type(real_scalar_gauss_legendre_field), intent(in)  :: ice2
    type(real_scalar_gauss_legendre_field), intent(inout) ::  sl2


    
    
    return
  end subroutine fingerprint


  !=================================================================!
  !                     general utility routines                    !
  !=================================================================!


  subroutine ocean_function(sl,ice,ofun)
    type(real_scalar_gauss_legendre_field), intent(in)  :: sl
    type(real_scalar_gauss_legendre_field), intent(in)  :: ice
    type(real_scalar_gauss_legendre_field), intent(inout)  :: ofun
    where(rho_water*sl%rdata - rho_ice*ice%rdata > 0.0_dp)
       ofun%rdata = 1.0_dp
    elsewhere
       ofun%rdata = 0.0_dp
    end where    
    return
  end subroutine ocean_function

  subroutine ocean_function_mask(sl,ice,fun,comp)
    type(real_scalar_gauss_legendre_field), intent(in)  :: sl
    type(real_scalar_gauss_legendre_field), intent(in)  :: ice
    type(real_scalar_gauss_legendre_field), intent(inout)  :: fun
    logical, intent(in), optional :: comp
    logical :: complement,test
    if(present(comp)) then
       complement = comp
    else
       complement  = .true.
    end if
    where(rho_water*sl%rdata - rho_ice*ice%rdata <= 0.0_dp .eqv. .not.complement)
       fun%rdata = 0.0_dp
    end where    
    return
  end subroutine ocean_function_mask


  real(dp) function ocean_area(grid,sl,ice) result(area)
    type(gauss_legendre_grid), intent(in) :: grid
    type(real_scalar_gauss_legendre_field), intent(in)  :: sl
    type(real_scalar_gauss_legendre_field), intent(in)  :: ice
    integer(i4b) :: ith,iph,nth,nph
    real(dp) :: atmp,fac
    nth = grid%nth
    nph = grid%nph
    area = 0.0_dp
    do ith = 1,nth
       atmp = 0.0_dp
       do iph = 1,nph          
          if(rho_water*sl%get(iph,ith) - rho_ice*ice%get(iph,ith) > 0.0_dp) then
             atmp = atmp + 1.0_dp             
          end if          
       end do
       fac = grid%w(ith)*twopi/nph
       area = area + atmp*fac
    end do
    return
  end function ocean_area
  
  subroutine load(sl1,ice1,sl2,ice2,sigma)
    type(real_scalar_gauss_legendre_field), intent(in)  :: sl1
    type(real_scalar_gauss_legendre_field), intent(in)  :: ice1
    type(real_scalar_gauss_legendre_field), intent(in)  :: sl2
    type(real_scalar_gauss_legendre_field), intent(in)  :: ice2
    type(real_scalar_gauss_legendre_field), intent(inout)  :: sigma

    integer(i4b) :: ith,iph
    real(dp) :: sig,f1,f2


    do ith = 1,sl1%nth
       do iph = 1,sl1%nph
          
          f1  = rho_water*sl1%get(iph,ith)
          f2  = rho_ice*ice1%get(iph,ith)
          if(f1 > f2) then
             sig = f1
          else
             sig = f2
          end if

          f1 = rho_water*sl2%get(iph,ith)
          f2 = rho_ice*ice2%get(iph,ith)
          if(f1 > f2) then
             sig = f1-sig
          else
             sig = f2-sig
          end if
          call sigma%set(iph,ith,sig)
       end do
    end do
    
    return
  end subroutine load

  
  
end module module_sea_level
