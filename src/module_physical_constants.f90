module module_physical_constants

  use module_constants
  implicit none


  ! constants for non-dimensionalisation
  real(dp), parameter :: length_norm   = 6371000.0_dp
  real(dp), parameter :: mass_norm     = 5.972e24_dp
  real(dp), parameter :: time_norm     = 3600.0_dp
  real(dp), parameter :: frequency_norm = 1.0_dp/time_norm
  real(dp), parameter :: velocity_norm = length_norm/time_norm
  real(dp), parameter :: acceleration_norm = length_norm/time_norm**2
  real(dp), parameter :: density_norm = mass_norm/length_norm**3
  real(dp), parameter :: force_norm = mass_norm*length_norm/time_norm**2
  real(dp), parameter :: gravitational_potential_norm = length_norm**2/time_norm**2
  real(dp), parameter :: gravitational_constant_norm = length_norm**3/(mass_norm*time_norm**2)
  real(dp), parameter :: modulus_norm = mass_norm/(time_norm**2*length_norm)
  real(dp), parameter :: viscosity_norm = mass_norm/(time_norm*length_norm)
  real(dp), parameter :: action_norm = length_norm**2*mass_norm/time_norm
  real(dp), parameter :: load_norm = mass_norm/(length_norm**2)
  
  ! values for some useful physical constants in SI units
  real(dp), parameter :: bigg = 6.6743e-11_dp/gravitational_constant_norm

  
end module module_physical_constants
