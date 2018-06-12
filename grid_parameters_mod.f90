!!$ This module defines the available parameters relevant for
!!$ the grid setup

module grid_parameters_mod

  use real_type_mod

!!$ Define the available parameters and their default values.

  real(kind=wp) :: &
       xmin = -1.0_wp, xmax = 1.0_wp, &
       last_time = 1.0_wp, &
       courant = 0.25_wp, &
       amr_error = 0.001_wp

  integer :: &
       base_n_points = 100, &
       ghostzones = 1, &
       n_convergence = 1, &
       max_refinement_levels = 1, &
       amr_buffer_width = 4
  
!!$  Derived parameters

  real(kind=wp) :: base_dx, base_dt, dtodx, crossing_time

  logical :: output_exact_and_errors

end module grid_parameters_mod
