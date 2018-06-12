!!$ This module defines the available parameters relevant for
!!$ the evolution method

module method_parameters_mod

  use real_type_mod

!!$ Define the available parameters and their default values.

  integer :: &
       centred_differencing_order = 2
  character(len=99) :: boundary = "periodic", &
                       rkmethod = "RK2", &
                       spatial_differencing = "TVD", &
                       limiter = "minmod", &
                       reconstruct_primitive = "yes", &
                       reconstruction_method = "TVD", &
                       method_type = "FVolume", &
                       riemann_solver = "HLLE", &
                       variable_timestep = "no"

  logical :: periodic, outflow, reflection, &
             RK2, RK4, &
             centred_differencing, TVD_differencing, &
             limiter_minmod, limiter_mc, &
             reconstruct_primitive_vars, &
             TVD_Reconstruct, WENO_Reconstruct, &
             Finite_Volume, Finite_Difference, &
             rs_HLLE, rs_HLLC, &
             variable_dt

end module method_parameters_mod
