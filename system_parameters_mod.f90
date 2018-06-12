!!$ This module defines the available parameters relevant for
!!$ the evolution system

module system_parameters_mod

  use real_type_mod

!!$ Define the available parameters and their default values.

  real(kind=wp) :: &
       a = 1.0_wp, &
       burgers_xl = -1.0_wp / 3.0_wp, burgers_xr = 1.0_wp / 3.0_wp, &
       burgers_ql = -1.0_wp, burgers_qc = 1.0_wp, burgers_qr = -1.0_wp, &
       scalarwave_sigma = 0.3_wp, &
       euler_gamma = 1.4_wp, &
       euler_x_l   = 0.0_wp, euler_x_r = 0.0_wp, &
       euler_rho_l =   1.0_wp, euler_v_l = 0.0_wp, euler_p_l = 1.0_wp, &
       euler_rho_m =   1.0_wp, euler_v_m = 0.0_wp, euler_p_m = 1.0_wp, &
       euler_rho_r = 0.125_wp, euler_v_r = 0.0_wp, euler_p_r = 0.1_wp
  character(len=50) :: system = "advection", &
                       advection_initial = "sine", &
                       scalarwave_initial = "tsgaussian"
  
!!$  Derived parameters

  logical :: advection, burgers, scalarwave, euler, &
             advection_sine, advection_squarewave, advection_bell, &
             scalarwave_tsgaussian, scalarwave_separatinggaussian

  integer :: nvariables

end module system_parameters_mod
