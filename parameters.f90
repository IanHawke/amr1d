!!$ Read the parameters from input.par

subroutine read_parameters

  use real_type_mod
  use grid_parameters_mod
  use system_parameters_mod
  use method_parameters_mod

  implicit none

!!$ Define the namelist for the available parameters.
  namelist /grid_params/ &
       xmin, xmax, &
       last_time, &
       base_n_points, &
       ghostzones, &
       n_convergence, &
       max_refinement_levels, &
       amr_error, &
       amr_buffer_width
  
  namelist /system_params/ &
       system, &
       a, &
       advection_initial, &
       burgers_xl, burgers_xr, &
       burgers_ql, burgers_qc, burgers_qr, &
       scalarwave_sigma, &
       scalarwave_initial, &
       euler_gamma, &
       euler_x_l, euler_x_r, &
       euler_rho_l, euler_v_l, euler_p_l, &
       euler_rho_m, euler_v_m, euler_p_m, &
       euler_rho_r, euler_v_r, euler_p_r, &
       nvariables
  
  namelist /method_params/ &
       courant, &
       boundary, &
       rkmethod, &
       spatial_differencing, &
       centred_differencing_order, &
       limiter, &
       reconstruct_primitive, &
       reconstruction_method, &
       method_type, &
       riemann_solver, &
       variable_timestep
  
  
!!$ Open the parameter file.
  open (3, file='input.par', status = 'old' )
  
!!$ Read the namelist from the file. NOTE: not all parameters have to
!!$ be in the parameter file.
  read (3, nml = system_params)
  read (3, nml = grid_params)
  read (3, nml = method_params)
  
!!$ Close the file.
  close(3)
  
!!$ Set derived parameters
  
  base_dx = (xmax - xmin) / dble(base_n_points)
  base_dt = courant * base_dx
  dtodx = courant
  crossing_time = (xmax - xmin) / a
  
  output_exact_and_errors = .true.
  advection = .false.
  burgers = .false.
  scalarwave = .false.
  nvariables = -1
  if (system == 'advection') then
    advection = .true.
    nvariables = 1
  else if (system == 'burgers') then
    burgers = .true.
    nvariables = 1
  else if (system == 'scalarwave') then
    scalarwave = .true.
    nvariables = 3
  else if (system == 'euler') then
    euler = .true.
    nvariables = 6
    output_exact_and_errors = .false.
  else
    write(*,*) 'System to be evolved not recognized!'
    stop
  end if
  
  periodic = .false.
  outflow = .false.
  reflection = .false.
  if (boundary == 'periodic') then
    periodic = .true.
  else if (boundary == 'outflow') then
    outflow = .true.
  else if (boundary == 'reflection') then
    reflection = .true.
  else
    write(*,*) 'Boundary condition not recognized!'
    stop
  end if
  
  RK2 = .false.
  RK4 = .false.
  if (rkmethod == 'RK2') then
    RK2 = .true.
  else if (rkmethod == 'RK4') then
    RK4 = .true.
  else
    write(*,*) 'Time evolution (ODE) method not recognized!'
    stop
  end if
  
  advection_sine = .false.
  advection_squarewave = .false.
  advection_bell = .false.
  if (advection_initial == 'sine') then
    advection_sine = .true.
  else if (advection_initial == 'squarewave') then
    advection_squarewave = .true.
  else if (advection_initial == 'bell') then
    advection_bell = .true.
  else
    write(*,*) 'Initial data type for advection equation not recognized!'
    stop
  end if
  
  scalarwave_tsgaussian = .false.
  scalarwave_separatinggaussian = .false.
  if (scalarwave_initial == 'tsgaussian') then
    scalarwave_tsgaussian = .true.
  else if (scalarwave_initial == 'separatinggaussian') then
    scalarwave_separatinggaussian = .true.
  else
    write(*,*) 'Initial data type for scalar wave equation not recognized!'
    stop
  end if
  
  centred_differencing = .false.
  TVD_differencing = .false.
  if (spatial_differencing == 'centred') then
    centred_differencing = .true.
  else if (spatial_differencing == 'TVD') then
    TVD_differencing = .true.
  else
    write(*,*) 'Spatial differencing method not recognized!'
    stop
  end if
  
  TVD_Reconstruct  = .false.
  WENO_Reconstruct = .false.
  if (reconstruction_method == 'TVD') then
    TVD_Reconstruct = .true.
  else if (reconstruction_method == 'WENO') then
    WENO_Reconstruct = .true.
  else
    write(*,*) 'Reconstruction method not recognized!'
    stop
  end if
  
  Finite_Volume  = .false.
  Finite_Difference = .false.
  if (method_type == 'FVolume') then
    Finite_Volume = .true.
  else if (method_type == 'FDifference') then
    Finite_Difference = .true.
  else
    write(*,*) &
      'Finite volume or differece? method_type not recognized!'
    stop
  end if
  
  limiter_minmod = .false.
  limiter_mc = .false.
  if (limiter == 'minmod') then
    limiter_minmod = .true.
  else if (limiter == 'mc') then
    limiter_mc = .true.
  else
    write(*,*) 'TVD limiter not recognized!'
    stop
  end if

  rs_HLLE = .false.
  rs_HLLC = .false.
  if (riemann_solver == 'HLLE') then
    rs_HLLE = .true.
  else if (riemann_solver == 'HLLC') then
    rs_HLLC = .true.
  else
    write(*,*) 'Riemann solver not recognized!'
    stop
  end if

  reconstruct_primitive_vars = .true.
  if (reconstruct_primitive == 'no') then
    reconstruct_primitive_vars = .false.
  end if

  variable_dt = .false.
  if (variable_timestep == 'yes') then
    variable_dt = .true.
  end if
  
  write(*,'(a11,es8.1,a1,es8.1,a1)') 'Domain is [',xmin,',',xmax,']'
  write(*,'(a10,i6,a20,es8.1)') 'There are ', base_n_points, &
       ' cells with spacing ', base_dx
  if (advection) then
    write(*,'(a16,es8.1,a17,es8.1)') 'Advection speed ', a, &
         ', crossing time ', crossing_time
    write(*, '(a16)', advance='NO') 'Initial data is '
    if (advection_sine) then
      write(*, '(a11)') 'a sine wave'
    else if (advection_squarewave) then
      write(*, '(a13)') 'a square wave'
    else if (advection_bell) then
      write(*, '(a17)') 'the inverted bell'
    else
      write(*, '(a10)') 'unknown?!?'
    end if
  else if (burgers) then
    write(*,*) 'Evolving Burger''s equation'
  end if
  if (centred_differencing) then
    write(*,'(a30,i2)') 'Centred differencing, order = ', &
         centred_differencing_order
  else if (TVD_differencing) then
    if (TVD_Reconstruct) then
      write(*,'(a26)', advance='NO') 'TVD differencing, limiter '
      if (limiter_minmod) then
        write(*,'(a6)') 'minmod'
      else if (limiter_mc) then
        write(*,'(a11)') 'Van Leer MC'
      else
        write(*,'(a7)') 'Unknown'
      end if
    else if (WENO_Reconstruct) then
      write(*, '(a19)') 'WENO Reconstruction'
    end if
  end if
  if (periodic) then
    write(*,'(a27)') 'Boundary condition periodic'
  else if (outflow) then
    write(*,'(a26)') 'Boundary condition outflow'
  else if (reflection) then
    write(*,'(a26)') 'Boundary condition reflection'
  else
    write(*,'(a24)') 'What boundary condition?'
  end if
  if (RK2) then
    write(*,'(a25)') 'Using RK2 time integrator'
  else if (RK4) then
    write(*,'(a25)') 'Using RK4 time integrator'
  end if
  
end subroutine read_parameters


