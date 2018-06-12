!!$ This module defines the available parameters and defines a subroutine to
!!$ read the parameters from input.par

module parameters_mod

  use real_type_mod

!!$ Define the available parameters and their default values.

  real(kind=wp) :: &
       xmin = -1.0_wp, xmax = 1.0_wp, &
       a = 1.0_wp, &
       courant = 0.25_wp, &
       last_time = 1.0_wp, &
       burgers_xl = -1.0_wp / 3.0_wp, burgers_xr = 1.0_wp / 3.0_wp, &
       burgers_ql = -1.0_wp, burgers_qc = 1.0_wp, burgers_qr = -1.0_wp, &
       euler_gamma = 1.4_wp
  integer :: &
       n_points = 100, &
       ghostzones = 1, &
       n_convergence = 1, &
       centred_differencing_order = 2
  character(len=10) :: system = "advection", &
                      boundary = "periodic", &
                      rkmethod = "RK2", &
                      advection_initial = "sine", &
                      spatial_differencing = "TVD"
  
!!$  Derived parameters

  real(kind=wp) :: dx, dt, dtodx, crossing_time
  logical :: advection, burgers, &
             periodic, outflow, &
             RK2, RK4, &
             advection_sine, advection_squarewave, advection_bell, &
             centred_differencing, TVD_differencing

contains

  subroutine read_parameters

    implicit none

!!$   Define the namelist for the available parameters.
    namelist /grid_params/ &
         xmin, xmax, &
         last_time, &
         n_points, &
         ghostzones, &
         n_convergence

    namelist /system_params/ &
         system, &
         a, &
         advection_initial, &
         burgers_xl, burgers_xr, &
         burgers_ql, burgers_qc, burgers_qr
         
    namelist /method_params/ &
         courant, &
         boundary, &
         rkmethod, &
         spatial_differencing, &
         centred_differencing_order


!!$   Open the parameter file.
    open (3, file='input.par', status = 'old' )

!!$   Read the namelist from the file. NOTE: not all parameters have to be
!!$   in the parameter file.
    read (3, nml = system_params)
    read (3, nml = grid_params)
    read (3, nml = method_params)

!!$   Close the file.
    close(3)

!!$    Set derived parameters

    dx = (xmax - xmin) / dble(n_points)
    dt = courant * dx
    dtodx = courant
    crossing_time = (xmax - xmin) / a

    advection = .false.
    burgers   = .false.
    if (system == 'advection') then
      advection = .true.
    else if (system == 'burgers') then
      burgers = .true.
    else
      write(*,*) 'System to be evolved not recognized!'
      stop
    end if

    periodic  = .false.
    outflow   = .false.
    if (boundary == 'periodic') then
      periodic = .true.
    else if (boundary == 'outflow') then
      outflow = .true.
    else
      write(*,*) 'Boundary condition not recognized!'
      stop
    end if

    RK2  = .false.
    RK4  = .false.
    if (rkmethod == 'RK2') then
      RK2 = .true.
    else if (rkmethod == 'RK4') then
      RK4 = .true.
    else
      write(*,*) 'Time evolution (ODE) method not recognized!'
      stop
    end if

    advection_sine         = .false.
    advection_squarewave   = .false.
    advection_bell         = .false.
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

    centred_differencing  = .false.
    TVD_differencing   = .false.
    if (spatial_differencing == 'centred') then
      centred_differencing = .true.
    else if (spatial_differencing == 'TVD') then
      TVD_differencing = .true.
    else
      write(*,*) 'Spatial differencing method not recognized!'
      stop
    end if
    
    write(*,'(a11,es8.1,a1,es8.1,a1)') 'Domain is [',xmin,',',xmax,']'
    write(*,'(a10,i6,a20,es8.1)') 'There are ', n_points, &
         ' cells with spacing ', dx
    if (advection) then
      write(*,'(a16,es8.1,a17,es8.1)') 'Advection speed ', a, &
           ', crossing time ', crossing_time
      write(*, '(a16)',advance='NO') 'Initial data is '
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
      write(*,*) 'TVD differencing'
    end if
    if (RK2) then
      write(*,'(a25)') 'Using RK2 time integrator'
    else if (RK4) then
      write(*,'(a25)') 'Using RK4 time integrator'
    end if
    
  end subroutine read_parameters

end module parameters_mod

