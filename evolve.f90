!!$Evolve using RK2 - Heuns

subroutine evolve_one_step(one_grid)

  use real_type_mod
  use method_parameters_mod
  use grid_mod
  use rhs_mod

  implicit none

  type(grid) :: one_grid

  real(kind=wp), dimension(:, :), pointer :: u, u_p
  real(kind=wp), dimension(:, :), pointer :: u_int1, u_int2, u_int3, rhs

  integer :: ierr

  integer :: nv, np, ngz
  real(kind=wp) :: dt, dx, char_max, tplus

  integer :: i, k

  nv  = one_grid%nvars
  np  = one_grid%n_points
  ngz = one_grid%ghostzones
  dx  = one_grid%dx
  
  dt  = one_grid%dt

!!$  RK2

  if (RK2) then

    u      => one_grid%u
    u_p    => one_grid%u_p

    u_int1 => one_grid%grid_data%u_int1

    rhs    => one_grid%grid_data%rhs

    call compute_RHS(ngz, u_p, one_grid)
    
    u_int1 = u_p + dt * rhs
    
    call boundary_conditions(nv, np, ngz, u_int1)

    if (associated(one_grid%parent)) then
      tplus = 1.0_wp
      call prolongate_boundaries(one_grid, u_int1, tplus)
    end if
    
    call conservative2primitive(nv, np, ngz, u_int1)
    
    call compute_RHS(ngz, u_int1, one_grid)
    
    u = 0.5_wp * (u_p + u_int1 + dt * rhs)

  else if (RK4) then

    u      => one_grid%u
    u_p    => one_grid%u_p

    u_int1 => one_grid%grid_data%u_int1
    u_int2 => one_grid%grid_data%u_int2
    u_int3 => one_grid%grid_data%u_int3

    rhs    => one_grid%grid_data%rhs

    call compute_RHS(ngz, u_p, one_grid)
    
    u_int1 = u_p + 0.5_wp * dt * rhs
    
    call boundary_conditions(nv, np, ngz, u_int1)

    if (associated(one_grid%parent)) then
      tplus = 0.5_wp
      call prolongate_boundaries(one_grid, u_int1, tplus)
    end if

    call conservative2primitive(nv, np, ngz, u_int1)
    
    call compute_RHS(ngz, u_int1, one_grid)
    
    u_int2 = u_p + 0.5_wp * dt * rhs
    
    call boundary_conditions(nv, np, ngz, u_int2)

    if (associated(one_grid%parent)) then
      tplus = 0.5_wp
      call prolongate_boundaries(one_grid, u_int2, tplus)
    end if

    call conservative2primitive(nv, np, ngz, u_int2)
    
    call compute_RHS(ngz, u_int2, one_grid)
    
    u_int3 = u_p + dt * rhs
    
    call boundary_conditions(nv, np, ngz, u_int3)
    
    if (associated(one_grid%parent)) then
      tplus = 1.0_wp
      call prolongate_boundaries(one_grid, u_int3, tplus)
    end if

    call conservative2primitive(nv, np, ngz, u_int3)

    call compute_RHS(ngz, u_int3, one_grid)

    u = (-2.0_wp * u_p + &
          2.0_wp * u_int1 + 4.0_wp * u_int2 + 2.0_wp * u_int3 + &
          dt * rhs) / 6.0_wp

  end if
    
  call boundary_conditions(nv, np, ngz, u)

  if (associated(one_grid%parent)) then
    tplus = 1.0_wp
    call prolongate_boundaries(one_grid, u, tplus)
  end if

  call conservative2primitive(nv, np, ngz, u)

  one_grid%iteration = one_grid%iteration + one_grid%iteration_step

end subroutine evolve_one_step

subroutine set_dt(one_grid)

  use real_type_mod
  use method_parameters_mod
  use grid_mod
  use rhs_mod

  implicit none

  type(grid) :: one_grid

  integer :: nv, np, ngz
  real(kind=wp) :: dt, dx, char_max

  nv  = one_grid%nvars
  np  = one_grid%n_points
  ngz = one_grid%ghostzones
  dx  = one_grid%dx

!!$  Set the timestep if varying it

  if (variable_dt) then

    call maxval_char(nv, np, ngz, one_grid%u, char_max)

    one_grid%dt = courant * dx / char_max

  end if

end subroutine set_dt
