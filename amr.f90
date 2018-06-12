!!$Driver program for the slope limited TVD solver using AMR

program amr

  use real_type_mod
  use grid_parameters_mod
  use level_mod
    
  implicit none

!!$  The base level

  type(level), pointer :: base, child, finest

!!$  The current time

  real(kind=wp) :: t

!!$  Indices

  integer :: i, k, n

  integer :: iteration

!!$  Error check

  integer :: ierr

!!$  Child size and shape

  type(box), pointer  :: child_box1, child_box2
  type(list_type) :: child_box_list

!!$  The base grid is temporarily required to initialize the child

  type(grid), pointer :: base_grid

  real(kind=wp) :: max_char

!!$  Read parameters

  call read_parameters

!!$  Check

  if (ghostzones < 2) then
    write(*,*) 'TVD requires 2 ghostzones!'
    stop
  end if

!!$  Now we actually do things

  t = 0.0_wp
  
  allocate(base)
  
  call setup_base_level(base, base_dx, base_dt)

  call initialize_base_level(base)
  
!!$    Initialize child

  allocate(child)

  call setup_box_list(child_box_list)

  allocate(child_box1, child_box2)

!!$  The box should cover the entire base because of the shadow heirarchy

  child_box1%lbnd = 1
  child_box1%ubnd = base_n_points
  child_box1%npoints = 2 * base_n_points
  child_box1%nghosts = ghostzones
  child_box1%bbox = 1
  call add_box_to_list(child_box_list, child_box1)

  call setup_level(child, base)
  call get_base_grid(base, base_grid)
  call initialize_grid_list(child%grid_list, child_box_list, base_grid)

!!$  Now look for a sensible initial grid setup

  finest => base%child
  do while (associated(finest))
    call set_initial_error_on_level(finest)
    call regrid_level(finest)
    finest => finest%child
  end do
  
!!$  Having initialized, run until the end time 

  iteration = 0

!!$    Output

  call output_all_levels(base, iteration)

  do while ( last_time - t > 0.000001_wp * base_dt)

!!$    Set timestep

!!$    call set_dt(base_grid)
    max_char = 0.0_wp
    call find_max_char_level(base, max_char)
    base_grid%dt = courant * base_grid%dx / max_char
!!$    This hack does not help the CW case
!!$    It seems that the oscillations occur at collision
    if (iteration < 3) then
      base_grid%dt = base_grid%dt / 10.0_wp
    end if
    base%dt = base_grid%dt
    call set_dt_from_parent(base%child)

!!$    Then evolve
      
    t = t + base%dt

!!$      Single evolution step
    
    call evolve_level(base)

    iteration = iteration + 1

!!$    Output

    call output_all_levels(base, iteration)
    
  end do

  write(*,*) 'Did ', iteration, ' iterations.'

  call output_level(base, t)
      
!!$    call delete_grid(child)
!!$    call delete_grid(base)

end program amr
