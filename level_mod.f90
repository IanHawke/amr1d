!!$This is a single refinement level

module level_mod

  use grid_list_mod
  use box_list_mod

  type level

     type(list_type), pointer :: grid_list

     type(level), pointer :: parent, child
     
     real(kind=wp) :: dx, dt

     integer :: depth

!!$iteration is the current iteration that this level is on.
!!$iteration_step is the number of iterations that a single evolution step
!!$will increase iteration by, at this refinement level.

     integer :: iteration
     integer :: iteration_step

  end type level
  
contains

  subroutine setup_base_level(base, dx, dt)

    implicit none

    type(level), intent(out) :: base
    real(kind=wp) :: dx, dt

    allocate(base%grid_list)
    call setup_base_grid_list(base%grid_list)

    base%parent => NULL()
    base%child  => NULL()

    base%dx = dx
    base%dt = dt

    base%depth = 1

    base%iteration      = 0
    base%iteration_step = 2**(max_refinement_levels - 1)

  end subroutine setup_base_level

  subroutine setup_level(base, parent)

    implicit none

    type(level), pointer :: base
    type(level), pointer :: parent

    allocate(base%grid_list)

    base%parent => parent
    base%child  => NULL()
    parent%child => base
    base%dx = parent%dx / 2.0_wp
    base%dt = parent%dt / 2.0_wp
    base%depth = parent%depth + 1

    base%iteration      = parent%iteration
    if (parent%iteration_step == 1) then
      write(*,*) 'Tried to allocate a level finer than max_refinement_levels allows!'
      STOP
    end if
    base%iteration_step = parent%iteration_step / 2

    call setup_grid_list(base%grid_list)
  
  end subroutine setup_level
  
  subroutine delete_level(base)

    implicit none

    type(level) :: base

    if (associated(base%child)) then
      write(*,*) 'You cannot delete a level with a child!'
      STOP
    end if
    
    call delete_grid_list(base%grid_list)
    deallocate(base%grid_list)
    nullify(base%grid_list)

  end subroutine delete_level

  subroutine initialize_base_level(base)

    implicit none

    type(level) :: base

    call initialize_base_grid_list(base%grid_list)
  
  end subroutine initialize_base_level
  
  subroutine get_base_grid(base, base_grid)

    implicit none

    type(level) :: base
    type(grid), pointer :: base_grid

    call get_head_grid_list(base%grid_list, base_grid)

  end subroutine get_base_grid

  recursive subroutine evolve_level(base)

    implicit none

    type(level), pointer :: base

    call swap_timelevels_grid_list(base%grid_list)
    call evolve_grid_list(base%grid_list)

    base%iteration = base%iteration + base%iteration_step

    if (associated(base%child)) then
      call evolve_level(base%child)
      call evolve_level(base%child)
    end if

    if (associated(base%parent)) then
      if (base%iteration == base%parent%iteration) then
        call restrict_grid_list(base%grid_list)
      end if
    end if

!!$    Physical and symmetry boundary conditions may need
!!$    reapplying after restriction
!!$    This requires more thought

!!$    call apply_physical_boundaries_level(base)

!!$    Check for regridding
!!$    The recursion means that the regrid_level call should naturally
!!$    unwind fine -> coarse without any recursion inside that routine

    if (associated(base%parent)) then
      if (base%iteration == base%parent%iteration) then
        call regrid_level(base)
      end if
    end if
    
  end subroutine evolve_level

  recursive subroutine output_level(base, t)

    implicit none

    type(level) :: base
    real(kind=wp), intent(in) :: t

    call output_grid_list(base%grid_list, base%depth, t)

    if (associated(base%child)) then
      call output_level(base%child, t)
    end if
    
  end subroutine output_level

  recursive subroutine regrid_level(base)

    implicit none

    type(level), pointer :: base

    type(level), pointer :: tmp_level

    type(list_type), pointer :: new_child_list
    
    if (base%depth < max_refinement_levels) then

      allocate(new_child_list)
      call LI_Init_List(new_child_list)

      if (associated(base%child)) then
        if (associated(base%child%child)) then
          call add_children_to_errors(base%child%child%grid_list)
        end if
      end if
      
      call regrid_grid_list(base%grid_list, new_child_list)

      if (associated(base%child) .or. &
          (LI_Get_Len(new_child_list) > 0) )  then
        if (.not.associated(base%child)) then
          allocate(base%child)
          call setup_level(base%child, base)
        end if

!!$        Copy over the data from the old grids

        call copy_data_from_lists(new_child_list, base%child%grid_list)

!!$        Get rid of the old list, install the new

        call delete_grid_list(base%child%grid_list)
        if (LI_Get_Len(new_child_list) > 0) then
          base%child%grid_list => new_child_list
          call reconnect_child_list(base%grid_list, base%child%grid_list)
          if (associated(base%child%child)) then
            call reconnect_child_list(base%child%grid_list, base%child%child%grid_list)
          end if
      
        else
          if (associated(base%child)) then
            if (associated(base%child%child)) then
              write(*,*) 'We are about to orphan a level!'
              STOP
            end if
          end if
          base%child => NULL()
        end if
      end if

    end if
    
  end subroutine regrid_level

  subroutine output_all_levels(base, iteration)

    type(level), pointer :: base
    integer, intent(in) :: iteration

    type(grid), pointer :: base_grid

    integer :: fileunit

    integer :: total_grid_points

    fileunit = 28

    total_grid_points = 0

    if (iteration == 0) then
      open(fileunit, file='amr_all.dat', status='replace')
    else
      open(fileunit, file='amr_all.dat', status='old', access='append')
    end if
    
    call get_head_grid_list(base%grid_list, base_grid)
    
    call output_all_grids(base, base_grid, fileunit, total_grid_points)
    write(fileunit, *)
    write(fileunit, *)

    close(fileunit)

    write(*,*) 'Outputting iteration', iteration, &
               ' Total grid points were', total_grid_points

  end subroutine output_all_levels

  recursive subroutine output_all_grids(base, base_grid, fileunit, &
                                        total_grid_points)

    type(level), pointer :: base
    type(grid), pointer  :: base_grid
    integer, intent(in) :: fileunit
    integer, intent(inout) :: total_grid_points

    type(Link_Ptr_Type) :: base_link, child_link
    type(grid_node_ptr), pointer :: base_node_ptr, child_node_ptr

    type(grid), pointer :: child_grid
    integer :: istart, iend

    allocate(base_node_ptr, child_node_ptr)

    base_link = LI_Get_Head(base%grid_list)

    do while (LI_Associated(base_link))
      
      base_node_ptr = TRANSFER(base_link, base_node_ptr)
      if (associated(base_grid, base_node_ptr%node_ptr%mygrid)) then

        istart = 1
        iend = base_grid%n_points

        total_grid_points = total_grid_points + base_grid%n_points

        if (associated(base%child)) then

          child_link = LI_Get_Head(base%child%grid_list)
          do while (LI_Associated(child_link))
            
            child_node_ptr = TRANSFER(child_link, child_node_ptr)
            child_grid => child_node_ptr%node_ptr%mygrid
            
            if (associated(child_grid%parent, base_grid)) then

              do i = istart, child_grid%mybox%lbnd - 1
                write(fileunit, '(es21.12E3,a1)',advance='NO') &
                     base_grid%x(i), ' '
                do k = 1, base_grid%nvars
                  write(fileunit, '(es21.12E3,a1)',advance='NO') &
                       base_grid%u(k, i), ' '
                end do
                write(fileunit, *)
              end do
              istart = child_grid%mybox%lbnd + child_grid%n_points / 2 
              call output_all_grids(base%child, child_grid, fileunit, &
                                    total_grid_points)
            end if
            
            child_link = LI_Get_Next(child_link)
          end do

        end if
        
        do i = istart, iend
          write(fileunit, '(es21.12E3,a1)',advance='NO') base_grid%x(i), ' '
          do k = 1, base_grid%nvars
            write(fileunit, '(es21.12E3,a1)',advance='NO') &
                 base_grid%u(k, i), ' '
          end do
          write(fileunit, *)
        end do
        
      end if
  
      base_link = LI_Get_Next(base_link)
      
    end do
    
    deallocate(base_node_ptr, child_node_ptr)

  end subroutine output_all_grids

  subroutine set_initial_error_on_level(base)

    implicit none

    type(level), pointer :: base
  
    call set_initial_error_on_grid_list(base%grid_list)
  
  end subroutine set_initial_error_on_level

  subroutine apply_physical_boundaries_level(base)

    implicit none

    type(level), pointer :: base

    call apply_physical_boundaries_grid_list(base%grid_list)
  
  end subroutine apply_physical_boundaries_level

  recursive subroutine set_dt_from_parent(base)

    implicit none

    type(level) :: base

    base%dt = base%parent%dt / 2.0_wp
    call set_dt_from_parent_grid_list(base%grid_list)

    if (associated(base%child)) then
      call set_dt_from_parent(base%child)
    end if
    
  end subroutine set_dt_from_parent

  recursive subroutine find_max_char_level(base, max_char)

    implicit none

    type(level) :: base
    real(kind=wp), intent(inout) :: max_char

    call find_max_char_grid_list(base%grid_list, max_char)

    if (associated(base%child)) then
      call find_max_char_level(base%child, max_char)
    end if
    
  end subroutine find_max_char_level
  
end module level_mod
