!!$A list of grids

module grid_list_mod

  use generic_list
  use grid_mod

  type grid_node
     type(link_type) :: link
     type(grid), pointer :: mygrid
  end type grid_node

  type grid_node_ptr
     type(grid_node), pointer :: node_ptr
  end type grid_node_ptr
  
contains

  subroutine setup_base_grid_list(grid_list)

    implicit none

    type(list_type), intent(out) :: grid_list

    type(Link_Ptr_Type) :: link
    type(grid_node_ptr) :: base_grid_node_ptr

    allocate(base_grid_node_ptr%node_ptr)
    allocate(base_grid_node_ptr%node_ptr%mygrid)

    call setup_base_grid(base_grid_node_ptr%node_ptr%mygrid)

    call LI_Init_List(grid_list)
  
    link = TRANSFER(base_grid_node_ptr, link)
    call LI_Add_To_Head(link, grid_list)

  end subroutine setup_base_grid_list
  
  subroutine setup_grid_list(grid_list)

    implicit none

    type(list_type), intent(out) :: grid_list

    call LI_Init_List(grid_list)
  
  end subroutine setup_grid_list
  
  subroutine add_grid_to_list(grid_list, add_grid)

    implicit none

    type(list_type), intent(out) :: grid_list
    type(grid), pointer :: add_grid

    type(Link_Ptr_Type) :: link
    type(grid_node_ptr) :: add_grid_node_ptr

    allocate(add_grid_node_ptr%node_ptr)
    add_grid_node_ptr%node_ptr%mygrid => add_grid

    link = TRANSFER(add_grid_node_ptr, link)
    call LI_Add_To_Tail(link, grid_list)

  end subroutine add_grid_to_list

  subroutine delete_grid_list(grid_list)

    implicit none

    type(list_type) :: grid_list

    type(Link_Ptr_Type) :: link
    type(grid_node_ptr), pointer :: delete_node_ptr

    link = LI_Remove_Head(grid_list)
    
    allocate(delete_node_ptr)

    do while (LI_associated(link))

      delete_node_ptr = TRANSFER(link, delete_node_ptr)
      call delete_grid(delete_node_ptr%node_ptr%mygrid)
      deallocate(delete_node_ptr%node_ptr%mygrid)
      deallocate(delete_node_ptr%node_ptr)
      link = LI_Remove_Head(grid_list)

    end do
        
    deallocate(delete_node_ptr)

  end subroutine delete_grid_list

  subroutine initialize_base_grid_list(grid_list)

    implicit none

    type(list_type) :: grid_list

    type(Link_Ptr_Type) :: link
    type(grid_node_ptr), pointer :: node_ptr

    link = LI_Get_Head(grid_list)
    allocate(node_ptr)
    node_ptr = TRANSFER(link, node_ptr)

    call initialize_grid(node_ptr%node_ptr%mygrid)
    deallocate(node_ptr)

  end subroutine initialize_base_grid_list

  subroutine get_head_grid_list(grid_list, base_grid)

    implicit none

    type(list_type) :: grid_list
    type(grid), pointer :: base_grid

    type(Link_Ptr_Type) :: link
    type(grid_node_ptr), pointer :: base_node_ptr

    link = LI_Get_Head(grid_list)

    allocate(base_node_ptr)

    base_node_ptr = TRANSFER(link, base_node_ptr)

    if (associated(base_node_ptr%node_ptr)) then
      base_grid => base_node_ptr%node_ptr%mygrid
    else
      base_grid => NULL()
    end if
    
    deallocate(base_node_ptr)

  end subroutine get_head_grid_list
  
  subroutine initialize_grid_list(grid_list, boxes, parent)

    use box_list_mod

    implicit none

    type(list_type) :: grid_list
    type(list_type)  :: boxes
    type(grid), pointer :: parent

    type(Link_Ptr_Type) :: link
    type(box_node_ptr),  pointer :: boxnode_ptr

    type(grid), pointer :: new_grid

    link = LI_Get_Head(boxes)
    allocate(boxnode_ptr)

    do while (LI_Associated(link))

      boxnode_ptr = TRANSFER(link, boxnode_ptr)
      allocate(new_grid)
      call setup_grid(new_grid, boxnode_ptr%node_ptr%mybox, parent)
      call prolongate_grid(new_grid)
      call add_grid_to_list(grid_list, new_grid)
      nullify(new_grid)

      link = LI_Get_Next(link)

    end do

    deallocate(boxnode_ptr)
    
  end subroutine initialize_grid_list

  subroutine swap_timelevels_grid_list(grid_list)

    implicit none

    type(list_type) :: grid_list

    type(Link_Ptr_Type) :: link
    type(grid_node_ptr), pointer :: node_ptr

    link = LI_Get_Head(grid_list)

    allocate(node_ptr)

    do while (LI_associated(link))

      node_ptr = TRANSFER(link, node_ptr)
      call swap_timelevels(node_ptr%node_ptr%mygrid)
      link = LI_Get_Next(link)

    end do

    deallocate(node_ptr)

  end subroutine swap_timelevels_grid_list
  
  subroutine evolve_grid_list(grid_list)

    implicit none

    type(list_type) :: grid_list

    type(Link_Ptr_Type) :: link
    type(grid_node_ptr), pointer :: node_ptr

    link = LI_Get_Head(grid_list)

    allocate(node_ptr)

    do while (LI_Associated(link))

      node_ptr = TRANSFER(link, node_ptr)
      call evolve_one_step(node_ptr%node_ptr%mygrid)
      link = LI_Get_Next(link)

    end do

    deallocate(node_ptr)

  end subroutine evolve_grid_list

  subroutine restrict_grid_list(grid_list)

    implicit none

    type(list_type) :: grid_list

    type(Link_Ptr_Type) :: link
    type(grid_node_ptr), pointer :: node_ptr

    link = LI_Get_Head(grid_list)

    allocate(node_ptr)

    do while (LI_Associated(link))

      node_ptr = TRANSFER(link, node_ptr)
      call restrict_grid(node_ptr%node_ptr%mygrid)
      link = LI_Get_Next(link)

    end do

    deallocate(node_ptr)

  end subroutine restrict_grid_list
    
  subroutine output_grid_list(grid_list, depth, t)

    implicit none

    type(list_type) :: grid_list
    integer, intent(in) :: depth
    real(kind=wp), intent(in) :: t

    type(Link_Ptr_Type) :: link
    type(grid_node_ptr), pointer :: node_ptr
    integer :: ngrid

    link = LI_Get_Head(grid_list)
    ngrid = 0

    allocate(node_ptr)

    do while (LI_Associated(link))

      ngrid = ngrid + 1
      node_ptr = TRANSFER(link, node_ptr)
      call output_grid(node_ptr%node_ptr%mygrid, depth, ngrid, t)
      link = LI_Get_Next(link)

    end do

    deallocate(node_ptr)

  end subroutine output_grid_list

  subroutine regrid_grid_list(grid_list, child_list)

    implicit none

    type(list_type) :: grid_list
    type(list_type), pointer :: child_list

    type(Link_Ptr_Type) :: link
    type(grid_node_ptr), pointer :: node_ptr

    type(list_type), pointer :: grid_child_list

    link = LI_Get_Head(grid_list)

    allocate(node_ptr)

    do while (LI_Associated(link))

      allocate(grid_child_list)
      call LI_Init_List(grid_child_list)

      node_ptr = TRANSFER(link, node_ptr)
      call regrid_one_grid(node_ptr%node_ptr%mygrid, grid_child_list)
      link = LI_Get_Next(link)

      call add_grid_list_to_grid_list(grid_child_list, child_list)

      grid_child_list => NULL()

    end do

    deallocate(node_ptr)

  end subroutine regrid_grid_list

  subroutine regrid_one_grid(base_grid, new_child_list)

    use box_list_mod

    implicit none

    type(grid), pointer :: base_grid
    type(list_type), pointer :: new_child_list

    type(list_type), pointer :: new_box_list

    type(Link_Ptr_Type) :: link
    type(box_node_ptr), pointer :: node_ptr

    type(grid), pointer :: new_grid

    allocate(new_box_list)
    call LI_Init_List(new_box_list)

    call find_error_boxes(base_grid, new_box_list)

!!$    We now have to set up each of the child grids and put them in
!!$    the list

    allocate(node_ptr)
    link = LI_Remove_Head(new_box_list)
    do while (LI_Associated(link))

      node_ptr = TRANSFER(link, node_ptr)
      allocate(new_grid)
      call setup_grid(new_grid, node_ptr%node_ptr%mybox, base_grid)
      call prolongate_grid(new_grid)
      call add_grid_to_list(new_child_list, new_grid)
      new_grid => NULL()
      deallocate(node_ptr%node_ptr%mybox)
      deallocate(node_ptr%node_ptr)
      link = LI_Remove_Head(new_box_list)

    end do
    
  end subroutine regrid_one_grid
  
  subroutine add_grid_list_to_grid_list(new_list, base_list)

    implicit none

    type(list_type) :: new_list, base_list

    type(Link_Ptr_Type) :: link
    type(grid_node_ptr), pointer :: node_ptr

    link = LI_Get_Head(new_list)

    allocate(node_ptr)

    do while (LI_Associated(link))

      node_ptr = TRANSFER(link, node_ptr)
      call add_grid_to_list(base_list, node_ptr%node_ptr%mygrid)
      link = LI_Get_Next(link)

    end do

    deallocate(node_ptr)

  end subroutine add_grid_list_to_grid_list

  subroutine reconnect_child_list(parent_list, child_list)

    implicit none

    type(list_type) :: parent_list, child_list

    type(Link_Ptr_Type) :: parent_link, child_link
    type(grid_node_ptr), pointer :: parent_node_ptr, child_node_ptr

    allocate(parent_node_ptr, child_node_ptr)

    child_link = LI_Get_Head(child_list)

    do while (LI_Associated(child_link))

      child_node_ptr = TRANSFER(child_link, child_node_ptr)
      
      parent_link = LI_Get_Head(parent_list)
      do while (LI_Associated(parent_link))
        parent_node_ptr = TRANSFER(parent_link, parent_node_ptr)
        if (grid_child_of_parent(child_node_ptr%node_ptr%mygrid, &
                                 parent_node_ptr%node_ptr%mygrid)) exit
        parent_link = LI_Get_Next(parent_link)
      end do
      if (.not.LI_Associated(parent_link)) then
        write(*,*) 'Failed to reconnect a child grid!'
        STOP
      end if
      child_node_ptr%node_ptr%mygrid%parent => parent_node_ptr%node_ptr%mygrid
      call correct_box(child_node_ptr%node_ptr%mygrid)

      child_link = LI_Get_Next(child_link)

    end do

    deallocate(parent_node_ptr, child_node_ptr)

  end subroutine reconnect_child_list

  subroutine add_children_to_errors(child_list)

    implicit none

    type(list_type) :: child_list

    type(Link_Ptr_Type) :: child_link
    type(grid_node_ptr), pointer :: child_node_ptr

    type(grid), pointer :: parent, grandparent

    integer :: parent_i_start, parent_i_end, &
         grandparent_i_start, grandparent_i_end, &
         grid_n_points, parent_n_points

    allocate(child_node_ptr)

    child_link = LI_Get_Head(child_list)

    do while (LI_Associated(child_link))

      child_node_ptr = TRANSFER(child_link, child_node_ptr)

      parent => child_node_ptr%node_ptr%mygrid%parent
      grandparent => child_node_ptr%node_ptr%mygrid%parent%parent
      if (.not.associated(grandparent)) then
        write(*,*) 'Grandparent not associated. Should not be here!!!'
        STOP
      end if
      
      parent_i_start  = child_node_ptr%node_ptr%mygrid%mybox%lbnd
      grid_n_points   = child_node_ptr%node_ptr%mygrid%n_points / 2 
      parent_i_end    = parent_i_start + grid_n_points - 1

      grandparent_i_start  = parent%mybox%lbnd + parent_i_start / 2
      parent_n_points      = grid_n_points / 2 + mod(grid_n_points, 2)
      grandparent_i_end    = grandparent_i_start + parent_n_points - 1 
      
      grandparent%local_error(:, grandparent_i_start:grandparent_i_end) = &
           1.0_wp + 10.0_wp * amr_error

      child_link = LI_Get_Next(child_link)

    end do
    
  end subroutine add_children_to_errors

  subroutine output_grid_structure(grid_list, depth)

    implicit none

    type(list_type), pointer :: grid_list
    integer, intent(in) :: depth

    type(Link_Ptr_Type) :: link
    type(grid_node_ptr), pointer :: node_ptr

    allocate(node_ptr)

    link = LI_Get_Head(grid_list)

    write(*,*) 'Grid structure for level', depth

    do while (LI_Associated(link))

      node_ptr = TRANSFER(link, node_ptr)

      write(*,*) 'Grid from', node_ptr%node_ptr%mygrid%x0, 'to', &
           node_ptr%node_ptr%mygrid%x0 + &
           node_ptr%node_ptr%mygrid%n_points * node_ptr%node_ptr%mygrid%dx

      link = LI_Get_Next(link)

    end do
    
    deallocate(node_ptr)

  end subroutine output_grid_structure
  
  subroutine set_initial_error_on_grid_list(grid_list)

    implicit none

    type(list_type) :: grid_list

    type(Link_Ptr_Type) :: link
    type(grid_node_ptr), pointer :: node_ptr

    allocate(node_ptr)

    link = LI_Get_Head(grid_list)

    do while (LI_Associated(link))

      node_ptr = TRANSFER(link, node_ptr)
      call set_initial_error_on_grid(node_ptr%node_ptr%mygrid)

      link = LI_Get_Next(link)

    end do

    deallocate(node_ptr)
    
  end subroutine set_initial_error_on_grid_list
  
  subroutine apply_physical_boundaries_grid_list(grid_list)

    implicit none

    type(list_type) :: grid_list

    type(Link_Ptr_Type) :: link
    type(grid_node_ptr), pointer :: node_ptr

    allocate(node_ptr)

    link = LI_Get_Head(grid_list)

    do while (LI_Associated(link))

      node_ptr = TRANSFER(link, node_ptr)
      call apply_physical_boundaries_grid(node_ptr%node_ptr%mygrid)

      link = LI_Get_Next(link)

    end do

    deallocate(node_ptr)
    
  end subroutine apply_physical_boundaries_grid_list
  
  subroutine copy_data_from_lists(new_list, old_list)

    implicit none

    type(list_type), pointer :: new_list, old_list

    type(Link_Ptr_Type) :: new_link, old_link
    type(grid_node_ptr), pointer :: new_node_ptr, old_node_ptr

    type(grid), pointer :: new_grid, old_grid

    integer :: istart, iend, offset

    allocate(new_node_ptr, old_node_ptr)

    new_link = LI_Get_Head(new_list)

    do while (LI_Associated(new_link))

      new_node_ptr = TRANSFER(new_link, new_node_ptr)
      new_grid => new_node_ptr%node_ptr%mygrid

      old_link = LI_Get_Head(old_list)

      do while (LI_Associated(old_link))

        istart = -10
        iend   = -10
        offset = 0

        old_node_ptr = TRANSFER(old_link, old_node_ptr)
        old_grid => old_node_ptr%node_ptr%mygrid

        if (grid_intersection(new_grid, &
                              old_grid, &
                              istart, iend, offset)) then
          new_grid%u(:, istart:iend) = old_grid%u(:, istart+offset:iend+offset)
        end if
        
        old_link = LI_Get_Next(old_link)

      end do

      new_link = LI_Get_Next(new_link)

    end do
    
    deallocate(new_node_ptr, old_node_ptr)

  end subroutine copy_data_from_lists

  subroutine set_dt_from_parent_grid_list(grid_list)

    implicit none

    type(list_type), pointer :: grid_list

    type(Link_Ptr_Type) :: link
    type(grid_node_ptr), pointer :: node_ptr

    allocate(node_ptr)

    link = LI_Get_Head(grid_list)

    do while (LI_Associated(link))

      node_ptr = TRANSFER(link, node_ptr)
      call set_dt_from_parent_grid(node_ptr%node_ptr%mygrid)

      link = LI_Get_Next(link)

    end do

    deallocate(node_ptr)

  end subroutine set_dt_from_parent_grid_list

  subroutine find_max_char_grid_list(grid_list, max_char)

    implicit none

    type(list_type), pointer :: grid_list
    real(kind=wp), intent(inout) :: max_char

    type(Link_Ptr_Type) :: link
    type(grid_node_ptr), pointer :: node_ptr

    allocate(node_ptr)

    link = LI_Get_Head(grid_list)

    do while (LI_Associated(link))

      node_ptr = TRANSFER(link, node_ptr)
      call find_max_char_grid(node_ptr%node_ptr%mygrid, max_char)

      link = LI_Get_Next(link)

    end do

    deallocate(node_ptr)

  end subroutine find_max_char_grid_list
  
end module grid_list_mod
