!!$This is a linked list of gnodes

module glist_mod

  use gnode_mod
  use boxlist_mod

  type glist

     type(gnode), pointer :: head

  end type glist

contains

  subroutine setup_base_glist(list)

    implicit none

    type(glist), intent(out) :: list
    
    allocate(list%head)
    call setup_base_gnode(list%head)

  end subroutine setup_base_glist

  subroutine copy_glist(list, other)

    implicit none

    type(glist) :: list, other

    list = other

  end subroutine copy_glist
  
  subroutine setup_glist(list, boxes, nvars, dx, dt, parent)

    implicit none

    type(glist), intent(out) :: list
    type(boxlist), intent(in) :: boxes
    type(boxnode), pointer :: bnode
    type(gnode), pointer :: new_gnode
    integer, intent(in)       :: nvars
    real(kind=wp), intent(in) :: dx, dt
    type(gnode), pointer :: parent

    bnode => boxes%head
    do while (associated(bnode%next_boxnode))
      allocate(new_gnode)
      call setup_gnode(new_gnode, bnode, nvars, dx, dt, parent)
      call add_node_to_glist(list, new_gnode)
      new_gnode => NULL()
      bnode => bnode%next_boxnode
    end do
    
    if (associated(bnode)) then
      allocate(new_gnode)
      call setup_gnode(new_gnode, bnode, nvars, dx, dt, parent)
      call add_node_to_glist(list, new_gnode)
      new_gnode => NULL()
    end if
    
  end subroutine setup_glist
  
  subroutine add_node_to_glist(list, node)

    implicit none

    type(glist) :: list
    type(gnode), pointer :: node
    type(gnode), pointer :: current_node

    current_node => list%head

    do while (associated(current_node%next_gnode)) 
      current_node => current_node%next_gnode
    end do
    
    current_node%next_gnode => node
    nullify(node%next_gnode)

  end subroutine add_node_to_glist

  subroutine delete_glist(list)

    implicit none

    type(glist) :: list
    type(gnode), pointer :: current_node, next_node
  
    current_node => list%head
    next_node => list%head

    do while (associated(current_node%next_gnode))
      next_node => current_node%next_gnode
      call delete_gnode(current_node)
      current_node => next_node
    end do
    
    call delete_gnode(next_node)
    nullify(list%head)

  end subroutine delete_glist
  
end module glist_mod

    
    
