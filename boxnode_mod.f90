!!$This is a "node" for a single box.

module boxnode_mod

  use box_mod

  type boxnode

     type(box),     pointer :: this_box
     type(boxnode), pointer :: next_boxnode

  end type boxnode
  
contains

  subroutine setup_base_boxnode(node, npoints, nghosts)

    implicit none

    type(boxnode), intent(out) :: node
    integer :: npoints, nghosts

    allocate(node%this_box)
    call setup_base_box(node%this_box, npoints, nghosts)
    nullify(node%next_boxnode)

  end subroutine setup_base_boxnode
  
  subroutine setup_boxnode(node, zero_box)

    implicit none

    type(boxnode), intent(out) :: node
    type(box), intent(in)   :: zero_box

    if (associated(node%this_box)) then
      write(*,*) 'This node already has the box allocated!'
      STOP
    end if
    
    allocate(node%this_box)
    call setup_box(node%this_box, zero_box)
    nullify(node%next_boxnode)

  end subroutine setup_boxnode

end module boxnode_mod
