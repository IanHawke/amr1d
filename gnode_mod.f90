!!$This module is a "node" for a single grid.

module gnode_mod

  use grid_mod
  use boxnode_mod

  type gnode

     type(grid),  pointer :: this_grid
     type(gnode), pointer :: next_gnode

  end type gnode

contains

  subroutine setup_base_gnode(node)

    implicit none

    type(gnode), intent(out) :: node

!!$    if (associated(node%this_grid)) then
!!$      write(*,*) 'This node already has the grid allocated!'
!!$      STOP
!!$    end if
    
    allocate(node%this_grid)
    call setup_base_grid(node%this_grid)
    nullify(node%next_gnode)

  end subroutine setup_base_gnode

  subroutine setup_gnode(node, bnode, nvars, dx, dt, parent)

    implicit none

    type(gnode), intent(out)    :: node
    type(boxnode), intent(in)   :: bnode
    integer, intent(in)         :: nvars
    real(kind=wp), intent(in)   :: dx, dt
    type(gnode), pointer        :: parent

    if (associated(node%this_grid)) then
      write(*,*) 'This node already has the grid allocated!'
      STOP
    end if
    
    allocate(node%this_grid)
    call setup_grid(node%this_grid, bnode%this_box, nvars, dx, dt, &
         parent%this_grid)
    nullify(node%next_gnode)

  end subroutine setup_gnode

  subroutine delete_gnode(node)

    implicit none

    type(gnode) :: node

    call delete_grid(node%this_grid)

  end subroutine delete_gnode

end module gnode_mod
