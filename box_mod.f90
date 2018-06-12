!!$This defines the shape of a grid or an array or...

module box_mod

!!$  The upper and lower bounds give the location of this box
!!$  with respect to a parent grid
!!$  The actual local numbering runs [1 - nghosts : npoints + nghosts]
!!$  The bbox uses Cactus syntax.
!!$  That is, 1 means a "physical" boundary and 0 a code (processor or
!!$  MR boundary
!!$  Of course, "physical" is with respect to the current map / domain

  type box

     integer :: lbnd, ubnd
     integer :: npoints
     integer :: nghosts
     integer, dimension(2) :: bbox

  end type box
  
contains

  subroutine setup_base_box(mbox, npoints, nghosts)

    implicit none

    type(box), intent(out) :: mbox
    integer, intent(in) :: npoints, nghosts

!!$    There is no parent

    mbox%lbnd    = -1
    mbox%ubnd    = -1

    mbox%npoints = npoints
    mbox%nghosts = nghosts

!!$    All boundaries are physical

    mbox%bbox    = 1

  end subroutine setup_base_box
  
  subroutine setup_box(mbox, zero_box)

    implicit none

    type(box), intent(out) :: mbox
    type(box), intent(in) :: zero_box

    mbox = zero_box

  end subroutine setup_box
  
end module box_mod
