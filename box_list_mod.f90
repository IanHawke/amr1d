!!$A list of boxes

module box_list_mod

  use generic_list
  use box_mod

  type box_node
     type(link_type) :: link
     type(box), pointer :: mybox
  end type box_node

  type box_node_ptr
     type(box_node), pointer :: node_ptr
  end type box_node_ptr
  
contains

  subroutine setup_box_list(box_list)

    implicit none

    type(list_type), intent(out) :: box_list

    call LI_Init_List(box_list)
  
  end subroutine setup_box_list

  subroutine add_box_to_list(box_list, add_box)

    implicit none

    type(list_type), intent(out) :: box_list
    type(box), pointer :: add_box

    type(Link_Ptr_Type) :: link
    type(box_node_ptr), pointer :: add_box_node_ptr

    allocate(add_box_node_ptr)
    allocate(add_box_node_ptr%node_ptr)
    add_box_node_ptr%node_ptr%mybox => add_box

    link = TRANSFER(add_box_node_ptr, link)
    call LI_Add_To_Tail(link, box_list)

  end subroutine add_box_to_list

end module box_list_mod
