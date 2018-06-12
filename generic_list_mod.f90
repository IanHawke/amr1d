!!$A slight modification of the generic_list module of Young and McGavin
!!$to include doubly linked lists
!
!
! The problem is how to treat generic lists in Fortran 90. 
! Lists can be (1) homogeneous (elements all of the same type) or 
! (2) heterogeneous.
! 
! (1) Often there is a need to work with many (homogeneous) lists which
! may be of different types. Patrice Lignelet has shown how generic list
! properties may be logically treated in Fortran 90. However it still
! appears that the list operations (initialization, addition/removal of
! elements etc) must be separately defined for each list type, which
! leads to considerable duplication of code.
! 
! (2) Jean Vezina has shown how to handle a heterogeneous list by employing the 
! F90 TRANSFER() function. 
! 
! Peter McGavin at Industrial Research Limited (p.mcgavin@irl.cri.nz)
! has constructed a species of generic list for Fortran 90. The method
! is based on 2 ideas: the properties of the TRANSFER() function, and
! the fact that a pointer to a derived data type also points to the
! *first field* within the data type (and conversely). 
! 
! Since the TRANSFER() function does not accept pointer arguments the
! method requires the introduction of 2 auxilliary data types which
! contain the pointers, one in the generic list module and one for each
! list type in the calling program (the same construction is adopted
! when defining "arrays of pointers"). In spite of this complication the
! method represents a big saving, both conceptually and practically, when
! many lists of different types are involved.
! 
! To make clear the method we present a simple generic list module
! together with a calling program. The list module defines a
! uni-directional linked list with a few sample operations, but
! obviously more complicated generic lists could be substituted in its
! place (eg include back pointers).
! 
! Roger Young
! Peter McGavin
! 
! .........................................................................

MODULE Generic_List
! Defines a generic uni-directional linked list with a few sample operations

IMPLICIT NONE

PRIVATE

PUBLIC :: &
     Link_Type,        &! Put a Link_Type field first in your structure
     Link_Ptr_Type,    &! Mold this to and from your type ptr with TRANSFER
     List_Type          ! You should declare a List_Type variable

PUBLIC :: & 
     LI_Init_List,        &! Initialise the List_Type variable before use
     LI_Get_Head,         &! Returns the first Link in the list
     LI_Get_Tail,         &! Returns the last Link in the list
     LI_Get_Next,         &! Return the next Link after a given one
     LI_Get_Prev,         &! Return the previous Link before a given one
     LI_Add_To_Head,      &! Add a Link to the head of the list
     LI_Add_To_Tail,      &! Add a Link to the tail of the list
     LI_Remove_Head,      &! Remove the first Link and return it
     LI_Remove_Tail,      &! Remove the last Link and return it
     LI_Get_Len,          &! Compute list length
     LI_Associated,       &! Check if list member is associated
     LI_Check_List         ! Aborts program if list is invalid or corrupt

TYPE Link_Type
  PRIVATE
  TYPE(Link_Type), POINTER :: Next
  TYPE(Link_Type), POINTER :: Prev
END TYPE Link_Type

! Auxilliary data type required for the transfer function 
TYPE Link_Ptr_Type       ! Use TRANSFER() function to mold Link_Ptr_Type
  PRIVATE                ! to your pointer type and vice versa
  TYPE(Link_Type), POINTER :: P
END TYPE Link_Ptr_Type

TYPE List_Type             
  PRIVATE
  TYPE(Link_Type) :: Head   ! Dummy Link always at head of list
  TYPE(Link_Type) :: Tail   ! Dummy Link always at tail of list
END TYPE List_Type

CONTAINS

!-----------------------------------------------------------------------
SUBROUTINE Abort(Message)
IMPLICIT NONE
CHARACTER *(*) Message

WRITE(6,*) Message
WRITE(6,*) 'Program aborted'
STOP

END SUBROUTINE Abort

!-----------------------------------------------------------------------
SUBROUTINE LI_Check_List(List,Message)
IMPLICIT NONE
TYPE(List_Type) List
CHARACTER *(*) Message

IF(.NOT.ASSOCIATED(List%Head%Next))THEN
   WRITE(6,*) Message
   CALL Abort('List is not initialised in call to LI_Check_List()')
ENDIF

END SUBROUTINE LI_Check_List

!-----------------------------------------------------------------------
SUBROUTINE LI_Init_List(List)
  IMPLICIT NONE
  TYPE(List_Type),INTENT(INOUT),TARGET :: List

  NULLIFY(List%Head%Next)
  NULLIFY(List%Tail%Next)
  NULLIFY(List%Head%Prev)
  NULLIFY(List%Tail%Prev)

  RETURN
END SUBROUTINE LI_Init_List

!-----------------------------------------------------------------------
SUBROUTINE LI_Add_To_Head(Link,List)
  IMPLICIT NONE
  TYPE(List_Type),INTENT(INOUT)     :: List
  TYPE(Link_Ptr_Type),INTENT(INOUT) :: Link

  Link%P%Next => List%Head%Next
  List%Head%Next => Link%P

  NULLIFY(Link%P%Prev)

  if (ASSOCIATED(Link%P%Next)) then ! This is not the first node
    Link%P%Next%Prev => Link%P
  else                              ! This is the first node
    List%Tail%Prev => Link%P
  end if
  
  RETURN
END SUBROUTINE LI_Add_To_Head

!-----------------------------------------------------------------------
SUBROUTINE LI_Add_To_Tail(Link,List)
  IMPLICIT NONE
  TYPE(List_Type),INTENT(INOUT)     :: List
  TYPE(Link_Ptr_Type),INTENT(INOUT) :: Link

  Link%P%Prev => List%Tail%Prev
  List%Tail%Prev => Link%P

  NULLIFY(Link%P%Next)

  if (ASSOCIATED(Link%P%Prev)) then ! This is not the first node
    Link%P%Prev%Next => Link%P
  else                              ! This is the first node
    List%Head%Next => Link%P
  end if
  
  RETURN
END SUBROUTINE LI_Add_To_Tail

!-----------------------------------------------------------------------
INTEGER FUNCTION LI_Get_Len(List)
  IMPLICIT NONE
  TYPE(List_Type), INTENT(IN),TARGET :: LIST
  TYPE(Link_Ptr_Type) :: Link
  INTEGER N

  Link%P => List%Head
  N = 0
  DO WHILE(ASSOCIATED(Link%P%Next))
     Link%P => Link%P%Next
     N = N+1
  ENDDO
  LI_Get_Len = N

  RETURN
END FUNCTION LI_Get_Len

!-----------------------------------------------------------------------
FUNCTION LI_Associated(Link)
  IMPLICIT NONE
  LOGICAL :: LI_Associated
  TYPE(Link_Ptr_Type),INTENT(IN) :: Link

  LI_Associated = .FALSE.
  IF(ASSOCIATED(Link%P))LI_Associated=.TRUE.

  RETURN
END FUNCTION LI_Associated

!-----------------------------------------------------------------------
FUNCTION LI_Get_Next(Link)
  IMPLICIT NONE
  Type(Link_Ptr_Type)           :: LI_Get_Next
  TYPE(Link_Ptr_Type),INTENT(IN) :: Link

  IF(.NOT.ASSOCIATED(Link%P%Next))THEN
     NULLIFY(LI_Get_Next%P)
  ELSE   
     LI_Get_Next%P => Link%P%Next
  ENDIF

  RETURN
END FUNCTION LI_Get_Next

!-----------------------------------------------------------------------
FUNCTION LI_Get_Prev(Link)
  IMPLICIT NONE
  Type(Link_Ptr_Type)           :: LI_Get_Prev
  TYPE(Link_Ptr_Type),INTENT(IN) :: Link

  IF(.NOT.ASSOCIATED(Link%P%Prev))THEN
     NULLIFY(LI_Get_Prev%P)
  ELSE   
     LI_Get_Prev%P => Link%P%Prev
  ENDIF

  RETURN
END FUNCTION LI_Get_Prev

!-----------------------------------------------------------------------
FUNCTION LI_Get_Head(List)
  IMPLICIT NONE
  TYPE(Link_Ptr_Type)               :: LI_Get_Head
  TYPE(List_Type),INTENT(IN),TARGET :: List

  LI_Get_Head%P => List%Head%Next

  RETURN
END FUNCTION LI_Get_Head

!-----------------------------------------------------------------------
FUNCTION LI_Get_Tail(List)
  IMPLICIT NONE
  TYPE(Link_Ptr_Type)               :: LI_Get_Tail
  TYPE(List_Type),INTENT(IN),TARGET :: List

  LI_Get_Tail%P => List%Tail%Prev

  RETURN
END FUNCTION LI_Get_Tail

!-----------------------------------------------------------------------
FUNCTION LI_Remove_Head(List)
  IMPLICIT NONE
  TYPE(Link_Ptr_Type)                  :: LI_Remove_Head
  TYPE(List_Type),INTENT(INOUT),TARGET :: List
  TYPE(Link_Ptr_Type) :: Link

  Link%P => List%Head%Next
  IF(ASSOCIATED(Link%P))THEN
     List%Head%Next => Link%P%Next
     if (associated(link%p, list%tail%prev)) then
       list%tail%prev => NULL()
     end if
     NULLIFY(Link%P%Next)
     NULLIFY(Link%P%Prev)
  ENDIF

  LI_Remove_Head%P => Link%P

  RETURN
END FUNCTION LI_Remove_Head

!-----------------------------------------------------------------------
FUNCTION LI_Remove_Tail(List)
  IMPLICIT NONE
  TYPE(Link_Ptr_Type)                  :: LI_Remove_Tail
  TYPE(List_Type),INTENT(INOUT),TARGET :: List
  TYPE(Link_Ptr_Type) :: Link

  Link%P => List%Tail%Prev
  IF(ASSOCIATED(Link%P))THEN
     List%Tail%Prev => Link%P%Prev
     if (associated(link%p, list%head%next)) then
       list%head%next => NULL()
     end if
     NULLIFY(Link%P%Next)
     NULLIFY(Link%P%Prev)
  ENDIF

  LI_Remove_Tail%P => Link%P

  RETURN
END FUNCTION LI_Remove_Tail

!-----------------------------------------------------------------------
END MODULE Generic_List
