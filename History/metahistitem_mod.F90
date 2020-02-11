!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: metahistitem_mod.F90
!
! !DESCRIPTION: Contains types and methods to create a METAHISTORY ITEM
!  object, which is a container for a HISTORY ITEM.  In other words,
!  a METAHISTORY ITEM represents a single node of a linked list that is
!  used to contain HISTORY ITEMS.
!\\
!\\
!  In practice, we can think of a METAHISTORY ITEM as a list of HISTORY ITEMS
!  that will be archived to netCDF output at a specified frequency (e.g.
!  instantaneous, daily, hourly, etc.)
!\\
!\\
! !INTERFACE:
!
MODULE MetaHistItem_Mod
!
! !USES:
!
  USE HistItem_Mod, ONLY : HistItem
  USE Precision_Mod

  IMPLICIT NONE
  PRIVATE
!
! !PRIVATE MEMBER FUNCTIONS:
!
  PRIVATE :: MetaHistItem_Create
  PRIVATE :: MetaHistItem_Insert
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC  :: MetaHistItem_AddNew
  PUBLIC  :: MetaHistItem_Count
  PUBLIC  :: MetaHistItem_Destroy
  PUBLIC  :: MetaHistItem_Print
!
! !PUBLIC TYPES:
!
  !=========================================================================
  ! This is the derived type for a METAHISTORY ITEM object, which represents
  ! a SINGLE NODE OF A LINKED LIST consisting of HISTORY ITEMS.
  !
  ! As such, the METAHISTORY ITEM does not contain any data itself,
  ! but is a wrapper for a single HISTORY ITEM object, plus a pointer
  ! to another METAHISTORY ITEM (i.e. the next node in the list).
  !=========================================================================
  TYPE, PUBLIC ::  MetaHistItem

     ! Pointer to the next METAHISTORY ITEM object
     ! (i.e. the next node in the linked list)
     TYPE(MetaHistItem), POINTER :: Next => NULL()

     ! The HISTORY ITEM object (which represents a diagnostic
     ! quantity that will be archived to netCDF file format)
     TYPE(HistItem),     POINTER :: Item => NULL()

  END TYPE MetaHistItem
!
! !REMARKS:
!  As described above, a METAHISTORY ITEM can be thought of as a SINGLE NODE
!  OF A LINKED LIST INTENDED TO HOLD HISTORY ITEMS.  It looks like this:
!
!      +-------------------------+   +-------------------------+
!      | METAHISTORY ITEM n      |   | METAHISTORY ITEM n+1    |
!      | (aka NODE n of list)    |   | (aka NODE n+1 of list)  |
!      |                         |   |                         |
!      | Contains:               |   | Contains:               |
!      |                         |   |                         |
!      |   HISTORY ITEM n        |   |   HISTORY ITEM n+1      |
!      |                         |   |                         |
! =======> Pointer to next    =========> Pointer to next    ========> etc ...
!      |    METAHISTORY ITEM     |   |    METAHISTORY ITEM     |
!      +-------------------------+   +-------------------------+
!
!  Linked list routines taken from original code (linkedlist.f90)
!  by Arjen Markus; http://flibs.sourceforge.net/linked_list.html

! !REVISION HISTORY:
!  14 Jun 2017 - R. Yantosca - Initial version, based on code by Arjen Markus
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: MetaHistItem_AddNew
!
! !DESCRIPTION: Wrapper for methods MetaHistItem\_Create and
!  MetaHistItem\_Insert.  Will create a METAHISTORY ITEM (containing a
!  HISTORY ITEM) and (1) set it as the head node of a new linked list, or
!  (2) append it to an existing linked list.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE MetaHistItem_AddNew( Input_Opt, Node, Item, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE HistItem_Mod,  ONLY : HistItem
    USE Input_Opt_Mod, ONLY : OptInput
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput),     INTENT(IN)  :: Input_Opt  ! Input Options object
    TYPE(HistItem),     POINTER     :: Item       ! HISTORY ITEM object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(MetaHistItem), POINTER     :: Node       ! METAHISTORY ITEM object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,            INTENT(OUT) :: RC         ! Success or failure
!
! !REVISION HISTORY:
!  13 Jun 2017 - R. Yantosca - Initial version, based on code by Arjen Markus
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Strings
    CHARACTER(LEN=255) :: ErrMsg, ThisLoc

    !=======================================================================
    ! Initialize
    !=======================================================================

    ! Assume success
    RC     = GC_SUCCESS

    ! For error output
    ErrMsg  = ''
    ThisLoc = ' -> at MetaHistItem_Add (in History/metahistitem_mod.F90)'

    !=======================================================================
    ! Test if the METAHISTORY ITEM (aka "Node") has been allocated memory
    ! and is therefore part of an existing linked list
    !=======================================================================
    IF ( .not. ASSOCIATED( Node ) ) THEN

       !--------------------------------------------------------------------
       ! If not, then create a new METAHISTORY ITEM (named "Node"),
       ! and set it at the head of a new linked list
       !--------------------------------------------------------------------
       CALL MetaHistItem_Create( Input_Opt, Node, Item, RC )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Could not create "Node" as the head node of a list!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

    ELSE

       !--------------------------------------------------------------------
       ! Otherwise, create a new METAHISTORY ITEM (named "Node"),
       ! and append it to the list, immediately following the head node
       !--------------------------------------------------------------------
       CALL MetaHistItem_Insert( Input_Opt, Node, Item, RC )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Could not insert "Node" into an existing linked list!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

    ENDIF

  END SUBROUTINE MetaHistItem_AddNew
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: MetaHistItem_Create
!
! !DESCRIPTION: This method creates a new METAHISTORY ITEM (to contain the
!  supplied HISTORY ITEM) and sets it as the head node of a linked list.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE MetaHistItem_Create( Input_Opt, Node, Item, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE HistItem_Mod,  ONLY : HistItem
    USE Input_Opt_Mod, ONLY : OptInput
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput),     INTENT(IN)  :: Input_Opt  ! Input Options object
    TYPE(HistItem),     POINTER     :: Item       ! HISTORY ITEM object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(MetaHistItem), POINTER     :: Node       ! METAHISTORY ITEM object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,            INTENT(OUT) :: RC         ! Success or failure
!
! !REMARKS:
!  This method is not intended to be called directly, but is rather
!  wrapped by the MetaHistItem_AddNew method.
!
! !REVISION HISTORY:
!  13 Jun 2017 - R. Yantosca - Initial version, based on code by Arjen Markus
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Strings
    CHARACTER(LEN=255) :: ErrMsg, ThisLoc

    !=======================================================================
    ! Initialize
    !=======================================================================

    ! Assume success
    RC      = GC_SUCCESS

    ! For error output
    ErrMsg  = ''
    ThisLoc = ' -> at MetaHistItem_Create (in History/metahistitem_mod.F90)'

    !=======================================================================
    ! Initialize the METAHISTORY ITEM itself
    !=======================================================================

    ! Allocate memory
    ALLOCATE( Node, STAT=RC )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Could not allocate "Node"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    ! Nullify the pointer to the next node (it doesn't exist yet)
    Node%Next => NULL()

    !=======================================================================
    ! Initialize the field that will store the HISTORY ITEM
    !=======================================================================

    ! Because this is the first METAHISTORY ITEM that is being created,
    ! we can consider this to be the head node of a linked list.
    IF ( .not. ASSOCIATED( Node%Item ) ) THEN
       ALLOCATE( Node%Item, STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Could not allocate "Node%Item"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDIF

    ! Attach the given HISTORY ITEM to the METAHISTORY ITEM
    ! (i.e. place it into the head node of a linked list)
    Node%Item = Item

  END SUBROUTINE MetaHistItem_Create
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: MetaHistItem_Insert
!
! !DESCRIPTION: Creates a new METAHISTORY ITEM (to contain the supplied
!  HISTORY ITEM), and pops it into an existing linked list, immediately
!  following the head node.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE MetaHistItem_Insert( Input_Opt, Node, Item, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE HistItem_Mod,  ONLY : HistItem
    USE Input_Opt_Mod, ONLY : OptInput
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput),     INTENT(IN)  :: Input_Opt  ! Input Options object
    TYPE(HistItem),     POINTER     :: Item      ! HISTORY ITEM object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(MetaHistItem), POINTER     :: Node      ! METAHISTORY ITEM object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,            INTENT(OUT) :: RC        ! Success or failure
!
! !REMARKS:
!  This method is not intended to be called directly, but is rather
!  wrapped by the MetaHistItem_AddNew method.
!
! !REVISION HISTORY:
!  13 Jun 2017 - R. Yantosca - Initial version, based on code by Arjen Markus
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Strings
    CHARACTER(LEN=255)          :: ErrMsg, ThisLoc

    ! Objects
    TYPE(MetaHistItem), POINTER :: Head

    !=======================================================================
    ! Initialize
    !=======================================================================

    ! Assume success
    RC      = GC_SUCCESS

    ! For error output
    ErrMsg  = ''
    ThisLoc = ' -> at MetaHistItem_Insert (in History/metahistitem_mod.F90)'

    !=======================================================================
    ! Initialize a METAHISTORY ITEM named "Head", which will become the
    ! head of the existing list.  "Head" will contain a new HISTORY ITEM.
    !=======================================================================

    ! Allocate the "Head" object
    ALLOCATE( Head, STAT=RC )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Could not allocate "Head"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    ! Allocate the "Head%Item" field, which will hold the HISTORY ITEM
    ALLOCATE( Head%Item, STAT=RC )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Could not allocate "Head%Item"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !=======================================================================
    ! Insert "Head" at the start of the existing linked list
    !=======================================================================

    ! Save the HISTORY ITEM argument in the "Item" field of "Head"
    Head%Item  =  Item

    ! The "Next" field of "Head" points to the current head of the list
    Head%Next  => Node

    ! Set "Head" as the new head of the linked list
    Node       => Head

  END SUBROUTINE MetaHistItem_Insert
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: MetaHistItem_Count
!
! !DESCRIPTION: Counts the number of METAHISTORY ITEMS stored in a linked
!  list.  By extension, this is also the number of HISTORY ITEMS stored in
!  the list, because each METAHISTORY ITEM contains only one HISTORY ITEM.
!\\
!\\
! !INTERFACE:
!
  FUNCTION MetaHistItem_Count( List ) RESULT( nNodes )

!
! !INPUT PARAMETERS:
!
    TYPE(MetaHistItem), POINTER :: List    ! Linked list of METAHISTORY ITEMS
!
! !RETURN VALUE:
!
    INTEGER                     :: nNodes  ! Number of METAHISTORY ITEMS
!
! !REVISION HISTORY:
!  14 Jun 2017 - R. Yantosca - Initial version, based on code by Arjen Markus
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Objects
    TYPE(MetaHistItem), POINTER  :: Current

    !=======================================================================
    ! Initialize
    !=======================================================================

    ! Initialize
    nNodes  =  0
    Current => NULL()

    !=======================================================================
    ! Count the number of METAHISTORY ITEMS in the list
    !=======================================================================

    ! If the list does not point to NULL() ...
    IF ( ASSOCIATED( List ) ) THEN

       ! ... then there is at least 1 node (the head node)
       nNodes  =  1

       ! Set the CURRENT pointer to the head node
       Current => List

       ! As long as the following node doesn't point to NULL()
       DO WHILE ( ASSOCIATED( Current%Next ) )

          ! Set CURRENT to the following node
          Current => Current%Next

          ! Increment the node count
          nNodes  =  nNodes + 1

       ENDDO

       ! Free pointers
       Current => NULL()

    ENDIF

  END FUNCTION MetaHistItem_Count
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: MetaHistItem_Print
!
! !DESCRIPTION: This method will print information about the HISTORY ITEM
!  belonging to each METAHISTORY ITEM (aka node) of a linked list.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE MetaHistItem_Print( Input_Opt, List, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE HistItem_Mod,  ONLY : HistItem, HistItem_Print
    USE Input_Opt_Mod, ONLY : OptInput
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput),     INTENT(IN)  :: Input_Opt  ! Input Options object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(MetaHistItem), POINTER     :: List        ! List of history items
!
! !OUTPUT PARAMETERS:
!
    INTEGER,            INTENT(OUT) :: RC          ! Success or failure
!
! !REMARKS:
!
! !REVISION HISTORY:
!  14 Jun 2017 - R. Yantosca - Initial version, based on code by Arjen Markus
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Strings
    CHARACTER(LEN=255)          :: ErrMsg, ThisLoc

    ! Objects
    TYPE(MetaHistitem), POINTER :: Current

    !=======================================================================
    ! Initialize
    !=======================================================================

    ! Assume success
    RC      =  GC_SUCCESS

    ! Free pointers
    Current => NULL()

    ! For error output
    ErrMsg  =  ''
    ThisLoc =  ' -> at MetaHistItem_Print (in History/metahistitem_mod.F90)'

    !=======================================================================
    ! Print information about each METAHISTORY ITEM (aka node)
    ! of the linked list, only if we are on the root CPU
    !=======================================================================
    IF ( Input_Opt%amIRoot ) THEN

       ! Point CURRENT to the head node of the list
       Current => List

       ! As long as the current node is valid
       DO WHILE( ASSOCIATED( Current ) )

          ! Print info about the history item corresponding to this node
          CALL HistItem_Print( Input_Opt, Current%Item, RC )
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = 'Could not print info for "Current%Item"!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF

          ! Point to the next node for the next iteration
          Current => Current%Next
       ENDDO

       ! Free pointers
       Current => NULL()
    ENDIF

  END SUBROUTINE MetaHistItem_Print
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: MetaHistItem_Destroy
!
! !DESCRIPTION:  This method will destroy the HISTORY ITEM belonging to
!  each METAHISTORY ITEM (aka node) of a linked list.  It will then destroy
!  each METAHISTORY ITEM itself.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE MetaHistItem_Destroy( List, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE HistItem_Mod,  ONLY : HistItem, HistItem_Destroy
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(MetaHistItem), POINTER     :: List       ! List of METAHISTORY ITEMS
!
! !OUTPUT PARAMETERS:
!
    INTEGER,            INTENT(OUT) :: RC         ! Success or failure?
!
! !REMARKS:
!
! !REVISION HISTORY:
!  14 Jun 2017 - R. Yantosca - Initial version, based on code by Arjen Markus
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Strings
    CHARACTER(LEN=255)           :: ErrMsg,  ThisLoc

    ! Objects
    TYPE(MetaHistItem), POINTER  :: Current, Node

    !=======================================================================
    ! Initialize
    !=======================================================================

    ! Assume success
    RC      = GC_SUCCESS

    ! Free pointers
    Current => NULL()
    Node    => NULL()

    ! For error output
    ErrMsg  =  ''
    ThisLoc =  ' -> at MetaHistItem_Destroy (in History/metahistitem_mod.F90)'

    !=======================================================================
    ! Destroy each METAHISTORY ITEM in the linked list
    !=======================================================================

    ! Point to the first METAHISTORY ITEM (aka node) in the list
    Node => List

    ! As long as it doesn't point to NULL()
    DO WHILE ( ASSOCIATED( Node ) )

       ! Set the CURRENT pointer to the current METAHISTORY ITEM
       Current => Node

       ! Destroy the HISTORY ITEM contained within this METAHISTORY ITEM
       CALL HistItem_Destroy( Current%Item, RC )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Cannot deallocate the "Current%Item" HISTORY ITEM!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       ! Point to the next METAHISTORY ITEM for the next iteration
       Node => Current%Next

       ! And destroy the current METAHISTORY ITEM in the list
       DEALLOCATE( Current, STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Cannot deallocate the "Current" META HISTORY ITEM!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

    ENDDO

    ! Free pointers
    Current => NULL()
    Node    => NULL()

  END SUBROUTINE MetaHistItem_Destroy
!EOC
END MODULE MetaHistItem_Mod

