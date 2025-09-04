!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: metahistcontainer_mod.F90
!
! !DESCRIPTION: Contains types and methods to create a METAHISTORY CONTAINER
!  object, which is a container for a HISTORY CONTAINER.  In other words,
!  a METAHISTORY CONTAINER represents a single node of a linked list that is
!  used to contain HISTORY CONTAINERS.
!\\
!\\
!  In practice, we can think of a METAHISTORY CONTAINER as a list of
!  diagnostic collections, each of which contains a list of HISTORY ITEMS
!  to be archived to netCDF output at a specified frequency (e.g.
!  instantaneous, daily, hourly, etc.)
!!\\
!\\
! !INTERFACE:
!
MODULE MetaHistContainer_Mod
!
! !USES:
!
  USE HistContainer_Mod, ONLY : HistContainer
  USE Precision_Mod

  IMPLICIT NONE
  PRIVATE
!
! !PRIVATE MEMBER FUNCTIONS:
!
  PRIVATE :: MetaHistContainer_Create
  PRIVATE :: MetaHistContainer_Insert
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC  :: MetaHistContainer_AddNew
  PUBLIC  :: MetaHistContainer_Count
  PUBLIC  :: MetaHistContainer_Destroy
  PUBLIC  :: MetaHistContainer_Print
!
! !PUBLIC TYPES:
!
  !=========================================================================
  ! This is the derived type for a METAHISTORY CONTAINER object, which
  ! represents a SINGLE NODE OF A LINKED LIST consisting of HISTORY
  ! CONTAINERS.
  !
  ! As such, the METAHISTORY CONTAINER does not contain any data itself,
  ! but is a wrapper for a single HISTORY CONTAINER object, plus a pointer
  ! to another METAHISTORY CONTAINER (i.e. the next node in the list).
  !=========================================================================
  TYPE, PUBLIC ::  MetaHistContainer

     ! Pointer to the next METAHISTORY CONTAINER object
     ! (i.e. the next node in the linked list)
     TYPE(MetaHistContainer), POINTER :: Next      => NULL()

     ! The HISTORY CONTAINER object (which represents a diagnostic
     ! quantity that will be archived to netCDF file format)
     TYPE(HistContainer),     POINTER :: Container => NULL()

  END TYPE MetaHistContainer
!
! !REMARKS:
!  As described above, a METAHISTORY CONTAINER can be thought of as a SINGLE
!  NODE OF A LINKED LIST INTENDED TO HOLD HISTORY CONTAINERS.  It looks like
!  this:
!
!      +----------------------------+   +----------------------------+
!      | METAHISTORY CONTAINER n    |   | METAHISTORY CONTAINER n+1  |
!      | (aka NODE n of list)       |   | (aka NODE n+1 of list)     |
!      |                            |   |                            |
!      | Contains:                  |   | Contains:                  |
!      |                            |   |                            |
!      |   HISTORY CONTAINER n      |   |   HISTORY CONTAINER n+1    |
!      |                            |   |                            |
! =======> Pointer to next    ============> Pointer to next    =========> etc
!      |    METAHISTORY CONTAINER   |   |    METAHISTORY CONTAINER   |
!      +----------------------------+   +----------------------------+
!
!  Linked list routines taken from original code (linkedlist.f90)
!  by Arjen Markus; http://flibs.sourceforge.net/linked_list.html

! !REVISION HISTORY:
!  16 Jun 2017 - R. Yantosca - Initial version, based on code by Arjen Markus
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
! !IROUTINE: MetaHistContainer_AddNew
!
! !DESCRIPTION: Wrapper for methods MetaHistContainer\_Create and
!  MetaHistContainer\_Insert.  Will create a METAHISTORY CONTAINER (containing
!  a HISTORY CONTAINER) and (1) set it as the head node of a new linked list,
!  or (2) append it to an existing linked list.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE MetaHistContainer_AddNew( Input_Opt, Node, Container, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE HistContainer_Mod, ONLY : HistContainer
    USE Input_Opt_Mod,     ONLY : OptInput
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput),          INTENT(IN)  :: Input_Opt ! Input Options object
    TYPE(HistContainer),     POINTER     :: Container ! HISTORY CONTAINER
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(MetaHistContainer), POINTER     :: Node      ! METAHISTORY CONTAINER
!
! !OUTPUT PARAMETERS:
!
    INTEGER,                 INTENT(OUT) :: RC        ! Success or failure
!
! !REVISION HISTORY:
!  16 Jun 2017 - R. Yantosca - Initial version, based on code by Arjen Markus
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
    RC      = GC_SUCCESS
    ErrMsg  = ''
    ThisLoc = &
      ' -> at MetaHistContainer_Add (in History/metahistcontainer_mod.F90)'

    !=======================================================================
    ! Test if the METAHISTORY CONTAINER (aka "Node") has been allocated
    ! memory  and is therefore part of an existing linked list
    !=======================================================================
    IF ( .not. ASSOCIATED( Node ) ) THEN

       !--------------------------------------------------------------------
       ! If not, then create a new METAHISTORY CONTAINER (named "Node"),
       ! and set it at the head of a new linked list
       !--------------------------------------------------------------------
       CALL MetaHistContainer_Create( Input_Opt, Node, Container, RC )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Could not create "Node" as the head node of a list!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

    ELSE

       !--------------------------------------------------------------------
       ! Otherwise, create a new METAHISTORY CONTAINER (named "Node"),
       ! and append it to the list, immediately following the head node
       !--------------------------------------------------------------------
       CALL MetaHistContainer_Insert( Input_Opt, Node, Container, RC )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Could not insert "Node" into an existing linked list!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

    ENDIF

  END SUBROUTINE MetaHistContainer_AddNew
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: MetaHistContainer_Create
!
! !DESCRIPTION: This method creates a new METAHISTORY CONTAINER (to contain the
!  supplied HISTORY CONTAINER) and sets it as the head node of a linked list.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE MetaHistContainer_Create( Input_Opt, Node, Container, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE HistContainer_Mod, ONLY : HistContainer
    USE Input_Opt_Mod,     ONLY : OptInput
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput),          INTENT(IN)  :: Input_Opt ! Input Options object
    TYPE(HistContainer),     POINTER     :: Container ! HISTORY CONTAINER
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(MetaHistContainer), POINTER     :: Node      ! METAHISTORY CONTAINER
!
! !OUTPUT PARAMETERS:
!
    INTEGER,                 INTENT(OUT) :: RC        ! Success or failure
!
! !REMARKS:
!  This method is not intended to be called directly, but is rather
!  wrapped by the MetaHistContainer_AddNew method.
!
! !REVISION HISTORY:
!  16 Jun 2017 - R. Yantosca - Initial version, based on code by Arjen Markus
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
    RC      = GC_SUCCESS
    ErrMsg  = ''
    ThisLoc = &
      ' -> at MetaHistContainer_Create (in History/metahistcontainer_mod.F90)'

    !=======================================================================
    ! Initialize the METAHISTORY CONTAINER itself
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
    ! Initialize the field that will store the HISTORY CONTAINER
    !=======================================================================

    ! Because this is the first METAHISTORY CONTAINER that is being created,
    ! we can consider this to be the head node of a linked list.
    IF ( .not. ASSOCIATED( Node%Container ) ) THEN
       ALLOCATE( Node%Container, STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Could not allocate "Node%Container"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDIF

    ! Attach the given HISTORY CONTAINER to the METAHISTORY CONTAINER
    ! (i.e. place it into the head node of a linked list)
    Node%Container = Container

  END SUBROUTINE MetaHistContainer_Create
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: MetaHistContainer_Insert
!
! !DESCRIPTION: Creates a new METAHISTORY CONTAINER (to contain the supplied
!  HISTORY CONTAINER), and pops it into an existing linked list, immediately
!  following the head node.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE MetaHistContainer_Insert( Input_Opt, Node, Container, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE HistContainer_Mod, ONLY : HistContainer
    USE Input_Opt_Mod,     ONLY : OptInput
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput),          INTENT(IN)  :: Input_Opt ! Input Options object
    TYPE(HistContainer),     POINTER     :: Container ! HISTORY CONTAINER
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(MetaHistContainer), POINTER     :: Node      ! METAHISTORY CONTAINER
!
! !OUTPUT PARAMETERS:
!
    INTEGER,                 INTENT(OUT) :: RC        ! Success or failure
!
! !REMARKS:
!  This method is not intended to be called directly, but is rather
!  wrapped by the MetaHistContainer_AddNew method.
!
! !REVISION HISTORY:
!  16 Jun 2017 - R. Yantosca - Initial version, based on code by Arjen Markus
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Strings
    CHARACTER(LEN=255)               :: ErrMsg, ThisLoc

    ! Objects
    TYPE(MetaHistContainer), POINTER :: Head

    !=======================================================================
    ! Initialize
    !=======================================================================
    RC      = GC_SUCCESS
    ErrMsg  = ''
    ThisLoc = &
      ' -> at MetaHistContainer_Insert (in History/metahistcontainer_mod.F90)'

    !=======================================================================
    ! Initialize a METAHISTORY CONTAINER named "Head", which will
    ! become  the head of the existing list.  "Head" will contain
    ! a new HISTORY CONTAINER.
    !=======================================================================

    ! Allocate the "Head" object
    ALLOCATE( Head, STAT=RC )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Could not allocate "Head"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    ! Allocate the "HeadContainer" field,
    ! which will hold the HISTORY CONTAINER
    ALLOCATE( Head%Container, STAT=RC )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Could not allocate "Head%Container"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !=======================================================================
    ! Insert "Head" at the start of the existing linked list
    !=======================================================================

    ! Save the HISTORY CONTAINER argument in the "Container" field of "Head"
    Head%Container =  Container

    ! The "Next" field of "Head" points to the current head of the list
    Head%Next      => Node

    ! Set "Head" as the new head of the linked list
    Node           => Head

  END SUBROUTINE MetaHistContainer_Insert
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: MetaHistContainer_Count
!
! !DESCRIPTION: Counts the number of METAHISTORY CONTAINERS stored in a linked
!  list.  By extension, this is also the number of HISTORY CONTAINERS stored in
!  the list, because each METAHISTORY CONTAINER contains only one HISTORY !
!  CONTAINER.
!\\
!\\
! !INTERFACE:
!
  FUNCTION MetaHistContainer_Count( List ) RESULT( nNodes )

!
! !INPUT PARAMETERS:
!
    TYPE(MetaHistContainer), POINTER :: List   ! List of METAHISTORY CONTAINERS
!
! !RETURN VALUE:
!
    INTEGER                          :: nNodes ! # of METAHISTORY CONTAINERS
                                               ! (aka nodes) in the list
!
! !REVISION HISTORY:
!  16 Jun 2017 - R. Yantosca - Initial version, based on code by Arjen Markus
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Objects
    TYPE(MetaHistContainer), POINTER  :: Current

    !=======================================================================
    ! Initialize
    !=======================================================================
    nNodes  =  0
    Current => NULL()

    !=======================================================================
    ! Count the number of METAHISTORY CONTAINERS in the list
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

  END FUNCTION MetaHistContainer_Count
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: MetaHistContainer_Print
!
! !DESCRIPTION: This method will print information about the HISTORY CONTAINER
!  belonging to each METAHISTORY CONTAINER (aka node) of a linked list.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE MetaHistContainer_Print( Input_Opt, List, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE HistContainer_Mod, ONLY : HistContainer, HistContainer_Print
    USE Input_Opt_Mod,     ONLY : OptInput
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput),          INTENT(IN)  :: Input_Opt ! Input Options object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(MetaHistContainer), POINTER     :: List      ! List of METAHISTORY
                                                      !  CONTAINERS
! !OUTPUT PARAMETERS:
!
    INTEGER,                 INTENT(OUT) :: RC        ! Success or failure
!
! !REMARKS:
!
! !REVISION HISTORY:
!  16 Jun 2017 - R. Yantosca - Initial version, based on code by Arjen Markus
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Strings
    CHARACTER(LEN=255)               :: ErrMsg, ThisLoc

    ! Objects
    TYPE(MetaHistcontainer), POINTER :: Current

    !=======================================================================
    ! Initialize
    !=======================================================================
    RC      =  GC_SUCCESS
    Current => NULL()
    ErrMsg  =  ''
    ThisLoc =  &
      ' -> at MetaHistContainer_Print (in History/metahistcontainer_mod.F90)'

    !=======================================================================
    ! Print information about each METAHISTORY CONTAINER (aka node)
    ! of the linked list, only if we are on the root CPU.
    !=======================================================================
    IF ( Input_Opt%amIRoot ) THEN

       ! Point CURRENT to the head node of the list
       Current => List

       ! As long as the current node is valid
       DO WHILE( ASSOCIATED( Current ) )

          ! Print info about the history container corresponding to this node
          CALL HistContainer_Print( Input_Opt, Current%Container, RC )
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = 'Could not print info for "Current%Container"!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF

          ! Point to the next node for the next iteration
          Current => Current%Next
       ENDDO

       ! Free pointers
       Current => NULL()
    ENDIF

  END SUBROUTINE MetaHistContainer_Print
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: MetaHistContainer_Destroy
!
! !DESCRIPTION:  This method will destroy the HISTORY CONTAINER belonging to
!  each METAHISTORY CONTAINER (aka node) of a linked list.  It will then
!  destroy each METAHISTORY CONTAINER in the list.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE MetaHistContainer_Destroy( List, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE HistContainer_Mod, ONLY : HistContainer, HistContainer_Destroy
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(MetaHistContainer), POINTER     :: List      ! List of METAHISTORY
                                                      !  CONTAINERS
!
! !OUTPUT PARAMETERS:
!
    INTEGER,                 INTENT(OUT) :: RC        ! Success or failure?
!
! !REVISION HISTORY:
!  16 Jun 2017 - R. Yantosca - Initial version, based on code by Arjen Markus
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Strings
    CHARACTER(LEN=255)                :: ErrMsg,  ThisLoc

    ! Objects
    TYPE(MetaHistContainer), POINTER  :: Current, Node

    !=======================================================================
    ! Initialize
    !=======================================================================
    RC      =  GC_SUCCESS
    Current => NULL()
    Node    => NULL()
    ErrMsg  =  ''
    ThisLoc =  &
      ' -> at MetaHistContainer_Destroy (in History/metahistcontainer_mod.F90)'

    !=======================================================================
    ! Destroy each METAHISTORY CONTAINER in the linked list
    !=======================================================================

    ! Point to the first METAHISTORY CONTAINER (aka node) in the list
    Node => List

    ! As long as it doesn't point to NULL()
    DO WHILE ( ASSOCIATED( Node ) )

       ! Set the CURRENT pointer to the current METAHISTORY CONTAINER
       Current => Node

       ! Destroy the HISTORY CONTAINER contained within
       ! this METAHISTORY CONTAINER
       CALL HistContainer_Destroy( Current%Container, RC )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = &
            'Cannot deallocate the "Current%Container" HISTORY CONTAINER!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       ! Point to the next METAHISTORY CONTAINER for the next iteration
       Node => Current%Next

       ! And destroy the current METAHISTORY CONTAINER in the list
       DEALLOCATE( Current, STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Cannot deallocate the "Current" METAHISTORY CONTAINER!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

    ENDDO

    ! Free pointers
    Current => NULL()
    Node    => NULL()

  END SUBROUTINE MetaHistContainer_Destroy
!EOC
END MODULE MetaHistContainer_Mod

