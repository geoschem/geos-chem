!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: taggeddiaglist_mod.F90
!
! !DESCRIPTION: Module taggeddiaglist\_mod contains type definitions and
!  routines to define link list with detailed diagnostic information for each
!  State_Diag diagnostic contained in the DiagList object.
!\\
!\\
! !INTERFACE:
!
MODULE TaggedDiagList_Mod
!
! !USES:
!
  USE DiagList_Mod
  USE ErrCode_Mod
  USE Precision_Mod

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC  :: Cleanup_TaggedDiagList
  PUBLIC  :: Init_TaggedDiagList
  PUBLIC  :: Print_TaggedDiagItem
  PUBLIC  :: Print_TaggedDiagList
  PUBLIC  :: Print_TagList
  PUBLIC  :: Query_Tag_in_TagList
  PUBLIC  :: Query_TaggedDiagList
!
! !PRIVATE MEMBER FUNCTIONS:
!
  PRIVATE :: Cleanup_TagList
  PRIVATE :: Init_TagItem
  PRIVATE :: Init_TaggedDiagItem
  PRIVATE :: InsertBeginning_TagList
  PRIVATE :: InsertBeginning_TaggedDiagList
  PRIVATE :: Update_TaggedDiagList
!
! !PUBLIC DATA TYPES:
!
  !=========================================================================
  ! Derived type for tag items (can be used for any linked list of strings)
  !=========================================================================
  TYPE, PUBLIC :: DgnTagItem
     CHARACTER(LEN=63)             :: name         ! Tag or wildcard name
     INTEGER                       :: index        ! Position in the list
     TYPE(DgnTagItem), POINTER     :: next         ! Points to next tag/WC
  END TYPE DgnTagItem

  !=========================================================================
  ! Derived type for tag list (can be used for any linked list of strings)
  !=========================================================================
  TYPE, PUBLIC :: DgnTagList
     INTEGER                       :: count        ! # of tags/WCs in list
     TYPE(DgnTagItem), POINTER     :: head         ! The start of the list
  END TYPE DgnTagList

  !=========================================================================
  ! Derived type for tagged diagnostic items, e.g. DryDep
  !=========================================================================
  TYPE, PUBLIC :: TaggedDgnItem
     CHARACTER(LEN=63)             :: metadataID   ! Diagnostic name
     LOGICAL                       :: isWildcard   ! Is it a wildcard?
     TYPE(DgnTagList)              :: tagList      ! Tags for this diagnostic
     TYPE(DgnTagList)              :: wildcardList ! WCs for this diagnostic
     TYPE(TaggedDgnItem), POINTER  :: next         ! Points to next diagnostic
  END TYPE TaggedDgnItem

  !=========================================================================
  ! Derived type for tagged diagnostic list, e.g. DryDep, SpeciesConc, etc
  !=========================================================================
  TYPE, PUBLIC :: TaggedDgnList
     TYPE(TaggedDgnItem), POINTER  :: head         ! Start of the list
  END TYPE TaggedDgnList
!
! !REVISION HISTORY:
!  18 Nov 2019 - E. Lundgren - Initial version
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
! !IROUTINE: Init_TaggedDiagList
!
! !DESCRIPTION: Subroutine Init\_TaggedDiagList initializes the TaggedDiagList
!  corresponding to each diagnostic collection.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Init_TaggedDiagList( am_I_Root, DiagList, TaggedDiagList, RC )
!
! !USES:
!
    USE Charpak_Mod,   ONLY : To_UpperCase
!
! !INPUT PARAMETERS:
!
    LOGICAL,              INTENT(IN)    :: am_I_Root       ! On root thread?
    TYPE(DgnList),        INTENT(IN)    :: DiagList        ! Diagnostics List
!
! !INPUT AND OUTPUT PARAMETERS:
!
    TYPE(TaggedDgnList),  INTENT(INOUT) :: TaggedDiagList  ! Tagged Diag List
!
! !OUTPUT PARAMETERS:
!
    INTEGER,              INTENT(OUT)   :: RC              ! Success/failure?
!
! !REVISION HISTORY:
!  18 Nov 2019 - E. Lundgren - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    LOGICAL                      :: isWildCard
    LOGICAL                      :: taggedDiagListExists
    INTEGER                      :: listIndex
    INTEGER                      :: numTags
    INTEGER                      :: numWildCards

    ! Strings
    CHARACTER(LEN=63)            :: tagName
    CHARACTER(LEN=255)           :: errMsg
    CHARACTER(LEN=255)           :: thisLoc

    ! Objects
    TYPE(DgnItem),       POINTER :: diagnostic
    TYPE(DgnTagItem),    POINTER :: current
    TYPE(TaggedDgnItem), POINTER :: TaggedDiagItem
    TYPE(TaggedDgnItem), POINTER :: NewTaggedDiagItem

    !=======================================================================
    ! Init_TaggedDiagList begins here!
    !=======================================================================

    ! Initialize
    RC      = GC_SUCCESS
    errMsg  = ''
    thisLoc = &
     '-> at Init_TaggedDiagList (in module Headers/taggeddiaglist_mod.F90)'

    ! Initialize pointers
    diagnostic          => DiagList%head
    current             => NULL()
    TaggedDiagList%head => NULL()

    !=======================================================================
    ! For each State_Diag diagnostic containd in DiagList,
    ! initialize the corresponding TaggedDiagList.
    !=======================================================================
    DO WHILE ( ASSOCIATED( diagnostic ) )

       ! Only proceed for State_Diag diagnostics
       IF ( diagnostic%state == 'DIAG' .AND.                                 &
            ( diagnostic%isTagged .OR. diagnostic%isWildcard ) ) THEN

          !--------------------------------------------------------------
          ! First check if the the TaggedDiagList corresponding
          ! to this diagnostic exists or not
          !--------------------------------------------------------------
          CALL Query_TaggedDiagList(                                         &
               TaggedDiagList = TaggedDiagList,                              &
               diagName       = diagnostic%metadataID,                       &
               Found          = taggedDiagListExists,                        &
               RC             = RC                                          )

          ! Trap potential errors
          IF ( RC /= GC_SUCCESS ) THEN
             errMsg = 'Error encountered in "Query_TaggedDiagList" (#1)!'
             CALL GC_Error( errMsg, RC, thisLoc )
             RETURN
          ENDIF

          ! Name of each tag or wildcard in the current State_Diag diagnostic
          IF ( diagnostic%isTagged ) THEN
             tagName = diagnostic%tag
          ELSE
             tagName = diagnostic%wildcard
          ENDIF

          ! Check if the list exists
          IF ( taggedDiagListExists ) THEN

             !--------------------------------------------------------------
             ! If the TaggedDiagList (of wildcards or tags) already exists:
             ! (1) Add a new TaggedDiagItem into it
             ! (2) Set the index to a placeholder.
             !     Indices will be updated in the following section.
             !--------------------------------------------------------------
             CALL Update_TaggedDiagList(                                     &
                  metadataID        = diagnostic%metadataID,                 &
                  isWildCard        = diagnostic%isWildCard,                 &
                  tagName           = tagName,                               &
                  index             = 0,                                     &
                  TaggedDiagList    = TaggedDiagList,                        &
                  RC                = RC                                    )

             ! Trap potential errors
             IF ( RC /= GC_SUCCESS ) THEN
                errMsg = 'Error encountered in "Update_TaggedDiagList"!'
                CALL GC_Error( errMsg, RC, thisLoc )
                RETURN
             ENDIF

          ELSE

             !--------------------------------------------------------------
             ! If the TaggedDiagList does not exist:
             ! (1) Create a new TaggedDiagItem
             ! (2) Set the index to a placeholder.
             !     Indices will be updated in the following section.
             !--------------------------------------------------------------
             CALL Init_TaggedDiagItem(                                       &
                  NewTaggedDiagItem = NewTaggedDiagItem,                     &
                  metadataID        = diagnostic%metadataID,                 &
                  isWildcard        = diagnostic%isWildCard,                 &
                  tagName           = tagName,                               &
                  index             = 0,                                     &
                  RC                = RC                                    )

             ! Trap potential errors
             IF ( RC /= GC_SUCCESS ) THEN
                errMsg = 'Error encountered in "Init_TaggedDiagItem"!'
                CALL GC_Error( errMsg, RC, thisLoc )
                RETURN
             ENDIF

             !--------------------------------------------------------------
             ! (3) Create the TaggedDiagList
             ! (4) Add the TaggedDiagItem to the head of the TaggedDiagList
             !--------------------------------------------------------------
             CALL InsertBeginning_TaggedDiagList(                            &
                  TaggedDiagItem    = NewTaggedDiagItem,                     &
                  TaggedDiagList    = TaggedDiagList,                        &
                  RC                = RC                                    )

             ! Trap potential errors
             IF ( RC /= GC_SUCCESS ) THEN
                errMsg = &
                   'Error encountered in "InsertBeginning_TaggedDiagList"!'
                CALL GC_Error( errMsg, RC, thisLoc )
                RETURN
             ENDIF
          ENDIF
       ENDIF

       ! Advance to next item in DiagList
       diagnostic => diagnostic%next
    ENDDO

    !=======================================================================
    ! Now cycle through each diagnostic in DiagList again and reverse
    ! the indices of the tags or wildcards into proper order.  We have to
    ! do this as a second pass because in the first loop over diagnostics
    ! above, we didn't yet know the total number of tags or wildcards.
    !=======================================================================

    ! Re-initialize pointers
    diagnostic => DiagList%head

    DO WHILE ( ASSOCIATED( diagnostic ) )

       ! Initialize
       TaggedDiagItem => NULL()
       current        => NULL()
       numTags        =  0
       numWildCards   =  0

       ! Only proceed for State_Diag diagnostics
       IF ( diagnostic%state == 'DIAG' .AND.                                 &
            ( diagnostic%isTagged .OR. diagnostic%isWildcard ) ) THEN

          !-----------------------------------------------------------------
          ! Get the number of tags and wildcards in the TaggedDiagList
          ! that corresponds to this State_Diag diagnostic
          !-----------------------------------------------------------------
          CALL Query_TaggedDiagList(                                         &
               TaggedDiagList = TaggedDiagList,                              &
               diagName       = diagnostic%metadataID,                       &
               Found          = taggedDiagListExists,                        &
               isWildCard     = isWildCard,                                  &
               numTags        = numTags,                                     &
               numWildCards   = numWildCards,                                &
               RC             = RC                                          )

          ! Trap potential errors
          IF ( RC /= GC_SUCCESS ) THEN
             errMsg = 'Error encountered in "Query_TaggedDiagList" (#2)!'
             CALL GC_Error( errMsg, RC, thisLoc )
             RETURN
          ENDIF

          ! Skip if the TaggedDiagList for this diagnostic does not exist
          ! (but this should never happen since we created a TaggedDiagList
          ! for each diagnostic in the above section!)
          IF ( .not. taggedDiagListExists ) CYCLE

          !-----------------------------------------------------------------
          ! Loop over all tags of wildcards belonging to this diagnostic
          !-----------------------------------------------------------------
          TaggedDiagItem => TaggedDiagList%head

          DO WHILE ( ASSOCIATED( TaggedDiagItem ) )

             ! Make sure we the diagnostic field name (metadataId) in the
             ! TaggedDiagItem matches that of the current diagnostic
             IF ( TRIM( TaggedDiagItem%metaDataId ) ==                      &
                  TRIM( diagnostic%metaDataId     )     ) THEN

                !-----------------------------------------------------------
                ! Reset the indices of entries in WildCardList
                !-----------------------------------------------------------
                IF ( isWildCard ) THEN
                   current   => TaggedDiagItem%wildCardList%head
                   listIndex =  0
                   DO WHILE ( ASSOCIATED( current ) )
                      listIndex     =  listIndex + 1
                      current%index =  listIndex
                      current       => current%next
                   ENDDO
                   current => NULL()

                !-----------------------------------------------------------
                ! Reset the indices of entries in TagList
                !-----------------------------------------------------------
                ELSE
                    current   => TaggedDiagItem%tagList%head
                    listIndex =  0
                    DO WHILE ( ASSOCIATED( current ) )
                       listIndex     =  listIndex + 1
                       current%index =  listIndex
                       current       => current%next
                    ENDDO
                    current => NULL()

                ENDIF
             ENDIF

             ! Advance to the next item in TaggedDiagList
             TaggedDiagItem => TaggedDiagItem%next

          ENDDO

          ! Free pointers
          TaggedDiagItem => NULL()
       ENDIF

       ! Advance to next diagnostic in DiagList
       diagnostic => diagnostic%next
    ENDDO

    !-----------------------------------------------------------------------
    ! Cleanup and quit
    !-----------------------------------------------------------------------
    diagnostic => NULL()
    current    => NULL()

  END SUBROUTINE Init_TaggedDiagList
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Print_TaggedDiagItem
!
! !DESCRIPTION: Prints information contained in a single TaggedDiagItem object.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Print_TaggedDiagItem( am_I_Root, TaggedDiagItem, RC )
!
! !INPUT PARAMETERS:
!
    LOGICAL,             INTENT(IN)  :: am_I_Root
    TYPE(TaggedDgnItem), INTENT(IN)  :: TaggedDiagItem
!
! !OUTPUT PARAMETERS:
!
    INTEGER,             INTENT(OUT) :: RC
!
! !REVISION HISTORY:
!  24 Mar 2020 - R. Yantosca - Initial version
!  See the subsequent Git history with the gitk browser!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Return success
    RC = GC_SUCCESS

    ! Proceed to print info if we are on the root CPU
    IF ( am_I_Root ) THEN

       ! Print name of diagnostic to which this TaggedDataItem belongs
       WRITE( 6, 120 ) TRIM( TaggedDiagItem%metadataID )

       IF ( TaggedDiagItem%isWildCard ) THEN

          !---------------------------------
          ! Print info about wildcards ...
          !---------------------------------
          WRITE( 6, 130 ) ADJUSTL( 'numWildCards:' ),                        &
                          TaggedDiagItem%WildCardList%count
          CALL Print_TagList( am_I_Root, TaggedDiagItem%WildCardList, RC )

       ELSE

          !----------------------------------
          ! ... or print info about tags
          !----------------------------------
          WRITE( 6, 130 )  ADJUSTL( 'numTags:' ),                            &
                           TaggedDiagItem%TagList%count
          CALL Print_TagList( am_I_Root, TaggedDiagItem%TagList, RC )

       ENDIF

       ! Print spacer
       WRITE( 6, 120 ) ""

       ! FORMAT statemetns
 120   FORMAT( A       )
 130   FORMAT( A15, I5 )

    ENDIF

  END SUBROUTINE Print_TaggedDiagItem
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Print_TaggedDiagList
!
! !DESCRIPTION: Prints information for all TaggedDiagItem members within
!  a TaggedDiagList linked list.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Print_TaggedDiagList( am_I_Root, TaggedDiagList, RC )
!
! !INPUT PARAMETERS:
!
    LOGICAL,             INTENT(IN)  :: am_I_Root
    TYPE(TaggedDgnList), INTENT(IN)  :: TaggedDiagList
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,             INTENT(OUT) :: RC
!
! !REVISION HISTORY:
!  18 Nov 2019 - E. Lundgren - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    TYPE(TaggedDgnItem), POINTER :: current

    ! ================================================================
    ! Print_TaggedDiagList begins here
    ! ================================================================

    ! Initialize
    RC      =  GC_SUCCESS
    current => NULL()

    ! Print tagged diagnostic list only if we are on the root core
    IF ( am_I_Root ) THEN
       WRITE( 6, 100 ) REPEAT( '=', 30 )
       WRITE( 6, 110 ) 'Summary of tagged diagnostics'

       ! Point to the first item in the TaggedDiagList
       current => TaggedDiagList%head

       ! Keep looping over all items in TaggedDiagList
       DO WHILE ( ASSOCIATED( current ) )

          ! Print wildcard or tag info
          CALL Print_TaggedDiagItem( am_I_Root, current, RC )

          ! Advance to next item in TaggedDiagList
          current => current%next
       ENDDO

       ! Spacer
       WRITE( 6, 120 ) ""

       ! Free pointer
       current => NULL()

       ! FORMAT statements
 100   FORMAT( /,   A  )
 110   FORMAT( A,   /  )
 120   FORMAT( A       )
 130   FORMAT( A15, I5 )

    ENDIF

  END SUBROUTINE Print_TaggedDiagList
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Query_TaggedDiagList
!
! !DESCRIPTION: Returns information about a TaggedDiagItem within a
!  TaggedDiagList option linked list, given the name of the corresponding
!  diagnostic.  The information that is returned can include if the
!  diagnostic is a wildcard, a list of wildcards, number of wildcards, a
!  list of non-wildcard tags, or number of non-wildcard tags.   The entire
!  TaggedDiagItem itself may also be returned if so desired.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Query_TaggedDiagList( TaggedDiagList,                           &
                                   diagName,      RC,                        &
                                   found,         isWildcard,                &
                                   numWildcards,  numTags,                   &
                                   WildCardList,  TagList,                   &
                                   TaggedDiagItem                           )
!
! !USES:
!
    USE Charpak_Mod,   ONLY : To_UpperCase
!
! !INPUT PARAMETERS:
!
    TYPE(TaggedDgnList), INTENT(IN)  :: TaggedDiagList  ! Tagged diag list
    CHARACTER(LEN=*),    INTENT(IN)  :: diagName        ! Name of diagnostic
!
! !OUTPUT PARAMETERS:
!
    LOGICAL,             INTENT(OUT) :: found           ! Was name in list?
    INTEGER,             INTENT(OUT) :: RC              ! Success or failure?
    LOGICAL,             OPTIONAL    :: isWildcard      ! Is diag a wildcard?
    INTEGER,             OPTIONAL    :: numWildcards    ! # of wildcards (WCs)
    TYPE(DgnTagList),    OPTIONAL    :: WildCardList    ! List of wildcards
    INTEGER,             OPTIONAL    :: numTags         ! # of non-WC tags
    TYPE(DgnTagList),    OPTIONAL    :: TagList         ! List of non-WC tags
    TYPE(TaggedDgnItem), OPTIONAL    :: TaggedDiagItem  ! Item in linked list
!
! !REVISION HISTORY:
!  18 Nov 2019 - E. Lundgren - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    LOGICAL                      :: P_isWildCard
    LOGICAL                      :: P_numWildCards
    LOGICAL                      :: P_numTags
    LOGICAL                      :: P_WildCardList
    LOGICAL                      :: P_TagList
    LOGICAL                      :: P_TagDiagItem

    ! Strings
    CHARACTER(LEN=255)           :: thisDiagName

    ! Objects
    TYPE(TaggedDgnItem), POINTER :: current

    !=======================================================================
    ! Query_TaggedDiagList begins here!
    !=======================================================================

    ! Initialize
    RC    = GC_SUCCESS
    found = .FALSE.

    ! Test if each optional arguments is present outside of the loop
    P_isWildCard   = PRESENT( isWildcard      )
    P_numWildCards = PRESENT( numWildcards    )
    P_numTags      = PRESENT( numTags         )
    P_WildCardList = PRESENT( WildCardList    )
    P_TagList      = PRESENT( TagList         )
    P_TagDiagItem  = PRESENT( TaggedDiagItem  )

    ! Search for name in list and return optional arguments
    current => TaggedDiagList%head
    DO WHILE ( ASSOCIATED( current ) )
       thisDiagName = To_Uppercase( current%metadataID )
       IF ( TRIM( ThisDiagName ) == TRIM( To_Uppercase( diagName ) ) ) THEN
          found = .TRUE.
          IF ( P_isWildcard   ) isWildcard      = current%isWildcard
          IF ( P_numWildcards ) numWildcards    = current%WildCardList%count
          IF ( P_numTags      ) numTags         = current%TagList%count
          IF ( P_WildCardList ) WildCardList    = current%WildCardList
          IF ( P_TagList      ) TagList         = current%TagList
          IF ( P_TagDiagItem  ) TaggedDiagItem  = current
          EXIT
       ENDIF
       current => current%next
    ENDDO
    current => NULL()

  END SUBROUTINE Query_TaggedDiagList
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Cleanup_TaggedDiagList
!
! !DESCRIPTION: Subroutine Cleanup\_TaggedDiagList deallocates a TaggedDiagList
!  object and all of its member objects including the linked list of
!  TaggedDiagItem objects.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Cleanup_TaggedDiagList( TaggedDiagList, RC )
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(TaggedDgnList), INTENT(INOUT) :: TaggedDiagList
!
! !OUTPUT PARAMETERS:
!
    INTEGER,             INTENT(OUT)   :: RC
!
! !REVISION HISTORY:
!  18 Nov 2019 - E. Lundgren - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Strings
    CHARACTER(LEN=255)           :: errMsg
    CHARACTER(LEN=255)           :: thisLoc

    ! Pointers
    TYPE(TaggedDgnItem), POINTER :: current
    TYPE(TaggedDgnItem), POINTER :: next

    !=======================================================================
    ! Cleanup_TaggedDiagList begins here
    !=======================================================================

    ! Initialize
    RC      = GC_SUCCESS
    errMsg  = ''
    thisLoc = &
    ' -> at Cleanup_TaggedDiagList (in module Headers/taggeddiaglist_mod.F90)'

    ! Deallocate each item in the linked list of DiagExport objects
    current => TaggedDiagList%head
    IF ( ASSOCIATED( current ) ) next => current%next

    DO WHILE ( ASSOCIATED( current ) )

       ! First, free the list of tags
       CALL Cleanup_TagList( current%taglist, RC )
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = 'Error encountered in "Cleanup_TagList" (for tags)!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       ! Then free the list of wildcards
       CALL Cleanup_TagList( current%wildcardlist, RC )
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = 'Error encountered in "Cleanup_TagList" (for wildcards)!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       ! Next, free the current TaggedDgnItem
       DEALLOCATE( current, STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = 'Could not deallocate the "current" TaggedDgnItem object!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       ! Exit if we are at the end of the TaggedDiagList ...
       IF ( .NOT. ASSOCIATED ( next ) ) EXIT

       ! ...or if not, advance to next item in list
       current => next
       next => current%next
    ENDDO

    ! Free pointers
    current => NULL()
    next    => NULL()

  END SUBROUTINE Cleanup_TaggedDiagList
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Init_TaggedDiagItem
!
! !DESCRIPTION: Initializes a TaggedDiagItem object, which contains metadata
!  as well as the lists of tags and wildcards for each State_Diag diagnostic.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Init_TaggedDiagItem( NewTaggedDiagItem, metadataID, isWildcard,  &
                                  tagName,           index,      RC          )
!
! !INPUT PARAMETERS:
!
    CHARACTER(LEN=*),    INTENT(IN)  :: metadataID   ! Collection name
    LOGICAL,             INTENT(IN)  :: isWildcard   ! Is it a wildcard?
    CHARACTER(LEN=*),    INTENT(IN)  :: tagName      ! Name of each tag
    INTEGER,             INTENT(IN)  :: index        ! Position of each tag
!
! !OUTPUT PARAMETERS:
!
    TYPE(TaggedDgnItem), POINTER     :: NewTaggedDiagItem  ! New TagItem
    INTEGER,             INTENT(OUT) :: RC                 ! Success/failure?
!
! !REVISION HISTORY:
!  18 Nov 2019 - E. Lundgren - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Strings
    CHARACTER(LEN=255)        :: diagID
    CHARACTER(LEN=255)        :: errMsg
    CHARACTER(LEN=255)        :: thisLoc

    ! Objects
    TYPE(DgnTagItem), POINTER :: NewTagItem

    !=======================================================================
    ! Init_TaggedDiagItem begins here
    !=======================================================================

    ! Initialize
    RC      = GC_SUCCESS
    errMsg  = ''
    thisLoc = &
     ' -> at Init_TaggedDiagItem (in module Headers/taggeddiaglist_mod.F90)'

    ! Create a new entry for TaggedDiagList
    ALLOCATE( NewTaggedDiagItem, STAT=RC )
    diagId = 'taggeddiaglist_mod.F90:NewTaggedDiagItem'
    CALL GC_CheckVar( diagId, 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    !-----------------------------------------------------------------------
    ! Initialize fields of the TaggedDiagItem
    !-----------------------------------------------------------------------
    NewTaggedDiagItem%metadataID         =  TRIM(metadataID)
    NewTaggedDiagItem%isWildcard         =  isWildcard
    NewTaggedDiagItem%wildcardList%head  => NULL()
    NewTaggedDiagItem%wildcardList%count =  0
    NewTaggedDiagItem%tagList%head       => NULL()
    NewTaggedDiagItem%tagList%count      =  0

    !-----------------------------------------------------------------------
    ! Create a new DgnTagItem object, which represents an individual
    ! tag or wildcard for a given State_Diag diagnostic.
    !-----------------------------------------------------------------------
    CALL Init_TagItem(                                                       &
         NewTagItem   = NewTagItem,                                          &
         name         = tagName,                                             &
         index        = index,                                               &
         RC           = RC                                                  )

    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error encountered in "Init_TagItem"!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    IF ( isWildcard ) THEN

       !--------------------------------------------------------------------
       ! Add the DgnTagItem object to the list of wildcards
       ! belonging to the NewTaggedDiagItem ...
       !--------------------------------------------------------------------
       CALL InsertBeginning_TagList(                                         &
            TagItem   = NewTagItem,                                          &
            TagList   = NewTaggedDiagItem%wildcardList,                      &
            RC        = RC                                                  )

       ! Trap potential errors
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = 'Error encountered in "InsertBeginning_TagList" (wildcard)!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

    ELSE

       !--------------------------------------------------------------------
       ! ...or, if it is a tag and not a wildcard, add it to the list of
       ! tags belonging to the NewTaggedDiagItem object.
       !--------------------------------------------------------------------
       CALL InsertBeginning_TagList(                                         &
            TagItem   = NewTagItem,                                          &
            TagList   = NewTaggedDiagItem%tagList,                           &
            RC        = RC                                                  )

       ! Trap potential errors
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = 'Error encountered in "InsertBeginning_TagList" (wildcard)!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

    ENDIF

  END SUBROUTINE Init_TaggedDiagItem
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Init_TagItem
!
! !DESCRIPTION: Initializes a TagItem object, which represents a single tag
!  or wildcard belonging to a State_Diag diagnostic.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Init_TagItem( NewTagItem, name, index, RC )
!
! !INPUT PARAMETERS:
!
    TYPE(DgnTagItem), POINTER     :: NewTagItem    ! TagItem object
    CHARACTER(LEN=*), INTENT(IN)  :: name          ! Name of quantity
    INTEGER,          INTENT(IN)  :: index         ! Position of quantity
!
! !OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(OUT) :: RC            ! Success or failure?
!
! !REVISION HISTORY:
!  18 Nov 2019 - E. Lundgren - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Strigns
    CHARACTER(LEN=255) :: diagId

    !=================================================================
    ! Init_TagItem begins here
    !=================================================================

    ! Initialize
    RC  = GC_SUCCESS

    ! Allocate the NewTagItem Object
    ALLOCATE( NewTagItem, STAT=RC )
    diagId = 'taggeddiaglist_mod.F90:NewTagItem'
    CALL GC_CheckVar( diagId, 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    ! Define fields
    NewTagItem%name  = TRIM(name)
    NewTagItem%index = index

  END SUBROUTINE Init_TagItem
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: InsertBeginning_TaggedDiagList
!
! !DESCRIPTION: Inserts a new TaggedDiagItem (containing metadata and
!  lists of tags or wildcards for a single State_Diag diagnostic) to the
!  beginning of the TaggedDiagList linked list object.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE InsertBeginning_TaggedDiagList( TaggedDiagItem,    &
                                             TaggedDiagList, RC  )
!
! !INPUT PARAMETERS:
!
    TYPE(TaggedDgnItem), POINTER       :: TaggedDiagItem
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(TaggedDgnList), INTENT(INOUT) :: TaggedDiagList
!
! !OUTPUT PARAMETERS:
!
    INTEGER,             INTENT(OUT)   :: RC
!
! !REVISION HISTORY:
!  18 Nov 2019 - E. Lundgren - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Objects
    TYPE(TaggedDgnItem), POINTER :: NewTaggedDiagItem

    ! ================================================================
    ! InsertBeginning_TaggedDiagList begins here
    ! ================================================================

    ! Initialize
    RC = GC_SUCCESS

    ! Add new object to the beginning of the linked list
    TaggedDiagItem%next => TaggedDiagList%head
    TaggedDiagList%head => TaggedDiagItem

  END SUBROUTINE InsertBeginning_TaggedDiagList
!EOC
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: InsertBeginning_TagList
!
! !DESCRIPTION: Inserts a new node at the beginning of the TagList linked
!  list object.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE InsertBeginning_TagList( TagItem, TagList, RC )
!
! !INPUT PARAMETERS:
!
    TYPE(DgnTagItem), POINTER       :: TagItem     ! Tag or wildcard
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(DgnTagList), INTENT(INOUT) :: TagList     ! Tag list or wildcard list
!
! !OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(OUT)   :: RC          ! Success or failure
!
! !REVISION HISTORY:
!  18 Nov 2019 - E. Lundgren - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Objects
    TYPE(DgnTagItem), POINTER :: NewTagItem

    ! ================================================================
    ! InsertBeginning_TagList begins here
    ! ================================================================

    ! Initialize
    RC = GC_SUCCESS

    ! Add new object to the beginning of the linked list
    TagItem%next => TagList%head
    TagList%head => TagItem
    TagList%count = TagList%count + 1

  END SUBROUTINE InsertBeginning_TagList
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Query_Tag_in_TagList
!
! !DESCRIPTION: Searches for a given tag within a list of tags, or
!  a wildcard within a list of wildcards.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Query_Tag_in_TagList( TagList, name, found, RC )
!
! !USES:
!
    USE Charpak_Mod,   ONLY : To_UpperCase
!
! !INPUT PARAMETERS:
!
    TYPE(DgnTagList), INTENT(IN)  :: TagList     ! List of tags or wildcards
    CHARACTER(LEN=*), INTENT(IN)  :: name        ! Search string
!
! !OUTPUT PARAMETERS:
!
    LOGICAL,          INTENT(OUT) :: found       ! Is the tag/WC? in the list?
    INTEGER,          INTENT(OUT) :: RC          ! Success or failure
!
! !REVISION HISTORY:
!  18 Nov 2019 - E. Lundgren - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Pointers
    TYPE(DgnTagItem), POINTER :: current

    ! Strings
    CHARACTER(LEN=255)        :: thisTagName

    !=======================================================================
    ! Query_Tag_in_Taglist begins here!
    !=======================================================================

    ! Initialize
    RC       = GC_SUCCESS
    found    = .FALSE.

    ! Search for name in list
    current => TagList%head
    DO WHILE ( ASSOCIATED( current ) )
       thisTagName = To_Uppercase(current%name)
       IF ( TRIM(thisTagName) == TRIM(To_Uppercase(name)) ) THEN
          found = .TRUE.
          EXIT
       ENDIF
       current => current%next
    ENDDO

    ! Free pointer
    current => NULL()

  END SUBROUTINE Query_Tag_in_TagList
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Print_TagList
!
! !DESCRIPTION: Subroutine Print\_TagList prints information for all
!  TagItem members in a TagList linked list.  This represents a list of
!  tags or wildcards for a single State_Diag diagnostic.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Print_TagList( am_I_Root, TagList, RC )
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN)    :: am_I_Root
    TYPE(DgnTagList), INTENT(IN)    :: TagList    ! List of tags or wildcards
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(INOUT) :: RC         ! Success or failure?
!
! !REVISION HISTORY:
!  18 Nov 2019 - E. Lundgren - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    TYPE(DgnTagItem), POINTER :: current
    CHARACTER(LEN=255)        :: thisLoc

    ! ================================================================
    ! Print_TagList begins here
    ! ================================================================

    ! Initialize
    RC      =  GC_SUCCESS
    current => TagList%head

    ! Only print taglist if we are on the root core
    IF ( am_I_Root ) THEN
       DO WHILE ( ASSOCIATED( current ) )
          WRITE( 6, 100 ) ADJUSTL( TRIM( current%name ) ), current%index
 100      FORMAT( 21x, A, I5 )
          current => current%next
       ENDDO
    ENDIF
    current => NULL()

  END SUBROUTINE Print_TagList
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Update_TaggedDiagList
!
! !DESCRIPTION: Updates a TaggedDiagList object with a new TaggedDiagItem.
!  This represents adding a new tag or wildcard to an existing tag list or
!  wildcard list for a State_Diag diagnostic quantity.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Update_TaggedDiagList( metadataID, isWildcard,     tagName, &
                                    index,      TaggedDiagList, RC        )
!
! !USES:
!
    USE Charpak_Mod,   ONLY : To_UpperCase
!
! !INPUT PARAMETERS:
!
    CHARACTER(LEN=*),    INTENT(IN)    :: metadataID     ! Diagnostic name
    LOGICAL,             INTENT(IN)    :: isWildcard     ! Does this diagnostic
                                                         !  use a wildcard?
    CHARACTER(LEN=*),    INTENT(IN)    :: tagName        ! Tag or WC name
    INTEGER,             INTENT(IN)    :: index          ! Position of the
                                                         !  tag or WC in
                                                         !  the taglist/WClist
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(TaggedDgnList), INTENT(INOUT) :: TaggedDiagList ! List containing
                                                         !  metadata and
                                                         !  list of tags and
                                                         !  WCs per diagnostic
!
! !OUTPUT PARAMETERS:
!
    INTEGER,             INTENT(OUT)   :: RC             ! Success or failure?
!
! !REVISION HISTORY:
!  18 Nov 2019 - E. Lundgren - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Objects
    TYPE(TaggedDgnItem), POINTER :: current
    TYPE(DgnTagItem),    POINTER :: NewTagItem

    ! Strings
    CHARACTER(LEN=255)           :: errMsg
    CHARACTER(LEN=255)           :: thisLoc
    CHARACTER(LEN=255)           :: thisDiagName

    !=======================================================================
    ! Update_TaggedDiagList begins here
    !=======================================================================

    ! Initialize
    RC      = GC_SUCCESS
    errMsg  = ''
    thisLoc = ' -> at Update_TaggedDiagItem (in module taggeddiaglist_mod.F90)'

    ! Search for item in list and update if found
    current => TaggedDiagList%head
    DO WHILE ( ASSOCIATED( current ) )
       thisDiagName = To_Uppercase(current%metadataID)

       ! If the diagnostic (collection) name matches
       IF ( TRIM(thisDiagName) == TRIM(To_Uppercase(metadataID)) ) THEN

          !-----------------------------------------------------------------
          ! Create a new DgnTagList item.  This represents a single
          ! tag or wildcard belonging to a State_Diag diagnostic.
          !-----------------------------------------------------------------
          CALL Init_TagItem(                                                 &
               NewTagItem   = NewTagItem,                                    &
               name         = tagName,                                       &
               index        = index,                                         &
               RC           = RC                                            )

          ! Trap potential errors
          IF ( RC /= GC_SUCCESS ) THEN
             errMsg = 'Error encountered in "Init_TagItem!'
             CALL GC_Error( errMsg, RC, thisLoc )
             RETURN
          ENDIF

          IF ( isWildcard ) THEN

             !--------------------------------------------------------------
             ! Add the the DgnTagItem to the list of wildcards ...
             !--------------------------------------------------------------
             current%isWildcard = .TRUE.
             CALL InsertBeginning_TagList(                                   &
                  TagItem   = NewTagItem,                                    &
                  TagList   = current%wildcardList,                          &
                  RC        = RC                                            )

             ! Trap potential errors
             IF ( RC /= GC_SUCCESS ) THEN
                errMsg = &
                  'Error encountered in "InsertBeginning_tagList" (wildcard)!'
                CALL GC_Error( errMsg, RC, thisLoc )
                RETURN
             ENDIF

          ELSE

             !--------------------------------------------------------------
             ! ... or to the list of tags for the State_Diag diagnostic.
             !--------------------------------------------------------------
             current%isWildCard = .FALSE.
             CALL InsertBeginning_TagList(                                   &
                  TagItem   = NewTagItem,                                    &
                  TagList   = current%tagList,                               &
                  RC        = RC                                            )

             ! Trap potential errors
             IF ( RC /= GC_SUCCESS ) THEN
                errMsg = &
                  'Error encountered in "InsertBeginning_tagList" (tag)!'
                CALL GC_Error( errMsg, RC, thisLoc )
                RETURN
             ENDIF

          ENDIF
       ENDIF

       ! Go to the next item in DiagList
       current => current%next
    ENDDO

    ! Free pointers
    current => NULL()

  END SUBROUTINE Update_TaggedDiagList
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Cleanup_TagList
!
! !DESCRIPTION: Subroutine Cleanup\_TagList deallocates a TagList object
!  and all of its member objects including the linked list of TagItem objects.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Cleanup_TagList( TagList, RC )
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(DgnTagList), INTENT(INOUT) :: TagList
!
! !OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(OUT)   :: RC
!
! !REVISION HISTORY:
!  18 Sep 2019 - E. Lundgren - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Objects
    TYPE(DgnTagItem), POINTER :: current
    TYPE(DgnTagItem), POINTER :: next

    ! Strings
    CHARACTER(LEN=255)        :: diagId

    ! ================================================================
    ! Cleanup_TagList begins here
    ! ================================================================

    ! Initialize
    RC     = GC_SUCCESS
    diagId = 'taggeddiaglist_mod.F90:DgnTagItem'

    ! Deallocate each DgnTagItem in the DgnTagList, which is
    ! a list of tags or wildcards for each diagnostic
    current => TagList%head
    IF ( ASSOCIATED( current ) ) next => current%next

    DO WHILE ( ASSOCIATED( current ) )

       ! Free the DgnTagItem
       DEALLOCATE( current, STAT=RC )
       CALL GC_CheckVar( diagId, 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN

       ! Exit if we are at the end of the DgnTagList ...
       IF ( .NOT. ASSOCIATED ( next ) ) EXIT

       ! ...or if not, advance to the next DgnTagItem in the DgnTagList
       current => next
       next => current%next

    ENDDO

    ! Free pointers
    current => NULL()
    next    => NULL()

  END SUBROUTINE Cleanup_TagList
!EOC
END MODULE TaggedDiagList_Mod
