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
  PUBLIC  :: Init_TaggedDiagList
  PUBLIC  :: Print_TaggedDiagList
  PUBLIC  :: Query_TaggedDiagList
  PUBLIC  :: Query_Tag_in_TagList
  PUBLIC  :: Cleanup_TaggedDiagList
  PUBLIC  :: Print_TagList
!
! !PRIVATE MEMBER FUNCTIONS
!
  PRIVATE :: Init_TaggedDiagItem
  PRIVATE :: InsertBeginning_TaggedDiagList
  PRIVATE :: Update_TaggedDiagList
  PRIVATE :: Init_TagItem
  PRIVATE :: InsertBeginning_TagList
  PRIVATE :: Cleanup_TagList
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
  SUBROUTINE Init_TaggedDiagList( Input_Opt, DiagList, TaggedDiagList, RC )
!
! !USES:
!
    USE CharPak_Mod,   ONLY : To_UpperCase
    USE Input_Opt_Mod, ONLY : OptInput
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput),       INTENT(IN)    :: Input_Opt       ! Input Options
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
    LOGICAL                      :: TaggedDiagListExists
    INTEGER                      :: numTags
    INTEGER                      :: numWildCards
    INTEGER                      :: reverseIndex

    ! Strings
    CHARACTER(LEN=63)            :: tagName
    CHARACTER(LEN=63)            :: tagNameUpper
    CHARACTER(LEN=255)           :: ErrMsg
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
    thisLoc = '-> at Init_TaggedDiagList (located in module '            // &
              'Headers/taggeddiaglist_mod.F90)'

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
               Input_Opt      = Input_Opt,                                   &
               TaggedDiagList = TaggedDiagList,                              &
               name           = diagnostic%metadataID,                       &
               Found          = TaggedDiagListExists,                        &
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
          IF ( TaggedDiagListExists ) THEN

             !--------------------------------------------------------------
             ! If the TaggedDiagList (of wildcards or tags) already exists:
             ! (1) Add a new TaggedDiagItem into it
             ! (2) Update the index (which is actually in reverse order)
             !--------------------------------------------------------------
             reverseIndex = reverseIndex + 1
             CALL Update_TaggedDiagList(                                     &
                  Input_Opt         = Input_Opt,                             &
                  metadataID        = diagnostic%metadataID,                 &
                  isWildCard        = diagnostic%isWildCard,                 &
                  tagName           = tagName,                               &
                  index             = reverseIndex,                          &
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
             ! (2) Set the starting index of the TaggedDiagItem to 1.
             !     (NOTE: This index is actually in reverse order)
             !--------------------------------------------------------------
             reverseIndex = 1
             CALL Init_TaggedDiagItem(                                       &
                  Input_Opt         = Input_Opt,                             &
                  NewTaggedDiagItem = NewTaggedDiagItem,                     &
                  metadataID        = diagnostic%metadataID,                 &
                  isWildcard        = diagnostic%isWildCard,                 &
                  tagName           = tagName,                               &
                  index             = reverseIndex,                          &
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
                  Input_Opt         = Input_Opt,                             &
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
    TaggedDiagItem => NULL()
    current        => NULL()
    diagnostic     => DiagList%head

    DO WHILE ( ASSOCIATED( diagnostic ) )

       ! Only proceed for State_Diag diagnostics
       IF ( diagnostic%state == 'DIAG' .AND.                                 &
            ( diagnostic%isTagged .OR. diagnostic%isWildcard ) ) THEN

          !-----------------------------------------------------------------
          ! Get the number of tags and wildcards in the TaggedDiagList
          ! that corresponds to this State_Diag diagnostic
          !-----------------------------------------------------------------
          CALL Query_TaggedDiagList(                                         &
               Input_Opt      = Input_Opt,                                   &
               TaggedDiagList = TaggedDiagList,                              &
               name           = diagnostic%metadataID,                       &
               Found          = TaggedDiagListExists,                        &
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
          IF ( .not. TaggedDiagListExists ) CYCLE

          !-----------------------------------------------------------------
          ! Get the name of each tag or wildcard
          !-----------------------------------------------------------------
          IF ( diagnostic%isTagged ) THEN
             tagName = diagnostic%tag
          ELSE
             tagName = diagnostic%wildcard
          ENDIF

          ! Convert tag name to upper case
          tagNameUpper = To_UpperCase( tagName )

          !-----------------------------------------------------------------
          ! Loop over all tags of wildcards belonging to this diagnostic
          !-----------------------------------------------------------------
          TaggedDiagItem => TaggedDiagList%head

          DO WHILE ( ASSOCIATED( TaggedDiagItem ) )

             !--------------------------------------------------------------
             ! Reverse the indices of tags -- skip if there is only one
             !--------------------------------------------------------------
             IF ( numTags > 1 ) THEN
                current => TaggedDiagItem%tagList%head
                DO WHILE ( ASSOCIATED( current ) )
                   IF ( To_UpperCase(current%name) == tagNameUpper ) THEN
                      current%index = numtags - current%index + 1
                      EXIT
                   ENDIF
                   current => current%next
                ENDDO
                current => NULL()

             !--------------------------------------------------------------
             ! Reverse indices of wildcards -- skip if there is only one
             !--------------------------------------------------------------
             ELSE IF ( numWildCards > 1 ) THEN
                current => TaggedDiagItem%wildCardList%head
                DO WHILE ( ASSOCIATED( current ) )
                   IF ( To_UpperCase(current%name) == tagNameUpper ) THEN
                      current%index = numWildCards - current%index + 1
                      EXIT
                   ENDIF
                   current => current%next
                ENDDO
                current => NULL()
             ENDIF

             ! Advance to the next item in TaggedDiagList
             TaggedDiagItem => TaggedDiagItem%next

          ENDDO
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
! !IROUTINE: Print_TaggedDiagList
!
! !DESCRIPTION: Subroutine Print\_TaggedDiagList prints information for all
!  TaggedDiagItem members in a TaggedDiagList linked list.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Print_TaggedDiagList( Input_Opt, TaggedDiagList, RC )
!
! !USES:
!
    USE Input_Opt_Mod, ONLY : OptInput
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput),      INTENT(IN)  :: Input_Opt       ! Input Options object
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
    IF ( Input_Opt%amIRoot ) THEN
       WRITE( 6, 100 ) REPEAT( '=', 30 )
       WRITE( 6, 110 ) 'Summary of tagged diagnostics'

       ! Point to the first item in the TaggedDiagList
       current => TaggedDiagList%head

       ! Keep looping over all items in TaggedDiagList
       DO WHILE ( ASSOCIATED( current ) )

          ! Print name of diagnostic and if it is a wildcard
          WRITE( 6, 120 ) TRIM(current%metadataID)

          ! Print info about wildcards or tags
          IF ( current%isWildCard ) THEN
             WRITE( 6, 130 ) ADJUSTL('numWildcards:'),                       &
                             current%wildcardList%count
             CALL Print_TagList( Input_Opt, current%wildcardList, RC )
          ELSE
             WRITE( 6, 130 )  ADJUSTL('numTags:'), current%tagList%count
             CALL Print_TagList( Input_Opt, current%tagList, RC )
          ENDIF
          WRITE( 6, 120 ) ""

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
! !DESCRIPTION: Return information about a TaggedDiagList option, such
!  if name is found and information about that named object. Info includes
!  if is wildcard, a list of wildcards, number of wildcards, a list of
!  non-wildcard tags, or number of non-wildcard tags.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Query_TaggedDiagList( Input_Opt,    TaggedDiagList,             &
                                   name,         RC,                         &
                                   found,        isWildcard,                 &
                                   numWildcards, numTags,                    &
                                   wildcardList, tagList                    )
!
! !USES:
!
    USE Charpak_Mod,   ONLY : To_UpperCase
    USE Input_Opt_Mod, ONLY : OptInput
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput),      INTENT(IN)  :: Input_Opt
    TYPE(TaggedDgnList), INTENT(IN)  :: TaggedDiagList
    CHARACTER(LEN=*),    INTENT(IN)  :: name
!
! !OUTPUT PARAMETERS:
!
    INTEGER,             INTENT(OUT) :: RC
    LOGICAL,             OPTIONAL    :: found
    LOGICAL,             OPTIONAL    :: isWildcard
    INTEGER,             OPTIONAL    :: numWildcards
    INTEGER,             OPTIONAL    :: numTags
    TYPE(DgnTagList),    OPTIONAL    :: wildcardList
    TYPE(DgnTagList),    OPTIONAL    :: tagList
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
    CHARACTER(LEN=255)           :: thisDiagName

    ! Objects
    TYPE(TaggedDgnItem), POINTER :: current

    !=======================================================================
    ! Query_TaggedDiagList begins here!
    !=======================================================================

    ! Initialize
    RC    = GC_SUCCESS
    found = .FALSE.

    ! Search for name in list
    current => TaggedDiagList%head
    DO WHILE ( ASSOCIATED( current ) )
       thisDiagName = To_Uppercase(current%metadataID)
       IF ( TRIM(ThisDiagName) == TRIM(To_Uppercase(name)) ) THEN
          found = .TRUE.
          IF (PRESENT(isWildcard  )) isWildcard   = current%isWildcard
          IF (PRESENT(numWildcards)) numWildcards = current%wildcardList%count
          IF (PRESENT(numTags     )) numTags      = current%tagList%count
          IF (PRESENT(wildcardList)) wildcardList = current%wildcardList
          IF (PRESENT(tagList     )) tagList      = current%tagList
          current=> NULL()
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


    ! ================================================================
    ! Cleanup_TaggedDiagList begins here
    ! ================================================================

    ! Initialize
    RC      = GC_SUCCESS
    errMsg  = ''
    thisLoc = ' -> at Cleanup_TaggedDiagList (located in module '          // &
              'Headers/taggeddiaglist_mod.F90)'

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
  SUBROUTINE Init_TaggedDiagItem( Input_Opt,  NewTaggedDiagItem,             &
                                  metadataID, isWildcard,                    &
                                  tagName,    index,                         &
                                  RC                                        )
!
! !USES:
!
    USE Input_Opt_Mod, ONLY : OptInput
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput),      INTENT(IN)  :: Input_Opt    ! Input Options object
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

    ! ======================================================================
    ! Init_TaggedDiagList begins here
    ! ======================================================================

    ! Initialize
    RC      = GC_SUCCESS
    errMsg  = ''
    thisLoc = ' ->at Init_TaggedDiagItem (located in module '             // &
              'Headers/taggeddiaglist_mod.F90)'

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
         Input_Opt    = Input_Opt,                                           &
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
            Input_Opt = Input_Opt,                                           &
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
            Input_Opt = Input_Opt,                                           &
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
  SUBROUTINE Init_TagItem( Input_Opt, NewTagItem, name, index, RC )
!
! !USES:
!
    USE Input_Opt_Mod, ONLY : OptInput
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput),   INTENT(IN)  :: Input_Opt     ! Input Options object
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
  SUBROUTINE InsertBeginning_TaggedDiagList( Input_Opt,      TaggedDiagItem, &
                                             TaggedDiagList, RC             )
!
! USES:
!
    USE Input_Opt_Mod, ONLY : OptInput
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput),      INTENT(IN)    :: Input_Opt
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
  SUBROUTINE InsertBeginning_TagList( Input_Opt, TagItem, TagList, RC )
!
! !USES:
!
    USE Input_Opt_Mod, ONLY : OptInput
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput),   INTENT(IN)    :: Input_Opt   ! Input Options object
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
  SUBROUTINE Query_Tag_in_TagList( Input_Opt, TagList, name, found, RC )
!
! !USES:
!
    USE Charpak_Mod,   ONLY : To_UpperCase
    USE Input_Opt_Mod, ONLY : OptInput
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput),   INTENT(IN)  :: Input_Opt   ! Input Options object
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
  SUBROUTINE Print_TagList( Input_Opt, TagList, RC )
!
! !USES:
!
    USE Input_Opt_Mod, ONLY : OptInput
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput),   INTENT(IN)    :: Input_Opt  ! Input Options object
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
    RC      = GC_SUCCESS
    current => TagList%head

    ! Only print taglist if we are on the root core
    IF ( Input_Opt%amIRoot ) THEN
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
  SUBROUTINE Update_TaggedDiagList( Input_Opt, metadataID, isWildcard,       &
                                    tagName,   index,      TaggedDiagList,   &
                                    RC                                      )
!
! !USES:
!
    USE Charpak_Mod,   ONLY : To_UpperCase
    USE Input_Opt_Mod, ONLY : OptInput
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput),      INTENT(IN)    :: Input_Opt      ! Input Options object
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
    thisLoc = ' -> at Update_TaggedDiagItem (located in module '          // &
              'taggeddiaglist_mod.F90)'

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
               Input_Opt    = Input_Opt,                                     &
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
                  Input_Opt = Input_Opt,                                     &
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
                  Input_Opt = Input_Opt,                                     &
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
