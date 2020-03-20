!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: taggeddiaglist_mod.F90
!
! !DESCRIPTION: Module taggeddiaglist\_mod
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
!
! !PRIVATE MEMBER FUNCTIONS
!
  PRIVATE :: Init_TaggedDiagItem
  PRIVATE :: InsertBeginning_TaggedDiagList
  PRIVATE :: Update_TaggedDiagList
  PRIVATE :: Init_TagItem
  PRIVATE :: InsertBeginning_TagList
  PRIVATE :: Print_TagList
  PRIVATE :: Cleanup_TagList
!
! !PUBLIC DATA TYPES:
!
  !=========================================================================
  ! Derived type for tag items (can be used for any linked list of strings)
  !=========================================================================
  TYPE, PUBLIC :: DgnTagItem
     CHARACTER(LEN=63)         :: name
     TYPE(DgnTagItem), POINTER :: next
  END TYPE DgnTagItem

  !=========================================================================
  ! Derived type for tag list (can be used for any linked list of strings)
  !=========================================================================
  TYPE, PUBLIC :: DgnTagList
     INTEGER                    :: count
     TYPE(DgnTagItem), POINTER  :: head
  END TYPE DgnTagList

  !=========================================================================
  ! Derived type for tagged diagnostic items, e.g. DryDep
  !=========================================================================
  TYPE, PUBLIC :: TaggedDgnItem
     CHARACTER(LEN=63)             :: metadataID
     LOGICAL                       :: isWildcard
     TYPE(DgnTagList)              :: tagList
     TYPE(DgnTagList)              :: wildcardList
     TYPE(TaggedDgnItem), POINTER  :: next
  END TYPE TaggedDgnItem

  !=========================================================================
  ! Derived type for tagged diagnostic list, e.g. DryDep, SpeciesConc, etc
  !=========================================================================
  TYPE, PUBLIC :: TaggedDgnList
     TYPE(TaggedDgnItem), POINTER  :: head
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
! !DESCRIPTION: Subroutine Init\_TaggedDiagList
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Init_TaggedDiagList( am_I_Root, DiagList, TaggedDiagList, RC )
!
! !INPUT PARAMETERS:
!
    LOGICAL,              INTENT(IN)    :: am_I_Root
    TYPE(DgnList),        INTENT(IN)    :: DiagList
!
! !INPUT AND OUTPUT PARAMETERS:
!
    TYPE(TaggedDgnList),  INTENT(INOUT) :: TaggedDiagList
!
! !OUTPUT PARAMETERS:
!
    INTEGER,              INTENT(OUT)   :: RC
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
    LOGICAL                      :: FOUND
    TYPE(DgnItem),       POINTER :: current
    TYPE(TaggedDgnItem), POINTER :: NewTaggedDiagItem
    CHARACTER(LEN=63)            :: tagName

    !=======================================================================
    ! Init_TaggedDiagList begins here
    !=======================================================================

    ! Initialize list of tagged diagnostics
    TaggedDiagList%head => NULL()

    ! Loop over all diagnostics listed in diagnostics config file
    current => DiagList%head
    DO WHILE ( ASSOCIATED( current ) )
       IF ( current%state == 'DIAG' .AND. &
            ( current%isTagged .OR. current%isWildcard ) ) THEN
          CALL Query_TaggedDiagList( am_I_Root, TaggedDiagList, &
                                     name=current%metadataID,   &
                                     FOUND=FOUND,               &
                                     RC=RC)
          If ( current%isTagged ) THEN
             tagName = current%tag
          ELSE
             tagName = current%wildcard
          ENDIF
          IF ( FOUND ) THEN
             CALL Update_TaggedDiagList( am_I_Root,                 &
                                         current%metadataID,        &
                                         current%isWildcard,        &
                                         tagName,                   &
                                         TaggedDiagList,            &
                                         RC)
          ELSE
             CALL Init_TaggedDiagItem( am_I_Root,                     &
                                       NewTaggedDiagItem,             &
                                       metadataID=current%metadataID, &
                                       isWildcard=current%isWildcard, &
                                       tagName=tagName,               &
                                       RC=RC)
             CALL InsertBeginning_TaggedDiagList( am_I_Root,         &
                                                  NewTaggedDiagItem, &
                                                  TaggedDiagList, RC )
          ENDIF
       ENDIF
       current => current%next
    ENDDO
    current => NULL()

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
  SUBROUTINE Print_TaggedDiagList( am_I_Root, TaggedDiagList, RC )
!
! !INPUT PARAMETERS:
!
    LOGICAL,             INTENT(IN)    :: am_I_Root
    TYPE(TaggedDgnList), INTENT(IN)    :: TaggedDiagList
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,             INTENT(INOUT) :: RC
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
    CHARACTER(LEN=255)           :: thisLoc

    ! ================================================================
    ! Print_TaggedDiagList begins here
    ! ================================================================
    thisLoc = 'Print_TaggedDiagList (taggeddiaglist_mod.F90)'
    current => TaggedDiagList%head

    IF ( am_I_Root ) THEN
       PRINT *, " "
       PRINT *, "============================="
       PRINT *, "Summary of tagged diagnostics"
       PRINT *, " "
    ENDIF
    DO WHILE ( ASSOCIATED( current ) )
       IF ( am_I_Root ) THEN
          PRINT *, TRIM(current%metadataID)
          PRINT '(A15,L3)', ADJUSTL('isWildcard:'), current%isWildcard
          PRINT '(A15,I3)', ADJUSTL('numWildcards:'), &
                              current%wildcardList%count
          CALL Print_TagList( am_I_Root, current%wildcardList, RC )
          PRINT '(A15,I3)', ADJUSTL('numTags:'), current%tagList%count
          CALL Print_TagList( am_I_Root, current%tagList, RC )
          PRINT *, " "
       ENDIF
       current => current%next
    ENDDO
    current => NULL()
    IF ( am_I_Root ) PRINT *, " "

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
  SUBROUTINE Query_TaggedDiagList ( am_I_Root, TaggedDiagList, name, &
                                    found, isWildcard, numWildcards, &
                                    numTags, wildcardList, tagList, RC )
!
! !USES:
!
    USE Charpak_Mod, ONLY : To_UpperCase
!
! !INPUT PARAMETERS:
!
    LOGICAL,             INTENT(IN) :: am_I_Root
    TYPE(TaggedDgnList), INTENT(IN) :: TaggedDiagList
    CHARACTER(LEN=*),    INTENT(IN) :: name
!
! !OUTPUT PARAMETERS:
!
    LOGICAL,             OPTIONAL   :: found
    LOGICAL,             OPTIONAL   :: isWildcard
    INTEGER,             OPTIONAL   :: numWildcards
    INTEGER,             OPTIONAL   :: numTags
    TYPE(DgnTagList),    OPTIONAL   :: wildcardList
    TYPE(DgnTagList),    OPTIONAL   :: tagList
    INTEGER,             OPTIONAL   :: RC
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
    CHARACTER(LEN=255)           :: thisLoc, thisDiagName

    ! Initialize
    thisLoc = 'Query_Name_in_TaggedDiagList (taggeddiaglist_mod.F90)'
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
    TYPE(TaggedDgnItem), POINTER :: current
    TYPE(TaggedDgnItem), POINTER :: next
    CHARACTER(LEN=255)           :: thisLoc

    ! ================================================================
    ! Cleanup_TaggedDiagList begins here
    ! ================================================================
    thisLoc = 'Cleanup_TaggedDiagList (taggeddiaglist_mod.F90)'

    ! Deallocate each item in the linked list of DiagExport objects
    current => TaggedDiagList%head
    IF ( ASSOCIATED( current ) ) next => current%next
    DO WHILE ( ASSOCIATED( current ) )
       CALL Cleanup_TagList( current%taglist, RC )
       CALL Cleanup_TagList( current%wildcardlist, RC )
       DEALLOCATE( current, STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN
          PRINT *, "Error in ", trim(thisLoc)
          RETURN
       ENDIF
       IF ( .NOT. ASSOCIATED ( next ) ) EXIT
       current => next
       next => current%next
    ENDDO
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
! !DESCRIPTION: Initializes a taggedDiagItem object
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Init_TaggedDiagItem ( am_I_Root, NewTaggedDiagItem, metadataID, &
                                   isWildcard, tagName, RC  )
!
! !INPUT PARAMETERS:
!
    LOGICAL,             INTENT(IN) :: am_I_Root
    CHARACTER(LEN=*),    OPTIONAL   :: metadataID
    LOGICAL,             OPTIONAL   :: isWildcard
    CHARACTER(LEN=*),    OPTIONAL   :: tagName
!
! !OUTPUT PARAMETERS:
!
    TYPE(TaggedDgnItem), POINTER    :: NewTaggedDiagItem
    INTEGER,             OPTIONAL   :: RC
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
    CHARACTER(LEN=255)        :: thisLoc
    TYPE(DgnTagItem), POINTER :: NewTagItem

    ! ================================================================
    ! Init_TaggedDiagList begins here
    ! ================================================================
    thisLoc = 'Init_TaggedDiagItem (taggeddiaglist_mod.F90)'

    ALLOCATE(NewTaggedDiagItem)
    NewTaggedDiagItem%metadataID = ' '
    NewTaggedDiagItem%isWildcard = .FALSE.

    ! Initialize wildcard list
    NewTaggedDiagItem%wildcardList%head => NULL()
    NewTaggedDiagItem%wildcardList%count = 0

    ! Initialize non-wildcard list of tags
    NewTaggedDiagItem%tagList%head => NULL()
    NewTaggedDiagItem%tagList%count = 0

    ! Set true values if passed
    IF (PRESENT(metadataID) ) NewTaggedDiagItem%metadataID = TRIM(metadataID)
    IF (PRESENT(isWildcard) ) NewTaggedDiagItem%isWildcard = isWildcard
    IF (PRESENT(tagName) .AND. PRESENT(isWildcard) ) THEN
       CALL Init_TagItem( am_I_Root, NewTagItem, tagName, RC )
       IF ( isWildcard ) THEN
          CALL InsertBeginning_TagList( am_I_Root, NewTagItem, &
                                        NewTaggedDiagItem%wildcardList, RC )
       ELSE
          CALL InsertBeginning_TagList( am_I_Root, NewTagItem, &
                                        NewTaggedDiagItem%tagList, RC )
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
! !DESCRIPTION: Initializes a TagItem object, which contains information
!  about a single GEOS-Chem diagnostic tag.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Init_TagItem ( am_I_Root, NewTagItem, name, RC  )
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN) :: am_I_Root
    TYPE(DgnTagItem), POINTER    :: NewTagItem
    CHARACTER(LEN=*)             :: name
!
! !OUTPUT PARAMETERS:
!
    INTEGER,          OPTIONAL   :: RC
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
    CHARACTER(LEN=255) :: thisLoc

    ! ================================================================
    ! Init_TagItem begins here
    ! ================================================================
    thisLoc = 'Init_TagItem (taggeddiaglist_mod.F90)'

    ALLOCATE(NewTagItem)
    NewTagItem%name      = TRIM(name)

  END SUBROUTINE Init_TagItem
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: InsertBeginning_TaggedDiagList
!
! !DESCRIPTION: Inserts a new node at the beginning of the TaggedDiagList linked
!  list object.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE InsertBeginning_TaggedDiagList ( am_I_Root, TaggedDiagItem, &
                                              TaggedDiagList, RC )
!
! !INPUT PARAMETERS:
!
    LOGICAL,             INTENT(IN)    :: am_I_Root
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
    TYPE(TaggedDgnItem), POINTER :: NewTaggedDiagItem
    CHARACTER(LEN=255)           :: thisLoc

    ! ================================================================
    ! InsertBeginning_TaggedDiagList begins here
    ! ================================================================
    thisLoc = 'InsertBeginning_TaggedDiagList (taggeddiaglist_mod.F90)'

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
  SUBROUTINE InsertBeginning_TagList ( am_I_Root, TagItem, TagList, RC )
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN)    :: am_I_Root
    TYPE(DgnTagItem), POINTER       :: TagItem
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
!  18 Nov 2019 - E. Lundgren - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    TYPE(DgnTagItem), POINTER :: NewTagItem
    CHARACTER(LEN=255)        :: thisLoc

    ! ================================================================
    ! InsertBeginning_TagList begins here
    ! ================================================================
    thisLoc = 'InsertBeginning_TagList (taggeddiaglist_mod.F90)'

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
! !DESCRIPTION: Searches for a given tag name within TagList object.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Query_Tag_in_TagList( am_I_Root, TagList, name, found, RC )
!
! !USES:
!
    USE Charpak_Mod, ONLY : To_UpperCase
!
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN)  :: am_I_Root
    TYPE(DgnTagList), INTENT(IN)  :: TagList
    CHARACTER(LEN=*), INTENT(IN)  :: name
!
! !OUTPUT PARAMETERS:
!
    LOGICAL,          INTENT(OUT) :: found
    INTEGER,          INTENT(OUT) :: RC
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
    CHARACTER(LEN=255)        :: thisLoc, thisTagName

    ! Initialize
    thisLoc = 'Query_Tag_in_TagList (taggeddiaglist_mod.F90)'
    found = .FALSE.

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
!  TagItem members in a TagList linked list.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Print_TagList( am_I_Root, TagList, RC )
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN)    :: am_I_Root
    TYPE(DgnTagList), INTENT(IN)    :: TagList
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(INOUT) :: RC
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
    thisLoc = 'Print_TagList (taggeddiaglist_mod.F90)'
    current => TagList%head

    DO WHILE ( ASSOCIATED( current ) )
       IF ( am_I_Root ) THEN
          PRINT *, '                  ', TRIM(current%name)
       ENDIF
       current => current%next
    ENDDO
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
! !DESCRIPTION: Update the taggedDiagList object with new information about
!  an item
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Update_TaggedDiagList ( am_I_Root, metadataID, isWildcard, &
                                     tagName, TaggedDiagList, RC  )
!
! !USES:
!
    USE Charpak_Mod, ONLY : To_UpperCase
!
! !INPUT PARAMETERS:
!
    LOGICAL,             INTENT(IN)    :: am_I_Root
    CHARACTER(LEN=*),    INTENT(IN)    :: metadataID
    LOGICAL,             INTENT(IN)    :: isWildcard
    CHARACTER(LEN=*),    INTENT(IN)    :: tagName
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
    TYPE(TaggedDgnItem), POINTER :: current
    TYPE(DgnTagItem),    POINTER :: NewTagItem
    CHARACTER(LEN=255)           :: thisLoc, thisDiagName

    ! ================================================================
    ! Update_TaggedDiagList begins here
    ! ================================================================
    thisLoc = 'Update_TaggedDiagItem (taggeddiaglist_mod.F90)'

    ! Search for item in list and update if found
    current => TaggedDiagList%head
    DO WHILE ( ASSOCIATED( current ) )
       thisDiagName = To_Uppercase(current%metadataID)
       IF ( TRIM(ThisDiagName) == TRIM(To_Uppercase(metadataID)) ) THEN
          CALL Init_TagItem( am_I_Root, NewTagItem, tagName, RC )
          IF ( isWildcard ) THEN
             current%isWildcard = .TRUE.
             CALL InsertBeginning_TagList( am_I_Root, NewTagItem, &
                                           current%wildcardList, RC )
          ELSE
             CALL InsertBeginning_TagList( am_I_Root, NewTagItem, &
                                           current%tagList, RC )
          ENDIF
       ENDIF
       current => current%next
    ENDDO
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
    TYPE(DgnTagItem), POINTER :: current
    TYPE(DgnTagItem), POINTER :: next
    CHARACTER(LEN=255)        :: thisLoc

    ! ================================================================
    ! Cleanup_TagList begins here
    ! ================================================================
    thisLoc = 'Cleanup_TagList (taggeddiaglist_mod.F90)'

    ! Deallocate each item in the linked list of collection objects
    current => TagList%head
    IF ( ASSOCIATED( current ) ) next => current%next
    DO WHILE ( ASSOCIATED( current ) )
       DEALLOCATE( current, STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN
          PRINT *, "Error in ", trim(thisLoc)
          RETURN
       ENDIF
       IF ( .NOT. ASSOCIATED ( next ) ) EXIT
       current => next
       next => current%next
    ENDDO
    current => NULL()
    next    => NULL()

  END SUBROUTINE Cleanup_TagList
!EOC
END MODULE TaggedDiagList_Mod
