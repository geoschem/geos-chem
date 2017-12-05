!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: diagnostics_mod.F90
!
! !DESCRIPTION: Module diagnostics\_mod.F90 contains the derived types
!  and subroutines used for reading and storing user-configured diagnostic
!  information from the history configuration file. The diagnostics list is
!  used to allocate diagnostics stored in container State_Diag and to 
!  declare exports in GCHP. It does not store collection information.
!
! TODO: Also read input.geos to get the wavelengths in the radiation menu
!       Store in module-level variables RadWL1, RadWL2, and RadWL3 (strings).
!       If wavelength not present, store as 'WL2' etc. This must be done
!       during init_diaglist. Will use these in state_diag_mod, in both
!       init_state_diag and in get_metdata_state_diag.
! 
! !INTERFACE:
!
MODULE Diagnostics_Mod
!
! !USES:
!
  USE ErrCode_Mod
  USE Precision_Mod

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: Init_DiagList
  PUBLIC :: Print_DiagList
  PUBLIC :: Check_DiagList
  PUBLIC :: Cleanup_DiagList
!
! !PRIVATE MEMBER FUNCTIONS
!
  PRIVATE :: Init_DiagItem
  PRIVATE :: InsertBeginning_DiagList
  PRIVATE :: Search_DiagList
!
! !PUBLIC DATA TYPES:
!
  !=========================================================================
  ! Derived type for Diagnostics List
  !=========================================================================
  TYPE, PUBLIC :: DgnList
     TYPE(DgnItem), POINTER  :: head
  END TYPE DgnList

  !=========================================================================
  ! Derived type for Diagnostics Item (unique name in HISTORY.rc)
  !=========================================================================
  TYPE, PUBLIC :: DgnItem
     CHARACTER(LEN=63)      :: name 
     CHARACTER(LEN=63)      :: state
     CHARACTER(LEN=63)      :: metadataID
     CHARACTER(LEN=63)      :: registryID
     LOGICAL                :: isWildcard
     CHARACTER(LEN=7)       :: wildcard
     LOGICAL                :: isTagged
     CHARACTER(LEN=63)      :: tag
     TYPE(DgnItem), POINTER :: next
  END TYPE DgnItem

  !=========================================================================
  ! Configurable Settings Used for Diagnostic Names at Run-time
  !=========================================================================
  CHARACTER(LEN=5), PUBLIC :: RadWL(3) ! Wavelengths configured in rad menu
!
! !REVISION HISTORY:
!  22 Sep 2017 - E. Lundgren - Initial version
!  01 Nov 2017 - R. Yantosca - Moved ReadOneLine, CleanText to charpak_mod.F90
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
! !IROUTINE: Init_DiagList 
!
! !DESCRIPTION: 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Init_DiagList ( am_I_Root, historyConfigFile, DiagList, RC )
!
! !USES:
!
    USE Charpak_Mod
    USE InquireMod,       ONLY: findFreeLun
!
! !INPUT PARAMETERS:
!
    LOGICAL,              INTENT(IN)    :: am_I_Root
    CHARACTER(LEN=*),     INTENT(IN)    :: historyConfigFile
!
! !INPUT AND OUTPUT PARAMETERS:
!
    TYPE(DgnList),        INTENT(INOUT) :: DiagList
!
! !OUTPUT PARAMETERS:
!
    INTEGER,              INTENT(OUT)   :: RC
!
! !REVISION HISTORY:
!  22 Sep 2017 - E. Lundgren - initial version
!  28 Nov 2017 - C. J. Lee   - Allow state MET variables to have underscores
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER                 :: fId, IOS, N, I, J
    INTEGER                 :: numSpcWords, numIDWords
    CHARACTER(LEN=255)      :: errMsg, thisLoc, nameAllCaps
    CHARACTER(LEN=255)      :: line, SubStrs(500), SubStr
    CHARACTER(LEN=255)      :: wildcard, tag, name, state
    CHARACTER(LEN=255)      :: metadataID, registryID
    LOGICAL                 :: EOF, found, isWildcard, isTagged
    TYPE(DgnItem),  POINTER :: NewDiagItem

    ! ================================================================
    ! Init_DiagList begins here
    ! ================================================================
    thisLoc = 'Init_DiagList (diagnostics_mod.F90)'

    ! Init
    EOF = .FALSE.
    found = .FALSE.
    NewDiagItem => NULL()
    RadWL(:) = ['WL1  ','WL2  ','WL3  ']

    ! Create DiagList object
    DiagList%head => NULL()

    ! Open the history config file
    fId = FindFreeLun()
    OPEN( fId, FILE=TRIM(historyConfigFile), STATUS='OLD', IOSTAT=RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Could not open "' //TRIM(historyConfigFile) // '"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    ! Read data from the file
    DO
    
       ! Read line and strip leading/trailing spaces
       Line = ReadOneLine( fId, EOF, IOS, Squeeze=.TRUE. )
       IF ( EOF ) EXIT
       IF ( IOS > 0 ) THEN
          ErrMsg = 'Unexpected end-of-file in "'       // &
                    TRIM( historyConfigFile ) // '"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       ! Skip line if GIGCchem not present
       IF ( INDEX( line, 'GIGCchem' ) .le. 0 ) CYCLE
    
       ! Remove commas, spaces, and tabs
       Line = CleanText( Line )    

       ! Skip if the line is commented out
       IF ( Line(1:1) == "#"  ) CYCLE

       ! Get the item name
       CALL StrSplit( Line, " ", SubStrs, N )
       IF ( INDEX(line, '.fields') > 0 .AND. N > 1 ) THEN
          name = TRIM(SubStrs(2))
       ELSE
          name = TRIM(SubStrs(1))
       ENDIF
       nameAllCaps = To_Uppercase( TRIM(name) )

       ! Set GC state
       IF ( nameAllCaps(1:4) == 'MET_' ) THEN
          state = 'MET'
       ELSEIF ( nameAllCaps(1:5) == 'CHEM_' ) THEN
          state = 'CHEM'
#if defined( ESMF_ )
       ! Emissions diagnostics are included in HISTORY.rc in GCHP only
       ELSEIF ( nameAllCaps(1:4) == 'EMIS' ) THEN
          state = 'EMISSIONS'
       ELSEIF ( nameAllCaps(1:4) == 'SPC_' ) THEN
          state = 'INTERNAL'
#endif
       ELSE
          state = 'DIAG'
       ENDIF

       ! Get wildcard, if any 
       ! NOTE: Must be prefaced with single underscore in HISTORY.rc!
       isWildcard = .FALSE.
       wildcard   = ''
       IF ( INDEX( name, '?' ) > 0 ) THEN
#if defined( ESMF_ )
          ! Exit with an error if using GCHP and wildcard is present
          CALL GC_Error( 'Wildcards not allowed in GCHP', RC, ThisLoc )
          RETURN
#endif
          isWildcard = .TRUE.
          CALL StrSplit( name, '?', SubStrs, N )
          wildcard = SubStrs(N-1)
       ENDIF

       ! Get tag, if any
       isTagged  = .FALSE.
       tag = ''
       IF ( .NOT. isWildcard ) THEN
          CALL StrSplit( name, '_', SubStrs, N )
          IF ( TRIM(state) == 'DIAG' .AND. N == 2 ) THEN
             isTagged = .TRUE.
             tag = SubStrs(2)
          ELSEIF ( TRIM(state) == 'CHEM' &
                   .AND. N == 3 ) THEN
             isTagged = .TRUE.
             tag = SubStrs(3)
          ENDIF
       ENDIF

       ! Get registryID - start with the full name in HISTORY.rc
       registryID = TRIM(nameAllCaps)
       ! Strip off the state prefix, if any
       IF ( TRIM(state) == 'MET' ) THEN
          registryID = registryID(5:)
       ELSE IF ( TRIM(state) == 'CHEM' ) THEN
          registryID = registryID(6:)
       ENDIF
       ! Strip off the wildcard, if any
       IF ( isWildcard ) THEN
          I = INDEX( TRIM(registryID), '_' ) 
          IF ( I .le. 0 ) THEN
             ErrMsg = 'Error setting registryID. Single underscore must' &
                      // ' precede wildcard in HISTORY.rc!'
             CALL GC_ERROR( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF
          registryID = registryID(1:I-1)
       ENDIF

       ! Get metadataID - start with the registry ID
       metadataID = registryID
       ! Strip off the tag suffix, if any
       IF ( isTagged ) THEN
          I = INDEX( TRIM(metadataID), '_' ) 
          metadataID = metadataID(1:I-1)
       ENDIF
       
       ! Skip if already encountered entry with same full name
       ! TODO: make this into a function that returns true or false
       ! Have a series of functions that do quality checks too
       CALL Search_DiagList( am_I_Root, DiagList, name, found, RC )
       IF ( found ) CYCLE

       ! Create a new DiagItem object
       CALL Init_DiagItem( am_I_Root, NewDiagItem, &
                           name=name,              &
                           state=state,            &
                           metadataID=metadataID,  &
                           registryID=registryID,  &
                           isWildcard=isWildcard,  &
                           wildcard=wildcard,      &
                           isTagged=isTagged,     &
                           tag=tag,        &
                           RC=RC  )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Error initializing DiagItem ' // TRIM(name)
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    
       ! Add new DiagItem to linked list
       CALL InsertBeginning_DiagList( am_I_Root, NewDiagItem, DiagList, RC )
       IF ( RC /= GC_SUCCESS ) RETURN

    ENDDO
    
    ! Close the file
    CLOSE( fId )

! ewl new:
    ! Open the input.geos config file
    fId = FindFreeLun()
    OPEN( fId, FILE=TRIM('input.geos'), STATUS='OLD', IOSTAT=RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Could not open "' //TRIM('input.geos') // '"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    
    ! Read data from the file
    DO    
       ! Read line and strip leading/trailing spaces
       Line = ReadOneLine( fId, EOF, IOS, Squeeze=.TRUE. )
       IF ( EOF ) EXIT
       IF ( IOS > 0 ) THEN
          ErrMsg = 'Unexpected end-of-file in "'       // &
                    TRIM('input.geos') // '"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    
       ! Skip line if not relevant
       IF ( INDEX( Line, 'AOD Wavelength' ) .le. 0 ) CYCLE
    
       ! Update wavelength(s) with string in file
       I = INDEX( Line, ':' )
       CALL StrSplit( Line(I:), ' ', SubStrs, N )
       DO J = 1, N-1
          WRITE ( RadWL(J), "(a5)" ) SubStrs(J+1)
       ENDDO
    
       ! End the loop
       EXIT
    ENDDO
    
    ! Close the file
    CLOSE( fId )

  END SUBROUTINE Init_DiagList
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Init_DiagItem
!
! !DESCRIPTION: 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Init_DiagItem ( am_I_Root,  NewDiagItem, name,       state,     &
                             metadataID, registryID,  isWildcard, wildcard,  &
                             isTagged,   tag,         RC  )
!
! !INPUT PARAMETERS:
!
    LOGICAL,             INTENT(IN) :: am_I_Root
    CHARACTER(LEN=*),    OPTIONAL   :: name
    CHARACTER(LEN=*),    OPTIONAL   :: state
    CHARACTER(LEN=*),    OPTIONAL   :: metadataID
    CHARACTER(LEN=*),    OPTIONAL   :: registryID
    LOGICAL,             OPTIONAL   :: isWildcard
    CHARACTER(LEN=*),    OPTIONAL   :: wildcard
    LOGICAL,             OPTIONAL   :: isTagged
    CHARACTER(LEN=*),    OPTIONAL   :: tag
!
! !OUTPUT PARAMETERS:
!
    TYPE(DgnItem), POINTER  :: NewDiagItem
    INTEGER,       OPTIONAL :: RC
!
! !REVISION HISTORY:
!  21 Sep 2017 - E. Lundgren - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    CHARACTER(LEN=255) :: thisLoc

    ! ================================================================
    ! Init_DiagList begins here
    ! ================================================================
    thisLoc = 'Init_DiagItem (diagnostics_mod.F90)'

    ALLOCATE(NewDiagItem)
    NewDiagItem%name       = TRIM(name)
    NewDiagItem%state      = TRIM(state)
    NewDiagItem%metadataID = TRIM(metadataID)
    NewDiagItem%registryID = TRIM(registryID)
    NewDiagItem%isWildcard = isWildcard
    NewDiagItem%wildcard   = TRIM(wildcard)
    NewDiagItem%isTagged   = isTagged
    NewDiagItem%tag        = TRIM(tag)

  END SUBROUTINE Init_DiagItem
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: InsertBeginning_DiagList 
!
! !DESCRIPTION: 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE InsertBeginning_DiagList ( am_I_Root, DiagItem, DiagList, RC )
!
! !USES:
!
!
! !INPUT PARAMETERS:
!
    LOGICAL,         INTENT(IN)    :: am_I_Root
    TYPE(DgnItem),   POINTER       :: DiagItem
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(DgnList),   INTENT(INOUT) :: DiagList
!
! !OUTPUT PARAMETERS:
!
    INTEGER,         INTENT(OUT)   :: RC
!
! !REVISION HISTORY:
!  22 Sep 2017 - E. Lundgren - initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    TYPE(DgnItem),  POINTER :: NewDiagItem
    CHARACTER(LEN=255)      :: thisLoc

    ! ================================================================
    ! InsertBeginning_DiagList begins here
    ! ================================================================
    thisLoc = 'InsertBeginning_DiagList (diagnostics_mod.F90)'

    ! Add new object to the beginning of the linked list
    DiagItem%next => DiagList%head
    DiagList%head => DiagItem

  END SUBROUTINE InsertBeginning_DiagList
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Search_DiagList
!
! !DESCRIPTION:
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Search_DiagList ( am_I_Root, DiagList, name, found, RC )
!
! !INPUT PARAMETERS:
!
    LOGICAL,           INTENT(IN) :: am_I_Root
    TYPE(DgnList),     INTENT(IN) :: DiagList
    CHARACTER(LEN=*),  INTENT(IN) :: name
!
! !OUTPUT PARAMETERS:
!
    LOGICAL,           INTENT(OUT) :: found
    INTEGER,           INTENT(OUT) :: RC 
!
! !REVISION HISTORY:
!  22 Sep 2017 - E. Lundgren - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    TYPE(DgnItem), POINTER :: current
    CHARACTER(LEN=255)     :: thisLoc

    ! Initialize
    thisLoc = 'Search_DiagList (diagnostics_mod.F90)'
    found = .FALSE.

    ! Search for name in list
    current => DiagList%head
    DO WHILE ( ASSOCIATED( current ) )
       IF ( current%name == name ) THEN
          found = .TRUE.
          EXIT
       ENDIF
       current => current%next    
    ENDDO
    current => NULL()

  END SUBROUTINE Search_DiagList
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Check_DiagList
!
! !DESCRIPTION:
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Check_DiagList( am_I_Root, DiagList, substr, found, RC )
!
! !USES:
!
    USE Charpak_Mod, ONLY : To_UpperCase
!
! !INPUT PARAMETERS:
!
    LOGICAL,           INTENT(IN)  :: am_I_Root   ! Are we on the root CPU?
    TYPE(DgnList),     INTENT(IN)  :: DiagList    ! Diagnostic list object
    CHARACTER(LEN=*),  INTENT(IN)  :: substr      ! Substring
!
! !OUTPUT PARAMETERS:
!
    LOGICAL,           INTENT(OUT) :: found       ! Was a match found (T/F)?
    INTEGER,           INTENT(OUT) :: RC          ! Success or failure?
!
! !REVISION HISTORY:
!  22 Sep 2017 - E. Lundgren - Initial version
!  01 Nov 2017 - R. Yantosca - Now use To_UpperCase from charpak_mod.F90
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    TYPE(DgnItem), POINTER :: current
    CHARACTER(LEN=255)     :: thisLoc, substr_AllCaps, currentName_AllCaps

    ! Initialize
    thisLoc = 'Check_DiagList (diagnostics_mod.F90)'
    found = .FALSE.

    ! Convert strings to uppercase for comparison
    substr_AllCaps = To_Uppercase( TRIM( substr ) )

    ! Search for name in list
    current => DiagList%head
    DO WHILE ( ASSOCIATED( current ) )
       currentName_AllCaps = To_Uppercase( current%name )
!------------------------------------------------------------------------------
! Prior to 11/16/17:
! Don't just restrict to State_Diag (bmy, 11/16/17)
!       IF ( ( current%state == 'DIAG' ) .AND. &
!            INDEX( currentName_AllCaps, TRIM( substr_AllCaps ) ) > 0 ) THEN
!------------------------------------------------------------------------------
       IF ( INDEX( currentName_AllCaps, TRIM( substr_AllCaps ) ) > 0 ) THEN
          found = .TRUE.
          EXIT
       ENDIF
       current => current%next    
    ENDDO
    current => NULL()

  END SUBROUTINE Check_DiagList
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Print_DiagList 
!
! !DESCRIPTION: Subroutine Print_DiagList prints information for all
!  DiagItem members in a DiagList linked list.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Print_DiagList( am_I_Root, DiagList, RC )
!
! !INPUT PARAMETERS:
!
    LOGICAL,           INTENT(IN)    :: am_I_Root     ! root CPU?
    TYPE(DgnList),     INTENT(IN)    :: DiagList
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,           INTENT(INOUT) :: RC            ! Success?
!
! !REVISION HISTORY:
!  22 Sep 2017 - E. Lundgren - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    TYPE(DgnItem), POINTER :: current
    CHARACTER(LEN=255)     :: thisLoc

    ! ================================================================
    ! Print_DiagList begins here (come back to replace with write instead)
    ! ================================================================
    thisLoc = 'Print_DiagList (diagnostics_mod.F90)'
    current => DiagList%head

    IF ( am_I_Root ) THEN
       PRINT *, " "
       PRINT *, "===================="
       PRINT *, "Contents of DiagList"
       PRINT *, " "
    ENDIF
    DO WHILE ( ASSOCIATED( current ) )

       ! Print info
       IF ( am_I_Root ) THEN
          PRINT *, TRIM(current%name)
          PRINT *, "   state:      ", TRIM(current%state)
          PRINT *, "   metadataID: ", TRIM(current%metadataID)
          PRINT *, "   registryID: ", TRIM(current%registryID)
          IF ( current%isWildcard ) THEN
             PRINT *, "   isWildcard: ", current%isWildcard
             PRINT *, "   wildcard:   ", TRIM(current%wildcard)
          ENDIF
          IF ( current%isTagged ) THEN
             PRINT *, "   isTagged:  ", current%isTagged
             PRINT *, "   tag:    ", TRIM(current%tag)
          ENDIF
          PRINT *, " "
       ENDIF

       ! Set up for next item
       current => current%next    
    ENDDO

    ! cleanup
    current => NULL()
    IF ( am_I_Root ) THEN
       PRINT *, "===================="
       PRINT *, " "
    ENDIF

  END SUBROUTINE Print_DiagList
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Cleanup_DiagList 
!
! !DESCRIPTION: Subroutine Cleanup_DiagList deallocates a DiagList
!  object and all of its member objects including the linked list of
!  DiagItem objects.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Cleanup_DiagList ( am_I_Root, DiagList, RC )
!
! !INPUT PARAMETERS:
!
    LOGICAL,           INTENT(IN)    :: am_I_Root     ! root CPU?
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(DgnList),     INTENT(INOUT) :: DiagList
!
! !OUTPUT PARAMETERS:
!
    INTEGER,           INTENT(OUT)   :: RC            ! Success?
!
! !REVISION HISTORY:
!  21 Sep 2017 - E. Lundgren - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    TYPE(DgnItem), POINTER :: current
    TYPE(DgnItem), POINTER :: next
    CHARACTER(LEN=255)     :: thisLoc

    ! ================================================================
    ! Cleanup_DiagList begins here
    ! ================================================================
    thisLoc = 'Cleanup_DiagList (diagnostics_mod.F90)'

    ! Deallocate each item in the linked list of DiagExport objects
    current => DiagList%head
    IF ( ASSOCIATED( current ) ) next => current%next
    DO WHILE ( ASSOCIATED( current ) )
       DEALLOCATE( current, STAT=RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       IF ( .NOT. ASSOCIATED ( next ) ) EXIT
       current => next
       next => current%next
    ENDDO

    ! Final cleanup
    current => NULL()
    next    => NULL()

  END SUBROUTINE Cleanup_DiagList
!EOC
END MODULE Diagnostics_Mod
