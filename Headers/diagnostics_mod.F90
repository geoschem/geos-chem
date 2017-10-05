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
  PRIVATE :: ReadOneLine ! copied from history_mod. Consider consolidating.
  PRIVATE :: CleanText   ! copied from history_mod. Consider consolidating.
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
     CHARACTER(LEN=7)       :: state
     CHARACTER(LEN=63)      :: metadataID
     CHARACTER(LEN=63)      :: registryID
     LOGICAL                :: isWildcard
     CHARACTER(LEN=7)       :: wildcard
     LOGICAL                :: isSpecies
     CHARACTER(LEN=63)      :: species
     TYPE(DgnItem), POINTER :: next
  END TYPE DgnItem
!
! !REVISION HISTORY:
!  22 Sep 2017 - E. Lundgren - Initial version
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
    USE Charpak_Mod,      ONLY: StrSplit
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
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER                 :: fId, IOS, N, I
    INTEGER                 :: numSpcWords, numIDWords
    CHARACTER(LEN=255)      :: errMsg, thisLoc
    CHARACTER(LEN=255)      :: line, SubStrs(500), SubStr
    CHARACTER(LEN=255)      :: wildcard, species, name, state
    CHARACTER(LEN=255)      :: metadataID, registryID
    LOGICAL                 :: EOF, found, isWildcard, isSpecies
    TYPE(DgnItem),  POINTER :: NewDiagItem

    ! ================================================================
    ! Init_DiagList begins here
    ! ================================================================
    thisLoc = 'Init_DiagList (diagnostics_mod.F90)'

    ! Init
    EOF = .FALSE.
    found = .FALSE.
    NewDiagItem => NULL()

    ! Create DiagList object
    DiagList%head => NULL()

    ! Open the file
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
       ! TODO: Beware that this currently does not strip tabs!!!
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

       ! Get GC state
       IF ( INDEX( name, '_') > 0 ) THEN
          CALL StrSplit( name, '_', SubStrs, N )
          IF ( SubStrs(1) == 'MET' .or. SubStrs(1) == 'CHEM' ) THEN
             state = TRIM( SubStrs(1) )
          ELSE
             state = 'DIAG'
          ENDIF
       ELSE
          ! TODO: need to not classify GCHP internal state as DIAG!
          state = 'DIAG'
       ENDIF

       ! Get wildcard or species, if any 
       ! NOTE: Must be prefaced with double underscore in HISTORY.rc!
       isWildcard = .FALSE.
       isSpecies  = .FALSE.
       wildcard   = ''
       species    = ''
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
       IF ( .NOT. isWildcard ) THEN
          I = INDEX( TRIM(name), '__' ) 
          IF ( I > 0 ) THEN
             isSpecies = .TRUE.
             species = name(I+2:)
          ENDIF
       ENDIF

       ! Get registryID (name minus state prefix and wildcard suffix, but 
       ! keep species)
       registryID = TRIM(name)
       IF ( TRIM(state) == 'MET' ) THEN
          registryID = registryID(5:)
       ELSE IF ( TRIM(state) == 'CHEM' ) THEN
          registryID = registryID(6:)
       ENDIF
       IF ( isWildcard ) THEN
          I = INDEX( TRIM(registryID), '__' ) 
          IF ( I .le. 0 ) THEN
             ErrMsg = 'Error setting registryID. Double underscore must' &
                      // ' precede wildcard in HISTORY.rc!'
             CALL GC_ERROR( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF
          registryID = registryID(1:I-1)
       ENDIF

       ! Get metadataID (name minus state prefix and species/wildcard suffix)
       ! Start with registryID since that already has state prefix and wilcard
       ! (if any) stripped. It only needs species stripped, if any.
       metadataID = registryID
       IF ( isSpecies ) THEN
          I = INDEX( TRIM(metadataID), '__' ) 
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
                           isSpecies=isSpecies,     &
                           species=species,        &
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
                             isSpecies,  species,     RC  )
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
    LOGICAL,             OPTIONAL   :: isSpecies
    CHARACTER(LEN=*),    OPTIONAL   :: species
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
    NewDiagItem%isSpecies  = isSpecies
    NewDiagItem%species    = TRIM(species)

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
  SUBROUTINE Check_DiagList ( am_I_Root, DiagList, substr, found, RC )
!
! !USES:
!
    USE Registry_Mod, ONLY : To_Uppercase
!
! !INPUT PARAMETERS:
!
    LOGICAL,           INTENT(IN) :: am_I_Root
    TYPE(DgnList),     INTENT(IN) :: DiagList
    CHARACTER(LEN=*),  INTENT(IN) :: substr
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
       IF ( ( current%state == 'DIAG' ) .AND. &
            INDEX( currentName_AllCaps, TRIM( substr_AllCaps ) ) > 0 ) THEN
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
  SUBROUTINE Print_DiagList ( am_I_Root, DiagList, RC )
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
          PRINT *, "   registryID (GCHP-only): ", TRIM(current%registryID)
          PRINT *, "   isWildcard: ", current%isWildcard
          IF ( current%isWildcard ) THEN
             PRINT *, "   wildcard:   ", TRIM(current%wildcard)
          ENDIF
          PRINT *, "   isSpecies:  ", current%isSpecies
          IF ( current%isSpecies ) THEN
             PRINT *, "   species:    ", TRIM(current%species)
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
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!BOP
!
! !IROUTINE: ReadOneLine
!
! !DESCRIPTION: Subroutine READ\_ONE\_LINE reads a line from the input file.  
!  If the global variable VERBOSE is set, the line will be printed to stdout.  
!  READ\_ONE\_LINE can trap an unexpected EOF if LOCATION is passed.  
!  Otherwise, it will pass a logical flag back to the calling routine, 
!  where the error trapping will be done.
!\\
!\\
! !INTERFACE:
!
  FUNCTION ReadOneLine( fId, EndOfFile, IoStatus, Squeeze ) RESULT( Line )
!
! !USES:
!
    USE Charpak_Mod,  ONLY : StrSqueeze, CStrip
!
! !INPUT PARAMETERS:
!
    INTEGER, INTENT(IN)  :: fId        ! File unit number
    LOGICAL, OPTIONAL    :: Squeeze    ! Call Strsqueeze?
!
! !OUTPUT PARAMETERS:
!
    LOGICAL, INTENT(OUT) :: EndOfFile  ! Denotes EOF condition
    INTEGER, INTENT(OUT) :: IoStatus   ! I/O status code
!
! !RETURN VALUE:
!
    CHARACTER(LEN=255)   :: Line       ! Single line from the input file
! 
! !REVISION HISTORY: 
!  16 Jun 2017 - R. Yantosca - Initial version, based on GEOS-Chem
!  22 Sep 2017 - E. Lundgren - Copied from history_mod to allow use in
!                              Headers subdirectory. Consider placing
!                              elsewhere for common use!!!
!EOP
!------------------------------------------------------------------------------
!BOC

    !=================================================================
    ! Initialize
    !=================================================================
    EndOfFile = .FALSE.
    IoStatus  = 0
    Line      = ''

    !=================================================================
    ! Read data from the file
    !=================================================================

    ! Read a line from the file
    READ( fId, '(a)', IOSTAT=IoStatus ) Line

    ! IO Status < 0: EOF condition
    IF ( IoStatus < 0 ) THEN 
       EndOfFile = .TRUE.
       RETURN
    ENDIF

    ! Remove blanks and null characters
    CALL CSTRIP( Line ) 

    ! If desired, call StrSqueeze to strip leading and trailing blanks
    IF ( PRESENT( Squeeze ) ) THEN
       IF ( Squeeze ) THEN
          CALL StrSqueeze( Line )
       ENDIF
    ENDIF

  END FUNCTION ReadOneLine
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: CleanText
!
! !DESCRIPTION: Strips commas, apostrophes, spaces, and tabs from a string.
!\\
!\\
! !INTERFACE:
!
  FUNCTION CleanText( Str ) RESULT( CleanStr )
!
! !USES:
!
    USE Charpak_Mod, ONLY : CStrip, StrRepl, StrSqueeze
!
! !INPUT PARAMETERS: 
!
    CHARACTER(LEN=*), INTENT(IN) :: Str        ! Original string
!
! !RETURN VALUE
!
    CHARACTER(LEN=255)           :: CleanStr   ! Cleaned-up string
!
! !REMARKS:
!
! !REVISION HISTORY:
!  06 Jan 2015 - R. Yantosca - Initial version
!  21 Jun 2017 - R. Yantosca - Now call CSTRIP to remove tabs etc.
!  05 Oct 2017 - E. Lundgren - Copied from history_mod to allow use in
!                              Headers subdirectory. Consider placing
!                              elsewhere for common use!!!
!EOP
!------------------------------------------------------------------------------
!BOC

    ! Initialize
    CleanStr = Str

    ! Strip out non-printing characters (e.g. tabs)
    CALL CStrip    ( CleanStr           )

    ! Remove commas and quotes
    CALL StrRepl   ( CleanStr, ",", " " )
    CALL StrRepl   ( CleanStr, "'", " " )
    
    ! Remove leading and trailing spaces
    CALL StrSqueeze( CleanStr           ) 

  END FUNCTION CleanText
!EOC
END MODULE Diagnostics_Mod
