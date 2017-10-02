!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: diagnostics_mod.F90
!
! !DESCRIPTION: Module diagnostics\_mod.F90 contains the derived types
!  and subroutines used for reading and storing user-configured diagnostic
!  information.
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
  PRIVATE :: Append_DiagList
  PRIVATE :: Search_DiagList
  PRIVATE :: ReadOneLine ! copied from history_mod. Consider consolidating.
!
! !PUBLIC DATA TYPES:
!
  !=========================================================================
  ! Derived type for Diagnostics List
  !=========================================================================
  TYPE, PUBLIC :: DgnList
     TYPE(DgnItem), POINTER  :: head
     INTEGER                 :: numItems
  END TYPE DgnList

  !=========================================================================
  ! Derived type for Diagnostics Item (unique item in HISTORY.rc)
  !=========================================================================
  TYPE, PUBLIC :: DgnItem
     CHARACTER(LEN=63)      :: name 
     CHARACTER(LEN=7)       :: state
     CHARACTER(LEN=63)      :: metadataID
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
!              
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
    USE Charpak_Mod,      ONLY: STRSPLIT
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
    CHARACTER(LEN=255)      :: errMsg, thisLoc
    CHARACTER(LEN=255)      :: line, SubStrs(500), SubStr
    CHARACTER(LEN=255)      :: wildcard, species, name, state, metadataID
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
    DiagList%numItems = 0

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
    
       ! Skip line if certain conditions met
       IF ( line(1:1) == "#" ) CYCLE
       IF ( INDEX( line, 'GIGCchem' ) .le. 0 ) CYCLE
    
       ! Get name in history config file
       CALL StrSplit( line, ' ', SubStrs, N )
       IF ( INDEX(line, '.fields') .le. 0 ) THEN
          SubStr = TRIM(SubStrs(1))
       ELSE 
          SubStr = TRIM(SubStrs(2))
       ENDIF
       CALL StrSplit( SubStr, "'", SubStrs, N )
       name = TRIM( SubStrs(1) )

       ! Get state prefix
       CALL StrSplit( name, '_', SubStrs, N )
       IF ( SubStrs(1) == 'MET' .or. SubStrs(1) == 'CHEM' ) THEN
          state = TRIM( SubStrs(1) )
       ELSE
          state = 'DIAG'
       ENDIF

       ! Get wildcard or species, if any
       isWildcard = .FALSE.
       isSpecies  = .FALSE.
       wildcard   = ''
       species    = ''
       IF ( INDEX ( name, '__' ) > 0 ) THEN
          IF ( INDEX( name, '?' ) > 0 ) THEN
#if defined( ESMF_ )
             ! Exit with an error if using GCHP and wildcard is present
             CALL GC_Error( 'Wildcards not allowed in GCHP', RC, ThisLoc )
             RETURN
#endif
             isWildcard = .TRUE.
             CALL StrSplit( name, '?', SubStrs, N )
             wildcard = SubStrs(N-1)
          ELSE 
             CALL StrSplit( name, '__', SubStrs, N )
             isSpecies  = .TRUE.
             species = SubStrs(N)
          ENDIF
       ENDIF
       
       ! Get metadataID which is the name w/out state prefix and spc suffix
       ! First, trim state prefix
       IF ( TRIM( state ) == 'DIAG' ) THEN
          metadataID = name
       ELSEIF ( TRIM( state ) == 'MET' ) THEN
          metadataID = name(5:)
       ELSEIF ( TRIM( state ) == 'CHEM' ) THEN
          metadataID = name(6:)
       ENDIF
       ! Next, remove wildcard or species name if present
       IF ( isWildcard .OR. isSpecies ) THEN
          CALL StrSplit( metadataID, '__', SubStrs, N )
          metadataID = SubStrs(1)
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
       CALL Append_DiagList( am_I_Root, NewDiagItem, DiagList, RC )
       IF ( RC /= GC_SUCCESS ) RETURN

    ENDDO
    
    ! Close the file
    CLOSE( fId )

  END SUBROUTINE Init_DiagList
!EOC
!------------------------------------------------------------------------------
!              
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
  SUBROUTINE Init_DiagItem ( am_I_Root,  NewDiagItem, name,     state,     &
                             metadataID, isWildcard,  wildcard, isSpecies, &
                             species,    RC  )
!
! !INPUT PARAMETERS:
!
    LOGICAL,             INTENT(IN) :: am_I_Root
    CHARACTER(LEN=*),    OPTIONAL   :: name
    CHARACTER(LEN=*),    OPTIONAL   :: state
    CHARACTER(LEN=*),    OPTIONAL   :: metadataID
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
    NewDiagItem%isWildcard = isWildcard
    NewDiagItem%wildcard   = TRIM(wildcard)
    NewDiagItem%isSpecies  = isSpecies
    NewDiagItem%species    = TRIM(species)

  END SUBROUTINE Init_DiagItem
!EOC
!------------------------------------------------------------------------------
!              
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Append_DiagList 
!
! !DESCRIPTION: 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Append_DiagList ( am_I_Root, DiagItem, DiagList, RC )
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
    ! Append_DiagList begins here
    ! ================================================================
    thisLoc = 'Append_DiagList (diagnostics_mod.F90)'

    ! Add new object to the beginning of the linked list
    DiagItem%next => DiagList%head
    DiagList%head => DiagItem

    ! Update # of list items
    DiagList%numItems = DiagList%numItems + 1

  END SUBROUTINE Append_DiagList
!EOC
!------------------------------------------------------------------------------
!
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
!
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
!
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
!
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

    current => NULL()
    next => NULL()

    ! Cleanup each item in the linked list of DiagExport objects
    current => DiagList%head
    IF ( ASSOCIATED( current ) ) next => current%next
    DO WHILE ( ASSOCIATED( current ) )
       DEALLOCATE( current, STAT=RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       IF ( .NOT. ASSOCIATED ( next ) ) EXIT
       current => next
       next => current%next
    ENDDO

    ! Deallocate the DiagList object
    !DEALLOCATE( DiagList, STAT=RC )
    !IF ( RC /= GC_SUCCESS ) RETURN

    ! Final cleanup
    current => NULL()
    next    => NULL()

  END SUBROUTINE Cleanup_DiagList
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
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
    USE Charpak_Mod,  ONLY : StrSqueeze
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

    ! If desired, call StrSqueeze to strip leading and trailing blanks
    IF ( PRESENT( Squeeze ) ) THEN
       IF ( Squeeze ) THEN
          CALL StrSqueeze( Line )
       ENDIF
    ENDIF

  END FUNCTION ReadOneLine
!EOC

END MODULE Diagnostics_Mod
