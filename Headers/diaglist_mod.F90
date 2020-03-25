!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: diaglist_mod.F90
!
! !DESCRIPTION: Module diaglist\_mod.F90 contains the derived types
!  and subroutines used for reading and storing user-configured diagnostic
!  information from the history configuration file, specifically names
!  and information derived from the names. The diagnostics list is
!  used to allocate diagnostics stored in container State\_Diag and to
!  declare exports in GCHP. It does not store collection information. A
!  module-level collection list containing names all collections that
!  are declared in HISTORY.rc with names not commented out is also in
!  this module. This is used to prevent adding diagnostics to the
!  diagnostic list that are in collections not turned on in HISTORY.rc,
!  thereby preventing their analogous State\_Diag array initialization
!  in GEOS-Chem.
!\\
!\\
! !INTERFACE:
!
MODULE DiagList_Mod
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
  PUBLIC  :: Init_DiagList
  PUBLIC  :: Print_DiagList
  PUBLIC  :: Check_DiagList
  PUBLIC  :: Cleanup_DiagList
  PUBLIC  :: Search_CollList
!
! !PRIVATE MEMBER FUNCTIONS
!
  PRIVATE :: Init_DiagItem
  PRIVATE :: InsertBeginning_DiagList
  PRIVATE :: Search_DiagList

  ! Private collection list in Init_DiagList
  ! Make public for more widespread use?
  PRIVATE :: Init_ColItem
  PRIVATE :: Set_ColItem
  PRIVATE :: InsertBeginning_ColList
  PRIVATE :: Print_ColList
  PRIVATE :: Cleanup_ColList
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
     ! could also add a list of collections this diagnostic is part of
     TYPE(DgnItem), POINTER :: next
  END TYPE DgnItem

  !=========================================================================
  ! Configurable Settings Used for Diagnostic Names at Run-time
  !=========================================================================
  CHARACTER(LEN=5),  PUBLIC  :: RadWL(3)     ! Wavelengths in radiation menu
  CHARACTER(LEN=2),  PUBLIC  :: RadFlux(12)  ! Names of RRTMG flux outputs
  INTEGER,           PUBLIC  :: nRadFlux     ! # of selected RRTMG flux outputs
  LOGICAL,           PUBLIC  :: IsFullChem   ! Is this a fullchem simulation?
  CHARACTER(LEN=10), PUBLIC  :: AltAboveSfc  ! Alt for O3, HNO3 diagnostics

  !=========================================================================
  ! Derived type for Collections List
  !=========================================================================
  TYPE, PUBLIC :: ColList
     TYPE(ColItem), POINTER :: head
  END TYPE ColList

  !=========================================================================
  ! Derived type for Collections Item (uncommented in HISTORY.rc)
  !=========================================================================
  TYPE, PUBLIC :: ColItem
     CHARACTER(LEN=63)      :: cname
     TYPE(ColItem), POINTER :: next
  END TYPE ColItem
!
! !PUBLIC DATA MEMBERS:
!
  TYPE(ColList),    PUBLIC  :: CollList      ! Collection list object
!
! !REVISION HISTORY:
!  22 Sep 2017 - E. Lundgren - Initial version
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
! !IROUTINE: Init_DiagList
!
! !DESCRIPTION: Reads the HISTORY.rc and input.geos input files to get
!  determine which GEOS-Chem diagnostics have been requested.  Then it
!  uses this information to initialize the master list of diagnostics,
!  aka, the DiagList object.
!\\
!\\
!  NOTE: This routine has to be called before any of the GEOS-Chem objects
!  Input\_Opt, State\_Chm, State\_Met, and State\_Diag are created. When using
!  GCHP, we must create the ESMF/MAPL export objects for the diagnostics
!  in the Set\_Services routine.  Set\_Services is called before GEOS-Chem
!  is initialized.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Init_DiagList ( am_I_Root, historyConfigFile, DiagList, RC )
!
! !USES:
!
    USE Charpak_Mod
    USE InquireMod,       ONLY : findFreeLun
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
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    LOGICAL                  :: EOF, found, isWildcard, isTagged
    LOGICAL                  :: InDefSection, InFieldsSection
    INTEGER                  :: QMatch, CMatch, LineNum
    INTEGER                  :: fId, IOS, N, N1, N2, N3, I, J, C, F, G
    INTEGER                  :: IWLMAX, IWLMAXLOC(1), IWL(3)
    INTEGER                  :: numSpcWords, numIDWords

    ! Strings
    CHARACTER(LEN=80 )       :: ErrorLine
    CHARACTER(LEN=255)       :: errMsg, thisLoc, nameAllCaps
    CHARACTER(LEN=255)       :: line, SubStrs(500), SubStr
    CHARACTER(LEN=255)       :: wildcard, tag, name, state
    CHARACTER(LEN=255)       :: metadataID, registryID, registryIDprefix
    CHARACTER(LEN=255)       :: collname, AttName, AttValue
    CHARACTER(LEN=255)       :: AttComp,  FieldName

    ! SAVEd variables
    CHARACTER(LEN=255), SAVE :: LastCollName

    ! Pointers
    TYPE(DgnItem),   POINTER :: NewDiagItem
    TYPE(ColItem),   POINTER :: NewCollItem

    !=======================================================================
    ! Init_DiagList begins here
    !=======================================================================

    ! Initialize
    RC              =  GC_SUCCESS
    ErrMsg          =  ''
    ThisLoc         =  ' -> at Init_DiagList (Headers/diaglist_mod.F90)'
    EOF             = .FALSE.
    found           = .FALSE.
    NewDiagItem     => NULL()
    RadWL           =  ''
    RadFlux         =  ''
    nRadFlux        =  0
    IsFullChem      = .FALSE.
    InDefSection    = .FALSE.
    InFieldsSection = .FALSE.
    LastCollName    =  ''

    ! Create DiagList object
    DiagList%head   => NULL()

    ! Create ColList object
    CollList%head   => NULL()

    !=======================================================================
    ! Read the input.geos configuration file to find out:
    ! (1) Which wavelength has been selected for optical depth diag output
    ! (2) If this is a fullchem simulation
    !=======================================================================

    ! Open input.geos file
    fId = FindFreeLun()
    OPEN( fId, FILE=TRIM('input.geos'), STATUS='OLD', IOSTAT=RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( 'Could not open input.geos!', RC, ThisLoc )
       RETURN
    ENDIF

    ! Read data from the input.geos file
    DO
       ! Read line and strip leading/trailing spaces
       Line = ReadOneLine( fId, EOF, IOS, Squeeze=.TRUE. )
       IF ( EOF ) EXIT
       IF ( IOS > 0 ) THEN
          ErrMsg = 'Unexpected end-of-file input.geos!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       ! Find out if this is a full-chemistry simulation
       IF ( INDEX( Line, 'Simulation name' ) > 0 ) THEN
          CALL StrSplit( Line, ':', SubStrs, N )
          SELECT CASE( To_UpperCase( ADJUSTL( SubStrs(2) ) ) )
             CASE( 'ACIDUPTAKE', 'APM',              'BENCHMARK',            &
                   'COMPLEXSOA', 'COMPLEXSOA_SVPOA', 'HEMCO',                &
                   'MARINEPOA',  'RRTMG',            'STANDARD',             &
                   'TROPCHEM',   'TOMAS12',          'TOMAS15',              &
                   'TOMAS30',    'TOMAS40'                                  )
                IsFullChem = .TRUE.
             CASE DEFAULT
                IsFullChem = .FALSE.
          END SELECT
       ENDIF

       ! Find the altitude above thes surface for O3, HNO3 diagnostics
       IF ( INDEX( Line, 'Diag alt above sfc [m]' ) > 0 ) THEN
          CALL StrSplit( Line, ':', SubStrs, N )
          AltAboveSfc = TRIM( ADJUSTL( SubStrs(2) ) ) // 'm'
          found       = .TRUE.
       ENDIF

       ! Update wavelength(s) with string in file
       IF ( INDEX( Line, 'AOD Wavelength' ) > 0 ) THEN
          I = INDEX( Line, ':' )
          CALL StrSplit( Line(I:), ' ', SubStrs, N )
          DO J = 1, N-1
             WRITE ( RadWL(J), "(a5)" ) SubStrs(J+1)
             RadWL(J) = ADJUSTL( RadWL(J) )
          ENDDO

          ! Exit the search
          IF ( found ) EXIT
       ENDIF
    ENDDO

    ! Close the file
    CLOSE( fId )
    found = .FALSE.

    !=======================================================================
    ! Read data from the HISTORY.rc configuration file
    !=======================================================================

    ! Open the history config file
    fId = FindFreeLun()
    OPEN( fId, FILE=TRIM(historyConfigFile), STATUS='OLD', IOSTAT=RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Could not open "' // TRIM(historyConfigFile) // '"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    ! Zero the line counter
    LineNum = 0

    !====================================================================
    ! Read history config file line by line in a loop
    !====================================================================
    DO

       ! Read line and strip leading/trailing spaces. Skip if commented out
       Line    = ReadOneLine( fId, EOF, IOS, Squeeze=.TRUE. )
       LineNum = LineNum + 1
       IF ( EOF ) EXIT
       IF ( IOS > 0 ) THEN
          ErrMsg = 'Unexpected end-of-file in "'       // &
                    TRIM( historyConfigFile ) // '" (1)!'
          WRITE( ErrorLine, '(i6)' ) LineNum
 250      FORMAT( ' -> ERROR occurred at (or near) line ', i6,               &
                      ' of the HISTORY.rc file' )
          CALL GC_Error( ErrMsg, RC, ThisLoc, ErrorLine )
          RETURN
       ENDIF
       IF ( Line(1:1) == '#' ) CYCLE

       !====================================================================
       ! Set collection name list (uncommented names only)
       !====================================================================
       IF ( INDEX( Line, 'COLLECTIONS:' ) .gt. 0 ) THEN

          ! Get the first collection name; remove commas, apost, and whitespace
          CALL CStrip( Line, KeepSpaces=.TRUE. )
          CALL StrSplit( Line, ":", SubStrs, N )
          collname = CleanText( SubStrs(2) )

          ! Read through file until end of COLLECTIONS section
          DO WHILE ( INDEX( Line, '::' ) .le. 0 )

             ! Add name to collection list if not commented out
             IF ( collname(1:1) /= "#"  ) THEN
                CALL Init_ColItem( am_I_Root, NewCollItem, collname )
                CALL InsertBeginning_ColList( am_I_Root, NewCollItem, &
                                              CollList, RC )
             ENDIF

             ! Read the next line and strip leading/trailing spaces
             Line    = ReadOneLine( fId, EOF, IOS, Squeeze=.TRUE. )
             LineNum = LineNum + 1

             IF ( IOS > 0 .OR. EOF ) THEN
                ErrMsg = 'Unexpected end-of-file in "'       // &
                          TRIM( historyConfigFile ) // '" (2)!'
                WRITE( ErrorLine, 250 ) LineNum
                CALL GC_Error( ErrMsg, RC, ThisLoc, ErrorLine )
                RETURN
             ENDIF

             ! Get next collection name
             collname = CleanText( Line )

          ENDDO
          CALL Print_ColList( am_I_Root, CollList, RC )
          CYCLE
       ENDIF

       !====================================================================
       ! Skip collection section if not in collection name list
       !====================================================================
       IF ( INDEX( Line, '.template' ) .gt. 0 .or. &
            INDEX( Line, '.filename' ) .gt. 0 ) THEN

          ! Check if collection was uncommented in the COLLECTIONS section
          CALL CStrip( Line, KeepSpaces=.TRUE. )
          CALL StrSplit( Line, ".", SubStrs, N )
          collname = CleanText( SubStrs(1) )
          CALL Search_CollList( am_I_Root, CollList, collname, Found, RC )

          ! Skip this collection if not found in list
          IF ( .NOT. Found ) THEN
             DO WHILE ( INDEX( Line, '::' ) .le. 0 )
                Line    = ReadOneLine( fId, EOF, IOS, Squeeze=.TRUE. )
                LineNum = LineNum + 1
                IF ( IOS > 0 .OR. EOF ) THEN
                   ErrMsg = 'Unexpected end-of-file in "'       // &
                             TRIM( historyConfigFile ) // '" (4)!'
                   WRITE( ErrorLine, 250 ) LineNum
                   CALL GC_Error( ErrMsg, RC, ThisLoc, ErrorLine )
                   RETURN
                ENDIF
             ENDDO
             CYCLE
          ENDIF
       ENDIF

#ifdef MODEL_CLASSIC
       !====================================================================
       ! Add some extra error checks for collections that are in the
       ! collection name list (and therefore will be archived)
       !====================================================================

       ! The double-colon indicates the end of a collection definition
       IF ( INDEX( Line, '::' ) > 0 ) THEN
          InDefSection    = .FALSE.
          InFieldsSection = .FALSE.
          LastCollName    = ''
       ENDIF

       !--------------------------------------------------------------------
       ! If the line has a "." character, then this denotes that we are
       ! in the section where collection attributes are defined
       !--------------------------------------------------------------------
       IF ( INDEX( Line, '.' ) > 0 ) THEN

          ! Denote that we are in the definition section
          InDefSection = .TRUE.

          ! Split the line into substrings
          CALL StrSplit( Line, " ", SubStrs, N )
          AttName  = SubStrs(1)           ! Attribute name
          AttValue = SubStrs(2)           ! Attribute value
          AttComp  = Substrs(3)           ! Gridded component name (for GCHP)

          ! If the .fields tag is found, then denote that we are
          ! in the section where collection fields are defined
          ! we have entered into th
          IF ( INDEX( Line, '.fields' ) > 0 ) THEN
             InFieldsSection = .TRUE.
             LastCollName    = collName
          ENDIF

          ! We expect at least 2 substrings
          IF ( LEN_TRIM( AttValue ) == 0 ) THEN
             ErrMsg = 'The value of attribute "'// TRIM( AttName )        // &
                      '" is missing! Please check the HISTORY.rc file.'
             WRITE( ErrorLine, 250 ) LineNum
             CALL GC_Error( ErrMsg, RC, ThisLoc, ErrorLine )
             RETURN
          ENDIF

          ! Each collection attribute definition needs to end with a colon
          C = LEN_TRIM( AttName )
          IF ( AttName(C:C) /= ':' ) THEN
             ErrMsg = 'The "' // TRIM( AttName ) // '" '                  // &
                      'collection attribute did not end with a ":" '      // &
                      'character!  Please check the HISTORY.rc file.'
             WRITE( ErrorLine, 250 ) LineNum
             CALL GC_Error( ErrMsg, RC, ThisLoc, ErrorLine )
             RETURN
          ENDIF

          !-----------------------------------------------------------------
          ! The "template", "filename", and "format", attribute values
          ! must be single-quoted strings that are followed by a comma,
          ! or else GCHP will die with an error.
          !-----------------------------------------------------------------
          IF ( INDEX( AttName, '.template' ) > 0   .or.                      &
               INDEX( AttName, '.filename' ) > 0   .or.                      &
               INDEX( AttName, '.format'   ) > 0  ) THEN

             ! Make sure that the value starts with a single quote
             IF ( AttValue(1:1) /= "'" ) THEN
                ErrMsg = 'The value of attribute "'// TRIM( AttName )     // &
                          '" does not begin with a single quote '         // &
                          'character! Please check the HISTORY.rc file.'
                WRITE( ErrorLine, 250 ) LineNum
                CALL GC_Error( ErrMsg, RC, ThisLoc, ErrorLine )
                RETURN
             ENDIF

             ! Make sure that the value ends with a single quote
             ! (comma is optional)
             C = LEN_TRIM( AttValue )
             IF ( AttValue(C-1:C) /= "'," ) THEN
                ErrMsg = 'The value of attribute "'// TRIM( AttName )     // &
                         '" must end with a single quote character, '     // &
                         'followed by a comma. '                          // &
                         'Please check the HISTORY.rc file.'
                WRITE( ErrorLine, 250 ) LineNum
                CALL GC_Error( ErrMsg, RC, ThisLoc, ErrorLine )
                RETURN
             ENDIF
          ENDIF

          !-----------------------------------------------------------------
          ! The "mode" attribute value must be a quoted string.  The comma
          ! is optional, or at least omitting it isn't fatal for GCHP.
          !-----------------------------------------------------------------
          IF ( INDEX( AttName, '.mode' ) > 0 ) THEN

             ! Make sure that the value starts with a single quote
             IF ( AttValue(1:1) /= "'" ) THEN
                ErrMsg = 'The value of attribute "'// TRIM( AttName )     // &
                          '" does not begin with a single quote '         // &
                          'character! Please check the HISTORY.rc file.'
                WRITE( ErrorLine, 250 ) LineNum
                CALL GC_Error( ErrMsg, RC, ThisLoc, ErrorLine )
                RETURN
             ENDIF

             ! Make sure that the value ends with a single quote
             ! (comma is optional)
             C = LEN_TRIM( AttValue )
             IF ( AttValue(C:C) == ',' ) C = C -1
             IF ( AttValue(C:C) /= "'" ) THEN
                ErrMsg = 'The value of attribute "'// TRIM( AttName )     // &
                         '" must end with a single quote character. '     // &
                         'Please check the HISTORY.rc file.'
                WRITE( ErrorLine, 250 ) LineNum
                CALL GC_Error( ErrMsg, RC, ThisLoc, ErrorLine )
                RETURN
             ENDIF
          ENDIF

       ENDIF

       !-----------------------------------------------------------------
       ! Further error checks for the fields attribute
       !-----------------------------------------------------------------
       IF ( InFieldsSection ) THEN

          ! Extract the current collection name, which forms the first
          ! part of the attribute name (up to the "." character.  If this
          ! does not match the expected collection name, then we have a
          ! missing separator ("::") somewhere.  Stop with an error.
          C = INDEX( AttName, '.' )
          IF ( AttName(1:C-1) /= TRIM( LastCollName ) ) THEN
             ErrMsg = 'Attribute "' // TRIM( AttName ) // ' specifies a ' // &
                      'value for collection "' // TRIM( AttName(1:C-1) )  // &
                      '", but the expected collection name is "'          // &
                      TRIM( LastCollName ) // '".  This indicates that '  // &
                      'the end-of-collection delimiter (i.e. "::") is '   // &
                      'missing.  Please check the HISTORY.rc file.'
             WRITE( ErrorLine, 250 ) LineNum
             CALL GC_Error( ErrMsg, RC, ThisLoc, ErrorLine )
             RETURN
          ENDIF

          ! Throw an error if we cannot find the gridcomp name
          ! (e.g. "'GIGCchem',").  GCHP will choke if this isn't found.
          G = INDEX( Line, "'GIGCchem'," )
          IF ( G == 0 ) THEN
             ErrMsg = 'The name of the GCHP gridded component '           // &
                      "(e.g. 'GIGCchem') for attribute "  // '" '         // &
                      TRIM( AttName ) // '" must be enclosed in '         // &
                      'single quotes and be followed by a comma. '        // &
                      'Please check the HISTORY.rc file.'
             WRITE( ErrorLine, 250 ) LineNum
             CALL GC_Error( ErrMsg, RC, ThisLoc, ErrorLine )
             RETURN
          ENDIF

          ! Save into LineSq the text of the line, skipping over
          ! the attribute name (if we are on the first line),
          ! as well as the gridcomp name
          F = INDEX( Line, '.fields' )
          IF ( F > 0 ) THEN
             C = INDEX( Line, ':' )
             FieldName = Line(C+1:G-1)
          ELSE
             FieldName = Line(1:G-1)
          ENDIF

          ! Pack all whitespace in LineSq
          CALL StrSqueeze( FieldName )
          CALL CStrip( FieldName )

          ! Make sure that the value starts with a single quote
          IF ( FieldName(1:1) /= "'" ) THEN
             ErrMsg = 'The diagnostic field name "' // TRIM( FieldName )  // &
                      '" does not begin with a single quote '             // &
                      'character! Please check the HISTORY.rc file.'
             WRITE( ErrorLine, 250 ) LineNum
             CALL GC_Error( ErrMsg, RC, ThisLoc, ErrorLine )
             RETURN
          ENDIF

          ! Make sure that the value ends with a single quote
          C = LEN_TRIM( FieldName )
          IF ( FieldName(C-1:C) /= "'," ) THEN
             ErrMsg = 'The diagnostic field name "' // TRIM( FieldName )  // &
                      '" must end with a single quote character '         // &
                      'followed by a comma, '                             // &
                      'Please check the HISTORY.rc file.'
             WRITE( ErrorLine, 250 ) LineNum
             CALL GC_Error( ErrMsg, RC, ThisLoc, ErrorLine )
             RETURN
          ENDIF
       ENDIF
#endif

       !====================================================================
       ! Add unique diagnostic names to diag list
       !====================================================================

       ! Skip line if GIGCchem not present
       ! GEOS-Chem is names 'GEOSCHEMCHEM' on NCCS discover,
       ! scan accordingly (ckeller, 12/29/17)
#ifdef MODEL_GEOS
       IF ( INDEX( Line, 'GEOSCHEMCHEM' ) .le. 0 ) CYCLE
#else
       IF ( INDEX( Line, 'GIGCchem' ) .le. 0 ) CYCLE
#endif

       ! Get diagnostic name
       CALL CStrip( Line, KeepSpaces=.TRUE. )
       CALL StrSplit( Line, " ", SubStrs, N )
       IF ( INDEX(Line, '.fields') > 0 .AND. N > 1 ) THEN
          name = CleanText( SubStrs(2) )
       ELSE
          name = CleanText( SubStrs(1) )
       ENDIF

       ! Skip if diagnostic name is commented out
       IF ( name(1:1) == '#' ) CYCLE

       ! Skip if name is already in diag list
       CALL Search_DiagList( am_I_Root, DiagList, name, Found, RC )
       IF ( Found ) CYCLE

       ! Set GC state
       nameAllCaps = To_Uppercase( TRIM(name) )
       IF ( nameAllCaps(1:4) == 'MET_' ) THEN
          state = 'MET'
       ELSEIF ( nameAllCaps(1:5) == 'CHEM_' ) THEN
          state = 'CHEM'
#ifdef ESMF_
       ! Emissions diagnostics are included in HISTORY.rc in GCHP only
       ELSEIF ( nameAllCaps(1:4) == 'EMIS' ) THEN
          state = 'EMISSIONS'
       ! Emissions inventory diagnostics are included in HISTORY.rc in GCHP only
       ELSEIF ( nameAllCaps(1:3) == 'INV' ) THEN
          state = 'EMISSIONS'
#ifdef MODEL_GEOS
       ! GEOS uses a different internal state prefix than GCHP and
       ! and also can have custom diagnostics
       ELSEIF ( nameAllCaps(1:5) == 'GEOS_' ) THEN
          state = 'GEOS'
       ELSEIF ( nameAllCaps(1:4) == 'TRC_' ) THEN
#else
       ELSEIF ( nameAllCaps(1:4) == 'SPC_' ) THEN
#endif
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
#ifdef ESMF_
          ! Exit with an error if using GCHP and wildcard is present
          ErrMsg = 'ERROR: HISTORY.rc wildcard handling is not ' // &
                   'implemented in GCHP: ' // TRIM(name) // '. Replace ' // &
                   'wildcard with a specific tag.'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
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
       ! Then strip off the state prefix, if any
       IF ( TRIM(state) == 'MET' ) THEN
          registryID = registryID(5:)
       ELSE IF ( TRIM(state) == 'CHEM' ) THEN
          registryID = registryID(6:)
       ENDIF
       ! Then strip off the wildcard, if any
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
       ! Then strip off the tag suffix, if any
       IF ( isTagged ) THEN
          I = INDEX( TRIM(metadataID), '_' )
          metadataID = metadataID(1:I-1)
       ENDIF

       ! For registryID and metdataID, handle special case of AOD wavelength
       ! Update registryID
       IWL(1) = INDEX( TRIM(registryID), 'WL1' )
       IWL(2) = INDEX( TRIM(registryID), 'WL2' )
       IWL(3) = INDEX( TRIM(registryID), 'WL3' )
       IWLMAX = MAX(IWL(1),IWL(2),IWL(3))
       IF ( IWLMAX > 0 ) THEN
          IWLMAXLOC = MAXLOC(IWL)
          registryIDprefix = registryID(1:IWL(IWLMAXLOC(1))-1) // &
                             TRIM(RadWL(IWLMAXLOC(1))) // 'NM'
          I = INDEX( TRIM(registryID), '_' )
          IF ( I > 0 ) THEN
             registryID = TRIM(registryIDprefix) // registryID(I:)
          ELSE
             registryID = registryIDprefix
          ENDIF
       ENDIF
       ! Update metadataID
       IWL(1) = INDEX( TRIM(metadataID), 'WL1' )
       IWL(2) = INDEX( TRIM(metadataID), 'WL2' )
       IWL(3) = INDEX( TRIM(metadataID), 'WL3' )
       IWLMAX = MAX(IWL(1),IWL(2),IWL(3))
       IF ( IWLMAX > 0 ) THEN
          IWLMAXLOC = MAXLOC(IWL(:))
          metadataID = metadataID(1:IWL(IWLMAXLOC(1))-1) //  &
                       TRIM(RadWL(IWLMAXLOC(1))) // 'NM'
       ENDIF

       ! Special handling for the RRTMG diagnostic flux outputs
       ! Store the list of the requested fluxes in RadFlux
       IWL(1) = INDEX( TRIM(metadataID), 'RADCLRSKY' )
       IWL(2) = INDEX( TRIM(metadataID), 'RADALLSKY' )
       IWLMAX = MAX( IWL(1), IWL(2) )
       IF ( IWLMAX > 0 ) THEN

          IF ( LEN_TRIM( Tag ) > 0 ) THEN

             ! If a tag is specified explicitly, then add each
             ! tag name to the RadFlux array and update nRadFlux
             IF ( .not. ANY( RadFlux == TRIM(Tag) ) ) THEN
                nRadFlux          = nRadFlux + 1
                RadFlux(nRadFlux) = TRIM( Tag )
             ENDIF

          ELSE

             ! If a tag is not explicity stated, then manually define all
             ! slots of the RadFlux array (and update nRadFlux accordingly)
             ! See: wiki.geos-chem.org/Coupling_GEOS-Chem_with_RRTMG
             RadFlux(1 ) = 'O3'
             RadFlux(2 ) = 'ME'
             RadFlux(3 ) = 'SU'
             RadFlux(4 ) = 'NI'
             RadFlux(5 ) = 'AM'
             RadFlux(6 ) = 'BC'
             RadFlux(7 ) = 'OA'
             RadFlux(8 ) = 'SS'
             RadFlux(9 ) = 'DU'
             RadFlux(10) = 'PM'
             RadFlux(11) = 'ST'
             nRadFlux    = 11

          ENDIF
       ENDIF

       ! Special handling for diagnostics at a specific height
       ! (e.g. rename O3CONCATALT --> O3CONCAT10M)
       IWL(1) = INDEX( TRIM(registryID), 'ALT1' )
       IF ( IWL(1) > 0 ) THEN
          registryIDprefix = registryID(1:IWL(1)-1) // TRIM( AltAboveSfc )
          I = INDEX( TRIM(registryID), '_' )
          IF ( I > 0 ) THEN
             registryID = TRIM(registryIDprefix) // registryID(I:)
          ELSE
             registryID = registryIDprefix
          ENDIF
       ENDIF
       IWL(2) = INDEX( TRIM(metadataID), 'ALT1' )
       IF ( IWL(2) > 0 ) THEN
          metadataID = metadataID(1:IWL(2)-1) // TRIM( AltAboveSfc )
       ENDIF

       !====================================================================
       ! Create a new DiagItem object
       !====================================================================
       CALL Init_DiagItem( am_I_Root,              &
                           NewDiagItem,            &
                           name=name,              &
                           state=state,            &
                           metadataID=metadataID,  &
                           registryID=registryID,  &
                           isWildcard=isWildcard,  &
                           wildcard=wildcard,      &
                           isTagged=isTagged,      &
                           tag=tag,                &
                           RC=RC  )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Error initializing DiagItem ' // TRIM(name)
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       !====================================================================
       ! Add new DiagItem to linked list
       !====================================================================
       CALL InsertBeginning_DiagList( am_I_Root, NewDiagItem, DiagList, RC )
       IF ( RC /= GC_SUCCESS ) RETURN

    ENDDO

    !====================================================================
    ! Close the file
    !====================================================================
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
! !DESCRIPTION: Initializes a DiagItem object, which contains information
!  about a single GEOS-Chem diagnostic.  Several DiagItem objects will be
!  linked together in the master diagnostics list (DiagList).
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
!  See https://github.com/geoschem/geos-chem for complete history
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
    thisLoc = 'Init_DiagItem (diaglist_mod.F90)'

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
! !IROUTINE: Init_ColItem
!
! !DESCRIPTION: Initializes a ColItem object, which contains information
!  about a single GEOS-Chem collection.  Several ColItem objects will be
!  linked together in the master collections list (ColList).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Init_ColItem ( am_I_Root, NewCollItem, cname, RC  )
!
! !INPUT PARAMETERS:
!
    LOGICAL,        INTENT(IN) :: am_I_Root
    TYPE(ColItem),  POINTER    :: NewCollItem
    CHARACTER(LEN=*)           :: cname
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        OPTIONAL   :: RC
!
! !REVISION HISTORY:
!  25 Jan 2018 - E. Lundgren - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    CHARACTER(LEN=255) :: thisLoc

    ! ================================================================
    ! Init_ColItem begins here
    ! ================================================================
    thisLoc = 'Init_ColItem (diaglist_mod.F90)'

    ALLOCATE(NewCollItem)
    NewCollItem%cname      = TRIM(cname)

  END SUBROUTINE Init_ColItem
!EOC
!!------------------------------------------------------------------------------
!!                  GEOS-Chem Global Chemical Transport Model                  !
!!------------------------------------------------------------------------------
!!BOP
!!
!! !IROUTINE: Get_ColItem
!!
!! !DESCRIPTION: Gets a pointer to a collection (ColItem) object
!!\\
!!\\
!! !INTERFACE:
!!
!  SUBROUTINE Get_ColItem ( am_I_Root, CollName, CollList, CollItem, Found, RC )
!!
!! !INPUT PARAMETERS:
!!
!    LOGICAL,          INTENT(IN) :: am_I_Root
!    CHARACTER(LEN=*), INTENT(IN) :: CollName
!    TYPE(ColList),    INTENT(IN) :: CollList
!!
!! !OUTPUT PARAMETERS:
!!
!    TYPE(ColItem),    POINTER    :: CollItem
!    LOGICAL,          OPTIONAL   :: Found
!    INTEGER,          OPTIONAL   :: RC
!!
!! !REVISION HISTORY:
!!  25 Jan 2018 - E. Lundgren - Initial version
!!  See https://github.com/geoschem/geos-chem for complete history
!!EOP
!!------------------------------------------------------------------------------
!!BOC
!!
!! !LOCAL VARIABLES:
!!
!    TYPE(ColItem), POINTER :: current
!    CHARACTER(LEN=255)     :: thisLoc
!
!    ! ================================================================
!    ! Get_ColList begins here
!    ! ================================================================
!
!    ! Initialize
!    thisLoc = 'Get_ColItem (diaglist_mod.F90)'
!    IF ( PRESENT( Found ) ) Found = .FALSE.
!
!    ! Search for name in list
!    current => CollList%head
!    DO WHILE ( ASSOCIATED( current ) )
!       IF ( current%cname == CollName ) THEN
!          IF ( PRESENT( Found ) ) Found = .TRUE.
!          CollItem = current
!          EXIT
!       ENDIF
!       current => current%next
!    ENDDO
!    current => NULL()
!
!  END SUBROUTINE Get_ColItem
!!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Set_ColItem
!
! !DESCRIPTION: Sets a ColItem object, which contains information
!  about a single GEOS-Chem collection.  Several ColItem objects will be
!  linked together in the master collections list (ColList).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Set_ColItem ( am_I_Root, Collname, CollList, Found, RC  )
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN) :: am_I_Root
    CHARACTER(LEN=*)             :: CollName
    TYPE(ColList)                :: CollList
!
! !OUTPUT PARAMETERS:
!
    LOGICAL,          OPTIONAL   :: Found
    INTEGER,          OPTIONAL   :: RC
!
! !REVISION HISTORY:
!  25 Jan 2018 - E. Lundgren - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    CHARACTER(LEN=255)     :: thisLoc
    TYPE(ColItem), POINTER :: current

    ! ================================================================
    ! Set_ColList begins here
    ! ================================================================

    ! Initialize
    thisLoc = 'Set_ColItem (diaglist_mod.F90)'
    IF ( PRESENT( Found ) ) Found = .FALSE.

    ! Search for name in list
    current => CollList%head
    DO WHILE ( ASSOCIATED( current ) )
       IF ( current%cname == CollName ) THEN
          IF ( PRESENT( Found ) ) Found = .TRUE.
          EXIT
       ENDIF
       current => current%next
    ENDDO

    ! Exit with error if no collection matches the input name
    IF ( .NOT. FOUND ) THEN
       CALL GC_ERROR("Error setting collection item", RC, ThisLoc)
       RETURN
    ENDIF

    ! Null pointer
    current => NULL()

  END SUBROUTINE Set_ColItem
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: InsertBeginning_DiagList
!
! !DESCRIPTION: Inserts a new node at the beginning of the DiagList linked
!  list object.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE InsertBeginning_DiagList ( am_I_Root, DiagItem, DiagList, RC )
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
!  See https://github.com/geoschem/geos-chem for complete history
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
    thisLoc = 'InsertBeginning_DiagList (diaglist_mod.F90)'

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
! !IROUTINE: InsertBeginning_ColList
!
! !DESCRIPTION: Inserts a new node at the beginning of the ColList linked
!  list object.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE InsertBeginning_ColList ( am_I_Root, CollItem, CollList, RC )
!
! !INPUT PARAMETERS:
!
    LOGICAL,         INTENT(IN)    :: am_I_Root
    TYPE(ColItem),   POINTER       :: CollItem
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ColList),   INTENT(INOUT) :: CollList
!
! !OUTPUT PARAMETERS:
!
    INTEGER,         INTENT(OUT)   :: RC
!
! !REVISION HISTORY:
!  25 Jan 2018 - E. Lundgren - initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    TYPE(ColItem),  POINTER :: NewCollItem
    CHARACTER(LEN=255)      :: thisLoc

    ! ================================================================
    ! InsertBeginning_ColList begins here
    ! ================================================================
    thisLoc = 'InsertBeginning_ColList (diaglist_mod.F90)'

    ! Add new object to the beginning of the linked list
    CollItem%next => CollList%head
    CollList%head => CollItem

  END SUBROUTINE InsertBeginning_ColList
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Search_DiagList
!
! !DESCRIPTION: Searches for a given diagnostic name within the DiagList
!  diagnostic list object.
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
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    TYPE(DgnItem), POINTER :: current
    CHARACTER(LEN=255)     :: thisLoc

    ! Initialize
    thisLoc = 'Search_DiagList (diaglist_mod.F90)'
    found = .FALSE.

    ! Search for name in list
    current => DiagList%head
    DO WHILE ( ASSOCIATED( current ) )
       IF ( TRIM(current%name) == TRIM(name) ) THEN
          found = .TRUE.
          current=> NULL()
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
! !IROUTINE: Search_CollList
!
! !DESCRIPTION: Searches for a given collection name within the ColList
!  collection list object.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Search_CollList ( am_I_Root, CollList, name, found, RC )
!
! !INPUT PARAMETERS:
!
    LOGICAL,           INTENT(IN) :: am_I_Root
    TYPE(ColList),     INTENT(IN) :: CollList
    CHARACTER(LEN=*),  INTENT(IN) :: name
!
! !OUTPUT PARAMETERS:
!
    LOGICAL,           INTENT(OUT) :: found
    INTEGER,           INTENT(OUT) :: RC
!
! !REVISION HISTORY:
!  25 Jan 2018 - E. Lundgren - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    TYPE(ColItem), POINTER :: current
    CHARACTER(LEN=255)     :: thisLoc

    ! Initialize
    thisLoc = 'Search_CollList (diaglist_mod.F90)'
    found = .FALSE.

    ! Search for name in list
    current => CollList%head
    DO WHILE ( ASSOCIATED( current ) )
       IF ( TRIM(current%cname) == TRIM(name) ) THEN
          found = .TRUE.
          EXIT
       ENDIF
       current => current%next
    ENDDO
    current => NULL()

  END SUBROUTINE Search_CollList
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Check_DiagList
!
! !DESCRIPTION: Returns TRUE if a string matches the metadataID for at
!  least one diagnostic in the passed DiagList object. An exact match is
!  required (case-insensitive) unless optional partial logical argument
!  is passed.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Check_DiagList( am_I_Root, DiagList, name, found, RC, partial )
!
! !USES:
!
    USE Charpak_Mod, ONLY : To_UpperCase
!
! !INPUT PARAMETERS:
!
    LOGICAL,           INTENT(IN)  :: am_I_Root   ! Are we on the root CPU?
    TYPE(DgnList),     INTENT(IN)  :: DiagList    ! Diagnostic list object
    CHARACTER(LEN=*),  INTENT(IN)  :: name        ! Diagnostic metadata name
    LOGICAL,           OPTIONAL    :: partial     ! Allow partial name match?
!
! !OUTPUT PARAMETERS:
!
    LOGICAL,           INTENT(OUT) :: found       ! Was a match found (T/F)?
    INTEGER,           INTENT(OUT) :: RC          ! Success or failure?
!
! !REVISION HISTORY:
!  22 Sep 2017 - E. Lundgren - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    LOGICAL                :: doPartialMatch
    INTEGER                :: matchInd
    INTEGER                :: matchLen

    ! Strings
    CHARACTER(LEN=255)     :: thisLoc
    CHARACTER(LEN=255)     :: inName_AllCaps
    CHARACTER(LEN=255)     :: currentName_AllCaps

    ! Pointers
    TYPE(DgnItem), POINTER :: current

    ! Initialize
    RC      = GC_SUCCESS
    thisLoc = ' -> at Check_DiagList (in module Headers/diaglist_mod.F90)'
    found   = .FALSE.

    ! Get the optional exactMatch argument, whichg determines
    ! if we should force an exact name match or not (bmy, 10/29/18)
    IF ( PRESENT( partial ) ) THEN
       doPartialMatch = partial
    ELSE
       doPartialMatch = .FALSE.
    ENDIF

    ! Convert strings to uppercase for comparison
    inName_AllCaps = To_Uppercase( TRIM( name ) )

    ! Get the length of inName_AllCaps excluding whitespace
    matchLen = LEN_TRIM( inName_AllCaps )

    ! Search for name in list
    current => DiagList%head
    DO WHILE ( ASSOCIATED( current ) )

       ! Name of the diagnostic metadata at this point in the linked list
       currentName_AllCaps = To_Uppercase( current%metadataID )

       ! Test if the substring matches all or part of the diagnostic name
       matchInd = INDEX( currentName_AllCaps, TRIM( inName_AllCaps ) )

       ! Determine if we need to have an exact or partial match
       IF ( .NOT. doPartialMatch ) THEN

          ! Exact match: inName_AllCaps matches a sequence of characters
          ! starting with the first character of currentName_AllCaps.
          ! AND has the same trimmed length as currentName_AllCaps
          IF ( ( matchInd == 1                               )   .and.       &
               ( matchLen == LEN_TRIM( currentName_AllCaps ) ) ) THEN
             found = .TRUE.
             EXIT
          ENDIF

       ELSE

          ! Partial match: inName_AllCaps matches a sequence of characters
          ! somewhere within currentName_AllCaps (but not necessarily
          ! starting from the beginning).
          IF ( matchInd > 0 ) THEN
             found = .TRUE.
             EXIT
          ENDIF

       ENDIF

       ! Move to next diagnostic in the linked list
       current => current%next
    ENDDO

    ! Free pointer
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
! !DESCRIPTION: Subroutine Print\_DiagList prints information for all
!  DiagItem members in a DiagList linked list.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Print_DiagList( am_I_Root, DiagList, RC )
!
! !INPUT PARAMETERS:
!
    LOGICAL,           INTENT(IN)    :: am_I_Root
    TYPE(DgnList),     INTENT(IN)    :: DiagList
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,           INTENT(INOUT) :: RC            ! Success?
!
! !REVISION HISTORY:
!  22 Sep 2017 - E. Lundgren - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
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
    thisLoc = 'Print_DiagList (diaglist_mod.F90)'
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
    IF ( am_I_Root ) PRINT *, " "

  END SUBROUTINE Print_DiagList
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Print_ColList
!
! !DESCRIPTION: Subroutine Print\_ColList prints information for all
!  ColItem members in a ColList linked list.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Print_ColList( am_I_Root, CollList, RC )
!
! !INPUT PARAMETERS:
!
    LOGICAL,           INTENT(IN)    :: am_I_Root
    TYPE(ColList),     INTENT(IN)    :: CollList
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,           INTENT(INOUT) :: RC            ! Success?
!
! !REVISION HISTORY:
!  25 Jan 2018 - E. Lundgren - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    TYPE(ColItem), POINTER :: current
    CHARACTER(LEN=255)     :: thisLoc

    ! ================================================================
    ! Print_ColList begins here (come back to replace with write instead)
    ! ================================================================
    thisLoc = 'Print_ColList (diaglist_mod.F90)'
    current => CollList%head

    IF ( am_I_Root ) THEN
       PRINT *, " "
       PRINT *, "======================"
       PRINT *, "Contents of  CollList"
       PRINT *, " "
    ENDIF
    DO WHILE ( ASSOCIATED( current ) )

       ! Print info
       IF ( am_I_Root ) THEN
          PRINT *, TRIM(current%cname)
       ENDIF

       ! Set up for next item
       current => current%next
    ENDDO

    ! cleanup
    current => NULL()
    IF ( am_I_Root ) PRINT *, " "

  END SUBROUTINE Print_ColList
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Cleanup_DiagList
!
! !DESCRIPTION: Subroutine Cleanup\_DiagList deallocates a DiagList
!  object and all of its member objects including the linked list of
!  DiagItem objects.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Cleanup_DiagList( DiagList, RC )
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
!  See https://github.com/geoschem/geos-chem for complete history
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
    thisLoc = 'Cleanup_DiagList (diaglist_mod.F90)'

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

    ! Also get rid of module-level collection list when cleanup up diaglist
    CALL Cleanup_ColList( CollList, RC )

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
! !IROUTINE: Cleanup_ColList
!
! !DESCRIPTION: Subroutine Cleanup\_ColList deallocates a ColList
!  object and all of its member objects including the linked list of
!  ColItem objects.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Cleanup_ColList ( CollList, RC )
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ColList),     INTENT(INOUT) :: CollList
!
! !OUTPUT PARAMETERS:
!
    INTEGER,           INTENT(OUT)   :: RC            ! Success?
!
! !REVISION HISTORY:
!  21 Sep 2017 - E. Lundgren - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    TYPE(ColItem), POINTER :: current
    TYPE(ColItem), POINTER :: next
    CHARACTER(LEN=255)     :: thisLoc

    ! ================================================================
    ! Cleanup_ColList begins here
    ! ================================================================
    thisLoc = 'Cleanup_ColList (diaglist_mod.F90)'

    ! Deallocate each item in the linked list of collection objects
    current => CollList%head
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

  END SUBROUTINE Cleanup_ColList
!EOC
END MODULE DiagList_Mod
