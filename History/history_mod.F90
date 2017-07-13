!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: history_mod.F90
!
! !DESCRIPTION: Driver module for the new netCDF diagnostics package.
!\\
!\\
! !INTERFACE:
!
MODULE History_Mod
!
! !USES:
!
  USE History_Params_Mod
  USE Precision_Mod
  USE HistContainer_MOd,     ONLY : HistContainer
  USE MetaHistContainer_Mod, ONLY : MetaHistContainer

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
! 
  PUBLIC  :: History_Init
  PUBLIC  :: History_Cleanup
!
! PRIVATE MEMBER FUNCTIONS:

  PRIVATE :: History_ReadCollectionNames
  PRIVATE :: History_ReadCollectionData
  PRIVATE :: CleanText
  PRIVATE :: ReadOneLine
!
! !REMARKS:
!
! !REVISION HISTORY:
!  06 Jan 2015 - R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
  ! Scalars
  INTEGER                          :: CollectionCount

  ! Strings
  CHARACTER(LEN=255), ALLOCATABLE  :: CollectionName      (:)
  CHARACTER(LEN=255), ALLOCATABLE  :: CollectionTemplate  (:)
  CHARACTER(LEN=255), ALLOCATABLE  :: CollectionSubsetDims(:)
  CHARACTER(LEN=255), ALLOCATABLE  :: CollectionFormat    (:)
  CHARACTER(LEN=255), ALLOCATABLE  :: CollectionFrequency (:)
  CHARACTER(LEN=255), ALLOCATABLE  :: CollectionDuration  (:)
  CHARACTER(LEN=255), ALLOCATABLE  :: CollectionMode      (:)

  ! Objects
  TYPE(MetaHistContainer), POINTER :: CollectionList

CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: History_Init
!
! !DESCRIPTION: Reads the HISTORY.rc file and creates the linked list of
!  collections (i.e. netCDF diagnostic files containing several data fields
!  with a specified archival frequency).  The list of fields belonging to
!  each collection is also determined.  
!\\
!\\
!  Each collection is described by a HISTORY CONTAINER object, which also
!  contains a linked list of diagnostic quantities (i.e. a METAHISTORY ITEM)
!  that will be archived to netCDF format.  The list of diagnostic quantities
!  is determined here by parsing the HISTORY.rc file.
!\\
!\\
!  NOTE: The HISTORY.rc file is read twice.  The first (done by method
!  History\_ReadCollectionNames) reads the list of all collections.  Then,
!  for each defined collection, the list of diagnostic quantities belonging
!  to that collection is determined by routine History\_ReadCollectionData.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE History_Init( am_I_root,  Input_Opt, State_Chm,  &
                           State_Diag, State_Met, RC         )
!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod,  ONLY : OptInput
    USE State_Chm_Mod , ONLY : ChmState
    USE State_Diag_Mod, ONLY : DgnState
    USE State_Met_Mod,  ONLY : MetState
!
! !INPUT PARAMETERS: 
!
    LOGICAL,          INTENT(IN)  :: am_I_Root
    TYPE(OptInput),   INTENT(IN)  :: Input_Opt
    TYPE(ChmState),   INTENT(IN)  :: State_Chm
    TYPE(DgnState),   INTENT(IN)  :: State_Diag
    TYPE(MetState),   INTENT(IN)  :: State_Met
!
! !OUTPUT PARAMETERS: 
!
    INTEGER,          INTENT(OUT) :: RC
!
! !REMARKS:
!  Calls internal routines History_ReadCollectionNames and
!  History_Read_CollectionData
!
! !REVISION HISTORY:
!  06 Jan 2015 - R. Yantosca - Initial version
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
     ' -> at History_Init (in module History/history_mod.F90)'

    !=======================================================================
    ! First initialize the list of collections
    ! ("collection" = a netCDF file with a specific archival frequency)
    !=======================================================================
    CALL History_ReadCollectionNames( am_I_root,  Input_Opt, State_Chm,  &
                                      State_Diag, State_Met, RC         )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered in "History_ReadCollectionNames"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !=======================================================================
    ! Then determine the fields that will be saved to each collection
    !=======================================================================
    CALL History_ReadCollectionData( am_I_root,  Input_Opt, State_Chm,  &
                                     State_Diag, State_Met, RC         )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered in "History_ReadCollectionNames"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

  END SUBROUTINE History_Init
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Read_Collection_Names
!
! !DESCRIPTION: Reads the History input file (e.g. HISTORY.rc) and determines
!  the names of each individual diagnostic collection.  It stores this 
!  information in module variables for use in the next step.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE History_ReadCollectionNames( am_I_Root,  Input_Opt, State_Chm,  &
                                          State_Diag, State_Met, RC         )
!
! !USES:
!
    USE Charpak_Mod
    USE ErrCode_Mod
    USE Input_Opt_Mod,  ONLY : OptInput
    USE InquireMod,     ONLY : FindFreeLun
    USE State_Chm_Mod , ONLY : ChmState
    USE State_Diag_Mod, ONLY : DgnState
    USE State_Met_Mod,  ONLY : MetState
!
! !INPUT PARAMETERS: 
!
    LOGICAL,          INTENT(IN)  :: am_I_Root    ! Are we on the root CPU?
    TYPE(OptInput),   INTENT(IN)  :: Input_Opt    ! Input Options object
    TYPE(ChmState),   INTENT(IN)  :: State_Chm    ! Chemistry State object
    TYPE(DgnState),   INTENT(IN)  :: State_Diag   ! Diagnostic State object
    TYPE(MetState),   INTENT(IN)  :: State_Met    ! Meteorology State object
!
! !OUTPUT PARAMETERS: 
!
    INTEGER,          INTENT(OUT) :: RC           ! Success or failure?
!
! !REMARKS:
!  Called from routine History_Init.
!
! !REVISION HISTORY:
!  16 Jun 2017 - R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
    ! Ridiculously big number
    INTEGER, PARAMETER :: MAX_COLLECTIONS = 500
!
! !LOCAL VARIABLES:
!
    ! Scalars
    LOGICAL            :: EOF
    INTEGER            :: fId,    IOS      
    INTEGER            :: N,      nSubs1,  nSubs2

    ! Strings
    CHARACTER(LEN=255) :: ErrMsg, ThisLoc, Line,  Line2
    
    ! String arrays
    CHARACTER(LEN=255) :: Subs1(255)
    CHARACTER(LEN=255) :: Subs2(255)
    CHARACTER(LEN=255) :: TmpCollectionName(MAX_COLLECTIONS)

    !=======================================================================
    ! Initialize
    !=======================================================================

    ! Assume success
    RC      = GC_SUCCESS

    ! For error output
    ErrMsg  = ''
    ThisLoc = &
     ' -> at History_ReadCollectionNames (in module History/history_mod.F90)'

    ! Zero global variables
    CollectionCount   = 0
    TmpCollectionName = ''

    !=======================================================================
    ! Open the file containing the list of requested diagnostics
    !=======================================================================

    ! Find a free file unit
    fId     = FindFreeLun()

    ! Open the file
    OPEN( fId, FILE=TRIM(Input_Opt%HistoryInputFile), STATUS='OLD', IOSTAT=RC )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Could not open "' // TRIM(Input_Opt%HistoryInputFile) // '"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !=======================================================================
    ! Read data from the file
    !=======================================================================
    DO

       ! Read a single line, and strip leading/trailing spaces
       Line = ReadOneLine( fId, EOF, IOS, Squeeze=.TRUE. )

       ! Exit the loop if it's the end of the file
       IF ( EOF ) EXIT

       ! If it's a real I/O error, quit w/ error message
       IF ( IOS > 0 ) THEN
          ErrMsg = 'Unexpected end-of-file in "'       // &
                    TRIM( Input_Opt%HistoryInputFile ) // '"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       ! Skip if the line is commented out
       IF ( Line(1:1) == "#" ) CYCLE
      
       !-------------------------------------------------------------------
       ! Get the list of collections
       !-------------------------------------------------------------------
       IF ( Line(1:12) == 'COLLECTIONS:' ) THEN

          ! Start
          CollectionCount                    = CollectionCount + 1
          TmpCollectionName(CollectionCount) = CleanText( Line(14:) )
          
          ! Loop over all collections
          DO

             ! Read next line and strip leading/trailing blanks
             Line2 = ReadOneLine( fId, EOF, IOS, Squeeze=.TRUE. )

             ! Exit the loop if it's the end of the file
             IF ( EOF ) EXIT

             ! If it's a real I/O error, quit w/ error message
             IF ( IOS > 0 ) THEN
                ErrMsg = 'Unexpected error in "'        // &
                     TRIM( Input_Opt%HistoryInputFile ) // '"!'
                CALL GC_Error( ErrMsg, RC, ThisLoc )
                RETURN
             ENDIF

             ! Skip if the line is commented out
             IF ( Line2(1:1) == "#" ) CYCLE
      
             ! 2 colons signal the end of the collecto
             IF ( Line2(1:2) == '::' ) EXIT

             CollectionCount                    = CollectionCount + 1
             TmpCollectionName(CollectionCount) = CleanText( Line2 )
          ENDDO
       ENDIF
    ENDDO

    !=======================================================================
    ! Cleanup and quit
    !=======================================================================

    ! Close the file
    CLOSE( fId )

    ! Now that we now the number of diagnostic collections, we can allocate
    ! the array that will hold the name of each diagnostic collection.
    IF ( .not. ALLOCATED( CollectionName ) ) THEN
       ALLOCATE( CollectionName( CollectionCount ), STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN 
          ErrMsg = 'Could not allocate "CollectionName"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDIF

    ! Copy the collection names from the temporary array
    DO N = 1, CollectionCount
       CollectionName(N) = TmpCollectionName(N)
       WRITE( 6, '(i3, ":", a )' ) N, TRIM( CollectionName(N) )
    ENDDO
    CALL FLUSH( 6 )

    ! Allocate CollectionTemplate
    IF ( .not. ALLOCATED( CollectionTemplate ) ) THEN
       ALLOCATE( CollectionTemplate( CollectionCount ), STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN 
          ErrMsg = 'Could not allocate "CollectionTemplate"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
       CollectionTemplate = ''
    ENDIF

    ! Allocate CollectionFormat
    IF ( .not. ALLOCATED( CollectionFormat ) ) THEN
       ALLOCATE( CollectionFormat( CollectionCount ), STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN 
          ErrMsg = 'Could not allocate "CollectionFormat"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
       CollectionFormat = ''
    ENDIF

    ! Allocate CollectionFrequency
    IF ( .not. ALLOCATED( CollectionFrequency ) ) THEN
       ALLOCATE( CollectionFrequency( CollectionCount ), STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN 
          ErrMsg = 'Could not allocate "CollectionFrequency"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
       CollectionFrequency = ''
    ENDIF

    ! Allocate CollectionDuration
    IF ( .not. ALLOCATED( CollectionDuration ) ) THEN
       ALLOCATE( CollectionDuration( CollectionCount ), STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN 
          ErrMsg = 'Could not allocate "CollectionDuration"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
       CollectionDuration = ''
    ENDIF

    ! Allocate CollectionSubsetDims
    IF ( .not. ALLOCATED( CollectionSubsetDims ) ) THEN
       ALLOCATE( CollectionSubSetDims( CollectionCount ), STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN 
          ErrMsg = 'Could not allocate "CollectionSubsetDims"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
       CollectionSubsetDims = ''
    ENDIF

    ! Allocate CollectionMode
    IF ( .not. ALLOCATED( CollectionMode ) ) THEN
       ALLOCATE( CollectionMode( CollectionCount ), STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN 
          ErrMsg = 'Could not allocate "CollectionMode"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
       CollectionMode = ''
    ENDIF

  END SUBROUTINE History_ReadCollectionNames
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Read_Collection_Data
!
! !DESCRIPTION: Parses the History input file (e.g. HISTORY.rc) and compiles
!  the list of diagnostic quantities belonging to each collection.  In other
!  words, this is the list of individual fields that will be archived to a
!  particular netCDF file with a given archival frequency.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE History_ReadCollectionData( am_I_Root,  Input_Opt, State_Chm,  &
                                         State_Diag, State_Met, RC         )
!
! !USES:
!
    USE Charpak_Mod
    USE ErrCode_Mod
    USE HistContainer_Mod
    USE HistItem_Mod
    USE Input_Opt_Mod,         ONLY : OptInput
    USE InquireMod,            ONLY : FindFreeLun
    USE MetaHistContainer_Mod
    USE MetaHistItem_Mod
    USE Species_Mod,           ONLY : Species
    USE State_Chm_Mod
    USE State_Diag_Mod
    USE State_Met_Mod
!
! !INPUT PARAMETERS: 
!
    LOGICAL,          INTENT(IN)  :: am_I_Root    ! Are we on the root CPU?
    TYPE(OptInput),   INTENT(IN)  :: Input_Opt    ! Input Options object
    TYPE(ChmState),   INTENT(IN)  :: State_Chm    ! Chemistry State object
    TYPE(DgnState),   INTENT(IN)  :: State_Diag   ! Diagnostic State object
    TYPE(MetState),   INTENT(IN)  :: State_Met
!
! !OUTPUT PARAMETERS: 
!
    INTEGER,          INTENT(OUT) :: RC           ! Success or failure?
!
! !REMARKS:
!  Called from routine History_Init.
!
! !REVISION HISTORY:
!  16 Jun 2017 - R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
     ! Scalars
    LOGICAL                      :: EOF   
    INTEGER                      :: C,           N,           W    
    INTEGER                      :: fId,         IOS      
    INTEGER                      :: nSubs1,      nSubs2
    INTEGER                      :: Ind1,        Ind2
    INTEGER                      :: ArchivalYmd, ArchivalHms
    INTEGER                      :: ItemCount,   SpaceDim
    INTEGER                      :: Ind_All,     Ind_Adv,     Ind_Aer
    INTEGER                      :: Ind_Dry,     Ind_Fix,     Ind_Gas
    INTEGER                      :: Ind_Kpp,     Ind_Pho,     Ind_Rst
    INTEGER                      :: Ind_Var,     Ind_Wet,     Ind

    ! Strings
    CHARACTER(LEN=255)           :: Line,        FileName
    CHARACTER(LEN=255)           :: ErrMsg,      ThisLoc
    CHARACTER(LEN=255)           :: MetaData,    Reference
    CHARACTER(LEN=255)           :: Title,       Units
    CHARACTER(LEN=255)           :: ItemName,    ItemTemplate
    CHARACTER(LEN=255)           :: Description
    
    ! Arrays
    INTEGER                      :: SubsetDims(3)
    CHARACTER(LEN=255)           :: Subs1(255)
    CHARACTER(LEN=255)           :: Subs2(255)

    ! Objects
    TYPE(HistContainer), POINTER :: Collection
    TYPE(HistItem),      POINTER :: Item
    TYPE(Species),       POINTER :: ThisSpc

    ! Pointer arrays
    REAL(fp),            POINTER :: Ptr3d  (:,:,:)
    REAL(f4),            POINTER :: Ptr3d_4(:,:,:)

    !=======================================================================
    ! Initialize
    !=======================================================================

    ! Assume success
    RC          =  GC_SUCCESS

    ! Initialize variables
    ArchivalYmd =  0 
    ArchivalHms =  0
    SpaceDim    =  0
    SubsetDims  =  0 

    ! Initialize objects and pointers
    Collection  => NULL()
    Item        => NULL()
    Ptr3d       => NULL()
    Ptr3d_4     => NULL()
    ThisSpc     => NULL()

    ! Initialize Strings
    Description =  ''
    ErrMsg      =  ''
    Reference   =  'www.geos-chem.org; wiki.geos-chem.org'
    ThisLoc     =  &
     ' -> at History_ReadCollectionData (in module History/history_mod.F90)'
    Units       =  ''

    !=======================================================================
    ! Open the file containing the list of requested diagnostics
    !=======================================================================

    ! Find a free file unit
    fId     = FindFreeLun()

    ! Open the file
    OPEN( fId, FILE=TRIM(Input_Opt%HistoryInputFile), STATUS='OLD', IOSTAT=RC )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Could not open "' //TRIM(Input_Opt%HistoryInputFile) // '"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !=======================================================================
    ! Read data from the file
    !=======================================================================
    DO

       ! Read a single line, and strip leading/trailing spaces
       Line = ReadOneLine( fId, EOF, IOS, Squeeze=.TRUE. )

       ! Exit the loop if it's the end of the file
       IF ( EOF ) EXIT

       ! If it's a real I/O error, quit w/ error message
       IF ( IOS > 0 ) THEN
          ErrMsg = 'Unexpected end-of-file in "'       // &
                    TRIM( Input_Opt%HistoryInputFile ) // '"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       ! Skip if the line is commented out
       IF ( Line(1:1) == "#" ) CYCLE

       !====================================================================
       ! The HISTORY.rc file specifies collection metadata as:
       !
       !   instantaneous.template:  '%y4%m2%d2.nc4',
       !   instantaneous.format:    'CFIO',
       !   instantaneous.frequency:  010000,
       !   instantaneous.duration:   240000
       !   etc.
       !
       ! where in this example, "instantaneous" is the collection name
       ! and "template", "format", "frequency", "duration" are the
       ! metadata fields.  
       !
       ! Get the metadata belonging to each collection and store them
       ! in the proper arrays for later use.  NOTE: this method does not
       ! assume that the collections are in the same order as they
       ! are listed under the COLLECTIONS section.
       !====================================================================
       IF ( INDEX( Line, 'template' ) > 0 ) THEN 
          CALL GetCollectionMetaData( Line, 'template',   MetaData, C )
          CollectionTemplate(C) = Metadata
       ENDIF

       IF ( INDEX( Line, 'format' ) > 0 ) THEN
          CALL GetCollectionMetaData( Line, 'format',     MetaData, C )
          CollectionFormat(C) = Metadata
       ENDIF

       IF ( INDEX( Line, 'frequency' ) > 0 ) THEN
          CALL GetCollectionMetaData( Line, 'frequency',  MetaData, C )
          CollectionFrequency(C) = Metadata
          READ( CollectionFrequency(C), '(i6.6)' ) ArchivalHms
       ENDIF

       IF ( INDEX( Line, 'duration' ) > 0 ) THEN
          CALL GetCollectionMetaData( Line, 'format',     MetaData, C )
          CollectionDuration(C) = Metadata
       ENDIF

       IF ( INDEX( Line, 'mode' ) > 0 ) THEN
          CALL GetCollectionMetaData( Line, 'mode',       MetaData, C )
          CollectionMode(C) = Metadata
       ENDIF

       IF ( INDEX( Line, 'subsetdims' ) > 0 ) THEN
          
          ! First split the line by colon
          CALL StrSplit( Line, ":", Subs1, nSubs1 )
          CollectionSubsetDims(C) = Subs1(2)

          ! Then split by spaces and convert to INTEGER
          CALL StrSplit( CollectionSubsetDims(C), " ", Subs2, nSubs2 )          
          IF ( nSubs2 > 0 ) THEN
             DO N = 1, nSubs2
                READ( Subs2(N), '(i10)' ) SubsetDims(N)
             ENDDO

             ! Define the number of dimensions
             SpaceDim = nSubs2
             IF ( SpaceDim == 2 ) SubsetDims(3) = 1
          ENDIF
       ENDIF

       !====================================================================
       ! NOTE: We assume FIELDS is the last metadata tag for the
       ! collection.  We need to create the collection object and
       ! an object for each history item stored in the collection.
       !====================================================================
       IF ( INDEX( Line, 'fields' ) > 0 ) THEN

          ! Zero the counter of items
          ItemCount = 0

          ! Create title string for collection
          Title = 'GEOS-Chem diagnostic collection: ' //                    &
                   TRIM( CollectionName(C) )
          
          !=================================================================
          ! Create a HISTORY CONTAINER object for this collection
          !=================================================================
          CALL HistContainer_Create( am_I_Root    = am_I_Root,              &
                                     Container    = Collection,             &
                                     Name         = CollectionName(C),      &
                                     Id           = C,                      &
                                     ArchivalYmd  = ArchivalYmd,            &
                                     ArchivalHms  = ArchivalHms,            &
                                     Conventions  = 'COARDS',               &
                                     FileTemplate = CollectionTemplate(C),  & 
                                     NcFormat     = 'netCDF-4',             &
                                     Reference    = Reference,              &
                                     Title        = Title,                  &
                                     RC           = RC                     )

          ! Trap potential error
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = 'Could not create Collection: ' // &
                      TRIM( CollectionName(C) ) 
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF
         
          !=================================================================
          ! Create a list of HISTORY ITEMS that will be contained in this
          ! collection, for each entry under the "fields" tag.
          !=================================================================
          DO 

             IF ( ItemCount == 0 ) THEN

                !----------------------------------------------------------
                ! If we are on the same line as the "fields" tag, then
                ! the name of the HISTORY ITEM will be the first substring
                ! of MetaData (split on spaces).
                !----------------------------------------------------------
                CALL GetCollectionMetaData( Line, 'fields', MetaData, C )
                CALL StrSplit( MetaData, " ", Subs1, nSubs1 )
                ItemName = Subs1(1)
             
             ELSE

                !----------------------------------------------------------
                ! Otherwise, read the next line to get the name for
                ! each subsequent HISTORY ITEM.  The name will be the
                ! first substring of the line.
                !----------------------------------------------------------

                ! Read a single line, and strip leading/trailing spaces
                Line = ReadOneLine( fId, EOF, IOS, Squeeze=.TRUE. )

                ! IF we have hit the end of file then 
                iF ( EOF ) GOTO 999

                ! If it's a real I/O error, quit w/ error message
                IF ( IOS > 0 ) THEN
                   ErrMsg = 'Unexpected end-of-file in '        // &
                             TRIM( Input_Opt%HistoryInputFile ) //'!'
                   CALL GC_Error( ErrMsg, RC, ThisLoc )
                   RETURN
                ENDIF
                             
                ! Remove commas, spaces, and tabs
                Line = CleanText( Line )    

                ! Skip if the line is commented out
                IF ( Line(1:1) == "#"  ) CYCLE

                ! 2 colons denotes the end of the "fields" list
                IF ( Line(1:2) == '::' ) EXIT

                ! The HISTORY ITEM name will be the 1st substring of Line
                CALL StrSplit( Line, " ", Subs1, nSubs1 )
                ItemName = Subs1(1)

             ENDIF

             !--------------------------------------------------------------
             ! Create the
             !--------------------------------------------------------------

             ! Save a temporary 
             ItemTemplate = ItemName

             ! Test if there are wild cards present, otherwise skip
             IF ( INDEX( ItemTemplate, '?' ) >  0 ) THEN 

                ! Test for further wild cards
                Ind_Adv = INDEX( ItemTemplate, '?ADV?' ) ! Advected species
                Ind_All = INDEX( ItemTemplate, '?ALL?' ) ! All species
                Ind_Aer = INDEX( ItemTemplate, '?AER?' ) ! Aerosol species
                Ind_Dry = INDEX( ItemTemplate, '?DRY?' ) ! Drydep species
                Ind_Fix = INDEX( ItemTemplate, '?FIX?' ) ! KPP fixed species
                Ind_Gas = INDEX( ItemTemplate, '?GAS?' ) ! Gas-phase species
                Ind_Kpp = INDEX( ItemTemplate, '?KPP?' ) ! KPP species
                Ind_Pho = INDEX( ItemTemplate, '?PHO?' ) ! Photolysis species
                Ind_Var = INDEX( ItemTemplate, '?VAR?' ) ! KPP active species
                Ind_Wet = INDEX( ItemTemplate, '?WET?' ) ! Wetdep species

                ! Loop over all species
                DO N = 1, State_Chm%nSpecies

                   ! Point to this entry of the species database 
                   ThisSpc => State_Chm%SpcData(N)%Info

                   ! Skip to next entry if necessary
                   IF ( Ind_All > 0 ) THEN
                      Ind = Ind_All
                   ELSE IF ( Ind_Adv > 0 ) THEN
                      IF ( .not. ThisSpc%Is_Advected   ) CYCLE
                      Ind = Ind_Adv
                   ELSE IF ( Ind_Aer > 0 ) THEN
                      IF ( ThisSpc%Is_Gas              ) CYCLE
                      Ind = Ind_Aer
                   ELSE IF ( Ind_Dry > 0 ) THEN
                      IF ( .not. ThisSpc%Is_DryDep     ) CYCLE
                      Ind = Ind_Dry
                   ELSE IF ( Ind_Fix > 0 ) THEN
                      IF ( .not. ThisSpc%Is_FixedChem  ) CYCLE
                      Ind = Ind_Fix
                   ELSE IF ( Ind_Gas > 0 ) THEN
                      IF ( .not. ThisSpc%Is_Gas        ) CYCLE
                      Ind = Ind_Gas
                   ELSE IF ( Ind_Kpp > 0 ) THEN
                      IF ( .not. ThisSpc%Is_Kpp        ) CYCLE
                      Ind = Ind_Kpp
                   ELSE IF ( Ind_Pho > 0 ) THEN
                      IF ( .not. ThisSpc%Is_Photolysis ) CYCLE
                      Ind = Ind_Pho
                   ELSE IF ( Ind_Var > 0 ) THEN
                      IF ( .not. ThisSpc%Is_ActiveChem ) CYCLE
                      Ind = Ind_Var
                   ELSE IF ( Ind_Wet > 0 ) THEN
                      IF ( .not. ThisSpc%Is_WetDep     ) CYCLE
                      Ind = Ind_Wet
                   ELSE
                      Ind = -1
                   ENDIF

                   IF ( Ind <= 0 ) THEN
                      ErrMsg = 'Could not find wild card!'
                      CALL GC_Error( ErrMsg, RC, ThisLoc )
                      RETURN
                   ENDIF 

                   ! Construct the item name
                   ItemName    = ItemTemplate( 1:Ind-1 ) // &
                                 TRIM( ThisSpc%Name    ) // &
                                 ItemTemplate( Ind+5:  ) 

                   ! Increment the item count
                   ItemCount   = ItemCount + 1

                   !----------------------------------------------------------
                   ! Get a pointer to the data
                   !--------------------------------------------------------
                   IF ( INDEX( ItemName, State_Chm%State ) > 0 ) THEN
                      CALL Lookup_State_Chm( am_I_Root   = am_I_Root,        &
                                             State_Chm   = State_Chm,        &
                                             Variable    = ItemName,         &
                                             Description = Description,      &
                                             Units       = Units,            &
                                             Ptr3d       = Ptr3d,            &
                                             Ptr3d_4     = Ptr3d_4,          &
                                             RC          = RC               )

                      ! Trap potential error
                      IF ( RC /= GC_SUCCESS ) THEN
                         ErrMsg = 'Error in "Lookup_State_Chm" for ' // &
                                  ' diagnostic ' // TRIM( ItemName )
                         CALL GC_Error( ErrMsg, RC, ThisLoc )
                         RETURN
                      ENDIF
                   ENDIF

                   !--------------------------------------------------------
                   ! After all this setup, create the HISTORY ITEM object 
                   !--------------------------------------------------------
                   CALL HistItem_Create( am_I_Root   = am_I_Root,            &
                                         Item        = Item,                 &
                                         Id          = ItemCount,            & 
                                         ContainerId = C,                    &
                                         Name        = ItemName,             &
                                         LongName    = Description,          &
                                         Units       = Units,                &
                                         SpaceDim    = SpaceDim,             &
                                         NX          = SubsetDims(1),        &
                                         NY          = SubsetDims(2),        &
                                         NZ          = SubsetDims(3),        & 
                                         RC          = RC                   )


                   ! Trap potential error
                   IF ( RC /= GC_SUCCESS ) THEN
                      ErrMsg = 'Could not create Item: ' // TRIM( ItemName )
                      CALL GC_Error( ErrMsg, RC, ThisLoc )
                      RETURN
                   ENDIF
                   print*, '### after create!'
                   call flush(6)

                   !--------------------------------------------------------
                   ! Attach this HISTORY ITEM to the METAHISTORY ITEM (aka 
                   ! list of HISTORY ITEMS) belonging to the HISTORY 
                   ! CONTAINER object for the given collection.
                   !
                   ! In other words, this denotes the list of fields that 
                   ! will be written to the netCDF file, and with the 
                   ! archiving frequency, specified by this given diagnostic 
                   ! collection. 
                   !--------------------------------------------------------
                   CALL MetaHistItem_AddNew( am_I_Root = am_I_Root,          &
                                             Node = Collection%HistItems,    &
                                             Item      = Item,               &
                                             RC        = RC                 )

                   ! Trap potential error
                   IF ( RC /= GC_SUCCESS ) THEN
                      ErrMsg = 'Could not add Item' // TRIM( ItemName ) //  &
                           ' to ' // TRIM( CollectionName(C) ) // '%HistItems!'
                      CALL GC_Error( ErrMsg, RC, ThisLoc )
                      RETURN
                   ENDIF

                   ! Free the species database pointer
                   ThisSpc => NULL()

                ENDDO

             ELSE
             
                !=========================
                ! No wild cards present
                !==========================

                ! Increment the number of HISTORY items
                ItemCount = ItemCount + 1

                 !--------------------------------------------------------------
                 ! After all this setup, create the HISTORY ITEM object itself
                 !--------------------------------------------------------------
                CALL HistItem_Create( am_I_Root   = am_I_Root,               &
                                      Item        = Item,                    &
                                      Id          = ItemCount,               &
                                      ContainerId = C,                       &
                                      Name        = ItemName,                &
                                      LongName    = ItemName,                &
                                      Units       = 'TBD',                   &
                                      SpaceDim    = SpaceDim,                &
                                      NX          = SubsetDims(1),           &
                                      NY          = SubsetDims(2),           &
                                      NZ          = SubsetDims(3),           & 
                                      RC          = RC                      )


                ! Trap potential error
                IF ( RC /= GC_SUCCESS ) THEN
                   ErrMsg = 'Could not create Item: ' // TRIM( ItemName )
                   CALL GC_Error( ErrMsg, RC, ThisLoc )
                   RETURN
                ENDIF

             !--------------------------------------------------------------
             ! Attach this HISTORY ITEM to the METAHISTORY ITEM (aka list
             ! of HISTORY ITEMS) belonging to the HISTORY CONTAINER object
             ! for the given collection.
             !
             ! In other words, this denotes the list of fields that will
             ! be written to the netCDF file, and with the archiving
             ! frequency, specified by this given diagnostic collection.
             !--------------------------------------------------------------
                CALL MetaHistItem_AddNew( am_I_Root = am_I_Root,             &
                                          Node      = Collection%HistItems,  &
                                          Item      = Item,                  &
                                          RC        = RC                    )

                ! Trap potential error
                IF ( RC /= GC_SUCCESS ) THEN
                   ErrMsg = 'Could not add Item' // TRIM( ItemName ) //      &
                     ' to ' // TRIM( CollectionName(C) ) // '%HistItems!'
                   CALL GC_Error( ErrMsg, RC, ThisLoc )
                   RETURN
                ENDIF

             ENDIF
          ENDDO

          !=================================================================
          ! Add this HISTORY CONTAINER object (i.e. this collection) into
          ! the METAHISTORY OBJECT (i.e. the master list of collections).
          !=================================================================
          CALL MetaHistContainer_AddNew( am_I_Root   = am_I_Root,            &
                                         Node        = CollectionList,       &
                                         Container   = Collection,           &
                                         RC          = RC                   )

          ! Trap potential error
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = 'Could not add Container' //                           &
                      TRIM( CollectionName(C) ) //                           &
                      ' to the list of collections!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF
       ENDIF

    ENDDO
       
    !=======================================================================
    ! Cleanup and quit
    !=======================================================================
999 CONTINUE

    ! Close the file
    CLOSE( fId )

    DO C = 1, CollectionCount
       print*, 'Collection       ', TRIM( CollectionName      (C) )
       print*, '  -> Template    ', TRIM( CollectionTemplate  (C) )
       print*, '  -> Format      ', TRIM( CollectionFormat    (C) )
       print*, '  -> Frequency   ', TRIM( CollectionFrequency (C) )
       print*, '  -> Duration    ', TRIM( CollectionDuration  (C) )
       print*, '  -> Subset Dims ', TRIM( CollectionSubsetDims(C) )
       print*, '  -> Mode        ', TRIM( CollectionMode      (C) )
    ENDDO

    ! Print information about each diagnostic collection
    CALL MetaHistContainer_Print( am_I_Root, CollectionList, RC )

  END SUBROUTINE History_ReadCollectionData
!EOC
!!------------------------------------------------------------------------------
!!                  GEOS-Chem Global Chemical Transport Model                  !
!!------------------------------------------------------------------------------
!!BOP
!!
!! !IROUTINE: 
!!
!! !DESCRIPTION: 
!!\\
!!\\
!! !INTERFACE:
!!
!  SUBROUTINE History_AddItemToCollection( )
!!
!! !USES:
!!
!    USE Charpak_Mod
!    USE ErrCode_Mod
!    USE HistContainer_Mod
!    USE HistItem_Mod
!    USE Input_Opt_Mod,         ONLY : OptInput
!    USE MetaHistContainer_Mod
!    USE MetaHistItem_Mod
!    USE State_Chm_Mod
!    USE State_Diag_Mod
!    USE State_Met_Mod
!!
!! !INPUT PARAMETERS: 
!!
!    LOGICAL,             INTENT(IN)  :: am_I_Root     ! Are we on the root CPU?
!    TYPE(OptInput),      INTENT(IN)  :: Input_Opt     ! Input Options object
!    TYPE(ChmState),      INTENT(IN)  :: State_Chm     ! Chemistry State object
!    TYPE(DgnState),      INTENT(IN)  :: State_Diag    ! Diagnostic State object
!    TYPE(MetState),      INTENT(IN)  :: State_Met     ! Meteorology State object
!    CHARACTER(LEN=255),  INTENT(IN)  :: ItemName      ! Name of HISTORY ITEM 
!    INTEGER,             INTENT(IN)  :: ItemCount     ! Index of HISTORY ITEM
!    INTEGER,             INTENT(IN)  :: CollectionID  ! Collection ID number
!    INTEGER,             INTENT(IN)  :: SpaceDim      ! Number of dimensions
!    INTEGER,             INTENT(IN)  :: NX            ! # of longitudes
!    INTEGER,             INTENT(IN)  :: NY            ! # of latitudes
!    INTEGER,             INTENT(IN)  :: NZ            ! # of levels
!!
!! !INPUT/OUTPUT PARAMETERS:
!!
!    TYPE(HistContainer), POINTER     :: Collection    ! 
!!
!! !OUTPUT PARAMETERS: 
!!
!    INTEGER,             INTENT(OUT) :: RC            ! Success or failure?
!!
!! !RETURN VALUE:
!!
!!
!! !REMARKS:
!!
!! !REVISION HISTORY:
!!  06 Jan 2015 - R. Yantosca - Initial version
!!EOP
!!------------------------------------------------------------------------------
!!BOC
!!
!! !LOCAL VARIABLES:
!!  Called from routine History_Init.
!!
!! !REVISION HISTORY:
!!  16 Jun 2017 - R. Yantosca - Initial version
!!EOP
!!------------------------------------------------------------------------------
!!BOC
!!
!! !LOCAL VARIABLES:
!!
!     ! Scalars
!    LOGICAL                      :: EOF   
!    INTEGER                      :: C,           N,           W    
!    INTEGER                      :: fId,         IOS      
!    INTEGER                      :: nSubs1,      nSubs2
!    INTEGER                      :: Ind1,        Ind2
!    INTEGER                      :: ArchivalYmd, ArchivalHms
!    INTEGER                      :: ItemCount,   SpaceDim
!
!    ! Strings
!    CHARACTER(LEN=255)           :: Line,        FileName
!    CHARACTER(LEN=255)           :: ErrMsg,      ThisLoc
!    CHARACTER(LEN=255)           :: MetaData,    Reference
!    CHARACTER(LEN=255)           :: Title,       Units
!    CHARACTER(LEN=255)           :: ItemName,    ItemTemplate
!    CHARACTER(LEN=255)           :: Description
!    
!    ! Arrays
!    INTEGER                      :: SubsetDims(3)
!    CHARACTER(LEN=255)           :: Subs1(255)
!    CHARACTER(LEN=255)           :: Subs2(255)
!
!    ! Objects
!    TYPE(HistItem),      POINTER :: Item
!
!    ! Pointer arrays
!    REAL(fp),            POINTER :: Ptr2d  (:,:  )
!    REAL(f4),            POINTER :: Ptr2d_4(:,:  )
!    REAL(fp),            POINTER :: Ptr3d  (:,:,:)
!    REAL(f4),            POINTER :: Ptr3d_4(:,:,:)
!
!    !=======================================================================
!    ! Initialize
!    !=======================================================================
!
!
!
!    !=======================================================================
!    ! For each HISTORY ITEM, find the matching entry in the relevant
!    ! registry (in State_Chm, State_Diag, State_Met) and get a pointer
!    ! to the data source 
!    !=======================================================================
!    IF ( INDEX( ItemName, State_Chm%State ) > 0 ) THEN
!
!       !-------------------------
!       ! Chemistry State
!       !-------------------------
!       CALL Lookup_State_Chm(  am_I_Root   = am_I_Root,                      &
!                               State_Chm   = State_Chm,                      &
!                               Variable    = ItemName,                       &
!                               Description = Description,                    &
!                               Units       = Units,                          &
!                               Ptr3d       = Ptr3d,                          &
!                               Ptr3d_4     = Ptr3d_4,                        &
!                               RC          = RC                             )
!
!       ! Trap potential error
!       IF ( RC /= GC_SUCCESS ) THEN
!          ErrMsg = 'Error in "Lookup_State_Chm" for diagnostic ' //          &
!                   TRIM( ItemName )
!          CALL GC_Error( ErrMsg, RC, ThisLoc )
!          RETURN
!       ENDIF
!
!    ELSE IF ( INDEX( ItemName, State_Diag%State ) > 0 ) THEN
!
!       !-------------------------
!       ! Diagnostic State 
!       !-------------------------
!       CALL Lookup_State_Diag( am_I_Root   = am_I_Root,                      &
!                               State_Diag  = State_Met,                      &
!                               Variable    = ItemName,                       &
!                               Description = Description,                    &
!                               Units       = Units,                          &
!                               Ptr3d       = Ptr3d,                          &
!                               Ptr3d_4     = Ptr3d_4,                        &
!                               RC          = RC                             )
!
!       ! Trap potential error
!       IF ( RC /= GC_SUCCESS ) THEN
!          ErrMsg = 'Error in "Lookup_State_Diag for diagnostic ' //          &
!                   TRIM( ItemName )
!          CALL GC_Error( ErrMsg, RC, ThisLoc )
!          RETURN
!       ENDIF
!
!    ELSE IF ( INDEX( ItemName, State_Met%State ) > 0 ) THEN
!
!       !-------------------------
!       ! Meteorology State
!       !-------------------------
!       CALL Lookup_State_Met(  am_I_Root   = am_I_Root,                      &
!                               State_Met   = State_Met,                      &
!                               Variable    = ItemName,                       &
!                               Description = Description,                    &
!                               Units       = Units,                          &
!                               Ptr2d       = Ptr2d,                          &
!                               Ptr2d_4     = Ptr2d_4,                        &
!                               Ptr3d       = Ptr3d,                          &
!                               Ptr3d_4     = Ptr3d_4,                        &
!                               RC          = RC                             )
!
!       ! Trap potential error
!       IF ( RC /= GC_SUCCESS ) THEN
!          ErrMsg = 'Error in "Lookup_State_Met for diagnostic ' //           &
!                   TRIM( ItemName )
!          CALL GC_Error( ErrMsg, RC, ThisLoc )
!          RETURN
!       ENDIF
!
!
!
!    ENDIF
!
!    !=======================================================================
!    ! After all this setup, create the HISTORY ITEM object 
!    !=======================================================================
!    CALL HistItem_Create( am_I_Root   = am_I_Root,                           &
!                          Item        = Item,                                &
!                          Id          = ItemCount,                           & 
!                          ContainerId = CollectionId,                        &
!                          Name        = ItemName,                            &
!                          LongName    = Description,                         &
!                          Units       = Units,                               &
!                          SpaceDim    = SpaceDim,                            &
!                          NX          = NX,          ,                       &
!                          NY          = NY,                                  &
!                          NZ          = NZ,                                  & 
!                          RC          = RC                                  )
!
!
!    ! Trap potential error
!    IF ( RC /= GC_SUCCESS ) THEN
!       ErrMsg = 'Could not create Item: ' // TRIM( ItemName )
!       CALL GC_Error( ErrMsg, RC, ThisLoc )
!       RETURN
!    ENDIF
!
!    !=======================================================================
!    ! Attach this HISTORY ITEM to the METAHISTORY ITEM (aka 
!    ! list of HISTORY ITEMS) belonging to the HISTORY 
!    ! CONTAINER object for the given collection.
!    !
!    ! In other words, this denotes the list of fields that 
!    ! will be written to the netCDF file, and with the 
!    ! archiving frequency, specified by this given diagnostic 
!    ! collection. 
!    !=======================================================================
!    CALL MetaHistItem_AddNew( am_I_Root = am_I_Root,                         &
!                              Node      = Collection%HistItems,              &
!                              Item      = Item,                              &
!                              RC        = RC                                )
!
!    ! Trap potential error
!    IF ( RC /= GC_SUCCESS ) THEN
!       ErrMsg = 'Could not add Item' // TRIM( ItemName ) //  &
!            ' to ' // TRIM( CollectionName(C) ) // '%HistItems!'
!       CALL GC_Error( ErrMsg, RC, ThisLoc )
!       RETURN
!    ENDIF
!
!  END SUBROUTINE History_AddItemToCollection
!                 
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
    CHARACTER(LEN=255)      :: CleanStr   ! Cleaned-up string
!
! !REMARKS:
!
! !REVISION HISTORY:
!  06 Jan 2015 - R. Yantosca - Initial version
!  21 Jun 2017 - R. Yantosca - Now call CSTRIP to remove tabs etc.
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
    USE Charpak_Mod, ONLY : StrSqueeze
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
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GetCollectionMetaData
!
! !DESCRIPTION: tbd
!\\
!\\
! !INTERFACE:
!
   SUBROUTINE GetCollectionMetaData( Line, Pattern, MetaData, nCollection )
!
! !USES:
!
     USE Charpak_Mod, ONLY: StrSplit
!
! !INPUT PARAMETERS: 
!
     CHARACTER(LEN=*),   INTENT(IN)  :: Line
     CHARACTER(LEN=*),   INTENT(IN)  :: Pattern
!
! !OUTPUT PARAMETERS:
!
     CHARACTER(LEN=255), INTENT(OUT) :: MetaData
     INTEGER,            INTENT(OUT) :: nCollection
!
! !REMARKS:
!
! !REVISION HISTORY:
!  06 Jan 2015 - R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
     INTEGER                 :: C, Ind, nSubStr
     CHARACTER(LEN=255) :: SubStr(255)
     
     !======================================================================
     ! Find the metadata for the given collection
     !======================================================================

     ! Loop over all collection names
     DO C = 1, CollectionCount

        ! Then check to see which collection this is in
        Ind = INDEX( Line, TRIM( CollectionName(C) ) )

        ! If we find a match, then return it
        IF ( Ind > 0 ) THEN
           CALL StrSplit( Line, ':', SubStr, nSubStr )
           IF ( nSubStr == 2 ) THEN
              nCollection = C
              MetaData    = CleanText( SubStr(2) )
              EXIT
           ENDIF
        ENDIF
     ENDDO

 END SUBROUTINE GetCollectionMetaData
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: History_Cleanup
!
! !DESCRIPTION: Deallocates all module variables and objects.
!\\
!\\
! !INTERFACE:
!
   SUBROUTINE History_Cleanup( am_I_Root, RC )
!
! !USES:
!
     USE ErrCode_Mod
!
! !INPUT PARAMETERS: 
!
     LOGICAL, INTENT(IN)  :: am_I_Root
!
! !OUTPUT PARAMETERS: 
!
     INTEGER, INTENT(OUT) :: RC
!
! !REVISION HISTORY:
!  16 Jun 2017 - R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
     ! Strings
     CHARACTER(LEN=255) :: ErrMsg, ThisLoc

     !======================================================================
     ! Initialize
     !======================================================================

     ! Assume success
     RC      = GC_SUCCESS

     ErrMsg  = ''
     ThisLoc = ' -> History_Cleanup (in module History/history_mod.F90'

     !======================================================================
     ! Deallocate module variables
     !======================================================================

     IF ( ASSOCIATED( CollectionList ) ) THEN
        DEALLOCATE( CollectionList, STAT=RC )
        IF ( RC /= GC_SUCCESS ) THEN
           ErrMsg = 'Could not deallocate "CollectionList"!'
           CALL GC_Error( ErrMsg, RC, ThisLoc )
           RETURN
        ENDIF        
     ENDIF

     IF ( ALLOCATED( CollectionName ) ) THEN
        DEALLOCATE( CollectionName, STAT=RC )
        IF ( RC /= GC_SUCCESS ) THEN
           ErrMsg = 'Could not deallocate "CollectionName"!'
           CALL GC_Error( ErrMsg, RC, ThisLoc )
           RETURN
        ENDIF
     ENDIF
     
     IF ( ALLOCATED( CollectionTemplate ) ) THEN
        DEALLOCATE( CollectionTemplate, STAT=RC )
        IF ( RC /= GC_SUCCESS ) THEN
           ErrMsg = 'Could not deallocate "CollectionTemplate"!'
           CALL GC_Error( ErrMsg, RC, ThisLoc )
           RETURN
        ENDIF
     ENDIF

     IF ( ALLOCATED( CollectionSubsetDims ) ) THEN
        DEALLOCATE( CollectionSubsetDims, STAT=RC )
        IF ( RC /= GC_SUCCESS ) THEN
           ErrMsg = 'Could not deallocate "CollectionSubsetDims"!'
           CALL GC_Error( ErrMsg, RC, ThisLoc )
           RETURN
        ENDIF
     ENDIF

     IF ( ALLOCATED( CollectionFormat ) ) THEN
        DEALLOCATE( CollectionFormat, STAT=RC )
        IF ( RC /= GC_SUCCESS ) THEN
           ErrMsg = 'Could not deallocate "CollectionFormat"!'
           CALL GC_Error( ErrMsg, RC, ThisLoc )
           RETURN
        ENDIF
     ENDIF

     IF ( ALLOCATED( CollectionFrequency ) ) THEN
        DEALLOCATE( CollectionFrequency, STAT=RC )
        IF ( RC /= GC_SUCCESS ) THEN
           ErrMsg = 'Could not deallocate "CollectionFrequency"!'
           CALL GC_Error( ErrMsg, RC, ThisLoc )
           RETURN
        ENDIF
     ENDIF

     IF ( ALLOCATED( CollectionDuration ) ) THEN
        DEALLOCATE( CollectionDuration, STAT=RC )
        IF ( RC /= GC_SUCCESS ) THEN
           ErrMsg = 'Could not deallocate "CollectionDuration"!'
           CALL GC_Error( ErrMsg, RC, ThisLoc )
           RETURN
        ENDIF
     ENDIF

     IF ( ALLOCATED( CollectionMode ) ) THEN
        DEALLOCATE( CollectionMode, STAT=RC )
        IF ( RC /= GC_SUCCESS ) THEN
           ErrMsg = 'Could not deallocate "CollectionMode"!'
           CALL GC_Error( ErrMsg, RC, ThisLoc )
           RETURN
        ENDIF
     ENDIF

   END SUBROUTINE History_Cleanup
!EOC
END MODULE History_Mod
