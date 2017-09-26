!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: history_mod.F90
!
! !DESCRIPTION: Driver module for GEOS-Chem's netCDF diagnostics package, aka
!  the "History Component".
!\\
!\\
! !INTERFACE:
!
MODULE History_Mod
!
! !USES:
!
  USE Precision_Mod
  USE HistContainer_Mod,     ONLY : HistContainer
  USE MetaHistContainer_Mod, ONLY : MetaHistContainer

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
! 
  PUBLIC  :: History_Init
  PUBLIC  :: History_SetTime
  PUBLIC  :: History_Update
  PUBLIC  :: History_Write
  PUBLIC  :: History_Close_AllFiles
  PUBLIC  :: History_Cleanup
!
! PRIVATE MEMBER FUNCTIONS:

  PRIVATE :: History_ReadCollectionNames
  PRIVATE :: History_ReadCollectionData
  PRIVATE :: History_AddItemToCollection
  PRIVATE :: CleanText
  PRIVATE :: ReadOneLine
!
! !REMARKS:
!  
!
! !REVISION HISTORY:
!  06 Jan 2015 - R. Yantosca - Initial version
!  02 Aug 2017 - R. Yantosca - Added History_Update routine
!  14 Aug 2017 - R. Yantosca - Now read the "acc_interval" field for
!                              time-averaged data collections
!  16 Aug 2017 - R. Yantosca - Add subroutine TestTimeForAction to avoid
!                              duplicating similar code
!  16 Aug 2017 - R. Yantosca - Now close all netCDF files in routine
!                              History_Close_AllFiles
!  18 Aug 2017 - R. Yantosca - Added routine History_SetTime
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
  ! Scalars
  INTEGER                          :: CollectionCount

  ! Strings
  CHARACTER(LEN=255), ALLOCATABLE  :: CollectionName       (:)
  CHARACTER(LEN=255), ALLOCATABLE  :: CollectionTemplate   (:)
  CHARACTER(LEN=255), ALLOCATABLE  :: CollectionSubsetDims (:)
  CHARACTER(LEN=255), ALLOCATABLE  :: CollectionFormat     (:)
  CHARACTER(LEN=255), ALLOCATABLE  :: CollectionFrequency  (:)
  CHARACTER(LEN=255), ALLOCATABLE  :: CollectionAccInterval(:)
  CHARACTER(LEN=255), ALLOCATABLE  :: CollectionDuration   (:)
  CHARACTER(LEN=255), ALLOCATABLE  :: CollectionMode       (:)

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
!  with a specified update frequency).  The list of fields belonging to
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
  SUBROUTINE History_Init( am_I_root, Input_Opt, State_Chm, State_Diag,      &
                           State_Met, yyyymmdd,  hhmmss,    RC              )
!
! !USES:
!
    USE ErrCode_Mod
    USE History_Netcdf_Mod, ONLY : History_Netcdf_Init
    USE History_Util_Mod
    USE Input_Opt_Mod,      ONLY : OptInput
    USE State_Chm_Mod ,     ONLY : ChmState
    USE State_Diag_Mod,     ONLY : DgnState
    USE State_Met_Mod,      ONLY : MetState
!
! !INPUT PARAMETERS: 
!
    LOGICAL,          INTENT(IN)  :: am_I_Root
    TYPE(OptInput),   INTENT(IN)  :: Input_Opt
    TYPE(ChmState),   INTENT(IN)  :: State_Chm
    TYPE(DgnState),   INTENT(IN)  :: State_Diag
    TYPE(MetState),   INTENT(IN)  :: State_Met
    INTEGER,          INTENT(IN)  :: yyyymmdd
    INTEGER,          INTENT(IN)  :: hhmmss
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
    ! Initialize the history_netcdf_mod.F90 module
    !=======================================================================
    CALL History_Netcdf_Init( am_I_Root, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered in "HistoryNetCdfInit"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !=======================================================================
    ! First initialize the list of collections
    ! ("collection" = a netCDF file with a specific update frequency)
    !=======================================================================
    CALL History_ReadCollectionNames( am_I_root,  Input_Opt, State_Chm,      &
                                      State_Diag, State_Met, RC             )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered in "History_ReadCollectionNames"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !=======================================================================
    ! Then determine the fields that will be saved to each collection
    !=======================================================================
    CALL History_ReadCollectionData( am_I_root,  Input_Opt, State_Chm,       &
                                     State_Diag, State_Met, yyyymmdd,        &
                                     hhmmss,     RC                         )
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
! !IROUTINE: History_Read_Collection_Names
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
    USE History_Util_Mod
    USE Input_Opt_Mod,     ONLY : OptInput
    USE InquireMod,        ONLY : FindFreeLun
    USE State_Chm_Mod ,    ONLY : ChmState
    USE State_Diag_Mod,    ONLY : DgnState
    USE State_Met_Mod,     ONLY : MetState
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
!  Private routine, called from routine History_Init.
!
! !REVISION HISTORY:
!  16 Jun 2017 - R. Yantosca - Initial version
!  15 Aug 2017 - R. Yantosca - Now initialize string arrays to UNDEFINED_STR
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
    ENDDO

    ! Allocate CollectionTemplate
    IF ( .not. ALLOCATED( CollectionTemplate ) ) THEN
       ALLOCATE( CollectionTemplate( CollectionCount ), STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN 
          ErrMsg = 'Could not allocate "CollectionTemplate"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
       CollectionTemplate = UNDEFINED_STR
    ENDIF

    ! Allocate CollectionFormat
    IF ( .not. ALLOCATED( CollectionFormat ) ) THEN
       ALLOCATE( CollectionFormat( CollectionCount ), STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN 
          ErrMsg = 'Could not allocate "CollectionFormat"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
       CollectionFormat = UNDEFINED_STR
    ENDIF

    ! Allocate CollectionFrequency
    IF ( .not. ALLOCATED( CollectionFrequency ) ) THEN
       ALLOCATE( CollectionFrequency( CollectionCount ), STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN 
          ErrMsg = 'Could not allocate "CollectionFrequency"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
       CollectionFrequency = UNDEFINED_STR
    ENDIF

    ! Allocate CollectionAccInterval
    IF ( .not. ALLOCATED( CollectionAccInterval ) ) THEN
       ALLOCATE( CollectionAccInterval( CollectionCount ), STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN 
          ErrMsg = 'Could not allocate "CollectionAccInterval"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
       CollectionAccInterval = UNDEFINED_STR
    ENDIF

    ! Allocate CollectionDuration
    IF ( .not. ALLOCATED( CollectionDuration ) ) THEN
       ALLOCATE( CollectionDuration( CollectionCount ), STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN 
          ErrMsg = 'Could not allocate "CollectionDuration"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
       CollectionDuration = UNDEFINED_STR
    ENDIF

    ! Allocate CollectionSubsetDims
    IF ( .not. ALLOCATED( CollectionSubsetDims ) ) THEN
       ALLOCATE( CollectionSubSetDims( CollectionCount ), STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN 
          ErrMsg = 'Could not allocate "CollectionSubsetDims"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
       CollectionSubsetDims = UNDEFINED_STR
    ENDIF

    ! Allocate CollectionMode
    IF ( .not. ALLOCATED( CollectionMode ) ) THEN
       ALLOCATE( CollectionMode( CollectionCount ), STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN 
          ErrMsg = 'Could not allocate "CollectionMode"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
       CollectionMode = UNDEFINED_STR
    ENDIF

  END SUBROUTINE History_ReadCollectionNames
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: History_Read_Collection_Data
!
! !DESCRIPTION: Parses the History input file (e.g. HISTORY.rc) and compiles
!  the list of diagnostic quantities belonging to each collection.  In other
!  words, this is the list of individual fields that will be archived to a
!  particular netCDF file with a given update frequency.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE History_ReadCollectionData( am_I_Root,  Input_Opt, State_Chm,  &
                                         State_Diag, State_Met, yyyymmdd,   &
                                         hhmmss,     RC                    )
!
! !USES:
!
    USE Charpak_Mod
    USE CMN_Size_Mod,          ONLY : IIPAR, JJPAR, LLPAR
    USE ErrCode_Mod
    USE HistContainer_Mod
    USE HistItem_Mod
    USE History_Util_Mod
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
    TYPE(MetState),   INTENT(IN)  :: State_Met    ! Meteorology State object
    INTEGER,          INTENT(IN)  :: yyyymmdd     ! Current date in YMD
    INTEGER,          INTENT(IN)  :: hhmmss       ! Current time in hms
!
! !OUTPUT PARAMETERS: 
!
    INTEGER,          INTENT(OUT) :: RC           ! Success or failure?
!
! !REMARKS:
!  Private routine, called from History_Init.
!
! !REVISION HISTORY:
!  16 Jun 2017 - R. Yantosca - Initial version
!  03 Aug 2017 - R. Yantosca - Pass OPERATION to History_AddItemToCollection
!  14 Aug 2017 - R. Yantosca - FileWrite{Ymd,Hms} and FileClose{Ymd,Hms} are
!                              now computed properly, w/r/t acc_interval
!  30 Aug 2017 - R. Yantosca - Now write collection info only on the root CPU
!  18 Sep 2017 - R. Yantosca - Don't allow acc_interval for inst collections
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
     ! Scalars
    LOGICAL                      :: EOF   
    INTEGER                      :: C,              N,             W
    INTEGER                      :: nX,             nY,            nZ
    INTEGER                      :: fId,            IOS
    INTEGER                      :: nSubs1,         nSubs2
    INTEGER                      :: Ind1,           Ind2
    INTEGER                      :: UpdateYmd,      UpdateHms    
    INTEGER                      :: FileCloseYmd,   FileCloseHms
    INTEGER                      :: FileWriteYmd,   FileWriteHms
    INTEGER                      :: ItemCount,      SpaceDim,      Operation
    INTEGER                      :: Ind_All,        Ind_Adv,       Ind_Aer
    INTEGER                      :: Ind_Dry,        Ind_Fix,       Ind_Gas
    INTEGER                      :: Ind_Kpp,        Ind_Pho,       Ind_Rst
    INTEGER                      :: Ind_Var,        Ind_Wet,       Ind
    REAL(f8)                     :: UpdateAlarm,    JulianDate
    REAL(f8)                     :: FileWriteAlarm, FileCloseAlarm
    REAL(f8)                     :: HeartBeatDtSec

    ! Strings
    CHARACTER(LEN=255)           :: Line,           FileName
    CHARACTER(LEN=255)           :: ErrMsg,         ThisLoc
    CHARACTER(LEN=255)           :: MetaData,       Reference
    CHARACTER(LEN=255)           :: Title,          Units
    CHARACTER(LEN=255)           :: ItemName,       ItemTemplate
    CHARACTER(LEN=255)           :: Description,    TmpMode
    CHARACTER(LEN=255)           :: Contact,        Pattern

    ! Arrays
    INTEGER                      :: SubsetDims(3)
    CHARACTER(LEN=255)           :: Subs1(255)
    CHARACTER(LEN=255)           :: Subs2(255)

    ! Objects
    TYPE(HistContainer), POINTER :: Container
    TYPE(HistItem),      POINTER :: Item
    TYPE(Species),       POINTER :: ThisSpc

    ! Pointer arrays
    REAL(fp),            POINTER :: Ptr3d  (:,:,:)
    REAL(f4),            POINTER :: Ptr3d_4(:,:,:)

    !=======================================================================
    ! Initialize
    !=======================================================================

    ! Assume success
    RC             =  GC_SUCCESS

    ! Initialize variables
    UpdateYmd      =  0 
    UpdateHms      =  0
    FileCloseYmd   =  0
    FileCloseHms   =  0
    FileWriteYmd   =  0
    FileWriteHms   =  0
    SpaceDim       =  0
    SubsetDims     =  0 
    HeartBeatDtSec =  DBLE( Input_Opt%TS_DYN ) * SECONDS_PER_MINUTE

    ! Initialize objects and pointers
    Container      => NULL()
    Item           => NULL()
    Ptr3d          => NULL()
    Ptr3d_4        => NULL()
    ThisSpc        => NULL()

    ! Initialize Strings
    Description    =  ''
    ErrMsg         =  ''
    Contact        =  &
     'GEOS-Chem Support Team (geos-chem-support@as.harvard.edu)'
    Reference      =  'www.geos-chem.org; wiki.geos-chem.org'
    ThisLoc        =  &
     ' -> at History_ReadCollectionData (in module History/history_mod.F90)'
    Units          =  ''

    !=======================================================================
    ! Open the file containing the list of requested diagnostics
    !=======================================================================

    ! Find a free file unit
    fId     = FindFreeLun()

    ! Open the filei
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

       ! "template": Specifies the year/month/day/hr/min/sec in filenames
       Pattern = 'template'
       IF ( INDEX( TRIM( Line ), TRIM( Pattern ) ) > 0 ) THEN 
          CALL GetCollectionMetaData( Line, Pattern, MetaData, C )
          IF ( C > 0 ) CollectionTemplate(C) = Metadata
       ENDIF       

       ! "format": Specifies the file output format (e.g. netCDF-4, CFIO)
       Pattern = 'format'
       IF ( INDEX( TRIM( Line ), TRIM( Pattern ) ) > 0 ) THEN
          CALL GetCollectionMetaData( Line, Pattern, MetaData, C )
          IF ( C > 0 ) CollectionFormat(C) = Metadata
       ENDIF

       ! "frequency": Specifies how often diagnostics are updated
       Pattern = 'frequency'
       IF ( INDEX( TRIM( Line ), TRIM( Pattern ) ) > 0 ) THEN
          CALL GetCollectionMetaData( Line, Pattern, MetaData, C )
          IF ( C > 0 ) CollectionFrequency(C) = Metadata
       ENDIF

       ! "acc_interval": Specifies how often time-averaged diagnostics
       ! are updated.  Defaults to "frequncy"
       Pattern = 'acc_interval'
       IF ( INDEX( TRIM( Line ), TRIM( Pattern ) ) > 0 ) THEN
          CALL GetCollectionMetaData( Line, Pattern, MetaData, C )
          IF ( C > 0 ) CollectionAccInterval(C) = Metadata
       ENDIF

       ! "duration:: Specifies how often files will be written
       ! Defaults to "frequency"
       Pattern = 'duration'
       IF ( INDEX( TRIM( Line ), TRIM( Pattern ) ) > 0 ) THEN
          CALL GetCollectionMetaData( Line, Pattern, MetaData, C )
          IF ( C > 0 ) CollectionDuration(C) = Metadata
       ENDIF

       ! "mode": Specifies instantaneous or time-averaged archiving
       Pattern = 'mode'
       IF ( INDEX( TRIM( Line ), TRIM( Pattern ) ) > 0 ) THEN
          CALL GetCollectionMetaData( Line, Pattern, MetaData, C )
          IF ( C > 0 ) CollectionMode(C) = Metadata
       ENDIF

       ! "subsetdims": Specifies a subset of the data grid
       ! NOTE: Currently not used at the present time
       Pattern = 'subsetdims'
       IF ( INDEX( TRIM( Line ), TRIM( Pattern ) ) > 0 ) THEN
          
          ! First split the line by colon
          CALL StrSplit( Line, ":", Subs1, nSubs1 )
          IF ( C > 0 ) THEN
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
       ENDIF

       !====================================================================
       ! NOTE: We assume FIELDS is the last metadata tag for the
       ! collection.  We need to create the collection object and
       ! an object for each history item stored in the collection.
       !====================================================================
       IF ( INDEX( TRIM( Line ), 'fields' ) > 0 ) THEN

          ! If we can't find the metadata for the collection in HISTORY.rc,
          ! then this might point to a mismatch between names under the
          ! "COLLECTIONS:" list and the corresponding metadata section.
          ! Quit with error if this occurs.
          IF ( C == UNDEFINED_INT ) THEN
             ErrMsg = 'Mismatch between listed collection names and the ' // &
                      'metadata!  Please check HISTORY.rc for typos.'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF

          ! Zero the counter of items
          ItemCount = 0

          ! Create title string for collection
          Title = 'GEOS-Chem diagnostic collection: ' //                     &
                   TRIM( CollectionName(C) )

          !-----------------------------------------------------------------
          ! Determine the operation code (i.e. copy or accumulate from the 
          ! source pointer to Item's data array for further analysis), 
          ! based on the value of CollectionMode.
          !-----------------------------------------------------------------
          TmpMode = CollectionMode(C)
          CALL TranUc( TmpMode )
          SELECT CASE( TmpMode )
             CASE( 'TIME-AVERAGED', 'TIMEAVERAGED' )
                Operation = ACCUM_FROM_SOURCE
             CASE DEFAULT
                Operation = COPY_FROM_SOURCE

                ! Throw an error if the "acc_interval" is defined,
                ! but the collection is instantaneous.
                IF ( .not. TRIM( CollectionAccInterval(C) ) ==               &
                                 UNDEFINED_STR                 ) THEN
                   ErrMsg = 'Acc_interval cannot be defined for ' //         &
                            'instantaneous collection: "'         //         &
                            TRIM( CollectionName(C) )             // '"!' 
                   CALL GC_Error( ErrMsg, RC, ThisLoc )
                   RETURN
                ENDIF

          END SELECT

          !-----------------------------------------------------------------
          ! Define UpdateHms and UpdateYmd
          ! 
          ! NOTE: If CollectionFrequency is 6 digits long, then assume that
          ! to be UpdateHms.  If longer, then assume that it is specifying
          ! both UpdateYmd and UpdateHms. (sde, bmy, 8/4/17)
          !-----------------------------------------------------------------
          IF ( LEN_TRIM( CollectionFrequency(C) ) == 6 ) THEN
             READ( CollectionFrequency(C), '(i6.6)'  ) UpdateHms
          ELSE
             READ( CollectionFrequency(C), '(2i6.6)' ) UpdateYmd,            &
                                                       UpdateHms
          ENDIF

          ! SPECIAL CASE: If UpdateHms is 240000 then set
          ! and set UpdateYmd=000001 and UpdateHms=000000
          IF ( UpdateHms == 240000 ) THEN
             UpdateYmd = 000001
             UpdateHms = 000000
          ENDIF

          !-----------------------------------------------------------------
          ! Define FileWriteYmd and FileWriteHms
          !-----------------------------------------------------------------
          IF ( Operation == COPY_FROM_SOURCE ) THEN

             ! Instantaneous data: this defaults to the frequency attribute
             ! in the HISTORY.rc file, which is UpdateYmd, UpdateHms
             FileWriteYmd = UpdateYmd
             FileWriteHms = UpdateHms

          ELSE
             
             ! Time-averaged data: Use the acc_interval field to denote
             ! when it is time to write to disk.  If acc_interval is not
             ! specified, then default to UpdateYmd and UpdateHms,
             ! which are specified by the "frequency" attribute.
             IF ( TRIM( CollectionAccInterval(C) ) == UNDEFINED_STR ) THEN
                FileWriteYmd = UpdateYmd
                FileWriteHms = UpdateHms
             ELSE IF ( LEN_TRIM( CollectionAccInterval(C) ) == 6 ) THEN
                READ( CollectionAccInterval(C), '(i6.6)'  ) FileWriteHms
             ELSE
                READ( CollectionAccInterval(C), '(2i6.6)' ) FileWriteYmd,    &
                                                            FileWriteHms
             ENDIF
          ENDIF
          
          ! SPECIAL CASE: If FileWriteYmd is 240000 then set
          ! and set FileWriteYmd=000001 and FileWriteHms=000000
          IF ( FileWriteHms == 240000 ) THEN
             FileWriteYmd = 000001
             FileWriteHms = 000000
          ENDIF
   
          !-----------------------------------------------------------------
          ! Define FileCloseHms and FileCloseYmd
          !
          ! NOTE: If CollectionDuration is 6 digits long, then assume that
          ! to be UpdateHms.  If longer, then assume that it is specifying
          ! both FileWriteYmd and FileWriteHms. (sde, bmy, 8/4/17)
          !
          ! ALSO NOTE: If "duration" is not found, then default to th
          ! same values specified by the "frequency" attribute.
          !-----------------------------------------------------------------
          IF ( TRIM( CollectionDuration(C) ) == UNDEFINED_STR ) THEN
             FileCloseYmd = UpdateYmd
             FileCloseHms = UpdateHms
          ELSE IF ( LEN_TRIM( CollectionDuration(C) ) == 6 ) THEN
             READ( CollectionDuration(C), '(i6.6)'  ) FileCloseHms
          ELSE
             READ( CollectionDuration(C), '(2i6.6)' ) FileCloseYmd,          &
                                                      FileCloseHms
          ENDIF

          ! SPECIAL CASE: If FileCloseHms is 240000 then set
          ! and set FileCloseYmd=000001 and FileCloseYmd=000000
          IF ( FileCloseHms == 240000 ) THEN
             FileCloseYmd = 000001
             FileCloseHms = 000000
          ENDIF
          
          !=================================================================
          ! Create a HISTORY CONTAINER object for this collection
          !=================================================================

          ! Compute the Astronomical Julian Date corresponding to 
          ! the initial yyyymmdd and hhmmss of the simulation.
          ! This is needed to set the EpochJd field.
          CALL Compute_Julian_Date( yyyymmdd, hhmmss, JulianDate )
          
          ! Create the HISTORY CONTAINER object itself.
          ! This will also define the alarm intervals and initial alarm times
          CALL HistContainer_Create( am_I_Root      = am_I_Root,             &
                                     Container      = Container,             &
                                     Id             = C,                     &
                                     Name           = CollectionName(C),     &
                                     EpochJd        = JulianDate,            &
                                     CurrentYmd     = yyyymmdd,              &
                                     CurrentHms     = hhmmss,                &
                                     UpdateMode     = CollectionMode(C),     &
                                     UpdateYmd      = UpdateYmd,             &
                                     UpdateHms      = UpdateHms,             &
                                     Operation      = Operation,             &
                                     HeartBeatDtSec = HeartBeatDtSec,        &
                                     FileWriteYmd   = FileWriteYmd,          &
                                     FileWriteHms   = FileWriteHms,          &
                                     FileCloseYmd   = FileCloseYmd,          &
                                     FileCloseHms   = FileCloseHms,          &
                                     Conventions    = 'COARDS',              &
                                     FileTemplate   = CollectionTemplate(C), & 
                                     NcFormat       = CollectionFormat(C),   &
                                     Reference      = Reference,             &
                                     Title          = Title,                 &
                                     Contact        = Contact,               &
                                     RC             = RC                    )

          ! Trap potential error
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = 'Could not create Collection: ' // &
                      TRIM( CollectionName(C) ) 
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF

          ! Set elapsed time quantities in the HISTORY CONTAINER object
          CALL HistContainer_SetTime( am_I_Root   = am_I_Root,               &
                                      Container   = Container,               &
                                      HeartBeatDt = 0.0_f8,                  &
                                      RC          = RC                      )


          ! Trap potential error
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = 'Error encountered in "HistContainer_SetTime"'   //    &
                      ' for collection: ' // TRIM( CollectionName(C) ) 
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
             ! Create the a HISTORY ITEM object for each diagnostic 
             ! entry read from HISTORY.rc and add to the given COLLECTION
             !--------------------------------------------------------------

             ! Save the item name in a temporary variable
             ! so that we can parse it for wild cards
             ItemTemplate = ItemName

             ! Test if there are wild cards present, otherwise skip
             IF ( INDEX( ItemTemplate, '?' ) >  0 ) THEN 

                !-----------------------------------------------------------
                ! Item name contains wild cards; fill in species names
                !-----------------------------------------------------------
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

                   ! For a given wild card, skip the unnecessary species
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

                   ! Create the a HISTORY ITEM object for this diagnostic
                   ! and add it to the given DIAGNOSTIC COLLECTION
                   CALL History_AddItemToCollection(                         &
                            am_I_Root    = am_I_Root,                        &
                            Input_Opt    = Input_Opt,                        &
                            State_Chm    = State_Chm,                        &
                            State_Diag   = State_Diag,                       &
                            State_Met    = State_Met,                        &
                            Collection   = Container,                        &
                            CollectionId = C,                                &
                           !SubsetDims   = CollectionSubsetDims(C),          &
                            ItemName     = ItemName,                         &
                            ItemCount    = ItemCount,                        &
                            RC           = RC                               )

                   ! Trap potential error
                   IF ( RC /= GC_SUCCESS ) THEN
                      ErrMsg = 'Could not add diagnostic :'            //     &
                               TRIM( ItemName ) // '" to collection: ' //    &
                               TRIM( CollectionName(C) ) 
                      CALL GC_Error( ErrMsg, RC, ThisLoc )
                      RETURN
                   ENDIF
         
                   ! Free the species database pointer
                   ThisSpc => NULL()
                ENDDO

             ELSE

                !-----------------------------------------------------------
                ! Item name does not have wild cards; no special handling
                !-----------------------------------------------------------

                ! Increment the number of HISTORY items
                ItemCount = ItemCount + 1

                ! Create the a HISTORY ITEM object for this diagnostic
                ! and add it to the given DIAGNOSTIC COLLECTION
                CALL History_AddItemToCollection(                            &
                         am_I_Root    = am_I_Root,                           &
                         Input_Opt    = Input_Opt,                           &
                         State_Chm    = State_Chm,                           &
                         State_Diag   = State_Diag,                          &
                         State_Met    = State_Met,                           &
                         Collection   = Container,                           &
                         CollectionId = C,                                   &
                        !SubsetDims   = CollectionSubsetDims(C),             &
                         ItemName     = ItemName,                            &
                         ItemCount    = ItemCount,                           &
                         RC           = RC                                  )

                ! Trap potential error
                IF ( RC /= GC_SUCCESS ) THEN
                   ErrMsg = 'Could not create add diagnostic :'     //       &
                            TRIM( ItemName ) // '" to collection: ' //       &
                            TRIM( CollectionName(C) ) 
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
                                         Container   = Container,            &
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

    ! Write output
    IF ( am_I_Root ) THEN
       WRITE( 6, '(/,a)' ) REPEAT( '=', 79 )
       WRITE( 6, '(a  )' ) 'DEFINED DIAGNOSTIC COLLECTIONS:'
       WRITE( 6, '(a  )' ) REPEAT( '=', 79 )

       DO C = 1, CollectionCount
          print*, 'Collection        ', TRIM( CollectionName       (C) )
          print*, '  -> Template     ', TRIM( CollectionTemplate   (C) )
          print*, '  -> Format       ', TRIM( CollectionFormat     (C) )
          print*, '  -> Frequency    ', TRIM( CollectionFrequency  (C) )
          print*, '  -> Acc_Interval ', TRIM( CollectionAccInterval(C) )
          print*, '  -> Duration     ', TRIM( CollectionDuration   (C) )
!       print*, '  -> Subset Dims  ', TRIM( CollectionSubsetDims (C) )
          print*, '  -> Mode         ', TRIM( CollectionMode       (C) )

          ! Trap error if the collection name is undefined
          ! This indicates an error in parsing the file
          IF ( TRIM( CollectionTemplate(C) ) == UNDEFINED_STR ) THEN
             ErrMsg = 'Collection: ' // TRIM( CollectionName(C) ) //         &
                      ' is undefined!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF
       ENDDO
    ENDIF

    ! Print information about each diagnostic collection
    CALL MetaHistContainer_Print( am_I_Root, CollectionList, RC )

    ! Write spacer
    WRITE( 6, '(a,/)' ) REPEAT( '=', 79 )   

  END SUBROUTINE History_ReadCollectionData
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: History_AddItemToCollection
!
! !DESCRIPTION: Creates a HISTORY ITEM object for a given diagnostic quantity,
!  and then attaches it to a given diagnostic collection.  Given the name
!  of the diagnostic quantity, it will obtain metadata (and pointers to the
!  data array) via the appropriate state registry.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE History_AddItemToCollection( am_I_Root,    Input_Opt,           &
                                          State_Chm,    State_Diag,          &
                                          State_Met,    Collection,          &
                                          CollectionId, ItemName,            &
                                          ItemCount,    SubsetDims,          &
                                          RC                                )
!
! !USES:
!
    USE Charpak_Mod,           ONLY : TranUc
    USE ErrCode_Mod
    USE HistContainer_Mod
    USE HistItem_Mod
    USE History_Util_Mod,      ONLY : UNDEFINED_INT
    USE Input_Opt_Mod,         ONLY : OptInput
    USE MetaHistContainer_Mod
    USE MetaHistItem_Mod
    USE State_Chm_Mod
    USE State_Diag_Mod
    USE State_Met_Mod
!
! !INPUT PARAMETERS: 
!
    ! Required arguments
    LOGICAL,             INTENT(IN)  :: am_I_Root      ! Are we on the root CPU?
    TYPE(OptInput),      INTENT(IN)  :: Input_Opt      ! Input Options State
    TYPE(ChmState),      INTENT(IN)  :: State_Chm      ! Chemistry State
    TYPE(DgnState),      INTENT(IN)  :: State_Diag     ! Diagnostic State
    TYPE(MetState),      INTENT(IN)  :: State_Met      ! Meteorology State
    INTEGER,             INTENT(IN)  :: CollectionID   ! Collection ID number
    CHARACTER(LEN=255),  INTENT(IN)  :: ItemName       ! Name of HISTORY ITEM 
    INTEGER,             INTENT(IN)  :: ItemCount      ! Index of HISTORY ITEM

    ! Optional arguments
    INTEGER,             OPTIONAL    :: SubsetDims(3)  ! Dimensions specified
                                                       !  by the collection
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(HistContainer), POINTER     :: Collection     ! Diagnostic Collection
!
! !OUTPUT PARAMETERS: 
!
    INTEGER,             INTENT(OUT) :: RC             ! Success or failure?
!
! !REMARKS:
!  Private routine, called from History_Init.
!
! !REVISION HISTORY:
!  06 Jan 2015 - R. Yantosca - Initial version
!  03 Aug 2017 - R. Yantosca - Inherit operation code from the Collection
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    LOGICAL                      :: OnLevelEdges
    INTEGER                      :: KindVal
    INTEGER                      :: Rank
    INTEGER                      :: NX, NY, NZ

    ! Arrays
    INTEGER                      :: Dimensions(3)
    INTEGER                      :: ItemDims(3)

    ! Strings
    CHARACTER(LEN=255)           :: Description
    CHARACTER(LEN=255)           :: ErrMsg
    CHARACTER(LEN=255)           :: ThisLoc
    CHARACTER(LEN=255)           :: Units

    ! Objects
    TYPE(HistItem),      POINTER :: Item

    ! Pointer arrays
    REAL(fp),            POINTER :: Ptr2d  (:,:  )
    REAL(f4),            POINTER :: Ptr2d_4(:,:  )
    INTEGER,             POINTER :: Ptr2d_I(:,:  )
    REAL(fp),            POINTER :: Ptr3d  (:,:,:)
    REAL(f4),            POINTER :: Ptr3d_4(:,:,:)
    INTEGER,             POINTER :: Ptr3d_I(:,:,:)

    !=======================================================================
    ! Initialize
    !=======================================================================
    RC          =  GC_SUCCESS
    Description =  ''
    Dimensions  =  0
    KindVal     =  0
    Rank        =  0
    Units       =  ''
    ErrMsg      =  ''
    ThisLoc     =  &
                ' -> History_AddItemToCollection (in History/history_mod.F90)'
    Ptr2d       => NULL()
    Ptr2d_4     => NULL()
    Ptr2d_I     => NULL()
    Ptr3d       => NULL()
    Ptr3d_4     => NULL()
    Ptr3d_I     => NULL()

    !=======================================================================
    ! For each HISTORY ITEM, find the matching entry in the relevant
    ! registry (in State_Chm, State_Diag, State_Met) and get a pointer
    ! to the data source 
    !=======================================================================
    IF ( INDEX( ItemName, State_Chm%State ) > 0 ) THEN

       !--------------------------------------------------------------------
       ! Chemistry State
       !--------------------------------------------------------------------
       CALL Lookup_State_Chm(  am_I_Root    = am_I_Root,                     &
                               State_Chm    = State_Chm,                     &
                               Variable     = ItemName,                      &
                               Description  = Description,                   &
                               Dimensions   = Dimensions,                    &
                               KindVal      = KindVal,                       &
                               Units        = Units,                         &
                               OnLevelEdges = OnLevelEdges,                  &
                               Rank         = Rank,                          &
                               Ptr3d        = Ptr3d,                         &
                               Ptr3d_4      = Ptr3d_4,                       &
                               RC           = RC                            )

       ! Trap potential error
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Error in "Lookup_State_Chm" for diagnostic ' //          &
                   TRIM( ItemName )
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

    ELSE IF ( INDEX( ItemName, State_Diag%State ) > 0 ) THEN

       !--------------------------------------------------------------------
       ! Diagnostic State 
       !--------------------------------------------------------------------
       CALL Lookup_State_Diag( am_I_Root    = am_I_Root,                     &
                               State_Diag   = State_Diag,                    &
                               Variable     = ItemName,                      &
                               Description  = Description,                   &
                               Dimensions   = Dimensions,                    &
                               KindVal      = KindVal,                       &
                               Rank         = Rank,                          &
                               Units        = Units,                         &
                               OnLevelEdges = OnLevelEdges,                  &
                               Ptr2d_4      = Ptr2d_4,                       &
                               Ptr3d_4      = Ptr3d_4,                       &
                               RC           = RC                            )

       ! Trap potential error
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Error in "Lookup_State_Diag for diagnostic ' //          &
                   TRIM( ItemName )
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

    ELSE IF ( INDEX( ItemName, State_Met%State ) > 0 ) THEN

       !--------------------------------------------------------------------
       ! Meteorology State
       !--------------------------------------------------------------------
       CALL Lookup_State_Met(  am_I_Root    = am_I_Root,                     &
                               State_Met    = State_Met,                     &
                               Variable     = ItemName,                      &
                               Description  = Description,                   &
                               Dimensions   = Dimensions,                    &
                               KindVal      = KindVal,                       &
                               Rank         = Rank,                          &
                               Units        = Units,                         &
                               OnLevelEdges = OnLevelEdges,                  &
                               Ptr2d        = Ptr2d,                         &
                               Ptr2d_I      = Ptr2d_I,                       &
                               Ptr3d        = Ptr3d,                         &
                               Ptr3d_I      = Ptr3d_I,                       &
                               RC           = RC                            )

       ! Trap potential error
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Error in "Lookup_State_Met for diagnostic ' //           &
                   TRIM( ItemName )
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF


    ELSE

       !--------------------------------------------------------------------
       ! ERROR! Item name could not be found in a registry!
       !-------------------------------------------------------------------- 
       ErrMsg = 'Could not locate ' // TRIM( ItemName )  // '. Try prefacing the state name ("CHEM_", "MET_", "DIAG_") to the item name.'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN

    ENDIF

    !=======================================================================
    ! If the optional SUBSETDIMS argument is passed, then use that to
    ! size the data arrays.  Otherwise assume the data arrays will be
    ! the same size as the pointer to the data source (as will be true
    ! in most cases.) 
    !=======================================================================
    IF ( PRESENT( SubsetDims ) ) THEN
       NX = SubsetDims(1)
       NY = SubsetDims(2)
       NZ = SubsetDims(3)
    ELSE
       NX = Dimensions(1)
       NY = Dimensions(2)
       NZ = Dimensions(3)
    ENDIF

    !=======================================================================
    ! Now that we have obtained information (and pointers to the data)
    ! corresponding to the given diagnostic quantity, use that to create
    ! a HISTORY ITEM object for that diagnostic quantity.
    !=======================================================================
    CALL HistItem_Create( am_I_Root      = am_I_Root,                        &
                          Item           = Item,                             &
                          Id             = ItemCount,                        & 
                          ContainerId    = CollectionId,                     &
                          Name           = ItemName,                         &
                          LongName       = Description,                      &
                          Units          = Units,                            &
                          OnLevelEdges   = OnLevelEdges,                     &
                          SpaceDim       = Rank,                             &
                          Operation      = Collection%Operation,             &
                          Source_KindVal = KindVal,                          &
                          Source_2d      = Ptr2d,                            &
                          Source_2d_4    = Ptr2d_4,                          &
                          Source_2d_I    = Ptr2d_I,                          &
                          Source_3d      = Ptr3d,                            &
                          Source_3d_4    = Ptr3d_4,                          &
                          Source_3d_I    = Ptr3d_I,                          &
                          Dimensions     = ItemDims,                         &
                          RC             = RC                               )

    ! Trap potential error
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Could not create Item: "' // TRIM( ItemName ) // '"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !=======================================================================
    ! Attach this HISTORY ITEM to the METAHISTORY ITEM (aka list of HISTORY 
    ! ITEMS) belonging to the HISTORY CONTAINER object for the given 
    ! diagnostic collection.
    !
    ! In other words, we are adding this diagnostic quantity to the list
    ! of diagnostic quantities that belong to this diagnostic collection.
    ! These quantities will be written to the netCDF file described by
    ! the collection, with the specified update frequency.
    !=======================================================================
    CALL MetaHistItem_AddNew( am_I_Root = am_I_Root,                         &
                              Node      = Collection%HistItems,              &
                              Item      = Item,                              &
                              RC        = RC                                )

    ! Trap potential error
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Could not add Item "' //                                    &
                TRIM( ItemName )       //  '" to '   //                      &
                TRIM( CollectionName(CollectionId) ) // '%HistItems!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !=======================================================================
    ! Define the dimensions of the collection from the HISTORY ITEMS,
    ! and also whether vertical data is on the level centers or edges
    !=======================================================================

    ! Define the X dimension of the collection from the 
    ! first HISTORY ITEM that has an X dimension
    IF ( Collection%NX == UNDEFINED_INT ) THEN
       SELECT CASE( Item%DimNames )
          CASE( 'xyz', 'xz', 'xy', 'x' )
             Collection%NX = ItemDims(1)
          CASE DEFAULT
             ! Nothing
       END SELECT
    ENDIF

    ! Define the Y dimension of the collection from the
    ! first HISTORY ITEM that has a Y dimension
    IF ( Collection%NY == UNDEFINED_INT ) THEN
       SELECT CASE( Item%DimNames )
          CASE( 'xyz', 'xy' )
             Collection%NY = ItemDims(2)
          CASE( 'yz', 'y' )
             Collection%NY = ItemDims(1)
          CASE DEFAULT
             ! Nothing
       END SELECT
    ENDIF

    ! Define the Z dimension of the collection from the first HISTORY ITEM 
    ! that has a Z dimension. Also define whether the collection will
    ! contain data that is centered or edged on vertical levels.
    IF ( Collection%NZ == UNDEFINED_INT ) THEN
       SELECT CASE( Item%DimNames )
          CASE( 'xyz' )
             Collection%NZ = ItemDims(3)
          CASE( 'xz', 'yz' )
             Collection%NZ = ItemDims(2)
          CASE( 'z' )
             Collection%NZ = ItemDims(1)
          CASE DEFAULT
             ! Nothing
       END SELECT

       Collection%OnLevelEdges = Item%OnLevelEdges
    ENDIF

    !=======================================================================
    ! Make sure that all the HISTORY ITEMS in this collection are 
    ! placed on the level centers or edges, but not both.  The netCDF
    ! COARDS/CF conventions do not allow for data on more than one 
    ! vertical dimension per file.
    !=======================================================================
    IF ( Item%SpaceDim == 3 ) THEN
       IF ( Collection%OnLevelEdges .neqv. Item%OnLevelEdges ) THEN
          ErrMsg = TRIM( Item%Name )                                    //   &
                   ' has the wrong vertical alignment for collection> ' //   &
                   TRIM( Collection%Name )
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDIF

    !=======================================================================
    ! Cleanup and quit
    !=======================================================================

    ! Increment number of HISTORY ITEMS in this collection
    Collection%nHistItems = Collection%nHistItems + 1

    ! Free pointers
    Ptr2d   => NULL()
    Ptr2d_4 => NULL()
    Ptr2d_I => NULL()
    Ptr3d   => NULL()
    Ptr3d_4 => NULL()
    Ptr3d_I => NULL()

  END SUBROUTINE History_AddItemToCollection
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: History_SetTime
!
! !DESCRIPTION: Sets the time values for each HISTORY CONTAINER object
!  that specifies a diagnostic collection. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE History_SetTime( am_I_Root, RC )
! 
! !USES:
!
    USE ErrCode_Mod
    USE HistContainer_Mod,     ONLY : HistContainer_SetTime
    USE History_Util_Mod
    USE MetaHistContainer_Mod, ONLY : MetaHistContainer
!
! !INPUT PARAMETERS: 
!
    LOGICAL,  INTENT(IN)  :: am_I_Root        ! Are we on the root CPU?
!
! !OUTPUT PARAMETERS: 
!
    INTEGER,  INTENT(OUT) :: RC               ! Success or failure
!
! !REMARKS:
!  This routine is meant to be called after History_Update() but before
!  History_Write().
!
! !REVISION HISTORY:
!  18 Aug 2017 - R. Yantosca - Initial version
!  29 Aug 2017 - R. Yantosca - Remove HeartBeatDtMin as an argument; now the
!                              Container object contains heartbeat timesteps
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Strings
    CHARACTER(LEN=255)               :: ErrMsg
    CHARACTER(LEN=255)               :: ThisLoc

    ! Objects
    TYPE(MetaHistContainer), POINTER :: Collection

    !=======================================================================
    ! Initialize
    !=======================================================================
    RC         =  GC_SUCCESS
    Collection => NULL()
    ErrMsg     =  ''
    ThisLoc    =  ' -> at History_SetTime (in History/history_mod.F90)' 

    !=======================================================================
    ! Loop through each DIAGNOSTIC COLLECTION in the master list
    !=======================================================================

    ! Point to the first COLLECTION in the master collection list
    Collection => CollectionList
    
    ! As long as this current COLLECTION is valid ...
    DO WHILE( ASSOCIATED( Collection ) ) 

       ! Update the time settings for the next timestep
       CALL HistContainer_SetTime( am_I_Root   = am_I_Root,                  &
                                   Container   = Collection%Container,       &
                                   RC          = RC                         )

       ! Trap error
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Error encountered in "HistContainer_SetTime" ' //        &
                   ' for container : ' // TRIM( Collection%Container%Name )
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       ! Go to next collection
       Collection => Collection%Next
    ENDDO

    !=======================================================================
    ! Cleanup and quit
    !=======================================================================

    ! Free pointers
    Collection => NULL()

  END SUBROUTINE History_SetTime
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: History_Update
!
! !DESCRIPTION: For each HISTORY ITEM belonging to a diagnostic COLLECTION,
!  the data from the target variable is copied or accumulated into the 
!  HISTORY ITEM's data field for further analysis.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE History_Update( am_I_Root, RC )
! 
! !USES:
!
    USE ErrCode_Mod
    USE HistItem_Mod,          ONLY : HistItem
    USE HistContainer_Mod,     ONLY : HistContainer
    USE History_Util_Mod
    USE MetaHistContainer_Mod, ONLY : MetaHistContainer
    USE MetaHistItem_Mod,      ONLY : MetaHistItem
    USE Registry_Params_Mod
!
! !INPUT PARAMETERS: 
!
    LOGICAL, INTENT(IN)  :: am_I_Root  ! Are we on the root CPU?
!
! !OUTPUT PARAMETERS: 
!
    INTEGER, INTENT(OUT) :: RC         ! Success or failure
!
! !REMARKS:
!  This routine is called from the main program at the end of each 
!  "heartbeat" timestep.
!
! !REVISION HISTORY:
!  03 Aug 2017 - R. Yantosca - Initial version
!  11 Aug 2017 - R. Yantosca - Remove references to 0d pointers, data arrays
!  16 Aug 2017 - R. Yantosca - Now call TestTimeForAction to test if it is
!                              time to update the diagnostic collection.
!  21 Aug 2017 - R. Yantosca - Now get yyyymmdd, hhmmss from the container
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    LOGICAL                          :: DoUpdate
    
    ! Strings
    CHARACTER(LEN=255)               :: ErrMsg
    CHARACTER(LEN=255)               :: ThisLoc

    ! Objects
    TYPE(MetaHistContainer), POINTER :: Collection
    TYPE(HistContainer),     POINTER :: Container
    TYPE(MetaHistItem),      POINTER :: Current
    TYPE(HistItem),          POINTER :: Item

    !=======================================================================
    ! Initialize
    !=======================================================================
    RC         =  GC_SUCCESS
    DoUpdate   = .FALSE.
    Collection => NULL()
    Container  => NULL()
    Current    => NULL()
    Item       => NULL()
    ErrMsg     =  ''
    ThisLoc    =  ' -> at History_Update (in History/history_mod.F90)' 

    !=======================================================================
    ! Loop through each DIAGNOSTIC COLLECTION in the master list, and
    ! then loop through the HISTORY ITEMS belonnging to each COLLECTION.
    ! Update each HISTORY ITEM if it is the proper time.
    !=======================================================================

    ! Point to the first COLLECTION in the master collection list
    Collection => CollectionList
    
    ! As long as this current COLLECTION is valid ...
    DO WHILE( ASSOCIATED( Collection ) ) 

       ! Point to the HISTORY CONTAINER object in this COLLECTION
       Container => Collection%Container

       !--------------------------------------------------------------------
       ! Test if it is time to update this collection
       !--------------------------------------------------------------------

       ! Test if the "UpdateAlarm" is ringing
       DoUpdate = ( ( Container%UpdateAlarm - Container%ElapsedSec ) < EPS )

       ! Skip to next collection if it isn't
       IF ( .not. DoUpdate ) THEN
          Container  => NULL()
          Collection => Collection%Next
          CYCLE
       ENDIF

#if defined( DEBUG )
       ! Debug output
       IF ( am_I_Root ) THEN
          WRITE( 6, 100 ) Container%Name
 100      FORMAT( '     - Updating collection: ', a20 ) 
       ENDIF
#endif
       
       !--------------------------------------------------------------------
       ! If it is time to update the collection, then loop through all of
       ! the associated HISTORY ITEMS and either copy or accumulate the
       ! data from the source pointer into the HISTORY ITEM's data array.
       !--------------------------------------------------------------------

       ! Point to the first HISTORY ITEM belonging to the 
       ! HISTORY CONTAINER object for the current COLLECTION
       Current => Container%HistItems

       ! As long as this HISTORY ITEM is valid ...
       DO WHILE ( ASSOCIATED( Current ) )

          ! Get the HISTORY ITEM object contained in this node
          ! of the linked list of HISTORY ITEMS for this COLLECTION
          Item => Current%Item

          ! Test the rank of the data
          SELECT CASE( Item%SpaceDim ) 

             !--------------------------------------------------------------
             ! Update 3-D data field
             !--------------------------------------------------------------
             CASE( 3 )

                ! Flex-precision floating point
                IF ( Item%Source_KindVal == KINDVAL_FP ) THEN

                   IF ( Item%Operation == COPY_FROM_SOURCE ) THEN
                      Item%Data_3d  = Item%Source_3d
                      Item%nUpdates = 1.0_f8
                   ELSE
                      Item%Data_3d  = Item%Data_3d  + Item%Source_3d
                      Item%nUpdates = Item%nUpdates + 1.0_f8
                   ENDIF

                ! 4-byte floating point
                ELSE IF ( Item%Source_KindVal == KINDVAL_F4 ) THEN

                   IF ( Item%Operation == COPY_FROM_SOURCE ) THEN
                      Item%Data_3d  = Item%Source_3d_4
                      Item%nUpdates = 1.0_f8
                   ELSE
                      Item%Data_3d  = Item%Data_3d  + Item%Source_3d_4
                      Item%nUpdates = Item%nUpdates + 1.0_f8
                   ENDIF

                ! Integer
                ELSE IF ( Item%Source_KindVal == KINDVAL_I4 ) THEN

                   IF ( Item%Operation == COPY_FROM_SOURCE ) THEN
                      Item%Data_3d  = Item%Source_3d_I
                      Item%nUpdates = 1.0_f8
                   ELSE
                      Item%Data_3d  = Item%Data_3d  + Item%Source_3d_I
                      Item%nUpdates = Item%nUpdates + 1.0_f8
                   ENDIF

                ENDIF

             !--------------------------------------------------------------
             ! Update 2-D data field
             !--------------------------------------------------------------
             CASE( 2 )

                ! Flex-precision floating point
                IF ( Item%Source_KindVal == KINDVAL_FP ) THEN

                   IF ( Item%Operation == COPY_FROM_SOURCE ) THEN
                      Item%Data_2d  = Item%Source_2d
                      Item%nUpdates = 1.0_f8
                   ELSE 
                      Item%Data_2d  = Item%Data_2d  + Item%Source_2d
                      Item%nUpdates = Item%nUpdates + 1.0_f8
                   ENDIF

                ! 4-byte floating point
                ELSE IF ( Item%Source_KindVal == KINDVAL_F4 ) THEN

                   IF ( Item%Operation == COPY_FROM_SOURCE ) THEN
                      Item%Data_2d  = Item%Source_2d_4
                      Item%nUpdates = 1.0_f8
                   ELSE
                      Item%Data_2d  = Item%Data_2d + Item%Source_2d_4
                      Item%nUpdates = Item%nUpdates + 1.0_f8 
                   ENDIF

                ! Integer
                ELSE IF ( Item%Source_KindVal == KINDVAL_I4 ) THEN

                   IF ( Item%Operation == COPY_FROM_SOURCE ) THEN
                      Item%Data_2d  = Item%Source_2d_I
                      Item%nUpdates = 1.0_f8
                   ELSE 
                      Item%Data_2d  = Item%Data_2d  + Item%Source_2d_I
                      Item%nUpdates = Item%nUpdates + 1.0_f8 
                   ENDIF

                ENDIF

             !--------------------------------------------------------------
             ! Update 1-D data field
             !--------------------------------------------------------------
             CASE( 1 )

                ! Flex-precision floating point
                IF ( Item%Source_KindVal == KINDVAL_FP ) THEN

                   IF ( Item%Operation == COPY_FROM_SOURCE ) THEN
                      Item%Data_1d  = Item%Source_1d
                      Item%nUpdates = 1.0_f8
                   ELSE 
                      Item%Data_1d  = Item%Data_1d  + Item%Source_1d
                      Item%nUpdates = Item%nUpdates + 1.0_f8 
                   ENDIF

                ! 4-byte floating point
                ELSE IF ( Item%Source_KindVal == KINDVAL_F4 ) THEN

                   IF ( Item%Operation == COPY_FROM_SOURCE ) THEN
                      Item%Data_1d  = Item%Source_1d_4
                      Item%nUpdates = 1.0_f8
                   ELSE 
                      Item%Data_1d  = Item%Data_1d  + Item%Source_1d_4
                      Item%nUpdates = Item%nUpdates + 1.0_f8  
                   ENDIF

                ! Integer
                ELSE IF ( Item%Source_KindVal == KINDVAL_I4 ) THEN

                   IF ( Item%Operation == COPY_FROM_SOURCE ) THEN
                      Item%Data_1d  = Item%Source_1d_I
                      Item%nUpdates = 1.0_f8 
                   ELSE
                      Item%Data_1d  = Item%Data_1d  + Item%Source_1d_I
                      Item%nUpdates = Item%nUpdates + 1.0_f8 
                   ENDIF

                ENDIF

          END SELECT

! Uncomment more detailed debug output if you need it!
!#if defined( DEBUG ) 
!          ! Debug output
!          IF ( am_I_Root ) THEN
!             WRITE( 6, 110 ) TRIM(Container%Name),                        &
!                             TRIM(Item%Name),       Item%nUpdates
! 110         FORMAT( a20, 1x, a20, 1x, f7.1 )
!          ENDIF
!#endif

          ! Free pointer
          Item => NULL()

          ! Go to the next HISTORY item
          Current => Current%Next
       ENDDO

       !------------------------------------------------------------------
       ! Prepare to go to the next collection
       !------------------------------------------------------------------ 

       ! Update the "UpdateAlarm" time for the next updating interval.
       Container%UpdateAlarm = Container%UpdateAlarm +                    &
                               Container%UpdateIvalSec

       ! Free pointers
       Current    => NULL()
       Container  => NULL()

       ! Go to the next entry in the list of collections
       Collection => Collection%Next
    ENDDO

    !=======================================================================   
    ! Cleanup and quit
    !=======================================================================

    ! Free pointers
    Collection => NULL()
    Container  => NULL()
    Current    => NULL()
    Item       => NULL()

  END SUBROUTINE History_Update
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: History_Write
!
! !DESCRIPTION: For each HISTORY ITEM belonging to a diagnostic COLLECTION,
!  the data from the target variable is copied or accumulated into the 
!  HISTORY ITEM's data field for further analysis.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE History_Write( am_I_Root, Spc_Units, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE HistContainer_Mod
    USE HistItem_Mod,          ONLY : HistItem
    USE History_Netcdf_Mod
    USE History_Util_Mod
    USE MetaHistContainer_Mod, ONLY : MetaHistContainer
    USE MetaHistItem_Mod,      ONLY : MetaHistItem
    USE Registry_Params_Mod
!
! !INPUT PARAMETERS: 
!
    LOGICAL,          INTENT(IN)  :: am_I_Root   ! Are we on the root CPU?
    CHARACTER(LEN=*), INTENT(IN)  :: Spc_Units   ! Units of SC%Species array 
!
! !OUTPUT PARAMETERS: 
!
    INTEGER,          INTENT(OUT) :: RC          ! Success or failure
!
! !REMARKS:
!  This routine is called from the main program at the end of each 
!  "heartbeat" timestep.
!
! !REVISION HISTORY:
!  03 Aug 2017 - R. Yantosca - Initial version
!  21 Aug 2017 - R. Yantosca - Now get yyyymmdd, hhmmss from the container
!  28 Aug 2017 - R. Yantosca - Now save the species units to the container
!  06 Sep 2017 - R. Yantosca - Now recompute the file write and file close
!                               intervals, if they are 1 month or longer
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    LOGICAL                          :: DoClose
    LOGICAL                          :: DoWrite

    ! Strings
    CHARACTER(LEN=20)                :: TmpUnits
    CHARACTER(LEN=255)               :: ErrMsg
    CHARACTER(LEN=255)               :: ThisLoc

    ! Objects
    TYPE(MetaHistContainer), POINTER :: Collection
    TYPE(HistContainer),     POINTER :: Container
    TYPE(MetaHistItem),      POINTER :: Current

    !=======================================================================
    ! Initialize
    !=======================================================================
    RC         =  GC_SUCCESS
    DoClose    = .FALSE.
    DoWrite    = .FALSE.
    Collection => NULL()
    Current    => NULL()
    ErrMsg     =  ''
    ThisLoc    =  ' -> at History_Write (in History/history_mod.F90)' 

    !=======================================================================
    ! Loop through each DIAGNOSTIC COLLECTION in the master list, and
    ! then loop through the HISTORY ITEMS belonnging to each COLLECTION.
    !=======================================================================

    ! Point to the first COLLECTION in the master collection list
    Collection => CollectionList
    
    ! As long as this current COLLECTION is valid ...
    DO WHILE( ASSOCIATED( Collection ) ) 

       ! Point to the HISTORY CONTAINER object in this COLLECTION
       Container => Collection%Container

       !====================================================================
       ! Test if it is time to close/repopen the file or to write data
       !====================================================================

       ! Test if the "FileCloseAlarm" is ringing
       DoClose = ( ( Container%FileCloseAlarm - Container%ElapsedSec ) < EPS )

       ! Test if the "FileWriteAlarm" is ringing
       DoWrite = ( ( Container%FileWriteAlarm - Container%ElapsedSec ) < EPS )

       !====================================================================
       ! %%% GEOS-Chem "Classic" %%%
       !
       ! It is time to create a new netCDF file (closing the old one)
       !====================================================================
       IF ( DoClose ) THEN

          ! Save the units of State_Chm%Species in the container, 
          ! so that we can redefine the unit string from "TBD".
          ! Copy into a temp variable so that Gfortran won't choke.
          TmpUnits            = Spc_Units
          Container%Spc_Units = TmpUnits
         
          !-----------------------------------------------------------------
          ! If the netCDF file specified by this collection is open, 
          ! then close it and undefine all relevant object fields.
          !-----------------------------------------------------------------
          CALL History_Netcdf_Close( am_I_Root = am_I_Root,                  &
                                     Container = Container,                  &
                                     RC        = RC                         )

          ! Trap error
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = 'Error returned from "History_Netcdf_Close"!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF

          !-----------------------------------------------------------------
          ! Create the netCDF file for this HISTORY CONTAINER object,
          ! Defines each variable, saves global attributes, and writes 
          ! the index variable data to the file.
          !-----------------------------------------------------------------
          CALL History_Netcdf_Define( am_I_Root  = am_I_Root,                &
                                      Container  = Container,                &
                                      RC         = RC                       )

          ! Trap error
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = 'Error returend from "History_Netcdf_Define"!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF

          !-----------------------------------------------------------------
          ! Update "FileClose" alarm for next interval
          !-----------------------------------------------------------------

          ! Recompute the file close alarm interval if it 1 month or longer,
          ! as we will have to take into account leap years, etc.
          IF ( Container%FileCloseYmd >= 000100 ) THEN
             CALL HistContainer_FileCloseIvalSet( am_I_Root, Container, RC )
          ENDIF

          ! Update the alarm
          Container%FileCloseAlarm = Container%FileCloseAlarm                &
                                   + Container%FileCloseIvalSec

       ENDIF

       !=================================================================
       ! %%% GEOS-Chem "Classic" %%%
       !
       ! It is time to write data to the netCDF file
       !=================================================================
       IF ( DoWrite ) THEN

          !-----------------------------------------------------------------
          ! Write the HISTORY ITEMS for this collection to the netCDF file.  
          !-----------------------------------------------------------------
          CALL History_Netcdf_Write( am_I_Root = am_I_Root,                  &
                                     Container = Container,                  &
                                     RC        = RC                         )

          ! Trap error
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = 'Error returned from "History_Netcdf_Write"!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF

          !-----------------------------------------------------------------
          ! Update "FileWrite" alarm for next interval
          !-----------------------------------------------------------------

          ! Recompute the file write alarm interval if it 1 month or longer,
          ! as we will have to take into account leap years, etc.
          IF ( Container%FileWriteYmd >= 000100 ) THEN
             CALL HistContainer_FileWriteIvalSet( am_I_Root, Container, RC )
          ENDIF

          ! Update the alarm
          Container%FileWriteAlarm = Container%FileWriteAlarm                &
                                   + Container%FileWriteIvalSec

       ENDIF

       ! Skip to the next collection
       Container  => NULL()
       Collection => Collection%Next

    ENDDO

    ! Free pointers
    Container  => NULL()
    Collection => NULL()

  END SUBROUTINE History_Write
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
! !DESCRIPTION: Parses a line of the HISTORY.rc file and returns metadata
!  for a given attribute (e.g. "frequency", "template", etc.)
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GetCollectionMetaData( Line, Pattern, MetaData, nCollection )
!
! !USES:
!
    USE Charpak_Mod,       ONLY: StrSplit
    USE History_Util_Mod
!
! !INPUT PARAMETERS: 
!
    CHARACTER(LEN=*),   INTENT(IN)  :: Line          ! Line to be searched
    CHARACTER(LEN=*),   INTENT(IN)  :: Pattern       ! Search pattern
!
! !OUTPUT PARAMETERS:
!
    CHARACTER(LEN=255), INTENT(OUT) :: MetaData      ! Metadata value
    INTEGER,            INTENT(OUT) :: nCollection   ! Collection Id
!
!
! !REVISION HISTORY:
!  06 Jan 2015 - R. Yantosca - Initial version
!  03 Aug 2017 - R. Yantosca - Make search algorithm more robust
!  14 Aug 2017 - R. Yantosca - Initialize MetaData and nCollection
!  15 Aug 2017 - R. Yantosca - Bug fix: TRIM string arguments to INDEX, and
!                              initialize output arguments to undefined values
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER                  :: C, Ind, nSubStr, N, P

    ! Strings
    CHARACTER(LEN=255)       :: Name
    CHARACTER(LEN=255)       :: SubStr(255)

    ! SAVEd variables
    INTEGER,            SAVE :: CollectionStart = 0
    CHARACTER(LEN=255), SAVE :: LastName        = ''
     
    !=======================================================================
    ! Initialize
    !=======================================================================
    nCollection = UNDEFINED_INT
    MetaData    = UNDEFINED_STR

    !=======================================================================
    ! Find the metadata for the given collection
    !=======================================================================

    ! The collection name is between column 1 and the first "." character
    Ind  = INDEX( TRIM( Line ), '.' )
    Name = Line(1:Ind-1) 

    ! Non-white-space lengths of the collection cname and search pattern
    N = LEN_TRIM( Name    )
    P = LEN_TRIM( Pattern )

    ! Increment the number of the collection that we'll search from
    ! if the current collection name is different from the prior one.
    IF ( TRIM( Name ) /= TRIM( LastName ) ) THEN
       LastName        = Name
       CollectionStart = CollectionStart + 1
    ENDIF

    ! Loop over all collection names
    DO C = CollectionStart, CollectionCount

       ! Check to see if the current line matches the collection name
       ! Then check to see which collection this is in
       IF ( Name(1:30) == CollectionName(C)(1:30) ) THEN
          Ind = 1
       ELSE
          Ind = 0
       ENDIF

       ! If the we match the current collection, then ...
       IF ( Ind > 0 ) THEN

          ! Split the line on the the colon
          CALL StrSplit( Line, ':', SubStr, nSubStr )

          ! If there are 2 substrings ...
          IF ( nSubStr == 2 ) THEN

             ! Make sure the first substring matches the name 
             ! of the metadata field we would like to obtain.
             ! if it does, then we have found a match, and so return
             IF ( SubStr(1)(N+2:P+N+1) == Pattern(1:P) ) THEN
                nCollection = C
                MetaData    = CleanText( SubStr(2) )
                EXIT
             ENDIF
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
! !IROUTINE: History_Close_AllFiles
!
! !DESCRIPTION: Closes the netCDF file described by each HISTORY CONTAINER
!  object in the master list of diagnostic collections.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE History_Close_AllFiles( am_I_Root, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE History_Netcdf_Mod,    ONLY : History_Netcdf_Close
    USE MetaHistContainer_Mod, ONLY : MetaHistContainer
!
! !INPUT PARAMETERS: 
!
    LOGICAL, INTENT(IN)  :: am_I_Root   ! Are we on the root CPU?
!
! !OUTPUT PARAMETERS: 
!
    INTEGER, INTENT(OUT) :: RC          ! Success or failure
!
! !REMARKS:
!  This is called from History_Cleanup, but may also be called in other
!  locations (e.g. when processing abnormal exits)
!
! !REVISION HISTORY:
!  16 Aug 2017 - R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Strings
    CHARACTER(LEN=255)               :: ErrMsg, ThisLoc

    ! Objects
    TYPE(MetaHistContainer), POINTER :: Current

    !=======================================================================
    ! Initialize
    !=======================================================================
    RC      =  GC_SUCCESS
    Current => NULL()
    ErrMsg  =  ''
    ThisLoc =  &
       ' -> at History_Close_AllFiles (in module History/history_mod.F90)'

    !=======================================================================
    ! Close the netCDF file for each diagnostic collection
    !=======================================================================

    ! Set CURRENT to the first entry in the list of HISTORY CONTAINERS
    Current => CollectionList

    ! If this entry is not null ...
    DO WHILE ( ASSOCIATED( Current ) )

       ! Close the file (if it's open) and reset all relevant fields
       ! in the HISTORY CONTAINER object
       CALL History_Netcdf_Close( am_I_Root = am_I_Root,                    &
                                  Container = Current%Container,            &
                                  RC        = RC                           )

       ! Trap error
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Error returned from "History_Netcdf_Close"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          Current => NULL()
          RETURN
       ENDIF
       
       ! Go to the next entry in the list of HISTORY CONTAINERS
       Current => Current%Next
    ENDDO

    !=======================================================================
    ! Cleanup and quit
    !=======================================================================

    ! Free pointer
    Current => NULL()

  END SUBROUTINE History_Close_AllFiles
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: History_Cleanup
!
! !DESCRIPTION: Deallocates all module variables and objects.  Also closes
!  any remaining open netCDF files.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE History_Cleanup( am_I_Root, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE History_Netcdf_Mod,    ONLY : History_Netcdf_Cleanup
    USE MetaHistContainer_Mod, ONLY : MetaHistContainer_Destroy
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
!  14 Aug 2017 - R. Yantosca - Call History_Netcdf_Close to close open files
!  16 Aug 2017 - R. Yantosca - Move netCDF close code to History_Close_AllFiles
!  26 Sep 2017 - R. Yantosca - Now call MetaHistItem_Destroy to finalize the
!                              ContainerList object, instead of DEALLOCATE
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
     ! Strings
     CHARACTER(LEN=255)               :: ErrMsg, ThisLoc

     ! Objects
     TYPE(MetaHistContainer), POINTER :: Current

     !======================================================================
     ! Initialize
     !======================================================================
     RC      =  GC_SUCCESS
     Current => NULL()
     ErrMsg  =  ''
     ThisLoc =  ' -> at History_Cleanup (in module History/history_mod.F90)'

     !======================================================================
     ! Close all remanining netCDF files
     !======================================================================
     CALL History_Close_AllFiles( am_I_Root, RC )
     IF ( RC /= GC_SUCCESS ) THEN
        ErrMsg = 'Error returned from "History_Close_AllFiles"!'
        CALL GC_Error( ErrMsg, RC, ThisLoc )
        RETURN
     ENDIF

     !======================================================================
     ! Then finalize the history_netcdf_mod.F90 module
     !======================================================================
     CALL History_Netcdf_Cleanup( am_I_Root, RC )
     IF ( RC /= GC_SUCCESS ) THEN
        ErrMsg = 'Error returned from "History_Netcdf_Cleanup"!'
        CALL GC_Error( ErrMsg, RC, ThisLoc )
        RETURN
     ENDIF

     !======================================================================
     ! And deallocate variables belonging to history_mod.F90
     !======================================================================
     IF ( ASSOCIATED( CollectionList ) ) THEN
        CALL MetaHistContainer_Destroy( am_I_Root, CollectionList, RC )
        IF ( RC /= GC_SUCCESS ) THEN
           ErrMsg = 'Could not destroy "CollectionList"!'
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

     IF ( ALLOCATED( CollectionAccInterval ) ) THEN
        DEALLOCATE( CollectionAccInterval, STAT=RC )
        IF ( RC /= GC_SUCCESS ) THEN
           ErrMsg = 'Could not deallocate "CollectionAccInterval"!'
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
