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
  PUBLIC  :: History_Cleanup
!
! PRIVATE MEMBER FUNCTIONS:
!
  PRIVATE :: History_ReadCollectionNames
  PRIVATE :: History_ReadCollectionData
  PRIVATE :: History_AddItemToCollection
  PRIVATE :: History_Close_AllFiles
!
! !REMARKS:
!
!
! !REVISION HISTORY:
!  06 Jan 2015 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
  ! Scalars
  INTEGER                              :: CollectionCount

  ! Strings
  CHARACTER(LEN=255),      ALLOCATABLE :: CollectionName       (:  )
  CHARACTER(LEN=255),      ALLOCATABLE :: CollectionFileName   (:  )
  CHARACTER(LEN=255),      ALLOCATABLE :: CollectionTemplate   (:  )
  CHARACTER(LEN=255),      ALLOCATABLE :: CollectionFormat     (:  )
  CHARACTER(LEN=255),      ALLOCATABLE :: CollectionFrequency  (:  )
  CHARACTER(LEN=255),      ALLOCATABLE :: CollectionAccInterval(:  )
  CHARACTER(LEN=255),      ALLOCATABLE :: CollectionDuration   (:  )
  CHARACTER(LEN=255),      ALLOCATABLE :: CollectionMode       (:  )
  CHARACTER(LEN=255),      ALLOCATABLE :: CollectionLonRange   (:  )
  CHARACTER(LEN=255),      ALLOCATABLE :: CollectionLatRange   (:  )
  INTEGER,                 ALLOCATABLE :: CollectionSubsetInd  (:,:)
  CHARACTER(LEN=255),      ALLOCATABLE :: CollectionLevels     (:  )
  INTEGER,                 ALLOCATABLE :: CollectionLevelInd   (:,:)

  ! Objects
  TYPE(MetaHistContainer), POINTER     :: CollectionList
!
! !DEFINED PARAMETERS:
!
  ! Maximum number of collections (set to a ridiculously big number)
  INTEGER,                 PARAMETER   :: MAX_COLLECTIONS = 500

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
  SUBROUTINE History_Init( Input_Opt, State_Met, State_Chm, State_Diag, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE History_Util_Mod
    USE Input_Opt_Mod,      ONLY : OptInput
    USE State_Chm_Mod ,     ONLY : ChmState
    USE State_Diag_Mod,     ONLY : DgnState
    USE State_Met_Mod,      ONLY : MetState
!
! !INPUT PARAMETERS:
!
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
!  History_ReadCollectionData
!
! !REVISION HISTORY:
!  06 Jan 2015 - R. Yantosca - Initial version
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
     ' -> at History_Init (in module History/history_mod.F90)'

    !=======================================================================
    ! First initialize the list of collections
    ! ("collection" = a netCDF file with a specific update frequency)
    !=======================================================================
    IF ( .not. Input_Opt%DryRun ) THEN
       CALL History_ReadCollectionNames( Input_Opt,  State_Chm, &
                                         State_Diag, State_Met, RC )

       ! Trap potential errors
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Error encountered in "History_ReadCollectionNames"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDIF

    !=======================================================================
    ! Then determine the fields that will be saved to each collection
    ! NOTE: For dry-run, enter to print out file name & status
    !=======================================================================
    CALL History_ReadCollectionData( Input_Opt,  State_Chm, &
                                     State_Diag, State_Met, RC )

    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered in "History_ReadCollectionData"!'
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
  SUBROUTINE History_ReadCollectionNames( Input_Opt, State_Chm,  &
                                          State_Diag, State_Met, RC )
!
! !USES:
!
    USE Charpak_Mod
    USe DiagList_Mod,      ONLY : CollList, ColItem
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
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    LOGICAL                :: EOF
    INTEGER                :: fId,    IOS
    INTEGER                :: N,      nSubs1,  nSubs2

    ! Strings
    CHARACTER(LEN=255)     :: ErrMsg, ThisLoc, Line,  Line2

    ! String arrays
    CHARACTER(LEN=255)     :: Subs1(255)
    CHARACTER(LEN=255)     :: Subs2(255)
    CHARACTER(LEN=255)     :: TmpCollectionName(MAX_COLLECTIONS)

    ! Objects
    TYPE(ColItem), POINTER :: Current

    !=======================================================================
    ! Initialize
    !=======================================================================

    ! Zero local variables
    EOF     = .FALSE.
    IOS     = 0
    RC      = GC_SUCCESS
    ErrMsg  = ''
    ThisLoc = &
     ' -> at History_ReadCollectionNames (in module History/history_mod.F90)'

    ! Zero global variables
    CollectionCount   = 0
    TmpCollectionName = ''

    !=======================================================================
    ! Get the number of collections and list of collection names by
    ! querying the collection list object (CollList, from diaglist_mod.F90).
    !
    ! NOTE: We are importing CollList from diaglist_mod.F90 via a USE
    ! association.  This might not be the best way to share data (it
    ! violates data encapsulation).  But it works for now.  Maybe figure
    ! out a more elegant method later. (bmy, 2/28/18)
    !=======================================================================

    ! Initialize
    CollectionCount = 0

    ! Point to head of collection list
    Current => CollList%Head

    ! While we are not at the end of the collection list
    DO WHILE ( ASSOCIATED( Current ) )

       ! Increment the collection count
       CollectionCount = CollectionCount + 1

       ! Save the collection name in a temporary arrayu
       TmpCollectionName(CollectionCount) = TRIM( Current%CName )

       ! Point to next collection
       Current => Current%Next

    ENDDO

    ! Free pointer
    Current => NULL()

    !=======================================================================
    ! Now that we now the number of diagnostic collections, we can
    ! allocate the arrays that will hold various collection attributes
    !=======================================================================

    ! Allocate CollectionName
    IF ( .not. ALLOCATED( CollectionName ) ) THEN
       ALLOCATE( CollectionName( CollectionCount ), STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Could not allocate "CollectionName"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDIF

    ! Copy the collection names from the temporary array
    ! NOTE: The linked-list stors in reverse order, so revert here
    DO N = 1, CollectionCount
       CollectionName(N) = TmpCollectionName(CollectionCount-N+1)
    ENDDO

    ! Allocate CollectionFileName
    IF ( .not. ALLOCATED( CollectionFileName ) ) THEN
       ALLOCATE( CollectionFileName( CollectionCount ), STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Could not allocate "CollectionFileName"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
       CollectionFileName = UNDEFINED_STR
    ENDIF

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

    ! Allocate CollectionLonRange
    IF ( .not. ALLOCATED( CollectionLonRange ) ) THEN
       ALLOCATE( CollectionLonRange( CollectionCount ), STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Could not allocate "CollectionLonRange"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
       CollectionLonRange = UNDEFINED_STR
    ENDIF

    ! Allocate CollectionLatRange
    IF ( .not. ALLOCATED( CollectionLatRange ) ) THEN
       ALLOCATE( CollectionLatRange( CollectionCount ), STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Could not allocate "CollectionLatRange"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
       CollectionLatRange = UNDEFINED_STR
    ENDIF

    ! Allocate CollectionSubsetInd
    IF ( .not. ALLOCATED( CollectionSubsetInd ) ) THEN
       ALLOCATE( CollectionSubsetInd( 4, CollectionCount ), STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Could not allocate "CollectionSubsetInd"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
       CollectionSubsetInd = UNDEFINED_INT
    ENDIF

    ! Allocate CollectionLevels
    IF ( .not. ALLOCATED( CollectionLevels ) ) THEN
       ALLOCATE( CollectionLevels( CollectionCount ), STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Could not allocate "CollectionLevels"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
       CollectionLevels = UNDEFINED_STR
    ENDIF

    ! Allocate CollectionLevelInd
    IF ( .not. ALLOCATED( CollectionLevelInd ) ) THEN
       ALLOCATE( CollectionLevelInd( 2, CollectionCount ), STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Could not allocate "CollectionLevelInt"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
       CollectionLevelInd = UNDEFINED_INT
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
  SUBROUTINE History_ReadCollectionData( Input_Opt, State_Chm,  &
                                         State_Diag, State_Met, RC )
!
! !USES:
!
    USE Charpak_Mod
    USE DiagList_Mod,          ONLY : CollList, Search_CollList
    USE ErrCode_Mod
    USE Grid_Registry_Mod,     ONLY : Lookup_Grid
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
!  Private routine, called from History_Init.
!
! !REVISION HISTORY:
!  16 Jun 2017 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
     ! Scalars
    LOGICAL                      :: EOF,            Found
    LOGICAL                      :: FileExists
    INTEGER                      :: yyyymmdd,       hhmmss
    INTEGER                      :: yyyymmdd_end,   hhmmss_end
    INTEGER                      :: DeltaYMD,       DeltaHMS
    INTEGER                      :: X,              Y
    INTEGER                      :: C,              N,             W
    INTEGER                      :: nX,             nY,            nZ
    INTEGER                      :: fId,            IOS,           LineNum
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
    INTEGER                      :: HbHrs,          HbMin,         HbSec
    INTEGER                      :: HeartBeatHms,   nTags
    REAL(f8)                     :: UpdateAlarm,    HeartBeatDtSec
    REAL(f8)                     :: FileWriteAlarm, FileCloseAlarm
    REAL(f8)                     :: JulianDate,     JulianDateEnd
    REAL(f8)                     :: UpdateCheck,    FileWriteCheck
    REAL(f8)                     :: SimLengthSec

    ! Strings
    CHARACTER(LEN=6  )           :: TStr
    CHARACTER(LEN=8  )           :: DStr
    CHARACTER(LEN=20 )           :: StartTimeStamp, EndTimeStamp
    CHARACTER(LEN=63 )           :: CName
    CHARACTER(LEN=80 )           :: ErrorLine
    CHARACTER(LEN=255)           :: FileExpId
    CHARACTER(LEN=255)           :: Line,           FileName
    CHARACTER(LEN=255)           :: OutputName,     ThisLoc
    CHARACTER(LEN=255)           :: MetaData,       Reference
    CHARACTER(LEN=255)           :: Title,          Units
    CHARACTER(LEN=255)           :: ItemTemplate,   ItemTemplateUC
    CHARACTER(LEN=255)           :: ItemName,       Description
    CHARACTER(LEN=255)           :: TmpMode,        Contact
    CHARACTER(LEN=255)           :: Pattern,        ItemPrefix
    CHARACTER(LEN=255)           :: tagId,          tagName
    CHARACTER(LEN=512)           :: ErrMsg,         FileMsg

    ! Arrays
    REAL(f8)                     :: Subset(2)
    INTEGER                      :: Levels(200)
    CHARACTER(LEN=255)           :: Subs1(255)
    CHARACTER(LEN=255)           :: Subs2(255)
    CHARACTER(LEN=255)           :: SubStrs(255)

    ! Objects
    TYPE(HistContainer), POINTER :: Container
    TYPE(HistItem),      POINTER :: Item
    TYPE(Species),       POINTER :: ThisSpc

    ! Pointer arrays
    REAL(f8),            POINTER :: Grid_Lat (:    )
    REAL(f8),            POINTER :: Grid_LatE(:    )
    REAL(f8),            POINTER :: Grid_Lon (:    )
    REAL(f8),            POINTER :: Grid_LonE(:    )
    REAL(fp),            POINTER :: Ptr3d    (:,:,:)
    REAL(f4),            POINTER :: Ptr3d_4  (:,:,:)

    !=======================================================================
    ! Initialize
    !=======================================================================

    ! Assume success
    RC             =  GC_SUCCESS

    ! Skip for GEOS-Chem dry-run simulations
    IF ( .not. Input_Opt%DryRun ) THEN

       ! Initialize variables
       EOF            =  .FALSE.
       IOS            =  0
       UpdateYmd      =  0
       UpdateHms      =  0
       FileCloseYmd   =  0
       FileCloseHms   =  0
       FileWriteYmd   =  0
       FileWriteHms   =  0
       LineNum        =  0
       SpaceDim       =  0
       HeartBeatDtSec =  DBLE( Input_Opt%TS_DYN )
       yyyymmdd       =  Input_Opt%NymdB
       hhmmss         =  Input_Opt%NhmsB
       yyyymmdd_end   =  Input_Opt%NymdE
       hhmmss_end     =  Input_Opt%NhmsE
       Subset         =  UNDEFINED_DBL
       Levels         =  UNDEFINED_INT

       ! Compute the YMD and HMS intervals for collections specified with "End",
       ! such as for restart files.  NOTE: This algorithm should work with most
       ! common model simulation intervals, but there might be some edge cases
       ! that will cause it to fail.  It is still an improvement. (bmy, 2/26/19)
       CALL Compute_DeltaYmdHms_For_End( yyyymmdd,     hhmmss,               &
                                         yyyymmdd_end, hhmmss_end,           &
                                         deltaYMD,     deltaHMS             )

       ! Convert the HeartBeatDtSec into hours:minutes:seconds
       ! for defining the Update interval for time-averaged collections
       HbMin          = HeartBeatDtSec / 60
       HbHrs          = HbMin / 60
       HbSec          = HeartBeatDtSec - ( HbMin * 60 ) - ( HbHrs * 3600 )
       HeartBeatHms   = ( HbHrs * 10000 ) + ( HbMin * 100 ) + HbSec

       ! Initialize objects and pointers
       Container      => NULL()
       Item           => NULL()
       Ptr3d          => NULL()
       Ptr3d_4        => NULL()
       ThisSpc        => NULL()
       Grid_Lat       => NULL()
       Grid_LatE      => NULL()
       Grid_Lon       => NULL()
       Grid_LonE      => NULL()

       ! Initialize Strings
       Description    =  ''
       ErrMsg         =  ''
       Contact        =  &
         'GEOS-Chem Support Team (geos-chem-support@as.harvard.edu)'
       Reference      =  'www.geos-chem.org; wiki.geos-chem.org'
       ThisLoc        =  &
         ' -> at History_ReadCollectionData (in module History/history_mod.F90)'
       Units          =  ''
       FileExpId      =  ''

       ! Create the timestamp at the start of the simulation
       WRITE( DStr,          '(i8.8)' ) yyyymmdd
       WRITE( TStr,          '(i6.6)' ) hhmmss
       WRITE( StartTimeStamp, 300     ) DStr(1:4), DStr(5:6), DStr(7:8),       &
                                        TStr(1:2), TStr(3:4), TStr(5:6)

       ! Create the timestamp at the end of the simulation
       WRITE( DStr,        '(i8.8)' ) yyyymmdd_end
       WRITE( TStr,        '(i6.6)' ) hhmmss_end
       WRITE( EndTimeStamp, 300     ) DStr(1:4), DStr(5:6), DStr(7:8),         &
                                      TStr(1:2), TStr(3:4), TStr(5:6)

       ! Format string
 300   FORMAT( a4, '-', a2, '-', a2, ' ', a2, ':', a2, ':', a2, 'z' )

       ! Compute the Astronomical Julian Date corresponding to the yyyymmdd
       ! and hhmmss values at the start and end of the simulation, which are
       ! needed below.  This can be done outside of the DO loop below.
       CALL Compute_Julian_Date( yyyymmdd,     hhmmss,     JulianDate    )
       CALL Compute_Julian_Date( yyyymmdd_end, hhmmss_end, JulianDateEnd )

       ! Compute the length of the simulation, in elapsed seconds
       SimLengthSec   = NINT( ( JulianDateEnd - JulianDate ) * SECONDS_PER_DAY )

       !====================================================================
       ! Get pointers to the grid longitudes and latitudes
       !====================================================================

       ! Lookup latitude centers
       CALL Lookup_Grid( Input_Opt = Input_Opt,  &
                         Variable  = 'GRID_LAT', &
                         Ptr1d_8   = Grid_Lat,   &
                         RC        = RC )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Could not get pointer to latitudes (aka GRID_LAT)!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       ! Lookup latitude edges
       CALL Lookup_Grid( Input_Opt = Input_Opt,   &
                         Variable  = 'GRID_LATE', &
                         Ptr1d_8   = Grid_LatE,   &
                         RC        = RC )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Could not get pointer to latitude edges (aka GRID_LATE)!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       ! Lookup longitude centers
       CALL Lookup_Grid( Input_Opt = Input_Opt,  &
                         Variable  = 'GRID_LON', &
                         Ptr1d_8   = Grid_Lon,   &
                         RC        = RC )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Could not get pointer to longitudes (aka GRID_LON)!'
          CALL GC_Error( ErrMsg, RC, ThisLoc)
          RETURN
       ENDIF

       ! Lookup longitude edges
       CALL Lookup_Grid( Input_Opt = Input_Opt,   &
                         Variable  = 'GRID_LONE', &
                         Ptr1d_8   = Grid_LonE,   &
                         RC        = RC )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Could not get pointer to longitude edges (aka GRID_LONE)!'
          CALL GC_Error( ErrMsg, RC, ThisLoc)
          RETURN
       ENDIF
    ENDIF

    !=======================================================================
    ! Open the file containing the list of requested diagnostics
    !=======================================================================

    ! Test if the file exists
    INQUIRE( FILE=TRIM( Input_Opt%HistoryInputFile ), EXIST=FileExists )

    ! Test if the file exists and define an output string
    IF ( FileExists ) THEN
       FileMsg = 'HISTORY (INIT): Opening'
    ELSE
       FileMsg = 'HISTORY (INIT): REQUIRED FILE NOT FOUND'
    ENDIF

    ! Write message to stdout for both regular and dry-run simulations
    IF ( Input_Opt%AmIRoot ) THEN
       WRITE( 6, 350 ) TRIM( FileMsg ), TRIM( Input_Opt%HistoryInputFile )
 350   FORMAT( a, ' ', a )
    ENDIF

    ! For dry-run simulations, return to calling program.
    ! For regular simulations, throw an error if we can't find the file.
    IF ( Input_Opt%DryRun ) THEN
       RETURN
    ELSE
       IF ( .not. FileExists ) THEN
          WRITE( ErrMsg, 350 ) TRIM( FileMsg                    ),           &
                               TRIM( Input_Opt%HistoryInputFile )
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDIF

    ! Find a free file unit
    fId     = FindFreeLun()

    ! Open the file
    OPEN( fId, FILE=TRIM(Input_Opt%HistoryInputFile), STATUS='OLD', IOSTAT=RC )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error opening "' //TRIM(Input_Opt%HistoryInputFile) // '"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !=======================================================================
    ! Read data from the file
    !=======================================================================
    DO

       ! Read a single line, and strip leading/trailing spaces
       ! and keep track of the line number for error output
 500   CONTINUE
       Line    = ReadOneLine( fId, EOF, IOS, Squeeze=.TRUE. )
       LineNum = LineNum + 1

       ! Exit the loop if it's the end of the file
       IF ( EOF ) EXIT

       ! If it's a real I/O error, quit w/ error message
       IF ( IOS > 0 ) THEN
          ErrMsg = 'Unexpected end-of-file in "'                          // &
                    TRIM( Input_Opt%HistoryInputFile )
          WRITE( ErrorLine, 250 ) LineNum
 250      FORMAT( ' -> ERROR occurred at (or near) line ', i6,               &
                      ' of the HISTORY.rc file' )
          CALL GC_Error( ErrMsg, RC, ThisLoc, ErrorLine )
          RETURN
       ENDIF

       ! Skip if the line is commented out
       IF ( Line(1:1) == "#" ) CYCLE

       ! Zero variables
       FileCloseYmd   = 0
       FileCloseHms   = 0
       FileWriteYmd   = 0
       FileWriteHms   = 0
       FileWriteCheck = 0.0_f8
       UpdateYmd      = 0
       UpdateHms      = 0
       UpdateCheck    = 0.0_f8

       !====================================================================
       ! Get the EXPID string.  This is the "front part" of the netCDF
       ! file path for each collection.  In other words, if EXPID is
       ! "OutputDir/GEOSChem", then the default SpeciesConc collection file
       ! names will be "OutputDir/GEOSChem.SpeciesConc_YYYYMMDD_hhmmz.nc4"
       !====================================================================
       IF ( INDEX( Line, 'EXPID' ) > 0  ) THEN

          ! Split the line on the colon
          CALL StrSplit( Line, ":", Subs1, nSubs1 )

          ! Stop with error if there are more than 2 substrings
          IF ( nSubs1 /= 2 ) THEN
             ErrMsg = 'Error in extracting the EXPID value from the '     // &
                      'HISTORY.rc file.  This forms the start of the '    // &
                      'netCDF file name for each collection.  Please '    // &
                      'check the HISTORY.rc file for typos.'
             WRITE( ErrorLine, 250 ) LineNum
             CALL GC_Error( ErrMsg, RC, ThisLoc, ErrorLine )
          ENDIF

          ! Save the EXPID parameter
          FileExpId = Subs1(2)
          CALL CStrip( FileExpId )
       ENDIF

       !====================================================================
       ! The HISTORY.rc file specifies collection metadata as:
       !
       !   SpeciesConc.filename:  './output/GEOSChem.inst.%y4%m2%d2.nc4'
       !   SpeciesConc.template:  '%y4%m2%d2_%h2%n2.nc4',
       !   SpeciesConc.format:    'CFIO',
       !   SpeciesConc.frequency:  010000,
       !   SpeciesConc.duration:   240000
       !   etc.
       !
       ! where in this example, "instantaneous" is the collection name
       ! and "filename', "template", "format", "frequency", "duration"
       ! are the metadata fields.
       !
       ! Get the metadata belonging to each collection and store them
       ! in the proper arrays for later use.  NOTE: this method does not
       ! assume that the collections are in the same order as they
       ! are listed under the COLLECTIONS section.
       !====================================================================

       ! "filename": Specifies the full filename path
       ! Can be omitted if "template" is specified
       Pattern = 'filename'
       IF ( INDEX( TRIM( Line ), TRIM( Pattern ) ) > 0 ) THEN
          CALL GetCollectionMetaData( Input_Opt, Line, Pattern, MetaData, C )
          IF ( C > 0 ) CollectionFileName(C) = Metadata
       ENDIF

       ! "template": Specifies the year/month/day/hr/min/sec in filenames
       ! Can be omitted if "filename" is specified
       Pattern = 'template'
       IF ( INDEX( TRIM( Line ), TRIM( Pattern ) ) > 0 ) THEN
          CALL GetCollectionMetaData( Input_Opt, Line, Pattern, MetaData, C )
          IF ( C > 0 ) CollectionTemplate(C) = Metadata
       ENDIF

       ! "format": Specifies the file output format (e.g. netCDF-4, CFIO)
       Pattern = 'format'
       IF ( INDEX( TRIM( Line ), TRIM( Pattern ) ) > 0 ) THEN
          CALL GetCollectionMetaData( Input_Opt, Line, Pattern, MetaData, C )
          IF ( C > 0 ) CollectionFormat(C) = Metadata
       ENDIF

       ! "frequency": Specifies how often diagnostics are updated,
       ! Must be either in "YYYYMMDD hhmmss" or "hhmmss" format.
       Pattern = 'frequency'
       IF ( INDEX( TRIM( Line ), TRIM( Pattern ) ) > 0 ) THEN
          CALL GetCollectionMetaData( Input_Opt, Line, Pattern, MetaData, C )
          IF ( C > 0 ) THEN
             IF ( LEN_TRIM( MetaData ) == 6     .or.                         &
                  LEN_TRIM( MetaData ) == 14    .or.                         &
                  TRIM(     MetaData ) == 'End' .or.                         &
                  TRIM(     MetaData ) == 'end' ) THEN
                CollectionFrequency(C) = Metadata
             ELSE
                ErrMsg = 'Error in defining "frequency" for collection "' // &
                         TRIM( CollectionName(C) ) // '"!  This field '   // &
                         'must either be of the format "YYYYMMDD '        // &
                         'hhmmss", "hhmmss", or "End".  Please check the '// &
                         '"frequency" setting in the HISTORY.rc file.'
                WRITE( ErrorLine, 250 ) LineNum
                CALL GC_Error( ErrMsg, RC, ThisLoc, ErrorLine )
                RETURN
             ENDIF
          ENDIF
       ENDIF

       ! "acc_interval": Specifies how often time-averaged diagnostics
       ! are updated.   ! Must be either in "YYYYMMDD hhmmss" or "hhmmss"
       ! format.   If omitted, "acc_interval" will be set from "frequency"
       !%%%%% NOTE: The "acc_interval" attribute is not really needed;
       !%%%%% we only really need "frequency" and "duration".  We will
       !%%%%% leave this as an "undocumented feature". (bmy, 3/26/18)
       Pattern = 'acc_interval'
       IF ( INDEX( TRIM( Line ), TRIM( Pattern ) ) > 0 ) THEN
          CALL GetCollectionMetaData( Input_Opt, Line, Pattern, MetaData, C )
          IF ( C > 0 ) THEN
             IF ( LEN_TRIM( MetaData ) == 6   .or.                           &
                  LEN_TRIM( MetaData ) == 14 ) THEN
                CollectionAccInterval(C) = Metadata
             ELSE
                ErrMsg = 'Error in defining "acc_interval" for '          // &
                         'collection "'// TRIM( CollectionName(C) )       // &
                         '"!  This field must either be of the format '   // &
                         '"YYYYMMDD hhmmss" or "hhmmss".  Please check '  // &
                         'the "acc_interval" setting in the HISTORY.rc '  // &
                         'file.'
                WRITE( ErrorLine, 250 ) LineNum
                CALL GC_Error( ErrMsg, RC, ThisLoc, ErrorLine )
                RETURN
              ENDIF
          ENDIF
       ENDIF

       ! "duration:: Specifies how often files will be written.
       ! Must be either in "YYYYMMDD hhmmss" or "hhmmss" format.
       ! If omitted, "duration" will be set from "frequency"
       Pattern = 'duration'
       IF ( INDEX( TRIM( Line ), TRIM( Pattern ) ) > 0 ) THEN
          CALL GetCollectionMetaData( Input_Opt, Line, Pattern, MetaData, C )
          IF ( C > 0 ) THEN
             IF ( LEN_TRIM( MetaData ) == 6     .or.                         &
                  LEN_TRIM( MetaData ) == 14    .or.                         &
                  TRIM(     MetaData ) == 'End' .or.                         &
                  TRIM(     MetaData ) == 'end' ) THEN
                CollectionDuration(C) = Metadata
             ELSE
                ErrMsg = 'Error in defining "duration" for collection "'  // &
                         TRIM( CollectionName(C) ) // '"!  This field '   // &
                         'must either be of the format "YYYYMMDD '        // &
                         'hhmmss", "hhmmss", or "End".  Please check the '// &
                         '"duration" setting in the HISTORY.rc file.'
                WRITE( ErrorLine, 250 ) LineNum
                CALL GC_Error( ErrMsg, RC, ThisLoc, ErrorLine )
                RETURN
             ENDIF
          ENDIF
       ENDIF

       ! "mode": Specifies instantaneous or time-averaged archiving.
       ! Throw an error if anything else is specified
       Pattern = 'mode'
       IF ( INDEX( TRIM( Line ), TRIM( Pattern ) ) > 0 ) THEN
          CALL GetCollectionMetaData( Input_Opt, Line, Pattern, MetaData, C )
          IF ( C > 0 ) THEN
             TmpMode = Metadata
             CALL TranUc( TmpMode )
             SELECT CASE( TmpMode )
                CASE( 'INSTANTANEOUS', 'TIME-AVERAGED', 'TIMEAVERAGED' )
                   CollectionMode(C) = Metadata
                CASE DEFAULT
                   ErrMsg = 'Error in defining "mode" for collection "'   // &
                             TRIM( CollectionName(C) ) // '"!  The mode ' // &
                            'value can either be "instantaneous" or '     // &
                            '"time-averaged".  Please check the "mode" '  // &
                            'setting in the HISTORY.rc file.'
                   WRITE( ErrorLine, 250 ) LineNum
                   CALL GC_Error( ErrMsg, RC, ThisLoc, ErrorLine )
                   RETURN
             END SELECT
          ENDIF
       ENDIF

       ! "LON_RANGE": Specifies a longitude range for subsetting
       ! the data grid. The required order is: lonMin, lonMax
       Pattern = 'LON_RANGE'
       Subset  =  UNDEFINED_DBL
       IF ( INDEX( TRIM( Line ), TRIM( Pattern ) ) > 0 ) THEN

          ! First split the line by colon
          CALL StrSplit( Line, ":", Subs1, nSubs1 )
          IF ( C > 0 ) THEN

             ! Replace any commas with spaces
             CALL StrRepl( Subs1(2), ",", " " )
             CollectionLonRange(C) = Subs1(2)

             ! Then split by spaces and convert to INTEGER
             CALL StrSplit( CollectionLonRange(C), " ", Subs2, nSubs2 )
             IF ( nSubs2 == 2 ) THEN
                DO N = 1, nSubs2
                   READ( Subs2(N), '(f13.6)' ) Subset(N)
                ENDDO
             ELSE
                ErrMsg = 'Subsets must be specified as: lonmin, lonmax!'
                WRITE( ErrorLine, 250 ) LineNum
                CALL GC_Error( ErrMsg, RC, ThisLoc, ErrorLine )
                RETURN
             ENDIF

             ! Find the longitude indices for lonMin and lonMax values
             DO X = 1, SIZE( Grid_LonE )-1
                IF ( Grid_LonE(X  ) <= Subset(1)  .and.                      &
                     Grid_LonE(X+1) >  Subset(1) ) THEN
                   CollectionSubsetInd(1,C) = X
                ENDIF
                IF ( Grid_LonE(X  ) <= Subset(2)  .and.                      &
                     Grid_LonE(X+1) >  Subset(2) ) THEN
                   CollectionSubsetInd(2,C) = X
                ENDIF
             ENDDO

             ! Error check longitudes
             DO N = 1, 2
                IF ( CollectionSubsetInd(N,C) < -180.0_f8  .or.              &
                     CollectionSubsetInd(N,C) >  180.0_f8 ) THEN
                   ErrMsg = 'Invalid longitude subset values for '   //      &
                            'collection "'// TRIM(CollectionName(C)) // '"!'
                   WRITE( ErrorLine, 250 ) LineNum
                   CALL GC_Error( ErrMsg, RC, ThisLoc, ErrorLine )
                   RETURN
                ENDIF
             ENDDO
          ENDIF
       ENDIF

       ! "LON_RANGE": Specifies a latitude range for subsetting
       ! the data grid. The required order is: latMin, latMax
       Pattern = 'LAT_RANGE'
       Subset  =  UNDEFINED_DBL
       IF ( INDEX( TRIM( Line ), TRIM( Pattern ) ) > 0 ) THEN

          ! First split the line by colon
          CALL StrSplit( Line, ":", Subs1, nSubs1 )
          IF ( C > 0 ) THEN

             ! Replace any commas with spaces
             CALL StrRepl( Subs1(2), ",", " " )
             CollectionLatRange(C) = Subs1(2)

             ! Then split by spaces and convert to INTEGER
             CALL StrSplit( CollectionLatRange(C), " ", Subs2, nSubs2 )
             IF ( nSubs2 == 2 ) THEN
                DO N = 1, nSubs2
                   READ( Subs2(N), '(f13.6)' ) Subset(N)
                ENDDO
             ELSE
                ErrMsg = 'Subsets must be specified as: latMin, latMax!'
                WRITE( ErrorLine, 250 ) LineNum
                CALL GC_Error( ErrMsg, RC, ThisLoc, ErrorLine )
                RETURN
             ENDIF

             ! Find the latitude indices for latMin and latMax values
             DO Y = 1, SIZE( Grid_LatE )-1
                IF ( Grid_LatE(Y  ) <= Subset(1)  .and.                      &
                     Grid_LatE(Y+1) >  Subset(1) ) THEN
                   CollectionSubsetInd(3,C) = Y
                ENDIF
                IF ( Grid_LatE(Y  ) <= Subset(2)  .and.                      &
                     Grid_LatE(Y+1) >  Subset(2) ) THEN
                   CollectionSubsetInd(4,C) = Y
                ENDIF
             ENDDO

             ! Error check latitudes
             DO N = 3, 4
                IF ( CollectionSubsetInd(N,C) < -90.0_f8  .or.               &
                     CollectionSubsetInd(N,C) >  90.0_f8 ) THEN
                   ErrMsg = 'Invalid latitude subset values for '     //     &
                            'collection " '// TRIM(CollectionName(C)) // '"!'
                   WRITE( ErrorLine, 250 ) LineNum
                   CALL GC_Error( ErrMsg, RC, ThisLoc, ErrorLine )
                   RETURN
                ENDIF
             ENDDO
          ENDIF
       ENDIF

       ! "levels: Specifies a vertical subset of the data grid
       Pattern  = 'levels'
       IF ( INDEX( TRIM( Line ), TRIM( Pattern ) ) > 0 ) THEN

          ! First split the line by colon
          CALL StrSplit( Line, ":", Subs1, nSubs1 )
          IF ( C > 0 ) THEN

             ! Replace any commas with spaces
             CALL StrRepl( Subs1(2), ",", " " )
             CollectionLevels(C) = Subs1(2)

             ! Then split by spaces and convert to INTEGER
             ! Also compute the min and max level
             CALL StrSplit( CollectionLevels(C), " ", Subs2, nSubs2 )
             IF ( nSubs2 <= SIZE( Levels ) ) THEN
                DO N = 1, nSubs2
                   READ( Subs2(N), '(i10)' ) Levels(N)
                   IF ( Levels(N) < 0 ) THEN
                      ErrMsg = TRIM( CollectionName(C) ) // '.levels '    // &
                               'must not have any negative values!'
                      WRITE( ErrorLine, 250 ) LineNum
                      CALL GC_Error( ErrMsg, RC, ThisLoc, ErrorLine )
                      RETURN
                   ENDIF
                ENDDO
             ELSE
                ErrMsg = 'Too many levels specified for collection "'     // &
                          TRIM( CollectionName(C) ) // '" Must be <= 200.'
                WRITE( ErrorLine, 250 ) LineNum
                CALL GC_Error( ErrMsg, RC, ThisLoc, ErrorLine )
                RETURN
             ENDIF

             ! Of all the levels that are specified, store the min and max
             ! in the CollectionLevelInd array.  We will save out all of
             ! the levels between the min and max.
             ! NOTE: GCHP HISTORY can archive out individual levels, but
             ! this is trickier to implement in GC "Classic".  It is easier
             ! to point to a contiguous array subslice, so we will just
             ! archive everything between the min and max level for the
             ! time being. (bmy, 7/18/19)
             CollectionLevelInd(1,C) = MINVAL( Levels(1:nSubs2) )
             CollectionLevelInd(2,C) = MAXVAL( Levels(1:nSubs2) )
          ENDIF
       ENDIF

       !====================================================================
       ! NOTE: We assume FIELDS is the last metadata tag for the
       ! collection.  We need to create the collection object and
       ! an object for each history item stored in the collection.
       !====================================================================
       IF ( INDEX( TRIM( Line ), 'fields' ) > 0 ) THEN

          !-----------------------------------------------------------------
          ! If we can't find the metadata for the collection in HISTORY.rc,
          ! then this might point to a mismatch between names under the
          ! "COLLECTIONS:" list and the corresponding metadata section.
          ! Do some further error checking.
          !-----------------------------------------------------------------
          IF ( C == UNDEFINED_INT ) THEN

             !--------------------------------------------------------------
             ! If the collection corresponding to this ".fields" tag is
             ! not active, then keep reading lines from HISTORY.rc
             ! until we reach the next collection definition section.
             ! then cycle back to the top of the loop.
             !--------------------------------------------------------------

             ! Get the collection name (its to the left of the first ".")
             N     = INDEX( TRIM( Line ), '.' )
             CName = Line(1:N-1)

             ! This means skipping over all of the fields listed under
             ! this collection until we get to the :: separator
             CALL Search_CollList( Input_Opt%amIRoot, CollList, CName, Found, RC )
             IF ( .not. Found ) THEN
                DO
                   Line    = ReadOneLine( fId, EOF, IOS, Squeeze=.TRUE. )
                   LineNum = LineNum + 1
                   CALL StrSqueeze( Line )
                   IF ( TRIM( Line ) == '::' ) GOTO 500
                ENDDO
             ENDIF

             !--------------------------------------------------------------
             ! If we get to this point, then there is a true error
             ! condition.  Print an error message asking the user to
             ! check the HISTORY.rc file for inconsistencies.
             !--------------------------------------------------------------

             ! List the defined collections
             WRITE( 6, '(/,a)' ) REPEAT( '=', 79 )
             WRITE( 6, 200   )
 200         FORMAT( 'GEOS-Chem ERROR: One or more collection ',             &
                     'attributes do not correspond', /                       &
                     'to any of these defined collection names '             &
                     'in the "HISTORY.rc" input file:', /                   )

             DO N = 1, CollectionCount
                WRITE( 6, 210 ) N, TRIM( CollectionName(N) )
 210            FORMAT( i3, ') ', a )
             ENDDO

             WRITE( 6, 220 )
 220         FORMAT( /, 'Please check the HISTORY.rc file for any ',         &
                     'missing ":" or "," characters', /,                     &
                     'in the collection attributes.'                        )
             WRITE( 6, '(a,/)' ) REPEAT( '=', 79 )

             ! Write error message and then return
             ErrMsg = 'Inconsistency in collection names and attributes!' // &
                      ' Please check "HISTORY.rc" for typos.'
             WRITE( ErrorLine, 250 ) LineNum
             CALL GC_Error( ErrMsg, RC, ThisLoc, ErrorLine )
             RETURN
          ENDIF

          !=================================================================
          ! At this point we are sure that the collection has been
          ! defined, so we can continue to populate it with fields.
          !=================================================================

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
                   ErrMsg = 'Acc_interval cannot be defined for '         // &
                            'instantaneous collection: "'                 // &
                            TRIM( CollectionName(C) )                     // &
                            '"!'
                   WRITE( ErrorLine, 250 ) LineNum
                   CALL GC_Error( ErrMsg, RC, ThisLoc, ErrorLine )
                   RETURN
                ENDIF

          END SELECT

          !-----------------------------------------------------------------
          ! %%%%% INSTANTANEOUS AND TIME-AVERAGED COLLECTIONS %%%%%
          !
          ! Define the "File Write" interval
          !
          ! The ".frequency" tag in HISTORY.rc specifies the interval at
          ! which data will be written to the netCDF file.  Thus, we can
          ! set FileWriteYmd and FileWriteHms from CollectionFrequency.
          !
          ! NOTE: If CollectionFrequency is 6 digits long, then assume
          ! that to be FileWriteHms.  If longer, then assume that it is
          ! both FileWriteYmd and FileWriteHms.  This is a hack that we
          ! introduced for GEOS-Chem "Classic" only, as this feature is
          ! not yet supported in MAPL.  (sde, bmy, 8/4/17, 10/26/17)
          !
          ! Add capability to set frequency to 'End'. In that case, the
          ! data will be written to the netCDF file at the end of the
          ! simulation. This is especially useful for the Restart collection
          ! for saving fields needed for subsequent GEOS-Chem runs.
          ! FileCloseYmd and FileCloseHms will be computed as the amount
          ! of time between the start and end of the simulation.
          ! (mps, 10/12/18)
          !-----------------------------------------------------------------
          IF ( LEN_TRIM( CollectionFrequency(C) ) == 6 ) THEN
             READ( CollectionFrequency(C), '(i6.6)'  ) FileWriteHms
          ELSE IF ( LEN_TRIM( CollectionFrequency(C) ) == 14 ) THEN
             READ( CollectionFrequency(C), '(i8,i6)' ) FileWriteYmd,         &
                                                       FileWriteHms
          ELSE IF ( TRIM( CollectionFrequency(C) ) == 'End'   .or.         &
                    TRIM( CollectionFrequency(C) ) == 'end' ) THEN
             FileWriteYmd = DeltaYMD
             FileWriteHms = DeltaHMS
          ENDIF

          ! SPECIAL CASE: If FileWriteHms is 240000, set
          ! FileWriteYmd=000001 and FileWriteHms=000000
          IF ( FileWriteHms == 240000 ) THEN
             FileWriteYmd = 00000001
             FileWriteHms = 000000
          ENDIF

          !-----------------------------------------------------------------
          ! %%%%% INSTANTANEOUS AND TIME-AVERAGED COLLECTIONS %%%%%
          !
          ! Define the "File Close" interval
          !
          ! The ".duration" tag in HISTORY.rc denotes the interval when
          ! each new netCDF file will be produced.  Thus, we can set
          ! FileCloseYmd and FileCloseHms from CollectionDuration.
          !
          ! If ".duration" is not specified in HISTORY.rc, then both
          ! FileCloseYmd and FileCloseHms will both be defined from
          ! the ".frequency" tag (stored in CollectionFrequency).
          !
          ! NOTE: If CollectionDuration is 6 digits long, then assume
          ! that to be FileCloseHms.  If longer, then assume that it is
          ! both FileCloseYmd and FileCloseHms.  This is a hack that we
          ! introduced for GEOS-Chem "Classic" only, as this feature is
          ! not yet supported in MAPL.  (sde, bmy, 8/4/17, 10/26/17)
          !
          ! Add capability to set duration to 'End'. In that case, the
          ! netCDF file will be closed at the end of the simulation.
          ! This is especially useful for the Restart collection
          ! for saving fields needed for subsequent GEOS-Chem runs.
          ! FileCloseYmd and FileCloseHms will be computed as the amount
          ! of time between the start and end of the simulation.
          ! (mps, 10/12/18)
          !-----------------------------------------------------------------
          IF ( TRIM( CollectionDuration(C) ) == UNDEFINED_STR ) THEN
             FileCloseYmd = FileWriteYmd
             FileCloseHms = FileWriteHms
          ELSE IF ( LEN_TRIM( CollectionDuration(C) ) == 6 ) THEN
             READ( CollectionDuration(C), '(i6.6)'  ) FileCloseHms
          ELSE IF ( LEN_TRIM( CollectionDuration(C) ) == 14 ) THEN
             READ( CollectionDuration(C), '(i8,i6)' ) FileCloseYmd,          &
                                                      FileCloseHms
          ELSE IF ( TRIM( CollectionDuration(C) ) == 'End'   .or.          &
                    TRIM( CollectionDuration(C) ) == 'end' ) THEN
             FileCloseYmd = DeltaYMD
             FileCloseHms = DeltaHMS
          ENDIF

          ! SPECIAL CASE: If FileCloseHms is 240000, set
          ! FileCloseYmd=000001 and FileCloseYmd=000000
          IF ( FileCloseHms == 240000 ) THEN
             FileCloseYmd = 00000001
             FileCloseHms = 000000
          ENDIF

          IF ( Operation == COPY_FROM_SOURCE ) THEN

             !--------------------------------------------------------------
             ! %%%%% INSTANTANEOUS COLLECTION %%%%%
             !
             ! Define the "Update" interval
             !
             ! Because there is no time-averaging, each field is written to
             ! the netCDF file as soon as it is updated.  Thus, we can set
             ! UpdateYmd and UpdateHms from the ".frequency" tag in
             ! HISTORY.rc (stored in CollectionFrequency).
             !--------------------------------------------------------------
             UpdateYmd = FileWriteYmd
             UpdateHms = FileWriteHms

          ELSE

             !--------------------------------------------------------------
             ! %%%% TIME-AVERAGED COLLECTION %%%%
             !
             ! Define the "Update" interval
             !
             ! Normally, we will set UpdateYmd and UpdateHms directly from
             ! the "heartbeat" timestep of the simulation in seconds.
             !
             ! If the ".acc_interval" tag is specified in HISTORY.rc,
             ! then we will set UpdateYmd and UpdateHms from
             ! CollectionAccInterval.  But if using this option, note
             ! that the ".acc_interval" tag must not specify an interval
             ! that is longer than the interval specified by ".frequency".
             !
             ! NOTE: If CollectionAccInterval is 6 digits long, then assume
             ! that to be UpdateHms.  If longer, then assume that it is
             ! both UpdateYmd and UpdateHms.  This is a hack that we
             ! introduced for GEOS-Chem "Classic" only, as this feature is
             ! not supported in MAPL.  (sde, bmy, 8/4/17, 10/26/17)
             !--------------------------------------------------------------
             IF ( TRIM( CollectionAccInterval(C) ) == UNDEFINED_STR ) THEN

                ! Set UpdateYmd and UpdateHms from the HeartBeat timestep
                UpdateYmd = 00000000
                UpdateHms = HeartBeatHms

                ! SPECIAL CASE: If FileWriteYmd is 240000 then set
                ! and set FileWriteYmd=000001 and FileWriteHms=000000
                IF ( UpdateHms == 240000 ) THEN
                   UpdateYmd = 00000001
                   UpdateHms = 000000
                ENDIF

             ELSE

                ! Set UpdateYmd and UpdateHms from the ".acc_interval" tag
                IF ( LEN_TRIM( CollectionAccInterval(C) ) == 6 ) THEN
                   READ( CollectionAccInterval(C), '(i6.6)'  ) UpdateHms
                ELSE IF ( LEN_TRIM( CollectionAccInterval(C) ) == 14 ) THEN
                   READ( CollectionAccInterval(C), '(i8,i6)' ) UpdateYmd,    &
                                                               UpdateHms
                ENDIF

                ! SPECIAL CASE: If FileWriteYmd is 240000 then set
                ! and set FileWriteYmd=000001 and FileWriteHms=000000
                IF ( UpdateHms == 240000 ) THEN
                   UpdateYmd = 00000001
                   UpdateHms = 000000
                ENDIF

                ! Combine UpdateYmd and UpdateHms
                UpdateCheck    = ( DBLE( UpdateYmd    ) * 1.0e6_f8 )         &
                               + ( DBLE( UpdateHms    )            )

                ! Combine FileWriteYmd and FileWriteHms
                FileWriteCheck = ( DBLE( FileWriteYmd ) * 1.0e6_f8 )         &
                               + ( DBLE( FileWriteHMs )            )

                ! Error check: If using acc_interval, then the Update interval
                ! has to be smaller or equal to the File Write interval
                IF ( UpdateCheck > FileWriteCheck ) THEN
                   ErrMsg = 'Update interval is greater than File Write ' // &
                            'interval for collection: '                   // &
                            TRIM( CollectionName(C) )
                   WRITE( ErrorLine, 250 ) LineNum
                   CALL GC_Error( ErrMsg, RC, ThisLoc, ErrorLine )
                   RETURN
                ENDIF

             ENDIF

          ENDIF

          !=================================================================
          ! Create a HISTORY CONTAINER object for this collection
          !=================================================================

          ! Create the HISTORY CONTAINER object itself.
          ! This will also define the alarm intervals and initial alarm times
          CALL HistContainer_Create( Input_Opt      = Input_Opt,             &
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
                                     FileExpId      = FileExpId,             &
                                     FileWriteYmd   = FileWriteYmd,          &
                                     FileWriteHms   = FileWriteHms,          &
                                     FileCloseYmd   = FileCloseYmd,          &
                                     FileCloseHms   = FileCloseHms,          &
                                     Conventions    = 'COARDS',              &
                                     FileName       = CollectionFileName(C), &
                                     FileTemplate   = CollectionTemplate(C), &
                                     NcFormat       = CollectionFormat(C),   &
                                     Reference      = Reference,             &
                                     Title          = Title,                 &
                                     Contact        = Contact,               &
                                     StartTimeStamp = StartTimeStamp,        &
                                     EndTimeStamp   = EndTimeStamp,          &
                                     RC             = RC                    )

          ! Update CollectionFileName
          CollectionFileName(C) = TRIM( Container%FileName )

          ! Trap potential error
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = 'Could not create Collection: ' // &
                      TRIM( CollectionName(C) )
             WRITE( ErrorLine, 250 ) LineNum
             CALL GC_Error( ErrMsg, RC, ThisLoc, ErrorLine )
             RETURN
          ENDIF

          ! Set elapsed time quantities in the HISTORY CONTAINER object
          CALL HistContainer_SetTime( Input_Opt   = Input_Opt,               &
                                      Container   = Container,               &
                                      HeartBeatDt = 0.0_f8,                  &
                                      RC          = RC                      )


          ! Trap potential error
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = 'Error encountered in "HistContainer_SetTime"'      // &
                      ' for collection: ' // TRIM( CollectionName(C) )
             WRITE( ErrorLine, 250 ) LineNum
             CALL GC_Error( ErrMsg, RC, ThisLoc, ErrorLine )
             RETURN
          ENDIF

          !-----------------------------------------------------------------
          ! ERROR CHECK: Make sure that the length of the simulation is
          ! not shorter than the requested "File Write" interval.  This
          ! will prevent simulations without diagnostic output.
          !-----------------------------------------------------------------
          IF ( SimLengthSec < Container%FileWriteAlarm ) THEN

             ! Construct error message
             ErrMsg =                                                        &
                'No diagnostic output will be created for collection: "'  // &
                 TRIM( CollectionName(C) ) // '"!  Make sure that the '   // &
                'length of the simulation as specified in "input.geos" '  // &
                '(check the start and end times) is not shorter than '    // &
                'the frequency setting in HISTORY.rc!  For example, if '  // &
                'the frequency is 010000 (1 hour) but the simulation '    // &
                'is set up to run for only 20 minutes, then this error '  // &
                'will occur.'

             ! Return error
             WRITE( ErrorLine, 250 ) LineNum
             CALL GC_Error( ErrMsg, RC, ThisLoc, ErrorLine )
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
                CALL GetCollectionMetaData( Input_Opt, Line, 'fields', &
                                            MetaData, C )
                CALL StrSplit( MetaData, " ", Subs1, nSubs1 )
                ItemName = Subs1(1)

             ELSE

                !----------------------------------------------------------
                ! Otherwise, read the next line to get the name for
                ! each subsequent HISTORY ITEM.  The name will be the
                ! first substring of the line.
                !----------------------------------------------------------

                ! Read a single line, and strip leading/trailing spaces
                Line    = ReadOneLine( fId, EOF, IOS, Squeeze=.TRUE. )
                LineNum = LineNum + 1

                ! IF we have hit the end of file then
                iF ( EOF ) GOTO 999

                ! If it's a real I/O error, quit w/ error message
                IF ( IOS > 0 ) THEN
                   ErrMsg = 'Unexpected end-of-file in '        // &
                             TRIM( Input_Opt%HistoryInputFile ) //'!'
                   WRITE( ErrorLine, 250 ) LineNum
                   CALL GC_Error( ErrMsg, RC, ThisLoc, ErrorLine )
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

             ! Save the item name in temporary variables
             ! so that we can parse for wildcards
             ItemTemplate   = ItemName
             ItemTemplateUC = To_UpperCase( ItemTemplate )

             ! Test if there are wildcards present, otherwise skip
             IF ( INDEX( ItemTemplate, '?' ) >  0 ) THEN

                ! Split the name to get wildcard and string prior to wildcard
                CALL StrSplit( ItemTemplate, '?', SubStrs, N )
                tagId = SubStrs(N-1)
                ItemPrefix = SubStrs(1)

                ! Get number of tags for this wildcard
                CALL Get_TagInfo( Input_Opt, tagId, State_Chm, Found, RC, &
                                  nTags=nTags )
                IF ( RC /= GC_SUCCESS ) THEN
                   ErrMsg = 'Error retrieving # of tags for' //              &
                            ' wildcard ' // TRIM(tagId)
                   WRITE( ErrorLine, 250 ) LineNum
                   CALL GC_Error( ErrMsg, RC, ThisLoc, ErrorLine )
                   RETURN
                ENDIF

                ! Add each tagged name as a separate item in the collection
                DO N = 1, nTags
                   ! Construct the item name

                   ! Get tag, if any
                   CALL Get_TagInfo( Input_Opt, tagId, State_Chm, Found, RC, &
                                     N=N, tagName=tagName )
                   IF ( RC /= GC_SUCCESS ) THEN
                      ErrMsg = 'Error retrieving tag name for' //            &
                               ' wildcard ' // TRIM(tagId)
                      WRITE( ErrorLine, 250 ) LineNum
                      CALL GC_Error( ErrMsg, RC, ThisLoc, ErrorLine )
                      RETURN
                   ENDIF

                   ! Append the tag name to the output name
                   ItemName = TRIM( ItemPrefix ) // TRIM( tagName )

                   ! Update the ItemName if dependent on input parameters
                   CALL Get_NameInfo( Input_Opt, ItemName, OutputName, RC )

                   ! Increment the item count
                   ItemCount   = ItemCount + 1

                   ! Create the a HISTORY ITEM object for this diagnostic
                   ! and add it to the given DIAGNOSTIC COLLECTION
                   CALL History_AddItemToCollection(                         &
                            Input_Opt    = Input_Opt,                        &
                            State_Chm    = State_Chm,                        &
                            State_Diag   = State_Diag,                       &
                            State_Met    = State_Met,                        &
                            Collection   = Container,                        &
                            CollectionId = C,                                &
                            SubsetInd    = CollectionSubsetInd(:,C),         &
                            LevelInd     = CollectionLevelInd(:,C),          &
                            ItemName     = OutputName,                       &
                            ItemCount    = ItemCount,                        &
                            RC           = RC                               )

                   ! Error checking
                   IF ( RC /= GC_SUCCESS ) THEN
                      ErrMsg = 'Could not add diagnostic "'               // &
                               TRIM( OutputName ) // '" to collection: '  // &
                               TRIM( CollectionName(C) )
                      WRITE( ErrorLine, 250 ) LineNum
                      CALL GC_Error( ErrMsg, RC, ThisLoc, ErrorLine )
                      RETURN
                   ENDIF
                ENDDO

             ELSE

                !-----------------------------------------------------------
                ! Item name does not have wildcards; no special handling
                !-----------------------------------------------------------

                ! Update the ItemName if dependent on input parameters
                CALL Get_NameInfo( Input_Opt, ItemName, OutputName, RC )

                ! Increment the number of HISTORY items
                ItemCount = ItemCount + 1

                ! Create the a HISTORY ITEM object for this diagnostic
                ! and add it to the given DIAGNOSTIC COLLECTION
                CALL History_AddItemToCollection(                            &
                         Input_Opt    = Input_Opt,                           &
                         State_Chm    = State_Chm,                           &
                         State_Diag   = State_Diag,                          &
                         State_Met    = State_Met,                           &
                         Collection   = Container,                           &
                         CollectionId = C,                                   &
                         SubsetInd    = CollectionSubsetInd(:,C),            &
                         LevelInd     = CollectionLevelInd(:,C),             &
                         ItemName     = OutputName,                          &
                         ItemCount    = ItemCount,                           &
                         RC           = RC                                  )

                ! Trap potential error
                IF ( RC /= GC_SUCCESS ) THEN
                   ErrMsg = 'Could not add diagnostic "' // TRIM(OutputName) &
                            // '" to collection: ' // TRIM( CollectionName(C) )
                   WRITE( ErrorLine, 250 ) LineNum
                   CALL GC_Error( ErrMsg, RC, ThisLoc, ErrorLine )
                   RETURN
                ENDIF
             ENDIF
          ENDDO

          !=================================================================
          ! Add this HISTORY CONTAINER object (i.e. this collection) into
          ! the METAHISTORY OBJECT (i.e. the master list of collections).
          !=================================================================
          CALL MetaHistContainer_AddNew( Input_Opt   = Input_Opt,            &
                                         Node        = CollectionList,       &
                                         Container   = Container,            &
                                         RC          = RC                   )

          ! Trap potential error
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = 'Could not add Container' //                           &
                      TRIM( CollectionName(C) ) //                           &
                      ' to the list of collections!'
             WRITE( ErrorLine, 250 ) LineNum
             CALL GC_Error( ErrMsg, RC, ThisLoc, ErrorLine )
             RETURN
          ENDIF

       ELSE

          !=================================================================
          ! If the we have gotten this far tthrough a collection definition
          ! section, but still haven't found the ".fields" tag, then we
          ! need to do some further error checking.
          !=================================================================
          IF ( C == UNDEFINED_INT ) THEN

             !--------------------------------------------------------------
             ! First, check if the collection isn't activated.  If that
             ! is the case, then skip over all of the lines in the
             ! collection definition section until we hit the "::"
             ! termination character.  Then cycle up to the top of the
             ! loop to read the next collection definition section.
             !--------------------------------------------------------------

             ! Get the collection name (its to the left of the first ".")
             N     = INDEX( TRIM( Line ), '.' )
             CName = Line(1:N-1)

             ! Search for this collection in the list of active collections
             ! and skip to the next collection if not found
             CALL Search_CollList( Input_Opt%amIRoot, CollList, CName, Found, RC )
             IF ( .not. Found ) THEN
                DO
                   Line    = ReadOneLine( fId, EOF, IOS, Squeeze=.TRUE. )
                   LineNum = LineNum + 1
                   CALL StrSqueeze( Line )
                   IF ( TRIM( Line ) == '::' ) GOTO 500
                ENDDO
             ENDIF

             !--------------------------------------------------------------
             ! If we have gotten down to this point, then a true error
             ! condition exists.  Print an error message asking the user
             ! to check the HISTORY.rc file for inconsistencies.
             !--------------------------------------------------------------

             ! List the defined collections
             WRITE( 6, '(/,a)' ) REPEAT( '=', 79 )
             WRITE( 6, 200   )
             DO N = 1, CollectionCount
                WRITE( 6, 210 ) N, TRIM( CollectionName(N) )
             ENDDO
             WRITE( 6, 220 )
             WRITE( 6, '(a,/)' ) REPEAT( '=', 79 )

             ! Write error message and then return
             ErrMsg = 'Inconsistency in collection names and attributes!' // &
                      ' Please check "HISTORY.rc" for typos.'
             WRITE( ErrorLine, 250 ) LineNum
             CALL GC_Error( ErrMsg, RC, ThisLoc, ErrorLine )
             RETURN
          ENDIF
       ENDIF

    ENDDO

    !=======================================================================
    ! Cleanup and quit
    !=======================================================================
999 CONTINUE

    ! Free pointers
    Grid_Lat  => NULL()
    Grid_LatE => NULL()
    Grid_Lon  => NULL()
    Grid_LonE => NULL()

    ! Close the file
    CLOSE( fId )

    ! Write output
    IF ( Input_Opt%amIRoot ) THEN
       WRITE( 6, '(/,a)' ) REPEAT( '=', 79 )
       WRITE( 6, '(a  )' ) 'DEFINED DIAGNOSTIC COLLECTIONS:'
       WRITE( 6, '(a  )' ) REPEAT( '=', 79 )

       DO C = 1, CollectionCount
          print*, 'Collection        ', TRIM( CollectionName       (C) )
          print*, '  -> FileName     ', TRIM( CollectionFileName   (C) )
          print*, '  -> Format       ', TRIM( CollectionFormat     (C) )
          print*, '  -> Frequency    ', TRIM( CollectionFrequency  (C) )
          IF ( CollectionAccInterval(C) /= UNDEFINED_STR ) THEN
             print*, '  -> Acc_Interval ', TRIM( CollectionAccInterval(C) )
          ENDIF
          print*, '  -> Duration     ', TRIM( CollectionDuration   (C) )
          print*, '  -> Mode         ', TRIM( CollectionMode       (C) )
          IF ( CollectionLonRange(C) /= UNDEFINED_STR ) THEN
             print*, '  -> LON_RANGE    ',                                   &
                  TRIM(ADJUSTL(ADJUSTR( CollectionLonRange(C) )))
             print*, '     -> X0 X1  ', ((CollectionSubsetInd(N,C)), N=1,2)
          ENDIF
          IF ( CollectionLatRange(C) /= UNDEFINED_STR ) THEN
             print*, '  -> LAT_RANGE    ',                                   &
                  TRIM(ADJUSTL(ADJUSTR( CollectionLatRange(C) )))
             print*, '     -> Y0 Y1  ', ((CollectionSubsetInd(N,C)), N=3,4)
          ENDIF
          IF ( CollectionLevels(C) /= UNDEFINED_STR ) THEN
             print*, '  -> Levels    ' , TRIM( CollectionLevels(C) )
             print*, '     -> Z0 Z1  ', ((CollectionLevelInd(N,C)), N=1,2)
          ENDIF

          ! Trap error if the collection frequency is undefined
          ! This indicates an error in parsing the file
          IF ( TRIM( CollectionFrequency(C) ) == UNDEFINED_STR ) THEN
             ErrMsg = 'Collection: ' // TRIM( CollectionName(C) ) //         &
                      ' is undefined!'
             WRITE( ErrorLine, 250 ) LineNum
             CALL GC_Error( ErrMsg, RC, ThisLoc, ErrorLine )
             RETURN
          ENDIF
       ENDDO
    ENDIF

    IF ( Input_Opt%LPRT .and. Input_Opt%amIRoot ) THEN
       ! Print information about each diagnostic collection
       CALL MetaHistContainer_Print( Input_Opt, CollectionList, RC )
    ENDIF

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
  SUBROUTINE History_AddItemToCollection( Input_Opt,                         &
                                          State_Chm,    State_Diag,          &
                                          State_Met,    Collection,          &
                                          CollectionId, ItemName,            &
                                          ItemCount,    SubsetInd,           &
                                          LevelInd,     RC                  )
!
! !USES:
!
    USE Charpak_Mod,           ONLY : To_UpperCase
    USE ErrCode_Mod
    USE HistContainer_Mod
    USE HistItem_Mod
    USE History_Util_Mod,      ONLY : UNDEFINED_INT
    USE Input_Opt_Mod,         ONLY : OptInput
    USE MetaHistContainer_Mod
    USE MetaHistItem_Mod
    USE Registry_Mod,          ONLY : Registry_Lookup
    USE State_Chm_Mod
    USE State_Diag_Mod
    USE State_Met_Mod
!
! !INPUT PARAMETERS:
!
    ! Required arguments
    TYPE(OptInput),      INTENT(IN)  :: Input_Opt      ! Input Options State
    TYPE(ChmState),      INTENT(IN)  :: State_Chm      ! Chemistry State
    TYPE(DgnState),      INTENT(IN)  :: State_Diag     ! Diagnostic State
    TYPE(MetState),      INTENT(IN)  :: State_Met      ! Meteorology State
    INTEGER,             INTENT(IN)  :: CollectionID   ! Collection ID number
    CHARACTER(LEN=255),  INTENT(IN)  :: ItemName       ! Name of HISTORY ITEM
    INTEGER,             INTENT(IN)  :: ItemCount      ! Index of HISTORY ITEM

    ! Optional arguments
    INTEGER,             OPTIONAL    :: SubsetInd(4)    ! X0,X1,Y0,Y1 indices
    INTEGER,             OPTIONAL    :: LevelInd(2)     ! Z0,Z1 indices
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
!  See https://github.com/geoschem/geos-chem for complete history
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
    INTEGER                      :: NX, X0, X1
    INTEGER                      :: NY, Y0, Y1
    INTEGER                      :: NZ, Z0, Z1

    ! Arrays
    INTEGER                      :: Dimensions(3)
    INTEGER                      :: ItemDims(3)
    INTEGER                      :: Subset_X(2)
    INTEGER                      :: Subset_Y(2)
    INTEGER                      :: Subset_Z(2)

    ! Strings
    CHARACTER(LEN=4  )           :: StateMetUC
    CHARACTER(LEN=5  )           :: StateChmUC
    CHARACTER(LEN=255)           :: ItemNameUC
    CHARACTER(LEN=255)           :: Description
    CHARACTER(LEN=255)           :: ThisLoc
    CHARACTER(LEN=255)           :: Units
    CHARACTER(LEN=512)           :: ErrMsg

    ! Objects
    TYPE(HistItem),      POINTER :: Item

    ! Pointer arrays
    REAL(fp),            POINTER :: Ptr0d
    REAL(f8),            POINTER :: Ptr0d_8
    REAL(f4),            POINTER :: Ptr0d_4
    INTEGER,             POINTER :: Ptr0d_I
    REAL(fp),            POINTER :: Ptr1d  (:    )
    REAL(f8),            POINTER :: Ptr1d_8(:    )
    REAL(f4),            POINTER :: Ptr1d_4(:    )
    INTEGER,             POINTER :: Ptr1d_I(:    )
    REAL(fp),            POINTER :: Ptr2d  (:,:  )
    REAL(f8),            POINTER :: Ptr2d_8(:,:  )
    REAL(f4),            POINTER :: Ptr2d_4(:,:  )
    INTEGER,             POINTER :: Ptr2d_I(:,:  )
    REAL(fp),            POINTER :: Ptr3d  (:,:,:)
    REAL(f8),            POINTER :: Ptr3d_8(:,:,:)
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
    ItemNameUC  = To_UpperCase( ItemName )
    StateMetUC  = State_Met%State // '_'   ! State_Met%State is uppercase
    StateChmUC  = State_Chm%State // '_'   ! State_Chm%State is uppercase
    Ptr0d       => NULL()
    Ptr0d_8     => NULL()
    Ptr0d_4     => NULL()
    Ptr0d_I     => NULL()
    Ptr1d       => NULL()
    Ptr1d_8     => NULL()
    Ptr1d_4     => NULL()
    Ptr1d_I     => NULL()
    Ptr2d       => NULL()
    Ptr2d_8     => NULL()
    Ptr2d_4     => NULL()
    Ptr2d_I     => NULL()
    Ptr3d       => NULL()
    Ptr3d_8     => NULL()
    Ptr3d_4     => NULL()
    Ptr3d_I     => NULL()

    !=======================================================================
    ! For each HISTORY ITEM, find the matching entry in the relevant
    ! registry (in State_Chm, State_Diag, State_Met) and get a pointer
    ! to the data source
    !=======================================================================
    IF ( ItemNameUC(1:5) == StateChmUC ) THEN

       !--------------------------------------------------------------------
       ! Chemistry State
       !--------------------------------------------------------------------
       CALL Registry_Lookup( am_I_Root    = Input_Opt%amIRoot,               &
                             Registry     = State_Chm%Registry,              &
                             State        = State_Chm%State,                 &
                             Variable     = ItemName,                        &
                             Description  = Description,                     &
                             Dimensions   = Dimensions,                      &
                             KindVal      = KindVal,                         &
                             Rank         = Rank,                            &
                             Units        = Units,                           &
                             OnLevelEdges = OnLevelEdges,                    &
                             Ptr0d        = Ptr0d,                           &
                             Ptr1d        = Ptr1d,                           &
                             Ptr2d        = Ptr2d,                           &
                             Ptr3d        = Ptr3d,                           &
                             Ptr0d_8      = Ptr0d_8,                         &
                             Ptr1d_8      = Ptr1d_8,                         &
                             Ptr2d_8      = Ptr2d_8,                         &
                             Ptr3d_8      = Ptr3d_8,                         &
                             Ptr0d_4      = Ptr0d_4,                         &
                             Ptr1d_4      = Ptr1d_4,                         &
                             Ptr2d_4      = Ptr2d_4,                         &
                             Ptr3d_4      = Ptr3d_4,                         &
                             Ptr0d_I      = Ptr0d_I,                         &
                             Ptr1d_I      = Ptr1d_I,                         &
                             Ptr2d_I      = Ptr2d_I,                         &
                             Ptr3d_I      = Ptr3d_I,                         &
                             RC           = RC                                 )

       ! Trap potential not found error
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Could not locate ' // TRIM( ItemName )  // &
                   ' chemistry state registry.'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

    ELSE IF ( ItemNameUC(1:4) == StateMetUC ) THEN

       !--------------------------------------------------------------------
       ! Meteorology State
       !--------------------------------------------------------------------
       CALL Registry_Lookup( am_I_Root    = Input_Opt%amIRoot,               &
                             Registry     = State_Met%Registry,              &
                             State        = State_Met%State,                 &
                             Variable     = ItemName,                        &
                             Description  = Description,                     &
                             Dimensions   = Dimensions,                      &
                             KindVal      = KindVal,                         &
                             Rank         = Rank,                            &
                             Units        = Units,                           &
                             OnLevelEdges = OnLevelEdges,                    &
                             Ptr0d        = Ptr0d,                           &
                             Ptr1d        = Ptr1d,                           &
                             Ptr2d        = Ptr2d,                           &
                             Ptr3d        = Ptr3d,                           &
                             Ptr0d_8      = Ptr0d_8,                         &
                             Ptr1d_8      = Ptr1d_8,                         &
                             Ptr2d_8      = Ptr2d_8,                         &
                             Ptr3d_8      = Ptr3d_8,                         &
                             Ptr0d_4      = Ptr0d_4,                         &
                             Ptr1d_4      = Ptr1d_4,                         &
                             Ptr2d_4      = Ptr2d_4,                         &
                             Ptr3d_4      = Ptr3d_4,                         &
                             Ptr0d_I      = Ptr0d_I,                         &
                             Ptr1d_I      = Ptr1d_I,                         &
                             Ptr2d_I      = Ptr2d_I,                         &
                             Ptr3d_I      = Ptr3d_I,                         &
                             RC           = RC                                 )

       ! Trap potential not found error
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Could not locate ' // TRIM( ItemName )  // &
                   ' meteorology state registry.'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

    ELSE

       !--------------------------------------------------------------------
       ! Diagnostic State
       !--------------------------------------------------------------------
       CALL Registry_Lookup( am_I_Root    = Input_Opt%amIRoot,               &
                             Registry     = State_Diag%Registry,             &
                             State        = State_Diag%State,                &
                             Variable     = ItemName,                        &
                             Description  = Description,                     &
                             Dimensions   = Dimensions,                      &
                             KindVal      = KindVal,                         &
                             Rank         = Rank,                            &
                             Units        = Units,                           &
                             OnLevelEdges = OnLevelEdges,                    &
                             Ptr0d        = Ptr0d,                           &
                             Ptr1d        = Ptr1d,                           &
                             Ptr2d        = Ptr2d,                           &
                             Ptr3d        = Ptr3d,                           &
                             Ptr0d_8      = Ptr0d_8,                         &
                             Ptr1d_8      = Ptr1d_8,                         &
                             Ptr2d_8      = Ptr2d_8,                         &
                             Ptr3d_8      = Ptr3d_8,                         &
                             Ptr0d_4      = Ptr0d_4,                         &
                             Ptr1d_4      = Ptr1d_4,                         &
                             Ptr2d_4      = Ptr2d_4,                         &
                             Ptr3d_4      = Ptr3d_4,                         &
                             Ptr0d_I      = Ptr0d_I,                         &
                             Ptr1d_I      = Ptr1d_I,                         &
                             Ptr2d_I      = Ptr2d_I,                         &
                             Ptr3d_I      = Ptr3d_I,                         &
                             RC           = RC                                 )

       ! Trap potential not found error
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Could not locate ' // TRIM( ItemName )  // &
                   ' diagnostics state registry.'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDIF

    !=======================================================================
    ! If the optional SUBSETS and/or LEVELS arguments are passed, then use
    ! these to size the data arrays.  Otherwise assume the data arrays will
    ! be the same size as the pointer to the data source (as will be true
    ! in most cases.)
    !=======================================================================

    !-------------------------
    ! Default values
    !-------------------------

    ! By default, use the size of the data array to define the
    ! X0, Y0, X1, Y1, Z0, and Z1 indices for the subset region.
    X0 = 1
    X1 = MAX( Dimensions(1), 1 )
    Y0 = 1
    Y1 = MAX( Dimensions(2), 1 )
    Z0 = 1
    Z1 = MAX( Dimensions(3), 1 )

    !-------------------------
    ! Horizontal subsetting
    !-------------------------

    ! If SubsetInd has valid values, use them to redefine X0, Y0, X1, and Y1.
    IF ( PRESENT( SubsetInd ) ) THEN
       IF ( SubsetInd(1) /= UNDEFINED_INT ) X0 = SubsetInd(1)
       IF ( SubsetInd(2) /= UNDEFINED_INT ) X1 = SubsetInd(2)
       IF ( SubsetInd(3) /= UNDEFINED_INT ) Y0 = SubsetInd(3)
       IF ( SubsetInd(4) /= UNDEFINED_INT ) Y1 = SubsetInd(4)
    ENDIF

    !-------------------------
    ! Vertical subsetting
    !-------------------------

    ! If LevelInd has valid values, use them to redefine Z0 and Z1.
    IF ( PRESENT( LevelInd ) ) THEN
       IF ( LevelInd(1) /= UNDEFINED_INT ) Z0 = LevelInd(1)
       IF ( LevelInd(2) /= UNDEFINED_INT ) Z1 = LevelInd(2)
    ENDIF

    ! Error check X-dimension indices
    IF ( X1 < X0 ) THEN
       WRITE( ErrMsg, 100 ) X0, X1, TRIM( Collection%Name )
 100   FORMAT(  'Invalid X-dimension indices: ', 2i6, ' for collection', a )
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    ! Error-check Y-dimension indices
    IF ( Y1 < Y0 ) THEN
       WRITE( ErrMsg, 110 ) Y0, Y1, TRIM( Collection%Name )
 110   FORMAT(  'Invalid Y-dimension indices: ', 2i6, ' for collection', a )
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    ! Error-check Z-dimension indices
    IF ( Z1 < Z0 ) THEN
       WRITE( ErrMsg, 120 ) Z0, Z1, TRIM( Collection%Name )
 120   FORMAT(  'Invalid Y-dimension indices: ', 2i6, ' for collection', a )
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    ! Compute dimension extent
    NX = X1 - X0 + 1
    NY = Y1 - Y0 + 1
    NZ = Z1 - Z0 + 1

    ! Indices for subsetting the data
    Subset_X = (/ X0, X1 /)
    Subset_Y = (/ Y0, Y1 /)
    Subset_Z = (/ Z0, Z1 /)

    ! Save the subsets to the collection
    Collection%X0 = X0
    Collection%X1 = X1
    Collection%Y0 = Y0
    Collection%Y1 = Y1
    Collection%Z0 = Z0
    !NOTE: Z1 is not needed, we compute that later!

    !=======================================================================
    ! Now that we have obtained information (and pointers to the data)
    ! corresponding to the given diagnostic quantity, use that to create
    ! a HISTORY ITEM object for that diagnostic quantity.
    !=======================================================================
    CALL HistItem_Create( Input_Opt      = Input_Opt,                       &
                          Item           = Item,                             &
                          Id             = ItemCount,                        &
                          ContainerId    = CollectionId,                     &
                          Name           = ItemName,                         &
                          LongName       = Description,                      &
                          Units          = Units,                            &
                          OnLevelEdges   = OnLevelEdges,                     &
                          SpaceDim       = Rank,                             &
                          Operation      = Collection%Operation,             &
                          Subset_X       = Subset_X,                         &
                          Subset_Y       = Subset_Y,                         &
                          Subset_Z       = Subset_Z,                         &
                          Source_KindVal = KindVal,                          &
                          Source_0d_8    = Ptr0d_8,                          &
                          Source_1d      = Ptr1d,                            &
                          Source_1d_8    = Ptr1d_8,                          &
                          Source_1d_4    = Ptr1d_4,                          &
                          Source_1d_I    = Ptr1d_I,                          &
                          Source_2d      = Ptr2d,                            &
                          Source_2d_8    = Ptr2d_8,                          &
                          Source_2d_4    = Ptr2d_4,                          &
                          Source_2d_I    = Ptr2d_I,                          &
                          Source_3d      = Ptr3d,                            &
                          Source_3d_8    = Ptr3d_8,                          &
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
    CALL MetaHistItem_AddNew( Input_Opt = Input_Opt,                         &
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
          ErrMsg = TRIM( Item%Name )                                      // &
                   ' has the wrong vertical alignment for collection: "'  // &
                   TRIM( Collection%Name )  // '".  Please check your '   // &
                   'HISTORY.rc file to make sure that this collection '   // &
                   'only contains 3-D diagnostics with the same vertical '// &
                   'alignment.  You cannot add diagnostics that are '     // &
                   'defined on level centers and diagnostics that are '   // &
                   'defined on level edges in the same collection, as '   // &
                   'per netCDF conventions.'
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
    Ptr0d   => NULL()
    Ptr0d_8 => NULL()
    Ptr0d_4 => NULL()
    Ptr0d_I => NULL()
    Ptr1d   => NULL()
    Ptr1d_8 => NULL()
    Ptr1d_4 => NULL()
    Ptr1d_I => NULL()
    Ptr2d   => NULL()
    Ptr2d_8 => NULL()
    Ptr2d_4 => NULL()
    Ptr2d_I => NULL()
    Ptr3d   => NULL()
    Ptr3d_8 => NULL()
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
  SUBROUTINE History_SetTime( Input_Opt, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE HistContainer_Mod,     ONLY : HistContainer_SetTime
    USE History_Util_Mod
    USE Input_Opt_Mod,         ONLY : OptInput
    USE MetaHistContainer_Mod, ONLY : MetaHistContainer
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)  :: Input_Opt        ! Input Options object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT) :: RC               ! Success or failure
!
! !REMARKS:
!  This routine is meant to be called after History_Update() but before
!  History_Write().
!
! !REVISION HISTORY:
!  18 Aug 2017 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
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
       CALL HistContainer_SetTime( Input_Opt   = Input_Opt,                  &
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
  SUBROUTINE History_Update( Input_Opt, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE HistItem_Mod,          ONLY : HistItem
    USE HistContainer_Mod,     ONLY : HistContainer
    USE HistContainer_Mod,     ONLY : HistContainer_UpdateIvalSet
    USE History_Util_Mod
    USE Input_Opt_Mod,         ONLY : OptInput
    USE MetaHistContainer_Mod, ONLY : MetaHistContainer
    USE MetaHistItem_Mod,      ONLY : MetaHistItem
    USE Registry_Params_Mod
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
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
!  See https://github.com/geoschem/geos-chem for complete history
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

       ! Debug output
       IF ( Input_Opt%LPRT .and. Input_Opt%amIRoot ) THEN
          WRITE( 6, 100 ) Container%Name
 100      FORMAT( '     - Updating collection: ', a20 )
       ENDIF

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

                ! 8-byte floating point
                ELSE IF ( Item%Source_KindVal == KINDVAL_F8 ) THEN

                   IF ( Item%Operation == COPY_FROM_SOURCE ) THEN
                      Item%Data_3d = Item%Source_3d_8
                      Item%nUpdates = 1.0_f8
                   ELSE
                      Item%Data_3d  = Item%Data_3d  + Item%Source_3d_8
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

                ! 8-byte floating point
                ELSE IF ( Item%Source_KindVal == KINDVAL_F8 ) THEN

                   IF ( Item%Operation == COPY_FROM_SOURCE ) THEN
                      Item%Data_2d  = Item%Source_2d_8
                      Item%nUpdates = 1.0_f8
                   ELSE
                      Item%Data_2d  = Item%Data_2d + Item%Source_2d_8
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

                ! 8-byte floating point
                ELSE IF ( Item%Source_KindVal == KINDVAL_F8 ) THEN

                   IF ( Item%Operation == COPY_FROM_SOURCE ) THEN
                      Item%Data_1d  = Item%Source_1d_8
                      Item%nUpdates = 1.0_f8
                   ELSE
                      Item%Data_1d  = Item%Data_1d  + Item%Source_1d_8
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

             !--------------------------------------------------------------
             ! Update 0D data field
             !--------------------------------------------------------------
             CASE( 0 )

                ! Flex-precision floating point
                IF ( Item%Source_KindVal == KINDVAL_F8 ) THEN

                   IF ( Item%Operation == COPY_FROM_SOURCE ) THEN
                      Item%Data_0d  = Item%Source_0d_8
                      Item%nUpdates = 1.0_f8
                   ELSE
                      Item%Data_0d  = Item%Data_0d  + Item%Source_0d_8
                      Item%nUpdates = Item%nUpdates + 1.0_f8
                   ENDIF

                ENDIF

          END SELECT

! Uncomment more detailed debug output if you need it!
!          ! Debug output
!          IF ( Input_Opt%LPRT .and. Input_Opt%amIRoot ) THEN
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

       ! Recompute the update alarm interval if it 1 month or longer,
       ! as we will have to take into account leap years, etc.
       IF ( Container%UpdateYmd >= 000100 ) THEN
          CALL HistContainer_UpdateIvalSet( Input_Opt, Container, RC )
       ENDIF

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
  SUBROUTINE History_Write( Input_Opt, Spc_Units, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE HistContainer_Mod
    USE HistItem_Mod,          ONLY : HistItem
    USE History_Netcdf_Mod
    USE History_Util_Mod
    USE Input_Opt_Mod,         ONLY : OptInput
    USE MetaHistContainer_Mod, ONLY : MetaHistContainer
    USE MetaHistItem_Mod,      ONLY : MetaHistItem
    USE Registry_Params_Mod
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput),   INTENT(IN)  :: Input_Opt   ! Input Options object
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
!  See https://github.com/geoschem/geos-chem for complete history
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
          CALL History_Netcdf_Close( Container = Container,                  &
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
          CALL History_Netcdf_Define( Input_Opt = Input_Opt,                &
                                      Container = Container,                &
                                      RC        = RC                       )

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
             CALL HistContainer_FileCloseIvalSet( Input_Opt, Container, RC )
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
          CALL History_Netcdf_Write( Input_Opt = Input_Opt,                  &
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
             CALL HistContainer_FileWriteIvalSet( Input_Opt, Container, RC )
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
! !IROUTINE: GetCollectionMetaData
!
! !DESCRIPTION: Parses a line of the HISTORY.rc file and returns metadata
!  for a given attribute (e.g. "frequency", "template", etc.)
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GetCollectionMetaData( Input_Opt, Line, Pattern, MetaData, &
                                    nCollection )
!
! !USES:
!
    USE Charpak_Mod,      ONLY : CleanText, StrSplit
    USE DiagList_Mod,     ONLY : CollList,  Search_CollList
    USE History_Util_Mod
    USE Input_Opt_Mod,    ONLY : OptInput
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput),     INTENT(IN)  :: Input_Opt     ! Input Options object
    CHARACTER(LEN=*),   INTENT(IN)  :: Line          ! Line to be searched
    CHARACTER(LEN=*),   INTENT(IN)  :: Pattern       ! Search pattern
!
! !OUTPUT PARAMETERS:
!
    CHARACTER(LEN=255), INTENT(OUT) :: MetaData      ! Metadata value
    INTEGER,            INTENT(OUT) :: nCollection   ! Collection Id
!
! !REVISION HISTORY:
!  16 Aug 2017 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    LOGICAL                  :: Found
    INTEGER                  :: C, Ind, nSubStr, N, P, RC

    ! Strings
    CHARACTER(LEN=255)       :: Name
    CHARACTER(LEN=255)       :: SubStr(255)

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

    ! Exit if the collection name is not in the list of active collections
    CALL Search_CollList( Input_Opt%amIRoot, CollList, Name, Found, RC )
    IF ( .not. Found ) RETURN

    ! Non-white-space lengths of the collection cname and search pattern
    N = LEN_TRIM( Name    )
    P = LEN_TRIM( Pattern )

    ! Loop over all collection names
    ! NOTE: This algorithm may not be the most efficient, as it will
    ! not skip collections that we have already encountered.  But it
    ! only gets done during the init phase, so it might not be a huge
    ! expenditure of time anyway.  Worry about this later. (bmy, 1/18/18)
    DO C = 1, CollectionCount

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
  SUBROUTINE History_Close_AllFiles( RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE History_Netcdf_Mod,    ONLY : History_Netcdf_Close
    USE MetaHistContainer_Mod, ONLY : MetaHistContainer
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT) :: RC          ! Success or failure
!
! !REMARKS:
!  This is called from History_Cleanup, but may also be called in other
!  locations (e.g. when processing abnormal exits)
!
! !REVISION HISTORY:
!  16 Aug 2017 - R. Yantosca - Initial version
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
       CALL History_Netcdf_Close( Container = Current%Container, &
                                  RC        = RC                 )

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
  SUBROUTINE History_Cleanup( RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE MetaHistContainer_Mod, ONLY : MetaHistContainer_Destroy
!
! !OUTPUT PARAMETERS:
!
     INTEGER,       INTENT(OUT) :: RC
!
! !REVISION HISTORY:
!  16 Jun 2017 - R. Yantosca - Initial version
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
     CALL History_Close_AllFiles( RC )
     IF ( RC /= GC_SUCCESS ) THEN
        ErrMsg = 'Error returned from "History_Close_AllFiles"!'
        CALL GC_Error( ErrMsg, RC, ThisLoc )
        RETURN
     ENDIF

     !======================================================================
     ! And deallocate variables belonging to history_mod.F90
     !======================================================================
     IF ( ASSOCIATED( CollectionList ) ) THEN
        CALL MetaHistContainer_Destroy( CollectionList, RC )
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

     IF ( ALLOCATED( CollectionFileName ) ) THEN
        DEALLOCATE( CollectionFileName, STAT=RC )
        IF ( RC /= GC_SUCCESS ) THEN
           ErrMsg = 'Could not deallocate "CollectionFileName"!'
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

     IF ( ALLOCATED( CollectionLonRange ) ) THEN
        DEALLOCATE( CollectionLonRange, STAT=RC )
        IF ( RC /= GC_SUCCESS ) THEN
           ErrMsg = 'Could not deallocate "CollectionLonRange"!'
           CALL GC_Error( ErrMsg, RC, ThisLoc )
           RETURN
        ENDIF
     ENDIF

     IF ( ALLOCATED( CollectionLatRange ) ) THEN
        DEALLOCATE( CollectionLatRange, STAT=RC )
        IF ( RC /= GC_SUCCESS ) THEN
           ErrMsg = 'Could not deallocate "CollectionLatRange"!'
           CALL GC_Error( ErrMsg, RC, ThisLoc )
           RETURN
        ENDIF
     ENDIF

     IF ( ALLOCATED( CollectionSubsetInd ) ) THEN
        DEALLOCATE( CollectionSubsetInd, STAT=RC )
        IF ( RC /= GC_SUCCESS ) THEN
           ErrMsg = 'Could not deallocate "CollectionSubsetInd"!'
           CALL GC_Error( ErrMsg, RC, ThisLoc )
           RETURN
        ENDIF
     ENDIF

     IF ( ALLOCATED( CollectionLevels ) ) THEN
        DEALLOCATE( CollectionLevels, STAT=RC )
        IF ( RC /= GC_SUCCESS ) THEN
           ErrMsg = 'Could not deallocate "CollectionLevels"!'
           CALL GC_Error( ErrMsg, RC, ThisLoc )
           RETURN
        ENDIF
     ENDIF

     IF ( ALLOCATED( CollectionLevelInd ) ) THEN
        DEALLOCATE( CollectionLevelInd, STAT=RC )
        IF ( RC /= GC_SUCCESS ) THEN
           ErrMsg = 'Could not deallocate "CollectionLevelInd"!'
           CALL GC_Error( ErrMsg, RC, ThisLoc )
           RETURN
        ENDIF
     ENDIF

   END SUBROUTINE History_Cleanup
!EOC
END MODULE History_Mod
