!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: histcontainer_mod.F90
!
! !DESCRIPTION: Contains methods to create a HISTORY CONTAINER object.
!  A HISTORY CONTAINER represents a collection of HISTORY ITEMS that will
!  be archived to a netCDF file at a specific temporal frequencly (e.g.
!  instantaneous, hourly, daily, monthly, end-of-run, etc.)
!\\
!\\
!  In other words, the HISTORY CONTAINER provides metadata for the
!  netCDF file, and the HISTORY ITEMS belonging to the HISTORY CONTAINER
!  contains the data and attributes for each variable that will be
!  saved to the netCDF file.
!
! !INTERFACE:
!
MODULE HistContainer_Mod
!
! !USES:
!
  USE MetaHistItem_Mod,  ONLY: MetaHistItem
  USE Precision_Mod

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC  :: HistContainer_Create
  PUBLIC  :: HistContainer_Print
  PUBLIC  :: HistContainer_Destroy
  PUBLIC  :: HistContainer_SetTime
  PUBLIC  :: HistContainer_UpdateIvalSet
  PUBLIC  :: HistContainer_FileCloseIvalSet
  PUBLIC  :: HistContainer_FileWriteIvalSet
!
! !PRIVATE MEMBER FUNCTIONS:
!
  PRIVATE :: AlarmIncrementMonths
  PRIVATE :: AlarmIncrementYears
!
! !PUBLIC TYPES:
!
  !=========================================================================
  ! This is the derived type for a single HISTORY CONTAINER OBJECT, which
  ! contains several HISTORY ITEMS that will be archived at a specified
  ! frequency (e.g. instantaneous, hourly, daily, etc) to netCDF output.
  !=========================================================================
  TYPE, PUBLIC :: HistContainer

     !----------------------------------------------------------------------
     ! Identifying information
     !----------------------------------------------------------------------
     CHARACTER(LEN=255)          :: Name                ! Container name
     INTEGER                     :: Id                  !  and ID number
     INTEGER                     :: nX                  ! X (or lon) dim size
     INTEGER                     :: nY                  ! Y (or lat) dim size
     INTEGER                     :: nZ                  ! Z (or lev) dim size
     INTEGER                     :: X0, X1              ! X (or lon) indices
     INTEGER                     :: Y0, Y1              ! Y (or lon) indices
     INTEGER                     :: Z0                  ! Z (or lev) indices
     LOGICAL                     :: OnLevelEdges        ! =T if data is defined
                                                        !    on level edges;
                                                        ! =F if on centers

     !----------------------------------------------------------------------
     ! List of history items in this collection
     !----------------------------------------------------------------------
     TYPE(MetaHistItem), POINTER :: HistItems => NULL() ! List and # of
     INTEGER                     :: nHistItems          !  HISTORY ITEMS
                                                        !  in this container

     !----------------------------------------------------------------------
     ! Time quantities measured since start of simulation
     !----------------------------------------------------------------------
     REAL(f8)                    :: EpochJd             ! Astronomical Julian
                                                        !  date @ start of sim
     REAL(f8)                    :: EpochJsec           ! Astronomical Julian
                                                        !  secs @ start of sim
     INTEGER                     :: CurrentYmd          ! Current YMD date
     INTEGER                     :: CurrentHms          ! Current hms time
     REAL(f8)                    :: CurrentJd           ! Astronomical Julian
                                                        !  date @ current time
     REAL(f8)                    :: CurrentJsec         ! Astronomical Julian
                                                        !  secs @ current time
     REAL(f8)                    :: ElapsedSec          ! Elapsed seconds
                                                        !  since start of sim
     REAL(f8)                    :: UpdateAlarm         ! Alarm (elapsed sec)
                                                        !  for data updating
     REAL(f8)                    :: FileCloseAlarm      ! Alarm (elapsed sec)
                                                        !  for file close/open
     REAL(f8)                    :: FileWriteAlarm      ! Alarm (elapsed sec)
                                                        !  for file write

     !----------------------------------------------------------------------
     ! Time quantities measured since the time of netCDF file creation
     !----------------------------------------------------------------------
     INTEGER                     :: ReferenceYmd        ! Reference YMD & hms
     INTEGER                     :: ReferenceHms        !  for the "time" dim
     REAL(f8)                    :: ReferenceJd         ! Julian Date at the
                                                        !  reference YMD & hms
     REAL(f8)                    :: ReferenceJsec       ! Julian Seconds @ the
                                                        !  reference YMD & hms
     INTEGER                     :: CurrTimeSlice       ! Current time slice
                                                        !  for the "time" dim
     REAL(f8)                    :: TimeStamp           ! Elapsed minutes w/r/t
                                                        !  reference YMD & hms

     !----------------------------------------------------------------------
     ! Quantities that govern the updating/time averaging of data
     !----------------------------------------------------------------------
     CHARACTER(LEN=255)          :: UpdateMode          ! e.g. inst or time-avg
     INTEGER                     :: UpdateYmd           ! Update frequency
     INTEGER                     :: UpdateHms           !  in YMD and hms
     REAL(f8)                    :: UpdateIvalSec       ! Update interval [sec]
     INTEGER                     :: Operation           ! Operation code
                                                        !  0=copy from source
                                                        !  1=accum from source
     REAL(f8)                    :: HeartBeatDtSec      ! The "heartbeat"
                                                        !  timestep [sec]

     !----------------------------------------------------------------------
     ! Quantities for file creation, writing, and I/O status
     !----------------------------------------------------------------------
     INTEGER                     :: FileWriteYmd        ! File write frequency
     INTEGER                     :: FileWriteHms        !  in YMD and hms
     REAL(f8)                    :: FileWriteIvalSec    ! File write interval
                                                        !  in seconds

     INTEGER                     :: FileCloseYmd        ! File closing time
     INTEGER                     :: FileCloseHms        !  in YMD and hms
     REAL(f8)                    :: FileCloseIvalSec    ! File close interval
                                                        !  in seconds

     LOGICAL                     :: IsFileDefined       ! Have we done netCDF
                                                        !  define mode yet?
     LOGICAL                     :: IsFileOpen          ! Is the netCDF file
                                                        !  currently open?

     !----------------------------------------------------------------------
     ! netCDF file identifiers and attributes
     !----------------------------------------------------------------------
     LOGICAL                     :: FirstInst           ! 1st inst file write?
     INTEGER                     :: FileId              ! netCDF file ID
     INTEGER                     :: xDimId              ! X (or lon ) dim ID
     INTEGER                     :: yDimId              ! Y (or lat ) dim ID
     INTEGER                     :: zDimId              ! Z (or lev ) dim ID
     INTEGER                     :: iDimId              ! I (or ilev) dim ID
     INTEGER                     :: tDimId              ! T (or time) dim ID
     CHARACTER(LEN=20)           :: StartTimeStamp      ! Timestamps at start
     CHARACTER(LEN=20)           :: EndTimeStamp        !  and end of sim
     CHARACTER(LEN=20)           :: Spc_Units           ! Units of SC%Species
     CHARACTER(LEN=255)          :: FileExpId           ! Filename ExpId
     CHARACTER(LEN=255)          :: FilePrefix          ! Filename prefix
     CHARACTER(LEN=255)          :: FileTemplate        ! YMDhms template
     CHARACTER(LEN=255)          :: FileName            ! Name of nc file
     CHARACTER(LEN=255)          :: Conventions         ! e.g. "COARDS"
     CHARACTER(LEN=255)          :: NcFormat            ! e.g. "netCDF-4"
     CHARACTER(LEN=255)          :: History             ! History
     CHARACTER(LEN=255)          :: ProdDateTime        ! When produced
     CHARACTER(LEN=255)          :: Reference           ! Reference string
     CHARACTER(LEN=255)          :: Contact             ! Contact string
     CHARACTER(LEN=255)          :: Title               ! Title string

  END TYPE HistContainer
!
! !REMARKS:
!  Linked list routines taken from original code (linkedlist.f90)
!  by Arjen Markus; http://flibs.sourceforge.net/linked_list.html
!
! !REVISION HISTORY:
!  12 Jun 2017 - R. Yantosca - Initial version, based on history_list_mod.F90
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
! !IROUTINE: HistContainer_Create
!
! !DESCRIPTION: Initializes a single HISTORY CONTAINER object, which
!  will hold a METAHISTORY ITEM (which is a list of HISTORY ITEMS), to
!  archive to netCDF output.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HistContainer_Create( Input_Opt,      Container,                &
                                   Id,             Name,                     &
                                   RC,             EpochJd,                  &
                                   CurrentYmd,     CurrentHms,               &
                                   UpdateMode,     UpdateYmd,                &
                                   UpdateHms,      UpdateAlarm,              &
                                   Operation,      HeartBeatDtSec,           &
                                   FileWriteYmd,   FileWriteHms,             &
                                   FileWriteAlarm, FileCloseYmd,             &
                                   FileCloseHms,   FileCloseAlarm,           &
                                   FileId,         FileExpId,                &
                                   FilePrefix,     FileName,                 &
                                   FileTemplate,   Conventions,              &
                                   NcFormat,       History,                  &
                                   ProdDateTime,   Reference,                &
                                   Title,          Contact,                  &
                                   StartTimeStamp, EndTimeStamp             )
!
! !USES:
!
    USE ErrCode_Mod
    USE History_Util_Mod
    USE Input_Opt_Mod,    ONLY : OptInput
    USE MetaHistItem_Mod, ONLY : MetaHistItem
!
! !INPUT PARAMETERS:
!
    !-----------------------------------------------------------------------
    ! REQUIRED INPUTS
    !-----------------------------------------------------------------------
    TYPE(OptInput),      INTENT(IN)  :: Input_Opt      ! Input Options object
    INTEGER,             INTENT(IN)  :: Id             ! Container Id #
    CHARACTER(LEN=*),    INTENT(IN)  :: Name           ! Container name

    !-----------------------------------------------------------------------
    ! OPTIONAL INPUTS: Time and date quantities
    !-----------------------------------------------------------------------
    REAL(f8),            OPTIONAL    :: EpochJd        ! Astronomical Julian
                                                       !  date @ start of sim
    INTEGER,             OPTIONAL    :: CurrentYmd     ! Current YMD date
    INTEGER,             OPTIONAL    :: CurrentHms     ! Current hms time

    !-----------------------------------------------------------------------
    ! OPTIONAL INPUTS: quantities controlling data updates
    !-----------------------------------------------------------------------
    CHARACTER(LEN=*),    OPTIONAL    :: UpdateMode     ! e.g. inst or time-avg
    INTEGER,             OPTIONAL    :: UpdateYmd      ! Update frequency
    INTEGER,             OPTIONAL    :: UpdateHms      !  in both YMD and hms
    REAL(f8),            OPTIONAL    :: UpdateAlarm    ! JD for data update
    INTEGER,             OPTIONAL    :: Operation      ! Operation code:
                                                       !  0=copy  from source
                                                       !  1=accum from source
    REAL(f8),            OPTIONAL    :: HeartBeatDtSec ! Model "heartbeat"
                                                       !  timestep [sec]

    !-----------------------------------------------------------------------
    ! OPTIONAL INPUTS: quantities controlling file write and close/reopen
    !-----------------------------------------------------------------------
    INTEGER,             OPTIONAL    :: FileWriteYmd   ! File write frequency
    INTEGER,             OPTIONAL    :: FileWriteHms   !  in both YMD and hms
    REAL(f8),            OPTIONAL    :: FileWriteAlarm ! JD for file write

    INTEGER,             OPTIONAL    :: FileCloseYmd   ! File close/open freq
    INTEGER,             OPTIONAL    :: FileCloseHms   !  in both YMD and hm
    REAL(f8),            OPTIONAL    :: FileCloseAlarm ! JD for file close

    !-----------------------------------------------------------------------
    ! OPTIONAL INPUTS: netCDF file identifiers and metadata
    !-----------------------------------------------------------------------
    INTEGER,             OPTIONAL    :: FileId         ! netCDF file ID
    CHARACTER(LEN=*),    OPTIONAL    :: FileExpId      ! Dir name + file string
    CHARACTER(LEN=*),    OPTIONAL    :: FilePrefix     ! Filename prefix
    CHARACTER(LEN=*),    OPTIONAL    :: FileTemplate   ! YMDhms template
    CHARACTER(LEN=*),    OPTIONAL    :: Conventions    ! e.g. "COARDS"
    CHARACTER(LEN=*),    OPTIONAL    :: Filename       ! Name of nc file
    CHARACTER(LEN=*),    OPTIONAL    :: NcFormat       ! e.g. "netCDF-4"
    CHARACTER(LEN=*),    OPTIONAL    :: History        ! History
    CHARACTER(LEN=*),    OPTIONAL    :: ProdDateTime   ! When produced
    CHARACTER(LEN=*),    OPTIONAL    :: Reference      ! Reference string
    CHARACTER(LEN=*),    OPTIONAL    :: Title          ! Title string
    CHARACTER(LEN=*),    OPTIONAL    :: Contact        ! Contact string
    CHARACTER(LEN=*),    OPTIONAL    :: StartTimeStamp ! Timestamps at start
    CHARACTER(LEN=*),    OPTIONAL    :: EndTimeStamp   !  & end of simulation
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(HistContainer), POINTER     :: Container      ! Collection object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,             INTENT(OUT) :: RC             ! Success or failure
!
! !REMARKS:
!  (1) We need to copy string data to a temporary string of length 255
!       characters, or else Gfortran will choke.
!
! !REVISION HISTORY:
!  16 Jun 2017 - R. Yantosca - Initial version, based on history_list_mod.F90
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER            :: ThisId, C

    ! Strings
    CHARACTER(LEN=255) :: ErrMsg, ThisLoc, TempStr

    !========================================================================
    ! Initialize
    !========================================================================

    ! Set initial values
    RC      = GC_SUCCESS
    ErrMsg  = ''
    ThisLoc = ' -> at HistContainer_Create (in History/histcontainer_mod.F90)'

    ! Allocate the Container object
    ALLOCATE( Container, STAT=RC )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Cannot allocate the "Container" object!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    ! Local value for ID
    ThisId = 0

    !========================================================================
    ! Required inputs, handle these first
    !========================================================================

    !---------------------------------
    ! Container ID
    !---------------------------------
    IF ( Id >= 0 ) THEN
       Container%Id = Id
    ELSE
       ErrMsg = 'History Container ID # cannot be negative!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !-----------------------
    ! Name
    !-----------------------
    IF ( LEN_TRIM( Name ) > 0 ) THEN
       TempStr   = Name
       Container%Name = TempStr
    ELSE
       ErrMsg = 'Must specify a name for this collection!!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !========================================================================
    ! Optional inputs, handle these next
    !========================================================================

    !----------------------------------
    ! EpochJd (Julian date @ start)
    !----------------------------------
    IF ( PRESENT( EpochJd ) ) THEN
       Container%EpochJd = EpochJd
    ELSE
       Container%EpochJd = UNDEFINED_DBL
    ENDIF

    !----------------------------------
    ! Current date in YYYMMDD
    !----------------------------------
    IF ( PRESENT( CurrentYmd ) ) THEN
       Container%CurrentYmd = CurrentYmd
    ELSE
       Container%CurrentYmd = 0
    ENDIF

    !----------------------------------
    ! Current time in hh:mm:ss
    !----------------------------------
    IF ( PRESENT( CurrentHms ) ) THEN
       Container%CurrentHms = CurrentHms
    ELSE
       Container%CurrentHms = 0
    ENDIF

    !----------------------------------
    ! Update mode
    !----------------------------------
    IF ( PRESENT( UpdateMode ) ) THEN
       Container%UpdateMode = UpdateMode
    ELSE
       Container%UpdateMode = ''
    ENDIF

    !----------------------------------
    ! Update frequency in YYYYMMDD
    !----------------------------------
    IF ( PRESENT( UpdateYmd ) ) THEN
       Container%UpdateYmd = UpdateYmd
    ELSE
       Container%UpdateYmd = 0
    ENDIF

    !----------------------------------
    ! Update frequency in hhmmss
    !----------------------------------
    IF ( PRESENT( UpdateHms ) ) THEN
       Container%UpdateHms = UpdateHms
    ELSE
       Container%UpdateHms = 0
    ENDIF

    !----------------------------------
    ! Update alarm (Julian date)
    !----------------------------------
    IF ( PRESENT( UpdateAlarm ) ) THEN
       Container%UpdateAlarm = UpdateAlarm
    ELSE
       Container%UpdateAlarm = UNDEFINED_DBL
    ENDIF

    !----------------------------------
    ! Operation code
    !----------------------------------
    IF ( PRESENT( Operation ) ) THEN
       Container%Operation = Operation
    ELSE
       Container%Operation = COPY_FROM_SOURCE
    ENDIF

    !----------------------------------
    ! Heartbeat timestep [min]
    !----------------------------------
    IF ( PRESENT( HeartBeatDtSec ) ) THEN
       Container%HeartBeatDtSec = HeartBeatDtSec
    ELSE
       Container%HeartBeatDtSec = UNDEFINED_DBL
    ENDIF

    !----------------------------------
    ! File write frequency in YYYYMMDD
    !----------------------------------
    IF ( PRESENT( FileWriteYmd ) ) THEN
       Container%FileWriteYmd = FileWriteYmd
    ELSE
       Container%FileWriteYmd = 0
    ENDIF

    !----------------------------------
    ! File write frequency in hhmmss
    !----------------------------------
    IF ( PRESENT( FileWriteHms ) ) THEN
       Container%FileWriteHms = FileWriteHms
    ELSE
       Container%FileWriteHms = 0
    ENDIF

    !----------------------------------
    ! File write alarm (Julian date)
    !----------------------------------
    IF ( PRESENT( FileWriteAlarm ) ) THEN
       Container%FileWriteAlarm = FileWriteAlarm
    ELSE
       Container%FileWriteAlarm = UNDEFINED_DBL
    ENDIF

    !----------------------------------
    ! File close frequency in YYYYMMDD
    !----------------------------------
    IF ( PRESENT( FileCloseYmd ) ) THEN
       Container%FileCloseYmd = FileCloseYmd
    ELSE
       Container%FileCloseYmd = 0
    ENDIF

    !----------------------------------
    ! File close frequency in hhmmss
    !----------------------------------
    IF ( PRESENT( FileCloseHms ) ) THEN
       Container%FileCloseHms = FileCloseHms
    ELSE
       Container%FileCloseHms = 0
    ENDIF

    !----------------------------------
    ! File close alarm (Julian date)
    !----------------------------------
    IF ( PRESENT( FileCloseAlarm ) ) THEN
       Container%FileCloseAlarm = FileCloseAlarm
    ELSE
       Container%FileCloseAlarm = UNDEFINED_DBL
    ENDIF

    !----------------------------------
    ! File ExpId (the dir name plus
    ! beginning of file name)
    !----------------------------------
    IF ( LEN_TRIM( FileExpId ) > 0 ) THEN
       TempStr             = FileExpId
       Container%FileExpId = TempStr
    ELSE
       Container%FileExpId = 'GEOSChem'
    ENDIF

    ! Add an error check.  The netCDF routines apparently cannot write
    ! files with "./" in the file path.  Strip out such occurrences.
    C = INDEX( Container%FileExpId, './' )
    IF ( C > 0 ) THEN
       Container%FileExpId = Container%FileExpId(C+2:)
    ENDIF

    !----------------------------------
    ! File Prefix
    !----------------------------------
    IF ( LEN_TRIM( FilePrefix ) > 0 ) THEN
       TempStr              = FilePrefix
       Container%FilePrefix = TempStr
    ELSE
       Container%FilePrefix = TRIM( Container%FileExpId ) // '.' //         &
                              TRIM( Name                ) // '.'
    ENDIF

    !----------------------------------
    ! File Template
    !----------------------------------
    IF ( LEN_TRIM( FileTemplate ) > 0 ) THEN

       ! If the FILETEMPLATE argument is passed (and not the undefined
       ! string) then use it.  Otherwise, construct a default template.
       IF ( TRIM( FileTemplate ) /= UNDEFINED_STR ) THEN
          TempStr                = FileTemplate
          Container%FileTemplate = TempStr
       ELSE
          Container%FileTemplate = '%y4%m2%d2_%h2%n2z.nc4'
       ENDIF

    ELSE

       ! If the FILETEMPLATE argument isn't passed,
       ! then construct a default template
       Container%FileTemplate = '%y4%m2%d2_%h2%n2z.nc4'

    ENDIF

    !----------------------------------
    ! File Name
    !----------------------------------
    IF ( LEN_TRIM( FileName ) > 0 ) THEN

       ! If the FILENAME argument is passed, then use it,
       ! otherwise, construct a default file name
       IF ( TRIM( FileName ) /= UNDEFINED_STR ) THEN
          TempStr                = FileName
          Container%FileName     = TempStr
          Container%FilePrefix   = UNDEFINED_STR
          Container%FileTemplate = UNDEFINED_STR
       ELSE
          Container%FileName = TRIM( Container%FilePrefix   ) // &
                               TRIM( Container%FileTemplate )
       ENDIF

    ELSE

       ! If the FILENAME argument isn't passed,
       ! construct a default file name
       Container%FileName = TRIM( Container%FilePrefix   ) // &
                            TRIM( Container%FileTemplate )
    ENDIF

    !----------------------------------
    ! Conventions
    !----------------------------------
    IF ( PRESENT( Conventions ) ) THEN
       TempStr               = Conventions
       Container%Conventions = TempStr
    ELSE
       Container%Conventions = ''
    ENDIF

    !----------------------------------
    ! NcFormat
    !----------------------------------
    IF ( PRESENT( NcFormat ) ) THEN
       TempStr            = NcFormat
       Container%NcFormat = TempStr
    ELSE
       Container%NcFormat = ''
    ENDIF

#if !defined( ESMF_ ) && !defined( NC_HAS_COMPRESSION )

    ! For GEOS-Chem Classic simulations compiled with either DEBUG=y or
    ! NC_NODEFLATE=y, set NcFormat to "NetCDF-3 with large file support",
    ! in order to denote that compression and chunking are disabled.
    Container%NcFormat = 'NetCDF-3 with large file support'

#endif

    !----------------------------------
    ! History
    !----------------------------------
    IF ( PRESENT( History ) ) THEN
       TempStr           = History
       Container%History = TempStr
    ELSE
       Container%History = ''
    ENDIF

    !----------------------------------
    ! ProdDateTime
    !----------------------------------
    IF ( PRESENT( ProdDateTime ) ) THEN
       TempStr                = ProdDateTime
       Container%ProdDateTime = TempStr
    ELSE
       Container%ProdDateTime = ''
    ENDIF

    !----------------------------------
    ! Reference
    !----------------------------------
    IF ( PRESENT( Reference ) ) THEN
       TempStr             = Reference
       Container%Reference = TempStr
    ELSE
       Container%Reference = ''
    ENDIF

    !----------------------------------
    ! Title
    !----------------------------------
    IF ( PRESENT( Title ) ) THEN
       TempStr         = Title
       Container%Title = TempStr
    ELSE
       Container%Title = ''
    ENDIF

    !----------------------------------
    ! Contact
    !----------------------------------
    IF ( PRESENT( Contact ) ) THEN
       TempStr           = Contact
       Container%Contact = TempStr
    ELSE
       Container%Contact = ''
    ENDIF

    !----------------------------------
    ! StartTimeStamp
    !----------------------------------
    IF ( PRESENT( StartTimeStamp ) ) THEN
       TempStr                  = StartTimeStamp
       Container%StartTimeStamp = TempStr
    ELSE
       Container%StartTimeStamp = ''
    ENDIF

    !----------------------------------
    ! EndTimeStamp
    !----------------------------------
    IF ( PRESENT( EndTimeStamp ) ) THEN
       TempStr                = EndTimeStamp
       Container%EndTimeStamp = TempStr
    ELSE
       Container%EndTimeStamp = ''
    ENDIF

    !=======================================================================
    ! Set other fields to initial or undefined values
    !=======================================================================

    ! These fields won't get defined until we open/write the netCDF file
    Container%IsFileDefined   = .FALSE.
    Container%IsFileOpen      = .FALSE.
    Container%FileId          = UNDEFINED_INT
    Container%xDimId          = UNDEFINED_INT
    Container%yDimId          = UNDEFINED_INT
    Container%zDimId          = UNDEFINED_INT
    Container%iDimId          = UNDEFINED_INT
    Container%tDimId          = UNDEFINED_INT
    Container%Spc_Units       = ''

    ! Set the other time/date fields from EpochJd, CurrentYmd, CurrentHms, etc.
    Container%EpochJsec       = Container%EpochJd * SECONDS_PER_DAY
    Container%CurrentJsec     = Container%EpochJSec
    Container%CurrentJd       = Container%EpochJd
    Container%ReferenceJsec   = Container%EpochJsec
    Container%ReferenceJd     = Container%EpochJd
    Container%ReferenceYmd    = Container%CurrentYmd
    Container%ReferenceHms    = Container%CurrentHms

    ! These other time fields will be defined later
    Container%ElapsedSec      = 0.0_f8
    Container%CurrTimeSlice   = UNDEFINED_INT
    Container%TimeStamp       = 0.0_f8

    ! Spatial information fields will be defined according to the
    ! dimensions of the HISTORY ITEMS belonging to the collection
    Container%NX              = UNDEFINED_INT
    Container%NY              = UNDEFINED_INT
    Container%NZ              = UNDEFINED_INT
    Container%X0              = UNDEFINED_INT
    Container%X1              = UNDEFINED_INT
    Container%Y0              = UNDEFINED_INT
    Container%Y1              = UNDEFINED_INT
    Container%Z0              = UNDEFINED_INT
    Container%OnLevelEdges    = .FALSE.

    ! If the collection is instantaneous, then set a flag to denote that
    ! first the netCDF file reference date/time should be the start-of-the-
    ! simulation time.  This will ensure that all timestamps and filenames
    ! for instantaneous collections are consistent.
    IF ( Container%Operation == COPY_FROM_SOURCE ) THEN
       Container%FirstInst    = .TRUE.
    ELSE
       Container%FirstInst    = .FALSE.
    ENDIF

    !=======================================================================
    ! Initialize the alarms (elapsed seconds since start of run)
    !=======================================================================

    !----------------------------------
    ! Initial UpdateAlarm interval
    !----------------------------------
    CALL HistContainer_UpdateIvalSet( Input_Opt, Container, RC )

    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered in "HistContainer_UpdateIvalSet"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !----------------------------------
    ! Initial FileCloseAlarm interval
    !----------------------------------
    CALL HistContainer_FileCloseIvalSet( Input_Opt, Container, RC )

    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered in "HistContainer_FileCloseIvalSet"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !----------------------------------
    ! Initial FileWriteAlarm interval
    !----------------------------------
    CALL HistContainer_FileWriteIvalSet( Input_Opt, Container, RC )

    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered in "HistContainer_FileWriteIvalSet"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !--------------------------------------
    ! Initial "UpdateAlarm" setting
    !--------------------------------------

    !-------------------------------------------------------------------------
    ! Prior to 3/5/19:
    ! Subtract the "heartbeat" timestep in seconds from UpdateAlarm.
    ! This will ensure that (1) Instantaneous file collections will be
    ! updated just before the file write, (2) Time-averaged collections
    ! will be averaged on the same timestep as the "historical" GEOS-Chem
    ! diagnostics, thus allowing for a direct comparison.
    !Container%UpdateAlarm = Container%UpdateIvalSec - Container%HeartBeatDtSec
    !-------------------------------------------------------------------------

    ! Set the initial UpdateAlarm value to the update interval in seconds
    ! NOTE: We no longer have to subtract the heartbeat timestep, because
    ! in the main program, we now call History_SetTime to advance the
    ! clock before calling History_Update to update the diagnostics.
    ! This will now allow us to recompute monthly or yearly intervals
    ! that span leap year days properly. (bmy, 3/5/19)
    Container%UpdateAlarm = Container%UpdateIvalSec

    ! Trap error if negative
    IF ( Container%UpdateAlarm < 0 ) THEN
       ErrMsg = 'UpdateAlarm for collection ' //                            &
                TRIM( Container%Name )        // ' is negative!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !--------------------------------------
    ! Initial "FileWriteAlarm" setting
    !--------------------------------------

    ! Set the file write alarm to its computed interval
    Container%FileWriteAlarm = Container%FileWriteIvalSec

    ! Trap error if negative
    IF ( Container%FileWriteAlarm < 0 ) THEN
       ErrMsg = 'FileWriteAlarm for collection ' //                         &
            TRIM( Container%Name )               // ' is negative!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !--------------------------------------
    ! Initial "FileCloseAlarm" setting
    !--------------------------------------

    IF ( Container%Operation == COPY_FROM_SOURCE ) THEN

       ! %%% INSTANTANEOUS %%%
       ! Create a new file ASAP so that we can start writing data to it
       Container%FileCloseAlarm = 0.0_f8

    ELSE

       ! %%% TIME-AVERAGED %%%
       ! Set the initial file close/reopen time to the first write time.
       ! (We will subtract this off later, when computing the reference
       ! date and time for the netCDF file.)
       Container%FileCloseAlarm = Container%FileWriteIvalSec

    ENDIF

    ! Trap error if negative
    IF ( Container%FileCloseAlarm < 0 ) THEN
       ErrMsg = 'FileCloseAlarm for collection ' //                         &
            TRIM( Container%Name )               // ' is negative!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

  END SUBROUTINE HistContainer_Create
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HistContainer_Print
!
! !DESCRIPTION: Prints information stored in a single HISTORY CONTAINER object.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HistContainer_Print( Input_Opt, Container, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod, ONLY : OptInput
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput),      INTENT(IN)  :: Input_Opt  ! Input Options object
    TYPE(HistContainer), POINTER     :: Container  ! HISTORY CONTAINER object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,             INTENT(OUT) :: RC         ! Success or failure
!
! !REVISION HISTORY:
!  16 Jun 2017 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOcAL VARIABLES:
!
    ! Strings
    CHARACTER(LEN=255)          :: DimStr

    ! String arrays
    CHARACTER(LEN=22)           :: OpCode(0:1) =                  &
                                     (/ 'Copy from source      ',  &
                                        'Accumulate from source' /)

    ! Objects
    TYPE(MetaHistItem), POINTER :: Current

    !=======================================================================
    ! Initialize
    !=======================================================================

    ! Assume success
    RC = GC_SUCCESS

    ! Free pointers
    Current => NULL()

    !=======================================================================
    ! Print information about this HISTORY CONTAINER
    ! only if we are on the root CPU
    !=======================================================================
    IF ( ASSOCIATED( Container ) .and. Input_Opt%amIRoot ) THEN
       WRITE( 6, 110 ) REPEAT( '-', 78 )
       WRITE( 6, 110 ) REPEAT( '-', 78 )
       WRITE( 6, 120 ) 'Container Name   : ', TRIM( Container%Name  )
       WRITE( 6, 130 ) 'Container Id #   : ', Container%Id
       WRITE( 6, 130 ) 'nX               : ', Container%nX
       WRITE( 6, 130 ) 'nY               : ', Container%nY
       WRITE( 6, 130 ) 'nZ               : ', Container%nZ
       WRITE( 6, 160 ) 'EpochJsec        : ', Container%EpochJsec
       WRITE( 6, 160 ) 'EpochJd          : ', Container%EpochJd
       WRITE( 6, 160 ) 'CurrentJSec      : ', Container%CurrentJSec
       WRITE( 6, 160 ) 'CurrentJd        : ', Container%CurrentJd
       WRITE( 6, 135 ) 'CurrentYmd       : ', Container%CurrentYmd
       WRITE( 6, 145 ) 'CurrentHms       : ', Container%CurrentHms
       WRITE( 6, 160 ) 'ElapsedSec       : ', Container%ElapsedSec
       WRITE( 6, 120 ) 'UpdateMode       : ', TRIM( Container%UpdateMode )
       WRITE( 6, 135 ) 'UpdateYmd        : ', Container%UpdateYmd
       WRITE( 6, 145 ) 'UpdateHms        : ', Container%UpdateHms
       WRITE( 6, 160 ) 'UpdateIvalSec    : ', Container%UpdateIvalSec
       WRITE( 6, 160 ) 'UpdateAlarm      : ', Container%UpdateAlarm
       WRITE( 6, 120 ) 'Operation        : ', OpCode( Container%Operation )
       WRITE( 6, 160 ) 'HeartBeatDtSec   : ', Container%HeartBeatDtSec
       WRITE( 6, 135 ) 'ReferenceYmd     : ', Container%ReferenceYmd
       WRITE( 6, 145 ) 'ReferenceHms     : ', Container%ReferenceHms
       WRITE( 6, 160 ) 'ReferenceJsec    : ', Container%ReferenceJd
       WRITE( 6, 160 ) 'ReferenceJd      : ', Container%ReferenceJd
       WRITE( 6, 135 ) 'FileWriteYmd     : ', Container%FileWriteYmd
       WRITE( 6, 145 ) 'FileWriteHms     : ', Container%FileWriteHms
       WRITE( 6, 160 ) 'FileWriteIvalSec : ', Container%FileWriteIvalSec
       WRITE( 6, 160 ) 'FileWriteAlarm   : ', Container%FileWriteAlarm
       WRITE( 6, 135 ) 'FileCloseYmd     : ', Container%FileCloseYmd
       WRITE( 6, 145 ) 'FileCloseHms     : ', Container%FileCloseHms
       WRITE( 6, 160 ) 'FileCloseIvalSec : ', Container%FileCloseIvalSec
       WRITE( 6, 160 ) 'FileCloseAlarm   : ', Container%FileCloseAlarm
       WRITE( 6, 130 ) 'CurrTimeSlice    : ', Container%CurrTimeSlice
       WRITE( 6, 150 ) 'IsFileOpen       : ', Container%IsFileOpen
       WRITE( 6, 150 ) 'IsFileDefined    : ', Container%IsFileDefined
       WRITE( 6, 130 ) 'FileId           : ', Container%FileId
       WRITE( 6, 130 ) 'xDimId           : ', Container%xDimId
       WRITE( 6, 130 ) 'yDimId           : ', Container%yDimId
       WRITE( 6, 130 ) 'zDimId           : ', Container%zDimId
       WRITE( 6, 130 ) 'tDimId           : ', Container%tDimId
       WRITE( 6, 120 ) 'FileExpId        : ', TRIM( Container%FileExpId    )
       WRITE( 6, 120 ) 'FilePrefix       : ', TRIM( Container%FilePrefix   )
       WRITE( 6, 120 ) 'FileTemplate     : ', TRIM( Container%FileTemplate )
       WRITE( 6, 120 ) 'Filename         : ', TRIM( Container%FileName     )
       WRITE( 6, 120 ) 'Conventions      : ', TRIM( Container%Conventions  )
       WRITE( 6, 120 ) 'NcFormat         : ', TRIM( Container%NcFormat     )
       WRITE( 6, 120 ) 'History          : ', TRIM( Container%History      )
       WRITE( 6, 120 ) 'ProdDateTime     : ', TRIM( Container%ProdDateTime )
       WRITE( 6, 120 ) 'Reference        : ', TRIM( Container%Reference    )
       WRITE( 6, 120 ) 'Title            : ', TRIM( Container%Title        )
       WRITE( 6, 120 ) 'Contact          : ', TRIM( Container%Contact      )
       WRITE( 6, 120 ) 'StartTimeStamp   : ', Container%StartTimeStamp
       WRITE( 6, 120 ) 'EndTimeStamp     : ', Container%EndTimeStamp
       WRITE( 6, 110 ) ''
       WRITE( 6, 110 ) 'Items archived in this collection:'

       ! FORMAT statements
 110   FORMAT( 1x, a           )
 120   FORMAT( 1x, a, a        )
 130   FORMAT( 1x, a, 7x, i8   )
 135   FORMAT( 1x, a, 7x, i8.8 )
 140   FORMAT( 1x, a, i6       )
 145   FORMAT( 1x, a, 9x, i6.6 )
 150   FORMAT( 1x, a, L15      )
 160   FORMAT( 1x, a, f17.1    )

       ! If there are HISTORY ITEMS belonging to this container ...
       IF ( ASSOCIATED( Container%HistItems ) ) THEN

          ! Point to the start of the list of HISTORY ITEMS
          Current => Container%HistItems

          ! As long as this HISTORY ITEM is valid ...
          DO WHILE ( ASSOCIATED( Current ) )

             ! Print the name, long-name, and units of each HISTORY ITEM
             ! that is stored in the METAHISTORY ITEM belonging to this
             ! HISTORY CONTAINER.  In other words, these are the diagnostic
             ! quantities that will get archived to the netCDF file.
             WRITE( 6, 100 ) Current%Item%Name,        &
                             Current%Item%LongName,    &
                             Current%Item%DimNames,    &
                             TRIM( Current%Item%Units )
 100         FORMAT( 2x, a20, ' | ', a35, ' | ', a3, ' | ', a )

             ! Skip to net item
             Current => Current%Next
          ENDDO

          ! Free pointers
          Current => NULL()
       ENDIF
    ENDIF

  END SUBROUTINE HistContainer_Print
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HistContainer_Destroy
!
! !DESCRIPTION: This method will destroy the METAHISTORY ITEM belonging to
!  a HISTORY CONTAINER.  It will then destroy the HISTORY CONTAINER itself.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HistContainer_Destroy( Container, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE MetaHistItem_Mod, ONLY : MetaHistItem_Destroy
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(HistContainer), POINTER     :: Container  ! HISTORY CONTAINER object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,             INTENT(OUT) :: RC         ! Success or failure
!
! !REVISION HISTORY:
!  16 Jun 2017 - R. Yantosca - Initial version, based on code by Arjen Markus
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
    ThisLoc = ' -> at HistContainer_Destroy (in History/histcontainer_mod.F90)'

    !=======================================================================
    ! Destroy the METAHISTORY ITEM belonging to this HISTORY CONTAINER
    !=======================================================================
    IF ( ASSOCIATED( Container%HistItems ) ) THEN
       CALL MetaHistItem_Destroy( Container%HistItems, RC )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'MetaHistItem_Destroy returned with error!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDIF

    !=======================================================================
    ! Then destroy the HISTORY CONTAINER itself
    !=======================================================================
    IF ( ASSOCIATED( Container ) ) THEN
       DEALLOCATE( Container, STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Could not destroy the "Container" HISTORY CONTAINER!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDIF

  END SUBROUTINE HistContainer_Destroy
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HistContainer_UpdateIvalSet
!
! !DESCRIPTION: Defines the alarm interval for the UPDATE operation.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HistContainer_UpdateIvalSet( Input_Opt, Container, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE History_Util_Mod
    USE Input_Opt_Mod,    ONLY : OptInput
    USE Time_Mod,         ONLY : Its_A_Leapyear, Ymd_Extract
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput),      INTENT(IN)  :: Input_Opt  ! Input Options object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(HistContainer), POINTER     :: Container  ! HISTORY CONTAINER object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,             INTENT(OUT) :: RC         ! Success or failure
!
!
! !REMARKS:
!  Assume that we will always update data more frequently than 1 month.
!  This means that we only have to compute this interval at initialization.
!
! !REVISION HISTORY:
!  06 Sep 2017 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER            :: Year,   Month,   Day
    INTEGER            :: Hour,   Minute,  Second

    ! Strings
    CHARACTER(LEN=255) :: ErrMsg, ThisLoc
    logical :: debug
    !=======================================================================
    ! Initialize
    !=======================================================================
    RC      = GC_SUCCESS
    ErrMsg  = ''
    ThisLoc = &
     ' -> at HistContainer_UpdateIvalSet (in History/histcontainer_mod.F90)'

    !=======================================================================
    ! Compute the interval for the "UpdateAlarm"
    !=======================================================================
    IF ( Container%UpdateYmd >= 010000 ) THEN

       !--------------------------------------------------------------------
       ! Update interval is one year or greater
       !--------------------------------------------------------------------

       ! Split the current date & time into its constituent values
       CALL Ymd_Extract( Container%CurrentYmd, Year, Month, Day )

       ! Update the alarm increment
       CALL AlarmIncrementYears( IntervalYmd = Container%UpdateYmd,          &
                                 Year        = Year,                         &
                                 Increment   = Container%UpdateIvalSec      )

    ELSE IF ( Container%UpdateYmd <  001200  .and.                           &
              Container%UpdateYmd >= 000100 ) THEN

       !--------------------------------------------------------------------
       ! Update interval is between 1 month and 1 year
       !--------------------------------------------------------------------

       ! Split the current date & time into its constituent values
       CALL Ymd_Extract( Container%CurrentYmd, Year, Month, Day )

       ! Update the alarm increment
       CALL AlarmIncrementMonths( IntervalYmd = Container%UpdateYmd,         &
                                  Year        = Year,                        &
                                  Month       = Month,                       &
                                  Increment   = Container%UpdateIvalSec     )

    ELSE

       !--------------------------------------------------------------------
       ! Update interval is less than 1 month
       !--------------------------------------------------------------------

       ! Split the file close interval date/time into its constituent values
       CALL Ymd_Extract( Container%UpdateYmd, Year, Month,  Day    )
       CALL Ymd_Extract( Container%UpdateHms, Hour, Minute, Second )

       ! "Update" interval in seconds
       Container%UpdateIvalSec = ( DBLE( Day    ) * SECONDS_PER_DAY    ) +   &
                                 ( DBLE( Hour   ) * SECONDS_PER_HOUR   ) +   &
                                 ( DBLE( Minute ) * SECONDS_PER_MINUTE ) +   &
                                 ( DBLE( Second )                      )
    ENDIF

  END SUBROUTINE HistContainer_UpdateIvalSet
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HistContainer_FileCloseIvalSet
!
! !DESCRIPTION: Defines the alarm interval for the FILE CLOSE/REOPEN operation.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HistContainer_FileCloseIvalSet( Input_Opt, Container, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE History_Util_Mod
    USE Input_Opt_Mod,    ONLY : OptInput
    USE Time_Mod,         ONLY : Its_A_Leapyear, Ymd_Extract
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput),      INTENT(IN)  :: Input_Opt  ! Input Options object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(HistContainer), POINTER     :: Container  ! HISTORY CONTAINER object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,             INTENT(OUT) :: RC         ! Success or failure
!
!
! !REMARKS:
!  The algorithm may not be as robust when straddling leap-year months, so we
!  would recommend selecting an interval of 1 month or 1 year at a time.
!
! !REVISION HISTORY:
!  06 Sep 2017 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER            :: Year,   Month,   Day
    INTEGER            :: Hour,   Minute,  Second

    ! Strings
    CHARACTER(LEN=255) :: ErrMsg, ThisLoc

    !=======================================================================
    ! Initialize
    !=======================================================================
    RC      = GC_SUCCESS
    ErrMsg  = ''
    ThisLoc = &
    ' -> at HistContainer_FileCloseIvalSet (in History/histcontainer_mod.F90)'

    !=======================================================================
    ! Compute the interval for the "FileCloseAlarm"
    !=======================================================================
    IF ( Container%FileCloseYmd >= 010000 ) THEN

       !--------------------------------------------------------------------
       ! File close interval is 1 year or greater
       !--------------------------------------------------------------------

       ! Split the current date & time into its constituent values
       CALL Ymd_Extract( Container%CurrentYmd, Year, Month, Day )

       ! Update the alarm increment
       CALL AlarmIncrementYears( IntervalYmd = Container%FileCloseYmd,       &
                                 Year        = Year,                         &
                                 Increment   = Container%FileCloseIvalSec   )

    ELSE IF ( Container%FileCloseYmd <  001200  .and.                        &
              Container%FileCloseYmd >= 000100 ) THEN

       !--------------------------------------------------------------------
       ! File close interval is between 1 month and 1 year
       !--------------------------------------------------------------------

       ! Split the current date & time into its constituent values
       CALL Ymd_Extract( Container%CurrentYmd, Year, Month, Day )

       ! Update the alarm increment
       CALL AlarmIncrementMonths( IntervalYmd = Container%FileCloseYmd,      &
                                  Year        = Year,                        &
                                  Month       = Month,                       &
                                  Increment   = Container%FileCloseIvalSec  )

    ELSE

       !--------------------------------------------------------------------
       ! File close interval is less than 1 month
       !--------------------------------------------------------------------

       ! Split the file close interval date/time into its constituent values
       CALL Ymd_Extract( Container%FileCloseYmd, Year, Month,  Day    )
       CALL Ymd_Extract( Container%FileCloseHms, Hour, Minute, Second )

       ! "FileClose" interval in seconds
       Container%FileCloseIvalSec = ( DBLE(Day   ) * SECONDS_PER_DAY    ) +  &
                                    ( DBLE(Hour  ) * SECONDS_PER_HOUR   ) +  &
                                    ( DBLE(Minute) * SECONDS_PER_MINUTE ) +  &
                                    ( DBLE(Second)                      )
    ENDIF

  END SUBROUTINE HistContainer_FileCloseIvalSet
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HistContainer_FileWriteIvalSet
!
! !DESCRIPTION: Defines the alarm intervals for the FILE WRITE operation.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HistContainer_FileWriteIvalSet( Input_Opt, Container, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE History_Util_Mod
    USE Input_Opt_Mod,    ONLY : OptInput
    USE Time_Mod,         ONLY : Its_A_Leapyear, Ymd_Extract
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput),      INTENT(IN)  :: Input_Opt  ! Input Options object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(HistContainer), POINTER     :: Container  ! HISTORY CONTAINER object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,             INTENT(OUT) :: RC         ! Success or failure
!
!
! !REMARKS:
!  The algorithm may not be as robust when straddling leap-year months, so we
!  would recommend selecting an interval of 1 month or 1 year at a time.
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
    INTEGER            :: Year,   Month,   Day
    INTEGER            :: Hour,   Minute,  Second

    ! Strings
    CHARACTER(LEN=255) :: ErrMsg, ThisLoc
    !=======================================================================
    ! Initialize
    !=======================================================================
    RC      = GC_SUCCESS
    ErrMsg  = ''
    ThisLoc = &
    ' -> at HistContainer_FileWriteIvalSet (in History/histcontainer_mod.F90)'

    !=======================================================================
    ! Compute the interval for the "FileWriteAlarm"
    !=======================================================================
    IF ( Container%FileWriteYmd >= 010000 ) THEN

       !--------------------------------------------------------------------
       ! File write interval is one year or greater
       !--------------------------------------------------------------------

       ! Split the current date & time into its constituent values
       CALL Ymd_Extract( Container%CurrentYmd, Year, Month, Day )

       ! Update the alarm increment
       CALL AlarmIncrementYears( IntervalYmd = Container%FileWriteYmd,       &
                                 Year        = Year,                         &
                                 Increment   = Container%FileWriteIvalSec   )

    ELSE IF ( Container%FileWriteYmd <  001200  .and.                        &
              Container%FileWriteYmd >= 000100 ) THEN

       !--------------------------------------------------------------------
       ! File write interval is one or more months but less than a year
       !
       ! This will probably be the most common option.
       ! Now accounts properly for leap year days. (bmy, 2/26/19)
       !--------------------------------------------------------------------

       ! Split the current date & time into its constituent values
       CALL Ymd_Extract( Container%CurrentYmd, Year, Month, Day )

       ! Update the alarm increment
       CALL AlarmIncrementMonths( IntervalYmd = Container%FileWriteYmd,      &
                                  Year        = Year,                        &
                                  Month       = Month,                       &
                                  Increment   = Container%FileWriteIvalSec  )

    ELSE

       !--------------------------------------------------------------------
       ! File write interval is less than a month
       !--------------------------------------------------------------------

       ! Split the file write interval date/time into its constituent values
       CALL Ymd_Extract( Container%FileWriteYmd, Year, Month,  Day    )
       CALL Ymd_Extract( Container%FileWriteHms, Hour, Minute, Second )

       ! "FileWrite" interval in seconds
       Container%FileWriteIvalSec = ( DBLE(Day   ) * SECONDS_PER_DAY    ) +  &
                                    ( DBLE(Hour  ) * SECONDS_PER_HOUR   ) +  &
                                    ( DBLE(Minute) * SECONDS_PER_MINUTE ) +  &
                                    ( DBLE(Second)                      )
    ENDIF

  END SUBROUTINE HistContainer_FileWriteIvalSet
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HistContainer_SetTime
!
! !DESCRIPTION: Increments the current astronomical Julian Date of a HISTORY
!  CONTAINER object by the HeartBeat interval (in fractional days).  Then it
!  recomputes the corresponding date/time and elapsed minutes.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HistContainer_SetTime( Input_Opt, Container, HeartBeatDt, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE History_Util_Mod
    USE Input_Opt_Mod,   ONLY : OptInput
    USE Julday_Mod,      ONLY :CALDATE
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput),      INTENT(IN)  :: Input_Opt   ! Input Options object
    TYPE(HistContainer), POINTER     :: Container   ! HISTORY CONTAINER object
    REAL(f8),            OPTIONAL    :: HeartBeatDt ! Heartbeat increment for
                                                    !  for timestepping [days]
!
! !OUTPUT PARAMETERS:
!
    INTEGER,             INTENT(OUT) :: RC          ! Success or failure
!
! !REMARKS:
!  This routine is called after the initial creation of the HISTORY
!  CONTAINER object.  It is also called from History_SetTime, which is
!  placed after the call to History_Update but before History_Write.
!
! !REVISION HISTORY:
!  21 Aug 2017 - R. Yantosca - Initial version
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
    RC      =  GC_SUCCESS
    ErrMsg  =  ''
    ThisLoc =  &
         ' -> at HistContainer_SetTime (in History/history_mod.F90)'

    !========================================================================
    ! Update the current time by the heartbeat time (in seconds)
    !========================================================================

    ! Update the Astronomical Julian seconds value by the heartbeat interval.
    ! Increment in seconds instead of days to avoid roundoff errors.
    IF ( PRESENT( HeartBeatDt ) ) THEN
       Container%CurrentJsec = Container%CurrentJsec +                       &
                               HeartBeatDt
    ELSE
       Container%CurrentJsec = Container%CurrentJsec +                       &
                               Container%HeartBeatDtSec
    ENDIF

    ! Convert Astronomical Julian Seconds to Astronomical Julian Date,
    ! for the conversion to calendar date and time. (bmy, 7/11/18)
    Container%CurrentJd = Container%CurrentJsec / SECONDS_PER_DAY

    ! Convert the Astronomical Julian Date to calendar date and time
    CALL CalDate( JulianDay = Container%CurrentJd,                           &
                  yyyymmdd  = Container%CurrentYmd,                          &
                  hhmmss    = Container%CurrentHms                          )

    !========================================================================
    ! Compute elapsed time quantities
    !========================================================================

   ! Compute the elapsed time in seconds since the start of the run
   CALL Compute_Elapsed_Time( CurrentJsec  = Container%CurrentJsec,          &
                              TimeBaseJsec = Container%EpochJsec,            &
                              ElapsedSec   = Container%ElapsedSec           )



  END SUBROUTINE HistContainer_SetTime
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ComputeAlarmIncrementYears
!
! !DESCRIPTION: Given an interval, computes the number of seconds to add
!  to an alarm, properly accounting for leap years.  This is for the case
!  when the update frequency is 1 year or greater.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE AlarmIncrementYears( IntervalYmd, Year, Increment )
!
! !USES:
!
    USE History_Util_Mod, ONLY : SECONDS_PER_DAY
    USE Time_Mod,         ONLY : Its_A_Leapyear
!
! !INPUT PARAMETERS:
!
    INTEGER,  INTENT(IN)  :: IntervalYmd  ! Update frequency in YYYYMMDD format
    INTEGER,  INTENT(IN)  :: Year         ! Current year
!
! !OUTPUT PARAMETERS:
!
    REAL(f8), INTENT(OUT) :: Increment    ! Number of seconds to add to alarm
!
! !REVISION HISTORY:
!  26 Feb 2019 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: nYears, T, YYYY

    !=======================================================================
    ! AlarmIncrementYears begins here!
    !=======================================================================

    ! Initialize
    Increment = 0.0_fp
    nYears    = IntervalYmd / 10000

    ! Loop over the requested # of years
    DO T = 0, nYears-1

       ! Increment the year from the starting year
       YYYY = Year + T

       ! Compute the increment, accounting for leap years
       IF ( Its_A_LeapYear( YYYY ) ) THEN
          Increment = Increment + ( 366.0_f8 * SECONDS_PER_DAY )
       ELSE
          Increment = Increment + ( 365.0_f8 * SECONDS_PER_DAY )
       ENDIF

    ENDDO

  END SUBROUTINE AlarmIncrementYears
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: AlarmIncremenmtMonths
!
! !DESCRIPTION: Given an interval, computes the number of seconds to add
!  to an alarm, properly accounting for leap years.  This is for the case
!  when the update frequency is between 1 month and 1 year.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE AlarmIncrementMonths( IntervalYmd, Year, Month, Increment )
!
! !USES:
!
    USE History_Util_Mod, ONLY : SECONDS_PER_DAY
    USE Time_Mod,         ONLY : Its_A_Leapyear
!
! !INPUT PARAMETERS:
!
    INTEGER,  INTENT(IN)  :: IntervalYmd  ! Update frequency in YYYYMMDD format
    INTEGER,  INTENT(IN)  :: Year         ! Current year
    INTEGER,  INTENT(IN)  :: Month        ! Current month
!
! !OUTPUT PARAMETERS:
!
    REAL(f8), INTENT(OUT) :: Increment    ! Number of seconds to add to alarm
!
! !REVISION HISTORY:
!  26 Feb 2019 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: MM, nDays, nMonths, T, YYYY
!
! !DEFINED PARAMETERS:
!
    ! Days in non-leap-year months:     J  F  M  A  M  J  J  A  S  O  N  D
    INTEGER, PARAMETER :: DpM(12) = (/ 31,28,31,30,31,30,31,31,30,31,30,31 /)

    !=======================================================================
    ! AlarmIncrememntMonths begins here!
    !=======================================================================

    ! Initialize
    Increment = 0.0_fp
    MM        = Month
    nDays     = 0
    nMonths   = IntervalYmd / 100
    YYYY      = Year

    ! Loop over the requested # of months
    DO T = 0, nMonths-1

       ! Keep a running total of the number of days in the interval
       nDays = nDays + DpM(MM)

       ! Add the leap year day if necessary
       IF ( Its_A_LeapYear( YYYY ) .and. MM == 2 ) THEN
          nDays = nDays + 1
       ENDIF

       ! Increment the month for next iteration
       MM = MM + 1

       ! Also increment the year if we straddle New Year's Day
       IF ( MM > 12 ) THEN
          MM   = 1
          YYYY = YYYY + 1
       ENDIF

    ENDDO

    ! Convert from days to seconds
    Increment = DBLE( nDays ) * SECONDS_PER_DAY

  END SUBROUTINE AlarmIncrementMonths
!EOC
END MODULE HistContainer_Mod

