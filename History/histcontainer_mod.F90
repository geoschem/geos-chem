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
  PUBLIC :: HistContainer_Create
  PUBLIC :: HistContainer_Print
  PUBLIC :: HistContainer_Destroy
  PUBLIC :: HistContainer_SetTime
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
                                                        !  1=accum from source
     INTEGER                     :: CurrentYmd          ! Current YMD date
     INTEGER                     :: CurrentHms          ! Current hms time
     REAL(f8)                    :: CurrentJd           ! Astronomical Julian
                                                        !  date @ current time
     REAL(f8)                    :: ElapsedMin          ! Elapsed minutes
                                                        !  since start of sim
     REAL(f8)                    :: UpdateAlarm         ! Alarm (elapsed min)
                                                        !  for data updating
     REAL(f8)                    :: FileCloseAlarm      ! Alarm (elapsed min)
                                                        !  for file close/open
     REAL(f8)                    :: FileWriteAlarm      ! Alarm (elapsed min)
                                                        !  for file write

     !----------------------------------------------------------------------
     ! Time quantities measured since the time of netCDF file creation          
     !----------------------------------------------------------------------
     INTEGER                     :: ReferenceYmd        ! Reference YMD & hms 
     INTEGER                     :: ReferenceHms        !  for the "time" dim
     REAL(f8)                    :: ReferenceJD         ! Julian Date at the
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
     REAL(f8)                    :: UpdateIvalMin       ! Update interval
                                                        !  in minutes
     INTEGER                     :: Operation           ! Operation code
                                                        !  0=copy from source

     !----------------------------------------------------------------------
     ! Quantities for file creation, writing, and I/O status
     !----------------------------------------------------------------------
     INTEGER                     :: FileWriteYmd        ! File write frequency
     INTEGER                     :: FileWriteHms        !  in YMD and hms
     REAL(f8)                    :: FileWriteIvalMin    ! File write interval
                                                        ! in minutes

     INTEGER                     :: FileCloseYmd        ! File closing time
     INTEGER                     :: FileCloseHms        !  in YMD and hms
     REAL(f8)                    :: FileCloseIvalMin    ! File close interval
                                                        !  in minutes

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
     CHARACTER(LEN=20)           :: Spc_Units           ! Units of SC%Species
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
!  07 Aug 2017 - R. Yantosca - Add FileWriteYmd, FileWriteHms
!  08 Aug 2017 - R. Yantosca - Add IsFileDefined, IsFileOpen, nX, nY, nZ and 
!                              the ouptuts xDimId, yIDimd, zDimId, tDimId
!  16 Aug 2017 - R. Yantosca - Rename Archival* variables to Update*
!  17 Aug 2017 - R. Yantosca - Added the *Alarm variables
!  18 Aug 2017 - R. Yantosca - Added EpochJd so that we can compute Julian
!                              dates as relative to the start of the run
!  18 Aug 2017 - R. Yantosca - Added ElapsedMin
!  18 Aug 2017 - R. Yantosca - Add HistContainer_ElapsedTime routine
!  21 Aug 2017 - R. Yantosca - Removed *_AlarmCheck, *_AlarmSet routines
!  24 Aug 2017 - R. Yantosca - Added iDimId as the dimension ID for ilev,
!                               which is the vertical dimension on interfaces
!  28 Aug 2017 - R. Yantosca - Added SpcUnits, FirstInst to type HistContainer
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
  SUBROUTINE HistContainer_Create( am_I_Root,      Container,                &
                                   Id,             Name,                     &
                                   RC,             EpochJd,                  &
                                   CurrentYmd,     CurrentHms,               &
                                   UpdateMode,     UpdateYmd,                &
                                   UpdateHms,      UpdateAlarm,              &
                                   Operation,      FileWriteYmd,             &
                                   FileWriteHms,   FileWriteAlarm,           &
                                   FileCloseYmd,   FileCloseHms,             &
                                   FileCloseAlarm, FileId,                   &
                                   FilePrefix,     FileName,                 &
                                   FileTemplate,   Conventions,              &
                                   NcFormat,       History,                  &
                                   ProdDateTime,   Reference,                &
                                   Title,          Contact                  )
!
! !USES:
!
  USE ErrCode_Mod
  USE History_Util_Mod
  USE MetaHistItem_Mod, ONLY : MetaHistItem
!
! !INPUT PARAMETERS: 
!
    !-----------------------------------------------------------------------
    ! REQUIRED INPUTS
    !-----------------------------------------------------------------------
    LOGICAL,             INTENT(IN)  :: am_I_Root      ! Root CPU?
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
!  07 Aug 2017 - R. Yantosca - Add FileWriteYmd and FileWriteHms
!  09 Aug 2017 - R. Yantosca - Add nX, ny, and, nZ
!  11 Aug 2017 - R. Yantosca - Add FileCloseYmd, FileCloseHms, ReferenceYmd,
!                              ReferenceHms, and CurrTimeSlice 
!  14 Aug 2017 - R. Yantosca - Add FileCloseYmd and FileCloseHms arguments
!  16 Aug 2017 - R. Yantosca - Renamed Archival* variables to Update*
!  17 Aug 2017 - R. Yantosca - Add *Alarm and Reference* arguments
!  18 Aug 2017 - R. Yantosca - Now initialize CurrentJd with EpochJd
!  21 Aug 2017 - R. Yantosca - Reorganize arguments, now define several time
!                              fields from EpochJd, CurrentYmd, CurrentHms
!  21 Aug 2017 - R. Yantosca - Now define initial alarm intervals and alarms
!  28 Aug 2017 - R. Yantosca - Now initialize Container%Spc_Units to null str
!  29 Aug 2017 - R. Yantosca - Reset NcFormat if netCDF compression is off
!                              for GEOS-Chem "Classic" simulations.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER            :: ThisId

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
    ! File Prefix
    !----------------------------------
    IF ( LEN_TRIM( FilePrefix ) > 0 ) THEN
       TempStr              = FilePrefix
       Container%FilePrefix = TempStr
    ELSE
       Container%FilePrefix = 'GEOSChem.History.' // TRIM( Name ) // '.'
    ENDIF

    !----------------------------------
    ! File Template
    !----------------------------------
    IF ( LEN_TRIM( FileTemplate) > 0 ) THEN
       TempStr                = FileTemplate
       Container%FileTemplate = TempStr
    ELSE
       Container%FileTemplate = ''
    ENDIF

    !----------------------------------
    ! File Name
    !----------------------------------
    IF ( LEN_TRIM( FileName ) > 0 ) THEN
       TempStr            = FileName
       Container%FileName = TempStr
    ELSE
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

    !=======================================================================
    ! Set other fields to initial or undefined values
    !=======================================================================
    
    ! These fields won't get defined until we open/write the netCDF file
    Container%IsFileDefined = .FALSE.
    Container%IsFileOpen    = .FALSE.
    Container%FileId        = UNDEFINED_INT
    Container%xDimId        = UNDEFINED_INT
    Container%yDimId        = UNDEFINED_INT
    Container%zDimId        = UNDEFINED_INT
    Container%iDimId        = UNDEFINED_INT
    Container%tDimId        = UNDEFINED_INT
    Container%Spc_Units     = ''

    ! Set the other time/date fields from EpochJd, CurrentYmd, CurrentHms
    Container%CurrentJd     = Container%EpochJd
    Container%ReferenceJd   = Container%EpochJd
    Container%ReferenceYmd  = Container%CurrentYmd
    Container%ReferenceHms  = Container%CurrentHms

    ! These other time fields will be defined later
    Container%ElapsedMin    = 0.0_f8
    Container%CurrTimeSlice = UNDEFINED_INT
    Container%TimeStamp     = 0.0_f8

    ! Spatial information fields will be defined according to the
    ! dimensions of the HISTORY TTEMS belonging to the collection
    Container%NX            = UNDEFINED_INT
    Container%NY            = UNDEFINED_INT
    Container%NZ            = UNDEFINED_INT
    Container%OnLevelEdges  = .FALSE.
    
    ! If the collection is instantaneous, then set a flag to denote that
    ! first the netCDF file reference date/time should be the start-of-the- 
    ! simulation time.  This will ensure that all timestamps and filenames
    ! for instantaneous collections are consistent.
    IF ( Container%Operation == COPY_FROM_SOURCE ) THEN
       Container%FirstInst  = .TRUE.
    ELSE
       Container%FirstInst  = .FALSE.
    ENDIF

    !=======================================================================
    ! Initialize the alarms
    !=======================================================================

    ! Get the initial alarm inervals (in elapsed min since start of run)
    CALL HistContainer_AlarmIntervalSet( am_I_Root, Container, RC )

    ! Trap potential error
    IF ( RC /= GC_SUCCESS ) THEN 
       ErrMsg = 'Error encountered in "HistContainer_AlarmIntervalSet"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    ! Set the initial update and file write alarms to the increment 
    ! values (measured in elapsed minutes since the start of the simulation)
    Container%UpdateAlarm    = Container%UpdateIvalMin
    Container%FileWriteAlarm = Container%FileWriteIvalMin

    IF ( Container%Operation == COPY_FROM_SOURCE ) THEN

       ! Open the file as soon as possible
       Container%FileCloseAlarm = 0.0_f8

    ELSE

       ! Set the initial file close/reopen to the first write time
       ! (We will subtract this off when computing the reference time)
       Container%FileCloseAlarm = Container%FileWriteIvalMin

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
  SUBROUTINE HistContainer_Print( am_I_Root, Container, RC )
!
! !USES:
! 
    USE ErrCode_Mod
!
! !INPUT PARAMETERS: 
!
    LOGICAL,             INTENT(IN)  :: am_I_Root  ! Are we on the root CPU?
    TYPE(HistContainer), POINTER     :: Container  ! HISTORY CONTAINER object
!
! !OUTPUT PARAMETERS: 
!
    INTEGER,             INTENT(OUT) :: RC         ! Success or failure 
!
! !REVISION HISTORY:
!  16 Jun 2017 - R. Yantosca - Initial version
!  16 Aug 2017 - R. Yantosca - Renamed Archival* variables to Update*
!  17 Aug 2017 - R. Yantosca - Now print *Alarm variables.  Also use the
!                              Item%DimNames field to print the data dims
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
    IF ( ASSOCIATED( Container ) .and. am_I_Root ) THEN
       WRITE( 6, 110 ) REPEAT( '-', 78 )
       WRITE( 6, 110 ) REPEAT( '-', 78 )
       WRITE( 6, 120 ) 'Container Name   : ', TRIM( Container%Name  )
       WRITE( 6, 130 ) 'Container Id #   : ', Container%Id
       WRITE( 6, 130 ) 'nX               : ', Container%nX
       WRITE( 6, 130 ) 'nY               : ', Container%nY
       WRITE( 6, 130 ) 'nZ               : ', Container%nZ
       WRITE( 6, 160 ) 'EpochJd          : ', Container%EpochJd
       WRITE( 6, 160 ) 'CurrentJd        : ', Container%CurrentJd
       WRITE( 6, 135 ) 'CurrentYmd       : ', Container%CurrentYmd
       WRITE( 6, 145 ) 'CurrentHms       : ', Container%CurrentHms
       WRITE( 6, 160 ) 'ElapsedMin       : ', Container%ElapsedMin
       WRITE( 6, 120 ) 'UpdateMode       : ', TRIM( Container%UpdateMode )
       WRITE( 6, 135 ) 'UpdateYmd        : ', Container%UpdateYmd
       WRITE( 6, 145 ) 'UpdateHms        : ', Container%UpdateHms
       WRITE( 6, 160 ) 'UpdateIvalMin    : ', Container%UpdateIvalMin
       WRITE( 6, 160 ) 'UpdateAlarm      : ', Container%UpdateAlarm
       WRITE( 6, 120 ) 'Operation        : ', OpCode( Container%Operation )
       WRITE( 6, 135 ) 'ReferenceYmd     : ', Container%ReferenceYmd
       WRITE( 6, 145 ) 'ReferenceHms     : ', Container%ReferenceHms
       WRITE( 6, 160 ) 'ReferenceJd      : ', Container%ReferenceJd
       WRITE( 6, 135 ) 'FileWriteYmd     : ', Container%FileWriteYmd
       WRITE( 6, 145 ) 'FileWriteHms     : ', Container%FileWriteHms
       WRITE( 6, 160 ) 'FileWriteIvalMin : ', Container%FileWriteIvalMin
       WRITE( 6, 160 ) 'FileWriteAlarm   : ', Container%FileWriteAlarm
       WRITE( 6, 135 ) 'FileCloseYmd     : ', Container%FileCloseYmd
       WRITE( 6, 145 ) 'FileCloseHms     : ', Container%FileCloseHms
       WRITE( 6, 160 ) 'FileCloseIvalMin : ', Container%FileCloseIvalMin
       WRITE( 6, 160 ) 'FileCloseAlarm   : ', Container%FileCloseAlarm
       WRITE( 6, 130 ) 'CurrTimeSlice    : ', Container%CurrTimeSlice
       WRITE( 6, 150 ) 'IsFileOpen       : ', Container%IsFileOpen
       WRITE( 6, 150 ) 'IsFileDefined    : ', Container%IsFileDefined
       WRITE( 6, 130 ) 'FileId           : ', Container%FileId
       WRITE( 6, 130 ) 'xDimId           : ', Container%xDimId 
       WRITE( 6, 130 ) 'yDimId           : ', Container%yDimId 
       WRITE( 6, 130 ) 'zDimId           : ', Container%zDimId 
       WRITE( 6, 130 ) 'tDimId           : ', Container%tDimId 
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
       WRITE( 6, 110 ) ''
       WRITE( 6, 110 ) 'Items archived in this collection:'

       ! FORMAT statements
 110   FORMAT( 1x, a        )
 120   FORMAT( 1x, a, a     )
 130   FORMAT( 1x, a, i8    )
 135   FORMAT( 1x, a, i8.8  )
 140   FORMAT( 1x, a, i6    )
 145   FORMAT( 1x, a, i6.6  )
 150   FORMAT( 1x, a, L8    )
 160   FORMAT( 1x, a, f13.4 )

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
  SUBROUTINE HistContainer_Destroy( am_I_Root, Container, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE MetaHistItem_Mod, ONLY : MetaHistItem_Destroy 
!
! !INPUT PARAMETERS:
!
    LOGICAL,             INTENT(IN)  :: am_I_Root  ! Are we on the root CPU?
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
       CALL MetaHistItem_Destroy( am_I_Root, Container%HistItems, RC )
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
! !IROUTINE: HistContainer_AlarmIntervalSet
!
! !DESCRIPTION: Defines the alarm intervals for the update, file write, and
!  file close/reopen operations.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HistContainer_AlarmIntervalSet( am_I_root, Container, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE History_Util_Mod
    USE Time_Mod,         ONLY : Ymd_Extract
!
! !INPUT PARAMETERS:
!
    LOGICAL,             INTENT(IN)  :: am_I_Root  ! Are we on the root CPU?
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
    INTEGER            :: Year, Month, Day, Hour, Minute, Second

    ! Strings
    CHARACTER(LEN=255) :: ErrMsg, ThisLoc

    !=======================================================================
    ! Initialize
    !=======================================================================
    RC      = GC_SUCCESS
    ErrMsg  = ''
    ThisLoc = ' -> at HistContainer_Destroy (in History/histcontainer_mod.F90)'

    !=======================================================================
    ! Compute the interval for the "UpdateAlarm" 
    !=======================================================================

    ! Split 
    CALL Ymd_Extract( Container%UpdateYmd, Year, Month,  Day    )
    CALL Ymd_Extract( Container%UpdateHms, Hour, Minute, Second )

    IF ( Container%UpdateYmd >= 000100 ) THEN

       ! Need to handle a way of updating the months

    ELSE
       
       ! "Update" interval in minutes (interval < 1 month)
       Container%UpdateIvalMin    = ( DBLE(Day   ) * MINUTES_PER_DAY    ) +  & 
                                    ( DBLE(Hour  ) * MINUTES_PER_HOUR   ) +  & 
                                    ( DBLE(Minute)                      ) +  &
                                    ( DBLE(Second) / SECONDS_PER_MINUTE )
    ENDIF
    
    !=======================================================================
    ! Compute the interval for the "FileWriteAlarm"
    !=======================================================================

    CALL Ymd_Extract( Container%FileWriteYmd, Year, Month,  Day    )
    CALL Ymd_Extract( Container%FileWriteHms, Hour, Minute, Second )

    IF ( Container%FileWriteYmd >= 000100 ) THEN
       
       ! Need to handle a way of doing the months

    ELSE
       

       ! "FileWrite" interval in minutes (interval < 1 month)
       Container%FileWriteIvalMin = ( DBLE(Day   ) * MINUTES_PER_DAY    ) +  &
                                    ( DBLE(Hour  ) * MINUTES_PER_HOUR   ) +  &
                                    ( DBLE(Minute)                      ) +  &
                                    ( DBLE(Second) / SECONDS_PER_MINUTE )
    ENDIF

    !=======================================================================
    ! Compute the interval for the "FileWriteAlarm"
    !=======================================================================

    CALL Ymd_Extract( Container%FileCloseYmd, Year, Month,  Day    )
    CALL Ymd_Extract( Container%FileCloseHms, Hour, Minute, Second )

    IF ( Container%FileCloseYmd >= 000100 ) THEN
       
    ELSE

       ! "FileClose" interval in minutes (interval < 1 month)
       Container%FileCloseIvalMin = ( DBLE(Day   ) * MINUTES_PER_DAY    ) +  &
                                    ( DBLE(Hour  ) * MINUTES_PER_HOUR   ) +  &
                                    ( DBLE(Minute)                      ) +  &
                                    ( DBLE(Second) / SECONDS_PER_MINUTE )
    ENDIF

  END SUBROUTINE HistContainer_AlarmIntervalSet
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
  SUBROUTINE HistContainer_SetTime( am_I_Root, Container, HeartBeatDt, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE History_Util_Mod
    USE Julday_Mod,      ONLY :CALDATE
!
! !INPUT PARAMETERS: 
!
    LOGICAL,             INTENT(IN)  :: am_I_Root   ! Are we on the root CPU?
    TYPE(HistContainer), POINTER     :: Container   ! HISTORY CONTAINER object
    REAL(f8),            INTENT(IN)  :: HeartBeatDt ! Heartbeat increment for
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
    ! Update the current Julian date by the heartbeat time (in days)
    !========================================================================

    ! Update the Julian date by the heart beat interval in decimal days
    Container%CurrentJd = Container%CurrentJd + HeartBeatDt

    ! Convert the Julian Day to year/month/day and hour/minutes/seconds
    CALL CalDate( JulianDay = Container%CurrentJd,                           &
                  yyyymmdd  = Container%CurrentYmd,                          &
                  hhmmss    = Container%CurrentHms                          )

    !========================================================================
    ! Compute elapsed time quantities
    !========================================================================

    ! Compute the elapsed time in minutes since the start of the run
    CALL Compute_Elapsed_Time( CurrentJd  = Container%CurrentJd,             &
                               TimeBaseJd = Container%EpochJd,               &
                               ElapsedMin = Container%ElapsedMin            )

  END SUBROUTINE HistContainer_SetTime
!EOC
END MODULE HistContainer_Mod

