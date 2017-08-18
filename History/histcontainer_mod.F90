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
  PUBLIC :: HistContainer_AlarmSet
  PUBLIC :: HistContainer_AlarmCheck
  PUBLIC :: HistContainer_ElapsedTime
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
     INTEGER                     :: Operation           ! Operation code
                                                        !  0=copy from source

     !----------------------------------------------------------------------
     ! Quantities for file creation, writing, and I/O status
     !----------------------------------------------------------------------
     INTEGER                     :: FileWriteYmd        ! File write frequency
     INTEGER                     :: FileWriteHms        !  in YMD and hms
     INTEGER                     :: FileCloseYmd        ! File closing time
     INTEGER                     :: FileCloseHms        !  in YMD and hms
     LOGICAL                     :: IsFileDefined       ! Have we done netCDF
                                                        !  define mode yet?
     LOGICAL                     :: IsFileOpen          ! Is the netCDF file
                                                        !  currently open?

     !----------------------------------------------------------------------
     ! netCDF file identifiers and attributes
     !----------------------------------------------------------------------
     INTEGER                     :: FileId              ! netCDF file ID
     INTEGER                     :: xDimId              ! X (or lon ) dim ID
     INTEGER                     :: yDimId              ! Y (or lat ) dim ID
     INTEGER                     :: zDimId              ! Z (or lev ) dim ID
     INTEGER                     :: tDimId              ! T (or time) dim ID
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
                                   nX,             nY,                       &
                                   nZ,             RC,                       &
                                   UpdateMode,     UpdateYmd,                &
                                   UpdateHms,      UpdateAlarm,              &
                                   Operation,      FileWriteYmd,             & 
                                   FileWriteHms,   FileWriteAlarm,           &
                                   FileCloseYmd,   FileCloseHms,             &
                                   FileCloseAlarm, EpochJd,                  &
                                   ReferenceYmd,   ReferenceHms,             &
                                   ReferenceJd,    FileId,                   &
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
    INTEGER,             INTENT(IN)  :: nX             ! X (or lon) dim size
    INTEGER,             INTENT(IN)  :: nY             ! Y (or lat) dim size
    INTEGER,             INTENT(IN)  :: nZ             ! Z (or lev) dim size

    !-----------------------------------------------------------------------
    ! OPTIONAL INPUTS: data update parameters
    !-----------------------------------------------------------------------
    CHARACTER(LEN=*),    OPTIONAL    :: UpdateMode     ! e.g. inst or time-avg
    INTEGER,             OPTIONAL    :: UpdateYmd      ! Update frequency
    INTEGER,             OPTIONAL    :: UpdateHms      !  in both YMD and hms
    REAL(f8),            OPTIONAL    :: UpdateAlarm    ! JD for data update
    INTEGER,             OPTIONAL    :: Operation      ! Operation code:
                                                       !  0=copy  from source
                                                       !  1=accum from source

    !-----------------------------------------------------------------------
    ! OPTIONAL INPUTS: File definition and I/O info
    !-----------------------------------------------------------------------
    REAL(f8),            OPTIONAL    :: EpochJd        ! Astronomical Julian
                                                       !  date @ start of sim
    INTEGER,             OPTIONAL    :: ReferenceYmd   ! Ref YMD for "time" var
    INTEGER,             OPTIONAL    :: ReferenceHms   ! Ref hms for "time" var
    REAL(f8),            OPTIONAL    :: ReferenceJd    ! Relative Julian date
                                                       !  @ ref date/time

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

    !-----------------------
    ! nX, nY, nZ
    !-----------------------
    Container%nX = nX
    Container%nY = nY
    Container%nZ = nZ

    !========================================================================
    ! Optional inputs, handle these next
    !========================================================================

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
    ! EpochJd (Julian date @ start)
    !----------------------------------
    IF ( PRESENT( EpochJd ) ) THEN
       Container%EpochJd = EpochJd
    ELSE
       Container%EpochJd = UNDEFINED_DBL
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
    ! Reference date in YYYMMDD
    !----------------------------------
    IF ( PRESENT( ReferenceYmd ) ) THEN
       Container%ReferenceYmd = ReferenceYmd
    ELSE
       Container%ReferenceYmd = 0
    ENDIF

    !----------------------------------
    ! ReferenceHms in hh:mm:ss
    !----------------------------------
    IF ( PRESENT( ReferenceHms ) ) THEN
       Container%ReferenceHms = ReferenceHms
    ELSE
       Container%ReferenceHms = 0
    ENDIF

    !----------------------------------
    ! Reference Julian Date
    !----------------------------------
    IF ( PRESENT( ReferenceJd ) ) THEN
       Container%ReferenceJd = ReferenceJd
    ELSE
       Container%ReferenceJd = UNDEFINED_DBL
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
    ! Set these fields to initial or undefined values
    ! Many of these won't be used until we open the netCDF file.
    !=======================================================================
    Container%IsFileDefined  = .FALSE.
    Container%IsFileOpen     = .FALSE.
    Container%FileId         = UNDEFINED_INT
    Container%xDimId         = UNDEFINED_INT
    Container%yDimId         = UNDEFINED_INT
    Container%zDimId         = UNDEFINED_INT
    Container%tDimId         = UNDEFINED_INT
    Container%CurrTimeSlice  = UNDEFINED_INT
    Container%CurrentJd      = Container%EpochJd
    Container%ElapsedMin     = 0.0_f8
    Container%TimeStamp      = 0.0_f8

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
       WRITE( 6, 120 ) 'Container Name : ', TRIM( Container%Name  )
       WRITE( 6, 130 ) 'Container Id # : ', Container%Id
       WRITE( 6, 130 ) 'nX             : ', Container%nX
       WRITE( 6, 130 ) 'nY             : ', Container%nY
       WRITE( 6, 130 ) 'nZ             : ', Container%nZ
       WRITE( 6, 160 ) 'EpochJd        : ', Container%EpochJd
       WRITE( 6, 160 ) 'CurrentJd      : ', Container%CurrentJd
       WRITE( 6, 160 ) 'ElapsedMin     : ', Container%ElapsedMin
       WRITE( 6, 120 ) 'UpdateMode     : ', TRIM( Container%UpdateMode )
       WRITE( 6, 135 ) 'UpdateYmd      : ', Container%UpdateYmd
       WRITE( 6, 145 ) 'UpdateHms      : ', Container%UpdateHms
       WRITE( 6, 160 ) 'UpdateAlarm    : ', Container%UpdateAlarm
       WRITE( 6, 120 ) 'Operation      : ', OpCode( Container%Operation )
       WRITE( 6, 135 ) 'ReferenceYmd   : ', Container%ReferenceYmd
       WRITE( 6, 145 ) 'ReferenceHms   : ', Container%ReferenceHms
       WRITE( 6, 160 ) 'ReferenceJd    : ', Container%ReferenceJd
       WRITE( 6, 135 ) 'FileWriteYmd   : ', Container%FileWriteYmd
       WRITE( 6, 145 ) 'FileWriteHms   : ', Container%FileWriteHms
       WRITE( 6, 160 ) 'FileWriteAlarm : ', Container%FileWriteAlarm
       WRITE( 6, 135 ) 'FileCloseYmd   : ', Container%FileCloseYmd
       WRITE( 6, 145 ) 'FileCloseHms   : ', Container%FileCloseHms
       WRITE( 6, 160 ) 'FileCloseAlarm : ', Container%FileCloseAlarm
       WRITE( 6, 130 ) 'CurrTimeSlice  : ', Container%CurrTimeSlice
       WRITE( 6, 150 ) 'IsFileOpen     : ', Container%IsFileOpen
       WRITE( 6, 150 ) 'IsFileDefined  : ', Container%IsFileDefined
       WRITE( 6, 130 ) 'FileId         : ', Container%FileId
       WRITE( 6, 130 ) 'xDimId         : ', Container%xDimId 
       WRITE( 6, 130 ) 'yDimId         : ', Container%yDimId 
       WRITE( 6, 130 ) 'zDimId         : ', Container%zDimId 
       WRITE( 6, 130 ) 'tDimId         : ', Container%tDimId 
       WRITE( 6, 120 ) 'FilePrefix     : ', TRIM( Container%FilePrefix   )
       WRITE( 6, 120 ) 'FileTemplate   : ', TRIM( Container%FileTemplate )
       WRITE( 6, 120 ) 'Filename       : ', TRIM( Container%FileName     )
       WRITE( 6, 120 ) 'Conventions    : ', TRIM( Container%Conventions  )
       WRITE( 6, 120 ) 'NcFormat       : ', TRIM( Container%NcFormat     )
       WRITE( 6, 120 ) 'History        : ', TRIM( Container%History      )
       WRITE( 6, 120 ) 'ProdDateTime   : ', TRIM( Container%ProdDateTime )
       WRITE( 6, 120 ) 'Reference      : ', TRIM( Container%Reference    )
       WRITE( 6, 120 ) 'Title          : ', TRIM( Container%Title        )
       WRITE( 6, 120 ) 'Contact        : ', TRIM( Container%Contact      )
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
! !IROUTINE: HistContainer_AlarmSet
!
! !DESCRIPTION: Given the current date and time, as well as in interval,
!  this routine computes the "alarm time" for data updating, file writing,
!  and file closing.  An "alarm time", which is represented in elapsed minutes
!  since the start of the simulation, represents the next date/time for which 
!  a given operation will occur.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HistContainer_AlarmSet( am_I_Root, Container, yyyymmdd,        &
                                     hhmmss,    AlarmType, RC              )
!
! !USES:
!
    USE ErrCode_Mod
    USE History_Util_Mod
    USE Roundoff_Mod,     ONLY : Roundoff
    USE Time_Mod,         ONLY : Its_A_LeapYear, Ymd_Extract
!
! !INPUT PARAMETERS: 
!
    LOGICAL,             INTENT(IN)  :: am_I_Root  ! Are we on the root CPU?
    TYPE(HistContainer), POINTER     :: Container  ! HISTORY CONTAINER object
    INTEGER,             INTENT(IN)  :: yyyymmdd   ! Current date in YMD
    INTEGER,             INTENT(IN)  :: hhmmss     ! Current time in hms
    INTEGER,             INTENT(IN)  :: AlarmType  ! 0=update
                                                   ! 1=file write
                                                   ! 2=file close/reopen
!
! !OUTPUT PARAMETERS: 
!
    INTEGER,             INTENT(OUT) :: RC         ! Success or failure
!
! !REMARKS:
!
! !REVISION HISTORY:
!  18 Aug 2017 - R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER            :: Year,      Month,      Day
    INTEGER            :: Hour,      Minute,     Second
    INTEGER            :: AlarmDate, AlarmTime
    REAL(f8)           :: AlarmJd,   AlarmMins

    ! Strings
    CHARACTER(LEN=255) :: ErrMsg,    ThisLoc

    !=======================================================================
    ! Initialize
    !=======================================================================
    RC        = GC_SUCCESS
    AlarmJd   = 0.0_f8
    AlarmMins = 0.0_f8
    ErrMsg    = ''
    ThisLoc   = &
       ' -> at HistContainer_AlarmSet (in History/histcontainer_mod.F90)'

    !=======================================================================
    ! Add the relevant time interval to the current date & time
    !=======================================================================
    IF ( AlarmType == ALARM_UPDATE ) THEN

       ! We are setting the UPDATE alarm
       ! Add the date & time interval to the current date & time
       AlarmDate = yyyymmdd + Container%UpdateYmd
       AlarmTime = hhmmss   + Container%UpdateHms

    ELSE IF ( AlarmType == ALARM_FILE_WRITE ) THEN

       ! We are setting the FILE WRITE alarm
       ! Add the date & time interval to the current date & time
       AlarmDate = yyyymmdd + Container%FileWriteYmd
       AlarmTime = hhmmss   + Container%FileWriteHms

    ELSE IF ( AlarmType == ALARM_FILE_CLOSE ) THEN

       ! We are setting the FILE CLOSE/REOPEN alarm
       ! Add the date & time interval to the current date & time
       AlarmDate = yyyymmdd + Container%FileCloseYmd
       AlarmTime = hhmmss   + Container%FileCloseHms

    ELSE

       ! Trap potential error
       ErrMsg = 'Invalid value for AlarmType!"'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN

    ENDIF

    !=======================================================================
    ! Compute the alarm time as elapsed minutes since start of simulation.
    !=======================================================================

    ! Split the current date/time into constituent fields
    CALL Ymd_Extract( AlarmDate, Year, Month,  Day    )
    CALL Ymd_Extract( AlarmTime, Hour, Minute, Second )

    ! Check seconds
    IF ( Second > 59 ) THEN
       Second = Second - 60
       Minute = Minute + 1
    ENDIF

    ! Check minutes
    IF ( Minute > 59 ) THEN
       Minute = Minute - 60
       Hour   = Hour   + 1
    ENDIF

    ! Check hours
    IF ( Hour > 23 ) THEN
       Hour   = Hour   - 24
       Day    = Day    + 1
    ENDIF

    ! Check days
    SELECT CASE ( Month )

       ! 30 days hath September, April, June, and November.
       CASE( 4, 6, 9, 11 )
          IF ( Day > 30 ) THEN
             Day   = Day   - 30
             Month = Month + 1
          ENDIF

       ! All the rest have 31.
       CASE( 1, 3, 5, 7, 8, 10, 12 )
          IF ( Day > 31 ) THEN
             Day   = Day   - 31
             Month = Month + 1
          ENDIF

       ! February has 28 alone.
       ! Except for leap years, that's the time
       ! when February's days are 29.
       CASE( 2 )
           IF ( Its_A_LeapYear( Year ) ) THEN
              IF ( Day > 29 ) THEN
                 Day   = Day - 29
                 Month = 3
              ENDIF
           ELSE
              IF ( Day > 28 ) THEN
                 Day   = Day - 28
                 Month = 3
              ENDIF
           ENDIF
        
    END SELECT

    ! Check months
    IF ( Month > 12 ) THEN
       Month = Month - 12
       Year  = Year  + 1
    ENDIF
           
    !=======================================================================
    ! Compute the alarm time and update the HISTORY CONTAINER object
    !=======================================================================

    ! Reconstruct the alarm date and alarm time, after having corrected
    ! the individual year, month, day, hour, minute, second values
    AlarmDate = ( Year * 10000 ) + ( Month  * 100 ) + Day
    AlarmTime = ( Hour * 10000 ) + ( Minute * 100 ) + Second
    
    ! Compute the elapsed minutes since the start of the 
    ! simulation, corresponding to AlarmDate and AlarmTime
    CALL HistContainer_ElapsedTime( am_I_Root  = am_I_Root,                  &
                                    Container  = Container,                  &
                                    yyyymmdd   = AlarmDate,                  &
                                    hhmmss     = AlarmTime,                  &
                                    TimeBase   = FROM_START_OF_SIM,          &
                                    ElapsedMin = AlarmMins,                  &
                                    RC         = RC                         )

    ! Trap potential error
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered in "HistContainer_ElapsedTime"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    ! Update the relevant alarm
    IF ( AlarmType == ALARM_UPDATE ) THEN

       ! Alarm for updating data values
       Container%UpdateAlarm = AlarmMins
       
    ELSE IF ( AlarmType == ALARM_FILE_WRITE ) THEN

       ! Alarm for writing to the netCDF file
       Container%FileWriteAlarm = AlarmMins

    ELSE IF ( AlarmType == ALARM_FILE_CLOSE ) THEN

       ! Alarm for closing the current netCDF file and 
       ! reopening it for the next diagnostic interval
       Container%FileCloseAlarm = AlarmMins

    ENDIF
    
  END SUBROUTINE HistContainer_AlarmSet
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HistContainer_AlarmCheck
!
! !DESCRIPTION: Checks if is time to: (1) Update the HISTORY ITEMS belonging
!  to the given HISTORY CONTAINER object; (2) Write the HISTORY ITEMS to the
!  netCDF file specified by this HISTORY CONTAINER object, or; (3) Close
!  the netCDF file and reopen it for the next diagnostic interval.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HistContainer_AlarmCheck( am_I_Root, Container, yyyymmdd,      &
                                       hhmmss,    RC,        DoUpdate,      &
                                       DoWrite,   DoClose                  )
!
! !USES:
!
    USE ErrCode_Mod
    USE History_Util_Mod
    USE Roundoff_Mod
!
! !INPUT PARAMETERS: 
!
    LOGICAL,             INTENT(IN)  :: am_I_Root  ! Are we on the root CPU?
    TYPE(HistContainer), POINTER     :: Container  ! HISTORY CONTAINER object
    INTEGER,             INTENT(IN)  :: yyyymmdd   ! Current date in YMD
    INTEGER,             INTENT(IN)  :: hhmmss     ! Current time in hms
!
! !OUTPUT PARAMETERS: 
!
    LOGICAL,             OPTIONAL    :: DoUpdate   ! Is it time for update?
    LOGICAL,             OPTIONAL    :: DoWrite    ! Is it time for file write? 
    LOGICAL,             OPTIONAL    :: DoClose    ! Is it time for file close?
    INTEGER,             INTENT(OUT) :: RC         ! Success or failure
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
    REAL(f8)           :: DiffMin, ElapsedMin

    ! Strings
    CHARACTER(LEN=255) :: ErrMsg,  ThisLoc

    !=======================================================================
    ! Initialize
    !=======================================================================
    RC         = GC_SUCCESS
    DiffMin    = 0.0_f8
    ElapsedMin = 0.0_f8
    ErrMsg     = ''
    ThisLoc    = &
       ' -> at HistContainer_ElapsedTime (in History/histcontainer_mod.F90)'

    !=======================================================================
    ! First, compute the elapsed time
    !=======================================================================

    ! Compute the elapsed time since the start of the simulation
    CALL HistContainer_ElapsedTime( am_I_Root  = am_I_Root,                  &
                                    Container  = Container,                  &
                                    yyyymmdd   = yyyymmdd,                   &
                                    hhmmss     = hhmmss,                     &
                                    TimeBase   = FROM_START_OF_SIM,          &
                                    ElapsedMin = Container%ElapsedMin,       &
                                    RC         = RC                         )

    ! Trap potential error
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered in "HistContainer_ElapsedTime"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
       
    !=======================================================================
    ! Epsilon-test the current time against the UpdateAlarm time to see
    ! if it is time to update the HISTORY ITEMS in this collection
    !=======================================================================
    IF ( PRESENT( DoUpdate ) ) THEN
       DiffMin  = ABS( Container%UpdateAlarm - Container%ElapsedMin )
       DoUpdate = ( DiffMin < EPSILON )  
!       print*, '@@@@@: ', TRIM( Container%Name ), DiffMin, DoUpdate
    ENDIF

    !=======================================================================
    ! Epsilon-test the current time against the FileWriteAlarm time to see
    ! if it is time to write the HISTORY ITEMS to the netCDF file
    !=======================================================================
    IF ( PRESENT( DoWrite ) ) THEN
       DiffMin  = ABS( Container%FileWriteAlarm - Container%ElapsedMin )
       DoWrite  = ( DiffMin < EPSILON )  
    ENDIF

    !=======================================================================
    ! Epsilon-test the current time against the FileCloseAlarm time to see
    ! if it is time to close this file and reopen it for the next interval
    !=======================================================================
    IF ( PRESENT( DoClose ) ) THEN
       DiffMin  = ABS( Container%FileCloseAlarm - Container%ElapsedMin )
       DoClose  = ( DiffMin < EPSILON )
    ENDIF

  END SUBROUTINE HistContainer_AlarmCheck
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HistContainer_ElapsedTime
!
! !DESCRIPTION: Computes elapsed time in minutes, with respect to either
!  the start of the simulation or the netCDF file reference date/time.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HistContainer_ElapsedTime( am_I_Root,  Container, TimeBase,      &
                                        ElapsedMin, RC,        yyyymmdd,      &
                                        hhmmss                               )
!
! !USES:
!
    USE ErrCode_Mod
    USE History_Util_Mod
    USE Roundoff_Mod
!
! !INPUT PARAMETERS: 
!
    LOGICAL,             INTENT(IN)  :: am_I_Root  ! Are we on the root CPU?
    TYPE(HistContainer), POINTER     :: Container  ! HISTORY CONTAINER object
    INTEGER,             OPTIONAL    :: yyyymmdd   ! Current date in YMD
    INTEGER,             OPTIONAL    :: hhmmss     ! Current time in hms
    INTEGER,             INTENT(IN)  :: TimeBase   ! 0=From start of simulation
                                                   ! 1=From netCDF ref date
!
! !OUTPUT PARAMETERS: 
!
    REAL(f8),            INTENT(OUT) :: ElapsedMin ! Time in elapsed minutes
    INTEGER,             INTENT(OUT) :: RC         ! Success or failure
!
! !REMARKS:
!  If yyyymmdd and hhmmss are passed, then routine Compute_Julian_Date
!  will be called to obtain the corresponding Astronomical Julian Date.
!  Otherwise, the Julian Date will be taken from the Container%CurrentJd
!  field.
!
!  The netCDF file reference date and time are given by the ReferenceYmd,
!  ReferenceHms, and ReferenceJd fields of the Container object.  This
!  denotes the simulation date & time when the netCDF file was created.
!
! !REVISION HISTORY:
!  18 Aug 2017 - R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    REAL(f8)           :: JulianDate

    ! Strings
    CHARACTER(LEN=255) :: ErrMsg, ThisLoc

    !=======================================================================
    ! Initialize
    !=======================================================================
    RC         = GC_SUCCESS
    ElapsedMin = 0.0_f8
    JulianDate = 0.0_f8
    ErrMsg     = ''
    ThisLoc    = &
       ' -> at HistContainer_ElapsedTime (in History/histcontainer_mod.F90)'

    !=======================================================================
    ! Get the Astronomical Julian Date for the given time
    !=======================================================================
    IF ( PRESENT( yyyymmdd ) .and. PRESENT( hhmmss ) ) THEN
       
       ! If yyyymmdd and hhhmms arguments are present, then
       ! compute the corresponding Astronomical Julian Date
       CALL Compute_Julian_Date( yyyymmdd, hhmmss, JulianDate )

    ELSE

       ! Otherwise, take the current Julian Date value 
       ! from this HISTORY CONTAINER object.
       JulianDate = Container%CurrentJd
          
    ENDIF

    !=======================================================================
    ! Compute the elapsed time
    !=======================================================================
    IF ( TimeBase == FROM_START_OF_SIM ) THEN
       
       ! Compute elapsed minutes since start of simulation
       ElapsedMin = ( JulianDate - Container%EpochJd ) * 1440.0_f8

       ! Round off to a few places to avoid numerical noise 
       ElapsedMin = Roundoff( ElapsedMin, ROUNDOFF_DECIMALS )
       
    ELSE IF ( TimeBase == FROM_FILE_CREATE ) THEN

       ! Reference time is the netCDF file reference time/date
       ! (which is the model date & time when the file was created)
       ElapsedMin = ( JulianDate - Container%ReferenceJd ) * 1440.0_f8

       ! Round off to a few places to avoid numerical noise 
       ElapsedMin = Roundoff( ElapsedMin, ROUNDOFF_DECIMALS )

    ELSE

       ! Trap potential error
       ErrMsg = 'Invalid value for TimeBase, must be 0 or 1!"'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN

    ENDIF

  END SUBROUTINE HistContainer_ElapsedTime
!EOC
END MODULE HistContainer_Mod

