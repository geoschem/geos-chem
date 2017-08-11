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
     ! Time-averaging information                     
     !----------------------------------------------------------------------
     CHARACTER(LEN=255)          :: ArchivalMode        ! Archival mode
     INTEGER                     :: ArchivalYmd         ! Archival freq
     INTEGER                     :: ArchivalHms         !  in YMD and hms
     INTEGER                     :: Operation           ! Operation code
                                                        !  0=copy from source
                                                        !  1=accum from source
     INTEGER                     :: ReferenceYmd        ! Reference YMD & hms 
     INTEGER                     :: ReferenceHms        !  for the "time" dim
     INTEGER                     :: CurrTimeSlice       ! Current time slice
                                                        !  for the "time" dim

     !----------------------------------------------------------------------
     ! File definition and I/O status information
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
  SUBROUTINE HistContainer_Create( am_I_Root,    Container,    Id,           &
                                   Name,         nX,           nY,           &
                                   nZ,           RC,           ArchivalMode, &
                                   ArchivalYmd,  ArchivalHms,  Operation,    &
                                   FileWriteYmd, FileWriteHms, FileId,       &
                                   FilePrefix,   FileName,     FileTemplate, &
                                   Conventions,  NcFormat,     History,      & 
                                   ProdDateTime, Reference,    Title,        &
                                   Contact                                   )
!
! !USES:
!
  USE ErrCode_Mod
  USE History_Params_Mod
  USE MetaHistItem_Mod, ONLY : MetaHistItem
!
! !INPUT PARAMETERS: 
!
    !-----------------------------------------------------------------------
    ! REQUIRED INPUTS
    !-----------------------------------------------------------------------
    LOGICAL,             INTENT(IN)  :: am_I_Root     ! Root CPU?
    INTEGER,             INTENT(IN)  :: Id            ! Container Id #
    CHARACTER(LEN=*),    INTENT(IN)  :: Name          ! Container name
    INTEGER,             INTENT(IN)  :: nX            ! X (or lon) dim size
    INTEGER,             INTENT(IN)  :: nY            ! Y (or lat) dim size
    INTEGER,             INTENT(IN)  :: nZ            ! Z (or lev) dim size
    
    !-----------------------------------------------------------------------
    ! OPTIONAL INPUTS: data archival
    !-----------------------------------------------------------------------
    CHARACTER(LEN=*),    OPTIONAL    :: ArchivalMode  ! Archival mode
    INTEGER,             OPTIONAL    :: ArchivalYmd   ! Archival frequency
    INTEGER,             OPTIONAL    :: ArchivalHms   !  in both YMD and hms
    INTEGER,             OPTIONAL    :: Operation     ! Operation code:
                                                      !  0=copy  from source
                                                      !  1=accum from source
    !-----------------------------------------------------------------------
    ! OPTIONAL INPUTS: File definition and I/O info
    !-----------------------------------------------------------------------
    INTEGER,             OPTIONAL    :: FileWriteYmd  ! File write frequency
    INTEGER,             OPTIONAL    :: FileWriteHms  !  in both YMD and hms

    !-----------------------------------------------------------------------
    ! OPTIONAL INPUTS: netCDF file identifiers and metadata
    !-----------------------------------------------------------------------
    INTEGER,             OPTIONAL    :: FileId        ! netCDF file ID
    CHARACTER(LEN=*),    OPTIONAL    :: FilePrefix    ! Filename prefix
    CHARACTER(LEN=*),    OPTIONAL    :: FileTemplate  ! YMDhms template
    CHARACTER(LEN=*),    OPTIONAL    :: Conventions   ! e.g. "COARDS"
    CHARACTER(LEN=*),    OPTIONAL    :: Filename      ! Name of nc file
    CHARACTER(LEN=*),    OPTIONAL    :: NcFormat      ! e.g. "netCDF-4"
    CHARACTER(LEN=*),    OPTIONAL    :: History       ! History
    CHARACTER(LEN=*),    OPTIONAL    :: ProdDateTime  ! When produced
    CHARACTER(LEN=*),    OPTIONAL    :: Reference     ! Reference string
    CHARACTER(LEN=*),    OPTIONAL    :: Title         ! Title string
    CHARACTER(LEN=*),    OPTIONAL    :: Contact       ! Contact string
!
! !INPUT/OUTPUT PARAMETERS: 
!
    TYPE(HistContainer), POINTER     :: Container     ! Collection object
!
! !OUTPUT PARAMETERS: 
!
    INTEGER,             INTENT(OUT) :: RC            ! Success or failure
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
    ! Archival mode
    !----------------------------------
    IF ( PRESENT( ArchivalMode ) ) THEN
       Container%ArchivalMode = ArchivalMode
    ELSE
       Container%ArchivalMode = ''
    ENDIF

    !----------------------------------
    ! Archival frequency in YY/MM/DD
    !----------------------------------
    IF ( PRESENT( ArchivalYmd ) ) THEN
       Container%ArchivalYmd = ArchivalYmd 
    ELSE
       Container%ArchivalYmd = 0
    ENDIF

    !----------------------------------
    ! Archival frequency in hh:mm:ss
    !----------------------------------
    IF ( PRESENT( ArchivalHms ) ) THEN
       Container%ArchivalHms = ArchivalHms 
    ELSE
       Container%ArchivalHms = 0
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
    ! File write frequency in YY/MM/DD
    !----------------------------------
    IF ( PRESENT( FileWriteYmd ) ) THEN
       Container%FileWriteYmd = FileWriteYmd 
    ELSE
       Container%FileWriteYmd = 0
    ENDIF

    !----------------------------------
    ! File write frequency in hh:mm:ss
    !----------------------------------
    IF ( PRESENT( FileWriteHms ) ) THEN
       Container%FileWriteHms = FileWriteHms 
    ELSE
       Container%FileWriteHms = 0
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
    ! These fields will not be set until we create the netCDF file,
    ! so for now set them to UNDEFINED values.
    !=======================================================================
    Container%IsFileDefined = .FALSE.
    Container%IsFileOpen    = .FALSE.
    Container%FileId        = UNDEFINED_INT
    Container%xDimId        = UNDEFINED_INT
    Container%yDimId        = UNDEFINED_INT
    Container%zDimId        = UNDEFINED_INT
    Container%tDimId        = UNDEFINED_INT
    Container%FileCloseYmd  = UNDEFINED_INT
    Container%FileCloseHms  = UNDEFINED_INT
    Container%ReferenceYmd  = UNDEFINED_INT
    Container%ReferenceHms  = UNDEFINED_INT
    Container%CurrTimeSlice = UNDEFINED_INT

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
       WRITE( 6, 120 ) 'ArchivalMode   : ', TRIM( Container%ArchivalMode )
       WRITE( 6, 135 ) 'ArchivalYmd    : ', Container%ArchivalYmd
       WRITE( 6, 145 ) 'ArchivalHms    : ', Container%ArchivalHms
       WRITE( 6, 120 ) 'Operation      : ', OpCode( Container%Operation )
       WRITE( 6, 150 ) 'IsFileDefined  : ', Container%IsFileDefined
       WRITE( 6, 150 ) 'IsFileOpen     : ', Container%IsFileOpen
       WRITE( 6, 135 ) 'ReferenceYmd   : ', Container%ReferenceYmd
       WRITE( 6, 145 ) 'ReferenceHms   : ', Container%ReferenceHms
       WRITE( 6, 135 ) 'FileWriteYmd   : ', Container%FileWriteYmd
       WRITE( 6, 145 ) 'FileWriteHms   : ', Container%FileWriteHms
       WRITE( 6, 135 ) 'FileCloseYmd   : ', Container%FileCloseYmd
       WRITE( 6, 145 ) 'FileCloseHms   : ', Container%FileCloseHms
       WRITE( 6, 130 ) 'CurrTimeSlice  : ', Container%CurrTimeSlice
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
 110   FORMAT( 1x, a       )
 120   FORMAT( 1x, a, a    )
 130   FORMAT( 1x, a, i8   )
 135   FORMAT( 1x, a, i8.8 )
 140   FORMAT( 1x, a, i6   )
 145   FORMAT( 1x, a, i6.6 )
 150   FORMAT( 1x, a, L8   )

       ! If there are HISTORY ITEMS belonging to this container ...
       IF ( ASSOCIATED( Container%HistItems ) ) THEN 
          
          ! Point to the start of the list of HISTORY ITEMS
          Current => Container%HistItems

          ! Print the name, long-name, and units of each HISTORY ITEM 
          ! that is stored in the METAHISTORY ITEM belonging to this 
          ! HISTORY CONTAINER.  In other words, these are the diagnostic
          ! quantities that will get archived to the netCDF file.
          DO WHILE ( ASSOCIATED( Current ) ) 

             DimStr = ''
             
             ! Pick a string for the # of dimensions
             SELECT CASE( Current%Item%SpaceDim )
                CASE( 3 )
                   DimStr = 'xyz'
                CASE( 2 )
                   Dimstr = 'xy'
                CASE( 1 )
                   Dimstr = 'x'
             END SELECT

             WRITE( 6, 100 ) Current%Item%Name,      &
                             Current%Item%LongName,  &
                             DimStr,                 &
                             TRIM( Current%Item%Units )
 100         FORMAT( 2x, a20, ' | ', a35, ' | ', a3, ' | ', a )

             ! Skip to net item
             Current => Current%Next
          ENDDO

          ! Free pointers
          Current => NULL()
       ENDIF
    ENDIF

    ! Format statements
 
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
END MODULE HistContainer_Mod

