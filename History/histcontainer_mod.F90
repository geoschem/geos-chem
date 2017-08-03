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
                                                        !  0=nothing
                                                        !  1=add
                                                        !  2=multiply
                                                      
     !----------------------------------------------------------------------
     ! netCDF file info and attributes                
     !----------------------------------------------------------------------
     INTEGER                     :: FileId              ! netCDF file ID
     CHARACTER(LEN=255)          :: FilePrefix          ! Filename prefix
     CHARACTER(LEN=255)          :: FileTemplate        ! YMDhms template
     CHARACTER(LEN=255)          :: FileName            ! Name of nc file
     CHARACTER(LEN=255)          :: Conventions         ! e.g. "COARDS"
     CHARACTER(LEN=255)          :: NcFormat            ! e.g. "netCDF-4"
     CHARACTER(LEN=255)          :: History             ! History
     CHARACTER(LEN=255)          :: ProdDateTime        ! When produced
     CHARACTER(LEN=255)          :: Reference           ! Reference string
     CHARACTER(LEN=255)          :: Title               ! Title string

  END TYPE HistContainer
!
! !REMARKS:
!  Linked list routines taken from original code (linkedlist.f90) 
!  by Arjen Markus; http://flibs.sourceforge.net/linked_list.html
!
! !REVISION HISTORY:
!  12 Jun 2017 - R. Yantosca - Initial version, based on history_list_mod.F90
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
                                   Name,         RC,           ArchivalMode, &
                                   ArchivalYmd,  ArchivalHms,  Operation,    &
                                   FileId,       FilePrefix,   FileName,     &
                                   FileTemplate, Conventions,  NcFormat,     &
                                   History,      ProdDateTime, Reference,    &
                                   Title                                    )
!
! !USES:
!
  USE ErrCode_Mod
  USE History_Params_Mod
  USE MetaHistItem_Mod, ONLY : MetaHistItem
!
! !INPUT PARAMETERS: 
!
    ! Required inputs
    LOGICAL,             INTENT(IN)  :: am_I_Root     ! Root CPU?
    INTEGER,             INTENT(IN)  :: Id            ! Container Id #
    CHARACTER(LEN=*),    INTENT(IN)  :: Name          ! Container name

    ! Optional inputs: data archival
    CHARACTER(LEN=*),    OPTIONAL    :: ArchivalMode  ! Archival mode
    INTEGER,             OPTIONAL    :: ArchivalYmd   ! Frequency in YY/MM/DD
    INTEGER,             OPTIONAL    :: ArchivalHms   ! Frequency in hh:mm:ss
    INTEGER,             OPTIONAL    :: Operation     ! Operation code
                                                      !  0=copy
                                                      !  1=add
                                                      !  2=multiply

    ! Optional inputs: filename and metadata
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

    !---------------------------------
    ! Archival mode
    !---------------------------------
    IF ( PRESENT( ArchivalMode ) ) THEN
       Container%ArchivalMode = ArchivalMode
    ELSE
       Container%ArchivalMode = ''
    ENDIF

    !---------------------------------
    ! Archival frequency in YY/MM/DD
    !---------------------------------
    IF ( PRESENT( ArchivalYmd ) ) THEN
       Container%ArchivalYmd = ArchivalYmd 
    ELSE
       Container%ArchivalYmd = 0
    ENDIF

    !---------------------------------
    ! Archival frequency in hh:mm:ss
    !---------------------------------
    IF ( PRESENT( ArchivalHms ) ) THEN
       Container%ArchivalHms = ArchivalHms 
    ELSE
       Container%ArchivalHms = 0
    ENDIF

    !---------------------------------
    ! Operation code
    !---------------------------------
    IF ( PRESENT( Operation ) ) THEN
       Container%Operation = Operation
    ELSE
       Container%Operation = COPY_SOURCE
    ENDIF

    !---------------------------------
    ! File Id
    !---------------------------------
    IF ( PRESENT( FileId ) ) THEN
       Container%FileId = FileId
    ELSE
       Container%FileId = -1
    ENDIF

    !---------------------------------
    ! File Prefix
    !---------------------------------
    IF ( LEN_TRIM( FilePrefix ) > 0 ) THEN
       TempStr              = FilePrefix
       Container%FilePrefix = TempStr
    ELSE
       Container%FilePrefix = 'GEOSChem.History.' // TRIM( Name ) // '.'
    ENDIF

    !---------------------------------
    ! File Template
    !---------------------------------
    IF ( LEN_TRIM( FileTemplate) > 0 ) THEN
       TempStr                = FileTemplate
       Container%FileTemplate = TempStr
    ELSE
       Container%FileTemplate = ''
    ENDIF

    !---------------------------------
    ! File Name
    !---------------------------------
    IF ( LEN_TRIM( FileName ) > 0 ) THEN
       TempStr            = FileName
       Container%FileName = TempStr
    ELSE
       Container%FileName = TRIM( Container%FilePrefix   ) // &
                            TRIM( Container%FileTemplate )
    ENDIF

    !---------------------------------
    ! Conventions
    !---------------------------------
    IF ( PRESENT( Conventions ) ) THEN
       TempStr               = Conventions
       Container%Conventions = TempStr
    ELSE
       Container%Conventions = ''
    ENDIF

    !---------------------------------
    ! NcFormat
    !---------------------------------
    IF ( PRESENT( NcFormat ) ) THEN
       TempStr            = NcFormat
       Container%NcFormat = TempStr
    ELSE
       Container%NcFormat = ''
    ENDIF

    !---------------------------------
    ! History
    !---------------------------------
    IF ( PRESENT( History ) ) THEN
       TempStr           = History
       Container%History = TempStr
    ELSE
       Container%History = ''
    ENDIF

    !---------------------------------
    ! ProdDateTime
    !---------------------------------
    IF ( PRESENT( ProdDateTime ) ) THEN
       TempStr                = ProdDateTime
       Container%ProdDateTime = TempStr
    ELSE
       Container%ProdDateTime = ''
    ENDIF

    !---------------------------------
    ! Reference
    !---------------------------------
    IF ( PRESENT( Reference ) ) THEN
       TempStr             = Reference
       Container%Reference = TempStr
    ELSE
       Container%Reference = ''
    ENDIF

    !---------------------------------
    ! Title
    !---------------------------------
    IF ( PRESENT( Title ) ) THEN
       TempStr         = Title
       Container%Title = TempStr
    ELSE
       Container%Title = ''
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
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOcAL VARIABLES:
!
    ! Strings
    CHARACTER(LEN=255)          :: DimStr

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
       PRINT*, REPEAT( '-', 78 )
       PRINT*, REPEAT( '-', 78 )
       PRINT*, 'Container Name : ', TRIM( Container%Name         )
       PRINT*, 'Container Id # : ', Container%Id
       PRINT*, 'ArchivalMode   : ', TRIM( Container%ArchivalMode )
       PRINT*, 'ArchivalYmd    : ', Container%ArchivalYmd
       PRINT*, 'ArchivalHms    : ', Container%ArchivalHms
       PRINT*, 'Operation code : ', Container%Operation
       PRINT*, 'FileId         : ', Container%FileId
       PRINT*, 'FilePrefix     : ', TRIM( Container%FilePrefix   )
       PRINT*, 'FileTemplate   : ', TRIM( Container%FileTemplate )
       PRINT*, 'Filename       : ', TRIM( Container%FileName     )
       PRINT*, 'Conventions    : ', TRIM( Container%Conventions  )
       PRINT*, 'NcFormat       : ', TRIM( Container%NcFormat     )
       PRINT*, 'History        : ', TRIM( Container%History      )
       PRINT*, 'ProdDateTime   : ', TRIM( Container%ProdDateTime )
       PRINT*, 'Reference      : ', TRIM( Container%Reference    )
       PRINT*, 'Title          : ', TRIM( Container%Title        )
       PRINT*
       PRINT*, 'Items archived in this collection:'

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
                CASE( 0 )
                   DimStr = '0'
             END SELECT

             WRITE( 6, 100 ) Current%Item%Name,      &
                             Current%Item%LongName,  &
                             DimStr,                 &
                             Current%Item%Units
  100        FORMAT( 2x, a20, ' | ', a35, ' | ', a3, ' | ', a15 )
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
END MODULE HistContainer_Mod

