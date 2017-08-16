!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: history_netcdf_mod.F90 
!
! !DESCRIPTION: Contains routines to create a netCDF file for each GEOS-Chem
!  diagnostic collection (as specified by each HISTORY CONTAINER in the 
!  master collection list located within in history_mod.F90).
!\\
!\\
! !INTERFACE:
!
MODULE History_Netcdf_Mod
!
! !USES:
!
  USE Precision_Mod
  USE MetaHistItem_Mod, ONLY : MetaHistItem
  
  IMPLICIT NONE
  PRIVATE

# include "netcdf.inc"
!
! !PUBLIC MEMBER FUNCTIONS
!
  PUBLIC  :: History_Netcdf_Close
  PUBLIC  :: History_Netcdf_Define
  PUBLIC  :: History_Netcdf_Write
  PUBLIC  :: History_Netcdf_Init
  PUBLIC  :: History_Netcdf_Cleanup
!
! !PRIVATE MEMBER FUNCTIONS:
!
  PRIVATE :: Compute_Julian_Date
  PRIVATE :: Compute_TimeAvgOffset
  PRIVATE :: Expand_Date_Time
  PRIVATE :: Get_Var_DimIds
  PRIVATE :: Set_RefDateTime
!
! !REMARKS:
!
! !REVISION HISTORY:
!  10 Aug 2017 - R. Yantosca - Initial version
!  16 Aug 2017 - R. Yantosca - Reorder placement of routines
!  16 Aug 2017 - R. Yantosca - Rename History_Expand_Date to Expand_Date_Time
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE TYPES:
!
  ! Linked list of HISTORY ITEMS for netCDF index varaibles (lon, lat, etc)
  TYPE(MetaHistItem), POINTER :: IndexVarList => NULL()

CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: History_Netcdf_Close
!
! !DESCRIPTION: Closes the netCDF file specified by the the given HISTORY 
!  CONTAINER object.  Also resets the relevant fields of the HISTORY CONTAINER 
!  object (as well as the fields in each HISTORY ITEM contained within the
!  HISTORY CONTAINER) to undefined values.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE History_Netcdf_Close( am_I_Root, Container, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE HistContainer_Mod,     ONLY : HistContainer
    USE History_Params_Mod
    USE MetaHistContainer_Mod, ONLY : MetaHistContainer
    USE Ncdf_Mod
!
! !INPUT PARAMETERS: 
!
    LOGICAL,             INTENT(IN)  :: am_I_Root   ! Are we on the root CPU?
!
! !INPUT/OUTPUT PARAMETERS: 
!
    TYPE(HistContainer), POINTER     :: Container   ! HISTORY CONTAINER obj
!
! !OUTPUT PARAMETERS: 
!
    INTEGER,             INTENT(OUT) :: RC          ! Success or failure?
!
! !REMARKS:
!  This had been part of History_Netcdf_Define, but is now its own routine.
!
! !REVISION HISTORY:
!  14 Aug 2017 - R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    TYPE(MetaHIstItem), POINTER :: Current

    !=======================================================================
    ! Initialize
    !=======================================================================
    RC      =  GC_SUCCESS
    Current => NULL()

    !=======================================================================
    ! Close the previous file in this collection if it is still open
    !=======================================================================
    IF ( Container%IsFileOpen .or. Container%IsFileDefined ) THEN

       !--------------------------------------------------------------------
       ! Close file and undefine fields of the HISTORY CONTAINER object
       !--------------------------------------------------------------------

       ! Close the netCDF file
       CALL Nc_Close( Container%FileId )

       ! Undefine fields
       Container%IsFileOpen    = .FALSE.
       Container%IsFileDefined = .FALSE.
       Container%ReferenceYmd  = 0
       Container%ReferenceHms  = 0
       Container%ReferenceJd   = 0.0_f8
       Container%CurrTimeSlice = UNDEFINED_INT

       !--------------------------------------------------------------------
       ! Undefine relevant fields of each HISTORY ITEM object
       ! belonging to this HISTORY CONTAINER object
       !--------------------------------------------------------------------

       ! Set CURRENT to the first entry in the list of 
       ! HISTORY ITEMS belonging to this collection
       Current => Container%HistItems

       ! As long as this node of the list is valid ...
       DO WHILE( ASSOCIATED( Current ) )

          ! Undefine quantities for the file we just closed
          Current%Item%NcXDimId  = UNDEFINED_INT
          Current%Item%NcYDimId  = UNDEFINED_INT
          Current%Item%NcZDimId  = UNDEFINED_INT
          Current%Item%NcTDimId  = UNDEFINED_INT
          Current%Item%NcVarId   = UNDEFINED_INT
          
          ! Go to the next entry in the list of HISTORY ITEMS
          Current => Current%Next
       ENDDO

       ! Free pointer
       Current => NULL()
    ENDIF

  END SUBROUTINE History_Netcdf_Close
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: History_Netcdf_Define
!
! !DESCRIPTION: Creates the netCDF file specified by each HISTORY CONTAINER
!  object, and defines the variables specified by the HISTORY ITEMS beloinging
!  to the HISTORY CONTAINER.  Index variables lon, lat, lev, time, as well
!  as the AREA variable, are written to the netCDF file with the proper
!  metadata.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE History_Netcdf_Define( am_I_Root, Container,                    &
                                    yyyymmdd,  hhmmss,    RC                )
!
! !USES:
!
    USE ErrCode_Mod
    USE HistContainer_Mod,  ONLY : HistContainer, HistContainer_Print
    USE HistItem_Mod,       ONLY : HistItem,      HistItem_Print
    USE History_Params_Mod
    USE MetaHistItem_Mod,   ONLY : MetaHistItem
    USE Ncdf_Mod
!
! !INPUT PARAMETERS: 
!
    LOGICAL,             INTENT(IN)  :: am_I_Root ! Are we on the root CPU?
    INTEGER,             INTENT(IN)  :: yyyymmdd  ! Current Year/month/day
    INTEGER,             INTENT(IN)  :: hhmmss    ! Current hour/minute/second
!
! !INPUT/OUTPUT PARAMETERS           ::
!
    TYPE(HistContainer), POINTER     :: Container ! Diagnostic collection obj
!
! !OUTPUT PARAMETERS: 
!
    INTEGER,             INTENT(OUT) :: RC        ! Success or failure
!
! !REMARKS:
!
! !REVISION HISTORY:
!  03 Aug 2017 - R. Yantosca - Initial version
!  14 Aug 2017 - R. Yantosca - Call History_Netcdf_Close from 1 level higher
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER                     :: N,         VarCt
    INTEGER                     :: VarXDimId, VarYDimId
    INTEGER                     :: VarZDimId, VarTDimId

    ! Strings                   
    CHARACTER(LEN=5)            :: Z
    CHARACTER(LEN=8)            :: D
    CHARACTER(LEN=10)           :: T
    CHARACTER(LEN=255)          :: FileName
    CHARACTER(LEN=255)          :: ErrMsg,    ThisLoc,     VarUnits
    CHARACTER(LEN=255)          :: VarAxis,   VarPositive, VarCalendar

    ! Arrays
    INTEGER                     :: V(8)

    ! Pointers
    REAL(fp),           POINTER :: Data1d(:)
    REAL(fp),           POINTER :: Data2d(:)

    ! Objects
    TYPE(MetaHistItem), POINTER :: Current

    !=======================================================================
    ! Initialize
    !=======================================================================
    RC          =  GC_SUCCESS
    Current     => NULL()
    ErrMsg      =  ''
    ThisLoc     =  ' -> at History_Netcdf_Define (in History/history_mod.F90)'
    FileName    =  ''
    VarAxis     =  ''
    VarCalendar =  ''
    VarPositive =  ''
    VarUnits    =  ''

    !=======================================================================
    ! Create the netCDF file with global attributes
    ! Do not exit netCDF define mode just yet
    !=======================================================================
    IF ( .not. Container%IsFileOpen ) THEN

       !--------------------------------------------------------------------
       ! Replace time and date tokens in the netCDF file name
       !--------------------------------------------------------------------

       ! Save the collection's file name in a temporary variable
       FileName = TRIM( Container%FileName )

       ! Replace date and time tokens in the file name
       CALL Expand_Date_Time( FileName, yyyymmdd, hhmmss, MAPL_Style=.TRUE. )
   
       !--------------------------------------------------------------------
       ! Compute reference date and time fields in the HISTORY CONTAINER
       ! These are needed to compute the time stamps for each data field
       ! that is written to the netCDF file.  Also resets the current
       ! time slice index.
       !--------------------------------------------------------------------
       CALL Set_RefDateTime( Container, yyyymmdd, hhmmss )

       !--------------------------------------------------------------------
       ! Create the timestamp for the History and ProdDateTime attributes
       !--------------------------------------------------------------------

       ! Call F90 intrinsic DATE_AND_TIME Function
       D = 'ccyymmdd'
       T = 'hhmmss.sss'
       CALL Date_And_Time( Date=D, Time=T, Zone=Z, Values=V )  ! GMT time
     
       ! Create timestamp strings
       WRITE( Container%History,      10 ) V(1),V(2),V(3),V(5),V(6),V(7),Z
       WRITE( Container%ProdDateTime, 10 ) V(1),V(2),V(3),V(5),V(6),V(7),Z
 10    FORMAT( 'Produced on ', i4.4, '/', i2.2, '/', i2.2, 1x,               &
                               i2.2, ':', i2.2, ':', i2.2, ' UTC', a        )

       !--------------------------------------------------------------------
       ! Create the file and add global attributes
       ! Remain in netCDF define mode upon exiting this routine
       !
       ! NOTE: Container%Reference is a global attribute lists the GEOS-Chem
       ! web and wiki page.  It has nothing to do with the reference date
       ! and time fields that are computed by History_Set_RefDateTime.
       !--------------------------------------------------------------------
       CALL Nc_Create( Create_Nc4   = .TRUE.,                                &
                       NcFile       = FileName,                              &
                       nLon         = Container%nX,                          &
                       nLat         = Container%nY,                          &
                       nLev         = Container%nZ,                          &
                       nTime        = NF_UNLIMITED,                          &
                       NcFormat     = Container%NcFormat,                    &
                       Conventions  = Container%Conventions,                 &
                       History      = Container%History,                     &
                       ProdDateTime = Container%ProdDateTime,                &
                       Reference    = Container%Reference,                   &
                       Title        = Container%Title,                       &
                       Contact      = Container%Contact,                     &
                       fId          = Container%FileId,                      &
                       TimeId       = Container%tDimId,                      &
                       LevId        = Container%zDimId,                      &
                       LatId        = Container%yDimId,                      &
                       LonId        = Container%xDimId,                      &
                       KeepDefMode  = .TRUE.,                                &
                       Varct        = VarCt                                 )
         
       !--------------------------------------------------------------------
       ! Denote that the file has been created and is open
       !--------------------------------------------------------------------
       Container%IsFileOpen = .TRUE.
    ENDIF

    !=======================================================================
    ! Define all of the Create the netCDF file with global attributes
    !=======================================================================
    IF ( .not. Container%IsFileDefined ) THEN

       !--------------------------------------------------------------------
       ! Add the index dimension data
       !--------------------------------------------------------------------

       ! Set CURRENT to the first node in the list of HISTORY ITEMS
       Current => IndexVarList

       ! As long as this node of the list is valid ...
       DO WHILE( ASSOCIATED( Current ) )

          ! Get the dimension ID's that are relevant to each HISTORY ITEM
          ! Also get Axis, Calendar, and Positive attributes for index vars,
          ! and replace time & date tokens in the units string for "time".
          CALL Get_Var_DimIds( Item        = Current%Item,                   & 
                               xDimId      = Container%xDimId,               &
                               yDimId      = Container%yDimId,               &
                               zDimId      = Container%zDimId,               &
                               tDimId      = Container%tDimId,               &
                               RefDate     = Container%ReferenceYmd,         &
                               RefTime     = Container%ReferenceHms,         &
                               VarAxis     = VarAxis,                        &
                               VarPositive = VarPositive,                    &
                               VarCalendar = VarCalendar,                    &
                               VarUnits    = VarUnits                       )

          ! Define each HISTORY ITEM in this collection to the netCDF file
          CALL Nc_Var_Def( DefMode      = .TRUE.,                            &
                           Compress     = .TRUE.,                            &
                           fId          = Container%FileId,                  &
                           DataType     = 8,                                 &
                           VarName      = Current%Item%Name,                 &
                           VarCt        = Current%Item%NcVarId,              &
                           timeId       = Current%Item%NcTDimId,             &
                           levId        = Current%Item%NcZDimId,             &
                           latId        = Current%Item%NcYDimId,             &
                           lonId        = Current%Item%NcXDimId,             &
                           VarLongName  = Current%Item%LongName,             &
                           VarUnit      = VarUnits,                          &
                           Axis         = VarAxis,                           &
                           Calendar     = VarCalendar,                       &
                           Positive     = VarPositive                       )

          ! Debug print
          !CALL HistItem_Print( am_I_Root, Current%Item, RC )

          ! Go to next entry in the list of HISTORY ITEMS
          Current => Current%Next
       ENDDO

       ! Free pointers
       Current => NULL()

       !--------------------------------------------------------------------
       ! Then define each HISTORY ITEM belonging to this collection
       !--------------------------------------------------------------------

       ! Set CURRENT to the first node in the list of HISTORY ITEMS
       Current => Container%HistItems

       ! As long as this node of the list is valid ...
       DO WHILE( ASSOCIATED( Current ) )

          ! Get the dimension ID's that are relevant to each HISTORY ITEM
          ! and save them in fields of the HISTORY ITEM
          CALL Get_Var_DimIds( Item   = Current%Item,                        & 
                               xDimId = Container%xDimId,                    &
                               yDimId = Container%yDimId,                    &
                               zDimId = Container%zDimId,                    &
                               tDimId = Container%tDimId                    )

          ! Define each HISTORY ITEM in this collection to the netCDF file
          CALL Nc_Var_Def( DefMode      = .TRUE.,                            &
                           Compress     = .TRUE.,                            &
                           fId          = Container%FileId,                  &
                           DataType     = 4,                                 &
                           VarName      = Current%Item%Name,                 &
                           VarCt        = Current%Item%NcVarId,              &
                           timeId       = Current%Item%NcTDimId,             &
                           levId        = Current%Item%NcZDimId,             &
                           latId        = Current%Item%NcYDimId,             &
                           lonId        = Current%Item%NcXDimId,             &
                           VarLongName  = Current%Item%LongName,             &
                           VarUnit      = Current%Item%Units,                &
                           MissingValue = Current%Item%MissingValue,         &
                           AvgMethod    = Current%Item%AvgMethod            )
           
          ! Go to next entry in the list of HISTORY ITEMS
          Current => Current%Next
       ENDDO

       ! Free pointers
       Current => NULL()
 
       !--------------------------------------------------------------------
       ! Write the index variable data
       !--------------------------------------------------------------------

       ! Close definition section
       CALL Nc_Set_DefMode( Container%FileId, Off=.TRUE. )

       ! Set CURRENT to the first node in the list of HISTORY ITEMS
       Current => IndexVarList

       ! As long as this node of the list is valid ...
       DO WHILE( ASSOCIATED( Current ) )

          ! Write the data to the netCDF file
          SELECT CASE( TRIM( Current%Item%Name ) ) 
             CASE( 'AREA' ) 
                CALL Nc_Var_Write( fId     = Container%FileId,               &
                                   VarName = Current%Item%Name,              &
                                   Arr2d   = Current%Item%Source_2d         )
             CASE DEFAULT 
                CALL Nc_Var_Write( fId     = Container%FileId,               &
                                   VarName = Current%Item%Name,              &
                                   Arr1d   = Current%Item%Source_1d         )
           END SELECT

          ! Go to next entry in the list of HISTORY ITEMS
          Current => Current%Next
       ENDDO

       ! Free pointers
       Current => NULL()

       !--------------------------------------------------------------------
       ! We can now consider this collection to have been "defined"
       !--------------------------------------------------------------------
       Container%IsFileDefined = .TRUE.
    ENDIF

  END SUBROUTINE History_Netcdf_Define

!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: History_Netcdf_Write
!
! !DESCRIPTION: Writes the data contained in each HISTORY ITEM to the netCDF
!  file specified by a given HISTORY CONTAINER object.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE History_Netcdf_Write( am_I_Root, Container,                     &
                                   yyyymmdd,  hhmmss,    RC                 )
!
! !USES:
!
    USE ErrCode_Mod
    USE Gc_Grid_Mod,        ONLY : RoundOff
    USE HistItem_Mod,       ONLY : HistItem
    USE HistContainer_Mod,  ONLY : HistContainer
    USE History_Params_Mod
    USE M_NetCdf_Io_Write,  ONLY : NcWr
    USE MetaHistItem_Mod,   ONLY : MetaHistItem
!
! !INPUT PARAMETERS: 
!
    LOGICAL,             INTENT(IN)  :: am_I_Root ! Are we on the root CPU?
    INTEGER,             INTENT(IN)  :: yyyymmdd  ! Current Year/month/day
    INTEGER,             INTENT(IN)  :: hhmmss    ! Current hour/minute/second
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(HistContainer), POINTER     :: Container ! Diagnostic collection obj
!
! !OUTPUT PARAMETERS: 
!
    INTEGER,             INTENT(OUT) :: RC        ! Success or failure
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
    INTEGER                     :: NcFileId,         NcVarId
    INTEGER                     :: Dim1,             Dim2,       Dim3
    REAL(f8)                    :: ElapsedMin,       Jd,         OffsetMin

    ! Strings
    CHARACTER(LEN=255)          :: ErrMsg,           ThisLoc

    ! Arrays
    INTEGER                     :: St1d(1),          Ct1d(1)
    INTEGER                     :: St2d(2),          Ct2d(2)
    INTEGER                     :: St3d(3),          Ct3d(3)
    INTEGER                     :: St4d(4),          Ct4d(4)
    REAL(f4),       ALLOCATABLE :: NcData_1d(:    )
    REAL(f4),       ALLOCATABLE :: NcData_2d(:,:  )
    REAL(f4),       ALLOCATABLE :: NcData_3d(:,:,:)
    REAL(f8)                    :: NcTimeVal(1    )

    ! Objects
    TYPE(MetaHistItem), POINTER :: Current 
    TYPE(HistItem),     POINTER :: Item

    !=======================================================================
    ! Make sure the netCDF file is open and defined
    !=======================================================================
    IF ( ( .not. Container%IsFileOpen   )    .and.                          &
         ( .not. Container%IsFileDefined ) ) THEN 
       RC     = GC_FAILURE
       ErrMsg = 'NetCDF file is not open or defined for collection: ' // &
                 TRIM( Container%Name ) 
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !=======================================================================
    ! Initialize 
    !=======================================================================
    RC        =  GC_SUCCESS
    Dim1      =  UNDEFINED_INT
    Dim2      =  UNDEFINED_INT
    Dim3      =  UNDEFINED_INT
    NcFileId  =  Container%FileId
    Current   => NULL()
    Item      => NULL()
    ErrMsg    =  ''
    ThisLoc   =  &
         ' -> at History_Netcdf_Write (in History/history_netcdf_mod.F90)'

    !=======================================================================
    ! Compute time elapsed since the reference time
    !=======================================================================

    ! Increment the time index for the netCDF file
    Container%CurrTimeSlice = Container%CurrTimeSlice + 1

    ! Compute the Astronomical Julian Date at this time
    CALL Compute_Julian_Date( yyyymmdd, hhmmss, Jd )

    ! Compute the time in minutes elapsed since the reference time
    ElapsedMin = ( ( Jd - Container%ReferenceJd ) * 1440.0_f8 )

    ! Round off to 5 decimal places, this should be OK even
    ! if we have fractional minutes as time values.
    ElapsedMin = RoundOff( ElapsedMin, 5 )

    !=======================================================================
    ! Compute the time stamp value for the current time slice
    !=======================================================================

    ! Compute the timestamp offset
    CALL Compute_TimeAvgOffset( Container, yyyymmdd, hhmmss, OffsetMin )

    ! Time stamp for current time slice
    Container%TimeStamp = ElapsedMin + OffsetMin

    !=======================================================================
    ! Write the time stamp to the netCDF File
    !=======================================================================

    ! netCDF start and count arrays
    St1d      = (/ Container%CurrTimeSlice /)
    Ct1d      = (/ 1                       /)

    ! Time stamp value
    NcTimeVal = (/ Container%TimeStamp     /)

    ! Write the time stamp to the file
    CALL NcWr( NcTimeVal, NcFileId, 'time', St1d, Ct1d )

    !=======================================================================
    ! Loop over all of the HISTORY ITEMS belonging to this collection
    !=======================================================================

    ! Set CURRENT to the first entry in the list of HISTORY ITEMS
    Current => Container%HistItems

    ! As long as this entry of the list is valid ...
    DO WHILE( ASSOCIATED( Current ) )

       ! Point to the HISTORY ITEM object in this entry
       Item => Current%Item

       !--------------------------------------------------------------------
       ! For instantaneous diagnostic quantities:
       ! (1) Copy the Item's data array to the 4-byte local array
       ! (2) Zero the Item's data array
       ! (3) Zero the Item's update counter
       !
       ! For time-averaged diagnostic quantities:
       ! (1) Divide the Item's data array by the number diagnostic updates
       ! (2) Copy the Item's data array to the 4-byte local array
       ! (3) Zero the Item's data array
       ! (4) Zero the Item's update counter
       !--------------------------------------------------------------------
       SELECT CASE( Item%SpaceDim )
          
          !------------
          ! 3-D data
          !------------
          CASE( 3 )   
             
             ! Get dimensions of data
             Dim1 = SIZE( Item%Data_3d, 1 )
             Dim2 = SIZE( Item%Data_3d, 2 )
             Dim3 = SIZE( Item%Data_3d, 3 )

             ! Allocate the REAL*4 output array
             ALLOCATE( NcData_3d( Dim1, Dim2, Dim3 ), STAT=RC )

             ! Copy or average the data and store in a REAL*4 array
             IF ( Item%Operation == COPY_FROM_SOURCE ) THEN
                NcData_3d     = Item%Data_3d
                Item%Data_3d  = 0.0_f8
                Item%nUpdates = 0
             ELSE
                Item%Data_3d  = Item%Data_3d / DBLE( Item%nUpdates )
                NcData_3d     = Item%Data_3d
                Item%Data_3d  = 0.0_f8
                Item%nUpdates = 0
             ENDIF

             ! Compute start and count fields
             St4d = (/ 1,    1,    1,    Container%CurrTimeSlice /)
             Ct4d = (/ Dim1, Dim2, Dim3, 1                       /)

             ! Write data to disk
             CALL NcWr( NcData_3d, NcFileId, Item%Name, St4d, Ct4d )
             
             ! Deallocate output array
             DEALLOCATE( NcData_3d, STAT=RC )

          !------------
          ! 2-D data
          !------------
          CASE( 2 )

             ! Get dimensions of data
             Dim1 = SIZE( Item%Data_3d, 1 )
             Dim2 = SIZE( Item%Data_3d, 2 )

             ! Allocate the REAL*4 output array
             ALLOCATE( NcData_2d( Dim1, Dim2 ), STAT=RC )
 
             ! Copy or average the data and store in a REAL*4 array
             IF ( Item%Operation == COPY_FROM_SOURCE ) THEN
                NcData_2d     = Item%Data_2d
                Item%Data_2d  = 0.0_f8
                Item%nUpdates = 0
             ELSE
                Item%Data_2d  = Item%Data_2d / DBLE( Item%nUpdates )
                NcData_2d     = Item%Data_2d
                Item%Data_2d  = 0.0_f8
                Item%nUpdates = 0
             ENDIF

             ! Compute start and count fields
             St3d = (/ 1,    1,    Container%CurrTimeSlice /)
             Ct3d = (/ Dim1, Dim2, 1                       /)

             ! Write data to disk
             CALL NcWr( NcData_3d, NcFileId, Item%Name, St3d, Ct3d )
             
             ! Deallocate output array
             DEALLOCATE( NcData_2d, STAT=RC )

          !------------
          ! 1-D data
          !------------
          CASE( 1 )

             ! Get dimensions of data
             Dim1 = SIZE( Item%Data_3d, 1 )

             ! Allocate the REAL*4 output array
             ALLOCATE( NcData_1d( Dim1 ), STAT=RC )

             ! Copy or average the data and store in a REAL*4 array
             IF ( Item%Operation == COPY_FROM_SOURCE ) THEN
                NcData_1d     = Item%Data_1d
                Item%Data_1d  = 0.0_f8
                Item%nUpdates = 0
             ELSE
                Item%Data_1d  = Item%Data_1d / DBLE( Item%nUpdates )
                NcData_1d     = Item%Data_1d
                Item%Data_1d  = 0.0_f8
                Item%nUpdates = 0
             ENDIF


             ! Compute start and count fields
             St2d = (/ 1,    Container%CurrTimeSlice /)
             Ct2d = (/ Dim1, 1                       /)

             ! Write data to disk
             CALL NcWr( NcData_1d, NcFileId, Item%Name, St2d, Ct2d )
             
             ! Deallocate output array
             DEALLOCATE( NcData_2d, STAT=RC )

       END SELECT

       !--------------------------------------------------------------------
       ! Go to next entry in the list of HISTORY ITEMS
       !-------------------------------------------------------------------- 
       Current => Current%Next
       Item    => NULL()
    ENDDO

    !========================================================================
    ! Cleanup and quit
    !========================================================================

    ! Free pointers
    Current => NULL()
    Item    => NULL()

  END SUBROUTINE History_NetCdf_Write
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Compute_Julian_Date
!
! !DESCRIPTION: Computes the Astronomical Julian Date corresponding to a 
!  given date and time.  This is useful for computing elapsed times.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Compute_Julian_Date( yyyymmdd, hhmmss, Jd )
!
! !USES:
!
    USE Julday_Mod, ONLY : Julday
    USE Time_Mod,   ONLY : Ymd_Extract
!
! !INPUT PARAMETERS: 
!
    INTEGER,  INTENT(IN)  :: yyyymmdd  ! Current Year/month/day
    INTEGER,  INTENT(IN)  :: hhmmss    ! Current hour/minute/second
!
! !OUTPUT PARAMETERS: 
!
    REAL(f8), INTENT(OUT) :: Jd        ! Astronomical Julian date
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
    INTEGER  :: Year, Month, Day, Hour, Minute, Second
    REAL(f8) :: FracDay

    ! Extract year/month/day and hour/minute/seconds from the time
    CALL Ymd_Extract( yyyymmdd, Year, Month,  Day    )
    CALL Ymd_Extract( hhmmss,   Hour, Minute, Second )

    ! Compute the fractional day
    FracDay = DBLE( Day ) + ( DBLE( Hour   ) / 24.0_f8    )  +               & 
                            ( DBLE( Minute ) / 1440.0_f8  )  +               &
                            ( DBLE( Second ) / 86400.0_f8 ) 

    ! Return the Astronomical Julian Date
    Jd = JulDay( Year, Month, FracDay )

  END SUBROUTINE Compute_Julian_Date
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Compute_TimeAvgOffset
!
! !DESCRIPTION: Computes the offset (in minutes) for the netCDF timestamp. 
!  Instantaneous data will have an offset of zero.  Time-averaged data will
!  have an offset of 1/2 the accumulation interval.
!\\
!\\
! !INTERFACE:
!
 SUBROUTINE Compute_TimeAvgOffset( Container, yyyymmdd, hhmmss, OffsetMin )
!
! !USES:
!
   USE HistContainer_Mod,  ONLY : HistContainer
   USE History_Params_Mod
   USE Time_Mod,           ONLY : Ymd_Extract
!
! !INPUT PARAMETERS:
!
   TYPE(HistContainer), POINTER     :: Container   ! HISTORY CONTAINER object
   INTEGER,             INTENT(IN)  :: yyyymmdd    ! Current date in YMD
   INTEGER,             INTENT(IN)  :: hhmmss      ! Current time in hms
!
! !OUTPUT PARAMETERS:
!
   REAL(f8),            INTENT(OUT) :: OffsetMin   ! Timestamp offset in min
!
! !REMARKS:
!  Instantaneous data has an offset of 0.0 minutes.
!
! !REVISION HISTORY:
!  06 Jan 2015 - R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
   INTEGER :: Year, Month, Day, Hour, Minute, Second

   IF ( Container%Operation == COPY_FROM_SOURCE ) THEN

      !------------------------------------------
      ! Instantaneous data: offset is zero.
      !------------------------------------------
      OffsetMin = 0.0_f8
      RETURN

   ELSE

      !------------------------------------------
      ! Time averaged data: compute offset
      !------------------------------------------
      OffsetMin = 0.0_f8

      ! Split file writing date & times into individual constitutents
      CALL Ymd_Extract( Container%FileWriteYmd, Year, Month,  Day    )
      CALL Ymd_Extract( Container%FileWriteHms, Hour, Minute, Second )

      IF ( Day > 0 ) THEN
         OffsetMin = OffsetMin + ( DBLE( Day    ) * 1440.0_f8 )
      ENDIF

      IF ( Hour > 0 ) THEN
         OffsetMin = OffsetMin + ( DBLE( Hour   ) * 60.0_f8   )
      ENDIF

      IF ( Minute > 0 ) THEN
         OffsetMin = OffsetMin + ( DBLE( Minute )             )
      ENDIF

      IF ( Second > 0 ) THEN
         OffsetMin = OffsetMin + ( DBLE( Second ) / 60.0_f8   )
      ENDIF

      ! Take the midpoint of the interval
      OffsetMin = OffsetMin * 0.5_f8

   ENDIF

 END SUBROUTINE Compute_TimeAvgOffset
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Expand_Date_Time
!
! !DESCRIPTION: Replaces date and time tokens in a string with actual 
!  date and time values.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Expand_Date_Time( DateStr, yyyymmdd, hhmmss, MAPL_Style )
!
! !USES:
!
    USE Charpak_Mod, ONLY : StrRepl
    USE Time_Mod,    ONLY : Ymd_Extract
!
! !INPUT PARAMETERS: 
!
    INTEGER,          INTENT(IN)    :: yyyymmdd    ! Date in YYYYMMDD format
    INTEGER,          INTENT(IN)    :: hhmmss      ! Time in hhmmss format
    LOGICAL,          OPTIONAL      :: MAPL_Style  ! Use MAPL-style tokens
!
! !INPUT/OUTPUT PARAMETERS: 
!
    CHARACTER(LEN=*), INTENT(INOUT) :: DateStr     ! String with date tokens
!
! !REMARKS:
!  Based on EXPAND_DATE from GeosUtil/time_mod.F.
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
    LOGICAL          :: Is_Mapl_Style 
    INTEGER          :: Year,          Month,      Day
    INTEGER          :: Hour,          Minute,     Second

    ! Strings
    CHARACTER(LEN=2) :: MonthStr,      DayStr 
    CHARACTER(LEN=2) :: HourStr,       MinuteStr,  SecondStr
    CHARACTER(LEN=4) :: YearStr

    !=======================================================================
    ! Initialize
    !=======================================================================
    IF ( PRESENT( MAPL_Style ) ) THEN
       Is_Mapl_Style = MAPL_Style
    ELSE
       Is_Mapl_Style = .FALSE.
    ENDIF

    !=======================================================================
    ! Split the date and time into individal variables
    !=======================================================================

    ! Extract year/month/day and hour/minute/seconds from the time
    CALL Ymd_Extract( yyyymmdd, Year, Month,  Day    )
    CALL Ymd_Extract( hhmmss,   Hour, Minute, Second )

    ! Convert to strings
    WRITE( YearStr,   '(i4.4)' ) Year
    WRITE( MonthStr,  '(i2.2)' ) Month
    WRITE( DayStr,    '(i2.2)' ) Day
    WRITE( HourStr,   '(i2.2)' ) Hour
    WRITE( MinuteStr, '(i2.2)' ) Minute
    WRITE( SecondStr, '(i2.2)' ) Second

    !=======================================================================
    ! Replace the date and time tokens in the string
    !=======================================================================

    IF ( Is_Mapl_Style ) THEN
       
       ! Use MAPL-style tokens
       CALL StrRepl( DateStr, '%y4',  YearStr   )
       CALL StrRepl( DateStr, '%m2',  MonthStr  )
       CALL StrRepl( DateStr, '%d2',  DayStr    )
       CALL StrRepl( DateStr, '%h2',  HourStr   )
       CALL StrRepl( DateStr, '%n2',  MinuteStr )
       CALL StrRepl( DateStr, '%s2',  SecondStr )

    ELSE
       
       ! Use GEOS-Chem style tokens
       CALL StrRepl( DateStr, 'YYYY', YearStr   )
       CALL StrRepl( DateStr, 'MM',   MonthStr  )
       CALL StrRepl( DateStr, 'DD',   DayStr    )
       CALL StrRepl( DateStr, 'hh',   HourStr   )
       CALL StrRepl( DateStr, 'mm',   MinuteStr )
       CALL StrRepl( DateStr, 'ss',   SecondStr )

    ENDIF

  END SUBROUTINE Expand_Date_Time
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get_Var_DimIds
!
! !DESCRIPTION: For a given HISTORY ITEM, returns the name and the netCDF
!  dimension ID's pertaining to the data array.  Dimension ID's that do not
!  pertain to the data will be set to UNDEFINED_INT.  Certain metadata for
!  netCDF index variables will also be returned.  In particular, the unit 
!  string for the "time" index variable will be updated with the reference 
!  date and time.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Get_Var_DimIds( Item,        xDimId,      yDimId,  zDimId,      &
                             tDimID,      RefDate,     RefTime, VarAxis,     &
                             VarCalendar, VarPositive, VarUnits             )
!
! !USES:
!
    USE History_Params_Mod
    USE HistItem_Mod,       ONLY : HistItem
!
! !INPUT PARAMETERS: 
!
    TYPE(HistItem),     POINTER     :: Item        ! HISTORY ITEM object
    INTEGER,            INTENT(IN)  :: xDimId      ! Id # of netCDF X dimension 
    INTEGER,            INTENT(IN)  :: yDimId      ! Id # of netCDF Y dimension 
    INTEGER,            INTENT(IN)  :: zDimId      ! Id # of netCDF Z dimension 
    INTEGER,            INTENT(IN)  :: tDimId      ! Id # of netCDF T dimension 
    INTEGER,            OPTIONAL    :: RefDate     ! Ref YMD for "time" variable
    INTEGER,            OPTIONAL    :: RefTime     ! Ref hms for "time" variable
!
! !OUTPUT PARAMETERS
!
    CHARACTER(LEN=255), OPTIONAL    :: VarAxis     ! Axis attr for index vars
    CHARACTER(LEN=255), OPTIONAL    :: VarCalendar ! Calendar attr for "time"
    CHARACTER(LEN=255), OPTIONAL    :: VarPositive ! Positive attr for "lev"
    CHARACTER(LEN=255), OPTIONAL    :: VarUnits    ! Unit string
!
! !REMARKS:
!  Call this routine before calling NC_VAR_DEF.
!
! !REVISION HISTORY:
!  10 Aug 2017 - R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER            :: ReferenceYmd, ReferenceHms  
    INTEGER            :: Year,         Month,        Day
    INTEGER            :: Hour,         Minute,       Second

    ! Strings
    CHARACTER(LEN=2)   :: MonthStr,     DayStr 
    CHARACTER(LEN=2)   :: HourStr,      MinuteStr,    SecondStr
    CHARACTER(LEN=4)   :: YearStr
    CHARACTER(LEN=255) :: TmpAxis,      TmpCalendar
    CHARACTER(LEN=255) :: TmpPositive,  TmpUnits

    !=======================================================================
    ! Initialize
    !=======================================================================
    TmpAxis         = ''
    TmpCalendar     = ''
    TmpPositive     = ''
    TmpUnits        = Item%Units
    Item%NcXDimId   = UNDEFINED_INT
    Item%NcYDimId   = UNDEFINED_INT
    Item%NcZDimId   = UNDEFINED_INT
    Item%NcTDimId   = UNDEFINED_INT
    
    IF ( PRESENT( RefDate ) ) THEN
       ReferenceYmd = RefDate
    ELSE
       ReferenceYmd = UNDEFINED_INT
    ENDIF

    IF ( PRESENT( RefTime ) ) THEN
       ReferenceHms = RefTime
    ELSE
       ReferenceHms = UNDEFINED_INT
    ENDIF

    !=======================================================================
    ! Return relevant dim ID's and metadata for the HISTORY ITEM
    !=======================================================================
    SELECT CASE( TRIM( Item%Name ) )

       ! lon
       CASE( 'lon' )
          Item%NcXDimId = xDimId
          TmpAxis       = 'X'

       ! lat
       CASE( 'lat' )
          Item%NcYDimId = yDimId
          TmpAxis       = 'Y'

       ! lev
       CASE( 'lev' )
          Item%NcZDimId = zDimId
          TmpAxis       = 'Z'
          TmpPositive   = 'up'

       ! time
       CASE( 'time' )
          Item%NcTDimId = tDimId
          TmpAxis       = 'T'
          TmpCalendar   = 'gregorian'

          ! Replace date and time tokens in the unit string
          ! with the netCDF file's reference date and time
          IF ( ReferenceYmd > 0 ) THEN
             CALL Expand_Date_Time( TmpUnits, ReferenceYmd, ReferenceHms )
          ENDIF

       ! area
       CASE( 'AREA' )
          Item%NcXDimId = xDimId
          Item%NcYDimId = yDimId

       ! All other variable names
       CASE DEFAULT

          ! Pick the 
          SELECT CASE( TRIM( Item%DimNames ) )
             CASE( 'xyz' )
                Item%NcXDimId = xDimId
                Item%NcYDimId = yDimId
                Item%NcZDimId = ZDimId
                Item%NcTDimId = tDimId
             CASE( 'xy'  )
                Item%NcXDimId = xDimId
                Item%NcYDimId = yDimId
                Item%NcTDimId = tDimId
             CASE( 'yz'  )
                Item%NcYDimId = yDimId
                Item%NcZDimId = zDimId
                Item%NcTDimId = tDimId
             CASE( 'xz'  )
                Item%NcXDimId = xDimId
                Item%NcZDimId = zDimId
                Item%NcTDimId = tDimId
             CASE( 'x'   )
                Item%NcXDimId = xDimId
                Item%NcTDimId = tDimId
             CASE( 'y'   )
                Item%NcYDimId = yDimId
                Item%NcTDimId = tDimId
             CASE( 'z'   )
                Item%NcZDimId = zDimId
                Item%NcTDimId = tDimId
          END SELECT

    END SELECT      
    
    ! Return optional attributes for index variables: axis and calendar
    IF ( PRESENT( VarAxis     ) ) VarAxis     = TmpAxis
    IF ( PRESENT( VarCalendar ) ) VarCalendar = TmpCalendar
    IF ( PRESENT( VarPositive ) ) VarPositive = TmpPositive
    IF ( PRESENT( VarUnits    ) ) VarUnits    = TmpUnits

  END SUBROUTINE Get_Var_DimIds
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Set_RefDateTime
!
! !DESCRIPTION: Defines the reference date and time fields of the HISTORY
!  CONTAINER object.  These are needed to compute the timestamps along the
!  time dimension of each HISTORY ITEM.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Set_RefDateTime( Container, yyyymmdd, hhmmss )
!
! !USES:
!
    USE HistContainer_Mod, ONLY : HistContainer
    USE Time_Mod,          ONLY :
!
! !INPUT PARAMETERS: 
!
    INTEGER,             INTENT(IN) :: yyyymmdd   ! Current Year/month/day
    INTEGER,             INTENT(IN) :: hhmmss     ! Current hour/minute/second
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(HistContainer), POINTER    :: Container  ! Diagnostic collection obj
!
! !REMARKS:
!  Also initializes the Container%CurrTimeSlice field.
!
! !REVISION HISTORY:
!  06 Jan 2015 - R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC

    !=======================================================================
    ! Compute the reference date/time quantities
    !=======================================================================

    ! Current time slice (will increment this in History_Netcdf_Write)
    Container%CurrTimeSlice = 0

    ! Reference date
    Container%ReferenceYmd  = yyyymmdd

    ! Reference time
    Container%ReferenceHms  = hhmmss

    ! Corresponding astronomical Julian date
    CALL Compute_Julian_Date( yyyymmdd, hhmmss, Container%ReferenceJd )

  END SUBROUTINE Set_RefDateTime
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: History_Netcdf_Init
!
! !DESCRIPTION: Creates a HISTORY ITEM for each netCDF index variable (e.g.
!  lon, lat, lev, time, area) and adds it to the METAHISTORY ITEM IndexVarList.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE History_Netcdf_Init( am_I_Root, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE GC_Grid_Mod,      ONLY : Lookup_Grid
    USE HistItem_Mod
    USE MetaHistItem_Mod
!
! !INPUT PARAMETERS: 
!
    LOGICAL, INTENT(IN)  :: am_I_Root
!
! !OUTPUT PARAMETERS: 
!
    INTEGER, INTENT(OUT) :: RC
!
! !REMARKS:
!
! !REVISION HISTORY:
!  10 Aug 2017 - R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER                  :: N
    INTEGER                  :: KindVal
    INTEGER                  :: Rank
    INTEGER                  :: Dimensions(3)

    ! Strings
    CHARACTER(LEN=3)         :: ItemDimNames(5)
    CHARACTER(LEN=4)         :: ItemName(5)
    CHARACTER(LEN=9)         :: RegistryName(5)
    CHARACTER(LEN=255)       :: Description
    CHARACTER(LEN=255)       :: ErrMsg
    CHARACTER(LEN=255)       :: ThisLoc
    CHARACTER(LEN=255)       :: Units

    ! Pointer arrays
    REAL(fp),        POINTER :: Ptr1d(:    )
    REAL(fp),        POINTER :: Ptr2d(:,:  )

    ! Objects
    TYPE(HistItem),  POINTER :: Item

    !=======================================================================
    ! Initialize
    !=======================================================================
    RC           =  GC_SUCCESS
    Description  =  ''
    Dimensions   =  0
    KindVal      =  0
    Rank         =  0
    Units        =  ''
    ErrMsg       =  ''
    ThisLoc      =  &
                ' -> at History_Netcdf_Init (in History/history_mod.F90)'
    RegistryName = (/ 'GRID_AREA', 'GRID_TIME', 'GRID_LEV ',                 &
                      'GRID_LAT ', 'GRID_LON '                             /)
    ItemName     = (/ 'AREA', 'time', 'lev ', 'lat ', 'lon '               /)
    ItemDimNames = (/ 'xy ' , 't  ' , 'z  ' , 'y  ' , 'x  '                /)
    Ptr1d        => NULL()
    Ptr2d        => NULL()

    !=======================================================================
    ! Create a HISTORY ITEM for each of the index fields (lon, lat, area)
    ! of gc_grid_mod.F90 and add them to a METAHISTORY ITEM list
    !=======================================================================
    DO N = 1, SIZE( ItemName )

       !---------------------------------------------------------------------
       ! Look up one of the index fields from gc_grid_mod.F90
       !---------------------------------------------------------------------
       CALL Lookup_Grid( am_I_Root   = am_I_Root,                            &
                         Variable    = RegistryName(N),                      &
                         Description = Description,                          &
                         Dimensions  = Dimensions,                           &
                         KindVal     = KindVal,                              &
                         Rank        = Rank,                                 &
                         Units       = Units,                                &
                         Ptr1d       = Ptr1d,                                &
                         Ptr2d       = Ptr2d,                                &
                         RC          = RC                                   )

       ! Trap potential error
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Error in "Lookup_Grid" for diagnostic ' //               &
                   TRIM( RegistryName(N) )
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       !---------------------------------------------------------------------
       ! Create a HISTORY ITEM for this index field
       !---------------------------------------------------------------------
       CALL HistItem_Create( am_I_Root      = am_I_Root,                     &
                             Item           = Item,                          &
                             Id             = N,                             & 
                             ContainerId    = 0,                             &
                             Name           = ItemName(N),                   &
                             LongName       = Description,                   &
                             Units          = Units,                         &
                             SpaceDim       = Rank,                          &
                             NX             = Dimensions(1),                 &
                             NY             = Dimensions(2),                 &
                             NZ             = Dimensions(3),                 &
                             DimNames       = ItemDimNames(N),               &
                             Operation      = 0,                             &
                             Source_KindVal = KindVal,                       &
                             Source_1d      = Ptr1d,                         &
                             Source_2d      = Ptr2d,                         &
                             RC             = RC                            )

       ! Trap potential error
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Could not create Item: "' // TRIM( ItemName(N) ) // '"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       !CALL HistItem_Print( am_I_Root, Item, RC )

       !---------------------------------------------------------------------
       ! Add this item to the Dimension list
       !---------------------------------------------------------------------
       CALL MetaHistItem_AddNew( am_I_Root = am_I_Root,                      &
                                 Node      = IndexVarList,                   &
                                 Item      = Item,                           &
                                 RC        = RC                             )

       ! Trap potential error
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Could not add Item "' // TRIM( ItemName(N) ) //          &
                   '" to Dimensionlist!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDDO

  END SUBROUTINE History_Netcdf_Init
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: History_Netcdf_Cleanup
!
! !DESCRIPTION: Finalizes all module variables (e.g. IndexVarList).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE History_Netcdf_Cleanup( am_I_Root, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE MetaHistItem_Mod,  ONLY: MetaHistItem_Destroy
!
! !INPUT PARAMETERS: 
!
    LOGICAL, INTENT(IN)  :: am_I_Root
!
! !OUTPUT PARAMETERS: 
!
    INTEGER, INTENT(OUT) :: RC
!
! !REMARKS:
!
! !REVISION HISTORY:
!  10 Aug 2017 - R. Yantosca - Initial version
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
    ErrMsg  =  ''
    ThisLoc =  ' -> at MetaHistItem_Destroy (in History/metahistitem_mod.F90)'

    !=======================================================================
    ! Destroy the METAHISTORY ITEM list of index variables for netCDF
    !=======================================================================
    CALL MetaHistItem_Destroy( am_I_Root, IndexVarList, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Cannot deallocate the "IndexVarListt" META HISTORY ITEM!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

  END SUBROUTINE History_Netcdf_Cleanup
!EOC
END MODULE History_Netcdf_Mod

