!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: history_netcdf_mod.F90 
!
! !DESCRIPTION: 
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
  PUBLIC  :: History_Netcdf_Define
  PUBLIC  :: History_Netcdf_Init
  PUBLIC  :: History_Netcdf_Cleanup
!
! !PRIVATE MEMBER FUNCTIONS:
!
!  PRIVATE :: History_Netcdf_DefineDims
!
! !REMARKS:
!
! !REVISION HISTORY:
!  06 Jan 2015 - R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE TYPES:
!
  TYPE(MetaHistItem), POINTER :: IndexVarList => NULL()

CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: History_Netcdf_Define
!
! !DESCRIPTION: 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE History_NetCdf_Define( am_I_Root, Container,                    &
                                    yyyymmdd,  hhmmss,    RC                )
!
! !USES:
!
    USE Charpak_Mod,        ONLY : StrRepl
    USE ErrCode_Mod
    USE HistContainer_Mod,  ONLY : HistContainer, HistContainer_Print
    USE HistItem_Mod,       ONLY : HistItem,      HistItem_Print
    USE MetaHistItem_Mod,   ONLY : MetaHistItem
    USE Ncdf_Mod
    USE Time_Mod,           ONLY : Ymd_Extract
!
! !INPUT PARAMETERS: 
!
    LOGICAL,             INTENT(IN)  :: am_I_Root ! Are we on the root CPU?
    TYPE(HistContainer), POINTER     :: Container ! Diagnostic collection obj
    INTEGER,             INTENT(IN)  :: yyyymmdd  ! Current Year/month/day
    INTEGER,             INTENT(IN)  :: hhmmss    ! Current hour/minute/second
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
    INTEGER                     :: N,           VarCt
    INTEGER                     :: Year,        Month,       Day
    INTEGER                     :: Hour,        Minute,      Second
    INTEGER                     :: VarXDimId,   VarYDimId
    INTEGER                     :: VarZDimId,   VarTDimId

    ! Strings                   
    CHARACTER(LEN=2)            :: MonthStr,    DayStr 
    CHARACTER(LEN=2)            :: HourStr,     MinuteStr,   SecondStr
    CHARACTER(LEN=4)            :: YearStr
    CHARACTER(LEN=5)            :: Z
    CHARACTER(LEN=8)            :: D
    CHARACTER(LEN=10)           :: T
    CHARACTER(LEN=255)          :: FileName,    TimeString
    CHARACTER(LEN=255)          :: ErrMsg,      ThisLoc
    CHARACTER(LEN=255)          :: VarName,     VarAxis    
    CHARACTER(LEN=255)          :: VarPositive, VarCalendar

    ! Arrays
    INTEGER                     :: V(8)

    ! Pointers
    REAL(fp),           POINTER :: Data1d(:)
    REAL(fp),           POINTER :: Data2d(:)

    ! Objects
    TYPE(MetaHistItem), POINTER :: Current
    TYPE(HistItem),     POINTER :: Item
    
    !=======================================================================
    ! Initialize
    !=======================================================================
    RC      =  GC_SUCCESS
    Current => NULL()
    Item    => NULL()
    ErrMsg  =  ''
    ThisLoc =  ' -> at History_Netcdf_Define (in History/history_mod.F90)'

    !=======================================================================
    ! Construct the file name using the passed date and time
    !=======================================================================

    ! Save the collection's file name in a temporary variable
    FileName = TRIM( Container%FileName )

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

    ! Replace
    CALL StrRepl( FileName, '%y4', YearStr   )
    CALL StrRepl( FileName, '%m2', MonthStr  )
    CALL StrRepl( FileName, '%d2', DayStr    )
    CALL StrRepl( FileName, '%h2', HourStr   )
    CALL StrRepl( FileName, '%n2', MinuteStr )
    CALL StrRepl( FileName, '%s2', SecondStr )
   
#if defined( DEBUG )
    !###  Debug output for development
    WRITE(6,'(a     )') REPEAT( '#', 79 )
    WRITE(6,'(a,a   )') '### Time to write : ', TRIM( Container%Name )
    WRITE(6,'(a,i8.8)') '### YYYYMMDD      : ', yyyymmdd
    WRITE(6,'(a,i6.6)') '### hhmmss        : ', hhmmss
    WRITE(6,'(a,L3  )') '### IsFileDefined : ', Container%IsFileDefined
    WRITE(6,'(a,i8.8)') '### FileWriteYmd  : ', Container%FileWriteYmd
    WRITE(6,'(a,i6.6)') '### FileWriteHms  : ', Container%FileWriteHms
    WRITE(6,'(a,a   )') '### FileName      : ', TRIM( FileName )
    WRITE(6,'(a     )') REPEAT( '#', 79 )
#endif

    !=======================================================================
    ! Create the timestamp for the 
    !=======================================================================

    ! Call F90 intrinsic DATE_AND_TIME Function
    D = 'ccyymmdd'
    T = 'hhmmss.sss'
    CALL Date_And_Time( Date=D, Time=T, Zone=Z, Values=V )  ! GMT time
     
    ! Create timestamp strings
    WRITE( Container%History,      10 ) V(1), V(2), V(3), V(5), V(6), V(7), Z
    WRITE( Container%ProdDateTime, 10 ) V(1), V(2), V(3), V(5), V(6), V(7), Z
 10 FORMAT( 'Produced on ', i4.4, '/', i2.2, '/', i2.2, 1x,                  &
                            i2.2, ':', i2.2, ':', i2.2, ' UTC', a           )

    !=======================================================================
    ! Close the previous file in this collection if it is still open,
    ! so that we can generate the file for the next iteration.
    !=======================================================================
    IF ( Container%IsFileOpen ) THEN
       CALL Nc_Close( Container%FileId )
       Container%IsFileOpen    = .FALSE.
       Container%IsFileDefined = .FALSE.
    ENDIF

    !=======================================================================
    ! Create the netCDF file with global attributes
    ! Do not exit netCDF define mode just yet
    !=======================================================================
    IF ( .not. Container%IsFileOpen ) THEN

       !--------------------------------------------------------------------
       ! Create the file and add global attributes
       ! Remain in netCDF define mode upon exiting this routine
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
          ! (non-relevant dimensions will be given the value UNDEFINED_INT).
          ! Also get the Axis and Calendar attributes for index variables
          CALL Get_Var_DimIds( Item        = Current%Item,                   & 
                               xDimId      = Container%xDimId,               &
                               yDimId      = Container%yDimId,               &
                               zDimId      = Container%zDimId,               &
                               tDimId      = Container%tDimId,               &
                               VarName     = VarName,                        &
                               VarXDimId   = VarXDimId,                      &
                               VarYDimId   = VarYDimId,                      &
                               VarZDimId   = VarZDimId,                      &
                               VarTDimId   = VarTDimId,                      &
                               VarAxis     = VarAxis,                        &
                               VarPositive = VarPositive,                    &
                               VarCalendar = VarCalendar                    )
          

          ! Define each HISTORY ITEM in this collection to the netCDF file
          CALL Nc_Var_Def( DefMode      = .TRUE.,                            &
                           Compress     = .TRUE.,                            &
                           fId          = Container%FileId,                  &
                           DataType     = 8,                                 &
                           VarName      = VarName,                           &
                           timeId       = VarTDimId,                         &
                           levId        = VarZDimId,                         &
                           latId        = VarYDimId,                         &
                           lonId        = VarXDimId,                         &
                           VarLongName  = Current%Item%LongName,             &
                           VarUnit      = Current%Item%Units,                &
                           Axis         = VarAxis,                           &
                           Calendar     = VarCalendar,                       &
                           Positive     = VarPositive,                       &
                           VarCt        = Current%Item%NcVarId              )

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
          ! (non-defined dimensions will be given -1)
          CALL Get_Var_DimIds( Item      = Current%Item,                     & 
                               xDimId    = Container%xDimId,                 &
                               yDimId    = Container%yDimId,                 &
                               zDimId    = Container%zDimId,                 &
                               tDimId    = Container%tDimId,                 &
                               VarName   = VarName,                          &
                               VarXDimId = VarXDimId,                        &
                               VarYDimId = VarYDimId,                        &
                               VarZDimId = VarZDimId,                        &
                               VarTDimId = VarTDimId                        )

          ! Define each HISTORY ITEM in this collection to the netCDF file
          CALL Nc_Var_Def( DefMode      = .TRUE.,                            &
                           Compress     = .TRUE.,                            &
                           fId          = Container%FileId,                  &
                           DataType     = 4,                                 &
                           VarName      = VarName,                           &
                           timeId       = VarTDimId,                         &
                           levId        = VarZDimId,                         &
                           latId        = VarYDimId,                         &
                           lonId        = VarXDimId,                         &
                           VarLongName  = Current%Item%LongName,             &
                           VarUnit      = Current%Item%Units,                &
                           MissingValue = Current%Item%MissingValue,         &
                           AvgMethod    = Current%Item%AvgMethod,            &
                           VarCt        = Current%Item%NcVarId              )
           
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
             CASE( 'GRID_AREA' ) 
                CALL Nc_Var_Write( fId     = Container%FileId,               &
                                   VarName = 'AREA',                         &
                                   Arr2d   = Current%Item%Source_2d         )
             CASE( 'GRID_TIME' ) 
                CALL Nc_Var_Write( fId     = Container%FileId,               &
                                   VarName = 'time',                         &
                                   Arr1d   = Current%Item%Source_1d         )
             CASE( 'GRID_LEV' ) 
                CALL Nc_Var_Write( fId     = Container%FileId,               &
                                   VarName = 'lev',                          &
                                   Arr1d   = Current%Item%Source_1d         )
             CASE( 'GRID_LAT' ) 
                CALL Nc_Var_Write( fId     = Container%FileId,               &
                                   VarName = 'lat',                          &
                                   Arr1d   = Current%Item%Source_1d         )
             CASE( 'GRID_LON' ) 
                CALL Nc_Var_Write( fId     = Container%FileId,               &
                                   VarName = 'lon',                          &
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
! !IROUTINE: Get_Var_DimIds
!
! !DESCRIPTION: For a given HISTORY ITEM, returns the name and the netCDF
!  dimension ID's pertaining to the data array.  Dimension ID's that do not
!  pertain to the data will be set to UNDEFINED_INT.  Also, the name of
!  index variables will be adjusted for COARDS compliance (i.e. "GRID_LON"
!  will be changed to "lon", etc).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Get_Var_DimIds( Item,       xDimId,    yDimId,    zDimId,       &
                             tDimID,     VarName,   VarXDimId, VarYDimId,    &
                             VarZDimId,  VarTDimId, VarAxis,   VarCalendar,  &
                             VarPositive                                    )
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
!
! !OUTPUT PARAMETERS
!
    CHARACTER(LEN=255), INTENT(OUT) :: VarName     ! Adjusted name for netCDF
    INTEGER,            INTENT(OUT) :: VarXDimId   ! Adjusted X dimension ID
    INTEGER,            INTENT(OUT) :: VarYDimId   ! Adjusted Y dimension ID
    INTEGER,            INTENT(OUT) :: VarZDimId   ! Adjusted Z dimension ID
    INTEGER,            INTENT(OUT) :: VarTDimId   ! Adjusted T dimension ID
    CHARACTER(LEN=255), OPTIONAL    :: VarAxis     ! Axis attr for index vars
    CHARACTER(LEN=255), OPTIONAL    :: VarCalendar ! Calendar attr for time
    CHARACTER(LEN=255), OPTIONAL    :: VarPositive ! Positive attr for lev
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
    ! Strings
    CHARACTER(LEN=255) :: TmpAxis, TmpCalendar, TmpPositive

    !=======================================================================
    ! Initialize
    !=======================================================================
    VarXDimId   = UNDEFINED_INT
    VarYDimId   = UNDEFINED_INT
    VarZDimId   = UNDEFINED_INT
    VarTDimId   = UNDEFINED_INT
    TmpAxis     = ''
    TmpCalendar = ''
    TmpPositive = ''

    !=======================================================================
    ! Return relevant dim ID's and metadata for the HISTORY ITEM
    !=======================================================================
    SELECT CASE( TRIM( Item%Name ) )

       ! lon
       CASE( 'GRID_LON' )
          VarName     = 'lon'
          VarXDimId   = xDimId
          TmpAxis     = 'X'

       ! lat
       CASE( 'GRID_LAT' )
          VarName     = 'lat'
          VarYDimId   = yDimId
          TmpAxis     = 'Y'

       ! lev
       CASE( 'GRID_LEV' )
          VarName     = 'lev'
          VarZDimId   = zDimId
          TmpAxis     = 'Z'
          TmpPositive = 'up'

       ! time
       CASE( 'GRID_TIME' )
          VarName     = 'time'
          VarTDimId   = tDimId
          TmpAxis     = 'T'
          TmpCalendar = 'gregorian'

       ! area
       CASE( 'GRID_AREA' )
          VarName     = 'AREA'
          VarXDimId   = xDimId
          VarYDimId   = yDimId

       ! All other variable names
       CASE DEFAULT
          VarName = Item%Name

          SELECT CASE( Item%SpaceDim )
             CASE( 3 )
                VarXDimId = xDimId
                VarYDimId = yDimId
                VarZDimId = ZDimId
                VartDimId = tDimId
             CASE( 2 )
                VarXDimId = xDimId
                VarYDimId = yDimId
                VartDimId = tDimId
             CASE( 1 )
                ! NOTE: Need a way to resolve X Y Z
                VarXDimId = xDimId
                VartDimId = tDimId
             CASE DEFAULT
                ! Nothing
          END SELECT

    END SELECT      
    
    ! Return optional attributes for index variables: axis and calendar
    IF ( PRESENT( VarAxis     ) ) VarAxis     = TmpAxis
    IF ( PRESENT( VarCalendar ) ) VarCalendar = TmpCalendar
    IF ( PRESENT( VarPositive ) ) VarPositive = TmpPositive

  END SUBROUTINE Get_Var_DimIds
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
    CHARACTER(LEN=255)       :: Description
    CHARACTER(LEN=255)       :: ErrMsg
    CHARACTER(LEN=255)       :: ThisLoc
    CHARACTER(LEN=255)       :: Units
    CHARACTER(LEN=9)         :: ItemName(5)

    ! Pointer arrays
    REAL(fp),        POINTER :: Ptr1d  (:    )
    REAL(fp),        POINTER :: Ptr2d  (:,:  )

    ! Objects
    TYPE(HistItem),  POINTER :: Item

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
                ' -> History_Netcdf_Init (in History/history_mod.F90)'
    ItemName    = (/ 'GRID_AREA', 'GRID_TIME', 'GRID_LEV ',                  &
                     'GRID_LAT ', 'GRID_LON '                              /)
    Ptr1d       => NULL()
    Ptr2d       => NULL()

    !=======================================================================
    ! Create a HISTORY ITEM for each of the index fields (lon, lat, area)
    ! of gc_grid_mod.F90 and add them to a METAHISTORY ITEM list
    !=======================================================================
    DO N = 1, SIZE( ItemName )

       !---------------------------------------------------------------------
       ! Look up one of the index fields from gc_grid_mod.F90
       !---------------------------------------------------------------------
       CALL Lookup_Grid( am_I_Root   = am_I_Root,                            &
                         Variable    = TRIM( ItemName(N) ),                  &
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
          ErrMsg = 'Error in "Lookup_State_Chm" for diagnostic ' //          &
                   TRIM( ItemName(N) )
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
    USE MetaHistItem_Mod, ONLY: MetaHistItem_Destroy
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

