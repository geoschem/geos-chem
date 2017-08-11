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
  PRIVATE :: History_Expand_Date
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
    VarAxis     =  ''
    VarCalendar =  ''
    VarPositive =  ''
    VarUnits    =  ''
    
    !=======================================================================
    ! Construct the file name using the passed date and time
    !=======================================================================

    ! Save the collection's file name in a temporary variable
    FileName = TRIM( Container%FileName )

    ! Replace date and time tokens in the file name
    CALL History_Expand_Date( FileName, yyyymmdd, hhmmss, MAPL_Style=.TRUE. )
   
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
    ! Create the timestamp for the History and ProdDateTime attributes
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

       ! Close the netCDF file
       CALL Nc_Close( Container%FileId )

       ! Set fields to denotet his file is closed
       Container%IsFileOpen    = .FALSE.
       Container%IsFileDefined = .FALSE.

       ! Set CURRENT to the first node in the list of 
       ! HISTORY ITEMS belonging to this collection
       Current => Container%HistItems

       ! As long as this node of the list is valid ...
       DO WHILE( ASSOCIATED( Current ) )

          ! Undefine quantities for the file we just closed
          Container%ReferenceYmd = 0
          Container%ReferenceHms = 0
          Current%Item%NcXDimId  = UNDEFINED_INT
          Current%Item%NcYDimId  = UNDEFINED_INT
          Current%Item%NcZDimId  = UNDEFINED_INT
          Current%Item%NcTDimId  = UNDEFINED_INT
          Current%Item%NcVarId   = UNDEFINED_INT
          
          ! Go to the next HISTORY ITEM
          Current => Current%Next
       ENDDO

       ! Free pointer
       Current => NULL()
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

       ! Set the reference time
       Container%ReferenceYmd = yyyymmdd
       Container%ReferenceHms = hhmmss

       ! Create the file
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
       ! Also set the reference time and date accordingly
       !--------------------------------------------------------------------
       Container%IsFileOpen   = .TRUE.
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
             CALL History_Expand_Date( TmpUnits, ReferenceYmd, ReferenceHms )
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
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: History_Expand_Date
!
! !DESCRIPTION: Replaces date and time tokens in a string with actual 
!  date and time values.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE History_Expand_Date( DateStr, yyyymmdd, hhmmss, MAPL_Style )
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

  END SUBROUTINE History_Expand_Date
!EOC
END MODULE History_Netcdf_Mod

