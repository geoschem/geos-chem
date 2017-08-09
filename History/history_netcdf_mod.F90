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
  
  IMPLICIT NONE
  PRIVATE

# include "netcdf.inc"
!
! !PUBLIC MEMBER FUNCTIONS
!
  PUBLIC :: HIstory_Netcdf_Define
!
! PRIVATE MEMBER FUNCTIONS:
!
!
! !REMARKS:
!
! !REVISION HISTORY:
!  06 Jan 2015 - R. Yantosca - Initial version
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
    INTEGER                     :: N,        VarCt
    INTEGER                     :: Year,     Month,      Day
    INTEGER                     :: Hour,     Minute,     Second

    ! Strings                   
    CHARACTER(LEN=2)            :: MonthStr, DayStr 
    CHARACTER(LEN=2)            :: HourStr,  MinuteStr,  SecondStr
    CHARACTER(LEN=4)            :: YearStr
    CHARACTER(LEN=5)            :: Z
    CHARACTER(LEN=8)            :: D
    CHARACTER(LEN=10)           :: T
    CHARACTER(LEN=255)          :: FileName, TimeString
    CHARACTER(LEN=255)          :: ErrMsg,   ThisLoc

    ! Arrays
    INTEGER                     :: V(8)

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

       ! Create the file and add global attributes
       ! Remain in netCDF define mode upon exiting this routine
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
         
!       ! Handle potential error (and exit define mode)
!       IF ( RC /= GC_SUCCESS ) THEN
!          ErrMsg = 'Error in Nc_Create!'
!          CALL Gc_Error( ErrMsg, Rc, ThisLoc )
!          CALL Nc_Set_DefMode( Container%FileId, Off=.TRUE. )
!          RETURN
!       ENDIF

       ! Denote that the file has been created and is open
       Container%IsFileOpen = .TRUE.
    ENDIF

    !=======================================================================
    ! Define all of the Create the netCDF file with global attributes
    !=======================================================================
    IF ( .not. Container%IsFileDefined ) THEN

       !--------------------------------------------------------------------
       ! Add the index dimension data
       !--------------------------------------------------------------------


       !--------------------------------------------------------------------
       ! Then define each HISTORY ITEM belonging to this collection
       !--------------------------------------------------------------------

       ! Set CURRENT to the first node in the list of HISTORY ITEMS
       Current => Container%HistItems

       ! As long as this node of the list is valid ...
       DO WHILE( ASSOCIATED( Current ) )

          ! Define each HISTORY ITEM in this collection to the netCDF file
          CALL Nc_Var_Def( fId          = Container%FileId,                  &
                           timeId       = Container%tDimId,                  &
                           levId        = Container%zDimId,                  &
                           latId        = Container%yDimId,                  &
                           lonId        = Container%xDimId,                  &
                           VarName      = Current%Item%Name,                 &
                           VarLongName  = Current%Item%LongName,             &
                           VarUnit      = Current%Item%Units,                &
                           AddOffset    = Current%Item%AddOffset,            &
                           MissingValue = Current%Item%MissingValue,         &
                           ScaleFactor  = Current%Item%ScaleFactor,          &
                           DataType     = 4,                                 &
                           VarCt        = VarCt,                             &
                           DefMode      = .TRUE.,                            &
                           Compress     = .TRUE.                            ) 

          ! Save the returned netCDF variable Id to this HISTORY ITEM
          Current%Item%NcVarId = VarCt

          ! Debug print
          !CALL HistItem_Print( am_I_Root, Current%Item, RC )

          ! Go to next entry in the list of HISTORY ITEMS
          Current => Current%Next
       ENDDO

       ! Free pointers
       Current => NULL()
    ENDIF

    !=======================================================================
    ! Cleanup and quit
    !=======================================================================

    ! Close definition section
    CALL Nc_Set_DefMode( Container%FileId, Off=.TRUE. )

    ! We can now consider this collection to have been "defined"
    Container%IsFileDefined = .TRUE.

    !CALL HistContainer_Print( am_I_Root, Container, RC )

    !STOP

  END SUBROUTINE History_NetCdf_Define
!EOC
END MODULE History_Netcdf_Mod

