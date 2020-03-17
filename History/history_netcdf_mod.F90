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
!
! !PRIVATE MEMBER FUNCTIONS:
!
  PRIVATE :: Expand_Date_Time
  PRIVATE :: Get_Number_Of_Levels
  PRIVATE :: Get_Var_DimIds
  PRIVATE :: IndexVarList_Create
  PRIVATE :: IndexVarList_Destroy
!
! !REMARKS:
!
! !REVISION HISTORY:
!  10 Aug 2017 - R. Yantosca - Initial version
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
  SUBROUTINE History_Netcdf_Close( Container, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE HistContainer_Mod,     ONLY : HistContainer
    USE History_Util_Mod
    USE MetaHistContainer_Mod, ONLY : MetaHistContainer
    USE Ncdf_Mod
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
!  See https://github.com/geoschem/geos-chem for complete history
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
       Container%ReferenceYmd  = UNDEFINED_INT
       Container%ReferenceHms  = UNDEFINED_INT
!       Container%ReferenceJd   = UNDEFINED_DBL
       Container%ReferenceJsec = UNDEFINED_DBL
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
  SUBROUTINE History_Netcdf_Define( Input_Opt, Container, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE HistContainer_Mod,   ONLY : HistContainer, HistContainer_Print
    USE HistItem_Mod,        ONLY : HistItem,      HistItem_Print
    USE History_Util_Mod
    USE Input_Opt_Mod,       ONLY : OptInput
    USE JulDay_Mod,          ONLY : CalDate
    USE MetaHistItem_Mod,    ONLY : MetaHistItem
    USE Ncdf_Mod
    USE Registry_Params_Mod, ONLY : KINDVAL_F4
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput),      INTENT(IN)  :: Input_Opt   ! Input options
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(HistContainer), POINTER     :: Container  ! Diagnostic collection obj
!
! !OUTPUT PARAMETERS:
!
    INTEGER,             INTENT(OUT) :: RC         ! Success or failure
!
! !REMARKS:
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
    INTEGER                     :: N,          VarCt
    INTEGER                     :: VarXDimId,  VarYDimId
    INTEGER                     :: VarZDimId,  VarTDimId
    INTEGER                     :: yyyymmdd,   hhmmss
    INTEGER                     :: nLev,       nILev
    INTEGER                     :: DataType

    ! Strings
    CHARACTER(LEN=5)            :: Z
    CHARACTER(LEN=8)            :: D
    CHARACTER(LEN=10)           :: T
    CHARACTER(LEN=255)          :: FileName
    CHARACTER(LEN=255)          :: ErrMsg,     ThisLoc,     VarUnits
    CHARACTER(LEN=255)          :: VarAxis,    VarPositive, VarCalendar
    CHARACTER(LEN=255)          :: VarStdName, VarFormula

    ! Arrays
    INTEGER                     :: V(8)

    ! Pointers
    REAL(fp),           POINTER :: Data1d(:)
    REAL(fp),           POINTER :: Data2d(:)

    ! Objects
    TYPE(MetaHistItem), POINTER :: Current
    TYPE(MetaHistItem), POINTER :: IndexVarList

    !=======================================================================
    ! Initialize
    !=======================================================================
    RC          =  GC_SUCCESS
    Current     => NULL()
    IndexVarList => NULL()
    ErrMsg      =  ''
    ThisLoc     =  ' -> at History_Netcdf_Define (in History/history_mod.F90)'
    FileName    =  ''
    VarAxis     =  ''
    VarCalendar =  ''
    VarPositive =  ''
    VarUnits    =  ''
    yyyymmdd    =  Container%CurrentYmd
    hhmmss      =  Container%CurrentHms

    !=======================================================================
    ! Create the netCDF file with global attributes
    ! Do not exit netCDF define mode just yet
    !=======================================================================
    IF ( .not. Container%IsFileOpen ) THEN

       !--------------------------------------------------------------------
       ! Compute reference date and time fields in the HISTORY CONTAINER
       ! These are needed to compute the time stamps for each data field
       ! that is written to the netCDF file.  Also resets the current
       ! time slice index.
       !--------------------------------------------------------------------

       ! Reset the current time slice index
       Container%CurrTimeSlice = 0

       IF ( Container%Operation == ACCUM_FROM_SOURCE ) THEN

          !-----------------------------
          ! TIME-AVERAGED COLLECTIONS
          !-----------------------------

          ! Subtract the file write alarm interval that we added to
          ! the current date/time (CurrentJd) field at initialization
          Container%ReferenceJsec = Container%CurrentJsec                    &
                                  - Container%FileWriteIvalSec

       ELSE

          !-----------------------------
          ! INSTANTANEOUS COLLECTIONS
          !-----------------------------

          ! If this is the first time we are writing a file, then set
          ! the reference date/time to the date/time at the start of
          ! the simulation (EpochJd).  This will make sure the file names
          ! and timestamps in all instantaneous files will be consistent.
          ! For all future file writes, set the reference date/time to the
          ! current date/time (CurrentJd).
          IF ( Container%FirstInst ) THEN
             Container%ReferenceJsec = Container%EpochJsec
             Container%FirstInst   = .FALSE.
          ELSE
             Container%ReferenceJsec = Container%CurrentJsec
          ENDIF
       ENDIF

       ! Convert the reference time from Astronomical Julian Seconds
       ! to Astronomical Julian Date.  This is needed for the conversion
       ! to calendar date and time via routine CALDATE.
       Container%ReferenceJd = Container%ReferenceJsec / SECONDS_PER_DAY

       ! Recompute the ReferenceYmd and ReferenceHms fields
       CALL CalDate( JulianDay = Container%ReferenceJd,                      &
                     yyyymmdd  = Container%ReferenceYmd,                     &
                     hhmmss    = Container%ReferenceHms                     )

       !--------------------------------------------------------------------
       ! Replace time and date tokens in the netCDF file name
       !--------------------------------------------------------------------

       ! Save the collection's file name in a temporary variable
       FileName = TRIM( Container%FileName )

       ! Replace date and time tokens in the file name
       CALL Expand_Date_Time( DateStr    = FileName,                         &
                              yyyymmdd   = Container%ReferenceYmd,           &
                              hhmmss     = Container%ReferenceHms,           &
                              MAPL_Style = .TRUE. )

!------------------------------------------------------------------------------
! TEMPORARY FIX (bmy, 9/20/17)
! NOTE: The different timestamps will cause the binary diff in the unit
! tests and difference tests to fail, so comment these out for now.
! We will look into a better way to check netCDF files soon.
!       !--------------------------------------------------------------------
!       ! Create the timestamp for the History and ProdDateTime attributes
!       !--------------------------------------------------------------------
!
!       ! Call F90 intrinsic DATE_AND_TIME Function
!       D = 'ccyymmdd'
!       T = 'hhmmss.sss'
!       CALL Date_And_Time( Date=D, Time=T, Zone=Z, Values=V )  ! GMT time
!
!       ! Create timestamp strings
!       WRITE( Container%History,      10 ) V(1),V(2),V(3),V(5),V(6),V(7),Z
!       WRITE( Container%ProdDateTime, 10 ) V(1),V(2),V(3),V(5),V(6),V(7),Z
! 10    FORMAT( 'Produced on ', i4.4, '/', i2.2, '/', i2.2, 1x,               &
!                               i2.2, ':', i2.2, ':', i2.2, ' UTC', a        )
!
       ! For now, just set History and ProdDateTime to blanks
       ! to get binary file diffs to pass.
       Container%History      = ''
       Container%ProdDateTime = ''
!------------------------------------------------------------------------------

       ! Get the number of levels (nLev) and level interfaces (nIlev)
       CALL Get_Number_Of_Levels( Container, nLev, nIlev )

       ! Do not create netCDF file on first timestep if instantaneous collection
       ! and frequency = duration. This will avoid creation of a netCDF file
       ! containing all missing values.
       IF ( TRIM( Container%UpdateMode ) == 'instantaneous'       .and. &
            Container%UpdateIvalSec == Container%FileCloseIvalSec .and. &
            Container%FileCloseAlarm  == 0.0  ) THEN
          RETURN
       ELSE

          ! Echo info about the file we are creating
          IF ( Input_Opt%amIRoot ) THEN
             WRITE( 6, 100 ) TRIM( Container%Name ),                         &
                             Container%ReferenceYmd,                         &
                             Container%ReferenceHms
             WRITE( 6, 110 ) TRIM( FileName       )
          ENDIF
100       FORMAT( '     - Creating file for ', a, '; reference = ',i8.8,1x,i6.6)
110       FORMAT( '        with filename = ', a                                )

          !--------------------------------------------------------------------
          ! Create the file and add global attributes
          ! Remain in netCDF define mode upon exiting this routine
          !
          ! NOTE: Container%Reference is a global attribute lists the GEOS-Chem
          ! web and wiki page.  It has nothing to do with the reference date
          ! and time fields that are computed by History_Set_RefDateTime.
          !--------------------------------------------------------------------
          CALL Nc_Create( Create_Nc4     = .TRUE.,                           &
                          NcFile         = FileName,                         &
                          nLon           = Container%nX,                     &
                          nLat           = Container%nY,                     &
                          nLev           = nLev,                             &
                          nIlev          = nILev,                            &
                          nTime          = NF_UNLIMITED,                     &
                          NcFormat       = Container%NcFormat,               &
                          Conventions    = Container%Conventions,            &
                          History        = Container%History,                &
                          ProdDateTime   = Container%ProdDateTime,           &
                          Reference      = Container%Reference,              &
                          Title          = Container%Title,                  &
                          Contact        = Container%Contact,                &
                          StartTimeStamp = Container%StartTimeStamp,         &
                          EndTimeStamp   = Container%EndTimeStamp,           &
                          fId            = Container%FileId,                 &
                          TimeId         = Container%tDimId,                 &
                          LevId          = Container%zDimId,                 &
                          ILevId         = Container%iDimId,                 &
                          LatId          = Container%yDimId,                 &
                          LonId          = Container%xDimId,                 &
                          KeepDefMode    = .TRUE.,                           &
                          Varct          = VarCt                               )

          !--------------------------------------------------------------------
          ! Denote that the file has been created and is open
          !--------------------------------------------------------------------
          Container%IsFileOpen = .TRUE.

       ENDIF

    ENDIF

    !=======================================================================
    ! Define all of the Create the netCDF file with global attributes
    !=======================================================================
    IF ( .not. Container%IsFileDefined ) THEN

       !--------------------------------------------------------------------
       ! Define the index variables
       !--------------------------------------------------------------------

       ! Define the linked list of index variables (IndexVarList) that
       ! have the same dimension subsets as the current container
       CALL IndexVarList_Create( Input_Opt, Container, IndexVarList, RC )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Error encountered in "History_NetCdf_Init"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       ! Set CURRENT to the first node in the list of HISTORY ITEMS
       Current => IndexVarList

       ! As long as this node of the list is valid ...
       DO WHILE( ASSOCIATED( Current ) )

          ! Get the dimension ID's that are relevant to each HISTORY ITEM
          ! Also get Axis, Calendar, and Positive attributes for index vars,
          ! and replace time & date tokens in the units string for "time".
          CALL Get_Var_DimIds( xDimId       = Container%xDimId,              &
                               yDimId       = Container%yDimId,              &
                               zDimId       = Container%zDimId,              &
                               iDimId       = Container%iDimId,              &
                               tDimId       = Container%tDimId,              &
                               RefDate      = Container%ReferenceYmd,        &
                               RefTime      = Container%ReferenceHms,        &
                               OnLevelEdges = Container%OnLevelEdges,        &
                               Item         = Current%Item,                  &
                               VarAxis      = VarAxis,                       &
                               VarPositive  = VarPositive,                   &
                               VarCalendar  = VarCalendar,                   &
                               VarUnits     = VarUnits,                      &
                               VarStdName   = VarStdName,                    &
                               VarFormula   = VarFormula                    )

          ! Set a flag for the precision of the data
          IF ( Current%Item%Source_KindVal == KINDVAL_F4 ) THEN
             DataType = 4
          ELSE
             DataType = 8
          ENDIF

          ! Define each HISTORY ITEM in this collection to the netCDF file
          CALL Nc_Var_Def( DefMode      = .TRUE.,                            &
                           Compress     = .TRUE.,                            &
                           fId          = Container%FileId,                  &
                           DataType     = DataType,                          &
                           VarName      = Current%Item%Name,                 &
                           VarCt        = Current%Item%NcVarId,              &
                           timeId       = Current%Item%NcTDimId,             &
                           levId        = Current%Item%NcZDimId,             &
                           iLevId       = Current%Item%NcIDimId,             &
                           latId        = Current%Item%NcYDimId,             &
                           lonId        = Current%Item%NcXDimId,             &
                           VarLongName  = Current%Item%LongName,             &
                           VarUnit      = VarUnits,                          &
                           Axis         = VarAxis,                           &
                           Calendar     = VarCalendar,                       &
                           Positive     = VarPositive,                       &
                           StandardName = VarStdName,                        &
                           FormulaTerms = VarFormula                        )

          ! Debug print
          !CALL HistItem_Print( Current%Item, RC )

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
          CALL Get_Var_DimIds( xDimId   = Container%xDimId,                  &
                               yDimId   = Container%yDimId,                  &
                               zDimId   = Container%zDimId,                  &
                               iDimId   = Container%iDimId,                  &
                               tDimId   = Container%tDimId,                  &
                               Item     = Current%Item,                      &
                               VarUnits = VarUnits                          )

          ! Replace "TBD"  with the current units of State_Chm%Species
          IF ( TRIM( VarUnits ) == 'TBD' ) THEN
             VarUnits = Container%Spc_Units
          ENDIF

          ! Define each HISTORY ITEM in this collection
          ! as a variable in the netCDF file
          CALL Nc_Var_Def( DefMode      = .TRUE.,                            &
                           Compress     = .TRUE.,                            &
                           fId          = Container%FileId,                  &
                           DataType     = 4,                                 &
                           VarName      = Current%Item%Name,                 &
                           VarCt        = Current%Item%NcVarId,              &
                           timeId       = Current%Item%NcTDimId,             &
                           levId        = Current%Item%NcZDimId,             &
                           iLevId       = Current%Item%NcIDimId,             &
                           latId        = Current%Item%NcYDimId,             &
                           lonId        = Current%Item%NcXDimId,             &
                           VarLongName  = Current%Item%LongName,             &
                           VarUnit      = VarUnits,                          &
                           MissingValue = Current%Item%MissingValue,         &
                           AvgMethod    = Current%Item%AvgMethod            )

#if defined( NC_HAS_COMPRESSION )
          ! Turn on netCDF chunking for this HISTORY ITEM
          ! NOTE: This will only work if the netCDF library supports netCDF-4
          ! files with compression.  Also note: file compression is turned off
          ! by default when using DEBUG=y, because otherwise the compression
          ! makes it difficult to compare files generated by difference tests.
          IF ( ASSOCIATED( Current%Item%NcChunkSizes ) ) THEN

             ! Apply the chunk sizes to this variable
             CALL Nc_Var_Chunk( fId        = Container%FileId,               &
                                vId        = Current%Item%NcVarId,           &
                                ChunkSizes = Current%Item%NcChunkSizes,      &
                                RC         = RC                             )

             ! Trap potential error
             IF ( RC /= GC_SUCCESS ) THEN
                ErrMsg = 'Error encountered in "Nc_Var_Chunk"!'           // &
                         ' Try DEBUG=n and/or check if your netCDF-4 '    // &
                         ' library supports compression'
                CALL GC_Error( ErrMsg, RC, ThisLoc )
                RETURN
             ENDIF

          ENDIF
#endif

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

          ! Write data for index variables to the netCDF file
          IF ( Current%Item%SpaceDim == 2 ) THEN

             ! AREA is the only 2-D array (4-byte precison)
             CALL Nc_Var_Write( fId     = Container%FileId,                  &
                                VarName = Current%Item%Name,                 &
                                Arr2d   = Current%Item%Source_2d_4          )

          ELSE IF ( Current%Item%SpaceDim == 1 ) THEN

             ! All other index fields are 1-D (8-byte precision) ...
             CALL Nc_Var_Write( fId     = Container%FileId,                  &
                                VarName = Current%Item%Name,                 &
                                Arr1d   = Current%Item%Source_1d_8          )


          ELSE

             ! ... except P0, which is a scalar (8-byte precision)
             CALL Nc_Var_Write( fId     = Container%FileId,                  &
                                VarName = Current%Item%Name,                 &
                                Var     = Current%Item%Source_0d_8          )

          ENDIF

          ! Go to next entry in the list of HISTORY ITEMS
          Current => Current%Next
       ENDDO

       ! Free pointers
       Current => NULL()

       !--------------------------------------------------------------------
       ! We can now consider this collection to have been "defined"
       !--------------------------------------------------------------------
       Container%IsFileDefined = .TRUE.

       CALL IndexVarList_Destroy( IndexVarList, RC )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Error returned from "History_Netcdf_Cleanup"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

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
  SUBROUTINE History_Netcdf_Write( Input_Opt, Container, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE HistItem_Mod,       ONLY : HistItem
    USE HistContainer_Mod,  ONLY : HistContainer
    USE History_Util_Mod
    USE Input_Opt_Mod,      ONLY : OptInput
    USE M_Netcdf_Io_Write,  ONLY : NcWr
    USE MetaHistItem_Mod,   ONLY : MetaHistItem
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput),      INTENT(IN)  :: Input_Opt ! Input Options object
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
!  Use routine NcWr from NcdfUtil/m_netcdf_io_write.F90 instead of the
!  NC_VAR_WRITE routine from NcdfUtil/netcdf_mod.F90, because this gives us
!  better control of the start and count values.
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
    INTEGER                     :: NcFileId,         NcVarId
    INTEGER                     :: Dim1,             Dim2,       Dim3

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

    !=======================================================================
    ! Compute the time stamp value for the current time slice
    !=======================================================================

    ! Compute the elapsed time in seconds since the file creation
    CALL Compute_Elapsed_Time( CurrentJsec  = Container%CurrentJsec,         &
                               TimeBaseJsec = Container%ReferenceJsec,       &
                               ElapsedSec   = Container%TimeStamp           )

    ! For time-averaged collections, offset the timestamp
    ! by 1/2 of the file averaging interval in minutes
    IF ( Container%Operation == ACCUM_FROM_SOURCE ) THEN
       Container%TimeStamp = Container%TimeStamp -                           &
                             ( Container%FileWriteIvalSec * 0.5_f8 )
    ENDIF

    ! Convert to minutes since the reference time
    Container%TimeStamp = Container%TimeStamp / SECONDS_PER_MINUTE

    ! Debug output
    IF ( Input_Opt%LPRT .and. Input_opt%amIRoot ) THEN
       WRITE( 6, 110 ) TRIM( Container%name ), Container%TimeStamp
110    FORMAT( '     - Writing data to ', a, '; timestamp = ', f13.4 )
    ENDIF

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
                Item%Data_3d  = Item%Data_3d / Item%nUpdates
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
             Dim1 = SIZE( Item%Data_2d, 1 )
             Dim2 = SIZE( Item%Data_2d, 2 )

             ! Allocate the REAL*4 output array
             ALLOCATE( NcData_2d( Dim1, Dim2 ), STAT=RC )

             ! Copy or average the data and store in a REAL*4 array
             IF ( Item%Operation == COPY_FROM_SOURCE ) THEN
                NcData_2d     = Item%Data_2d
                Item%Data_2d  = 0.0_f8
                Item%nUpdates = 0
             ELSE
                Item%Data_2d  = Item%Data_2d / Item%nUpdates
                NcData_2d     = Item%Data_2d
                Item%Data_2d  = 0.0_f8
                Item%nUpdates = 0
             ENDIF

             ! Compute start and count fields
             St3d = (/ 1,    1,    Container%CurrTimeSlice /)
             Ct3d = (/ Dim1, Dim2, 1                       /)

             ! Write data to disk
             CALL NcWr( NcData_2d, NcFileId, Item%Name, St3d, Ct3d )

             ! Deallocate output array
             DEALLOCATE( NcData_2d, STAT=RC )

          !------------
          ! 1-D data
          !------------
          CASE( 1 )

             ! Get dimensions of data
             Dim1 = SIZE( Item%Data_1d, 1 )

             ! Allocate the REAL*4 output array
             ALLOCATE( NcData_1d( Dim1 ), STAT=RC )

             ! Copy or average the data and store in a REAL*4 array
             IF ( Item%Operation == COPY_FROM_SOURCE ) THEN
                NcData_1d     = Item%Data_1d
                Item%Data_1d  = 0.0_f8
                Item%nUpdates = 0
             ELSE
                Item%Data_1d  = Item%Data_1d / Item%nUpdates
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
             DEALLOCATE( NcData_1d, STAT=RC )

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
!  See https://github.com/geoschem/geos-chem for complete history
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
  SUBROUTINE Get_Var_DimIds( Item,     xDimId,      yDimId,                  &
                             zDimId,   iDimId,      tDimID,                  &
                             RefDate,  RefTime,     OnLevelEdges,            &
                             VarAxis,  VarCalendar, VarPositive,             &
                             VarUnits, VarStdName,  VarFormula              )
!
! !USES:
!
    USE History_Util_Mod
    USE HistItem_Mod,       ONLY : HistItem
!
! !INPUT PARAMETERS:
!
    INTEGER,            INTENT(IN)  :: xDimId       ! Id # of X (lon     ) dim
    INTEGER,            INTENT(IN)  :: yDimId       ! Id # of Y (lat     ) dim
    INTEGER,            INTENT(IN)  :: zDimId       ! Id # of Z (lev cntr) dim
    INTEGER,            INTENT(IN)  :: iDimId       ! Id # of I (lev edge) dim
    INTEGER,            INTENT(IN)  :: tDimId       ! Id # of T (time    ) dim
    INTEGER,            OPTIONAL    :: RefDate      ! Ref YMD for "time" var
    INTEGER,            OPTIONAL    :: RefTime      ! Ref hms for "time" var
    LOGICAL,            OPTIONAL    :: OnLevelEdges ! Is 3D data on lvl edges?
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(HistItem),     POINTER     :: Item         ! HISTORY ITEM object
!
! !OUTPUT PARAMETERS
!
    CHARACTER(LEN=255), OPTIONAL    :: VarAxis      ! Axis attr for index vars
    CHARACTER(LEN=255), OPTIONAL    :: VarCalendar  ! Calendar attr for "time"
    CHARACTER(LEN=255), OPTIONAL    :: VarPositive  ! Positive attr for "lev"
    CHARACTER(LEN=255), OPTIONAL    :: VarUnits     ! Unit string
    CHARACTER(LEN=255), OPTIONAL    :: VarStdName   ! Standard name
    CHARACTER(LEN=255), OPTIONAL    :: VarFormula   ! Formula terms

!
! !REMARKS:
!  Call this routine before calling NC_VAR_DEF.
!
! !REVISION HISTORY:
!  10 Aug 2017 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    LOGICAL            :: IsOnLevelEdges
    INTEGER            :: ReferenceYmd,   ReferenceHms
    INTEGER            :: Year,           Month,        Day
    INTEGER            :: Hour,           Minute,       Second

    ! Strings
    CHARACTER(LEN=2)   :: MonthStr,       DayStr
    CHARACTER(LEN=2)   :: HourStr,        MinuteStr,    SecondStr
    CHARACTER(LEN=4)   :: YearStr
    CHARACTER(LEN=255) :: TmpAxis,        TmpCalendar,  TmpStdName
    CHARACTER(LEN=255) :: TmpPositive,    TmpUnits,     TmpFormula

    !=======================================================================
    ! Initialize
    !=======================================================================
    TmpAxis         = ''
    TmpCalendar     = ''
    TmpPositive     = ''
    TmpUnits        = Item%Units
    TmpStdName      = ''
    TmpFormula      = ''
    Item%NcXDimId   = UNDEFINED_INT
    Item%NcYDimId   = UNDEFINED_INT
    Item%NcZDimId   = UNDEFINED_INT
    Item%NcIDimId   = UNDEFINED_INT
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

    IF ( PRESENT( OnLevelEdges ) ) THEN
       IsOnLevelEdges = OnLevelEdges
    ELSE
       IsOnLevelEdges = .FALSE.
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
          TmpAxis       = ''
          TmpPositive   = 'up'
          TmpStdName    = 'atmosphere_hybrid_sigma_pressure_coordinate'
          TmpFormula    = 'a: hyam b: hybm p0: P0 ps: PS'

          ! If the collection contains is level-centered data
          ! then "lev" (and not "ilev") is the "Z" axis
          IF ( .not. IsOnLevelEdges ) THEN
             TmpAxis    = 'Z'
          ENDIF

       ! ilev
       CASE( 'ilev' )
          Item%NcZDimId = iDimId
          TmpAxis       = ''
          TmpPositive   = 'up'
          TmpStdName    = 'atmosphere_hybrid_sigma_pressure_coordinate'
          TmpFormula    = 'a: hyai b: hybi p0: P0 ps: PS'

          ! If the collection contains is level-centered data
          ! then "ilev" (and not "lev") is the "Z" axis
          IF ( IsOnLevelEdges ) THEN
             TmpAxis    = 'Z'
          ENDIF

       ! hybrid coordinates, level centers
       CASE( 'hyam', 'hybm' )
          Item%NcZDimId = zDimId

       ! hybrid coordinates, level edges
       CASE( 'hyai', 'hybi' )
          Item%NcZDimId = iDimId

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

          ! Set the various netCDF dimension variables that will be passed
          ! to NC_CREATE.  If the data is defined on vertical level edges
          ! (aka "interfaces), then use iDimId instead of zDimId.
          SELECT CASE( TRIM( Item%DimNames ) )

             CASE( 'xyz' )
                Item%NcXDimId = xDimId
                Item%NcYDimId = yDimId
                Item%NcTDimId = tDimId

                IF ( Item%OnLevelEdges ) THEN
                   Item%NcIDimId = iDimId
                ELSE
                   Item%NcZDimId = zDimId
                ENDIF

             CASE( 'xy'  )
                Item%NcXDimId = xDimId
                Item%NcYDimId = yDimId
                Item%NcTDimId = tDimId

             CASE( 'yz'  )
                Item%NcYDimId = yDimId
                Item%NcTDimId = tDimId

                IF ( Item%OnLevelEdges ) THEN
                   Item%NcIDimId = iDimId
                ELSE
                   Item%NcZDimId = zDimId
                ENDIF

             CASE( 'xz'  )
                Item%NcXDimId = xDimId
                Item%NcTDimId = tDimId

                IF ( Item%OnLevelEdges ) THEN
                   Item%NcIDimId = iDimId
                ELSE
                   Item%NcZDimId = zDimId
                ENDIF

             CASE( 'x'   )
                Item%NcXDimId = xDimId
                Item%NcTDimId = tDimId

             CASE( 'y'   )
                Item%NcYDimId = yDimId
                Item%NcTDimId = tDimId

             CASE( 'z'   )
                Item%NcTDimId = tDimId

                IF ( Item%OnLevelEdges ) THEN
                   Item%NcIDimId = iDimId
                ELSE
                   Item%NcZDimId = zDimId
                ENDIF

             CASE DEFAULT
                ! Nothing

          END SELECT

    END SELECT

    ! Return optional attributes for index variables: axis and calendar
    IF ( PRESENT( VarAxis     ) ) VarAxis     = TmpAxis
    IF ( PRESENT( VarCalendar ) ) VarCalendar = TmpCalendar
    IF ( PRESENT( VarPositive ) ) VarPositive = TmpPositive
    IF ( PRESENT( VarUnits    ) ) VarUnits    = TmpUnits
    IF ( PRESENT( VarStdName  ) ) VarStdName  = TmpStdName
    IF ( PRESENT( VarFormula  ) ) VarFormula  = TmpFormula

  END SUBROUTINE Get_Var_DimIds
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: IndexVarList_Create
!
! !DESCRIPTION: Creates a HISTORY ITEM for each netCDF index variable (e.g.
!  lon, lat, lev, time, area) and adds it to the METAHISTORY ITEM IndexVarList.
!  Subsets each index variable according to the subset indices from the
!  given collection (passed via the Container argument).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE IndexVarList_Create( Input_Opt, Container, IndexVarList, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE Grid_Registry_Mod, ONLY : Lookup_Grid
    USE HistContainer_Mod, ONLY : HistContainer
    USE HistItem_Mod
    USE Input_Opt_Mod,     ONLY : OptInput
    USE MetaHistItem_Mod
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput),      INTENT(IN)  :: Input_Opt    ! Input Options object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(HistContainer), POINTER     :: Container    ! Collection object
!
! !OUTPUT PARAMETERS:
!
    TYPE(MetaHistItem),  POINTER     :: IndexVarList ! Linked list of index
                                                     !  variables for netCDF
    INTEGER,             INTENT(OUT) :: RC           ! Success or failure?
!
! !REMARKS:
!
! !REVISION HISTORY:
!  10 Aug 2017 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    LOGICAL                  :: OnLevelEdges
    INTEGER                  :: N
    INTEGER                  :: KindVal
    INTEGER                  :: Rank
    INTEGER                  :: nILev
    INTEGER                  :: nLev

    ! Arrays
    INTEGER                  :: Dimensions(3)
    INTEGER                  :: Subset_X(2)
    INTEGER                  :: Subset_Y(2)
    INTEGER                  :: Subset_Z(2)
    INTEGER                  :: Subset_Zc(2)
    INTEGER                  :: Subset_Ze(2)

    ! Strings
    CHARACTER(LEN=20)        :: ItemDimName(11)
    CHARACTER(LEN=20)        :: ItemName(11)
    CHARACTER(LEN=20)        :: RegistryName(11)
    CHARACTER(LEN=255)       :: Description
    CHARACTER(LEN=255)       :: ErrMsg
    CHARACTER(LEN=255)       :: ThisLoc
    CHARACTER(LEN=255)       :: Units

    ! Pointer arrays
    REAL(f8),        POINTER :: Ptr0d_8
    REAL(f8),        POINTER :: Ptr1d_8(:    )
    REAL(f4),        POINTER :: Ptr2d_4(:,:  )

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
    Ptr0d_8      => NULL()
    Ptr1d_8      => NULL()
    Ptr2d_4      => NULL()

    !=======================================================================
    ! Define the names that will be used to create the HISTORY ITEMS
    ! for fields to be used as netCDF metadata
    !=======================================================================

    ! Fields saved in the Registry object in GeosUtil/grid_registry_mod.F90
    RegistryName(1 ) = 'GRID_AREA'
    RegistryName(2 ) = 'GRID_P0'
    RegistryName(3 ) = 'GRID_HYBI'
    RegistryName(4 ) = 'GRID_HYAI'
    RegistryName(5 ) = 'GRID_HYBM'
    RegistryName(6 ) = 'GRID_HYAM'
    RegistryName(7 ) = 'GRID_LON'
    RegistryName(8 ) = 'GRID_LAT'
    RegistryName(9 ) = 'GRID_ILEV'
    RegistryName(10) = 'GRID_LEV'
    RegistryName(11) = 'GRID_TIME'

    ! Name for each HISTORY ITEM
    ItemName(1 )     = 'AREA'
    ItemName(2 )     = 'P0'
    ItemName(3 )     = 'hybi'
    ItemName(4 )     = 'hyai'
    ItemName(5 )     = 'hybm'
    ItemName(6 )     = 'hyam'
    ItemName(7 )     = 'lon'
    ItemName(8 )     = 'lat'
    ItemName(9 )     = 'ilev'
    ItemName(10)     = 'lev'
    ItemName(11)     = 'time'

    ! Dimensions for each HISTORY ITEM
    ItemDimName(1 )  = 'xy'
    ItemDimName(2 )  = '-'
    ItemDimName(3 )  = 'z'
    ItemDimName(4 )  = 'z'
    ItemDimName(5 )  = 'z'
    ItemDimName(6 )  = 'z'
    ItemDimName(7 )  = 'x'
    ItemDimName(8 )  = 'y'
    ItemDimName(9 )  = 'z'
    ItemDimName(10)  = 'z'
    ItemDimName(11)  = 't'

    !=======================================================================
    ! Pick the dimensions of the lev and ilev variables properly
    !=======================================================================

    ! Get the number of levels (nLev) and level interfaces (nIlev)
    CALL Get_Number_Of_Levels( Container, nLev, nIlev )

    ! Subset indices
    Subset_X  = (/ Container%X0, Container%X1 /)
    Subset_Y  = (/ Container%Y0, Container%Y1 /)
    Subset_Zc = (/ Container%Z0, nLev         /)
    Subset_Ze = (/ Container%Z0, nILev        /)

    !=======================================================================
    ! Create a HISTORY ITEM for each of the index fields (lon, lat, area)
    ! of grid_registry_mod.F90 and add them to a METAHISTORY ITEM list
    !=======================================================================
    DO N = 1, SIZE( RegistryName )

       !---------------------------------------------------------------------
       ! Look up one of the index fields from gc_grid_mod.F90
       !---------------------------------------------------------------------
       CALL Lookup_Grid( Input_Opt    = Input_Opt,                           &
                         Variable     = RegistryName(N),                     &
                         Description  = Description,                         &
                         Dimensions   = Dimensions,                          &
                         KindVal      = KindVal,                             &
                         Rank         = Rank,                                &
                         Units        = Units,                               &
                         OnLevelEdges = OnLevelEdges,                        &
                         Ptr0d_8      = Ptr0d_8,                             &
                         Ptr1d_8      = Ptr1d_8,                             &
                         Ptr2d_4      = Ptr2d_4,                             &
                         RC           = RC                                  )

       ! Trap potential error
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Error in "Lookup_Grid" for diagnostic ' //               &
                   TRIM( RegistryName(N) )
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       ! Pick proper subset indices for index variables placed on
       ! edges  (hyai, hybi, ilev) or centers (everything else)
       SELECT CASE( N )
          CASE( 3, 4, 9 )
             Subset_Z = Subset_Ze
          CASE DEFAULT
             Subset_Z = Subset_Zc
       END SELECT

       !---------------------------------------------------------------------
       ! Create a HISTORY ITEM for this index field
       !---------------------------------------------------------------------
       CALL HistItem_Create( Input_Opt      = Input_Opt,                     &
                             Item           = Item,                          &
                             Id             = N,                             &
                             ContainerId    = 0,                             &
                             Name           = ItemName(N),                   &
                             LongName       = Description,                   &
                             Units          = Units,                         &
                             SpaceDim       = Rank,                          &
                             OnLevelEdges   = OnLevelEdges,                  &
                             DimNames       = ItemDimName(N),                &
                             Operation      = 0,                             &
                             Subset_X       = Subset_X,                      &
                             Subset_Y       = Subset_Y,                      &
                             Subset_Z       = Subset_Z,                      &
                             Source_KindVal = KindVal,                       &
                             Source_0d_8    = Ptr0d_8,                       &
                             Source_1d_8    = Ptr1d_8,                       &
                             Source_2d_4    = Ptr2d_4,                       &
                             RC             = RC                            )

       ! Trap potential error
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Could not create Item: "' // TRIM( ItemName(N) ) // '"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       !### Debug: Print full info about this HISTORY ITEM
       !### You can leave this commented out unless you are debugging
       !CALL HistItem_Print( Item, RC )

       !---------------------------------------------------------------------
       ! Add this item to the Dimension list
       !---------------------------------------------------------------------
       CALL MetaHistItem_AddNew( Input_Opt    = Input_Opt,                   &
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

  END SUBROUTINE IndexVarList_Create
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: IndexVarList_Destroy
!
! !DESCRIPTION: Finalizes the IndexVarList linked list of netCDF
!  index variables.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE IndexVarList_Destroy( IndexVarList, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE MetaHistItem_Mod,  ONLY : MetaHistItem_Destroy
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(MetaHistItem),  POINTER     :: IndexVarList ! Linked list of index
                                                     !  variables for netCDF
!
! !OUTPUT PARAMETERS:
!
    INTEGER,             INTENT(OUT) :: RC           ! Success or failure
!
! !REVISION HISTORY:
!  10 Aug 2017 - R. Yantosca - Initial version
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
    ErrMsg  =  ''
    ThisLoc =  ' -> at MetaHistItem_Destroy (in History/metahistitem_mod.F90)'

    !=======================================================================
    ! Destroy the METAHISTORY ITEM list of index variables for netCDF
    !=======================================================================
    CALL MetaHistItem_Destroy( IndexVarList, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Cannot deallocate the "IndexVarList" META HISTORY ITEM!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

  END SUBROUTINE IndexVarList_Destroy
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get_Number_Of_Levels
!
! !DESCRIPTION: Given the vertical dimension of the container, returns the
!  values NLEV (number of levels) and NILEV (number of level interfaces).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Get_Number_Of_Levels( Container, nLev, nILev )
!
! !USES:
!
    USE HistContainer_Mod, ONLY : HistContainer
!
! !INPUT PARAMETERS:
!
    TYPE(HistContainer), POINTER     :: Container ! Diagnostic collection obj
!
! !OUTPUT PARAMETERS:
!
    INTEGER,             INTENT(OUT) :: nLev      ! Number of levels
    INTEGER,             INTENT(OUT) :: nIlev     ! Number of level interfaces
!
! !REVISION HISTORY:
!  05 Jun 2019 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC

    !=======================================================================
    ! Pick the dimensions of the lev and ilev variables properly
    ! so that we can use that for writing to then netCDF files.
    !
    ! If the vertical dimension (Container%NZ) is undefined, then
    ! this indicates that there is only 2-D data in the collection.
    ! Thus, there will be 1 level (the surface) and 2 level edges.
    !=======================================================================
    IF ( Container%OnLevelEdges ) THEN
       nILev = MAX( Container%NZ, 2 )
       nLev  = nILev - 1
    ELSE
       nLev  = MAX( Container%NZ, 1 )
       nILev = nLev  + 1
    ENDIF

  END SUBROUTINE Get_Number_Of_Levels
!EOC
END MODULE History_Netcdf_Mod
