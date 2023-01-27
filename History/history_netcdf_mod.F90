!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: history_netcdf_mod.F90
!
! !DESCRIPTION: Contains routines to create a netCDF file for each GEOS-Chem
!  diagnostic collection (as specified by each HISTORY CONTAINER in the
!  main collection list located within in history_mod.F90).
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
    TYPE(MetaHistItem), POINTER :: Current

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
      !Container%ReferenceJd   = UNDEFINED_DBL
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
          Current%Item%NcBDimId  = UNDEFINED_INT
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
    USE CharPak_Mod,            ONLY : To_UpperCase
    USE ErrCode_Mod
    USE HistContainer_Mod,      ONLY : HistContainer, HistContainer_Print
    USE HistItem_Mod,           ONLY : HistItem,      HistItem_Print
    USE History_Util_Mod
    USE Input_Opt_Mod,          ONLY : OptInput
    USE JulDay_Mod,             ONLY : CalDate
    USE MetaHistItem_Mod,       ONLY : MetaHistItem
    USE Ncdf_Mod
    USE Registry_Params_Mod,    ONLY : KINDVAL_F4
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
!  For instantaneous file collections, if a file already exists (e.g. at 0h
!  on the day when a run ended), then we will append into that file instead
!  of opening a new file. This will prevent clobbering of existing data.
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
    LOGICAL                     :: appendToFile, isRestart
    INTEGER                     :: N,            VarCt
    INTEGER                     :: VarXDimId,    VarYDimId
    INTEGER                     :: VarZDimId,    VarTDimId
    INTEGER                     :: yyyymmdd,     hhmmss
    INTEGER                     :: nLev,         nILev
    INTEGER                     :: DataType,     RC2
    REAL(f8)                    :: offset

    ! Strings
    CHARACTER(LEN=5)            :: Z
    CHARACTER(LEN=8)            :: D
    CHARACTER(LEN=10)           :: T
    CHARACTER(LEN=255)          :: FileName
    CHARACTER(LEN=255)          :: ErrMsg,       ThisLoc,     VarUnits
    CHARACTER(LEN=255)          :: VarAxis,      VarPositive, VarCalendar
    CHARACTER(LEN=255)          :: VarStdName,   VarFormula,  VarBounds

    ! Arrays
    INTEGER                     :: V(8)

    ! Pointers
    REAL(fp),           POINTER :: Data1d(:)
    REAL(fp),           POINTER :: Data2d(:)

    ! Objects
    TYPE(MetaHistItem), POINTER :: Current
    TYPE(MetaHistItem), POINTER :: IndexVarList

    !========================================================================
    ! Initialize
    !========================================================================
    RC           =  GC_SUCCESS
    RC2          =  GC_SUCCESS
    Current      => NULL()
    IndexVarList => NULL()
    ErrMsg       =  ''
    ThisLoc      =  ' -> at History_Netcdf_Define (in History/history_mod.F90)'
    FileName     =  ''
    VarAxis      =  ''
    VarCalendar  =  ''
    VarPositive  =  ''
    VarUnits     =  ''
    appendToFile = .FALSE.
    yyyymmdd     =  Container%CurrentYmd
    hhmmss       =  Container%CurrentHms

    ! Test if this collection is a restart file (for which we will always
    ! want to create a new file instead of appending to an open file).
    isRestart = ( INDEX( To_UpperCase(TRIM(Container%Name)), 'RESTART' ) > 0 )

    !========================================================================
    ! Create the netCDF file with global attributes (or append)
    ! Do not exit netCDF define mode just yet
    !========================================================================
    IF ( .not. Container%IsFileOpen ) THEN

       !---------------------------------------------------------------------
       ! Compute reference date and time fields in the HISTORY CONTAINER
       ! These are needed to compute the time stamps for each data field
       ! that is written to the netCDF file.  Also resets the current
       ! time slice index.
       !---------------------------------------------------------------------

       ! Reset the current time slice index
       Container%CurrTimeSlice = 0

       IF ( Container%Operation == ACCUM_FROM_SOURCE ) THEN

          !------------------------------------------------------------------
          ! %%% TIME-AVERAGED COLLECTIONS %%%
          !------------------------------------------------------------------

          ! REFERENCE TIMESTAMP:
          ! Subtract the file write alarm interval that we added to
          ! the current date/time (CurrentJd) field at initialization
          Container%ReferenceJsec = Container%CurrentJsec                    &
                                  - Container%FileWriteIvalSec

       ELSE

          !------------------------------------------------------------------
          ! %%% INSTANTANEOUS COLLECTIONS %%%
          !------------------------------------------------------------------

          ! REFERENCE TIMESTAMP: Use the date/time when the file is created.
          Container%ReferenceJsec = Container%CurrentJsec

       ENDIF

       ! Convert reference time from Astronomical Julian Seconds to date/time
       Container%ReferenceJd = Container%ReferenceJsec / SECONDS_PER_DAY
       CALL CalDate( JulianDay = Container%ReferenceJd,                      &
                     yyyymmdd  = Container%ReferenceYmd,                     &
                     hhmmss    = Container%ReferenceHms                     )

       ! Replace time and date tokens in the netCDF file name
       ! (Save the collection's file name in a temporary variable)
       FileName = TRIM( Container%FileName )
       CALL Expand_Date_Time( DateStr    = FileName,                         &
                              yyyymmdd   = Container%ReferenceYmd,           &
                              hhmmss     = Container%ReferenceHms,           &
                              MAPL_Style = .TRUE.                           )

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

       !---------------------------------------------------------------------
       ! Do not create netCDF file on first timestep if instantaneous
       ! collection and frequency = duration. This will avoid creating
       ! of a netCDF file containing all missing values.
       !---------------------------------------------------------------------
       IF ( Container%Operation      == COPY_FROM_SOURCE            .and.    &
            Container%UpdateIvalSec  == Container%FileCloseIvalSec  .and.    &
            Container%FileCloseAlarm == 0.0                         .and.    &
            TRIM(Container%Name) .ne. 'BoundaryConditions' )  THEN
          RETURN

       ELSE

          !------------------------------------------------------------------
          ! Create the file and add global attributes
          ! Remain in netCDF define mode upon exiting this routine
          !
          ! NOTE: Container%Reference is a global attribute lists the
          ! GEOS-Chem web and wiki page.  It has nothing to do with the
          ! reference data/time computed by History_Set_RefDateTime.
          !------------------------------------------------------------------

          ! Echo info about the file we are creating
          IF ( Input_Opt%amIRoot ) THEN
             WRITE( 6, 100 ) TRIM( Container%Name ),                         &
                             Container%ReferenceYmd,                         &
                             Container%ReferenceHms
             WRITE( 6, 110 ) TRIM( FileName       )
          ENDIF

          ! Create the file
          CALL Nc_Create( Create_Nc4     = .TRUE.,                           &
                          NcFile         = fileName,                         &
                          nLon           = Container%nX,                     &
                          nLat           = Container%nY,                     &
                          nLev           = nLev,                             &
                          nIlev          = nILev,                            &
                          nTime          = NF_UNLIMITED,                     &
                          nBounds        = 2,                                &
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
                          boundsId       = Container%bDimId,                 &
                          KeepDefMode    = .TRUE.,                           &
                          Varct          = VarCt                            )

       ENDIF

       !---------------------------------------------------------------------
       ! Denote that the file has been created and is open
       !---------------------------------------------------------------------
       Container%IsFileOpen = .TRUE.

    ENDIF

    ! Format strings for use above
100 FORMAT( '     - Creating file for ',  a, '; reference = ',i8.8,1x,i6.6 )
110 FORMAT( '        with filename = ', a                                  )

    !========================================================================
    ! Define all of the Create the netCDF file with global attributes
    ! Skip if we are appending to an existing file
    !========================================================================
    IF ( .not. Container%IsFileDefined ) THEN

       !---------------------------------------------------------------------
       ! Define the index variables
       !---------------------------------------------------------------------

       ! Define the linked list of index variables (IndexVarList) that
       ! have the same dimension subsets as the current container
       CALL IndexVarList_Create( Input_Opt, Container, IndexVarList, RC )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Error encountered in "IndexVarList_Create"!'
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
                               bDimId       = Container%bDimId,              &
                               RefDate      = Container%ReferenceYmd,        &
                               RefTime      = Container%ReferenceHms,        &
                               OnLevelEdges = Container%OnLevelEdges,        &
                               Item         = Current%Item,                  &
                               VarAxis      = VarAxis,                       &
                               VarPositive  = VarPositive,                   &
                               VarCalendar  = VarCalendar,                   &
                               VarUnits     = VarUnits,                      &
                               VarStdName   = VarStdName,                    &
                               VarFormula   = VarFormula,                    &
                               VarBounds    = VarBounds                     )

          ! Set a flag for the precision of the data
          IF ( Current%Item%Output_KindVal == KINDVAL_F4 ) THEN
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
                           boundsId     = Current%Item%NcBDimId,             &
                           VarLongName  = Current%Item%LongName,             &
                           VarUnit      = VarUnits,                          &
                           Axis         = VarAxis,                           &
                           Calendar     = VarCalendar,                       &
                           Positive     = VarPositive,                       &
                           StandardName = VarStdName,                        &
                           FormulaTerms = VarFormula,                        &
                           Bounds       = VarBounds                         )

          ! Debug print
          !CALL HistItem_Print( Current%Item, RC )

          ! Go to next entry in the list of HISTORY ITEMS
          Current => Current%Next
       ENDDO

       ! Free pointers
       Current => NULL()

       !---------------------------------------------------------------------
       ! Then define each HISTORY ITEM belonging to this collection
       !---------------------------------------------------------------------

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
                               bDimId   = Container%bDimId,                  &
                               Item     = Current%Item,                      &
                               VarUnits = VarUnits                          )

          ! Replace "TBD"  with the current units of State_Chm%Species(:)%Conc
          IF ( TRIM( VarUnits ) == 'TBD' ) THEN
             VarUnits = Container%Spc_Units
          ENDIF

          ! Set a flag for the precision of the data
          IF ( Current%Item%Output_KindVal == KINDVAL_F4 ) THEN
             DataType = 4
          ELSE
             DataType = 8
          ENDIF

          !---------------------------------------------------------------
          ! Define a HISTORY ITEM in this collection as a 4-byte real
          ! or 8-byte real data variable for the netCDF file output
          !---------------------------------------------------------------
          CALL Nc_Var_Def( DefMode      = .TRUE.,                         &
                           Compress     = .TRUE.,                         &
                           fId          = Container%FileId,               &
                           DataType     = DataType,                       &
                           VarName      = Current%Item%Name,              &
                           VarCt        = Current%Item%NcVarId,           &
                           timeId       = Current%Item%NcTDimId,          &
                           levId        = Current%Item%NcZDimId,          &
                           iLevId       = Current%Item%NcIDimId,          &
                           latId        = Current%Item%NcYDimId,          &
                           lonId        = Current%Item%NcXDimId,          &
                           boundsId     = Current%Item%NcBDimId,          &
                           varLongName  = Current%Item%LongName,          &
                           varUnit      = VarUnits,                       &
                          !missingValue = Current%Item%MissingValue,      &
                           avgMethod    = Current%Item%AvgMethod,         &
                           bounds       = VarBounds                      )

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

       !---------------------------------------------------------------------
       ! Write the index variable data
       !---------------------------------------------------------------------

       ! Close definition section
       CALL Nc_Set_DefMode( Container%FileId, Off=.TRUE. )

       ! Set CURRENT to the first node in the list of HISTORY ITEMS
       Current => IndexVarList

       ! As long as this node of the list is valid ...
       DO WHILE( ASSOCIATED( Current ) )

          ! Write data for index variables to the netCDF file
          IF ( Current%Item%SpaceDim == 2 ) THEN

             ! Check the dimension names
             SELECT CASE( Current%Item%DimNames )

                ! lon_bnds or lat_bnds
                CASE( 'bx', 'by' )
                   CALL Nc_Var_Write( fId     = Container%FileId,            &
                                      VarName = Current%Item%Name,           &
                                      Arr2d   = Current%Item%Source_2d_8    )

                ! AREA
                CASE DEFAULT
                   CALL Nc_Var_Write( fId     = Container%FileId,            &
                                      VarName = Current%Item%Name,           &
                                      Arr2d   = Current%Item%Source_2d_4    )
             END SELECT

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

       !---------------------------------------------------------------------
       ! We can now consider this collection to have been "defined"
       !---------------------------------------------------------------------
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
  SUBROUTINE History_Netcdf_Write( Input_Opt, State_Diag, Container, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE HistItem_Mod,        ONLY : HistItem
    USE HistContainer_Mod,   ONLY : HistContainer
    USE History_Util_Mod
    USE Input_Opt_Mod,       ONLY : OptInput
    USE State_Diag_Mod,      ONLY : DgnState
    USE M_Netcdf_Io_Write,   ONLY : NcWr
    USE MetaHistItem_Mod,    ONLY : MetaHistItem
    USE Registry_Params_Mod, ONLY : KINDVAL_F4

!
! !INPUT PARAMETERS:
!
    TYPE(OptInput),      INTENT(IN)  :: Input_Opt  ! Input Options object
    TYPE(DgnState),      INTENT(IN)  :: State_Diag ! Diagnostics state obj
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
    LOGICAL                     :: output4Bytes
    INTEGER                     :: NcFileId,         NcVarId
    INTEGER                     :: Dim1,             Dim2,       Dim3

    ! Strings
    CHARACTER(LEN=255)          :: ErrMsg,           ThisLoc

    ! Arrays
    INTEGER                     :: St1d(1),          Ct1d(1)
    INTEGER                     :: St2d(2),          Ct2d(2)
    INTEGER                     :: St3d(3),          Ct3d(3)
    INTEGER                     :: St4d(4),          Ct4d(4)
    REAL(f4),       ALLOCATABLE :: NcData_1d4(:    )
    REAL(f4),       ALLOCATABLE :: NcData_2d4(:,:  )
    REAL(f4),       ALLOCATABLE :: NcData_3d4(:,:,:)
    REAL(f8),       ALLOCATABLE :: NcData_1d8(:    )
    REAL(f8),       ALLOCATABLE :: NcData_2d8(:,:  )
    REAL(f8),       ALLOCATABLE :: NcData_3d8(:,:,:)
    REAL(f8)                    :: NcTimeVal (1    )

    ! Objects
    TYPE(MetaHistItem), POINTER :: Current
    TYPE(HistItem),     POINTER :: Item

    !========================================================================
    ! Make sure the netCDF file is open and defined
    !========================================================================
    IF ( ( .not. Container%IsFileOpen   )    .and.                           &
         ( .not. Container%IsFileDefined ) ) THEN
       RC     = GC_FAILURE
       ErrMsg = 'NetCDF file is not open or defined for collection: '     // &
                 TRIM( Container%Name )
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !========================================================================
    ! Initialize
    !========================================================================
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

    !========================================================================
    ! Compute time elapsed since the reference time
    !========================================================================

    ! Increment the time index for the netCDF file
    Container%CurrTimeSlice = Container%CurrTimeSlice + 1

    !========================================================================
    ! Compute the time stamp value for the current time slice
    !========================================================================

    ! Compute the elapsed time in seconds since the file creation
    CALL Compute_Elapsed_Time( CurrentJsec  = Container%CurrentJsec,         &
                               TimeBaseJsec = Container%ReferenceJsec,       &
                               ElapsedSec   = Container%TimeStamp           )

    ! For time-averaged collections, we need to subtract the file write
    ! interval from the current time.  This will make sure that time[0] = 0,
    ! or in other words, that the first time slice matches up to the
    ! reference datetime.  -- Bob Yantosca (13 Jan 2023)
    IF ( Container%Operation == ACCUM_FROM_SOURCE ) THEN
       Container%TimeStamp = Container%TimeStamp - Container%FileWriteIvalSec
    ENDIF

    ! Convert to minutes since the reference time
    Container%TimeStamp = Container%TimeStamp / SECONDS_PER_MINUTE

    ! Debug output
    IF ( Input_Opt%LPRT .and. Input_opt%amIRoot ) THEN
       WRITE( 6, 110 ) TRIM( Container%name ), Container%TimeStamp
110    FORMAT( '     - Writing data to ', a, '; timestamp = ', f13.4 )
    ENDIF

    !========================================================================
    ! Write the time stamp to the netCDF File
    !========================================================================

    ! netCDF start and count arrays
    St1d      = (/ Container%CurrTimeSlice /)
    Ct1d      = (/ 1                       /)

    ! Time stamp value
    NcTimeVal = (/ Container%TimeStamp     /)

    ! Write the time stamp to the file
    CALL NcWr( NcTimeVal, NcFileId, 'time', St1d, Ct1d )

    !========================================================================
    ! Loop over all of the HISTORY ITEMS belonging to this collection
    !========================================================================

    ! Set CURRENT to the first entry in the list of HISTORY ITEMS
    Current => Container%HistItems

    ! As long as this entry of the list is valid ...
    DO WHILE( ASSOCIATED( Current ) )

       ! Point to the HISTORY ITEM object in this entry
       Item => Current%Item

       ! Does this HISTORY ITEM request output as 4-byte reals?
       ! If not, we will assume output will be 8-byte reals.
       output4Bytes = ( Item%Output_KindVal == KINDVAL_F4 )

       !---------------------------------------------------------------------
       ! For instantaneous diagnostic quantities:
       ! (1) Copy the Item's data array to the 4-byte or 8-byte local array
       ! (2) Zero the Item's data array
       ! (3) Zero the Item's update counter
       !
       ! For time-averaged diagnostic quantities:
       ! (1) Divide the Item's data array by the number diagnostic updates
       ! (2) Copy the Item's data array to the 4-byte or 8-byte local array
       ! (3) Zero the Item's data array
       ! (4) Zero the Item's update counter
       !---------------------------------------------------------------------
       SELECT CASE( Item%SpaceDim )

          !------------------------------------------------------------------
          ! 3-D data
          !------------------------------------------------------------------
          CASE( 3 )

             ! Get dimensions of data
             Dim1 = SIZE( Item%Data_3d, 1 )
             Dim2 = SIZE( Item%Data_3d, 2 )
             Dim3 = SIZE( Item%Data_3d, 3 )

             ! Get average for satellite diagnostic:
             IF ( Container%name == 'SatDiagn' ) THEN
                Item%Data_3d = Item%Data_3d / State_Diag%SatDiagnCount
                Item%nUpdates = 1.0
             ENDIF

             ! Allocate the 4-byte or 8-byte output array
             IF ( output4bytes ) THEN
                ALLOCATE( NcData_3d4( Dim1, Dim2, Dim3 ), STAT=RC )
             ELSE
                ALLOCATE( NcData_3d8( Dim1, Dim2, Dim3 ), STAT=RC )
             ENDIF

             ! Copy or average the data and store in a 4-byte or 8-byte array
             IF ( Item%Operation == COPY_FROM_SOURCE ) THEN

                !%%% Instantaneous output %%%
                IF ( output4Bytes ) THEN
                   NcData_3d4 = Item%Data_3d
                ELSE
                   NcData_3d8 = Item%Data_3d
                ENDIF
                Item%Data_3d  = 0.0_f8
                Item%nUpdates = 0.0_f8

             ELSE

                !%%% Time-averaged output %%%
                Item%Data_3d  = Item%Data_3d / Item%nUpdates
                IF ( output4Bytes ) THEN
                   NcData_3d4 = Item%Data_3d
                ELSE
                   NcData_3d8 = Item%Data_3d
                ENDIF
                Item%Data_3d  = 0.0_f8
                Item%nUpdates = 0.0_f8

             ENDIF

             ! Compute start and count fields
             St4d = (/ 1,    1,    1,    Container%CurrTimeSlice /)
             Ct4d = (/ Dim1, Dim2, Dim3, 1                       /)

             ! Write data to disk and deallocate output array
             IF ( output4bytes ) THEN
                CALL NcWr( NcData_3d4, NcFileId, Item%Name, St4d, Ct4d )
                DEALLOCATE( NcData_3d4, STAT=RC )
             ELSE
                CALL NcWr( NcData_3d8, NcFileId, Item%Name, St4d, Ct4d )
                DEALLOCATE( NcData_3d8, STAT=RC )
             ENDIF

          !------------------------------------------------------------------
          ! 2-D data
          !------------------------------------------------------------------
          CASE( 2 )

             ! Get dimensions of data
             Dim1 = SIZE( Item%Data_2d, 1 )
             Dim2 = SIZE( Item%Data_2d, 2 )

             ! Get average for satellite diagnostic:
             IF ( Container%name == 'SatDiagn' ) THEN
                Item%Data_2d = Item%Data_2d / State_Diag%SatDiagnCount(:,:,1)
                Item%nUpdates = 1.0
             ENDIF

             ! Allocate the 4-byte or 8-byte output array
             IF ( output4bytes ) THEN
                ALLOCATE( NcData_2d4( Dim1, Dim2 ), STAT=RC )
             ELSE
                ALLOCATE( NcData_2d8( Dim1, Dim2 ), STAT=RC )
             ENDIF

             ! Copy or average the data and store in a 4-byte or 8-byte array
             IF ( Item%Operation == COPY_FROM_SOURCE ) THEN

                !%%% Instantaneous output %%%
                IF ( output4bytes ) THEN
                   NcData_2d4 = Item%Data_2d
                ELSE
                   NcData_2d8 = Item%Data_2d
                ENDIF
                Item%Data_2d  = 0.0_f8
                Item%nUpdates = 0.0_f8

             ELSE

                !%%% Time-averaged output %%%
                Item%Data_2d  = Item%Data_2d / Item%nUpdates
                IF ( output4bytes ) THEN
                   NcData_2d4 = Item%Data_2d
                ELSE
                   NcData_2d8 = Item%Data_2d
                ENDIF
                Item%Data_2d  = 0.0_f8
                Item%nUpdates = 0.0_f8

             ENDIF

             ! Compute start and count fields
             St3d = (/ 1,    1,    Container%CurrTimeSlice /)
             Ct3d = (/ Dim1, Dim2, 1                       /)

             ! Write data to disk
             IF ( output4bytes ) THEN
                CALL NcWr( NcData_2d4, NcFileId, Item%Name, St3d, Ct3d )
                DEALLOCATE( NcData_2d4, STAT=RC )
             ELSE
                CALL NcWr( NcData_2d8, NcFileId, Item%Name, St3d, Ct3d )
                DEALLOCATE( NcData_2d8, STAT=RC )
             ENDIF

          !------------------------------------------------------------------
          ! 1-D data
          !------------------------------------------------------------------
          CASE( 1 )

             ! Get dimensions of data
             Dim1 = SIZE( Item%Data_1d, 1 )

             ! Allocate the 4-byte or 8-byte output array
             IF ( output4bytes ) THEN
                ALLOCATE( NcData_1d4( Dim1 ), STAT=RC )
             ELSE
                ALLOCATE( NcData_1d8( Dim1 ), STAT=RC )
             ENDIF

             ! Copy or average the data and store in a 4-byte or 8-byte array
             IF ( Item%Operation == COPY_FROM_SOURCE ) THEN

                !%%% Instantaneous output %%%
                IF ( output4bytes ) THEN
                   NcData_1d4 = Item%Data_1d
                ELSE
                   NcData_1d8 = Item%Data_1d
                ENDIF
                Item%Data_1d  = 0.0_f8
                Item%nUpdates = 0.0_f8

             ELSE

                ! %%% Time-averaged output %%%
                Item%Data_1d  = Item%Data_1d / Item%nUpdates
                IF ( output4bytes ) THEN
                   NcData_1d4 = Item%Data_1d
                ELSE
                   NcData_1d8 = Item%Data_1d
                ENDIF
                Item%Data_1d  = 0.0_f8
                Item%nUpdates = 0.0_f8

             ENDIF

             ! Compute start and count fields
             St2d = (/ 1,    Container%CurrTimeSlice /)
             Ct2d = (/ Dim1, 1                       /)

             ! Write data to disk
             IF ( output4bytes ) THEN
                CALL NcWr( NcData_1d4, NcFileId, Item%Name, St2d, Ct2d )
                DEALLOCATE( NcData_1d4, STAT=RC )
             ELSE
                CALL NcWr( NcData_1d8, NcFileId, Item%Name, St2d, Ct2d )
                DEALLOCATE( NcData_1d8, STAT=RC )
             ENDIF

       END SELECT

       !---------------------------------------------------------------------
       ! Go to next entry in the list of HISTORY ITEMS
       !---------------------------------------------------------------------
       Current => Current%Next
       Item    => NULL()
    ENDDO

    !---------------------------------------------------------------------
    ! Set count for satellite diagnostic to zero:
    !---------------------------------------------------------------------
    IF ( State_Diag%Archive_SatDiagnCount ) THEN
       State_Diag%SatDiagnCount = 0.0e+0_fp
    ENDIF

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

    !========================================================================
    ! Initialize
    !========================================================================
    IF ( PRESENT( MAPL_Style ) ) THEN
       Is_Mapl_Style = MAPL_Style
    ELSE
       Is_Mapl_Style = .FALSE.
    ENDIF

    !========================================================================
    ! Split the date and time into individal variables
    !========================================================================

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

    !========================================================================
    ! Replace the date and time tokens in the string
    !========================================================================

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
  SUBROUTINE Get_Var_DimIds( Item,         xDimId,   yDimId,                 &
                             zDimId,       iDimId,   tDimID,                 &
                             bDimId,       RefDate,  RefTime,                &
                             OnLevelEdges, VarAxis,  VarCalendar,            &
                             VarPositive,  VarUnits, VarStdName,             &
                             VarFormula,   VarBounds                        )
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
    INTEGER,            INTENT(IN)  :: bDimId       ! Id # of B (bounds  ) dim
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
    CHARACTER(LEN=255), OPTIONAL    :: VarBounds    ! X or Y bounds var name
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
    CHARACTER(LEN=255) :: TmpBounds

    !========================================================================
    ! Initialize
    !========================================================================
    TmpAxis         = ''
    TmpCalendar     = ''
    TmpPositive     = ''
    TmpUnits        = Item%Units
    TmpStdName      = ''
    TmpFormula      = ''
    TmpBounds       = ''
    Item%NcXDimId   = UNDEFINED_INT
    Item%NcYDimId   = UNDEFINED_INT
    Item%NcZDimId   = UNDEFINED_INT
    Item%NcIDimId   = UNDEFINED_INT
    Item%NcTDimId   = UNDEFINED_INT
    Item%NcBDimId   = UNDEFINED_INT

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

    !========================================================================
    ! Return relevant dim ID's and metadata for the HISTORY ITEM
    !========================================================================
    SELECT CASE( TRIM( Item%Name ) )

       ! lon
       CASE( 'lon' )
          Item%NcXDimId = xDimId
          TmpAxis       = 'X'
          TmpBounds     = 'lon_bnds'

       ! lat
       CASE( 'lat' )
          Item%NcYDimId = yDimId
          TmpAxis       = 'Y'
          TmpBounds     = 'lat_bnds'

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

             ! For longitude bounds
             CASE( 'bx' )
                Item%NcBDimId = bDimId
                Item%NcXDimId = xDimId

             ! For latitude bounds
             CASE( 'by' )
                Item%NcBDimId = bDimId
                Item%NcYDimId = yDimId

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
    IF ( PRESENT( VarBounds   ) ) VarBounds   = TmpBounds

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
    INTEGER                  :: Output_KindVal
    INTEGER                  :: Source_KindVal
    INTEGER                  :: Rank
    INTEGER                  :: nILev
    INTEGER                  :: nLev

    ! Arrays
    INTEGER                  :: Dimensions(3)
    INTEGER                  :: Subset_X(2), Subset_Xc(2)
    INTEGER                  :: Subset_Y(2), Subset_Yc(2)
    INTEGER                  :: Subset_Z(2), Subset_Zc(2), Subset_Ze(2)

    ! Strings
    CHARACTER(LEN=20)        :: ItemDimName(13)
    CHARACTER(LEN=20)        :: ItemName(13)
    CHARACTER(LEN=20)        :: RegistryName(13)
    CHARACTER(LEN=255)       :: Description
    CHARACTER(LEN=255)       :: ErrMsg
    CHARACTER(LEN=255)       :: ThisLoc
    CHARACTER(LEN=255)       :: Units

    ! Pointer arrays
    REAL(f8),        POINTER :: Ptr0d_8
    REAL(f8),        POINTER :: Ptr1d_8(:  )
    REAL(f8),        POINTER :: Ptr2d_8(:,:)
    REAL(f4),        POINTER :: Ptr2d_4(:,:)

    ! Objects
    TYPE(HistItem),  POINTER :: Item

    !========================================================================
    ! Initialize
    !========================================================================
    RC             =  GC_SUCCESS
    Description    =  ''
    Dimensions     =  0
    Source_KindVal =  0
    Output_KindVal =  0
    Rank           =  0
    Units          =  ''
    ErrMsg         =  ''
    ThisLoc        =  &
         ' -> at History_Netcdf_Init (in History/history_mod.F90)'
    Ptr0d_8      => NULL()
    Ptr1d_8      => NULL()
    Ptr2d_4      => NULL()

    !========================================================================
    ! Define the names that will be used to create the HISTORY ITEMS
    ! for fields to be used as netCDF metadata
    !========================================================================

    ! Fields saved in the Registry object in GeosUtil/grid_registry_mod.F90
    RegistryName(1 ) = 'GRID_AREA'
    RegistryName(2 ) = 'GRID_P0'
    RegistryName(3 ) = 'GRID_HYBI'
    RegistryName(4 ) = 'GRID_HYAI'
    RegistryName(5 ) = 'GRID_HYBM'
    RegistryName(6 ) = 'GRID_HYAM'
    RegistryName(7 ) = 'GRID_LON'
    RegistryName(8 ) = 'GRID_LONBND'
    RegistryName(9 ) = 'GRID_LAT'
    RegistryName(10) = 'GRID_LATBND'
    RegistryName(11) = 'GRID_ILEV'
    RegistryName(12) = 'GRID_LEV'
    RegistryName(13) = 'GRID_TIME'

    ! Name for each HISTORY ITEM
    ItemName(1 )     = 'AREA'
    ItemName(2 )     = 'P0'
    ItemName(3 )     = 'hybi'
    ItemName(4 )     = 'hyai'
    ItemName(5 )     = 'hybm'
    ItemName(6 )     = 'hyam'
    ItemName(7 )     = 'lon'
    ItemName(8 )     = 'lon_bnds'
    ItemName(9 )     = 'lat'
    ItemName(10)     = 'lat_bnds'
    ItemName(11)     = 'ilev'
    ItemName(12)     = 'lev'
    ItemName(13)     = 'time'

    ! Dimensions for each HISTORY ITEM
    ItemDimName(1 )  = 'xy'
    ItemDimName(2 )  = '-'
    ItemDimName(3 )  = 'z'
    ItemDimName(4 )  = 'z'
    ItemDimName(5 )  = 'z'
    ItemDimName(6 )  = 'z'
    ItemDimName(7 )  = 'x'
    ItemDimName(8 )  = 'bx'
    ItemDimName(9 )  = 'y'
    ItemDimName(10)  = 'by'
    ItemDimName(11)  = 'z'
    ItemDimName(12)  = 'z'
    ItemDimName(13)  = 't'

    !========================================================================
    ! Pick the dimensions of the lev and ilev variables properly
    !========================================================================

    ! Get the number of levels (nLev) and level interfaces (nIlev)
    CALL Get_Number_Of_Levels( Container, nLev, nIlev )

    ! Subset indices
    Subset_Xc = (/ Container%X0, Container%X1 /)
    Subset_Yc = (/ Container%Y0, Container%Y1 /)
    Subset_Zc = (/ Container%Z0, nLev         /)
    Subset_Ze = (/ Container%Z0, nILev        /)

    !========================================================================
    ! Create a HISTORY ITEM for each of the index fields (lon, lat, area)
    ! of grid_registry_mod.F90 and add them to a METAHISTORY ITEM list
    !========================================================================
    DO N = 1, SIZE( RegistryName )

       !---------------------------------------------------------------------
       ! Look up one of the index fields from grid_registry_mod.F90
       !---------------------------------------------------------------------
       CALL Lookup_Grid( Input_Opt      = Input_Opt,                         &
                         Variable       = RegistryName(N),                   &
                         Description    = Description,                       &
                         Dimensions     = Dimensions,                        &
                         Source_KindVal = Source_KindVal,                    &
                         Output_KindVal = Output_KindVal,                    &
                         Rank           = Rank,                              &
                         Units          = Units,                             &
                         OnLevelEdges   = OnLevelEdges,                      &
                         Ptr0d_8        = Ptr0d_8,                           &
                         Ptr1d_8        = Ptr1d_8,                           &
                         Ptr2d_8        = Ptr2d_8,                           &
                         Ptr2d_4        = Ptr2d_4,                           &
                         RC             = RC                                )

       ! Trap potential error
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Error in "Lookup_Grid" for diagnostic ' //               &
                   TRIM( RegistryName(N) )
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       ! Pick proper subset indices for index variables placed on
       ! edges  (hyai, hybi, ilev) or centers (everything else)
       ! NOTE: special handling for lat_bnds (12) and lon_bnds (13)
       SELECT CASE( N )
          CASE( 2 )                   ! P0
             Subset_X = (/ 0, 0 /)
             Subset_Y = (/ 0, 0 /)
             Subset_Z = (/ 0, 0 /)
          CASE( 3, 4, 11 )            ! hybi, hyai, ilev
             Subset_X = Subset_Xc
             Subset_Y = Subset_Yc
             Subset_Z = Subset_Ze
          CASE( 8  )                  ! lon_bnds
             Subset_X = (/ 1, 2 /)
             Subset_Y = Subset_Xc
             Subset_Z = (/ 0, 0 /)
          CASE( 10 )                  ! lat_bnds
             Subset_X = (/ 1, 2 /)
             Subset_Y = Subset_Yc
             Subset_Z = (/ 0, 0 /)
          CASE DEFAULT                ! everything else
             Subset_X = Subset_Xc
             Subset_Y = Subset_Yc
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
                             Source_KindVal = Source_KindVal,                &
                             Output_KindVal = Output_KindVal,                &
                             Source_0d_8    = Ptr0d_8,                       &
                             Source_1d_8    = Ptr1d_8,                       &
                             Source_2d_4    = Ptr2d_4,                       &
                             Source_2d_8    = Ptr2d_8,                       &
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
       CALL MetaHistItem_AddNew( Input_Opt = Input_Opt,                      &
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

       ! Free Item now that we have added it to IndexVarList
       DEALLOCATE( Item )
       Item => NULL()
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

    !========================================================================
    ! Initialize
    !========================================================================
    RC      = GC_SUCCESS
    ErrMsg  =  ''
    ThisLoc =  ' -> at MetaHistItem_Destroy (in History/metahistitem_mod.F90)'

    !========================================================================
    ! Destroy the METAHISTORY ITEM list of index variables for netCDF
    !========================================================================
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

    !========================================================================
    ! Pick the dimensions of the lev and ilev variables properly
    ! so that we can use that for writing to then netCDF files.
    !
    ! If the vertical dimension (Container%NZ) is undefined, then
    ! this indicates that there is only 2-D data in the collection.
    ! Thus, there will be 1 level (the surface) and 2 level edges.
    !========================================================================
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
