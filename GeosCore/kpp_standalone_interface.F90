!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: kpp_standalone_interface.F90
!
! !DESCRIPTION: Contains routines to print the full chemical state
!  which can be used as input to the KPP Standalone.
!\\
!\\
! !INTERFACE:
!
MODULE KPP_Standalone_Interface
!
! !USES:
!
  USE Precision_Mod
  USE HCO_Error_Mod, ONLY : hp

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBERS:
!
  PUBLIC  :: Check_Domain
  PUBLIC  :: Check_ActiveCell
  PUBLIC  :: Config_KPP_Standalone
  PUBLIC  :: Write_Samples
  PUBLIC  :: Cleanup_KPP_Standalone
!
! !DERIVED TYPES:
!
  ! Type to hold information read from the YAML config file
  TYPE, PRIVATE :: KPP_Standalone_Interface_Type
     INTEGER                         :: NLOC
     LOGICAL                         :: SkipIt
     CHARACTER(LEN=255)              :: Output_Directory
     CHARACTER(LEN=255), ALLOCATABLE :: LocationName(:)
     REAL(hp),           ALLOCATABLE :: LocationLons(:)
     REAL(hp),           ALLOCATABLE :: LocationLats(:)
     INTEGER,            ALLOCATABLE :: IDX(:)
     INTEGER,            ALLOCATABLE :: JDX(:)
     INTEGER,            ALLOCATABLE :: Levels(:)
  END TYPE KPP_Standalone_Interface_Type

  ! Type to denote active cells
  TYPE, PRIVATE :: KPP_Standalone_ActiveCell_Type
     LOGICAL                         :: Active_Cell
     CHARACTER(LEN=255)              :: Active_Cell_Name
  END TYPE KPP_Standalone_ActiveCell_Type
!
! !PRIVATE DATA MEMBERS:
!
  TYPE(KPP_Standalone_Interface_Type),  PRIVATE :: KPP_Standalone_YAML
  TYPE(KPP_Standalone_ActiveCell_Type), PRIVATE :: KPP_Standalone_ActiveCell
  !$OMP THREADPRIVATE( KPP_Standalone_ActiveCell )
!
! !AUTHORS:
!  P. Obin Sturm (psturm@usc.edu)
!
! !REVISION HISTORY:
!  11 Mar 2024 - P. Obin Sturm - Initial version
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
! !IROUTINE: check_domain
!
! !DESCRIPTION: Subroutine Check_Domain is used to identify if a
! specified latitude and longitude falls within a grid cell on the
! current CPU. Multiple lat/lon pairs can be checked simultaneously.
!\\
!\\
! !INTERFACE:
!
   SUBROUTINE Check_Domain( RC )
!
! !USES:
!
     USE HCO_GeoTools_Mod, ONLY : HCO_GetHorzIJIndex
     USE HCO_State_GC_Mod, ONLY : HcoState
!
! !OUTPUT PARAMETERS
!
     INTEGER, INTENT(out) :: RC
!
! !REVISION HISTORY:
!  11 Mar 2024 - P. Obin Sturm - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
    ! Early exit if no locations
    IF ( KPP_Standalone_YAML%SkipIt ) RETURN

    CALL HCO_GetHorzIJIndex( HcoState,                      &
                             KPP_Standalone_YAML%NLOC,         &
                             KPP_Standalone_YAML%LocationLons, &
                             KPP_Standalone_YAML%LocationLats, &
                             KPP_Standalone_YAML%IDX,          &
                             KPP_Standalone_YAML%JDX,          &
                             RC)

   END SUBROUTINE Check_Domain
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: check_activecell
!
! !DESCRIPTION: Subroutine Check_ActiveCell is used to identify if a grid cell
! is within a specified latitude and longitude to print the full chemical state
! (all concentrations, reaction rates, rate constants, and meteo metadata).
!\\
!\\
! !INTERFACE:
!
   SUBROUTINE Check_ActiveCell( I, J, L )
!
! !INPUT PARAMETERS:
!
     INTEGER, INTENT(IN) :: I, J, L   ! Grid Indices
!
! !REVISION HISTORY:
!  11 Mar 2024 - P. Obin Sturm - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES
!
    INTEGER :: K

    ! Early exit if there was no YAML file or no active cells
    IF ( KPP_Standalone_YAML%SkipIt ) RETURN

    ! Initialize
    KPP_Standalone_ActiveCell%Active_Cell      = .FALSE.
    KPP_Standalone_ActiveCell%Active_Cell_Name = ''

    IF ( ANY( L == KPP_Standalone_YAML%Levels ) ) THEN
       DO K = 1, KPP_Standalone_YAML%NLOC
          IF ( KPP_Standalone_YAML%IDX(K) == I  .AND.                        &
               KPP_Standalone_YAML%JDX(K) == J ) THEN
             KPP_Standalone_ActiveCell%Active_Cell = .TRUE.
             KPP_Standalone_ActiveCell%Active_Cell_Name =                    &
                  KPP_Standalone_YAML%LocationName(K)
          ENDIF
       ENDDO
    ENDIF

  END SUBROUTINE Check_ActiveCell
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Config_KPP_Standalone
!
! !DESCRIPTION: Subroutine Config_KPP_Standalone reads a set of gridcells
!  to be sampled and the full chemical state printed.
!\\
!\\
! !INTERFACE:
!
   SUBROUTINE Config_KPP_Standalone( Input_Opt, RC )
      USE QfYaml_Mod
      USE ErrCode_Mod
      USE Input_Opt_Mod, ONLY : OptInput
      USE RoundOff_Mod,  ONLY : Cast_and_RoundOff
      USE inquireMod,    ONLY : findFreeLUN
! !INPUT PARAMETERS:
!
      TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT)   :: RC          ! Success or failure
!
! !REVISION HISTORY:
!  11 Mar 2024 - P. Obin Sturm - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Scalars
      INTEGER                      :: I, N
      INTEGER                      :: IU_FILE     ! Available unit for writing
      INTEGER                      :: path_exists
      LOGICAL                      :: file_exists
      LOGICAL                      :: v_bool

      ! Strings
      CHARACTER(LEN=255)           :: thisLoc
      CHARACTER(LEN=512)           :: errMsg
      CHARACTER(LEN=QFYAML_NamLen) :: key
      CHARACTER(LEN=QFYAML_StrLen) :: v_str

      ! Objects
      TYPE(QFYAML_t)               :: Config, ConfigAnchored

      ! Arrays
      INTEGER                      :: a_int(QFYAML_MaxArr)

      ! String arrays
      CHARACTER(LEN=QFYAML_NamLen) :: a_str(QFYAML_MaxArr)

      ! YAML configuration file name to be read
      CHARACTER(LEN=30), PARAMETER :: configFile = &
           './kpp_standalone_interface.yml'

      ! Inquire if YAML interface exists -- if not, skip initializing
      KPP_Standalone_YAML%SkipIt = .FALSE.
      INQUIRE( FILE=configFile, EXIST=file_exists )
      IF ( .NOT. file_exists ) THEN
         KPP_Standalone_YAML%SkipIt = .TRUE.
         IF ( Input_Opt%amIRoot ) THEN
            WRITE( 6, 100 ) TRIM( configFile )
 100        FORMAT( "Config file ", a ", not found, ",                       &
                    "skipping KPP standalone interface" )
            RETURN
         ENDIF
      ENDIF

      ! Assume success
      RC      = GC_SUCCESS
      errMsg  = ''
      thisLoc = ' -> at Config_KPP_Standalone (in module GeosCore/kpp_standalone_interface.F90)'

      !========================================================================
      ! Read the YAML file into the Config object
      !========================================================================
      CALL QFYAML_Init( configFile, Config, ConfigAnchored, RC )
      IF ( RC /= GC_SUCCESS ) THEN
         errMsg = 'Error reading configuration file: ' // TRIM( configFile )
         CALL GC_Error( errMsg, RC, thisLoc )
         RETURN
      ENDIF

      !========================================================================
      ! Read the main on/off switch; Exit if the switch is turned off
      !========================================================================
      key = "settings%activate"
      v_bool = MISSING_BOOL
      CALL QFYAML_Add_Get( Config, TRIM( key ), v_bool, "", RC )
      IF ( RC /= GC_SUCCESS ) THEN
         errMsg = 'Error parsing ' // TRIM( key ) // '!'
         CALL GC_Error( errMsg, RC, thisLoc )
         RETURN
      ENDIF
      KPP_Standalone_YAML%SkipIt = ( .not. v_bool )
      IF ( KPP_Standalone_YAML%SkipIt ) THEN
         WRITE( 6, 110 )
 110     FORMAT( "KPP standalone interface was manually disabled" )
         RETURN
      ENDIF

      !========================================================================
      ! Read the list of active cells
      !========================================================================
      key = "active_cells"
      a_str = MISSING_STR
      CALL QFYAML_Add_Get( Config, key, a_str, "", RC, dynamic_size=.TRUE. )
      IF ( RC /= GC_SUCCESS ) THEN
         errMsg = 'Error parsing ' // TRIM( key ) // '!'
         CALL GC_Error( errMsg, RC, thisLoc )
         RETURN
      ENDIF

      !========================================================================
      ! Get the number of active cells (if 0, return) and the list of names
      !========================================================================
      KPP_Standalone_YAML%NLOC = Find_Number_of_Locations( a_str )
      IF ( KPP_Standalone_YAML%NLOC .eq. 0 ) THEN
         ! Set SkipIt flag to short circuit other subroutines
         KPP_Standalone_YAML%SkipIt = .TRUE.
         IF ( Input_Opt%amIRoot ) THEN
            WRITE( 6, 120 )
 120        FORMAT( "No active cells for box modeling ",                      &
                    "in kpp_standalone_interface.yml")
            RETURN
         ENDIF
      ENDIF
      ALLOCATE( KPP_Standalone_YAML%LocationName( KPP_Standalone_YAML%NLOC ), STAT=RC )
      CALL GC_CheckVar( 'KPP_Standalone_YAML%LocationName', 0, RC )
      IF ( RC /= GC_SUCCESS ) RETURN
      DO I = 1,KPP_Standalone_YAML%NLOC
         KPP_Standalone_YAML%LocationName(I) = TRIM( a_str(I) )
      END DO

      !========================================================================
      ! Read latitude and longitude of active cells
      !========================================================================

      ! Allocate number of locations for lats and lons
      ALLOCATE( KPP_Standalone_YAML%LocationLons( KPP_Standalone_YAML%NLOC ), STAT=RC )
      CALL GC_CheckVar( 'KPP_Standalone_YAML%LocationLons', 0, RC )
      IF ( RC /= GC_SUCCESS ) RETURN

      ALLOCATE( KPP_Standalone_YAML%LocationLats( KPP_Standalone_YAML%NLOC ), STAT=RC )
      CALL GC_CheckVar( 'KPP_Standalone_YAML%LocationLats', 0, RC )
      IF ( RC /= GC_SUCCESS ) RETURN

      ! Read coordinates
      DO I = 1,KPP_Standalone_YAML%NLOC
      ! Read longitudes
         key = "locations%"//TRIM( KPP_Standalone_YAML%LocationName(I) )//"%longitude"
         v_str = MISSING_STR
         CALL QFYAML_Add_Get( Config, TRIM( key ), v_str, "", RC )
         IF ( RC /= GC_SUCCESS ) THEN
            errMsg = 'Error parsing ' // TRIM( key ) // '!'
            CALL GC_Error( errMsg, RC, thisLoc )
            RETURN
         ENDIF
         KPP_Standalone_YAML%LocationLons( I ) = Cast_and_RoundOff( TRIM( v_str ), places=-1 )
      ! Read latitudes
         key = "locations%"//TRIM( KPP_Standalone_YAML%LocationName(I) )//"%latitude"
         v_str = MISSING_STR
         CALL QFYAML_Add_Get( Config, TRIM( key ), v_str, "", RC )
         IF ( RC /= GC_SUCCESS ) THEN
            errMsg = 'Error parsing ' // TRIM( key ) // '!'
            CALL GC_Error( errMsg, RC, thisLoc )
            RETURN
         ENDIF
         KPP_Standalone_YAML%LocationLats( I ) = Cast_and_RoundOff( TRIM( v_str ), places=-1 )
      END DO

      ! Allocate IDX and JDX (masks for whether a location is on the CPU)
      ALLOCATE( KPP_Standalone_YAML%IDX( KPP_Standalone_YAML%NLOC ), STAT=RC )
      CALL GC_CheckVar( 'KPP_Standalone_YAML%IDX', 0, RC )
      IF ( RC /= GC_SUCCESS ) RETURN

      ALLOCATE( KPP_Standalone_YAML%JDX( KPP_Standalone_YAML%NLOC ), STAT=RC )
      CALL GC_CheckVar( 'KPP_Standalone_YAML%JDX', 0, RC )
      IF ( RC /= GC_SUCCESS ) RETURN

      KPP_Standalone_YAML%IDX(:) = -1
      KPP_Standalone_YAML%JDX(:) = -1

      !========================================================================
      ! Get the list of levels and number of levels
      !========================================================================
      ! TODO: could add capability for location specific levels
      key = "settings%levels"
      a_int = MISSING_INT
      CALL QFYAML_Add_Get( Config, key, a_int, "", RC, dynamic_size=.TRUE. )
      IF ( RC /= GC_SUCCESS ) THEN
         errMsg = 'Error parsing ' // TRIM( key ) // '!'
         CALL GC_Error( errMsg, RC, thisLoc )
         RETURN
      ENDIF
      N = Find_Number_of_Levels( a_int )
      ! if no specified levels, print the surface
      IF ( N .eq. 0 ) THEN
         N = 1
         a_int(1) = 1
      END IF
      ALLOCATE( KPP_Standalone_YAML%Levels( N ), STAT=RC )
      CALL GC_CheckVar( 'KPP_Standalone_YAML%Levels', 0, RC )
      IF ( RC /= GC_SUCCESS ) RETURN
      DO I = 1,N
         KPP_Standalone_YAML%Levels(I) = a_int(I)
      END DO

      !========================================================================
      ! Set the output directory
      !========================================================================
      ! Get that value
      key = "settings%output_directory"
      v_str = MISSING_STR
      CALL QFYAML_Add_Get( Config, TRIM( key ), v_str, "", RC )
      IF ( RC /= GC_SUCCESS ) THEN
         errMsg = 'Error parsing ' // TRIM( key ) // '!'
         CALL GC_Error( errMsg, RC, thisLoc )
         RETURN
      ENDIF
      ! Check to see if the directory exists
      ! Do this in a portable way that works across compilers
      ! The directory specifier in inquire might be specific to ifort
      ! So instead try to open a test file within the output directory
      ! Check ./OutputDir (which exists for GEOS-Chem and GCHP) as backup
      IU_FILE = findFreeLUN()
      OPEN( IU_FILE, FILE   = trim(v_str)//'/.test_directory_existence',     &
                     ACTION = "WRITE",                                       &
                     IOSTAT = path_exists,                                   &
                     ACCESS = 'SEQUENTIAL'                                  )

      ! If the specified folder doesn't exist, try OutputDir
      IF ( path_exists /= 0 ) THEN
         OPEN( IU_FILE, FILE   = './OutputDir/.test_directory_existence',    &
                        ACTION = "WRITE",                                    &
                        IOSTAT = path_exists,                                &
                        ACCESS ='SEQUENTIAL'                                )
         KPP_Standalone_YAML%Output_Directory = "./OutputDir"
         IF ( Input_Opt%amIRoot ) THEN
            WRITE( 6, '(a)' )                                                &
             "KPP Standalone Interface warning: Specified output directory ",&
             trim(v_str),                                                    &
             " was not found, trying default output path './OutputDir' "
         ENDIF

         ! If OutputDir doesn't exist, write to the current directory
         IF ( path_exists /= 0 ) THEN
            IF ( Input_Opt%amIRoot ) THEN
               WRITE( 6, '(a)' )                                             &
            "KPP Standalone Interface warning: Specified output directory ", &
               trim(v_str),                                                  &
               " and default output directory './OutputDir' "            //  &
               "were not found, writing output to the current directory './'"
               KPP_Standalone_YAML%Output_Directory = "./"
            ENDIF
         ENDIF
      ELSE
         KPP_Standalone_YAML%Output_Directory = trim(v_str)
         close(IU_FILE)
      END IF

      !=======================================================================
      ! Print information about sites that will be archived
      !=======================================================================
      IF ( Input_Opt%amIRoot ) THEN
         WRITE( 6, '(a)'   ) REPEAT( "=", 79 )
         WRITE( 6, '(a,/)' ) "KPP STANDALONE INTERFACE"
         WRITE( 6, '(a,/)' ) "Model state will be archived at these sites:"
         DO I = 1, KPP_Standalone_YAML%NLOC
            WRITE( 6, 150 ) KPP_Standalone_YAML%LocationName(I),             &
                            KPP_Standalone_YAML%LocationLons(I),             &
                            KPP_Standalone_YAML%LocationLats(I)
 150        FORMAT( a25, "( ", f9.4, ", ", f9.4, " )")
         ENDDO
         WRITE( 6, '(/,a)'   ) "For GEOS-Chem vertical levels:"
         WRITE( 6, '(100i4)' ) KPP_Standalone_YAML%Levels
         WRITE( 6, '(a)'   ) REPEAT( "=", 79 )
      ENDIF

   END SUBROUTINE Config_KPP_Standalone
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Write_Samples
!
! !DESCRIPTION: Subroutine Write_Samples writes the full chemical state
! (concentrations, reaction rates and rate constants, meteorological
! conditions).
!\\
!\\
! !INTERFACE:
!
   SUBROUTINE Write_Samples( I,           J,           L,                    &
                             initC,       localRCONST, initHvalue,           &
                             exitHvalue,  State_Grid,  State_Chm,            &
                             State_Met,   Input_Opt,   KPP_TotSteps,         &
                             RC,          FORCE_WRITE, CELL_NAME            )
!
! !USES:
!
      USE ErrCode_Mod
      USE State_Grid_Mod,           ONLY : GrdState
      USE State_Chm_Mod,            ONLY : ChmState
      USE State_Met_Mod,            ONLY : MetState
      USE Input_Opt_Mod,            ONLY : OptInput
      USE GcKpp_Function
      USE GcKpp_Parameters,         ONLY : NSPEC, NREACT, NVAR
      USE TIME_MOD,                 ONLY : GET_TS_CHEM
      USE TIME_MOD,                 ONLY : TIMESTAMP_STRING
      USE TIME_MOD,                 ONLY : Get_Minute
      USE TIME_MOD,                 ONLY : Get_Hour
      USE TIME_MOD,                 ONLY : Get_Day
      USE TIME_MOD,                 ONLY : Get_Month
      USE TIME_MOD,                 ONLY : Get_Year
      USE Pressure_Mod,             ONLY : Get_Pcenter
      USE inquireMod,               ONLY : findFreeLUN
!
! !INPUT PARAMETERS:
!
      INTEGER,        INTENT(IN)    :: I                   ! Longitude index
      INTEGER,        INTENT(IN)    :: J                   ! Latitude index
      INTEGER,        INTENT(IN)    :: L                   ! Vertical level
      INTEGER,        INTENT(IN)    :: KPP_TotSteps        ! Total integr. steps
      TYPE(GrdState), INTENT(IN)    :: State_Grid          ! Grid State object
      TYPE(ChmState), INTENT(IN)    :: State_Chm           ! Chem State obj
      TYPE(MetState), INTENT(IN)    :: State_Met           ! Met State obj
      TYPE(OptInput), INTENT(IN)    :: Input_Opt           ! Input Options obj
      REAL(dp), INTENT(IN)          :: initC(NSPEC)        ! Initial conc.
      REAL(dp), INTENT(IN)          :: localRCONST(NREACT) ! Rate constants
      REAL(dp)                      :: initHvalue          ! Initial timestep
      REAL(dp)                      :: exitHvalue          ! Final timestep:
                                                           !  RSTATE(Nhexit)
      LOGICAL, OPTIONAL             :: FORCE_WRITE         ! Write even if not
                                                           ! in an active cell
      CHARACTER(LEN=*), OPTIONAL    :: CELL_NAME           ! Customize name of
                                                           !  this file
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT)   :: RC          ! Success or failure
!
! !REVISION HISTORY:
!  11 Mar 2024 - P. Obin Sturm - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Integers
      INTEGER            :: N
      INTEGER            :: IU_FILE
      INTEGER            :: SpcID
      REAL(fp)           :: DT
      LOGICAL            :: FORCE_WRITE_AUX
      CHARACTER(LEN=255) :: CELL_NAME_AUX

      ! Strings
      CHARACTER(LEN=255) :: YYYYMMDD_hhmmz
      CHARACTER(LEN=255) :: level_string
      CHARACTER(LEN=512) :: errMsg, filename

      ! Arrays
      REAL(dp)           :: Aout(NREACT)
      REAL(dp)           :: Vloc(NVAR)

      !======================================================================
      ! Write_Samples begins here!
      !======================================================================

      ! Did a user want to write the chemical state
      ! even if not in an active cell?
      FORCE_WRITE_AUX = .FALSE.
      IF ( PRESENT( FORCE_WRITE ) ) FORCE_WRITE_AUX = FORCE_WRITE

      ! Quit early if there's no writing to be done
      IF ( .not. KPP_Standalone_ActiveCell%Active_Cell  .AND.               &
           .not. FORCE_WRITE_AUX                       ) THEN
         RETURN
      END IF

      ! Did the call include an optional cell name?
      CELL_NAME_AUX = ''
      IF ( PRESENT( CELL_NAME ) ) CELL_NAME_AUX = CELL_NAME

      ! Get KPP state
      CALL Fun( V    = initC(1:NVAR),                                        &
                F    = initC(NVAR+1:NSPEC),                                  &
                RCT  = localRCONST,                                          &
                Vdot = Vloc,                                                 &
                Aout = Aout                                                 )

      ! Chemistry timestep (seconds)
      DT = GET_TS_CHEM()

      !======================================================================
      ! Write the file.  We need to place this into an !$OMP CRITICAL
      ! block to ensure that only one thread can open & write to the file
      ! at a time.  Otherwise we will get corrupted files
      !======================================================================
      !$OMP CRITICAL

      ! Find a free file LUN
      IU_FILE = findFreeLUN()
      WRITE(level_string,'(I0)') L
      WRITE( YYYYMMDD_hhmmz,'(I0.4,I0.2,I0.2,a,I0.2,I0.2)' )                 &
         Get_Year(), Get_Month(), Get_Day(), '_', Get_Hour(), Get_Minute()

      ! Filename for output
      filename = TRIM( KPP_Standalone_YAML%Output_Directory             ) // &
                 '/'                                                      // &
                 TRIM( Cell_Name_Aux                                    ) // &
                 TRIM( KPP_Standalone_ActiveCell%Active_Cell_Name       ) // &
                 '_L'                                                     // &
                 trim( level_string                                     ) // &
                 '_'                                                      // &
                 TRIM( YYYYMMDD_hhmmz                                   ) // &
                 '.txt'

      ! Open the file
      OPEN( IU_FILE,  FILE=TRIM(filename), ACTION="WRITE",                   &
                      IOSTAT=RC,           ACCESS='SEQUENTIAL')

      ! NOTE: Cannot exit from an !$OMP CRITICAL block, so comment out
      ! for now
      !IF ( RC /= GC_SUCCESS ) THEN
      !   errMsg = 'Error writing chemical state to KPP Standalone file'
      !   CALL GC_Error( errMsg, RC, '' )
      !   RETURN
      !ENDIF

      ! Write header to file
      WRITE( IU_FILE, '(a)' ) '48'
      WRITE( IU_FILE, '(a)' ) REPEAT("=", 76 )
      WRITE( IU_FILE, '(a)' ) ''
      WRITE( IU_FILE, '(a)' ) &
       '                  KPP Standalone Atmospheric Chemical State'
      WRITE( IU_FILE, '(a)' ) 'File Description:'
      WRITE( IU_FILE, '(a)' ) &
       'This file contains model output of the atmospheric chemical state'
      WRITE( IU_FILE, '(a)' ) &
       'as simulated by the GEOS-Chem chemistry module in a 3D setting.'
      WRITE( IU_FILE, '(a)' ) &
       'Each grid cell represents the chemical state of an individual location,'
      WRITE( IU_FILE, '(a)' ) &
       'suitable for input into a separate KPP Standalone program which will'
      WRITE( IU_FILE, '(a)' ) &
       'replicate the chemical evolution of that grid cell for mechanism analysis.'
      WRITE( IU_FILE, '(a)' ) &
       'Note that the KPP Standalone will only use concentrations, rate constants,'
      WRITE( IU_FILE, '(a)' ) &
       'and KPP-specific fields. All other fields are for reference. The first line'
      WRITE( IU_FILE, '(a)' ) &
       'contains the number of lines in this header. If wanting to use this output'
      WRITE( IU_FILE, '(a)' ) &
       'for other analysis, a Python class to read these fields is available by'
      WRITE( IU_FILE, '(a)' ) &
       'request, contact Obin Sturm (psturm@usc.edu).'
      WRITE( IU_FILE, '(a)' ) ''
      WRITE( IU_FILE, '(a)' ) 'Generated by the GEOS-Chem Model'
      WRITE( IU_FILE, '(a)' ) '       (https://geos-chem.org/)'
      WRITE( IU_FILE, '(a)' ) 'Using the KPP Standalone Interface'
      WRITE( IU_FILE, '(a)' ) 'github.com/GEOS-ESM/geos-chem/tree/feature/psturm/kpp_standalone_interface'
      WRITE( IU_FILE, '(a)' ) '     With contributions from:'
      WRITE( IU_FILE, '(a)' ) '        Obin Sturm (psturm@usc.edu)'
      WRITE( IU_FILE, '(a)' ) '        Christoph Keller'
      WRITE( IU_FILE, '(a)' ) '        Michael Long'
      WRITE( IU_FILE, '(a)' ) '        Sam Silva'
      WRITE( IU_FILE, '(a)' ) ''

      ! Write the grid cell metadata as part of the header
      WRITE( IU_FILE, '(a)'       )                                          &
         'Meteorological and general grid cell metadata    '
      WRITE( IU_FILE, '(a,a)'     )                                          &
         'Location:                                        '              // &
          TRIM( CELL_NAME_AUX                              )              // &
          TRIM( KPP_Standalone_ActiveCell%ACTIVE_CELL_NAME )
      WRITE( IU_FILE, '(a,a)'     )                                          &
         'Timestamp:                                       ',                &
          TIMESTAMP_STRING()
      WRITE( IU_FILE, '(a,f11.4)' )                                          &
         'Longitude (degrees):                             ',                &
          State_Grid%XMid(I,J)
      WRITE( IU_FILE, '(a,f11.4)' )                                          &
         'Latitude (degrees):                              ',                &
          State_Grid%YMid(I,J)
      WRITE( IU_FILE, '(a,i6)'    )                                          &
         'GEOS-Chem Vertical Level:                        ',                &
          L
      WRITE( IU_FILE, '(a,f11.4)' )                                          &
         'Pressure (hPa):                                  ',                &
          Get_Pcenter( I, J, L )
      WRITE( IU_FILE, '(a,f11.2)' )                                          &
         'Temperature (K):                                 ',                &
          State_Met%T(I,J,L)
      WRITE( IU_FILE, '(a,e11.4)' )                                          &
         'Dry air density (molec/cm3):                     ',                &
          State_Met%AIRNUMDEN(I,J,L)
      WRITE( IU_FILE, '(a,e11.4)' )                                          &
         'Water vapor mixing ratio (vol H2O/vol dry air):  ',                &
          State_Met%AVGW(I,J,L)
      WRITE( IU_FILE, '(a,e11.4)' )                                          &
         'Cloud fraction:                                  ',                &
          State_Met%CLDF(I,J,L)
      WRITE( IU_FILE, '(a,e11.4)' )                                          &
         'Cosine of solar zenith angle:                    ',                &
          State_Met%SUNCOSmid(I,J)
      WRITE( IU_FILE, '(a)'       )                                          &
         'KPP Integrator-specific parameters               '
      WRITE( IU_FILE, '(a,f11.4)' )                                          &
         'Init KPP Timestep (seconds):                     ',                &
         initHvalue
      WRITE( IU_FILE, '(a,f11.4)' )                                          &
         'Exit KPP Timestep (seconds):                     ',                &
          exitHvalue
      WRITE( IU_FILE, '(a,f11.4)' )                                          &
         'Chemistry operator timestep (seconds):           ',                &
          DT
      WRITE( IU_FILE, '(a,i6)'    )                                          &
         'Number of internal timesteps:                    ',                &
          KPP_TotSteps
      WRITE( IU_FILE, '(a)'       )                                          &
         'CSV data of full chemical state, including species concentrations,'
      WRITE( IU_FILE, '(a)'       )                                          &
         'rate constants (R) and instantaneous reaction rates (A).'
      WRITE( IU_FILE, '(a)'       )                                          &
         'All concentration units are in molecules/cc and rates in molec/cc/s.'
      WRITE( IU_FILE, '(a)'       ) ''
      WRITE( IU_FILE, '(a)'       ) REPEAT("=", 76 )
      WRITE( IU_FILE, '(a)'       ) 'Name,   Value'

      ! Write species concentrations
      DO N = 1, NSPEC
         SpcID = State_Chm%Map_KppSpc(N)
         IF ( SpcID <= 0 ) THEN
            WRITE( IU_FILE, '(a,i0,a,e25.16e3)' ) "C", N, ",", initC(N)
            CYCLE
         ENDIF
         WRITE( IU_FILE, '(a,a,e25.16e3)' )                                   &
            TRIM(State_Chm%SpcData(SpcID)%Info%Name), ',', initC(N)
      ENDDO

      ! Write reaction rates
      DO N = 1, NREACT
         WRITE( IU_FILE,'(a,I0,a,e25.16e3)' ) 'R', N, ',', localRCONST(N)
      ENDDO

      ! Write instantaneous reaction rates
      DO N = 1, NREACT
         WRITE( IU_FILE,'(A,I0,A,E25.16E3)' ) 'A', N, ',', Aout(N)
      ENDDO

      ! Close file
      CLOSE( IU_FILE )
      !$OMP END CRITICAL

   END SUBROUTINE Write_Samples
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cleanup_kpp_standalone
!
! !DESCRIPTION: Deallocates module variables that may have been allocated
!               at run time and unnecessary files required during the process
!\\
!\\
! !INTERFACE:
!
   SUBROUTINE Cleanup_KPP_Standalone( RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE inquireMod,      ONLY : findFreeLUN
!
! !OUTPUT PARAMETERS:
!
    INTEGER, INTENT(OUT) :: RC          ! Success or failure?
!
! !REVISION HISTORY:
!  11 Mar 2024 - P. Obin Sturm - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES::
!
    ! Strings
    CHARACTER(LEN=255) :: arrayId

    ! Assume success
    RC = GC_SUCCESS

    IF ( ALLOCATED( KPP_Standalone_YAML%LocationName ) ) THEN
       arrayId = 'kpp_standalone_interface.F90:KPP_Standalone_YAML%LocationName'
       DEALLOCATE( KPP_Standalone_YAML%LocationName, STAT=RC  )
       CALL GC_CheckVar( arrayId, 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ALLOCATED( KPP_Standalone_YAML%LocationLons ) ) THEN
       arrayId = 'kpp_standalone_interface.F90:KPP_Standalone_YAML%LocationLons'
       DEALLOCATE( KPP_Standalone_YAML%LocationLons, STAT=RC  )
       CALL GC_CheckVar( arrayId, 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ALLOCATED( KPP_Standalone_YAML%LocationLats ) ) THEN
       arrayId = 'kpp_standalone_interface.F90:KPP_Standalone_YAML%LocationLats'
       DEALLOCATE( KPP_Standalone_YAML%LocationLats, STAT=RC  )
       CALL GC_CheckVar( arrayId, 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ALLOCATED( KPP_Standalone_YAML%IDX ) ) THEN
       arrayId = 'kpp_standalone_interface.F90:KPP_Standalone_YAML%IDX'
       DEALLOCATE( KPP_Standalone_YAML%IDX, STAT=RC  )
       CALL GC_CheckVar( arrayId, 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ALLOCATED( KPP_Standalone_YAML%JDX ) ) THEN
       arrayId = 'kpp_standalone_interface.F90:KPP_Standalone_YAML%JDX'
       DEALLOCATE( KPP_Standalone_YAML%JDX, STAT=RC  )
       CALL GC_CheckVar( arrayId, 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ALLOCATED( KPP_Standalone_YAML%Levels ) ) THEN
       arrayId = 'kpp_standalone_interface.F90:KPP_Standalone_YAML%Levels'
       DEALLOCATE( KPP_Standalone_YAML%Levels, STAT=RC  )
       CALL GC_CheckVar( arrayId, 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

   END SUBROUTINE Cleanup_KPP_Standalone
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Find_Number_of_Locations
!
! !DESCRIPTION: Searches a string array containing location names and returns
!  the number of valid locations (i.e. char that do not match MISSING_STR).
!  Assumes all the valid locations will be listed contiguously at the front
!  of the array. Taken from Get_Number_of_Species from input_mod.F90
!\\
!\\
! !INTERFACE:
   FUNCTION Find_Number_of_Locations( a_str ) RESULT( n_valid )
!
! !INPUT PARAMETERS:
!
    CHARACTER(LEN=*), INTENT(IN) :: a_str(:)
!
! !RETURN VALUE:
!
    INTEGER                      :: n_valid
!
! !REVISION HISTORY:
!  11 Mar 2024 - P. Obin Sturm - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: N

    ! Return the number of valid locations
    n_valid = 0
    DO N = 1, SIZE( a_str )
       IF ( TRIM( a_str(N) ) == MISSING_STR ) EXIT
       n_valid = n_valid + 1
    ENDDO

  END FUNCTION Find_Number_of_Locations
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Find_Number_of_Levels
!
! !DESCRIPTION: Searches an integer array containing location names and returns
!  the number of valid levels (i.e. int that do not match MISSING_INT).
!  Assumes all the valid levels will be listed contiguously at the front
!  of the array. Taken from Get_Number_of_Species from input_mod.F90
!\\
!\\
! !INTERFACE:
   FUNCTION Find_Number_of_Levels( a_int ) RESULT( n_valid )
!
! !INPUT PARAMETERS:
!
    INTEGER, INTENT(IN) :: a_int(:)
!
! !RETURN VALUE:
!
    INTEGER             :: n_valid
!  11 Mar 2024 - P. Obin Sturm - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: N

    ! Return the number of valid locations
    n_valid = 0
    DO N = 1, SIZE( a_int )
       IF ( a_int(N) == MISSING_INT ) EXIT
       n_valid = n_valid + 1
    ENDDO

  END FUNCTION Find_Number_of_Levels
!EOC
END MODULE KPP_Standalone_Interface
