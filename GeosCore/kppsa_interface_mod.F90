!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: kppsa_interface_mod.F90
!
! !DESCRIPTION: Contains routines to print the full chemical state
!  which can be used as input to the KPP Standalone.
!\\
!\\
! !INTERFACE:
!
MODULE KppSa_Interface_Mod
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
  PUBLIC :: KppSa_Check_ActiveCell
  PUBLIC :: KppSa_Check_Domain
  PUBLIC :: KppSa_Check_Time
  PUBLIC :: KppSa_Cleanup
  PUBLIC :: KppSa_Config
  PUBLIC :: KppSa_Write_Samples
!
! !DERIVED TYPES:
!
  ! Type to hold information read from the YAML config file
  TYPE, PRIVATE :: KppSa_Interface_Type
     INTEGER                           :: NLOC
     INTEGER                           :: Start_Output(2)
     INTEGER                           :: Stop_Output(2)
     LOGICAL                           :: SkipIt
     LOGICAL                           :: SkipWriteAtThisTime
     CHARACTER(LEN=255)                :: Output_Directory
     CHARACTER(LEN=255),   ALLOCATABLE :: LocationName(:)
     REAL(hp),             ALLOCATABLE :: LocationLons(:)
     REAL(hp),             ALLOCATABLE :: LocationLats(:)
     INTEGER,              ALLOCATABLE :: IDX(:)
     INTEGER,              ALLOCATABLE :: JDX(:)
     INTEGER,              ALLOCATABLE :: Levels(:)
  END TYPE KppSa_Interface_Type

  ! Type to denote active cells
  TYPE, PRIVATE :: KppSa_ActiveCell_Type
     LOGICAL                           :: Active_Cell
     CHARACTER(LEN=255)                :: Active_Cell_Name
  END TYPE KppSa_ActiveCell_Type
!
! !PRIVATE DATA MEMBERS:
!
  TYPE(KppSa_Interface_Type),  PRIVATE :: KppSa_State
  TYPE(KppSa_ActiveCell_Type), PRIVATE :: KppSa_ActiveCell
  !$OMP THREADPRIVATE( KppSa_ActiveCell )
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
! !IROUTINE: kppsa_check_domain
!
! !DESCRIPTION: Subroutine Check_Domain is used to identify if a
!  specified latitude and longitude falls within a grid cell on the
!  current CPU. Multiple lat/lon pairs can be checked simultaneously.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE KppSa_Check_Domain( RC )
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
    IF ( KppSa_State%SkipIt ) RETURN

    ! Compute (I,J) indices of grid boxes
    CALL HCO_GetHorzIJIndex( HcoState,                                       &
                             KppSa_State%NLOC,                               &
                             KppSa_State%LocationLons,                       &
                             KppSa_State%LocationLats,                       &
                             KppSa_State%IDX,                                &
                             KppSa_State%JDX,                                &
                             RC                                             )

  END SUBROUTINE KppSa_Check_Domain
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: kppsa_check_time
!
! !DESCRIPTION: Subroutine Check_Domain is used to identify if a
! specified latitude and longitude falls within a grid cell on the
! current CPU. Multiple lat/lon pairs can be checked simultaneously.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE KppSa_Check_Time( RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE Time_Mod,   ONLY : Get_Nymd, Get_Nhms
!
! !OUTPUT PARAMETERS
!
    INTEGER, INTENT(OUT) :: RC   ! Success or failure?
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
    INTEGER :: yyyymmdd, hhmmss

    ! Initialize
    RC = GC_SUCCESS

    ! Early exit if no locations
    IF ( KppSa_State%SkipIt ) RETURN

    ! Assume we will not write to disk at this date/time
    KppSa_State%SkipWriteAtThisTime  = .TRUE.

    ! Get current date & time
    yyyymmdd = Get_Nymd()
    hhmmss   = Get_Nhms()

    ! Exit if we are outside the window for archiving model state
    IF ( yyyymmdd < KppSa_State%Start_Output(1) ) RETURN
    IF ( yyyymmdd > KppSa_State%Stop_Output(1)  ) RETURN
    IF ( hhmmss   < KppSa_State%Start_Output(2) ) RETURN
    IF ( hhmmss   > KppSa_State%Stop_Output(2)  ) RETURN

    ! If we get this far, we're in the time window where we
    ! archive the chemical state for the KPP standalone
    KppSa_State%SkipWriteAtThisTime = .FALSE.

   END SUBROUTINE KppSa_Check_Time
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: kppsa_check_activecell
!
! !DESCRIPTION: Identifies if a grid cell is within a specified latitude
!  and longitude to print the full chemical state (all concentrations,
!  reaction rates, rate constants, and meteo metadata).
!\\
!\\
! !INTERFACE:
!
   SUBROUTINE KppSa_Check_ActiveCell( I, J, L )
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

    ! Early exit if KPP standalone interface is disabled
    IF ( KppSa_State%SkipIt ) RETURN

    ! Initialize
    KppSa_ActiveCell%Active_Cell      = .FALSE.
    KppSa_ActiveCell%Active_Cell_Name = ''

    ! Skip if we are outside the time interval
    IF ( KppSa_State%SkipWriteAtThisTime ) RETURN

    ! Flag active cells
    IF ( ANY( L == KppSa_State%Levels ) ) THEN
       DO K = 1, KppSa_State%NLOC
          IF ( KppSa_State%IDX(K) == I  .AND.                                &
               KppSa_State%JDX(K) == J ) THEN
             KppSa_ActiveCell%Active_Cell = .TRUE.
             KppSa_ActiveCell%Active_Cell_Name =                             &
                  KppSa_State%LocationName(K)
          ENDIF
       ENDDO
    ENDIF

  END SUBROUTINE KppSa_Check_ActiveCell
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: kppsa_config
!
! !DESCRIPTION: Subroutine Config_KPP_Standalone reads a set of gridcells
!  to be sampled and the full chemical state printed.
!\\
!\\
! !INTERFACE:
!
   SUBROUTINE KppSa_Config( Input_Opt, RC )
!
! !USES:
!
      USE QfYaml_Mod
      USE ErrCode_Mod
      USE Input_Opt_Mod, ONLY : OptInput
      USE RoundOff_Mod,  ONLY : Cast_and_RoundOff
      USE inquireMod,    ONLY : findFreeLUN
!
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
      KppSa_State%SkipIt = .FALSE.
      INQUIRE( FILE=configFile, EXIST=file_exists )
      IF ( .NOT. file_exists ) THEN
         KppSa_State%SkipIt = .TRUE.
         IF ( Input_Opt%amIRoot ) THEN
            WRITE( 6, 100 ) TRIM( configFile )
 100        FORMAT( "Config file ", a ", not found, ",                       &
                    "skipping KPP standalone interface" )
         ENDIF
         RETURN
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
      KppSa_State%SkipIt = ( .not. v_bool )
      IF ( KppSa_State%SkipIt ) THEN
         IF ( Input_Opt%amIRoot ) WRITE( 6, 110 )
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
      KppSa_State%NLOC = Find_Number_of_Locations( a_str )
      IF ( KppSa_State%NLOC .eq. 0 ) THEN
         ! Set SkipIt flag to short circuit other subroutines
         KppSa_State%SkipIt = .TRUE.
         IF ( Input_Opt%amIRoot ) THEN
            WRITE( 6, 120 )
 120        FORMAT( "No active cells for box modeling ",                      &
                    "in kpp_standalone_interface.yml")
            RETURN
         ENDIF
      ENDIF
      ALLOCATE( KppSa_State%LocationName( KppSa_State%NLOC ), STAT=RC )
      CALL GC_CheckVar( 'KppSa_State%LocationName', 0, RC )
      IF ( RC /= GC_SUCCESS ) RETURN
      DO I = 1,KppSa_State%NLOC
         KppSa_State%LocationName(I) = TRIM( a_str(I) )
      END DO

      !========================================================================
      ! Read latitude and longitude of active cells
      !========================================================================

      ! Allocate number of locations for lats and lons
      ALLOCATE( KppSa_State%LocationLons( KppSa_State%NLOC ), STAT=RC )
      CALL GC_CheckVar( 'KppSa_State%LocationLons', 0, RC )
      IF ( RC /= GC_SUCCESS ) RETURN

      ALLOCATE( KppSa_State%LocationLats( KppSa_State%NLOC ), STAT=RC )
      CALL GC_CheckVar( 'KppSa_State%LocationLats', 0, RC )
      IF ( RC /= GC_SUCCESS ) RETURN

      ! Read coordinates
      DO I = 1,KppSa_State%NLOC
      ! Read longitudes
         key = "locations%"//TRIM( KppSa_State%LocationName(I) )//"%longitude"
         v_str = MISSING_STR
         CALL QFYAML_Add_Get( Config, TRIM( key ), v_str, "", RC )
         IF ( RC /= GC_SUCCESS ) THEN
            errMsg = 'Error parsing ' // TRIM( key ) // '!'
            CALL GC_Error( errMsg, RC, thisLoc )
            RETURN
         ENDIF
         KppSa_State%LocationLons( I ) = Cast_and_RoundOff( TRIM( v_str ), places=-1 )
      ! Read latitudes
         key = "locations%"//TRIM( KppSa_State%LocationName(I) )//"%latitude"
         v_str = MISSING_STR
         CALL QFYAML_Add_Get( Config, TRIM( key ), v_str, "", RC )
         IF ( RC /= GC_SUCCESS ) THEN
            errMsg = 'Error parsing ' // TRIM( key ) // '!'
            CALL GC_Error( errMsg, RC, thisLoc )
            RETURN
         ENDIF
         KppSa_State%LocationLats( I ) = Cast_and_RoundOff( TRIM( v_str ), places=-1 )
      END DO

      ! Allocate IDX and JDX (masks for whether a location is on the CPU)
      ALLOCATE( KppSa_State%IDX( KppSa_State%NLOC ), STAT=RC )
      CALL GC_CheckVar( 'KppSa_State%IDX', 0, RC )
      IF ( RC /= GC_SUCCESS ) RETURN

      ALLOCATE( KppSa_State%JDX( KppSa_State%NLOC ), STAT=RC )
      CALL GC_CheckVar( 'KppSa_State%JDX', 0, RC )
      IF ( RC /= GC_SUCCESS ) RETURN

      KppSa_State%IDX(:) = -1
      KppSa_State%JDX(:) = -1

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
      ALLOCATE( KppSa_State%Levels( N ), STAT=RC )
      CALL GC_CheckVar( 'KppSa_State%Levels', 0, RC )
      IF ( RC /= GC_SUCCESS ) RETURN
      DO I = 1,N
         KppSa_State%Levels(I) = a_int(I)
      END DO

      !========================================================================
      ! Get the start & stop date/time for which output will be printed
      !========================================================================
      key = "settings%start_output_at"
      a_int = MISSING_INT
      CALL QFYAML_Add_Get( Config, key, a_int(1:2), "", RC )
      IF ( RC /= GC_SUCCESS ) THEN
         errMsg = 'Error parsing ' // TRIM( key ) // '!'
         CALL GC_Error( errMsg, RC, thisLoc )
         RETURN
      ENDIF
      KppSa_State%Start_Output = a_int(1:2)

      key = "settings%stop_output_at"
      a_int = MISSING_INT
      CALL QFYAML_Add_Get( Config, key, a_int(1:2), "", RC )
      IF ( RC /= GC_SUCCESS ) THEN
         errMsg = 'Error parsing ' // TRIM( key ) // '!'
         CALL GC_Error( errMsg, RC, thisLoc )
         RETURN
      ENDIF
      KppSa_State%Stop_Output = a_int(1:2)

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
         KppSa_State%Output_Directory = "./OutputDir"
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
               KppSa_State%Output_Directory = "./"
            ENDIF
         ENDIF
      ELSE
         KppSa_State%Output_Directory = trim(v_str)
         close(IU_FILE)
      END IF

      !=======================================================================
      ! Print information about sites that will be archived
      !=======================================================================
      IF ( Input_Opt%amIRoot ) THEN
         WRITE( 6, '(a)'   ) REPEAT( "=", 79 )
         WRITE( 6, '(a,/)' ) "KPP STANDALONE INTERFACE"
         WRITE( 6, '(a,/)' ) "Model state will be archived at these sites:"
         DO I = 1, KppSa_State%NLOC
            WRITE( 6, 150 ) KppSa_State%LocationName(I),             &
                            KppSa_State%LocationLons(I),             &
                            KppSa_State%LocationLats(I)
 150        FORMAT( a25, "( ", f9.4, ", ", f9.4, " )")
         ENDDO
         WRITE( 6, '(/,a)'   ) "For GEOS-Chem vertical levels:"
         WRITE( 6, '(100i4)' ) KppSa_State%Levels
         WRITE( 6, 160       ) KppSa_State%Start_Output
 160     FORMAT( "Starting at ", i8.8, 1x, i6.6 )
         WRITE( 6, 170       ) KppSa_State%Stop_Output
 170     FORMAT( "Ending at   ", i8.8, 1x, i6.6 )
         WRITE( 6, '(a)'     ) REPEAT( "=", 79 )
      ENDIF

   END SUBROUTINE KppSa_Config
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: kppsa_write_samples
!
! !DESCRIPTION: Subroutine Write_Samples writes the full chemical state
! (concentrations, reaction rates and rate constants, meteorological
! conditions).
!\\
!\\
! !INTERFACE:
!
   SUBROUTINE KppSa_Write_Samples( I,           J,            L,             &
                                   initC,       localRCONST,  initHvalue,    &
                                   exitHvalue,  ICNTRL,       RCNTRL,        &
                                   State_Grid,  State_Chm,    State_Met,     &
                                   Input_Opt,   KPP_TotSteps, RC,            &
                                   FORCE_WRITE, CELL_NAME                   )
!
! !USES:
!
      USE ErrCode_Mod
      USE State_Grid_Mod,   ONLY : GrdState
      USE State_Chm_Mod,    ONLY : ChmState
      USE State_Met_Mod,    ONLY : MetState
      USE Input_Opt_Mod,    ONLY : OptInput
      USE GcKpp_Global,     ONLY : ATOL
      USE GcKpp_Function
      USE GcKpp_Parameters, ONLY : NSPEC, NREACT, NVAR
      USE TIME_MOD,         ONLY : GET_TS_CHEM
      USE TIME_MOD,         ONLY : TIMESTAMP_STRING
      USE TIME_MOD,         ONLY : Get_Minute
      USE TIME_MOD,         ONLY : Get_Hour
      USE TIME_MOD,         ONLY : Get_Day
      USE TIME_MOD,         ONLY : Get_Month
      USE TIME_MOD,         ONLY : Get_Year
      USE Pressure_Mod,     ONLY : Get_Pcenter
      USE inquireMod,       ONLY : findFreeLUN
!
! !INPUT PARAMETERS:
!
      INTEGER,          INTENT(IN)  :: I                   ! Longitude index
      INTEGER,          INTENT(IN)  :: J                   ! Latitude index
      INTEGER,          INTENT(IN)  :: L                   ! Vertical level
      INTEGER,          INTENT(IN)  :: KPP_TotSteps        ! Total integr. steps
      INTEGER,          INTENT(IN)  :: ICNTRL(20)          ! Integrator options
      TYPE(GrdState),   INTENT(IN)  :: State_Grid          ! Grid State object
      TYPE(ChmState),   INTENT(IN)  :: State_Chm           ! Chem State obj
      TYPE(MetState),   INTENT(IN)  :: State_Met           ! Met State obj
      TYPE(OptInput),   INTENT(IN)  :: Input_Opt           ! Input Options obj
      REAL(dp),         INTENT(IN)  :: initC(NSPEC)        ! Initial conc.
      REAL(dp),         INTENT(IN)  :: localRCONST(NREACT) ! Rate constants
      REAL(dp)                      :: initHvalue          ! Initial timestep
      REAL(dp)                      :: exitHvalue          ! Final timestep:
                                                           !  RSTATE(Nhexit)
      REAL(dp),         INTENT(IN)  :: RCNTRL(20)          ! Integrator options
      LOGICAL,          OPTIONAL    :: FORCE_WRITE         ! Write even if not
                                                           ! in an active cell
      CHARACTER(LEN=*), OPTIONAL    :: CELL_NAME           ! Customize name of
                                                           !  this file
!
! !OUTPUT PARAMETERS:
!
      INTEGER,          INTENT(OUT) :: RC          ! Success or failure
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
      IF ( .not. KppSa_ActiveCell%Active_Cell  .AND.                         &
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
      filename = TRIM( KppSa_State%Output_Directory                     ) // &
                 '/'                                                      // &
                 TRIM( Cell_Name_Aux                                    ) // &
                 TRIM( KppSa_ActiveCell%Active_Cell_Name       ) // &
                 '_L'                                                     // &
                 trim( level_string                                     ) // &
                 '_'                                                      // &
                 TRIM( YYYYMMDD_hhmmz                                   ) // &
                 '.txt'

      ! Open the file
      ! NOTE: We cannot exit from within an !$OMP CRITICAL block
      OPEN( IU_FILE,  FILE=TRIM(filename), ACTION="WRITE",                   &
                      IOSTAT=RC,           ACCESS='SEQUENTIAL')

      ! Write header to file
      WRITE( IU_FILE, '(a)' ) '60'
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
      WRITE( IU_FILE, '(a,/)'     )                                          &
         'Meteorological and general grid cell metadata    '
      WRITE( IU_FILE, '(a,a)'     )                                          &
         'Location:                                        '              // &
          TRIM( CELL_NAME_AUX                     )                       // &
          TRIM( KppSa_ActiveCell%ACTIVE_CELL_NAME )
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
      WRITE( IU_FILE, '(/,a,/)'   )                                          &
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
      WRITE( IU_FILE, '(a)'       ) 'ICNTRL integrator options used:'
      WRITE( IU_FILE, '(10i6)'    ) ICNTRL( 1:10)
      WRITE( IU_FILE, '(10i6)'    ) ICNTRL(11:20)
      WRITE( IU_FILE, '(a)'       ) 'RCNTRL integrator options used:'
      WRITE( IU_FILE, '(5F13.6)'  ) RCNTRL( 1: 5)
      WRITE( IU_FILE, '(5F13.6)'  ) RCNTRL( 6:10)
      WRITE( IU_FILE, '(5F13.6)'  ) RCNTRL(11:15)
      WRITE( IU_FILE, '(5F13.6)'  ) RCNTRL(16:20)
      WRITE( IU_FILE, '(/,a)'     )                                          &
         'CSV data of full chemical state, including species concentrations,'
      WRITE( IU_FILE, '(a)'       )                                          &
         'rate constants (R) and instantaneous reaction rates (A).'
      WRITE( IU_FILE, '(a)'       )                                          &
         'All concentration units are in molec/cm3 and rates in molec/cm3/s.'
      WRITE( IU_FILE, '(a)'       ) ''
      WRITE( IU_FILE, '(a)'       ) REPEAT("=", 76 )
      WRITE( IU_FILE, '(a)'       ) 'Name,   Value,   Absolute Tolerance'

      ! Write species concentrations and absolute tolerances
      DO N = 1, NSPEC
         SpcID = State_Chm%Map_KppSpc(N)
         IF ( SpcID <= 0 ) THEN
            WRITE( IU_FILE, 120 ) N, initC(N), ATOL(N)
 120        FORMAT( "C", i0, ",", es25.16e3, ",", es10.2e2 )
            CYCLE
         ENDIF
         WRITE( IU_FILE, 130 )                                               &
            TRIM(State_Chm%SpcData(SpcID)%Info%Name), initC(N), ATOL(N)
 130     FORMAT( a, ",", es25.16e3, ",", es10.2e2 )
      ENDDO

      ! Write reaction rates
      DO N = 1, NREACT
         WRITE( IU_FILE, 140 ) N, localRCONST(N)
 140     FORMAT( "R", i0, ",", es25.16e3 )
      ENDDO

      ! Write instantaneous reaction rates
      DO N = 1, NREACT
         WRITE( IU_FILE, 150 ) N, Aout(N)
 150     FORMAT( "A", i0, ",", es25.16e3 )
      ENDDO

      ! Close file
      CLOSE( IU_FILE )
      !$OMP END CRITICAL

   END SUBROUTINE KppSa_Write_Samples
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: kppsa_cleanup
!
! !DESCRIPTION: Deallocates module variables that may have been allocated
!               at run time and unnecessary files required during the process
!\\
!\\
! !INTERFACE:
!
   SUBROUTINE KppSa_Cleanup( RC )
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

    IF ( ALLOCATED( KppSa_State%LocationName ) ) THEN
       arrayId = 'kpp_standalone_interface.F90:KppSa_State%LocationName'
       DEALLOCATE( KppSa_State%LocationName, STAT=RC  )
       CALL GC_CheckVar( arrayId, 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ALLOCATED( KppSa_State%LocationLons ) ) THEN
       arrayId = 'kpp_standalone_interface.F90:KppSa_State%LocationLons'
       DEALLOCATE( KppSa_State%LocationLons, STAT=RC  )
       CALL GC_CheckVar( arrayId, 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ALLOCATED( KppSa_State%LocationLats ) ) THEN
       arrayId = 'kpp_standalone_interface.F90:KppSa_State%LocationLats'
       DEALLOCATE( KppSa_State%LocationLats, STAT=RC  )
       CALL GC_CheckVar( arrayId, 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ALLOCATED( KppSa_State%IDX ) ) THEN
       arrayId = 'kpp_standalone_interface.F90:KppSa_State%IDX'
       DEALLOCATE( KppSa_State%IDX, STAT=RC  )
       CALL GC_CheckVar( arrayId, 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ALLOCATED( KppSa_State%JDX ) ) THEN
       arrayId = 'kpp_standalone_interface.F90:KppSa_State%JDX'
       DEALLOCATE( KppSa_State%JDX, STAT=RC  )
       CALL GC_CheckVar( arrayId, 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ALLOCATED( KppSa_State%Levels ) ) THEN
       arrayId = 'kpp_standalone_interface.F90:KppSa_State%Levels'
       DEALLOCATE( KppSa_State%Levels, STAT=RC  )
       CALL GC_CheckVar( arrayId, 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

   END SUBROUTINE KppSa_Cleanup
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
    DO N = 1, SIZE( a_int )
       IF ( a_int(N) == MISSING_INT ) EXIT
       n_valid = n_valid + 1
    ENDDO

  END FUNCTION Find_Number_of_Levels
!EOC
END MODULE KppSa_Interface_Mod
