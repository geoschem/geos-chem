!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: input_mod.F90
!
! !DESCRIPTION: Module INPUT\_MOD contains routines that read the GEOS-Chem
!  input file at the start of the run and pass the information to several
!  other GEOS-Chem F90 modules.
!\\
!\\
! !INTERFACE:
!
MODULE Input_Mod
!
! !USES:
!
  USE CharPak_Mod, ONLY : MaxDim  => MaxStrLen
  USE QfYaml_Mod
  USE Precision_Mod

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC  :: Read_Input_File
  PUBLIC  :: Do_Error_Checks
  PUBLIC  :: Validate_Directories
!
! !REVISION HISTORY:
!  20 Jul 2004 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE TYPES:
!
  ! QFYAML missing values
  LOGICAL,           PARAMETER, PRIVATE :: MISSING_BOOL = .FALSE.
  INTEGER,           PARAMETER, PRIVATE :: MISSING_INT  = -999
  REAL(yp),          PARAMETER, PRIVATE :: MISSING_MW   = 1.0_yp
  REAL(yp),          PARAMETER, PRIVATE :: MISSING_REAL = -999e+0_yp
  CHARACTER(LEN=7),  PARAMETER, PRIVATE :: MISSING_STR  = 'UNKNOWN'

  ! YAML configuration file name to be read
  CHARACTER(LEN=19), PARAMETER, PRIVATE :: configFile   = './input_options.yml'

CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_input_file
!
! !DESCRIPTION: Subroutine READ\_INPUT\_FILE is the driver program for
!  reading the GEOS-Chem input file "input.geos" from disk.
!\\
!\\
! In an ESMF environment, all time steps (chemistry, convection, emissions,
! dynamics) must be specified externally before calling this routine.
! The time steps specified in input.geos are ignored.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Read_Input_File( Input_Opt, State_Grid, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod,  ONLY : OptInput
    USE State_Grid_Mod, ONLY : GrdState
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(INOUT) :: Input_Opt   ! Input options
    TYPE(GrdState), INTENT(INOUT) :: State_Grid  ! Grid State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure
!
! !REVISION HISTORY:
!  20 Jul 2004 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Strings
    CHARACTER(LEN=255) :: thisLoc
    CHARACTER(LEN=512) :: errMsg

    ! Objects
    TYPE(QFYAML_t)     :: Config, ConfigAnchored

    !========================================================================
    ! Read_Input_File begins here!
    !========================================================================

    ! Echo output
    IF ( Input_Opt%amIRoot ) THEN
       WRITE( 6, '(a  )' ) REPEAT( '=', 79 )
       WRITE( 6, '(a,/)' ) 'G E O S - C H E M   U S E R   I N P U T'
       WRITE( 6, 100   ) TRIM( configFile )
100    FORMAT( 'READ_INPUT_FILE: Opening ', a )
    ENDIF


    ! Assume success
    RC      = GC_SUCCESS
    errMsg  = ''
    thisLoc = ' -> at Read_Input_File (in module GeosCore/input_mod.F90)'

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
    ! Get basic simulation settings from the YAML Config object
    !========================================================================

    ! Simulation config settings
    CALL Config_Simulation( Config, Input_Opt, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error in "Config_Simulation"!'
       CALL GC_Error( errMsg, RC, thisLoc  )
       CALL QFYAML_CleanUp( Config         )
       CALL QFYAML_CleanUp( ConfigAnchored )
       RETURN
    ENDIF

#if !(defined( EXTERNAL_GRID ) || defined( EXTERNAL_FORCING ))
    ! Grid config settings
    ! Skip if we are gettomg the grid from an external model
    CALL Config_Grid( Config, Input_Opt, State_Grid, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error in "Config_Grid"!'
       CALL GC_Error( errMsg, RC, thisLoc  )
       CALL QFYAML_CleanUp( Config         )
       CALL QFYAML_CleanUp( ConfigAnchored )
       RETURN
    ENDIF
#endif

    ! Timesteps config settings
    CALL Config_Timesteps( Config, Input_Opt, State_Grid, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error in "Config_Grid"!'
       CALL GC_Error( errMsg, RC, thisLoc  )
       CALL QFYAML_CleanUp( Config         )
       CALL QFYAML_CleanUp( ConfigAnchored )
       RETURN
    ENDIF

    !========================================================================
    ! Get settings for GEOS-Chem operations from the YAML Config object
    ! (read in the same order as we did from input.geos)
    !========================================================================

    ! Transport settings
    CALL Config_Transport( Config, Input_Opt, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error in "Config_Transport"!'
       CALL GC_Error( errMsg, RC, thisLoc )
       CALL QFYAML_CleanUp( Config         )
       CALL QFYAML_CleanUp( ConfigAnchored )
       RETURN
    ENDIF

    ! Convection and PBL mixing settings
    CALL Config_Convection_Mixing( Config, Input_Opt, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error in "Config_Chemistry"!'
       CALL GC_Error( errMsg, RC, thisLoc  )
       CALL QFYAML_CleanUp( Config         )
       CALL QFYAML_CleanUp( ConfigAnchored )
       RETURN
    ENDIF

    ! Aerosol settings
    CALL Config_Aerosol( Config, Input_Opt, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error in "Config_Aerosol"!'
       CALL GC_Error( errMsg, RC, thisLoc  )
       CALL QFYAML_CleanUp( Config         )
       CALL QFYAML_CleanUp( ConfigAnchored )
       RETURN
    ENDIF

    ! Dry deposition and wet deposition settings
    CALL Config_DryDep_WetDep( Config, Input_Opt, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error in "Config_DryDep_WetDep"!'
       CALL GC_Error( errMsg, RC, thisLoc  )
       CALL QFYAML_CleanUp( Config         )
       CALL QFYAML_CleanUp( ConfigAnchored )
       RETURN
    ENDIF

    ! Chemistry settings
    CALL Config_Chemistry( Config, Input_Opt, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error in "Config_Chemistry"!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    ! Photolysis settings
    CALL Config_Photolysis( Config, Input_Opt, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error in "Config_Photolysis"!'
       CALL GC_Error( errMsg, RC, thisLoc  )
       CALL QFYAML_CleanUp( Config         )
       CALL QFYAML_CleanUp( ConfigAnchored )
       RETURN
    ENDIF

    ! RRTMG (radiative transfer model) settings
    CALL Config_RRTMG( Config, Input_Opt, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error in "Config_RRTMG"!'
       CALL GC_Error( errMsg, RC, thisLoc  )
       CALL QFYAML_CleanUp( Config         )
       CALL QFYAML_CleanUp( ConfigAnchored )
       RETURN
    ENDIF

    !========================================================================
    ! Get settings for extra diagnostics from the YAML Config object
    !========================================================================
    ! Obspack diagnostic settings
    CALL Config_ObsPack( Config, Input_Opt, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error in "Config_ObsPack"!'
       CALL GC_Error( errMsg, RC, thisLoc )
       CALL QFYAML_CleanUp( Config         )
       CALL QFYAML_CleanUp( ConfigAnchored )
       RETURN
    ENDIF

#if !(defined( EXTERNAL_GRID ) || defined( EXTERNAL_FORCING ))

    ! Planeflight diagnostic settings
    ! (Skip if we are connecting to an external model)
    CALL Config_PlaneFlight( Config, Input_Opt, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error in "Config_PlaneFlight"!'
       CALL GC_Error( errMsg, RC, thisLoc  )
       CALL QFYAML_CleanUp( Config         )
       CALL QFYAML_CleanUp( ConfigAnchored )
       RETURN
    ENDIF

#ifdef BPCH_DIAG
    ! Planeflight diagnostic settings
    ! (Skip if we are connecting to an external model)
    CALL Config_Gamap( Config, Input_Opt, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error in "Config_Gamap"!'
       CALL GC_Error( errMsg, RC, thisLoc  )
       CALL QFYAML_CleanUp( Config         )
       CALL QFYAML_CleanUp( ConfigAnchored )
       RETURN
    ENDIF
#endif

#endif

!       ELSE IF ( INDEX( LINE, 'CO SIM MENU'     ) > 0 ) THEN
!          CALL READ_CO_SIM_MENU( Input_Opt, RC )
!          IF ( RC /= GC_SUCCESS ) THEN
!             errMsg = 'Error in "Read_CO_Sim_Menu"!'
!             CALL GC_Error( errMsg, RC, thisLoc )
!             RETURN
!          ENDIF
!
!       ELSE IF ( INDEX( LINE, 'CO2 SIM MENU'     ) > 0 ) THEN
!          CALL READ_CO2_SIM_MENU( Input_Opt, RC )
!          IF ( RC /= GC_SUCCESS ) THEN
!             errMsg = 'Error in "Read_CO2_Sim_Menu"!'
!             CALL GC_Error( errMsg, RC, thisLoc )
!             RETURN
!          ENDIF
!
!       ELSE IF ( INDEX( LINE, 'POPS MENU'        ) > 0 ) THEN
!          CALL READ_POPS_MENU( Input_Opt, RC )
!          IF ( RC /= GC_SUCCESS ) THEN
!             errMsg = 'Error in "Read_POPS_Menu"!'
!             CALL GC_Error( errMsg, RC, thisLoc )
!             RETURN
!          ENDIF
!
!       ELSE IF ( INDEX( LINE, 'MERCURY MENU'     ) > 0 ) THEN
!          CALL READ_MERCURY_MENU( Input_Opt, RC  )
!          IF ( RC /= GC_SUCCESS ) THEN
!             errMsg = 'Error in "Read_Mercury_Menu"!'
!             CALL GC_Error( errMsg, RC, thisLoc )
!             RETURN
!          ENDIF
!
!       ELSE IF ( INDEX( LINE, 'CH4 MENU'         ) > 0 ) THEN
!          CALL READ_CH4_MENU( Input_Opt, RC )
!          IF ( RC /= GC_SUCCESS ) THEN
!             errMsg = 'Error in "Read_CH4_Menu"!'
!             CALL GC_Error( errMsg, RC, thisLoc )
!             RETURN
!          ENDIF
!
!       ELSE IF ( INDEX( LINE, 'PASSIVE SPECIES' ) > 0 ) THEN
!          CALL READ_PASSIVE_SPECIES_MENU( Input_Opt, RC )
!          IF ( RC /= GC_SUCCESS ) THEN
!             errMsg = 'Error in "Read_Passive_Species_Menu"!'
!             CALL GC_Error( errMsg, RC, thisLoc )
!             RETURN
!          ENDIF
!
!
!#ifdef BPCH_DIAG
!       !=====================================================================
!       ! Skip BPCH-related menus unless compiled with BPCH_DIAG=y
!       ! Always skip BPCH-related menus for external ESMs
!       !=====================================================================
!       ELSE IF ( INDEX( LINE, 'OUTPUT MENU'      ) > 0 ) THEN
!          CALL READ_OUTPUT_MENU( Input_Opt, RC )
!          IF ( RC /= GC_SUCCESS ) THEN
!             errMsg = 'Error in "Read_Output_Menu"!'
!             CALL GC_Error( errMsg, RC, thisLoc )
!             RETURN
!          ENDIF
!
!       ELSE IF ( INDEX( LINE, 'DIAGNOSTIC MENU'  ) > 0 ) THEN
!          CALL READ_DIAGNOSTIC_MENU( Input_Opt, RC )
!          IF ( RC /= GC_SUCCESS ) THEN
!             errMsg = 'Error in "Read_Diagnostic_Menu"!'
!             CALL GC_Error( errMsg, RC, thisLoc )
!             RETURN
!          ENDIF
!
!       ELSE IF ( INDEX( LINE, 'ND51 MENU'        ) > 0 ) THEN
!          CALL READ_ND51_MENU( Input_Opt, RC )
!          IF ( RC /= GC_SUCCESS ) THEN
!             errMsg = 'Error in "Read_ND51_Menu"!'
!             CALL GC_Error( errMsg, RC, thisLoc )
!             RETURN
!          ENDIF
!
!       ELSE IF ( INDEX( LINE, 'ND51b MENU'       ) > 0 ) THEN
!          CALL READ_ND51b_MENU( Input_Opt, RC )
!          IF ( RC /= GC_SUCCESS ) THEN
!             errMsg = 'Error in "Read_ND51b_Menu"!'
!             CALL GC_Error( errMsg, RC, thisLoc )
!             RETURN
!          ENDIF
!
!#if defined TOMAS
!       ELSE IF ( INDEX( LINE, 'PROD & LOSS MENU' ) > 0 ) THEN
!          CALL READ_PROD_LOSS_MENU( Input_Opt, RC )
!          IF ( RC /= GC_SUCCESS ) THEN
!             errMsg = 'Error in "Read_Prod_Loss_Menu"!'
!             CALL GC_Error( errMsg, RC, thisLoc )
!             RETURN
!          ENDIF
!#endif
!#endif
!#endif
!
!       ELSE IF ( INDEX( LINE, 'END OF FILE'      ) > 0 ) THEN
!          EXIT
!
!       ENDIF
!    ENDDO
!
!    ! Close input file
!    CLOSE( IU_GEOS )
!

    !========================================================================
    ! Check GEOS-CHEM timesteps
    !========================================================================
    CALL Check_Time_Steps( Input_Opt, State_Grid, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error in "Check_Time_Steps"!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !========================================================================
    ! Further error-checking and initialization
    !========================================================================
    CALL QFYAML_CleanUp( Config         )
    CALL QFYAML_CleanUp( ConfigAnchored )

  END SUBROUTINE Read_Input_File
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: config_simulation
!
! !DESCRIPTION: Copies simulation information from the Config object
!  to Input_Opt, and does necessary checks.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Config_Simulation( Config, Input_Opt, RC )
!
! !USES:
!
    USE Charpak_Mod,   ONLY : To_UpperCase
    USE ErrCode_Mod
    USE Input_Opt_Mod, ONLY : OptInput
    USE Time_Mod
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(QFYAML_t), INTENT(INOUT) :: Config      ! YAML Config object
    TYPE(OptInput), INTENT(INOUT) :: Input_Opt   ! Input Options object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    LOGICAL                      :: v_bool
    INTEGER                      :: N,                C
    REAL(fp)                     :: JulianDateStart,  JulianDateEnd
#if defined( ESMF_ ) || defined( MODEL_ )
    INTEGER                      :: H,       M,       S
    REAL(f4)                     :: init_UTC
#endif

    ! Arrays
    INTEGER                      :: a_int(2)

    ! Strings
    CHARACTER(LEN=6)             :: timeStr
    CHARACTER(LEN=8)             :: dateStr
    CHARACTER(LEN=12)            :: met
    CHARACTER(LEN=24)            :: sim
    CHARACTER(LEN=255)           :: thisLoc
    CHARACTER(LEN=512)           :: errMsg
    CHARACTER(LEN=QFYAML_NamLen) :: key
    CHARACTER(LEN=QFYAML_StrLen) :: v_str

    !========================================================================
    ! Parse_Simulation begins here!
    !========================================================================

    ! Initialize
    RC      = GC_SUCCESS
    errMsg  = ''
    thisLoc = &
     ' -> at Read_Simulation_Menu (in module GeosCore/input_mod.F90)'

    !------------------------------------------------------------------------
    ! Simulation start date
    !------------------------------------------------------------------------
    key   = "simulation%start_date"
    a_int = MISSING_INT
    CALL QFYAML_Add_Get( Config, key, a_int, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Input_Opt%NYMDb = a_int(1)
    Input_Opt%NHMSb = a_int(2)

    ! Make sure the starting date NYMDb is valid
    IF ( .not. Valid_Date( Input_Opt%NYMDb ) ) THEN
       WRITE( DateStr, '(i8.8)' ) Input_Opt%NYMDb
       errMsg = 'Input%Opt%NYMDb = ' // DateStr        // &
                ' is not a valid calendar date!'       // &
                ' Please check your "input.geos" file.'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    ! Make sure the starting time NHMSb is valid
    IF ( .not. Valid_Time( Input_Opt%NHMSb ) ) THEN
       WRITE( TimeStr, '(i6.6)' ) Input_Opt%NHMSb
       errMsg = 'Input%Opt%NHMSb = ' // TimeStr        // &
                ' is not a valid clock time!'          // &
                ' Please check your "input.geos" file.'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! Simulation end date
    !------------------------------------------------------------------------
    key   = "simulation%end_date"
    a_int = MISSING_INT
    CALL QFYAML_Add_Get( Config, key, a_int, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Input_Opt%NYMDe = a_int(1)
    Input_Opt%NHMSe = a_int(2)

    ! Make sure the starting date NYMDb is valid
    IF ( .not. Valid_Date( Input_Opt%NYMDe ) ) THEN
       WRITE( DateStr, '(i8.8)' ) Input_Opt%NYMDe
       errMsg = 'Input%Opt%NYMDe = ' // DateStr        // &
                ' is not a valid calendar date!'       // &
                ' Please check your "input.geos" file.'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    ! Make sure the ending time NHMSe is valid
    IF ( .not. Valid_Time( Input_Opt%NHMSe ) ) THEN
       WRITE( TimeStr, '(i6.6)' ) Input_Opt%NHMSe
       errMsg = 'Input%Opt%NHMSe = ' // TimeStr        // &
                ' is not a valid clock time!'          // &
                ' Please check your "input.geos" file.'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    ! Compute the length of the simulation, in elapsed seconds
    JulianDateStart = GET_JD( Input_Opt%NymdB, Input_Opt%NhmsB )
    JulianDateEnd   = GET_JD( Input_Opt%NymdE, Input_Opt%NhmsE )
    Input_Opt%SimLengthSec = NINT( ( JulianDateEnd - JulianDateStart ) &
                             * 86400_f8)

    !------------------------------------------------------------------------
    ! Root data directory
    !------------------------------------------------------------------------
    key   = "simulation%root_data_dir"
    v_str = MISSING_STR
    CALL QFYAML_Add_Get( Config, key, v_str, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Input_Opt%Data_Dir = TRIM( v_str )

    ! Make sure DATA-DIR ends with a "/" character
    C = LEN_TRIM( Input_Opt%DATA_DIR )
    IF ( Input_Opt%DATA_DIR(C:C) /= '/' ) THEN
       Input_Opt%DATA_DIR = TRIM( Input_Opt%DATA_DIR ) // '/'
    ENDIF

    ! Create CHEM_INPUTS directory
    Input_Opt%CHEM_INPUTS_DIR = TRIM( Input_Opt%DATA_DIR ) // &
                                'CHEM_INPUTS/'

    !------------------------------------------------------------------------
    ! Meteorology field
    !------------------------------------------------------------------------
    key   = "simulation%met_field"
    v_str = MISSING_STR
    CALL QFYAML_Add_Get( Config, key, v_str, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Input_Opt%MetField = TRIM( v_str )

    ! Make sure a valid met field is specified
    Met = To_UpperCase( TRIM( Input_Opt%MetField ) )
    SELECT CASE( TRIM( Met ) )
    CASE( 'GEOS-FP', 'GEOSFP' )
       Input_Opt%MetField = 'GEOSFP'
    CASE( 'MERRA-2', 'MERRA2' )
       Input_Opt%MetField = 'MERRA2'
    CASE( 'MODELE2.1' )
       Input_Opt%MetField = 'MODELE2.1'
    CASE( 'MODELE2.2' )
       Input_Opt%MetField = 'MODELE2.2'
    CASE DEFAULT
       errMsg = Trim( Input_Opt%MetField) // ' is not a valid '  // &
                ' met field. Supported met fields are GEOS-FP, '   // &
                ' MERRA-2 and ModelE2.1. Please check your "input.geos" file.'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    END SELECT

    !------------------------------------------------------------------------
    ! Simulation type
    !------------------------------------------------------------------------
    key   = "simulation%name"
    v_str = MISSING_STR
    CALL QFYAML_Add_Get( Config, key, v_str, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Input_Opt%SimulationName = TRIM( v_str )

    ! Error check simulation name
    Sim = To_UpperCase( TRIM( Input_Opt%SimulationName ) )
    IF ( TRIM(Sim) /= 'AEROSOL' .and. TRIM(Sim) /= 'CH4'               .and. &
         TRIM(Sim) /= 'CO2'     .and. TRIM(Sim) /= 'FULLCHEM'          .and. &
         TRIM(Sim) /= 'HG'      .and. TRIM(Sim) /= 'METALS'            .and. &
         TRIM(Sim) /= 'POPS'    .and. TRIM(Sim) /= 'TRANSPORTTRACERS'  .and. &
         TRIM(Sim) /= 'TAGCO'   .and. TRIM(Sim) /= 'TAGCH4'            .and. &
         TRIM(Sim) /= 'TAGHG'   .and. TRIM(Sim) /= 'TAGO3'           ) THEN
       errMsg = Trim( Input_Opt%SimulationName) // ' is not a'            // &
                ' valid simulation. Supported simulations are:'           // &
                ' aerosol, CH4, CO2, fullchem, Hg, POPs,'                 // &
                ' TransportTracers, TagCO, TagCH4, TagHg, or TagO3.'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    ! Set simulation type flags in Input_Opt
    Input_Opt%ITS_A_CH4_SIM        = ( TRIM(Sim) == 'CH4'              .or.  &
                                       TRIM(Sim) == 'TAGCH4'                )
    Input_Opt%ITS_A_CO2_SIM        = ( TRIM(Sim) == 'CO2'                   )
    Input_Opt%ITS_A_FULLCHEM_SIM   = ( TRIM(Sim) == 'FULLCHEM'              )
    Input_Opt%ITS_A_MERCURY_SIM    = ( TRIM(Sim) == 'HG'               .or.  &
                                       TRIM(Sim) == 'TAGHG'                 )
    Input_Opt%ITS_A_POPS_SIM       = ( TRIM(Sim) == 'POPS'                  )
    Input_Opt%ITS_A_RnPbBe_SIM     = ( TRIM(Sim) == 'TRANSPORTTRACERS'      )
    Input_Opt%ITS_A_TAGO3_SIM      = ( TRIM(Sim) == 'TAGO3'                 )
    Input_Opt%ITS_A_TAGCO_SIM      = ( TRIM(Sim) == 'TAGCO'                 )
    Input_Opt%ITS_AN_AEROSOL_SIM   = ( TRIM(Sim) == 'AEROSOL'               )
    Input_Opt%ITS_A_TRACEMETAL_SIM = ( TRIM(SIM) == 'METALS'                )

    !------------------------------------------------------------------------
    ! Species database file
    !------------------------------------------------------------------------
    key   = "simulation%species_database_file"
    v_str = MISSING_STR
    CALL QFYAML_Add_Get( Config, key, v_str, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Input_Opt%SpcDataBaseFile = TRIM( v_str )

    !------------------------------------------------------------------------
    ! Turn on debug output
    !------------------------------------------------------------------------
    key    = "simulation%debug_printout"
    v_bool = MISSING_BOOL
    CALL QFYAML_Add_Get( Config, key, v_bool, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Input_Opt%LPRT = v_bool

    !------------------------------------------------------------------------
    ! Turn on GEOS-Chem timers
    !------------------------------------------------------------------------
    key    = "simulation%use_gcclassic_timers"
    v_bool = MISSING_BOOL
    CALL QFYAML_Add_Get( Config, key, v_bool, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Input_Opt%UseTimers = v_bool

    ! Return success
    RC = GC_SUCCESS

    !========================================================================
    ! Print to screen
    !========================================================================
    IF ( Input_Opt%amIRoot ) THEN
       WRITE( 6, '(/,a)' ) 'SIMULATION SETTINGS'
       WRITE( 6, '(  a)' ) '-------------------'
       WRITE( 6, 100 ) 'Start time of run           : ', &
                        Input_Opt%NYMDb, Input_Opt%NHMSb
       WRITE( 6, 100 ) 'End time of run             : ', &
                        Input_Opt%NYMDe, Input_Opt%NHMSe
       WRITE( 6, 110 ) 'Run directory               : ', &
                        TRIM( Input_Opt%RUN_DIR )
       WRITE( 6, 110 ) 'Data Directory              : ', &
                        TRIM( Input_Opt%DATA_DIR )
       WRITE( 6, 110 ) 'CHEM_INPUTS directory       : ', &
                        TRIM( Input_Opt%CHEM_INPUTS_DIR )
       WRITE( 6, 110 ) 'Meteorology field           : ', &
                        TRIM( Input_Opt%MetField )
       WRITE( 6, 110 ) 'Simulation name             : ', &
                        TRIM( Input_Opt%SimulationName )
       WRITE( 6, 110 ) 'Species database file       : ', &
                        TRIM( Input_Opt%SpcDatabaseFile )
       WRITE( 6, 120 ) 'Turn on debug output        : ', &
                        Input_Opt%LPRT
       WRITE( 6, 120 ) 'Turn on GEOS-Chem timers    : ', &
                        Input_Opt%useTimers
    ENDIF

    ! Format statements
100 FORMAT( A, I8.8, 1X, I6.6 )
110 FORMAT( A, A              )
120 FORMAT( A, L5             )

    !========================================================================
    ! Call setup routines from other GEOS-CHEM modules
    !========================================================================

#if defined( EXTERNAL_GRID ) || defined( EXTERNAL_FORCING )
    !-----------------------------------------------------------------
    !         %%%%%%% GEOS-Chem HP (with ESMF & MPI) %%%%%%%
    !
    ! If we are connecting to the ESMF interface, we need to take
    ! the start & end dates as defined in the ESMF resource file.
    ! (i.e. GEOSCHEMchem_GridComp_mod.rc) instead of those in
    ! input.geos.  This is because the ESMF Clock object needs to be
    ! defined at the highest level (in the driver routine), before
    ! input.geos is ever read.
    !
    ! Therefore, we will assign the start & end date fields (i.e.
    ! Input_Opt%NYMDb, Input_Opt%NYMDe, Input_Opt%NHMSb, and
    ! Input_Opt%NHMSe) in the Gridded Component module file
    ! GEOSCHEMchem_GridComp_Mod.F90 (i.e. two levels higher
    ! in the code).  We don'need to define those fields here, so
    ! we have bracketed this with an #ifdef.
    !-----------------------------------------------------------------

    ! Get the starting UTC time from Input_Opt%NHMSb for use below
    CALL YMD_Extract( Input_Opt%NHMSb, H, M, S )
    init_UTC = ( H + ( M / 60 ) + ( S / 3600 ) )

    ! Pass the values for the start & end times of the simulation directly
    ! to GeosUtil/time_mod.F90 via subroutine ACCEPT_EXTERNAL_DATE_TIME.
    ! (bmy, 12/6/12)
    CALL Accept_External_Date_Time( value_NYMDb = Input_Opt%NYMDb,   &
                                    value_NHMSb = Input_Opt%NHMSb,   &
                                    value_NYMDe = Input_Opt%NYMDe,   &
                                    value_NHMSe = Input_Opt%NHMSe,   &
                                    value_NYMD  = Input_Opt%NYMDb,   &
                                    value_NHMS  = Input_Opt%NHMSb,   &
                                    value_UTC   = init_UTC,          &
                                    RC          = RC                 )
#else
    !------------------------------------------------------------------------
    !         %%%%%%% GEOS-Chem CLASSIC (with OpenMP) %%%%%%%
    !
    ! If we are not using ESMF, then call the traditional GEOS-Chem
    ! timing routines (from GeosUtil/time_mod.F90) to set the start &
    ! end times of the simulation, as well as the current time.
    ! (bmy, 12/6/12)
    !------------------------------------------------------------------------

    ! Set start time of run in "time_mod.F90"
    CALL Set_Begin_Time( Input_Opt%NYMDb, Input_Opt%NHMSb, RC  )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error encountered in "Set_Begin_Time"!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    ! Set end time of run in "time_mod.F90"
    errMsg = 'Error encountered in "Set_Begin_Time"!'
    CALL Set_End_Time( Input_Opt%NYMDe, Input_Opt%NHMSe, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    ! Set the current time
    CALL Set_Current_Time()
#endif

#ifdef BPCH_DIAG
    ! Set the start of the 1st diagnostic interval
    CALL Set_DiagB( GET_TAU() )
#endif

  END SUBROUTINE Config_Simulation
!EOC
#if !(defined( EXTERNAL_GRID ) || defined( EXTERNAL_FORCING ))
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: config_grid
!
! !DESCRIPTION: Copies grid information from the Config object
!  to Input_Opt, and does necessary checks.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Config_Grid( Config, Input_Opt, State_Grid, RC )
!
! !USES:
!
    USE CharPak_Mod,    ONLY : StrSplit
    USE ErrCode_Mod
    USE Input_Opt_Mod,  ONLY : OptInput
    USE RoundOff_Mod,   ONLY : RoundOff
    USE State_Grid_Mod, ONLY : GrdState
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(QFYAML_t), INTENT(INOUT) :: Config      ! YAML Config object
    TYPE(OptInput), INTENT(INOUT) :: Input_Opt   ! Input options
    TYPE(GrdState), INTENT(INOUT) :: State_Grid  ! Grid State
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure
!
! !REVISION HISTORY:
!  20 Oct 2018 - M. Sulprizio- Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    LOGICAL                      :: v_bool
    INTEGER                      :: v_int
    INTEGER                      :: nSubStrs
    INTEGER                      :: N
    INTEGER                      :: C

    ! Arrays
    INTEGER                      :: a_int(4)

    ! Strings
    CHARACTER(LEN=10)            :: xMin_Str, xMax_Str
    CHARACTER(LEN=10)            :: yMin_Str, yMax_Str
    CHARACTER(LEN=255)           :: thisLoc,  nLev
    CHARACTER(LEN=512)           :: errMsg
    CHARACTER(LEN=QFYAML_StrLen) :: key
    CHARACTER(LEN=QFYAML_StrLen) :: v_str

    ! String arrays
    CHARACTER(LEN=255)           :: subStrs(MAXDIM)
    CHARACTER(LEN=QFYAML_StrLen) :: a_str(2)

    !========================================================================
    ! Config_Grid begins here!
    !========================================================================

    ! Initialize
    RC      = GC_SUCCESS
    errMsg  = ''
    thisLoc = ' -> at Read_Grid_Menu (in GeosCore/input_mod.F90)'

    !------------------------------------------------------------------------
    ! Grid resolution
    !------------------------------------------------------------------------
    key   = "grid%resolution"
    v_str = MISSING_STR
    CALL QFYAML_Add_Get( Config, key, v_str, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    State_Grid%GridRes = TRIM( v_str )

    ! Split into two values, separated by 'x'
    CALL StrSplit( TRIM( State_Grid%GridRes ) , 'x', SubStrs, nSubStrs )

    ! Stop with error if there are more than two substrings
    IF ( nSubStrs /= 2 ) THEN
       errMsg = 'Error in extracting delta X and Y values from'    // &
                ' State_Grid%GridRes. Values must be separated by' // &
                ' an x. Please check your "input.geos" file.'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    ! Save the delta X and Y values
    READ( SubStrs(1), '(f10.4)' ) State_Grid%DY
    READ( SubStrs(2), '(f10.4)' ) State_Grid%DX

    !------------------------------------------------------------------------
    ! Level range
    !------------------------------------------------------------------------
    key   = "grid%number_of_levels"
    v_int = MISSING_INT
    CALL QFYAML_Add_Get( Config, key, v_int, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    State_Grid%NZ = v_int

    !------------------------------------------------------------------------
    ! Longitude range
    !------------------------------------------------------------------------
    key   = "grid%longitude%range"
    a_str = MISSING_STR
    CALL QFYAML_Add_Get( Config, key, a_str, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    READ( a_str(1), * ) State_Grid%XMin
    READ( a_str(2), * ) State_Grid%XMax
    State_Grid%XMin = Roundoff( State_Grid%XMin, 4 )
    State_Grid%XMax = Roundoff( State_Grid%XMax, 4 )

    ! Make sure values are in valid rangre
    IF ( State_Grid%XMin >= State_Grid%XMax ) THEN
       WRITE( XMin_Str, '(i10)' ) State_Grid%XMin
       WRITE( XMax_Str, '(i10)' ) State_Grid%XMax
       errMsg = 'Lower lon must be smaller than upper lon: ' // &
                TRIM( XMin_Str ) // ' ' // TRIM( XMax_Str )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! Center longitude on International Date Line?Longitude range
    !------------------------------------------------------------------------
    key    = "grid%longitude%center_at_180"
    v_bool = MISSING_BOOL
    CALL QFYAML_Add_Get( Config, key, v_bool, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    State_Grid%Center180 = v_bool

    !------------------------------------------------------------------------
    ! Latitude range
    !------------------------------------------------------------------------
    key   = "grid%latitude%range"
    a_str = MISSING_STR
    CALL QFYAML_Add_Get( Config, key, a_str, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    READ( a_str(1), * ) State_Grid%YMin
    READ( a_str(2), * ) State_Grid%YMax
    State_Grid%YMin = Roundoff( State_Grid%YMin, 4 )
    State_Grid%YMax = Roundoff( State_Grid%YMax, 4 )

    ! Make sure values are in valid range
    IF ( State_Grid%YMin >= State_Grid%YMax ) THEN
       WRITE( YMin_Str, '(i10)' ) State_Grid%YMin
       WRITE( YMax_Str, '(i10)' ) State_Grid%YMax
       errMsg = 'Lower lat must be smaller than upper lat: ' // &
                TRIM( YMin_Str ) // ' ' // TRIM( YMax_Str )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    ! Restrict latitude values to -90.0 and 90.0
    IF ( State_Grid%YMin < -90.0_fp ) THEN
       WRITE( YMin_Str, '(i10)' ) State_Grid%YMin
       errMsg = 'Lower latitude must be between -90 and 90 degN: ' // &
                TRIM( YMin_Str )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    IF ( State_Grid%YMax > 90.0_fp ) THEN
       WRITE( YMax_Str, '(i10)' ) State_Grid%YMax
       errMsg = 'Upper latitude must be between -90 and 90 degN: ' // &
                TRIM( YMax_Str )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! Use half-sized polar boxes in latitude?
    !------------------------------------------------------------------------
    key    = "grid%latitude%half_size_polar_boxes"
    v_bool = MISSING_BOOL
    CALL QFYAML_Add_Get( Config, key, v_bool, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    State_Grid%HalfPolar = v_bool

    !------------------------------------------------------------------------
    ! Nested grid settings
    !------------------------------------------------------------------------
    key    = "grid%nested_grid_simulation%activate"
    v_bool = MISSING_BOOL
    CALL QFYAML_Add_Get( Config, key, v_bool, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    State_Grid%NestedGrid = v_bool

    IF ( State_Grid%NestedGrid ) THEN
       ! Increase NX by 1
       State_Grid%NX        = State_Grid%NX + 1

       ! For now hardcode HalfPolar to false when using a nested grid
       State_Grid%HalfPolar = .FALSE.
    ENDIF

    !------------------------------------------------------------------------
    ! Nested grid transport offsets
    !------------------------------------------------------------------------
    key   = "grid%nested_grid_simulations%buffer_zone_NSEW"
    a_int = MISSING_INT
    CALL QFYAML_Add_Get( Config, key, a_int, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    State_Grid%NorthBuffer = a_int(1)
    State_Grid%SouthBuffer = a_int(2)
    State_Grid%EastBuffer  = a_int(3)
    State_Grid%WestBuffer  = a_int(4)

    ! Set buffers to zero for global grids
    IF ( .not. State_Grid%NestedGrid ) THEN
       State_Grid%NorthBuffer = 0
       State_Grid%SouthBuffer = 0
       State_Grid%EastBuffer  = 0
       State_Grid%WestBuffer  = 0
    ENDIF

    ! Compute grid horizontal dimensions
    State_Grid%NX =                                                          &
       FLOOR( ( State_Grid%XMax - State_Grid%XMin ) / State_Grid%DX )
    IF ( State_Grid%HalfPolar .and. .not. State_Grid%NestedGrid ) THEN
       State_Grid%NY =                                                       &
          FLOOR( ( State_Grid%YMax - State_Grid%YMin ) / State_Grid%DY ) + 1
    ELSE
       State_Grid%NY =                                                       &
          FLOOR( ( State_Grid%YMax - State_Grid%YMin ) / State_Grid%DY )
    ENDIF

    ! Return success
    RC = GC_SUCCESS

    !========================================================================
    ! Print to screen
    !========================================================================
    IF ( Input_Opt%amIRoot ) THEN
       WRITE( 6, '(/,a)' ) 'GRID SETTINGS'
       WRITE( 6, '(  a)' ) '-------------'
       WRITE( 6, 100 ) 'Grid resolution             : ',                     &
                        TRIM( State_Grid%GridRes )
       WRITE( 6, 110 ) 'Min/max longitude           : ',                     &
                        State_Grid%XMin, State_Grid%XMax
       WRITE( 6, 110 ) 'Min/max latitude            : ',                     &
                        State_Grid%YMin, State_Grid%YMax
       WRITE( 6, 120 ) 'X grid dimension            : ',                     &
                        State_Grid%NX
       WRITE( 6, 120 ) 'Y grid dimension            : ',                     &
                        State_Grid%NY
       WRITE( 6, 120 ) 'Z grid dimension            : ',                     &
                        State_Grid%NZ
       WRITE( 6, 130 ) 'Use half-sized polar boxes? : ',                     &
                        State_Grid%HalfPolar
       WRITE( 6, 130 ) 'Center on Intl Date Line?   : ',                     &
                        State_Grid%Center180
       WRITE( 6, 130 ) 'Is this a nested-grid sim?  : ',                     &
                        State_Grid%NestedGrid
       WRITE( 6, 140 ) ' --> Buffer zone (N S E W ) : ',                     &
                        State_Grid%NorthBuffer,                              &
                        State_Grid%SouthBuffer,                              &
                        State_Grid%EastBuffer,                               &
                        State_Grid%WestBuffer
    ENDIF

    ! Format statements
100 FORMAT( A, A                )
110 FORMAT( A, F10.4, 1X, F10.4 )
120 FORMAT( A, I5               )
130 FORMAT( A, L5               )
140 FORMAT( A, 4I5              )

  END SUBROUTINE Config_Grid
!EOC
#endif
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: config_timesteps
!
! !DESCRIPTION: Copies timestep information from the Config object
!  to Input_Opt, and does necessary checks.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Config_Timesteps( Config, Input_Opt, State_Grid, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod,      ONLY : OptInput
    USE State_Grid_Mod,     ONLY : GrdState
!
! !INPUT PARAMETERS:
!
    TYPE(GrdState), INTENT(IN)    :: State_Grid  ! Grid State object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(QFYAML_t), INTENT(INOUT) :: Config      ! YAML Config object
    TYPE(OptInput), INTENT(INOUT) :: Input_Opt   ! Input options
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER                      :: v_int

    ! Strings
    CHARACTER(LEN=255)           :: thisLoc
    CHARACTER(LEN=512)           :: errMsg
    CHARACTER(LEN=QFYAML_StrLen) :: key

    !========================================================================
    ! Config_Timesteps begins here!
    !========================================================================

    ! Initialize
    RC      = GC_SUCCESS
    errMsg  = 'Error reading the "input.geos" file!'
    thisLoc = ' -> at Read_Timestep_Menu (in module GeosCore/input_mod.F90)'

    !------------------------------------------------------------------------
    ! Transport/convection timestep
    !------------------------------------------------------------------------
    key   = "timesteps%transport_timestep_in_s"
    v_int = MISSING_INT
    CALL QFYAML_Add_Get( Config, key, v_int, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Input_Opt%TS_DYN  = v_int
    Input_Opt%TS_CONV = v_int

    !------------------------------------------------------------------------
    ! Chemistry/emissions timestep
    !------------------------------------------------------------------------
    key   = "timesteps%chemistry_timestep_in_s"
    v_int = MISSING_INT
    CALL QFYAML_Add_Get( Config, key, v_int, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Input_Opt%TS_CHEM = v_int
    Input_Opt%TS_EMIS = v_int

    !------------------------------------------------------------------------
    ! Chemistry/emissions timestep
    !------------------------------------------------------------------------
    key   = "timesteps%radiation_timestep_in_s"
    v_int = MISSING_INT
    CALL QFYAML_Add_Get( Config, key, v_int, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Input_Opt%TS_RAD = v_int

    !========================================================================
    ! Error checks
    !========================================================================

    IF ( Input_Opt%SimLengthSec < Input_Opt%TS_DYN                      .or. &
         Input_Opt%SimLengthSec < Input_Opt%TS_CHEM )                   THEN
       IF ( Input_Opt%amIRoot )                                         THEN
          WRITE( 6, '(a)' ) ''
          WRITE( 6, '(a)' ) 'The length of the simulation is shorter '
          WRITE( 6, '(a)' ) 'than the transport and/or chemistry '
          WRITE( 6, '(a)' ) 'timesteps. Check the settings in '
          WRITE( 6, '(a)' ) 'the "input.geos" file.'
          WRITE( 6, '(a)' ) ''
          WRITE( 6, 100   ) 'Transport/Convection [sec]  : ', &
                             Input_Opt%TS_DYN
          WRITE( 6, 100   ) 'Chemistry/Emissions  [sec]  : ', &
                             Input_Opt%TS_CHEM
          WRITE( 6, 100   ) 'Simulation duration  [sec]  : ', &
                             Input_Opt%SimLengthSec
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
    ENDIF

    IF ( TRIM( Input_Opt%MetField ) == 'MERRA2'                        .and. &
         TRIM( State_Grid%GridRes ) == '0.5x0.625' )                   THEN
       IF ( Input_Opt%ITS_A_CH4_SIM .or. Input_Opt%ITS_A_CO2_SIM )     THEN
          IF ( Input_Opt%TS_DYN > 300 .or. Input_Opt%TS_CHEM > 600 )   THEN
             IF ( Input_Opt%amIRoot ) THEN
                WRITE( 6,'(a)' ) ''
                WRITE( 6,'(a)' ) 'It has been noted that MERRA-2 nested grid'
                WRITE( 6,'(a)' ) ' simulations can have very high species'
                WRITE( 6,'(a)' ) ' concentrations in the stratosphere caused'
                WRITE( 6,'(a)' ) ' by a violation of the CFL condition due to'
                WRITE( 6,'(a)' ) ' strong stratospheric winds. This is'
                WRITE( 6,'(a)' ) ' especially problematic when using total'
                WRITE( 6,'(a)' ) ' column concentrations. To avoid the issue,'
                WRITE( 6,'(a)' ) ' a timestep of 5/10 instead of 10/20 is'
                WRITE( 6,'(a)' ) ' recommended for CH4 and CO2 simulations.'
                WRITE( 6,'(a)' ) ''
                WRITE( 6,'(a)' ) 'You may remove this trap at your own peril,'
                WRITE( 6,'(a)' ) ' by commenting out the call to GC_ERROR in'
                WRITE( 6,'(a)' ) ' GeosCore/input_mod.F90. '
                WRITE( 6,'(a)' ) ''
                WRITE( 6,'(a)' ) 'See the MERRA-2 implementation details page'
                WRITE( 6,'(a)' ) ' on the GEOS-Chem wiki for details'
                CALL GC_Error( errMsg, RC, thisLoc )
                RETURN
             ENDIF
          ENDIF
       ENDIF
    ENDIF

    ! Return success
    RC = GC_SUCCESS

    !========================================================================
    ! Print to screen
    !========================================================================
    IF ( Input_Opt%amIRoot ) THEN
       WRITE( 6, '(/,a)' ) 'TIMESTEP MENU'
       WRITE( 6, '(  a)' ) '---------------'
       WRITE( 6, 100     ) 'Transport/Convection [sec]  : ',                 &
                            Input_Opt%TS_DYN
       WRITE( 6, 100     ) 'Chemistry/Emissions  [sec]  : ',                 &
                            Input_Opt%TS_CHEM
       WRITE( 6, 100     ) 'RRTMG rad transfer   [sec]  : ',                 &
                            Input_Opt%TS_RAD
    ENDIF

    ! Format statements
100 FORMAT( A, I5  )

  END SUBROUTINE Config_Timesteps
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: config_transport
!
! !DESCRIPTION: Copies grid information from the Config object
!  to Input_Opt, and does necessary checks.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Config_Transport( Config, Input_Opt, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod, ONLY : OptInput, Set_Input_Opt_Advect
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(QFYAML_t), INTENT(INOUT) :: Config      ! YAML Config object
    TYPE(OptInput), INTENT(INOUT) :: Input_Opt   ! Input Options object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure
!
! !REVISION HISTORY:
!  20 Jul 2004 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    LOGICAL                      :: v_bool
    INTEGER                      :: N

    ! Arrays
    INTEGER                      :: a_int(3)

    ! Strings
    CHARACTER(LEN=255)           :: thisLoc
    CHARACTER(LEN=512)           :: errMsg
    CHARACTER(LEN=QFYAML_NamLen) :: key

    ! String arrays
    CHARACTER(LEN=14)            :: a_str(QFYAML_MaxArr)

    !========================================================================
    ! Config_Transport begins here!
    !========================================================================

    ! Initialize
    RC      = GC_SUCCESS
    errMsg  = 'Error reading the "input.geos" file!'
    thisLoc = ' -> at Config_Transport (in GeosCore/input_mod.F90)'

    !------------------------------------------------------------------------
    ! Turn on transport?
    !------------------------------------------------------------------------
    key    = "operations%transport%activate"
    v_bool = MISSING_BOOL
    CALL QFYAML_Add_Get( Config, key, v_bool, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Input_Opt%LTRAN = v_bool

    !------------------------------------------------------------------------
    ! Fill negative values generated by TPCORE?
    !------------------------------------------------------------------------
    key    = "operations%transport%fill_negative_values"
    v_bool = MISSING_BOOL
    CALL QFYAML_Add_Get( Config, key, v_bool, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Input_Opt%LFILL = v_bool

    !------------------------------------------------------------------------
    ! IORD, JORD, KORD settings for TPCORE
    !------------------------------------------------------------------------
    key   = "operations%transport%iord_jord_kord"
    a_int = MISSING_INT
    CALL QFYAML_Add_Get( Config, key, a_int, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Input_Opt%TPCORE_IORD = a_int(1)
    Input_Opt%TPCORE_JORD = a_int(2)
    Input_Opt%TPCORE_KORD = a_int(3)

    !------------------------------------------------------------------------
    ! Transported (advected) species list
    !------------------------------------------------------------------------
    key   = "operations%transport%transported_species"
    a_str = MISSING_STR
    CALL QFYAML_Add_Get( Config, key, a_str, "", RC, dynamic_size=.TRUE. )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    ! Find number of transported (advected) species by counting the
    ! number of elements that are not equal to the missing value
    ! Find the number of species
    Input_Opt%N_Advect = 0
    DO N = 1, SIZE( a_str )
       IF ( TRIM( a_str(N) ) == MISSING_STR ) EXIT
       Input_Opt%N_Advect = Input_Opt%N_Advect + 1
    ENDDO

    ! Throw an error if
    IF ( Input_Opt%N_Advect > Input_Opt%Max_AdvectSpc ) THEN
       errMsg = 'Number of advected species exceeds maximum. ' // &
            'This value can be modified in input_opt_mod.F90.'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    ! Copy the advected species names into Input_Opt%Advect_Spc_Name
    Input_Opt%AdvectSpc_Name(1:Input_Opt%N_Advect) = &
         a_str(1:Input_Opt%N_Advect)

    !=================================================================
    ! Print to screen
    !=================================================================
    IF ( Input_Opt%amIRoot ) THEN
       WRITE( 6, '(/,a)' ) 'TRANSPORT SETTINGS'
       WRITE( 6, '(  a)' ) '------------------'
       WRITE( 6, 100     ) 'Turn on transport?          : ', &
                            Input_Opt%LTRAN
       WRITE( 6, 100     ) 'Let TPCORE Fill negatives?  : ', &
                            Input_Opt%LFILL
       WRITE( 6, 110     ) 'IORD, JORD, KORD for TPCORE?: ', &
                            Input_Opt%TPCORE_IORD,           &
                            Input_Opt%TPCORE_JORD,           &
                            Input_Opt%TPCORE_KORD
    ENDIF

    ! FORMAT statements
100 FORMAT( A, L5  )
110 FORMAT( A, 5I5 )

    !=================================================================
    ! Call setup routines from other F90 modules
    !=================================================================

    ! Split into tagged species (turn off for full-chemistry)
    IF ( Input_Opt%ITS_A_FULLCHEM_SIM ) THEN

       ! There are no tagged species for fullchem
       Input_Opt%LSPLIT = .FALSE.

    ELSE

       ! All other simulations: tagged if more than 1 species listed
       Input_Opt%LSPLIT = ( Input_Opt%N_ADVECT > 1 )

    ENDIF

    ! Initialize arrays in Input_Opt that depend on N_ADVECT
    CALL Set_Input_Opt_Advect( Input_Opt, RC )

  END SUBROUTINE Config_Transport
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: config_aerosol
!
! !DESCRIPTION: Copies aerosol information from the Config object
!  to Input_Opt, and does necessary checks.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Config_Aerosol( Config, Input_Opt, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod, ONLY : OptInput
    USE RoundOff_Mod,  ONLY : RoundOff
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(QFYAML_t), INTENT(INOUT) :: Config      ! YAML Config object
    TYPE(OptInput), INTENT(INOUT) :: Input_Opt   ! Input Options object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure
!
! !REMARKS:
!  Move error checks that depend on species indices to the subroutine
!  DO_ERROR_CHECKS.  This is now called from GC_INIT_EXTRA, after the
!  initialization of the species database.
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER                      :: N, T
    INTEGER                      :: v_int
    LOGICAL                      :: v_bool
    REAL(yp)                     :: v_real

    ! Strings
    CHARACTER(LEN=255)           :: thisLoc
    CHARACTER(LEN=255)           :: errMsg
    CHARACTER(LEN=QFYAML_NamLen) :: key
    CHARACTER(LEN=QFYAML_StrLen) :: v_str
    CHARACTER(LEN=QFYAML_StrLen) :: a_str(2)

    !========================================================================
    ! Config_Aerosol begins here!
    !========================================================================

    ! Initialize
    RC      = GC_SUCCESS
    errMsg  = ''
    thisLoc = ' -> at Read_Aerosol_Menu (in module GeosCore/input_mod.F90)'

    !------------------------------------------------------------------------
    ! Use online carbon aerosols?
    !------------------------------------------------------------------------
    key    = "aerosols%carbon%activate"
    v_bool = MISSING_BOOL
    CALL QFYAML_Add_Get( Config, key, v_bool, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Input_Opt%LCARB = v_bool

    !------------------------------------------------------------------------
    ! Use brown carbon aerosols?
    !------------------------------------------------------------------------
    key    = "aerosols%carbon%use_brown_carbon"
    v_bool = MISSING_BOOL
    CALL QFYAML_Add_Get( Config, key, v_bool, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Input_Opt%LBRC = v_bool

    !------------------------------------------------------------------------
    ! Include BC absorption enhancement due to coating?
    !------------------------------------------------------------------------
    key    = "aerosols%carbon%enhance_black_carbon_absorption%activate"
    v_bool = MISSING_BOOL
    CALL QFYAML_Add_Get( Config, key, v_bool, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Input_Opt%LBCAE = v_bool

    !------------------------------------------------------------------------
    ! Define BC absorption enhancement
    !------------------------------------------------------------------------
    key   = "aerosols%carbon%enhance_black_carbon_absorption%hydrophilic"
    v_str = MISSING_STR
    CALL QFYAML_Add_Get( Config, key, v_str, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    READ( v_str, * ) Input_Opt%BCAE_1
    Input_Opt%BCAE_1 = Roundoff( Input_Opt%BCAE_1, 2 )

    !------------------------------------------------------------------------
    ! Define BC absorption enhancement (xnw, 8/24/15)
    !------------------------------------------------------------------------
    key   = "aerosols%carbon%enhance_black_carbon_absorption%hydrophobic"
    v_str = MISSING_STR
    CALL QFYAML_Add_Get( Config, key, v_str, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    READ( v_str, * ) Input_Opt%BCAE_2
    Input_Opt%BCAE_2 = Roundoff( Input_Opt%BCAE_2, 2 )

    !------------------------------------------------------------------------
    ! Use secondary organic aerosols?
    !------------------------------------------------------------------------
    key    = "aerosols%complex_SOA%activate"
    v_bool = MISSING_BOOL
    CALL QFYAML_Add_Get( Config, key, v_bool, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Input_Opt%LSOA = v_bool

    !------------------------------------------------------------------------
    ! Use semi-volatile POA?
    !------------------------------------------------------------------------
    key    = "aerosols%complex_SOA%semivolatile_POA"
    v_bool = MISSING_BOOL
    CALL QFYAML_Add_Get( Config, key, v_bool, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Input_Opt%LSVPOA = v_bool

    !------------------------------------------------------------------------
    ! Use online dust aerosols ?
    !------------------------------------------------------------------------
    key    = "aerosols%dust%activate"
    v_bool = MISSING_BOOL
    CALL QFYAML_Add_Get( Config, key, v_bool, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Input_Opt%LDUST = v_bool

    !------------------------------------------------------------------------
    ! Use SO2 and HNO3 uptake on dust aerosols
    !------------------------------------------------------------------------
    key    = "aerosols%dust%acid_uptake_on_dust"
    v_bool = MISSING_BOOL
    CALL QFYAML_Add_Get( Config, key, v_bool, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Input_Opt%LDSTUP = v_bool

    !------------------------------------------------------------------------
    ! Use online sea-salt aerosols?
    !------------------------------------------------------------------------
    key    = "aerosols%sea_salt%activate"
    v_bool = MISSING_BOOL
    CALL QFYAML_Add_Get( Config, key, v_bool, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Input_Opt%LSSALT = v_bool

    !------------------------------------------------------------------------
    ! Accum mode seasalt radii bin edges [um]
    !------------------------------------------------------------------------
    key   = "aerosols%sea_salt%SALA_radius_bin_in_um"
    a_str = MISSING_STR
    CALL QFYAML_Add_Get( Config, key, a_str, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    READ( a_str(1), * ) Input_Opt%SALA_Redge_um(1)
    READ( a_str(2), * ) Input_Opt%SALA_Redge_um(2)
    Input_Opt%SALA_Redge_um(1) = RoundOff( Input_Opt%SALA_Redge_um(1), 2 )
    Input_Opt%SALA_Redge_um(2) = RoundOff( Input_Opt%SALA_Redge_um(2), 2 )

    !------------------------------------------------------------------------
    ! Coarse mode seasalt radii bin edges [um]
    !------------------------------------------------------------------------
    key   = "aerosols%sea_salt%SALC_radius_bin_in_um"
    a_str = MISSING_STR
    CALL QFYAML_Add_Get( Config, key, a_str, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    READ( a_str(1), * ) Input_Opt%SALC_Redge_um(1)
    READ( a_str(2), * ) Input_Opt%SALC_Redge_um(2)

    !------------------------------------------------------------------------
    ! Use marine organic aerosols?
    !------------------------------------------------------------------------
    key    = "aerosols%sea_salt%marine_organic_aerosols"
    v_bool = MISSING_BOOL
    CALL QFYAML_Add_Get( Config, key, v_bool, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Input_Opt%LMPOA = v_bool

    !------------------------------------------------------------------------
    ! Apply gravitational settling in stratosphere?
    !------------------------------------------------------------------------
    key    = "aerosols%stratosphere%settle_strat_aerosol"
    v_bool = MISSING_BOOL
    CALL QFYAML_Add_Get( Config, key, v_bool, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Input_Opt%LGRAVSTRAT = v_bool

    !------------------------------------------------------------------------
    ! Use solid polar stratospheric clouds (PSCs)?
    !------------------------------------------------------------------------
    key    = "aerosols%stratosphere%polar_strat_clouds%activate"
    v_bool = MISSING_BOOL
    CALL QFYAML_Add_Get( Config, key, v_bool, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Input_Opt%LSOLIDPSC = v_bool

    !------------------------------------------------------------------------
    ! Perform heterogeneous chemistry on PSCs?
    !------------------------------------------------------------------------
    key    = "aerosols%stratosphere%polar_strat_clouds%het_chem"
    v_bool = MISSING_BOOL
    CALL QFYAML_Add_Get( Config, key, v_bool, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Input_Opt%LPSCCHEM = v_bool

    !------------------------------------------------------------------------
    ! Allow homogeneous NAT?
    !------------------------------------------------------------------------
    key    = "aerosols%stratosphere%allow_homogeneous_NAT"
    v_bool = MISSING_BOOL
    CALL QFYAML_Add_Get( Config, key, v_bool, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Input_Opt%LHOMNUCNAT = v_bool

    !------------------------------------------------------------------------
    ! NAT supercooling requirement (K)
    !------------------------------------------------------------------------
    key   = "aerosols%stratosphere%NAT_supercooling_req_in_K"
    v_str = MISSING_STR
    CALL QFYAML_Add_Get( Config, key, v_str, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    READ( v_str, * ) Input_Opt%T_NAT_SUPERCOOL
    Input_Opt%T_NAT_SUPERCOOL = RoundOff( Input_Opt%T_NAT_SUPERCOOL, 2 )

    !------------------------------------------------------------------------
    ! Ice supersaturation ratio requirement
    !------------------------------------------------------------------------
    key   = "aerosols%stratosphere%supersat_factor_req_for_ice_nucl"
    v_str = MISSING_STR
    CALL QFYAML_Add_Get( Config, key, v_str, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    READ( v_str, * ) Input_Opt%P_ICE_SUPERSAT
    Input_Opt%P_ICE_SUPERSAT = RoundOff( Input_Opt%P_ICE_SUPERSAT, 2 )

    !------------------------------------------------------------------------
    ! Include stratospheric aerosols optical depths?
    !------------------------------------------------------------------------
    key    = "aerosols%stratosphere%calc_strat_aero_optdepth"
    v_bool = MISSING_BOOL
    CALL QFYAML_Add_Get( Config, key, v_bool, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Input_Opt%LSTRATOD = v_bool

    !------------------------------------------------------------------------
    ! Use online sulfate aerosols?
    !------------------------------------------------------------------------
    key    = "aerosols%sulfate%activate"
    v_bool = MISSING_BOOL
    CALL QFYAML_Add_Get( Config, key, v_bool, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Input_Opt%LSULF = v_bool

    !------------------------------------------------------------------------
    ! Use metal catalyzed oxidation of SO2?
    !------------------------------------------------------------------------
    key    = "aerosols%sulfate%metal_cat_SO2_oxidation"
    v_bool = MISSING_BOOL
    CALL QFYAML_Add_Get( Config, key, v_bool, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Input_Opt%LMETALCATSO2 = v_bool

    !=================================================================
    ! Error checks
    !=================================================================

    ! Make sure that SALA, SALC bins are contiguous
    IF ( Input_Opt%SALA_REDGE_um(2) /= &
         Input_Opt%SALC_REDGE_um(1)     ) THEN
       errMsg = 'SALA and SALC bin edges are not contiguous!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    ! Turn off switches for simulations that don't use aerosols
    IF ( ( .not. Input_Opt%ITS_A_FULLCHEM_SIM )  .and. &
         ( .not. Input_Opt%ITS_AN_AEROSOL_SIM ) ) THEN
       Input_Opt%LSULF        = .FALSE.
       Input_Opt%LMETALCATSO2 = .FALSE.
       Input_Opt%LCARB        = .FALSE.
       Input_Opt%LBRC         = .FALSE.
       Input_Opt%LSOA         = .FALSE.
       Input_Opt%LDUST        = .FALSE.
       Input_Opt%LSSALT       = .FALSE.
       Input_Opt%LMPOA        = .FALSE.
       Input_Opt%LSVPOA       = .FALSE.
       Input_Opt%LBCAE        = .FALSE.
    ENDIF

    ! Return success
    RC = GC_SUCCESS

    !=================================================================
    ! Print to screen
    !=================================================================
    IF ( Input_Opt%amIRoot ) THEN
       WRITE( 6, '(/,a)' ) 'AEROSOL SETTINGS'
       WRITE( 6, '(  a)' ) '----------------'
       WRITE( 6, 100     ) 'Online SULFATE AEROSOLS?    : ', &
                            Input_Opt%LSULF
       WRITE( 6, 100     ) 'Metal catalyzed SO2 ox.?    : ', &
                            Input_Opt%LMETALCATSO2
       WRITE( 6, 100     ) 'Online CARBON AEROSOLS?     : ', &
                            Input_Opt%LCARB
       WRITE( 6, 100     ) 'Brown Carbon Aerosol?       : ', &
                            Input_Opt%LBRC
       WRITE( 6, 100     ) 'BC Absorption Enhancement?  : ', &
                            Input_Opt%LBCAE
       WRITE( 6, 105     ) 'Hydrophilic BC AE factor    : ', &
                            Input_Opt%BCAE_1
       WRITE( 6, 105     ) 'Hydrophobic BC AE factor    : ', &
                            Input_Opt%BCAE_2
       WRITE( 6, 100     ) 'Online COMPLEX SOA?         : ', &
                            Input_Opt%LSOA
       WRITE( 6, 100     ) 'Semivolatile POA?           : ', &
                            Input_Opt%LSVPOA
       WRITE( 6, 100     ) 'Online DUST AEROSOLS?       : ', &
                            Input_Opt%LDUST
       WRITE( 6, 100     ) 'Acid uptake on dust?        : ', &
                            Input_Opt%LDSTUP
       WRITE( 6, 100     ) 'Online SEA SALT AEROSOLS?   : ', &
                            Input_Opt%LSSALT
       WRITE( 6, 110     ) 'Accum  SEA SALT radii [um]  : ', &
                            Input_Opt%SALA_REDGE_um(1),      &
                            Input_Opt%SALA_REDGE_um(2)
       WRITE( 6, 110     ) 'Coarse SEA SALT radii [um]  : ', &
                            Input_Opt%SALC_REDGE_um(1),      &
                            Input_Opt%SALC_REDGE_um(2)
       WRITE( 6, 100     ) 'MARINE ORGANIC AEROSOLS?    : ', &
                            Input_Opt%LMPOA
       WRITE( 6, 100     ) 'Settle strat. aerosols?     : ', &
                            Input_Opt%LGRAVSTRAT
       WRITE( 6, 100     ) 'Online SOLID PSC aerosols?  : ', &
                            Input_Opt%LSOLIDPSC
       WRITE( 6, 100     ) 'Allow hom. NAT nucleation?  : ', &
                            Input_Opt%LHOMNUCNAT
       WRITE( 6, 120     ) 'NAT supercooling requirement: ', &
                            Input_Opt%T_NAT_SUPERCOOL
       WRITE( 6, 120     ) 'Ice supersaturation req.    : ', &
                            ((Input_Opt%P_ICE_SUPERSAT-1)*1.e+2_fp)
       WRITE( 6, 100     ) 'Perform PSC het. chemistry? : ', &
                            Input_Opt%LPSCCHEM
       WRITE( 6, 100     ) 'Use strat. aerosol OD?      : ', &
                            Input_Opt%LSTRATOD
    ENDIF

100 FORMAT( A, L5                )
105 FORMAT( A, f8.2              )
110 FORMAT( A, f8.2, ' - ', f8.2 )
120 FORMAT( A, f8.2, 'K'         )

  END SUBROUTINE Config_Aerosol
!EOC
!!------------------------------------------------------------------------------
!!                  GEOS-Chem Global Chemical Transport Model                  !
!!------------------------------------------------------------------------------
!!BOP
!!
!! !IROUTINE: read_co_sim_menu
!!
!! !DESCRIPTION: Subroutine READ\_CO\_SIM\_MENU reads the CO SIM MENU
!!  section of the GEOS-Chem input file.
!!\\
!!\\
!! !INTERFACE:
!!
!  SUBROUTINE READ_CO_SIM_MENU( Input_Opt, RC )
!!
!! !USES:
!!
!    USE ErrCode_Mod
!    USE Input_Opt_Mod, ONLY : OptInput
!!
!! !INPUT/OUTPUT PARAMETERS:
!!
!    TYPE(OptInput), INTENT(INOUT) :: Input_Opt   ! Input options
!!
!! !OUTPUT PARAMETERS:
!!
!    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure
!!
!! !REVISION HISTORY:
!!  24 Mar 2017 - J. Fisher - Initial version
!!  See https://github.com/geoschem/geos-chem for complete history
!!EOP
!!------------------------------------------------------------------------------
!!BOC
!!
!! !LOCAL VARIABLES:
!!
!    ! Scalars
!    INTEGER            :: N
!    CHARACTER(LEN=255) :: MSG, LOC
!
!    ! Arrays
!    CHARACTER(LEN=255) :: SUBSTRS(MAXDIM)
!
!    !=================================================================
!    ! READ_CO_SIM_MENU begins here!
!    !=================================================================
!
!    ! Assume success
!    RC = GC_SUCCESS
!
!    ! Location for error messages
!    LOC = ' -> at READ_CO_SIM_MENU (in GeosCore/input_mod.F90)'
!
!    ! Error check
!    IF ( CT1 /= 2 ) THEN
!       MSG = 'SIMULATION MENU & ADVECTED SPECIES MENU ' // &
!             'must be read in first!'
!       CALL GC_Error( Msg, RC, Loc )
!       RETURN
!    ENDIF
!
!    ! Use P(CO) from CH4 from full chemistry simulation?
!    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'LPCO_CH4', RC )
!    READ( SUBSTRS(1:N), * ) Input_Opt%LPCO_CH4
!
!    ! Use P(CO) from CH4 from full chemistry simulation?
!    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'LPCO_NMVOC', RC )
!    READ( SUBSTRS(1:N), * ) Input_Opt%LPCO_NMVOC
!
!    !=================================================================
!    ! Print to screen
!    !=================================================================
!    IF ( Input_Opt%ITS_A_TAGCO_SIM .and. Input_Opt%amIRoot ) THEN
!       WRITE(6,'(/,a)') 'CO SIMULATION MENU ' // &
!             '(overwrites any other settings related to CO)'
!       WRITE(6,'(  a)') '-------------------------------------'
!       WRITE(6,100    ) 'Use full chem. P(CO) from CH4?:', &
!                        Input_Opt%LPCO_CH4
!       WRITE(6,100    ) '                   from NMVOC?:', &
!                        Input_Opt%LPCO_NMVOC
!       WRITE(6,'(  a)') '-------------------------------------'
!    ENDIF
!
!    ! FORMAT statements
!90  FORMAT( A )
!100 FORMAT( A, L5 )
!110 FORMAT( A, L5, A )
!
!  END SUBROUTINE READ_CO_SIM_MENU
!!EOC
!!------------------------------------------------------------------------------
!!                  GEOS-Chem Global Chemical Transport Model                  !
!!------------------------------------------------------------------------------
!!BOP
!!
!! !IROUTINE: read_co2_sim_menu
!!
!! !DESCRIPTION: Subroutine READ\_CO2\_SIM\_MENU reads the CO2 SIM MENU
!!  section of the GEOS-Chem input file.
!!\\
!!\\
!! !INTERFACE:
!!
!  SUBROUTINE READ_CO2_SIM_MENU( Input_Opt, RC )
!!
!! !USES:
!!
!    USE ErrCode_Mod
!    USE Input_Opt_Mod,      ONLY : OptInput
!!
!! !INPUT/OUTPUT PARAMETERS:
!!
!    TYPE(OptInput), INTENT(INOUT) :: Input_Opt   ! Input options
!!
!! !OUTPUT PARAMETERS:
!!
!    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure
!!
!! !REVISION HISTORY:
!!  02 Mar 2009 - R. Nassar   - Initial version
!!  See https://github.com/geoschem/geos-chem for complete history
!!EOP
!!------------------------------------------------------------------------------
!!BOC
!!
!! !LOCAL VARIABLES:
!!
!    ! Scalars
!    INTEGER            :: N
!
!    ! Strings
!    CHARACTER(LEN=255) :: errMsg, thisLoc
!
!    ! Arrays
!    CHARACTER(LEN=255) :: SUBSTRS(MAXDIM)
!
!    !=================================================================
!    ! READ_CO2_SIM_MENU begins here!
!    !=================================================================
!
!    ! Initialize
!    RC      = GC_SUCCESS
!    errMsg  = 'Error reading the "input.geos" file!'
!    thisLoc = ' -> at Read_CO2_SIM_Menu (in module GeosCore/input_mod.F90)'
!
!    ! Error check
!    IF ( CT1 /= 2 ) THEN
!       errMsg = 'SIMULATION MENU & ADVECTED SPECIES MENU ' // &
!                'must be read in first!'
!       CALL GC_Error( errMsg, RC, thisLoc )
!       RETURN
!    ENDIF
!
!    ! Use Fossil Fuel emissions?
!    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'LFOSSIL', RC )
!    IF ( RC /= GC_SUCCESS ) THEN
!       CALL GC_Error( errMsg, RC, thisLoc )
!       RETURN
!    ENDIF
!    READ( SUBSTRS(1:N), * ) Input_Opt%LFOSSIL
!
!    ! Use Ocean Exchange?
!    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'LOCEAN', RC )
!    IF ( RC /= GC_SUCCESS ) THEN
!       CALL GC_Error( errMsg, RC, thisLoc )
!       RETURN
!    ENDIF
!    READ( SUBSTRS(1:N), * ) Input_Opt%LOCEAN
!
!    ! Turn on (balanced) biosphere with diurnal cycle?
!    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'LBIODIURNAL', RC )
!    IF ( RC /= GC_SUCCESS ) THEN
!       CALL GC_Error( errMsg, RC, thisLoc )
!       RETURN
!    ENDIF
!    READ( SUBSTRS(1:N), * ) Input_Opt%LBIODIURNAL
!
!    ! Use Net Terrestrial Exchange Climatology?
!    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'LBIONETCLIM', RC )
!    IF ( RC /= GC_SUCCESS ) THEN
!       CALL GC_Error( errMsg, RC, thisLoc )
!       RETURN
!    ENDIF
!    READ( SUBSTRS(1:N), * ) Input_Opt%LBIONETCLIM
!
!    ! Turn on Ship emissions?
!    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'LSHIP', RC )
!    IF ( RC /= GC_SUCCESS ) THEN
!       CALL GC_Error( errMsg, RC, thisLoc )
!       RETURN
!    ENDIF
!    READ( SUBSTRS(1:N), * ) Input_Opt%LSHIP
!
!    ! Turn on Aviation emissions?
!    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'LPLANE', RC )
!    IF ( RC /= GC_SUCCESS ) THEN
!       CALL GC_Error( errMsg, RC, thisLoc )
!       RETURN
!    ENDIF
!    READ( SUBSTRS(1:N), * ) Input_Opt%LPLANE
!
!    ! Turn on CO2 3D chemical source and surface correction?
!    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'LCHEMCO2', RC )
!    IF ( RC /= GC_SUCCESS ) THEN
!       CALL GC_Error( errMsg, RC, thisLoc )
!       RETURN
!    ENDIF
!    READ( SUBSTRS(1:N), * ) Input_Opt%LCHEMCO2
!
!    ! Tagged CO2 Separator line
!    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'separator 1', RC )
!    IF ( RC /= GC_SUCCESS ) THEN
!       CALL GC_Error( errMsg, RC, thisLoc )
!       RETURN
!    ENDIF
!
!    ! Background CO2 (no emissions or exchange) for Tagged-CO2 runs
!    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'LFFBKGRD',RC )
!    IF ( RC /= GC_SUCCESS ) THEN
!       CALL GC_Error( errMsg, RC, thisLoc )
!       RETURN
!    ENDIF
!    READ( SUBSTRS(1:N), * ) Input_Opt%LFFBKGRD
!
!    ! Turn on biosphere and ocean exchange region tagged species?
!    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'LBIOSPHTAG', RC )
!    IF ( RC /= GC_SUCCESS ) THEN
!       CALL GC_Error( errMsg, RC, thisLoc )
!       RETURN
!    ENDIF
!    READ( SUBSTRS(1:N), * ) Input_Opt%LBIOSPHTAG
!
!    ! Turn on fossil fuel emission region tagged species?
!    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'LFOSSILTAG', RC )
!    IF ( RC /= GC_SUCCESS ) THEN
!       CALL GC_Error( errMsg, RC, thisLoc )
!       RETURN
!    ENDIF
!    READ( SUBSTRS(1:N), * ) Input_Opt%LFOSSILTAG
!
!    ! Turn on global ship emissions tagged species?
!    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'LSHIPTAG', RC )
!    IF ( RC /= GC_SUCCESS ) THEN
!       CALL GC_Error( errMsg, RC, thisLoc )
!       RETURN
!    ENDIF
!    READ( SUBSTRS(1:N), * ) Input_Opt%LSHIPTAG
!
!    ! Turn on global aircraft emissions tagged species?
!    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'LPLANETAG', RC )
!    IF ( RC /= GC_SUCCESS ) THEN
!       CALL GC_Error( errMsg, RC, thisLoc )
!       RETURN
!    ENDIF
!    READ( SUBSTRS(1:N), * ) Input_Opt%LPLANETAG
!
!    !=================================================================
!    ! Print to screen
!    !=================================================================
!    IF ( Input_Opt%ITS_A_CO2_SIM .and. Input_Opt%amIRoot ) THEN
!       WRITE(6,'(/,a)') 'CO2 SIMULATION MENU ' // &
!               '(overwrites any other settings related to CO2)'
!       WRITE(6,'(  a)') '-------------------------------------'
!       WRITE(6,100    ) 'National Fossil Fuel Emission :', &
!                         Input_Opt%LFOSSIL
!       WRITE(6,100    ) 'Ocean CO2 Uptake/Emission     :', &
!                         Input_Opt%LOCEAN
!       WRITE(6,100    ) 'Biosphere seas/diurnal cycle  :', &
!                         Input_Opt%LBIODIURNAL
!       WRITE(6,100    ) 'Net Terr Exch - Climatology   :', &
!                         Input_Opt%LBIONETCLIM
!       WRITE(6,100    ) 'Intl/Domestic Ship emissions  :', &
!                         Input_Opt%LSHIP
!       WRITE(6,100    ) 'Intl/Domestic Aviation emiss  :', &
!                         Input_Opt%LPLANE
!       WRITE(6,100    ) 'CO2 from oxidation (CO,CH4,..):', &
!                         Input_Opt%LCHEMCO2
!       WRITE(6, 90    ) 'Tagged CO2 settings'
!       WRITE(6,100    ) '  Save Fossil CO2 in Bckgrnd  :', &
!                         Input_Opt%LFFBKGRD
!       WRITE(6,100    ) '  Tag Biosphere/Ocean CO2     :', &
!                           Input_Opt%LBIOSPHTAG
!       WRITE(6,100    ) '  Tag Fossil Fuel CO2         :', &
!                           Input_Opt%LFOSSILTAG
!       WRITE(6,100    ) '  Tag Global Ship CO2         :', &
!                           Input_Opt%LSHIPTAG
!       WRITE(6,100    ) '  Tag Global Aviation CO2     :', &
!                           Input_Opt%LPLANETAG
!       WRITE(6,'(  a)') '-------------------------------------'
!    ENDIF
!
!    ! FORMAT statements
!90  FORMAT( A )
!100 FORMAT( A, L5 )
!
!  END SUBROUTINE READ_CO2_SIM_MENU
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: config_chemistry
!
! !DESCRIPTION: Copies chemistry information from the Config object
!  to Input_Opt, and does necessary checks.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Config_Chemistry( Config, Input_Opt, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod, ONLY : OptInput
    USE RoundOff_Mod,  ONLY : RoundOff
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(QFYAML_t), INTENT(INOUT) :: Config      ! YAML Config object
    TYPE(OptInput), INTENT(INOUT) :: Input_Opt   ! Input Options object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER                      :: N
    LOGICAL                      :: v_bool

    ! Strings
    CHARACTER(LEN=255)           :: thisLoc
    CHARACTER(LEN=512)           :: errMsg
    CHARACTER(LEN=QFYAML_NamLen) :: key
    CHARACTER(LEN=QFYAML_StrLen) :: v_str

    !========================================================================
    ! Config_Chemistry begins here!
    !========================================================================

    ! Initialize
    RC      = GC_SUCCESS
    errMsg  = 'Error reading the "input.geos" file!'
    thisLoc = ' -> at Config_Chemistry (in module GeosCore/input_mod.F90)'

    !------------------------------------------------------------------------
    ! Turn on chemistry?
    !------------------------------------------------------------------------
    key    = "operations%chemistry%activate"
    v_bool = MISSING_BOOL
    CALL QFYAML_Add_Get( Config, key, v_bool, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Input_Opt%LCHEM = v_bool

    !------------------------------------------------------------------------
    ! Turn on linearized chemistry above chemistry grid?
    !------------------------------------------------------------------------
    key    = "operations%chemistry%linear_chemistry_aloft%activate"
    v_bool = MISSING_BOOL
    CALL QFYAML_Add_Get( Config, key, v_bool, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Input_Opt%LINEAR_CHEM = v_bool

    !------------------------------------------------------------------------
    ! Use Linoz for ozone above chemistry grid? (Otherwise, Synoz is used)
    !------------------------------------------------------------------------
    key    = "operations%chemistry%linear_chemistry_aloft%use_linoz_for_O3"
    v_bool = MISSING_BOOL
    CALL QFYAML_Add_Get( Config, key, v_bool, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Input_Opt%LLINOZ = v_bool
    IF ( .not. Input_Opt%LLINOZ ) Input_Opt%LSYNOZ = .TRUE.

    !------------------------------------------------------------------------
    ! Turn on online stratospheric H2O?
    !------------------------------------------------------------------------
    key    = "operations%chemistry%active_strat_H2O%activate"
    v_bool = MISSING_BOOL
    CALL QFYAML_Add_Get( Config, key, v_bool, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Input_Opt%LACTIVEH2O = v_bool

    !------------------------------------------------------------------------
    ! Turn on online stratospheric H2O?
    !------------------------------------------------------------------------
    key    = "operations%chemistry%active_strat_H2O%use_static_bnd_cond"
    v_bool = MISSING_BOOL
    CALL QFYAML_Add_Get( Config, key, v_bool, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Input_Opt%LStaticH2OBC = v_bool

    !------------------------------------------------------------------------
    ! GAMMA HO2 ?
    !------------------------------------------------------------------------
    key   = "operations%chemistry%gamma_HO2"
    v_str = MISSING_STR
    CALL QFYAML_Add_Get( Config, key, v_str, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    READ( v_str, * ) Input_Opt%GAMMA_HO2
    Input_Opt%GAMMA_HO2 = RoundOff( Input_Opt%GAMMA_HO2, 2 )

    ! Return success
    RC = GC_SUCCESS

    !=================================================================
    ! Print to screen
    !=================================================================
    IF ( Input_Opt%amIRoot ) THEN
       WRITE( 6, '(/,a)' ) 'CHEMISTRY SETTINGS'
       WRITE( 6, '(  a)' ) '------------------'
       WRITE( 6, 100     ) 'Turn on chemistry?          : ', &
                            Input_Opt%LCHEM
       WRITE( 6, 100     ) 'Use linear. mesospheric chem: ', &
                            Input_Opt%LINEAR_CHEM
       WRITE( 6, 100     ) ' => Use Linoz for O3?       : ', &
                            Input_Opt%LLINOZ
       WRITE( 6, 100     ) 'Online strat. H2O?          : ', &
                            Input_Opt%LACTIVEH2O
       WRITE( 6, 100     ) 'Use robust strat H2O BC?    : ', &
                            Input_Opt%LStaticH2OBC
       WRITE( 6, 110     ) 'GAMMA HO2                   : ', &
                            Input_Opt%GAMMA_HO2
    ENDIF

    ! FORMAT statements
100 FORMAT( A, L5  )
110 FORMAT( A, F4.2 )

  END SUBROUTINE Config_Chemistry
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Config_RRTMG
!
! !DESCRIPTION: Copies RRTMG information from the Config object
!  to Input_Opt, and does necessary checks.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Config_RRTMG( Config, Input_Opt, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod, ONLY : OptInput
    USE RoundOff_Mod,  ONLY : RoundOff
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(QFYAML_t), INTENT(INOUT) :: Config      ! YAML Config object
    TYPE(OptInput), INTENT(INOUT) :: Input_Opt   ! Input Options object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure
!
! !REMARKS:
!  Flux outputs are now scheduled in the HISTORY.rc file, and the relevant
!  fields of Input_Opt will be populated in the RRTMG module routine
!  Init_RRTMG_Indices (called at startup).
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER                      :: I
    INTEGER                      :: N
    LOGICAL                      :: v_bool

    ! Strings
    CHARACTER(LEN=20)            :: str
    CHARACTER(LEN=255)           :: thisLoc
    CHARACTER(LEN=512)           :: errMsg
    CHARACTER(LEN=QFYAML_NamLen) :: key
    CHARACTER(LEN=QFYAML_NamLen) :: a_str(3)

    !========================================================================
    ! Config_RRTMG begins here!
    !========================================================================

    ! Initialize
    RC      = GC_SUCCESS
    errMsg  = ''
    thisLoc = ' -> at Read_Radiation_Menu (in module GeosCore/input_mod.F90)'

    !------------------------------------------------------------------------
    ! Turn on RRTMG?
    !------------------------------------------------------------------------
    key    = "operations%rrtmg%activate"
    v_bool = MISSING_BOOL
    CALL QFYAML_Add_Get( Config, key, v_bool, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Input_Opt%LRAD = v_bool

    !------------------------------------------------------------------------
    ! AOD wavelength selection? (You can have up to 3)
    !------------------------------------------------------------------------
    key   = "operations%rrtmg%aod_wavelengths_in_nm"
    a_str = MISSING_STR
    CALL QFYAML_Add_Get( Config, key, a_str, "", RC, dynamic_size=.TRUE. )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    ! Copy values into Input_Opt
    I = 0
    DO N = 1, SIZE( a_str )
       IF ( a_str(N) == MISSING_STR ) EXIT
       I = I + 1
       Input_Opt%nWvSelect      = I
       Input_Opt%StrWvSelect(I) = TRIM( ADJUSTL( a_str(N) ) )
       READ( a_str(N), *      )   Input_Opt%WvSelect(I)
       Input_Opt%WvSelect(I)    = RoundOff( Input_Opt%WvSelect(I), 2 )
    ENDDO

    !------------------------------------------------------------------------
    ! Turn on LW radiation calculation?
    !------------------------------------------------------------------------
    key    = "operations%rrtmg%longwave_fluxes"
    v_bool = MISSING_BOOL
    CALL QFYAML_Add_Get( Config, key, v_bool, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Input_Opt%LLWRAD = v_bool

    !------------------------------------------------------------------------
    ! Turn on SW radiation calculation?
    !------------------------------------------------------------------------
    key    = "operations%rrtmg%shortwave_fluxes"
    v_bool = MISSING_BOOL
    CALL QFYAML_Add_Get( Config, key, v_bool, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Input_Opt%LSWRAD = v_bool

    !------------------------------------------------------------------------
    ! Calculate for clear-sky?
    !------------------------------------------------------------------------
    key    = "operations%rrtmg%clear_sky_flux"
    v_bool = MISSING_BOOL
    CALL QFYAML_Add_Get( Config, key, v_bool, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Input_Opt%LSKYRAD(1) = v_bool

    !------------------------------------------------------------------------
    ! Calculate for all-sky?
    !------------------------------------------------------------------------
    key    = "operations%rrtmg%all_sky_flux"
    v_bool = MISSING_BOOL
    CALL QFYAML_Add_Get( Config, key, v_bool, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Input_Opt%LSKYRAD(2) = v_bool

    !=================================================================
    ! Error check settings
    !=================================================================

    ! Use of RRTMG necessitates recompilation
#if !defined( RRTMG )
    IF ( Input_Opt%LRAD ) THEN
       errMsg = 'LRAD=T but RRTMG is not defined at compile time!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
#endif

    ! Make sure radiation switches are turned off if RRTMG is off
    IF ( ( .not. Input_Opt%LRAD ) .and. Input_Opt%LLWRAD ) THEN
       errMsg = 'Cannot have LW fluxes turned on without RRTMG'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    IF ( ( .not. Input_Opt%LRAD ) .and. Input_Opt%LSWRAD ) THEN
       errMsg = 'Cannot have SW fluxes turned on without RRTMG'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    IF ( ( .not. Input_Opt%LRAD ) .and. Input_Opt%LSKYRAD(1) ) THEN
       errMsg = 'Cannot have clear-sky flux turned on without RRTMG'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    IF ( ( .not. Input_Opt%LRAD ) .and. Input_Opt%LSKYRAD(2) ) THEN
       errMsg = 'Cannot have all-sky flux turned on without RRTMG'
       CALL GC_Error( errMsg, RC, thisLoc )
    ENDIF

    !=================================================================
    ! Print to screen
    !=================================================================
    IF ( Input_Opt%amIRoot ) THEN
       WRITE( 6, '(/,a)' ) 'RRTMG SETTINGS'
       WRITE( 6, '(  a)' ) '--------------'
       DO N=1, Input_Opt%NWVSELECT
          WRITE( 6, 115     ) 'AOD output wavelength (nm)  : ', &
                               Input_Opt%WVSELECT(N)
       ENDDO
       WRITE( 6, 100 ) 'Turn on radiation?          : ', &
                        Input_Opt%LRAD
       WRITE( 6, 100 ) 'Consider LW                 : ', &
                        Input_Opt%LLWRAD
       WRITE( 6, 100 ) 'Consider SW                 : ', &
                        Input_Opt%LSWRAD
       WRITE( 6, 125 ) 'Clear-sky/All-sky           : ', &
                        Input_Opt%LSKYRAD(1), '/',       &
                        Input_Opt%LSKYRAD(2)
    ENDIF

    ! FORMAT statements
100 FORMAT( A, L5           )
110 FORMAT( A, I5           )
115 FORMAT( A, F7.1         )
120 FORMAT( A, 11I1         )
125 FORMAT( A, L5, A, L5    )
130 FORMAT( A, 12( A2, 1x ) )

  END SUBROUTINE Config_RRTMG
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Config_Photolysis
!
! !DESCRIPTION: Copies photolysis information from the Config object
!  to Input_Opt, and does necessary checks.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Config_Photolysis( Config, Input_Opt, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod, ONLY : OptInput
    USE RoundOff_Mod,  ONLY : RoundOff
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(QFYAML_t), INTENT(INOUT) :: Config      ! YAML Config object
    TYPE(OptInput), INTENT(INOUT) :: Input_Opt   ! Input Options object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    LOGICAL                      :: v_bool
    REAL(yp)                     :: v_real

    ! Strings
    CHARACTER(LEN=255)           :: thisLoc
    CHARACTER(LEN=512)           :: errMsg
    CHARACTER(LEN=QFYAML_NamLen) :: key
    CHARACTER(LEN=QFYAML_StrLen) :: v_str

    !=================================================================
    ! READ_PHOTOLYSIS_MENU begins here!
    !=================================================================

    ! Initialize
    RC      = GC_SUCCESS
    errMsg  = ''
    thisLoc = &
      ' -> at Config_Photolysis (in module GeosCore/input_mod.F90)'

    !-----------------------------------------------------------------
    ! Directory with photolysis input files
    !-----------------------------------------------------------------
    key   = "operations%photolysis%input_directory"
    v_str = MISSING_STR
    CALL QFYAML_Add_Get( Config, key, v_str, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Input_Opt%FAST_JX_DIR = TRIM( v_str )

    !------------------------------------------------------------------------
    ! Use online ozone in extinction calculations for FAST-JX?
    !------------------------------------------------------------------------
    key    = "operations%photolysis%overhead_O3%use_online_O3_from_model"
    v_bool = MISSING_BOOL
    CALL QFYAML_Add_Get( Config, key, v_bool, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Input_Opt%USE_ONLINE_O3 = v_bool

    !------------------------------------------------------------------------
    ! Use ozone columns from met fields?
    !------------------------------------------------------------------------
    key    = "operations%photolysis%overhead_O3%use_column_O3_from_met"
    v_bool = MISSING_BOOL
    CALL QFYAML_Add_Get( Config, key, v_bool, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Input_Opt%USE_O3_FROM_MET = v_bool

    !------------------------------------------------------------------------
    ! Use ozone columns from TOMS?
    !------------------------------------------------------------------------
    key    = "operations%photolysis%overhead_O3%use_TOMS_SBUV_O3"
    v_bool = MISSING_BOOL
    CALL QFYAML_Add_Get( Config, key, v_bool, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Input_Opt%USE_TOMS_O3 = v_bool

    !------------------------------------------------------------------------
    ! Photoylse nitrate aerosol?
    !------------------------------------------------------------------------
    key    = "operations%photolysis%photolyze_nitrate_aerosol%activate"
    v_bool = MISSING_BOOL
    CALL QFYAML_Add_Get( Config, key, v_bool, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Input_Opt%hvAerNIT = v_bool

    !------------------------------------------------------------------------
    ! Scalar for JHNO3 for photoylsing NITs aerosol
    !------------------------------------------------------------------------
    key    = &
     "operations%photolysis%photolyze_nitrate_aerosol%NITs_Jscale_JHNO3"
    v_str = MISSING_STR
    CALL QFYAML_Add_Get( Config, key, v_str, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    READ( v_str, * ) Input_Opt%hvAerNIT_JNITs
    Input_Opt%hvAerNIT_JNITs = Roundoff( Input_Opt%hvAerNIT_JNITs, 3 )

    !------------------------------------------------------------------------
    ! scalar for JHNO3 for photoylsing NIT aerosol (TMS, 23/08/18)
    !------------------------------------------------------------------------
    key    = &
     "operations%photolysis%photolyze_nitrate_aerosol%NIT_Jscale_JHNO2"
    v_str = MISSING_STR
    CALL QFYAML_Add_Get( Config, key, v_str, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    READ( v_str, * ) Input_Opt%hvAerNIT_JNIT
    Input_Opt%hvAerNIT_JNIT = Roundoff( Input_Opt%hvAerNIT_JNIT, 3 )

    !------------------------------------------------------------------------
    ! Fraction for JNITS/NIT channel A (HNO2) for NITs photoylsis
    !------------------------------------------------------------------------
    key   = &
     "operations%photolysis%photolyze_nitrate_aerosol%percent_channel_A_HONO"
    v_str = MISSING_STR
    CALL QFYAML_Add_Get( Config, key, v_str, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    READ( v_str, * ) Input_Opt%JNITChanA
    Input_Opt%JNITChanA = RoundOff( Input_Opt%JNITChanA, 3 )

    !------------------------------------------------------------------------
    ! Fraction for JNITs/NIT channel B (NO2) for NITs photoylsis
    !------------------------------------------------------------------------
    key    = &
     "operations%photolysis%photolyze_nitrate_aerosol%percent_channel_B_NO2"
    v_str = MISSING_STR
    CALL QFYAML_Add_Get( Config, key, v_str, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    READ( v_str, * ) Input_Opt%JNITChanB
    Input_Opt%JNITChanB= RoundOff( Input_Opt%JNITChanB, 3 )

    !========================================================================
    ! Error check settings
    !========================================================================

#ifndef MODEL_GEOS
    ! Cannot use Synoz with linearized mesospheric chemistry
    IF ( Input_Opt%ITS_A_FULLCHEM_SIM .and. Input_Opt%LINEAR_CHEM ) THEN
       IF (.not.Input_Opt%LLINOZ) THEN
          errMsg = 'Cannot use Synoz with linearized mesospheric chem.!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
    ENDIF
#endif

    ! FAST-JX is only used for fullchem and offline aerosol
    IF ( Input_Opt%ITS_A_FULLCHEM_SIM  .or. &
         Input_Opt%ITS_AN_AEROSOL_SIM  .or. &
         Input_Opt%ITS_A_MERCURY_SIM  ) THEN

       ! Make sure either O3 from met or TOMS is selected
       IF ( .not. Input_Opt%USE_O3_FROM_MET .and. &
            .not. Input_Opt%USE_TOMS_O3 ) THEN
          errMsg = 'Must select either O3 from met or TOMS/SBUV O3' &
                // 'for O3 values above the chemistry grid!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
       IF ( Input_Opt%USE_O3_FROM_MET .and. &
            Input_Opt%USE_TOMS_O3 ) THEN
          errMsg = 'Must select either O3 from met or TOMS/SBUV O3' &
                // 'for O3 values above the chemistry grid!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       ! Make sure specialty simulations select O3 from met or TOMS
       IF ( Input_Opt%ITS_AN_AEROSOL_SIM ) THEN
          IF ( Input_Opt%USE_ONLINE_O3 ) THEN
             errMsg= 'Cannot use online O3 for specialty simulations! ' &
                  // 'Select O3 from met or TOMS O3 instead.'
             CALL GC_Error( errMsg, RC, thisLoc )
             RETURN
          ENDIF
       ENDIF

    ELSE

       Input_Opt%USE_ONLINE_O3   = .FALSE.
       Input_Opt%USE_O3_FROM_MET = .FALSE.
       Input_Opt%USE_TOMS_O3     = .FALSE.

    ENDIF

    ! Turn off switches for simulations that don't use aerosols
    IF ( ( .not. Input_Opt%ITS_A_FULLCHEM_SIM )                        .and. &
         ( .not. Input_OPt%ITS_AN_AEROSOL_SIM ) ) THEN
       Input_Opt%hvAerNIT       = .FALSE.
       Input_Opt%hvAerNIT_JNITs = MISSING_REAL
       Input_Opt%hvAerNIT_JNIT  = MISSING_REAL
       Input_Opt%JNITChanA      = MISSING_REAL
       Input_Opt%JNITChanB      = MISSING_REAL
    ENDIF

    ! Return success
    RC = GC_SUCCESS

    !=================================================================
    ! Print to screen
    !=================================================================
    IF ( Input_Opt%amIRoot ) THEN
       WRITE( 6, '(/,a)' ) 'PHOTOLYSIS SETTINGS'
       WRITE( 6, '(  a)' ) '-------------------'
       WRITE( 6, 120     ) 'FAST-JX input directory     : ', &
                            TRIM( Input_Opt%FAST_JX_DIR )
       WRITE( 6, 100     ) 'Online ozone for FAST-JX?   : ', &
                            Input_Opt%USE_ONLINE_O3
       WRITE( 6, 100     ) 'Ozone from met for FAST-JX? : ', &
                            Input_Opt%USE_O3_FROM_MET
       WRITE( 6, 100     ) 'TOMS/SBUV ozone for FAST-JX?: ', &
                            Input_Opt%USE_TOMS_O3
       WRITE( 6, 100     ) 'Photolyse nitrate aerosol?  : ', &
                            Input_Opt%hvAerNIT
       WRITE( 6, 105     ) 'JNITs scaling of JHNO3      : ', &
                            Input_Opt%hvAerNIT_JNITs
       WRITE( 6, 105     ) 'JNIT scaling of JHNO3       : ', &
                            Input_Opt%hvAerNIT_JNIT
       WRITE( 6, 105     ) 'JNIT(s) channel A (HONO)    : ', &
                            Input_Opt%JNITChanA
       WRITE( 6, 105     ) 'JNIT(s) channel B (NO2)     : ', &
                            Input_Opt%JNITChanB
       ! Write more info
       IF ( Input_Opt%USE_ONLINE_O3 ) THEN
          WRITE( 6, '(a)' ) ''
          WRITE( 6, '(a)' ) 'NOTE ABOUT OVERHEAD O3 FOR FAST-JX:'
          WRITE( 6, '(a)' ) ' Online O3 from GEOS-Chem will be used'
          WRITE( 6, '(a)' ) ' to weight the O3 column within the'
          WRITE( 6, '(a)' ) ' chemistry grid and O3 from met or TOMS'
          WRITE( 6, '(a)' ) ' will be used outside the chemistry grid.'
       ENDIF
    ENDIF

    ! FORMAT statements
100 FORMAT( A, L5   )
105 FORMAT( A, F8.3 )
110 FORMAT( A, F4.2 )
120 FORMAT( A, A    )

    END SUBROUTINE Config_Photolysis
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: config_convection_mixing
!
! !DESCRIPTION: Copies convection & PBL mixing information from the Config
!  object to Input_Opt, and does necessary checks.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Config_Convection_Mixing( Config, Input_Opt, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod, ONLY : OptInput
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(QFYAML_t), INTENT(INOUT) :: Config      ! YAML Config object
    TYPE(OptInput), INTENT(INOUT) :: Input_Opt   ! Input options
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    LOGICAl                      :: v_bool

    ! Strings
    CHARACTER(LEN=255)           :: thisLoc
    CHARACTER(LEN=512)           :: errMsg
    CHARACTER(LEN=QFYAML_NamLen) :: key

    !========================================================================
    ! Config_Convection_Mixing begins here!
    !========================================================================

    ! Initialize
    RC      = GC_SUCCESS
    errMsg  = ""
    thisLoc = &
      ' -> at Config_Convection_Mixing (in module GeosCore/input_mod.F90)'

    !------------------------------------------------------------------------
    ! Turn on convection?
    !------------------------------------------------------------------------
    key    = "operations%convection%activate"
    v_bool = MISSING_BOOL
    CALL QFYAML_Add_Get( Config, key, v_bool, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Input_Opt%LCONV = v_bool

    !------------------------------------------------------------------------
    ! Turn on PBL mixing
    !------------------------------------------------------------------------
    key    = "operations%pbl_mixing%activate"
    v_bool = MISSING_BOOL
    CALL QFYAML_Add_Get( Config, key, v_bool, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Input_Opt%LTURB = v_bool

    !------------------------------------------------------------------------
    ! Use non-local PBL mixing?
    !------------------------------------------------------------------------
    key    = "operations%pbl_mixing%use_non_local_pbl"
    v_bool = MISSING_BOOL
    CALL QFYAML_Add_Get( Config, key, v_bool, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Input_Opt%LNLPBL = v_bool

    ! Set the PBL drydep flag. This determines if dry deposition is
    ! applied (and drydep frequencies are calculated) over the entire
    ! PBL or the first model layer only. For now, set this value
    ! automatically based upon the selected PBL scheme: 1st model layer
    ! for the non-local PBL scheme, full PBL for the full-mixing scheme.
    Input_Opt%PBL_DRYDEP = ( .not. Input_Opt%LNLPBL )

    ! Return success
    RC = GC_SUCCESS

    !=================================================================
    ! Print to screen
    !=================================================================
    IF ( Input_Opt%amIRoot ) THEN
       WRITE( 6, '(/,a)' ) 'CONVECTION SETTINGS'
       WRITE( 6, '(  a)' ) '-------------------'
       WRITE( 6, 100     ) 'Turn on cloud convection?   : ', &
                            Input_Opt%LCONV

       WRITE( 6, '(/,a)' ) 'PBL MIXING SETTINGS'
       WRITE( 6, '(  a)' ) '-------------------'
       WRITE( 6, 100     ) 'Turn on PBL mixing?         : ', &
                            Input_Opt%LTURB
       WRITE( 6, 100     ) 'Turn on non-local PBL?      : ', &
                            Input_Opt%LNLPBL
    ENDIF

    ! FORMAT statements
100 FORMAT( A, L5 )

  END SUBROUTINE Config_Convection_Mixing
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: config_drydep_wetdep
!
! !DESCRIPTION: Subroutine READ\_DEPOSITION\_MENU reads the DEPOSITION MENU
!  section of the GEOS-Chem input file.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Config_DryDep_WetDep( Config, Input_Opt, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod, ONLY : OptInput
    USE RoundOff_Mod,  ONLY : RoundOff
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(QFYAML_t), INTENT(INOUT) :: Config      ! YAML Config object
    TYPE(OptInput), INTENT(INOUT) :: Input_Opt   ! Input Options object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    LOGICAL                      :: v_bool
    INTEGER                      :: v_int

    ! Strings
    CHARACTER(LEN=255)           :: thisLoc
    CHARACTER(LEN=512)           :: errMsg
    CHARACTER(LEN=QFYAML_NamLen) :: key
    CHARACTER(LEN=QFYAML_StrLen) :: v_str

    !========================================================================
    ! Config_DryDep_WetDep begins here!
    !========================================================================

    ! Initialize
    RC      = GC_SUCCESS
    errMsg  = ''
    thisLoc = ' -> at Read_Deposition_Menu (in module GeosCore/input_mod.F90)'

    !------------------------------------------------------------------------
    ! Turn on drydep?
    !------------------------------------------------------------------------
    key    = "operations%dry_deposition%activate"
    v_bool = MISSING_BOOL
    CALL QFYAML_Add_Get( Config, key, v_bool, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Input_Opt%LDRYD = v_bool

    !------------------------------------------------------------------------
    ! Turn on CO2 effect on drydep?
    !------------------------------------------------------------------------
    key    = "operations%dry_deposition%CO2_effect%activate"
    v_bool = MISSING_BOOL
    CALL QFYAML_Add_Get( Config, key, v_bool, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Input_Opt%CO2_EFFECT = v_bool

    !------------------------------------------------------------------------
    ! CO2 level at simulation
    !------------------------------------------------------------------------
    key   = "operations%dry_deposition%CO2_effect%CO2_level"
    v_str = MISSING_STR
    CALL QFYAML_Add_Get( Config, key, v_str, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    READ( v_str, * ) Input_Opt%CO2_LEVEL
    Input_Opt%CO2_LEVEL = RoundOff( Input_Opt%CO2_LEVEL, 2 )

    !------------------------------------------------------------------------
    ! Reference CO2 level
    !------------------------------------------------------------------------
    key   = "operations%dry_deposition%CO2_effect%reference_CO2_level"
    v_str = MISSING_STR
    CALL QFYAML_Add_Get( Config, key, v_str, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    READ( v_str, * ) Input_Opt%CO2_REF
    Input_Opt%CO2_REF = RoundOff( Input_Opt%CO2_REF, 2 )

    !------------------------------------------------------------------------
    ! Diag for RA_alt above surface in meters
    !------------------------------------------------------------------------
    key   = "operations%dry_deposition%diag_alt_above_sfc_in_m"
    v_int = MISSING_INT
    CALL QFYAML_Add_Get( Config, key, v_int, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Input_Opt%RA_Alt_Above_Sfc = v_int

    !------------------------------------------------------------------------
    ! Turn on wetdep?
    !------------------------------------------------------------------------
    key    = "operations%wet_deposition%activate"
    v_bool = MISSING_BOOL
    CALL QFYAML_Add_Get( Config, key, v_bool, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Input_Opt%LWETD = v_bool

    !========================================================================
    ! Error check settings
    !========================================================================

    ! Turn off drydep for simulations that don't need it
    IF ( Input_Opt%ITS_A_TAGCO_SIM   ) Input_Opt%LDRYD = .FALSE.

    ! Turn off wetdep for simulations that don't need it
    IF ( Input_Opt%ITS_A_TAGO3_SIM   ) Input_Opt%LWETD = .FALSE.
    IF ( Input_Opt%ITS_A_TAGCO_SIM   ) Input_Opt%LWETD = .FALSE.
    IF ( Input_Opt%ITS_A_CH4_SIM     ) Input_Opt%LWETD = .FALSE.

    ! If CO2 effect on RS in turned on, calculate the scaling factor
    ! on Rs based on Franks et al. (2013) (ayhwong, 6/25/2019)
    If (Input_Opt%CO2_EFFECT) THEN
       Input_Opt%RS_SCALE = Input_Opt%CO2_LEVEL / Input_Opt%CO2_REF * &
                           (Input_Opt%CO2_LEVEL + 80.0_fp) *          &
                           (Input_Opt%CO2_REF   - 40.0_fp) /          &
                           (Input_Opt%CO2_LEVEL - 40.0_fp) /          &
                           (Input_Opt%CO2_REF   + 80.0_fp)
    ENDIF

    !========================================================================
    ! Print to screen
    !========================================================================
    IF ( Input_Opt%amIRoot ) THEN
       WRITE( 6, '(/,a)' ) 'DRY DEPOSITION SETTINGS'
       WRITE( 6, '(  a)' ) '-----------------------'
       WRITE( 6, 100     ) 'Turn on dry deposition?     : ', &
                            Input_Opt%LDRYD
       WRITE( 6, 100     ) 'Dry dep over full PBL?      : ', &
                            Input_Opt%PBL_DRYDEP
       WRITE( 6, 100     ) 'Turn on CO2 effect?         : ', &
                            Input_Opt%CO2_EFFECT
       WRITE( 6, 110     ) 'CO2 level                   : ', &
                            Input_Opt%CO2_LEVEL
       WRITE( 6, 110     ) 'CO2 reference level         : ', &
                            Input_Opt%CO2_REF
       WRITE( 6, 110     ) 'RIX scaling factor          : ', &
                            Input_Opt%RS_SCALE

       WRITE( 6, '(/,a)' ) 'WET DEPOSITION SETTINGS'
       WRITE( 6, '(  a)' ) '----------------------'
       WRITE( 6, 100     ) 'Turn on wet deposition?     : ', &
                            Input_Opt%LWETD
    ENDIF

    ! FORMAT statements
100 FORMAT( A, L5 )
110 FORMAT( A, f8.2 )

  END SUBROUTINE Config_DryDep_WetDep
!!EOC
#if !(defined( EXTERNAL_GRID ) || defined( EXTERNAL_FORCING ))
#ifdef BPCH_DIAG
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_gamap_menu
!
! !DESCRIPTION: Subroutine READ\_GAMAP\_MENU reads the GAMAP MENU section
!  of the GEOS-Chem input file.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Config_Gamap( Config, Input_Opt, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod, ONLY : OptInput
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(QFYAML_t), INTENT(INOUT) :: Config      ! YAML Config object
    TYPE(OptInput), INTENT(INOUT) :: Input_Opt   ! Input Options object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure
!
! !REVISION HISTORY:
!  25 Apr 2005 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Strings
    CHARACTER(LEN=255)           :: thisLoc
    CHARACTER(LEN=512)           :: errMsg
    CHARACTER(LEN=QFYAML_NamLen) :: key
    CHARACTER(LEN=QFYAML_StrLen) :: v_str

    !========================================================================
    ! Config_Gamap begins here!
    !========================================================================

    ! Initialize
    RC      = GC_SUCCESS
    errMsg  = ''
    thisLoc = ' -> at Read_Gamap_Menu (in module GeosCore/input_mod.F90)'

    !------------------------------------------------------------------------
    ! diaginfo.dat
    !------------------------------------------------------------------------
    key   = "extra_diagnostics%gamap%diaginfo_dat_file"
    v_str = MISSING_STR
    CALL QFYAML_Add_Get( Config, key, v_str, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Input_Opt%GAMAP_DIAGINFO = TRIM( ADJUSTL( v_str ) )

    !------------------------------------------------------------------------
    ! tracerinfo.dat
    !------------------------------------------------------------------------
    key   = "extra_diagnostics%gamap%tracerinfo_dat_file"
    v_str = MISSING_STR
    CALL QFYAML_Add_Get( Config, key, v_str, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Input_Opt%GAMAP_TRACERINFO = TRIM( ADJUSTL( v_str ) )

    !========================================================================
    ! Print to screen
    !========================================================================
    IF ( Input_Opt%amIRoot ) THEN
       WRITE( 6, '(/,a)' ) 'GAMAP SETTINGS (when -DBPCH_DIAG=y)'
       WRITE( 6, '(  a)' ) '-----------------------------------'
       WRITE( 6, '(a,a)' ) 'GAMAP "diaginfo.dat"   file : ',                 &
                            TRIM( Input_Opt%GAMAP_DIAGINFO   )
       WRITE( 6, '(a,a)' ) 'GAMAP "tracerinfo.dat" file : ',                 &
                            TRIM( Input_Opt%GAMAP_TRACERINFO )
    ENDIF

  END SUBROUTINE Config_Gamap
!EOC
!!------------------------------------------------------------------------------
!!                  GEOS-Chem Global Chemical Transport Model                  !
!!------------------------------------------------------------------------------
!!BOP
!!
!! !IROUTINE: read_output_menu
!!
!! !DESCRIPTION: Subroutine READ\_OUTPUT\_MENU reads the OUTPUT MENU section of
!!  the GEOS-Chem input file.
!!\\
!!\\
!! !INTERFACE:
!!
!  SUBROUTINE READ_OUTPUT_MENU( Input_Opt, RC )
!!
!! !USES:
!!
!    USE ErrCode_Mod
!    USE Input_Opt_Mod,  ONLY : OptInput
!!
!! !INPUT/OUTPUT PARAMETERS:
!!
!    TYPE(OptInput), INTENT(INOUT) :: Input_Opt   ! Input options
!!
!! !OUTPUT PARAMETERS:
!!
!    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure
!!
!! !REVISION HISTORY:
!!  20 Jul 2004 - R. Yantosca - Initial version
!!  See https://github.com/geoschem/geos-chem for complete history
!!EOP
!!------------------------------------------------------------------------------
!!BOC
!!
!! !LOCAL VARIABLES:
!!
!    ! Scalars
!    INTEGER            :: IOS
!
!    ! Strings
!    CHARACTER(LEN=255) :: errMsg, thisLoc
!
!    !=================================================================
!    ! READ_OUTPUT_MENU begins here!
!    !=================================================================
!
!    ! Initialize
!    RC      = GC_SUCCESS
!    errMsg  = 'Error reading the "input.geos" file!'
!    thisLoc = ' -> at Read_Output_Menu (in module GeosCore/input_mod.F90)'
!
!    ! Read info
!    READ( IU_GEOS, 100, IOSTAT=IOS ) Input_Opt%NJDAY
!100 FORMAT( 26x, 31i1, /  26x, 29i1, /, 26x, 31i1, /, 26x, 30i1, /, &
!            26x, 31i1, /, 26x, 30i1, /, 26x, 31i1, /, 26x, 31i1, /, &
!            26x, 30i1, /  26x, 31i1, /, 26x, 30i1, /, 26x, 31i1 )
!
!    ! Trap potential errors
!    IF ( IOS /= GC_SUCCESS ) THEN
!       CALL GC_Error( errMsg, RC, thisLoc )
!       RETURN
!    ENDIF
!
!    !=================================================================
!    ! Print to screen
!    !=================================================================
!    IF( Input_Opt%amIRoot ) THEN
!       WRITE( 6, '(/,a)' ) 'OUTPUT MENU'
!       WRITE( 6, '(  a)' ) '-----------'
!       WRITE( 6, 110     )
!       WRITE( 6, 120     )
!       WRITE( 6, 130     )
!       WRITE( 6, 140     ) Input_Opt%NJDAY
!    ENDIF
!
!    ! FORMAT statements
!110 FORMAT( '              1111111111222222222233' )
!120 FORMAT( '     1234567890123456789012345678901' )
!130 FORMAT( '     -------------------------------' )
!140 FORMAT( 'JAN--', 31i1, /, 'FEB--', 29i1, /, 'MAR--', 31i1, /, &
!            'APR--', 30i1, /, 'MAY--', 31i1, /, 'JUN--', 30i1, /, &
!            'JUL--', 31i1, /, 'AUG--', 31i1, /, 'SEP--', 30i1, /, &
!            'OCT--', 31i1, /, 'NOV--', 30i1, /, 'DEC--', 31i1 )
!
!    ! Make sure we have output at end of run
!    CALL IS_LAST_DAY_GOOD( Input_Opt, RC )
!
!    ! Trap potential errors
!    IF ( RC /= GC_SUCCESS ) THEN
!       CALL GC_Error( errMsg, RC, thisLoc )
!       RETURN
!    ENDIF
!
!  END SUBROUTINE READ_OUTPUT_MENU
!!EOC
!!------------------------------------------------------------------------------
!!                  GEOS-Chem Global Chemical Transport Model                  !
!!------------------------------------------------------------------------------
!!BOP
!!
!! !IROUTINE: read_diagnostic_menu
!!
!! !DESCRIPTION: Subroutine READ\_DIAGNOSTIC\_MENU reads the DIAGNOSTIC MENU
!!  section of the GEOS-Chem input file.
!!\\
!!\\
!! !INTERFACE:
!!
!  SUBROUTINE READ_DIAGNOSTIC_MENU( Input_Opt, RC )
!!
!! !USES:
!!
!    USE CMN_DIAG_MOD        ! Needed for timeseries diags (binary only)
!    USE CMN_SIZE_MOD,  ONLY : NDSTBIN
!    USE BPCH2_MOD,     ONLY : OPEN_BPCH2_FOR_WRITE
!    USE DIAG03_MOD,    ONLY : ND03,      PD03,      PD03_PL
!    USE DIAG53_MOD,    ONLY : ND53,      PD53
!    USE DRYDEP_MOD,    ONLY : NUMDEP
!    USE ErrCode_Mod
!    USE FILE_MOD,      ONLY : IU_BPCH
!    USE Input_Opt_Mod, ONLY : OptInput
!    USE TIME_MOD,      ONLY : GET_NYMDb, GET_NHMSb, EXPAND_DATE
!!
!! !INPUT/OUTPUT PARAMETERS:
!!
!    TYPE(OptInput), INTENT(INOUT) :: Input_Opt   ! Input options
!!
!! !OUTPUT PARAMETERS:
!!
!    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure
!!
!! !REVISION HISTORY:
!!  20 Jul 2004 - R. Yantosca - Initial version
!!  See https://github.com/geoschem/geos-chem for complete history
!!EOP
!!------------------------------------------------------------------------------
!!BOC
!!
!! !LOCAL VARIABLES:
!!
!    ! Scalars
!    INTEGER            :: M, N, N_TMP
!
!    ! Strings
!    CHARACTER(LEN=255) :: errMsg, thisLoc, BPCH_FILE
!
!    ! Arrays
!    CHARACTER(LEN=255) :: SUBSTRS(MAXDIM)
!
!    !=================================================================
!    ! READ_DIAGNOSTIC_MENU begins here!
!    !=================================================================
!
!    ! Initialize
!    RC      = GC_SUCCESS
!
!    errMsg  = 'Error reading the "input.geos" file!'
!    thisLoc = ' -> at Read_Diagnostic_Menu (in module GeosCore/input_mod.F90)'
!
!    ! Error check
!    IF ( CT1 /= 2 ) THEN
!       errMsg= 'SIMULATION MENU & ADVECTED SPECIES MENU ' // &
!               'must be read in first!'
!       CALL GC_Error( errMsg, RC, thisLoc )
!       RETURN
!    ENDIF
!
!    ! Binary punch file name
!    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1,  'BPCH_FILE', RC )
!    IF ( RC /= GC_SUCCESS ) THEN
!       CALL GC_Error( errMsg, RC, thisLoc )
!       RETURN
!    ENDIF
!    READ( SUBSTRS(1:N), '(a)' ) BPCH_FILE
!
!    ! Separator line
!    CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'separator', RC )
!    IF ( RC /= GC_SUCCESS ) THEN
!       CALL GC_Error( errMsg, RC, thisLoc )
!       RETURN
!    ENDIF
!
!    !--------------------------
!    ! ND03: Hg diagnostics
!    !--------------------------
!    CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'ND03', RC )
!    IF ( RC /= GC_SUCCESS ) THEN
!       CALL GC_Error( errMsg, RC, thisLoc )
!       RETURN
!    ENDIF
!    READ( SUBSTRS(1), * ) ND03
!    IF ( .not. Input_Opt%ITS_A_MERCURY_SIM ) ND03 = 0
!    CALL SET_TINDEX( Input_Opt, 03, ND03, SUBSTRS(2:N), N-1, PD03 )
!
!    !--------------------------
!    ! ND06: Dust emissions
!    !--------------------------
!    CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'ND06', RC )
!    IF ( RC /= GC_SUCCESS ) THEN
!       CALL GC_Error( errMsg, RC, thisLoc )
!       RETURN
!    ENDIF
!#ifdef TOMAS
!    READ( SUBSTRS(1), * ) ND06
!    IF ( .not. Input_Opt%LDUST ) ND06 = 0
!    CALL SET_TINDEX( Input_Opt, 06, ND06, SUBSTRS(2:N), N-1, NDSTBIN)
!#endif
!
!    !--------------------------
!    ! ND44 drydep vel & flux
!    !--------------------------
!
!    ! Number of species depends on simulation type
!    IF ( Input_Opt%ITS_A_TAGO3_SIM ) THEN
!       N_TMP = Input_Opt%N_ADVECT
!    ELSE
!       N_TMP = NUMDEP
!    ENDIF
!
!    CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'ND44', RC )
!    IF ( RC /= GC_SUCCESS ) THEN
!       CALL GC_Error( errMsg, RC, thisLoc )
!       RETURN
!    ENDIF
!#ifdef TOMAS
!    READ( SUBSTRS(1), * ) ND44
!    IF ( .not. Input_Opt%LDRYD ) ND44 = 0
!    CALL SET_TINDEX( Input_Opt, 44, ND44, SUBSTRS(2:N), N-1, N_TMP )
!#endif
!
!    !--------------------------
!    ! ND53: POPS
!    !--------------------------
!    CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'ND53', RC )
!    IF ( RC /= GC_SUCCESS ) THEN
!       CALL GC_Error( errMsg, RC, thisLoc )
!       RETURN
!    ENDIF
!    READ( SUBSTRS(1), * ) ND53
!    IF ( .not. Input_Opt%ITS_A_POPS_SIM ) ND53 = 0
!    CALL SET_TINDEX( Input_Opt, 53, ND53, SUBSTRS(2:N), N-1, PD53 )
!
!    !--------------------------
!    ! ND59: TOMAS aerosol emiss
!    !--------------------------
!    CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'ND59', RC )
!    IF ( RC /= GC_SUCCESS ) THEN
!       CALL GC_Error( errMsg, RC, thisLoc )
!       RETURN
!    ENDIF
!#ifdef TOMAS
!    READ( SUBSTRS(1), * ) ND59
!    CALL SET_TINDEX( Input_Opt, 59, ND59, SUBSTRS(2:N), N-1, PD59 )
!#endif
!
!    !--------------------------
!    ! ND60: TOMAS rate
!    !--------------------------
!    CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'ND60', RC )
!    IF ( RC /= GC_SUCCESS ) THEN
!       CALL GC_Error( errMsg, RC, thisLoc )
!       RETURN
!    ENDIF
!#ifdef TOMAS
!    READ( SUBSTRS(1), * ) ND60
!    CALL SET_TINDEX( Input_Opt, 60, ND60, SUBSTRS(2:N), N-1, PD60 )
!#endif
!
!    !--------------------------
!    ! ND61: 3-D TOMAS rate
!    !--------------------------
!    CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'ND61', RC )
!    IF ( RC /= GC_SUCCESS ) THEN
!       CALL GC_Error( errMsg, RC, thisLoc )
!       RETURN
!    ENDIF
!#ifdef TOMAS
!    READ( SUBSTRS(1), * ) ND61
!    CALL SET_TINDEX( Input_Opt, 61, ND61, SUBSTRS(2:N), N-1, PD61 )
!#endif
!
!    !=================================================================
!    ! %%%%% IF BPCH DIAGNOSTICS ARE ACTIVATED (BPCH_DIAG=y) %%%%%
!    !
!    ! Copy shadow variables to diagnostic variables and
!    ! call various bpch diagnostic setup routines
!    !=================================================================
!    Input_Opt%ND53  = ND53
!
!#ifdef TOMAS
!    Input_Opt%ND06  = ND06
!    Input_Opt%ND44  = ND44
!    Input_Opt%ND59  = ND59
!    Input_Opt%ND60  = ND60
!    Input_Opt%ND61  = ND61
!#endif
!
!    ! Loop over # of diagnostics
!    DO M = 1, Input_Opt%Max_BPCH_Diag
!       Input_Opt%TCOUNT(M)       = TCOUNT(M)
!       Input_Opt%TMAX(M)         = TMAX(M)
!
!       ! Loop over tracers per diagnostic
!       DO N = 1, Input_Opt%N_ADVECT
!          Input_Opt%TINDEX(M,N)  = TINDEX(M,N)
!       ENDDO
!    ENDDO
!
!    !=================================================================
!    ! Call other setup routines
!    !================================================================
!
!    ! Expand YYYYMMDD tokens in the bpch file name
!    CALL EXPAND_DATE( BPCH_FILE, GET_NYMDb(), GET_NHMSb() )
!
!    ! Find a free file LUN
!    IU_BPCH = findFreeLUN()
!
!    ! Open the binary punch file for output
!    CALL OPEN_BPCH2_FOR_WRITE( IU_BPCH, BPCH_FILE )
!
!    ! Return success
!    RC = GC_SUCCESS
!
!  END SUBROUTINE READ_DIAGNOSTIC_MENU
!!EOC
!!------------------------------------------------------------------------------
!!                  GEOS-Chem Global Chemical Transport Model                  !
!!------------------------------------------------------------------------------
!!BOP
!!
!! !IROUTINE: set_tindex
!!
!! !DESCRIPTION: Subroutine SET\_TINDEX sets the TINDEX and TMAX arrays,
!!  which determine how many tracers to print to the punch file.
!!\\
!!\\
!! !INTERFACE:
!!
!  SUBROUTINE SET_TINDEX(Input_Opt, N_DIAG, L_DIAG, SUBSTRS, N, NMAX)
!!
!! !USES:
!!
!#ifdef TOMAS
!    USE CHARPAK_MOD,   ONLY : TXTEXT   ! (win, 7/14/09)
!#endif
!    USE CMN_DIAG_MOD        ! TMAX, TINDEX
!    USE Input_Opt_Mod, ONLY : OptInput
!!
!! !INPUT PARAMETERS:
!!
!    TYPE(OptInput),     INTENT(IN) :: Input_Opt   ! Input Options object
!    INTEGER,            INTENT(IN) :: N_DIAG      ! GEOS-Chem diagnostic #
!    INTEGER,            INTENT(IN) :: N           ! # of valid substrs passed
!    INTEGER,            INTENT(IN) :: NMAX        ! Max # of tracers allowed
!    INTEGER,            INTENT(IN) :: L_DIAG      ! # of levels to save
!    CHARACTER(LEN=255), INTENT(IN) :: SUBSTRS(N)  ! Substrs passed from
!                                                  !  READ_DIAGNOSTIC_MENU
!!
!! !REMARKS:
!!  NOTE: This routine converts to a stub when BPCH_DIAG=n (bmy, 1/16/18)
!!
!! !REVISION HISTORY:
!!  20 Jul 2004 - R. Yantosca - Initial version
!!  See https://github.com/geoschem/geos-chem for complete history
!!EOP
!!------------------------------------------------------------------------------
!!BOC
!!
!! !LOCAL VARIABLES:
!!
!    LOGICAL, SAVE       :: FIRST = .TRUE.
!    LOGICAL             :: IS_ALL
!    INTEGER             :: M
!#ifdef TOMAS
!    INTEGER             :: NN,     COL,     IFLAG, TC     ! (win, 7/14/09)
!    CHARACTER (LEN=255) :: WORD,   SUBWORD, TMP1,  TMP2   ! (win, 7/14/09)
!    INTEGER             :: MINTMP, MAXTMP                 ! (win, 7/14/09)
!#endif
!
!    !=================================================================
!    ! SET_TINDEX begins here!
!    !=================================================================
!
!    ! Error check
!    IF ( N < 1 ) THEN
!       IF ( Input_Opt%amIRoot ) THEN
!          WRITE( 6, '(a)' ) 'ERROR: N must be 1 or greater!'
!          WRITE( 6, '(a)' ) 'STOP in SET_TINDEX (input_mod.F90)'
!          WRITE( 6, '(a)' ) REPEAT( '=', 79 )
!       ENDIF
!       STOP
!    ENDIF
!
!    !=================================================================
!    ! If the word "all" is present, then set TMAX, TINDEX to all
!    ! available tracers for the given diagnostic.  Otherwise, just
!    ! use the tracers that were read in from the line
!    !=================================================================
!    IF ( TRIM( SUBSTRS(1) ) == 'all'  .or. &
!         TRIM( SUBSTRS(1) ) == 'ALL' ) THEN
!
!       ! TMAX is the max # of tracers to print out
!       TMAX(N_DIAG) = NMAX
!
!       ! Fill TINDEX with all possible diagnostic tracer numbers
!       DO M = 1, TMAX(N_DIAG)
!          TINDEX(N_DIAG,M) = M
!       ENDDO
!
!       ! Set flag
!       IS_ALL = .TRUE.
!
!    ELSE
!
!#ifdef TOMAS
!!(win, 7/14/09)  use TXTEXT and split the read in characters by -
!
!       COL      = 1
!       NN       = 0
!       SUBWORD  = ''
!       IFLAG    = 0
!
!       ! Use explicit DO-loop
!       DO M = 1, N
!          WORD = SUBSTRS(M)
!
!          ! Check if the characters are a range with - in the middle
!          CALL TXTEXT ( '-', WORD, COL, SUBWORD, IFLAG )
!
!          ! Found a dash!  Get the numbers on both sides of the dash
!          ! since these the min and max of the tracer range
!          IF ( IFLAG == 0 ) THEN
!             TMP1 = TRIM( WORD(      1:COL-1      ) )
!             TMP2 = TRIM( WORD( COL+1:LEN_TRIM( WORD ) ) )
!
!             READ( TMP1, * ) MINTMP
!             READ( TMP2, * ) MAXTMP
!
!             DO TC = MINTMP, MAXTMP
!                NN = NN + 1
!                TINDEX( N_DIAG, NN ) = TC
!             ENDDO
!
!          ! If we haven't found a dash, then there is only one number,
!          ! so that number is both the min and max of the tracer range
!          ELSE IF ( IFLAG == -1 ) THEN
!             NN = NN + 1
!             TMP1 = TRIM( WORD )
!             READ( TMP1, * ) TINDEX( N_DIAG, NN )
!          ENDIF
!
!       ENDDO
!
!       ! Set TMAX to the counted # of tracers
!       TMAX( N_DIAG ) = NN
!#else
!       ! Otherwise, set TMAX, TINDEX to the # of tracers
!       ! listed in "input.ctm" -- need some error checks too
!       TMAX(N_DIAG) = N
!
!       ! Use explicit DO-loop
!       DO M = 1, N
!          READ( SUBSTRS(M:M), * ) TINDEX(N_DIAG,M)
!       ENDDO
!#endif
!       ! Set flag
!       IS_ALL = .FALSE.
!
!    ENDIF
!
!    !=================================================================
!    ! Print to screen
!    !=================================================================
!
!    ! First-time printing only
!    IF ( FIRST ) THEN
!       IF( Input_Opt%amIRoot ) THEN
!          WRITE( 6,'(/,a)' ) 'DIAGNOSTIC MENU'
!          WRITE( 6,'(  a)' ) '---------------'
!          WRITE( 6,'(  a)' ) 'Diag    L   Tracers being saved to disk'
!       ENDIF
!       FIRST = .FALSE.
!    ENDIF
!
!    ! Test if all tracers are being printed out
!    IF ( IS_ALL ) THEN
!
!       ! Print abbreviated output string
!       IF ( L_DIAG > 0 ) THEN
!          IF ( Input_Opt%amIRoot ) THEN
!             WRITE( 6, 100 ) N_DIAG, L_DIAG, 1, TMAX(N_DIAG)
!          ENDIF
!100       FORMAT( 'ND', i2.2, 2x, i3, 1x, i3, ' -', i3 )
!       ENDIF
!
!    ELSE
!
!       ! Or just list each tracer
!       ! Print each diagnostic and # of tracers that will print out
!       IF ( L_DIAG > 0 ) THEN
!          IF ( Input_Opt%amIRoot ) THEN
!             WRITE( 6, 110 ) N_DIAG, L_DIAG, &
!                             ( TINDEX(N_DIAG,M), M=1,TMAX(N_DIAG) )
!          ENDIF
!110       FORMAT( 'ND', i2, 2x, i3, 1x, 100i3 )
!       ENDIF
!
!    ENDIF
!
!  END SUBROUTINE SET_TINDEX
!!EOC
#endif
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: config_planeflight
!
! !DESCRIPTION: Copies grid information from the Config object
!  to Input_Opt, and does necessary checks.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Config_PlaneFlight( Config, Input_Opt, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod,   ONLY : OptInput
    USE PlaneFlight_Mod, ONLY : SET_PLANEFLIGHT
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(QFYAML_t), INTENT(INOUT) :: Config      ! YAML Config object
    TYPE(OptInput), INTENT(INOUT) :: Input_Opt   ! Input Options Object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    LOGICAL                      :: v_bool

    ! Strings
    CHARACTER(LEN=255)           :: key
    CHARACTER(LEN=255)           :: thisLoc
    CHARACTER(LEN=255)           :: errMsg
    CHARACTER(LEN=QFYAML_StrLen) :: v_str

    !========================================================================
    ! READ_PLANEFLIGHT_MENU begins here!
    !========================================================================

    ! Initialize
    RC      = GC_SUCCESS
    errMsg  = ''
    thisLoc = ' -> at Config_PlaneFlight (in module GeosCore/input_mod.F90)'

    !------------------------------------------------------------------------
    ! Turn on planeflight diagnostic?
    !------------------------------------------------------------------------
    key    = "extra_diagnostics%planeflight%activate"
    v_bool = MISSING_BOOL
    CALL QFYAML_Add_Get( Config, key, v_bool, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Input_Opt%Do_Planeflight = v_bool

    !------------------------------------------------------------------------
    ! Input file name (w/ flight track data points)
    !------------------------------------------------------------------------
    key   = "extra_diagnostics%planeflight%flight_track_file"
    v_str = MISSING_STR
    CALL QFYAML_Add_Get( Config, key, v_str, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Input_Opt%Planeflight_InFile = TRIM( v_str )

    !------------------------------------------------------------------------
    ! Output file name
    !------------------------------------------------------------------------
    key   = "extra_diagnostics%planeflight%output_file"
    v_str = MISSING_STR
    CALL QFYAML_Add_Get( Config, key, v_str, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Input_Opt%Planeflight_OutFile = TRIM( v_str )

    !========================================================================
    ! Print to screen
    !========================================================================
    IF( Input_Opt%amIRoot ) THEN
       WRITE( 6, '(/,a)' ) 'PLANEFLIGHT DIAGNOSTIC SETTINGS'
       WRITE( 6, '(  a)' ) '-------------------------------'
       WRITE( 6, 100 ) 'Turn on planeflight diag?   : ', &
                        Input_Opt%Do_Planeflight
       WRITE( 6, 110 ) 'Flight track input file     : ', &
                        TRIM( Input_Opt%Planeflight_InFile )
       WRITE( 6, 110 ) 'Output file name            : ', &
                        TRIM( Input_Opt%Planeflight_OutFile )
    ENDIF

    ! FORMAT statements
100 FORMAT( A, L5    )
110 FORMAT( A, A     )

    !========================================================================
    ! Call setup routines from other F90 modules
    !========================================================================

    ! Pass variables to "planeflight_mod.F90"
    CALL Set_PlaneFlight( PF       = Input_Opt%Do_Planeflight,               &
                          In_File  = Input_Opt%Planeflight_InFile,           &
                          Out_File = Input_Opt%Planeflight_OutFile          )

  END SUBROUTINE Config_PlaneFlight
#endif
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_obspack_menu
!
! !DESCRIPTION: Subroutine READ\_OBSPACK\_MENU reads the OBSPACK MENU
!  section of the GEOS-Chem input file.  This turns on the ObsPack diagnostic.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Config_ObsPack( Config, Input_Opt, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod, ONLY : OptInput
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(QFYAML_t), INTENT(INOUT) :: Config      ! YAML Config object
    TYPE(OptInput), INTENT(INOUT) :: Input_Opt   ! Input Options object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER                      :: N
    LOGICAL                      :: v_bool

    ! Strings
    CHARACTER(LEN=255)           :: thisLoc
    CHARACTER(LEN=255)           :: errMsg
    CHARACTER(LEN=QFYAML_NamLen) :: key
    CHARACTER(LEN=QFYAML_StrLen) :: v_str

    ! String arrays
    CHARACTER(LEN=QFYAML_StrLen) :: a_str(QFYAML_MaxArr)

    !========================================================================
    ! READ_OBSPACK_MENU begins here!
    !========================================================================

    ! Initialize
    RC      = GC_SUCCESS
    errMsg  = 'Error reading the "input.geos" file!'
    thisLoc = ' -> at Read_ObsPack_Menu (in module GeosCore/input_mod.F90)'

    !------------------------------------------------------------------------
    ! Turn on ObsPack diagnostic?
    !------------------------------------------------------------------------
    key    = "extra_diagnostics%obspack%activate"
    v_bool = MISSING_BOOL
    CALL QFYAML_Add_Get( Config, key, v_bool, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Input_Opt%Do_ObsPack = v_bool

    !------------------------------------------------------------------------
    ! ObsPack quiet output?
    !------------------------------------------------------------------------
    key    = "extra_diagnostics%obspack%quiet_logfile_output"
    v_bool = MISSING_BOOL
    CALL QFYAML_Add_Get( Config, key, v_bool, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Input_Opt%ObsPack_Quiet = v_bool

    !------------------------------------------------------------------------
    ! Input file name (w/ coordinates and sampling strategy)
    !------------------------------------------------------------------------
    key   = "extra_diagnostics%obspack%input_file"
    v_str = MISSING_STR
    CALL QFYAML_Add_Get( Config, key, v_str, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Input_Opt%ObsPack_InputFile = TRIM( v_str )

    !------------------------------------------------------------------------
    ! Output file name
    !------------------------------------------------------------------------
    key   = "extra_diagnostics%obspack%output_file"
    v_str = MISSING_STR
    CALL QFYAML_Add_Get( Config, key, v_str, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Input_Opt%ObsPack_OutputFile = TRIM( v_str )

    !------------------------------------------------------------------------
    ! Species names
    !------------------------------------------------------------------------
    key   = "extra_diagnostics%obspack%output_species"
    a_str = MISSING_STR
    CALL QFYAML_Add_Get( Config, key, v_str, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! Copy species names into Input_Opt
    !------------------------------------------------------------------------
    Input_Opt%ObsPack_nSpc = 0
    DO N = 1, SIZE( a_str )

       ! Stop iterationg when we find a missing value
       IF ( TRIM( a_str(N) ) == MISSING_STR ) EXIT

       ! If wildcard for all species is requested then update the
       ! list of species to track to be the list of advected species
       ! and exit from further
       IF ( N==1 .AND. INDEX( a_str(1) , '?ALL' ) >  0)  THEN
          Input_Opt%ObsPack_SpcName = Input_Opt%AdvectSpc_Name
          Input_Opt%ObsPack_nSpc   = Input_Opt%N_Advect
          EXIT
       ENDIF

       ! Otherwise, increment the count and copy the Obspack species
       ! name into the Input_Opt object
       Input_Opt%ObsPack_nSpc = Input_Opt%ObsPack_nSpc + 1
       Input_Opt%ObsPack_SpcName(Input_Opt%ObsPack_nSpc) = TRIM( a_str(N) )
    ENDDO

    !=================================================================
    ! Print to screen
    !=================================================================
    IF( Input_Opt%amIRoot ) THEN
       WRITE( 6, '(/,a)' ) 'OBSPACK_MENU'
       WRITE( 6, '(  a)' ) '----------------'
       WRITE( 6, 100 ) 'Turn on ObsPack diagnostic? : ', &
                        Input_Opt%Do_ObsPack
       WRITE( 6, 100 ) 'Suppress logfile output?    : ', &
                        Input_Opt%ObsPack_Quiet
       WRITE( 6, 110 ) 'ObsPack input file          : ', &
                        TRIM( Input_Opt%ObsPack_InputFile  )
       WRITE( 6, 110 ) 'ObsPack output file         : ', &
                        TRIM( Input_Opt%ObsPack_OutputFile )
    ENDIF

    ! FORMAT statements
100 FORMAT( A, L5    )
110 FORMAT( A, A     )

  END SUBROUTINE Config_ObsPack
!EOC
!#if !(defined( EXTERNAL_GRID ) || defined( EXTERNAL_FORCING ))
!#ifdef BPCH_DIAG
!!EOC
!!------------------------------------------------------------------------------
!!                  GEOS-Chem Global Chemical Transport Model                  !
!!------------------------------------------------------------------------------
!!BOP
!!
!! !IROUTINE: read_nd51_menu
!!
!! !DESCRIPTION: Subroutine READ\_ND51\_MENU reads the ND51 MENU section of
!!  the GEOS-Chem input file.
!!\\
!!\\
!! !INTERFACE:
!!
!  SUBROUTINE READ_ND51_MENU( Input_Opt, RC )
!!
!! !USES:
!!
!    USE ErrCode_Mod
!    USE Input_Opt_Mod,      ONLY : OptInput
!!
!! !INPUT/OUTPUT PARAMETERS:
!!
!    TYPE(OptInput), INTENT(INOUT) :: Input_Opt   ! Input options
!!
!! !OUTPUT PARAMETERS:
!!
!    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure
!!
!! !REVISION HISTORY:
!!  20 Jul 2004 - R. Yantosca - Initial version
!!  See https://github.com/geoschem/geos-chem for complete history
!!EOP
!!------------------------------------------------------------------------------
!!BOC
!!
!! !LOCAL VARIABLES:
!!
!    ! Scalars
!    INTEGER            :: N
!
!    ! Strings
!    CHARACTER(LEN=255) :: errMsg, thisLoc
!
!    ! Arrays
!    CHARACTER(LEN=255) :: SUBSTRS(MAXDIM)
!
!    !=================================================================
!    ! READ_ND51_MENU begins here!
!    !=================================================================
!
!    ! Initialize
!    RC      = GC_SUCCESS
!
!    errMsg  = 'Error reading the "input.geos" file!'
!    thisLoc = ' -> at Read_ND51_Menu (in module GeosCore/input_mod.F90)'
!
!    ! Error check
!    IF ( CT1 /= 2 ) THEN
!       errMsg = 'SIMULATION MENU & ADVECTED SPECIES MENU ' // &
!                'must be read in first!'
!       CALL GC_Error( errMsg, RC, thisLoc )
!       RETURN
!    ENDIF
!
!    ! Turn on ND51 diagnostic
!    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1,  'read_nd51_menu:1', RC )
!    IF ( RC /= GC_SUCCESS ) THEN
!       CALL GC_Error( errMsg, RC, thisLoc )
!       RETURN
!    ENDIF
!    READ( SUBSTRS(1:N), * ) Input_Opt%DO_ND51
!
!    ! Instantaneous 3-D timeseries file
!    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1,  'read_nd51_menu:2', RC )
!    IF ( RC /= GC_SUCCESS ) THEN
!       CALL GC_Error( errMsg, RC, thisLoc )
!       RETURN
!    ENDIF
!    READ( SUBSTRS(1:N), '(a)' ) Input_Opt%ND51_FILE
!
!    ! Output as hdf
!    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1,  'read_nd51_menu:3', RC )
!    IF ( RC /= GC_SUCCESS ) THEN
!       CALL GC_Error( errMsg, RC, thisLoc )
!       RETURN
!    ENDIF
!    READ( SUBSTRS(1:N), * ) Input_Opt%LND51_HDF
!
!    ! (bmy, 12/21/09)
!    Input_Opt%LND51_HDF = .FALSE.
!
!    ! Tracers to include
!    CALL SPLIT_ONE_LINE( SUBSTRS, Input_Opt%N_ND51, -1, &
!                         'read_nd51_menu:4', RC )
!    IF ( RC /= GC_SUCCESS ) THEN
!       CALL GC_Error( errMsg, RC, thisLoc )
!       RETURN
!    ENDIF
!    DO N = 1, Input_Opt%N_ND51
!       READ( SUBSTRS(N), * ) Input_Opt%ND51_TRACERS(N)
!    ENDDO
!
!    ! NHMS_WRITE
!    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1,  'read_nd51_menu:6', RC )
!    IF ( RC /= GC_SUCCESS ) THEN
!       CALL GC_Error( errMsg, RC, thisLoc )
!       RETURN
!    ENDIF
!    READ( SUBSTRS(1:N), * ) Input_Opt%ND51_HR_WRITE
!
!    ! Make sure ND51_HR_WRITE is in the range 0-23.999 hrs
!    Input_Opt%ND51_HR_WRITE = MOD( Input_Opt%ND51_HR_WRITE, 24e+0_fp )
!
!    ! HR1, HR2
!    CALL SPLIT_ONE_LINE( SUBSTRS, N, 2,  'read_nd51_menu:7', RC )
!    IF ( RC /= GC_SUCCESS ) THEN
!       CALL GC_Error( errMsg, RC, thisLoc )
!       RETURN
!    ENDIF
!    READ( SUBSTRS(1:N), * ) Input_Opt%ND51_HR1, Input_Opt%ND51_HR2
!
!    ! IMIN, IMAX
!    CALL SPLIT_ONE_LINE( SUBSTRS, N, 2,  'read_nd51_menu:8', RC )
!    IF ( RC /= GC_SUCCESS ) THEN
!       CALL GC_Error( errMsg, RC, thisLoc )
!       RETURN
!    ENDIF
!    READ( SUBSTRS(1:N), * ) Input_Opt%ND51_IMIN, Input_Opt%ND51_IMAX
!
!    ! JMIN, JMAX
!    CALL SPLIT_ONE_LINE( SUBSTRS, N, 2,  'read_nd51_menu:9', RC )
!    IF ( RC /= GC_SUCCESS ) THEN
!       CALL GC_Error( errMsg, RC, thisLoc )
!       RETURN
!    ENDIF
!    READ( SUBSTRS(1:N), * ) Input_Opt%ND51_JMIN, Input_Opt%ND51_JMAX
!
!    ! LMIN, LMAX
!    CALL SPLIT_ONE_LINE( SUBSTRS, N, 2,  'read_nd51_menu:10', RC )
!    IF ( RC /= GC_SUCCESS ) THEN
!       CALL GC_Error( errMsg, RC, thisLoc )
!       RETURN
!    ENDIF
!    READ( SUBSTRS(1:N), * ) Input_Opt%ND51_LMIN, Input_Opt%ND51_LMAX
!
!    ! Separator line
!    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1,  'read_nd51_menu:11', RC )
!    IF ( RC /= GC_SUCCESS ) THEN
!       CALL GC_Error( errMsg, RC, thisLoc )
!       RETURN
!    ENDIF
!
!    !=================================================================
!    ! Print to screen
!    !=================================================================
!    IF ( Input_Opt%amIRoot ) THEN
!       WRITE( 6, '(/,a)' ) 'ND51 MORNING OR AFTERNOON TIMESERIES MENU'
!       WRITE( 6, '(  a)' ) '-----------------------------------------'
!       WRITE( 6, 100 ) 'Turn on ND51 timeseries?    : ', &
!                        Input_Opt%DO_ND51
!       WRITE( 6, 110 ) 'ND51 timeseries file name   : ', &
!                        TRIM( Input_Opt%ND51_FILE )
!       WRITE( 6, 100 ) 'Output as HDF?              : ', &
!                        Input_Opt%LND51_HDF
!       WRITE( 6, 120 ) 'ND51 timeseries tracers     : ',  &
!                        ( Input_Opt%ND51_TRACERS(N), N=1, &
!                          Input_Opt%N_ND51 )
!       WRITE( 6, 140 ) 'ND51 hour to write to disk  : ', &
!                        Input_Opt%ND51_HR_WRITE
!       WRITE( 6, 140 ) 'ND51 averaging period [GMT] : ', &
!                        Input_Opt%ND51_HR1,  Input_Opt%ND51_HR2
!       WRITE( 6, 130 ) 'ND51 longitude limits       : ', &
!                        Input_Opt%ND51_IMIN, Input_Opt%ND51_IMAX
!       WRITE( 6, 130 ) 'ND51 latitude  limits       : ', &
!                        Input_Opt%ND51_JMIN, Input_Opt%ND51_JMAX
!       WRITE( 6, 130 ) 'ND51 altitude  limits       : ', &
!                        Input_Opt%ND51_LMIN, Input_Opt%ND51_LMAX
!    ENDIF
!
!    ! FORMAT statements
!100 FORMAT( A, L5    )
!110 FORMAT( A, A     )
!120 FORMAT( A, 100I4 )
!130 FORMAT( A, 2I5   )
!140 FORMAT( A, 2F5.1 )
!
!  END SUBROUTINE READ_ND51_MENU
!!EOC
!!------------------------------------------------------------------------------
!!                  GEOS-Chem Global Chemical Transport Model                  !
!!------------------------------------------------------------------------------
!!BOP
!!
!! !IROUTINE: read_nd51b_menu
!!
!! !DESCRIPTION: Subroutine READ\_ND51b\_MENU reads the ND51 MENU section
!!  of the GEOS-Chem input file.
!!\\
!!\\
!! !INTERFACE:
!!
!  SUBROUTINE READ_ND51b_MENU( Input_Opt, RC )
!!
!! !USES:
!!
!    USE ErrCode_Mod
!    USE Input_Opt_Mod, ONLY : OptInput
!!
!! !INPUT/OUTPUT PARAMETERS:
!!
!    TYPE(OptInput), INTENT(INOUT) :: Input_Opt   ! Input options
!!
!! !OUTPUT PARAMETERS:
!!
!    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure
!!
!! !REVISION HISTORY:
!!  21 Dec 2009 - Aaron van D - Initial version
!!  See https://github.com/geoschem/geos-chem for complete history
!!EOP
!!------------------------------------------------------------------------------
!!BOC
!!
!! !LOCAL VARIABLES:
!!
!    ! Scalars
!    INTEGER              :: N
!
!    ! Strings
!    CHARACTER(LEN=255)   :: errMsg, thisLoc
!
!    ! Arrays
!    CHARACTER(LEN=255)   :: SUBSTRS(MAXDIM)
!
!    !=================================================================
!    ! READ_ND51b_MENU begins here!
!    !=================================================================
!
!    ! Initialize
!    RC      = GC_SUCCESS
!
!    errMsg  = 'Error reading the "input.geos" file!'
!    thisLoc = ' -> at Read_ND51b_MENU  (in module GeosCore/input_mod.F90)'
!
!    ! Error check
!    IF ( CT1 /= 2 ) THEN
!       errMsg = 'SIMULATION MENU & ADVECTED SPECIES MENU ' // &
!                'must be read in first!'
!       CALL GC_Error( errMsg, RC, thisLoc )
!       RETURN
!    ENDIF
!
!    ! Turn on ND51b diagnostic
!    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1,  'read_nd51b_menu:1', RC )
!    IF ( RC /= GC_SUCCESS ) THEN
!       CALL GC_Error( errMsg, RC, thisLoc )
!       RETURN
!    ENDIF
!    READ( SUBSTRS(1:N), * ) Input_Opt%DO_ND51b
!
!    ! Instantaneous 3-D timeseries file
!    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1,  'read_nd51b_menu:2', RC )
!    IF ( RC /= GC_SUCCESS ) THEN
!       CALL GC_Error( errMsg, RC, thisLoc )
!       RETURN
!    ENDIF
!    READ( SUBSTRS(1:N), '(a)' ) Input_Opt%ND51b_FILE
!
!    ! Output as hdf
!    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1,  'read_nd51b_menu:3', RC )
!    IF ( RC /= GC_SUCCESS ) THEN
!       CALL GC_Error( errMsg, RC, thisLoc )
!       RETURN
!    ENDIF
!    READ( SUBSTRS(1:N), * ) Input_Opt%LND51b_HDF
!
!    Input_Opt%LND51b_HDF = .FALSE.
!
!    ! Tracers to include
!    CALL SPLIT_ONE_LINE( SUBSTRS, Input_Opt%N_ND51b, -1, &
!                         'read_nd51b_menu:4',RC )
!    IF ( RC /= GC_SUCCESS ) THEN
!       CALL GC_Error( errMsg, RC, thisLoc )
!       RETURN
!    ENDIF
!    DO N = 1, Input_Opt%N_ND51b
!       READ( SUBSTRS(N), * ) Input_Opt%ND51b_TRACERS(N)
!    ENDDO
!
!    ! NHMS_WRITE
!    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1,  'read_nd51b_menu:5', RC )
!    IF ( RC /= GC_SUCCESS ) THEN
!       CALL GC_Error( errMsg, RC, thisLoc )
!       RETURN
!    ENDIF
!    READ( SUBSTRS(1:N), * ) Input_Opt%ND51b_HR_WRITE
!
!    ! Make sure ND51b_HR_WRITE is in the range 0-23.999 hrs
!    Input_Opt%ND51b_HR_WRITE = MOD(Input_Opt%ND51b_HR_WRITE,24e+0_fp)
!
!    ! HR1, HR2
!    CALL SPLIT_ONE_LINE( SUBSTRS, N, 2,  'read_nd51b_menu:6', RC )
!    IF ( RC /= GC_SUCCESS ) THEN
!       CALL GC_Error( errMsg, RC, thisLoc )
!       RETURN
!    ENDIF
!    READ( SUBSTRS(1:N), * ) Input_Opt%ND51b_HR1, Input_Opt%ND51b_HR2
!
!    ! IMIN, IMAX
!    CALL SPLIT_ONE_LINE( SUBSTRS, N, 2,  'read_nd51b_menu:7', RC )
!    IF ( RC /= GC_SUCCESS ) THEN
!       CALL GC_Error( errMsg, RC, thisLoc )
!       RETURN
!    ENDIF
!    READ( SUBSTRS(1:N), * ) Input_Opt%ND51b_IMIN, Input_Opt%ND51b_IMAX
!
!    ! JMIN, JMAX
!    CALL SPLIT_ONE_LINE( SUBSTRS, N, 2,  'read_nd51b_menu:8', RC )
!    IF ( RC /= GC_SUCCESS ) THEN
!       CALL GC_Error( errMsg, RC, thisLoc )
!       RETURN
!    ENDIF
!    READ( SUBSTRS(1:N), * ) Input_Opt%ND51b_JMIN, Input_Opt%ND51b_JMAX
!
!    ! LMIN, LMAX
!    CALL SPLIT_ONE_LINE( SUBSTRS, N, 2,  'read_nd51b_menu:9', RC )
!    IF ( RC /= GC_SUCCESS ) THEN
!       CALL GC_Error( errMsg, RC, thisLoc )
!       RETURN
!    ENDIF
!    READ( SUBSTRS(1:N), * ) Input_Opt%ND51b_LMIN, Input_Opt%ND51b_LMAX
!
!    ! Separator line
!    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1,  'read_nd51b_menu:10', RC )
!    IF ( RC /= GC_SUCCESS ) THEN
!       CALL GC_Error( errMsg, RC, thisLoc )
!       RETURN
!    ENDIF
!
!    ! Return success
!    RC = GC_SUCCESS
!
!    !=================================================================
!    ! Print to screen
!    !=================================================================
!    IF ( Input_Opt%amIRoot ) THEN
!       WRITE( 6, '(/,a)' )'ND51b MORNING OR AFTERNOON TIMESERIES MENU'
!       WRITE( 6, '(  a)' )'-----------------------------------------'
!       WRITE( 6, 100 ) 'Turn on ND51b timeseries?    : ', &
!                        Input_Opt%DO_ND51
!       WRITE( 6, 110 ) 'ND51b timeseries file name   : ', &
!                        TRIM( Input_Opt%ND51b_FILE )
!       WRITE( 6, 100 ) 'Output as HDF?               : ', &
!                        Input_Opt%LND51b_HDF
!       WRITE( 6, 120 ) 'ND51b timeseries tracers      : ',  &
!                        ( Input_Opt%ND51b_TRACERS(N), N=1, &
!                          Input_Opt%N_ND51b )
!       WRITE( 6, 140 ) 'ND51b hour to write to disk  : ', &
!                        Input_Opt%ND51b_HR_WRITE
!       WRITE( 6, 140 ) 'ND51b averaging period [GMT] : ', &
!                        Input_Opt%ND51b_HR1,  Input_Opt%ND51b_HR2
!       WRITE( 6, 130 ) 'ND51b longitude limits       : ', &
!                        Input_Opt%ND51b_IMIN, Input_Opt%ND51b_IMAX
!       WRITE( 6, 130 ) 'ND51b latitude  limits       : ', &
!                        Input_Opt%ND51b_JMIN, Input_Opt%ND51b_JMAX
!       WRITE( 6, 130 ) 'ND51b altitude  limits       : ', &
!                        Input_Opt%ND51b_LMIN, Input_Opt%ND51b_LMAX
!    ENDIF
!
!    ! FORMAT statements
!100 FORMAT( A, L5    )
!110 FORMAT( A, A     )
!120 FORMAT( A, 100I3 )
!130 FORMAT( A, 2I5   )
!140 FORMAT( A, 2F5.1 )
!
!  END SUBROUTINE READ_ND51b_MENU
!#ifdef TOMAS
!!EOC
!!------------------------------------------------------------------------------
!!                  GEOS-Chem Global Chemical Transport Model                  !
!!------------------------------------------------------------------------------
!!BOP
!!
!! !IROUTINE: read_prod_loss_menu
!!
!! !DESCRIPTION: Subroutine READ\_PROD\_LOSS\_MENU reads the PROD AND LOSS MENU
!!  section of the GEOS-Chem input file
!!\\
!!\\
!! !INTERFACE:
!!
!  SUBROUTINE READ_PROD_LOSS_MENU( Input_Opt, RC )
!!
!! !USES:
!!
!    USE CMN_DIAG_MOD
!    USE ErrCode_Mod
!    USE gckpp_Parameters,   ONLY : NFAM
!    USE gckpp_Monitor,      ONLY : FAM_NAMES
!    USE Input_Opt_Mod,      ONLY : OptInput
!!
!! !INPUT/OUTPUT PARAMETERS:
!!
!    TYPE(OptInput), INTENT(INOUT) :: Input_Opt   ! Input options
!!
!! !OUTPUT PARAMETERS:
!!
!    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure
!!
!! !REVISION HISTORY:
!!  20 Jul 2004 - R. Yantosca - Initial version
!!  See https://github.com/geoschem/geos-chem for complete history
!!EOP
!!------------------------------------------------------------------------------
!!BOC
!!
!! !LOCAL VARIABLES:
!!
!    INTEGER            :: F, N, N_ADVECT
!
!    ! Strings
!    CHARACTER(LEN=255) :: errMsg, thisLoc
!
!    ! Arrays
!    CHARACTER(LEN=255) :: SUBSTRS(MAXDIM)
!
!    !=================================================================
!    ! READ_PROD_LOSS_MENU begins here!
!    !=================================================================
!
!    ! Initialize
!    RC      = GC_SUCCESS
!
!    errMsg  = 'Error reading the "input.geos" file!'
!    thisLoc = ' -> at Read_Prod_Loss_Menu (in module GeosCore/input_mod.F90)'
!
!    ! Error check
!    IF ( CT1 /= 2 ) THEN
!       errMsg = 'SIMULATION MENU & ADVECTED SPECIES MENU ' // &
!                'must be read in first!'
!       CALL GC_Error( errMsg, RC, thisLoc )
!       RETURN
!    ENDIF
!
!    !=================================================================
!    ! Read info about prod & loss families
!    !=================================================================
!
!    ! Turn on production & loss diagnostic (e.g. ND65 diagnostic)
!    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'DO_SAVE_PL', RC )
!    IF ( RC /= GC_SUCCESS ) THEN
!       CALL GC_Error( errMsg, RC, thisLoc )
!       RETURN
!    ENDIF
!    READ( SUBSTRS(1:N), * ) Input_Opt%DO_SAVE_PL
!
!    ! Read number of levels for ND65 diagnostic
!    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'ND65', RC )
!    IF ( RC /= GC_SUCCESS ) THEN
!       CALL GC_Error( errMsg, RC, thisLoc )
!       RETURN
!    ENDIF
!    READ( SUBSTRS(1:N), * ) Input_Opt%ND65
!
!    ! Copy field to variable in CMN_DIAG
!    ND65 = Input_Opt%ND65
!
!    !=================================================================
!    ! Error check families for certain types of simulations
!    !=================================================================
!
!    ! Offline aerosol -- turn off DO_SAVE_PL, since we use ND05,
!    ! ND06, ND07, ND08, ND13 etc diagnostics instead of ND65
!    IF ( Input_Opt%ITS_AN_AEROSOL_SIM ) THEN
!       Input_Opt%DO_SAVE_PL    = .FALSE.
!       Input_Opt%ND65          = 0
!    ENDIF
!
!    !=================================================================
!    ! Set fields of Input Options object
!    !=================================================================
!
!    ! Number of advected species
!    N_ADVECT = Input_Opt%N_ADVECT
!
!    IF ( Input_Opt%DO_SAVE_PL ) THEN
!       IF ( Input_Opt%ITS_A_FULLCHEM_SIM ) THEN
!          ! Fullchem - Obtain NFAM from KPP
!          Input_Opt%NFAM = NFAM
!       ELSEIF ( Input_Opt%ITS_A_TAGO3_SIM ) THEN
!          ! Tagged O3
!          Input_Opt%NFAM = 2*N_ADVECT
!       ELSEIF ( Input_Opt%ITS_A_TAGCO_SIM ) THEN
!          ! Tagged CO
!          IF ( Input_Opt%LPCO_NMVOC ) THEN
!             Input_Opt%NFAM = N_ADVECT+2
!          ELSE
!             Input_Opt%NFAM = N_ADVECT+6
!          ENDIF
!       ENDIF
!    ENDIF
!
!    ! Return if there are no prod/loss families
!    ! or if we have turned off this diagnostic
!    IF ( .not. ( Input_Opt%DO_SAVE_PL .and. Input_Opt%NFAM > 0 )) THEN
!       Input_Opt%DO_SAVE_PL = .FALSE.
!       Input_Opt%ND65       = 0
!    ENDIF
!
!    ! Loop over families
!    DO F = 1, Input_Opt%NFAM
!
!       IF ( Input_Opt%ITS_A_FULLCHEM_SIM ) THEN
!
!          ! Fullchem - Obtain FAM_NAME from KPP
!          Input_Opt%FAM_NAME(F) = FAM_NAMES(F)
!
!       ELSEIF ( Input_Opt%ITS_A_TAGO3_SIM ) THEN
!
!          ! Tagged O3
!          IF ( F <= N_ADVECT ) THEN
!             Input_Opt%FAM_NAME(F) = &
!                  'P' // TRIM(Input_Opt%AdvectSpc_Name(F))
!          ELSE
!             Input_Opt%FAM_NAME(F) = &
!                  'L' // TRIM(Input_Opt%AdvectSpc_Name(F-N_ADVECT))
!          ENDIF
!
!       ELSEIF ( Input_Opt%ITS_A_TAGCO_SIM ) THEN
!
!          ! Tagged CO
!          IF ( F <= N_ADVECT ) THEN
!             Input_Opt%FAM_NAME(F) = 'L'//Input_Opt%AdvectSpc_Name(F)
!          ELSEIF ( F == N_ADVECT+1 ) THEN
!             Input_Opt%FAM_NAME(F) = 'PCO_CH4'
!          ELSEIF ( F == N_ADVECT+2 ) THEN
!             Input_Opt%FAM_NAME(F) = 'PCO_NMVOC'
!          ELSEIF ( F == N_ADVECT+3 ) THEN
!             Input_Opt%FAM_NAME(F) = 'PCO_ISOP'
!          ELSEIF ( F == N_ADVECT+4 ) THEN
!             Input_Opt%FAM_NAME(F) = 'PCO_CH3OH'
!          ELSEIF ( F == N_ADVECT+5 ) THEN
!             Input_Opt%FAM_NAME(F) = 'PCO_MONO'
!          ELSEIF ( F == N_ADVECT+6 ) THEN
!             Input_Opt%FAM_NAME(F) = 'PCO_ACET'
!          ENDIF
!
!       ENDIF
!
!       ! Get family type as prod or loss
!       IF ( Input_Opt%FAM_NAME(F)(1:1) == 'P'   .or. &
!            Input_Opt%FAM_NAME(F)(1:1) == 'p' ) THEN
!          Input_Opt%FAM_TYPE(F) = 'prod'
!       ELSE
!          Input_Opt%FAM_TYPE(F) = 'loss'
!       ENDIF
!
!    ENDDO
!
!    !=================================================================
!    ! Print to screen
!    !=================================================================
!    IF ( Input_Opt%amIRoot ) THEN
!       WRITE( 6, '(/,a)' ) 'PROD & LOSS DIAGNOSTIC MENU'
!       WRITE( 6, '(  a)' ) '---------------------------'
!       WRITE( 6, 100 ) 'Turn on prod & loss diag?   : ', &
!                        Input_Opt%DO_SAVE_PL
!       WRITE( 6, 110 ) '# of levels for P/L diag    : ', &
!                        Input_Opt%ND65
!
!       ! Loop over families
!       DO F = 1, Input_Opt%NFAM
!
!          ! Write family name and type
!          WRITE( 6, 120 ) TRIM(Input_Opt%FAM_NAME(F)), &
!                          TRIM(Input_Opt%FAM_TYPE(F))
!
!       ENDDO
!
!    ENDIF
!
!    ! FORMAT statements
!100 FORMAT( A, L5 )
!110 FORMAT( A, I5 )
!120 FORMAT( /, 'Family=', A10, '  Type=', A4 )
!
!  END SUBROUTINE READ_PROD_LOSS_MENU
!!EOC
!#endif
!#endif
!#endif
!!------------------------------------------------------------------------------
!!                  GEOS-Chem Global Chemical Transport Model                  !
!!------------------------------------------------------------------------------
!!BOP
!!
!! !IROUTINE: read_mercury_menu
!!
!! !DESCRIPTION: Subroutine READ\_MERCURY\_MENU reads the BENCHMARK MENU
!!  section of the GEOS-Chem input file.
!!\\
!!\\
!! !INTERFACE:
!!
!  SUBROUTINE READ_MERCURY_MENU( Input_Opt, RC )
!!
!! !USES:
!!
!    USE ErrCode_Mod
!    USE Input_Opt_Mod, ONLY : OptInput
!!
!! !INPUT/OUTPUT PARAMETERS:
!!
!    TYPE(OptInput), INTENT(INOUT) :: Input_Opt   ! Input options
!!
!! !OUTPUT PARAMETERS:
!!
!    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure
!!
!! !REVISION HISTORY:
!!  24 Feb 2006 - R. Yantosca - Initial version
!!  See https://github.com/geoschem/geos-chem for complete history
!!EOP
!!------------------------------------------------------------------------------
!!BOC
!!
!! !LOCAL VARIABLES:
!!
!    ! Scalars
!    LOGICAL            :: LDYNOCEAN,      LPREINDHG
!    LOGICAL            :: LGTMM,          USE_CHECKS
!    LOGICAL            :: LARCTICRIV,     LKRedUV
!    INTEGER            :: N
!    CHARACTER(LEN=255) :: GTMM_RST_FILE
!
!    ! Strings
!    CHARACTER(LEN=255) :: errMsg,         thisLoc
!
!    ! Arrays
!    CHARACTER(LEN=255) :: SUBSTRS(MAXDIM)
!
!    !=================================================================
!    ! READ_MERCURY_MENU begins here!
!    !=================================================================
!
!    ! Initialize
!    RC      = GC_SUCCESS
!    errMsg  = 'Error reading the "input.geos" file!'
!    thisLoc = ' -> at Read_Mercury_Menu (in module GeosCore/input_mod.F90)'
!
!    ! Use error check for tag/tot Hg?
!    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'USE_CHECKS', RC )
!    IF ( RC /= GC_SUCCESS ) THEN
!       CALL GC_Error( errMsg, RC, thisLoc )
!       RETURN
!    ENDIF
!    READ( SUBSTRS(1:N), * ) Input_Opt%USE_CHECKS
!
!    ! Use dynamic ocean model?
!    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'LDYNOCEAN', RC )
!    IF ( RC /= GC_SUCCESS ) THEN
!       CALL GC_Error( errMsg, RC, thisLoc )
!       RETURN
!    ENDIF
!    READ( SUBSTRS(1:N), * ) Input_Opt%LDYNOCEAN
!
!    ! Use preindustrial simulation?
!    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'LPREINDHG', RC )
!    IF ( RC /= GC_SUCCESS ) THEN
!       CALL GC_Error( errMsg, RC, thisLoc )
!       RETURN
!    ENDIF
!    READ( SUBSTRS(1:N), * ) Input_Opt%LPREINDHG
!
!    ! Use GTMM ?
!    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'LGTMM', RC )
!    IF ( RC /= GC_SUCCESS ) THEN
!       CALL GC_Error( errMsg, RC, thisLoc )
!       RETURN
!    ENDIF
!    READ( SUBSTRS(1:N), * ) Input_Opt%LGTMM
!
!    ! Name of GTMM restart file
!    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'GTMM_RST_FILE', RC )
!    IF ( RC /= GC_SUCCESS ) THEN
!       CALL GC_Error( errMsg, RC, thisLoc )
!       RETURN
!    ENDIF
!    READ( SUBSTRS(1:N), * ) Input_Opt%GTMM_RST_FILE
!
!    ! Use Arctic river Hg?
!    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'LARCTICRIV', RC )
!    IF ( RC /= GC_SUCCESS ) THEN
!       CALL GC_Error( errMsg, RC, thisLoc )
!       RETURN
!    ENDIF
!    READ( SUBSTRS(1:N), * ) Input_Opt%LARCTICRIV
!
!    ! Tie reducible HgII(aq) to UV-B?
!    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'LKREDUV', RC )
!    IF ( RC /= GC_SUCCESS ) THEN
!       CALL GC_Error( errMsg, RC, thisLoc )
!       RETURN
!    ENDIF
!    READ( SUBSTRS(1:N), * ) Input_Opt%LKRedUV
!
!    ! Separator line
!    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'separator', RC )
!    IF ( RC /= GC_SUCCESS ) THEN
!       CALL GC_Error( errMsg, RC, thisLoc )
!       RETURN
!    ENDIF
!
!    ! Check on logical
!    IF (.NOT.( Input_Opt%ITS_A_MERCURY_SIM ) ) THEN
!       Input_Opt%LGTMM      = .FALSE.
!       Input_Opt%LDYNOCEAN  = .FALSE.
!       Input_Opt%LARCTICRIV = .FALSE.
!       Input_Opt%LKRedUV    = .FALSE.
!    ENDIF
!
!    !=================================================================
!    ! Print to screen
!    !=================================================================
!    IF ( Input_Opt%amIRoot ) THEN
!       WRITE( 6, '(/,a)' ) 'MERCURY MENU'
!       WRITE( 6, '(  a)' ) '------------'
!       WRITE( 6, 110 ) 'Error check tag & total Hg? : ', &
!                        Input_Opt%USE_CHECKS
!       WRITE( 6, 110 ) 'Use dynamic ocean Hg model? : ', &
!                        Input_Opt%LDYNOCEAN
!       WRITE( 6, 110 ) 'Preindustrial simulation?   : ', &
!                        Input_Opt%LPREINDHG
!       WRITE( 6, 110 ) 'Use GTMM ?                  : ', &
!                        Input_Opt%LGTMM
!       WRITE( 6, 120 ) '=> GTMM restart file        : ', &
!                        TRIM( Input_Opt%GTMM_RST_FILE )
!       WRITE( 6, 110 ) 'Use Arctic river Hg ?       : ', &
!                        Input_Opt%LARCTICRIV
!       WRITE( 6, 110 ) 'Tie HgII(aq) red. to UV-B?  : ', &
!                        Input_Opt%LKRedUV
!    ENDIF
!
!    ! FORMAT statements
!100 FORMAT( A, I4  )
!110 FORMAT( A, L5  )
!120 FORMAT( A, A   )
!
!  END SUBROUTINE READ_MERCURY_MENU
!!EOC
!!------------------------------------------------------------------------------
!!                  GEOS-Chem Global Chemical Transport Model                  !
!!------------------------------------------------------------------------------
!!BOP
!!
!! !IROUTINE: read_ch4_menu
!!
!! !DESCRIPTION: Subroutine READ\_CH4\_MENU reads the CH4 MENU section of
!!  the GEOS-Chem input file; this defines options for CH4 simulations.
!!\\
!!\\
!! !INTERFACE:
!!
!  SUBROUTINE READ_CH4_MENU( Input_Opt, RC )
!!
!! !USES:
!!
!    USE ErrCode_Mod
!    USE Input_Opt_Mod, ONLY : OptInput
!!
!! !INPUT/OUTPUT PARAMETERS:
!!
!    TYPE(OptInput), INTENT(INOUT) :: Input_Opt   ! Input options
!!
!! !OUTPUT PARAMETERS:
!!
!    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!!
!! !REVISION HISTORY:
!!  03 Aug 2009 - K. Wecht, C. Pickett-Heaps - Initial version
!!  See https://github.com/geoschem/geos-chem for complete history
!!EOP
!!------------------------------------------------------------------------------
!!BOC
!!
!! !LOCAL VARIABLES:
!!
!    ! Scalars
!    INTEGER            :: N
!
!    ! Strings
!    CHARACTER(LEN=255) :: errMsg, thisLoc
!
!    ! Arrays
!    CHARACTER(LEN=255) :: SUBSTRS(MAXDIM)
!
!    !=================================================================
!    ! READ_CH4_MENU begins here!
!    !=================================================================
!
!    ! Initialize
!    RC      = GC_SUCCESS
!    errMsg  = 'Error reading the "input.geos" file!'
!    thisLoc = ' -> at Read_CH4_Menu (in module GeosCore/input_mod.F90)'
!
!    ! Use GOSAT CH4 observation operator?
!    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'GOSAT_CH4_OBS', RC )
!    IF ( RC /= GC_SUCCESS ) THEN
!       CALL GC_Error( errMsg, RC, thisLoc )
!       RETURN
!    ENDIF
!    READ( SUBSTRS(1:N), * ) Input_Opt%GOSAT_CH4_OBS
!
!    ! Use TCCON CH4 observation operator?
!    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'TCCON_CH4_OBS', RC )
!    IF ( RC /= GC_SUCCESS ) THEN
!       CALL GC_Error( errMsg, RC, thisLoc )
!       RETURN
!    ENDIF
!    READ( SUBSTRS(1:N), * ) Input_Opt%TCCON_CH4_OBS
!
!    ! Separator line
!    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'separator', RC )
!    IF ( RC /= GC_SUCCESS ) THEN
!       CALL GC_Error( errMsg, RC, thisLoc )
!       RETURN
!    ENDIF
!
!    !=================================================================
!    ! Print to screen
!    !=================================================================
!    IF ( Input_Opt%amIRoot ) THEN
!       WRITE( 6, '(/,a)' ) 'CH4 MENU'
!       WRITE( 6, '(  a)' ) '-----------'
!       WRITE( 6, 100     ) 'Use GOSAT obs operator: ', &
!                            Input_Opt%GOSAT_CH4_OBS
!       WRITE( 6, 100     ) 'Use TCCON obs operator: ', &
!            Input_Opt%TCCON_CH4_OBS
!    ENDIF
!
!    ! FORMAT statements
!100 FORMAT( A, L5  )
!
!    ! Return success
!    RC = GC_SUCCESS
!
!  END SUBROUTINE READ_CH4_MENU
!!EOC
!!------------------------------------------------------------------------------
!!                  GEOS-Chem Global Chemical Transport Model                  !
!!------------------------------------------------------------------------------
!!BOP
!!
!! !IROUTINE: read_pops_menu
!!
!! !DESCRIPTION: Subroutine READ\_POPS\_MENU reads the POPS MENU section of
!!  the GEOS-Chem input file; this defines options for POPs simulations.
!!\\
!!\\
!! !INTERFACE:
!!
!  SUBROUTINE READ_POPS_MENU( Input_Opt, RC )
!!
!! !USES:
!!
!    USE ErrCode_Mod
!    USE Input_Opt_Mod, ONLY : OptInput
!!
!! !INPUT/OUTPUT PARAMETERS:
!!
!    TYPE(OptInput), INTENT(INOUT) :: Input_Opt   ! Input options
!!
!! !OUTPUT PARAMETERS:
!!
!    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!!
!! !REVISION HISTORY:
!!  01 Oct 2012 - C. Friedman - Initial version
!!  See https://github.com/geoschem/geos-chem for complete history
!!EOP
!!------------------------------------------------------------------------------
!!BOC
!!
!! !LOCAL VARIABLES:
!!
!    ! Scalars
!    INTEGER            :: N
!
!    ! Strings
!    CHARACTER(LEN=255) :: errMsg, thisLoc
!
!    ! Arrays
!    CHARACTER(LEN=255) :: SUBSTRS(MAXDIM)
!
!    !=================================================================
!    ! READ_POPS_MENU begins here!
!    !=================================================================
!
!    ! Initialize
!    RC      = GC_SUCCESS
!    errMsg  = 'Error reading the "input.geos" file!'
!    thisLoc = ' -> at Read_POPs_Menu (in module GeosCore/input_mod.F90)'
!
!    ! POP species
!    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'POP_TYPE', RC )
!    IF ( RC /= GC_SUCCESS ) THEN
!       CALL GC_Error( errMsg, RC, thisLoc )
!       RETURN
!    ENDIF
!    READ( SUBSTRS(1:N), * ) Input_Opt%POP_TYPE
!
!    ! Dummy for future process logical switches
!    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'CHEM_PROCESS', RC )
!    IF ( RC /= GC_SUCCESS ) THEN
!       CALL GC_Error( errMsg, RC, thisLoc )
!       RETURN
!    ENDIF
!    READ( SUBSTRS(1:N), * ) Input_Opt%CHEM_PROCESS
!
!    ! Molecular weight
!    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'POP_XMW', RC )
!    IF ( RC /= GC_SUCCESS ) THEN
!       CALL GC_Error( errMsg, RC, thisLoc )
!       RETURN
!    ENDIF
!    READ( SUBSTRS(1:N), * ) Input_Opt%POP_XMW
!
!    ! KOA
!    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'POP_KOA', RC )
!    IF ( RC /= GC_SUCCESS ) THEN
!       CALL GC_Error( errMsg, RC, thisLoc )
!       RETURN
!    ENDIF
!    READ( SUBSTRS(1:N), * ) Input_Opt%POP_KOA
!
!    ! KBC
!    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'POP_KBC', RC )
!    IF ( RC /= GC_SUCCESS ) THEN
!       CALL GC_Error( errMsg, RC, thisLoc )
!       RETURN
!    ENDIF
!    READ( SUBSTRS(1:N), * ) Input_Opt%POP_KBC
!
!    ! OH oxidation
!    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, ' POP_K_POPG_OH', RC )
!    IF ( RC /= GC_SUCCESS ) THEN
!       CALL GC_Error( errMsg, RC, thisLoc )
!       RETURN
!    ENDIF
!    READ( SUBSTRS(1:N), * ) Input_Opt%POP_K_POPG_OH
!
!    ! O3 oxidation 1
!    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, ' POP_K_POPP_O3A', RC )
!    IF ( RC /= GC_SUCCESS ) THEN
!       CALL GC_Error( errMsg, RC, thisLoc )
!       RETURN
!    ENDIF
!    READ( SUBSTRS(1:N), * ) Input_Opt%POP_K_POPP_O3A
!
!    ! O3 oxidation 2
!    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'POP_K_POPP_O3B', RC )
!    IF ( RC /= GC_SUCCESS ) THEN
!       CALL GC_Error( errMsg, RC, thisLoc )
!       RETURN
!    ENDIF
!    READ( SUBSTRS(1:N), * ) Input_Opt%POP_K_POPP_O3B
!
!    ! H*
!    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'POP_HSTAR', RC )
!    IF ( RC /= GC_SUCCESS ) THEN
!       CALL GC_Error( errMsg, RC, thisLoc )
!       RETURN
!    ENDIF
!    READ( SUBSTRS(1:N), * ) Input_Opt%POP_HSTAR
!
!    ! DEL_H
!    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'POP_DEL_H', RC )
!    IF ( RC /= GC_SUCCESS ) THEN
!       CALL GC_Error( errMsg, RC, thisLoc )
!       RETURN
!    ENDIF
!    READ( SUBSTRS(1:N), * ) Input_Opt%POP_DEL_H
!
!    ! DEL_Hw
!    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'POP_DEL_Hw', RC )
!    IF ( RC /= GC_SUCCESS ) THEN
!       CALL GC_Error( errMsg, RC, thisLoc )
!       RETURN
!    ENDIF
!    READ( SUBSTRS(1:N), * ) Input_Opt%POP_DEL_Hw
!
!    ! Separator line
!    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'separator', RC )
!    IF ( RC /= GC_SUCCESS ) THEN
!       CALL GC_Error( errMsg, RC, thisLoc )
!       RETURN
!    ENDIF
!
!    !=================================================================
!    ! Print to screen
!    !=================================================================
!    IF ( Input_Opt%amIRoot ) THEN
!       WRITE( 6, '(/,a)' ) 'POPS MENU'
!       WRITE( 6, '(  a)' ) '------------'
!       WRITE( 6, 120     ) 'Species of POP        : ', &
!                            Input_Opt%POP_TYPE
!       WRITE( 6, 110     ) 'Chemistry on?         : ', &
!                            Input_Opt%CHEM_PROCESS
!       WRITE( 6, 130     ) 'POP_XMW               : ', &
!                            Input_Opt%POP_XMW
!       WRITE( 6, 130     ) 'POP_KOA               : ', &
!                            Input_Opt%POP_KOA
!       WRITE( 6, 130     ) 'POP_KBC               : ', &
!                            Input_Opt%POP_KBC
!       WRITE( 6, 130     ) 'POP_K_POPG_OH         : ', &
!                            Input_Opt%POP_K_POPG_OH
!       WRITE( 6, 130     ) 'POP_K_POPP_O3A        : ', &
!                            Input_Opt%POP_K_POPP_O3A
!       WRITE( 6, 130     ) 'POP_K_POPP_O3B        : ', &
!                            Input_Opt%POP_K_POPP_O3B
!       WRITE( 6, 130     ) 'POP_HSTAR             : ', &
!                            Input_Opt%POP_HSTAR
!       WRITE( 6, 130     ) 'POP_DEL_H             : ', &
!                            Input_Opt%POP_DEL_H
!       WRITE( 6, 130     ) 'POP_DEL_Hw            : ', &
!                            Input_Opt%POP_DEL_Hw
!    ENDIF
!
!    ! FORMAT statements
!110 FORMAT( A, L5  )
!120 FORMAT( A, A   )
!130 FORMAT( A, ES10.2 )
!
!    ! Return success
!    RC = GC_SUCCESS
!
!  END SUBROUTINE READ_POPS_MENU
!!EOC
!!------------------------------------------------------------------------------
!!                  GEOS-Chem Global Chemical Transport Model                  !
!!------------------------------------------------------------------------------
!!BOP
!!
!! !IROUTINE: read_passive_species_menu
!!
!! !DESCRIPTION: Subroutine READ\_PASSIVE\_SPECIES\_MENU reads the passive
!!  species menu section of the GEOS-Chem input file; this defines passive
!!  species to be used for this simulation.
!!\\
!!\\
!! !INTERFACE:
!!
!  SUBROUTINE READ_PASSIVE_SPECIES_MENU( Input_Opt, RC )
!!
!! !USES:
!!
!    USE ErrCode_Mod
!    USE Input_Opt_Mod, ONLY : OptInput
!!
!! !INPUT/OUTPUT PARAMETERS:
!!
!    TYPE(OptInput), INTENT(INOUT) :: Input_Opt   ! Input options
!!
!! !OUTPUT PARAMETERS:
!!
!    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!!
!! !REVISION HISTORY:
!!  04 Sep 2015 - C. Keller   - Initial version
!!  See https://github.com/geoschem/geos-chem for complete history
!!EOP
!!------------------------------------------------------------------------------
!!BOC
!!
!! !LOCAL VARIABLES:
!!
!    ! Scalars
!    INTEGER            :: N, P, D
!
!    ! Strings
!    CHARACTER(LEN=255) :: errMsg, thisLoc
!
!    ! Arrays
!    CHARACTER(LEN=255) :: SUBSTRS(MAXDIM)
!
!    !=================================================================
!    ! READ_PASSIVE_SPECIES_MENU begins here!
!    !=================================================================
!
!    ! Initialize
!    RC      = GC_SUCCESS
!    errMsg  = 'Error reading the "input.geos" file!'
!    thisLoc = ' -> at Read_Passive_Species_Menu (in GeosCore/input_mod.F90)'
!
!    ! Initialize
!    P = 0
!    D = 0
!
!    IF ( Input_Opt%amIRoot ) THEN
!       WRITE( 6, '(/,a)' ) 'PASSIVE SPECIES MENU'
!       WRITE( 6, '(  a)' ) '---------------------'
!    ENDIF
!
!    ! Do until exit
!    DO
!
!       ! Read passive species information for each passive species
!       ! Every passive species line is expected to have 4 entries:
!       ! - Species name
!       ! - Species molecular weight
!       ! - Atmospheric lifetime (s)
!       ! - Initial atmospheric concentration (v/v)
!       CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'Passive Species', RC )
!       IF ( RC /= GC_SUCCESS ) THEN
!          CALL GC_Error( errMsg, RC, thisLoc )
!          RETURN
!       ENDIF
!
!       ! Exit here if separator line is encountered
!       IF ( INDEX ( TRIM(SUBSTRS(1)), '-----' ) > 0 ) EXIT
!
!       ! Make sure there are at least 4 entries
!       IF ( N < 4 ) THEN
!          errMsg = 'Each passive species is expected to have '      // &
!                   'at least four entries: Name, MW, TAU, initial ' // &
!                   'concentration, and long name (optional).'
!          CALL GC_Error( errMsg, RC, thisLoc )
!          RETURN
!       ENDIF
!
!       ! Increase number of passive species by one
!       P = P + 1
!
!       ! Stop simulation and print warning if we exceed maximum number
!       ! of passive species hardcoded in input_opt_mod.F90.
!       IF ( P > Input_Opt%Max_PassiveSpc ) THEN
!          errMsg = 'Number of passive species exceeds maximum. ' // &
!                   'This value can be modified in input_opt_mod.F90.'
!          CALL GC_Error( errMsg, RC, thisLoc )
!          RETURN
!       ENDIF
!
!       ! Read and store species information
!       Input_Opt%PASSIVE_ID(P)   = P
!       Input_Opt%PASSIVE_NAME(P) = TRIM( SUBSTRS(1) )
!       READ( SUBSTRS(2), * ) Input_Opt%PASSIVE_MW(P)
!       READ( SUBSTRS(3), * ) Input_Opt%PASSIVE_TAU(P)
!       READ( SUBSTRS(4), * ) Input_Opt%PASSIVE_INITCONC(P)
!
!       ! Check if optional entry for long name is included
!       IF ( N > 4 ) THEN
!          READ( SUBSTRS(5), * ) Input_Opt%PASSIVE_LONGNAME(P)
!       ENDIF
!
!       ! Determine if the passive species decays (i.e. if it has
!       ! an atmospheric lifetime that is not -1).  This will allow
!       ! us to skip those passive species that do not decay in
!       ! routine CHEM_PASSIVE_SPECIES, to speed up execution.
!       IF ( Input_Opt%PASSIVE_TAU(P) > 0.0_fp ) THEN
!          D                            = D + 1
!          Input_Opt%NPASSIVE_DECAY     = D
!          Input_Opt%PASSIVE_DECAYID(D) = P
!       ENDIF
!
!       ! Verbose
!       IF ( Input_Opt%amIRoot ) THEN
!          WRITE( 6, '(a)' ) 'Added passive species: '
!          WRITE( 6, 110   ) ' - Species name                : ', &
!                TRIM( Input_Opt%PASSIVE_NAME(P) )
!          WRITE( 6, 120   ) ' - Molec. weight [g/mol]       : ', &
!                Input_Opt%PASSIVE_MW(P)
!          WRITE( 6, 130   ) ' - Lifetime [s]                : ', &
!                Input_Opt%PASSIVE_TAU(P)
!          WRITE( 6, 130   ) ' - Initial concentration [v/v] : ', &
!                Input_Opt%PASSIVE_INITCONC(P)
!          WRITE( 6, 110   ) ' - Species long name           : ', &
!                TRIM( Input_Opt%PASSIVE_LONGNAME(P) )
!       ENDIF
!
!    ENDDO
!
!    IF ( P < 0 ) THEN
!       errMsg = 'Cannot add passive species '     // &
!                TRIM(Input_Opt%PASSIVE_NAME(P) ) // &
!                ': # of passive species is smaller than 1!'
!       CALL GC_Error( errMsg, RC, thisLoc )
!       RETURN
!    ENDIF
!
!    ! Total number of passive species
!    Input_Opt%NPASSIVE = P
!
!110 FORMAT( A, A )
!120 FORMAT( A, F10.2  )
!130 FORMAT( A, ES10.2 )
!
!    ! Return success
!    RC = GC_SUCCESS
!
!  END SUBROUTINE READ_PASSIVE_SPECIES_MENU
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: validate_directories
!
! !DESCRIPTION: Subroutine VALIDATE\_DIRECTORIES makes sure that each of the
!  directories that we have read from the GEOS-Chem input file are valid.
!  Also, trailing separator characters will be added.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE VALIDATE_DIRECTORIES( Input_Opt, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod, ONLY : OptInput
    USE Time_Mod,      ONLY : Expand_Date
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(INOUT) :: Input_Opt   ! Input Options object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REVISION HISTORY:
!  20 Jul 2004 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Strings
    CHARACTER(LEN=255) :: errMsg, thisLoc, Dir

    !=================================================================
    ! VALIDATE_DIRECTORIES begins here!
    !=================================================================

    ! Initialize
    RC      = GC_SUCCESS
    errMsg  = 'Invalid directory encountered!'
    thisLoc = ' -> at Validate_Directories (in module GeosCore/input_mod.F90)'

    ! Skip for dry-runs
    IF ( Input_Opt%DryRun ) RETURN

    ! Check directories
    CALL Check_Directory( Input_Opt, Input_Opt%DATA_DIR, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    CALL Check_Directory( Input_Opt, Input_Opt%CHEM_INPUTS_DIR, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    CALL Check_Directory( Input_Opt, Input_Opt%RUN_DIR, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

  END SUBROUTINE Validate_Directories
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: check_directory
!
! !DESCRIPTION: Subroutine CHECK\_DIRECTORY makes sure that the given
!  directory is valid.  Also a trailing slash character will be added if
!  necessary.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Check_Directory( Input_Opt, dir, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE FILE_MOD,      ONLY : File_Exists
    USE Input_Opt_Mod, ONLY : OptInput
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(OptInput),   INTENT(INOUT) :: Input_Opt   ! Input Options object
    CHARACTER(LEN=*), INTENT(INOUT) :: dir         ! Dir to be checked
!
! !OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(OUT)   :: RC          ! Success or failure
!
! !REVISION HISTORY:
!  20 Mar 2003 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER            :: C

    ! Strings
    CHARACTER(LEN=255) :: errMsg, thisLoc

    !=================================================================
    ! CHECK_DIRECTORY begins here!
    !=================================================================

    ! Initialize
    RC      = GC_SUCCESS
    errMsg  = ''
    thisLoc = ' -> at Check_Directory (in module GeosCore/input_mod.F90)'

    ! Locate the last non-white-space character of NEWDIR
    C = LEN_TRIM( dir )

    ! Add the trailing directory separator if it is not present
    IF ( dir(C:C) /= '/' ) THEN
       dir(C+1:C+1) = TRIM( '/' )
    ENDIF

    !=================================================================
    ! Test if the directory actually exists
    !=================================================================

    ! If the directory does not exist then stop w/ an error message
    IF ( .not. File_Exists( dir ) ) THEN
       errMsg = 'Invalid directory: ' // TRIM( dir )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

  END SUBROUTINE Check_Directory
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: check_time_steps
!
! !DESCRIPTION: Subroutine CHECK\_TIME\_STEPS computes the smallest dynamic
!  time step for the model, based on which operation are turned on.  This
!  is called from routine READ\_INPUT\_FILE, after all of the timesteps and
!  logical flags have been read from "input.geos".
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Check_Time_Steps( Input_Opt, State_Grid, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod,  ONLY : OptInput
    USE State_Grid_Mod, ONLY : GrdState
    USE Time_Mod,       ONLY : Set_TimeSteps
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)  :: Input_Opt   ! Input Options object
    TYPE(GrdState), INTENT(IN)  :: State_Grid  ! Grid State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT) :: RC          ! Success or failure?
!
! !REVISION HISTORY:
!  20 Jul 2004 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    LOGICAL            :: LCONV, LCHEM,       LDRYD
    LOGICAL            :: LTRAN, LTURB
    INTEGER            :: I,     J,           K
    INTEGER            :: L,     TS_SMALLEST, TS_DIAG
    INTEGER            :: TS_CHEM, TS_EMIS, TS_CONV, TS_DYN
    INTEGER            :: TS_UNIT, TS_RAD,  MAX_DYN

    ! Strings
    CHARACTER(LEN=255) :: errMsg, thisLoc

    !=================================================================
    ! CHECK_TIME_STEPS begins here!
    !=================================================================

    ! Initialize
    RC      = GC_SUCCESS
    errMsg  = ''
    thisLoc = ' -> at Check_Time_Steps (in module GeosCore/input_mod.F90)'

    ! Copy fields from Input_Opt
    LCONV = Input_Opt%LCONV
    LCHEM = Input_Opt%LCHEM
    LDRYD = Input_Opt%LDRYD
    LTRAN = Input_Opt%LTRAN
    LTURB = Input_Opt%LTURB

    TS_CHEM = Input_Opt%TS_CHEM
    TS_EMIS = Input_Opt%TS_EMIS
    TS_CONV = Input_Opt%TS_CONV
    TS_DYN  = Input_Opt%TS_DYN
    TS_RAD  = Input_Opt%TS_RAD

    ! If we're doing the reverse integration
    ! multiply all the timesteps by -1 here
    if (TS_DYN < 0) THEN
       TS_CHEM = TS_CHEM * -1
       TS_EMIS = TS_EMIS * -1
       TS_CONV = TS_CONV * -1
       TS_DYN  = TS_DYN  * -1
       TS_RAD  = TS_RAD  * -1
    endif


    ! NUNIT is time step in minutes for unit conversion
    TS_UNIT = -1

    ! Define maximum timestep for transport
    IF ( TRIM(State_Grid%GridRes) == '4.0x5.0') THEN
       MAX_DYN = 1800
    ELSE IF ( TRIM(State_Grid%GridRes) == '2.0x2.5' ) THEN
       MAX_DYN = 900
    ELSE IF ( TRIM(State_Grid%GridRes) == '0.5x0.625' ) THEN
       MAX_DYN = 600
    ELSE IF ( TRIM(State_Grid%GridRes) == '0.25x0.3125' ) THEN
       MAX_DYN = 300
    ELSE
       MAX_DYN = 3600
    ENDIF

    ! If TS_DYN is greater than MAX_DYN, then stop w/ error
    IF ( .not. Input_Opt%isMPI ) THEN
       IF ( Input_Opt%TS_DYN > MAX_DYN .and. LTRAN ) THEN
          WRITE( errMsg, 300 ) 'Transport timestep exceeds max:', &
                                Input_Opt%TS_DYN, MAX_DYN
300       FORMAT( a, i8, ' >', i8 )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
    ENDIF

    ! Only do unit conversion if necessary
    IF ( LTRAN .or. LCONV .or. LTURB ) THEN
       TS_UNIT = MAX( TS_DYN, TS_CONV )
    ENDIF

    ! Compute NSMALLEST as the minimum of NDYN, NCONV, NSRCE, NCHEM
    I = TS_DYN
    J = TS_CONV
    K = TS_EMIS
    L = TS_CHEM

    ! SDE 2017-02-24: Always use LTRAN on the assumption that it will
    ! be used as a "heartbeat". This ensures that chemistry always
    ! takes place at the same time, regardless of whether or not
    ! transport is enabled.
    !IF ( .not. LTRAN                  ) I = 999999
    IF ( .not. LCONV .and..not. LTURB ) J = 999999
    IF ( .not. LDRYD                  ) K = 999999
    IF ( .not. LCHEM                  ) L = 999999

    ! Get the smallest of all of the above
    TS_SMALLEST = MIN( I, J, K, L )

    ! If all of the operators above are turned off,
    ! then set TS_SMALLEST to TS_DYN.
    IF ( TS_SMALLEST == 999999 ) THEN
       TS_SMALLEST = TS_DYN
    ENDIF

    IF ( LTRAN .and. TS_DYN /= TS_SMALLEST ) THEN
       errMsg = 'The transport time step should be the smallest one'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    ! If TS_DYN is smaller than TS_SMALLEST, reset TS_DYN
    ! to TS_SMALLEST.
    ! This is useful for runs where transport is turned off,
    ! but where chemistry is turned on.
    IF ( TS_DYN < TS_SMALLEST ) THEN
       TS_DYN = TS_SMALLEST
    ENDIF

    ! Define the largest time step, TS_DIAG, for diagnostics.
    ! Diagnostics should be incremented at the end of multiples of
    ! TS_DIAG, so that the system is at a physical state.
    ! (ccc, 5/13/09)
    IF ( .not. LTRAN                  ) I = -999999
    IF ( .not. LCONV .and..not. LTURB ) J = -999999
    IF ( .not. LDRYD                  ) K = -999999
    IF ( .not. LCHEM                  ) L = -999999

    TS_DIAG = MAX( I, J, K, L )

    ! If all the operators are turned off, then set TS_DIAG to TS_CHEM
    ! Usually the chemistry time step is large. (ccc, 5/13/09)
    IF ( TS_DIAG == -999999 ) THEN
       TS_DIAG = TS_CHEM
    ENDIF

    ! Check if all time steps are multiples of the smallest.
    ! (ccc, 5/13/09)
    IF ( L /= -999999 .and. MOD( TS_CHEM, TS_SMALLEST ) /= 0 ) THEN
       WRITE( errMsg, 100 ) 'Chemistry', TS_CHEM, TS_SMALLEST
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    IF ( K /= -999999 .and. MOD( TS_EMIS, TS_SMALLEST ) /= 0 ) THEN
       WRITE( ErrMSg, 100 ) 'Emission', TS_EMIS, TS_SMALLEST
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    IF ( J /= -999999 .and. MOD( TS_CONV, TS_SMALLEST ) /= 0 ) THEN
       WRITE( errMsg, 100 ) 'Convection', TS_CONV, TS_SMALLEST
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    IF ( I /= -999999 .and. MOD( TS_DYN, TS_SMALLEST ) /= 0 ) THEN
       WRITE( errMsg, 100 ) 'Transport', TS_DYN, TS_SMALLEST
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    ! Initialize timesteps in "time_mod.F90"
    CALL Set_Timesteps( Input_Opt,                                           &
                        CHEMISTRY  = TS_CHEM,                                &
                        EMISSION   = TS_EMIS,                                &
                        DYNAMICS   = TS_DYN,                                 &
                        UNIT_CONV  = TS_UNIT,                                &
                        CONVECTION = TS_CONV,                                &
                        DIAGNOS    = TS_DIAG,                                &
                        RADIATION  = TS_RAD )

100 FORMAT( A, ' time step must be a multiple of the smallest one:', i5, i5 )

  END SUBROUTINE Check_Time_Steps
!EOC
#if !(defined( EXTERNAL_GRID ) || defined( EXTERNAL_FORCING ))
#ifdef BPCH_DIAG
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: is_last_day_good
!
! !DESCRIPTION: Suborutine IS\_LAST\_DAY\_GOOD tests to see if there is
!  output scheduled on the last day of the run.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE IS_LAST_DAY_GOOD( Input_Opt, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod, ONLY : OptInput
    USE JULDAY_MOD,    ONLY : JULDAY
    USE TIME_MOD,      ONLY : GET_NYMDe, ITS_A_LEAPYEAR, YMD_EXTRACT
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(INOUT) :: Input_Opt  ! Input options
    INTEGER,        INTENT(OUT)   :: RC         ! Success or failure?
!
! !REVISION HISTORY:
!  20 Jul 2004 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    LOGICAL            :: IS_LEAPYEAR
    INTEGER            :: NYMDe, Y, M, D, LASTDAY
    REAL(fp)           :: JD, JD0

    ! Strings
    CHARACTER(LEN=255) :: errMsg, thisLoc

    !=================================================================
    ! IS_LAST_DAY_GOOD begins here!
    !=================================================================

    ! Initialize
    RC      = GC_SUCCESS
    errMsg  = ''
    thisLoc = ' -> at Is_Last_Day_Good (in module GeosCore/input_mod.F90)'

    ! Astronomical Julian Day corresponding to NYMDe
    NYMDe = GET_NYMDe()
    CALL YMD_EXTRACT( NYMDe, Y, M, D )
    JD = JULDAY( Y, M, DBLE( D ) )

    ! Astronomical Julian Day corresponding to the 1st of the year
    JD0 = JULDAY( Y, 1, 0d0 )

    ! LASTDAY is the day of year corresponding to NYMDe
    LASTDAY = JD - JD0

    ! Skip past the element of NJDAY for Feb 29, if necessary
    IF ( .not. ITS_A_LEAPYEAR( Y, .TRUE. ) .and. LASTDAY > 59 ) THEN
       LASTDAY = LASTDAY + 1
    ENDIF

    ! Exit w/ error if THIS_NJDAY = 0
    IF ( Input_Opt%NJDAY(LASTDAY) == 0 ) THEN
       errMsg = 'No output scheduled on last day of run!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

  END SUBROUTINE IS_LAST_DAY_GOOD
!EOC
#endif
#endif
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: do_error_checks
!
! !DESCRIPTION: Makes sure that certain species are defined in order to
!  proceed with a certain option.  Halts the simulation with an error message
!  if incorrect inputs  would have caused  a simulation to crash.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Do_Error_Checks( Input_Opt, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod, ONLY : OptInput
    USE State_Chm_Mod, ONLY : Ind_
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(INOUT) :: Input_Opt
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC
!
! !REMARKS:
!  These error checks were originally called when the various menus were
!  read in from disk.  However, in order to use the Ind_() function to look
!  up given species indices, we need to call this after the Species Database
!  (which is in State_Chm) is initialized.  Therefore, we have now moved these
!  error checks to this routine, which is now called from GC_Init_Extra.
!                                                                             .
!  The Ind_() function now defines all species ID's.  It returns -1 if
!  a species cannot be found.  The prior behavior was to return 0 if a
!  species wasn't found.  Therefore, in order to preserve the logic of the
!  error checks, we must force any -1's returned by Ind_() to 0's in
!  this subroutine.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER            :: I

    ! Strings
    CHARACTER(LEN=255) :: MSG, LOCATION

    !=================================================================
    ! Initialization
    !=================================================================

    ! Assume success
    RC       = GC_SUCCESS

    ! Define location string
    LOCATION = '-> at Do_Error_Checks (in GeosCore/input_mod.F90)'

    !=================================================================
    ! Error check SEASALT AEROSOLS
    !=================================================================
    I = MAX( Ind_('SALA','A'), 0 ) + MAX( Ind_('SALC','A'), 0 )

    IF ( Input_Opt%LSSALT ) THEN
       IF ( I == 0 ) THEN
          MSG = 'LSSALT=T but ONLINE SEASALT AEROSOLS are undefined!'
          CALL GC_Error( Msg, RC, Location )
          RETURN
       ENDIF
    ELSE
       IF ( I > 0 ) THEN
          MSG = 'Cannot use ONLINE SEASALT AEROSOLS if LSSALT=F!'
          CALL GC_Error( Msg, RC, Location )
          RETURN
       ENDIF
    ENDIF

    !=================================================================
    ! Error check MARINE ORGANIC AEROSOLS
    !=================================================================
    I = MAX( Ind_('MOPO','A'), 0 ) + MAX( Ind_('MOPI','A'), 0 )

    IF ( Input_Opt%LMPOA ) THEN
       IF ( .not. Input_Opt%LSSALT ) THEN
          MSG = 'LMPOA=T but LSSALT=F!'
          CALL GC_Error( Msg, RC, Location )
          RETURN
       ENDIF
       IF ( I == 0 ) THEN
          MSG = 'LMPOA=T but MARINE ORGANIC AEROSOLS are undefined!'
          CALL GC_Error( Msg, RC, Location )
          RETURN
       ENDIF
    ELSE
       IF ( I > 0 ) THEN
          MSG = 'Cannot use MARINE ORGANIC AEROSOLS if LMPOA=F!'
          CALL GC_Error( Msg, RC, Location )
          RETURN
       ENDIF
    ENDIF

    !=================================================================
    ! Error check SULFUR AEROSOLS
    !=================================================================
    I = MAX( Ind_('DMS' ,'A'), 0 ) + &
        MAX( Ind_('SO2' ,'A'), 0 ) + &
        MAX( Ind_('SO4' ,'A'), 0 ) + &
        MAX( Ind_('SO4s','A'), 0 ) + &
        MAX( Ind_('HMS' ,'A'), 0 ) + &! (jmm, 07/2/18)
        MAX( Ind_('MSA' ,'A'), 0 ) + &
        MAX( Ind_('NH3' ,'A'), 0 ) + &
        MAX( Ind_('NH4' ,'A'), 0 ) + &
        MAX( Ind_('NITs','A'), 0 )

    IF ( Input_Opt%LSULF ) THEN

       ! We now compute the production of SO4s and NITs, so when
       ! LSULF=T, then we must also have LSSALT=T (bec, bmy, 4/13/05)
       IF ( .not. Input_Opt%LSSALT ) THEN
          MSG = 'LSULF=T now also requires LSSALT=T!'
          CALL GC_Error( Msg, RC, Location )
          RETURN
       ENDIF

       ! Stop w/ error if everything is undefined
       IF ( I == 0 ) THEN
          MSG = 'LSULF=T but ONLINE SULFUR AEROSOLS are undefined!'
          CALL GC_Error( Msg, RC, Location )
          RETURN
       ENDIF

    ELSE

       ! If LSULF=F but we have defined species, stop w/ error
       IF ( I > 0 ) THEN
          MSG = 'Cannot use ONLINE SULFUR AEROSOLS if LSULF=F!'
          CALL GC_Error( Msg, RC, Location )
          RETURN
       ENDIF

    ENDIF

    !=================================================================
    ! Error check CARBON AEROSOLS
    !=================================================================

    ! SOAupdate: Add POA (hotp 10/11/09)
    I = MAX( Ind_('BCPO','A'), 0 ) + &
        MAX( Ind_('BCPI','A'), 0 ) + &
        MAX( Ind_('OCPO','A'), 0 ) + &
        MAX( Ind_('OCPI','A'), 0 ) + &
        MAX( Ind_('POA1','A'), 0 )

    IF ( Input_Opt%LCARB ) THEN
       IF ( I == 0 ) THEN
          MSG = 'LCARB=T but ONLINE CARBON AEROSOLS are undefined!'
          CALL GC_Error( Msg, RC, Location )
          RETURN
       ENDIF
    ELSE
       IF ( I > 0 ) THEN
          MSG = 'Cannot use ONLINE CARBON AEROSOLS if LCARB=F!'
          CALL GC_Error( Msg, RC, Location )
          RETURN
       ENDIF
    ENDIF

    IF ( Input_Opt%LSVPOA .and. ( .NOT. Input_Opt%LSOA ) ) THEN
       MSG = 'Semivolatile POA requires COMPLEX SOA (LSOA=T)'
       CALL GC_Error( Msg, RC, Location )
       RETURN
    ENDIF

    ! SOAupdate: Error check (hotp 8/24/09)
    ! OCPI and OCPO are the non-volatile POA species
    ! POA (along w/ POG, OPOA, and OPOG) are the semivol POA species
    ! You can't have both!
    I = MAX( Ind_('OCPI','A'), 0 ) + MAX( Ind_('OCPO','A'), 0 )

    IF ( Ind_('POA1') > 0 ) THEN
       IF ( I > 0 ) THEN
          MSG = 'Semivolatile POA species is defined in addition to ' // &
                 'Nonvolatile POA'
          CALL GC_Error( Msg, RC, Location )
          RETURN
       ENDIF
       IF ( ( .NOT. Input_Opt%LSOA   ) .or. &
            ( .NOT. Input_Opt%LSVPOA ) ) THEN
          MSG = 'Semivolatile POA requires LSOA=T and LSVPOA=T'
          CALL GC_Error( Msg, RC, Location )
          RETURN
       ENDIF
    ENDIF

    ! SOAupdate
    ! Options for organic aerosol species:
    ! IF LSOA = F: only OCPI and OCPO
    ! IF LSOA = T:
    !   OCPI OCPO SOA (non-vol + original traditional)
    !   POA POG OPOA OPOG SOA BTX NAP (semivol + orig trad + IVOC )
    ! NAP emissions are set in input.geos
    ! LSVPOA is just a check (doesn't do anything hotp 7/21/10)
    I = MAX( Ind_('POA1' ,'A'), 0 ) + &
        MAX( Ind_('POA2' ,'A'), 0 ) + &
        MAX( Ind_('POG1' ,'A'), 0 ) + &
        MAX( Ind_('POG2' ,'A'), 0 ) + &
        MAX( Ind_('OPOA1','A'), 0 ) + &
        MAX( Ind_('OPOA2','A'), 0 ) + &
        MAX( Ind_('OPOG1','A'), 0 ) + &
        MAX( Ind_('OPOG2','A'), 0 )

    IF ( Input_Opt%LSVPOA ) THEN
       IF ( I < 8 ) THEN
          MSG = 'Not enough semivolatile POA species!'
          CALL GC_Error( Msg, RC, Location )
          RETURN
       ENDIF
       IF ( Ind_('NAP','A') < 0 ) THEN
          MSG = 'Semivolatile POA requires IVOCs/NAP!'
          CALL GC_Error( Msg, RC, Location )
          RETURN
       ENDIF
    ENDIF

    !=================================================================
    ! Error check SECONDARY ORGANIC AEROSOLS
    !=================================================================

    ! Check for complex SOA species
    I = MAX( Ind_('TSOA1','A'), 0 ) + &
        MAX( Ind_('TSOA2','A'), 0 ) + &
        MAX( Ind_('TSOA3','A'), 0 ) + &
        MAX( Ind_('ASOA1','A'), 0 ) + &
        MAX( Ind_('ASOA2','A'), 0 ) + &
        MAX( Ind_('ASOA3','A'), 0 ) + &
        MAX( Ind_('ASOAN','A'), 0 )

    IF ( Input_Opt%LSOA ) THEN
       IF ( I == 0 ) THEN
          MSG = 'LSOA=T but COMPLEX SOA species are undefined!'
          CALL GC_Error( Msg, RC, Location )
          RETURN
       ENDIF
    ELSE
       IF ( I > 0 ) THEN
          MSG = 'Cannot use COMPLEX SOA species if LSOA=F!'
          CALL GC_Error( Msg, RC, Location )
          RETURN
       ENDIF
    ENDIF

    !=================================================================
    ! Error check DUST AEROSOLS
    !=================================================================

#ifdef TOMAS
    ! For TOMAS only: If DUST1 is present, the other dust species are too
    I = MAX( Ind_('DUST1','A'), 0 )
#else
    ! Non-TOMAS simulations: Need all DST1-DST4 species
    I = MAX( Ind_('DST1','A'), 0 ) + &
        MAX( Ind_('DST2','A'), 0 ) + &
        MAX( Ind_('DST3','A'), 0 ) + &
        MAX( Ind_('DST4','A'), 0 )
#endif

    IF ( Input_Opt%LDUST ) THEN
       IF ( I == 0 ) THEN
          MSG = 'LDUST=T but ONLINE DUST AEROSOLS are undefined!'
          CALL GC_Error( Msg, RC, Location )
          RETURN
       ENDIF
    ELSE
       IF ( I > 0 ) THEN
          MSG = 'Cannot use ONLINE DUST AEROSOLS if LDUST=F!'
          CALL GC_Error( Msg, RC, Location )
          RETURN
       ENDIF
    ENDIF

    !=================================================================
    ! Error check DUST NITRATE    AEROSOLS
    !             DUST SULFATE    AEROSOLS
    !             DUST ALKALINITY AEROSOLS
    !=================================================================
    I = MAX( Ind_('NITd1'  ,'A'), 0 ) + &
        MAX( Ind_('NITd2'  ,'A'), 0 ) + &
        MAX( Ind_('NITd3'  ,'A'), 0 ) + &
        MAX( Ind_('NITd4'  ,'A'), 0 ) + &
        MAX( Ind_('SO4d1'  ,'A'), 0 ) + &
        MAX( Ind_('SO4d2'  ,'A'), 0 ) + &
        MAX( Ind_('SO4d3'  ,'A'), 0 ) + &
        MAX( Ind_('SO4d4'  ,'A'), 0 ) + &
        MAX( Ind_('DSTAL1' ,'A'), 0 ) + &
        MAX( Ind_('DSTAL2' ,'A'), 0 ) + &
        MAX( Ind_('DSTAL3' ,'A'), 0 ) + &
        MAX( Ind_('DSTAL4' ,'A'), 0 )

    IF ( Input_Opt%LDSTUP ) THEN
       IF ( I < 12 ) THEN
          MSG = 'LDSTUP=T but COATED DUST AEROSOLS are undefined!'
          CALL GC_Error( Msg, RC, Location )
          RETURN
       ENDIF
    ELSE
       IF ( I > 0 ) THEN
          MSG = 'Cannot use COATED DUST AEROSOLS if LDSTUP=F!'
          CALL GC_Error( Msg, RC, Location )
          RETURN
       ENDIF
    ENDIF

    !=================================================================
    ! Error check SEASALT AEROSOLS
    !=================================================================
    I = MAX( Ind_('SALA','A'), 0 ) + MAX( Ind_('SALC','A'), 0 )

    IF ( Input_Opt%LSSALT ) THEN
       IF ( I == 0 ) THEN
          MSG = 'LSSALT=T but ONLINE SEASALT AEROSOLS are undefined!'
          CALL GC_Error( Msg, RC, Location )
          RETURN
       ENDIF
    ELSE
       IF ( I > 0 ) THEN
          MSG = 'Cannot use ONLINE SEASALT AEROSOLS if LSSALT=F!'
          CALL GC_Error( Msg, RC, Location )
          RETURN
       ENDIF
    ENDIF

    !=================================================================
    ! Error check stratospheric H2O
    !=================================================================
    IF ( Input_Opt%LSETH2O .and. Ind_('H2O') < 0 ) THEN
       WRITE( 6, '(a)'     ) REPEAT( '=', 79 )
       WRITE( 6, '(/,a,/)' ) 'Warning in input_mod.F90: ' &
            // 'H2O is set but H2O species is undefined.'
       Input_Opt%LSETH2O = .FALSE.
       WRITE( 6, '(a)'     ) REPEAT( '=', 79 )
    ENDIF

  END SUBROUTINE Do_Error_Checks
!EOC
END MODULE Input_Mod
