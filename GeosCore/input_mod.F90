!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: input_mod.F90
!
! !DESCRIPTION: Contains routines that read the GEOS-Chem configuration file at
!  the start of the run and pass the information to Input_Opt.
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
! !DEFINED PARAMETERS:
!
  ! YAML configuration file name to be read
  CHARACTER(LEN=21), PARAMETER, PRIVATE :: configFile ='./geoschem_config.yml'

CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_input_file
!
! !DESCRIPTION: Driver program for reading the GEOS-Chem input file.
!\\
!\\
! In an ESMF environment, all time steps (chemistry, convection, emissions,
! dynamics) must be specified externally before calling this routine.
! The time steps specified in the GEOS-Chem configuration file are ignored.
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
       errMsg = 'Error in "Config_Timesteps"!'
       CALL GC_Error( errMsg, RC, thisLoc  )
       CALL QFYAML_CleanUp( Config         )
       CALL QFYAML_CleanUp( ConfigAnchored )
       RETURN
    ENDIF

    !========================================================================
    ! Get settings for GEOS-Chem operations from the YAML Config object
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

    ! Passive Species settigns (if any)
    CALL Config_PassiveSpecies( Config, Input_Opt, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error in "Config_PassiveSpecies"!'
       CALL GC_Error( errMsg, RC, thisLoc )
       CALL QFYAML_CleanUp( Config         )
       CALL QFYAML_CleanUp( ConfigAnchored )
       RETURN
    ENDIF

    ! Convection and PBL mixing settings
    CALL Config_Convection_Mixing( Config, Input_Opt, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error in "Config_Convection_Mixing"!'
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
    ! Get settings for specialty simulations from the YAML Config object
    !========================================================================

    ! CH4/carbon simulation settings
    IF ( Input_Opt%Its_A_CH4_Sim .or. Input_Opt%Its_A_Carbon_Sim ) THEN
       CALL Config_CH4( Config, Input_Opt, RC )
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = 'Error in "Config_CH4"!'
          CALL GC_Error( errMsg, RC, thisLoc  )
          CALL QFYAML_CleanUp( Config         )
          CALL QFYAML_CleanUp( ConfigAnchored )
          RETURN
       ENDIF
    ENDIF

    ! CO simulation settings
    IF ( Input_Opt%Its_A_TagCO_Sim .or. Input_Opt%Its_A_Carbon_Sim ) THEN
       CALL Config_CO( Config, Input_Opt, RC )
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = 'Error in "Config_CO"!'
          CALL GC_Error( errMsg, RC, thisLoc  )
          CALL QFYAML_CleanUp( Config         )
          CALL QFYAML_CleanUp( ConfigAnchored )
          RETURN
       ENDIF
    ENDIF


    ! CO2/carbon simulation settings
    IF ( Input_Opt%Its_A_CO2_Sim .or. Input_Opt%Its_A_Carbon_Sim ) THEN
       CALL Config_CO2( Config, Input_Opt, RC )
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = 'Error in "Config_CO2"!'
          CALL GC_Error( errMsg, RC, thisLoc  )
          CALL QFYAML_CleanUp( Config         )
          CALL QFYAML_CleanUp( ConfigAnchored )
          RETURN
       ENDIF
    ENDIF

    ! Hg simulation settings
    IF ( Input_Opt%Its_A_Mercury_Sim ) THEN
       CALL Config_Hg( Config, Input_Opt, RC )
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = 'Error in "Config_Hg"!'
          CALL GC_Error( errMsg, RC, thisLoc  )
          CALL QFYAML_CleanUp( Config         )
          CALL QFYAML_CleanUp( ConfigAnchored )
          RETURN
       ENDIF
    ENDIF

    ! POPs simulation settings
    IF ( Input_Opt%Its_A_POPs_Sim ) THEN
       CALL Config_POPs( Config, Input_Opt, RC )
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = 'Error in "Config_POPs"!'
          CALL GC_Error( errMsg, RC, thisLoc  )
          CALL QFYAML_CleanUp( Config         )
          CALL QFYAML_CleanUp( ConfigAnchored )
          RETURN
       ENDIF
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
    ! GAMAP metadata files
    ! (Skip if we are connecting to an external model)
    CALL Config_Gamap( Config, Input_Opt, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error in "Config_Gamap"!'
       CALL GC_Error( errMsg, RC, thisLoc  )
       CALL QFYAML_CleanUp( Config         )
       CALL QFYAML_CleanUp( ConfigAnchored )
       RETURN
    ENDIF

    ! ND51 timeseries
    ! (Skip if we are connecting to an external model)
    CALL Config_ND51( Config, Input_Opt, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error in "Config_ND51"!'
       CALL GC_Error( errMsg, RC, thisLoc  )
       CALL QFYAML_CleanUp( Config         )
       CALL QFYAML_CleanUp( ConfigAnchored )
       RETURN
    ENDIF

    ! ND51b timeseries
    ! (Skip if we are connecting to an external model)
    CALL Config_ND51b( Config, Input_Opt, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error in "Config_ND51b"!'
       CALL GC_Error( errMsg, RC, thisLoc  )
       CALL QFYAML_CleanUp( Config         )
       CALL QFYAML_CleanUp( ConfigAnchored )
       RETURN
    ENDIF

    ! Bpch diagnostic output
    ! (Skip if we are connecting to an external model)
    CALL Config_Bpch_Output( Config, Input_Opt, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error in "Config_Bpch_Output"!'
       CALL GC_Error( errMsg, RC, thisLoc  )
       CALL QFYAML_CleanUp( Config         )
       CALL QFYAML_CleanUp( ConfigAnchored )
       RETURN
    ENDIF

#ifdef TOMAS
    ! Bpch prod & loss -- only needed for TOMAS
    ! (Skip if we are connecting to an external model)
    CALL Config_Bpch_ProdLoss( Config, Input_Opt, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error in "Config_Bpch_ProdLoss"!'
       CALL GC_Error( errMsg, RC, thisLoc  )
       CALL QFYAML_CleanUp( Config         )
       CALL QFYAML_CleanUp( ConfigAnchored )
       RETURN
    ENDIF
#endif
#endif

    !========================================================================
    ! Check GEOS-CHEM timesteps
    ! NOTE: Skip for GCHP/GEOS, as this is called from GCHP_Chunk_Run
    !========================================================================
    CALL Check_Time_Steps( Input_Opt, State_Grid, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error in "Check_Time_Steps"!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
#endif

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
    ! Config_Simulation begins here!
    !========================================================================

    ! Initialize
    RC      = GC_SUCCESS
    errMsg  = ''
    thisLoc = &
     ' -> at Config_Simulation (in module GeosCore/input_mod.F90)'

    !------------------------------------------------------------------------
    ! Simulation type
    !------------------------------------------------------------------------
    key   = "simulation%name"
    v_str = MISSING_STR
    CALL QFYAML_Add_Get( Config, TRIM( key ), v_str, "", RC )
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
         TRIM(Sim) /= 'TAGHG'   .and. TRIM(Sim) /= 'TAGO3'             .and. &
         TRIM(Sim) /= 'CARBON'                                       ) THEN
         
       errMsg = Trim( Input_Opt%SimulationName) // ' is not a'            // &
                ' valid simulation. Supported simulations are:'           // &
                ' aerosol, carbon, CH4, CO2, fullchem, Hg, POPs,'    // &
                ' TransportTracers, TagCO, TagCH4, or TagO3.'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    ! Set simulation type flags in Input_Opt
    Input_Opt%ITS_A_CH4_SIM        = ( TRIM(Sim) == 'CH4'              .or.  &
                                       TRIM(Sim) == 'TAGCH4'                )
    Input_Opt%ITS_A_CO2_SIM        = ( TRIM(Sim) == 'CO2'                   )
    Input_Opt%ITS_A_FULLCHEM_SIM   = ( TRIM(Sim) == 'FULLCHEM'              )
    Input_Opt%ITS_A_MERCURY_SIM    = ( TRIM(Sim) == 'HG'                    )
    Input_Opt%ITS_A_POPS_SIM       = ( TRIM(Sim) == 'POPS'                  )
    Input_Opt%ITS_A_RnPbBe_SIM     = ( TRIM(Sim) == 'TRANSPORTTRACERS'      )
    Input_Opt%ITS_A_TAGO3_SIM      = ( TRIM(Sim) == 'TAGO3'                 )
    Input_Opt%ITS_A_TAGCO_SIM      = ( TRIM(Sim) == 'TAGCO'                 )
    Input_Opt%ITS_AN_AEROSOL_SIM   = ( TRIM(Sim) == 'AEROSOL'               )
    Input_Opt%ITS_A_TRACEMETAL_SIM = ( TRIM(SIM) == 'METALS'                )
    Input_Opt%ITS_A_CARBON_SIM     = ( TRIM(Sim) == 'CARBON'                )

    !------------------------------------------------------------------------
    ! Species database file
    !------------------------------------------------------------------------
    key   = "simulation%species_database_file"
    v_str = MISSING_STR
    CALL QFYAML_Add_Get( Config, TRIM( key ), v_str, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Input_Opt%SpcDataBaseFile = TRIM( v_str )

    !------------------------------------------------------------------------
    ! Species metadata output file
    !------------------------------------------------------------------------
    key   = "simulation%species_metadata_output_file"
    v_str = MISSING_STR
    CALL QFYAML_Add_Get( Config, TRIM( key ), v_str, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Input_Opt%SpcMetaDataOutFile = TRIM( v_str )

    !------------------------------------------------------------------------
    ! Turn on debug output
    !------------------------------------------------------------------------
    key    = "simulation%debug_printout"
    v_bool = MISSING_BOOL
    CALL QFYAML_Add_Get( Config, TRIM( key ), v_bool, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Input_Opt%LPRT = v_bool

#if defined( MODEL_GCHP ) || defined( MODEL_GEOS )
    !========================================================================
    !          %%%%%%% GCHP and NASA/GEOS (with ESMF & MPI) %%%%%%%
    !
    ! Because GCHP and NASA/GEOS use ESMF, we need to take the start & end
    ! dates from the ESMF resource file (GEOSCHEMchem_GridComp_mod.rc)
    ! instead of those in geoschem_config.yml.  Therefore, the following
    ! fields have already been defined in the GCHP_Chunk_Init routine
    ! (located in module Interfaces/GCHP/Chem_GridCompMod.F90):
    !
    ! (1) Input_Opt%NYMDb
    ! (2) Input_Opt%NHMSb
    ! (3) Input_Opt%NYMDe
    ! (4) Input_Opt%NYMDe
    !
    ! Pass these fields of Input_Opt to routine Accept_External_Date_Time,
    ! as well as the starting UTC value in hours.
    !========================================================================

    ! Make sure the starting date NYMDb is valid
    IF ( .not. Valid_Date( Input_Opt%NYMDb ) ) THEN
       WRITE( DateStr, '(i8.8)' ) Input_Opt%NYMDb
       errMsg = 'Input%Opt%NYMDb = ' // DateStr // ' is not a valid '     // &
                'calendar date!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    ! Make sure the starting time NHMSb is valid
    IF ( .not. Valid_Time( Input_Opt%NHMSb ) ) THEN
       WRITE( TimeStr, '(i6.6)' ) Input_Opt%NHMSb
       errMsg = 'Input%Opt%NHMSb = ' // TimeStr // ' is not a valid '     // &
                'clock time!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    ! Make sure the ending date NYMDe is valid
    IF ( .not. Valid_Date( Input_Opt%NYMDe ) ) THEN
       WRITE( DateStr, '(i8.8)' ) Input_Opt%NYMDe
       errMsg = 'Input%Opt%NYMDe = ' // DateStr // ' is not a valid '     // &
                'calendar date!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    ! Make sure the ending time NHMSe is valid
    IF ( .not. Valid_Time( Input_Opt%NHMSe ) ) THEN
       WRITE( TimeStr, '(i6.6)' ) Input_Opt%NHMSe
       errMsg = 'Input%Opt%NHMSe = ' // TimeStr // ' is not a valid '     // &
                'clock time!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    ! Get the starting UTC time from Input_Opt%NHMSb for use below
    CALL YMD_Extract( Input_Opt%NHMSb, H, M, S )
    init_UTC = ( H + ( M / 60 ) + ( S / 3600 ) )

    ! Pass the values for the start & end times of the simulation directly
    ! to GeosUtil/time_mod.F90 via subroutine ACCEPT_EXTERNAL_DATE_TIME.
    ! (bmy, 12/6/12)
    CALL Accept_External_Date_Time( value_NYMDb = Input_Opt%NYMDb,           &
                                    value_NHMSb = Input_Opt%NHMSb,           &
                                    value_NYMDe = Input_Opt%NYMDe,           &
                                    value_NHMSe = Input_Opt%NHMSe,           &
                                    value_NYMD  = Input_Opt%NYMDb,           &
                                    value_NHMS  = Input_Opt%NHMSb,           &
                                    value_UTC   = init_UTC,                  &
                                    RC          = RC                        )

    !------------------------------------------------------------------------
    ! Chemistry inputs folder (GCHP/NASA-GEOS only)
    !------------------------------------------------------------------------

    key   = "simulation%chem_inputs_dir"
    v_str = MISSING_STR
    CALL QFYAML_Add_Get( Config, TRIM( key ), v_str, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Input_Opt%CHEM_INPUTS_DIR = TRIM( v_str )

    !------------------------------------------------------------------------
    ! Set other fields of Input_Opt accordingly
    !------------------------------------------------------------------------
    Input_Opt%MetField        = 'See ExtData.rc'
    Input_Opt%DATA_DIR        = 'N/A'
    Input_Opt%RUN_DIR         = 'N/A'

#else
    !========================================================================
    !             %%%%%%% GEOS-Chem CLASSIC (with OpenMP) %%%%%%%
    !
    ! If we aren't using ESMF, read extra settings from geoschem_config.yml
    !========================================================================

    !------------------------------------------------------------------------
    ! Simulation start date
    !------------------------------------------------------------------------
    key   = "simulation%start_date"
    a_int = MISSING_INT
    CALL QFYAML_Add_Get( Config, TRIM( key ), a_int, "", RC )
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
                ' Please check your "geoschem_config.yml" file.'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    ! Make sure the starting time NHMSb is valid
    IF ( .not. Valid_Time( Input_Opt%NHMSb ) ) THEN
       WRITE( TimeStr, '(i6.6)' ) Input_Opt%NHMSb
       errMsg = 'Input%Opt%NHMSb = ' // TimeStr        // &
                ' is not a valid clock time!'          // &
                ' Please check your "geoschem_config.yml" file.'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! Simulation end date
    !------------------------------------------------------------------------
    key   = "simulation%end_date"
    a_int = MISSING_INT
    CALL QFYAML_Add_Get( Config, TRIM( key ), a_int, "", RC )
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
                ' Please check your "geoschem_config.yml" file.'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    ! Make sure the ending time NHMSe is valid
    IF ( .not. Valid_Time( Input_Opt%NHMSe ) ) THEN
       WRITE( TimeStr, '(i6.6)' ) Input_Opt%NHMSe
       errMsg = 'Input%Opt%NHMSe = ' // TimeStr        // &
                ' is not a valid clock time!'          // &
                ' Please check your "geoschem_config.yml" file.'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! Root data directory
    !------------------------------------------------------------------------
    key   = "simulation%root_data_dir"
    v_str = MISSING_STR
    CALL QFYAML_Add_Get( Config, TRIM( key ), v_str, "", RC )
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
    CALL QFYAML_Add_Get( Config, TRIM( key ), v_str, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Input_Opt%MetField = TRIM( v_str )

#if !defined( MODEL_CESM )
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
          errMsg = Trim( Input_Opt%MetField ) // ' is not a valid '       // &
                ' met field. Supported met fields are GEOS-FP, '          // &
                ' MERRA-2 and ModelE2.1. Please check your '              // &
                '"geoschem_config.ymls" file.'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    END SELECT
#endif

    !------------------------------------------------------------------------
    ! Turn on timers
    !------------------------------------------------------------------------
    key    = "simulation%use_gcclassic_timers"
    v_bool = MISSING_BOOL
    CALL QFYAML_Add_Get( Config, TRIM( key ), v_bool, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Input_Opt%UseTimers = v_bool

    ! Error check simulation name
    Sim = To_UpperCase( TRIM( Input_Opt%SimulationName ) )
    IF ( TRIM(Sim) /= 'AEROSOL' .and. TRIM(Sim) /= 'CH4'               .and. &
         TRIM(Sim) /= 'CO2'     .and. TRIM(Sim) /= 'FULLCHEM'          .and. &
         TRIM(Sim) /= 'HG'      .and. TRIM(Sim) /= 'METALS'            .and. &
         TRIM(Sim) /= 'POPS'    .and. TRIM(Sim) /= 'TRANSPORTTRACERS'  .and. &
         TRIM(Sim) /= 'TAGCO'   .and. TRIM(Sim) /= 'TAGCH4'            .and. &
         TRIM(Sim) /= 'CARBON'                                         .and. &
         TRIM(Sim) /= 'TAGHG'   .and. TRIM(Sim) /= 'TAGO3'           ) THEN
       ErrMsg = Trim( Input_Opt%SimulationName) // ' is not a'            // &
                ' valid simulation. Supported simulations are:'           // &
                ' aerosol, CH4, CO2, fullchem, Hg, POPs,'                 // &
                ' TransportTracers, TagCO, TagCH4, TagHg, or TagO3.'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
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
    Input_Opt%ITS_A_CARBON_SIM     = ( TRIM(SIM) == 'CARBON'                )

    !------------------------------------------------------------------------
    ! Set start time of run in "time_mod.F90"
    !------------------------------------------------------------------------
    CALL Set_Begin_Time( Input_Opt%NYMDb, Input_Opt%NHMSb, RC  )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error encountered in "Set_Begin_Time"!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! Set end time of run in "time_mod.F90"
    !------------------------------------------------------------------------
    errMsg = 'Error encountered in "Set_Begin_Time"!'
    CALL Set_End_Time( Input_Opt%NYMDe, Input_Opt%NHMSe, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    ! Set the current time
    CALL Set_Current_Time()

#ifdef BPCH_DIAG
    !------------------------------------------------------------------------
    ! Set the start of the 1st diagnostic interval
    ! Only applies when GEOS-Chem is configured with -DBPCH_DIAG=y
    !------------------------------------------------------------------------
    CALL Set_DiagB( GET_TAU() )
#endif

#endif

    !========================================================================
    ! Compute the length of the simulation, in elapsed seconds
    !========================================================================
    JulianDateStart        = GET_JD( Input_Opt%NymdB, Input_Opt%NhmsB )
    JulianDateEnd          = GET_JD( Input_Opt%NymdE, Input_Opt%NhmsE )
    Input_Opt%SimLengthSec = NINT( ( JulianDateEnd - JulianDateStart  )      &
                           * 86400_f8)

    ! Return success
    RC = GC_SUCCESS

    !========================================================================
    ! Print to screen
    !========================================================================
    IF ( Input_Opt%amIRoot ) THEN
       WRITE( 6, 90  ) 'SIMULATION SETTINGS'
       WRITE( 6, 95  ) '-------------------'
       WRITE( 6, 110 ) 'Simulation name             : ',                     &
                        TRIM( Input_Opt%SimulationName )
       WRITE( 6, 110 ) 'CHEM_INPUTS directory       : ',                     &
                        TRIM( Input_Opt%CHEM_INPUTS_DIR )
       WRITE( 6, 110 ) 'Species database file       : ',                     &
                        TRIM( Input_Opt%SpcDatabaseFile )
       WRITE( 6, 120 ) 'Turn on debug output        : ',                     &
                        Input_Opt%LPRT
#ifdef MODEL_CLASSIC
       WRITE( 6, 100 ) 'Start time of run           : ',                     &
                        Input_Opt%NYMDb, Input_Opt%NHMSb
       WRITE( 6, 100 ) 'End time of run             : ',                     &
                        Input_Opt%NYMDe, Input_Opt%NHMSe
       WRITE( 6, 110 ) 'Data Directory              : ',                     &
                        TRIM( Input_Opt%DATA_DIR )
       WRITE( 6, 110 ) 'Meteorology field           : ',                     &
                        TRIM( Input_Opt%MetField )
       WRITE( 6, 120 ) 'Turn on GEOS-Chem timers    : ',                     &
                        Input_Opt%useTimers
#endif
    ENDIF

    ! Format statements
90  FORMAT( /, A              )
95  FORMAT( A                 )
100 FORMAT( A, I8.8, 1X, I6.6 )
110 FORMAT( A, A              )
120 FORMAT( A, L5             )

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
    USE RoundOff_Mod,   ONLY : Cast_and_RoundOff
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
    thisLoc = ' -> at Config_Grid (in GeosCore/input_mod.F90)'

    !------------------------------------------------------------------------
    ! Grid resolution
    !------------------------------------------------------------------------
    key   = "grid%resolution"
    v_str = MISSING_STR
    CALL QFYAML_Add_Get( Config, TRIM( key ), v_str, "", RC )
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
                ' an x. Please check your "geoschem_config.yml" file.'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    ! Save the delta X and Y values
    State_Grid%DY = Cast_and_RoundOff( subStrs(1), places=4 )
    State_Grid%DX = Cast_and_RoundOff( subStrs(2), places=4 )

    !------------------------------------------------------------------------
    ! Level range
    !------------------------------------------------------------------------
    key   = "grid%number_of_levels"
    v_int = MISSING_INT
    CALL QFYAML_Add_Get( Config, TRIM( key ), v_int, "", RC )
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
    CALL QFYAML_Add_Get( Config, TRIM( key ), a_str, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    State_Grid%XMin = Cast_and_RoundOff( a_str(1), places=4 )
    State_Grid%XMax = Cast_and_RoundOff( a_str(2), places=4 )

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
    CALL QFYAML_Add_Get( Config, TRIM( key ), v_bool, "", RC )
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
    CALL QFYAML_Add_Get( Config, TRIM( key ), a_str, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    State_Grid%YMin = Cast_and_RoundOff( a_str(1), places=4 )
    State_Grid%YMax = Cast_and_RoundOff( a_str(2), places=4 )

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
    CALL QFYAML_Add_Get( Config, TRIM( key ), v_bool, "", RC )
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
    CALL QFYAML_Add_Get( Config, TRIM( key ), v_bool, "", RC )
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
    key   = "grid%nested_grid_simulation%buffer_zone_NSEW"
    a_int = MISSING_INT
    CALL QFYAML_Add_Get( Config, TRIM( key ), a_int, "", RC )
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
       WRITE( 6, 90  ) 'GRID SETTINGS'
       WRITE( 6, 95  )  '------------'
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
90  FORMAT( /, A                )
95  FORMAT( A                   )
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
    USE Input_Opt_Mod,  ONLY : OptInput
    USE State_Grid_Mod, ONLY : GrdState
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
    errMsg  = 'Error reading the "geoschem_config.yml" file!'
    thisLoc = ' -> at Config_Timestep (in module GeosCore/input_mod.F90)'

#if !(defined( EXTERNAL_GRID ) || defined( EXTERNAL_FORCING ))
    !------------------------------------------------------------------------
    ! Transport/convection timestep
    !------------------------------------------------------------------------
    key   = "timesteps%transport_timestep_in_s"
    v_int = MISSING_INT
    CALL QFYAML_Add_Get( Config, TRIM( key ), v_int, "", RC )
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
    CALL QFYAML_Add_Get( Config, TRIM( key ), v_int, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Input_Opt%TS_CHEM = v_int
    Input_Opt%TS_EMIS = v_int

    !------------------------------------------------------------------------
    ! RRTMG radiation timestep
    !------------------------------------------------------------------------
    key   = "timesteps%radiation_timestep_in_s"
    v_int = 0
    CALL QFYAML_Add_Get( Config, TRIM( key ), v_int, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Input_Opt%TS_RAD = v_int
#endif

    !========================================================================
    ! Error checks
    !========================================================================

#ifdef MODEL_CLASSIC
    IF ( Input_Opt%SimLengthSec < Input_Opt%TS_DYN                      .or. &
         Input_Opt%SimLengthSec < Input_Opt%TS_CHEM )                   THEN
       IF ( Input_Opt%amIRoot )                                         THEN
          WRITE( 6,'(a)' ) ''
          WRITE( 6,'(a)' ) 'The length of the simulation is shorter '
          WRITE( 6,'(a)' ) 'than the transport and/or chemistry '
          WRITE( 6,'(a)' ) 'timesteps. Check the settings in '
          WRITE( 6,'(a)' ) 'the "geoschem_config.yml" file.'
          WRITE( 6,'(a)' ) ''
          WRITE( 6,100 ) 'Transport/Convection  [sec] : ', Input_Opt%TS_DYN
          WRITE( 6,100 ) 'Chemistry/Emissions   [sec] : ', Input_Opt%TS_CHEM
          WRITE( 6,100 ) 'Simulation duration   [sec] : ',                   &
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
#endif

    ! Return success
    RC = GC_SUCCESS

    !========================================================================
    ! Print to screen
    !========================================================================
    IF ( Input_Opt%amIRoot ) THEN
       WRITE( 6, 90  ) 'TIMESTEP SETTINGS'
       WRITE( 6, 95  ) '-----------------'
       WRITE( 6, 100 ) 'Transport/Convection [sec]  : ', Input_Opt%TS_DYN
       WRITE( 6, 100 ) 'Chemistry/Emissions  [sec]  : ', Input_Opt%TS_CHEM
       WRITE( 6, 100 ) 'RRTMG rad transfer   [sec]  : ', Input_Opt%TS_RAD
    ENDIF

    ! Format statements
90  FORMAT( /, A   )
95  FORMAT( A      )
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
! !DESCRIPTION: Copies transport information from the Config object
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
    INTEGER                      :: C

    ! Arrays
    INTEGER                      :: a_int(3)

    ! Strings
    CHARACTER(LEN=255)           :: thisLoc
    CHARACTER(LEN=512)           :: errMsg
    CHARACTER(LEN=QFYAML_NamLen) :: key

    ! String arrays
    CHARACTER(LEN=QFYAML_NamLen) :: a_str(QFYAML_MaxArr)

    !========================================================================
    ! Config_Transport begins here!
    !========================================================================

    ! Initialize
    RC      = GC_SUCCESS
    errMsg  = 'Error reading the "geoschem_config.yml" file!'
    thisLoc = ' -> at Config_Transport (in GeosCore/input_mod.F90)'

#ifdef MODEL_CLASSIC
    !========================================================================
    !             %%%%%%% GEOS-Chem CLASSIC (with OpenMP) %%%%%%%
    !========================================================================

    !------------------------------------------------------------------------
    ! Turn on TPCORE transport?
    !------------------------------------------------------------------------
    key    = "operations%transport%gcclassic_tpcore%activate"
    v_bool = MISSING_BOOL
    CALL QFYAML_Add_Get( Config, TRIM( key ), v_bool, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Input_Opt%LTRAN = v_bool

    !------------------------------------------------------------------------
    ! Fill negative values generated by TPCORE?
    !------------------------------------------------------------------------
    key    = "operations%transport%gcclassic_tpcore%fill_negative_values"
    v_bool = MISSING_BOOL
    CALL QFYAML_Add_Get( Config, TRIM( key ), v_bool, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Input_Opt%LFILL = v_bool

    !------------------------------------------------------------------------
    ! IORD, JORD, KORD settings for TPCORE
    !------------------------------------------------------------------------
    key   = "operations%transport%gcclassic_tpcore%iord_jord_kord"
    a_int = MISSING_INT
    CALL QFYAML_Add_Get( Config, TRIM( key ), a_int, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Input_Opt%TPCORE_IORD = a_int(1)
    Input_Opt%TPCORE_JORD = a_int(2)
    Input_Opt%TPCORE_KORD = a_int(3)

#else

    !========================================================================
    !             %%%%%%% OTHER INSTANCES OF GEOS-Chem %%%%%%%
    !========================================================================

    !------------------------------------------------------------------------
    ! Turn on transport?
    !------------------------------------------------------------------------
    key    = "operations%transport%activate"
    v_bool = MISSING_BOOL
    CALL QFYAML_Add_Get( Config, TRIM( key ), v_bool, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Input_Opt%LTRAN = v_bool

#endif

    !------------------------------------------------------------------------
    ! Transported (advected) species list
    !------------------------------------------------------------------------
    key   = "operations%transport%transported_species"
    a_str = MISSING_STR
    CALL QFYAML_Add_Get( Config, TRIM( key ), a_str,                         &
                         "",     RC,          dynamic_size=.TRUE.           )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    ! Find the number of valid species (not equal to MISSING_STR)
    Input_Opt%N_Advect = Get_Number_of_Species( a_str )

    ! Remove single quotes from species (i.e. 'NO' -> NO)
    DO N = 1, Input_Opt%N_Advect
       C = INDEX( a_str(N), "'" )
       IF ( C > 0 ) THEN
          a_str(N) = a_str(N)(C+1:)
          C = INDEX( a_str(N), "'" )
          IF ( C > 0 ) a_str(N) = a_str(N)(1:C-1)
       ENDIF
    ENDDO

    ! Throw an error if there are more species than the array length
    IF ( Input_Opt%N_Advect > Input_Opt%Max_AdvectSpc ) THEN
       errMsg = 'Number of advected species exceeds maximum. ' // &
            'This value can be modified in input_opt_mod.F90.'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    ! Copy the advected species names into Input_Opt%Advect_Spc_Name
    Input_Opt%AdvectSpc_Name(1:Input_Opt%N_Advect) = &
         a_str(1:Input_Opt%N_Advect)

    !========================================================================
    ! Print to screen
    !========================================================================
    IF ( Input_Opt%amIRoot ) THEN
       WRITE( 6,90  ) 'TRANSPORT SETTINGS'
       WRITE( 6,95  ) '------------------'
       WRITE( 6,100 ) 'Turn on transport?          : ',Input_Opt%LTRAN
#ifdef MODEL_CLASSIC
       WRITE( 6,100 ) 'Let TPCORE Fill negatives?  : ',Input_Opt%LFILL
       WRITE( 6,110 ) 'IORD, JORD, KORD for TPCORE?: ',Input_Opt%TPCORE_IORD,&
                                                       Input_Opt%TPCORE_JORD,&
                                                       Input_Opt%TPCORE_KORD
#endif
    ENDIF

    ! FORMAT statements
90  FORMAT( /, A   )
95  FORMAT( A      )
100 FORMAT( A, L5  )
110 FORMAT( A, 5I5 )

    !=================================================================
    ! Call setup routines from other F90 modules
    !=================================================================

    ! Split into tagged species (turn off for full-chemistry)
    IF ( Input_Opt%ITS_A_FULLCHEM_SIM ) THEN
       Input_Opt%LSPLIT = .FALSE.                     ! No tagged species
    ELSE IF ( Input_Opt%ITS_A_CARBON_SIM ) THEN
       Input_Opt%LSPLIT = ( Input_Opt%N_ADVECT > 4 )  ! Tags if > 4 species
    ELSE
       Input_Opt%LSPLIT = ( Input_Opt%N_ADVECT > 1 )  ! Tags if > 1 species
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
    USE RoundOff_Mod,  ONLY : Cast_and_RoundOff
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
    thisLoc = ' -> at Config_Aerosol (in module GeosCore/input_mod.F90)'

    !------------------------------------------------------------------------
    ! Use online carbon aerosols?
    !------------------------------------------------------------------------
    key    = "aerosols%carbon%activate"
    v_bool = MISSING_BOOL
    CALL QFYAML_Add_Get( Config, TRIM( key ), v_bool, "", RC )
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
    CALL QFYAML_Add_Get( Config, TRIM( key ), v_bool, "", RC )
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
    CALL QFYAML_Add_Get( Config, TRIM( key ), v_bool, "", RC )
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
    CALL QFYAML_Add_Get( Config, TRIM( key ), v_str, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Input_Opt%BCAE_1 = Cast_and_RoundOff( v_str, places=2 )

    !------------------------------------------------------------------------
    ! Define BC absorption enhancement (xnw, 8/24/15)
    !------------------------------------------------------------------------
    key   = "aerosols%carbon%enhance_black_carbon_absorption%hydrophobic"
    v_str = MISSING_STR
    CALL QFYAML_Add_Get( Config, TRIM( key ), v_str, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Input_Opt%BCAE_2 = Cast_and_RoundOff( v_str, places=2 )

    !------------------------------------------------------------------------
    ! Use secondary organic aerosols?
    !------------------------------------------------------------------------
    key    = "aerosols%complex_SOA%activate"
    v_bool = MISSING_BOOL
    CALL QFYAML_Add_Get( Config, TRIM( key ), v_bool, "", RC )
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
    CALL QFYAML_Add_Get( Config, TRIM( key ), v_bool, "", RC )
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
    CALL QFYAML_Add_Get( Config, TRIM( key ), v_bool, "", RC )
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
    CALL QFYAML_Add_Get( Config, TRIM( key ), v_bool, "", RC )
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
    CALL QFYAML_Add_Get( Config, TRIM( key ), v_bool, "", RC )
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
    CALL QFYAML_Add_Get( Config, TRIM( key ), a_str, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Input_Opt%SALA_Redge_um(1) = Cast_and_RoundOff( a_str(1), places=2 )
    Input_Opt%SALA_Redge_um(2) = Cast_and_RoundOff( a_str(2), places=2 )

    !------------------------------------------------------------------------
    ! Coarse mode seasalt radii bin edges [um]
    !------------------------------------------------------------------------
    key   = "aerosols%sea_salt%SALC_radius_bin_in_um"
    a_str = MISSING_STR
    CALL QFYAML_Add_Get( Config, TRIM( key ), a_str, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Input_Opt%SALC_Redge_um(1) = Cast_and_RoundOff( a_str(1), places=2 )
    Input_Opt%SALC_Redge_um(2) = Cast_and_RoundOff( a_str(2), places=2 )

    !------------------------------------------------------------------------
    ! Use marine organic aerosols?
    !------------------------------------------------------------------------
    key    = "aerosols%sea_salt%marine_organic_aerosols"
    v_bool = MISSING_BOOL
    CALL QFYAML_Add_Get( Config, TRIM( key ), v_bool, "", RC )
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
    CALL QFYAML_Add_Get( Config, TRIM( key ), v_bool, "", RC )
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
    CALL QFYAML_Add_Get( Config, TRIM( key ), v_bool, "", RC )
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
    CALL QFYAML_Add_Get( Config, TRIM( key ), v_bool, "", RC )
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
    CALL QFYAML_Add_Get( Config, TRIM( key ), v_bool, "", RC )
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
    CALL QFYAML_Add_Get( Config, TRIM( key ), v_str, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Input_Opt%T_NAT_SUPERCOOL = Cast_and_RoundOff( v_str, places=2 )

    !------------------------------------------------------------------------
    ! Ice supersaturation ratio requirement
    !------------------------------------------------------------------------
    key   = "aerosols%stratosphere%supersat_factor_req_for_ice_nucl"
    v_str = MISSING_STR
    CALL QFYAML_Add_Get( Config, TRIM( key ), v_str, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Input_Opt%P_ICE_SUPERSAT = Cast_and_RoundOff( v_str, places=2 )

    !------------------------------------------------------------------------
    ! Include stratospheric aerosols optical depths?
    !------------------------------------------------------------------------
    key    = "aerosols%stratosphere%calc_strat_aod"
    v_bool = MISSING_BOOL
    CALL QFYAML_Add_Get( Config, TRIM( key ), v_bool, "", RC )
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
    CALL QFYAML_Add_Get( Config, TRIM( key ), v_bool, "", RC )
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
    CALL QFYAML_Add_Get( Config, TRIM( key ), v_bool, "", RC )
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

    !========================================================================
    ! Print to screen
    !========================================================================
    IF ( Input_Opt%amIRoot ) THEN
       WRITE( 6, 90  ) 'AEROSOL SETTINGS'
       WRITE( 6, 95  ) '----------------'
       WRITE( 6, 100 ) 'Online SULFATE AEROSOLS?    : ', Input_Opt%LSULF
       WRITE( 6, 100 ) 'Metal catalyzed SO2 ox.?    : ', Input_Opt%LMETALCATSO2
       WRITE( 6, 100 ) 'Online CARBON AEROSOLS?     : ', Input_Opt%LCARB
       WRITE( 6, 100 ) 'Brown Carbon Aerosol?       : ', Input_Opt%LBRC
       WRITE( 6, 100 ) 'BC Absorption Enhancement?  : ', Input_Opt%LBCAE
       WRITE( 6, 105 ) 'Hydrophilic BC AE factor    : ', Input_Opt%BCAE_1
       WRITE( 6, 105 ) 'Hydrophobic BC AE factor    : ', Input_Opt%BCAE_2
       WRITE( 6, 100 ) 'Online COMPLEX SOA?         : ', Input_Opt%LSOA
       WRITE( 6, 100 ) 'Semivolatile POA?           : ', Input_Opt%LSVPOA
       WRITE( 6, 100 ) 'Online DUST AEROSOLS?       : ', Input_Opt%LDUST
       WRITE( 6, 100 ) 'Acid uptake on dust?        : ', Input_Opt%LDSTUP
       WRITE( 6, 100 ) 'Online SEA SALT AEROSOLS?   : ', Input_Opt%LSSALT
       WRITE( 6, 110 ) 'Accum  SEA SALT radii [um]  : ',                     &
                                                 Input_Opt%SALA_REDGE_um(1), &
                                                 Input_Opt%SALA_REDGE_um(2)
       WRITE( 6, 110 ) 'Coarse SEA SALT radii [um]  : ',                     &
                                                 Input_Opt%SALC_REDGE_um(1), &
                                                 Input_Opt%SALC_REDGE_um(2)
       WRITE( 6, 100 ) 'MARINE ORGANIC AEROSOLS?    : ', Input_Opt%LMPOA
       WRITE( 6, 100 ) 'Settle strat. aerosols?     : ', Input_Opt%LGRAVSTRAT
       WRITE( 6, 100 ) 'Online SOLID PSC aerosols?  : ', Input_Opt%LSOLIDPSC
       WRITE( 6, 100 ) 'Allow hom. NAT nucleation?  : ', Input_Opt%LHOMNUCNAT
       WRITE( 6, 120 ) 'NAT supercooling requirement: ',                     &
                                                   Input_Opt%T_NAT_SUPERCOOL
       WRITE( 6, 120 ) 'Ice supersaturation req.    : ',                     &
                                     ((Input_Opt%P_ICE_SUPERSAT-1)*1.e+2_fp)
       WRITE( 6, 100 ) 'Perform PSC het. chemistry? : ', Input_Opt%LPSCCHEM
       WRITE( 6, 100 ) 'Use strat. aerosol OD?      : ', Input_Opt%LSTRATOD
    ENDIF

90  FORMAT( /, A                 )
95  FORMAT( A                    )
100 FORMAT( A, L5                )
105 FORMAT( A, f8.2              )
110 FORMAT( A, f8.2, ' - ', f8.2 )
120 FORMAT( A, f8.2, 'K'         )

  END SUBROUTINE Config_Aerosol
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: config_co
!
! !DESCRIPTION: Copies CO simulation information from the Config object
!  to Input_Opt, and does necessary checks.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Config_CO( Config, Input_Opt, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod,      ONLY : OptInput
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
!
! !LOCAL VARIABLES:
!
    ! Scalars
    LOGICAL                      :: v_bool

    ! Strings
    CHARACTER(LEN=255)           :: thisLoc
    CHARACTER(LEN=512)           :: errMsg
    CHARACTER(LEN=QFYAML_NamLen) :: key

    !========================================================================
    ! Config_CO begins here!
    !========================================================================

    ! Initialize
    RC      = GC_SUCCESS
    errMsg  = ''
    thisLoc = ' -> at Config_CO2 (in module GeosCore/input_mod.F90)'

    !------------------------------------------------------------------------
    ! Use P(CO) from CH4 (archived from a fullchem simulation)?
    !------------------------------------------------------------------------
    key    = "CO_simulation_options%use_fullchem_PCO_from_CH4"
    v_bool = MISSING_BOOL
    CALL QFYAML_Add_Get( Config, key, v_bool, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Input_Opt%LPCO_CH4 = v_bool

    !------------------------------------------------------------------------
    ! Use P(CO) from NMVOC (archived from a fullchem simulation)?
    !------------------------------------------------------------------------
    key    = "CO_simulation_options%use_fullchem_PCO_from_NMVOC"
    v_bool = MISSING_BOOL
    CALL QFYAML_Add_Get( Config, key, v_bool, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Input_Opt%LPCO_NMVOC = v_bool

    !========================================================================
    ! Print to screen
    !========================================================================
    IF ( Input_Opt%ITS_A_TAGCO_SIM .and. Input_Opt%amIRoot ) THEN
       WRITE(6,90 ) 'TAGGED CO SIMULATION SETTINGS'
       WRITE(6,95 ) '(overwrites any other settings related to CO)'
       WRITE(6,95 ) '---------------------------------------------'
       WRITE(6,100) 'Use full chem. P(CO) from CH4?   :', Input_Opt%LPCO_CH4
       WRITE(6,100) 'Use full chem. P(CO) from NMVOC? :', Input_Opt%LPCO_NMVOC
    ENDIF

    ! FORMAT statements
90  FORMAT( /, A   )
95  FORMAT( A      )
100 FORMAT( A, L5  )

  END SUBROUTINE Config_CO
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: config_co2
!
! !DESCRIPTION: Copies CO2 simulation information from the Config object
!  to Input_Opt, and does necessary checks.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Config_CO2( Config, Input_Opt, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod,      ONLY : OptInput
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

    ! Strings
    CHARACTER(LEN=255)           :: thisLoc
    CHARACTER(LEN=512)           :: errMsg
    CHARACTER(LEN=QFYAML_NamLen) :: key

    !========================================================================
    ! Config_CO2 begins here!
    !========================================================================

    ! Initialize
    RC      = GC_SUCCESS
    errMsg  = ''
    thisLoc = ' -> at Config_CO2 (in module GeosCore/input_mod.F90)'

    !------------------------------------------------------------------------
    ! Use Fossil Fuel emissions?
    !------------------------------------------------------------------------
    key    = "CO2_simulation_options%sources%fossil_fuel_emissions"
    v_bool = MISSING_BOOL
    CALL QFYAML_Add_Get( Config, key, v_bool, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Input_Opt%LFOSSIL = v_bool

    !------------------------------------------------------------------------
    ! Use Ocean Exchange?
    !------------------------------------------------------------------------
    key    = "CO2_simulation_options%sources%ocean_exchange"
    v_bool = MISSING_BOOL
    CALL QFYAML_Add_Get( Config, key, v_bool, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Input_Opt%LOCEAN = v_bool

    !------------------------------------------------------------------------
    ! Turn on (balanced) biosphere with diurnal cycle?
    !------------------------------------------------------------------------
    key    = "CO2_simulation_options%sources%balanced_biosphere_exchange"
    v_bool = MISSING_BOOL
    CALL QFYAML_Add_Get( Config, key, v_bool, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Input_Opt%LBIODIURNAL = v_bool

    !------------------------------------------------------------------------
    ! Use Net Terrestrial Exchange Climatology?
    !------------------------------------------------------------------------
    key    = "CO2_simulation_options%sources%net_terrestrial_exchange"
    v_bool = MISSING_BOOL
    CALL QFYAML_Add_Get( Config, key, v_bool, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Input_Opt%LBIONETCLIM = v_bool

    !------------------------------------------------------------------------
    ! Turn on Ship emissions?
    !------------------------------------------------------------------------
    key    = "CO2_simulation_options%sources%ship_emissions"
    v_bool = MISSING_BOOL
    CALL QFYAML_Add_Get( Config, key, v_bool, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Input_Opt%LSHIP = v_bool

    !------------------------------------------------------------------------
    ! Turn on Aviation emissions?
    !------------------------------------------------------------------------
    key    = "CO2_simulation_options%sources%aviation_emissions"
    v_bool = MISSING_BOOL
    CALL QFYAML_Add_Get( Config, key, v_bool, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Input_Opt%LPLANE = v_bool

    !------------------------------------------------------------------------
    ! Turn on CO2 3D chemical source and surface correction?
    !------------------------------------------------------------------------
    key    = "CO2_simulation_options%sources%3D_chemical_oxidation_source"
    v_bool = MISSING_BOOL
    CALL QFYAML_Add_Get( Config, key, v_bool, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Input_Opt%LCHEMCO2 = v_bool

    !------------------------------------------------------------------------
    ! Background CO2 (no emissions or exchange) for Tagged-CO2 runs
    !------------------------------------------------------------------------
    key = "CO2_simulation_options%tagged_species%save_fossil_fuel_in_background"
    v_bool = MISSING_BOOL
    CALL QFYAML_Add_Get( Config, key, v_bool, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Input_Opt%LFFBKGRD = v_bool

    !------------------------------------------------------------------------
    ! Turn on biosphere and ocean exchange region tagged species?
    !------------------------------------------------------------------------
    key    = "CO2_simulation_options%tagged_species%tag_bio_and_ocean_CO2"
    v_bool = MISSING_BOOL
    CALL QFYAML_Add_Get( Config, key, v_bool, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Input_Opt%LBIOSPHTAG = v_bool

    !------------------------------------------------------------------------
    ! Turn on fossil fuel emission region tagged species?
    !------------------------------------------------------------------------
    key    = "CO2_simulation_options%tagged_species%tag_land_fossil_fuel_CO2"
    v_bool = MISSING_BOOL
    CALL QFYAML_Add_Get( Config, key, v_bool, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Input_Opt%LFOSSILTAG = v_bool

    !------------------------------------------------------------------------
    ! Turn on global ship emissions tagged species?
    !------------------------------------------------------------------------
    key    = "CO2_simulation_options%tagged_species%tag_global_ship_CO2"
    v_bool = MISSING_BOOL
    CALL QFYAML_Add_Get( Config, key, v_bool, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Input_Opt%LSHIPTAG = v_bool

    !------------------------------------------------------------------------
    ! Turn on global aircraft emissions tagged species?
    !------------------------------------------------------------------------
    key    = "CO2_simulation_options%tagged_species%tag_global_aircraft_CO2"
    v_bool = MISSING_BOOL
    CALL QFYAML_Add_Get( Config, key, v_bool, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Input_Opt%LPLANETAG = v_bool

    !=================================================================
    ! Print to screen
    !=================================================================
    IF ( Input_Opt%ITS_A_CO2_SIM .and. Input_Opt%amIRoot ) THEN
       WRITE( 6,90  ) 'CO2 SIMULATION SETTINGS'
       WRITE( 6,95  ) '(overwrites any other settings related to CO2)'
       WRITE( 6,95  ) '----------------------------------------------'
       WRITE( 6,100 ) 'National Fossil Fuel Emission :', Input_Opt%LFOSSIL
       WRITE( 6,100 ) 'Ocean CO2 Uptake/Emission     :', Input_Opt%LOCEAN
       WRITE( 6,100 ) 'Biosphere seas/diurnal cycle  :', Input_Opt%LBIODIURNAL
       WRITE( 6,100 ) 'Net Terr Exch - Climatology   :', Input_Opt%LBIONETCLIM
       WRITE( 6,100 ) 'Intl/Domestic Ship emissions  :', Input_Opt%LSHIP
       WRITE( 6,100 ) 'Intl/Domestic Aviation emiss  :', Input_Opt%LPLANE
       WRITE( 6,100 ) 'CO2 from oxidation (CO,CH4,..):', Input_Opt%LCHEMCO2
       WRITE( 6, 95 ) 'Tagged CO2 settings'
       WRITE( 6,100 ) '  Save Fossil CO2 in Bckgrnd  :', Input_Opt%LFFBKGRD
       WRITE( 6,100 ) '  Tag Biosphere/Ocean CO2     :', Input_Opt%LBIOSPHTAG
       WRITE( 6,100 ) '  Tag Fossil Fuel CO2         :', Input_Opt%LFOSSILTAG
       WRITE( 6,100 ) '  Tag Global Ship CO2         :', Input_Opt%LSHIPTAG
       WRITE( 6,100 ) '  Tag Global Aviation CO2     :', Input_Opt%LPLANETAG
    ENDIF

    ! FORMAT statements
90  FORMAT( /, A  )
95  FORMAT( A     )
100 FORMAT( A, L5 )

  END SUBROUTINE Config_CO2
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
    USE RoundOff_Mod,  ONLY : Cast_and_RoundOff
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
    errMsg  = ''
    thisLoc = ' -> at Config_Chemistry (in module GeosCore/input_mod.F90)'

    !------------------------------------------------------------------------
    ! Turn on chemistry?
    !------------------------------------------------------------------------
    key    = "operations%chemistry%activate"
    v_bool = MISSING_BOOL
    CALL QFYAML_Add_Get( Config, TRIM( key ), v_bool, "", RC )
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
    CALL QFYAML_Add_Get( Config, TRIM( key ), v_bool, "", RC )
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
    CALL QFYAML_Add_Get( Config, TRIM( key ), v_bool, "", RC )
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
    CALL QFYAML_Add_Get( Config, TRIM( key ), v_bool, "", RC )
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
    CALL QFYAML_Add_Get( Config, TRIM( key ), v_bool, "", RC )
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
    CALL QFYAML_Add_Get( Config, TRIM( key ), v_str, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Input_Opt%GAMMA_HO2 = Cast_and_RoundOff( v_str, places=2 )

    !------------------------------------------------------------------------
    ! Auto-reduce solver options (hplin, 10/3/22)
    ! autoreduce_solver:
    !   activate: false
    !   use_target_threshold:
    !     activate: true
    !     oh_tuning_factor: 0.00005
    !     no2_tuning_factor: 0.0001
    !   use_absolute_threshold:
    !     scale_by_pressure: true
    !     absolute_threshold: 100.0
    !   keep_halogens_active: false
    !   append_in_internal_timestep: false
    !------------------------------------------------------------------------
    key   = "operations%chemistry%autoreduce_solver%activate"
    v_bool = MISSING_BOOL
    CALL QFYAML_Add_Get( Config, TRIM( key ), v_bool, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Input_Opt%USE_AUTOREDUCE = v_bool

    ! Use target species (OH, NO2) based threshold?
    key   = "operations%chemistry%autoreduce_solver%use_target_threshold%activate"
    v_bool = MISSING_BOOL
    CALL QFYAML_Add_Get( Config, TRIM( key ), v_bool, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Input_Opt%AUTOREDUCE_IS_KEY_THRESHOLD = v_bool

    ! ... OH and NO2 tuning factors?
    key   = "operations%chemistry%autoreduce_solver%use_target_threshold%oh_tuning_factor"
    v_str = MISSING_STR
    CALL QFYAML_Add_Get( Config, TRIM( key ), v_str, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Input_Opt%AUTOREDUCE_TUNING_OH = Cast_and_RoundOff( v_str, places=0 )

    key   = "operations%chemistry%autoreduce_solver%use_target_threshold%no2_tuning_factor"
    v_str = MISSING_STR
    CALL QFYAML_Add_Get( Config, TRIM( key ), v_str, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Input_Opt%AUTOREDUCE_TUNING_NO2 = Cast_and_RoundOff( v_str, places=0 )

    ! If not target species, absolute rate threshold
    key   = "operations%chemistry%autoreduce_solver%use_absolute_threshold%absolute_threshold"
    v_str = MISSING_STR
    CALL QFYAML_Add_Get( Config, TRIM( key ), v_str, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Input_Opt%AUTOREDUCE_THRESHOLD = Cast_and_RoundOff( v_str, places=0 )

    ! Would this absolute threshold be scaled by pressure?
    key   = "operations%chemistry%autoreduce_solver%use_absolute_threshold%scale_by_pressure"
    v_bool = MISSING_BOOL
    CALL QFYAML_Add_Get( Config, TRIM( key ), v_bool, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Input_Opt%AUTOREDUCE_IS_PRS_THRESHOLD = v_bool

    ! Keep halogens active?
    key   = "operations%chemistry%autoreduce_solver%keep_halogens_active"
    v_bool = MISSING_BOOL
    CALL QFYAML_Add_Get( Config, TRIM( key ), v_bool, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Input_Opt%AUTOREDUCE_IS_KEEPACTIVE = v_bool

    ! Append species over the course of the external time step
    ! (aka. in internal timesteps?)
    key   = "operations%chemistry%autoreduce_solver%append_in_internal_timestep"
    v_bool = MISSING_BOOL
    CALL QFYAML_Add_Get( Config, TRIM( key ), v_bool, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Input_Opt%AUTOREDUCE_IS_APPEND = v_bool

    ! Return success
    RC = GC_SUCCESS

    !========================================================================
    ! Print to screen
    !========================================================================
    IF ( Input_Opt%amIRoot ) THEN
       WRITE( 6,90  ) 'CHEMISTRY SETTINGS'
       WRITE( 6,95  ) '------------------'
       WRITE( 6,100 ) 'Turn on chemistry?          : ', Input_Opt%LCHEM
       WRITE( 6,100 ) 'Use linear. mesospheric chem: ', Input_Opt%LINEAR_CHEM
       WRITE( 6,100 ) ' => Use Linoz for O3?       : ', Input_Opt%LLINOZ
       WRITE( 6,100 ) 'Online strat. H2O?          : ', Input_Opt%LACTIVEH2O
       WRITE( 6,100 ) 'Use robust strat H2O BC?    : ', Input_Opt%LStaticH2OBC
       WRITE( 6,110 ) 'GAMMA HO2                   : ', Input_Opt%GAMMA_HO2
       WRITE( 6,100 ) 'Use auto-reduce solver?     : ', Input_Opt%USE_AUTOREDUCE
       IF ( Input_Opt%AUTOREDUCE_IS_KEY_THRESHOLD ) THEN
         WRITE( 6,100 ) 'Use target species threshold: ', Input_Opt%AUTOREDUCE_IS_KEY_THRESHOLD
         WRITE( 6,130 ) 'OH tuning factor:             ', Input_Opt%AUTOREDUCE_TUNING_OH
         WRITE( 6,130 ) 'NO2 tuning factor:            ', Input_Opt%AUTOREDUCE_TUNING_NO2
       ELSE
         WRITE( 6,120 ) 'Absolute AR threshold       : ', Input_Opt%AUTOREDUCE_THRESHOLD
         WRITE( 6,100 ) 'Use prs. dependent thres?   : ', Input_Opt%AUTOREDUCE_IS_PRS_THRESHOLD
       ENDIF
       WRITE( 6,100 ) 'Keep halogen spec. active?  : ', Input_Opt%AUTOREDUCE_IS_KEEPACTIVE
       WRITE( 6,100 ) 'Use append in auto-reduce?  : ', Input_Opt%AUTOREDUCE_IS_APPEND
    ENDIF

    ! FORMAT statements
90  FORMAT( /, A    )
95  FORMAT( A       )
100 FORMAT( A, L5   )
110 FORMAT( A, F4.2 )
120 FORMAT( A, F5.1 )
130 FORMAT( A, ES7.1 )

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
    USE RoundOff_Mod,  ONLY : Cast_and_RoundOff
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
    thisLoc = ' -> at Config_RRTMG (in module GeosCore/input_mod.F90)'

    !------------------------------------------------------------------------
    ! Turn on RRTMG?
    !------------------------------------------------------------------------
    key    = "operations%rrtmg_rad_transfer_model%activate"
    v_bool = MISSING_BOOL
    CALL QFYAML_Add_Get( Config, TRIM( key ), v_bool, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Input_Opt%LRAD = v_bool

    !------------------------------------------------------------------------
    ! AOD wavelength selection? (You can have up to 3)
    !------------------------------------------------------------------------
    key   = "operations%rrtmg_rad_transfer_model%aod_wavelengths_in_nm"
    a_str = MISSING_STR
    CALL QFYAML_Add_Get( Config, TRIM( key ), a_str,                          &
                         "",     RC,          dynamic_size=.TRUE.            )
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
       Input_Opt%WvSelect(I)    = Cast_and_RoundOff( a_str(N), places=2 )
    ENDDO

    !------------------------------------------------------------------------
    ! Turn on LW radiation calculation?
    !------------------------------------------------------------------------
    key    = "operations%rrtmg_rad_transfer_model%longwave_fluxes"
    v_bool = MISSING_BOOL
    CALL QFYAML_Add_Get( Config, TRIM( key ), v_bool, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Input_Opt%LLWRAD = v_bool

    !------------------------------------------------------------------------
    ! Turn on SW radiation calculation?
    !------------------------------------------------------------------------
    key    = "operations%rrtmg_rad_transfer_model%shortwave_fluxes"
    v_bool = MISSING_BOOL
    CALL QFYAML_Add_Get( Config, TRIM( key ), v_bool, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Input_Opt%LSWRAD = v_bool

    !------------------------------------------------------------------------
    ! Calculate for clear-sky?
    !------------------------------------------------------------------------
    key    = "operations%rrtmg_rad_transfer_model%clear_sky_flux"
    v_bool = MISSING_BOOL
    CALL QFYAML_Add_Get( Config, TRIM( key ), v_bool, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Input_Opt%LSKYRAD(1) = v_bool

    !------------------------------------------------------------------------
    ! Calculate for all-sky?
    !------------------------------------------------------------------------
    key    = "operations%rrtmg_rad_transfer_model%all_sky_flux"
    v_bool = MISSING_BOOL
    CALL QFYAML_Add_Get( Config, TRIM( key ), v_bool, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Input_Opt%LSKYRAD(2) = v_bool

    !========================================================================
    ! Error check settings
    !========================================================================
#ifndef RRTMG
    ! Use of RRTMG necessitates recompilation
    IF ( Input_Opt%LRAD ) THEN
       errMsg = 'LRAD=T but RRTMG is not defined at compile time!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
#endif

    ! Make sure radiation switches are turned off if RRTMG is off
    IF ( .not. Input_Opt%LRAD ) THEN

       IF ( Input_Opt%LLWRAD ) THEN
          errMsg = 'Cannot have LW fluxes turned on without RRTMG!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       IF ( Input_Opt%LSWRAD ) THEN
          errMsg = 'Cannot have SW fluxes turned on without RRTMG!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       IF ( Input_Opt%LSKYRAD(1) ) THEN
          errMsg = 'Cannot have clear-sky flux turned on without RRTMG!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       IF ( Input_Opt%LSKYRAD(2) ) THEN
          errMsg = 'Cannot have all-sky flux turned on without RRTMG!'
          CALL GC_Error( errMsg, RC, thisLoc )
       ENDIF
    ENDIF

    !========================================================================
    ! Print to screen
    !========================================================================
    IF ( Input_Opt%amIRoot ) THEN
       WRITE( 6, 90 ) 'RRTMG SETTINGS'
       WRITE( 6, 95 ) '--------------'
       DO N=1, Input_Opt%NWVSELECT
          WRITE( 6, 115     ) 'AOD output wavelength (nm)  : ',              &
                               Input_Opt%WVSELECT(N)
       ENDDO
       WRITE( 6, 100 ) 'Turn on radiation?          : ', Input_Opt%LRAD
       WRITE( 6, 100 ) 'Consider longwave?          : ', Input_Opt%LLWRAD
       WRITE( 6, 100 ) 'Consider shortwave?         : ', Input_Opt%LSWRAD
       WRITE( 6, 100 ) 'Clear-sky flux?             : ', Input_Opt%LSKYRAD(1)
       WRITE( 6, 100 ) 'All-sky flux?               : ', Input_Opt%LSKYRAD(2)
    ENDIF

    ! FORMAT statements
90  FORMAT( /, A    )
95  FORMAT( A       )
100 FORMAT( A, L5   )
110 FORMAT( A, I5   )
115 FORMAT( A, F7.1 )
120 FORMAT( A, 11I1 )

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
    USE RoundOff_Mod,  ONLY : Cast_and_RoundOff
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

    !========================================================================
    ! Config_Photolysis begins here!
    !========================================================================

    ! Initialize
    RC      = GC_SUCCESS
    errMsg  = ''
    thisLoc = ' -> at Config_Photolysis (in module GeosCore/input_mod.F90)'

    !------------------------------------------------------------------------
    ! Directory with photolysis input files
    !------------------------------------------------------------------------
    key   = "operations%photolysis%input_dir"
    v_str = MISSING_STR
    CALL QFYAML_Add_Get( Config, TRIM( key ), v_str, "", RC )
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
    CALL QFYAML_Add_Get( Config, TRIM( key ), v_bool, "", RC )
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
    CALL QFYAML_Add_Get( Config, TRIM( key ), v_bool, "", RC )
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
    CALL QFYAML_Add_Get( Config, TRIM( key ), v_bool, "", RC )
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
    CALL QFYAML_Add_Get( Config, TRIM( key ), v_bool, "", RC )
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
    CALL QFYAML_Add_Get( Config, TRIM( key ), v_str, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Input_Opt%hvAerNIT_JNITs = Cast_and_RoundOff( v_str, places=3 )

    !------------------------------------------------------------------------
    ! scalar for JHNO3 for photoylsing NIT aerosol (TMS, 23/08/18)
    !------------------------------------------------------------------------
    key    = &
     "operations%photolysis%photolyze_nitrate_aerosol%NIT_Jscale_JHNO2"
    v_str = MISSING_STR
    CALL QFYAML_Add_Get( Config, TRIM( key ), v_str, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Input_Opt%hvAerNIT_JNIT = Cast_and_RoundOff( v_str, places=3 )

    !------------------------------------------------------------------------
    ! Fraction for JNITS/NIT channel A (HNO2) for NITs photoylsis
    !------------------------------------------------------------------------
    key   = &
     "operations%photolysis%photolyze_nitrate_aerosol%percent_channel_A_HONO"
    v_str = MISSING_STR
    CALL QFYAML_Add_Get( Config, TRIM( key ), v_str, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Input_Opt%JNITChanA = Cast_and_RoundOff( v_str, places=3 )

    !------------------------------------------------------------------------
    ! Fraction for JNITs/NIT channel B (NO2) for NITs photoylsis
    !------------------------------------------------------------------------
    key    = &
     "operations%photolysis%photolyze_nitrate_aerosol%percent_channel_B_NO2"
    v_str = MISSING_STR
    CALL QFYAML_Add_Get( Config, TRIM( key ), v_str, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Input_Opt%JNITChanB = Cast_and_RoundOff( v_str, places=3 )

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

    ! FAST-JX is only used for fullchem and mercury
    IF ( Input_Opt%ITS_A_FULLCHEM_SIM   .or.                                 &
         Input_Opt%ITS_AN_AEROSOL_SIM   .or.                                 &
         Input_Opt%ITS_A_MERCURY_SIM  ) THEN

       ! Make sure either O3 from met or TOMS is selected
       IF ( .not. Input_Opt%USE_O3_FROM_MET   .and.                          &
            .not. Input_Opt%USE_TOMS_O3     ) THEN
          errMsg = 'Must select either O3 from met or TOMS/SBUV O3'          &
                // 'for O3 values above the chemistry grid!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
       IF ( Input_Opt%USE_O3_FROM_MET .and. &
            Input_Opt%USE_TOMS_O3 ) THEN
          errMsg = 'Must select either O3 from met or TOMS/SBUV O3'          &
                // 'for O3 values above the chemistry grid!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       ! Make sure the aerosol-only simulation gets O3 from met or TOMS
       IF ( Input_Opt%ITS_AN_AEROSOL_SIM ) THEN
          IF ( Input_Opt%USE_ONLINE_O3 ) THEN
             errMsg= 'Cannot use online O3 for specialty simulations! '      &
                  // 'Select O3 from met or TOMS O3 instead.'
             CALL GC_Error( errMsg, RC, thisLoc )
             RETURN
          ENDIF
       ENDIF

    ELSE

       ! If not a simulation that uses photolysis, set options to FALSE
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

    !========================================================================
    ! Print to screen
    !========================================================================
    IF ( Input_Opt%amIRoot ) THEN
       WRITE( 6,90  ) 'PHOTOLYSIS SETTINGS'
       WRITE( 6,95  ) '-------------------'
       WRITE( 6,120 ) 'FAST-JX input directory     : ',                      &
                       TRIM( Input_Opt%FAST_JX_DIR )
       WRITE( 6,100 ) 'Online ozone for FAST-JX?   : ', Input_Opt%USE_ONLINE_O3
       WRITE( 6,100 ) 'Ozone from met for FAST-JX? : ',                      &
                       Input_Opt%USE_O3_FROM_MET
       WRITE( 6,100 ) 'TOMS/SBUV ozone for FAST-JX?: ', Input_Opt%USE_TOMS_O3
       WRITE( 6,100 ) 'Photolyse nitrate aerosol?  : ', Input_Opt%hvAerNIT
       WRITE( 6,105 ) 'JNITs scaling of JHNO3      : ', Input_Opt%hvAerNIT_JNITs
       WRITE( 6,105 ) 'JNIT scaling of JHNO3       : ', Input_Opt%hvAerNIT_JNIT
       WRITE( 6,105 ) 'JNIT(s) channel A (HONO)    : ', Input_Opt%JNITChanA
       WRITE( 6,105 ) 'JNIT(s) channel B (NO2)     : ', Input_Opt%JNITChanB
       ! Write more info
       IF ( Input_Opt%USE_ONLINE_O3 ) THEN
          WRITE( 6, 95 ) ''
          WRITE( 6, 95 ) 'NOTE ABOUT OVERHEAD O3 FOR FAST-JX:'
          WRITE( 6, 95 ) ' Online O3 from GEOS-Chem will be used'
          WRITE( 6, 95 ) ' to weight the O3 column within the'
          WRITE( 6, 95 ) ' chemistry grid and O3 from met or TOMS'
          WRITE( 6, 95 ) ' will be used outside the chemistry grid.'
       ENDIF
    ENDIF

    ! FORMAT statements
90  FORMAT ( /, A    )
95  FORMAT( A       )
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
    CALL QFYAML_Add_Get( Config, TRIM( key ), v_bool, "", RC )
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
    CALL QFYAML_Add_Get( Config, TRIM( key ), v_bool, "", RC )
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
    CALL QFYAML_Add_Get( Config, TRIM( key ), v_bool, "", RC )
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

    !=========================================================================
    ! Print to screen
    !=========================================================================
    IF ( Input_Opt%amIRoot ) THEN
       WRITE( 6, 90  ) 'CONVECTION SETTINGS'
       WRITE( 6, 95  ) '-------------------'
       WRITE( 6, 100 ) 'Turn on cloud convection?   : ', Input_Opt%LCONV

       WRITE( 6, 90  ) 'PBL MIXING SETTINGS'
       WRITE( 6, 95  ) '-------------------'
       WRITE( 6, 100 ) 'Turn on PBL mixing?         : ', Input_Opt%LTURB
       WRITE( 6, 100 ) 'Turn on non-local PBL?      : ', Input_Opt%LNLPBL
    ENDIF

    ! FORMAT statements
90  FORMAT( /, A  )
95  FORMAT( A     )
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
! !DESCRIPTION: Copies drydep and wetdep information from the Config object
!  to Input_Opt, and does necessary checks.
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
    USE RoundOff_Mod,  ONLY : Cast_and_RoundOff
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
    thisLoc = ' -> at Config_DryDep_WetDep (in module GeosCore/input_mod.F90)'

    !------------------------------------------------------------------------
    ! Turn on drydep?
    !------------------------------------------------------------------------
    key    = "operations%dry_deposition%activate"
    v_bool = MISSING_BOOL
    CALL QFYAML_Add_Get( Config, TRIM( key ), v_bool, "", RC )
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
    CALL QFYAML_Add_Get( Config, TRIM( key ), v_bool, "", RC )
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
    CALL QFYAML_Add_Get( Config, TRIM( key ), v_str, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Input_Opt%CO2_LEVEL = Cast_and_RoundOff( v_str, places=2 )

    !------------------------------------------------------------------------
    ! Reference CO2 level
    !------------------------------------------------------------------------
    key   = "operations%dry_deposition%CO2_effect%reference_CO2_level"
    v_str = MISSING_STR
    CALL QFYAML_Add_Get( Config, TRIM( key ), v_str, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Input_Opt%CO2_REF = Cast_and_RoundOff( v_str, places=2 )

    !------------------------------------------------------------------------
    ! Diag for RA_alt above surface in meters
    !------------------------------------------------------------------------
    key   = "operations%dry_deposition%diag_alt_above_sfc_in_m"
    v_int = MISSING_INT
    CALL QFYAML_Add_Get( Config, TRIM( key ), v_int, "", RC )
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
    CALL QFYAML_Add_Get( Config, TRIM( key ), v_bool, "", RC )
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
       WRITE( 6, 90  ) 'DRY DEPOSITION SETTINGS'
       WRITE( 6, 95  ) '-----------------------'
       WRITE( 6, 100 ) 'Turn on dry deposition?     : ', Input_Opt%LDRYD
       WRITE( 6, 100 ) 'Dry dep over full PBL?      : ', Input_Opt%PBL_DRYDEP
       WRITE( 6, 100 ) 'Turn on CO2 effect?         : ', Input_Opt%CO2_EFFECT
       WRITE( 6, 110 ) 'CO2 level                   : ', Input_Opt%CO2_LEVEL
       WRITE( 6, 110 ) 'CO2 reference level         : ', Input_Opt%CO2_REF
       WRITE( 6, 110 ) 'RIX scaling factor          : ', Input_Opt%RS_SCALE


       WRITE( 6, 90  ) 'WET DEPOSITION SETTINGS'
       WRITE( 6, 95  ) '-----------------------'
       WRITE( 6, 100 ) 'Turn on wet deposition?     : ', Input_Opt%LWETD
    ENDIF

    ! FORMAT statements
90  FORMAT( /, A    )
95  FORMAT( A       )
100 FORMAT( A, L5   )
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
! !IROUTINE: config_gamap
!
! !DESCRIPTION: Copies GAMAP information from the Config object
!  to Input_Opt, and does necessary checks.
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
    thisLoc = ' -> at Config_Gamap (in module GeosCore/input_mod.F90)'

    !------------------------------------------------------------------------
    ! diaginfo.dat
    !------------------------------------------------------------------------
    key   = "extra_diagnostics%legacy_bpch%gamap%diaginfo_dat_file"
    v_str = MISSING_STR
    CALL QFYAML_Add_Get( Config, TRIM( key ), v_str, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Input_Opt%GAMAP_DIAGINFO = TRIM( ADJUSTL( v_str ) )

    !------------------------------------------------------------------------
    ! tracerinfo.dat
    !------------------------------------------------------------------------
    key   = "extra_diagnostics%legacy_bpch%gamap%tracerinfo_dat_file"
    v_str = MISSING_STR
    CALL QFYAML_Add_Get( Config, TRIM( key ), v_str, "", RC )
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
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: config_bpch_output

! !DESCRIPTION: Copies bpch output information from the Config object
!  to Input_Opt, and does necessary checks.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Config_Bpch_Output( Config, Input_Opt, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod,  ONLY : OptInput
    USE QFYAML_Mod,     ONLY : QFYAML_t
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
    INTEGER                      :: N
    INTEGER                      :: v_int

    ! Strings
    CHARACTER(LEN=255)           :: thisLoc
    CHARACTER(LEN=512)           :: errMsg
    CHARACTER(LEN=QFYAML_NamLen) :: key
    CHARACTER(LEN=QFYAML_StrLen) :: v_str

    ! String arrays
    CHARACTER(LEN=3)             :: mon(12)=  (/'JAN', 'FEB', 'MAR', 'APR',  &
                                                'MAY', 'JUN', 'JUL', 'AUG',  &
                                                'SEP', 'OCT', 'NOV', 'DEC'/)

    !========================================================================
    ! Config_Bpch_Output begins here!
    !========================================================================

    ! Initialize
    RC      = GC_SUCCESS
    errMsg  = ''
    thisLoc = ' -> at Config_Bpch_Output (in module GeosCore/input_mod.F90)'

    !========================================================================
    ! Read data into the Input_Opt%NJDAY array
    !========================================================================
    DO N = 1, 12

       key = "extra_diagnostics%legacy_bpch%output_menu%" // &
             "schedule_output_for_"                       // mon(N)
       v_str = MISSING_STR
       CALL QFYAML_Add_Get( Config, TRIM( key ), v_str, "", RC )
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = 'Error parsing ' // TRIM( key ) // '!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       ! Parse string into NJDAY array by month
       SELECT CASE( N )
          CASE( 1 )
             READ( v_str, '(31i1)' ) Input_Opt%NJDAY(1:31)
          CASE( 2  )
             READ( v_str, '(29i1)' ) Input_Opt%NJDAY(32:60)
          CASE( 3 )
             READ( v_str, '(31i1)' ) Input_Opt%NJDAY(61:91)
          CASE( 4  )
             READ( v_str, '(30i1)' ) Input_Opt%NJDAY(92:121)
          CASE( 5  )
             READ( v_str, '(31i1)' ) Input_Opt%NJDAY(122:152)
          CASE( 6  )
             READ( v_str, '(30i1)' ) Input_Opt%NJDAY(153:182)
          CASE( 7  )
             READ( v_str, '(31i1)' ) Input_Opt%NJDAY(183:213)
          CASE( 8  )
             READ( v_str, '(31i1)' ) Input_Opt%NJDAY(214:244)
          CASE( 9  )
             READ( v_str, '(30i1)' ) Input_Opt%NJDAY(245:274)
          CASE( 10 )
             READ( v_str, '(31i1)' ) Input_Opt%NJDAY(275:305)
          CASE( 11 )
             READ( v_str, '(30i1)' ) Input_Opt%NJDAY(306:335)
          CASE( 12 )
             READ( v_str, '(31i1)' ) Input_Opt%NJDAY(336:366)
       END SELECT
    ENDDO

    !=================================================================
    ! Print to screen
    !=================================================================
    IF( Input_Opt%amIRoot ) THEN
       WRITE( 6, '(/,a)' ) 'BPCH OUTPUT SETTINGS'
       WRITE( 6, '(  a)' ) '--------------------'
       WRITE( 6, 110     )
       WRITE( 6, 120     )
       WRITE( 6, 130     )
       WRITE( 6, 140     ) Input_Opt%NJDAY
    ENDIF

    ! FORMAT statements
110 FORMAT( '              1111111111222222222233' )
120 FORMAT( '     1234567890123456789012345678901' )
130 FORMAT( '     -------------------------------' )
140 FORMAT( 'JAN--', 31i1, /, 'FEB--', 29i1, /, 'MAR--', 31i1, /, &
            'APR--', 30i1, /, 'MAY--', 31i1, /, 'JUN--', 30i1, /, &
            'JUL--', 31i1, /, 'AUG--', 31i1, /, 'SEP--', 30i1, /, &
            'OCT--', 31i1, /, 'NOV--', 30i1, /, 'DEC--', 31i1 )

    ! Make sure we have output at end of run
    CALL IS_LAST_DAY_GOOD( Input_Opt, RC )

    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

  END SUBROUTINE Config_Bpch_Output
!EOC
#endif
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: config_planeflight
!
! !DESCRIPTION: Copies PlaneFlight diagnostic information from the Config
!  object to Input_Opt, and does necessary checks.
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
    ! Config_PlaneFlight begins here!
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
    CALL QFYAML_Add_Get( Config, TRIM( key ), v_bool, "", RC )
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
    CALL QFYAML_Add_Get( Config, TRIM( key ), v_str, "", RC )
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
    CALL QFYAML_Add_Get( Config, TRIM( key ), v_str, "", RC )

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
       WRITE( 6, 90  ) 'PLANEFLIGHT DIAGNOSTIC SETTINGS'
       WRITE( 6, 95  ) '-------------------------------'
       WRITE( 6, 100 ) 'Turn on planeflight diag?   : ',                     &
                        Input_Opt%Do_Planeflight
       WRITE( 6, 110 ) 'Flight track input file     : ',                     &
                        TRIM( Input_Opt%Planeflight_InFile )
       WRITE( 6, 110 ) 'Output file name            : ',                     &
                        TRIM( Input_Opt%Planeflight_OutFile )
    ENDIF

    ! FORMAT statements
90  FORMAT( /, A   )
95  FORMAT( A      )
100 FORMAT( A, L5  )
110 FORMAT( A, A   )

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
! !IROUTINE: config_obspack
!
! !DESCRIPTION: Copies Obspack diagnostic information from the Config
!  object to Input_Opt, and does necessary checks.
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
    ! Config_ObsPack begins here!
    !========================================================================

    ! Initialize
    RC      = GC_SUCCESS
    errMsg  = 'Error reading the "geoschem_config.yml" file!'
    thisLoc = ' -> at Config_ObsPack (in module GeosCore/input_mod.F90)'

    !------------------------------------------------------------------------
    ! Turn on ObsPack diagnostic?
    !------------------------------------------------------------------------
    key    = "extra_diagnostics%obspack%activate"
    v_bool = MISSING_BOOL
    CALL QFYAML_Add_Get( Config, TRIM( key ), v_bool, "", RC )
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
    CALL QFYAML_Add_Get( Config, TRIM( key ), v_bool, "", RC )
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
    CALL QFYAML_Add_Get( Config, TRIM( key ), v_str, "", RC )
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
    CALL QFYAML_Add_Get( Config, TRIM( key ), v_str, "", RC )
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
    CALL QFYAML_Add_Get( Config, TRIM( key ), a_str,                         &
                         "",     RC,          dynamic_size=.TRUE.           )
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
          Input_Opt%ObsPack_nSpc    = Input_Opt%N_Advect
          EXIT
       ENDIF

       ! Otherwise, increment the count and copy the Obspack species
       ! name into the Input_Opt object
       Input_Opt%ObsPack_nSpc = Input_Opt%ObsPack_nSpc + 1
       Input_Opt%ObsPack_SpcName(Input_Opt%ObsPack_nSpc) = TRIM( a_str(N) )
    ENDDO

    !========================================================================
    ! Print to screen
    !========================================================================
    IF( Input_Opt%amIRoot ) THEN
       WRITE( 6, 90  ) 'OBSPACK SETTINGS'
       WRITE( 6, 95  ) '----------------'
       WRITE( 6, 100 ) 'Turn on ObsPack diagnostic? : ',                     &
                        Input_Opt%Do_ObsPack
       WRITE( 6, 100 ) 'Suppress logfile output?    : ',                     &
                        Input_Opt%ObsPack_Quiet
       WRITE( 6, 110 ) 'ObsPack input file          : ',                     &
                        TRIM( Input_Opt%ObsPack_InputFile  )
       WRITE( 6, 110 ) 'ObsPack output file         : ',                     &
                        TRIM( Input_Opt%ObsPack_OutputFile )
    ENDIF

    ! FORMAT statements
90  FORMAT( /, A   )
95  FORMAT( A      )
100 FORMAT( A, L5  )
110 FORMAT( A, A   )

  END SUBROUTINE Config_ObsPack
!EOC
#if !(defined( EXTERNAL_GRID ) || defined( EXTERNAL_FORCING ))
#ifdef BPCH_DIAG
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: config_nd51
!
! !DESCRIPTION: Copies ND51 satellite diagnostic information from the Config
!  object to Input_Opt, and does necessary checks.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Config_ND51( Config, Input_Opt, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod, ONLY : OptInput
    USE RoundOff_Mod,  ONLY : Cast_and_RoundOff
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
    INTEGER                      :: N

    ! Arrays
    INTEGER                      :: a_int(QFYAML_MaxArr)

    ! Strings
    CHARACTER(LEN=255)           :: thisLoc
    CHARACTER(LEN=255)           :: errMsg
    CHARACTER(LEN=QFYAML_NamLen) :: key
    CHARACTER(LEN=QFYAML_StrLen) :: v_str

    ! String arrays
    CHARACTER(LEN=QFYAML_StrLen) :: a_str(QFYAML_MaxArr)


    !========================================================================
    ! Config_ND51 begins here!
    !========================================================================

    ! Initialize
    RC      = GC_SUCCESS
    errMsg  = ''
    thisLoc = ' -> at Config_ND51 (in module GeosCore/input_mod.F90)'

    !------------------------------------------------------------------------
    ! Turn on ND51 diagnostic
    !------------------------------------------------------------------------
    key    = "extra_diagnostics%legacy_bpch%ND51_satellite%activate"
    v_bool = MISSING_BOOL
    CALL QFYAML_Add_Get( Config, key, v_bool, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Input_Opt%DO_ND51 = v_bool

    !------------------------------------------------------------------------
    ! Instantaneous 3-D timeseries file
    !------------------------------------------------------------------------
    key    = "extra_diagnostics%legacy_bpch%ND51_satellite%output_file"
    v_str = MISSING_STR
    CALL QFYAML_Add_Get( Config, key, v_str, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Input_Opt%ND51_FILE = TRIM( v_str )

    !------------------------------------------------------------------------
    ! Tracers to include
    !------------------------------------------------------------------------
    key   = "extra_diagnostics%legacy_bpch%ND51_satellite%tracers"
    a_int = MISSING_INT
    CALL QFYAML_Add_Get( Config, key, a_int, "", RC, dynamic_size=.TRUE. )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    DO N = 1, SIZE( a_int )
       IF ( a_int(N) == MISSING_INT ) EXIT
       Input_Opt%ND51_TRACERS(N) = a_int(N)
    ENDDO
    Input_Opt%N_ND51 = N - 1

    !------------------------------------------------------------------------
    ! NHMS_WRITE
    !------------------------------------------------------------------------
    key   = "extra_diagnostics%legacy_bpch%ND51_satellite%UTC_hour_for_write"
    v_str = MISSING_STR
    CALL QFYAML_Add_Get( Config, key, v_str, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Input_Opt%ND51_HR_WRITE = Cast_and_RoundOff( v_str, places=2 )

    ! Make sure ND51_HR_WRITE is in the range 0-23.999 hrs
    Input_Opt%ND51_HR_WRITE = MOD( Input_Opt%ND51_HR_WRITE, 24.0_fp )

    !------------------------------------------------------------------------
    ! HR1, HR2
    !------------------------------------------------------------------------
    key = "extra_diagnostics%legacy_bpch%ND51_satellite%averaging_period_in_LT"
    a_str = MISSING_STR
    CALL QFYAML_Add_Get( Config, key, a_str, "", RC, dynamic_size=.TRUE. )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Input_Opt%ND51_HR1 = Cast_and_RoundOff( a_str(1), places=2 )
    Input_Opt%ND51_HR2 = Cast_and_RoundOff( a_str(2), places=2 )

    !------------------------------------------------------------------------
    ! IMIN, IMAX
    !------------------------------------------------------------------------
    key = "extra_diagnostics%legacy_bpch%ND51_satellite%IMIN_and_IMAX_of_region"
    a_int = MISSING_INT
    CALL QFYAML_Add_Get( Config, key, a_int, "", RC, dynamic_size=.TRUE. )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Input_Opt%ND51_IMIN = a_int(1)
    Input_Opt%ND51_IMAX = a_int(2)

    !------------------------------------------------------------------------
    ! JMIN, JMAX
    !------------------------------------------------------------------------
    key = "extra_diagnostics%legacy_bpch%ND51_satellite%JMIN_and_JMAX_of_region"
    a_int = MISSING_INT
    CALL QFYAML_Add_Get( Config, key, a_int, "", RC, dynamic_size=.TRUE. )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Input_Opt%ND51_JMIN = a_int(1)
    Input_Opt%ND51_JMAX = a_int(2)

    !------------------------------------------------------------------------
    ! LMIN, LMAX
    !------------------------------------------------------------------------
    key = "extra_diagnostics%legacy_bpch%ND51_satellite%LMIN_and_LMAX_of_region"
    a_int = MISSING_INT
    CALL QFYAML_Add_Get( Config, key, a_int, "", RC, dynamic_size=.TRUE. )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Input_Opt%ND51_LMIN = a_int(1)
    Input_Opt%ND51_LMAX = a_int(2)

    !========================================================================
    ! Print to screen
    !========================================================================
    IF ( Input_Opt%amIRoot ) THEN
       WRITE( 6,90  ) 'ND51 TIMESERIES SETTINGS (when -DBPCH_DIAG=y)'
       WRITE( 6,95  ) '---------------------------------------------'
       WRITE( 6,100 ) 'Turn on ND51 timeseries?    : ', Input_Opt%DO_ND51
       WRITE( 6,110 ) 'ND51 timeseries file name   : ',                      &
                        TRIM( Input_Opt%ND51_FILE )
       WRITE( 6,120 ) 'ND51 timeseries tracers     : ',                      &
                       ( Input_Opt%ND51_TRACERS(N), N=1, Input_Opt%N_ND51 )
       WRITE( 6,140 ) 'ND51 hour to write to disk  : ', Input_Opt%ND51_HR_WRITE
       WRITE( 6,140 ) 'ND51 averaging period [GMT] : ', Input_Opt%ND51_HR1,  &
                                                        Input_Opt%ND51_HR2
       WRITE( 6,130 ) 'ND51 longitude limits       : ', Input_Opt%ND51_IMIN, &
                                                        Input_Opt%ND51_IMAX
       WRITE( 6,130 ) 'ND51 latitude  limits       : ', Input_Opt%ND51_JMIN, &
                                                        Input_Opt%ND51_JMAX
       WRITE( 6,130 ) 'ND51 altitude  limits       : ', Input_Opt%ND51_LMIN, &
                                                        Input_Opt%ND51_LMAX
    ENDIF

    ! FORMAT statements
90  FORMAT( /, A     )
95  FORMAT( A        )
100 FORMAT( A, L5    )
110 FORMAT( A, A     )
120 FORMAT( A, 100I4 )
130 FORMAT( A, 2I5   )
140 FORMAT( A, 2F5.1 )

  END SUBROUTINE Config_Nd51
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: config_nd51b
!
! !DESCRIPTION: Copies ND51b satellite diagnostic information from the Config
!  object to Input_Opt, and does necessary checks.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Config_ND51b( Config, Input_Opt, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod, ONLY : OptInput
    USE RoundOff_Mod,  ONLY : Cast_and_RoundOff
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
    INTEGER                      :: N

    ! Arrays
    INTEGER                      :: a_int(QFYAML_MaxArr)

    ! Strings
    CHARACTER(LEN=255)           :: thisLoc
    CHARACTER(LEN=255)           :: errMsg
    CHARACTER(LEN=QFYAML_NamLen) :: key
    CHARACTER(LEN=QFYAML_StrLen) :: v_str

    ! String arrays
    CHARACTER(LEN=QFYAML_StrLen) :: a_str(QFYAML_MaxArr)

    !========================================================================
    ! Config_ND51b begins here!
    !========================================================================

    ! Initialize
    RC      = GC_SUCCESS
    errMsg  = ''
    thisLoc = ' -> at Config_ND51b (in module GeosCore/input_mod.F90)'

    !------------------------------------------------------------------------
    ! Turn on ND51b diagnostic
    !------------------------------------------------------------------------
    key    = "extra_diagnostics%legacy_bpch%ND51b_satellite%activate"
    v_bool = MISSING_BOOL
    CALL QFYAML_Add_Get( Config, key, v_bool, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Input_Opt%DO_ND51b = v_bool

    !------------------------------------------------------------------------
    ! Instantaneous 3-D timeseries file
    !------------------------------------------------------------------------
    key    = "extra_diagnostics%legacy_bpch%ND51b_satellite%output_file"
    v_str = MISSING_STR
    CALL QFYAML_Add_Get( Config, key, v_str, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Input_Opt%ND51b_FILE = TRIM( v_str )

    !------------------------------------------------------------------------
    ! Tracers to include
    !------------------------------------------------------------------------
    key   = "extra_diagnostics%legacy_bpch%ND51b_satellite%tracers"
    a_int = MISSING_INT
    CALL QFYAML_Add_Get( Config, key, a_int, "", RC, dynamic_size=.TRUE. )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    DO N = 1, SIZE( a_int )
       IF ( a_int(N) == MISSING_INT ) EXIT
       Input_Opt%ND51b_TRACERS(N) = a_int(N)
    ENDDO
    Input_Opt%N_ND51b = N - 1

    !------------------------------------------------------------------------
    ! NHMS_WRITE
    !------------------------------------------------------------------------
    key   = "extra_diagnostics%legacy_bpch%ND51b_satellite%UTC_hour_for_write"
    v_str = MISSING_STR
    CALL QFYAML_Add_Get( Config, key, v_str, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Input_Opt%ND51b_HR_WRITE = Cast_and_RoundOff( v_str, places=2 )

    ! Make sure ND51b_HR_WRITE is in the range 0-23.999 hrs
    Input_Opt%ND51b_HR_WRITE = MOD( Input_Opt%ND51b_HR_WRITE, 24.0_fp )

    !------------------------------------------------------------------------
    ! HR1, HR2
    !------------------------------------------------------------------------
    key = "extra_diagnostics%legacy_bpch%ND51b_satellite%averaging_period_in_LT"
    a_str = MISSING_STR
    CALL QFYAML_Add_Get( Config, key, a_str, "", RC, dynamic_size=.TRUE. )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Input_Opt%ND51b_HR1 = Cast_and_RoundOff( a_str(1), places=2 )
    Input_Opt%ND51b_HR2 = Cast_and_RoundOff( a_str(2), places=2 )

    !------------------------------------------------------------------------
    ! IMIN, IMAX
    !------------------------------------------------------------------------
    key = "extra_diagnostics%legacy_bpch%ND51b_satellite%IMIN_and_IMAX_of_region"
    a_int = MISSING_INT
    CALL QFYAML_Add_Get( Config, key, a_int, "", RC, dynamic_size=.TRUE. )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Input_Opt%ND51b_IMIN = a_int(1)
    Input_Opt%ND51b_IMAX = a_int(2)

    !------------------------------------------------------------------------
    ! JMIN, JMAX
    !------------------------------------------------------------------------
    key = "extra_diagnostics%legacy_bpch%ND51b_satellite%JMIN_and_JMAX_of_region"
    a_int = MISSING_INT
    CALL QFYAML_Add_Get( Config, key, a_int, "", RC, dynamic_size=.TRUE. )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Input_Opt%ND51b_JMIN = a_int(1)
    Input_Opt%ND51b_JMAX = a_int(2)

    !------------------------------------------------------------------------
    ! LMIN, LMAX
    !------------------------------------------------------------------------
    key = "extra_diagnostics%legacy_bpch%ND51b_satellite%LMIN_and_LMAX_of_region"
    a_int = MISSING_INT
    CALL QFYAML_Add_Get( Config, key, a_int, "", RC, dynamic_size=.TRUE. )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Input_Opt%ND51b_LMIN = a_int(1)
    Input_Opt%ND51b_LMAX = a_int(2)

    !========================================================================
    ! Print to screen
    !========================================================================
    IF ( Input_Opt%amIRoot ) THEN
       WRITE( 6,90  ) 'ND51b TIMESERIES SETTINGS (when -DBPCH_DIAG=y)'
       WRITE( 6,95  ) '----------------------------------------------'
       WRITE( 6,100 ) 'Turn on ND51b timeseries?   : ', Input_Opt%DO_ND51b
       WRITE( 6,110 ) 'ND51b timeseries file name  : ',                      &
                        TRIM( Input_Opt%ND51b_FILE )
       WRITE( 6,120 ) 'ND51b timeseries tracers    : ',                      &
                       ( Input_Opt%ND51b_TRACERS(N), N=1, Input_Opt%N_ND51b )
       WRITE( 6,140 ) 'ND51b hour to write to disk : ', Input_Opt%ND51b_HR_WRITE
       WRITE( 6,140 ) 'ND51b averaging period [GMT]: ', Input_Opt%ND51b_HR1, &
                                                        Input_Opt%ND51b_HR2
       WRITE( 6,130 ) 'ND51b longitude limits      : ', Input_Opt%ND51b_IMIN,&
                                                        Input_Opt%ND51b_IMAX
       WRITE( 6,130 ) 'ND51b latitude  limits      : ', Input_Opt%ND51b_JMIN,&
                                                        Input_Opt%ND51b_JMAX
       WRITE( 6,130 ) 'ND51b altitude  limits      : ', Input_Opt%ND51b_LMIN,&
                                                        Input_Opt%ND51b_LMAX
    ENDIF

    ! FORMAT statements
90  FORMAT( /, A     )
95  FORMAT( A        )
100 FORMAT( A, L5    )
110 FORMAT( A, A     )
120 FORMAT( A, 100I4 )
130 FORMAT( A, 2I5   )
140 FORMAT( A, 2F5.1 )

  END SUBROUTINE Config_Nd51b
#ifdef TOMAS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: config_bpch_prodloss

! !DESCRIPTION: Copies bpch prodloss information from the Config object
!  to Input_Opt, and does necessary checks.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Config_Bpch_ProdLoss( Config, Input_Opt, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE CMN_DIAG_MOD,     ONLY : ND65
    USE gckpp_Monitor,    ONLY : Fam_Names
    USE gckpp_Parameters, ONLY : nFam
    USE Input_Opt_Mod,    ONLY : OptInput
    USE QFYAML_Mod,       ONLY : QFYAML_t
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
    INTEGER                      :: v_int, F, N, n_Advect

    ! Strings
    CHARACTER(LEN=255)           :: thisLoc
    CHARACTER(LEN=512)           :: errMsg
    CHARACTER(LEN=QFYAML_NamLen) :: key
    CHARACTER(LEN=QFYAML_StrLen) :: v_str

    ! Arrays
!    CHARACTER(LEN=255) :: SUBSTRS(MAXDIM)

    !=========================================================================
    ! Config_Bpch_ProdLoss begins here!
    !=========================================================================

    ! Initialize
    RC      = GC_SUCCESS
    errMsg  = 'Error reading the "geoschem_config.yml" file!'
    thisLoc = ' -> at Config_Bpch_ProdLoss (in module GeosCore/input_mod.F90)'

    !------------------------------------------------------------------------
    ! Activate bpch ND65 diagnostic?
    !------------------------------------------------------------------------
    key    = "extra_diagnostics%legacy_bpch%bpch_diagnostics%"            // &
             "ND65_prodloss%activate"
    v_bool = MISSING_BOOL
    CALL QFYAML_Add_Get( Config, key, v_bool, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Input_Opt%DO_SAVE_PL = v_bool

    !------------------------------------------------------------------------
    ! Number of levels for ND65 bpch diagnostic
    !------------------------------------------------------------------------
    key    = "extra_diagnostics%legacy_bpch%bpch_diagnostics%"            // &
             "ND65_prodloss%number_of_levels"
    v_int = MISSING_INT
    CALL QFYAML_Add_Get( Config, key, v_int, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Input_Opt%ND65 = v_int

    ! Copy field to variable in CMN_DIAG
    ND65 = Input_Opt%ND65

    !=================================================================
    ! Error check families for certain types of simulations
    !=================================================================

    ! Offline aerosol -- turn off DO_SAVE_PL
    IF ( Input_Opt%ITS_AN_AEROSOL_SIM ) THEN
       Input_Opt%DO_SAVE_PL    = .FALSE.
       Input_Opt%ND65          = 0
    ENDIF

    !=================================================================
    ! Set fields of Input Options object
    !=================================================================

    ! Number of advected species
    N_ADVECT = Input_Opt%N_ADVECT

    IF ( Input_Opt%DO_SAVE_PL ) THEN
       IF ( Input_Opt%ITS_A_FULLCHEM_SIM ) THEN
          ! Fullchem - Obtain NFAM from KPP
          Input_Opt%NFAM = NFAM
       ELSEIF ( Input_Opt%ITS_A_TAGO3_SIM ) THEN
          ! Tagged O3
          Input_Opt%NFAM = 2*N_ADVECT
       ELSEIF ( Input_Opt%ITS_A_TAGCO_SIM ) THEN
          ! Tagged CO
          IF ( Input_Opt%LPCO_NMVOC ) THEN
             Input_Opt%NFAM = N_ADVECT+2
          ELSE
             Input_Opt%NFAM = N_ADVECT+6
          ENDIF
       ENDIF
    ENDIF

    ! Return if there are no prod/loss families
    ! or if we have turned off this diagnostic
    IF ( .not. ( Input_Opt%DO_SAVE_PL .and. Input_Opt%NFAM > 0 )) THEN
       Input_Opt%DO_SAVE_PL = .FALSE.
       Input_Opt%ND65       = 0
    ENDIF

    ! Loop over families
    DO F = 1, Input_Opt%NFAM

       IF ( Input_Opt%ITS_A_FULLCHEM_SIM ) THEN

          ! Fullchem - Obtain FAM_NAME from KPP
          Input_Opt%FAM_NAME(F) = FAM_NAMES(F)

       ELSEIF ( Input_Opt%ITS_A_TAGO3_SIM ) THEN

          ! Tagged O3
          IF ( F <= N_ADVECT ) THEN
             Input_Opt%FAM_NAME(F) = &
                  'P' // TRIM(Input_Opt%AdvectSpc_Name(F))
          ELSE
             Input_Opt%FAM_NAME(F) = &
                  'L' // TRIM(Input_Opt%AdvectSpc_Name(F-N_ADVECT))
          ENDIF

       ELSEIF ( Input_Opt%ITS_A_TAGCO_SIM ) THEN

          ! Tagged CO
          IF ( F <= N_ADVECT ) THEN
             Input_Opt%FAM_NAME(F) = 'L'//Input_Opt%AdvectSpc_Name(F)
          ELSEIF ( F == N_ADVECT+1 ) THEN
             Input_Opt%FAM_NAME(F) = 'PCO_CH4'
          ELSEIF ( F == N_ADVECT+2 ) THEN
             Input_Opt%FAM_NAME(F) = 'PCO_NMVOC'
          ELSEIF ( F == N_ADVECT+3 ) THEN
             Input_Opt%FAM_NAME(F) = 'PCO_ISOP'
          ELSEIF ( F == N_ADVECT+4 ) THEN
             Input_Opt%FAM_NAME(F) = 'PCO_CH3OH'
          ELSEIF ( F == N_ADVECT+5 ) THEN
             Input_Opt%FAM_NAME(F) = 'PCO_MONO'
          ELSEIF ( F == N_ADVECT+6 ) THEN
             Input_Opt%FAM_NAME(F) = 'PCO_ACET'
          ENDIF

       ENDIF

       ! Get family type as prod or loss
       IF ( Input_Opt%FAM_NAME(F)(1:1) == 'P'   .or. &
            Input_Opt%FAM_NAME(F)(1:1) == 'p' ) THEN
          Input_Opt%FAM_TYPE(F) = 'prod'
       ELSE
          Input_Opt%FAM_TYPE(F) = 'loss'
       ENDIF

    ENDDO

    !=================================================================
    ! Print to screen
    !=================================================================
    IF ( Input_Opt%amIRoot ) THEN
       WRITE( 6, '(/,a)' ) 'PROD & LOSS DIAGNOSTIC SETTINGS'
       WRITE( 6, '(  a)' ) '-------------------------------'
       WRITE( 6, 100 ) 'Turn on prod & loss diag?   : ', &
                        Input_Opt%DO_SAVE_PL
       WRITE( 6, 110 ) '# of levels for P/L diag    : ', &
                        Input_Opt%ND65

       ! Loop over families
       DO F = 1, Input_Opt%NFAM

          ! Write family name and type
          WRITE( 6, 120 ) TRIM(Input_Opt%FAM_NAME(F)), &
                          TRIM(Input_Opt%FAM_TYPE(F))

       ENDDO

    ENDIF

    ! FORMAT statements
90  FORMAT( /, A                             )
95  FORMAT( A                                )
100 FORMAT( A, L5                            )
110 FORMAT( A, I5                            )
120 FORMAT( /, 'Family=', A10, '  Type=', A4 )

  END SUBROUTINE Config_Bpch_ProdLoss
!EOC
#endif
#endif
#endif
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: config_Hg
!
! !DESCRIPTION: Copies Hg simulation information from the Config object
!  to Input_Opt, and does necessary checks.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Config_Hg( Config, Input_Opt, RC )
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

    !=================================================================
    ! Config_Hg begins here!
    !=================================================================

    ! Initialize
    RC      = GC_SUCCESS
    errMsg  = ''
    thisLoc = ' -> at Config_Hg (in module GeosCore/input_mod.F90)'

    !------------------------------------------------------------------------
    ! Use dynamic ocean Hg?
    !------------------------------------------------------------------------
    key    = "Hg_simulation_options%sources%use_dynamic_ocean_Hg"
    v_bool = MISSING_BOOL
    CALL QFYAML_Add_Get( Config, key, v_bool, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Input_Opt%LDYNOCEAN = v_bool

    !------------------------------------------------------------------------
    ! Use preindustrial Hg?
    !------------------------------------------------------------------------
    key    = "Hg_simulation_options%sources%use_preindustrial_Hg"
    v_bool = MISSING_BOOL
    CALL QFYAML_Add_Get( Config, key, v_bool, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Input_Opt%LPREINDHG = v_bool

    !------------------------------------------------------------------------
    ! Use arctic river Hg?
    !------------------------------------------------------------------------
    key    = "Hg_simulation_options%sources%use_arctic_river_Hg"
    v_bool = MISSING_BOOL
    CALL QFYAML_Add_Get( Config, key, v_bool, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Input_Opt%LARCTICRIV = v_bool

    !------------------------------------------------------------------------
    ! Tie Hg2(aq) reduction to UV-B radiation?
    !------------------------------------------------------------------------
    key    = "Hg_simulation_options%chemistry%tie_HgIIaq_reduction_to_UVB"
    v_bool = MISSING_BOOL
    CALL QFYAML_Add_Get( Config, key, v_bool, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Input_Opt%LKRedUV = v_bool

    !------------------------------------------------------------------------
    ! Use GTMM soil model
    !
    ! NOTE: As of April 2022, GTMM is broken.  We look to the community
    ! to take the lead in restoring it.  Until that happens, these options
    ! will have no effect. -- Bob Yantosca (04 Apr 2022)
    !------------------------------------------------------------------------
    key    = "Hg_simulation_options%GTMM_soil_model%activate"
    v_bool = MISSING_BOOL
    CALL QFYAML_Add_Get( Config, key, v_bool, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Input_Opt%LGTMM = v_bool

    !------------------------------------------------------------------------
    ! GTMM restart file name
    !
    ! NOTE: As of April 2022, GTMM is broken.  We look to the community
    ! to take the lead in restoring it.  Until that happens, these options
    ! will have no effect. -- Bob Yantosca (04 Apr 2022)
    !------------------------------------------------------------------------
    key   = "Hg_simulation_options%GTMM_soil_model%restart_file"
    v_str = MISSING_STR
    CALL QFYAML_Add_Get( Config, key, v_str, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Input_Opt%GTMM_RST_FILE = TRIM( v_str )

    !------------------------------------------------------------------------
    ! Sanity checks
    !------------------------------------------------------------------------
    IF ( .not. Input_Opt%ITS_A_MERCURY_SIM ) THEN
       Input_Opt%LGTMM      = .FALSE.
       Input_Opt%LDYNOCEAN  = .FALSE.
       Input_Opt%LARCTICRIV = .FALSE.
       Input_Opt%LKRedUV    = .FALSE.
    ENDIF

    !========================================================================
    ! Print to screen
    !========================================================================
    IF ( Input_Opt%amIRoot ) THEN
       WRITE( 6, 90  ) 'MERCURY SIMULATION SETTINGS'
       WRITE( 6, 95  ) '---------------------------'
       WRITE( 6, 110 ) 'Use dynamic ocean Hg model? : ', Input_Opt%LDYNOCEAN
       WRITE( 6, 110 ) 'Preindustrial simulation?   : ', Input_Opt%LPREINDHG
       WRITE( 6, 110 ) 'Use Arctic river Hg ?       : ', Input_Opt%LARCTICRIV
       WRITE( 6, 110 ) 'Tie HgII(aq) red. to UV-B?  : ', Input_Opt%LKRedUV
       WRITE( 6, 110 ) 'Use GTMM ?                  : ', Input_Opt%LGTMM
       WRITE( 6, 120 ) '=> GTMM restart file        : ',                      &
                       TRIM( Input_Opt%GTMM_RST_FILE )
    ENDIF

    ! FORMAT statements
90  FORMAT( /, A  )
95  FORMAT( A     )
100 FORMAT( A, I4 )
110 FORMAT( A, L5 )
120 FORMAT( A, A  )

  END SUBROUTINE Config_Hg
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: config_ch4
!
! !DESCRIPTION: Copies CH4 simulation information from the Config
!  object to Input_Opt, and does necessary checks.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Config_CH4( Config, Input_Opt, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod, ONLY : OptInput
    USE RoundOff_Mod,  ONLY : Cast_and_RoundOff
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(QFYAML_t), INTENT(INOUT) :: Config      ! YAML Config object
    TYPE(OptInput), INTENT(INOUT) :: Input_Opt   ! Input options
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:

    ! Scalars
    INTEGER                      :: N
    INTEGER                      :: v_int
    LOGICAL                      :: v_bool

    ! Strings
    CHARACTER(LEN=255)           :: thisLoc
    CHARACTER(LEN=255)           :: errMsg
    CHARACTER(LEN=QFYAML_NamLen) :: key
    CHARACTER(LEN=QFYAML_StrLen) :: v_str

    !========================================================================
    ! Config_CH4 begins here!
    !========================================================================

    ! Initialize
    RC      = GC_SUCCESS
    errMsg  = 'Error reading the "geoschem_config.yml" file!'
    thisLoc = ' -> at Config_CH4 (in module GeosCore/input_mod.F90)'

    !------------------------------------------------------------------------
    ! Use AIRS observational operator?
    !------------------------------------------------------------------------
    key    = "CH4_simulation_options%use_observational_operators%AIRS"
    v_bool = MISSING_BOOL
    CALL QFYAML_Add_Get( Config, TRIM( key ), v_bool, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Input_Opt%AIRS_CH4_OBS = v_bool

    !------------------------------------------------------------------------
    ! Use GOSAT observational operator?
    !------------------------------------------------------------------------
    key    = "CH4_simulation_options%use_observational_operators%GOSAT"
    v_bool = MISSING_BOOL
    CALL QFYAML_Add_Get( Config, TRIM( key ), v_bool, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Input_Opt%GOSAT_CH4_OBS =  v_bool

    !------------------------------------------------------------------------
    ! Use TCCON observational operator?
    !------------------------------------------------------------------------
    key    = "CH4_simulation_options%use_observational_operators%TCCON"
    v_bool = MISSING_BOOL
    CALL QFYAML_Add_Get( Config, TRIM( key ), v_bool, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Input_Opt%TCCON_CH4_OBS = v_bool

    !------------------------------------------------------------------------
    ! Do an analytical inversion?
    !------------------------------------------------------------------------
    key    = "CH4_simulation_options%analytical_inversion%activate"
    v_bool = MISSING_BOOL
    CALL QFYAML_Add_Get( Config, TRIM( key ), v_bool, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Input_Opt%AnalyticalInv = v_bool

    !------------------------------------------------------------------------
    ! Emission perturbation
    !------------------------------------------------------------------------
    key   = "CH4_simulation_options%analytical_inversion%emission_perturbation"
    v_str = MISSING_STR
    CALL QFYAML_Add_Get( Config, TRIM( key ), v_str, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Input_Opt%PerturbEmis = Cast_and_RoundOff( v_str, places=4 )

    !------------------------------------------------------------------------
    ! Current state vector element number
    !------------------------------------------------------------------------
    key   = &
     "CH4_simulation_options%analytical_inversion%state_vector_element_number"
    v_int = MISSING_INT
    CALL QFYAML_Add_Get( Config, TRIM( key ), v_int, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Input_Opt%StateVectorElement = v_int

    !------------------------------------------------------------------------
    ! Use emission scale factor?
    !------------------------------------------------------------------------
    key = &
     "CH4_simulation_options%analytical_inversion%use_emission_scale_factor"
    v_bool = MISSING_BOOL
    CALL QFYAML_Add_Get( Config, TRIM( key ), v_bool, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Input_Opt%UseEmisSF = v_bool

    !------------------------------------------------------------------------
    ! Use OH scale factors?
    !------------------------------------------------------------------------
    key    = "CH4_simulation_options%analytical_inversion%use_OH_scale_factors"
    v_bool = MISSING_BOOL
    CALL QFYAML_Add_Get( Config, TRIM( key ), v_bool, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Input_Opt%UseOHSF = v_bool

    !========================================================================
    ! Print to screen
    !========================================================================
    IF ( Input_Opt%amIRoot ) THEN
       WRITE(6,90 ) 'CH4 SIMULATION SETTINGS'
       WRITE(6,95 ) '-----------------------'
       WRITE(6,100) 'Use AIRS obs operator    : ', Input_Opt%AIRS_CH4_OBS
       WRITE(6,100) 'Use GOSAT obs operator   : ', Input_Opt%GOSAT_CH4_OBS
       WRITE(6,100) 'Use TCCON obs operator   : ', Input_Opt%TCCON_CH4_OBS
       WRITE(6,100) 'Do analytical inversion  : ', Input_Opt%AnalyticalInv
       WRITE(6,110) 'Emission perturbation    : ', Input_Opt%PerturbEmis
       WRITE(6,120) 'Current state vector elem: ', Input_Opt%StateVectorElement
       WRITE(6,100) 'Use emis scale factors   : ', Input_Opt%UseEmisSF
       WRITE(6,100) 'Use OH scale factors     : ', Input_Opt%UseOHSF
    ENDIF

    ! FORMAT statements
90  FORMAT( /, A    )
95  FORMAT( A       )
100 FORMAT( A, L5   )
110 FORMAT( A, f6.2 )
120 FORMAT( A, I5   )

  END SUBROUTINE Config_CH4
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: config_pops
!
! !DESCRIPTION: Copies POPs simulation information from the Config
!  object to Input_Opt, and does necessary checks.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Config_POPs( Config, Input_Opt, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod, ONLY : OptInput
    USE RoundOff_Mod,  ONLY : Cast_and_RoundOff
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(QFYAML_t), INTENT(INOUT) :: Config      ! YAML Config object
    TYPE(OptInput), INTENT(INOUT) :: Input_Opt   ! Input options
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    LOGICAL                      :: v_bool

    ! Strings
    CHARACTER(LEN=255)           :: thisLoc
    CHARACTER(LEN=255)           :: errMsg
    CHARACTER(LEN=QFYAML_NamLen) :: key
    CHARACTER(LEN=QFYAML_StrLen) :: v_str

    !========================================================================
    ! Config_POPs begins here!
    !========================================================================

    ! Initialize
    RC      = GC_SUCCESS

    ! Continue initializing
    errMsg  = ''
    thisLoc = ' -> at Config_POPs (in module GeosCore/input_mod.F90)'

    !------------------------------------------------------------------------
    ! POP species
    !------------------------------------------------------------------------
    key   = "POPs_simulation_options%POP_type"
    v_str = MISSING_STR
    CALL QFYAML_Add_Get( Config, key, v_str, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Input_Opt%POP_TYPE = TRIM( v_str )

    !------------------------------------------------------------------------
    ! Dummy for future process logical switches
    !------------------------------------------------------------------------
    key    = "POPs_simulation_options%chemistry_processing"
    v_bool = MISSING_BOOL
    CALL QFYAML_Add_Get( Config, key, v_bool, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Input_Opt%CHEM_PROCESS = v_bool

    !------------------------------------------------------------------------
    ! Molecular weight
    !------------------------------------------------------------------------
    key   = "POPs_simulation_options%POP_XMW"
    v_str = MISSING_STR
    CALL QFYAML_Add_Get( Config, key, v_str, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Input_Opt%POP_XMW = Cast_and_RoundOff( v_str, places=0 )

    !------------------------------------------------------------------------
    ! KOA
    !------------------------------------------------------------------------
    key   = "POPs_simulation_options%POP_KOA"
    v_str = MISSING_STR
    CALL QFYAML_Add_Get( Config, key, v_str, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Input_Opt%POP_KOA = Cast_and_RoundOff( v_str, places=0 )

    !------------------------------------------------------------------------
    ! KBC
    !------------------------------------------------------------------------
    key   = "POPs_simulation_options%POP_KBC"
    v_str = MISSING_STR
    CALL QFYAML_Add_Get( Config, key, v_str, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Input_Opt%POP_KBC = Cast_and_RoundOff( v_str, places=0 )

    !------------------------------------------------------------------------
    ! OH oxidation
    !------------------------------------------------------------------------
    key   = "POPs_simulation_options%POP_K_POPG_OH"
    v_str = MISSING_STR
    CALL QFYAML_Add_Get( Config, key, v_str, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Input_Opt%POP_K_POPG_OH = Cast_and_RoundOff( v_str, places=0 )

    !------------------------------------------------------------------------
    ! O3 oxidation 1
    !------------------------------------------------------------------------
    key   = "POPs_simulation_options%POP_K_POPP_O3A"
    v_str = MISSING_STR
    CALL QFYAML_Add_Get( Config, key, v_str, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Input_Opt%POP_K_POPP_O3A = Cast_and_RoundOff( v_str, places=0 )

    !------------------------------------------------------------------------
    ! O3 oxidation 2
    !------------------------------------------------------------------------
    key   = "POPs_simulation_options%POP_K_POPP_O3B"
    v_str = MISSING_STR
    CALL QFYAML_Add_Get( Config, key, v_str, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Input_Opt%POP_K_POPP_O3B = Cast_and_RoundOff( v_str, places=0 )

    !------------------------------------------------------------------------
    ! H*
    !------------------------------------------------------------------------
    key   = "POPs_simulation_options%POP_HSTAR"
    v_str = MISSING_STR
    CALL QFYAML_Add_Get( Config, key, v_str, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Input_Opt%POP_HSTAR = Cast_and_RoundOff( v_str, places=0 )

    !------------------------------------------------------------------------
    ! DEL_H
    !------------------------------------------------------------------------
    key   = "POPs_simulation_options%POP_DEL_H"
    v_str = MISSING_STR
    CALL QFYAML_Add_Get( Config, key, v_str, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Input_Opt%POP_DEL_H = Cast_and_RoundOff( v_str, places=0 )

    !------------------------------------------------------------------------
    ! DEL_Hw
    !------------------------------------------------------------------------
    key   = "POPs_simulation_options%POP_DEL_Hw"
    v_str = MISSING_STR
    CALL QFYAML_Add_Get( Config, key, v_str, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Input_Opt%POP_DEL_Hw = Cast_and_RoundOff( v_str, places=0 )

    !=================================================================
    ! Print to screen
    !=================================================================
    IF ( Input_Opt%amIRoot ) THEN
       WRITE( 6, 90  ) 'POPs SIMULATION SETTINGS'
       WRITE( 6, 95  ) '------------------------'
       WRITE( 6, 120 ) 'Species of POP        : ', Input_Opt%POP_TYPE
       WRITE( 6, 110 ) 'Chemistry on?         : ', Input_Opt%CHEM_PROCESS
       WRITE( 6, 130 ) 'POP_XMW               : ', Input_Opt%POP_XMW
       WRITE( 6, 130 ) 'POP_KOA               : ', Input_Opt%POP_KOA
       WRITE( 6, 130 ) 'POP_KBC               : ', Input_Opt%POP_KBC
       WRITE( 6, 130 ) 'POP_K_POPG_OH         : ', Input_Opt%POP_K_POPG_OH
       WRITE( 6, 130 ) 'POP_K_POPP_O3A        : ', Input_Opt%POP_K_POPP_O3A
       WRITE( 6, 130 ) 'POP_K_POPP_O3B        : ', Input_Opt%POP_K_POPP_O3B
       WRITE( 6, 130 ) 'POP_HSTAR             : ', Input_Opt%POP_HSTAR
       WRITE( 6, 130 ) 'POP_DEL_H             : ', Input_Opt%POP_DEL_H
       WRITE( 6, 130 ) 'POP_DEL_Hw            : ', Input_Opt%POP_DEL_Hw
    ENDIF

    ! FORMAT statements
90  FORMAT( /, A      )
95  FORMAT( A         )
110 FORMAT( A, L5     )
120 FORMAT( A, A      )
130 FORMAT( A, ES10.2 )

    ! Return success
    RC = GC_SUCCESS

  END SUBROUTINE Config_POPs
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: config_passivespecies
!
! !DESCRIPTION: Copies passive species information from the Config
!  object to Input_Opt, and does necessary checks.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Config_PassiveSpecies( Config, Input_Opt, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod, ONLY : OptInput
    USE RoundOff_Mod,  ONLY : Cast_and_RoundOff
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(QFYAML_t), INTENT(INOUT) :: Config      ! YAML Config object
    TYPE(OptInput), INTENT(INOUT) :: Input_Opt   ! Input Options object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REVISION HISTORY:
!  04 Sep 2015 - C. Keller   - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER                      :: C
    INTEGER                      :: P
    REAL(yp)                     :: v_real

    ! Strings
    CHARACTER(LEN=255)           :: thisLoc
    CHARACTER(LEN=512)           :: errMsg
    CHARACTER(LEN=QFYAML_NamLen) :: key
    CHARACTER(LEN=QFYAML_StrLen) :: v_str

    ! String arrays
    CHARACTER(LEN=QFYAML_NamLen) :: passSpcName(Input_Opt%Max_PassiveSpc)

    !========================================================================
    ! Config_PassiveSpecies
    !========================================================================

    ! Initialize
    RC      = GC_SUCCESS
    errMsg  = ''
    thisLoc = &
     ' -> at Config_PassiveSpecies (in module GeosCore/input_mod.F90)'

    !========================================================================
    ! Check for passive species
    !========================================================================
    key = "operations%transport%passive_species%"
    CALL QFYAML_FindNextHigher( Config,             TRIM( key ),             &
                                Input_Opt%nPassive, passSpcName             )

    ! Return if there are no passive species
    IF ( Input_Opt%nPassive == 0 ) RETURN

    ! Throw an error if there are more than the max # of passive species
    IF ( Input_Opt%nPassive > Input_Opt%Max_PassiveSpc ) THEN
       errMsg = 'Number of passive species exceeds maximum. ' // &
            'This value can be modified in input_opt_mod.F90.'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !========================================================================
    ! Copy the passive species names & metadata into Input_Opt
    !========================================================================

    ! Write header
    IF ( Input_Opt%amIRoot ) THEN
       WRITE( 6, '(/,a)' ) 'PASSIVE SPECIES SETTINGS'
       WRITE( 6, '(  a)' ) '------------------------'
    ENDIF

    ! Loop over passive species
    DO P = 1, Input_Opt%nPassive

       !---------------------------------------------------------------------
       ! Passive species index and name
       !---------------------------------------------------------------------

       ! Index
       Input_Opt%Passive_Id(P)   = P

       ! NOTE: passSpcName is the full YAML variable name
       ! (e.g. operations%trasnport%passive_species%...",
       ! but we'll extract the last part of that for saving into
       ! Input_Opt%Passive_Name.
       C = INDEX( passSpcName(P), '%', back=.TRUE. )
       Input_Opt%Passive_Name(P) = TRIM( passSpcName(P)(C+1:) )

       !---------------------------------------------------------------------
       ! Passive species long name (if found)
       !---------------------------------------------------------------------
       key = TRIM( passSpcName(P) ) // '%long_name'
       v_str = MISSING_STR
       CALL QFYAML_Add_Get( Config, TRIM( key ), v_str, "", RC )
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = 'Error parsing ' // TRIM( key ) // '!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
       Input_Opt%Passive_LongName(P) = TRIM( v_str )

       !---------------------------------------------------------------------
       ! Passive species molecular weight (g)
       !---------------------------------------------------------------------
       key = TRIM( passSpcName(P) ) // '%mol_wt_in_g'
       v_str = MISSING_STR
       CALL QFYAML_Add_Get( Config, TRIM( key ), v_str, "", RC )
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = 'Error parsing ' // TRIM( key ) // '!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
       Input_Opt%Passive_MW(P) = Cast_and_RoundOff( v_str, places=2 )

       !---------------------------------------------------------------------
       ! Passive species lifetime (s); -1 means it never decays
       !---------------------------------------------------------------------
       key = TRIM( passSpcName(P) ) // '%lifetime_in_s'
       CALL QFYAML_Add_Get( Config, TRIM( key ), v_str, "", RC )
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = 'Error parsing ' // TRIM( key ) // '!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
       Input_Opt%Passive_Tau(P) = Cast_and_RoundOff( v_str, places=-1 )

       ! Determine if the passive species decays (i.e. if it has
       ! an atmospheric lifetime that is not -1).  This will allow
       ! us to skip those passive species that do not decay in
       ! routine CHEM_PASSIVE_SPECIES, to speed up execution.
       IF ( Input_Opt%Passive_Tau(P) > 0.0_fp ) THEN
          Input_Opt%nPassive_Decay = Input_Opt%nPassive_Decay + 1
          Input_Opt%Passive_DecayId(Input_Opt%nPassive_Decay) = P
       ENDIF

       !---------------------------------------------------------------------
       ! Default background concentration (v/v)
       !---------------------------------------------------------------------
       key = TRIM( passSpcName(P) ) // '%default_bkg_conc_in_vv'
       CALL QFYAML_Add_Get( Config, TRIM( key ), v_str, "", RC )
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = 'Error parsing ' // TRIM( key ) // '!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
       Input_Opt%Passive_Initconc(P) = Cast_and_RoundOff( v_str, places=-1 )

       !---------------------------------------------------------------------
       ! Logfile output
       !---------------------------------------------------------------------
       IF ( Input_Opt%amIRoot ) THEN
          WRITE( 6, 95  ) 'Added passive species: '
          WRITE( 6, 110 ) ' - Species name                : ',               &
                           TRIM( Input_Opt%Passive_Name(P) )
          WRITE( 6, 120 ) ' - Molec. weight [g/mol]       : ',               &
                          Input_Opt%Passive_MW(P)
          WRITE( 6, 130 ) ' - Lifetime [s]                : ',               &
                          Input_Opt%Passive_TAU(P)
          WRITE( 6, 130 ) ' - Initial concentration [v/v] : ',               &
                          Input_Opt%Passive_InitConc(P)
          WRITE( 6, 110 ) ' - Species long name           : ',               &
                          TRIM( Input_Opt%Passive_LongName(P) )
       ENDIF
    ENDDO

    ! Format statememnts
90  FORMAT( /, A      )
95  FORMAT( A         )
110 FORMAT( A, A      )
120 FORMAT( A, F10.2  )
130 FORMAT( A, ES10.2 )

  END SUBROUTINE Config_PassiveSpecies
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: validate_directories
!
! !DESCRIPTION: Makes sure that each of the directories that we have read from
!  the GEOS-Chem input file are valid. Also, trailing separator characters will
!  be added.
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
    ! Validate_Directories begins here!
    !=================================================================

    ! Initialize
    RC      = GC_SUCCESS
    errMsg  = 'Invalid directory encountered!'
    thisLoc = ' -> at Validate_Directories (in module GeosCore/input_mod.F90)'

    ! Skip for dry-runs
    IF ( Input_Opt%DryRun ) RETURN

#if !defined( MODEL_CESM )
    ! Check directories
    CALL Check_Directory( Input_Opt, Input_Opt%DATA_DIR, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
#endif

    CALL Check_Directory( Input_Opt, Input_Opt%CHEM_INPUTS_DIR, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

#if !defined( MODEL_CESM )
    CALL Check_Directory( Input_Opt, Input_Opt%RUN_DIR, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
#endif

  END SUBROUTINE Validate_Directories
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: check_directory
!
! !DESCRIPTION: Makes sure that the given directory is valid.  Also a trailing
!  slash character will be added if necessary.
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
    ! Check_Directory begins here!
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
! !DESCRIPTION: Computes the smallest dynamic time step for the model, based on
!  which operation are turned on.  This is called from routine
!  Read_Input_File, after all of the timesteps and logical flags have been
!  read from the configuration file.
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
    ! Check_Time_Steps begins here!
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
! !DESCRIPTION: Tests to see if there is output scheduled on the last day of
!  the run.
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
    ! Is_Last_Day_Good begins here!
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
    ! NAP emissions are set in HEMCO_Config.rc
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
        MAX( Ind_('ASOAN','A'), 0 ) + &
        MAX( Ind_('ASOG1','A'), 0 ) + &
        MAX( Ind_('ASOG2','A'), 0 ) + &
        MAX( Ind_('ASOG3','A'), 0 ) + &
        MAX( Ind_('TSOG0','A'), 0 ) + &
        MAX( Ind_('TSOG1','A'), 0 ) + &
        MAX( Ind_('TSOG2','A'), 0 ) + &
        MAX( Ind_('TSOG3','A'), 0 )

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
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Find_Number_of_Species
!
! !DESCRIPTION: Searches a string array containing species names and returns
!  the number of valid species (i.e. species that do not match MISSING_STR).
!  Assumes all the valid species will be listed contiguously at the front
!  of the array
!\\
!\\
! !INTERFACE:
!
   FUNCTION Get_Number_of_Species( a_str ) RESULT( n_valid )
!
! !INPUT PARAMETERS:
!
    CHARACTER(LEN=*), INTENT(IN) :: a_str(:)
!
! !RETURN VALUE:
!
    INTEGER                      :: n_valid
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: N

    ! Return the number of valid species
    n_valid = 0
    DO N = 1, SIZE( a_str )
       IF ( TRIM( a_str(N) ) == MISSING_STR ) EXIT
       n_valid = n_valid + 1
    ENDDO

  END FUNCTION Get_Number_of_Species
!EOC
END MODULE Input_Mod
