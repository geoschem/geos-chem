!BOC
#if defined( ESMF_ ) || defined( EXTERNAL_GRID ) || defined( MODEL_ )
!-----------------------------------------------------------------
!         %%%%%%% GEOS-Chem HP (with ESMF & MPI) %%%%%%%
!        %%%% GEOS-Chem Coupled with External Models %%%%
!
! When GEOS-Chem is connected to an external model or in GCHP,
! the GEOS-Chem classic main.F90 should not be built.
!-----------------------------------------------------------------
#else
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: main.F90
!
! !DESCRIPTION: Program GEOS\_CHEM is the main level driver program for the
!  GEOS-Chem model of atmospheric chemistry and composition.
!\\
!\\
! !INTERFACE:
!
PROGRAM GEOS_Chem
!
! !USES:
!
  !-----------------------------------------------------------------
  ! Parameters to define floating-point variables
  !-----------------------------------------------------------------
  USE PRECISION_MOD,   ONLY : fpp => fp  ! Flexible precision
  USE PRECISION_MOD,   ONLY : f4         ! 4-byte floating point
  USE PRECISION_MOD,   ONLY : f8         ! 8-byte floating point

  !-----------------------------------------------------------------
  ! Basic GEOS-Chem modules
  !-----------------------------------------------------------------
  USE DiagList_Mod          ! Derived type for diagnostics list
  USE TaggedDiagList_Mod    ! Derived type for tagged diagnostics list
  USE Diagnostics_Mod       ! Set select netcdf diagnostics
  USE ErrCode_Mod           ! Error codes for success or failure
  USE ERROR_MOD             ! For error checking
  USE FILE_MOD              ! For file I/O
  USE GC_Environment_Mod    ! For allocating derived type objects
  USE GC_GRID_MOD           ! For defining the lons/lats/areas of the grid
  USE Input_Opt_Mod         ! Derived type for Input Options
  USE INPUT_MOD             ! For reading settings from "input.geos"
  USE OLSON_LANDMAP_MOD     ! Computes IREG, ILAND, IUSE from Olson map
  USE PhysConstants         ! Physical constants
  USE PRESSURE_MOD          ! For computing pressure at grid boxes
  USE Grid_Registry_Mod     ! Registers horizontal/vertical grid metadata
  USE State_Chm_Mod         ! Derived type for Chemistry State object
  USE State_Diag_Mod        ! Derived type for Diagnostics State object
  USE State_Grid_Mod        ! Derived type for Grid State object
  USE State_Met_Mod         ! Derived type for Meteorology State object
  USE TIME_MOD              ! For computing date & time
  USE TIMERS_MOD            ! For GEOS-Chem timers (optional)
  USE UnitConv_Mod          ! For species conc unit conversions

  !-----------------------------------------------------------------
  ! GEOS-Chem chemistry modules
  !-----------------------------------------------------------------
  USE AEROSOL_MOD,     ONLY : Set_AerMass_Diagnostic
  USE CARBON_MOD            ! For SOA simulation
  USE CHEMISTRY_MOD         ! Driver routines for chemistry
  USE MERCURY_MOD           ! For offline Hg simulation (driver)
  USE OCEAN_MERCURY_MOD     ! For offline Hg simulation (ocean model)
  USE STRAT_CHEM_MOD        ! For linearized stratospheric chemistry
  USE TOMS_MOD              ! For overhead O3 columns (for FAST-J)
  USE UCX_MOD               ! For unified trop-strat chemistry (SDE)
  USE UVALBEDO_MOD          ! For reading UV albedoes (for FAST-J)
  USE SET_GLOBAL_CH4_MOD    ! For setting global CH4 concentrations

  !-----------------------------------------------------------------
  ! GEOS-Chem deposition modules
  !-----------------------------------------------------------------
  USE DEPO_MERCURY_MOD      ! Deposition for offline Hg simulation
  USE DRYDEP_MOD            ! For dry deposition
  USE WETSCAV_MOD           ! For wet deposition (rainout & washout)

  !-----------------------------------------------------------------
  ! GEOS-Chem diagnostics modules
  !-----------------------------------------------------------------
#ifdef BPCH_DIAG
  USE DIAG_MOD              ! G-C diagnostic arrays & counters
  USE CMN_DIAG_MOD          ! Logical switches for G-C diagnostics
  USE DIAG51_MOD            ! For ND51  (satellite timeseries) diag
  USE DIAG51b_MOD           ! For ND51b (satellite timeseries) diag
#endif
  USE PLANEFLIGHT_MOD       ! For planeflight track diag
  USE DIAG_OH_MOD           ! For OH,HO2,etc. prod diag
  USE HISTORY_MOD           ! Updated netCDF diagnostics
  USE OBSPACK_MOD           ! For ObsPack diagnostics
  USE GOSAT_CH4_MOD         ! For GOSAT observation operator
  USE TCCON_CH4_MOD         ! For TCCON observation operator
  USE HCOI_GC_MAIN_MOD      ! Writes out HEMCO diagnostics (C. Keller)

  !-----------------------------------------------------------------
  ! GEOS-Chem dynamics modules
  !-----------------------------------------------------------------
  USE CONVECTION_MOD        ! For deep cloud convection
  USE LINOZ_MOD             ! For LINOX linear strat chemistry
  USE PBL_MIX_MOD           ! To compute PBL height
  USE TRANSPORT_MOD         ! Driver routines for advection
  USE VDIFF_MOD             ! For non-local PBL mixing (J. Lin)

  !-----------------------------------------------------------------
  ! GEOS-Chem emissions modules
  !-----------------------------------------------------------------
  USE EMISSIONS_MOD         ! For interfacing with HEMCO emissions
  USE MIXING_MOD            ! performs tracer mixing
  USE MODIS_LAI_MOD         ! For MODIS leaf area indices (replacement)

  !-----------------------------------------------------------------
  ! GEOS-Chem meteorology field modules
  !-----------------------------------------------------------------
  USE Calc_Met_Mod          ! Met field calculations
  USE FLEXGRID_READ_MOD     ! For reading FLEXGRID data
#ifdef EXCHANGE
  USE EXCHANGE_MOD          ! For two-way coupled simulations
#endif

#ifdef RRTMG
  !-----------------------------------------------------------------
  ! Radiation modules (RRTMG)
  !-----------------------------------------------------------------
  USE RRTMG_RAD_TRANSFER_MOD, ONLY : Do_RRTMG_Rad_Transfer
  USE RRTMG_RAD_TRANSFER_MOD, ONLY : Init_RRTMG_Rad_Transfer
  USE RRTMG_RAD_TRANSFER_MOD, ONLY : Set_SpecMask
  USE RRTMG_LW_Init,          ONLY : RRTMG_LW_Ini
  USE RRTMG_SW_Init,          ONLY : RRTMG_SW_Ini
#endif

#ifdef APM
  !-----------------------------------------------------------------
  ! For APM aerosol microphysics simulation,
  ! see apm_driv_mod.f for more information
  !-----------------------------------------------------------------
  USE APM_INIT_MOD,      ONLY : APM_NTRACERS
  USE APM_INIT_MOD,      ONLY : APM_INIT
  USE APM_INIT_MOD,      ONLY : CLEANUP_APMARRAYS
  USE APM_DRIV_MOD,      ONLY : INIT_APM3D
  USE APM_DRIV_MOD,      ONLY : CLEANUP_APM3D
#endif

  IMPLICIT NONE
!
! !REMARKS:
!                                                                             .
!     GGGGGG  EEEEEEE  OOOOO  SSSSSSS       CCCCCC H     H EEEEEEE M     M
!    G        E       O     O S            C       H     H E       M M M M
!    G   GGG  EEEEEE  O     O SSSSSSS      C       HHHHHHH EEEEEE  M  M  M
!    G     G  E       O     O       S      C       H     H E       M     M
!     GGGGGG  EEEEEEE  OOOOO  SSSSSSS       CCCCCC H     H EEEEEEE M     M
!                                                                             .
!                                                                             .
!                 (formerly known as the Harvard-GEOS model)
!           for 4 x 5, 2 x 2.5 global grids and hi-res nested grids
!                                                                             .
!       Contact: GEOS-Chem Support Team (geos-chem-support@as.harvard.edu)
!
!                                                                             .
!  See the GEOS-Chem Web Site:
!                                                                             .
!     http://acmg.seas.harvard.edu/geos/
!                                                                             .
!  and the GEOS-Chem User's Guide:
!                                                                             .
!     http://acmg.seas.harvard.edu/geos/doc/man/
!                                                                             .
!  and the GEOS-Chem wiki:
!                                                                             .
!     http://wiki.seas.harvard.edu/geos-chem/
!                                                                             .
!  for the most up-to-date GEOS-Chem documentation on the following topics:
!                                                                             .
!     - installation, compilation, and execution
!     - coding practice and style
!     - input files and met field data files
!     - horizontal and vertical resolution
!     - modification history
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
  !-----------------------------
  ! Scalars
  !-----------------------------
  ! Logicals
  LOGICAL                  :: ITS_A_FULLCHEM_SIM
  LOGICAL                  :: ITS_A_MERCURY_SIM
  LOGICAL                  :: ITS_A_TAGCO_SIM
  LOGICAL                  :: ITS_AN_AEROSOL_SIM
  LOGICAL                  :: DO_DIAG_WRITE
  LOGICAL                  :: LCHEM
  LOGICAL                  :: LCONV
  LOGICAL                  :: LDRYD
  LOGICAL                  :: LDYNOCEAN
  LOGICAL                  :: LEMIS
  LOGICAL                  :: LGTMM
  LOGICAL                  :: LLINOZ
  LOGICAL                  :: LNLPBL
  LOGICAL                  :: LSTDRUN
  LOGICAL                  :: LSCHEM
  LOGICAL                  :: LSETH2O
  LOGICAL                  :: LTRAN
  LOGICAL                  :: LTURB
  LOGICAL                  :: LUCX
  LOGICAL                  :: LWETD
  LOGICAL                  :: prtDebug
  LOGICAL                  :: TimeForEmis
  LOGICAL                  :: notDryRun

  ! Integers
  INTEGER                  :: I,             IOS,         J
  INTEGER                  :: K,             L,           N
  INTEGER                  :: JDAY,          N_DYN
  INTEGER                  :: NNN,           N_DYN_STEPS, NSECb
  INTEGER                  :: N_STEP,        YEAR,        MONTH
  INTEGER                  :: DAY,           DAY_OF_YEAR
  INTEGER                  :: NYMD,          NYMDb,       NHMS
  INTEGER                  :: ELAPSED_SEC,   NHMSb,       RC
  INTEGER                  :: ELAPSED_TODAY, HOUR,        MINUTE
  INTEGER                  :: id_H2O,        id_CH4,      SECOND

  ! Reals
  REAL(f8)                 :: TAU,           TAUb

  ! Strings
  CHARACTER(LEN=255)       :: ThisLoc,       ZTYPE
  CHARACTER(LEN=255)       :: historyConfigFile
  CHARACTER(LEN=512)       :: ErrMsg
  CHARACTER(LEN=512)       :: Instr
  CHARACTER(LEN=255)       :: Argv

#ifdef RRTMG
  !-----------------------------
  ! Scalars specific to RRTMG
  !-----------------------------
  INTEGER                  :: iCld
  INTEGER                  :: iSeed
  INTEGER                  :: iSpecMenu
  INTEGER                  :: iNcDiag
  INTEGER                  :: RADSPEC
  LOGICAL, SAVE            :: FIRST_RT = .TRUE.
#endif

  !-----------------------------
  ! Derived type objects
  !-----------------------------
  TYPE(OptInput)           :: Input_Opt       ! Input Options object
  TYPE(ChmState)           :: State_Chm       ! Chemistry State object
  TYPE(DgnState)           :: State_Diag      ! Diagnostics State object
  TYPE(GrdState)           :: State_Grid      ! Grid State object
  TYPE(MetState)           :: State_Met       ! Meteorology State object
  TYPE(DgnList )           :: Diag_List       ! Diagnostics list object
  TYPE(TaggedDgnList )     :: TaggedDiag_List ! Tagged diagnostics list object
!
! !DEFINED PARAMETERS:
!
  ! When connecting G-C to an external GCM, we need to only write
  ! to stdout if we are on the root CPU.  Otherwise this will slow
  ! down the code.  This is why we introduced the am_I_Root logical
  ! variable.
  !
  ! However, if we are using the "traditional" G-C, then we don't
  ! need to restrict I/O to the root CPU.  (This is because each
  ! GEOS-Chem simulation starts out on a single CPU, with other
  ! CPUs joining only within parallel DO loops).  Therefore, we
  ! can just set am_I_Root = .true. here and then have it propagate
  ! down to all of the lower-level routines.  The main.F routine
  ! is not called when connecting G-C to an external GCM.
  ! (mlong, bmy, 7/30/12)
  LOGICAL, PARAMETER       :: am_I_Root = .TRUE.

  !=================================================================
  ! GEOS-CHEM starts here!
  !=================================================================

#ifdef TOMAS
  !(sfarina, 6/19/2013) It may seem strange, but this welcome message
  !                     fixes an issue where geoschem crashes with a
  !                     sigsegv immediately after starting.
  !                     This happens on ace-net's glooscap cluster with
  !                     ifort (IFORT) 11.1 20101201
  !                     this issue does not appear when running inside
  !                     a debugger, and is probably related to
  !                     some initialization garbage in memory
  !                     when using -O2 optimization
  !(bmy, 1/27/2014)   - Need to "CALL FLUSH(6).  FLUSH needs
  !                     an argument.  Unit 6 is Unix stdout.
  PRINT*, '%%%%% USING TOMAS MICROPHYSICS PACKAGE %%%%%'
  CALL FLUSH(6)
#endif

  ! Assume a successful return until otherwise
  RC      = GC_SUCCESS

  ! For error trapping
  ErrMsg  = ''
  ThisLoc = ' -> at GEOS-Chem (in GeosCore/main.F90)'

  ! Display model information
  CALL Display_Model_Info()

  !=================================================================
  !            ***** I N I T I A L I Z A T I O N *****
  !=================================================================

  !-----------------------------------------------------------------
  ! Read the user-defined options for the simulation, etc.
  !-----------------------------------------------------------------

  ! Initialize fields of the Input Options object
  CALL Set_Input_Opt( am_I_Root, Input_Opt, RC )
  IF ( RC /= GC_SUCCESS ) THEN
     ErrMsg = 'Error encountered within call to "Set_Input_Opt"!'
     CALL Error_Stop( ErrMsg, ThisLoc )
  ENDIF

  ! Initialize fields of the Grid State object
  CALL Init_State_Grid( Input_Opt, State_Grid, RC )
  IF ( RC /= GC_SUCCESS ) THEN
     ErrMsg = 'Error encountered within call to "Set_Grid_State"!'
     CALL Error_Stop( ErrMsg, ThisLoc )
  ENDIF

  ! Read GEOS-Chem input file at very beginning of simulation
  CALL Read_Input_File( Input_Opt, State_Grid, RC )
  IF ( RC /= GC_SUCCESS ) THEN
     ErrMsg = 'Error encountered in "Read_Input_File"!'
     CALL Error_Stop( ErrMsg, ThisLoc )
  ENDIF

  !-----------------------------------------------------------------
  ! Initialize GEOS-Chem timers
  !-----------------------------------------------------------------
  IF ( Input_Opt%useTimers ) THEN

     ! Call timer initilization
     CALL Timer_Setup( 1 )

     ! Add timers for various operations
     CALL Timer_Add( "GEOS-Chem",                    RC )
     CALL Timer_Add( "Initialization",               RC )
     CALL Timer_Add( "Timesteps",                    RC )
     CALL Timer_Add( "HEMCO",                        RC )
     CALL Timer_Add( "All chemistry",                RC )
     CALL Timer_Add( "=> Gas-phase chem",            RC )
     CALL Timer_Add( "=> FAST-JX photolysis",        RC )
     CALL Timer_Add( "=> All aerosol chem",          RC )
     CALL Timer_Add( "=> Strat chem",                RC )
     CALL Timer_Add( "=> Unit conversions",          RC )
     CALL Timer_Add( "Transport",                    RC )
     CALL Timer_Add( "Convection",                   RC )
     CALL Timer_Add( "Boundary layer mixing",        RC )
     CALL Timer_Add( "Dry deposition",               RC )
     CALL Timer_Add( "Wet deposition",               RC )
#ifdef RRTMG
     CALL Timer_Add( "RRTMG",                        RC )
#endif
     CALL Timer_Add( "All diagnostics",              RC )
     CALL Timer_Add( "=> HEMCO diagnostics",         RC )
#ifdef BPCH_DIAG
     CALL Timer_Add( "=> Binary punch diagnostics",  RC )
#endif
     CALL Timer_Add( "=> ObsPack diagnostics",       RC )
     CALL Timer_Add( "=> History (netCDF diags)",    RC )
     CALL Timer_Add( "Input",                        RC )
     CALL Timer_Add( "Output",                       RC )
     CALL Timer_Add( "Finalization",                 RC )

     ! Start running the main and initialization timer
     CALL Timer_Start( "GEOS-Chem",                  RC )
     CALL Timer_Start( "Initialization",             RC )
  ENDIF

  !-----------------------------------------------------------------
  ! Prepare the GEOS-Chem "dry run" option
  ! If in a "dry-run" mode, GEOS-Chem will simply check whether files
  ! are present (and possibly in the correct format) and go through
  ! time-steps to check met fields and other IO issues.
  ! No actual "compute" is performed.
  !
  ! The "dry-run" option is initialized using the command line extra
  ! argument ./geos --dry-run
  !
  ! A log file can be specified with --log FILENAME.
  ! If no log file is specified, the default logfile will be
  ! "GEOSChem.DryRun.log".
  !
  ! This option is currently only supported in GEOS-Chem Classic.
  !
  ! Additionally, this flag must be set after reading input file, or
  ! its value will be overwritten by READ_INPUT_FILE.
  ! (hplin, 11/1/19)
  !-----------------------------------------------------------------
  CALL Init_Dry_Run( Input_Opt, RC )
  IF ( RC /= GC_SUCCESS ) THEN
     ErrMsg = 'Error encountered in "Init_Dry_Run"!'
     CALL Error_Stop( ErrMsg, ThisLoc )
  ENDIF

  ! Define a local convenience variable for negating Input_Opt%DryRun
  notDryRun = ( .not. Input_Opt%DryRun )

  !-----------------------------------------------------------------
  ! Continue initializzation
  !-----------------------------------------------------------------

  ! Make sure all directories are valid
  CALL Validate_Directories( Input_Opt, RC )
  IF ( RC /= GC_SUCCESS ) THEN
     ErrMsg = 'Error encountered in  "Validate_Directories"!'
     CALL Error_Stop( ErrMsg, ThisLoc )
  ENDIF

  ! Initialize GEOS-Chem horizontal grid structure
  CALL GC_Init_Grid( Input_Opt, State_Grid, RC )
  IF ( RC /= GC_SUCCESS ) THEN
     ErrMsg = 'Error in "GC_Init_Grid"!'
     CALL Error_Stop( ErrMsg, ThisLoc )
  ENDIF

  ! Call the routine GC_Allocate_All (located in module file
  ! GeosCore/gc_environment_mod.F90) to allocate all lat/lon
  ! allocatable arrays used by GEOS-Chem.
  CALL GC_Allocate_All( Input_Opt, State_Grid, RC )
  IF ( RC /= GC_SUCCESS ) THEN
     ErrMsg = 'Error encountered in "GC_Allocate_All"!'
     CALL Error_Stop( ErrMsg, ThisLoc )
  ENDIF

  ! Store shadow copies of am_I_Root, Input_Opt in error_mod.F
  CALL Init_Error(Input_Opt, RC )
  IF ( RC /= GC_SUCCESS ) THEN
     ErrMsg = 'Error encountered in "Init_Error"!'
     CALL Error_Stop( ErrMsg, ThisLoc )
  ENDIF

  ! Copy values from Input_Opt.  These replace the variables
  ! from logical_mod.F and tracer_mod.F. (bmy, 3/29/13)
  ITS_A_FULLCHEM_SIM  =  Input_Opt%ITS_A_FULLCHEM_SIM
  ITS_A_MERCURY_SIM   =  Input_Opt%ITS_A_MERCURY_SIM
  ITS_A_TAGCO_SIM     =  Input_Opt%ITS_A_TAGCO_SIM
  ITS_AN_AEROSOL_SIM  =  Input_Opt%ITS_AN_AEROSOL_SIM
  DO_DIAG_WRITE       =  Input_Opt%DO_DIAG_WRITE
  LCHEM               =  Input_Opt%LCHEM
  LCONV               =  Input_Opt%LCONV
  LDRYD               =  Input_Opt%LDRYD
  LDYNOCEAN           =  Input_Opt%LDYNOCEAN
  LEMIS               =  Input_Opt%LEMIS
  LGTMM               =  Input_Opt%LGTMM
  LLINOZ              =  Input_Opt%LLINOZ
  LNLPBL              =  Input_Opt%LNLPBL
  LSCHEM              =  Input_Opt%LSCHEM
  LSETH2O             =  Input_Opt%LSETH2O
  LSTDRUN             =  Input_Opt%LSTDRUN
  LTRAN               =  Input_Opt%LTRAN
  LTURB               =  Input_Opt%LTURB
  LUCX                =  Input_Opt%LUCX
  LWETD               =  Input_Opt%LWETD

  ! Set a flag to denote if we should print ND70 debug output
  prtDebug            = ( Input_Opt%LPRT .and. am_I_Root )

  ! Turn off debug output for the dry-runs simulation
  prtDebug            = ( prtDebug .and. notDryRun )

  ! Debug output
  IF ( prtDebug ) CALL Debug_Msg( '### MAIN: a READ_INPUT_FILE' )

  !-----------------------------------------------------------------
  ! %%%% REPLICATING GCHP FUNCTIONALITY IN EXISTING GEOS-CHEM %%%%
  !
  ! Initialize the diagnostic list object which contains the
  ! unique entires in the history config file. Note that this is
  ! done in GCHP Set_Services and therefore must be done prior to
  ! initialization of the state objects. Also note that the diag_list
  ! obj may be stored in the HistoryConfig object in GCHP and we may
  ! want to replicate that behavior in GCC in the future.
  ! (ewl, 9/26/17)
  !-----------------------------------------------------------------
  IF ( notDryRun ) THEN
     IF ( Input_Opt%useTimers ) THEN
        CALL Timer_Start( "All diagnostics",           RC )
        CALL Timer_Start( "=> History (netCDF diags)", RC )
     ENDIF
     ! Don't initialize diagnostics when in dry-run mode
     historyConfigFile = 'HISTORY.rc' ! InputOpt not yet initialized
     CALL Init_DiagList( am_I_Root, historyConfigFile, Diag_List, RC )
     IF ( RC /= GC_SUCCESS ) THEN
        ErrMsg = 'Error encountered in "Init_DiagList"!'
        CALL Error_Stop( ErrMsg, ThisLoc )
     ENDIF
     CALL Init_TaggedDiagList( am_I_Root, Diag_List, TaggedDiag_List, RC )
     IF ( RC /= GC_SUCCESS ) THEN
        ErrMsg = 'Error encountered in "Init_TaggedDiagList"!'
        CALL Error_Stop( ErrMsg, ThisLoc )
     ENDIF

     !###  Print diagnostic lists if needed for debugging
     IF ( prtDebug ) CALL Print_DiagList( am_I_Root, Diag_List, RC )
     IF ( prtDebug ) CALL Print_TaggedDiagList( am_I_Root, TaggedDiag_List, RC )

     IF ( Input_Opt%useTimers ) THEN
        CALL Timer_End( "All diagnostics",           RC )
        CALL Timer_End( "=> History (netCDF diags)", RC )
     ENDIF
  ENDIF

  !-----------------------------------------------------------------
  ! %%%% REPLICATING GCHP FUNCTIONALITY IN EXISTING GEOS-CHEM %%%%
  !
  ! To replicate the functionality of the ESMF interface, we must
  ! initialize the Meteorology State (i.e. State_Met) and the
  ! Chemistry State (i.e. State_Chm) objects.  These objects hold
  ! several individual data fields that need to be passed as
  ! inputs to the chemistry routines.
  !
  ! The Meteorology State has replaced all of the individual
  ! met field arrays contained in module dao_mod.F. Likewise,
  ! the Chemistry State has replaced the STT tracer array
  ! and CSPEC chemical species array.
  !
  ! The Chemistry and Meteorology State objects facilitate using
  ! GEOS-Chem directly from the ESMF interface.  This is the main
  ! reason we are migrating towards used of these objects instead
  ! of the existing ALLOCATABLE module arrays. (bmy, 10/25/12)
  !-----------------------------------------------------------------

  ! Initialize State_Met, State_Chm, and State_Diag objects
  CALL GC_Init_StateObj( Diag_List,  Input_Opt,  State_Chm, &
                         State_Diag, State_Grid, State_Met, RC )
  IF ( RC /= GC_SUCCESS ) THEN
     ErrMsg = 'Error encountered in "GC_Init_All!"!'
     CALL Error_Stop( ErrMsg, ThisLoc )
  ENDIF

  !-----------------------------------------------------------------
  ! Skip these operations when in dry-run mode
  !-----------------------------------------------------------------
  IF ( notDryRun ) THEN
     ! Copy to State_Met%AREA_M2 to avoid breaking GCHP benchmarks,
     ! which require the AREA_M2 field saved out to the StateMet
     ! diagnostic collection for computing emission totals.
     State_Met%Area_M2 = State_Grid%Area_M2
  ENDIF

  !-----------------------------------------------------------------
  ! For dry-run simulations, call GC_Init_Extra, which will print
  ! out the Olson_drydep_inputs file name and exit.
  !
  ! For regular simulations, initialize various module arrays etc.
  ! This removes the init calls from the run-stage, which cannot
  ! happen when connecting GEOS-Chem to external ESMs.
  !-----------------------------------------------------------------     
  CALL GC_Init_Extra( Diag_List,  Input_Opt,  State_Chm, &
                      State_Diag, State_Grid, RC )
  IF ( RC /= GC_SUCCESS ) THEN
     ErrMsg = 'Error encountered in "GC_Init_Extra"!'
     CALL Error_Stop( ErrMsg, ThisLoc )
  ENDIF

  !-----------------------------------------------------------------
  ! Skip these operations when in dry-run mode
  !-----------------------------------------------------------------
  IF ( notDryRun ) THEN

     ! Initialize the regridding module by storing shadow copies
     CALL GC_Init_Regridding( Input_Opt, State_Grid, RC )
     IF ( RC /= GC_SUCCESS ) THEN
        ErrMsg = 'Error encountered in "Initialize_Regridding"!'
        CALL Error_Stop( ErrMsg, ThisLoc )
     ENDIF

     ! Define advected species ID flags for use below
     id_H2O   = Ind_('H2O', 'A')
     id_CH4   = Ind_('CH4', 'A')

     !-----------------------------------------------------------------
     ! OBSPACK Diagnostics: Get information from the species
     ! database for all requested ObsPack output species
     !-----------------------------------------------------------------
     IF ( Input_Opt%Do_ObsPack ) THEN
        CALL ObsPack_SpeciesMap_Init( Input_Opt, State_Chm, State_Diag, RC )
        IF ( RC /= GC_SUCCESS ) THEN
           ErrMsg = 'Error encountered in "ObsPack_SpeciesMap_Init"!'
           CALL Error_Stop( ErrMsg, ThisLoc )
        ENDIF
     ENDIF
  ENDIF

#ifdef RRTMG
  !-----------------------------------------------------------------
  ! Initializations for the RRTMG radiative transfer model
  !-----------------------------------------------------------------
  IF ( notDryRun ) THEN

     ! Initialize module variables
     CALL Init_RRTMG_Rad_Transfer( Input_Opt, State_Diag, State_Grid, RC )

     ! Trap potential errors
     IF ( RC /= GC_SUCCESS ) THEN
        ErrMsg = 'Error encountered in "Init_RRTMG_Rad_Transfer"!'
        CALL Error_Stop( ErrMsg, ThisLoc )
     ENDIF

     ! Initialize RRTMG code in the GeosRad folder
     CALL Rrtmg_Lw_Ini()
     CALL Rrtmg_Sw_Ini()

     ! Settings
     iCld  = 0
     iSeed = 10
  ENDIF
#endif

#ifdef APM
  !-----------------------------------------------------------------
  ! Initializations for the APM aerosol microphysics package
  !-----------------------------------------------------------------
  IF ( notDryRun ) THEN
     ! Initialize APM related variables, arrays
     CALL APM_NTRACERS( State_Chm )

     CALL APM_INIT(Input_Opt%CHEM_INPUTS_DIR)

     CALL Init_APM3D( Input_Opt, State_Grid, RC )
     IF ( RC /= GC_SUCCESS ) THEN
        ErrMsg = 'Error encountered in routine "Init_APM3D"!'
        CALL Error_Stop( ErrMsg, ThisLoc )
     ENDIF
  ENDIF
#endif

#ifdef BPCH_DIAG
  !-----------------------------------------------------------------
  ! Initialize bpch diagnostic arrays and counters
  !-----------------------------------------------------------------
  IF ( notDryRun ) THEN
     CALL Initialize( Input_Opt, State_Grid, 2, RC )
     CALL Initialize( Input_Opt, State_Grid, 3, RC )
     IF ( prtDebug ) CALL Debug_Msg( '### MAIN: a INITIALIZE' )
  ENDIF
#endif

  !-----------------------------------------------------------------
  ! Initialize the hybrid pressure module.  Define Ap and Bp.
  !-----------------------------------------------------------------
  CALL Init_Pressure( Input_Opt, State_Grid, RC )
  IF ( RC /= GC_SUCCESS ) THEN
     ErrMsg = 'Error encountered in "Init_Pressure"!'
     CALL Error_Stop( ErrMsg, ThisLoc )
  ENDIF
  IF ( prtDebug ) CALL Debug_Msg( '### MAIN: a INIT_PRESSURE' )

  !-----------------------------------------------------------------
  ! Register the horizontal and vertical grid information so that
  ! the History component can use it for netCDF metadata
  !-----------------------------------------------------------------
  IF ( notDryRun ) THEN
     CALL Init_Grid_Registry( Input_Opt, State_Grid, RC )
     IF ( RC /= GC_SUCCESS ) THEN
        ErrMsg = 'Error encountered in "Init_Grid_Registry"!'
        CALL Error_Stop( ErrMsg, ThisLoc )
     ENDIF
  ENDIF

  !-----------------------------------------------------------------
  ! Added to read input file for linoz strat chem
  !-----------------------------------------------------------------
  IF ( LLINOZ ) THEN
     CALL Linoz_Read( Input_Opt, RC )

     ! Trap potential errors
     IF ( RC /= GC_SUCCESS ) THEN
        ErrMsg = 'Error encountered in "Linoz_Read"!'
        CALL Error_Stop( ErrMsg, ThisLoc )
     ENDIF
  ENDIF

  ! Define time variables for use below
  NHMS  = GET_NHMS()
  NHMSb = GET_NHMSb()
  NYMD  = GET_NYMD()
  NYMDb = GET_NYMDb()
  TAU   = GET_TAU()
  TAUb  = GET_TAUb()

  !-----------------------------------------------------------------
  !    ***** H I S T O R Y   I N I T I A L I Z A T I O N *****
  !-----------------------------------------------------------------
  IF ( Input_Opt%useTimers .and. notDryRun ) THEN
     CALL Timer_Start( "All diagnostics",           RC )
     CALL Timer_Start( "=> History (netCDF diags)", RC )
  ENDIF

  ! For now, just hardwire the input file for the History component
  Input_Opt%HistoryInputFile = './HISTORY.rc'

  ! Initialize the GEOS-Chem history component
  ! (for dry-run, enter routine to print out HISTORY.rc status)
  CALL History_Init( Input_Opt,  State_Met, State_Chm, State_Diag, RC )

  ! Trap error
  IF ( RC /= GC_SUCCESS ) THEN
     ErrMsg = 'Error encountered in "History_Init"!'
     CALL Error_Stop( ErrMsg, ThisLoc )
  ENDIF

  IF ( Input_Opt%useTimers .and. notDryRun ) THEN
     CALL Timer_End( "All diagnostics",           RC )
     CALL Timer_End( "=> History (netCDF diags)", RC )
  ENDIF

  !-----------------------------------------------------------------
  !        ***** I N I T I A L I Z A T I O N  continued *****
  !-----------------------------------------------------------------

  ! To enable FlexGrid, need to initialize HEMCO and run phase 1
  ! before reading initial metfields.
  ! (Jiawei Zhuang 2017/6)

  ! Initialize HEMCO. This reads the HEMCO configuration file
  ! and creates entries for all data files needed for emission
  ! calculation. Also sets some logical switches in Input_Opt
  ! (e.g. LSOILNOX).
  ! Note: always call HEMCO, even if LEMIS is set to FALSE. This
  ! is to make sure that HEMCO can still be used to read
  ! non-emission data such as stratospheric Bry fields. If LEMIS
  ! is set to FALSE, the emissions driver routines will make sure
  ! that HEMCO does not calculate any emissions (ckeller, 1/12/15).
  !
  ! This call will also initialize the three built-in HEMCO
  ! diagnostics (default, manual, restart).
  IF ( Input_Opt%useTimers ) THEN
     CALL Timer_Start( "HEMCO", RC )
  ENDIF

  CALL Emissions_Init( Input_Opt, State_Chm, State_Grid, State_Met, RC )

  ! Trap potential errors
  IF ( RC /= GC_SUCCESS ) THEN
     ErrMsg = 'Error encountered in "Emissions_Init"!'
     CALL Error_Stop( ErrMsg, ThisLoc )
  ENDIF

  ! Run HEMCO phase 0 as simplfied phase 1 to get initial met fields
  ! and restart file fields
  CALL Emissions_Run( Input_Opt, State_Chm,   State_Diag, State_Grid,  &
                      State_Met, TimeForEmis, 0,          RC )

  ! Trap potential errors
  IF ( RC /= GC_SUCCESS ) THEN
     ErrMsg = 'Error encountered in "Emissions_Run", Phase 0'
     Instr  = 'This error can indicate a missing file. Please check '// &
              'the HEMCO log file for additional error messages! '
     CALL Error_Stop( ErrMsg, ThisLoc, Instr )
  ENDIF

  IF ( Input_Opt%useTimers ) THEN
     CALL Timer_End ( "HEMCO", RC )
  ENDIF

  ! Skip certain initializations
  IF ( notDryRun ) THEN

     ! Populate the State_Met%LandTypeFrac field with data from HEMCO
     CALL Init_LandTypeFrac( Input_Opt, State_Met, RC )
     IF ( RC /= GC_SUCCESS ) THEN
        ErrMsg = 'Error encountered in "Init_LandTypeFrac"!'
        CALL Error_Stop( ErrMsg, ThisLoc )
     ENDIF

     ! Compute the Olson landmap fields of State_Met
     ! (e.g. State_Met%IREG, State_Met%ILAND, etc.)
     CALL Compute_Olson_Landmap( Input_Opt, State_Grid, State_Met, RC )
     IF ( RC /= GC_SUCCESS ) THEN
        ErrMsg = 'Error encountered in "Compute_Olson_Landmap"!'
        CALL Error_Stop( ErrMsg, ThisLoc )
     ENDIF

     ! Initialize PBL quantities but do not do mixing
     ! Add option for non-local PBL (Lin, 03/31/09)
     CALL Init_Mixing( Input_Opt,  State_Chm, State_Diag, &
                       State_Grid, State_Met, RC )

     ! Trap potential errors
     IF ( RC /= GC_SUCCESS ) THEN
        ErrMsg = 'Error encountered in Init_Mixing!'
        CALL Error_Stop( ErrMsg, ThisLoc )
     ENDIF
  ENDIF

  ! Initialize chemistry
  ! Moved here because some of the variables are used for non-local
  ! PBL mixing BEFORE the first call of the chemistry routines
  ! (ckeller, 05/19/14).
  IF ( ITS_A_FULLCHEM_SIM .OR. ITS_AN_AEROSOL_SIM ) THEN
     CALL Init_Chemistry( Input_Opt,  State_Chm, State_Diag, State_Grid, RC )

     ! Trap potential errors
     IF ( RC /= GC_SUCCESS ) THEN
        ErrMsg = 'Error encountered in "Init_Chemistry"!'
        CALL Error_Stop( ErrMsg, ThisLoc )
     ENDIF
  ENDIF

  !=================================================================
  !       *****  I N I T I A L   C O N D I T I O N S *****
  !=================================================================

  ! Initialize the UCX module
  IF ( LUCX ) THEN
     CALL INIT_UCX( Input_Opt, State_Chm, State_Diag, State_Grid )
     IF ( prtDebug ) CALL DEBUG_MSG( '### MAIN: a INIT_UCX' )
  ENDIF

  ! Capture initial state of atmosphere for STE flux calc (ltm, 06/10/12)
  ! NOTE: Species concentrations enter the subroutine in [kg/kg dry air]
  ! and are converted locally to [kg] for chemistry (ewl, 9/18/15)
  IF ( LSCHEM .and. notDryRun ) THEN
     CALL Init_Strat_Chem( Input_Opt,  State_Chm, State_Met, State_Grid, RC )

     ! Trap potential errors
     IF ( RC /= GC_SUCCESS ) THEN
        ErrMsg = 'Error encountered in "Init_Strat_Chem"!'
        CALL Error_Stop( ErrMsg, ThisLoc )
     ENDIF
  ENDIF

  !-----------------------------------------------------------------------------
  ! TWO-WAY NESTING OPTION
  ! This is only invoked when compiling GEOS-Chem with COUPLE=y
  !
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! %%%%% NOTE: THIS OPTION WAS PROBABLY BROKEN WITH THE ADDITION OF HEMCO %%%%%
  ! %%%%% (v10-01), FLEXCHEM (v11-01), AND FLEXGRID (12.5.0). BUYER BEWARE.%%%%%
  ! %%%%%  --  Bob Yantosca (22 Jan 2018)                                  %%%%%
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !
#if defined( EXCHANGE )
  ! Initialize the two-way nesting.  This will only get done if you
  ! compile GEOS-Chem with EXCHANGE=yes. (yanyy, 03/28/14)
  CALL INIT_EXCHANGE( State_Grid )

  IF ( State_Grid%NestedGrid ) THEN

     ! Initialize exchange of nested boundary conditions
     CALL EXCHANGE_NESTED_INIT()

  ELSE

     ! Initialize exchange of global boundary conditions
     CALL EXCHANGE_GLOBAL_INIT()

  ENDIF
#endif
!-----------------------------------------------------------------------------

  IF ( Input_Opt%useTimers ) THEN
     CALL Timer_End( "Initialization", RC )
  ENDIF

  !=================================================================
  !        ***** O U T E R   T I M E S T E P   L O O P  *****
  !=================================================================      

  IF ( Input_Opt%useTimers ) THEN
     CALL Timer_Start( "Timesteps", RC )
  ENDIF

  ! Echo message before first timestep
  IF ( notDryRun ) THEN
     WRITE( 6, '(a)' )
     WRITE( 6, '(a)' ) REPEAT( '*', 44 )
     WRITE( 6, '(a)' ) '* B e g i n   T i m e   S t e p p i n g !! *'
     WRITE( 6, '(a)' ) REPEAT( '*', 44 )
     WRITE( 6, '(a)' )
  ENDIF

  ! NSTEP is the number of dynamic timesteps w/in the outer loop
  ! Timesteps are now retrieved in seconds (ewl, 2/6/2018)
  N_DYN_STEPS = 10800 / GET_TS_DYN()     ! 3hr interval

  ! Start a new outer loop
  DO

    ! Compute time parameters at start of 6-h loop
    CALL Set_Current_Time()

    ! NSECb is # of seconds (measured from 00 GMT today)
    ! at the start of this 6-hr timestepping loop.
    ! NOTE: Assume we start at the head of each minute (i.e. SECONDS=0)
    HOUR   = GET_HOUR()
    HOUR   = ( HOUR / 6 ) * 6
    MINUTE = GET_MINUTE()
    SECOND = GET_SECOND()
    NSECb  = ( HOUR * 3600 ) + ( MINUTE * 60 ) + SECOND

    ! Get dynamic timestep in seconds
    N_DYN  = GET_TS_DYN()

    !=================================================================
    !     ***** D Y N A M I C   T I M E S T E P   L O O P *****
    !     *****    a k a   H E A R T B E A T   L O O P    *****
    !=================================================================
    DO N_STEP = 1, N_DYN_STEPS

       ! Compute & print time quantities at start of dyn step
       CALL Set_Current_Time()
       IF ( notDryRun ) CALL Print_Current_Time()

       ! Set time variables for dynamic loop
       DAY_OF_YEAR   = GET_DAY_OF_YEAR()
       DAY           = GET_DAY()
       ELAPSED_SEC   = GET_ELAPSED_SEC()
       MONTH         = GET_MONTH()
       NHMS          = GET_NHMS()
       NYMD          = GET_NYMD()
       HOUR          = GET_HOUR()
       MINUTE        = GET_MINUTE()
       SECOND        = GET_SECOND()
       TAU           = GET_TAU()
       YEAR          = GET_YEAR()
       ELAPSED_TODAY = ( HOUR * 3600 ) + ( MINUTE * 60 ) + SECOND

       IF ( prtDebug ) THEN
          CALL Debug_Msg( '### MAIN: a SET_CURRENT_TIME' )
       ENDIF

       !--------------------------------------------------------------
       ! %%%%% HISTORY (netCDF diagnostics) %%%%%
       !
       ! Certain Hg simulation diagnostics (e.g. deposition of Hg2
       ! and HgP onto snow and ice) need to be zeroed out at the
       ! start each timestep, before operations like drydep, wetdep,
       ! and convection are executed.  Call a routine to do this.
       ! (bmy, 10/25/18)
       !--------------------------------------------------------------
       IF ( ITS_A_MERCURY_SIM .and. notDryRun ) THEN
          IF ( Input_Opt%useTimers ) THEN
             CALL Timer_Start( "All diagnostics",           RC )
             CALL Timer_Start( "=> History (netCDF diags)", RC )
          ENDIF

          CALL Reset_Hg_Diags( Input_Opt, State_Diag, RC )
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = 'Error encountered in "Reset_Hg_Diags!"'
             CALL Error_Stop( ErrMsg, ThisLoc )
          ENDIF

          IF ( Input_Opt%useTimers ) THEN
             CALL Timer_End( "All diagnostics",           RC )
             CALL Timer_End( "=> History (netCDF diags)", RC )
          ENDIF
       ENDIF

       !==============================================================
       !       ***** R U N   H E M C O   P H A S E   1 *****
       !
       !    Phase 1 updates the HEMCO clock and the content of the
       !    HEMCO data list. This should be done before writing the
       !    diagnostics organized in the HEMCO diagnostics structure,
       !    and before using any of the HEMCO data list fields.
       !    (ckeller, 4/1/15)
       !==============================================================
       IF ( Input_Opt%useTimers ) THEN
          CALL Timer_Start( "HEMCO", RC )
          CALL Timer_Start( "Input", RC )
       ENDIF

       ! Is it time for emissions?
       TimeForEmis = ITS_TIME_FOR_EMIS()

       ! Run HEMCO Phase 1
       CALL Emissions_Run( Input_Opt, State_Chm,   State_Diag, State_Grid, &
                           State_Met, TimeForEmis, 1,          RC )

       ! Trap potential errors
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Error encountered in "Emissions_Run"!'
          CALL Error_Stop( ErrMsg, ThisLoc )
       ENDIF

       IF ( Input_Opt%useTimers ) THEN
          CALL Timer_End( "HEMCO", RC )
          CALL Timer_End( "Input", RC )
       ENDIF

       IF ( prtDebug ) THEN
          CALL Debug_Msg( '### MAIN: a HEMCO PHASE 1' )
       ENDIf

       !==============================================================
       !  ***** W R I T E   H E M C O   D I A G N O S T I C S *****
       !==============================================================
       IF ( notDryRun ) THEN
          IF ( Input_Opt%useTimers ) THEN
             CALL Timer_Start( "All diagnostics",       RC )
             CALL Timer_Start( "=> HEMCO diagnostics",  RC )
             CALL Timer_Start( "Output",                RC )
          ENDIF

          ! Do not do actual output for dry-run
          ! Write HEMCO diagnostics (ckeller, 4/1/15)
          CALL HCOI_GC_WriteDiagn( Input_Opt, .FALSE., RC )

          ! Trap potential errors
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = 'Error encountered in "HCOI_GC_WriteDiagn" ' // &
                      '(writing HEMCO diagnostics)!'
             CALL Error_Stop( ErrMsg, ThisLoc )
          ENDIF

          IF ( Input_Opt%useTimers ) THEN
             CALL Timer_End( "All diagnostics",      RC )
             CALL Timer_End( "=> HEMCO diagnostics", RC )
             CALL Timer_End( "Output",               RC )
          ENDIF
       ENDIF

       !==============================================================
       !      ***** O B S P A C K   D I A G N O S T I C S *****
       !==============================================================
       IF ( Input_Opt%Do_ObsPack .and. &
          ( ELAPSED_TODAY == 0 ) .and. notDryRun ) THEN

          IF ( Input_Opt%useTimers ) THEN
             CALL Timer_Start( "All diagnostics",         RC )
             CALL Timer_Start( "=> ObsPack diagnostics",  RC )
          ENDIF

          ! Initialize Obspack for the new day
          CALL ObsPack_Init( NYMD,  NHMS, Input_Opt, State_Diag, RC )

          ! Trap potential errors
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = 'Error encountered in "ObsPack_Init"!'
             CALL Error_Stop( ErrMsg, ThisLoc )
          ENDIF

          IF ( Input_Opt%useTimers ) THEN
             CALL Timer_End( "All diagnostics",         RC )
             CALL Timer_End( "=> ObsPack diagnostics",  RC )
          ENDIF
       ENDIF

#ifdef BPCH_DIAG
       !===========================================================
       ! *****  W R I T E   B P C H   D I A G N O S T I C S *****
       !===========================================================
       IF ( ITS_TIME_FOR_BPCH( Input_Opt ) .and. notDryRun ) THEN

          IF ( Input_Opt%useTimers ) THEN
             CALL Timer_Start( "All diagnostics",              RC )
             CALL Timer_Start( "=> Binary punch diagnostics",  RC )
             CALL Timer_Start( "Output",                       RC )
          ENDIF

          ! Set time at end of diagnostic timestep
          CALL SET_DIAGe( TAU )

          ! Write bpch file
          IF ( DO_DIAG_WRITE ) THEN
             ! Write data to the "trac_avg." bpch file
             CALL DIAG3( Input_Opt, State_Chm, State_Grid, State_Met, RC )

             ! Flush file units
             CALL CTM_FLUSH()
          ENDIF

          ! Set time at beginning of next diagnostic timestep
          CALL SET_DIAGb( TAU )

          !===========================================================
          !   ***** Z E R O   B P C H   D I A G N O S T I C S *****
          !===========================================================
          CALL INITIALIZE( Input_Opt, State_Grid, 2, RC )
          CALL INITIALIZE( Input_Opt, State_Grid, 3, RC )

          IF ( Input_Opt%useTimers ) THEN
             CALL Timer_End( "All diagnostics",              RC )
             CALL Timer_End( "=> Binary punch diagnostics",  RC )
             CALL Timer_End( "Output",                       RC )
          ENDIF
       ENDIF
#endif

#ifdef GTMM_Hg
#ifdef BPCH_DIAG
       !=============================================================
       !   ***** W R I T E   G T M M   R E S T A R T   F I L E *****
       !     ***** MUST be done after call to diag3 ****
       !
       ! %%%%% NOTE: THIS MAY BE BROKEN %%%%
       !=============================================================
       ! Make land restart file: for GTMM runs only, beginning of each
       ! month but not start of the run.
       IF ( LGTMM .AND. ITS_A_NEW_MONTH() .AND. NYMD /= NYMDb ) THEN
          IF (.NOT.( ITS_TIME_FOR_BPCH( Input_Opt ) )) THEN

             ! Get the species ID (NNN) from the wetdep ID (N)
             N   = 1
             NNN = State_Chm%Map_Wetdep(N)

             DO
                ! Exit once we encounter Hg2
                If ( State_Chm%SpcData(NNN)%Info%Is_Hg2 ) THEN
                   EXIT
                ENDIF

                ! Get the species ID (NNN) from the wetdep ID (N)
                N   = N + 1
                NNN = State_Chm%Map_Wetdep(N)
             ENDDO
             CALL UPDATE_DEP( N )
          ENDIF
          CALL MAKE_GTMM_RESTART( Input_Opt, State_Grid, NYMD, NHMS, TAU, RC )
       ENDIF
#endif
#endif

       !==============================================================
       !       ***** T E S T   F O R   E N D   O F   R U N *****
       !==============================================================
       IF ( ITS_TIME_FOR_EXIT() ) GOTO 9999

       !==============================================================
       !        ***** L E A F   A R E A   I N D I C E S *****
       !==============================================================
       IF ( ITS_A_NEW_DAY() .and. notDryRun ) THEN

          ! Initialize the State_Met%XLAI_NATIVE field from HEMCO
          CALL Get_XlaiNative_from_HEMCO( Input_Opt, State_Met, RC )

          ! Trap potential errors
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = 'Error encountered in "Get_XlaiNative_from_HEMCO"!'
             CALL Error_Stop( ErrMsg, ThisLoc )
          ENDIF

          ! Compute State_Met%XLAI (for drydep) and State_Met%MODISLAI,
          ! which is the average LAI per grid box (for soil NOx emissions)
          CALL Compute_Xlai( Input_Opt, State_Grid, State_Met, RC )

          ! Trap potential errors
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = 'Error encountered in "Compute_Xlai"!'
             CALL Error_Stop( ErrMsg, ThisLoc )
          ENDIF
       ENDIF

       !==============================================================
       !   ***** I N T E R P O L A T E   Q U A N T I T I E S *****
       !==============================================================

       ! Do not compute any data in dry-run mode
       IF ( notDryRun ) THEN

          ! Interpolate I-3 fields to current dynamic timestep,
          ! based on their values at NSEC and NSEC+N_DYN
          CALL Interp( NSECb,     ELAPSED_TODAY, N_DYN, &
                       Input_Opt, State_Grid,    State_Met )

          ! If we are not doing transport, then make sure that
          ! the floating pressure is set to PSC2_WET (bdf, bmy, 8/22/02)
          ! Now also includes PSC2_DRY (ewl, 5/4/16)
          IF ( .not. LTRAN ) THEN
             CALL Set_Floating_Pressures( State_Grid, State_Met, RC )

             ! Trap potential errors
             IF ( RC /= GC_SUCCESS ) THEN
                ErrMsg = 'Error encountered in "Set_Floating_Pressures"!'
                CALL Error_Stop( ErrMsg, ThisLoc )
             ENDIF
          ENDIF

          ! Compute updated airmass quantities at each grid box
          ! and update tracer concentration to conserve tracer mass
          ! (ewl, 10/28/15)
          CALL AirQnt( Input_Opt, State_Chm, State_Grid, State_Met, &
                       RC, Update_Mixing_Ratio=.TRUE. )

          ! Trap potential errors
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = 'Error encountered in "AirQnt"!'
             CALL Error_Stop( ErrMsg, ThisLoc )
          ENDIF

          ! SDE 05/28/13: Set H2O to State_Chm tracer if relevant and,
          ! if LUCX=T and LSETH2O=F and LACTIVEH2O=T, update specific humidity
          ! in the stratosphere
          !
          ! NOTE: Specific humidity may change in SET_H2O_TRAC and
          ! therefore this routine may call AIRQNT again to update
          ! air quantities and tracer concentrations (ewl, 10/28/15)
          IF ( ITS_A_FULLCHEM_SIM .and. id_H2O > 0 ) THEN
             CALL Set_H2O_Trac( ( ( .not. LUCX ) .or. LSETH2O ), &
                                Input_Opt, State_Chm, State_Grid, &
                                State_Met, RC )

             ! Trap potential errors
             IF ( RC /= GC_SUCCESS ) THEN
                ErrMsg = 'Error encountered in "Set_H2O_Trac" #1!'
                CALL Error_Stop( ErrMsg, ThisLoc )
             ENDIF

             ! Only force strat once if using UCX
             IF (LSETH2O) LSETH2O = .FALSE.
          ENDIF

          ! Compute the cosine of the solar zenith angle array
          ! State_Met%SUNCOS     = at the current time
          ! State_Met%SUNCOSmid  = at the midpt of the chem timestep
          CALL Get_Cosine_SZA( Input_Opt, State_Grid, State_Met, RC )

          ! Trap potential errors
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = 'Error encountered in "Get_Cosine_SZA"!'
             CALL Error_Stop( ErrMsg, ThisLoc )
          ENDIF
       ENDIF

       IF ( prtDebug ) CALL Debug_Msg( '### MAIN: a INTERP, etc' )

       !----------------------------------------------------------
       ! %%% GET SOME NON-EMISSIONS DATA FIELDS VIA HEMCO %%%
       !
       ! HEMCO can track non-emission data fields for chemistry
       ! simulations.  Put these subroutine calls after the
       ! call to EMISSIONS_RUN, so that the HEMCO data structure
       ! will be initialized. (bmy, 3/20/15)
       !
       ! HEMCO data list is now updated further above, so can
       ! take these calls out of the emissions sequence.
       ! (ckeller, 4/01/15)
       !----------------------------------------------------------
       IF ( LCHEM .and. ITS_A_NEW_MONTH() .and. notDryRun ) THEN

          ! The following only apply when photolysis is used,
          ! that is for fullchem or aerosol simulations.
          IF ( ITS_A_FULLCHEM_SIM  .or. ITS_AN_AEROSOL_SIM ) THEN

             ! Copy UV Albedo data (for photolysis) into the
             ! State_Met%UVALBEDO field. (bmy, 3/20/15)
             CALL Get_UvAlbedo( Input_Opt, State_Met, RC )

             ! Trap potential errors
             IF ( RC /= GC_SUCCESS ) THEN
                ErrMsg = 'Error encountered in "Get_UvAlbedo"!'
                CALL Error_Stop( ErrMsg, ThisLoc )
             ENDIF

             IF ( Input_Opt%USE_TOMS_O3 ) THEN
                ! Get TOMS overhead O3 columns for photolysis from
                ! the HEMCO data structure (bmy, 3/20/15)
                CALL Read_TOMS( Input_Opt, RC )

                ! Trap potential errors
                IF ( RC /= GC_SUCCESS ) THEN
                   ErrMsg = 'Error encountered in "Read_TOMS"!'
                   CALL Error_Stop( ErrMsg, ThisLoc )
                ENDIF
             ENDIF

          ENDIF

          ! Read data required for Hg2 gas-particle partitioning
          ! (H Amos, 25 Oct 2011)
          IF ( ITS_A_MERCURY_SIM .and. notDryRun ) THEN
             CALL Read_Hg2_Partitioning( Input_Opt, State_Grid, State_Met, &
                                         MONTH,     RC )

             ! Trap potential errors
             IF ( RC /= GC_SUCCESS ) THEN
                ErrMsg = 'Error encountered in "Read_Hg2_Partitioning"!'
                CALL Error_Stop( ErrMsg, ThisLoc )
             ENDIF

             IF ( prtDebug ) THEN
                CALL Debug_Msg( '### MAIN: a READ_HG2_PARTITIONING')
             ENDIF
          ENDIF
       ENDIF

       ! Prescribe methane surface concentrations throughout PBL
       IF ( ITS_A_FULLCHEM_SIM .and. id_CH4 > 0 .and. notDryRun ) THEN

          IF ( prtDebug ) THEN
             CALL DEBUG_MSG( '### MAIN: Setting PBL CH4 conc')
          ENDIF

          ! Set CH4 concentrations
          CALL SET_CH4( Input_Opt, State_Chm, State_Diag, State_Grid, &
                        State_Met, RC )

          ! Trap potential errors
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = 'Error encountered in call to "SET_CH4"!'
             CALL Error_Stop( ErrMsg, ThisLoc )
          ENDIF
       ENDIF

       !==============================================================
       !              ***** T R A N S P O R T *****
       !==============================================================
       IF ( ITS_TIME_FOR_DYN() .and. notDryRun ) THEN

          IF ( Input_Opt%useTimers ) THEN
             CALL Timer_Start( "Transport", RC )
          ENDIF

          !--------------------------------------------------------------------
          ! TWO-WAY NESTING OPTION
          ! This is only invoked when compiling GEOS-Chem with COUPLE=y
          !
          ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          ! %%%%% NOTE: THIS OPTION WAS PROBABLY BROKEN WITH THE ADDITION OF
          ! %%%%% HEMCO (v10-01), FLEXCHEM (v11-01), and FLEXGRID (12.5.0).
          ! %%%%% BUYER BEWARE. --  Bob Yantosca (22 Jan 2018)
          ! %%%%%
          ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#if defined( EXCHANGE )
          IF ( State_Grid%NestedGrid ) THEN

             ! Exchange the position of POST (nested-grid simulations)
             IF ( ITS_TIME_FOR_EXCHANGE() ) THEN
                CALL EXCHANGE_NESTED_POST( Input_Opt, State_Chm, State_Grid, &
                                           State_Met, RC )
             ENDIF

          ELSE

             ! Exchange the position of POST (global simulations)
             IF ( ITS_TIME_FOR_EXCHANGE() ) THEN
                CALL EXCHANGE_GLOBAL_POST( Input_Opt, State_Chm, State_Grid, &
                                           State_Met, RC )
             ENDIF

          ENDIF
#endif
          !--------------------------------------------------------------------

          ! Call the appropriate version of TPCORE
          IF ( LTRAN ) THEN
             CALL Do_Transport( Input_Opt,  State_Chm, State_Diag, &
                                State_Grid, State_Met, RC )

             ! Trap potential error
             IF ( RC /= GC_SUCCESS ) THEN
                ErrMsg = 'Error encountered in "Do_Transport"!'
                CALL Error_Stop( ErrMsg, ThisLoc )
             ENDIF

             IF ( prtDebug ) THEN
                CALL Debug_Msg( '### MAIN: a DO_TRANSPORT' )
             ENDIF
          ENDIF

          ! Initialize wet scavenging and wetdep fields after
          ! the airmass quantities are reset after transport
#ifdef TOMAS
          ! ... TOMAS microphysics: Always call SETUP_WETSCAV ...
          CALL Setup_WetScav( Input_Opt, State_Chm, State_Grid, &
                              State_Met, RC )

          ! Trap potential error
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = 'Error encountered in "Setup_WetScav" (TOMAS)!'
             CALL Error_Stop( ErrMsg, ThisLoc )
          ENDIF
#else
          ! ... Standard GEOS-Chem: Call INIT_WETSCAV if   ...
          ! ... convection or wet scavenging or chemistry are turned on ...
          IF ( LCONV .or. LWETD .or. LCHEM ) THEN
             CALL Setup_WetScav( Input_Opt, State_Chm, State_Grid, &
                                 State_Met, RC )

             ! Trap potential error
             IF ( RC /= GC_SUCCESS ) THEN
                ErrMsg = 'Error encountered in "Setup_WetScav"!'
                CALL Error_Stop( ErrMsg, ThisLoc )
             ENDIF

             IF ( prtDebug ) THEN
                CALL Debug_Msg( '### MAIN: a SETUP_WETSCAV' )
             ENDIF
          ENDIF
#endif

          IF ( Input_Opt%useTimers ) THEN
             CALL Timer_End( "Transport", RC )
          ENDIF

       ENDIF

       ! Update age of air (skip if running in dry-run mode)
       IF ( notDryRun ) THEN
          CALL Set_Met_AgeOfAir( State_Grid, State_Met )
       ENDIF

       !==============================================================
       !     ***** C O M P U T E   P B L   H E I G H T  etc. *****
       !==============================================================
       IF ( notDryRun ) THEN
          IF ( Input_Opt%useTimers ) THEN
             CALL Timer_Start( "Boundary layer mixing", RC )
          ENDIF

          ! Move this call from the PBL mixing routines because the PBL
          ! height is used by drydep and some of the emissions routines.
          ! (ckeller, 3/5/15)
          CALL Compute_PBL_Height( State_Grid, State_Met, RC )

          ! Trap potential errors
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = 'Error encountered in "Compute_PBL_Height"!'
             CALL Error_Stop( ErrMsg, ThisLoc )
          ENDIF

          IF ( Input_Opt%useTimers ) THEN
             CALL Timer_End( "Boundary layer mixing", RC )
          ENDIF

          IF ( prtDebug ) THEN
             CALL Debug_Msg( '### MAIN: a COMPUTE_PBL_HEIGHT' )
          ENDIF
       ENDIF

       !--------------------------------------------------------------
       ! Test for emission timestep
       ! Now always do emissions here, even for full-mixing
       ! (ckeller, 3/5/15)
       !
       ! Emissions are ALWAYS done, even in dry-run mode. This is
       ! raison d'etre for --dry-run (hplin, 11/1/19)
       !--------------------------------------------------------------
       IF ( ITS_TIME_FOR_EMIS() ) THEN

#ifdef BPCH_DIAG
          ! Increment emission counter
          CALL Set_Ct_Emis( INCREMENT=.TRUE. )
#endif

          !===========================================================
          !         ***** D R Y   D E P O S I T I O N *****
          !===========================================================
          IF ( LDRYD .and. notDryRun ) THEN

             IF ( Input_Opt%useTimers ) THEN
                CALL Timer_Start( "Dry deposition", RC )
             ENDIF

             ! Compute drydep velocities
             CALL Do_Drydep( Input_Opt,  State_Chm, State_Diag, &
                             State_Grid, State_Met, RC )

             ! Trap potential errors
             IF ( RC /= GC_SUCCESS ) THEN
                ErrMsg = 'Error encountered in "Do_Drydep!"!'
                CALL Error_Stop( ErrMsg, ThisLoc )
             ENDIF

             IF ( Input_Opt%useTimers ) THEN
                CALL Timer_End ( "Dry deposition", RC )
             ENDIF

             IF ( prtDebug ) THEN
                CALL Debug_Msg( '### MAIN: a DO_DRYDEP' )
             ENDIF
          ENDIF

          !===========================================================
          !             ***** E M I S S I O N S *****
          !
          ! NOTE: For a complete description of how emissions from
          ! HEMCO are added into GEOS-Chem (and how they are mixed
          ! into the boundary layer), please see the wiki page:
          !
          ! http://wiki-geos-chem.org/Distributing_emissions_in_the_PBL
          !===========================================================
          IF ( Input_Opt%useTimers ) THEN
             CALL Timer_Start( "HEMCO", RC )
          ENDIF

          ! Is it time for emissions?
          TimeForEmis = ITS_TIME_FOR_EMIS()

          ! EMISSIONS_RUN will call HEMCO run phase 2. HEMCO run phase
          ! only calculates emissions. All data has been read to disk
          ! in phase 1 at the beginning of the time step.
          ! (ckeller, 4/1/15)
          CALL Emissions_Run( Input_Opt, State_Chm,   State_Diag, State_Grid, &
                              State_Met, TimeForEmis, 2,          RC )

          ! Trap potential errors
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = 'Error encountered in "Emissions_Run"! after drydep!'
             CALL Error_Stop( ErrMsg, ThisLoc )
          ENDIF

          IF ( prtDebug ) THEN
             CALL Debug_Msg( '### MAIN: a HEMCO PHASE 2' )
          ENDIF

          IF ( Input_Opt%useTimers ) THEN
             CALL Timer_End( "HEMCO", RC )
          ENDIF
       ENDIF

       !--------------------------------------------------------------
       ! Test for convection timestep
       !--------------------------------------------------------------
       IF ( ITS_TIME_FOR_CONV() .and. notDryRun ) THEN

#ifdef BPCH_DIAG
          ! Increment the convection timestep
          CALL Set_Ct_Conv( INCREMENT=.TRUE. )
#endif

          !===========================================================
          !      ***** M I X E D   L A Y E R   M I X I N G *****
          !===========================================================
          IF ( Input_Opt%useTimers ) THEN
             CALL Timer_Start( "Boundary layer mixing", RC )
          ENDIF

          ! Note: mixing routine expects tracers in v/v
          ! DO_MIXING applies the tracer tendencies (dry deposition,
          ! emission rates) to the tracer arrays and performs PBL
          ! mixing.
          ! In the non-local PBL scheme, dry deposition and emission
          ! fluxes below the PBL are handled within the PBL mixing
          ! routine. Otherwise, tracer concentrations are first updated
          ! and the full-mixing is then applied.
          ! (ckeller, 3/5/15)
          ! NOTE: Tracer concentration units are converted locally
          ! to [v/v dry air] for mixing. Eventually mixing should
          ! be updated to use [kg/kg total air] (ewl, 9/18/15)
          CALL Do_Mixing( Input_Opt,  State_Chm, State_Diag, &
                          State_Grid, State_Met, RC )

          ! Trap potential errors
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = 'Error encountered in "Do_Mixing"!'
             CALL ERror_Stop( ErrMsg, ThisLoc )
          ENDIF

          IF ( Input_Opt%useTimers ) THEN
             CALL Timer_End( "Boundary layer mixing", RC )
          ENDIF

          IF ( prtDebug ) CALL Debug_Msg( '### MAIN: a TURBDAY:2' )

          !===========================================================
          !        ***** C L O U D   C O N V E C T I O N *****
          !===========================================================
          IF ( LCONV ) THEN
             IF ( Input_Opt%useTimers ) THEN
                CALL Timer_Start( "Convection", RC )
             ENDIF

             ! Call the appropriate convection routine
             ! NOTE: Tracer concentration units are converted locally
             ! to [kg/kg total air] for convection (ewl, 9/18/15)
             CALL Do_Convection( Input_Opt,  State_Chm, State_Diag, &
                                 State_Grid, State_Met, RC )

             ! Trap potential errors
             IF ( RC /= GC_SUCCESS ) THEN
                ErrMsg = 'Error encountered in "Do_Convection"!'
                CALL Error_Stop( ErrMsg, ThisLoc )
             ENDIF

             IF ( prtDebug ) THEN
                CALL Debug_Msg( '### MAIN: a CONVECTION' )
             ENDIF

             IF ( Input_Opt%useTimers ) THEN
                CALL Timer_End( "Convection", RC )
             ENDIF
          ENDIF

       ENDIF

       !==============================================================
       !               ***** C H E M I S T R Y *****
       !============================================================== 
       IF ( notDryRun ) THEN
          IF ( Input_Opt%useTimers ) THEN
             CALL Timer_Start( "All chemistry", RC )
          ENDIF

          ! Get the overhead column O3 for use with FAST-J
          ! NOTE: Move to CHEMISTRY section.  This now has to come after
          ! the call to HEMCO emissions driver EMISSIONS_RUN. (bmy, 3/20/15)
          CALL Get_Overhead_O3_For_FastJ( Input_Opt, State_Grid, State_Met )

          ! Every chemistry timestep...
          IF ( ITS_TIME_FOR_CHEM() ) THEN

#ifdef BPCH_DIAG
             ! Increment chemistry timestep counter
             CALL Set_Ct_Chem( INCREMENT=.TRUE. )
#endif

             ! SDE 05/28/13: Set H2O to State_Chm tracer if relevant
             IF ( ITS_A_FULLCHEM_SIM .and. id_H2O > 0 ) THEN
                CALL Set_H2O_Trac( (.not. LUCX), &
                                   Input_Opt , State_Chm,    &
                                   State_Grid, State_Met, RC )

                ! Trap potential errors
                IF ( RC /= GC_SUCCESS ) THEN
                   ErrMsg = 'Error encountered in "Set_H2O_Trac" #2!'
                   CALL Error_Stop( ErrMsg, ThisLoc )
                ENDIF
             ENDIF

             ! Do GEOS-Chem chemistry
             ! NOTE: Tracer concentration units are converted locally
             ! to [kg] for all of chemistry. We will replace use of [kg]
             ! once FlexChem is implemented (ewl, 9/18/15)
             CALL Do_Chemistry( Input_Opt,  State_Chm, State_Diag, &
                                State_Grid, State_Met, RC )

             ! Trap potential errors
             IF ( RC /= GC_SUCCESS ) THEN
                ErrMsg = 'Error encountered in "Do_Chemistry"!'
                CALL Error_Stop( ErrMsg, ThisLoc )
             ENDIF

          ENDIF

          IF ( Input_Opt%useTimers ) THEN
             CALL Timer_End( "All chemistry", RC )
          ENDIF
       ENDIF

       !==============================================================
       ! ***** W E T   D E P O S I T I O N  (rainout + washout) *****
       !==============================================================
       IF ( LWETD .and. ITS_TIME_FOR_DYN() .and. notDryRun ) THEN

          IF ( Input_Opt%useTimers ) THEN
             CALL Timer_Start( "Wet deposition", RC )
          ENDIF

          ! Do wet deposition
          ! NOTE: Tracer concentration units are converted locally
          ! to [kg/m2] in wet deposition to enable calculations
          ! along the column (ewl, 9/18/15)
          CALL Do_WetDep( Input_Opt, State_Chm, State_Diag, State_Grid, &
                          State_Met, RC )

          ! Trap potential errors
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = 'Error encountered in "Do_WetDep"!'
             CALL Error_Stop( ErrMsg, ThisLoc )
          ENDIF

          IF ( Input_Opt%useTimers ) THEN
             CALL Timer_End( "Wet deposition", RC )
          ENDIF

       ENDIF

       !==============================================================
       !      ***** U P D A T E   O P T I C A L   D E P T H *****
       !==============================================================
       IF ( ITS_TIME_FOR_CHEM() .and. notDryRun ) THEN

          IF ( Input_Opt%useTimers ) THEN
             CALL Timer_Start( "All chemistry",       RC )
             CALL Timer_Start( "=> All aerosol chem", RC )
          ENDIF

          ! Recalculate the optical depth at the wavelength(s) specified
          ! in the Radiation Menu. This must be done before the call to any
          ! diagnostic and only on a chemistry timestep.
          ! (skim, 02/05/11)
          CALL Recompute_OD( Input_Opt,  State_Chm, State_Diag, &
                             State_Grid, State_Met, RC )

          ! Trap potential errors
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = 'Error encountered in "Recompute_OD"!'
             CALL Error_Stop( ErrMsg, ThisLoc )
          ENDIF

          IF ( Input_Opt%useTimers ) THEN
             CALL Timer_End( "All chemistry",       RC )
             CALL Timer_End( "=> All aerosol chem", RC )
          ENDIF
       ENDIF

#ifdef RRTMG
       !==============================================================
       !  ***** R R T M G   R A D I A T I V E   T R A N S F E R *****
       !
       ! NOTE: Tracer concentration units are converted locally to
       ! [kg] in RRTMG. Units should eventually be [kg/kg]
       ! (ewl, 9/18/15)
       !==============================================================
       IF ( Input_opt%LRAD  .and. ITS_TIME_FOR_RT() .and.  &
            notDryRun ) THEN

          IF ( Input_Opt%useTimers ) THEN
             CALL Timer_Start( "RRTMG", RC )
          ENDIF

          ! Splash page
          IF ( am_I_Root .and. FIRST_RT ) THEN
             WRITE( 6, '(a)' ) REPEAT( '#', 79 )
             WRITE( 6, 500 ) 'R R T M G : Radiative Transfer Model (by AER)'
500          FORMAT( '#####', 12x, a, 12x, '#####' )
             WRITE( 6, '(a)' ) REPEAT( '#', 79 )
             FIRST_RT = .FALSE.
          ENDIF

          iSeed = iSeed + 15

          !-----------------------------------------------------------
          ! Determine if we are doing clear-sky or all-sky
          !
          ! Clear-sky is output with all-sky, so we just need
          ! to run once regardless of whether both are required
          ! or just one.
          !-----------------------------------------------------------
          IF ( Input_Opt%LSKYRAD(2) ) THEN
             iCld = 1
          ELSE
             iCld = 0         !clouds are on
          ENDIF

          !-----------------------------------------------------------
          ! Calculation for each of the potential output types
          ! See: wiki.geos-chem.org/Coupling_GEOS-Chem_with_RRTMG
          !
          ! Flux outputs (scheduled in HISTORY.rc):
          !  0-BA  1=O3  2=ME  3=SU   4=NI  5=AM
          !  6=BC  7=OA  8=SS  9=DU  10=PM  11=ST (UCX only)
          !-----------------------------------------------------------
          DO N = 1, State_Diag%nRadFlux

             ! Index number for RRTMG (see list above)
             iSpecMenu = State_Diag%RadFluxInd(N)

             ! Slot # of netCDF diagnostic arrays to update
             iNcDiag = N

             ! Echo info
             WRITE( 6, 520 ) State_Diag%RadFluxName(N), iSpecMenu
520          FORMAT( 5x, '- Calling RRTMG to compute flux: ', &
                     a2, ' (Index = ', i2.2, ')' )

             ! Generate mask for species in RT
             CALL Set_SpecMask( iSpecMenu )

             ! Compute radiative fluxes for the given output
             CALL Do_RRTMG_Rad_Transfer( ThisDay    = Day,        &
                                         ThisMonth  = Month,      &
                                         iCld       = iCld,       &
                                         iSpecMenu  = iSpecMenu,  &
                                         iNcDiag    = iNcDiag,    &
                                         iSeed      = iSeed,      &
                                         Input_Opt  = Input_Opt,  &
                                         State_Chm  = State_Chm,  &
                                         State_Diag = State_Diag, &
                                         State_Grid = State_Grid, &
                                         State_Met  = State_Met,  &
                                         RC         = RC          )

             ! Trap potential errors
             IF ( RC /= GC_SUCCESS ) THEN
                ErrMsg = 'Error encountered in "Do_RRTMG_Rad_Transfer", ' // &
                         'for flux output = ' // &
                         TRIM( State_Diag%RadFluxName(N) )
                CALL Error_Stop( ErrMsg, ThisLoc )
             ENDIF
          ENDDO

#ifdef BPCH_DIAG
          ! Increment radiation timestep counter
          CALL Set_Ct_Rad( INCREMENT=.TRUE. )
#endif

          IF ( prtDebug ) THEN
             CALL Debug_Msg( '### MAIN: a DO_RRTMG_RAD_TRANSFER' )
          ENDIF

          IF ( Input_Opt%useTimers ) THEN
             CALL Timer_End( "RRTMG", RC )
          ENDIF
       ENDIF
#endif

       !==============================================================
       !      ***** D I A G N O S T I C S   A R C H I V A L *****
       !==============================================================
       IF ( notDryRun ) THEN

          !-----------------------------------------------------------
          !        ***** H I S T O R Y   U P D A T E  *****
          !-----------------------------------------------------------
          IF ( Input_Opt%useTimers ) THEN
             CALL Timer_Start( "All diagnostics",           RC )
             CALL Timer_Start( "Output",                    RC )
             CALL Timer_Start( "=> History (netCDF diags)", RC )
          ENDIF

          ! Set State_Diag arrays that rely on state at end of timestep
          CALL Set_Diagnostics_EndofTimestep( Input_Opt,  State_Chm,  &
                                              State_Diag, State_Grid, &
                                              State_Met,  RC )

          ! Trap potential errors
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = 'Error encountered in "Set_Diagnostics_EndOfTimestep"!'
             CALL Error_Stop( ErrMsg, ThisLoc )
          ENDIF

          ! Archive aerosol mass and PM2.5 diagnostics
          IF ( State_Diag%Archive_AerMass ) THEN
             CALL Set_AerMass_Diagnostic( Input_Opt,  State_Chm, State_Diag, &
                                          State_Grid, State_Met, RC )
          ENDIF

          ! Trap potential errors
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = 'Error encountered in "Set_AerMass_Diagnostic"!'
             CALL Error_Stop( ErrMsg, ThisLoc )
          ENDIF

          ! Increment the timestep values by the heartbeat time
          ! This is because we need to write out data with the timestamp
          ! at the end of the heartbeat timestep (i.e. at end of run)
          !
          ! NOTE: This should now go before HISTORY_UPDATE, so that we
          ! can recompute the update alarm interval properly for monthly
          ! or yearly intervals spanning leap years. (bmy, 3/5/19)
          CALL History_SetTime( Input_Opt, RC )

          ! Trap potential errors
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = 'Error encountered in "History_SetTime"!'
             CALL Error_Stop( ErrMsg, ThisLoc )
          ENDIF

          ! Update each HISTORY ITEM from its data source
          CALL History_Update( Input_Opt, RC )

          ! Trap potential errors
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = 'Error encountered in "History_Update"!'
             CALL Error_Stop( ErrMsg, ThisLoc )
          ENDIF

          IF ( Input_Opt%useTimers ) THEN
             CALL Timer_End( "All diagnostics",           RC )
             CALL Timer_End( "Output",                    RC )
             CALL Timer_End( "=> History (netCDF diags)", RC )
          ENDIF

          !-----------------------------------------------------------
          !     ***** O B S P A C K   D I A G N O S T I C S *****
          !-----------------------------------------------------------
          IF ( Input_Opt%Do_ObsPack ) THEN

             IF ( Input_Opt%useTimers ) THEN
                CALL Timer_Start( "All diagnostics",        RC )
                CALL Timer_Start( "=> ObsPack diagnostics", RC )
             ENDIF

             ! Sample the observations in today's ObsPack file
             CALL ObsPack_Sample( NYMD, NHMS, Input_Opt,  State_Chm, &
                                  State_Diag, State_Grid, State_Met, RC )

             IF ( Input_Opt%useTimers ) THEN
                CALL Timer_End( "All diagnostics",        RC )
                CALL Timer_End( "=> ObsPack diagnostics", RC )
             ENDIF

          ENDIF

          !--------------------------------------------------------------
          !   ***** P L A N E F L I G H T   D I A G   S E T U P  *****
          !--------------------------------------------------------------
          IF ( Input_Opt%Do_Planeflight .and. ITS_A_NEW_DAY() ) THEN

             IF ( Input_Opt%useTimers ) THEN
                CALL Timer_Start( "All diagnostics",             RC) 
                CALL Timer_Start( "Output",                      RC)
                CALL Timer_Start( "=> Binary punch diagnostics", RC)
             ENDIF

             ! Initialize planeflight diagnostic
             CALL Setup_PlaneFlight( Input_Opt, State_Chm, State_Grid, &
                                     State_Met, RC )

             ! Trap potential errors
             IF ( RC /= GC_SUCCESS ) THEN
                ErrMsg = 'Error encountered in "Setup_Planeflight"!'
                CALL Error_Stop( ErrMsg, ThisLoc )
             ENDIF

             IF ( Input_Opt%useTimers ) THEN
                CALL Timer_End( "All diagnostics",             RC )
                CALL Timer_End( "Output",                      RC )
                CALL Timer_End( "=> Binary punch diagnostics", RC )
             ENDIF

          ENDIF

          !-----------------------------------------------------------
          ! *** C H 4   S I M U L A T I O N   D I A G N O S I C S ***
          !
          ! CH4 columns from the GOSAT instrument (mps, 6/16/17)
          !-----------------------------------------------------------
          IF ( Input_Opt%GOSAT_CH4_OBS ) THEN
             IF ( ITS_A_NEW_HOUR() ) THEN
                CALL CALC_GOSAT_CH4_FORCE( Input_Opt, State_Chm, State_Grid, &
                                           State_Met )
             ENDIF
          ENDIF

          !--------------------------------------------------------------
          ! **** C H 4   S I M U L A T I O N   D I A G N O S I C S ****
          !
          ! CH4 columns from the TCCON instrument (mps, 8/17/17)
          !--------------------------------------------------------------
          IF ( Input_Opt%TCCON_CH4_OBS ) THEN
             IF ( ITS_A_NEW_HOUR() ) THEN
                CALL CALC_TCCON_CH4_FORCE( Input_Opt, State_Chm, State_Grid, &
                                           State_Met )
             ENDIF
          ENDIF
       ENDIF

       !==============================================================
       !   ***** I N C R E M E N T   E L A P S E D   T I M E *****
       !
       ! Moved before diagnostics to count the last timestep as done.
       ! Need to save timestamps for filenames. (ccc, 5/13/09)
       !==============================================================
       CALL Timestamp_Diag()
       CALL Set_Elapsed_Sec()
       CALL Set_Current_Time()
       IF ( prtDebug ) THEN
          CALL Debug_Msg( '### MAIN: after SET_ELAPSED_SEC' )
       ENDIF

       IF ( notDryRun ) THEN
          !===========================================================
          !   ***** D I A G N O S T I C S   A R C H I V A L *****
          !
          !             ***** C O N T I N U E D *****
          !===========================================================
          IF ( Input_Opt%useTimers ) THEN
             CALL Timer_Start( "All diagnostics",              RC )
             CALL Timer_Start( "Output",                       RC )
          ENDIF

          !-----------------------------------------------------------
          !  ***** P L A N E F L I G H T   D I A G N O S T I C  *****
          !-----------------------------------------------------------
          IF ( Input_Opt%Do_Planeflight ) THEN
             ! Archive data along the flight track
             CALL PLANEFLIGHT( Input_Opt,  State_Chm, State_Diag, &
                               State_Grid, State_Met, RC )

             ! Trap potential errors
             IF ( RC /= GC_SUCCESS ) THEN
                ErrMsg = 'Error encountered in "Planeflight"!'
                CALL Error_Stop( ErrMsg, ThisLoc )
             ENDIF
          ENDIF
          IF ( prtDebug ) CALL Debug_Msg( '### MAIN: after Planeflight' )

#ifdef BPCH_DIAG
          IF ( Input_Opt%useTimers ) THEN
             CALL Timer_Start( "=> Binary punch diagnostics",  RC )
          ENDIF
 
          !-----------------------------------------------------------
          !  *** A R C H I V E   B P C H   D I A G N O S T I C S ***
          !-----------------------------------------------------------
          IF ( ITS_TIME_FOR_DIAG() ) THEN

             IF ( prtDebug ) CALL Debug_Msg('### MAIN: b DIAGNOSTICS')

             ! Accumulate several diagnostic quantities
             CALL Diag1( Input_Opt, State_Chm, State_Grid, State_Met, RC )

             ! Trap potential errors
             IF ( RC /= GC_SUCCESS ) THEN
                ErrMsg = 'Error encountered in "Diag1"!'
                CALL Error_Stop( ErrMsg, ThisLoc )
             ENDIF
             IF ( prtDebug ) CALL Debug_Msg( '### MAIN: after DIAG1' )

             ! Increment diagnostic timestep counter. (ccc, 5/13/09)
             CALL Set_Ct_Diag( INCREMENT=.TRUE. )

             ! Planeflight diagnostic moved to be after chemistry, kyu
             IF ( prtDebug ) CALL Debug_Msg('### MAIN: a DIAGNOSTICS')
          ENDIF

          !------------------------------------------------------------
          !  ***** T I M E S E R I E S   D I A G N O S T I C S  *****
          !------------------------------------------------------------

          ! Morning or afternoon timeseries
          IF ( Input_Opt%DO_ND51 ) THEN
             CALL DIAG51( Input_Opt, State_Chm, State_Grid, State_Met, RC )

             ! Trap potential errors
             IF ( RC /= GC_SUCCESS ) THEN
                ErrMsg = 'Error encountered in "Diag51"!'
                CALL Error_Stop( ErrMsg, ThisLoc )
             ENDIF
          ENDIF

          IF ( Input_Opt%DO_ND51b ) THEN
             CALL DIAG51b( Input_Opt, State_Chm, State_Grid, State_Met, RC )

             ! Trap potential errors
             IF ( RC /= GC_SUCCESS ) THEN
                ErrMsg = 'Error encountered in "Diag51b"!'
                CALL Error_Stop( ErrMsg, ThisLoc )
             ENDIF
          ENDIF
          IF ( prtDebug ) CALL Debug_Msg( '### MAIN: after DIAG51' )

          IF ( prtDebug ) CALL Debug_Msg('### MAIN: after TIMESERIES')

          IF ( Input_Opt%useTimers ) THEN
             CALL Timer_End( "=> Binary punch diagnostics",  RC )
          ENDIF
#endif
          IF ( Input_Opt%useTimers ) THEN
             CALL Timer_End( "All diagnostics",              RC )
             CALL Timer_End( "Output",                       RC )
          ENDIF
       ENDIF

       !=================================================================
       !     ***** C O P Y   I - 3   F I E L D S *****
       !
       !     The I-3 fields at the end of an outer timestep (every 3 hours)
       !     become the fields at the beginning of the next timestep.
       !     This update must occur before writing History. (wbd1 12/03/19)   
       !=================================================================
       IF ( notDryRun ) THEN
          IF ((mod(get_hour(), 3) .eq. 0) .AND. (get_minute() .eq. 0)) THEN
             CALL Copy_I3_Fields( State_Met )
             IF ( prtDebug ) THEN
                CALL Debug_Msg( '### MAIN: after COPY_I3_FIELDS' )
             ENDIF
          ENDIF
       ENDIF

       !--------------------------------------------------------------
       !            ***** H I S T O R Y   W R I T E *****
       !--------------------------------------------------------------
       IF ( notDryRun ) THEN
          IF ( Input_Opt%useTimers ) THEN
             CALL Timer_Start( "All diagnostics",           RC )
             CALL Timer_Start( "Output",                    RC )
             CALL Timer_Start( "=> History (netCDF diags)", RC )
          ENDIF

          ! Write HISTORY ITEMS in each diagnostic collection to disk
          ! (or skip writing if it is not the proper output time.
          CALL History_Write( Input_Opt, State_Chm%Spc_Units, RC )

          ! Trap potential errors
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = 'Error encountered in "History_Write"!'
             CALL Error_Stop( ErrMsg, ThisLoc )
          ENDIF

          IF ( Input_Opt%useTimers ) THEN
             CALL Timer_End( "All diagnostics",           RC )
             CALL Timer_End( "Output",                    RC )
             CALL Timer_End( "=> History (netCDF diags)", RC )
          ENDIF
       ENDIF

       !==============================================================
       !  ***** E N D   O F   D Y N A M I C   T I M E S T E P *****
       !==============================================================
    ENDDO
  ENDDO

  !=================================================================
  !         ***** C L E A N U P   A N D   Q U I T *****
  !=================================================================
9999 CONTINUE

  ! Skip operations when running in dry-run mode
  IF ( notDryRun ) THEN

     IF ( Input_Opt%useTimers ) THEN
        CALL Timer_End ( "Timesteps", RC )
     ENDIF

    !--------------------------------------------------------------
    !     ***** W R I T E   H E M C O   R E S T A R T S *****
    !
    ! NOTE: If BPCH_DIAG=y, then this is done above whenever
    ! ITS_TIME_FOR_BPCH is TRUE.   If BPCH_DIAG=n, then we have to
    ! add this here to make sure we get HEMCO restart otuput.
    !--------------------------------------------------------------
     IF ( Input_Opt%useTimers ) THEN
        CALL Timer_Start( "HEMCO",  RC )
        CALL Timer_Start( "Output", RC )
     ENDIF

    ! Force the output of a HEMCO restart file (ckeller, 4/1/15)
    CALL HCOI_GC_WriteDiagn( Input_Opt, .TRUE., RC )

    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered in "HCOI_GC_WriteDiagn"!'
       CALL Error_Stop( ErrMsg, ThisLoc )
    ENDIF

    IF ( Input_Opt%useTimers ) THEN
       CALL Timer_End( "HEMCO",  RC )
       CALL Timer_End( "Output", RC )
    ENDIF

    !-----------------------------------------------------------------
    !         ***** O B S P A C K   D I A G N O S T I C S *****
    !
    ! Flush any unwritten ObsPack data to disk and finalize
    !-----------------------------------------------------------------
    IF ( Input_Opt%Do_ObsPack ) THEN

       IF ( Input_Opt%useTimers ) THEN
          CALL Timer_Start( "All diagnostics",        RC )
          CALL Timer_Start( "=> ObsPack diagnostics", RC )
       ENDIF

       ! Write any remaining ObsPack data to disk, and immediately
       ! thereafter free the ObsPack pointer fields of State_Diag
       IF ( ASSOCIATED( State_Diag%ObsPack_id ) ) THEN
          CALL ObsPack_Write_Output( Input_Opt, State_Diag, RC )
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = 'Error encountered in "ObsPack_Write_Output"!'
             CALL Error_Stop( ErrMsg, ThisLoc )
          ENDIF
       ENDIF

       ! Free the ObsPack species mapping fields of State_Diag
       CALL ObsPack_SpeciesMap_Cleanup( Input_Opt, State_Diag, RC )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Error encountered in "ObsPack_SpeciesMap_Cleanup"!'
          CALL Error_Stop( ErrMsg, ThisLoc )
       ENDIF

       IF ( Input_Opt%useTimers ) THEN
          CALL Timer_End( "All diagnostics",        RC )
          CALL Timer_End( "=> ObsPack diagnostics", RC )
       ENDIF
    ENDIF

    !--------------------------------------------------------------
    ! Print the mass-weighted mean OH concentration (if applicable)
    !--------------------------------------------------------------
    CALL Print_Diag_OH( Input_Opt, RC )

    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered in "Print_Diag_OH"!'
       CALL Error_Stop( ErrMsg, ThisLoc )
    ENDIF
  ENDIF

  !-----------------------------------------------------------------
  ! Finalize GEOS-Chem
  !-----------------------------------------------------------------
  IF ( Input_Opt%useTimers ) THEN
     CALL Timer_Start( "Finalization", RC )
  ENDIF

  ! Cleanup the dry-run simulation (if necessary)
  CALL Cleanup_Dry_Run( Input_Opt, RC )

  ! Close all files
  CALL CLOSE_FILES()
  IF ( prtDebug ) CALL Debug_Msg( '### MAIN: a CLOSE_FILES' )

  !%%% NOTE: Call HISTORY_CLEANUP from cleanup.F.  This will
  !%%% close all netCDF files upon both normal or abnormal exits.

  ! Deallocate fields of the Chemistry State object
  CALL Cleanup_State_Chm( State_Chm, RC )
  IF ( RC /= GC_SUCCESS ) THEN
     ErrMsg = 'Error encountered in "Cleanup_State_Chm"!'
     CALL Error_Stop( ErrMsg, ThisLoc )
  ENDIF
  IF ( prtDebug ) CALL Debug_Msg( '### MAIN: a cleanup State_Chm' )

  ! Deallocate fields of the Diagnostics State object
  CALL Cleanup_State_Diag( State_Diag, RC )
  IF ( RC /= GC_SUCCESS ) THEN
     ErrMsg = 'Error encountered in "Cleanup_State_Diag"!'
     CALL Error_Stop( ErrMsg, ThisLoc )
  ENDIF
  IF ( prtDebug ) CALL Debug_Msg( '### MAIN: a cleanup State_Diag' )

  ! Deallocate fields of the Meteorology State object
  CALL Cleanup_State_Met( State_Met, RC )
  IF ( RC /= GC_SUCCESS ) THEN
     ErrMsg = 'Error encountered in "Cleanup_State_Met"!'
     CALL Error_Stop( ErrMsg, ThisLoc )
  ENDIF
  IF ( prtDebug ) CALL Debug_Msg( '### MAIN: a cleanup State_Met' )

  ! Deallocate dynamic module arrays
  CALL CleanUp( Input_Opt, State_Grid, .FALSE., RC )
  IF ( RC /= GC_SUCCESS ) THEN
     ErrMsg = 'Error encountered in "Cleanup"!'
     CALL Error_Stop( ErrMsg, ThisLoc )
  ENDIF
  IF ( prtDebug ) CALL Debug_Msg( '### MAIN: a cleanup modules' )

  ! Deallocate fields of the Input Options object
  CALL Cleanup_Input_Opt( Input_Opt, RC )
  IF ( RC /= GC_SUCCESS ) THEN
     ErrMsg = 'Error encountered in "Cleanup_Input_Opt"!'
     CALL Error_Stop( ErrMsg, ThisLoc )
  ENDIF
  IF ( prtDebug ) CALL Debug_Msg( '### MAIN: a cleanup Input_Opt' )

  ! Deallocate fields of the Grid State object
  CALL Cleanup_State_Grid( State_Grid, RC )
  IF ( RC /= GC_SUCCESS ) THEN
     ErrMsg = 'Error encountered in "Cleanup_State_grid"!'
     CALL Error_Stop( ErrMsg, ThisLoc )
  ENDIF
  IF ( prtDebug ) CALL Debug_Msg( '### MAIN: a cleanup State_Grid' )

  ! Deallocate fields of the diagnostics list object
  IF ( notDryRun ) THEN
     CALL Cleanup_DiagList( Diag_List, RC )
     IF ( RC /= GC_SUCCESS ) THEN
        ErrMsg = 'Error encountered in "Cleanup_DiagList"!'
        CALL Error_Stop( ErrMsg, ThisLoc )
     ENDIF
     CALL Cleanup_TaggedDiagList( TaggedDiag_List, RC )
     IF ( RC /= GC_SUCCESS ) THEN
        ErrMsg = 'Error encountered in "Cleanup_TaggedDiagList"!'
        CALL Error_Stop( ErrMsg, ThisLoc )
     ENDIF
     IF ( prtDebug ) THEN
        CALL Debug_Msg( '### MAIN: a cleanup diag lists' )
     ENDIF
  ENDIF

#ifdef APM
  ! Clean up arrays for APM microphysics, etc.
  CALL CLEANUP_APMARRAYS()
  CALL CLEANUP_APM3D( Input_Opt, RC )
#endif

  !-----------------------------------------------------------------------------
  ! TWO-WAY NESTING OPTION
  ! This is only invoked when compiling GEOS-Chem with COUPLE=y
  !
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! %%%%% NOTE: THIS OPTION WAS PROBABLY BROKEN WITH THE ADDITION OF HEMCO %%%%%
  ! %%%%% AND FLEXCHEM INTO GEOS-CHEM V10-01 AND v11-01.  BUYER BEWARE.    %%%%%
  ! %%%%%  --  Bob Yantosca (22 Jan 2018)                                  %%%%%
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#ifdef EXCHANGE
  ! Finalize the two-way nesting
  CALL Cleanup_Exchange( Input_Opt )
#endif
!-----------------------------------------------------------------------------

  ! Free the shadow variables in error_mod.F
  CALL Cleanup_Error()

#ifdef GTMM_Hg
  ! Deallocate arrays from GTMM model for mercury simulation
  IF ( LGTMM ) CALL CleanupCASAarrays()
#endif

  IF ( prtDebug ) CALL Debug_Msg( '### MAIN: a CLEANUP' )

  IF ( Input_Opt%useTimers ) THEN
     ! Stop remaining timers
     CALL Timer_End( "Finalization", RC )
     CALL Timer_End( "GEOS-Chem",    RC )

     ! Print timer output (skip if a dry-run)
     IF ( notDryRun ) THEN
        CALL Timer_PrintAll( Input_Opt, RC )
     ENDIF
  ENDIF

  ! Print ending time of simulation
  CALL Display_End_Time()

  ! Flush the buffer to get output
  CALL Flush( 6 )

CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: display_model_info
!
! !DESCRIPTION: Internal Subroutine DISPLAY\_MODEL\_INFO displays the
!  appropriate messages for the given model and machine type.  It also
!  prints the starting time and date (local time) of the GEOS-Chem simulation.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Display_Model_Info()
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! For system time stamp
    CHARACTER(LEN=16) :: STAMP

#include "gc_classic_version.H"

    WRITE(  6, '(a)' ) &
         REPEAT( '*', 13 )                                          // &
         '   S T A R T I N G   G E O S - C H E M   '                // &
         REPEAT( '*', 13 )

    !-----------------------------------------------------------------
    ! Mode of simulation
    !-----------------------------------------------------------------
    WRITE( 6, 100 ) 'GEOS-Chem "Classic"'

    !-----------------------------------------------------------------
    ! Print model version
    !-----------------------------------------------------------------
    WRITE( 6, 110 ) TRIM( GC_CLASSIC_VERSION )

    !-----------------------------------------------------------------
    ! Print compiler
    !-----------------------------------------------------------------
#if defined( LINUX_IFORT )
    WRITE( 6, 120  ) 'Intel Fortran Compiler (aka ifort)'
#elif defined( LINUX_GFORTRAN )
    WRITE( 6, 120 ) 'GNU Fortran compiler (aka gfortran)'
#endif

    !-----------------------------------------------------------------
    ! Print status of OpenMP
    !-----------------------------------------------------------------
#ifdef NO_OMP
    WRITE( 6, 150 ) 'OFF'
#else
    WRITE( 6, 150 ) 'ON'
#endif

    !-----------------------------------------------------------------
    ! Print status of binary punch (bpch) diagnostics
    !-----------------------------------------------------------------
#ifdef BPCH_DIAG
    WRITE( 6, 160 ) 'ON'
#else
    WRITE( 6, 160 ) 'OFF'
#endif

    !-----------------------------------------------------------------
    ! Print status of netCDF diagnostics (aka History) - always on
    !-----------------------------------------------------------------
    WRITE( 6, 170 ) 'ON'

    !-----------------------------------------------------------------
    ! Print msg if netCDF compression is supported
    !-----------------------------------------------------------------
#ifdef NC_HAS_COMPRESSION
    WRITE( 6, 180 ) 'SUPPORTED'
#else
    WRITE( 6, 180 ) 'NOT SUPPORTED (or shut off w/ NC_NODEFLATE=y)'
#endif

    !-----------------------------------------------------------------
    ! Print msg if Luo et al wetdep scheme is supported
    !-----------------------------------------------------------------
#ifdef LUO_WETDEP
    WRITE( 6, 185 ) 'ON'
#else
    WRITE( 6, 185 ) 'OFF'
#endif

    !-----------------------------------------------------------------
    ! System time stamp
    !-----------------------------------------------------------------
    STAMP = SYSTEM_TIMESTAMP()
    WRITE( 6, 190 ) STAMP

    !-----------------------------------------------------------------
    ! Format strings
    !-----------------------------------------------------------------
100 FORMAT( /, '===> Mode of operation         : ', a             )
110 FORMAT(    '===> GEOS-Chem version         : ', a             )
120 FORMAT(    '===> Compiler                  : ', a             )
150 FORMAT(    '===> Parallelization w/ OpenMP : ', a             )
160 FORMAT(    '===> Binary punch diagnostics  : ', a             )
170 FORMAT(    '===> netCDF diagnostics        : ', a             )
180 FORMAT(    '===> netCDF file compression   : ', a             )
185 FORMAT(    '===> Luo et al (2019) wetdep?  : ', a             )
190 FORMAT( /, '===> SIMULATION START TIME: ',      a, ' <===', / )

  END SUBROUTINE Display_Model_Info
#ifdef BPCH_DIAG
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ctm_flush
!
! !DESCRIPTION: Internal subroutine CTM\_FLUSH flushes certain diagnostic
! file buffers to disk.
!\\
!\\
! CTM\_FLUSH should normally be called after each diagnostic output, so that 
! in case the run dies, the output files from the last diagnostic timestep 
! will not be lost.
!\\
!\\
! FLUSH is an intrinsic FORTRAN subroutine and takes as input the unit number
! of the file to be flushed to disk.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE CTM_Flush()
!
! !REVISION HISTORY:
!  31 Aug 2000 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC

    CALL FLUSH( IU_BPCH )

  END SUBROUTINE CTM_Flush
#endif
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: display_end_time
!
! !DESCRIPTION: Internal subroutine DISPLAY\_END\_TIME prints the ending
!  time of the GEOS-Chem simulation.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Display_End_Time()
!
! !REVISION HISTORY:
!  03 May 2005 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    CHARACTER(LEN=16) :: STAMP

    ! Print system time stamp
    STAMP = SYSTEM_TIMESTAMP()
    WRITE( 6, 100 ) STAMP
100 FORMAT( /, '===> SIMULATION END TIME: ', a, ' <===', / )

    ! Echo info
    WRITE ( 6, 3000 )
3000 FORMAT( /, '**************   E N D   O F   G E O S -- C H E M   ', &
                '**************' )

  END SUBROUTINE Display_End_Time
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_overhead_o3_for_fastj
!
! !DESCRIPTION: Internal subroutine GET\_OVERHEAD\_O3\_FOR\_FASTJ
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Get_Overhead_O3_For_FastJ( Input_Opt, State_Grid, State_Met )
!
! !USES:
!
    USE Input_Opt_Mod,      ONLY : OptInput
    USE State_Grid_Mod,     ONLY : GrdState
    USE State_Met_Mod,      ONLY : MetState
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt
    TYPE(GrdState), INTENT(IN)    :: State_Grid
    TYPE(MetState), INTENT(IN)    :: State_Met
!
! !REMARKS:
!  This routine makes use of variables declared in above in the main program
!  (which are visible in all sub-programs below the CONTAINS statement).
!                                                                             .
!  The original code was done in FAST-J routine "set_prof.F", but has been
!  split off to facilitate development of the grid-independent model.
!
! !REVISION HISTORY:
!  07 Mar 2012 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! FAST-J is only used for fullchem and offline aerosol, skip otherwise
    IF ( ITS_A_FULLCHEM_SIM .or. ITS_AN_AEROSOL_SIM  ) THEN

       ! Only execute this if we are doing chemistry
       ! and if it we are at a chemistry timestep
       IF ( LCHEM .and. ITS_TIME_FOR_CHEM() ) THEN

          ! Get the overhead O3 column for FAST-J.  Take either the
          ! TOMS O3 data or the column O3 directly from the met fields
          CALL Compute_Overhead_O3( Input_Opt, State_Grid, DAY, &
                                    Input_Opt%USE_O3_FROM_MET,  &
                                    State_Met%TO3 )
       ENDIF
    ENDIF

  END SUBROUTINE Get_Overhead_O3_For_FastJ
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Init_Dry_Run
!
! !DESCRIPTION: Looks at the input arguments to determine if the user
!  has selected to do a GEOS-Chem dry-run.  If so, then the proper
!  fields of Input\_Opt will be populated accordingly, and the dry-run
!  log file will be opened.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Init_Dry_Run( Input_Opt, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod,     ONLY : OptInput
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(INOUT) :: Input_Opt   ! Input Options Object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REMARKS:
!  If in a "dry-run" mode, GEOS-Chem will simply check whether files
!  are present (and possibly in the correct format) and go through
!  time-steps to check met fields and other IO issues.
!  No actual "compute" is performed.
!
! !REVISION HISTORY:
!  13 Nov 2019 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER            :: nArg,   ArgLen

    ! Strings
    CHARACTER(LEN=255) :: ArgVal, ErrMsg, ThisLoc

    !=================================================================
    ! Init_Dry_Run begins here!
    !=================================================================

    ! Initialize
    RC      = GC_SUCCESS
    nArg    = 0
    ErrMsg  = ''
    ThisLoc = ' -> at Init_Dry_Run (in GeosCore/main.F90)'

    !=================================================================
    ! Parse arguments to determine if the dry-run has been selected
    !=================================================================
    DO

       ! Initialize for next argument
       ArgLen  = 0
       ArgVal  = ''

       ! Get the next argument
       CALL Get_Command_Argument( nArg, ArgVal, ArgLen )
       IF ( ArgLen == 0 ) EXIT

       ! Parse the arguments
       SELECT CASE( TRIM( ArgVal ) )

       ! Test for the dry-run switch
       CASE( '--dryrun' )
          Input_Opt%DryRun  = .TRUE.

       ! Otherwise pass through
       CASE DEFAULT
          ! pass

       END SELECT

       ! Increment the argument counter
       nArg = nArg + 1
    ENDDO

    !=================================================================
    ! If GEOS-Chem is running in dry-run mode
    ! then print a warning to both to stdout and the HEMCO log
    !=================================================================
    IF ( Input_Opt%DryRun ) THEN
       CALL Print_Dry_Run_Warning( 6 )
    ENDIF

  END SUBROUTINE Init_Dry_Run
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Cleanup_Dry_Run
!
! !DESCRIPTION: Looks at the input arguments to determine if the user
!  has selected to do a GEOS-Chem dry-run.  If so, then the proper
!  fields of Input\_Opt will be populated accordingly.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Cleanup_Dry_Run( Input_Opt, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE HCO_Interface_Mod, ONLY : HcoState
    USE Input_Opt_Mod,     ONLY : OptInput
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)  :: Input_Opt   ! Input Options Object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT) :: RC          ! Success or failure?
!
! !REMARKS:
!  Uses the intrinsic F2003 function Get_Command_Argument.
!
! !REVISION HISTORY:
!  13 Nov 2019 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC

    ! Initialize
    RC = GC_SUCCESS

    ! For dry-run simulations, Write the dry-run header to
    ! stdout (aka GEOS-Chem log file) and the HEMCO log file.
    IF ( Input_Opt%DryRun ) THEN
       CALL Print_Dry_Run_Warning( 6 )
       CALL Print_Dry_Run_Warning( HcoState%Config%Err%LUN )
    ENDIF

  END SUBROUTINE Cleanup_Dry_Run
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Print_Dry_Run_Warning
!
! !DESCRIPTION: Prints the warning for the GEOS-Chem dry run to either
!  stdout (aka the GC log file) and the dry-run log file.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Print_Dry_Run_Warning( U )
!
! !INPUT PARAMETERS:
!
    INTEGER, INTENT(IN) :: U   ! Logical unit number
!
! !REVISION HISTORY:
!  06 Jan 2015 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
    !=================================================================
    ! Print warning info to the desired file
    !=================================================================
    WRITE( U, 100 )
    WRITE( U, 100 ) REPEAT( '!', 79 )
    WRITE( U, 100 ) '!!! GEOS-CHEM IS IN DRY-RUN MODE!'
    WRITE( U, 100 ) '!!!'
    WRITE( U, 100 ) '!!! You will NOT get output for this run!'
    WRITE( U, 100 ) '!!! Use this command to validate a '         // &
                    'GEOS-Chem run configuration:'
    WRITE( U, 100 ) '!!!   ./geos --dryrun > log'
    WRITE( U, 100 ) '!!!'
    WRITE( U, 100 ) '!!! REMOVE THE --dryrun ARGUMENT FROM THE '   // &
                    'COMMAND LINE'
    WRITE( U, 100 ) '!!! BEFORE RUNNING A GEOS-Chem PRODUCTION '   // &
                    'SIMULATION!'
    WRITE( U, 100 ) REPEAT( '!', 79 )
    WRITE( U, 120 ) '!!! Start Date       : ', &
                    Input_Opt%NYMDb, Input_Opt%NHMSb
    WRITE( U, 120 ) '!!! End Date         : ', &
                    Input_Opt%NYMDe, Input_Opt%NHMSe
    WRITE( U, 110 ) '!!! Simulation       : ', &
                     TRIM(Input_Opt%SimulationName)
    WRITE( U, 110 ) '!!! Meteorology      : ', &
                     TRIM(Input_Opt%MetField )
    WRITE( U, 110 ) '!!! Grid Resolution  : ', &
                    TRIM(State_Grid%GridRes )
    WRITE( U, 100 ) REPEAT( '!', 79 )
    WRITE( U, 100 )

    ! Format statements
100 FORMAT( a                 )
110 FORMAT( a, a              )
120 FORMAT( a, i8.8, 1x, i6.6 )

  END SUBROUTINE Print_Dry_Run_Warning
!EOC
END PROGRAM GEOS_Chem
#endif
