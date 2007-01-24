! $Id: main.f,v 1.45 2007/01/24 18:22:23 bmy Exp $
! $Log: main.f,v $
! Revision 1.45  2007/01/24 18:22:23  bmy
! GEOS-Chem v7-04-11, includes the following modifications:
! - Fixed near-land lightning; now scale 2x25, 4x5 to 6 Tg N/yr (ltm,bmy)
! - Now allow ND49, ND50, ND51 to save transects of a single lon or lat (bmy)
! - Add case for T > 293K to routine GET_LWC in "mercury_mod.f" (cdh,bmy)
! - Bug fix: correct indices for embedded chemistry in "transport_mod.f" (phs)
! - Bug fix: correct typo in SEASALT_CHEM routine in "sulfate_mod.f" (bmy)
! - Now prevent seg fault errors when LBIOMASS=F (bmy)
! - Now fixed minor bug that inverted TROPP1 and TROPP2 (phs)
! - CMN_SIZE: now define LLTROP_FIX for GCAP. (phs)
! - a3_read_mod.f : added SNOW and GETWETTOP fields for GCAP (phs)
! - main.f: remove duplicate call for unzip in GCAP case (phs)
! - time_mod.f: fix leap year problem in get_time_ahead for GCAP (phs)
! - extra fixes for the variable tropopause (phs)
! - minor diagnostic updates (phs)
!
! Revision 1.42  2006/10/17 17:51:14  bmy
! GEOS-Chem v7-04-10, includes the following modifications:
! - Includes variable tropopause with ND54 diagnostic
! - Added GFED2 biomass emissions for SO2, NH3, BC, OC, CO2
! - Rewrote default biomass emissions routines for clarity
! - Updates for GCAP: future emissions, met-field reading, TOMS-O3
! - Bug fix in planeflight_mod.f: set NCS variable correctly
! - Bug fix in SOA_LUMP; other minor bug fixes
!
! GEOS-Chem v7-04-09, includes the following modifications:
! - Updated CO for David Streets (2001 for China, 2000 elsewhere)
! - Now reset negative SPHU to a very small positive #
! - Remove use of TINY(1d0) to avoid NaN's on SUN platform
! - Minor bug fixes and deleted obsolete code
!
! Revision 1.38  2006/08/14 17:58:10  bmy
! GEOS-Chem v7-04-08, includes the following modifications:
! - Now add David Streets' emissions for China & SE Asia
! - Removed support for GEOS-1 and GEOS-STRAT met fields
! - Removed support for LINUX_IFC and LINUX_EFC compilers
!
! Revision 1.37  2006/06/28 17:26:52  bmy
! GEOS-Chem v7-04-06, includes the following modifications:
! - Now add BRAVO emissions (NOx, CO, SO2) over N. Mexico
! - Turn off HO2 uptake by aerosols in SMVGEAR mechanism
! - Bug fix: GEOS-4 convection now conserves mixing ratio
! - Other minor bug fixes & improvements
!
! Revision 1.36  2006/06/06 14:26:07  bmy
! GEOS-Chem v7-04-05, includes the following modifications:
! - Now gets ISOP that has reacted w/ OH from SMVGEAR (cf. D. Henze)
! - Incorporated IPCC future emission scale factors (cf. S. Wu)
! - Other minor bug fixes
!
! Revision 1.35  2006/05/26 17:45:24  bmy
! GEOS-Chem v7-04-04, includes the following modifications:
! - Now updated for SOA production from ISOP (cf D. Henze)
! - Now archive SOA concentrations in [ug/m3] ("diag42_mod.f")
! - Other minor bug fixes
!
! Revision 1.34  2006/05/15 17:52:52  bmy
! GEOS-Chem v7-04-03, includes the following modifications:
! - Added near-land formulation for lightning
! - Now can use CTH, MFLUX, PRECON params for lightning
!   (NOTE: new lightning is only applied for GEOS-4 for time being)
! - Added ND56 diagnostic for lightning flash rates
! - Fixed Feb 28 -> Mar 1 transition for GCAP (i.e. no leap years)
! - Other minor bug fixes
!
! Revision 1.33  2006/03/24 20:22:53  bmy
! GEOS-CHEM v7-04-01, includes the following modifications:
! - Updates to new Hg simulation (eck, cdh, sas)
! - Changed Reynold's # criterion for aerodyn smooth surfaces in drydep
! - Standardized several bug fixes implemented in v7-03-06 patch
! - Bug fix: Now call MAKE_RH from "main.f" to avoid problems in drydep
! - Removed obsolete code
!
      PROGRAM GEOS_CHEM
! 
!******************************************************************************
!                                                                       
!                                                                       
!     GGGGGG  EEEEEEE  OOOOO  SSSSSSS       CCCCCC H     H EEEEEEE M     M    
!    G        E       O     O S            C       H     H E       M M M M    
!    G   GGG  EEEEEE  O     O SSSSSSS      C       HHHHHHH EEEEEE  M  M  M    
!    G     G  E       O     O       S      C       H     H E       M     M    
!     GGGGGG  EEEEEEE  OOOOO  SSSSSSS       CCCCCC H     H EEEEEEE M     M    
!                                                                       
!                                                                       
!                 (formerly known as the Harvard-GEOS model)
!           for 4 x 5, 2 x 2.5 global grids and 1 x 1 nested grids
!
!       Contact: Bob Yantosca, Harvard University (bmy@io.as.harvard.edu)
!                                                                     
!******************************************************************************
!
!  See the GEOS-Chem Web Site:
!
!     http://www.as.harvard.edu/chemistry/trop/geos/
!
!  and  the GEOS-CHEM User's Guide:
!
!     http://www.as.harvard.edu/chemistry/trop/geos/doc/man/
!
!  for the most up-to-date GEOS-CHEM documentation on the following topics:
!
!     - installation, compilation, and execution
!     - coding practice and style
!     - input files and met field data files
!     - horizontal and vertical resolution
!     - modification history
!
!******************************************************************************
!
      ! References to F90 modules 
      USE A3_READ_MOD,       ONLY : GET_A3_FIELDS
      USE A3_READ_MOD,       ONLY : OPEN_A3_FIELDS
      USE A3_READ_MOD,       ONLY : UNZIP_A3_FIELDS
      USE A6_READ_MOD,       ONLY : GET_A6_FIELDS
      USE A6_READ_MOD,       ONLY : OPEN_A6_FIELDS
      USE A6_READ_MOD,       ONLY : UNZIP_A6_FIELDS
      USE BENCHMARK_MOD,     ONLY : STDRUN
      USE CHEMISTRY_MOD,     ONLY : DO_CHEMISTRY
      USE CONVECTION_MOD,    ONLY : DO_CONVECTION
      USE COMODE_MOD,        ONLY : INIT_COMODE
      USE DIAG_MOD,          ONLY : DIAGCHLORO
      USE DIAG41_MOD,        ONLY : DIAG41,          ND41
      USE DIAG42_MOD,        ONLY : DIAG42,          ND42
      USE DIAG48_MOD,        ONLY : DIAG48,          ITS_TIME_FOR_DIAG48
      USE DIAG49_MOD,        ONLY : DIAG49,          ITS_TIME_FOR_DIAG49
      USE DIAG50_MOD,        ONLY : DIAG50,          DO_SAVE_DIAG50
      USE DIAG51_MOD,        ONLY : DIAG51,          DO_SAVE_DIAG51
      USE DIAG_OH_MOD,       ONLY : PRINT_DIAG_OH
      USE DAO_MOD,           ONLY : AD,              AIRQNT  
      USE DAO_MOD,           ONLY : AVGPOLE,         CLDTOPS
      USE DAO_MOD,           ONLY : CONVERT_UNITS,   COPY_I6_FIELDS
      USE DAO_MOD,           ONLY : COSSZA,          INIT_DAO
      USE DAO_MOD,           ONLY : INTERP,          PS1
      USE DAO_MOD,           ONLY : PS2,             PSC2          
      USE DAO_MOD,           ONLY : T,               TS            
      USE DAO_MOD,           ONLY : SUNCOS,          SUNCOSB
      USE DAO_MOD,           ONLY : MAKE_RH
      USE DRYDEP_MOD,        ONLY : DO_DRYDEP
      USE EMISSIONS_MOD,     ONLY : DO_EMISSIONS
      USE ERROR_MOD,         ONLY : DEBUG_MSG
      USE FILE_MOD,          ONLY : IU_BPCH,         IU_DEBUG
      USE FILE_MOD,          ONLY : IU_ND48,         IU_SMV2LOG    
      USE FILE_MOD,          ONLY : CLOSE_FILES
      USE GLOBAL_CH4_MOD,    ONLY : INIT_GLOBAL_CH4, CH4_AVGTP
      USE GCAP_READ_MOD,     ONLY : GET_GCAP_FIELDS
      USE GCAP_READ_MOD,     ONLY : OPEN_GCAP_FIELDS
      USE GCAP_READ_MOD,     ONLY : UNZIP_GCAP_FIELDS
      USE GWET_READ_MOD,     ONLY : GET_GWET_FIELDS
      USE GWET_READ_MOD,     ONLY : OPEN_GWET_FIELDS
      USE GWET_READ_MOD,     ONLY : UNZIP_GWET_FIELDS
      USE I6_READ_MOD,       ONLY : GET_I6_FIELDS_1
      USE I6_READ_MOD,       ONLY : GET_I6_FIELDS_2
      USE I6_READ_MOD,       ONLY : OPEN_I6_FIELDS
      USE I6_READ_MOD,       ONLY : UNZIP_I6_FIELDS
      USE INPUT_MOD,         ONLY : READ_INPUT_FILE
      USE LAI_MOD,           ONLY : RDISOLAI
      USE LIGHTNING_NOX_MOD, ONLY : LIGHTNING
      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      !%%% NOTE: Temporary kludge: For GEOS-4 we want to use the new near-land
      !%%% lightning formulation.  But for the time being, we must keep the 
      !%%% existing lightning for other met field types. (ltm, bmy, 5/10/06)
      USE LIGHTNING_NOX_NL_MOD, ONLY : LIGHTNING_NL
      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      USE LOGICAL_MOD,       ONLY : LEMIS,     LCHEM, LUNZIP,  LDUST
      USE LOGICAL_MOD,       ONLY : LLIGHTNOX, LPRT,  LSTDRUN, LSVGLB
      USE LOGICAL_MOD,       ONLY : LWAIT,     LTRAN, LUPBD,   LCONV
      USE LOGICAL_MOD,       ONLY : LWETD,     LTURB, LDRYD,   LMEGAN  
      USE LOGICAL_MOD,       ONLY : LDYNOCEAN, LSOA,  LVARTROP
      USE MEGAN_MOD,         ONLY : INIT_MEGAN
      USE MEGAN_MOD,         ONLY : UPDATE_T_15_AVG
      USE MEGAN_MOD,         ONLY : UPDATE_T_DAY
      USE PBL_MIX_MOD,       ONLY : DO_PBL_MIX
      USE OCEAN_MERCURY_MOD, ONLY : MAKE_OCEAN_Hg_RESTART
      USE OCEAN_MERCURY_MOD, ONLY : READ_OCEAN_Hg_RESTART
      USE PLANEFLIGHT_MOD,   ONLY : PLANEFLIGHT
      USE PLANEFLIGHT_MOD,   ONLY : SETUP_PLANEFLIGHT 
      USE PRESSURE_MOD,      ONLY : INIT_PRESSURE
      USE PRESSURE_MOD,      ONLY : SET_FLOATING_PRESSURE
      USE TIME_MOD,          ONLY : GET_NYMDb,        GET_NHMSb
      USE TIME_MOD,          ONLY : GET_NYMD,         GET_NHMS
      USE TIME_MOD,          ONLY : GET_A3_TIME,      GET_FIRST_A3_TIME
      USE TIME_MOD,          ONLY : GET_A6_TIME,      GET_FIRST_A6_TIME
      USE TIME_MOD,          ONLY : GET_I6_TIME,      GET_MONTH
      USE TIME_MOD,          ONLY : GET_TAU,          GET_TAUb
      USE TIME_MOD,          ONLY : GET_TS_CHEM,      GET_TS_DYN
      USE TIME_MOD,          ONLY : GET_ELAPSED_SEC,  GET_TIME_AHEAD
      USE TIME_MOD,          ONLY : GET_DAY_OF_YEAR,  ITS_A_NEW_DAY
      USE TIME_MOD,          ONLY : ITS_A_NEW_SEASON, GET_SEASON
      USE TIME_MOD,          ONLY : ITS_A_NEW_MONTH,  GET_NDIAGTIME
      USE TIME_MOD,          ONLY : ITS_A_LEAPYEAR,   GET_YEAR
      USE TIME_MOD,          ONLY : ITS_TIME_FOR_A3,  ITS_TIME_FOR_A6
      USE TIME_MOD,          ONLY : ITS_TIME_FOR_I6,  ITS_TIME_FOR_CHEM
      USE TIME_MOD,          ONLY : ITS_TIME_FOR_CONV,ITS_TIME_FOR_DEL
      USE TIME_MOD,          ONLY : ITS_TIME_FOR_DIAG,ITS_TIME_FOR_DYN
      USE TIME_MOD,          ONLY : ITS_TIME_FOR_EMIS,ITS_TIME_FOR_EXIT
      USE TIME_MOD,          ONLY : ITS_TIME_FOR_UNIT,ITS_TIME_FOR_UNZIP
      USE TIME_MOD,          ONLY : SET_CT_CONV,      SET_CT_DYN
      USE TIME_MOD,          ONLY : SET_CT_EMIS,      SET_CT_CHEM
      USE TIME_MOD,          ONLY : SET_DIAGb,        SET_DIAGe
      USE TIME_MOD,          ONLY : SET_CURRENT_TIME, PRINT_CURRENT_TIME
      USE TIME_MOD,          ONLY : SET_ELAPSED_MIN,  SYSTEM_TIMESTAMP
      USE TRACER_MOD,        ONLY : CHECK_STT, N_TRACERS, STT, TCVV
      USE TRACER_MOD,        ONLY : ITS_AN_AEROSOL_SIM
      USE TRACER_MOD,        ONLY : ITS_A_CH4_SIM
      USE TRACER_MOD,        ONLY : ITS_A_FULLCHEM_SIM
      USE TRACER_MOD,        ONLY : ITS_A_MERCURY_SIM
      USE TRANSPORT_MOD,     ONLY : DO_TRANSPORT
      USE TROPOPAUSE_MOD,    ONLY : READ_TROPOPAUSE, CHECK_VAR_TROP
      USE RESTART_MOD,       ONLY : MAKE_RESTART_FILE, READ_RESTART_FILE
      USE UPBDFLX_MOD,       ONLY : DO_UPBDFLX,        UPBDFLX_NOY
      USE UVALBEDO_MOD,      ONLY : READ_UVALBEDO
      USE WETSCAV_MOD,       ONLY : INIT_WETSCAV,      DO_WETDEP
      USE XTRA_READ_MOD,     ONLY : GET_XTRA_FIELDS,   OPEN_XTRA_FIELDS
      USE XTRA_READ_MOD,     ONLY : UNZIP_XTRA_FIELDS

      ! Force all variables to be declared explicitly
      IMPLICIT NONE
      
      ! Header files 
#     include "CMN_SIZE"          ! Size parameters
#     include "CMN_DIAG"          ! Diagnostic switches, NJDAY
#     include "CMN_GCTM"          ! Physical constants

      ! Local variables
      LOGICAL            :: FIRST = .TRUE.
      LOGICAL            :: LXTRA 
      INTEGER            :: I,           IOS,   J,         K,      L
      INTEGER            :: N,           JDAY,  NDIAGTIME, N_DYN
      INTEGER            :: N_DYN_STEPS, NSECb, N_STEP,    DATE(2)
      INTEGER            :: YEAR,        MONTH, DAY,       DAY_OF_YEAR
      INTEGER            :: SEASON,      NYMD,  NYMDb,     NHMS
      INTEGER            :: ELAPSED_SEC, NHMSb
      REAL*8             :: TAU,         TAUb         
      CHARACTER(LEN=255) :: ZTYPE

      !=================================================================
      ! GEOS-CHEM starts here!                                            
      !=================================================================

      ! Display current grid resolution and data set type
      CALL DISPLAY_GRID_AND_MODEL

      !=================================================================
      !            ***** I N I T I A L I Z A T I O N *****
      !=================================================================

      ! Read input file and call init routines from other modules
      CALL READ_INPUT_FILE 
      IF ( LPRT ) CALL DEBUG_MSG( '### MAIN: a READ_INPUT_FILE' )

      ! Initialize met field arrays from "dao_mod.f"
      CALL INIT_DAO
      IF ( LPRT ) CALL DEBUG_MSG( '### MAIN: a INIT_DAO' )

      ! Initialize diagnostic arrays and counters
      CALL INITIALIZE( 2 )
      CALL INITIALIZE( 3 )
      IF ( LPRT ) CALL DEBUG_MSG( '### MAIN: a INITIALIZE' )

      ! Initialize the new hybrid pressure module.  Define Ap and Bp.
      CALL INIT_PRESSURE
      IF ( LPRT ) CALL DEBUG_MSG( '### MAIN: a INIT_PRESSURE' )

      ! Read annual mean tropopause if not a variable tropopause
      ! read_tropopause is obsolete with variable tropopause
      IF ( .not. LVARTROP ) THEN
         CALL READ_TROPOPAUSE
         IF ( LPRT ) CALL DEBUG_MSG( '### MAIN: a READ_TROPOPAUSE' )
      ENDIF

      ! Initialize allocatable SMVGEAR arrays
      IF ( LEMIS .or. LCHEM ) THEN
         IF ( ITS_A_FULLCHEM_SIM() ) CALL INIT_COMODE
         IF ( ITS_AN_AEROSOL_SIM() ) CALL INIT_COMODE
         IF ( LPRT ) CALL DEBUG_MSG( '### MAIN: a INIT_COMODE' )
      ENDIF
         
      ! Allocate arrays from "global_ch4_mod.f" for CH4 run 
      IF ( ITS_A_CH4_SIM() ) CALL INIT_GLOBAL_CH4

      ! Initialize MEGAN arrays, get 15-day avg temperatures
      IF ( LMEGAN ) THEN
         CALL INIT_MEGAN
         CALL INITIALIZE( 2 )
         CALL INITIALIZE( 3 )
         IF ( LPRT ) CALL DEBUG_MSG( '### MAIN: a INIT_MEGAN' )
      ENDIF

      ! Local flag for reading XTRA fields for GEOS-3
      !LXTRA = ( LDUST .or. LMEGAN )
      LXTRA = LMEGAN

      ! Define time variables for use below
      NHMS  = GET_NHMS()
      NHMSb = GET_NHMSb()
      NYMD  = GET_NYMD()
      NYMDb = GET_NYMDb()
      TAU   = GET_TAU()
      TAUb  = GET_TAUb()

      !=================================================================
      !   ***** U N Z I P   M E T   F I E L D S  @ start of run *****
      !=================================================================
      IF ( LUNZIP ) THEN

         !---------------------
         ! Remove all files
         !---------------------

         ! Type of unzip operation
         ZTYPE = 'remove all'
         
         ! Remove any leftover A-3, A-6, I-6, in temp dir
         CALL UNZIP_A3_FIELDS( ZTYPE )
         CALL UNZIP_A6_FIELDS( ZTYPE )
         CALL UNZIP_I6_FIELDS( ZTYPE )

#if   defined( GEOS_3 )
         ! Remove GEOS-3 GWET and XTRA files 
         IF ( LDUST ) CALL UNZIP_GWET_FIELDS( ZTYPE )
         IF ( LXTRA ) CALL UNZIP_XTRA_FIELDS( ZTYPE )
#endif

#if   defined( GCAP )
         ! Unzip GCAP PHIS field (if necessary)
         CALL UNZIP_GCAP_FIELDS( ZTYPE )
#endif

         !---------------------
         ! Unzip in foreground
         !---------------------

         ! Type of unzip operation
         ZTYPE = 'unzip foreground'

         ! Unzip A-3, A-6, I-6 files for START of run
         CALL UNZIP_A3_FIELDS( ZTYPE, NYMDb )
         CALL UNZIP_A6_FIELDS( ZTYPE, NYMDb )
         CALL UNZIP_I6_FIELDS( ZTYPE, NYMDb )

#if   defined( GEOS_3 )
         ! Unzip GEOS-3 GWET and XTRA fields for START of run
         IF ( LDUST ) CALL UNZIP_GWET_FIELDS( ZTYPE, NYMDb )
         IF ( LXTRA ) CALL UNZIP_XTRA_FIELDS( ZTYPE, NYMDb )
#endif

#if   defined( GCAP )
         ! Unzip GCAP PHIS field (if necessary)
         CALL UNZIP_GCAP_FIELDS( ZTYPE )
#endif

         !### Debug output
         IF ( LPRT ) CALL DEBUG_MSG( '### MAIN: a UNZIP' )
      ENDIF
      
      !=================================================================
      !    ***** R E A D   M E T   F I E L D S  @ start of run *****
      !=================================================================

      ! Open and read A-3 fields
      DATE = GET_FIRST_A3_TIME()
      CALL OPEN_A3_FIELDS( DATE(1), DATE(2) )
      CALL GET_A3_FIELDS(  DATE(1), DATE(2) )
      IF ( LPRT ) CALL DEBUG_MSG( '### MAIN: a 1st A3 TIME' )

      ! For MEGAN biogenics, update hourly temps w/in 15-day window
      IF ( LMEGAN ) THEN
         CALL UPDATE_T_DAY
         IF ( LPRT ) CALL DEBUG_MSG( '### MAIN: UPDATE T_DAY' )
      ENDIF

      ! Open & read A-6 fields
      DATE = GET_FIRST_A6_TIME()
      CALL OPEN_A6_FIELDS( DATE(1), DATE(2) ) 
      CALL GET_A6_FIELDS(  DATE(1), DATE(2) )      
      IF ( LPRT ) CALL DEBUG_MSG( '### MAIN: a 1st A6 TIME' )

      ! Open & read I-6 fields
      DATE = (/ NYMD, NHMS /)
      CALL OPEN_I6_FIELDS(  DATE(1), DATE(2) )
      CALL GET_I6_FIELDS_1( DATE(1), DATE(2) )
      IF ( LPRT ) CALL DEBUG_MSG( '### MAIN: a 1st I6 TIME' )
      
#if   defined( GEOS_3 )
      ! Open & read GEOS-3 GWET fields
      IF ( LDUST ) THEN
         DATE = GET_FIRST_A3_TIME()
         CALL OPEN_GWET_FIELDS( DATE(1), DATE(2) )
         CALL GET_GWET_FIELDS(  DATE(1), DATE(2) ) 
         IF ( LPRT ) CALL DEBUG_MSG( '### MAIN: a 1st GWET TIME' )
      ENDIF

      ! Open & read GEOS-3 XTRA fields
      IF ( LXTRA ) THEN
         DATE = GET_FIRST_A3_TIME()
         CALL OPEN_XTRA_FIELDS( DATE(1), DATE(2) )
         CALL GET_XTRA_FIELDS(  DATE(1), DATE(2) ) 
         IF ( LPRT ) CALL DEBUG_MSG( '### MAIN: a 1st XTRA TIME' )
      ENDIF
#endif

#if   defined( GCAP )
      ! Read GCAP PHIS and LWI fields (if necessary)
      CALL OPEN_GCAP_FIELDS
      CALL GET_GCAP_FIELDS

      ! Remove temporary file (if necessary)
      IF ( LUNZIP ) THEN
         CALL UNZIP_GCAP_FIELDS( 'remove date' )
      ENDIF
#endif

      ! Compute avg surface pressure near polar caps
      CALL AVGPOLE( PS1 )
      IF ( LPRT ) CALL DEBUG_MSG( '### MAIN: a AVGPOLE' )

      ! Call AIRQNT to compute air mass quantities from PS1
      CALL SET_FLOATING_PRESSURE( PS1 )
      CALL AIRQNT
      IF ( LPRT ) CALL DEBUG_MSG( '### MAIN: a AIRQNT' )

      ! Compute lightning NOx emissions [molec/box/6h]
      IF ( LLIGHTNOX ) THEN
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%% NOTE: Temporary kludge: For GEOS-4 we want to use the new near-land 
!%%% lightning formulation.  But for the time being, we must keep the existing
!%%% lightning for other met field types. (ltm, bmy, 5/10/06)
#if   defined( GEOS_4 )
         CALL LIGHTNING_NL
#else
         CALL LIGHTNING( T, CLDTOPS )
#endif
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

         IF ( LPRT ) CALL DEBUG_MSG( '### MAIN: a LIGHTNING' )
      ENDIF

      ! Read land types and fractions from "vegtype.global"
      CALL RDLAND   
      IF ( LPRT ) CALL DEBUG_MSG( '### MAIN: a RDLAND' )

      ! Initialize PBL quantities but do not do mixing
      CALL DO_PBL_MIX( .FALSE. )
      IF ( LPRT ) CALL DEBUG_MSG( '### MAIN: a TURBDAY:1' )

      !=================================================================
      !       *****  I N I T I A L   C O N D I T I O N S *****
      !=================================================================

      ! Read initial tracer conditions
      CALL READ_RESTART_FILE( NYMDb, NHMSb )
      IF ( LPRT ) CALL DEBUG_MSG( '### MAIN: a READ_RESTART_FILE' )

      ! Read ocean Hg initial conditions (if necessary)
      IF ( ITS_A_MERCURY_SIM() .and. LDYNOCEAN ) THEN
         CALL READ_OCEAN_Hg_RESTART( NYMDb, NHMSb )
         IF ( LPRT ) CALL DEBUG_MSG( '### MAIN: a READ_OCEAN_RESTART' )
      ENDIF

      ! Save initial tracer masses to disk for benchmark runs
      IF ( LSTDRUN ) CALL STDRUN( LBEGIN=.TRUE. )

      !=================================================================
      !      ***** 6 - H O U R   T I M E S T E P   L O O P  *****
      !=================================================================      

      ! Echo message before first timestep
      WRITE( 6, '(a)' )
      WRITE( 6, '(a)' ) REPEAT( '*', 44 )
      WRITE( 6, '(a)' ) '* B e g i n   T i m e   S t e p p i n g !! *'
      WRITE( 6, '(a)' ) REPEAT( '*', 44 )
      WRITE( 6, '(a)' ) 

      ! NSTEP is the number of dynamic timesteps w/in a 6-h interval
      N_DYN_STEPS = 360 / GET_TS_DYN()

      ! Start a new 6-h loop
      DO 

      ! Compute time parameters at start of 6-h loop
      CALL SET_CURRENT_TIME

      ! NSECb is # of seconds at the start of 6-h loop
      NSECb = GET_ELAPSED_SEC()

      ! Get dynamic timestep in seconds
      N_DYN = 60d0 * GET_TS_DYN()

      !=================================================================
      !     ***** D Y N A M I C   T I M E S T E P   L O O P *****
      !=================================================================
      DO N_STEP = 1, N_DYN_STEPS
    
         ! Compute & print time quantities at start of dyn step
         CALL SET_CURRENT_TIME
         CALL PRINT_CURRENT_TIME

         ! Set time variables for dynamic loop
         DAY_OF_YEAR = GET_DAY_OF_YEAR()
         ELAPSED_SEC = GET_ELAPSED_SEC()
         MONTH       = GET_MONTH()
         NHMS        = GET_NHMS()
         NYMD        = GET_NYMD()
         TAU         = GET_TAU()
         YEAR        = GET_YEAR()
         SEASON      = GET_SEASON()

         !### Debug
         IF ( LPRT ) CALL DEBUG_MSG( '### MAIN: a SET_CURRENT_TIME' )

         !==============================================================
         !   ***** W R I T E   D I A G N O S T I C   F I L E S *****
         !==============================================================
         IF ( ITS_TIME_FOR_BPCH() ) THEN
            
            ! Set time at end of diagnostic timestep
            CALL SET_DIAGe( TAU )

            ! Write bpch file
            CALL DIAG3  

            ! Flush file units
            CALL CTM_FLUSH

            !===========================================================
            !      ***** W R I T E   R E S T A R T   F I L E *****
            !===========================================================
            IF ( LSVGLB ) THEN

               ! Make atmospheric restart file
               CALL MAKE_RESTART_FILE( NYMD, NHMS, TAU )
                  
               ! Make ocean mercury restart file
               IF ( ITS_A_MERCURY_SIM() .and. LDYNOCEAN ) THEN
                  CALL MAKE_OCEAN_Hg_RESTART( NYMD, NHMS, TAU )
               ENDIF

               !### Debug
               IF ( LPRT ) THEN
                  CALL DEBUG_MSG( '### MAIN: a MAKE_RESTART_FILE' )
               ENDIF
            ENDIF

            ! Set time at beginning of next diagnostic timestep
            CALL SET_DIAGb( TAU )

            !===========================================================
            !        ***** Z E R O   D I A G N O S T I C S *****
            !===========================================================
            CALL INITIALIZE( 2 ) ! Zero arrays
            CALL INITIALIZE( 3 ) ! Zero counters
         ENDIF

         !==============================================================
         !       ***** T E S T   F O R   E N D   O F   R U N *****
         !==============================================================
         IF ( ITS_TIME_FOR_EXIT() ) GOTO 9999

         !===============================================================
         !        ***** U N Z I P   M E T   F I E L D S *****
         !===============================================================
         IF ( LUNZIP .and. ITS_TIME_FOR_UNZIP() ) THEN
            
            ! Get the date & time for 12h (720 mins) from now
            DATE = GET_TIME_AHEAD( 720 )

            ! If LWAIT=T then wait for the met fields to be
            ! fully unzipped before proceeding w/ the run.
            ! Otherwise, unzip fields in the background
            IF ( LWAIT ) THEN
               ZTYPE = 'unzip foreground'
            ELSE
               ZTYPE = 'unzip background'
            ENDIF
            
            ! Unzip A3, A6, I6 fields
            CALL UNZIP_A3_FIELDS( ZTYPE, DATE(1) )
            CALL UNZIP_A6_FIELDS( ZTYPE, DATE(1) )
            CALL UNZIP_I6_FIELDS( ZTYPE, DATE(1) )

#if   defined( GEOS_3 )
            ! Unzip GEOS-3 GWET & XTRA fields
            IF ( LDUST ) CALL UNZIP_GWET_FIELDS( ZTYPE, DATE(1) )
            IF ( LXTRA ) CALL UNZIP_XTRA_FIELDS( ZTYPE, DATE(1) )
#endif
         ENDIF     

         !===============================================================
         !        ***** R E M O V E   M E T   F I E L D S *****  
         !===============================================================
         IF ( LUNZIP .and. ITS_TIME_FOR_DEL() ) THEN

            ! Type of operation
            ZTYPE = 'remove date'

            ! Remove A-3, A-6, and I-6 files only for the current date
            CALL UNZIP_A3_FIELDS( ZTYPE, NYMD )
            CALL UNZIP_A6_FIELDS( ZTYPE, NYMD )
            CALL UNZIP_I6_FIELDS( ZTYPE, NYMD )

#if   defined( GEOS_3 )
            ! Remove GEOS-3 GWET & XTRA fields only for the current date
            IF ( LDUST ) CALL UNZIP_GWET_FIELDS( ZTYPE, NYMD )
            IF ( LXTRA ) CALL UNZIP_XTRA_FIELDS( ZTYPE, NYMD )
#endif
         ENDIF   

         !==============================================================
         !          ***** R E A D   A - 3   F I E L D S *****
         !==============================================================
         IF ( ITS_TIME_FOR_A3() ) THEN

            ! Get the date/time for the next A-3 data block
            DATE = GET_A3_TIME()

            ! Open & read A-3 fields
            CALL OPEN_A3_FIELDS( DATE(1), DATE(2) )
            CALL GET_A3_FIELDS(  DATE(1), DATE(2) )

            ! Update daily mean temperature archive for MEGAN biogenics
            IF ( LMEGAN ) CALL UPDATE_T_DAY 

#if   defined( GEOS_3 )
            ! Read GEOS-3 GWET fields
            IF ( LDUST ) THEN
               CALL OPEN_GWET_FIELDS( DATE(1), DATE(2) )
               CALL GET_GWET_FIELDS(  DATE(1), DATE(2) )           
            ENDIF
            
            ! Read GEOS-3 PARDF, PARDR, SNOW fields
            IF ( LXTRA ) THEN
               CALL OPEN_XTRA_FIELDS( DATE(1), DATE(2) )
               CALL GET_XTRA_FIELDS(  DATE(1), DATE(2) )           
            ENDIF
#endif
         ENDIF

         !==============================================================
         !          ***** R E A D   A - 6   F I E L D S *****   
         !==============================================================
         IF ( ITS_TIME_FOR_A6() ) THEN
            
            ! Get the date/time for the next A-6 data block
            DATE = GET_A6_TIME()

            ! Open and read A-6 fields
            CALL OPEN_A6_FIELDS( DATE(1), DATE(2) )
            CALL GET_A6_FIELDS(  DATE(1), DATE(2) )

            ! Since CLDTOPS is an A-6 field, update the
            ! lightning NOx emissions [molec/box/6h]
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%% NOTE: Temporary kludge: For GEOS-4 we want to use the new near-land 
!%%% lightning formulation.  But for the time being, we must keep the 
!%%% existing lightning for other met field types. (ltm, bmy, 5/10/06)
            IF ( LLIGHTNOX ) THEN
#if   defined( GEOS_4 )
               CALL LIGHTNING_NL
#else 
               CALL LIGHTNING( T, CLDTOPS )
#endif
            ENDIF
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         ENDIF

         !==============================================================
         !          ***** R E A D   I - 6   F I E L D S *****   
         !==============================================================
         IF ( ITS_TIME_FOR_I6() ) THEN

            ! Get the date/time for the next I-6 data block
            DATE = GET_I6_TIME()

            ! Open and read files
            CALL OPEN_I6_FIELDS(  DATE(1), DATE(2) )
            CALL GET_I6_FIELDS_2( DATE(1), DATE(2) )

            ! Compute avg pressure at polar caps 
            CALL AVGPOLE( PS2 )
         ENDIF

         !==============================================================
         ! ***** M O N T H L Y   O R   S E A S O N A L   D A T A *****
         !==============================================================

         ! UV albedoes
         IF ( LCHEM .and. ITS_A_NEW_MONTH() ) THEN
            CALL READ_UVALBEDO( MONTH )
         ENDIF

         ! Fossil fuel emissions (SMVGEAR)
         IF ( ITS_A_FULLCHEM_SIM() ) THEN
            IF ( LEMIS .and. ITS_A_NEW_SEASON() ) THEN
               CALL ANTHROEMS( SEASON )
            ENDIF
         ENDIF

         !==============================================================
         !              ***** D A I L Y   D A T A *****
         !==============================================================
         IF ( ITS_A_NEW_DAY() ) THEN 

            ! Read leaf-area index (needed for drydep)
            CALL RDLAI( DAY_OF_YEAR, MONTH )

            ! For MEGAN biogenics ...
            IF ( LMEGAN ) THEN

               ! Read AVHRR daily leaf-area-index
               CALL RDISOLAI( DAY_OF_YEAR, MONTH )

               ! Compute 15-day average temperature for MEGAN
               CALL UPDATE_T_15_AVG
            ENDIF
               
            ! Also read soil-type info for fullchem simulation
            IF ( ITS_A_FULLCHEM_SIM() ) CALL RDSOIL 

            !### Debug
            IF ( LPRT ) CALL DEBUG_MSG ( '### MAIN: a DAILY DATA' )
         ENDIF

         ! Get averaging intervals for local-time diagnostics
         ! (NOTE: maybe improve this later on)
         CALL DIAG_2PM
     
         !==============================================================
         !   ***** I N T E R P O L A T E   Q U A N T I T I E S *****   
         !==============================================================
         
         ! Interpolate I-6 fields to current dynamic timestep, 
         ! based on their values at NSEC and NSEC+N_DYN
         CALL INTERP( NSECb, ELAPSED_SEC, N_DYN )

         ! Case of variable tropopause:
         ! Check LLTROP and set LMIN, LMAX, and LPAUSE
         ! since this is not done with READ_TROPOPAUSE anymore.
         ! (Need to double-check that LMIN, Lmax are not used before-phs) 
         IF ( LVARTROP ) CALL CHECK_VAR_TROP
         
         ! If we are not doing transport, then make sure that
         ! the floating pressure is set to PSC2 (bdf, bmy, 8/22/02)
         IF ( .not. LTRAN ) CALL SET_FLOATING_PRESSURE( PSC2 )

         ! Compute airmass quantities at each grid box 
         CALL AIRQNT
         
         ! Compute the cosine of the solar zenith angle at each grid box
         CALL COSSZA( DAY_OF_YEAR, NHMSb, ELAPSED_SEC, SUNCOS )
         
         ! For SMVGEAR II, we also need to compute SUNCOS at
         ! the end of this chemistry timestep (bdf, bmy, 4/1/03)
         IF ( LCHEM .and. ITS_A_FULLCHEM_SIM() ) THEN
            CALL COSSZA( DAY_OF_YEAR,                  NHMSb, 
     &                   ELAPSED_SEC+GET_TS_CHEM()*60, SUNCOSB )
         ENDIF

         ! Compute tropopause height for ND55 diagnostic
         IF ( ND55 > 0 ) CALL TROPOPAUSE

#if   defined( GEOS_3 )

         ! 1998 GEOS-3 carries the ground temperature and not the air
         ! temperature -- thus TS will be 2-3 K too high.  As a quick fix, 
         ! copy the temperature at the first sigma level into TS. 
         ! (mje, bnd, bmy, 7/3/01)
         IF ( YEAR == 1998 ) TS(:,:) = T(:,:,1)
#endif

         ! Update dynamic timestep
         CALL SET_CT_DYN( INCREMENT=.TRUE. )

         !### Debug
         IF ( LPRT ) CALL DEBUG_MSG( '### MAIN: a INTERP, etc' )

         !==============================================================
         !   ***** U N I T   C O N V E R S I O N  ( kg -> v/v ) *****
         !==============================================================
         IF ( ITS_TIME_FOR_UNIT() ) THEN
            CALL CONVERT_UNITS( 1,  N_TRACERS, TCVV, AD, STT )

            !### Debug
            IF ( LPRT ) CALL DEBUG_MSG( '### MAIN: a CONVERT_UNITS:1' )
         ENDIF

         !==============================================================
         !     ***** S T R A T O S P H E R I C   F L U X E S *****
         !==============================================================
         IF ( LUPBD ) CALL DO_UPBDFLX

         !==============================================================
         !              ***** T R A N S P O R T *****
         !==============================================================
         IF ( ITS_TIME_FOR_DYN() ) THEN

            ! Call the appropritate version of TPCORE
            IF ( LTRAN ) CALL DO_TRANSPORT               

            ! Reset air mass quantities
            CALL AIRQNT

            ! Repartition [NOy] species after transport
            IF ( LUPBD .and. ITS_A_FULLCHEM_SIM() ) THEN
               CALL UPBDFLX_NOY( 2 )
            ENDIF

            ! Get relative humidity 
            ! (after recomputing pressure quantities)
            CALL MAKE_RH

            ! Initialize wet scavenging and wetdep fields after
            ! the airmass quantities are reset after transport
            IF ( LCONV .or. LWETD ) CALL INIT_WETSCAV
         ENDIF

         !-------------------------------
         ! Test for convection timestep
         !-------------------------------
         IF ( ITS_TIME_FOR_CONV() ) THEN

            ! Increment the convection timestep
            CALL SET_CT_CONV( INCREMENT=.TRUE. )

            !===========================================================
            !      ***** M I X E D   L A Y E R   M I X I N G *****
            !===========================================================
            CALL DO_PBL_MIX( LTURB )

            !### Debug
            IF ( LPRT ) CALL DEBUG_MSG( '### MAIN: a TURBDAY:2' )

            !===========================================================
            !        ***** C L O U D   C O N V E C T I O N *****
            !===========================================================
            IF ( LCONV ) THEN
               CALL DO_CONVECTION

               !### Debug
               IF ( LPRT ) CALL DEBUG_MSG( '### MAIN: a CONVECTION' )
            ENDIF 
         ENDIF                  

         !==============================================================
         !    ***** U N I T   C O N V E R S I O N  ( v/v -> kg ) *****
         !==============================================================
         IF ( ITS_TIME_FOR_UNIT() ) THEN 
            CALL CONVERT_UNITS( 2, N_TRACERS, TCVV, AD, STT )

            !### Debug
            IF ( LPRT ) CALL DEBUG_MSG( '### MAIN: a CONVERT_UNITS:2' )
         ENDIF

         !==============================================================
         !       ***** A R C H I V E   D I A G N O S T I C S *****
         !==============================================================
         IF ( ITS_TIME_FOR_DYN() ) THEN

            ! Accumulate several diagnostic quantities
            CALL DIAG1

            ! ND41: save PBL height in 1200-1600 LT (amf)
            ! (for comparison w/ Holzworth, 1967)
            IF ( ND41 > 0 ) CALL DIAG41

            ! ND42: SOA concentrations [ug/m3]
            IF ( ND42 > 0 ) CALL DIAG42

            !### Debug
            IF ( LPRT ) CALL DEBUG_MSG( '### MAIN: a DIAGNOSTICS' )
         ENDIF

         !-------------------------------
         ! Test for emission timestep
         !-------------------------------
         IF ( ITS_TIME_FOR_EMIS() ) THEN
               
            ! Increment emission counter
            CALL SET_CT_EMIS( INCREMENT=.TRUE. )

            !========================================================
            !         ***** D R Y   D E P O S I T I O N *****
            !========================================================
            IF ( LDRYD ) CALL DO_DRYDEP

            !========================================================
            !             ***** E M I S S I O N S *****
            !========================================================
            IF ( LEMIS ) CALL DO_EMISSIONS
         ENDIF    

         !===========================================================
         !               ***** C H E M I S T R Y *****
         !===========================================================    

         ! Also need to compute avg P, T for CH4 chemistry (bmy, 1/16/01)
         IF ( ITS_A_CH4_SIM() ) CALL CH4_AVGTP

         ! Every chemistry timestep...
         IF ( ITS_TIME_FOR_CHEM() ) THEN

            ! Increment chemistry timestep counter
            CALL SET_CT_CHEM( INCREMENT=.TRUE. )

            ! Call the appropriate chemistry routine
            CALL DO_CHEMISTRY
         ENDIF 
          
         !==============================================================
         ! ***** W E T   D E P O S I T I O N  (rainout + washout) *****
         !==============================================================
         IF ( LWETD .and. ITS_TIME_FOR_DYN() ) CALL DO_WETDEP

         ! Activate this here someday (bmy, 7/20/04)
         !!==============================================================
         !!       ***** A R C H I V E   D I A G N O S T I C S *****
         !!==============================================================
         !IF ( ITS_TIME_FOR_DYN() ) THEN
         !
         !   ! Accumulate several diagnostic quantities
         !   CALL DIAG1
         !
         !   ! ND41: save PBL height in 1200-1600 LT (amf)
         !   ! (for comparison w/ Holzworth, 1967
         !   IF ( ND41 > 0 ) CALL DIAG41
         !
         !   !### Debug
         !   IF ( LPRT ) CALL DEBUG_MSG( '### MAIN: a DIAGNOSTICS' )
         !ENDIF

         !==============================================================
         !   ***** T I M E S E R I E S   D I A G N O S T I C S  *****
         !
         ! NOTE: Since we are saving soluble tracers, we must move
         !       the ND40, ND49, and ND52 timeseries diagnostics
         !       to after the call to DO_WETDEP (bmy, 4/22/04)
         !============================================================== 

         ! Plane following diagnostic
         IF ( ND40 > 0 ) THEN 
         
            ! Call SETUP_PLANEFLIGHT routine if necessary
            IF ( ITS_A_NEW_DAY() ) THEN
               
               ! If it's a full-chemistry simulation but LCHEM=F,
               ! or if it's an offline simulation, call setup routine 
               IF ( ITS_A_FULLCHEM_SIM() ) THEN
                  IF ( .not. LCHEM ) CALL SETUP_PLANEFLIGHT
               ELSE
                  CALL SETUP_PLANEFLIGHT
               ENDIF
            ENDIF

            ! Archive data along the flight track
            CALL PLANEFLIGHT
         ENDIF
            
         ! Station timeseries
         IF ( ITS_TIME_FOR_DIAG48() ) CALL DIAG48

         ! 3-D timeseries
         IF ( ITS_TIME_FOR_DIAG49() ) CALL DIAG49

         ! 24-hr timeseries
         IF ( DO_SAVE_DIAG50 ) CALL DIAG50

         ! Morning or afternoon timeseries
         IF ( DO_SAVE_DIAG51 ) CALL DIAG51 

         ! Comment out for now 
         !! Column timeseries
         !IF ( ND52 > 0 .and. ITS_TIME_FOR_ND52() ) THEN
         !   CALL DIAG52
         !   IF ( LPRT ) CALL DEBUG_MSG( '### MAIN: a ND52' )
         !ENDIF

         !### After diagnostics
         IF ( LPRT ) CALL DEBUG_MSG( '### MAIN: after TIMESERIES' )

         !==============================================================
         !  ***** E N D   O F   D Y N A M I C   T I M E S T E P *****
         !==============================================================

         ! Check for NaN, Negatives, Infinities in STT once per hour
         IF ( ITS_TIME_FOR_DIAG() ) THEN
            CALL CHECK_STT( 'End of Dynamic Loop' )
         ENDIF

         ! Increment elapsed time
         CALL SET_ELAPSED_MIN
         IF ( LPRT ) CALL DEBUG_MSG( '### MAIN: after SET_ELAPSED_MIN' )
          
      ENDDO

      !=================================================================
      !            ***** C O P Y   I - 6   F I E L D S *****
      !
      !        The I-6 fields at the end of this timestep become
      !        the fields at the beginning of the next timestep
      !=================================================================
      CALL COPY_I6_FIELDS
      IF ( LPRT ) CALL DEBUG_MSG( '### MAIN: after COPY_I6_FIELDS' )

      ENDDO

      !=================================================================
      !         ***** C L E A N U P   A N D   Q U I T *****
      !=================================================================
 9999 CONTINUE

      ! Remove all files from temporary directory 
      IF ( LUNZIP ) THEN
         
         ! Type of operation
         ZTYPE = 'remove all'

         ! Remove A3, A6, I6 fields
         CALL UNZIP_A3_FIELDS( ZTYPE )
         CALL UNZIP_A6_FIELDS( ZTYPE )
         CALL UNZIP_I6_FIELDS( ZTYPE )

#if   defined( GEOS_3 )
         ! Remove GEOS-3 GWET & XTRA fields
         IF ( LDUST ) CALL UNZIP_GWET_FIELDS( ZTYPE )
         IF ( LXTRA ) CALL UNZIP_XTRA_FIELDS( ZTYPE )
#endif

#if   defined( GCAP )
         ! Remove GCAP PHIS field (if necessary)
         CALL UNZIP_GCAP_FIELDS( ZTYPE )
#endif

      ENDIF

      ! Print the mass-weighted mean OH concentration (if applicable)
      CALL PRINT_DIAG_OH

      ! For model benchmarking, save final masses of 
      ! Rn,Pb,Be or Ox to a binary punch file 
      IF ( LSTDRUN ) CALL STDRUN( LBEGIN=.FALSE. )

      ! Close all files
      CALL CLOSE_FILES
      IF ( LPRT ) CALL DEBUG_MSG( '### MAIN: a CLOSE_FILES' )

      ! Deallocate dynamic module arrays
      CALL CLEANUP
      IF ( LPRT ) CALL DEBUG_MSG( '### MAIN: a CLEANUP' )

      ! Print ending time of simulation
      CALL DISPLAY_END_TIME
!
!******************************************************************************
!  Internal procedures -- Use the F90 CONTAINS command to inline 
!  subroutines that only can be called from this main program. 
!
!  All variables referenced in the main program (local variables, F90 
!  module variables, or common block variables) also have scope within 
!  internal subroutines. 
!
!  List of Internal Procedures:
!  ============================================================================
!  (1 ) DISPLAY_GRID_AND_MODEL : Displays resolution, data set, & start time
!  (2 ) GET_NYMD_PHIS          : Gets YYYYMMDD for the PHIS data field
!  (3 ) DISPLAY_SIGMA_LAT_LON  : Displays sigma, lat, and lon information
!  (4 ) GET_WIND10M            : Wrapper for MAKE_WIND10M (from "dao_mod.f")
!  (5 ) CTM_FLUSH              : Flushes diagnostic files to disk
!  (6 ) DISPLAY_END_TIME       : Displays ending time of simulation
!  (7 ) MET_FIELD_DEBUG        : Prints min and max of met fields for debug
!******************************************************************************
!
      CONTAINS

!-----------------------------------------------------------------------------

      SUBROUTINE DISPLAY_GRID_AND_MODEL

      !=================================================================
      ! Internal Subroutine DISPLAY_GRID_AND_MODEL displays the 
      ! appropriate messages for the given model grid and machine type.
      ! It also prints the starting time and date (local time) of the
      ! GEOS-CHEM simulation. (bmy, 12/2/03, 10/18/05)
      !=================================================================

      ! For system time stamp
      CHARACTER(LEN=16) :: STAMP

      !-----------------------
      ! Print resolution info
      !-----------------------
#if   defined( GRID4x5   )
      WRITE( 6, '(a)' )                   
     &    REPEAT( '*', 13 )                                      //
     &    '   S T A R T I N G   4 x 5   G E O S--C H E M   '     //
     &    REPEAT( '*', 13 )

#elif defined( GRID2x25  )
      WRITE( 6, '(a)' ) 
     &    REPEAT( '*', 13 )                                      // 
     &    '   S T A R T I N G   2 x 2.5   G E O S--C H E M   '   //
     &    REPEAT( '*', 13 )

#elif defined( GRID1x125 )
      WRITE( 6, '(a)' ) 
     &    REPEAT( '*', 13 )                                      // 
     &    '   S T A R T I N G   1 x 1.25   G E O S--C H E M   '  //
     &    REPEAT( '*', 13 )

#elif defined( GRID1x1 )
      WRITE( 6, '(a)' ) 
     &    REPEAT( '*', 13 )                                      // 
     &    '   S T A R T I N G   1 x 1   G E O S -- C H E M   '     //
     &    REPEAT( '*', 13 )

#endif

      !-----------------------
      ! Print machine info
      !-----------------------

      ! Get the proper FORMAT statement for the model being used
#if   defined( COMPAQ    )
      WRITE( 6, '(a)' ) 'Created w/ HP/COMPAQ Alpha compiler'
#elif defined( IBM_AIX   )
      WRITE( 6, '(a)' ) 'Created w/ IBM-AIX compiler'
#elif defined( LINUX_PGI )
      WRITE( 6, '(a)' ) 'Created w/ LINUX/PGI compiler'
#elif defined( LINUX_IFORT )
      WRITE( 6, '(a)' ) 'Created w/ LINUX/IFORT (64-bit) compiler'
#elif defined( SGI_MIPS  )
      WRITE( 6, '(a)' ) 'Created w/ SGI MIPSpro compiler'
#elif defined( SPARC     )
      WRITE( 6, '(a)' ) 'Created w/ Sun/SPARC compiler'
#endif

      !-----------------------
      ! Print met field info
      !-----------------------
#if   defined( GEOS_3     )
      WRITE( 6, '(a)' ) 'Using GEOS-3 met fields'
#elif defined( GEOS_4     )
      WRITE( 6, '(a)' ) 'Using GEOS-4/fvDAS met fields'
#elif defined( GEOS_5     )
      WRITE( 6, '(a)' ) 'Using GEOS-5/fvDAS met fields'
#elif defined( GCAP       )
      WRITE( 6, '(a)' ) 'Using GCAP/GISS met fields'
#endif

      !-----------------------
      ! System time stamp
      !-----------------------
      STAMP = SYSTEM_TIMESTAMP()
      WRITE( 6, 100 ) STAMP
 100  FORMAT( /, '===> SIMULATION START TIME: ', a, ' <===', / )

      ! Return to MAIN program
      END SUBROUTINE DISPLAY_GRID_AND_MODEL

!-----------------------------------------------------------------------------

      FUNCTION ITS_TIME_FOR_BPCH() RESULT( DO_BPCH )

      !=================================================================
      ! Internal function ITS_TIME_FOR_BPCH returns TRUE if it is time
      ! to write to the binary punch file and FALSE otherwise.
      !=================================================================

      ! Local variables
      INTEGER :: TODAY, THIS_NJDAY, NHMS, NDIAGTIME
      
      ! Function value
      LOGICAL :: DO_BPCH

      !=================================================================
      ! ITS_TIME_FOR_BPCH begins here!
      !================================================================= 
      
      ! Return FALSE if it's the first timestep
      IF ( GET_TAU() == GET_TAUb() ) THEN
         DO_BPCH = .FALSE.
         RETURN
      ENDIF

      ! Current day of year
      TODAY     = GET_DAY_OF_YEAR()

      ! Current time of day
      NHMS      = GET_NHMS()

      ! Time of day to write bpch files to disk
      NDIAGTIME = GET_NDIAGTIME()

      ! Look up appropriate value of NJDAY array.  We may need to add a
      ! day to skip past the Feb 29 element of NJDAY for non-leap-years.
      IF ( .not. ITS_A_LEAPYEAR( FORCE=.TRUE. ) .and. TODAY > 59 ) THEN
         THIS_NJDAY = NJDAY( TODAY + 1 ) 
      ELSE
         THIS_NJDAY = NJDAY( TODAY )
      ENDIF

      ! Test if this is the day & time to write to the BPCH file!
      IF ( ( THIS_NJDAY > 0 ) .and. NHMS == NDIAGTIME ) THEN
         DO_BPCH = .TRUE.
      ELSE
         DO_BPCH = .FALSE.
      ENDIF

      ! Return to calling program
      END FUNCTION ITS_TIME_FOR_BPCH

!-----------------------------------------------------------------------------

      SUBROUTINE CTM_FLUSH

      !================================================================
      ! Internal subroutine CTM_FLUSH flushes certain diagnostic
      ! file buffers to disk. (bmy, 8/31/00, 7/1/02)
      !
      ! CTM_FLUSH should normally be called after each diagnostic 
      ! output, so that in case the run dies, the output files from 
      ! the last diagnostic timestep will not be lost.  
      !
      ! FLUSH is an intrinsic FORTRAN subroutine and takes as input 
      ! the unit number of the file to be flushed to disk.
      !================================================================
      CALL FLUSH( IU_ND48    )  
      CALL FLUSH( IU_BPCH    )  
      CALL FLUSH( IU_SMV2LOG )  
      CALL FLUSH( IU_DEBUG   ) 

      ! Return to MAIN program
      END SUBROUTINE CTM_FLUSH

!------------------------------------------------------------------------------

      SUBROUTINE DISPLAY_END_TIME

      !=================================================================
      ! Internal subroutine DISPLAY_END_TIME prints the ending time of
      ! the GEOS-CHEM simulation (bmy, 5/3/05)
      !=================================================================

      ! Local variables
      CHARACTER(LEN=16) :: STAMP

      ! Print system time stamp
      STAMP = SYSTEM_TIMESTAMP()
      WRITE( 6, 100 ) STAMP
 100  FORMAT( /, '===> SIMULATION END TIME: ', a, ' <===', / )

      ! Echo info
      WRITE ( 6, 3000 ) 
 3000 FORMAT
     &   ( /, '**************   E N D   O F   G E O S -- C H E M   ',
     &        '**************' )

      ! Return to MAIN program
      END SUBROUTINE DISPLAY_END_TIME

!------------------------------------------------------------------------------

      SUBROUTINE MET_FIELD_DEBUG

      !=================================================================
      ! Internal subroutine MET_FIELD_DEBUG prints out the maximum
      ! and minimum, and sum of DAO met fields for debugging 
      !=================================================================

      ! References to F90 modules
      USE DAO_MOD, ONLY : AD,       AIRDEN,  AIRVOL,   ALBD1,  ALBD2
      USE DAO_MOD, ONLY : ALBD,     AVGW,    BXHEIGHT, CLDFRC, CLDF     
      USE DAO_MOD, ONLY : CLDMAS,   CLDTOPS, DELP     
      USE DAO_MOD, ONLY : DTRAIN,   GWETTOP, HFLUX,    HKBETA, HKETA     
      USE DAO_MOD, ONLY : LWI,      MOISTQ,  OPTD,     OPTDEP, PBL      
      USE DAO_MOD, ONLY : PREACC,   PRECON,  PS1,      PS2,    PSC2     
      USE DAO_MOD, ONLY : RADLWG,   RADSWG,  RH,       SLP,    SNOW     
      USE DAO_MOD, ONLY : SPHU1,    SPHU2,   SPHU,     SUNCOS, SUNCOSB  
      USE DAO_MOD, ONLY : TMPU1,    TMPU2,   T,        TROPP,  TS       
      USE DAO_MOD, ONLY : TSKIN,    U10M,    USTAR,    UWND1,  UWND2     
      USE DAO_MOD, ONLY : UWND,     V10M,    VWND1,    VWND2,  VWND     
      USE DAO_MOD, ONLY : Z0,       ZMEU,    ZMMD,     ZMMU     

      ! Local variables
      INTEGER :: I, J, L, IJ

      !=================================================================
      ! MET_FIELD_DEBUG begins here!
      !=================================================================

      ! Define box to print out
      I  = 23
      J  = 34
      L  = 1
      IJ = ( ( J-1 ) * IIPAR ) + I

      !=================================================================
      ! Print out met fields at (I,J,L)
      !=================================================================
      IF ( ALLOCATED( AD       ) ) PRINT*, 'AD      : ', AD(I,J,L)        
      IF ( ALLOCATED( AIRDEN   ) ) PRINT*, 'AIRDEN  : ', AIRDEN(L,I,J) 
      IF ( ALLOCATED( AIRVOL   ) ) PRINT*, 'AIRVOL  : ', AIRVOL(I,J,L) 
      IF ( ALLOCATED( ALBD1    ) ) PRINT*, 'ALBD1   : ', ALBD1(I,J) 
      IF ( ALLOCATED( ALBD2    ) ) PRINT*, 'ALBD2   : ', ALBD2(I,J) 
      IF ( ALLOCATED( ALBD     ) ) PRINT*, 'ALBD    : ', ALBD(I,J) 
      IF ( ALLOCATED( AVGW     ) ) PRINT*, 'AVGW    : ', AVGW(I,J,L) 
      IF ( ALLOCATED( BXHEIGHT ) ) PRINT*, 'BXHEIGHT: ', BXHEIGHT(I,J,L) 
      IF ( ALLOCATED( CLDFRC   ) ) PRINT*, 'CLDFRC  : ', CLDFRC(I,J)
      IF ( ALLOCATED( CLDF     ) ) PRINT*, 'CLDF    : ', CLDF(L,I,J) 
      IF ( ALLOCATED( CLDMAS   ) ) PRINT*, 'CLDMAS  : ', CLDMAS(I,J,L) 
      IF ( ALLOCATED( CLDTOPS  ) ) PRINT*, 'CLDTOPS : ', CLDTOPS(I,J) 
      IF ( ALLOCATED( DELP     ) ) PRINT*, 'DELP    : ', DELP(L,I,J) 
      IF ( ALLOCATED( DTRAIN   ) ) PRINT*, 'DTRAIN  : ', DTRAIN(I,J,L) 
      IF ( ALLOCATED( GWETTOP  ) ) PRINT*, 'GWETTOP : ', GWETTOP(I,J) 
      IF ( ALLOCATED( HFLUX    ) ) PRINT*, 'HFLUX   : ', HFLUX(I,J) 
      IF ( ALLOCATED( HKBETA   ) ) PRINT*, 'HKBETA  : ', HKBETA(I,J,L) 
      IF ( ALLOCATED( HKETA    ) ) PRINT*, 'HKETA   : ', HKETA(I,J,L) 
      IF ( ALLOCATED( LWI      ) ) PRINT*, 'LWI     : ', LWI(I,J) 
      IF ( ALLOCATED( MOISTQ   ) ) PRINT*, 'MOISTQ  : ', MOISTQ(L,I,J) 
      IF ( ALLOCATED( OPTD     ) ) PRINT*, 'OPTD    : ', OPTD(L,I,J) 
      IF ( ALLOCATED( OPTDEP   ) ) PRINT*, 'OPTDEP  : ', OPTDEP(L,I,J) 
      IF ( ALLOCATED( PBL      ) ) PRINT*, 'PBL     : ', PBL(I,J) 
      IF ( ALLOCATED( PREACC   ) ) PRINT*, 'PREACC  : ', PREACC(I,J) 
      IF ( ALLOCATED( PRECON   ) ) PRINT*, 'PRECON  : ', PRECON(I,J) 
      IF ( ALLOCATED( PS1      ) ) PRINT*, 'PS1     : ', PS1(I,J) 
      IF ( ALLOCATED( PS2      ) ) PRINT*, 'PS2     : ', PS2(I,J) 
      IF ( ALLOCATED( PSC2     ) ) PRINT*, 'PSC2    : ', PSC2(I,J)
      IF ( ALLOCATED( RADLWG   ) ) PRINT*, 'RADLWG  : ', RADLWG(I,J)
      IF ( ALLOCATED( RADSWG   ) ) PRINT*, 'RADSWG  : ', RADSWG(I,J)
      IF ( ALLOCATED( RH       ) ) PRINT*, 'RH      : ', RH(I,J,L) 
      IF ( ALLOCATED( SLP      ) ) PRINT*, 'SLP     : ', SLP(I,J) 
      IF ( ALLOCATED( SNOW     ) ) PRINT*, 'SNOW    : ', SNOW(I,J) 
      IF ( ALLOCATED( SPHU1    ) ) PRINT*, 'SPHU1   : ', SPHU1(I,J,L) 
      IF ( ALLOCATED( SPHU2    ) ) PRINT*, 'SPHU2   : ', SPHU2(I,J,L) 
      IF ( ALLOCATED( SPHU     ) ) PRINT*, 'SPHU    : ', SPHU(I,J,L) 
      IF ( ALLOCATED( SUNCOS   ) ) PRINT*, 'SUNCOS  : ', SUNCOS(IJ) 
      IF ( ALLOCATED( SUNCOSB  ) ) PRINT*, 'SUNCOSB : ', SUNCOSB(IJ) 
      IF ( ALLOCATED( TMPU1    ) ) PRINT*, 'TMPU1   : ', TMPU1(I,J,L) 
      IF ( ALLOCATED( TMPU2    ) ) PRINT*, 'TMPU2   : ', TMPU2(I,J,L) 
      IF ( ALLOCATED( T        ) ) PRINT*, 'TMPU    : ', T(I,J,L)
      IF ( ALLOCATED( TROPP    ) ) PRINT*, 'TROPP   : ', TROPP(I,J) 
      IF ( ALLOCATED( TS       ) ) PRINT*, 'TS      : ', TS(I,J) 
      IF ( ALLOCATED( TSKIN    ) ) PRINT*, 'TSKIN   : ', TSKIN(I,J) 
      IF ( ALLOCATED( U10M     ) ) PRINT*, 'U10M    : ', U10M(I,J) 
      IF ( ALLOCATED( USTAR    ) ) PRINT*, 'USTAR   : ', USTAR(I,J) 
      IF ( ALLOCATED( UWND1    ) ) PRINT*, 'UWND1   : ', UWND1(I,J,L) 
      IF ( ALLOCATED( UWND2    ) ) PRINT*, 'UWND2   : ', UWND2(I,J,L) 
      IF ( ALLOCATED( UWND     ) ) PRINT*, 'UWND    : ', UWND(I,J,L) 
      IF ( ALLOCATED( V10M     ) ) PRINT*, 'V10M    : ', V10M(I,J)  
      IF ( ALLOCATED( VWND1    ) ) PRINT*, 'VWND1   : ', VWND1(I,J,L) 
      IF ( ALLOCATED( VWND2    ) ) PRINT*, 'VWND2   : ', VWND2(I,J,L) 
      IF ( ALLOCATED( VWND     ) ) PRINT*, 'VWND    : ', VWND(I,J,L) 
      IF ( ALLOCATED( Z0       ) ) PRINT*, 'Z0      : ', Z0(I,J) 
      IF ( ALLOCATED( ZMEU     ) ) PRINT*, 'ZMEU    : ', ZMEU(I,J,L) 
      IF ( ALLOCATED( ZMMD     ) ) PRINT*, 'ZMMD    : ', ZMMD(I,J,L) 
      IF ( ALLOCATED( ZMMU     ) ) PRINT*, 'ZMMU    : ', ZMMU(I,J,L) 

      ! Flush the output buffer
      CALL FLUSH( 6 )

      ! Return to MAIN program
      END SUBROUTINE MET_FIELD_DEBUG

!-----------------------------------------------------------------------------

      ! End of program
      END PROGRAM GEOS_CHEM

