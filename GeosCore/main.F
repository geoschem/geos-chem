!BOC
#if !defined(ESMF_)
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: geos_chem 
!
! !DESCRIPTION: Program GEOS\_CHEM is the main level driver program for the 
!  GEOS-Chem model of atmospheric chemistry and composition.
!\\
!\\
! !INTERFACE:
!
      PROGRAM GEOS_CHEM
!
! !USES:
!
      USE A3_READ_MOD,       ONLY : GET_A3_FIELDS
      USE A3_READ_MOD,       ONLY : OPEN_A3_FIELDS
      USE A3_READ_MOD,       ONLY : UNZIP_A3_FIELDS
      USE A6_READ_MOD,       ONLY : GET_A6_FIELDS
      USE A6_READ_MOD,       ONLY : OPEN_A6_FIELDS
      USE A6_READ_MOD,       ONLY : UNZIP_A6_FIELDS
      USE BENCHMARK_MOD,     ONLY : STDRUN
      ! (hotp 5/24/09) Modified for SOA from aroms
      !USE CARBON_MOD,        ONLY : WRITE_GPROD_APROD
      ! Add RECOMPUTE_OD to call AOD calculation twice (skim, 02/02/11)
      USE CHEMISTRY_MOD,     ONLY : DO_CHEMISTRY, RECOMPUTE_OD
      USE CONVECTION_MOD,    ONLY : DO_CONVECTION
      USE COMODE_MOD,        ONLY : INIT_COMODE
      USE GCKPP_COMODE_MOD,  ONLY : INIT_GCKPP_COMODE
      USE DIAG_MOD,          ONLY : DIAGCHLORO
      USE DIAG41_MOD,        ONLY : DIAG41,          ND41
      USE DIAG42_MOD,        ONLY : DIAG42,          ND42
      USE DIAG48_MOD,        ONLY : DIAG48,          ITS_TIME_FOR_DIAG48
      USE DIAG49_MOD,        ONLY : DIAG49,          ITS_TIME_FOR_DIAG49
      USE DIAG50_MOD,        ONLY : DIAG50,          DO_SAVE_DIAG50
      USE DIAG51_MOD,        ONLY : DIAG51,          DO_SAVE_DIAG51
      USE DIAG51b_MOD,       ONLY : DIAG51b,         DO_SAVE_DIAG51b
      USE DIAG63_MOD,        ONLY : DIAG63,          ITS_TIME_FOR_DIAG63
      USE DIAG_OH_MOD,       ONLY : PRINT_DIAG_OH
      USE DAO_MOD,           ONLY : AD,              AIRQNT  
      USE DAO_MOD,           ONLY : AVGPOLE,         CLDTOPS
      USE DAO_MOD,           ONLY : CONVERT_UNITS,   COPY_I6_FIELDS
      USE DAO_MOD,           ONLY : GET_COSINE_SZA,  INIT_DAO
      USE DAO_MOD,           ONLY : INTERP,          PS1
      USE DAO_MOD,           ONLY : PS2,             PSC2          
      USE DAO_MOD,           ONLY : T,               TS            
      USE DAO_MOD,           ONLY : SUNCOS,          SUNCOS_MID
      USE DAO_MOD,           ONLY : SUNCOS_MID_5hr
      USE DAO_MOD,           ONLY : MAKE_RH
      !Add MAKE_GTMM_RESTART for mercury simulation (ccc, 11/19/09)
      USE DEPO_MERCURY_MOD,  ONLY : MAKE_GTMM_RESTART, UPDATE_DEP
      USE DRYDEP_MOD,        ONLY : DO_DRYDEP
      USE EMISSIONS_MOD,     ONLY : DO_EMISSIONS
      USE ERROR_MOD,         ONLY : DEBUG_MSG,       ERROR_STOP
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
      USE LOGICAL_MOD,       ONLY : LEMIS,     LCHEM, LUNZIP,  LDUST
      USE LOGICAL_MOD,       ONLY : LLIGHTNOX, LPRT,  LSTDRUN, LSVGLB
      USE LOGICAL_MOD,       ONLY : LWAIT,     LTRAN, LUPBD,   LCONV
      USE LOGICAL_MOD,       ONLY : LWETD,     LTURB, LDRYD,   LMEGAN  
      USE LOGICAL_MOD,       ONLY : LDYNOCEAN, LSOA,  LVARTROP,LKPP
      USE LOGICAL_MOD,       ONLY : LLINOZ,    LWINDO
      ! Add LGTMM logical for mercury simulation (ccc, 11/19/09)
      USE LOGICAL_MOD,       ONLY : LGTMM
      USE MEGAN_MOD,         ONLY : INIT_MEGAN
      USE MEGAN_MOD,         ONLY : UPDATE_T_15_AVG
      USE MEGAN_MOD,         ONLY : UPDATE_T_DAY
      USE PBL_MIX_MOD,       ONLY : DO_PBL_MIX
      USE OCEAN_MERCURY_MOD, ONLY : MAKE_OCEAN_Hg_RESTART
      USE OCEAN_MERCURY_MOD, ONLY : READ_OCEAN_Hg_RESTART
      USE PLANEFLIGHT_MOD,   ONLY : PLANEFLIGHT
      USE PLANEFLIGHT_MOD,   ONLY : SETUP_PLANEFLIGHT 
      USE PRESSURE_MOD,      ONLY : INIT_PRESSURE
      USE PRESSURE_MOD,      ONLY : SET_FLOATING_PRESSURE, get_pedge
      ! add support for saving APROD, GPROD (dkh, 11/09/06)  
      USE SOAPROD_MOD,       ONLY : SET_SOAPROD, MAKE_SOAPROD_FILE
      USE SOAPROD_MOD,       ONLY : READ_SOAPROD_FILE
      ! hotp 5/25/09
      USE SOAPROD_MOD,       ONLY : FIRST_APRODGPROD
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
      USE TIME_MOD,          ONLY : ITS_TIME_FOR_BPCH
      USE TIME_MOD,          ONLY : SET_CT_CONV,      SET_CT_DYN
      USE TIME_MOD,          ONLY : SET_CT_EMIS,      SET_CT_CHEM
      USE TIME_MOD,          ONLY : SET_CT_DIAG
      USE TIME_MOD,          ONLY : SET_DIAGb,        SET_DIAGe
      USE TIME_MOD,          ONLY : SET_CURRENT_TIME, PRINT_CURRENT_TIME
      USE TIME_MOD,          ONLY : SET_ELAPSED_MIN,  SYSTEM_TIMESTAMP
      USE TIME_MOD,          ONLY : TIMESTAMP_DIAG
      USE TIME_MOD,          ONLY : GET_HOUR,         GET_MINUTE
      USE TIME_MOD,          ONLY : GET_FIRST_I6_TIME
      USE TPCORE_BC_MOD,     ONLY : SAVE_GLOBAL_TPCORE_BC
      USE TRACER_MOD,        ONLY : CHECK_STT, N_TRACERS, STT, TCVV
      USE TRACER_MOD,        ONLY : ITS_AN_AEROSOL_SIM
      USE TRACER_MOD,        ONLY : ITS_A_CH4_SIM
      USE TRACER_MOD,        ONLY : ITS_A_FULLCHEM_SIM
      USE TRACER_MOD,        ONLY : ITS_A_H2HD_SIM
      USE TRACER_MOD,        ONLY : ITS_A_MERCURY_SIM
      USE TRACER_MOD,        ONLY : ITS_A_TAGCO_SIM
      USE TRANSPORT_MOD,     ONLY : DO_TRANSPORT
      USE TROPOPAUSE_MOD,    ONLY : READ_TROPOPAUSE, CHECK_VAR_TROP
      USE TROPOPAUSE_MOD,    ONLY : DIAG_TROPOPAUSE
      USE RESTART_MOD,       ONLY : MAKE_RESTART_FILE, READ_RESTART_FILE
      USE UPBDFLX_MOD,       ONLY : DO_UPBDFLX,        UPBDFLX_NOY
      USE UVALBEDO_MOD,      ONLY : READ_UVALBEDO
      USE WETSCAV_MOD,       ONLY : INIT_WETSCAV,      DO_WETDEP
      USE XTRA_READ_MOD,     ONLY : GET_XTRA_FIELDS,   OPEN_XTRA_FIELDS
      USE XTRA_READ_MOD,     ONLY : UNZIP_XTRA_FIELDS
      USE ERROR_MOD,         ONLY : IT_IS_NAN, IT_IS_FINITE   !yxw
      USE ERROR_MOD,         ONLY : SAFE_DIV
      ! To save CSPEC_FULL restart (dkh, 02/12/09)
      USE LOGICAL_MOD,       ONLY : LSVCSPEC
      USE RESTART_MOD,       ONLY : MAKE_CSPEC_FILE
      ! Added (lin, 03/31/09)
      USE LOGICAL_MOD,       ONLY : LNLPBL
      USE VDIFF_MOD,         ONLY : DO_PBL_MIX_2
      USE LINOZ_MOD,         ONLY : LINOZ_READ
      USE TRACERID_MOD,      ONLY : IS_Hg2

      ! For GTMM for mercury simulations. (ccc, 6/7/10)
      USE WETSCAV_MOD,       ONLY : GET_WETDEP_IDWETD  
      USE MERCURY_MOD,       ONLY : PARTITIONHG

      ! For MERRA met fields (bmy, 8/19/10)
      USE MERRA_A1_MOD,      ONLY : GET_MERRA_A1_FIELDS
      USE MERRA_A1_MOD,      ONLY : OPEN_MERRA_A1_FIELDS
      USE MERRA_A3_MOD,      ONLY : GET_MERRA_A3_FIELDS
      USE MERRA_A3_MOD,      ONLY : OPEN_MERRA_A3_FIELDS
      USE MERRA_CN_MOD,      ONLY : GET_MERRA_CN_FIELDS
      USE MERRA_CN_MOD,      ONLY : OPEN_MERRA_CN_FIELDS
      USE MERRA_I6_MOD,      ONLY : GET_MERRA_I6_FIELDS_1
      USE MERRA_I6_MOD,      ONLY : GET_MERRA_I6_FIELDS_2
      USE MERRA_I6_MOD,      ONLY : OPEN_MERRA_I6_FIELDS
      USE TIME_MOD,          ONLY : GET_A1_TIME
      USE TIME_MOD,          ONLY : GET_FIRST_A1_TIME
      USE TIME_MOD,          ONLY : ITS_TIME_FOR_A1

      USE CMN_SIZE_MOD            ! Size parameters
      USE COMODE_LOOP_MOD         ! SMVGEAR common blocks
      USE CMN_DIAG_MOD            ! Diagnostic switches, NJDAY
      USE CMN_GCTM_MOD            ! Physical constants

      IMPLICIT NONE
#     include "define.h"
      
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
!           for 4 x 5, 2 x 2.5 global grids and 1 x 1 nested grids
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
!  13 Aug 2010 - R. Yantosca - Added ProTeX headers
!  13 Aug 2010 - R. Yantosca - Add modifications for MERRA (treat like GEOS-5)
!  19 Aug 2010 - R. Yantosca - Now call MERRA met field reader routines
!  02 Feb 2011 - S. Kim      - Call Compute_OD after wet deposition
!  05 Oct 2011 - R. Yantosca - Now get SUNCOS30 array from routine COSSZA
!  07 Oct 2011 - R. Yantosca - Rename SUNCOS30 to SUNCOS_MID, which is the
!                              cos(SZA) at the midpt of the chemistry timestep
!  07 Feb 2012 - M. Payer    - Replace call to COSSZA with GET_COSINE_SZA
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL            :: FIRST = .TRUE.
      LOGICAL            :: LXTRA 
      INTEGER            :: I,             IOS,   J,         K,      L
      INTEGER            :: N,             JDAY,  NDIAGTIME, N_DYN,  NN
      INTEGER            :: N_DYN_STEPS,   NSECb, N_STEP,    DATE(2)
      INTEGER            :: YEAR,          MONTH, DAY,       DAY_OF_YEAR
      INTEGER            :: SEASON,        NYMD,  NYMDb,     NHMS
      INTEGER            :: ELAPSED_SEC,   NHMSb, RC
      INTEGER            :: ELAPSED_TODAY, HOUR,  MINUTE
      REAL*8             :: TAU,           TAUb         
      REAL*8             :: HGPFRAC(IIPAR,JJPAR,LLPAR)
      CHARACTER(LEN=255) :: ZTYPE

      !=================================================================
      ! GEOS-CHEM starts here!                                            
      !=================================================================

      RC=0 ! Error flag. Default 0 is "no error"
      
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

      ! Initialize allocatable SMVGEAR/KPP arrays
      IF ( LEMIS .or. LCHEM ) THEN
         IF ( ITS_A_FULLCHEM_SIM() ) CALL INIT_COMODE
         IF ( ITS_AN_AEROSOL_SIM() ) CALL INIT_COMODE
         IF ( LKPP )  CALL INIT_GCKPP_COMODE( IIPAR, JJPAR, LLTROP,
     $        ITLOOP, NMTRATE, IGAS, RC )
         IF ( RC == 1 )
     $        CALL ERROR_STOP( "Alloc error", "INIT_GCKPP_COMODE" )
         IF ( LPRT ) CALL DEBUG_MSG( '### MAIN: a INIT_COMODE' )
      ENDIF
         
      ! Added to read input file for linoz strat (dbj, jliu, bmy, 10/16/09)
      IF ( LLINOZ ) CALL LINOZ_READ

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

#if   defined( MERRA )

      !=================================================================
      !    *****  R E A D   M E R R A   M E T   F I E L D S  *****
      !    *****  At the start of the GEOS-Chem simulation   *****
      !
      !    Handle MERRA met fields separately from other met products
      !=================================================================

      ! Open constant fields
      DATE = (/ 20000101, 000000 /)
      CALL OPEN_MERRA_CN_FIELDS( DATE(1), DATE(2) )
      CALL GET_MERRA_CN_FIELDS(  DATE(1), DATE(2) )
      IF ( LPRT ) CALL DEBUG_MSG( '### MAIN: a 1st MERRA CN TIME' )

      ! Open and read A-1 fields
      DATE = GET_FIRST_A1_TIME()
      CALL OPEN_MERRA_A1_FIELDS( DATE(1), DATE(2), RESET=.TRUE. )
      CALL GET_MERRA_A1_FIELDS(  DATE(1), DATE(2) )
      IF ( LPRT ) CALL DEBUG_MSG( '### MAIN: a 1st MERRA A1 TIME' )

      ! Open and read A-3 fields
      DATE = GET_FIRST_A3_TIME()
      CALL OPEN_MERRA_A3_FIELDS( DATE(1), DATE(2) )
      CALL GET_MERRA_A3_FIELDS(  DATE(1), DATE(2) )
      IF ( LPRT ) CALL DEBUG_MSG( '### MAIN: a 1st MERRA A3 TIME' )

      ! Open & read I-6 fields
      DATE = GET_FIRST_I6_TIME()
      CALL OPEN_MERRA_I6_FIELDS(  DATE(1), DATE(2) )
      CALL GET_MERRA_I6_FIELDS_1( DATE(1), DATE(2) )
      IF ( LPRT ) CALL DEBUG_MSG( '### MAIN: a 1st I6 TIME' )

#else

      !=================================================================
      !    *****      U N Z I P   M E T   F I E L D S        *****
      !    ***** At at the start of the GEOS-Chem simulation *****
      !
      !   Here we unzip the initial GEOS-3, GEOS-4, GEOS-5, GCAP data
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
      !      *****      R E A D   M E T   F I E L D S       *****
      !      ***** At the start of the GEOS-Chem simulation *****
      !
      !  Here we read in the initial GEOS-3, GEOS-4, GEOS-5, GCAP data
      !=================================================================

      ! Open and read A-3 fields
      DATE = GET_FIRST_A3_TIME()
      CALL OPEN_A3_FIELDS( DATE(1), DATE(2), RESET=.TRUE. )
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
      DATE = GET_FIRST_I6_TIME()
      CALL OPEN_I6_FIELDS(  DATE(1), DATE(2) )
      CALL GET_I6_FIELDS_1( DATE(1), DATE(2) )
      IF ( LPRT ) CALL DEBUG_MSG( '### MAIN: a 1st I6 TIME' )
      
#if   defined( GEOS_3 )

      !-----------------------------------------------------------------
      ! Read additional fields for GEOS-3 meteorology
      !-----------------------------------------------------------------

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

      !-----------------------------------------------------------------
      ! Read additional fields for GCAP meteorology
      !-----------------------------------------------------------------

      ! Read GCAP PHIS and LWI fields (if necessary)
      CALL OPEN_GCAP_FIELDS
      CALL GET_GCAP_FIELDS

      ! Remove temporary file (if necessary)
      IF ( LUNZIP ) THEN
         CALL UNZIP_GCAP_FIELDS( 'remove date' )
      ENDIF

#endif

#endif

      !=================================================================
      !        ***** I N I T I A L I Z A T I O N  continued *****
      !=================================================================

      ! Compute avg surface pressure near polar caps
      CALL AVGPOLE( PS1 )
      IF ( LPRT ) CALL DEBUG_MSG( '### MAIN: a AVGPOLE' )

      ! Call AIRQNT to compute air mass quantities from PS1
      CALL SET_FLOATING_PRESSURE( PS1 )
      IF ( LPRT ) CALL DEBUG_MSG( '### MAIN: a SET_FLT_PRS' )

      CALL AIRQNT
      IF ( LPRT ) CALL DEBUG_MSG( '### MAIN: a AIRQNT' )

      ! Compute lightning NOx emissions [molec/box/6h]
      IF ( LLIGHTNOX ) THEN
         CALL LIGHTNING
         IF ( LPRT ) CALL DEBUG_MSG( '### MAIN: a LIGHTNING' )
      ENDIF

      ! Read land types and fractions from "vegtype.global"
      CALL RDLAND   
      IF ( LPRT ) CALL DEBUG_MSG( '### MAIN: a RDLAND' )

      ! Initialize PBL quantities but do not do mixing
      ! Add option for non-local PBL (Lin, 03/31/09) 
      IF ( .NOT. LNLPBL ) THEN
        CALL DO_PBL_MIX( .FALSE. )
        IF ( LPRT ) CALL DEBUG_MSG( '### MAIN: a TURBDAY:1' )
      ELSE
        CALL DO_PBL_MIX_2( .FALSE. )
        IF ( LPRT ) CALL DEBUG_MSG( '### MAIN: a NLPBL 1' )
      ENDIF

      !=================================================================
      !       *****  I N I T I A L   C O N D I T I O N S *****
      !=================================================================

      ! Read initial tracer conditions
      CALL READ_RESTART_FILE( NYMDb, NHMSb )
      IF ( LPRT ) CALL DEBUG_MSG( '### MAIN: a READ_RESTART_FILE' )

      ! add support for making restart files of APROD and GPROD (dkh, 11/09/06)  
      IF ( LSOA ) THEN

         !! use this to make initial soaprod files  
         !CALL SET_SOAPROD
         !CALL FIRST_APRODGPROD()
         !CALL MAKE_SOAPROD_FILE( GET_NYMDb(), GET_NHMSb() )
         !goto 9999
         !!

         CALL SET_SOAPROD
         CALL READ_SOAPROD_FILE( GET_NYMDb(), GET_NHMSb() )

      ENDIF


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

      ! NSECb is # of seconds (measured from 00 GMT today) 
      ! at the start of this 6-hr timestepping loop.
      ! NOTE: Assume we start at the head of each minute (i.e. SECONDS=0)
      HOUR   = GET_HOUR()
      HOUR   = ( HOUR / 6 ) * 6
      MINUTE = GET_MINUTE()
      NSECb  = ( HOUR * 3600 ) + ( MINUTE * 60 )

      ! Get dynamic timestep in seconds
      N_DYN  = 60d0 * GET_TS_DYN()

      !=================================================================
      !     ***** D Y N A M I C   T I M E S T E P   L O O P *****
      !=================================================================
      DO N_STEP = 1, N_DYN_STEPS
    
         ! Compute & print time quantities at start of dyn step
         CALL SET_CURRENT_TIME
         CALL PRINT_CURRENT_TIME

         ! Set time variables for dynamic loop
         DAY_OF_YEAR   = GET_DAY_OF_YEAR()
         ELAPSED_SEC   = GET_ELAPSED_SEC()
         MONTH         = GET_MONTH()
         NHMS          = GET_NHMS()
         NYMD          = GET_NYMD()
         HOUR          = GET_HOUR()
         MINUTE        = GET_MINUTE()
         TAU           = GET_TAU()
         YEAR          = GET_YEAR()
         SEASON        = GET_SEASON()
         ELAPSED_TODAY = ( HOUR * 3600 ) + ( MINUTE * 60 )

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
            !    *****  W R I T E   R E S T A R T   F I L E S  *****
            !===========================================================
            IF ( LSVGLB ) THEN

               ! Make atmospheric restart file
               CALL MAKE_RESTART_FILE( NYMD, NHMS, TAU )
                  
               ! Make ocean mercury restart file
               IF ( ITS_A_MERCURY_SIM() .and. LDYNOCEAN ) THEN
                  CALL MAKE_OCEAN_Hg_RESTART( NYMD, NHMS, TAU )
               ENDIF

               ! (dkh, 11/09/06)  
!               ! Save SOA quantities GPROD & APROD
!               IF ( LSOA .and. LCHEM ) THEN 
!                  CALL WRITE_GPROD_APROD( NYMD, NHMS, TAU )
!               ENDIF

               IF ( LSOA .and. LCHEM ) THEN
                  CALL MAKE_SOAPROD_FILE( GET_NYMD(), GET_NHMS() )

                  !### Debug
                  IF ( LPRT ) THEN
                     CALL DEBUG_MSG( '### MAIN: a MAKE_SOAPROD_FILE' )
                  ENDIF
               ENDIF

               ! Save species concentrations (CSPEC_FULL). (dkh, 02/12/09)
               IF ( LCHEM .and. LSVCSPEC ) THEN 
                  CALL MAKE_CSPEC_FILE( NYMD, NHMS )
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

         !=============================================================
         !   ***** W R I T E   MERCURY RESTART  F I L E *****
         !     ***** MUST be done after call to diag3 *****
         !=============================================================
         ! Make land restart file: for GTMM runs only, beginning of each 
         ! month but not start of the run.
         IF ( LGTMM .AND. ITS_A_NEW_MONTH() .AND. NYMD /= NYMDb ) THEN
            IF (.NOT.( ITS_TIME_FOR_BPCH() )) THEN
               N = 1
               NN = GET_WETDEP_IDWETD( N )
               DO WHILE( .NOT.(IS_Hg2( NN )) )
               
                  N = N + 1
                  ! Tracer number
                  NN = GET_WETDEP_IDWETD( N )

               ENDDO

               CALL UPDATE_DEP( N )
            ENDIF
            CALL MAKE_GTMM_RESTART( NYMD, NHMS, TAU )
         ENDIF

         !==============================================================
         !       ***** T E S T   F O R   E N D   O F   R U N *****
         !==============================================================
         IF ( ITS_TIME_FOR_EXIT() ) GOTO 9999

#if   defined( MERRA )

         !==============================================================
         !    ***** R E A D   M E R R A   A - 1   F I E L D S *****
         !
         !    The MERRA archive contains hourly surface data fields.
         !==============================================================
         IF ( ITS_TIME_FOR_A1() ) THEN

            ! Get the date/time for the next A-3 data block
            DATE = GET_A1_TIME()

            ! Open & read A-3 fields
            CALL OPEN_MERRA_A1_FIELDS( DATE(1), DATE(2) )
            CALL GET_MERRA_A1_FIELDS ( DATE(1), DATE(2) )

            !%%% NEED TO UPDATE FOR MERRA %%%
            ! Update daily mean temperature archive for MEGAN biogenics
            !IF ( LMEGAN ) CALL UPDATE_T_DAY 
         ENDIf

         !==============================================================
         !    ***** R E A D   M E R R A   A - 3   F I E L D S *****
         !
         !     The MERRA archive contains 3-hourly 3-D data fields.
         !==============================================================
         IF ( ITS_TIME_FOR_A3() ) THEN
            
            ! Get the date/time for the next A-6 data block
            DATE = GET_A3_TIME()

            ! Open and read A-6 fields
            CALL OPEN_MERRA_A3_FIELDS( DATE(1), DATE(2) )
            CALL GET_MERRA_A3_FIELDS ( DATE(1), DATE(2) )

            ! Since CLDTOPS is an A-3 field, update the
            ! lightning NOx emissions [molec/box/6h]
            IF ( LLIGHTNOX ) CALL LIGHTNING
         ENDIF

         !==============================================================
         !    ***** R E A D   M E R R A   I - 6   F I E L D S *****
         !
         !    The MERRA archive contains 6-hourly instantaneous data.
         !==============================================================
         IF ( ITS_TIME_FOR_I6() ) THEN

            ! Get the date/time for the next I-6 data block
            DATE = GET_I6_TIME()

            ! Open and read files
            CALL OPEN_MERRA_I6_FIELDS ( DATE(1), DATE(2) )
            CALL GET_MERRA_I6_FIELDS_2( DATE(1), DATE(2) )

            ! Compute avg pressure at polar caps 
            CALL AVGPOLE( PS2 )
         ENDIF

#else

         !===============================================================
         !         ***** U N Z I P   M E T   F I E L D S *****
         !
         !      Some met data (except MERRA) are stored compressed.
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
         !
         !  All met data (except MERRA) contain 3-hourly surface data.
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
         !
         !      All other met fields contain 6-hourly 3-D data. 
         !==============================================================
         IF ( ITS_TIME_FOR_A6() ) THEN
            
            ! Get the date/time for the next A-6 data block
            DATE = GET_A6_TIME()

            ! Open and read A-6 fields
            CALL OPEN_A6_FIELDS( DATE(1), DATE(2) )
            CALL GET_A6_FIELDS(  DATE(1), DATE(2) )

            ! Since CLDTOPS is an A-6 field, update the
            ! lightning NOx emissions [molec/box/6h]
            IF ( LLIGHTNOX ) CALL LIGHTNING
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

#endif

         !==============================================================
         ! ***** M O N T H L Y   O R   S E A S O N A L   D A T A *****
         !==============================================================

         ! UV albedoes
         IF ( LCHEM .and. ITS_A_NEW_MONTH() ) THEN
            CALL READ_UVALBEDO( MONTH )
         ENDIF

         ! Fossil fuel emissions (SMVGEAR)
         ! NOTE: Only call ANTHROEMS if it's time to do emissions;
         ! this will function properly for the central chemistry time step
         ! (bmy, 10/14/11)
         IF ( ITS_A_FULLCHEM_SIM() .or. ITS_A_TAGCO_SIM() ) THEN
            IF ( ITS_TIME_FOR_EMIS() ) THEN 
               IF ( LEMIS .and. ITS_A_NEW_SEASON() ) THEN
                  CALL ANTHROEMS( SEASON )              
               ENDIF
            ENDIF
         ENDIF

         !==============================================================
         !              ***** D A I L Y   D A T A *****
         !==============================================================
         IF ( ITS_A_NEW_DAY() ) THEN 

            ! Read leaf-area index (needed for drydep)
            ! Now include the year (mpb,11/19/09)
            CALL RDLAI( DAY_OF_YEAR, MONTH, YEAR )

            ! For MEGAN biogenics ...
            IF ( LMEGAN ) THEN

               ! Read daily leaf-area-index
               ! Now include the year (mpb,11/19/09)
               CALL RDISOLAI( DAY_OF_YEAR, MONTH, YEAR )

               ! Compute 15-day average temperature for MEGAN
               CALL UPDATE_T_15_AVG
            ENDIF

            ! For mercury simulations ...
            IF ( ITS_A_MERCURY_SIM() ) THEN

               ! Read AVHRR daily leaf-area-index
               CALL RDISOLAI( DAY_OF_YEAR, MONTH, YEAR )

            ENDIF 
               
            ! Also read soil-type info for fullchem simulation
            IF ( ITS_A_FULLCHEM_SIM() .or. ITS_A_H2HD_SIM() ) THEN
               CALL RDSOIL 
            ENDIF

            !### Debug
            IF ( LPRT ) CALL DEBUG_MSG ( '### MAIN: a DAILY DATA' )
         ENDIF

         !==============================================================
         !   ***** I N T E R P O L A T E   Q U A N T I T I E S *****   
         !==============================================================
         
         ! Interpolate I-6 fields to current dynamic timestep, 
         ! based on their values at NSEC and NSEC+N_DYN
         CALL INTERP( NSECb, ELAPSED_TODAY, N_DYN )

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

         ! Compute the cosine of the solar zenith angle array
         ! SUNCOS         = cos(SZA) at the current time
         ! SUNCOS_MID     = cos(SZA) at the midpt of the chem timestep
         ! SUNCOS_MID_5hr = cos(SZA) at the midpt of the chem timestep 5hrs ago
!------------------------------------------------------------------------------
! Prior to 2/7/12:
!         CALL COSSZA( DAY_OF_YEAR, SUNCOS, SUNCOS_MID )
!------------------------------------------------------------------------------
         CALL GET_COSINE_SZA( SUNCOS, SUNCOS_MID, SUNCOS_MID_5hr )

         ! Compute tropopause height for ND55 diagnostic
         IF ( ND55 > 0 ) CALL DIAG_TROPOPAUSE

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

         ! Get averaging intervals for local-time diagnostics
         ! (NOTE: maybe improve this later on)
         ! Placed after interpolation to get correct value of TROPP. 
         ! (ccc, 12/9/08)
         CALL DIAG_2PM

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

            ! Output BC's
            ! Save boundary conditions (global grid) for future nested run
            IF ( LWINDO ) CALL SAVE_GLOBAL_TPCORE_BC

            ! Call the appropritate version of TPCORE
            IF ( LTRAN ) CALL DO_TRANSPORT               

            ! Reset air mass quantities
            CALL AIRQNT

            ! Repartition [NOy] species after transport
            IF ( LUPBD .and. ITS_A_FULLCHEM_SIM() ) THEN
               CALL UPBDFLX_NOY( 2 )
            ENDIF   

#if   !defined( GEOS_5 ) && !defined( MERRA )
            ! Get relative humidity (after recomputing pressures)
            ! NOTE: for GEOS-5 we'll read this from disk instead
            CALL MAKE_RH
#endif

            ! Initialize wet scavenging and wetdep fields after
            ! the airmass quantities are reset after transport
            IF ( LCONV .or. LWETD ) CALL INIT_WETSCAV
         ENDIF


         !-------------------------------
         ! Test for emission timestep
         !-------------------------------
         IF ( LNLPBL .AND. ITS_TIME_FOR_EMIS() ) THEN

            IF ( ITS_TIME_FOR_UNIT() )
     &         CALL CONVERT_UNITS( 2, N_TRACERS, TCVV, AD, STT ) ! v/v -> kg

            ! Increment emission counter
            CALL SET_CT_EMIS( INCREMENT=.TRUE. )

            !========================================================
            !         ***** D R Y   D E P O S I T I O N *****
            !========================================================
            IF ( LDRYD .and. ( .not. ITS_A_H2HD_SIM() ) ) CALL DO_DRYDEP

            !========================================================
            !             ***** E M I S S I O N S *****
            !========================================================
            IF ( LEMIS ) CALL DO_EMISSIONS

            IF ( ITS_TIME_FOR_UNIT() )
     &         CALL CONVERT_UNITS( 1, N_TRACERS, TCVV, AD, STT ) ! kg -> v/v

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
            ! Add option for non-local PBL. (Lin, 03/31/09)
            IF ( .NOT. LNLPBL ) THEN
               CALL DO_PBL_MIX( LTURB )
               IF ( LPRT ) CALL DEBUG_MSG( '### MAIN: a TURBDAY:1' )
            ELSE
               CALL DO_PBL_MIX_2( LTURB )
               IF ( LPRT ) CALL DEBUG_MSG( '### MAIN: a NLPBL 2' )
            ENDIF

            !### Debug
            IF ( LPRT ) CALL DEBUG_MSG( '### MAIN: a TURBDAY:2' )

            !===========================================================
            !        ***** C L O U D   C O N V E C T I O N *****
            !===========================================================
            IF ( LCONV ) THEN
               
               ! Partition Hg(II) between aerosol and gas
               IF ( ITS_A_MERCURY_SIM() ) THEN
                  CALL PARTITIONHG( 1, STT, HGPFRAC )
               ENDIF
      
               CALL DO_CONVECTION

               ! Return all reactive particulate Hg(II) to total Hg(II) tracer
               IF ( ITS_A_MERCURY_SIM() ) THEN
                  CALL PARTITIONHG( 2, STT, HGPFRAC )
               ENDIF
      
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

         !-------------------------------
         ! Test for emission timestep
         !-------------------------------
         IF ( ( .NOT. LNLPBL ) .AND. ITS_TIME_FOR_EMIS() ) THEN

            ! Increment emission counter
            CALL SET_CT_EMIS( INCREMENT=.TRUE. )

            !========================================================
            !         ***** D R Y   D E P O S I T I O N *****
            !========================================================
            IF ( LDRYD .and. ( .not. ITS_A_H2HD_SIM() ) ) CALL DO_DRYDEP

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
         IF ( LWETD .and. ITS_TIME_FOR_DYN() ) THEN

            ! Add partition Hg(II) between aerosol and gas
            IF ( ITS_A_MERCURY_SIM() ) THEN
               CALL PARTITIONHG( 1, STT, HGPFRAC )
            ENDIF            
     
            CALL DO_WETDEP
            
            ! Return all reactive particulate Hg(II) to total Hg(II) tracer
            IF ( ITS_A_MERCURY_SIM() ) THEN
               CALL PARTITIONHG( 2, STT, HGPFRAC )
            ENDIF 

         ENDIF

         !==============================================================
         !      ***** U P D A T E  O P T I C A L  D E P T H *****          
         !==============================================================
         ! Recalculate the optical depth at the wavelength specified by
         ! jv_spec_aod.dat This must be done before the call to any
         ! diagnostic and only on a chemistry timestep.
         ! (skim, 02/05/11)
         IF ( ITS_TIME_FOR_CHEM() ) THEN
            CALL RECOMPUTE_OD
         ENDIF

         !==============================================================
         !   ***** I N C R E M E N T   E L A P S E D   T I M E *****
         !============================================================== 
         ! Moved before diagnostics to count the last timestep as done.
         ! Need to save timestamps for filenames.
         ! (ccc, 5/13/09)
 
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
         ENDIF

         CALL TIMESTAMP_DIAG
         CALL SET_ELAPSED_MIN
         CALL SET_CURRENT_TIME
         IF ( LPRT ) CALL DEBUG_MSG( '### MAIN: after SET_ELAPSED_MIN' )

         !==============================================================
         !       ***** A R C H I V E   D I A G N O S T I C S *****
         !==============================================================
         IF ( ITS_TIME_FOR_DIAG() ) THEN

            !### Debug
            IF ( LPRT ) CALL DEBUG_MSG( '### MAIN: b DIAGNOSTICS' )

            ! Accumulate several diagnostic quantities
            CALL DIAG1
            IF ( LPRT ) CALL DEBUG_MSG( '### MAIN: after DIAG1' )

            ! ND41: save PBL height in 1200-1600 LT (amf)
            ! (for comparison w/ Holzworth, 1967)
            IF ( ND41 > 0 ) CALL DIAG41
            IF ( LPRT ) CALL DEBUG_MSG( '### MAIN: after DIAG41' )

            ! ND42: SOA concentrations [ug/m3]
            IF ( ND42 > 0 ) CALL DIAG42
            IF ( LPRT ) CALL DEBUG_MSG( '### MAIN: after DIAG42' )

            ! 24-hr timeseries
            IF ( DO_SAVE_DIAG50 ) CALL DIAG50

            ! Increment diagnostic timestep counter. (ccc, 5/13/09)
            CALL SET_CT_DIAG( INCREMENT=.TRUE. )

            ! Plane following diagnostic
            IF ( ND40 > 0 ) THEN 
               
               print*, 'Call planeflight'
               ! Archive data along the flight track
               CALL PLANEFLIGHT
            ENDIF
            IF ( LPRT ) CALL DEBUG_MSG( '### MAIN: after DIAG40' )

            !### Debug
            IF ( LPRT ) CALL DEBUG_MSG( '### MAIN: a DIAGNOSTICS' )
         ENDIF

         !==============================================================
         !   ***** T I M E S E R I E S   D I A G N O S T I C S  *****
         !
         ! NOTE: Since we are saving soluble tracers, we must move
         !       the ND40, ND49, and ND52 timeseries diagnostics
         !       to after the call to DO_WETDEP (bmy, 4/22/04)
         !============================================================== 

         IF ( LPRT ) CALL DEBUG_MSG( '### MAIN: before TIMESERIES' )

         ! Station timeseries
         IF ( ITS_TIME_FOR_DIAG48() ) CALL DIAG48
         IF ( LPRT ) CALL DEBUG_MSG( '### MAIN: after DIAG48' )

         ! 3-D timeseries
         IF ( ITS_TIME_FOR_DIAG49() ) CALL DIAG49
         IF ( LPRT ) CALL DEBUG_MSG( '### MAIN: after DIAG49' )

         ! Ship timeseries
         IF ( ITS_TIME_FOR_DIAG63() ) CALL DIAG63
         IF ( LPRT ) CALL DEBUG_MSG( '### MAIN: after DIAG63' )

         ! Morning or afternoon timeseries
         IF ( DO_SAVE_DIAG51 ) CALL DIAG51 
         IF ( DO_SAVE_DIAG51b ) CALL DIAG51b 
         IF ( LPRT ) CALL DEBUG_MSG( '### MAIN: after DIAG51' )

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

         ! Check for NaN, Negatives, Infinities in STT each time diag are
         ! saved. (ccc, 5/13/09)
         IF ( ITS_TIME_FOR_DIAG() ) THEN
            CALL CHECK_STT( 'End of Dynamic Loop' )
         ENDIF
          
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

#if defined( GTMM_Hg )
      ! Deallocate arrays from GTMM model for mercury simulation
      IF ( LGTMM ) CALL CleanupCASAarrays
#endif

      IF ( LPRT ) CALL DEBUG_MSG( '### MAIN: a CLEANUP' )

      ! Print ending time of simulation
      CALL DISPLAY_END_TIME

      CONTAINS
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: display_grid_and_model
!
! !DESCRIPTION: Internal Subroutine DISPLAY\_GRID\_AND\_MODEL displays the 
!  appropriate messages for the given model grid and machine type.  It also 
!  prints the starting time and date (local time) of the GEOS-Chem simulation.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE DISPLAY_GRID_AND_MODEL
! 
! !REVISION HISTORY: 
!  02 Dec 2003 - R. Yantosca - Initial version
!  13 Aug 2010 - R. Yantosca - Added ProTeX headers
!  13 Aug 2010 - R. Yantosca - Added extra output
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
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

#elif defined( GRID05x0666 )
      WRITE( 6, '(a)' ) 
     &    REPEAT( '*', 13 )                                          // 
     &    '   S T A R T I N G   0.5 x 0.666   G E O S -- C H E M   ' //
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
      WRITE( 6, '(a)' ) 'Created w/ LINUX/IFORT compiler'
      WRITE( 6, '(a)' ) 'Use ifort -V to print version information'
#elif defined( SGI_MIPS  )
      WRITE( 6, '(a)' ) 'Created w/ SGI MIPSpro compiler'
#elif defined( SPARC     )
      WRITE( 6, '(a)' ) 'Created w/ Sun/SPARC compiler'
#endif

      !-----------------------
      ! Print met field info
      !-----------------------
#if   defined( GEOS_3 )
      WRITE( 6, '(a)' ) 'Using GMAO GEOS-3 met fields'
#elif defined( GEOS_4 )
      WRITE( 6, '(a)' ) 'Using GMAO GEOS-4 met fields'
#elif defined( GEOS_5 )
      WRITE( 6, '(a)' ) 'Using GMAO GEOS-5 met fields'
#elif defined( MERRA )
      WRITE( 6, '(a)' ) 'Using GMAO MERRA met fields'
#elif defined( GCAP  )
      WRITE( 6, '(a)' ) 'Using GCAP/GISS met fields'
#endif

      !-----------------------
      ! System time stamp
      !-----------------------
      STAMP = SYSTEM_TIMESTAMP()
      WRITE( 6, 100 ) STAMP
 100  FORMAT( /, '===> SIMULATION START TIME: ', a, ' <===', / )

      END SUBROUTINE DISPLAY_GRID_AND_MODEL
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
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
      SUBROUTINE CTM_FLUSH
! 
! !REVISION HISTORY: 
!  31 Aug 2000 - R. Yantosca - Initial version
!  13 Aug 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
      CALL FLUSH( IU_ND48    )  
      CALL FLUSH( IU_BPCH    )  
      CALL FLUSH( IU_SMV2LOG )  
      CALL FLUSH( IU_DEBUG   ) 

      END SUBROUTINE CTM_FLUSH
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
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
      SUBROUTINE DISPLAY_END_TIME
! 
! !REVISION HISTORY: 
!  03 May 2005 - R. Yantosca - Initial version
!  13 Aug 2010 - R. Yantosca - Added ProTeX headers
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
 100  FORMAT( /, '===> SIMULATION END TIME: ', a, ' <===', / )

      ! Echo info
      WRITE ( 6, 3000 ) 
 3000 FORMAT
     &   ( /, '**************   E N D   O F   G E O S -- C H E M   ',
     &        '**************' )

      END SUBROUTINE DISPLAY_END_TIME
!EOC
      END PROGRAM GEOS_CHEM
#endif
