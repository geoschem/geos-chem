C $Id: main.f,v 1.8 2004/01/27 21:25:08 bmy Exp $
C $Log: main.f,v $
C Revision 1.8  2004/01/27 21:25:08  bmy
C GEOS-CHEM v6-01-05, includes the following modifications:
C - Updated deep & shallow convection code for GEOS-4 met fields
C - Bug fix in routine NFCLDMX to allow 2x25 runs on Altix
C - Updated code to read GEOS-4 met fields; now looks for ident strings
C - Bug fix in sulfate emissions code: convert GEOS-4 PBL from [m] to [hPa]
C - Updated comments, cosmetic changes
C
C Revision 1.7  2003/12/11 21:54:11  bmy
C GEOS-CHEM v6-01-04, includes the following modifications:
C - RADSWG, SNOW, USTAR, Z0 are now 2-D met field arrays in "dao_mod.f"
C - Renamed LRERD to LUNZIP in "input.geos"
C - Can now read either zipped or unzipped met field files
C - Now specify year/month info in the GEOS_1_DIR etc variables ("setup.f")
C
C Revision 1.6  2003/12/05 21:14:03  bmy
C GEOS-CHEM v6-01-03, includes the following modifications:
C - Added LINUX_IFC, LINUX_EFC C-preprocessor switches
C - Renamed SGI C-preprocessor switch to SGI_MIPS
C - Now includes Makefile.ifc and Makefile.efc for INTEL compiler
C - Updated "build" script to automatically select proper makefile2
C
C Revision 1.5  2003/11/06 21:07:18  bmy
C GEOS-CHEM v6-01-02, includes the following modifications:
C - Can now vertically regrid GEOS-4 55L to 30L on-the-fly
C - Bug fix in diag41 afternoon PBL diagnostic
C - Added File I/O error checks for IBM and INTEL_FC compilers
C
C Revision 1.4  2003/10/30 16:17:18  bmy
C GEOS-CHEM v6-01-01, includes the following modifications:
C - The fvDAS/TPCORE code conserves global airmass for GEOS-3/GEOS-4 winds
C - Improved 222Rn emissions so that they now reflect
C - Bug fix in "soilnoxems.f": now pass SUNCOS to "soilcrf.f"
C - Updated code to support for the IFC fortran compiler
C - Improved log file output formatting slightly
C - Removed obsolete code from several routines
C
C Revision 1.3  2003/10/01 20:32:22  bmy
C GEOS-CHEM v5.07.07, includes the following modifications:
C - GAMMA coeff for N2O5 hydrolysis is now a function of T, RH, aerosol type
C - Now apply drydep to entire PBL in single-tracer Ox routine chemo3.f
C - Added OpenMP parallelization commands to ND65 routine diagpl.f
C - Other minor bug fixes and updates
C
C Revision 1.2  2003/07/08 15:27:48  bmy
C GEOS-CHEM v5.07.01, now contains:
C - interannually varying latitudinal CH4 gradient full chemistry
C - removal of obsolete code; updates to documentation
C
C Revision 1.1.2.3  2003/06/30 18:32:06  bmy
C GEOS-CHEM v5.06, now contains:
C - Code for reading in fvDAS met fields
C - new (and yet untested) fvDAS convection code
C - Compilation patches for LINUX/PGI platform
C - Minor bug fixes in "rdust.f" and "tagged_co_mod.f"
C - Removed obsolete code; updated comments
C
C Revision 1.1.2.2  2003/06/24 15:28:23  bmy
C Synch up of NCCS code plus fvDAS met field modifications (bmy, 6/24/03)
C
C Revision 1.1.1.1  2003/06/13 16:35:02  bmy
C GEOS-CHEM v5.05.03, now contains:
C - new SMVGEAR II solver package for full chemistry
C - new TPCORE for GEOS-4/fvDAS wind fields
C - Phil Cameron-Smith P-Fixer to ensure new TPCORE conserves mass
C - lightning NOx for 1x1 China nested grid scaled to account for 1x1 gri
C - minor bug fixes in Rn-Pb-Be and Tagged CO
C - removed obsolete code, updated comments, cosmetic changes
C
C Revision 5.4  2003/04/08 20:27:32  bmy
C GEOS-CHEM v5.04, now contains:
C - new module for grid box lat, lon, surface area (grid_mod.f)
C - new module for date and time information (time_mod.f)
C - local time now computed consistently everywhere (time_mod.f)
C - new drivers for emission/chem (emission_mod.f, chemistry_mod.f)
C - new module for nested grid transport (tpcore_window_mod.f)
C - new module for nested grid boundary conditions (tpcore_bc_mod.f)
C - bug fix in FAST-J in solar constant computation
C - removed LOTS of historical baggage
C - updated look & feel of GEOS-CHEM output
C
C Revision 5.3  2003/02/03 15:04:12  bmy
C GEOS-CHEM v5.03, now contains:
C - new coupled fullchem/sulfate aerosol simulation
C - new offline sulfate simulation
C - drydep code now packaged into "drydep_mod.f"
C - wetdep code now packaged into "wetscav_mod.f"
C - new module "tracerid_mod.f" for tracer ID flags
C - new module "rpmares_mod.f" for aerosol thermodyn. equilib.
C - improved documentation and removed historical baggage
C
C Revision 5.2  2002/09/18 14:00:06  bmy
C GEOS-CHEM v5.02, now contains:
C - floating surface pressure (cf. Phil Cameron-Smith)
C - grid box pressure center/edges now computed for hybrid grid
C - ocean production of acetone for GEOS-3 now scaled by 0.83
C - removed lots of obsolete code and routines
C
C Revision 5.1  2002/08/20 13:52:09  bmy
C GEOS-CHEM version 5.01
C - coding changes to facilitate new TPCORE for fvDAS
C - added Mat Evans' plane flight diagnostic
C - consolidated all file unit numbers into "file_mod.f"
C - updated Rn-Pb-Be simulation w/ Hongyu Liu's most recent code
C - changes to facilitate benchmark runs w/ GEOS-3 winds
C - deleted a lot of obsolete code
C
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
!       Contact: Bob Yantosca, Harvard University (bmy@io.harvard.edu)
!                                                                     
!******************************************************************************
!
!  See the GEOS-CHEM User's Guide:
!
!    www-as.harvard.edu/chemistry/trop/geos/documentation/geos_chem_index.html
!
!  for the most up-to-date GEOS-CHEM documentation on the following topics:
!
!    - installation, compilation, and execution
!    - coding practice and style
!    - input files and met field data files
!    - horizontal and vertical resolution
!    - modification history
!
!******************************************************************************
!
!  Types of GEOS-CHEM Simulations:
!  ===========================================================================
!  (1 ) Rn-Pb-Be                        (6 ) Tagged Ox
!  (2 ) CH3I                            (7 ) Tagged CO
!  (3 ) NOx-Ox-HC (w/ or w/o aerosols)  (8 ) C2H6   
!  (4 ) HCN-CH3CN                       (9 ) CH4
!  (5 ) CO w/ parameterized OH          (10) Offline sulfate aerosols
!******************************************************************************
!
      ! References to F90 modules 
      USE A3_READ_MOD
      USE A6_READ_MOD
      USE BPCH2_MOD,         ONLY : BPCH2, GET_MODELNAME, 
     &                              OPEN_BPCH2_FOR_WRITE
      USE CHEMISTRY_MOD,     ONLY : DO_CHEMISTRY
      USE CONVECTION_MOD,    ONLY : DO_CONVECTION
      USE COMODE_MOD,        ONLY : INIT_COMODE
      USE DIAG_MOD,          ONLY : DIAGCHLORO
      USE DIAG51_MOD,        ONLY : DIAG51
      USE DAO_MOD
      USE DRYDEP_MOD,        ONLY : DO_DRYDEP
      USE EMISSIONS_MOD,     ONLY : DO_EMISSIONS
      USE ERROR_MOD
      USE FILE_MOD       
      USE GLOBAL_CH4_MOD,    ONLY : INIT_GLOBAL_CH4, CH4_AVGTP
      USE I6_READ_MOD
      USE PHIS_READ_MOD
      USE PLANEFLIGHT_MOD,   ONLY : SETUP_PLANEFLIGHT
      USE PRESSURE_MOD
      USE TIME_MOD
      USE RESTART_MOD,       ONLY : READ_RESTART_FILE, MAKE_RESTART_FILE
      USE UPBDFLX_MOD,       ONLY : DO_UPBDFLX, UPBDFLX_NOY
      USE UVALBEDO_MOD,      ONLY : READ_UVALBEDO
      USE TRACERID_MOD 
      USE TRANSPORT_MOD,     ONLY : DO_TRANSPORT
      USE TIME_MOD 
      USE WETSCAV_MOD,       ONLY : INIT_WETSCAV, DO_WETDEP

      ! Force all variables to be declared explicitly
      IMPLICIT NONE
      
      ! Header files 
#     include "CMN_SIZE"        ! Size parameters
#     include "CMN"             ! STT
#     include "CMN_DIAG"        ! Diagnostic switches
#     include "CMN_GCTM"        ! Physical constants
#     include "CMN_O3"          ! TCOBOX, FMOL, SAVEOH
#     include "CMN_SETUP"       ! GEOS-CTM logical flags and file I/O units
#     include "CMN_TIMES"       ! areas time series archival (diag49)
#     include "comode.h"        ! CSAVE, IDEMS

      ! Local variables
      LOGICAL :: FIRST = .TRUE.
      INTEGER :: I, IOS, J, K, L, N
      INTEGER :: JDAY, NSEASON_last, MONTH_last, NYMD_PHIS
      INTEGER :: N_DYN_STEPS, NSECb, N_STEP,  DATE(2)

      ! For setup.f (we will get rid of this soon!)
      INTEGER :: THISNYMDb,  THISNYMDe
      INTEGER :: THISNHMSb,  THISNHMSe
      INTEGER :: NDT,        NTDT,  IGD, J1, NDIAGTIME
      INTEGER :: IORD,       JORD,  KORD, KS
      REAL*8  :: PSTP, UMax, ALPHA_d, ALPHA_n

      !=================================================================
      ! GEOS-CHEM starts here!                                            
      !=================================================================

      ! Display current grid resolution and data set type
      CALL DISPLAY_GRID_AND_MODEL

      ! Get the YYYYMMDD value for the PHIS met field 
      NYMD_PHIS = GET_NYMD_PHIS()

      !=================================================================      
      !           ***** R E A D   I N P U T   F I L E *****
      !=================================================================
      CALL SETUP( THISNYMDb, THISNYMDe,  THISNHMSb, THISNHMSe,  
     &            NDT,       NTDT,       NDIAGTIME, ALPHA_d, 
     &            ALPHA_n,   IORD,       JORD,      KORD,    
     &            J1,        Umax,       IGD,       KS,  PSTP )

      NSEASON_last = 0           ! Saves the previous value of SEASON
      MONTH_last   = 0           ! Saves the previous value of MONTH

      !=================================================================
      !            ***** I N I T I A L I Z A T I O N *****
      !=================================================================

      ! Initialize quantities from "time_mod.f"
      CALL SET_BEGIN_TIME( THISNYMDb, THISNHMSb )
      CALL SET_END_TIME(   THISNYMDe, THISNHMSe )
      CALL SET_CURRENT_TIME
      CALL SET_DIAGb( GET_TAU() )
      
      ! Read settings from "input.ctm" and "inptr.ctm" 
      CALL INPUT
      CALL IS_LAST_DAY_GOOD

      ! Initialize met field arrays from "dao_mod.f"
      CALL INIT_DAO

      ! Initialize diagnostic arrays and countesr
      CALL INITIALIZE( 1 )
      CALL INITIALIZE( 2 )
      CALL INITIALIZE( 3 )

      ! Initialize the new hybrid pressure module.  Define Ap and Bp.
      CALL INIT_PRESSURE

      ! Read annual mean tropopause    
      CALL READ_TROPOPAUSE
 
      ! Initialize allocatable SMVGEAR arrays
      IF ( NSRCX == 3 ) THEN 
         IF ( LEMIS .or. LCHEM ) CALL INIT_COMODE
      ENDIF

      ! Allocate arrays from "global_ch4_mod.f" for CH4 run 
      IF ( NSRCX == 9 ) CALL INIT_GLOBAL_CH4

      ! Define flag for debug output
      LPRT = ( ND70 > 0 )
      IF ( LPRT ) CALL DEBUG_MSG( '### MAIN: a INPUT' )

      ! LSRCE is a synonym for LEMIS
      LSRCE = LEMIS

      !### Kludge: turn off wetdep for runs that don't require them
      !### Find a better solution to this later (bmy, 11/8/02)
      IF ( NSRCX /= 1 .AND. NSRCX /= 3 .AND. NSRCX /= 10 ) THEN
         WRITE( 6, '(a)' ) 'Setting LWETD = F!'
         LWETD = .FALSE.
      ENDIF

      ! Open station timeseries and binary punch files for output
      IF ( ND48 > 0 ) OPEN( IU_TS, FILE='ctm.ts', STATUS='UNKNOWN' )
      CALL OPEN_BPCH2_FOR_WRITE( IU_BPCH, 'ctm.bpch', XLABEL )

      !=================================================================
      !   ***** U N Z I P   M E T   F I E L D S  @ start of run *****
      !=================================================================
      IF ( LUNZIP ) THEN

         ! Remove any leftover A-3, A-6, I-6, and PHIS files in temp dir
         CALL UNZIP_A3_FIELDS(  'remove all' )
         CALL UNZIP_A6_FIELDS(  'remove all' )
         CALL UNZIP_I6_FIELDS(  'remove all' )
         CALL UNZIP_PHIS_FIELD( 'remove all' )

         ! Unzip A-3, A-6, I-6 files for START of run
         CALL UNZIP_A3_FIELDS(  'unzip foreground', GET_NYMDb() )
         CALL UNZIP_A6_FIELDS(  'unzip foreground', GET_NYMDb() )
         CALL UNZIP_I6_FIELDS(  'unzip foreground', GET_NYMDb() )
         
         ! Only unzip PHIS file for fullchem run
         IF ( NSRCX == 3 ) THEN
            CALL UNZIP_PHIS_FIELD( 'unzip foreground', NYMD_PHIS )
         ENDIF
      ENDIF
      
      !=================================================================
      !    ***** R E A D   M E T   F I E L D S  @ start of run *****
      !=================================================================

      ! Only read PHIS for fullchem runs
      IF ( NSRCX == 3 ) THEN
         CALL OPEN_PHIS_FIELD( NYMD_PHIS, 000000 )
         CALL GET_PHIS_FIELD(  NYMD_PHIS, 000000 ) 

         ! Remove PHIS field from temp dir
         IF ( LUNZIP ) THEN
            CALL UNZIP_PHIS_FIELD( 'remove date', NYMD_PHIS )
         ENDIF
      ENDIF

      ! Open and read A-3 fields
      DATE = GET_FIRST_A3_TIME()
      CALL OPEN_A3_FIELDS( DATE(1), DATE(2) )
      CALL GET_A3_FIELDS(  DATE(1), DATE(2) ) 

      ! Open & read A-6 fields
      DATE = GET_FIRST_A6_TIME()
      CALL OPEN_A6_FIELDS( DATE(1), DATE(2) ) 
      CALL GET_A6_FIELDS(  DATE(1), DATE(2) )      

      ! Open & read I-6 fields
      DATE = (/ GET_NYMD(), GET_NHMS() /)
      CALL OPEN_I6_FIELDS(  DATE(1), DATE(2) )
      CALL GET_I6_FIELDS_1( DATE(1), DATE(2) )
      
      ! Compute avg surface pressure near polar caps
      CALL AVGPOLE( PS1 )

      ! Call AIRQNT to compute air mass quantities from PS1
      CALL SET_FLOATING_PRESSURE( PS1 )
      CALL AIRQNT

      ! Compute lightning NOx emissions [molec/box/6h]
      IF ( LLIGHTNOX .and. NSRCX == 3 ) THEN
         CALL LIGHTNING( T, CLDTOPS )
      ENDIF

      ! Read land types and fractions from "vegtype.global"
      CALL RDLAND   

      ! Initialize XTRA2 array (PBL heights, # of layers)
      CALL TURBDAY( .FALSE.,  XTRA2,  NTRACE,          
     &               STT(:,:,:,1:NTRACE), TCVV(1:NTRACE) )

      IF ( LPRT ) CALL DEBUG_MSG( '### MAIN: a TURBDAY:1' )

      !=================================================================
      !    ***** R E A D   I N I T I A L   C O N D I T I O N S *****
      !=================================================================
      CALL READ_RESTART_FILE( GET_NYMDb(), GET_NHMSb() )

      ! Save initial tracer masses to disk for benchmark runs
      IF ( LSTDRUN ) CALL STDRUN( LBEGIN=.TRUE. )

      !=================================================================
      !     ***** 6 - H O U R   T I M E S T E P   L O O P  *****
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

      !=================================================================
      !     ***** D Y N A M I C   T I M E S T E P   L O O P *****
      !=================================================================
      DO N_STEP = 1, N_DYN_STEPS
    
         ! Compute & print time quantities at start of dyn step
         CALL SET_CURRENT_TIME
         CALL PRINT_CURRENT_TIME

         !==============================================================
         !   ***** W R I T E   D I A G N O S T I C   F I L E S *****
         !
         ! NOTE: If this is the first timestep, don't write diagnostics
         !       to the punch file since we haven't done anything yet.
         !==============================================================
         IF ( ITS_TIME_FOR_BPCH() .and. ( .not. FIRST ) ) THEN

            ! Set time at end of diagnostic timestep
            CALL SET_DIAGe( GET_TAU() )

            ! Write bpch and ND48 files
            CALL DIAG3          
            CALL DIAG5          

            ! Flush file units
            CALL CTM_FLUSH

            !===========================================================
            !      ***** W R I T E   R E S T A R T   F I L E *****
            !===========================================================
            IF ( LSVGLB ) THEN
               CALL MAKE_RESTART_FILE( GET_NYMD(), GET_NHMS() )
                  
               !### Debug
               IF ( LPRT ) THEN
                  CALL DEBUG_MSG( '### MAIN: a MAKE_RESTART_FILE' )
               ENDIF
            ENDIF

            ! Set time at beginning of next diagnostic timestep
            CALL SET_DIAGb( GET_TAU() )

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
            ! fully unzipped before proceeding w/ the run
            IF ( LWAIT ) THEN
               CALL UNZIP_A3_FIELDS( 'unzip foreground', DATE(1) )
               CALL UNZIP_A6_FIELDS( 'unzip foreground', DATE(1) )
               CALL UNZIP_I6_FIELDS( 'unzip foreground', DATE(1) )
            ELSE
               CALL UNZIP_A3_FIELDS( 'unzip background', DATE(1) )
               CALL UNZIP_A6_FIELDS( 'unzip background', DATE(1) )
               CALL UNZIP_I6_FIELDS( 'unzip background', DATE(1) )
            ENDIF
         ENDIF     

         !===============================================================
         !        ***** R E M O V E   M E T   F I E L D S *****  
         !===============================================================
         IF ( LUNZIP .and. ITS_TIME_FOR_DEL() ) THEN

            ! Get the current date
            DATE(1) = GET_NYMD()

            ! Remove A-3, A-6, and I-6 files only for the current date
            CALL UNZIP_A3_FIELDS( 'remove date', DATE(1) )
            CALL UNZIP_A6_FIELDS( 'remove date', DATE(1) )
            CALL UNZIP_I6_FIELDS( 'remove date', DATE(1) )
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
            IF ( LLIGHTNOX .and. NSRCX == 3 ) THEN 
               CALL LIGHTNING( T, CLDTOPS )
            ENDIF 
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
         IF ( GET_MONTH() /= MONTH_last .and. LCHEM ) THEN 
            CALL READ_UVALBEDO( GET_MONTH() )
            MONTH_last = GET_MONTH()
         ENDIF

         ! Fossil fuel emissions (SMVGEAR)
         IF ( GET_SEASON() /= NSEASON_last .and. 
     &        NSRCX == 3 .and. LEMIS ) THEN
            CALL ANTHROEMS( GET_SEASON() )
            NSEASON_last = GET_SEASON()
         ENDIF

         !==============================================================
         !              ***** D A I L Y   D A T A *****
         !
         ! RDLAI  returns today's leaf-area index
         ! RDSOIL returns today's soil type information
         !==============================================================
         IF ( FIRST .or. GET_NHMS() == 000000 ) THEN
            SELECT CASE ( NSRCX )
               CASE ( 1, 5, 6, 7, 10 ) 
                  CALL RDLAI( GET_DAY_OF_YEAR(), GET_MONTH() )

               CASE ( 3 )
                  CALL RDLAI( GET_DAY_OF_YEAR(), GET_MONTH() )
                  CALL RDSOIL

               CASE DEFAULT
                  ! Nothing

            END SELECT

            ! Reset first-time flag
            FIRST = .FALSE.

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
         ! based on their values at NSEC and NSEC+NTDT
         CALL INTERP( NSECb, GET_ELAPSED_SEC(), NTDT )         
         
         ! If we are not doing transport, then make sure that
         ! the floating pressure is set to PSC2 (bdf, bmy, 8/22/02)
         IF ( .not. LTRAN ) CALL SET_FLOATING_PRESSURE( PSC2 )
         
         ! Compute airmass quantities at each grid box 
         CALL AIRQNT
         
         ! Compute the cosine of the solar zenith angle at each grid box
         CALL COSSZA( GET_DAY_OF_YEAR(), GET_NHMSb(), 
     &                GET_ELAPSED_SEC(), SUNCOS )
         
         ! For SMVGEAR II, we also need to compute SUNCOS at
         ! the end of this chemistry timestep (bdf, bmy, 4/1/03)
         IF ( NSRCX == 3 .and. LCHEM ) THEN
            CALL COSSZA( GET_DAY_OF_YEAR(), GET_NHMSb(), 
     &                   GET_ELAPSED_SEC()+GET_TS_CHEM()*60, SUNCOSB )
         ENDIF

         ! Compute tropopause height for ND55 diagnostic
         IF ( ND55 > 0 ) CALL TROPOPAUSE

#if   defined( GEOS_STRAT )
         ! For GEOS-STRAT, if U10M and V10M are missing, compute 
         ! the resultant wind speed at 10 meters (bmy, 6/27/00)
         CALL GET_WIND10M( GET_NYMD() )

#elif defined( GEOS_3 )
         ! 1998 GEOS-3 carries the ground temperature and not the air
         ! temperature -- thus TS will be 2-3 K too high.  As a quick fix, 
         ! copy the temperature at the first sigma level into TS. 
         ! (mje, bnd, bmy, 7/3/01)
         IF ( GET_YEAR() == 1998 ) TS(:,:) = T(:,:,1)
#endif

         ! Update dynamic timestep
         CALL SET_CT_DYN( INCREMENT=.TRUE. )

         !### Debug
         IF ( LPRT ) CALL DEBUG_MSG( '### MAIN: a INTERP, etc' )

         !==============================================================
         !   ***** U N I T   C O N V E R S I O N  ( kg -> v/v ) *****
         !==============================================================
         IF ( ITS_TIME_FOR_UNIT() ) THEN
            CALL CONVERT_UNITS( 1, NTRACE, TCVV(1:NTRACE), 
     &                          AD,     STT(:,:,:,1:NTRACE) )

            !### Debug
            IF ( LPRT ) CALL DEBUG_MSG( '### MAIN: a CONVERT_UNITS:1' )
         ENDIF

         !==============================================================
         !     ***** S T R A T O S P H E R I C   F L U X E S *****
         !==============================================================
         IF ( LUPBD ) THEN
            CALL DO_UPBDFLX( IORD, JORD, KORD )
         ENDIF

         !==============================================================
         !              ***** T R A N S P O R T *****
         !==============================================================
         IF ( ITS_TIME_FOR_DYN() ) THEN

            ! Call the appropritate version of TPCORE
            IF ( LTRAN ) CALL DO_TRANSPORT( IORD, JORD, KORD )

            ! Reset air mass quantities
            CALL AIRQNT

            ! Repartition [NOy] species after transport
            IF ( LUPBD .and. NSRCX == 3 ) CALL UPBDFLX_NOY( 2 )

            ! Initialize wet scavenging and wetdep fields after
            ! the airmass quantities are reset after transport
            IF ( LCONV .or. LWETD ) CALL INIT_WETSCAV
         ENDIF

         ! Test for convection timestep
         IF ( ITS_TIME_FOR_CONV() ) THEN

            ! Increment the convection timestep
            CALL SET_CT_CONV( INCREMENT=.TRUE. )

            !===========================================================
            !      ***** M I X E D   L A Y E R   M I X I N G *****
            !===========================================================
            CALL TURBDAY( LTURB, XTRA2, NTRACE,         
     &                    STT(:,:,:,1:NTRACE), TCVV(1:NTRACE) )

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
         ENDIF                  ! mod(NMIN,NCONV)

         !==============================================================
         !    ***** U N I T   C O N V E R S I O N  ( v/v -> kg ) *****
         !==============================================================
         IF ( ITS_TIME_FOR_UNIT() ) THEN 
            CALL CONVERT_UNITS( 2, NTRACE, TCVV(1:NTRACE), 
     &                          AD,     STT(:,:,:,1:NTRACE) )

            IF ( LPRT ) CALL DEBUG_MSG( '### MAIN: a CONVERT_UNITS:2' )
         ENDIF

         !==============================================================
         !       ***** A R C H I V E   D I A G N O S T I C S *****
         !==============================================================
         IF ( ITS_TIME_FOR_DYN() ) THEN

            ! Accumulate several diagnostic quantities
            CALL DIAG1

            ! ND41: save PBL height in 1200-1600 LT (amf)
            ! (for comparison w/ Holzworth, 1967
            IF ( ND41 > 0 ) CALL DIAG41

            ! ND49: 3-D timeseries ("tsYYYYMMDD.bpch")
            IF ( ND49 > 0 ) THEN
               IF ( GET_TAU() >= DAT1_AREA .and. 
     &              GET_TAU() <= DAT2_AREA .and. 
     &              MOD( GET_ELAPSED_MIN(), FREQ_AREA ) == 0 ) THEN
                  CALL DIAG49
               ENDIF
            ENDIF

            ! ND50: 24-h avg timeseries ("ts24h.bpch")
            IF ( ND50 > 0 ) THEN
               IF ( GET_TAU() >= DAT1_AREA  .and. 
     &              GET_TAU() <= DAT2_AREA ) THEN
                  CALL DIAG50
               ENDIF
            ENDIF
    
            ! ND51: 10-12 LT or 13-16 LT timeseries ("ts1_4pm.bpch")
            IF ( ND51 > 0 ) THEN
               IF ( GET_TAU() >= DAT1_AREA  .and. 
     &              GET_TAU() <= DAT2_AREA ) THEN
                  CALL DIAG51
               ENDIF
            ENDIF
            
            !### Debug
            IF ( LPRT ) CALL DEBUG_MSG( '### MAIN: a DIAGNOSTICS' )
         ENDIF

         ! ND48: station timeseries data ("ctm.ts")
         IF ( ND48 > 0 .and. ITS_TIME_FOR_DIAG() ) THEN
            KDA48 = KDA48 + 1
            CALL DIAG48
         ENDIF

         ! Test for emission timestep
         IF ( ITS_TIME_FOR_EMIS() ) THEN
               
            ! Increment emission counter
            CALL SET_CT_EMIS( INCREMENT=.TRUE. )

            !========================================================
            !         ***** D R Y   D E P O S I T I O N *****
            !========================================================
            IF ( LDRYD ) THEN
               CALL DO_DRYDEP
            ENDIF

            !===========================================================
            !             ***** E M I S S I O N S *****
            !===========================================================
            IF ( LEMIS ) THEN 
               CALL DO_EMISSIONS
            ENDIF
         ENDIF    

         !===========================================================
         !               ***** C H E M I S T R Y *****
         !===========================================================    

         ! Setup fields for CO_OH parameterization
         CALL SETUP_CO_OH_CHEM

         ! Also need to compute avg P, T for CH4 chemistry (bmy, 1/16/01)
         IF ( NSRCX == 9 ) CALL CH4_AVGTP

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
         IF ( ITS_TIME_FOR_DYN() ) THEN
            IF ( LWETD ) CALL DO_WETDEP
         ENDIF
            
         !==============================================================
         !  ***** E N D   O F   D Y N A M I C   T I M E S T E P *****
         !==============================================================

         ! Check for NaN, Negatives, Infinities in STT once per hour
         IF ( ITS_TIME_FOR_DIAG() ) THEN
            CALL CHECK_STT( 'End of Dynamic Loop' )
         ENDIF

         ! Increment elapsed time
         CALL SET_ELAPSED_MIN
      ENDDO

      !=================================================================
      !            ***** C O P Y   I - 6   F I E L D S *****
      !
      !        The I-6 fields at the end of this timestep become
      !        the fields at the beginning of the next timestep
      !=================================================================
      PS1   = PS2 

#if   defined( GEOS_1 ) || defined( GEOS_STRAT ) || defined( GEOS_3 )
      ALBD1 = ALBD2    
      UWND1 = UWND2
      VWND1 = VWND2
      TMPU1 = TMPU2
      SPHU1 = SPHU2
#endif

      ENDDO

      !=================================================================
      !         ***** C L E A N U P   A N D   Q U I T *****
      !=================================================================
 9999 CONTINUE

      ! Remove all files from temporary directory 
      IF ( LUNZIP ) THEN
         CALL UNZIP_A3_FIELDS(  'remove all' )
         CALL UNZIP_A6_FIELDS(  'remove all' )
         CALL UNZIP_I6_FIELDS(  'remove all' )
         CALL UNZIP_PHIS_FIELD( 'remove all' )
      ENDIF

      ! Print the mass-weighted mean OH concentration (if applicable)
      CALL PRINT_MEAN_OH

      ! For model benchmarking, save final masses of 
      ! Rn,Pb,Be or Ox to a binary punch file 
      IF ( LSTDRUN ) CALL STDRUN( LBEGIN=.FALSE. )

      ! Close all files
      CALL CLOSE_FILES
      IF ( LPRT ) CALL DEBUG_MSG( '### MAIN: a CLOSE_FILES' )

      ! Deallocate dynamic module arrays
      CALL CLEANUP
      IF ( LPRT ) CALL DEBUG_MSG( '### MAIN: a CLEANUP' )

      ! Echo info
      WRITE ( 6, 3000 ) 
 3000 FORMAT
     &   ( /, /, '**************   E N D   O F   G E O S   C H E M   ',
     &    '***************')
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
!  (1 ) DISPLAY_GRID_AND_MODEL : Displays resolution & data set
!  (2 ) GET_NYMD_PHIS          : Gets YYYYMMDD for the PHIS data field
!  (3 ) DISPLAY_SIGMA_LAT_LON  : Displays sigma, lat, and lon information
!  (4 ) GET_WIND10M            : Wrapper for MAKE_WIND10M (from "dao_mod.f")
!  (5 ) CTM_FLUSH              : Flushes diagnostic files to disk
!  (6 ) SETUP_CO_OH_CHEM       : Calls setup routines for CO-OH param run
!  (7 ) CH3CCL3_LIFETIME       : Computes methyl chloroform lifetime
!  (8 ) MET_FIELD_DEBUG        : Prints min and max of met fields for debug
!  (9 ) STDRUN                 : Saves initial/final masses for benchmarks
!******************************************************************************
!
      CONTAINS

!-----------------------------------------------------------------------------

      SUBROUTINE DISPLAY_GRID_AND_MODEL

      !=================================================================
      ! Internal Subroutine DISPLAY_GRID_AND_MODEL displays the 
      ! appropriate messages for the given model grid (2 x 2.5, 4 x 5)
      ! and machine type (bmy, 12/2/03)
      !=================================================================

      !-----------------------
      ! Print resolution info
      !-----------------------
#if   defined( GRID4x5  )
      WRITE( 6, '(a)' )                   
     &    REPEAT( '*', 13 )                                    //
     &    '   S T A R T I N G   4 x 5   G E O S--C H E M   '   //
     &    REPEAT( '*', 13 )

#elif defined( GRID2x25 )
      WRITE( 6, '(a)' ) 
     &    REPEAT( '*', 13 )                                    // 
     &    '   S T A R T I N G   2 x 2.5   G E O S--C H E M   ' //
     &    REPEAT( '*', 13 )

#elif defined( GRID1x1 )
      WRITE( 6, '(a)' ) 
     &    REPEAT( '*', 13 )                                    // 
     &    '   S T A R T I N G   1 x 1   G E O S--C H E M   '   //
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
#elif defined( LINUX_IFC )
      WRITE( 6, '(a)' ) 'Created w/ LINUX/IFC (32-bit) compiler'
#elif defined( LINUX_EFC )
      WRITE( 6, '(a)' ) 'Created w/ LINUX/EFC (64-bit) compiler'
#elif defined( SGI_MIPS  )
      WRITE( 6, '(a)' ) 'Created w/ SGI MIPSpro compiler'
#elif defined( SPARC     )
      WRITE( 6, '(a)' ) 'Created w/ Sun/SPARC compiler'
#endif

      !-----------------------
      ! Print met field info
      !-----------------------
#if   defined( GEOS_1     )
      WRITE( 6, '(a)' ) 'Using GEOS-1 met fields' 
#elif defined( GEOS_STRAT )
      WRITE( 6, '(a)' ) 'Using GEOS-STRAT met fields'
#elif defined( GEOS_3     )
      WRITE( 6, '(a)' ) 'Using GEOS-3 met fields'
#elif defined( GEOS_4     )
      WRITE( 6, '(a)' ) 'Using GEOS-4/fvDAS met fields'
#endif

      ! Return to MAIN program
      END SUBROUTINE DISPLAY_GRID_AND_MODEL

!-----------------------------------------------------------------------------

      FUNCTION GET_NYMD_PHIS() RESULT( NYMD_PHIS )

      !=================================================================
      ! Internal Function GET_NYMD_PHIS selects the YYMMDD at which 
      ! there is PHIS data for the given model (GEOS-1, GEOS-STRAT, 
      ! GEOS-3, or fvDAS).  This date is Y2K compliant.
      !=================================================================

      ! Return value 
      INTEGER :: NYMD_PHIS

#if   defined( GEOS_1 )
      NYMD_PHIS = 19940101

#elif defined( GEOS_STRAT )
      NYMD_PHIS = 19960101

#elif defined( GEOS_3 )
      NYMD_PHIS = 20000101

#elif defined( GEOS_4 )
      NYMD_PHIS = 20030101

#endif

      ! Return to MAIN program
      END FUNCTION GET_NYMD_PHIS
      
!-----------------------------------------------------------------------------

      SUBROUTINE IS_LAST_DAY_GOOD

      !=================================================================
      ! Internal function IS_LAST_DAY_GOOD tests to see if there is
      ! output scheduled on the last day of the run.  (bmy, 9/25/03)
      !=================================================================

      ! References to F90 modules
      USE JULDAY_MOD, ONLY : JULDAY

      ! Local variables
      INTEGER :: NYMDe, Y, M, D, LASTDAY
      REAL*8  :: JD, JD0

      !=================================================================
      ! IS_LAST_DAY_GOOD begins here!
      !=================================================================

      ! Astronomical Julian Day corresponding to NYMDe
      NYMDe = GET_NYMDe()
      CALL YMD_EXTRACT( NYMDe, Y, M, D )
      JD = JULDAY( Y, M, DBLE( D ) )

      ! Astronomical Julian Day corresponding to the 1st of the year
      JD0 = JULDAY( Y, 1, 0d0 )

      ! LASTDAY is the day of year corresponding to NYMDe      
      LASTDAY = JD - JD0 

      ! Skip past the element of NJDAY for Feb 29, if necessary
      IF ( .not. ITS_A_LEAPYEAR( Y ) .and. LASTDAY > 59 ) THEN
         LASTDAY = LASTDAY + 1
      ENDIF

      ! Stop w/ error if THIS_NJDAY = 0 
      IF ( NJDAY(LASTDAY) == 0 ) THEN
         CALL ERROR_STOP( 'No output scheduled on last day of run!',
     &                    'IS_LAST_DAY_GOOD (main.f)' )
      ENDIF
     
      ! Return to calling program
      END SUBROUTINE IS_LAST_DAY_GOOD
      
!-----------------------------------------------------------------------------

      FUNCTION ITS_TIME_FOR_BPCH() RESULT( DO_BPCH )

      !=================================================================
      ! Internal function ITS_TIME_FOR_BPCH returns TRUE if it is time
      ! to write to the binary punch file and FALSE otherwise.
      !=================================================================

      ! Local variables
      INTEGER :: TODAY, THIS_NJDAY

      ! Function value
      LOGICAL :: DO_BPCH

      !=================================================================
      ! ITS_TIME_FOR_BPCH begins here!
      !=================================================================

      ! Current day of year
      TODAY = GET_DAY_OF_YEAR()

      ! Look up appropriate value of NJDAY array.  We may need to add a
      ! day to skip past the Feb 29 element of NJDAY for non-leap-years.
      IF ( .not. ITS_A_LEAPYEAR()  .and. TODAY > 59 ) THEN
         THIS_NJDAY = NJDAY( TODAY + 1 ) 
      ELSE
         THIS_NJDAY = NJDAY( TODAY )
      ENDIF

      ! Test if this is the day & time to write to the BPCH file!
      IF ( ( THIS_NJDAY > 0 ) .and. GET_NHMS() == NDIAGTIME ) THEN
         DO_BPCH = .TRUE.
      ELSE
         DO_BPCH = .FALSE.
      ENDIF

      ! Return to calling program
      END FUNCTION ITS_TIME_FOR_BPCH

!------------------------------------------------------------------------------

      SUBROUTINE GET_WIND10M( NYMD )

      !=================================================================
      ! Internal Subroutine GET_WIND10M is a wrapper for routine
      ! MAKE_WIND10M (from "dao_mod.f").  
      !=================================================================

      ! Arguments
      INTEGER, INTENT(IN) :: NYMD 

      !=================================================================
      ! GEOS-STRAT winds prior to 13 Feb 1997 do not contain U10M and
      ! V10M fields.  Call MAKE_WIND10M to compute these fields from 
      ! the existing values of UWND, VWND, BXHEIGHT, and Z0.
      !
      ! Also set logical flag USE_WIND_10M to let SFCWINDSQR know that 
      ! the U10M and V10M fields are missing and had to be constructed 
      ! by MAKE_WIND10M. (bmy, 6/27/00)
      !=================================================================
#if   defined( GEOS_STRAT )
      IF ( NYMD < 19970213 ) THEN
         CALL MAKE_WIND10M( UWND=UWND(:,:,1),         VWND=VWND(:,:,1), 
     &                      BXHEIGHT=BXHEIGHT(:,:,1), Z0=Z0 )
         USE_WIND_10M = .TRUE.
      ELSE
         USE_WIND_10M = .FALSE.
      ENDIF
#endif
      
      ! Return to MAIN program
      END SUBROUTINE GET_WIND10M

!-----------------------------------------------------------------------------

      SUBROUTINE SETUP_CO_OH_CHEM

      !=================================================================
      ! Internal Subroutine CO_OH_CHEM_SETUP calls subroutines which
      ! set up variables for the CO-OH parameterization simulation
      ! (when LGEOSCO is defined and when NSRCX == 5).
      !
      ! NOTE: These routines need to be called every dynamic timestep.
      !=================================================================
#if   defined( LGEOSCO )

      ! Get optical depth for the proper data set
#if   defined( GEOS_1 ) || defined( GEOS_STRAT )

      ! GEOS-1/GEOS-STRAT: Compute OPTD from T, CLMO, CLRO, DELP
      CALL OPTDEPTH( LM, CLMOSW, CLROSW, DELP, T, OPTD )

#elif defined( GEOS_2 ) || defined( GEOS_3 )

      ! GEOS-2/GEOS-3: Copy OPTDEP to OPTD, also save diagnostics
      CALL OPTDEPTH( LM, OPTDEP, OPTD )

#endif

      ! Compute 24h average reflectivities
      CALL AVGREFL( FIRSTCHEM, OPTD, NTDT, NHMSb, NSEC, XMID, NMIN )

      ! Compute 24h average temp & pressure
      CALL AVGTP( NTDT, NMIN )

      ! Compute 24h average water vapor
      CALL AVGAVGW( NTDT, NMIN, NCHEM )

#endif

      ! Return to MAIN program
      END SUBROUTINE SETUP_CO_OH_CHEM

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
      CALL FLUSH( IU_TS      )  
      CALL FLUSH( IU_BPCH    )  
      CALL FLUSH( IU_SMV2LOG )  
      CALL FLUSH( IU_DEBUG   ) 

      ! Return to calling program
      END SUBROUTINE CTM_FLUSH

!-----------------------------------------------------------------------------

      SUBROUTINE PRINT_MEAN_OH

      !=================================================================
      ! Internal Subroutine PRINT_MEAN_OH prints the mass-weighted OH
      ! concentration for selected simulation types. (bmy, 10/21/03)
      !=================================================================

      ! References to F90 modules
      USE DIAG_MOD, ONLY : DIAGCHLORO

      ! Local variables 
      LOGICAL :: DO_PRINT
      REAL*8  :: SUM_OHMASS, SUM_MASS, LIFETIME, OHCONC
     
      !=================================================================
      ! PRINT_MEAN_OH begins here!
      !=================================================================

      ! First figure out if this type of simulation has OH
      IF ( ( ND23 > 0   .and. LCHEM                      )  .and.
     &     ( NSRCX == 3 .or.  NSRCX == 5 .or. NSRCX == 9 ) ) THEN
         DO_PRINT = .TRUE.
      ELSE
         DO_PRINT = .FALSE.
      ENDIF
      
      ! Print out Mean OH if it's the appropriate simulation type
      IF ( DO_PRINT ) THEN
 
         ! Total Mass-weighted OH [molec OH/cm3] * [molec air]
         SUM_OHMASS = SUM( DIAGCHLORO(:,:,:,2) )

         ! Atmospheric air mass [molec air]
         SUM_MASS   = SUM( DIAGCHLORO(:,:,:,3) ) 
         
         ! Avoid divide-by-zero errors 
         IF ( SUM_MASS > 0d0 ) THEN 
            
            ! Divide OH by [molec air] and report as [1e5 molec/cm3]
            OHCONC = ( SUM_OHMASS / SUM_MASS ) / 1d5
         
            ! Write value to log file
            WRITE( 6, '(a)' ) 
            WRITE( 6, '(a)' ) REPEAT( '=', 79 ) 
            WRITE( 6, *     ) 'ND23: Mass-Weighted OH Concentration'
            WRITE( 6, *     ) 'Mean OH = ', OHCONC, ' [1e5 molec/cm3]' 
            WRITE( 6, '(a)' ) REPEAT( '=', 79 ) 

         ELSE

            ! Write error msg if SUM_MASS is zero
            WRITE( 6, '(a)' ) 
            WRITE( 6, '(a)' ) REPEAT( '=', 79 ) 
            WRITE( 6, '(a)' ) 'Could not print mass-weighted OH!'
            WRITE( 6, '(a)' ) 'Atmospheric air mass is zero!'
            WRITE( 6, '(a)' ) REPEAT( '=', 79 ) 
            
         ENDIF
      ENDIF

      ! Return to MAIN program
      END SUBROUTINE PRINT_MEAN_OH

!-----------------------------------------------------------------------------

      SUBROUTINE MET_FIELD_DEBUG

      !=================================================================
      ! Internal subroutine MET_FIELD_DEBUG prints out the maximum
      ! and minimum, and sum of DAO met fields for debugging 
      !=================================================================

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
      IF ( ALLOCATED( CLMOSW   ) ) PRINT*, 'CLMOSW  : ', CLMOSW(L,I,J) 
      IF ( ALLOCATED( CLROSW   ) ) PRINT*, 'CLROSW  : ', CLROSW(L,I,J) 
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
      IF ( ALLOCATED( PHIS     ) ) PRINT*, 'PHIS    : ', PHIS(I,J) 
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
      IF ( ALLOCATED( WIND_10M ) ) PRINT*, 'WIND_10M: ', WIND_10M(I,J) 
      IF ( ALLOCATED( Z0       ) ) PRINT*, 'Z0      : ', Z0(I,J) 
      IF ( ALLOCATED( ZMEU     ) ) PRINT*, 'ZMEU    : ', ZMEU(I,J,L) 
      IF ( ALLOCATED( ZMMD     ) ) PRINT*, 'ZMMD    : ', ZMMD(I,J,L) 
      IF ( ALLOCATED( ZMMU     ) ) PRINT*, 'ZMMU    : ', ZMMU(I,J,L) 

      ! Flush the output buffer
      CALL FLUSH( 6 )

      ! Return to MAIN program
      END SUBROUTINE MET_FIELD_DEBUG

!-----------------------------------------------------------------------------

      SUBROUTINE STDRUN( LBEGIN )

      !=================================================================
      ! Internal subroutine STDRUN dumps the mass of either Ox [kg] 
      ! or 222Rn, 210Pb, and 7Be [kg] at the start & end of each run, 
      ! for model benchmarking. (bmy, 8/12/02)
      !
      ! Arguments as Input:
      ! ----------------------------------------------------------------
      ! (1 ) LBEGIN (LOGICAL) : TRUE  denotes beginning of the run;
      !                         FALSE denotes the end of the run
      !
      ! NOTES:
      ! (1 ) Changed name from STDRUN_Ox to STDRUN, since we now can
      !      also save out Rn/Pb/Be for NSRCX==1.  Also deleted
      !      obsolete code from 6/02.  Added LBEGIN as an argument to
      !      determine if this is the start or end of the run.
      !      (bmy, 8/12/02)
      !=================================================================

      ! Arguments
      LOGICAL, INTENT(IN) :: LBEGIN 

      ! Local variables
      INTEGER             :: N
      INTEGER, PARAMETER  :: IFIRST=1, JFIRST=1, LFIRST=1
      INTEGER, PARAMETER  :: HALFPOLAR=1, CENTER180=1
      REAL*4              :: ARRAY(IIPAR,JJPAR,LLPAR)
      REAL*4              :: LONRES, LATRES
      CHARACTER(LEN=20)   :: MODELNAME 
      CHARACTER(LEN=40)   :: CATEGORY, RESERVED, UNIT
      CHARACTER(LEN=80)   :: TITLE
      CHARACTER(LEN=255)  :: FILENAME

      !=================================================================
      ! STDRUN begins here!
      !=================================================================

      ! Return if we are not doing either a radon or fullchem stdrun
      IF ( NSRCX /= 1 .and. NSRCX /= 3 ) RETURN

      ! Define variables for binary punch file
      MODELNAME = GET_MODELNAME()
      CATEGORY  = 'TCMASS-$'
      UNIT      = 'kg'
      RESERVED  = ''      
      LONRES    = DISIZE
      LATRES    = DJSIZE
    
      !=================================================================
      ! Save the mass of 222Rn, 210Pb, 7Be to a file
      !=================================================================
      IF ( NSRCX == 1 ) THEN

         ! Define filename for beginning or end
         IF ( LBEGIN ) THEN
            TITLE    = 'GEOS-CHEM Rn-Pb-Be benchmark: initial mass'
            FILENAME = 'Rn-Pb-Be.mass.initial'
         ELSE
            TITLE    = 'GEOS-CHEM Rn-Pb-Be benchmark: final mass'
            FILENAME = 'Rn-Pb-Be.mass.final'
         ENDIF
            
         ! Open binary punch file for writing
         CALL OPEN_BPCH2_FOR_WRITE( IU_FILE, FILENAME, TITLE )

         ! Loop over tracers
         DO N = 1, NTRACE

            ! Save Rn, Pb, Be as REAL*4
            ARRAY(:,:,:) = STT(:,:,:,N)

            ! Write Rn, Pb, Be to binary punch file
            CALL BPCH2( IU_FILE,   MODELNAME, LONRES,    LATRES,
     &                  HALFPOLAR, CENTER180, CATEGORY,  N+TRCOFFSET,    
     &                  UNIT,      GET_TAU(), GET_TAU(), RESERVED,   
     &                  IIPAR,     JJPAR,     LLPAR,     IFIRST,     
     &                  JFIRST,    LFIRST,    ARRAY(:,:,:) )

         ENDDO

      !=================================================================
      ! Save the mass of Ox to a file
      !=================================================================
      ELSE IF ( NSRCX == 3 .and. IDTOX > 0 ) THEN

         ! Define filename for beginning or end
         IF ( LBEGIN ) THEN
            TITLE    = 'GEOS-CHEM fullchem benchmark: Ox initial mass'
            FILENAME = 'Ox.mass.initial'
         ELSE
            TITLE    = 'GEOS-CHEM fullchem benchmark: Ox final mass'
            FILENAME = 'Ox.mass.final'
         ENDIF

         ! Open binary punch file for writing
         CALL OPEN_BPCH2_FOR_WRITE( IU_FILE, FILENAME, TITLE )
        
         ! Save Ox as REAL*4
         ARRAY(:,:,:) = STT(:,:,:,IDTOX)

         ! Write Ox to binary punch file
         CALL BPCH2( IU_FILE,   MODELNAME, LONRES,    LATRES,
     &               HALFPOLAR, CENTER180, CATEGORY,  IDTOX,    
     &               UNIT,      GET_TAU(), GET_TAU(), RESERVED,   
     &               IIPAR,     JJPAR,     LLPAR,     IFIRST,     
     &               JFIRST,    LFIRST,    ARRAY(:,:,:) )
               
      ENDIF

      ! Close file
      CLOSE( IU_FILE )

      !### Debug
      IF ( LPRT ) CALL DEBUG_MSG( '### MAIN: a STDRUN' )

      ! Return to MAIN program
      END SUBROUTINE STDRUN

!-----------------------------------------------------------------------------

      END PROGRAM GEOS_CHEM

