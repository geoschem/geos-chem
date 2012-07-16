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
      !-----------------------------------------------------------------
      ! Basic GEOS-Chem modules
      !-----------------------------------------------------------------
      USE CMN_SIZE_MOD          ! Size parameters
      USE CMN_GCTM_MOD          ! Physical constants
      USE ERROR_MOD             ! For error checking
      USE FILE_MOD              ! For file I/O
      USE INPUT_MOD             ! For reading settings from "input.geos"
      USE LOGICAL_MOD           ! Logical flags to toggle G-C options
      USE MAPPING_MOD           ! For regridding MODIS LAI
      USE OLSON_LANDMAP_MOD     ! Computes IREG, ILAND, IUSE from Olson map
      USE PRESSURE_MOD          ! For computing pressure at grid boxes
      USE RESTART_MOD           ! For restart file I/O
      USE TIME_MOD              ! For computing date & time 
      USE TRACER_MOD            ! Tracer array (STT) + related functions
      USE TRACERID_MOD          ! Flags for G-C tracers & chemical species

      !-----------------------------------------------------------------
      ! GEOS-Chem chemistry modules
      !-----------------------------------------------------------------
      USE CHEMISTRY_MOD         ! Driver routines for chemistry
      USE COMODE_MOD            ! Allocatable arrays for SMVGEAR solver
      USE COMODE_LOOP_MOD       ! Formerly common-block arrays for SMVGEAR
      USE GCKPP_COMODE_MOD      ! For the KPP chemical solver
      USE GLOBAL_CH4_MOD        ! For offline CH4 simulation
      USE MERCURY_MOD           ! For offline Hg simulation (driver)
      USE OCEAN_MERCURY_MOD     ! For offline Hg simulation (ocean model)
      USE SOAPROD_MOD           ! For SOA simulation
      USE STRAT_CHEM_MOD        ! For linearized stratospheric chemistry
      USE TOMS_MOD              ! For overhead O3 columns (for FAST-J)
      USE UVALBEDO_MOD          ! For reading UV albedoes (for FAST-J)

      !-----------------------------------------------------------------
      ! GEOS-Chem deposition modules
      !-----------------------------------------------------------------
      USE DEPO_MERCURY_MOD      ! Deposition for offline Hg simulation
      USE DRYDEP_MOD            ! For dry deposition
      USE WETSCAV_MOD           ! For wet deposition (rainout & washout)

      !-----------------------------------------------------------------
      ! GEOS-Chem diagnostics modules
      !-----------------------------------------------------------------
      USE BENCHMARK_MOD         ! For the 1-month benchmark simulations
      USE CMN_DIAG_MOD          ! Logical switches for G-C diagnostics
      USE DIAG_MOD              ! G-C diagnostic arrays & counters
      USE DIAG41_MOD            ! For ND41  (afternoon PBL       ) diag
      USE DIAG42_MOD            ! For ND42  (SOA products        ) diag
      USE DIAG48_MOD            ! For ND48  (station timeseries  ) diag
      USE DIAG49_MOD            ! For ND49  (inst. timeseries    ) diag
      USE DIAG50_MOD            ! For ND50  (24h avg timeseries  ) diag
      USE DIAG51_MOD            ! For ND51  (satellite timeseries) diag
      USE DIAG51b_MOD           ! For ND51b (satellite timeseries) diag
      USE DIAG63_MOD            ! For ND63  (PARANOX timeseries  ) diag
      USE DIAG_OH_MOD           ! For ND43  (OH,NO2,etc. prod    ) diag
      USE PLANEFLIGHT_MOD       ! For ND40  (plane flight track  ) diag

      !-----------------------------------------------------------------
      ! GEOS-Chem dynamics modules
      !-----------------------------------------------------------------
      USE CONVECTION_MOD        ! For deep cloud convection
      USE LINOZ_MOD             ! For LINOX linear strat chemistry
      USE PBL_MIX_MOD           ! For full PBL mixing (TURBDAY)
      USE TPCORE_BC_MOD         ! For nested-grid boundary conditions
      USE TRANSPORT_MOD         ! Driver routines for advection
      USE TROPOPAUSE_MOD        ! For the dynamic tropopause
      USE VDIFF_MOD             ! For non-local PBL mixing (J. Lin)

      !-----------------------------------------------------------------
      ! GEOS-Chem emissions modules
      !-----------------------------------------------------------------
      USE EMISSIONS_MOD         ! Driver routines for emissions
      USE MODIS_LAI_MOD         ! For MODIS leaf area indices (replacement)
      USE LIGHTNING_NOX_MOD     ! For lightning NOx emissions  
      USE MEGAN_MOD             ! For biogenic emissions
      USE BROMOCARB_MOD         ! For setting CH3Br concentrations in PBL, jpp

      !-----------------------------------------------------------------
      ! GEOS-Chem met field I/O modules
      !-----------------------------------------------------------------
      USE DAO_MOD               ! Met field definitions
      USE GCAP_READ_MOD         ! For reading GCAP met data
      USE GEOS57_READ_MOD       ! For reading GEOS-5.7.x data
      USE MERRA_A1_MOD          ! For reading MERRA A1 data
      USE MERRA_A3_MOD          ! For reading MERRA A3 data
      USE MERRA_CN_MOD          ! For reading MERRA CN data
      USE MERRA_I6_MOD          ! For reading MERRA I6 data
      USE A3_READ_MOD           ! For reading A3 data (all other met)
      USE A6_READ_MOD           ! For reading A6 data (all other met)
      USE I6_READ_MOD           ! For reading I6 data (all other met)

#if defined( DEVEL )
      USE LOGICAL_MOD,        ONLY : DO_DIAG_WRITE
      USE GC_ENVIRONMENT_MOD, ONLY : ALLOCATE_ALL, INIT_ALL
      USE GC_TYPE2_MOD,       ONLY : CHEMSTATE
      USE GC_TYPE_MOD,        ONLY : GC_MET_LOCAL
#endif

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
!  13 Aug 2010 - R. Yantosca - Added ProTeX headers
!  13 Aug 2010 - R. Yantosca - Add modifications for MERRA (treat like GEOS-5)
!  19 Aug 2010 - R. Yantosca - Now call MERRA met field reader routines
!  02 Feb 2011 - S. Kim      - Call Compute_OD after wet deposition
!  05 Oct 2011 - R. Yantosca - Now get SUNCOS30 array from routine COSSZA
!  07 Oct 2011 - R. Yantosca - Rename SUNCOS30 to SUNCOS_MID, which is the
!                              cos(SZA) at the midpt of the chemistry timestep
!  02 Feb 2012 - R. Yantosca - Added modifications for GEOS-5.7.x met fields
!  06 Feb 2012 - R. Yantosca - Reorganize USE statements for clarity
!  06 Feb 2012 - R. Yantosca - Renamed NN to NNN to avoid name confusion
!  07 Feb 2012 - R. Yantosca - Split off met field I/O into internal routines
!                              READ_INITIAL_MET_FIELDS and READ_MET_FIELDS
!  07 Feb 2012 - M. Payer    - Replace call to COSSZA with GET_COSINE_SZA
!  28 Feb 2012 - R. Yantosca - Removed support for GEOS-3
!  06 Mar 2012 - R. Yantosca - Now call READ_TOMS every month (this was
!                              formerly done within routine "fast_j.F")
!  06 Mar 2012 - R. Yantosca - Add subroutine GET_OVERHEAD_O3_FOR_FASTJ
!                              which calls COMPUTE_OVERHEAD_O3 (in toms_mod.F)
!                              to pre-compute the overhead O3 columsn for
!                              FAST-J photolysis.  This removes code from
!                              "set_prof.F" to facilitate the GI model.
!  19 Mar 2012 - R. Yantosca - Now call routines from olson_landmap_mod.F90
!                              to read the Olson land map data
!  04 Apr 2012 - R. Yantosca - Now call updated LAI routines from new module
!                              modis_lai_mod.F90.  Retire routine RDLAI.
!  05 Apr 2012 - R. Yantosca - Removed reference to LXTRA, it's obsolete
!  11 Apr 2012 - R. Yantosca - Replace lai_mod.F with modis_lai_mod.F90
!  11 Apr 2012 - R. Yantosca - Now call INIT_MODIS_LAI (in modis_lai_mod.F90)
!                              here so that we don't have to call it from 
!                              megan_mod.F and mercury_mod.F separately.
!  17 Apr 2012 - R. Yantosca - Need to set the mapping variable to NULL()
!  10 Jun 2012 - L. Murray   - Remove references to UPBDFLX_MOD.F
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Scalars
      LOGICAL                  :: FIRST = .TRUE.
      LOGICAL                  :: wasModisRead
      INTEGER                  :: I,             IOS,         J
      INTEGER                  :: K,             L,           N
      INTEGER                  :: JDAY,          NDIAGTIME,   N_DYN
      INTEGER                  :: NNN,           N_DYN_STEPS, NSECb 
      INTEGER                  :: N_STEP,        YEAR,        MONTH
      INTEGER                  :: DAY,           DAY_OF_YEAR, SEASON
      INTEGER                  :: NYMD,          NYMDb,       NHMS
      INTEGER                  :: ELAPSED_SEC,   NHMSb,       RC
      INTEGER                  :: ELAPSED_TODAY, HOUR,        MINUTE
      INTEGER                  :: laiYYYY,       laiMM
      REAL*8                   :: TAU,           TAUb         
      CHARACTER(LEN=255)       :: ZTYPE

      ! Arrays
      INTEGER                  :: DATE(2)

      ! Derived type objects
      TYPE(MapWeight), POINTER :: mapping(:,:) => NULL()

#if defined( DEVEL )
      TYPE(CHEMSTATE)          :: CHEM_STATE
      TYPE(GC_MET_LOCAL)       :: LOCAL_MET
#endif

      !=================================================================
      ! GEOS-CHEM starts here!                                            
      !=================================================================

      RC=0 ! Error flag. Default 0 is "no error"

#if defined( DEVEL )
      CALL ALLOCATE_ALL
#endif

      ! Display current grid resolution and data set type
      CALL DISPLAY_GRID_AND_MODEL

      !=================================================================
      !            ***** I N I T I A L I Z A T I O N *****
      !=================================================================

      ! Read input file and call init routines from other modules
      CALL READ_INPUT_FILE 
      IF ( LPRT ) CALL DEBUG_MSG( '### MAIN: a READ_INPUT_FILE' )

#if defined( DEVEL )
      ! Initialize all of CHEMSTATE
      CALL INIT_ALL(LOCAL_MET, CHEM_STATE)
#endif

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

      ! Define time variables for use below
      NHMS  = GET_NHMS()
      NHMSb = GET_NHMSb()
      NYMD  = GET_NYMD()
      NYMDb = GET_NYMDb()
      TAU   = GET_TAU()
      TAUb  = GET_TAUb()

      !=================================================================
      !        ***** I N I T I A L I Z A T I O N  continued *****
      !=================================================================

      ! Read the initial met fields from disk
      CALL READ_INITIAL_MET_FIELDS()

      ! Compute avg surface pressure near polar caps
      CALL AVGPOLE( PS1 )
      IF ( LPRT ) CALL DEBUG_MSG( '### MAIN: a AVGPOLE' )

      ! Call AIRQNT to compute air mass quantities from PS1
      CALL SET_FLOATING_PRESSURE( PS1 )
      IF ( LPRT ) CALL DEBUG_MSG( '### MAIN: a SET_FLT_PRS' )

#if defined( DEVEL )
      CALL AIRQNT(LOCAL_MET)
#else
      CALL AIRQNT
#endif
      IF ( LPRT ) CALL DEBUG_MSG( '### MAIN: a AIRQNT' )

      ! Compute lightning NOx emissions [molec/box/6h]
      IF ( LLIGHTNOX ) THEN
         CALL LIGHTNING
         IF ( LPRT ) CALL DEBUG_MSG( '### MAIN: a LIGHTNING' )
      ENDIF

      ! Initialize the derived type object containing
      ! mapping information for the MODIS LAI routines
      IF ( USE_OLSON_2001 ) THEN
         CALL Init_Mapping( 1440, 720, IIPAR, JJPAR, mapping )
      ELSE
         CALL Init_Mapping(  720, 360, IIPAR, JJPAR, mapping )
      ENDIF

      ! Compute the Olson land types that occur in each grid box
      ! (i.e. this is a replacement for rdland.F and vegtype.global)
      CALL Init_Olson_Landmap()
      CALL Compute_Olson_Landmap( mapping )
      CALL Cleanup_Olson_Landmap()

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

      ! Capture initial state of atmosphere for STE flux calc (ltm, 06/10/12)
      IF ( LSCHEM ) CALL INIT_STRAT_CHEM

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

      ! NSTEP is the number of dynamic timesteps w/in the outer loop
#if   defined( GEOS_57 ) 
      N_DYN_STEPS = 180 / GET_TS_DYN()     ! GEOS-5.7.x has a 3-hr interval
#else
      N_DYN_STEPS = 360 / GET_TS_DYN()     ! All other met has a 6hr interval
#endif

      ! Start a new outer loop
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
         DAY           = GET_DAY()
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
            IF ( DO_DIAG_WRITE ) THEN
               CALL DIAG3  

            ! Flush file units
               CALL CTM_FLUSH
            ENDIF
            !===========================================================
            !    *****  W R I T E   R E S T A R T   F I L E S  *****
            !===========================================================
            IF ( LSVGLB ) THEN

               ! Determine and print STE fluxes (ltm, 06/10/12)
               IF ( LSCHEM ) CALL CALC_STE

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
               NNN = GET_WETDEP_IDWETD( N )
               DO WHILE( .NOT.(IS_Hg2( NNN )) )
               
                  N = N + 1
                  ! Tracer number
                  NNN = GET_WETDEP_IDWETD( N )

               ENDDO

               CALL UPDATE_DEP( N )
            ENDIF
            CALL MAKE_GTMM_RESTART( NYMD, NHMS, TAU )
         ENDIF

         !==============================================================
         !       ***** T E S T   F O R   E N D   O F   R U N *****
         !==============================================================
         IF ( ITS_TIME_FOR_EXIT() ) GOTO 9999

         !==============================================================
         !          ****** R E A D   M E T   F I E L D S ******
         !
         !    For clarity, we have split off the code that reads met 
         !    fields into the internal subroutine READ_MET_FIELDS.
         !    Lightning NOx emissions are also computed at the same
         !    time when A3 or A6 met fields are read from disk.
         !==============================================================
         CALL READ_MET_FIELDS()

         !==============================================================
         ! ***** M O N T H L Y   O R   S E A S O N A L   D A T A *****
         !==============================================================
         IF ( LCHEM .and. ITS_A_NEW_MONTH() ) THEN

            ! Quantities for FAST-J photolysis
            IF ( ITS_A_FULLCHEM_SIM()  .or.
     &           ITS_AN_AEROSOL_SIM()  .or.
     &           ITS_A_CH3I_SIM()     ) THEN

               ! Read UV albedoes for FAST-J
               CALL READ_UVALBEDO( MONTH )

               ! Read TOMS overhead O3 columns for FAST-J
               CALL READ_TOMS( MONTH, YEAR )
            ENDIF
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

         ! Hg2 gas-particle partitioning (H Amos, 25 Oct 2011)
         IF ( ITS_A_MERCURY_SIM() .and. ITS_A_NEW_MONTH() ) THEN
            CALL PARTITION_Hg2( MONTH )
            IF ( LPRT ) CALL DEBUG_MSG( '### MAIN: a PARTITION_Hg2' )
         ENDIF

         !==============================================================
         !              ***** D A I L Y   D A T A *****
         !==============================================================
         IF ( ITS_A_NEW_DAY() ) THEN 

            !-----------------------------------------------------------
            !    %%%%%%% UPDATED LEAF AREA INDEX ALGORITHM %%%%%%%
            ! 
            ! Note that GEOS-Chem uses LAI data from two separate
            ! sources.  The dry deposition and soil NOx modules rely
            ! on the data from lai*.global ASCII files.  These files
            ! (which are pre-processed offline by IDL codes) are 
            ! generated for each specific GEOS-Chem grid configuration
            ! (e.g. 4x5, 2x25, 0.5x0.666 nested grids).  These files
            ! are read from disk by routine RDLAI, which saves
            ! the LAI data into the XLAI and XYLAI arrays.  XLAI and
            ! XYLAI store the leaf area index as a function of Olson 
            ! land type (cf Olson 1992 land map).
            !
            ! However, the MEGAN biogenic emissions code relies on 
            ! LAI data stored at 1x1 resolution stored in bpch format.
            ! These binary files are read by routine RDISOLAI (and
            ! other underlying routines in lai_mod.F), and are
            ! regridded on-the-fly to the current GEOS-Chem grid
            ! resolution.
            !
            ! Therefore, these two sources of LAI data present an
            ! inconsistency that should be resolved.  Also, for the 
            ! Grid-Indpendent GEOS-Chem project, we must move away
            ! from ASCII files (which prevent interfacing with 
            ! external GCMs).  We also cannot assume any particular
            ! horizontal grid, since that is now to be specified at
            ! the start of the simulation.
            !
            ! Also, to facilitate simulations at ultra-fine horizontal
            ! resolution, we will eventually adopt the Olson 2001 land 
            ! map, which has a native resolution of 0.25 x 0.25 
            ! degrees, and likewise use an updated version of the 
            ! MODIS LAI data at 0.25 x 0.25 resolution.
            !
            ! To resolve these issues, we have created a new module
            ! (modis_lai_mod.F90) which reads from the MODIS LAI data
            ! in netCDF format at the native resolution and then
            ! regrids the LAI data to GEOS-Chem resolution on-the-fly.
            ! The XLAI array is are populated for backwards
            ! compatibility with the existing legacy codes.  The LAI
            ! arrays used for MEGAN (ISOLAI, PMISOLAI, MISOLAI, and
            ! NMISOLAI) are now replaced by arrays GC_LAI, GC_LAI_PM, 
            ! GC_LAI_CM, and GC_LAI_NM) from modis_lai_mod.F.
            !
            ! We have validated that the new scheme generates identical 
            ! XLAI arrays w/r/t the old scheme.  The arrays GC_LAI etc. 
            ! differ from the ISOLAI etc. arrays slightly (but generally 
            ! agree to within 0.001).  This is due to the fact that the 
            ! ISOLAI arrays were regridded from 1 x 1 native resolution, 
            ! but now we are regridding from much finer resolution 
            ! (either 0.5 x 0.5 or 0.25 x 0.25).
            !
            ! NOTES:
            ! (1) At the present time, we have removed all references
            !     to the obsolete XYLAI array and its parent module
            !     CMN_VEL_mod.F.
            ! 
            !      -- Bob Yantosca (09 Apr 2012)
            !-----------------------------------------------------------

            ! Find the proper month and year for LAI  based on the current 
            ! month & year.  LAI months begin in the middle of the month.
            CALL Find_Lai_Month( day_of_year, month, year, 
     &                                        laiMM, laiYYYY )
            
            ! Read MODIS LAI from disk (if it is a new LAI month)
            ! Logical flag MODISREAD determines if data was just read in
            CALL Read_Modis_Lai( laiYYYY, laiMM, wasModisRead )

            ! Interpolate LAI quantities from native grid to GEOS-Chem grid
            CALL Compute_Modis_Lai( day_of_year, laiMM, 
     &                              mapping,     wasModisRead )

            ! For MEGAN biogenics, compute 15-day average temperature
            IF ( LMEGAN ) THEN
               CALL UPDATE_T_15_AVG
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
#if defined( DEVEL )
         CALL AIRQNT(LOCAL_MET)
#else
         CALL AIRQNT
#endif

         ! Compute the cosine of the solar zenith angle array
         ! SUNCOS         = cos(SZA) at the current time
         ! SUNCOS_MID     = cos(SZA) at the midpt of the chem timestep
         ! SUNCOS_MID_5hr = cos(SZA) at the midpt of the chem timestep 5hrs ago
         CALL GET_COSINE_SZA( SUNCOS, SUNCOS_MID, SUNCOS_MID_5hr )

         ! Compute tropopause height for ND55 diagnostic
         IF ( ND55 > 0 ) CALL DIAG_TROPOPAUSE

         ! Update dynamic timestep
         CALL SET_CT_DYN( INCREMENT=.TRUE. )

         !### Debug
         IF ( LPRT ) CALL DEBUG_MSG( '### MAIN: a INTERP, etc' )

         ! Get averaging intervals for local-time diagnostics
         ! (NOTE: maybe improve this later on)
         ! Placed after interpolation to get correct value of TROPP. 
         ! (ccc, 12/9/08)
         CALL DIAG_2PM

         ! Get the overhead column O3 for use with FAST-J
         CALL GET_OVERHEAD_O3_FOR_FASTJ()

         !==============================================================
         !   ***** U N I T   C O N V E R S I O N  ( kg -> v/v ) *****
         !==============================================================
         IF ( ITS_TIME_FOR_UNIT() ) THEN
            CALL CONVERT_UNITS( 1,  N_TRACERS, TCVV, AD, STT )

            !### Debug
            IF ( LPRT ) CALL DEBUG_MSG( '### MAIN: a CONVERT_UNITS:1' )
         ENDIF

!%%%!    ltm: As of v9-01-03, UPBDFLX and SCHEM are depreciated. These
!%%%!    subroutines and their children will be removed in the future.
!%%%!    !==============================================================
!%%%!    !     ***** S T R A T O S P H E R I C   F L U X E S *****
!%%%!    !==============================================================
!%%%!    IF ( LUPBD ) CALL DO_UPBDFLX

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
#if defined( DEVEL )
      CALL AIRQNT(LOCAL_MET)
#else
      CALL AIRQNT
#endif

!%%%!    ltm: As of v9-01-03, UPBDFLX and SCHEM are depreciated. These
!%%%!    subroutines and their children will be removed in the future.
!%%%!    ! Repartition [NOy] species after transport
!%%%!       IF ( LUPBD .and. ITS_A_FULLCHEM_SIM() ) THEN
!%%%!            CALL UPBDFLX_NOY( 2 )
!%%%!       ENDIF   

#if   !defined( GEOS_5 ) && !defined( MERRA ) && !defined( GEOS_57 )
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
#if defined( DEVEL )
            IF ( LEMIS ) CALL DO_EMISSIONS(CHEM_STATE)
#else
            IF ( LEMIS ) CALL DO_EMISSIONS
#endif

            IF ( ITS_TIME_FOR_UNIT() )
     &         CALL CONVERT_UNITS( 1, N_TRACERS, TCVV, AD, STT ) ! kg -> v/v


            IF ( ITS_A_FULLCHEM_SIM() ) THEN
               !========================================================
               !jpp, 2/12/08: putting a call to SET_CH3Br
               !              which is in bromocarb_mod.f
               !       ***** Fix CH3Br Concentration in PBL *****
               ! Kludge: eventually I want to keep the concentration
               !         entirely fixed! Ask around on how to...
               !========================================================
               IF ( LEMIS .and. ( IDTCH3Br > 0 ) )
     &            CALL SET_CH3Br( N_TRACERS, TCVV, AD, STT, 
     &                            ITS_TIME_FOR_UNIT() )

               ! ----------------------------------------------------
               ! If selected in input.geos, then set the MBL
               ! concentration of BrO equal to 1 pptv during daytime.
               ! ----------------------------------------------------
               IF ( LEMIS .and. ( IDTBrO > 0 ) )
     &            CALL SET_BRO( N_TRACERS, TCVV, AD, SUNCOS,
     &                          STT, ITS_TIME_FOR_UNIT() )

            ENDIF
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
                  CALL PARTITIONHG( 1, STT )
               ENDIF
      
               CALL DO_CONVECTION

               ! Return all reactive particulate Hg(II) to total Hg(II) tracer
               IF ( ITS_A_MERCURY_SIM() ) THEN
                  CALL PARTITIONHG( 2, STT )
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
#if defined( DEVEL )
            IF ( LEMIS ) CALL DO_EMISSIONS(CHEM_STATE)
#else
            IF ( LEMIS ) CALL DO_EMISSIONS
#endif

            IF ( ITS_A_FULLCHEM_SIM() ) THEN
               !========================================================
               !jpp, 2/12/08: putting a call to SET_CH3Br
               !              which is in bromocarb_mod.f
               !       ***** Fix CH3Br Concentration in PBL *****
               ! Kludge: eventually I want to keep the concentration
               !         entirely fixed! Ask around on how to...
               !========================================================
               IF ( LEMIS  .and. ( IDTCH3Br > 0 ) )
     &            CALL SET_CH3Br( N_TRACERS, TCVV, AD, STT,
     &                            ITS_TIME_FOR_UNIT() )

               ! ----------------------------------------------------
               ! If selected in input.geos, then set the MBL
               ! concentration of BrO equal to 1 pptv during daytime.
               ! ----------------------------------------------------
               IF ( LEMIS .and. ( IDTBrO > 0 ) )
     &            CALL SET_BRO( N_TRACERS, TCVV, AD, SUNCOS,
     &                          STT, ITS_TIME_FOR_UNIT() )

            ENDIF
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
#if defined( DEVEL )
! vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
! BELOW IS A TMP KLUGE TO TEST PASSING STT !
! TO "DO_CHEMISTRY" VIA THE "CHEM_STATE"   !
! DERIVED TYPE                             !
!
!    >-- STT SHOULD BE IN KG HERE --<
!
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
            CHEM_STATE%TRACERS(:,:,:,:N_TRACERS) = STT
! FIN DE KLUGE
            CALL DO_CHEMISTRY(CHEM_STATE)
#else
            CALL DO_CHEMISTRY
#endif

         ENDIF 

         ! Check STT (yxw)   
#if   defined( GEOS_5 ) && defined( GRID05x0666 )
         CALL CHECK_STT_05x0666( 'after chemistry' )
#endif 

         !==============================================================
         ! ***** W E T   D E P O S I T I O N  (rainout + washout) *****
         !==============================================================
         IF ( LWETD .and. ITS_TIME_FOR_DYN() ) THEN

            ! Add partition Hg(II) between aerosol and gas
            IF ( ITS_A_MERCURY_SIM() ) THEN
               CALL PARTITIONHG( 1, STT )
            ENDIF            
     
            ! Do wet deposition
            CALL DO_WETDEP
            
            ! Return all reactive particulate Hg(II) to total Hg(II) tracer
            IF ( ITS_A_MERCURY_SIM() ) THEN
               CALL PARTITIONHG( 2, STT )
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
         IF ( ND40 > 0 .and. DO_DIAG_WRITE ) THEN 
         
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

         ! Sometimes STT in the stratosphere can be negative at 
         ! the nested-grid domain edges. Force them to be zero before
         ! CHECK_STT (yxw)
#if   defined( GEOS_5 ) && defined( GRID05x0666 )
            CALL CHECK_STT_05x0666( 'after dynamics step' )
#endif

            CALL CHECK_STT( 'End of Dynamic Loop' )
         ENDIF
          
      ENDDO

      !=================================================================
      !     ***** C O P Y   I - 3   O R   I - 6   F I E L D S *****
      !
      !     The I-3 or I-6 fields at the end of this outer timestep 
      !     become the fields at the beginning of the next timestep
      !=================================================================
      CALL COPY_I3_I6_FIELDS
      IF ( LPRT ) CALL DEBUG_MSG( '### MAIN: after COPY_I3_I6_FIELDS' )

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

      ! Deallocate the derived type object w/ mapping info
      CALL CLEANUP_MAPPING( mapping )

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
!  02 Feb 2012 - R. Yantosca - Added output for GEOS-5.7.x met fields
!  28 Feb 2012 - R. Yantosca - Removed support for GEOS-3
!  19 Mar 2012 - R. Yantosca - Now echo info for 0.25 x 0.3125 runs
!  19 Mar 2012 - R. Yantosca - Now echo info if ISORROPIA is turned off
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! For system time stamp
      CHARACTER(LEN=16) :: STAMP

      !-------------------------------------------------
      ! Print resolution info
      !-------------------------------------------------
#if   defined( GRID4x5   )
      WRITE( 6, '(a)' )                   
     &    REPEAT( '*', 13 )                                          //
     &    '   S T A R T I N G   4 x 5   G E O S--C H E M   '         //
     &    REPEAT( '*', 13 )

#elif defined( GRID2x25  )
      WRITE( 6, '(a)' ) 
     &    REPEAT( '*', 13 )                                          // 
     &    '   S T A R T I N G   2 x 2.5   G E O S--C H E M   '       //
     &    REPEAT( '*', 13 )

#elif defined( GRID1x125 )
      WRITE( 6, '(a)' ) 
     &    REPEAT( '*', 13 )                                          // 
     &    '   S T A R T I N G   1 x 1.25   G E O S--C H E M   '      //
     &    REPEAT( '*', 13 )

#elif defined( GRID1x1 )
      WRITE( 6, '(a)' ) 
     &    REPEAT( '*', 13 )                                          // 
     &    '   S T A R T I N G   1 x 1   G E O S -- C H E M   '       //
     &    REPEAT( '*', 13 )

#elif defined( GRID05x0666 )
      WRITE( 6, '(a)' ) 
     &    REPEAT( '*', 13 )                                          // 
     &    '   S T A R T I N G   0.5 x 0.666   G E O S -- C H E M   ' //
     &    REPEAT( '*', 13 )

#elif defined( GRID025x03125 )
      WRITE( 6, '(a)' ) 
     &    REPEAT( '*', 13 )                                          // 
     &    '   S T A R T I N G  0.25 x 0.3125  G E O S -- C H E M   ' //
     &    REPEAT( '*', 13 )

#endif

      !-------------------------------------------------
      ! Print machine info
      !-------------------------------------------------

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

      !-------------------------------------------------
      ! Print met field info
      !-------------------------------------------------
#if   defined( GEOS_4 )
      WRITE( 6, '(a)' ) 'Using GMAO GEOS-4 met fields'
#elif defined( GEOS_5 )
      WRITE( 6, '(a)' ) 'Using GMAO GEOS-5 met fields'
#elif defined( GEOS_57 )
      WRITE( 6, '(a)' ) 'Using GMAO GEOS-5.7.x met fields'
#elif defined( MERRA )
      WRITE( 6, '(a)' ) 'Using GMAO MERRA met fields'
#elif defined( GCAP  )
      WRITE( 6, '(a)' ) 'Using GCAP/GISS met fields'
#endif

      !-------------------------------------------------
      ! Print msg if ISORROPIA is turned off for debug
      !-------------------------------------------------
#if   defined( NO_ISORROPIA ) 
      WRITE( 6, '(a)' ) REPEAT( '#', 39 )
      WRITE( 6, '(a)' ) '# ISORROPIA TURNED OFF FOR TESTING!!! #'
      WRITE( 6, '(a)' ) REPEAT( '#', 39 )
#endif
      
      !-------------------------------------------------
      ! System time stamp
      !-------------------------------------------------
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
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_initial_met_fields
!
! !DESCRIPTION: Internal subroutine READ\_INITIAL\_MET\_FIELDS calls the
!  various routines to read met fields at the beginning of a GEOS-Chem
!  simulation.  This code was moved out of the main routine for clarity,
!  due to the many \#if defined() blocks that are required.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE READ_INITIAL_MET_FIELDS()
! 
! !REMARKS:
!  All variables used in this routine are declared above in the main 
!  program, and as such, are visible here.
!                                                                             .
!  Also calls the following routines:
!  (1) AVGPOLE   (average pressure @ poles) when I3 or I6 fields are read
!  (2) LIGHTNING (lightning NOx emissions)  when A3 or A6 fields are read
!
! !REVISION HISTORY: 
!  07 Feb 2012 - R. Yantosca - Initial version
!  28 Feb 2012 - R. Yantosca - Removed support for GEOS-3
!EOP
!------------------------------------------------------------------------------
!BOC
      !=================================================================
      !    *****  R E A D   G E O S -- 5 . 7 , x  F I E L D S  *****
      !    *****  At the start of the GEOS-Chem simulation     *****
      !
      !    Handle GEOS-5.7.2 met fields separately from other met 
      !=================================================================
#if   defined( GEOS_57 )

      ! Read time-invariant data
      CALL GEOS57_READ_CN()

      ! Read 1-hr time-averaged data
      DATE = GET_FIRST_A1_TIME()
      CALL GEOS57_READ_A1( DATE(1), DATE(2) )
      
      ! Read 3-hr time averaged data
      DATE = GET_FIRST_A3_TIME()
      CALL GEOS57_READ_A3( DATE(1), DATE(2) )

      ! Read 3-hr time averaged data
      DATE = GET_FIRST_I3_TIME()
      CALL GEOS57_READ_I3_1( DATE(1), DATE(2) )

#elif defined( MERRA )

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
#if defined( DEVEL )
      CALL GET_A6_FIELDS(  DATE(1), DATE(2), LOCAL_MET )      
#else
      CALL GET_A6_FIELDS(  DATE(1), DATE(2) )      
#endif
      IF ( LPRT ) CALL DEBUG_MSG( '### MAIN: a 1st A6 TIME' )

      ! Open & read I-6 fields
      DATE = GET_FIRST_I6_TIME()
      CALL OPEN_I6_FIELDS(  DATE(1), DATE(2) )
      CALL GET_I6_FIELDS_1( DATE(1), DATE(2) )
      IF ( LPRT ) CALL DEBUG_MSG( '### MAIN: a 1st I6 TIME' )
      
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

      END SUBROUTINE READ_INITIAL_MET_FIELDS
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_met_fields
!
! !DESCRIPTION: Internal subroutine READ\_MET\_FIELDS calls the
!  various routines to read met fields in the main GEOS-Chem timestepping
!  loop.  This code was moved out of the main routine for clarity, due to
!  the many \#if defined() blocks that are required.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE READ_MET_FIELDS()
! 
! !REMARKS:
!  All variables used in this routine are declared above in the main 
!  program, and as such, are visible here.
!                                                                             .
!  Also calls the following routines:
!  (1) AVGPOLE   (average pressure @ poles) when I3 or I6 fields are read
!  (2) LIGHTNING (lightning NOx emissions)  when A3 or A6 fields are read
!
! !REVISION HISTORY: 
!  07 Feb 2012 - R. Yantosca - Initial version
!  28 Feb 2012 - R. Yantosca - Removed support for GEOS-3
!EOP
!------------------------------------------------------------------------------
!BOC

#if   defined( GEOS_57 )

      !==============================================================
      !  ****** R E A D   G E O S -- 5 . 7 . x   F I E L D S  *****
      !==============================================================

      !---------------------------------
      ! A-1 fields (1hr time averaged)
      !---------------------------------
      IF ( ITS_TIME_FOR_A1() ) THEN
         DATE = GET_A1_TIME()
         CALL GEOS57_READ_A1( DATE(1), DATE(2) )
      ENDIF

      !----------------------------------
      ! A-3 fields (3-hr time averaged)
      !----------------------------------
      IF ( ITS_TIME_FOR_A3() ) THEN
         DATE = GET_A3_TIME()
         CALL GEOS57_READ_A3( DATE(1), DATE(2) )

         ! Since CLDTOPS is an A-3 field, update the
         ! lightning NOx emissions [molec/box/6h]
         IF ( LLIGHTNOX ) CALL LIGHTNING
      ENDIF

      !----------------------------------
      ! I-3 fields (3-hr instantaneous
      !----------------------------------
      IF ( ITS_TIME_ FOR_I3() ) THEN
         DATE = GET_I3_TIME()
         CALL GEOS57_READ_I3_2( DATE(1), DATE(2) )

         ! Compute avg pressure at polar caps 
         CALL AVGPOLE( PS2 )
      ENDIF

#elif defined( MERRA )

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
#if defined( DEVEL )
         CALL GET_A6_FIELDS(  DATE(1), DATE(2), LOCAL_MET )      
#else
         CALL GET_A6_FIELDS(  DATE(1), DATE(2) )
#endif

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

      END SUBROUTINE READ_MET_FIELDS
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
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
      SUBROUTINE GET_OVERHEAD_O3_FOR_FASTJ()
! 
! !REMARKS:
!  This routine makes use of variables declared in above in the main program
!  (which are visible in all sub-programs below the CONTAINS statement).
!                                                                             .
!  The original code was done in FAST-J routine "set_prof.F", but has been
!  split off to facilitate development of the grid-independent model.
!
! !REVISION HISTORY: 
!   7 Mar 2012 - R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL :: USE_O3_FROM_MET

      ! FAST-J is only used for fullchem, offline aerosol
      ! and CH3I simulations, so skip otherwise
      IF ( ITS_A_FULLCHEM_SIM()  .or.
     &     ITS_AN_AEROSOL_SIM()  .or.
     &     ITS_A_CH3I_SIM()     ) THEN

         ! Only execute this if we are doing chemistry
         ! and if it we are at a chemistry timestep
         IF ( LCHEM .and. ITS_TIME_FOR_CHEM() ) THEN
             
#if   defined( GEOS_5 )

            ! When TOMS data is not available, then we will use
            ! the O3 columns directly from the met fields.  
            ! Currently, this is only done for GEOS-5 met.
            IF ( YEAR > LAST_TOMS_YEAR ) THEN
               USE_O3_FROM_MET = .TRUE.
            ELSE
               USE_O3_FROM_MET = .FALSE.
            ENDIF
#else
               
            ! When TOMS data is not available, then we will use
            ! the TOMS data from the last valid year.  This is
            ! done for all met field products other than GEOS-5.
            USE_O3_FROM_MET = .FALSE.

#endif

            ! Get the overhead O3 column for FAST-J.  Take either the
            ! TOMS O3 data or the column O3 directly from the met fields
            CALL COMPUTE_OVERHEAD_O3( DAY, USE_O3_FROM_MET, TO3 )
         ENDIF
      ENDIF

      END SUBROUTINE GET_OVERHEAD_O3_FOR_FASTJ
!EOC
      END PROGRAM GEOS_CHEM
#endif
