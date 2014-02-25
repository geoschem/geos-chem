!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: planeflight_mod
!
! !DESCRIPTION: Module PLANEFLIGHT\_MOD contains variables and routines which 
!  are used to "fly" a plane through the GEOS-Chem model simulation.  This is 
!  useful for comparing model results with aircraft observations. 
!\\
!\\
! !INTERFACE:
!
      MODULE PLANEFLIGHT_MOD
!
! !USES:
!
      USE inquireMod, ONLY : findFreeLUN

      IMPLICIT NONE
      PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC  :: ARCHIVE_RXNS_FOR_PF
      PUBLIC  :: CLEANUP_PLANEFLIGHT
      PUBLIC  :: PLANEFLIGHT
      PUBLIC  :: SETUP_PLANEFLIGHT
      PUBLIC  :: SET_PLANEFLIGHT
!
! !PRIVATE MEMBER FUNCTIONS:
!
      PRIVATE :: AN_SETUP
      PRIVATE :: INIT_PLANEFLIGHT
      PRIVATE :: NOY_SETUP
      PRIVATE :: READ_VARIABLES
      PRIVATE :: READ_POINTS
      PRIVATE :: RO2_SETUP
      PRIVATE :: TEST_VALID
      PRIVATE :: WRITE_VARS_TO_FILE
!
! !REMARKS:
!  The quantities that are saved to disk by the planeflight diagnostic were
!  requested by GEOS-Chem users.  If you would like to save out a new quantity,
!  then you will have to make your own modifications in this module.
!
! !REVISION HISTORY:
!  (1 ) Now references "pressure_mod.f" (dsa, bdf, bmy, 8/21/02)
!  (2 ) Now reference AD from "dao_mod.f".  Now also references "error_mod.f".
!        (bmy, 10/15/02)
!  (3 ) Bug fix: replace missing commas in FORMAT statement (bmy, 3/23/03)
!  (4 ) Now references "time_mod.f". (bmy, 3/27/03)
!  (5 ) Renamed PRATE to PRRATE to avoid conflict w/ SMVGEAR II (bmy, 4/1/03)
!  (6 ) Bug fix: use NAMEGAS instead of NAMESPEC (lyj, bmy, 7/9/03)
!  (7 ) Bug fix: avoid referencing JLOP for non-SMVGEAR runs (bmy, 7/18/03)
!  (8 ) Bug fix: Use T instead of T3 for GMAO temperature.  Also replace
!        NAMESPEC w/ NAMEGAS in RO2_SETUP.  Now locate reordered rxn 
!        numbers for SMVGEAR II.(tdf, mje, bmy, 8/1/03)
!  (9 ) Now print out N2O5 hydrolysis rxn as a special case.   Also rename
!        output file. (bmy, 8/8/03)
!  (10) Changed "DAO" to "GMAO" for met field variable names.  Now can save 
!        aerosol optical depths.  Bug fix in TEST_VALID. (bmy, 4/23/03)
!  (11) Now references "tracer_mod.f" (bmy, 7/20/04)
!  (12) Bug fix in READ_VARIABLES (1/7/05)
!  (13) Modified the plane flight diagnostic so that it writes output files
!        for each day where flight track files are defined. (bmy, 3/24/05)
!  (14) Minor bug fix in ARCHIVE_RXNS_FOR_PF (bmy, 5/20/05)
!  (15) Now split AOD's into column AOD's and AOD's below plane.  Also scale
!        AOD's to 400nm. (bmy, 10/25/05)
!  (16) Bug fixes in READ_VARIABLES (bmy, 10/16/06)
!  (17) Bug fix in PLANEFLIGHT (cdh, bmy, 12/12/06)
!  (18) Bug fix in RO2_SETUP (tmf, bmy, 4/23/07)
!  (19) Set very small values to zero.  (tmf, 1/7/09)
!  (20) Add new RO2 species according to 'globchem.dat' (tmf, 1/7/09) 
!  (21) Make sure we have 3 spaces in the exponential format (phs, 7/13/09)
!  (22) Output the grid cell indexes (kjw, 8/18/09)
!  (23) Add AN and NOy species. (fp, 3/10/10)
!  (24) Now scale AODs to wavelength specified in jv_spec_aod.dat(clh, 5/14/09)
!  29 Jul 2011 - R. Yantosca - Now also archive MERRA SEAICExx fields
!  29 Jul 2011 - R. Yantosca - Added ProTeX headers
!  03 Aug 2012 - R. Yantosca - Move calls to findFreeLUN out of DEVEL block
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE TYPES:
!
      !=================================================================
      ! MODULE VARIABLES:
      !
      ! DO_PF       : Turn on the planeflight diagnostic? (T/F)
      ! MAXVARS     : Maximum # of variables allowed
      ! MAXPOINTS   : Maximum # of flight track points allowed
      ! MAXREAC     : Maximum # of SMVGEAR reactions allowed
      ! MAXRO2      : Maximum # of RO2 constituents allowed
      ! NPOINTS     : Number of flight track points 
      ! PPOINT      : Pointer to last measured output
      ! PDATE       : Array of dates     at each flight point
      ! PTIME       : Array of times     at each flight point
      ! PTAU        : Array of TAU's     at each flight point
      ! PLAT        : Array of latitude  at each flight point
      ! PLON        : Array of longitude at each flight point
      ! PPRESS      : Array of pressure  at each flight point
      ! PTYPE       : Array of ID'#S     at each flight point
      ! NPVAR       : # of var's to be saved at each flight point
      ! PVAR        : Array of variable indices
      ! PNAME       : Array of variable names corresponding to PVAR
      ! NPREAC      : # of variables that are really SMVGEAR rxns
      ! PREAC       : Array of SMVGEAR rxn index numbers
      ! PRRATE      : Array of rxn rates for each entry in PREAC
      ! NRO2        : # number of RO2 constituents
      ! PRO2        : Array of SMVGEAR species that are RO2 const's
      ! INFILENAME  : Name of input file defining the flight track
      ! OUTFILENAME : Name of output file 
      !=================================================================

      ! Logicals
      LOGICAL                        :: DO_PF

      ! Parameters

#if   defined( TOMAS )
      INTEGER,           PARAMETER   :: MAXVARS   = 250  !(win, 7/28/09)
#else
      INTEGER,           PARAMETER   :: MAXVARS   = 200  !(sfarina, 1/28/13)
#endif
      INTEGER,           PARAMETER   :: MAXPOINTS = 10000
      INTEGER,           PARAMETER   :: MAXREAC   = 50
      INTEGER,           PARAMETER   :: MAXRO2    = 45
      INTEGER,           PARAMETER   :: MAXAN     = 10 ! FP 04/01/2010 (AN) 
      INTEGER,           PARAMETER   :: MAXNOY    = 12 ! FP 04/01/2010 (NOy)

      ! For specifying flight track points
      INTEGER                        :: NPOINTS           
      INTEGER                        :: PPOINT

      ! For specifying date/time
      INTEGER,           ALLOCATABLE :: PDATE(:)
      INTEGER,           ALLOCATABLE :: PTIME(:)
      ! Change to REAL*8 for more precise times (skim, 7/24/13)              
      REAL*8,            ALLOCATABLE :: PTAU(:)               

      ! For specifying lat/lon/alt and ID type
      REAL*4,            ALLOCATABLE :: PLAT(:)               
      REAL*4,            ALLOCATABLE :: PLON(:)               
      REAL*4,            ALLOCATABLE :: PPRESS(:)             
      CHARACTER(LEN=5),  ALLOCATABLE :: PTYPE(:)              

      ! For specifying variables to save at each flight point
      INTEGER                        :: NPVAR        
      INTEGER,           ALLOCATABLE :: PVAR(:) 
      CHARACTER(LEN=10), ALLOCATABLE :: PNAME(:)              
      
      ! For specifying SMVGEAR rxns to save at each flight point
      INTEGER                        :: NPREAC        
      INTEGER,           ALLOCATABLE :: PREAC(:) 
      REAL*8,            ALLOCATABLE :: PRRATE(:,:) 

      ! For specifying RO2 constituents at each flight point
      INTEGER                        :: NPRO2
      INTEGER                        :: PRO2(MAXRO2)

      ! For specifying NOY constituents at each flight point
      INTEGER                        :: NPNOY
      INTEGER                        :: PNOY(MAXNOY)

      ! For specifying AN constituents at each flight point
      INTEGER                        :: NPAN
      INTEGER                        :: P_AN(MAXAN)

      ! Input/output file names
      CHARACTER(LEN=255)             :: INFILENAME,  INF
      CHARACTER(LEN=255)             :: OUTFILENAME, OUTF

      ! Logical unit numbers (formerly in file_mod.F)
      INTEGER                        :: IU_FILE
      INTEGER                        :: IU_PLANE

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement 
      !=================================================================
      CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: setup_planeflight
!
! !DESCRIPTION: Subroutine SETUP\_PLANEFLIGHT reads information from the 
!  input file in order to initialize the planeflight diagnostic.  Also 
!  calls INIT\_PLANEFLIGHT to allocate and zero module arrays.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE SETUP_PLANEFLIGHT( am_I_Root )
!
! !USES:
!
      USE FILE_MOD,   ONLY : FILE_EXISTS
      USE FILE_MOD,   ONLY : IOERROR
      USE TIME_MOD,   ONLY : EXPAND_DATE
      USE TIME_MOD,   ONLY : GET_NYMD
      USE TIME_MOD,   ONLY : GET_NHMS
      USE TRACER_MOD, ONLY : ITS_A_FULLCHEM_SIM
!
! !INPUT PARAMETERS:
!
      LOGICAL, INTENT(IN) :: am_I_Root   ! Is this the root CPU?
!
! !REMARKS:
!  For SMVGEAR or KPP simulations, the call to SETUP_PLANEFLIGHT is made from 
!  routine "chemdr.f", after the "chem.dat" file is read.  This is necessary 
!  since we have to reference the SMVGEAR rxn rate and species numbers.
!                                                                             .
!  For offline simulations, the call to SETUP_PLANEFLIGHT can be made
!  at the start of the GEOS-Chem run (in "ndxx_setup.f" or similar routine).
! 
! !REVISION HISTORY: 
!  30 Jul 2002 - M. Evans    - Initial version
!  (1 ) Rename from "plane.dat" to "plane.log", since "*.dat" implies an input
!        file name. (bmy, 8/8/03)
!  (2 ) Add fancy output string (bmy, 4/26/04)
!  (3 ) Now references GET_NYMD, GET_NHMS, and EXPAND_DATE from "time_mod.f".
!        Now also replaces date & time tokens in the filenames. (bmy, 7/20/04)
!  (4 ) Now references FILE_EXISTS from "file_mod.f".  Modified so that we
!        check if a flight track file exists on each day.  Open file for 
!        output on each day and write header. (bmy, 3/25/05)
!  29 Jul 2011 - R. Yantosca - Added ProTeX headers
!  30 Jul 2012 - R. Yantosca - Now accept am_I_Root as an argument when
!                              running with the traditional driver main.F
!  06 Aug 2012 - R. Yantosca - Now use local IU_PLANE and not from file_mod.F
!  07 Aug 2012 - R. Yantosca - Now print LUN used to open file
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL, SAVE      :: FIRST = .TRUE.
      LOGICAL            :: IS_OPEN
      INTEGER            :: I,  IP,      N,   TEMP, LENGTH
      INTEGER            :: RN, COUNTER, IOS, NYMD, NHMS
      CHARACTER(LEN=7)   :: NAMES
      CHARACTER(LEN=20)  :: LINE
      CHARACTER(LEN=10)  :: TYPE

      !=================================================================
      ! SETUP_PLANEFLIGHT begins here!
      !=================================================================

      ! Assume that there is flight data for today
      DO_PF   = .TRUE.

      ! Find a free file LUN
      IU_FILE = findFreeLun()

      ! Get date & time
      NYMD    = GET_NYMD()
      NHMS    = GET_NHMS()

      ! Copy file names to local variables
      INF     = INFILENAME
      OUTF    = OUTFILENAME
      
      ! Replace any date & time tokens in the file names
      CALL EXPAND_DATE( INF,  NYMD, NHMS )
      CALL EXPAND_DATE( OUTF, NYMD, NHMS )

      ! If we can't find a flighttrack file for today's date, return
      IF ( .not. FILE_EXISTS( INF ) ) THEN 
         DO_PF = .FALSE.
         RETURN
      ENDIF

      ! Echo info
      IF ( am_I_Root ) THEN
         WRITE( 6, '(a)' ) REPEAT( '=', 79 )
         WRITE( 6, 99    )
         WRITE( 6, 100   ) TRIM( INF ), IU_FILE
 99      FORMAT( 'P L A N E   F L I G H T   D I A G N O S T I C' )
 100     FORMAT( /, 'SETUP_PLANEFLIGHT: Reading ',a, ' on unit ', i4 )
         WRITE( 6, '(a)' )  
      ENDIF

      ! Compute # of species and # of points & allocate arrays
      CALL INIT_PLANEFLIGHT( am_I_Root )

      ! Return if there are no flight track points for today
      IF ( NPOINTS == 0 ) THEN
         IF ( am_I_Root ) THEN
            WRITE( 6, '(a)' ) 'No flight track found for today!'
         ENDIF
         DO_PF = .FALSE.
         RETURN
      ENDIF

      !=================================================================
      ! Open file and read info
      !=================================================================

      ! Open input file
      OPEN( IU_FILE, FILE=TRIM( INF ), IOSTAT=IOS )
      IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'setup_planeflight:1')

      ! Read variables to be output -- sort into PVAR array by type
      CALL READ_VARIABLES( IU_FILE, am_I_Root )

      ! Read information about each point (date/time/lon/lat/alt)
      CALL READ_POINTS( IU_FILE, am_I_Root )

      ! Close the file
      CLOSE( IU_FILE )

      ! Set the pointer to the first record 
      PPOINT = 1

      !=================================================================
      ! Find the species # for all components of RO2 (SMVGEAR only)
      !=================================================================
      CALL RO2_SETUP

      !=================================================================
      ! Find the species # for all components of NOY (SMVGEAR only)
      !=================================================================
      CALL NOY_SETUP

      !=================================================================
      ! Find the species # for all components of AN (SMVGEAR only)
      !=================================================================
      CALL AN_SETUP
    
      ! Fancy output
      WRITE( 6, '(a)' ) REPEAT( '=', 79 )

      !=================================================================
      ! Open today's plane.log file and write file header
      !=================================================================

      ! Close previously opened plane.log file
      INQUIRE( IU_PLANE, OPENED=IS_OPEN )
      IF( IS_OPEN ) CLOSE( IU_PLANE )

      ! Close previously-opened file
      CLOSE( IU_PLANE )

      ! Find a free file LUN
      IU_PLANE = findFreeLUN()

      ! Open output file
      OPEN( IU_PLANE, FILE=TRIM( OUTF ), STATUS='UNKNOWN', IOSTAT=IOS )

      ! Error check
      IF ( IOS /= 0 ) THEN
         CALL IOERROR( IOS, IU_PLANE, 'setup_planeflight:1' )
      ENDIF

      ! Write header
      IF ( am_I_Root ) THEN
         WRITE( IU_PLANE, 110 ) 'POINT', 'TYPE', 'YYYYMMDD', 'HHMM',
     &        'LAT', 'LON', 'PRESS', ( PNAME(I), I=1,NPVAR )
      ENDIF

      ! FORMAT string
      ! Replace 95 by 250 (win)
      ! use 200 for non-Tomas sims (sfarina)

#if   defined( TOMAS )
 110  FORMAT( A5,X,A5,X,A8,X,A4,X,A7,X,A7,X,A7,X,250(a10,x) )
#else
 110  FORMAT( A5,X,A5,X,A8,X,A4,X,A7,X,A7,X,A7,X,200(a10,x) )
#endif

      END SUBROUTINE SETUP_PLANEFLIGHT
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_variables
!
! !DESCRIPTION: Subroutine READ\_VARIABLES reads the list of variables 
!  (SMVGEAR/KPP chemical species, SMVGEAR/KPP rxn rates, GMAO met fields, or 
!  GEOS-Chem tracers) to be printed out and sorts the information into the 
!  appropriate module variables.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE READ_VARIABLES( IU_FILE, am_I_Root )
!
! !USES:
!
      USE ERROR_MOD,  ONLY : GEOS_CHEM_STOP
      USE FILE_MOD,   ONLY : IOERROR
      USE TRACER_MOD, ONLY : N_TRACERS
      USE TRACER_MOD, ONLY : ITS_A_FULLCHEM_SIM
      USE CMN_SIZE_MOD                   ! Size parameters
      USE COMODE_LOOP_MOD                ! NAMEGAS, NSPEC
!
! !INPUT PARAMETERS:
!
      INTEGER, INTENT(IN) :: IU_FILE     ! Logical unit # for ASCII read
      LOGICAL, INTENT(IN) :: am_I_Root   ! Is this the root CPU?
! 
! !REVISION HISTORY:
!  30 Jul 2002 - M. Evans    - Initial version
!  (1 ) Now references GEOS_CHEM_STOP from "error_mod.f", which frees all
!        allocated memory before stopping the run. (bmy, 10/15/02)
!  (2 ) Bug fix: replace missing commas in FORMAT statement (bmy, 3/23/03)
!  (3 ) Bug fix: replace NAMESPEC w/ NAMEGAS for SMVGEAR II (lyj, bmy, 7/9/09)
!  (4 ) Now locate reordered rxn numbers for SMVGEAR II. (mje, bmy, 8/1/03)
!  (5 ) Now flag N2O5 hydrolysis rxn as a special case (bmy, 8/8/03)
!  (6 ) Changed variable name prefix "DAO" to "GMAO".  Also added aerosol
!        optical depths w/ tracer offset 2000. (bmy, 4/23/04)
!  (7 ) Now references N_TRACERS & ITS_A_FULLCHEM_SIM from "tracer_mod.f"
!        (bmy, 7/20/04)
!  (8 ) Bug fix: extract tracer # when reading rxn rates (bmy, 1/7/05)
!  (9 ) Now computes column AOD's and AOD's below plane (bmy, 10/24/05)
!  (10) We need to trim NAMEGAS before comparing to LINE so that comparisons 
!        for species like "O3" will work.  Also set NCS=NCSURBAN at the top
!        of the subroutine, to avoid out of bounds error. (dbm, bmy, 10/16/06)
!  (11) Add tracer TMS_?? for TOMAS microphysics rate diagnostic (win, 7/28/09)
!  29 Jul 2011 - R. Yantosca - Also search for MERRA SEAICExx met fields
!  29 Jul 2011 - R. Yantosca - Added ProTeX headers
!  30 Jul 2012 - R. Yantosca - Now accept am_I_Root as an argument when
!                              running with the traditional driver main.F
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL             :: IS_FULLCHEM
      INTEGER             :: M, N, NUM, R, IK, IOS
      CHARACTER(LEN=255)  :: LINE

      !=================================================================
      ! READ_VARIABLES begins here!
      !=================================================================

      ! Reset NCS to NCSURBAN for safety's sake (dbm, bmy, 10/16/06)
      NCS = NCSURBAN

      ! Test if this is a fullchem run
      IS_FULLCHEM = ITS_A_FULLCHEM_SIM()

      ! Read four lines of header
      DO N = 1, 4 
         READ( IU_FILE, '(a)', IOSTAT=IOS )
         IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'read_variables:1')
      ENDDO

      ! Read in the number of species to be output
      ! Read in as I4 now for the number of variables we save (skim, 7/24/13)
      READ( IU_FILE, '(i4)', IOSTAT=IOS ) NPVAR
      IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'read_variables:2' )

      ! Read in a separation line
      READ( IU_FILE, '(a)', IOSTAT=IOS )
      IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'read_variables:3' )

      ! Echo to stdout
      IF ( am_I_Root ) THEN
         WRITE( 6, '(a)' ) '   #    Species       PVAR'
         WRITE( 6, '(a)' ) '-----------------------------'
      ENDIF

      !=================================================================
      ! Sort variables by type; assign indices to PVAR, PREAC arrays
      ! NOTE: Variables for which PVAR(N) = 0 will be skipped!
      !=================================================================

      ! Zero reaction counter
      R = 0

      ! Loop over all variables
      DO N = 1, NPVAR

         ! Read each line
         READ( IU_FILE, '(a)', IOSTAT=IOS ) LINE
         IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'read_variables:4')

         ! Save the name of each variable into the global PNAME array
         PNAME(N) = LINE(1:10)

         ! We are searching for a ...
         SELECT CASE ( LINE(1:4) ) 

            !===========================================================
            ! GEOS-CHEM tracer: listed as "TRA_001", etc.
            ! PVAR offset: 100000
            !===========================================================
            CASE ( 'TRA_' )

               ! Extract tracer # from the string
               READ( LINE(5:14), '(i10)' ) NUM

               ! Make sure the tracer # is valid!
               IF ( NUM < 0 .or. NUM > N_TRACERS ) THEN                   
                  IF ( am_I_Root ) THEN
                     WRITE( 6, 100 ) TRIM( LINE )
 100                 FORMAT( 'TRACER ', i4, ' is out of range!' )
                     WRITE( 6, '(a)' ) 'STOP in SETUP_PLANEFLIGHT!'
                     WRITE( 6, '(a)' ) REPEAT( '=', 79 )
                  ENDIF
                  CALL GEOS_CHEM_STOP
               ENDIF

               ! Save in PVAR -- add offset of 100000
               PVAR(N) = 100000 + NUM

#if   defined( TOMAS )
            ! Add case matching for TOMAS microphysics rates (win, 7/28/09)
            !===========================================================
            ! GEOS-CHEM tracer: listed as "TMS_001", etc.
            ! PVAR offset: 200000
            !===========================================================
            CASE ( 'TMS_' )

               ! Extract tracer # from the string
               READ( LINE(5:14), '(i10)' ) NUM

               ! Make sure the tracer # is valid!
               IF ( NUM < 0 .or. NUM > N_TRACERS ) THEN                   
                  WRITE( 6, 100 ) TRIM( LINE )
                  WRITE( 6, '(a)' ) 'STOP in SETUP_PLANEFLIGHT!'
                  WRITE( 6, '(a)' ) REPEAT( '=', 79 )
                  CALL GEOS_CHEM_STOP
               ENDIF

               ! Save in PVAR -- add offset of 200000
               PVAR(N) = 200000 + NUM
#endif

            !===========================================================
            ! GMAO met field: listed as "GMAO_TEMP", etc.
            ! PVAR offset: 1000
            !===========================================================
            CASE ( 'GMAO' )
               
               IF ( LINE == 'GMAO_TEMP'  ) PVAR(N) = 1001
               IF ( LINE == 'GMAO_ABSH'  ) PVAR(N) = 1002
               IF ( LINE == 'GMAO_SURF'  ) PVAR(N) = 1003 
               IF ( LINE == 'GMAO_PSFC'  ) PVAR(N) = 1004
               IF ( LINE == 'GMAO_UWND'  ) PVAR(N) = 1005      
               IF ( LINE == 'GMAO_VWND'  ) PVAR(N) = 1006 
               IF ( LINE == 'GMAO_IIEV'  ) PVAR(N) = 1007
               IF ( LINE == 'GMAO_JJEV'  ) PVAR(N) = 1008
               IF ( LINE == 'GMAO_LLEV'  ) PVAR(N) = 1009
               IF ( LINE == 'GMAO_ICE00' ) PVAR(N) = 1100
               IF ( LINE == 'GMAO_ICE10' ) PVAR(N) = 1101
               IF ( LINE == 'GMAO_ICE20' ) PVAR(N) = 1102
               IF ( LINE == 'GMAO_ICE30' ) PVAR(N) = 1103
               IF ( LINE == 'GMAO_ICE40' ) PVAR(N) = 1104
               IF ( LINE == 'GMAO_ICE50' ) PVAR(N) = 1105
               IF ( LINE == 'GMAO_ICE60' ) PVAR(N) = 1106
               IF ( LINE == 'GMAO_ICE70' ) PVAR(N) = 1107
               IF ( LINE == 'GMAO_ICE80' ) PVAR(N) = 1108
               IF ( LINE == 'GMAO_ICE90' ) PVAR(N) = 1109

            !===========================================================
            ! Column aerosol optical depths (same order as for FAST-J)
            ! PVAR offset: 2000
            !===========================================================
            CASE ( 'AODC' )
            
               IF ( LINE == 'AODC_SULF'  ) PVAR(N) = 2001
               IF ( LINE == 'AODC_BLKC'  ) PVAR(N) = 2002
               IF ( LINE == 'AODC_ORGC'  ) PVAR(N) = 2003
               IF ( LINE == 'AODC_SALA'  ) PVAR(N) = 2004
               IF ( LINE == 'AODC_SALC'  ) PVAR(N) = 2005   

            !===========================================================
            ! Aerosol optical depths below the plane
            ! (same order as for FAST-J)  PVAR offset: 3000
            !===========================================================
            CASE ( 'AODB' )
            
               IF ( LINE == 'AODB_SULF'  ) PVAR(N) = 3001
               IF ( LINE == 'AODB_BLKC'  ) PVAR(N) = 3002
               IF ( LINE == 'AODB_ORGC'  ) PVAR(N) = 3003
               IF ( LINE == 'AODB_SALA'  ) PVAR(N) = 3004
               IF ( LINE == 'AODB_SALC'  ) PVAR(N) = 3005
 
            !===========================================================
            ! Hg(II) Partitioning - eds 10/27/11  PVAR offset: 4000
            !===========================================================
            CASE ( 'HG2_' )

               IF ( LINE == 'HG2_FRACG'  ) PVAR(N) = 4001
               IF ( LINE == 'HG2_FRACP'  ) PVAR(N) = 4002

            !===========================================================
            ! SMVGEAR rxn rate: listed as "REA_001", etc.
            ! PVAR offset: 10000
            !===========================================================
            CASE ( 'REA_' )

               ! Skip if not SMVGEAR!
               IF ( IS_FULLCHEM ) THEN 
               
                  ! Increment rxn counter
                  R = R + 1

                  IF ( TRIM( LINE ) == 'REA_O1D' ) THEN

                     ! O1D is a special rxn, give it offset of 20000
                     PVAR(N)  = 20000
                     PREAC(R) = 20000

                  ELSE IF ( TRIM( LINE ) == 'REA_N2O5' ) THEN

                     ! N2O5 hydrolysis is another special rxn
                     ! give it an offset of 21000
                     PVAR(N)  = 21000
                     PREAC(R) = 21000

                  ELSE
                     !==================================================
                     ! NOTE: the reaction numbers listed in smv2.log 
                     ! aren't really used to index SMVGEAR II rxns.  The 
                     ! rxns get reordered.  Find the right rxn number, 
                     ! which is stored in NOLDFNEW.  We assume only one 
                     ! chemistry scheme. (mje, bmy, 8/1/03)
                     !==================================================

                     ! Extract tracer # from the string
                     READ( LINE(5:14), '(i10)' ) NUM

                     ! Initialize
                     PVAR(N)  = -999
                     PREAC(R) = -999

                     ! Search for proper rxn number
                     DO IK = 1, NMTRATE 

                        ! Offset other reaction rates by 10000
                        IF ( NOLDFNEW(IK,1) == NUM ) THEN 
                           PVAR(N)  = 10000 + IK
                           PREAC(R) = 10000 + IK
                           EXIT
                        ENDIF
                     ENDDO

                     ! Stop w/ error 
                     IF ( PVAR(N) == -999 ) THEN 
                        IF ( am_I_Root ) THEN
                           WRITE (6,*) 'Cant match up reaction number'
                           WRITE (6,*) NUM
                           WRITE (6,*) 'Is it the second line of the'
                           WRITE (6,*) 'Three body reaction'
                           WRITE (6,*) 'Stopping'
                        ENDIF
                        CALL GEOS_CHEM_STOP
                     ENDIF
                  ENDIF
               ENDIF
               
            !===========================================================
            ! SMVGEAR chem species: listed as "O3", "C2H6", etc.
            ! PVAR offset: 0
            !===========================================================
            CASE DEFAULT

               ! Skip if not SMVGEAR!
               IF ( IS_FULLCHEM ) THEN

                  ! Loop over all SMVGEAR species -- 
                  ! match w/ species as read from disk
                  DO M = 1, NSPEC(NCS)
                     IF ( TRIM( NAMEGAS(M) ) == TRIM( LINE ) ) THEN
                        PVAR(N) = M
                        EXIT
                     ENDIF
                  ENDDO

                  ! Special flag for RO2 species
                  IF ( TRIM( LINE ) == 'RO2' ) PVAR(N) = 999

                  ! Special flag for AN species FP
                  IF ( TRIM( LINE ) == 'AN' ) PVAR(N) = 998

                  ! Special flag for NOy species FP
                  IF ( TRIM( LINE ) == 'NOy' ) PVAR(N) = 997

                  ! Error check
                  IF ( PVAR(N) == 0 ) THEN 
                     WRITE( 6, '(a)' ) 'ERROR: invalid species!'
                     WRITE( 6, 110   ) TRIM( LINE )
 110                 FORMAT( 'Species ', a, ' not found!' )
                     WRITE( 6, '(a)' ) 'STOP in PLANEFLIGHT!'
                     CALL GEOS_CHEM_STOP
                  ENDIF
               ENDIF

         END SELECT
      
         ! Echo species names/numbers to screen
         WRITE( 6, 120 ) N, TRIM( LINE ), PVAR(N)
 120     FORMAT( i4, 1x, a10, 1x, i10 )

      ENDDO

      END SUBROUTINE READ_VARIABLES
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_points
!
! !DESCRIPTION: Subroutine READ\_POINTS reads the information (ID, date, time, 
!  lat, lon, pressure) for each measurement listed in the input file, and 
!  sorts these into the appropriate module variables.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE READ_POINTS( IU_FILE, am_I_Root )
!
! !USES:
!
      USE BPCH2_MOD, ONLY : GET_TAU0
      USE ERROR_MOD, ONLY : GEOS_CHEM_STOP
      USE FILE_MOD,  ONLY : IOERROR
!
! !INPUT PARAMETERS:
!
      LOGICAL, INTENT(IN) :: am_I_Root   ! Is this the root CPU?
! 
! !REVISION HISTORY: 
!  30 Jul 2002 - M. Evans    - Initial version
!  (1 ) Now references GEOS_CHEM_STOP from "error_mod.f", which frees all
!        allocated memory before stopping the run. (bmy, 10/15/02)
!  29 Jul 2011 - R. Yantosca - Added ProTeX headers
!  30 Jul 2012 - R. Yantosca - Now accept am_I_Root as an argument when
!                              running with the traditional driver main.F
!EOP
!------------------------------------------------------------------------------
!BOC
      INTEGER, INTENT(IN) :: IU_FILE
!
! !LOCAL VARIABLES:
!
      INTEGER           :: N, IOS, QYY, QMM, QDD, QHH, QMN
      REAL*4            :: LAT, LON, PRES
      CHARACTER(LEN=10) :: TYPE

      !=================================================================
      ! READ_POINTS begins here!
      !=================================================================

      ! Read 4 header lines
      DO N = 1, 4
         READ( IU_FILE, '(a)', IOSTAT=IOS )
         IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'read_points:1' )
      ENDDO

      !=================================================================
      ! Read plane track points -- plane, lat/lon/alt, date/time
      ! We have previously computed NPOINTS in INIT_PLANEFLIGHT
      !=================================================================
      DO N = 1, NPOINTS
         
         ! Read a line from the file
         READ( IU_FILE, 100, IOSTAT=IOS ) 
     &        TYPE, QDD, QMM, QYY, QHH, QMN, LAT, LON, PRES
 100     FORMAT( 6x,a5,x,i2,x,i2,x,i4,x,i2,x,i2,x,f7.2,x,f7.2,x,f7.2 )

         ! Error check
         IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'read_points:2' )

         ! Exit if the word END is found
         IF ( INDEX( TYPE, 'END' ) > 0 ) EXIT

         !==============================================================
         ! Read date and time coordinates -- also do error checks
         !==============================================================

         ! Error check MONTH
         IF ( QMM < 1 .or. QMM > 12 ) THEN
            IF ( am_I_Root ) THEN
               WRITE( 6, 105 ) QMM
               WRITE( 6, 106 ) 
 105           FORMAT( 'ERROR: MONTH out of range: ', f8.3 )
 106           FORMAT( 'STOP in READ_POINTS (planeflight_mod.f)' )
            ENDIF
            CALL GEOS_CHEM_STOP
         ENDIF

         ! Error check DAY
         IF ( QDD < 1 .or. QDD > 31 ) THEN
            IF ( am_I_Root ) THEN
               WRITE( 6, 110 ) QDD
 111           WRITE( 6, 106 ) 
 110           FORMAT( 'ERROR: DAY out of range: ', f8.3 )
            ENDIF
            CALL GEOS_CHEM_STOP 
         ENDIF

         ! Error check HOUR
         IF ( QHH < 0 .or. QHH > 23 ) THEN
            IF ( am_I_Root ) THEN
               WRITE( 6, 115 ) QHH
               WRITE( 6, 106 )
 115           FORMAT( 'ERROR: HOUR out of range: ', f8.3 )
            ENDIF
            CALL GEOS_CHEM_STOP  
         ENDIF

         ! Error check MINUTES
         IF ( QMN < 0 .or. QMN > 59 ) THEN
            IF ( am_I_Root ) THEN
               WRITE( 6, 120 ) QMN
               WRITE( 6, 106 )
 120           FORMAT( 'ERROR: MINUTES out of range: ', f8.3 )
            ENDIF
            CALL GEOS_CHEM_STOP
         ENDIF

         ! Store type in the global PTYPE array
         PTYPE(N) = TYPE

         ! Store YYYYMMDD in the global PDATE array
         PDATE(N) = ( QYY * 10000 ) + ( QMM * 100 ) + QDD

         ! Store HHMMSS in the global PTIME array
         ! (actaully we read in just HHMM, assume seconds = 00)
         PTIME(N) = ( QHH * 100 ) + QMN

         ! Store TAU (hours since 1 Jan 1985) in the global PTAU array
         PTAU(N)  = GET_TAU0( QMM, QDD, QYY, QHH, QMN, 0 )

         !==============================================================
         ! Read lon/lat/alt coordinates -- also do error checks
         !==============================================================

         ! Put LONGITUDE in the range [-180...180]
         IF ( LON > 180.0 ) LON = LON - 360e0

         ! Error check LONGITUDE
         IF ( LON < -180 .OR. LON > 180 ) THEN 
            IF ( am_I_Root ) THEN
               WRITE( 6, 125 ) LON
               WRITE( 6, 106 ) 'STOP in READ_POINTS (planeflight_mod.f)'
 125           FORMAT( 'ERROR: Longitude out of range: ', f8.3 )
            ENDIF
            CALL GEOS_CHEM_STOP
         ENDIF

         ! Error check LATITUDE
         IF ( LAT < -90.0 .OR. LAT > 90.0 ) THEN 
            IF ( am_I_Root ) THEN
               WRITE( 6, 130 ) LAT
               WRITE( 6, 106 ) 'STOP in READ_POINTS (planeflight_mod.f)'
 130           FORMAT( 'ERROR: Latitude out of range: ', f8.3 )
            ENDIF
            CALL GEOS_CHEM_STOP
         ENDIF
        
         ! Assign LAT value into global PLAT array
         PLAT(N)   = LAT

         ! Assign LON value into global PLON array
         PLON(N)   = LON

         ! Assign PRES value into global PPRESS array
         PPRESS(N) = PRES

      ENDDO

      !=================================================================
      ! Echo number of points found and quit
      !=================================================================
      IF ( am_I_Root ) WRITE( 6, 135 ) NPOINTS
 135  FORMAT( /, 'Number of flight track points : ', i6 )

      END SUBROUTINE READ_POINTS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ro2_setup
!
! !DESCRIPTION: Subroutine RO2\_SETUP saves the SMVGEAR species indices of RO2 
!  constituents in the PRO2 array.  Also computes the count NPRO2. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE RO2_SETUP
!
! !USES:
!
      USE ERROR_MOD,  ONLY : GEOS_CHEM_STOP
      USE TRACER_MOD, ONLY : ITS_A_FULLCHEM_SIM

      USE CMN_SIZE_MOD     ! Size parameters
      USE COMODE_LOOP_MOD  ! NAMEGAS, NSPEC
!
! !REVISION HISTORY: 
!  01 Aug 2003 - M. Evans    - Initial version
!  (1 ) Now references GEOS_CHEM_STOP from "error_mod.f", which frees all
!        allocated memory before stopping the run. (bmy, 10/15/02)
!  (2 ) Now replace NAMESPEC w/ NAMEGAS for SMVGEAR II (bmy, 8/1/03)
!  (3 ) Now references ITS_A_FULLCHEM_SIM from "tracer_mod.f" (bmy, 7/20/04)
!  (4 ) Bug fix: PO3 should be PO2 (tmf, bmy, 4/23/07)
!  (5 ) NOTE: PO3 was a bug, that should have been PO2 (tmf, 2/10/09)
!  (6 ) Add new RO2 species according to 'globchem.dat' (tmf, 3/10/09)
!  29 Jul 2011 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER :: M

      !=================================================================
      ! RO2_SETUP begins here!
      !=================================================================
    
      ! Initialize 
      NPRO2 = 0

      ! We only need to proceed for SMVGEAR chemistry       
      IF ( .not. ITS_A_FULLCHEM_SIM() ) RETURN
      
      !=================================================================
      ! Loop over all SMVGEAR species, test for RO2 components
      !=================================================================
      DO M = 1, NSPEC(NCS)

         ! If we have found an RO2 compoent, add its species # to
         ! the PRO2 global array, and increment counter
         ! NOTE: PO3 was a bug, that should have been PO2 (tmf, 2/10/09) 
         SELECT CASE( TRIM( NAMEGAS(M) ) )

            CASE ( 'HO2',  'MO2',  'A3O2', 'ATO2', 'B3O2', 
     &             'ETO2', 'GCO3', 'IAO2', 'KO2',  'MAO3', 
     &             'MCO3', 'MRO2', 'PO2',  'RIO2', 'VRO2', 
     &             'ACO3', 'EO2', 'ENCO3', 'ENO2', 'GLCO3', 
     &             'IACO3', 'INO2', 'MACO3', 'NICO3', 'NIO2',
     &             'VOHRO2', 'RIO1', 'C59O2', 'MAPO2', 
     &             'BRO2', 'TRO2', 'XRO2', 'NRO2',
     &             'MOBAOO', 'ISOPNBO2', 'ISOPNDO2', 'IEPOXOO', 
     &             'PYPO2', 'MACRNO2', 'PMNO2', 'HC5OO', 'RCO3') 
               NPRO2       = NPRO2 + 1
               PRO2(NPRO2) = M

            CASE DEFAULT
               ! Nothing

         END SELECT

      ENDDO

      ! Error check
      IF ( NPRO2 > MAXRO2 ) THEN 
         WRITE( 6, '(a)' ) 'NPRO2 exceeds maximum allowed value!'
         WRITE( 6, '(a)' ) 'STOP in RO2_SETUP (planeflight_mod.f)'
         WRITE( 6, '(a)' ) REPEAT( '=', 79 )
         CALL GEOS_CHEM_STOP
      ENDIF

      !=================================================================
      ! Echo number of points found and quit
      !=================================================================
      WRITE( 6, 100 ) NPRO2
 100  FORMAT( 'Number of RO2 components      : ', i6 )

      END SUBROUTINE RO2_SETUP
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: noy_setup
!
! !DESCRIPTION: Subroutine NOY\_SETUP saves the SMVGEAR species indices of NOy 
!  constituents in the PNOY array.  Also computes the count NPNOY. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE NOY_SETUP
!
! !USES:
!
      USE ERROR_MOD,  ONLY : GEOS_CHEM_STOP
      USE TRACER_MOD, ONLY : ITS_A_FULLCHEM_SIM

      USE CMN_SIZE_MOD   ! Size parameters
      USE COMODE_LOOP_MOD       ! NSPEC, NAMEGAS, NCS
! 
! !REVISION HISTORY: 
!  01 Jun 2009 - F. Paulot   - Initial version
!  29 Jul 2011 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER :: M

      !=================================================================
      ! NOY_SETUP begins here!
      !=================================================================
    
      ! Initialize 
      NPNOY = 0

      ! We only need to proceed for SMVGEAR chemistry       
      IF ( .not. ITS_A_FULLCHEM_SIM() ) RETURN
      
      !=================================================================
      ! Loop over all SMVGEAR species, test for NOY components
      !=================================================================
      DO M = 1, NSPEC(NCS)

         SELECT CASE( TRIM( NAMEGAS(M) ) )

           CASE ( 'NO', 'NO2', 'NO3', 'HNO2', 'HNO4', 'HNO3',
     &            'PAN', 'PYPAN', 'PMN', 'PPN')

               NPNOY       = NPNOY + 1
               PNOY(NPNOY) = M

           CASE ( 'N2O5')

               NPNOY       = NPNOY + 1
               PNOY(NPNOY) = M

               NPNOY       = NPNOY + 1
               PNOY(NPNOY) = M

            CASE DEFAULT
               ! Nothing

         END SELECT

      ENDDO

      ! Error check
      IF ( NPNOY > MAXNOY ) THEN 
         WRITE( 6, '(a)' ) 'NPNOY exceeds maximum allowed value!'
         WRITE( 6, '(a)' ) 'STOP in NOY_SETUP (planeflight_mod.f)'
         WRITE( 6, '(a)' ) REPEAT( '=', 79 )
         CALL GEOS_CHEM_STOP
      ENDIF

      !=================================================================
      ! Echo number of points found and quit
      !=================================================================
      WRITE( 6, 100 ) NPNOY
 100  FORMAT( 'Number of NOY components      : ', i6 )

      END SUBROUTINE NOY_SETUP
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: an_setup
!
! !DESCRIPTION: Subroutine AN\_SETUP saves the SMVGEAR species indices of AN 
!  constituents in the P\_AN array.  Also computes the count NPAN.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE AN_SETUP
!
! !USES:
!
      USE ERROR_MOD,  ONLY : GEOS_CHEM_STOP
      USE TRACER_MOD, ONLY : ITS_A_FULLCHEM_SIM

      USE CMN_SIZE_MOD   ! Size parameters
      USE COMODE_LOOP_MOD       ! NSPEC, NAMEGAS, NCS
! 
! !REVISION HISTORY:
!  04 Jan 2010 - F. Paulot   - Initial version
!  29 Jul 2011 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER :: M

      !=================================================================
      ! AN_SETUP begins here!
      !=================================================================
    
      ! Initialize 
      NPAN = 0

      ! We only need to proceed for SMVGEAR chemistry       
      IF ( .not. ITS_A_FULLCHEM_SIM() ) RETURN
      
      !=================================================================
      ! Loop over all SMVGEAR species, test for AN components
      !=================================================================
      DO M = 1, NSPEC(NCS)

         ! If we have found an AN component, add its species # to
         ! the AN global array, and increment counter
         SELECT CASE( TRIM( NAMEGAS(M) ) )

           CASE ( 'ISOPNB', 'ISOPND', 'MVKN',  'MACRN',
     &            'ISN2',   'MNO3',   'ISNO3', 'R4N2',
     &            'PROPNN', 'ETHLN' ) 

               NPAN       = NPAN + 1
               P_AN(NPAN) = M

            CASE DEFAULT
               ! Nothing

         END SELECT

      ENDDO

      ! Error check
      IF ( NPAN > MAXAN ) THEN 
         WRITE( 6, '(a)' ) 'NPAN exceeds maximum allowed value!'
         WRITE( 6, '(a)' ) 'STOP in AN_SETUP (planeflight_mod.f)'
         WRITE( 6, '(a)' ) REPEAT( '=', 79 )
         CALL GEOS_CHEM_STOP
      ENDIF

      !=================================================================
      ! Echo number of points found and quit
      !=================================================================
      WRITE( 6, 100 ) NPAN
 100  FORMAT( 'Number of AN components      : ', i6 )

      END SUBROUTINE AN_SETUP
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: planeflight
!
! !DESCRIPTION: Subroutine PLANEFLIGHT saves concentrations to disk at 
!  locations corresponding to a flight track.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE PLANEFLIGHT( am_I_Root, Input_Opt,
     &                        State_Met, State_Chm, RC )
!
! !USES:
!
      USE COMODE_MOD,         ONLY : AIRDENS,        CSPEC,       JLOP
      USE COMODE_MOD,         ONLY : T3,             VOLUME,      ABSHUM
      USE COMODE_MOD,         ONLY : TAREA
#if   defined( TOMAS )
      USE DIAG_MOD,           ONLY : AD61_INST   ! (win, 7/28/09)
#endif
      USE ERROR_MOD,          ONLY : GEOS_CHEM_STOP  
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GIGC_State_Chm_Mod, ONLY : ChmState
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE PRESSURE_MOD,       ONLY : GET_PEDGE
      USE TIME_MOD,           ONLY : GET_TAU,        GET_TS_DIAG
      USE OCEAN_MERCURY_MOD,  ONLY : Fp, Fg !eds 10/27/11

      USE JV_CMN_MOD               ! ODAER, QAA, QAA_AOD
      USE COMODE_LOOP_MOD          ! CSPEC, etc.
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU?
      TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
      TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS:
!
      TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REVISION HISTORY: 
!  08 Jul 2002 - M. Evans    - Initial version
!  (1 ) Now reference AD from "dao_mod.f".  Now references GEOS_CHEM_STOP from
!        "error_mod.f", which frees memory before stopping. (bmy, 10/15/02)
!  (2 ) Now uses functions GET_TAU, GET_TS_CHEM from "time_mod.f".
!        (bmy, 3/27/03)
!  (3 ) Updated comments, cosmetic changes (bmy, 7/18/03)
!  (4 ) Now references T from "dao_mod.f", so that we can save out temperature
!        for non-SMVGEAR runs. (bmy, 8/1/03)
!  (5 ) Now references UWND and VWND from "dao_mod.f".  Now references
!        GET_PEDGE from "pressure_mod.f".  Added CASEs for surface pressure,
!        UWND, VWND to the CASE statement (bmy, 4/23/04)
!  (6 ) Now references STT & TCVV from "tracer_mod.f" (bmy, 7/20/04)
!  (7 ) Now return if DO_PF = .FALSE. (bmy, 3/24/05)
!  (8 ) Now compute column AOD's and AOD's below plane.  Also now scale
!        AOD's to 400nm. (bmy, 10/24/05)
!  (9 ) Bug fix: exit if PTAU(M) == PTAUE, so that we write out on the next !
!        planeflight timestep (cdh, bmy, 12/12/06)
!  (10) Change planeflight output time step. (ccc, 8/27/09)
!  (11) Add case matching for TOMAS rates (win, 7/28/09)
!  (12) Modify PTAUE calculation w/ ref to GET_TS_DYN  (win, 7/28/09)
!  (13) Now scale AOD's to jv_spec_aod.dat wavelength. (clh, 5/14/09)
!  29 Jul 2011 - R. Yantosca - Added ProTeX headers
!  09 Feb 2012 - R. Yantosca - Treat GEOS-5.7.x in the same way as MERRA
!  09 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!  25 Mar 2013 - R. Yantosca - Now accept am_I_Root, Input_Opt, State_Chm, RC
!  26 Sep 2013 - R. Yantosca - Renamed GEOS_57 Cpp switch to GEOS_FP
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL, SAVE       :: FIRST = .TRUE.
      LOGICAL             :: PCHEM
      INTEGER             :: I, IP, IRHN, J, L, JLOOP, M, N, R, RH, V
      INTEGER             :: LL     !eds
      INTEGER             :: LPLANE
      REAL*8              :: TK, PTAUS, PTAUE, CONSEXP, VPRESH2O, SAODnm
      REAL*8              :: VARI(NPVAR)

      ! Aerosol types: SULF, BLKC, ORGC, SALA, SALC 
      INTEGER             :: IND(5) = (/ 22, 29, 36, 43, 50 /)
!
! !DEFINED PARAMETERS:
!
      REAL*8,  PARAMETER  :: MISSING = -999.99999999d0  ! Missing data value
      REAL*8,  PARAMETER  :: TINY    = 1.d-36           ! arbitary small # to 
                                                        !  avoid faulty output

      ! Pointers
      ! We need to define local arrays to hold corresponding values 
      ! from the Chemistry State (State_Chm) object. (mpayer, 12/6/12)
      REAL*8, POINTER :: STT(:,:,:,:)

      !=================================================================
      ! PLANEFLIGHT begins here!
      !=================================================================

      ! Assume success
      RC = GIGC_SUCCESS

      ! Return if there is no flighttrack data for today
      IF ( .not. DO_PF ) RETURN

      ! Loop over all the locations that have not yet been found
      DO M = PPOINT, NPOINTS

         ! Starting & end times of chemistry interval
         PTAUE = GET_TAU()
         PTAUS = PTAUE - ( GET_TS_DIAG() / 60d0 )

         ! Initialize VARI to missing value for this point
         DO V = 1, NPVAR
            VARI(V) = MISSING
         ENDDO

         !==============================================================
         ! We haven't found the first plane point yet...
         !==============================================================
         IF ( PTAU(M) < PTAUS ) THEN

            PRINT*, 'BEGIN PLANEFLIGHT', PTYPE(M), PDATE(M), PTIME(M)
            PRINT*, PTAU(M), PTAUS

            ! Write all missing values to disk for point #M
            CALL WRITE_VARS_TO_FILE( M, VARI )      

            ! Increment pointer
            PPOINT = PPOINT + 1

         !==============================================================
         ! We have already found all of the plane points...
         !==============================================================
         ! Ensure that a model comparison is made for this point
         ! (skim, 7/24/13)
         ELSE IF ( PTAU(M) > PTAUE ) THEN

            PRINT*, 'END PLANEFLIGHT', PTYPE(M), PDATE(M), PTIME(M)
            PRINT*, PTAU(M), PTAUE

            ! Exit this loop and the subroutine
            EXIT

         !==============================================================
         ! We have found a plane point at the proper time & location!
         !==============================================================
         ELSE
            
            ! Print the flight track point number
            WRITE( 6, 100 ) PTYPE(M), PDATE(M), PTIME(M) 
 100        FORMAT( '     - PLANEFLIGHT: Archived ',a5,1x,i8.8,1x,i4.4 )

            ! Return grid box indices for the chemistry region
            ! NOTE: PCHEM and JLOOP are only defined for SMVGEAR runs!
            CALL TEST_VALID( M, PCHEM, JLOOP, I, J, L )

            ! Initialize SMVGEAR reaction counter
            R = 0

            ! Loop over all variables to save out
            DO V = 1, NPVAR

               ! Handle each variable
               SELECT CASE ( PVAR(V) )

                  !-------------------------
                  ! SMVGEAR species
                  !-------------------------
                  !997 and 998 are now special flags for NOy and AN
                  CASE ( 1:996) !( 1:998 )

                     ! Only archive where SMVGEAR chem is done
                     ! Save as mixing ratio [v/v]
                     IF ( JLOOP /= 0 ) THEN
                        VARI(V) = CSPEC(JLOOP,PVAR(V)) / AIRDENS(JLOOP)
                     ENDIF

                  ! FP 04/01/2010
                  !-------------------------
                  ! NOy family
                  !-------------------------   
                  CASE ( 997 )

                     ! Only archive where SMVGEAR chem is done
                     ! Sum all AN contributions, save as [v/v]
                     VARI(V) = 0d0

                     IF ( JLOOP /= 0 ) THEN 
                        DO N = 1, NPNOY
                           VARI(V) = VARI(V) + CSPEC(JLOOP,PNOY(N))
                        ENDDO

                        VARI(V) = VARI(V) / AIRDENS(JLOOP)
                     ENDIF

                  ! FP 04/01/2010
                  !-------------------------
                  ! AN family
                  !-------------------------   
                  CASE ( 998 )

                     ! Only archive where SMVGEAR chem is done
                     ! Sum all AN contributions, save as [v/v]
                     VARI(V) = 0d0
                     IF ( JLOOP /= 0 ) THEN 

                        DO N = 1, NPAN
                           VARI(V) = VARI(V) + CSPEC(JLOOP,P_AN(N))
                        ENDDO

                        VARI(V) = VARI(V) / AIRDENS(JLOOP)
                     ENDIF

                  !-------------------------
                  ! RO2 family
                  !-------------------------   
                  CASE ( 999 )

                     ! Only archive where SMVGEAR chem is done
                     ! Sum all RO2 contributions, save as [v/v]
                     VARI(V) = 0d0
                           
                        IF ( JLOOP /= 0 ) THEN
                           DO N = 1, NPRO2
                              VARI(V) = VARI(V) + CSPEC(JLOOP,PRO2(N))
                           ENDDO
                           
                           VARI(V) = VARI(V) / AIRDENS(JLOOP)
                        ENDIF

                  !--------------------------
                  ! GMAO temperature [K]
                  !--------------------------
                  CASE ( 1001 )
                     VARI(V) = State_Met%T(I,J,L)

                  !--------------------------
                  ! GMAO abs humidity [frac]
                  !--------------------------
                  CASE ( 1002 ) 
                     
                     ! Only archive where SMVGEAR chem is done
                     ! Code skalooched from "calcrate.f"
                        IF ( JLOOP /= 0 ) THEN
                           TK       = T3(JLOOP)
                           CONSEXP  = 17.2693882d0 * 
     &                                (TK - 273.16d0) / (TK - 35.86d0)
                           
                           VPRESH2O = CONSVAP * EXP(CONSEXP) * 1d0 / TK
                           
                           VARI(V)  = ABSHUM(JLOOP) * 
     &                                VPRESH2O      / AIRDENS(JLOOP)
                        ENDIF

                  !--------------------------
                  ! GMAO aerosol sfc area
                  !--------------------------
                  CASE ( 1003 )

                     ! Only archive where SMVGEAR chem is done
                        VARI(V) = 0d0

                        IF ( JLOOP /= 0 ) THEN
                           DO N = 1, NDUST + NAER
                              VARI(V) = VARI(V) + TAREA(JLOOP,N)
                           ENDDO
                        ENDIF

                  !--------------------------
                  ! GMAO sfc pressure [hPa]
                  !--------------------------
                  CASE ( 1004 )
                     VARI(V) = GET_PEDGE(I,J,1)

                  !-------------------------
                  ! GMAO U-wind [m/s]
                  !-------------------------
                  CASE ( 1005 )
                     VARI(V) = State_Met%U(I,J,L)

                  !--------------------------
                  ! GMAO V-wind [m/s]
                  !--------------------------
                  CASE ( 1006 )
                     VARI(V) = State_Met%V(I,J,L)

                  !--------------------------
                  ! GEOS-Chem Grid Box I
                  !--------------------------
                  CASE ( 1007 )
                     VARI(V) = I

                  !--------------------------
                  ! GEOS-Chem Grid Box J
                  !--------------------------
                  CASE ( 1008 )
                     VARI(V) = J

                  !--------------------------
                  ! GEOS-Chem Grid Box L
                  !--------------------------
                  CASE ( 1009 )
                     VARI(V) = L

#if   defined( MERRA ) || defined( GEOS_FP )
                  !--------------------------
                  ! GEOS-Chem SEAICE frac's
                  !--------------------------
                  CASE ( 1100 )
                     VARI(V) = State_Met%SEAICE00(I,J)
                  CASE ( 1101 )
                     VARI(V) = State_Met%SEAICE10(I,J)
                  CASE ( 1102 )
                     VARI(V) = State_Met%SEAICE20(I,J)
                  CASE ( 1103 )
                     VARI(V) = State_Met%SEAICE30(I,J)
                  CASE ( 1104 )
                     VARI(V) = State_Met%SEAICE40(I,J)
                  CASE ( 1105 )
                     VARI(V) = State_Met%SEAICE50(I,J)
                  CASE ( 1106 )
                     VARI(V) = State_Met%SEAICE60(I,J)
                  CASE ( 1107 )
                     VARI(V) = State_Met%SEAICE70(I,J)
                  CASE ( 1108 )
                     VARI(V) = State_Met%SEAICE80(I,J)
                  CASE ( 1109 )
                     VARI(V) = State_Met%SEAICE90(I,J)
#endif

                  !--------------------------
                  ! Column aerosol optical 
                  ! depths [unitless]
                  !--------------------------
                  CASE ( 2001:2005 )
                  
                     ! Only archive where SMVGEAR chem is done
                  
                        ! Remove MISSING flag
                        VARI(V) = 0d0
                  
                        ! Aerosol number
                        N = PVAR(V) - 2000
                  
                        ! Loop over RH bins
                        DO RH = 1, NRH
                        
                           ! Scaling factor for wavelength specified in
                           ! jv_spec_aod.dat
    !                       SAODnm  = QAA_AOD(IND(N)+RH-1) / 
    ! &                               QAA(4,IND(N)+RH-1)
                           ! We no longer need to scale by wavelength
                           ! (skim, 02/03/11)
                           SAODnm = 1.0

                           ! Index for type of aerosol and RH value
                           IRHN    = ( (N-1) * NRH ) + RH
                        
                           ! Sum AOD over all RH bins and store in VARI(V)
                           ! Sum over all vertical levels (bmy, 10/24/05)
                           VARI(V) = VARI(V) + 
     &                               SUM( SAODnm * ODAER(I,J,:,IRHN) )
                        ENDDO

                  !--------------------------
                  ! Aerosol optical depths
                  ! below plane [unitless]
                  !--------------------------
                  CASE ( 3001:3005 )
                  
                     IF ( PCHEM ) THEN

                       ! Remove MISSING flag
                        VARI(V) = 0d0
                  
                        ! Aerosol number
                        N = PVAR(V) - 3000
                  
                        ! Loop over RH bins
                        DO RH = 1, NRH
                        
                           ! Scaling factor for wavelength specified in
                           ! jv_spec_aod.dat
    !                       SAODnm  = QAA_AOD(IND(N)+RH-1) / 
    ! &                               QAA(4,IND(N)+RH-1)
                           ! We no longer need to scale by wavelength
                           ! (skim, 02/03/11)
                           SAODnm = 1.0 

                           ! Index for type of aerosol and RH value
                           IRHN    = ( (N-1) * NRH ) + RH
                          
                           ! Level of the plane.  AOD's are only computed
                           ! up to the tropopause, so if the plane goes into
                           ! the stratosphere, the AOD below plane will be
                           ! the same as the trop column at that point.
                           ! (bmy, 10/24/05)
                           LPLANE  = MIN( L, LLTROP )

                           ! Sum AOD over all RH bins and store in VARI(V)
                           ! Sum from surface to level where the plane is
                           VARI(V) = VARI(V) + 
     &                          SUM( SAODnm * ODAER(I,J,1:LPLANE,IRHN) )
                        ENDDO
                     ENDIF

                  !--------------------------
                  ! Hg(II) partitioning eds 10/27/11
                  !--------------------------
                  CASE ( 4001 )
                     VARI(V) = FG(I,J,LL) !L+1 sample 4/24/12

                  CASE ( 4002 )
                     VARI(V) = FP(I,J,LL) !L+1 sample 4/24/12

                  !--------------------------
                  ! SMVGEAR reaction rates
                  !--------------------------
                  CASE ( 10000:99999 )

                     ! Increment reaction count
                     R = R + 1

                     ! Only archive where SMVGEAR chem is done 
                     IF ( JLOOP /= 0 ) VARI(V) = PRRATE(JLOOP,R)

                  !--------------------------
                  ! GEOS-CHEM tracers [v/v]
                  !--------------------------
                  CASE( 100000:199999 )

                     ! Remove offset from PVAR
                     N = PVAR(V) - 100000

                     ! Initialize GEOS-Chem tracer array [kg] from
                     ! Chemistry State object (mpayer, 12/6/12)
                     STT => State_Chm%Tracers
                     
                     ! Convert from [kg] --> [v/v]
                     VARI(V) = STT(I,J,L,N) * Input_Opt%TCVV(N) /
     &                         State_Met%AD(I,J,L)

                     IF ( VARI(V) < TINY ) VARI(V) = 0.d0

                     ! Free pointer
                     NULLIFY( STT )

#if   defined( TOMAS )
                  !-------------------------------
                  ! TOMAS microphysics rate [kg/s] or [no./cm3/s]
                  !-------------------------------
                  CASE( 200000:299999 )

                     ! Remove offset from PVAR
                     N = PVAR(V) - 200000
                     
                     ! Archive the microphysics rate
                     VARI(V) = AD61_INST(I,J,L,N)
!debug
                     write (6,*) 'ARCHIVE TO PLANEFLIGHT DIAG',
     &                    'AD61_INST at',I,J,L,N,'=',AD61_INST(I,J,L,N)

#endif

                  !--------------------------
                  ! Otherwise it's an error!
                  !--------------------------
                  CASE DEFAULT
                     WRITE( 6, '(a)' ) REPEAT( '=', 79 )
                     WRITE( 6, '(a)' ) 'PLANEFLIGHT: Bad variable #!' 
                     WRITE( 6, '(a)' ) 'STOP in PLANEFLIGHT!'
                     WRITE( 6, '(a)' ) REPEAT( '=', 79 )
                     CALL GEOS_CHEM_STOP

               END SELECT
            ENDDO

            ! Write data for the Mth plane point out to disk
            CALL WRITE_VARS_TO_FILE( M, VARI )

            ! Increment the record pointer
            PPOINT = PPOINT + 1

         ENDIF
      ENDDO

      END SUBROUTINE PLANEFLIGHT
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: test_valid
!
! !DESCRIPTION: Subroutine TEST\_VALID tests to see if we are w/in the 
!  tropopause, which is where SMVGEAR chemistry is done.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE TEST_VALID( IND, PCHEM, JLOOP, I, J, L ) 
!
! !USES:
!
      USE COMODE_MOD,     ONLY : JLOP
      USE PRESSURE_MOD,   ONLY : GET_PEDGE
      USE TRACER_MOD,     ONLY : ITS_A_FULLCHEM_SIM
      USE GRID_MOD,       ONLY : GET_XOFFSET
      USE GRID_MOD,       ONLY : GET_YOFFSET

      USE CMN_SIZE_MOD                ! Size parameters
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN)  :: IND     ! # of the flight track point
!
! !OUTPUT PARAMETERS:
!
      LOGICAL, INTENT(OUT) :: PCHEM   ! =T if chemistry is done here
      INTEGER, INTENT(OUT) :: JLOOP   ! 1-D grid box index for SMVGEAR/KPP
      INTEGER, INTENT(OUT) :: I       ! GEOS-Chem longitude index
      INTEGER, INTENT(OUT) :: J       ! GEOS-Chem latitude index
      INTEGER, INTENT(OUT) :: L       ! GEOS-Chem level index
! 
! !REVISION HISTORY: 
!  08 Jul 2002 - M. Evans    - Initial version
!  (1 ) Now use GET_PEDGE of "pressure_mod.f" to return the pressure at the
!        bottom edge of box (I,J,L), for hybrid grid. (dsa, bdf, bmy, 8/21/02)
!  (2 ) Since JLOP is not allocated for non-SMVGEAR runs, set PCHEM=F and 
!        JLOOP=0 even if we are in the troposphere. (bmy, 7/18/03)
!  (3 ) Bug fix: add 0.5 in expression for I so that the rounding will
!        be done correctly.  Also make sure that I is computed correctly
!	 for points near the date line.  (bmy, 4/23/04)
!  (4 ) Now references ITS_A_FULLCHEM_SIM from "tracer_mod.f" (bmy, 7/20/04)
!  (5 ) Now references ITS_IN_THE_TROP from "tropopause_mod.f" (bmy, 8/22/05)
!  (6 ) Reference GET_XOFFSET and GET_YOFFSET from "grid_mod.f" and also 
!        add for the case of nested-grid simulation (win, 7/28/09)
!  29 Jul 2011 - R. Yantosca - Added ProTeX headers
!  08 Sep 2011 - L. Schiferl - Added correct definitions for I and J 
!                              based on nested regions
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER :: IL
      LOGICAL :: FOUND
      INTEGER :: I0, J0   ! (lds, 8/25/11)

      !=================================================================
      ! TEST_VALID begins here!
      !=================================================================

      ! We have not found a valid point
      FOUND = .FALSE.

      ! Added correct definitions for I and J based on nested regions 
      ! (lds, 8/25/11)
      I0 = GET_XOFFSET( GLOBAL=.TRUE. )
      J0 = GET_YOFFSET( GLOBAL=.TRUE. )

      ! Get I corresponding to PLON(IND)
      I = INT( ( PLON(IND) + 180d0 - (I0 * DISIZE) ) / DISIZE + 1.5d0 )

      ! Handle date line correctly (bmy, 4/23/04)
      IF ( I > IIPAR ) I = I - IIPAR

      ! Get J corresponding to PLAT(IND)
      J = INT( ( PLAT(IND) +  90d0 - (J0 * DJSIZE) ) / DJSIZE + 1.5d0 )

      ! Get L corresponding to PRESS(IND)
      L = 1
      DO IL = 1, LLPAR
         IF ( GET_PEDGE(I,J,IL) <= PPRESS(IND) .AND..NOT. FOUND ) THEN
            L     = IL-1
            FOUND =.TRUE.
            EXIT
         ENDIF          
      ENDDO

      ! Error check: L must be 1 or higher
      IF ( L == 0 ) L = 1

      !=================================================================
      ! We only do full-chemistry in the troposphere
      !=================================================================
      IF ( ITS_A_FULLCHEM_SIM() ) THEN
  
         ! JLOOP indicates if a box is in tropo (/=0) or not. 
         JLOOP = JLOP(I,J,L)
         PCHEM = .TRUE.

      ELSE

         ! For non-SMVGEAR runs, JLOOP has no meaning so we give the
         ! stratospheric value
         JLOOP = 0
         PCHEM = .FALSE.

      ENDIF

      END SUBROUTINE TEST_VALID
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: write_vars_to_file
!
! !DESCRIPTION: Subroutine WRITE\_VARS\_TO\_FILE writes the values of all 
!  the variables for a given flight track point to the output file.  
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE WRITE_VARS_TO_FILE( IND, VARI )
!
! !USES:
!
      USE FILE_MOD, ONLY : IOERROR
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN) :: IND           ! # of the flight track point
      REAL*8,  INTENT(IN) :: VARI(NPVAR)   ! Values to print to file
! 
! !REVISION HISTORY: 
!  08 Jul 2002 - M. Evans    - Initial version
!  (1 ) The max line length for output seems to be 1024 characters.  Adjust
!        MAXVARS accordingly so that we don't exceed this. (bmy, 7/8/02)
!  (2 ) Now do not write file header -- this is now done in subroutine
!        SETUP_PLANEFLIGHT at the start of each day (bmy, 3/25/05)
!  (3 ) Bug fix: make sure we have 3 spaces in exponential (phs, 7/13/09)
!  29 Jul 2011 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL, SAVE :: FIRST = .TRUE.
      INTEGER       :: I, IOS

      !=================================================================
      ! WRITE_VARS_TO_FILE begins here!
      !=================================================================

      ! Write data to file
      WRITE( IU_PLANE, 110, IOSTAT=IOS ) 
     &     IND, PTYPE(IND), INT( PDATE(IND) ), INT( PTIME(IND) ),
     &     PLAT(IND), PLON(IND), PPRESS(IND), ( VARI(I), I=1,NPVAR )
      
      ! Format string
      ! FP replace 95 by 200
! 110  FORMAT( I5,   X, A5,   X, I8.8, X, I4.4, X, 
!     &        F7.2, X, F7.2, X, F7.2, X, 95(es11.3e3,x) )
#if   defined( TOMAS )
 110  FORMAT(I5,X,A5,X,I8.8,X,I4.4,X,F7.2,X,F7.2,X,F7.2,X,
     &       250(es11.3e3,x))
#else
 110  FORMAT(I5,X,A5,X,I8.8,X,I4.4,X,F7.2,X,F7.2,X,F7.2,X,
     &       200(es11.3e3,x))
#endif

      ! Error check
      IF ( IOS /= 0 ) CALL IOERROR( IOS,IU_PLANE,'write_vars_to_file:1')

      ! Flush the file to disk
      CALL FLUSH( IU_PLANE )

      END SUBROUTINE WRITE_VARS_TO_FILE
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: archive_rxns_for_PF
!
! !DESCRIPTION: Subroutine ARCHIVE\_RXNS\_FOR\_PF is called from "calcrate.f" 
!  to pass reaction rates from the SMVGEAR solver for the planeflight 
!  diagnostic. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE ARCHIVE_RXNS_FOR_PF( JO1D, N2O5 )
!
! !USES:
!
      USE COMODE_MOD, ONLY : IXSAVE
      USE COMODE_MOD, ONLY : IYSAVE
      USE COMODE_MOD, ONLY : IZSAVE
      USE ERROR_MOD,  ONLY : GEOS_CHEM_STOP

      USE CMN_SIZE_MOD                     ! Size parameters
      USE COMODE_LOOP_MOD                  ! RRATE, JLOOPLO, KBLOOP
      USE CMN_DIAG_MOD                     ! ND40 switch
!
! !INPUT PARAMETERS: 
!
      REAL*8, INTENT(IN)  :: JO1D(KBLOOP)  ! JO1D photolysis rate [1/s]
      REAL*8, INTENT(IN)  :: N2O5(KBLOOP)  ! N2O5 hydrolysis rate [molec/cm3/s]
! 
! !REVISION HISTORY: 
!  08 Jul 2002 - M. Evans    - Initial version
!  (1 ) Now avoid overflow/underflow errors in PRATE (bmy, 7/8/02)
!  (2 ) Now reference GEOS_CHEM_STOP from "error_mod.f", which frees all
!        allocated memory before stopping the run (bmy, 10/15/02)
!  (3 ) Renamed PRATE to PRRATE to avoid conflict w/ SMVGEAR II (bmy, 4/1/03)
!  (4 ) Now also pass N2O5 hydrolysis rxn rate array via the arg list.  
!        Also bug fix: replace TMP with RATE in under/overflow checking
!        for JO1D and N2O5. (bmy, 8/8/03)      
!  (5 ) Bug fix: Replace with DO_PF since this variable is reset to either T 
!        or F each day depending on whether there is plane flight data 
!        available (bmy, 5/20/05)
!  29 Jul 2011 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER            :: KLOOP, JLOOP, V, R, I, J, L 
      REAL*8             :: RATE
!
! !DEFINED PARAMETERS:
!
      ! Smallest, largest REAL*4 #'s representable on this machine
      REAL*4, PARAMETER  :: SMALLEST=TINY(1e0), LARGEST=HUGE(1e0)

      !=================================================================
      ! ARCHIVE_RXNS_FOR_PF begins here!
      !=================================================================
      IF ( DO_PF ) THEN

         ! Loop over SMVGEAR reactions
         DO R = 1, NPREAC

            ! Test SMVGEAR rxn number
            SELECT CASE ( PREAC(R) ) 

               !-----------------------
               ! All except JO1D, N2O5
               !-----------------------
               CASE( 10000:19999 )
                         
                  ! Store rate in PRRATE
                  DO KLOOP = 1, KTLOOP
                     JLOOP = JLOOPLO + KLOOP
                     RATE  = RRATE(KLOOP,PREAC(R)-10000)

                     ! Avoid overflow/underflow
                     IF ( RATE < SMALLEST ) RATE = 0e0
                     IF ( RATE > LARGEST  ) RATE = LARGEST

                     PRRATE(JLOOP,R) = RATE
                  ENDDO

               !-----------------------
               ! JO1D photolysis rxn 
               !-----------------------
               CASE ( 20000 )

                  ! Store rate in PRATE
                  DO KLOOP = 1, KTLOOP
                     JLOOP = JLOOPLO + KLOOP
                     RATE  = JO1D(KLOOP)

                     ! Avoid overflow/underflow
                     IF ( RATE < SMALLEST ) RATE = 0e0
                     IF ( RATE > LARGEST  ) RATE = LARGEST 

                     PRRATE(JLOOP,R) = RATE
                  ENDDO
                 
               !-----------------------
               ! N2O5 hydrolysis rxn
               !-----------------------
               CASE ( 21000 ) 

                  ! Store rate in PRATE
                  DO KLOOP = 1, KTLOOP
                     JLOOP = JLOOPLO + KLOOP
                     RATE  = N2O5(KLOOP)

                     ! Avoid overflow/underflow
                     IF ( RATE < SMALLEST ) RATE = 0e0
                     IF ( RATE > LARGEST  ) RATE = LARGEST 

                     PRRATE(JLOOP,R) = RATE
                  ENDDO

               !-----------------------
               ! Error: invalid rxn
               !-----------------------
               CASE DEFAULT
                  WRITE( 6, '(a)' ) REPEAT( '=', 79 )
                  WRITE( 6, '(a)' ) 'ERROR -- Invalid SMVGEAR rxn #!'
                  WRITE( 6, '(a)' ) 'STOP in ARCHIVE_RXNS_FOR_PF!'
                  WRITE( 6, '(a)' ) REPEAT( '=', 79 )
                  CALL GEOS_CHEM_STOP

            END SELECT
         ENDDO
      ENDIF

      END SUBROUTINE ARCHIVE_RXNS_FOR_PF
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: set_planeflight
!
! !DESCRIPTION: Subroutine SET\_PLANEFLIGHT is used to pass values read in 
!  from the GEOS-Chem input file to "planeflight\_mod.f".
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE SET_PLANEFLIGHT( PF, IN_FILE, OUT_FILE )
!
! !INPUT PARAMETERS: 
!
      LOGICAL,            INTENT(IN) :: PF         ! Turn on planeflight diag?
      CHARACTER(LEN=255), INTENT(IN) :: IN_FILE    ! Input file to read
      CHARACTER(LEN=255), INTENT(IN) :: OUT_FILE   ! Output file to write
! 
! !REVISION HISTORY: 
!  20 Jul 2004 - R. Yantosca - Initial version
!  29 Jul 2011 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
      ! Save arguments to "shadow" module variables
      DO_PF       = PF
      INFILENAME  = TRIM( IN_FILE  )
      OUTFILENAME = TRIM( OUT_FILE )

      END SUBROUTINE SET_PLANEFLIGHT
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_planeflight
!
! !DESCRIPTION: Subroutine INIT\_PLANEFLIGHT reads the input file to compute 
!  the number of variables and flight track points to print out.  Also 
!  allocates all module arrays.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE INIT_PLANEFLIGHT( am_I_Root )
!
! !USES:
!
      USE ERROR_MOD, ONLY : ALLOC_ERR
      USE ERROR_MOD, ONLY : GEOS_CHEM_STOP
      USE FILE_MOD,  ONLY : IOERROR
      USE CMN_SIZE_MOD    ! Size Parameters
      USE COMODE_LOOP_MOD ! ITLOOP
!
! !INPUT PARAMETERS:
!
      LOGICAL, INTENT(IN) :: am_I_Root   ! Is this the root CPU?
! 
! !REVISION HISTORY: 
!  08 Jul 2002 - M. Evans    - Initial version
!  (1 ) Now reference GEOS_CHEM_STOP from "error_mod.f", which frees all
!        allocated memory before stopping the run.  Also reference ALLOC_ERR
!        from "error_mod.f" (bmy, 10/15/02)
!  (2 ) Renamed PRATE to PRRATE to avoid conflict w/ SMVGEAR II (bmy, 4/1/03)
!  (3 ) INIT_PLANEFLIGHT is now called each day but the arrays are only
!        allocated once.  Arrays are now allocated to the maximum size.
!        (bmy, 3/25/05)
!  29 Jul 2011 - R. Yantosca - Added ProTeX headers
!  30 Jul 2012 - R. Yantosca - Now accept am_I_Root as an argument when
!                              running with the traditional driver main.F
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL           :: IS_INIT = .FALSE.
      INTEGER           :: N, AS, IOS
      CHARACTER(LEN=20) :: LINE

      !=================================================================
      ! INIT_PLANEFLIGHT begins here!
      !=================================================================

      ! Find a free file LUN
      IU_FILE = findFreeLUN()

      ! Open file
      OPEN( IU_FILE, FILE=TRIM( INF ), IOSTAT=IOS )
      IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'init_planeflight:1' )

      ! Read four lines of header
      DO N = 1, 4 
         READ( IU_FILE, '(a)', IOSTAT=IOS )
         IF ( IOS /= 0 ) CALL IOERROR( IOS,IU_FILE,'init_planeflight:2')
      ENDDO

      !=================================================================
      ! Read in the number of variables to be output -- store in NPVAR
      !=================================================================
      ! Read in as I4 now for the number of variables we save (skim, 7/24/13)
      READ( IU_FILE, '(i4)', IOSTAT=IOS ) NPVAR
      IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'init_planeflight:3' )

      ! Make sure NPVAR is at least 1
      IF ( NPVAR < 1 ) THEN
         IF ( am_I_Root ) THEN
         WRITE( 6, '(a)') 'NPVAR cannot be zero or negative!'
         WRITE( 6, '(a)') 'STOP in INIT_PLANEFLIGHT (planeflight_mod.f)'
         WRITE( 6, '(a)') REPEAT( '=', 79 )
         ENDIF
         CALL GEOS_CHEM_STOP
      ENDIF

      ! Make sure NPVAR is less than MAXVARS
      IF ( NPVAR > MAXVARS ) THEN
         IF ( am_I_Root ) THEN
         WRITE( 6, '(a)') 'NPVAR exceeds maximum allowed value!'
         WRITE( 6, '(a)') 'STOP in INIT_PLANEFLIGHT (planeflight_mod.f)'
         WRITE( 6, '(a)') REPEAT( '=', 79 )
         ENDIF
         CALL GEOS_CHEM_STOP
      ENDIF

      ! Read in a separation line
      READ( IU_FILE, '(a)', IOSTAT=IOS )
      IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'init_planeflight:4' )

      ! Initialize SMVGEAR reaction counter
      NPREAC = 0

      ! Skip past the species declarations
      DO N = 1, NPVAR
         READ( IU_FILE, '(a)', IOSTAT=IOS ) LINE
         IF ( IOS /= 0 ) CALL IOERROR( IOS,IU_FILE,'init_planeflight:5')

         ! Increment number of SMVGEAR reactions found
         IF ( INDEX( LINE, 'REA_' ) > 0 ) NPREAC = NPREAC + 1
      ENDDO

      ! Read 4 header lines
      DO N = 1, 4
         READ( IU_FILE, '(a)', IOSTAT=IOS )
         IF ( IOS /= 0 ) CALL IOERROR( IOS,IU_FILE,'init_planeflight:6')
      ENDDO
  
      !=================================================================
      ! Read plane track points -- plane, lat/lon/alt, date/time
      !=================================================================
      NPOINTS = 0

      DO 
         
         ! Read a line from the file
         READ( IU_FILE, '(a)', IOSTAT=IOS ) LINE

         ! Exit at end of file
         IF ( IOS < 0 ) EXIT 
         IF ( IOS > 0 ) CALL IOERROR( IOS,IU_FILE,'init_planeflight:7' )
         
         ! Check for END
         IF ( INDEX( LINE, 'END' ) == 0 ) THEN
            NPOINTS = NPOINTS + 1
         ELSE
            EXIT
         ENDIF
      ENDDO

      ! Close file
      CLOSE( IU_FILE )
 
      ! If there are no flight-track points then just return
      IF ( NPOINTS < 1 ) THEN
         DO_PF = .FALSE.
         RETURN
      ENDIF

      ! Make sure NPOINTS is less than MAXPOINTS
      IF ( NPOINTS > MAXPOINTS ) THEN
         IF ( am_I_Root ) THEN
         WRITE( 6, '(a)') 'NPOINTS exceeds maximum allowed value!'
         WRITE( 6, '(a)') 'STOP in INIT_PLANEFLIGHT (planeflight_mod.f)'
         WRITE( 6, '(a)') REPEAT( '=', 79 )
         ENDIF
         CALL GEOS_CHEM_STOP
      ENDIF
         
      !=================================================================
      ! Allocate arrays to maximum sizes
      !
      ! NOTE: To save space, NPREAC is the actual number of reactions
      !       found.  We will worry about this later.  (bmy, 3/25/05)
      !=================================================================
      IF ( .not. IS_INIT ) THEN 

         !-------------------------
         ! Arrays of size NPREAC
         !-------------------------
         ALLOCATE( PREAC( MAX( NPREAC, 1 ) ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'PREAC' )

         ALLOCATE( PRRATE( ITLOOP, MAX( NPREAC, 1 ) ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'PRRATE' )

         !--------------------------
         ! Arrays of size MAXVARS
         !--------------------------
         ALLOCATE( PVAR( MAXVARS ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'PVAR' )

         ALLOCATE( PNAME( MAXVARS ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'PNAMES' )

         !---------------------------
         ! Arrays of size MAXPOINTS
         !---------------------------
         ALLOCATE( PTYPE( MAXPOINTS ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'PTYPE' )

         ALLOCATE( PDATE( MAXPOINTS ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'PDATE' )

         ALLOCATE( PTIME( MAXPOINTS ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'PTIME' )

         ALLOCATE( PTAU( MAXPOINTS ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'PTAU' )
         
         ALLOCATE( PLAT( MAXPOINTS ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'PLAT' )

         ALLOCATE( PLON( MAXPOINTS ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'PLON' )

         ALLOCATE( PPRESS( MAXPOINTS ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'PPRESS' )

         ! Reset IS_INIT flag
         IS_INIT = .TRUE.
      ENDIF

      !=================================================================
      ! Initialize arrays 
      !=================================================================
      PREAC  = 0
      PRRATE = 0e0
      PVAR   = 0
      PNAME  = ''
      PTYPE  = ''
      PDATE  = 0e0
      PTIME  = 0e0
      PTAU   = 0e0
      PLAT   = 0e0
      PLON   = 0e0
      PPRESS = 0e0

      END SUBROUTINE INIT_PLANEFLIGHT
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cleanup_planeflight
!
! !DESCRIPTION: Subroutine CLEANUP\_PLANEFLIGHT deallocates all allocatable 
!  module arrays.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CLEANUP_PLANEFLIGHT
! 
! !REVISION HISTORY: 
!  01 Jul 2001 - M. Evans    - Initial version
!  (1 ) Renamed PRATE to PRRATE to avoid conflict w/ SMVGEAR II (bmy, 4/1/03)
!  29 Jul 2011 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
      IF ( ALLOCATED( PVAR   ) ) DEALLOCATE( PVAR   )
      IF ( ALLOCATED( PREAC  ) ) DEALLOCATE( PREAC  )
      IF ( ALLOCATED( PNAME  ) ) DEALLOCATE( PNAME  )
      IF ( ALLOCATED( PRRATE ) ) DEALLOCATE( PRRATE )
      IF ( ALLOCATED( PTYPE  ) ) DEALLOCATE( PTYPE  ) 
      IF ( ALLOCATED( PDATE  ) ) DEALLOCATE( PDATE  ) 
      IF ( ALLOCATED( PTIME  ) ) DEALLOCATE( PTIME  ) 
      IF ( ALLOCATED( PTAU   ) ) DEALLOCATE( PTAU   ) 
      IF ( ALLOCATED( PLAT   ) ) DEALLOCATE( PLAT   )
      IF ( ALLOCATED( PLON   ) ) DEALLOCATE( PLON   )
      IF ( ALLOCATED( PPRESS ) ) DEALLOCATE( PPRESS )

      END SUBROUTINE CLEANUP_PLANEFLIGHT
!EOC
      END MODULE PLANEFLIGHT_MOD
