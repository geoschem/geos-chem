! $Id: input_mod.f,v 1.11 2005/02/10 19:53:26 bmy Exp $
      MODULE INPUT_MOD
!
!******************************************************************************
!  Module INPUT_MOD reads the GEOS_CHEM input file at the start of the run
!  and passes the information to several other GEOS-CHEM F90 modules.
!  (bmy, 7/20/04, 2/3/05)
! 
!  Module Variables:
!  ============================================================================
!  (1 ) VERBOSE   (LOGICAL )  : Turns on echo-back of lines read from disk.
!  (2 ) FIRSTCOL  (INTEGER )  : First column of the input file (default=26)
!  (3 ) MAXDIM    (INTEGER )  : Maximum number of substrings to read in
!  (4 ) TS_CHEM   (INTEGER )  : Placeholder for chemistry  timestep [min]
!  (5 ) TS_DYN    (INTEGER )  : Placeholder for dynamic    timestep [min]
!  (6 ) TS_CONV   (INTEGER )  : Placeholder for convection timestep [min]
!  (7 ) TS_EMIS   (INTEGER )  : Placeholder for emissions  timestep [min]
!  (8 ) TS_UNIT   (INTEGER )  : Placeholder for unit conv  timestep [min]
!  (9 ) FILENAME  (CHAR*255)  : GEOS-CHEM input file name
!  (10) TOPTITLE  (CHAR*255)  : Top line of input file
!
!  Module Routines:
!  ============================================================================
!  (1 ) READ_INPUT_FILE       : Driver routine for reading GEOS-CHEM input file
!  (2 ) READ_ONE_LINE         : Reads one line at a time
!  (3 ) SPLIT_ONE_LINE        : Splits one line into substrings (by spaces)
!  (4 ) READ_SIMULATION_MENU  : Reads the GEOS-CHEM simulation menu
!  (4 ) READ_TRACER_MENU      : Reads the GEOS-CHEM tracer menu
!  (6 ) READ_AEROSOL_MENU     : Reads the GEOS-CHEM aerosol menu
!  (7 ) READ_EMISSIONS_MENU   : Reads the GEOS-CHEM emission menu
!  (8 ) READ_CHEMISTRY_MENU   : Reads the GEOS-CHEM chemistry menu
!  (9 ) READ_TRANSPORT_MENU   : Reads the GEOS-CHEM transport menu
!  (10) READ_CONVECTION_MENU  : Reads the GEOS-CHEM convection menu
!  (11) READ_DEPOSITION_MENU  : Reads the GEOS-CHEM deposition menu
!  (12) READ_OUTPUT_MENU      : Reads the GEOS-CHEM output menu
!  (13) READ_DIAGNOSTIC_MENU  : Reads the GEOS-CHEM diagnostic menu
!  (14) SET_TINDEX            : Defines which tracers to print to the BPCH file
!  (15) READ_ND49_MENU        : Reads the GEOS-CHEM ND49 timeseries menu
!  (16) READ_ND50_MENU        : Reads the GEOS-CHEM ND50 timeseries menu
!  (17) READ_ND51_MENU        : Reads the GEOS-CHEM ND51 timeseries menu
!  (18) READ_PROD_LOSS_MENU   : Reads the GEOS-CHEM ND65 timeseries menu
!  (19) READ_UNIX_CMDS_MENU   : Reads the GEOS-CHEM unix commands menu
!  (20) READ_NESTED_GRID_MENU : Reads the GEOS-CHEM nested grid menu
!  (21) READ_ARCHIVED_OH_MENU : Reads the GEOS-CHEM archived OH menu
!  (22) READ_O3PL_MENU        : Reads the GEOS-CHEM O3 P/L menu
!  (23) READ_BENCHMARK_MENU   : Reads the GEOS-CHEM benchmark cmds menu
!  (24) VALIDATE_DIRECTORIES  : Makes sure all given directories are valid
!  (25) CHECK_DIRECTORY       : Checks a single directory for errors
!  (26) CHECK_TIME_STEPS      : Sets the GEOS_CHEM timesteps
!  (27) IS_LAST_DAY_GOOD      : Makes sure we have output on last day of run
!  (28) INIT_INPUT            : Initializes directory & logical variables
!
!  GEOS-CHEM modules referenced by "input_mod.f"
!  ============================================================================
!  (1 ) biofuel_mod.f         : Module w/ routines to read biofuel emissions
!  (2 ) biomass_mod.f         : Module w/ routines to read biomass emissions
!  (3 ) bpch2_mod.f           : Module w/ routines for binary punch file I/O
!  (4 ) charpak_mod.f         : Module w/ string handling routines
!  (5 ) dao_mod.f             : Module w/ arrays for DAO met fields
!  (6 ) diag_mod.f            : Module w/ GEOS-CHEM diagnostic arrays
!  (7 ) diag49_mod.f          : Module w/ routines for inst timeseries
!  (8 ) diag50_mod.f          : Module w/ routines for 24hr avg timeseries
!  (9 ) diag51_mod.f          : Module w/ routines for morning/aft t-series
!  (10) diag_oh_mod.f         : Module w/ arrays & routines for mean OH diag
!  (11) diag_pl_mod.f         : Module w/ routines for prod & loss diag's
!  (12) directory_mod.f       : Module w/ GEOS-CHEM data & met field dirs
!  (13) drydep_mod.f          : Module w/ GEOS-CHEM drydep routines
!  (14) error_mod.f           : Module w/ I/O error and NaN check routines
!  (15) file_mod.f            : Module w/ file unit numbers and error checks
!  (13) grid_mod.f            : Module w/ horizontal grid information
!  (14) logical_mod.f         : Module w/ GEOS-CHEM logical switches
!  (15) ocean_mercury_mod.f   : Module w/ routines for ocean flux of Hg0
!  (16) planeflight_mod.f     : Module w/ routines for flight track diag
!  (17) pressure_mod.f        : Module w/ routines to compute P(I,J,L)
!  (18) restart_mod.f         : Module w/ routines for restart file I/O
!  (19) time_mod.f            : Module w/ routines for computing time & date
!  (20) tpcore_bc_mod.f       : Module w/ routines to read/write TPCORE BC's
!  (21) tracer_mod.f          : Module w/ GEOS-CHEM tracer array STT etc.
!  (22) tracerid_mod.f        : Module w/ pointers to tracers & emissions  
!  (23) transport_mod.f       : Module w/ driver routine for TPCORE 
!  (24) unix_cmds_mod.f       : Module w/ Unix commands for unzipping etc
!  (25) upbdflx_mod.f         : Module w/ routines for strat O3, NOy BC's
!  (26) wetscav_mod.f         : Module w/ routines for wetdep/scavenging
!
!  NOTES:
!  (1 ) Now references LSOA in READ_AEROSOL_MENU (bmy, 9/28/04)
!  (2 ) Fixed error checks and assign LSPLIT for tagged Hg.  Also now 
!        refernces LAVHRRLAI from "logical_mod.f" (eck, bmy, 12/20/04)
!  (3 ) Updated for crystalline/aqueous aerosol tracers.  Also moved routine
!        IS_LAST_DAY_GOOD here from "main.f".  Also now references 
!        "ocean_mercury_mod.f".  Also now open the bpch file for output in
!        READ_DIAGNOSTIC_MENU instead of in "main.f".  (cas, sas, bmy, 2/3/05)
!******************************************************************************
!
      IMPLICIT NONE

      !=================================================================
      ! MODULE PRIVATE DECLARATIONS -- keep certain internal variables 
      ! and routines from being seen outside "input_mod.f"
      !=================================================================
     
      ! Make everything PRIVATE ...
      PRIVATE
 
      ! ... except these routines
      PUBLIC :: READ_INPUT_FILE
      
      !=================================================================
      ! MODULE VARIABLES 
      !=================================================================
      LOGICAL            :: VERBOSE  = .FALSE.
      INTEGER, PARAMETER :: FIRSTCOL = 26
      INTEGER, PARAMETER :: MAXDIM   = 255
      INTEGER            :: TS_CHEM
      INTEGER            :: TS_DYN
      INTEGER            :: TS_CONV
      INTEGER            :: TS_EMIS
      INTEGER            :: TS_UNIT
      INTEGER            :: CT1, CT2, CT3
      CHARACTER(LEN=255) :: FILENAME = 'input.geos'
      CHARACTER(LEN=255) :: TOPTITLE
      CHARACTER(LEN=255) :: BPCH_FILE
            
      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement 
      !=================================================================
      CONTAINS

!------------------------------------------------------------------------------

      SUBROUTINE READ_INPUT_FILE
!
!******************************************************************************
!  Subroutine READ_INPUT_FILE is the driver program for reading the GEOS_CHEM
!  input file "input.geos" from disk. (bmy, 7/20/04)
!
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE CHARPAK_MOD, ONLY : STRREPL
      USE FILE_MOD,    ONLY : IU_GEOS, IOERROR
      
      ! Local variables
      LOGICAL            :: EOF
      INTEGER            :: IOS
      CHARACTER(LEN=1)   :: TAB   = ACHAR(9)
      CHARACTER(LEN=1)   :: SPACE = ' '
      CHARACTER(LEN=255) :: LINE

      !=================================================================
      ! READ_INPUT_FILE begins here!
      !=================================================================  

      ! Echo output
      WRITE( 6, '(a  )' ) REPEAT( '=', 79 )
      WRITE( 6, '(a,/)' ) 'G E O S - C H E M   U S E R   I N P U T'
      WRITE( 6, 100   ) TRIM( FILENAME )
 100  FORMAT( 'READ_INPUT_FILE: Reading ', a )

      ! Initialize directory & logical variables
      CALL INIT_INPUT

      ! Open file
      OPEN( IU_GEOS, FILE=TRIM( FILENAME ), STATUS='OLD', IOSTAT=IOS )
      IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_GEOS, 'read_input_file:1' )

      ! Read TOPTITLE for binary punch file
      TOPTITLE = READ_ONE_LINE( EOF  )
      IF ( EOF ) RETURN

      ! Loop until EOF
      DO 
         
         ! Read a line from the file, exit if EOF
         LINE = READ_ONE_LINE( EOF ) 
         IF ( EOF ) EXIT
         
         ! Replace tab characters in LINE (if any) w/ spaces
         CALL STRREPL( LINE, TAB, SPACE )

         !=============================================================
         ! Call individual subroutines to read sections of the file
         ! 
         ! NOTE: You are pretty flexible in setting the order of the
         ! menus in the input file; however, a few guidelines apply:
         !
         ! (1) SIMULATION MENU should be listed first.
         ! (2) TRACER MENU should be listed second.
         ! (3) EMISSIONS, AEROSOL, CHEMISTRY, TRANSPORT, CONVECTION, 
         !      and DEPOSITION menus (in any order) should follow.
         ! (4) Diagnostic menus, including OUTPUT, DIAGNOSTIC,
         !      PLANEFLIGHT, ND48, ND49, ND50, ND51, and PROD_LOSS
         !      menus (in any order) should follow next.
         ! (5) The following menus have no other restriction and
         !      can be placed anywhere (but by convention we will
         !      place them after the diagnostic menu): NESTED GRID
         !      UNIX CMDS, ARCHIVED OH, and O3PL menus.
         !=============================================================
         IF      ( INDEX( LINE, 'SIMULATION MENU'  ) > 0 ) THEN
            CALL READ_SIMULATION_MENU             
                                     
         ELSE IF ( INDEX( LINE, 'TRACER MENU'      ) > 0 ) THEN
            CALL READ_TRACER_MENU                 

         ELSE IF ( INDEX( LINE, 'AEROSOL MENU'     ) > 0 ) THEN
            CALL READ_AEROSOL_MENU              

         ELSE IF ( INDEX( LINE, 'EMISSIONS MENU'   ) > 0 ) THEN
            CALL READ_EMISSIONS_MENU              
                                                  
         ELSE IF ( INDEX( LINE, 'CHEMISTRY MENU'   ) > 0 ) THEN
            CALL READ_CHEMISTRY_MENU              
                                                  
         ELSE IF ( INDEX( LINE, 'TRANSPORT MENU'   ) > 0 ) THEN
            CALL READ_TRANSPORT_MENU              
                                                  
         ELSE IF ( INDEX( LINE, 'CONVECTION MENU'  ) > 0 ) THEN
            CALL READ_CONVECTION_MENU             
                                                  
         ELSE IF ( INDEX( LINE, 'DEPOSITION MENU'  ) > 0 ) THEN
            CALL READ_DEPOSITION_MENU             
                                                  
         ELSE IF ( INDEX( LINE, 'OUTPUT MENU'      ) > 0 ) THEN
            CALL READ_OUTPUT_MENU                 
                                                  
         ELSE IF ( INDEX( LINE, 'DIAGNOSTIC MENU'  ) > 0 ) THEN
            CALL READ_DIAGNOSTIC_MENU             

         ELSE IF ( INDEX( LINE, 'PLANEFLIGHT MENU' ) > 0 ) THEN
            CALL READ_PLANEFLIGHT_MENU             
                                                  
         ELSE IF ( INDEX( LINE, 'ND48 MENU'        ) > 0 ) THEN
            CALL READ_ND48_MENU                   

         ELSE IF ( INDEX( LINE, 'ND49 MENU'        ) > 0 ) THEN
            CALL READ_ND49_MENU                   
                                                  
         ELSE IF ( INDEX( LINE, 'ND50 MENU'        ) > 0 ) THEN
            CALL READ_ND50_MENU                   
                                                  
         ELSE IF ( INDEX( LINE, 'ND51 MENU'        ) > 0 ) THEN
            CALL READ_ND51_MENU                   
                                                  
         ELSE IF ( INDEX( LINE, 'PROD & LOSS MENU' ) > 0 ) THEN
            CALL READ_PROD_LOSS_MENU               
                                                  
         ELSE IF ( INDEX( LINE, 'UNIX CMDS MENU'   ) > 0 ) THEN 
            CALL READ_UNIX_CMDS_MENU              

         ELSE IF ( INDEX( LINE, 'NESTED GRID MENU' ) > 0 ) THEN 
            CALL READ_NESTED_GRID_MENU

         ELSE IF ( INDEX( LINE, 'ARCHIVED OH MENU' ) > 0 ) THEN 
            CALL READ_ARCHIVED_OH_MENU

         ELSE IF ( INDEX( LINE, 'O3 P/L MENU'      ) > 0 ) THEN 
            CALL READ_O3PL_MENU

         ELSE IF ( INDEX( LINE, 'BENCHMARK MENU'   ) > 0 ) THEN 
            CALL READ_BENCHMARK_MENU
                                                  
         ELSE IF ( INDEX( LINE, 'END OF FILE'      ) > 0 ) THEN 
            EXIT

         ENDIF  
      ENDDO

      ! Close input file
      CLOSE( IU_GEOS )

      !=================================================================
      ! Further error-checking and initialization
      !=================================================================

      ! Make sure all directories are valid
      CALL VALIDATE_DIRECTORIES

      ! Check GEOS-CHEM timesteps
      CALL CHECK_TIME_STEPS

      ! Echo output
      WRITE( 6, '(a)' ) REPEAT( '=', 79 )

      ! Return to calling program
      END SUBROUTINE READ_INPUT_FILE

!------------------------------------------------------------------------------

      FUNCTION READ_ONE_LINE( EOF, LOCATION ) RESULT( LINE )
!
!******************************************************************************
!  Subroutine READ_ONE_LINE reads a line from the input file.  If the global 
!  variable VERBOSE is set, the line will be printed to stdout.  READ_ONE_LINE
!  can trap an unexpected EOF if LOCATION is passed.  Otherwise, it will pass
!  a logical flag back to the calling routine, where the error trapping will
!  be done. (bmy, 7/20/04, 7/20/04
! 
!  Arguments as Output:
!  ===========================================================================
!  (1 ) EOF      (CHARACTER) : Logical flag denoting EOF condition
!  (2 ) LOCATION (CHARACTER) : Name of calling routine; traps premature EOF
!
!  Function value:
!  ===========================================================================
!  (1 ) LINE     (CHARACTER) : A line of text as read from the file
!
!  NOTES:
!******************************************************************************
!      
      ! References to F90 modules
      USE FILE_MOD, ONLY : IU_GEOS, IOERROR

      ! Arguments
      LOGICAL,          INTENT(OUT)          :: EOF
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: LOCATION

      ! Local variables
      INTEGER                                :: IOS
      CHARACTER(LEN=255)                     :: LINE, MSG

      !=================================================================
      ! READ_ONE_LINE begins here!
      !=================================================================

      ! Initialize
      EOF = .FALSE.

      ! Read a line from the file
      READ( IU_GEOS, '(a)', IOSTAT=IOS ) LINE

      ! IO Status < 0: EOF condition
      IF ( IOS < 0 ) THEN
         EOF = .TRUE.

         ! Trap unexpected EOF -- stop w/ error msg if LOCATION is passed
         ! Otherwise, return EOF to the calling program
         IF ( PRESENT( LOCATION ) ) THEN
            MSG = 'READ_ONE_LINE: error at: ' // TRIM( LOCATION )
            WRITE( 6, '(a)' ) MSG
            WRITE( 6, '(a)' ) 'Unexpected end of file encountered!'
            WRITE( 6, '(a)' ) 'STOP in READ_ONE_LINE (input_mod.f)'
            WRITE( 6, '(a)' ) REPEAT( '=', 79 )
            STOP
         ELSE
            RETURN
         ENDIF
      ENDIF

      ! IO Status > 0: true I/O error condition
      IF ( IOS > 0 ) CALL IOERROR( IOS, IU_GEOS, 'read_one_line:1' )

      ! Print the line (if necessary)
      IF ( VERBOSE ) WRITE( 6, '(a)' ) TRIM( LINE )

      ! Return to calling program
      END FUNCTION READ_ONE_LINE

!------------------------------------------------------------------------------

      SUBROUTINE SPLIT_ONE_LINE( SUBSTRS, N_SUBSTRS, N_EXP, LOCATION ) 
!
!******************************************************************************
!  Subroutine SPLIT_ONE_LINE reads a line from the input file (via routine 
!  READ_ONE_LINE), and separates it into substrings. (bmy, 7/20/04)
!
!  SPLIT_ONE_LINE also checks to see if the number of substrings found is 
!  equal to the number of substrings that we expected to find.  However, if
!  you don't know a-priori how many substrings to expect a-priori, 
!  you can skip the error check.
! 
!  Arguments as Input:
!  ===========================================================================
!  (3 ) N_EXP     (INTEGER  ) : Number of substrings we expect to find
!                               (N_EXP < 0 will skip the error check!)
!  (4 ) LOCATION  (CHARACTER) : Name of routine that called SPLIT_ONE_LINE
!
!  Arguments as Output:
!  ===========================================================================
!  (1 ) SUBSTRS   (CHARACTER) : Array of substrings (separated by " ")
!  (2 ) N_SUBSTRS (INTEGER  ) : Number of substrings actually found
!
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE CHARPAK_MOD, ONLY: STRSPLIT
      
      ! Arguments
      CHARACTER(LEN=255), INTENT(OUT) :: SUBSTRS(MAXDIM)
      INTEGER,            INTENT(OUT) :: N_SUBSTRS
      INTEGER,            INTENT(IN)  :: N_EXP
      CHARACTER(LEN=*),   INTENT(IN)  :: LOCATION 

      ! Local varaibles
      LOGICAL                         :: EOF
      CHARACTER(LEN=255)              :: LINE, MSG

      !=================================================================
      ! SPLIT_ONE_LINE begins here!
      !=================================================================      

      ! Create error msg
      MSG = 'SPLIT_ONE_LINE: error at ' // TRIM( LOCATION )

      !=================================================================
      ! Read a line from disk
      !=================================================================
      LINE = READ_ONE_LINE( EOF )

      ! STOP on End-of-File w/ error msg
      IF ( EOF ) THEN
         WRITE( 6, '(a)' ) TRIM( MSG )
         WRITE( 6, '(a)' ) 'End of file encountered!' 
         WRITE( 6, '(a)' ) 'STOP in SPLIT_ONE_LINE (input_mod.f)!'
         WRITE( 6, '(a)' ) REPEAT( '=', 79 )
         STOP
      ENDIF

      !=================================================================
      ! Split the lines between spaces -- start at column FIRSTCOL
      !=================================================================
      CALL STRSPLIT( LINE(FIRSTCOL:), ' ', SUBSTRS, N_SUBSTRS )

      ! Sometimes we don't know how many substrings to expect,
      ! if N_EXP is greater than MAXDIM, then skip the error check
      IF ( N_EXP < 0 ) RETURN

      ! Stop if we found the wrong 
      IF ( N_EXP /= N_SUBSTRS ) THEN
         WRITE( 6, '(a)' ) TRIM( MSG )
         WRITE( 6, 100   ) N_EXP, N_SUBSTRS
         WRITE( 6, '(a)' ) 'STOP in SPLIT_ONE_LINE (input_mod.f)!'
         WRITE( 6, '(a)' ) REPEAT( '=', 79 )
         STOP
 100     FORMAT( 'Expected ',i2, ' substrs but found ',i3 )
      ENDIF
       
      ! Return to calling program
      END SUBROUTINE SPLIT_ONE_LINE

!------------------------------------------------------------------------------

      SUBROUTINE READ_SIMULATION_MENU
!
!******************************************************************************
!  Subroutine READ_SIMULATION_MENU reads the SIMULATION MENU section of 
!  the GEOS-CHEM input file (bmy, 7/20/04)
!
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE DIRECTORY_MOD
      USE GRID_MOD,    ONLY : SET_XOFFSET, SET_YOFFSET, COMPUTE_GRID
      USE LOGICAL_MOD, ONLY : LSVGLB,      LUNZIP,      LWAIT
      USE RESTART_MOD, ONLY : SET_RESTART
      USE TIME_MOD,    ONLY : SET_BEGIN_TIME,   SET_END_TIME, 
     &                        SET_CURRENT_TIME, SET_DIAGb,
     &                        SET_NDIAGTIME,    GET_TAU

      ! Local variables
      INTEGER            :: I0,    J0
      INTEGER            :: N,     NDIAGTIME
      INTEGER            :: NYMDb, NHMSb 
      INTEGER            :: NYMDe, NHMSe
      CHARACTER(LEN=255) :: SUBSTRS(MAXDIM)
      CHARACTER(LEN=255) :: IN_RST_FILE
      CHARACTER(LEN=255) :: OUT_RST_FILE

      !=================================================================
      ! READ_SIMULATION_MENU begins here!
      !=================================================================

      ! Start YYYYMMDD, HHMMSS
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 2, 'read_simulation_menu:1' )
      READ( SUBSTRS(1:N), * ) NYMDb, NHMSb

      ! End YYYYMMDD, HHMMSS
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 2, 'read_simulation_menu:2' )
      READ( SUBSTRS(1:N), * ) NYMDe, NHMSe

      ! Run directory
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_simulation_menu:3' )
      READ( SUBSTRS(1:N), '(a)' ) RUN_DIR

      ! Input restart file
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_simulation_menu:4' )
      READ( SUBSTRS(1:N), '(a)' ) IN_RST_FILE

      ! Make new restart file?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_simulation_menu:5' )
      READ( SUBSTRS(1:N), '(a)' ) LSVGLB

      ! Output restart file(s)
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_simulation_menu:6' )
      READ( SUBSTRS(1:N), '(a)' ) OUT_RST_FILE

      ! Root data dir
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_simulation_menu:7' )
      READ( SUBSTRS(1:N), '(a)' ) DATA_DIR

      ! GEOS-1 subdir
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_simulation_menu:8' )
      READ( SUBSTRS(1:N), '(a)' ) GEOS_1_DIR

      ! GEOS-STRAT subdir
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_simulation_menu:9' )
      READ( SUBSTRS(1:N), '(a)' ) GEOS_S_DIR

      ! GEOS-3 subdir
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_simulation_menu:10' )
      READ( SUBSTRS(1:N), '(a)' ) GEOS_3_DIR

      ! GEOS-4 subdir
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_simulation_menu:11' )
      READ( SUBSTRS(1:N), '(a)' ) GEOS_4_DIR

      ! Temp dir
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_simulation_menu:12' )
      READ( SUBSTRS(1:N), '(a)' ) TEMP_DIR

      ! Unzip met fields
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_simulation_menu:13' )
      READ( SUBSTRS(1:N), *     ) LUNZIP

      ! Wait for met fields?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_simulation_menu:14' )
      READ( SUBSTRS(1:N), *     ) LWAIT

      ! I0, J0
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 2, 'read_grid_size_menu:1' )
      READ( SUBSTRS(1:N), * ) I0, J0

      ! Separator line
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_simulation_menu:15' )

      !=================================================================
      ! Print to screen
      !=================================================================
      WRITE( 6, '(/,a)' ) 'SIMULATION MENU'
      WRITE( 6, '(  a)' ) '---------------'
      WRITE( 6, 100     ) 'Start time of run           : ', NYMDb, NHMSb
      WRITE( 6, 100     ) 'End time of run             : ', NYMDe, NHMSe
      WRITE( 6, 110     ) 'Run directory               : ',
     &                     TRIM( RUN_DIR )
      WRITE( 6, 110     ) 'Data Directory              : ',
     &                     TRIM( DATA_DIR )
      WRITE( 6, 110     ) 'GEOS-1     sub-directory    : ', 
     &                     TRIM( GEOS_1_DIR )
      WRITE( 6, 110     ) 'GEOS-STRAT sub-directory    : ', 
     &                     TRIM( GEOS_S_DIR )
      WRITE( 6, 110     ) 'GEOS-3     sub-directory    : ', 
     &                     TRIM( GEOS_3_DIR )
      WRITE( 6, 110     ) 'GEOS-4     sub-directory    : ', 
     &                     TRIM( GEOS_4_DIR )
      WRITE( 6, 110     ) 'Temporary Directory         : ', 
     &                     TRIM( TEMP_DIR )
      WRITE( 6, 110     ) 'Input restart file          : ', 
     &                     TRIM( IN_RST_FILE )
      WRITE( 6, 120     ) 'Create restart file?        : ', LSVGLB
      WRITE( 6, 110     ) 'Output restart file(s)      : ', 
     &                     TRIM( OUT_RST_FILE )
      WRITE( 6, 120     ) 'Unzip met fields?           : ', LUNZIP
      WRITE( 6, 120     ) 'Wait for met fields?        : ', LWAIT
      WRITE( 6, 130     ) 'Global offsets I0, J0       : ', I0, J0

      ! Format statements
 100  FORMAT( A, I8.8, 1X, I6.6 )
 110  FORMAT( A, A              )
 120  FORMAT( A, L5             )
 130  FORMAT( A, 2I5            )

      !=================================================================
      ! Call setup routines from other GEOS-CHEM modules
      !=================================================================

      ! Set start time of run in "time_mod.f"
      CALL SET_BEGIN_TIME( NYMDb, NHMSb )

      ! Set end time of run in "time_mod.f"
      CALL SET_END_TIME( NYMDe, NHMSe )

      ! Set the current time
      CALL SET_CURRENT_TIME()

      ! Set the time of day for writing bpch files
      NDIAGTIME = NHMSe !### test
      CALL SET_NDIAGTIME( NDIAGTIME )

      ! Set the start of the 1st diagnostic interval
      CALL SET_DIAGb( GET_TAU() )
     
      ! Set input & output restart file names
      CALL SET_RESTART( IN_RST_FILE, OUT_RST_FILE )

      ! Set global offsets
      CALL SET_XOFFSET( I0 )
      CALL SET_YOFFSET( J0 )

      ! Compute lat/lon/surface area variables
      CALL COMPUTE_GRID

      ! Set counter
      CT1 = CT1 + 1

      ! Return to calling program
      END SUBROUTINE READ_SIMULATION_MENU

!------------------------------------------------------------------------------

      SUBROUTINE READ_TRACER_MENU
!
!******************************************************************************
!  Subroutine READ_TRACER_MENU reads the TRACER MENU section of the 
!  GEOS-CHEM input file (bmy, 7/20/04, 1/20/05)
!
!  NOTES:
!  (1 ) Now set LSPLIT correctly for Tagged Hg simulation (eck, bmy, 12/13/04)
!  (2 ) Now initialize ocean mercury module (sas, bmy, 1/20/05)
!******************************************************************************
!
      ! References to F90 modules
      USE CHARPAK_MOD,       ONLY : ISDIGIT
      USE BIOMASS_MOD,       ONLY : SET_BIOTRCE
      USE BIOFUEL_MOD,       ONLY : SET_BFTRACE
      USE ERROR_MOD,         ONLY : ALLOC_ERR, ERROR_STOP
      USE LOGICAL_MOD,       ONLY : LSPLIT
      USE OCEAN_MERCURY_MOD, ONLY : INIT_OCEAN_MERCURY
      USE TRACER_MOD
      USE TRACERID_MOD,      ONLY : TRACERID

#     include "CMN_SIZE"  ! Size parameters
#     include "CMN_O3"    ! XNUMOL, XNUMOLAIR

      ! Local variables
      INTEGER            :: N,  T,  C
      CHARACTER(LEN=255) :: C1, C2, LINE, NAME
      CHARACTER(LEN=255) :: SUBSTRS(MAXDIM)

      !=================================================================
      ! READ_TRACER_MENU begins here!
      !=================================================================

      ! NSRCX
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_tracer_menu:1' )
      READ( SUBSTRS(1:N), * ) SIM_TYPE

      ! NTRACE
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_tracer_menu:2' )
      READ( SUBSTRS(1:N), * ) N_TRACERS

      ! Separator line
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_tracer_menu:3' )
      
      ! NTRACE cannot exceed NNPAR
      IF ( N_TRACERS > NNPAR ) THEN
         CALL ERROR_STOP( 'Error: N_TRACERS > NNPAR!', 
     &                    'read_tracer_menu (input_mod.f)' )
      ENDIF
      
      ! Initialize tracer arrays in tracer_mod.f
      CALL INIT_TRACER

      !=================================================================
      ! Read tracer ID, emission ID, biomass ID, biofuel ID, & mol wt
      !=================================================================
      DO T = 1, N_TRACERS

         ! Split line into substrings
         CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_tracer_menu:2' )

         ! Save tracer number
         READ( SUBSTRS(1), * ) ID_TRACER(T)

         ! Save tracer name
         TRACER_NAME(T)  = TRIM( SUBSTRS(2) )

         ! Save tracer molecular wt [g/mole]
         READ( SUBSTRS(3), * ) TRACER_MW_G(T) 

         ! Save tracer molecular wt [kg/mole]
         TRACER_MW_KG(T) = TRACER_MW_G(T) * 1d-3

         ! Ratio of MW air / MW tracer 
         TCVV(T)         = 28.97d0 / TRACER_MW_G(T)

         !---------------------------------
         ! Define quantities from "CMN_O3" 
         ! for backwards compatibility
         !---------------------------------

         ! Molecules air / kg air
         IF ( T == 1 ) THEN
            XNUMOLAIR    = 6.022d+23 / 28.9644d-3
         ENDIF

         ! Molecules tracer / kg tracer
         XNUMOL(T)       = 6.022d+23 / TRACER_MW_KG(T)
         
         ! Molecular weight [kg/mole]
         FMOL(T)         = TRACER_MW_KG(T)

         !----------------------------------------------------
         ! If this is a family tracer, get member species
         !---------------------------------------------------- 
         IF ( N > 3 ) THEN
  
            ! Number of species in this family
            TRACER_N_CONST(T) = N-3
               
            ! Loop over substrings
            DO C = 1, N-3
         
               ! Trim the name
               NAME = TRIM( SUBSTRS(C+3) )
         
               !--------------------------------
               ! Emitted species denoted by ()
               !--------------------------------
               IF ( NAME(1:1) == '(' ) THEN 

                  ! Pull out the text between parentheses
                  NAME = NAME( 2:LEN_TRIM( NAME )-1 )

                  ! Set coeff of emitted tracer to 1 for now
                  TRACER_COEFF(T,C) = 1d0

                  ! Index of which family member is the emitted species
                  ID_EMITTED(T) = C

                  ! 2nd character of name is a number
                  IF ( ISDIGIT( NAME(1:1) ) ) THEN

                     ! Coefficient of emitted tracer
                     READ( NAME(1:1), * ) TRACER_COEFF(T,C)

                     ! Get the rest of the name
                     NAME = NAME( 2:LEN_TRIM( NAME ) )
                  ENDIF

                  ! Tracer constituent name
                  TRACER_CONST(T,C) = NAME

               !--------------------------------  
               ! Non-emitted species 
               !--------------------------------
               ELSE 

                  ! Set coeff of emitted tracer to 1 for now
                  TRACER_COEFF(T,C) = 1d0

                  ! 2nd character of name is a number
                  IF ( ISDIGIT( NAME(1:1) ) ) THEN

                     ! Coefficient of emitted tracer
                     READ( NAME(1:1), * ) TRACER_COEFF(T,C)
                     
                     ! Tracer constituent name
                     NAME = NAME( 2:LEN_TRIM( NAME ) )

                  ENDIF

                  ! Tracer constituent name
                  TRACER_CONST(T,C) = NAME

               ENDIF

            ENDDO

         !----------------------------------------------------
         ! If not a family tracer, there is only one member
         !---------------------------------------------------- 
         ELSE

            ! Number of species in this family
            TRACER_N_CONST(T) = 1

            ! Index of which member is the emitted species
            ID_EMITTED(T)     = 0

            ! Coefficient of each tracer constituent
            TRACER_COEFF(T,1) = 1d0

            ! Tracer constituent name
            TRACER_CONST(T,1) = TRACER_NAME(T)

         ENDIF

      ENDDO

      ! Separator line
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_tracer_menu:5' )

      !=================================================================
      ! Print to screen
      !=================================================================
      WRITE( 6, '(/,a)' ) 
     &       'TRACER MENU (==> denotes SMVGEAR emitted species)' 
      WRITE( 6, '(  a)' ) REPEAT( '-', 48 )
      WRITE( 6, '(  a)' ) '  # Tracer          g/mole'

      ! Print info about each tracer
      DO T = 1, N_TRACERS

         ! Write tracer number, name, & mol wt
         WRITE( 6, 100 ) ID_TRACER(T), TRACER_NAME(T), TRACER_MW_G(T)
              
         ! If a family tracer (or just a tracer w/ emission)
         ! then also print info about species
         IF ( TRACER_N_CONST(T) > 1 .or. ID_EMITTED(T) > 0 ) THEN
            
            ! Loop over member species
            DO C = 1, TRACER_N_CONST(T)

               ! Also flag which is the emitted tracer
               IF ( ID_EMITTED(T) == C ) THEN
                  WRITE( 6, 110 ) TRACER_COEFF(T,C), TRACER_CONST(T,C) 
               ELSE
                  WRITE( 6, 120 ) TRACER_COEFF(T,C), TRACER_CONST(T,C)
               ENDIF
            ENDDO
         ENDIF
      ENDDO

      ! Format statement
 100  FORMAT( I3, 1x, A10, 6x, F6.1 )
 110  FORMAT( 5x, '===> ', f4.1, 1x, A4  )
 120  FORMAT( 5x, '---> ', f4.1, 1x, A4  )

      !=================================================================
      ! Call setup routines from other F90 modules
      !=================================================================

      ! Split into tagged tracers (turn off for full-chemistry)
      IF ( ITS_A_FULLCHEM_SIM() ) THEN

         ! There are no tagged tracers for fullchem
         LSPLIT = .FALSE.

      ELSE IF ( ITS_A_MERCURY_SIM() ) THEN

         ! Need Hg0, Hg2, HgP for tagged Mercury
         LSPLIT = ( N_TRACERS > 3 )

         ! Also initialize ocean flux module
         CALL INIT_OCEAN_MERCURY

      ELSE
         LSPLIT = ( N_TRACERS > 1 )

      ENDIF

      ! Set up tracer flags
      CALL TRACERID

      ! Set NBIOTRCE in "biomass_mod.f"
      CALL SET_BIOTRCE

      ! Set NBFTRACE in "biofuel_mod.f"
      CALL SET_BFTRACE

      ! Set counter
      CT1 = CT1 + 1

      ! Return to calling program
      END SUBROUTINE READ_TRACER_MENU

!------------------------------------------------------------------------------

      SUBROUTINE READ_AEROSOL_MENU
!
!******************************************************************************
!  Subroutine READ_AEROSOL_MENU reads the AEROSOL MENU section of 
!  the GEOS-CHEM input file. (bmy, 7/20/04, 1/14/05)
!
!  NOTES:
!  (1 ) Now reference LSOA (bmy, 9/28/04)
!  (2 ) Now stop run if LSOA=T and SOA tracers are undefined (bmy, 11/19/04)
!  (3 ) Now reference LCRYST from "logical_mod.f".  Also now check to make
!        prevent aerosol tracers from being undefined if the corresponding
!        logical switch is set. (cas, bmy, 1/14/05)
!******************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD,    ONLY : ERROR_STOP
      USE LOGICAL_MOD,  ONLY : LSULF, LCARB, LSOA, 
     &                         LDUST, LDEAD, LSSALT, LCRYST
      USE TRACER_MOD,   ONLY : N_TRACERS, 
     &                         SALA_REDGE_um,      SALC_REDGE_um,
     &                         ITS_AN_AEROSOL_SIM, ITS_A_FULLCHEM_SIM
      USE TRACERID_MOD

      ! Local variables
      INTEGER            :: N, T, I
      CHARACTER(LEN=255) :: SUBSTRS(MAXDIM), MSG, LOCATION

      !=================================================================
      ! READ_AEROSOL_MENU begins here!
      !=================================================================

      ! Location string for ERROR_STOP
      LOCATION = 'READ_AEROSOL_MENU ("input_mod.f")'

      ! Error check
      IF ( CT1 /= 2 ) THEN 
         MSG = 'SIMULATION MENU & TRACER MENU must be read in first!'
         CALL ERROR_STOP( MSG, LOCATION )
      ENDIF

      ! Use online sulfate aerosols 
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_aerosol_menu:1' )
      READ( SUBSTRS(1:N), * ) LSULF

      ! Use crystalline sulfate aerosols  
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_aerosol_menu:2' )
      READ( SUBSTRS(1:N), * ) LCRYST

      ! Use online carbon aerosols 
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_aerosol_menu:3' )
      READ( SUBSTRS(1:N), * ) LCARB

      ! Use secondary organic aerosols?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_aerosol_menu:4' )
      READ( SUBSTRS(1:N), * ) LSOA

      ! Use online dust aerosols 
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_aerosol_menu:5' )
      READ( SUBSTRS(1:N), * ) LDUST

      ! Use DEAD dust mobilization (=T) or GINOUX (=F)
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_aerosol_menu:6' )
      READ( SUBSTRS(1:N), * ) LDEAD      

      ! Use online sea-salt aerosols?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_aerosol_menu:7' )
      READ( SUBSTRS(1:N), * ) LSSALT      

      ! Accum mode seasalt radii bin edges [um]
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 2, 'read_aerosol_menu:8' )
      DO T = 1, N
         READ( SUBSTRS(T), * ) SALA_REDGE_um(T)
      ENDDO

      ! Coarse mode seasalt radii bin edges [um]
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 2, 'read_aerosol_menu:9' )
      DO T = 1, N
         READ( SUBSTRS(T), * ) SALC_REDGE_um(T)
      ENDDO

      ! Separator line
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_aerosol_menu:10' )

      !=================================================================
      ! Error checks
      !=================================================================

      ! Make sure that SALA, SALC bins are contiguous
      IF ( SALA_REDGE_um(2) /= SALC_REDGE_um(1) ) THEN
         MSG = 'SALA and SALC bin edges are not contiguous!'
         CALL ERROR_STOP( MSG, LOCATION )
      ENDIF

      ! Turn off switches for simulations that don't use aerosols
      IF ( ( .not. ITS_A_FULLCHEM_SIM() )  .and. 
     &     ( .not. ITS_AN_AEROSOL_SIM() ) ) THEN
         LSULF  = .FALSE.
         LCRYST = .FALSE.
         LCARB  = .FALSE.
         LSOA   = .FALSE.
         LDUST  = .FALSE.
         LDEAD  = .FALSE.
         LSSALT = .FALSE.
      ENDIF

      !---------------------------------
      ! Error check SULFUR AEROSOLS
      !---------------------------------
      I = IDTDMS + IDTSO2 + IDTSO4 + IDTMSA + IDTNH3 + IDTNH4 + IDTNIT

      IF ( LSULF ) THEN
         IF ( I == 0 ) THEN
            MSG = 'LSULF=T but ONLINE SULFUR AEROSOLS are undefined!'
            CALL ERROR_STOP( MSG, LOCATION )
         ENDIF
      ELSE
         IF ( I > 0 ) THEN
            MSG = 'Cannot use ONLINE SULFUR AEROSOLS if LSULF=F!'
            CALL ERROR_STOP( MSG, LOCATION )
         ENDIF
      ENDIF

      !---------------------------------
      ! Error check CRYST /AQ AEROSOLS
      !---------------------------------
      I = IDTAS + IDTAHS + IDTLET + IDTNH4aq + IDTSO4aq

      IF ( LCRYST ) THEN
         IF ( I == 0 ) THEN
            MSG = 'LCRYST=T but ONLINE CRYST/AQ AEROSOLS are undefined!'
            CALL ERROR_STOP( MSG, LOCATION )
         ENDIF
      ELSE
         IF ( I > 0 ) THEN
            MSG = 'Cannot use ONLINE CRYST/AQ AEROSOLS if LCRYST=F!'
            CALL ERROR_STOP( MSG, LOCATION )
         ENDIF
      ENDIF

      !---------------------------------
      ! Error check CARBON AEROSOLS
      !---------------------------------
      I = IDTBCPO + IDTBCPI + IDTOCPO + IDTOCPI

      IF ( LCARB ) THEN
         IF ( I == 0 ) THEN
            MSG = 'LCARB=T but ONLINE CARBON AEROSOLS are undefined!'
            CALL ERROR_STOP( MSG, LOCATION )
         ENDIF
      ELSE
         IF ( I > 0 ) THEN
            MSG = 'Cannot use ONLINE CARBON AEROSOLS if LCARB=F!'
            CALL ERROR_STOP( MSG, LOCATION )
         ENDIF
      ENDIF

      !---------------------------------
      ! Error check 2dy ORG AEROSOLS
      !---------------------------------
      I = IDTALPH + IDTLIMO + IDTALCO + IDTSOG1 + 
     &    IDTSOG2 + IDTSOG3 + IDTSOA1 + IDTSOA2 + IDTSOA3

      IF ( LSOA ) THEN
         IF ( I == 0 ) THEN
            MSG = 'LSOA=T but ONLINE 2dy ORG AEROSOLS are undefined!'
            CALL ERROR_STOP( MSG, LOCATION )
         ENDIF
      ELSE
         IF ( I > 0 ) THEN
            MSG = 'Cannot use ONLINE 2dy ORG AEROSOLS if LSOA=F!'
            CALL ERROR_STOP( MSG, LOCATION )
         ENDIF
      ENDIF

      !---------------------------------
      ! Error check DUST AEROSOLS
      !---------------------------------
      I = IDTDST1 + IDTDST2 + IDTDST3 + IDTDST4

      IF ( LDUST ) THEN
         IF ( I == 0 ) THEN
            MSG = 'LDUST=T but ONLINE DUST AEROSOLS are undefined!'
            CALL ERROR_STOP( MSG, LOCATION )
         ENDIF
      ELSE
         IF ( I > 0 ) THEN
            MSG = 'Cannot use ONLINE DUST AEROSOLS if LDUST=F!'
            CALL ERROR_STOP( MSG, LOCATION )
         ENDIF
      ENDIF

      !---------------------------------
      ! Error check SEASALT AEROSOLS
      !---------------------------------
      I = IDTSALA + IDTSALC

      IF ( LSSALT ) THEN
         IF ( I == 0 ) THEN
            MSG = 'LSSALT=T but ONLINE SEASALT AEROSOLS are undefined!'
            CALL ERROR_STOP( MSG, LOCATION )
         ENDIF
      ELSE
         IF ( I > 0 ) THEN
            MSG = 'Cannot use ONLINE SEASALT AEROSOLS if LSSALT=F!'
            CALL ERROR_STOP( MSG, LOCATION )
         ENDIF
      ENDIF

      !=================================================================
      ! Print to screen
      !=================================================================
      WRITE( 6, '(/,a)' ) 'AEROSOL MENU'
      WRITE( 6, '(  a)' ) '------------'
      WRITE( 6, 100     ) 'Online SULFATE AEROSOLS?    : ', LSULF
      WRITE( 6, 100     ) 'Online CRYST & AQ AEROSOLS? : ', LCRYST
      WRITE( 6, 100     ) 'Online CARBON AEROSOLS?     : ', LCARB
      WRITE( 6, 100     ) 'Online 2dy ORGANIC AEROSOLS?: ', LSOA
      WRITE( 6, 100     ) 'Online DUST AEROSOLS?       : ', LDUST
      WRITE( 6, 100     ) 'Use DEAD dust emissions?    : ', LDEAD
      WRITE( 6, 100     ) 'Online SEA SALT AEROSOLS?   : ', LSSALT
      WRITE( 6, 110     ) 'Accum  SEA SALT radii [um]  : ', 
     &                     SALA_REDGE_um(1), SALA_REDGE_um(2)
      WRITE( 6, 110     ) 'Coarse SEA SALT radii [um]  : ',
     &                     SALC_REDGE_um(1), SALC_REDGE_um(2)
 100  FORMAT( A, L5     )
 110  FORMAT( A, f6.2, ' - ', f6.2 )

      ! Return to calling program
      END SUBROUTINE READ_AEROSOL_MENU

!------------------------------------------------------------------------------

      SUBROUTINE READ_EMISSIONS_MENU
!
!******************************************************************************
!  Subroutine READ_EMISSIONS_MENU reads the EMISSIONS MENU section of 
!  the GEOS-CHEM input file. (bmy, 7/20/04, 11/5/04)
!
!  NOTES:
!  (1 ) Now read LNEI99 -- switch for EPA/NEI99 emissions (bmy, 11/5/04)
!  (2 ) Now read LAVHRRLAI-switch for using AVHRR-derived LAI (bmy, 12/20/04)
!******************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD,  ONLY : ERROR_STOP
      USE LOGICAL_MOD
      USE TRACER_MOD, ONLY : ITS_A_FULLCHEM_SIM

#     include "CMN_SIZE"  ! Size parameters
#     include "CMN_O3"    ! FSCALYR

      ! Local variables
      INTEGER            :: N
      CHARACTER(LEN=255) :: SUBSTRS(MAXDIM), MSG

      !=================================================================
      ! READ_EMISSIONS_MENU begins here!
      !=================================================================

      ! Error check
      IF ( CT1 /= 2 ) THEN 
         MSG = 'SIMULATION MENU & TRACER MENU must be read in first!'
         CALL ERROR_STOP( MSG, 'READ_EMISSIONS_MENU ("input_mod.f")' )
      ENDIF

      ! Turn on emissions?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_emissions_menu:1' )
      READ( SUBSTRS(1:N), * ) LEMIS

      ! Emissions timestep
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_emissions_menu:2' )
      READ( SUBSTRS(1:N), * ) TS_EMIS

      ! Include anthropogenic emissions?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_emissions_menu:3' )
      READ( SUBSTRS(1:N), * ) LANTHRO

      ! Scale 1985 to year ____ ?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_emissions_menu:4' )
      READ( SUBSTRS(1:N), * ) FSCALYR

      ! Include biofuel emissions?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_emissions_menu:5' )
      READ( SUBSTRS(1:N), * ) LBIOFUEL

      ! Include biogenic emissions?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_emissions_menu:6' )
      READ( SUBSTRS(1:N), * ) LBIOGENIC

      ! Scale Isoprene to Monoterpenes?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_emissions_menu:7' )
      READ( SUBSTRS(1:N), * ) LMONOT

      ! Include biomass emissions?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_emissions_menu:8' )
      READ( SUBSTRS(1:N), * ) LBIOMASS

      ! Seasonal biomass?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_emissions_menu:9' )
      READ( SUBSTRS(1:N), * ) LBBSEA

      ! Scaled to TOMSAI?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_emissions_menu:10' )
      READ( SUBSTRS(1:N), * ) LTOMSAI

      ! Separator line
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_emissions_menu:11' )

      ! Use aircraft NOx
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_emissions_menu:12' )
      READ( SUBSTRS(1:N), * ) LAIRNOX

      ! Use lightning NOx
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_emissions_menu:13' )
      READ( SUBSTRS(1:N), * ) LLIGHTNOX

      ! Use soil NOx
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_emissions_menu:14' )
      READ( SUBSTRS(1:N), * ) LSOILNOX

      ! Use ship SO2 emissions?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_emissions_menu:14' )
      READ( SUBSTRS(1:N), * ) LSHIPSO2

      ! Use EPA/NEI99 emissions?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_emissions_menu:15' )
      READ( SUBSTRS(1:N), * ) LNEI99

      ! Use AVHRR-derived LAI fields?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_emissions_menu:16' )
      READ( SUBSTRS(1:N), * ) LAVHRRLAI

      ! Separator line
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_emissions_menu:17' )

      !=================================================================
      ! Error check logical flags
      !=================================================================

      ! Define these flags for backwards compatibility
      LFOSSIL = LANTHRO
      LWOODCO = LBIOFUEL
      LBIONOX = LBIOMASS
      
      ! Turn off full-chem only switches 
      IF ( .not. ITS_A_FULLCHEM_SIM() ) THEN
         LLIGHTNOX = .FALSE.
         LAIRNOX   = .FALSE.
         LSOILNOX  = .FALSE.
      ENDIF
      
      !=================================================================
      ! Print to screen
      !=================================================================
      WRITE( 6, '(/,a)' ) 'EMISSIONS MENU'
      WRITE( 6, '(  a)' ) '--------------'
      WRITE( 6, 100     ) 'Turn on emissions?          : ', LEMIS
      WRITE( 6, 110     ) 'Emissions timestep [min]    : ', TS_EMIS
      WRITE( 6, 100     ) 'Turn on ANTHRO emissions    : ', LANTHRO
      WRITE( 6, 110     ) 'ANTHRO scale factor year    : ', FSCALYR
      WRITE( 6, 100     ) 'Turn on BIOFUEL emissions?  : ', LFOSSIL
      WRITE( 6, 100     ) 'Turn on BIOGENIC emissions? : ', LBIOGENIC
      WRITE( 6, 100     ) 'Scale ISOP to MONOTERPENES? : ', LMONOT
      WRITE( 6, 100     ) 'Turn on BIOMASS EMISSIONS   : ', LBIOMASS
      WRITE( 6, 100     ) 'Use seasonal BIOMASS emiss? : ', LBBSEA
      WRITE( 6, 100     ) 'Scale BIOMASS to TOMS-AI?   : ', LTOMSAI
      WRITE( 6, 100     ) 'Turn on LIGHTNING NOx?      : ', LLIGHTNOX
      WRITE( 6, 100     ) 'Turn on AIRCRAFT NOx?       : ', LAIRNOX
      WRITE( 6, 100     ) 'Turn on SOIL NOx?           : ', LSOILNOX
      WRITE( 6, 100     ) 'Turn on SHIP SO2 EMISSIONS? : ', LSHIPSO2
      WRITE( 6, 100     ) 'Turn on EPA/NEI99 EMISSIONS?: ', LNEI99
      WRITE( 6, 100     ) 'Turn on AVHRR-derived LAI?  : ', LAVHRRLAI
      
      ! FORMAT statements
 100  FORMAT( A, L5 )
 110  FORMAT( A, I5 )

      ! Return to calling program
      END SUBROUTINE READ_EMISSIONS_MENU

!------------------------------------------------------------------------------

      SUBROUTINE READ_CHEMISTRY_MENU
!
!******************************************************************************
!  Subroutine READ_CHEMISTRY_MENU reads the CHEMISTRY MENU section of 
!  the GEOS-CHEM input file. (bmy, 7/20/04)
!
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD,   ONLY : ERROR_STOP 
      USE LOGICAL_MOD, ONLY : LCHEM, LEMBED
      USE TIME_MOD,    ONLY : SET_CT_CHEM

#     include "CMN_SIZE"  ! Size parameters
#     include "CMN"       ! IEBD1 etc

      ! Local variables
      INTEGER            :: N
      CHARACTER(LEN=255) :: SUBSTRS(MAXDIM), MSG

      !=================================================================
      ! READ_CHEMISTRY_MENU begins here!
      !=================================================================

      ! Error check
      IF ( CT1 /= 2 ) THEN 
         MSG = 'SIMULATION MENU & TRACER MENU must be read in first!'
         CALL ERROR_STOP( MSG, 'READ_CHEMISTRY_MENU ("input_mod.f")' )
      ENDIF

      ! Turn on chemistry?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_chemistry_menu:1' )
      READ( SUBSTRS(1:N), * ) LCHEM

      ! Chemistry timestep
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_chemistry_menu:2' )
      READ( SUBSTRS(1:N), * ) TS_CHEM

      ! Use embedded chemistry?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_chemistry_menu:3' )
      READ( SUBSTRS(1:N), * ) LEMBED

      ! Embedded chemistry limits?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 4, 'read_chemistry_menu:4' )
      READ( SUBSTRS(1:N), * ) IEBD1, JEBD1, IEBD2, JEBD2

      ! Separator line
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_chemistry_menu:5' )

      !=================================================================
      ! Print to screen
      !=================================================================
      WRITE( 6, '(/,a)' ) 'CHEMISTRY MENU'
      WRITE( 6, '(  a)' ) '--------------'
      WRITE( 6, 100     ) 'Turn on chemistry?          : ', LCHEM
      WRITE( 6, 110     ) 'Chemistry timestep [min]    : ', TS_CHEM
      WRITE( 6, 100     ) 'Turn on EMBEDDED CHEMISTRY? : ', LEMBED
      WRITE( 6, 120     ) 'EMBEDDED CHEM lower L box:  : ', IEBD1, JEBD1
      WRITE( 6, 120     ) 'EMBEDDED CHEM upper R box   : ', IEBD2, JEBD2
      
      ! FORMAT statements
 100  FORMAT( A, L5  )
 110  FORMAT( A, I5  )
 120  FORMAT( A, 2I5 )

      ! Return to calling program
      END SUBROUTINE READ_CHEMISTRY_MENU

!------------------------------------------------------------------------------

      SUBROUTINE READ_TRANSPORT_MENU
!
!******************************************************************************
!  Subroutine READ_TRANSPORT_MENU reads the TRANSPORT MENU section of 
!  the GEOS-CHEM input file. (bmy, 7/20/04, 12/1/04)
!
!  NOTES:
!  (1 ) Now define MAX_DYN for 1 x 1.25 grid (bmy, 12/1/04)
!******************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD,     ONLY : ERROR_STOP
      USE LOGICAL_MOD,   ONLY : LTRAN, LUPBD, LMFCT, LFILL
      USE TRACER_MOD,    ONLY : ITS_A_FULLCHEM_SIM, ITS_A_TAGOX_SIM
      USE TRANSPORT_MOD, ONLY : SET_TRANSPORT
      USE UPBDFLX_MOD,   ONLY : INIT_UPBDFLX

      ! Local variables
      INTEGER            :: N, IORD, JORD, KORD, J1, KS, MAX_DYN
      CHARACTER(LEN=255) :: SUBSTRS(MAXDIM), MSG, LOCATION

      !=================================================================
      ! READ_TRANSPORT_MENU begins here!
      !=================================================================

      ! Location for err msg
      LOCATION = 'READ_TRANSPORT_MENU ("input_mod.f")'

      ! Error check
      IF ( CT1 /= 2 ) THEN 
         MSG = 'SIMULATION MENU & TRACER MENU must be read in first!'
         CALL ERROR_STOP( MSG, LOCATION )
      ENDIF

      ! Turn on transport?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_transport_menu:1' )
      READ( SUBSTRS(1:N), * ) LTRAN

      ! Do flux-corrected transport?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_transport_menu:2' )
      READ( SUBSTRS(1:N), * ) LMFCT

      ! Fill negative values
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_transport_menu:3' )
      READ( SUBSTRS(1:N), * ) LFILL

      ! IORD, JORD, KORD
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 3, 'read_transport_menu:4' )
      READ( SUBSTRS(1:N), * ) IORD, JORD, KORD

      ! Transport timestep
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_transport_menu:5' )
      READ( SUBSTRS(1:N), * ) TS_DYN

      ! Include strat O3/NOy
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_transport_menu:6' )
      READ( SUBSTRS(1:N), * ) LUPBD

      ! Separator line
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_transport_menu:7' )

      !=================================================================
      ! Error check settings
      !=================================================================      

      ! Turn off LUPBD if for simulations other than fullchem & tagox
      IF ( .not. ITS_A_FULLCHEM_SIM()  .and. 
     &     .not. ITS_A_TAGOX_SIM()    ) LUPBD = .FALSE.
      
      ! Define maximum timestep for transport
#if   defined( GRID4x5   ) 
      MAX_DYN = 30
#elif defined( GRID2x25  )
      MAX_DYN = 15
#elif defined( GRID1x125 )
      MAX_DYN = 10
#elif defined( GRID1x1   ) 
      MAX_DYN = 10
#endif

      ! If TS_DYN is greater than MAX_DYN, then stop w/ error
      IF ( TS_DYN > MAX_DYN ) THEN
         MSG = 'Transport timestep is too small!'
         CALL ERROR_STOP( MSG, LOCATION )
      ENDIF

      !=================================================================
      ! Print to screen
      !=================================================================
      WRITE( 6, '(/,a)' ) 'TRANSPORT MENU'
      WRITE( 6, '(  a)' ) '--------------'
      WRITE( 6, 100     ) 'Turn on transport?          : ', LTRAN
      WRITE( 6, 100     ) 'Use TPCORE flux-correction? : ', LMFCT
      WRITE( 6, 100     ) 'Let TPCORE Fill negatives?  : ', LFILL
      WRITE( 6, 110     ) 'IORD, JORD, KORD for TPCORE?: ', IORD, 
     &                                                      JORD, KORD
      WRITE( 6, 120     ) 'Transport timestep [min]    : ', TS_DYN
      WRITE( 6, 100     ) 'Use strat BC for O3 & NOy?  : ', LUPBD

      ! FORMAT statements
 100  FORMAT( A, L5  )
 110  FORMAT( A, 5I5 )
 120  FORMAT( A, I5  )
      
      !=================================================================
      ! Call setup routines from other F90 modules
      !=================================================================

      ! Pass parameters to "transport_mod.f"
      CALL SET_TRANSPORT( IORD, JORD, KORD )

      ! Pass parameters to "upbdflx_mod.f"
      CALL INIT_UPBDFLX( IORD, JORD, KORD )

      ! Return to calling program
      END SUBROUTINE READ_TRANSPORT_MENU

!------------------------------------------------------------------------------

      SUBROUTINE READ_CONVECTION_MENU
!
!******************************************************************************
!  Subroutine READ_CONVECTION_MENU reads the CONVECTION MENU section of 
!  the GEOS-CHEM input file. (bmy, 7/20/04)
!
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD,   ONLY : ERROR_STOP
      USE LOGICAL_MOD, ONLY : LCONV, LTURB

      ! Local variables
      INTEGER            :: N
      CHARACTER(LEN=255) :: SUBSTRS(MAXDIM), MSG

      !=================================================================
      ! READ_CONVECTION_MENU begins here!
      !=================================================================

      ! Error check
      IF ( CT1 /= 2 ) THEN 
         MSG = 'SIMULATION MENU & TRACER MENU must be read in first!'
         CALL ERROR_STOP( MSG, 'READ_TRANSPORT_MENU ("input_mod.f")' )
      ENDIF

      ! Turn on convection?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_convection_menu:1' )
      READ( SUBSTRS(1:N), * ) LCONV

      ! Turn on BL mixing
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_convection_menu:2' )
      READ( SUBSTRS(1:N), * ) LTURB

      ! Convection timestep
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_convection_menu:4' )
      READ( SUBSTRS(1:N), * ) TS_CONV

      ! Separator line
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_convection_menu:5' )

      !=================================================================
      ! Print to screen
      !=================================================================
      WRITE( 6, '(/,a)' ) 'CONVECTION MENU'
      WRITE( 6, '(  a)' ) '----------------'
      WRITE( 6, 100     ) 'Turn on cloud convection?   : ', LCONV
      WRITE( 6, 100     ) 'Turn on PBL mixing?         : ', LTURB
      WRITE( 6, 110     ) 'Convection timestep [min]   : ', TS_CONV

      ! FORMAT statements
 100  FORMAT( A, L5 )
 110  FORMAT( A, I5 )
      
      ! Return to calling program 
      END SUBROUTINE READ_CONVECTION_MENU

!------------------------------------------------------------------------------

      SUBROUTINE READ_DEPOSITION_MENU
!
!******************************************************************************
!  Subroutine READ_DEPOSITION_MENU reads the DEPOSITION MENU section of 
!  the GEOS-CHEM input file. (bmy, 7/20/04, 12/15/04)
!
!  NOTES:
!  (1 ) Now print an informational message for tagged Hg (bmy, 12/15/04)
!******************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD,   ONLY : ERROR_STOP
      USE DRYDEP_MOD,  ONLY : INIT_DRYDEP
      USE LOGICAL_MOD, ONLY : LDRYD, LWETD, LSPLIT
      USE TRACER_MOD
      USE WETSCAV_MOD, ONLY : WETDEPID

      ! Local variables
      INTEGER            :: N
      CHARACTER(LEN=255) :: SUBSTRS(MAXDIM), MSG

      !=================================================================
      ! READ_DEPOSITION_MENU begins here!
      !=================================================================

      ! Error check
      IF ( CT1 /= 2 ) THEN 
         MSG = 'SIMULATION MENU & TRACER MENU must be read in first!'
         CALL ERROR_STOP( MSG, 'READ_DEPOSITION_MENU ("input_mod.f")' )
      ENDIF

      ! Turn on drydep?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_deposition_menu:1' )
      READ( SUBSTRS(1:N), * ) LDRYD

      ! Turn on wetdep?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_deposition_menu:2' )
      READ( SUBSTRS(1:N), * ) LWETD

      ! Separator line
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_deposition_menu:3' )

      !=================================================================
      ! Error check settings
      !=================================================================

      ! Turn off drydep for simulations that don't need it
      IF ( ITS_A_CH3I_SIM()    ) LDRYD = .FALSE.
      IF ( ITS_A_COPARAM_SIM() ) LDRYD = .FALSE.
      IF ( ITS_A_TAGCO_SIM()   ) LDRYD = .FALSE.

      ! Turn off wetdep for simulations that don't need it
      IF ( ITS_A_CH3I_SIM()    ) LWETD = .FALSE.
      IF ( ITS_A_HCN_SIM()     ) LWETD = .FALSE.
      IF ( ITS_A_COPARAM_SIM() ) LWETD = .FALSE.
      IF ( ITS_A_TAGOX_SIM()   ) LWETD = .FALSE.
      IF ( ITS_A_TAGCO_SIM()   ) LWETD = .FALSE.
      IF ( ITS_A_C2H6_SIM()    ) LWETD = .FALSE.
      IF ( ITS_A_CH4_SIM()     ) LWETD = .FALSE.
      
      !=================================================================
      ! Print to screen
      !=================================================================
      WRITE( 6, '(/,a)' ) 'DEPOSITION MENU'
      WRITE( 6, '(  a)' ) '---------------'
      WRITE( 6, 100     ) 'Turn on dry deposition?     : ', LDRYD
      WRITE( 6, 100     ) 'Turn on wet deposition?     : ', LWETD

      ! FORMAT statements
 100  FORMAT( A, L5 )
 110  FORMAT( A, I5 )

      !=================================================================
      ! Call setup routines from other F90 modules
      !=================================================================

      ! Initialize dry deposition arrays
      IF ( LDRYD ) THEN

         ! Setup for dry deposition
         CALL INIT_DRYDEP
      
         ! Print extra info message for Hg simulation
         IF ( ITS_A_MERCURY_SIM() .and. LSPLIT ) THEN
            WRITE ( 6, 120 )
            WRITE ( 6, 121 )
         ENDIF
      ENDIF

      ! Initialize wet deposition tracers
      IF ( LWETD ) CALL WETDEPID

      ! FORMAT strings
 120  FORMAT( /, 'All tagged Hg2 tracers have the same dep velocity '
     &           'as the total Hg2 tracer.' )
 121  FORMAT( 'All tagged HgP tracers have the same dep velocity '
     &        'as the total HgP tracer.' )

      ! Return to calling program
      END SUBROUTINE READ_DEPOSITION_MENU

!------------------------------------------------------------------------------

      SUBROUTINE READ_OUTPUT_MENU
!
!******************************************************************************
!  Subroutine READ_OUTPUT_MENU reads the OUTPUT MENU section of 
!  the GEOS-CHEM input file. (bmy, 7/20/04)
!
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE FILE_MOD, ONLY : IU_GEOS, IOERROR
      
#     include "CMN_SIZE" ! Size parameters
#     include "CMN_DIAG" ! NJDAY

      ! Local variables
      INTEGER :: IOS

      !=================================================================
      ! READ_OUTPUT_MENU begins here!
      !=================================================================
      READ( IU_GEOS, 100, IOSTAT=IOS ) NJDAY
 100  FORMAT( 26x, 31i1, /  26x, 29i1, /, 26x, 31i1, /, 26x, 30i1, /, 
     &        26x, 31i1, /, 26x, 30i1, /, 26x, 31i1, /, 26x, 31i1, /,
     &        26x, 30i1, /  26x, 31i1, /, 26x, 30i1, /, 26x, 31i1 )

      ! Error check
      IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_GEOS, 'read_output_menu:1' )

      !=================================================================
      ! Print to screen
      !=================================================================

      WRITE( 6, '(/,a)' ) 'OUTPUT MENU'
      WRITE( 6, '(  a)' ) '-----------'
      WRITE( 6, 110     )
      WRITE( 6, 120     )
      WRITE( 6, 130     )
      WRITE( 6, 140     ) NJDAY

      ! FORMAT statements
 110  FORMAT( '              1111111111222222222233' )
 120  FORMAT( '     1234567890123456789012345678901' )
 130  FORMAT( '     -------------------------------' )
 140  FORMAT( 'JAN--', 31i1, /, 'FEB--', 29i1, /, 'MAR--', 31i1, /, 
     &        'APR--', 30i1, /, 'MAY--', 31i1, /, 'JUN--', 30i1, /, 
     &        'JUL--', 31i1, /, 'AUG--', 31i1, /, 'SEP--', 30i1, /,
     &        'OCT--', 31i1, /, 'NOV--', 30i1, /, 'DEC--', 31i1 )

      ! Make sure we have output at end of run
      CALL IS_LAST_DAY_GOOD

      ! Return to calling program
      END SUBROUTINE READ_OUTPUT_MENU

!------------------------------------------------------------------------------

      SUBROUTINE READ_DIAGNOSTIC_MENU
!
!******************************************************************************
!  Subroutine READ_DIAGNOSTIC_MENU reads the DIAGNOSTIC MENU section of 
!  the GEOS-CHEM input file. (bmy, 7/20/04, 2/3/05)
!
!  NOTES:
!  (1 ) Now reference IU_BPCH from "file_mod.f" and OPEN_BPCH2_FOR_WRITE
!        from "bpch2_mod.f".  Now opens the bpch file for output here
!        instead of w/in "main.f" (bmy, 2/3/05)
!******************************************************************************
!
      ! References to F90 modules
      USE BIOMASS_MOD,  ONLY : NBIOTRCE
      USE BIOFUEL_MOD,  ONLY : NBFTRACE
      USE BPCH2_MOD,    ONLY : OPEN_BPCH2_FOR_WRITE
      USE DIAG_MOD
      USE DIAG03_MOD,   ONLY : ND03, PD03, INIT_DIAG03
      USE DIAG_OH_MOD,  ONLY : INIT_DIAG_OH
      USE DRYDEP_MOD,   ONLY : NUMDEP
      USE ERROR_MOD,    ONLY : ERROR_STOP
      USE FILE_MOD,     ONLY : IU_BPCH
      USE LOGICAL_MOD
      USE TRACER_MOD
      USE TRACERID_MOD, ONLY : NEMANTHRO
      USE WETSCAV_MOD,  ONLY : GET_WETDEP_NMAX

#     include "CMN_SIZE" ! Size parameters
#     include "CMN_DIAG" 

      ! Local variables
      INTEGER            :: M, N, N_MAX, N_TMP
      CHARACTER(LEN=255) :: SUBSTRS(MAXDIM), MSG, LOCATION

      !=================================================================
      ! READ_DIAGNOSTIC_MENU begins here!
      !=================================================================

      ! Location for ERROR_STOP
      LOCATION = 'READ_DIAGNOSTIC_MENU ("input_mod.f")'

      ! Error check
      IF ( CT1 /= 2 ) THEN 
         MSG = 'SIMULATION MENU & TRACER MENU must be read in first!'
         CALL ERROR_STOP( MSG, LOCATION )
      ENDIF

      ! Get max number of tracers
      N_MAX = MIN( N_TRACERS, NNPAR )

      ! Binary punch file name
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1,  'read_diagnostic_menu:1' )
      READ( SUBSTRS(1:N), '(a)' ) BPCH_FILE

      ! Separator line
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:2' )

      !--------------------------
      ! ND01: Rn-Pb-Be source
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:3' )
      READ( SUBSTRS(1), * ) ND01
      IF ( .not. ITS_A_RnPbBe_SIM() ) ND01 = 0 
      CALL SET_TINDEX( 01, ND01, SUBSTRS(2:N), N-1, N_MAX )
      
      !--------------------------
      ! ND02: Rn-Pb-Be decay
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:4' )
      READ( SUBSTRS(1), * ) ND02
      IF ( .not. ITS_A_RnPbBe_SIM() ) ND02 = 0 
      CALL SET_TINDEX( 02, ND02, SUBSTRS(2:N), N-1, N_MAX )

      !--------------------------
      ! ND03: Hg diagnostics
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:5' )
      READ( SUBSTRS(1), * ) ND03
      IF ( .not. ITS_A_MERCURY_SIM() ) ND03 = 0
      CALL SET_TINDEX( 03, ND03, SUBSTRS(2:N), N-1, PD03 )

      !--------------------------
      ! ND04: Free diagnostic
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:6' )
      READ( SUBSTRS(1), * ) ND04
      CALL SET_TINDEX( 04, ND04, SUBSTRS(2:N), N-1, PD04 )

      !--------------------------
      ! ND05: Sulfate prod/loss
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:7' )
      READ( SUBSTRS(1), * ) ND05
      IF ( .not. LSULF ) ND05 = 0
      CALL SET_TINDEX( 05, ND05, SUBSTRS(2:N), N-1, PD05 )

      !--------------------------
      ! ND06: Dust emissions
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:8' )
      READ( SUBSTRS(1), * ) ND06
      IF ( .not. LDUST ) ND06 = 0
      CALL SET_TINDEX( 06, ND06, SUBSTRS(2:N), N-1, NDSTBIN )

      !--------------------------
      ! ND07: Carbon/SOA source
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:9' )
      READ( SUBSTRS(1), * ) ND07
      IF ( .not. LCARB ) ND07 = 0
      CALL SET_TINDEX( 07, ND07, SUBSTRS(2:N), N-1, PD07 )

      !--------------------------
      ! ND08: Sea salt source
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:10' )
      READ( SUBSTRS(1), * ) ND08
      IF ( .not. LSSALT ) ND08 = 0
      CALL SET_TINDEX( 08, ND08, SUBSTRS(2:N), N-1, PD08 )

      !--------------------------
      ! ND09: Free
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:11' )
      READ( SUBSTRS(1), * ) ND09
      CALL SET_TINDEX( 09, ND09, SUBSTRS(2:N), N-1, PD09 )

      !--------------------------
      ! ND10: Free
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:12' )
      READ( SUBSTRS(1), * ) ND10
      CALL SET_TINDEX( 10, ND10, SUBSTRS(2:N), N-1, PD10 )

      !--------------------------
      ! ND11: Acetone source
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:13' )
      READ( SUBSTRS(1), * ) ND11
      CALL SET_TINDEX( 11, ND11, SUBSTRS(2:N), N-1, PD11 )

      !--------------------------
      ! ND12: PBL distribution
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:14' )
      READ( SUBSTRS(1), * ) ND12
      CALL SET_TINDEX( 12, ND12, SUBSTRS(2:N), N-1, PD12 )

      !--------------------------
      ! ND13: Sulfur sources
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:15' )
      READ( SUBSTRS(1), * ) ND13
      IF ( .not. LSULF ) ND13 = 0
      CALL SET_TINDEX( 13, ND13, SUBSTRS(2:N), N-1, PD13 )

      !--------------------------
      ! ND14: Wet conv up flux
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:16' )
      READ( SUBSTRS(1), * ) ND14
      IF ( .not. LCONV ) ND14 = 0
      CALL SET_TINDEX( 14, ND14, SUBSTRS(2:N), N-1, N_MAX )

      !--------------------------
      ! ND15: Mass change in PBL
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:17' )
      READ( SUBSTRS(1), * ) ND15
      IF ( .not. LTURB ) ND15 = 0
      CALL SET_TINDEX( 15, ND15, SUBSTRS(2:N), N-1, N_MAX )

      !--------------------------
      ! ND16: Precip fractions
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:18' )
      READ( SUBSTRS(1), * ) ND16
      CALL SET_TINDEX( 16, ND16, SUBSTRS(2:N), N-1, N_MAX )

      !--------------------------
      ! ND17: Rainout losses
      !--------------------------
      N_TMP = GET_WETDEP_NMAX()
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:19' )
      READ( SUBSTRS(1), * ) ND17
      IF ( .not. LWETD ) ND17 = 0
      CALL SET_TINDEX( 17, ND17, SUBSTRS(2:N), N-1, N_TMP )

      !--------------------------
      ! ND18: Washout losses
      !--------------------------
      N_TMP = GET_WETDEP_NMAX()
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:20' )
      READ( SUBSTRS(1), * ) ND18
      IF ( .not. LWETD ) ND18 = 0
      CALL SET_TINDEX( 18, ND18, SUBSTRS(2:N), N-1, N_TMP )

      !--------------------------
      ! ND19: Free
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:21' )
      READ( SUBSTRS(1), * ) ND19
      CALL SET_TINDEX( 19, ND19, SUBSTRS(2:N), N-1, PD19 )

      !--------------------------
      ! ND21: Opt depths etc.
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:23' )
      READ( SUBSTRS(1), * ) ND21
      CALL SET_TINDEX( 21, ND21, SUBSTRS(2:N), N-1, PD21 )

      ! Error check 
      IF ( ND21 > 0 .and. SALA_REDGE_um(2) /= 0.5 ) THEN
         MSG = 'Cannot output seasalt AOD''s when radius bin' //
     &         ' is not split at 0.5 um!!'
         CALL ERROR_STOP( MSG, LOCATION )
      ENDIF

      !--------------------------
      ! ND22: J-values
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:24' )
      READ( SUBSTRS(1), * ) ND22
      CALL SET_TINDEX( 22, ND22, SUBSTRS(2:N), N-1, PD22 )

      CALL SPLIT_ONE_LINE( SUBSTRS, N, 2,  'read_diagnostic_menu:25' ) 
      READ( SUBSTRS(1:N), * ) HR1_JV, HR2_JV

      !--------------------------
      ! ND24: E/W transport flux
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:27' )
      READ( SUBSTRS(1), * ) ND24
      CALL SET_TINDEX( 24, ND24, SUBSTRS(2:N), N-1, N_MAX )

      !--------------------------
      ! ND25: N/S transport flux
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:28' )
      READ( SUBSTRS(1), * ) ND25
      CALL SET_TINDEX( 25, ND25, SUBSTRS(2:N), N-1, N_MAX )

      !--------------------------
      ! ND26: U/D transport flux
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:29' )
      READ( SUBSTRS(1), * ) ND26
      CALL SET_TINDEX( 26, ND26, SUBSTRS(2:N), N-1, N_MAX )

      !--------------------------
      ! ND27: Strat fluxes
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:30' )
      READ( SUBSTRS(1), * ) ND27
      CALL SET_TINDEX( 27, ND27, SUBSTRS(2:N), N-1, PD27 )

      !--------------------------
      ! ND28: Biomass emissions
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:31' )
      READ( SUBSTRS(1), * ) ND28
      IF ( .not. LBIOMASS ) ND28 = 0
      CALL SET_TINDEX( 28, ND28, SUBSTRS(2:N), N-1, NBIOTRCE )

      !--------------------------
      ! ND29: CO sources
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:32' )
      READ( SUBSTRS(1), * ) ND29
      CALL SET_TINDEX( 29, ND29, SUBSTRS(2:N), N-1, PD29 )

      !--------------------------
      ! ND30: Land map
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:33' )
      READ( SUBSTRS(1), * ) ND30
      CALL SET_TINDEX( 30, ND30, SUBSTRS(2:N), N-1, PD30 )

      !--------------------------
      ! ND31: Surface pressure
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:34' )
      READ( SUBSTRS(1), * ) ND31
      CALL SET_TINDEX( 31, ND31, SUBSTRS(2:N), N-1, PD31 )

      !--------------------------
      ! ND32: NOx sources
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:34' )
      READ( SUBSTRS(1), * ) ND32
      CALL SET_TINDEX( 32, ND32, SUBSTRS(2:N), N-1, PD32 )

      !--------------------------
      ! ND33: Column tracer
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:35' )
      READ( SUBSTRS(1), * ) ND33
      CALL SET_TINDEX( 33, ND33, SUBSTRS(2:N), N-1, N_MAX )

      !--------------------------
      ! ND34: Biofuel sources
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:36' )
      READ( SUBSTRS(1), * ) ND34
      IF ( .not. LBIOFUEL ) ND34 = 0
      CALL SET_TINDEX( 34, ND34, SUBSTRS(2:N), N-1, NBFTRACE )

      !--------------------------
      ! ND35: 500 hPa tracer
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:37' )
      READ( SUBSTRS(1), * ) ND35
      CALL SET_TINDEX( 35, ND35, SUBSTRS(2:N), N-1, N_MAX )

      !--------------------------
      ! ND36: Anthro emissions
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:38' )
      READ( SUBSTRS(1), * ) ND36
      CALL SET_TINDEX( 36, ND36, SUBSTRS(2:N), N-1, NEMANTHRO )

      !--------------------------
      ! ND37: Updraft scav frac
      !--------------------------
      N_TMP = GET_WETDEP_NMAX()
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:39' )
      READ( SUBSTRS(1), * ) ND37
      CALL SET_TINDEX( 37, ND37, SUBSTRS(2:N), N-1, N_TMP )

      !--------------------------
      ! ND38: Wet scav losses
      !--------------------------
      N_TMP = GET_WETDEP_NMAX()
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:40' )
      READ( SUBSTRS(1), * ) ND38
      IF ( .not. LWETD ) ND38 = 0
      CALL SET_TINDEX( 38, ND38, SUBSTRS(2:N), N-1, N_TMP )

      !--------------------------
      ! ND39: Wet scav losses
      !--------------------------
      N_TMP = GET_WETDEP_NMAX()
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:41' )
      READ( SUBSTRS(1), * ) ND39
      IF ( .not. LWETD ) ND39 = 0
      CALL SET_TINDEX( 39, ND39, SUBSTRS(2:N), N-1, N_TMP )

      !--------------------------
      ! ND41: Afternoon PBL ht
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:43' )
      READ( SUBSTRS(1), * ) ND41
      CALL SET_TINDEX( 41, ND41, SUBSTRS(2:N), N-1, PD41 )

      !--------------------------
      ! ND42: Free
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:44' )
      READ( SUBSTRS(1), * ) ND42
      CALL SET_TINDEX( 42, ND42, SUBSTRS(2:N), N-1, PD42 )

      !--------------------------
      ! ND43: OH, NO, etc.
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:45' )
      READ( SUBSTRS(1), * ) ND43
      IF ( .not. ITS_A_FULLCHEM_SIM() ) ND43 = 0
      CALL SET_TINDEX( 43, ND43, SUBSTRS(2:N), N-1, PD43 )

      CALL SPLIT_ONE_LINE( SUBSTRS, N, 2,  'read_diagnostic_menu:46' )
      READ( SUBSTRS(1:N), * ) HR1_OH, HR2_OH

      CALL SPLIT_ONE_LINE( SUBSTRS, N, 2,  'read_diagnostic_menu:47' )
      READ( SUBSTRS(1:N), * ) HR1_NO, HR2_NO

      !--------------------------
      ! ND44 drydep vel & flux
      !--------------------------

      ! Number of tracers depends on simulation type
      IF ( ITS_A_TAGOX_SIM() ) THEN
         N_TMP = N_TRACERS 
      ELSE
         N_TMP = NUMDEP
      ENDIF

      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:48' )
      READ( SUBSTRS(1), * ) ND44
      IF ( .not. LDRYD ) ND44 = 0
      CALL SET_TINDEX( 44, ND44, SUBSTRS(2:N), N-1, N_TMP )

      !--------------------------
      ! ND45: Tracer conc.
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:49' )
      READ( SUBSTRS(1), * ) ND45
      CALL SET_TINDEX( 45, ND45, SUBSTRS(2:N), N-1, N_TRACERS )

      CALL SPLIT_ONE_LINE( SUBSTRS, N, 2,  'read_diagnostic_menu:25' ) 
      READ( SUBSTRS(1:N), * ) HR1_OTH, HR2_OTH

      !--------------------------
      ! ND46: Biogenic sources
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:50' )
      READ( SUBSTRS(1), * ) ND46
      CALL SET_TINDEX( 46, ND46, SUBSTRS(2:N), N-1, PD46 )

      !--------------------------
      ! ND47: 24h avg tracer
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:51' )
      READ( SUBSTRS(1), * ) ND47
      CALL SET_TINDEX( 47, ND47, SUBSTRS(2:N), N-1, N_TRACERS )

      !--------------------------
      ! ND52: Free
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:52' )
      READ( SUBSTRS(1), * ) ND52
      CALL SET_TINDEX( 52, ND48, SUBSTRS(2:N), N-1, PD52 )

      !--------------------------
      ! ND53: Free
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:53' )
      READ( SUBSTRS(1), * ) ND53
      CALL SET_TINDEX( 53, ND53, SUBSTRS(2:N), N-1, N_TRACERS )

      !--------------------------
      ! ND54: Free
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:54' )
      READ( SUBSTRS(1), * ) ND54
      CALL SET_TINDEX( 54, ND54, SUBSTRS(2:N), N-1, N_TRACERS )

      !--------------------------
      ! ND55: Tropopause diags
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:55' )
      READ( SUBSTRS(1), * ) ND55
      CALL SET_TINDEX( 55, ND55, SUBSTRS(2:N), N-1, PD55 )

      !--------------------------
      ! ND56: Free
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:56' )
      READ( SUBSTRS(1), * ) ND56
      CALL SET_TINDEX( 56, ND56, SUBSTRS(2:N), N-1, PD56 )

      !--------------------------
      ! ND57: Free
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:57' )
      READ( SUBSTRS(1), * ) ND57
      CALL SET_TINDEX( 57, ND57, SUBSTRS(2:N), N-1, PD57 )

      !--------------------------
      ! ND58: Free
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:58' )
      READ( SUBSTRS(1), * ) ND58
      CALL SET_TINDEX( 58, ND58, SUBSTRS(2:N), N-1, PD58 )

      !--------------------------
      ! ND59: Free
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:59' )
      READ( SUBSTRS(1), * ) ND59
      CALL SET_TINDEX( 59, ND59, SUBSTRS(2:N), N-1, PD59 )

      !--------------------------
      ! ND60: Free
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:60' )
      READ( SUBSTRS(1), * ) ND60
      CALL SET_TINDEX( 60, ND60, SUBSTRS(2:N), N-1, PD60 )

      !--------------------------
      ! ND61: Free
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:61' )
      READ( SUBSTRS(1), * ) ND61
      CALL SET_TINDEX( 61, ND61, SUBSTRS(2:N), N-1, PD61 )

      !--------------------------
      ! ND62: Free
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:62' )
      READ( SUBSTRS(1), * ) ND62
      CALL SET_TINDEX( 62, ND62, SUBSTRS(2:N), N-1, PD62 )

      !--------------------------
      ! ND63: Free
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:63' )
      READ( SUBSTRS(1), * ) ND63
      CALL SET_TINDEX( 63, ND63, SUBSTRS(2:N), N-1, PD63 )

      !--------------------------
      ! ND64: Free
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:64' )
      READ( SUBSTRS(1), * ) ND64
      CALL SET_TINDEX( 64, ND64, SUBSTRS(2:N), N-1, PD64 )

      !--------------------------
      ! ND66: DAO 3-D fields
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:66' )
      READ( SUBSTRS(1), * ) ND66
      CALL SET_TINDEX( 66, ND66, SUBSTRS(2:N), N-1, PD66 )

      !--------------------------
      ! ND67: DAO 2-D fields
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:67' )
      READ( SUBSTRS(1), * ) ND67
      CALL SET_TINDEX( 67, ND67, SUBSTRS(2:N), N-1, PD67 )

      !--------------------------
      ! ND68: Air masses etc
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:68' )
      READ( SUBSTRS(1), * ) ND68
      CALL SET_TINDEX( 68, ND68, SUBSTRS(2:N), N-1, PD68 )

      !--------------------------
      ! ND69: Surface areas
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:69' )
      READ( SUBSTRS(1), * ) ND69
      CALL SET_TINDEX( 69, ND69, SUBSTRS(2:N), N-1, PD69 )

      !--------------------------
      ! ND70: Debug info
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:70' )
      READ( SUBSTRS(1), * ) ND70
      LPRT = ( ND70 > 0 )
      CALL SET_TINDEX( 70, ND70, SUBSTRS(2:N), N-1, PD70 )
     
      ! Separator line
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1,  'read_diagnostic_menu:71' )

      !=================================================================
      ! Call other setup routines
      !================================================================

      ! Allocate diagnostic arrays
      CALL INIT_DIAG03
      CALL NDXX_SETUP

      ! Enable Mean OH (or CH3CCl3) diag for runs which need it
      CALL INIT_DIAG_OH 

      ! Open the binary punch file for output 
      CALL OPEN_BPCH2_FOR_WRITE( IU_BPCH, BPCH_FILE )

      ! Return to calling program
      END SUBROUTINE READ_DIAGNOSTIC_MENU

!------------------------------------------------------------------------------

      SUBROUTINE SET_TINDEX( N_DIAG, L_DIAG, SUBSTRS, N, NMAX )
!
!******************************************************************************
!  Subroutine SET_TINDEX sets the TINDEX and TMAX arrays, which determine how 
!  many tracers to print to the punch file. (bmy, 7/20/04, 11/15/04)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) N_DIAG  (INTEGER  ) : Number of the GEOS-CHEM diagnostic
!  (2 ) L_DIAG  (INTEGER  ) : Number of levels that we are saving
!  (3 ) SUBSTRS (CHARACTER) : Substrings passed from READ_DIAGNOSTIC_MENU
!  (4 ) N       (INTEGER  ) : Number of valid substrings being passed
!  (5 ) NMAX    (INTEGER  ) : Maximum number of tracers for this diagnostic
!      
!  NOTES:
!  (1 ) Bug fix: now do not drop the last tracer number if "all" is not
!        explicitly specified (tmf, bmy, 11/15/04)
!******************************************************************************
!      
#     include "CMN_SIZE"  ! Size parameters
#     include "CMN_DIAG"  ! TMAX, TINDEX

      ! Arguments
      INTEGER,            INTENT(IN) :: N_DIAG, N, NMAX, L_DIAG
      CHARACTER(LEN=255), INTENT(IN) :: SUBSTRS(N)

      ! Local variables
      LOGICAL, SAVE                  :: FIRST = .TRUE.
      LOGICAL                        :: IS_ALL 
      INTEGER                        :: M

      !=================================================================
      ! SET_TINDEX begins here!
      !=================================================================     

      ! Error check
      IF ( N < 1 ) THEN
         WRITE( 6, '(a)' ) 'ERROR: N must be 1 or greater!'
         WRITE( 6, '(a)' ) 'STOP in SET_TINDEX (input_mod.f)'
         WRITE( 6, '(a)' ) REPEAT( '=', 79 )
         STOP
      ENDIF

      !=================================================================
      ! If the word "all" is present, then set TMAX, TINDEX to all
      ! available tracers for the given diagnostic.  Otherwise, just
      ! use the tracers that were read in from the line
      !=================================================================
      IF ( TRIM( SUBSTRS(1) ) == 'all'  .or. 
     &     TRIM( SUBSTRS(1) ) == 'ALL' ) THEN 

         ! TMAX is the max # of tracers to print out
         TMAX(N_DIAG) = NMAX 

         ! Fill TINDEX with all possible diagnostic tracer numbers
         DO M = 1, TMAX(N_DIAG)
            TINDEX(N_DIAG,M) = M
         ENDDO

         ! Set flag
         IS_ALL = .TRUE. 

      ELSE 

         ! Otherwise, set TMAX, TINDEX to the # of tracers
         ! listed in "input.ctm" -- need some error checks too
         TMAX(N_DIAG) = N

         ! Use explicit DO-loop
         DO M = 1, N
            READ( SUBSTRS(M:M), * ) TINDEX(N_DIAG,M)
         ENDDO

         ! Set flag
         IS_ALL = .FALSE.

      ENDIF

      !=================================================================
      ! Print to screen
      !=================================================================

      ! First-time printing only
      IF ( FIRST ) THEN 
         WRITE( 6, '(/,a)' ) 'DIAGNOSTIC MENU'
         WRITE( 6, '(  a)' ) '---------------'
         WRITE( 6, '(  a)' ) 'Diag    L   Tracers being saved to disk'
         FIRST = .FALSE.
      ENDIF

      ! Test if all tracers are being printed out
      IF ( IS_ALL ) THEN

         ! Print abbreviated output string
         IF ( L_DIAG > 0 ) THEN
            WRITE( 6, 100 ) N_DIAG, L_DIAG, 1, TMAX(N_DIAG)
 100        FORMAT( 'ND', i2.2, 2x, i3, 1x, i3, ' -', i3 ) 
         ENDIF

      ELSE

         ! Or just list each tracer
         ! Print each diagnostic and # of tracers that will print out
         IF ( L_DIAG > 0 ) THEN 
            WRITE( 6, 110 ) N_DIAG, L_DIAG, 
     &                      ( TINDEX(N_DIAG,M), M=1,TMAX(N_DIAG) )
 110        FORMAT( 'ND', i2, 2x, i3, 1x, 100i3 ) 
         ENDIF

      ENDIF

      ! Return to calling program
      END SUBROUTINE SET_TINDEX

!------------------------------------------------------------------------------

      SUBROUTINE READ_PLANEFLIGHT_MENU
!
!******************************************************************************
!  Subroutine READ_ND65_MENU reads the ND65 MENU section of the GEOS-CHEM
!  input file (bmy, 7/20/04)
!  
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD,       ONLY : ERROR_STOP
      USE PLANEFLIGHT_MOD, ONLY : SET_PLANEFLIGHT

#     include "CMN_SIZE"  ! MAXFAM
#     include "CMN_DIAG"  ! ND40

      ! Local variables
      LOGICAL            :: DO_PF
      INTEGER            :: N
      CHARACTER(LEN=255) :: IFILE
      CHARACTER(LEN=255) :: OFILE
      CHARACTER(LEN=255) :: SUBSTRS(MAXDIM), MSG

      !=================================================================
      ! READ_PLANEFLIOHT_MENU begins here!
      !=================================================================

      ! Error check
      IF ( CT1 /= 2 ) THEN 
         MSG = 'SIMULATION MENU & TRACER MENU must be read in first!'
         CALL ERROR_STOP( MSG, 'READ_PLANEFLIGHT_MENU ("input_mod.f")' )
      ENDIF

      ! Initialize
      ND40  = 0
      DO_PF = .FALSE.
      IFILE = ''
      OFILE = ''

      ! Turn on planeflight diagnostic?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1,  'read_planeflight_menu:1' )
      READ( SUBSTRS(1:N), * ) DO_PF

      ! Set ND40 flag from DO_PF
      IF ( DO_PF ) ND40 = 1

      ! Input file name (w/ flight track data points)
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1,  'read_planeflight_menu:2' )
      READ( SUBSTRS(1:N), '(a)' ) IFILE

      ! Output file name
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1,  'read_prod_loss_menu:3' )
      READ( SUBSTRS(1:N), '(a)' ) OFILE

      !=================================================================
      ! Print to screen
      !=================================================================
      WRITE( 6, '(/,a)' ) 'PLANEFLIGHT MENU'
      WRITE( 6, '(  a)' ) '----------------'
      WRITE( 6, 100     ) 'Turn on planeflight diag?   : ', DO_PF
      WRITE( 6, 110     ) 'Flight track input file     : ', TRIM(IFILE)
      WRITE( 6, 110     ) 'Output file name            : ', TRIM(OFILE)

      ! FORMAT statements
 100  FORMAT( A, L5    )
 110  FORMAT( A, A     )

      !=================================================================
      ! Call setup routines from other F90 modules
      !=================================================================

      ! Pass variables to "planeflight_mod.f"
      CALL SET_PLANEFLIGHT( DO_PF, IFILE, OFILE )
   
      ! Return to calling program
      END SUBROUTINE READ_PLANEFLIGHT_MENU

!-----------------------------------------------------------------------------

      SUBROUTINE READ_ND48_MENU
!
!******************************************************************************
!  Subroutine READ_ND48_MENU reads the ND48 MENU section of the GEOS-CHEM 
!  input file. (bmy, 7/20/04)
!
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE DIAG48_MOD, ONLY : INIT_DIAG48, ND48_MAX_STATIONS
      USE ERROR_MOD,  ONLY : ERROR_STOP

      ! Local variables
      LOGICAL            :: DO_ND48
      INTEGER            :: N
      INTEGER            :: FREQ
      INTEGER            :: N_STA
      INTEGER            :: IARR(ND48_MAX_STATIONS)
      INTEGER            :: JARR(ND48_MAX_STATIONS)
      INTEGER            :: LARR(ND48_MAX_STATIONS)
      INTEGER            :: NARR(ND48_MAX_STATIONS)
      CHARACTER(LEN=255) :: SUBSTRS(MAXDIM), MSG
      CHARACTER(LEN=255) :: FILE

      !=================================================================
      ! READ_ND48_MENU begins here!
      !=================================================================

      ! Error check
      IF ( CT1 /= 2 ) THEN 
         MSG = 'SIMULATION MENU & TRACER MENU must be read in first!'
         CALL ERROR_STOP( MSG, 'READ_ND48_MENU ("input_mod.f")' )
      ENDIF

      ! Turn on ND48 diagnostic
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_nd49_menu:1' )
      READ( SUBSTRS(1:N), * ) DO_ND48

      ! Timeseries file
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_nd49_menu:1' )
      READ( SUBSTRS(1:N), '(a)' ) FILE

      ! Frequency
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_nd49_menu:2' )
      READ( SUBSTRS(1:N), * ) FREQ

      ! Number of stations 
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_nd49_menu:3' )
      READ( SUBSTRS(1:N), * ) N_STA
      
      ! Read individual stations
      DO N = 1, N_STA
         CALL SPLIT_ONE_LINE( SUBSTRS, N, 4, 'read_nd49_menu:4' )
         READ( SUBSTRS(1), * ) IARR(N) 
         READ( SUBSTRS(2), * ) JARR(N) 
         READ( SUBSTRS(3), * ) LARR(N) 
         READ( SUBSTRS(4), * ) NARR(N) 
      ENDDO

      !=================================================================
      ! Print to screen
      !=================================================================
      WRITE( 6, '(/,a)' ) 'ND48 STATION TIMESERIES MENU'
      WRITE( 6, '(  a)' ) '----------------------------'
      WRITE( 6, 100     ) 'Turn on ND48 timeseries?    : ', DO_ND48
      WRITE( 6, 110     ) 'ND48 timeseries file name   : ', TRIM( FILE )
      WRITE( 6, 120     ) 'ND48 save frequency [min]   : ', FREQ

      DO N = 1, N_STA
         WRITE( 6, 130 ) N, IARR(N), JARR(N), LARR(N), NARR(N) 
      ENDDO

      ! FORMAT statements
 100  FORMAT( A, L5    )
 110  FORMAT( A, A     )
 120  FORMAT( A, I5    )                   
 130  FORMAT( 'ND48 timeseries station', i4, ' : ', 4i5 )

      !=================================================================
      ! Call setup routines from other F90 modules
      !=================================================================

      ! Initialize for ND48 timeseries
      CALL INIT_DIAG48( DO_ND48, FREQ, N_STA, IARR, 
     &                  JARR,    LARR, NARR,  FILE )

      ! Return to calling program
      END SUBROUTINE READ_ND48_MENU

!-----------------------------------------------------------------------------
    
      SUBROUTINE READ_ND49_MENU
!
!******************************************************************************
!  Subroutine READ_ND49_MENU reads the ND49 MENU section of the GEOS-CHEM 
!  input file. (bmy, 7/20/04)
!
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE DIAG49_MOD, ONLY : INIT_DIAG49
      USE ERROR_MOD,  ONLY : ERROR_STOP

#     include "CMN_SIZE"  ! Size parameters

      ! Local variables
      LOGICAL            :: DO_ND49
      INTEGER            :: N,    I,         AS
      INTEGER            :: ND49, N_TRACERS, TRACERS(100)
      INTEGER            :: IMIN, IMAX,      FREQ
      INTEGER            :: JMIN, JMAX,      N_ND49
      INTEGER            :: LMIN, LMAX
      CHARACTER(LEN=255) :: SUBSTRS(MAXDIM), MSG
      CHARACTER(LEN=255) :: FILE

      !=================================================================
      ! READ_ND49_MENU begins here!
      !=================================================================

      ! Error check
      IF ( CT1 /= 2 ) THEN 
         MSG = 'SIMULATION MENU & TRACER MENU must be read in first!'
         CALL ERROR_STOP( MSG, 'READ_ND49_MENU ("input_mod.f")' )
      ENDIF

      ! Initialize
      ND49       = 0
      TRACERS(:) = 0

      ! Turn on ND49 diagnostic
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1,  'read_nd49_menu:1' )
      READ( SUBSTRS(1:N), * ) DO_ND49

      ! Instantaneous 3-D timeseries file
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1,  'read_nd49_menu:2' )
      READ( SUBSTRS(1:N), '(a)' ) FILE

      ! Tracers to include
      CALL SPLIT_ONE_LINE( SUBSTRS, N_ND49, -1, 'read_nd49_menu:3' )
      DO N = 1, N_ND49
         READ( SUBSTRS(N), * ) TRACERS(N)
      ENDDO

      ! FREQ
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1,  'read_nd49_menu:4' )
      READ( SUBSTRS(1:N), * ) FREQ

      ! IMIN, IMAX
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 2,  'read_nd49_menu:5' )
      READ( SUBSTRS(1:N), * ) IMIN, IMAX

      ! JMIN, JMAX
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 2,  'read_nd49_menu:6' )
      READ( SUBSTRS(1:N), * ) JMIN, JMAX

      ! LMIN, LMAX
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 2,  'read_nd49_menu:7' )
      READ( SUBSTRS(1:N), * ) LMIN, LMAX

      ! Separator line
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1,  'read_nd49_menu:8' )

      !=================================================================
      ! Print to screen
      !=================================================================
      WRITE( 6, '(/,a)' ) 'ND49 3-D INSTANTANEOUS TIMESERIES MENU'
      WRITE( 6, '(  a)' ) '--------------------------------------'
      WRITE( 6, 100     ) 'Turn on ND49 timeseries?    : ', DO_ND49
      WRITE( 6, 110     ) 'ND49 timeseries file name   : ', TRIM( FILE )
      WRITE( 6, 120     ) 'ND49 timeseries tracers     : ', 
     &                     ( TRACERS(N), N=1, N_ND49 )
      WRITE( 6, 130     ) 'ND49 save frequency [min]   : ', FREQ
      WRITE( 6, 130     ) 'ND49 longitude limits       : ', IMIN, IMAX
      WRITE( 6, 130     ) 'ND49 latitude  limits       : ', JMIN, JMAX
      WRITE( 6, 130     ) 'ND49 level     limits       : ', LMIN, LMAX

      ! FORMAT statements
 100  FORMAT( A, L5    )
 110  FORMAT( A, A     )
 120  FORMAT( A, 100I3 )
 130  FORMAT( A, 2I5   )

      !=================================================================
      ! Call setup routines from other F90 modules
      !=================================================================

      ! Initialize for ND49 timeseries
      CALL INIT_DIAG49( DO_ND49, N_ND49, TRACERS, IMIN, 
     &                  IMAX,    JMIN,   JMAX,    LMIN,    
     &                  LMAX,    FREQ,   FILE )

      ! Return to calling program
      END SUBROUTINE READ_ND49_MENU

!------------------------------------------------------------------------------

      SUBROUTINE READ_ND50_MENU
!
!******************************************************************************
!  Subroutine READ_ND50_MENU reads the ND50 MENU section of the GEOS-CHEM 
!  input file. (bmy, 7/20/04)
!
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE DIAG50_MOD, ONLY : INIT_DIAG50
      USE ERROR_MOD,  ONLY : ERROR_STOP

      ! Local variables
      LOGICAL            :: DO_ND50
      INTEGER            :: N,      I,       AS  
      INTEGER            :: N_ND50, IMIN,    IMAX, TRACERS(100)   
      INTEGER            :: JMIN,   JMAX,    LMIN, LMAX
      CHARACTER(LEN=255) :: SUBSTRS(MAXDIM), MSG
      CHARACTER(LEN=255) :: FILE

      !=================================================================
      ! READ_ND50_MENU begins here!
      !=================================================================
      
      ! Error check
      IF ( CT1 /= 2 ) THEN 
         MSG = 'SIMULATION MENU & TRACER MENU must be read in first!'
         CALL ERROR_STOP( MSG, 'READ_ND50_MENU ("input_mod.f")' )
      ENDIF

      ! Initialize
      TRACERS(:) = 0

      ! Turn on ND49 diagnostic
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1,  'read_nd50_menu:1' )
      READ( SUBSTRS(1:N), * ) DO_ND50

      ! Instantaneous 3-D timeseries file
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1,  'read_nd50_menu:2' )
      READ( SUBSTRS(1:N), '(a)' ) FILE

      ! Tracers to include
      CALL SPLIT_ONE_LINE( SUBSTRS, N_ND50, -1, 'read_nd50_menu:3' )
      DO N = 1, N_ND50
         READ( SUBSTRS(N), * ) TRACERS(N)
      ENDDO

      ! IMIN, IMAX
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 2,  'read_nd50_menu:4' )
      READ( SUBSTRS(1:N), * ) IMIN, IMAX

      ! JMIN, JMAX
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 2,  'read_nd50_menu:5' )
      READ( SUBSTRS(1:N), * ) JMIN, JMAX

      ! LMIN, LMAX
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 2,  'read_nd50_menu:6' )
      READ( SUBSTRS(1:N), * ) LMIN, LMAX

      ! Separator line
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1,  'read_nd50_menu:8' )

      !=================================================================
      ! Print to screen
      !=================================================================
      WRITE( 6, '(/,a)' ) 'ND50 3-D 24hr AVG TIMESERIES MENU'
      WRITE( 6, '(  a)' ) '---------------------------------'
      WRITE( 6, 100     ) 'Turn on ND50 timeseries?    : ', DO_ND50
      WRITE( 6, 110     ) 'ND50 timeseries file name   : ', TRIM( FILE )
      WRITE( 6, 120     ) 'ND50 timeseries tracers     : ', 
     &                     ( TRACERS(N), N=1, N_ND50 )
      WRITE( 6, 130     ) 'ND50 longitude limits       : ', IMIN, IMAX
      WRITE( 6, 130     ) 'ND50 latitude  limits       : ', JMIN, JMAX
      WRITE( 6, 130     ) 'ND50 level     limits       : ', LMIN, LMAX

      ! FORMAT statements
 100  FORMAT( A, L5    )
 110  FORMAT( A, A     )
 120  FORMAT( A, 100I3 )
 130  FORMAT( A, 2I5   )

      !=================================================================
      ! Call setup routines from other F90 modules
      !=================================================================

      ! Initialize parameters for ND49 timeseries
      CALL INIT_DIAG50( DO_ND50, N_ND50, TRACERS, IMIN, IMAX, 
     &                  JMIN,    JMAX,   LMIN,    LMAX, FILE )

      ! Return to calling program
      END SUBROUTINE READ_ND50_MENU

!------------------------------------------------------------------------------

      SUBROUTINE READ_ND51_MENU
!
!******************************************************************************
!  Subroutine READ_ND51_MENU reads the ND51 MENU section of the GEOS-CHEM 
!  input file. (bmy, 7/20/04)
!
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE DIAG51_MOD, ONLY : INIT_DIAG51
      USE ERROR_MOD,  ONLY : ERROR_STOP

#     include "CMN_SIZE"
#     include "CMN_DIAG"

      ! Local variables
      LOGICAL            :: DO_ND51
      INTEGER            :: N,      I,       AS
      INTEGER            :: N_ND51, FREQ,    TRACERS(100)
      INTEGER            :: IMIN,   IMAX,    JMIN
      INTEGER            :: JMAX,   LMIN,    LMAX
      REAL*8             :: HR1,    HR2,     HR_WRITE
      CHARACTER(LEN=255) :: SUBSTRS(MAXDIM), MSG
      CHARACTER(LEN=255) :: FILE

      !=================================================================
      ! READ_ND51_MENU begins here!
      !=================================================================
      
      ! Error check
      IF ( CT1 /= 2 ) THEN 
         MSG = 'SIMULATION MENU & TRACER MENU must be read in first!'
         CALL ERROR_STOP( MSG, 'READ_ND51_MENU ("input_mod.f")' )
      ENDIF

      ! Initialize
      TRACERS(:) = 0

      ! Turn on ND51 diagnostic
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1,  'read_nd51_menu:1' )
      READ( SUBSTRS(1:N), * ) DO_ND51

      ! Instantaneous 3-D timeseries file
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1,  'read_nd51_menu:2' )
      READ( SUBSTRS(1:N), '(a)' ) FILE

      ! Tracers to include
      CALL SPLIT_ONE_LINE( SUBSTRS, N_ND51, -1, 'read_nd51_menu:3' )
      DO N = 1, N_ND51
         READ( SUBSTRS(N), * ) TRACERS(N)
      ENDDO

      ! NHMS_WRITE
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1,  'read_nd51_menu:4' )
      READ( SUBSTRS(1:N), * ) HR_WRITE

      ! HR1, HR2
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 2,  'read_nd51_menu:5' )
      READ( SUBSTRS(1:N), * ) HR1, HR2

      ! IMIN, IMAX
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 2,  'read_nd51_menu:5' )
      READ( SUBSTRS(1:N), * ) IMIN, IMAX

      ! JMIN, JMAX
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 2,  'read_nd51_menu:6' )
      READ( SUBSTRS(1:N), * ) JMIN, JMAX

      ! LMIN, LMAX
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 2,  'read_nd51_menu:7' )
      READ( SUBSTRS(1:N), * ) LMIN, LMAX

      ! Separator line
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1,  'read_nd51_menu:8' )

      !=================================================================
      ! Print to screen
      !=================================================================
      WRITE( 6, '(/,a)' ) 'ND51 MORNING OR AFTERNOON TIMESERIES MENU'
      WRITE( 6, '(  a)' ) '-----------------------------------------'
      WRITE( 6, 100     ) 'Turn on ND51 timeseries?    : ', DO_ND51
      WRITE( 6, 110     ) 'ND51 timeseries file name   : ', TRIM( FILE )
      WRITE( 6, 120     ) 'ND41 timeseries tracers     : ', 
     &                     ( TRACERS(N), N=1, N_ND51 )
      WRITE( 6, 140     ) 'ND51 hour to write to disk  : ', HR_WRITE
      WRITE( 6, 140     ) 'ND51 averaging period [GMT] : ', HR1,  HR2
      WRITE( 6, 130     ) 'ND51 longitude limits       : ', IMIN, IMAX
      WRITE( 6, 130     ) 'ND51 latitude  limits       : ', JMIN, JMAX
      WRITE( 6, 130     ) 'ND51 altitude  limits       : ', LMIN, LMAX

      ! FORMAT statements
 100  FORMAT( A, L5    )
 110  FORMAT( A, A     )
 120  FORMAT( A, 100I3 )
 130  FORMAT( A, 2I5   )
 140  FORMAT( A, 2F5.1 )

      !=================================================================
      ! Call setup routine from other F90 modules
      !=================================================================

      ! Initialize parameters for ND51 timeseries
      CALL INIT_DIAG51( DO_ND51, N_ND51, TRACERS, HR_WRITE, 
     &                  HR1,     HR2,    IMIN,    IMAX,   
     &                  JMIN,    JMAX,   LMIN,    LMAX,  FILE )

      ! Return to calling program
      END SUBROUTINE READ_ND51_MENU

!------------------------------------------------------------------------------

      SUBROUTINE READ_PROD_LOSS_MENU
!
!******************************************************************************
!  Subroutine READ_PROD_LOSS_MENU reads the PROD & LOSS MENU section of the 
!  GEOS-CHEM input file (bmy, 7/20/04, 10/29/04)
!  
!  NOTES:
!  (1 ) Bug fixes.  Only error check # of prod/loss families for TagOx and 
!        TagCO runs if DO_SAVE_PL=T.  Also turn off this diagnostic for
!        the offline aerosol run. (bmy, 10/29/04)
!******************************************************************************
!
      ! References to F90 modules
      USE CHARPAK_MOD, ONLY : ISDIGIT, STRSPLIT
      USE DIAG_PL_MOD, ONLY : INIT_DIAG_PL
      USE ERROR_MOD,   ONLY : ERROR_STOP
      USE TRACER_MOD,  ONLY : N_TRACERS,       ITS_A_TAGCO_SIM, 
     &                        ITS_A_TAGOX_SIM, ITS_AN_AEROSOL_SIM

#     include "CMN_SIZE"  ! MAXFAM
#     include "CMN_DIAG"  ! ND65

      ! Local variables
      LOGICAL            :: EOF, DO_SAVE_PL, DO_SAVE_O3
      INTEGER, PARAMETER :: MAXMEM=10
      INTEGER            :: F, M, N, NFAM
      INTEGER            :: FAM_NMEM(MAXFAM)
      REAL*8             :: FAM_COEF(MAXMEM,MAXFAM)
      CHARACTER(LEN=14 ) :: FAM_NAME(MAXFAM)
      CHARACTER(LEN=14 ) :: FAM_TYPE(MAXFAM)
      CHARACTER(LEN=14 ) :: FAM_MEMB(MAXMEM,MAXFAM)
      CHARACTER(LEN=255) :: LOCATION,        NAME         
      CHARACTER(LEN=255) :: SUBSTRS(MAXDIM), MSG 

      !=================================================================
      ! READ_PROD_LOSS_MENU begins here!
      !=================================================================

      ! Location string
      LOCATION = 'READ_PROD_LOSS_MENU ("input_mod.f")'

      ! Error check
      IF ( CT1 /= 2 ) THEN 
         MSG = 'SIMULATION MENU & TRACER MENU must be read in first!'
         CALL ERROR_STOP( MSG, LOCATION )
      ENDIF

      ! Initialize
      FAM_NAME(:)   = ''
      FAM_TYPE(:)   = ''
      FAM_NMEM(:)   = 0
      FAM_MEMB(:,:) = ''
      FAM_COEF(:,:) = 0d0

      !=================================================================
      ! Read info about prod & loss families
      !=================================================================

      ! Turn on production & loss diagnostic (e.g. ND65 diagnostic)
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1,  'read_prod_loss_menu:1' )
      READ( SUBSTRS(1:N), * ) DO_SAVE_PL

      ! Read number of levels for ND65 diagnostic 
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1,  'read_prod_loss_menu:2' )
      READ( SUBSTRS(1:N), * ) ND65

      ! Save P(O3) & L(O3) for tagged Ox run? (i.e. ND20 diagnostic)
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1,  'read_prod_loss_menu:3' )
      READ( SUBSTRS(1:N), * ) DO_SAVE_O3

      ! Read number of families 
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1,  'read_prod_loss_menu:4' )
      READ( SUBSTRS(1:N), * ) NFAM

      ! Loop over families
      DO F = 1, NFAM

         ! Get family members
         CALL SPLIT_ONE_LINE( SUBSTRS, N, -1,  'read_prod_loss_menu:5' )

         ! The first entry is the family name (strip colon at end)
         NAME        = SUBSTRS(1)
         FAM_NAME(F) = NAME( 1:LEN_TRIM( NAME )-1 )

         ! Get family type as prod or loss
         IF ( FAM_NAME(F)(1:1) == 'P'   .or. 
     &        FAM_NAME(F)(1:1) == 'p' ) THEN
            FAM_TYPE(F) = 'prod'
         ELSE
            FAM_TYPE(F) = 'loss'
         ENDIF

         ! Number of member species in this prodloss family
         FAM_NMEM(F) = N - 1

         ! Loop over substrings
         DO M = 1, N-1

            ! Family member name
            NAME          =  TRIM( SUBSTRS(M+1) )

            ! Family member coefficient (set to 1 for now)
            FAM_COEF(M,F) = 1d0

            ! If first char is a digit ...
            IF ( ISDIGIT( NAME(1:1) ) ) THEN
               
               ! Save new family coefficient 
               READ( NAME(1:1), * ) FAM_COEF(M,F)

               ! Get the rest of the member name (skip digit)
               NAME = NAME( 2:LEN_TRIM(NAME) )
            ENDIF

            ! Family member name
            FAM_MEMB(M,F) = NAME
         ENDDO
      ENDDO

      ! Separator line
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1,  'read_prod_loss_menu:7' )

      !=================================================================
      ! Error check families for certain types of simulations
      !=================================================================

      ! Tagged Ox
      IF ( DO_SAVE_PL .and. ITS_A_TAGOX_SIM() ) THEN
         IF ( NFAM /= 2*N_TRACERS ) THEN
            MSG = 'Wrong number of P/L families for Tagged Ox!'
            CALL ERROR_STOP( MSG, LOCATION )
         ENDIF
      ENDIF

      ! Tagged CO
      IF ( DO_SAVE_PL .and. ITS_A_TAGCO_SIM() ) THEN
         IF ( NFAM /= N_TRACERS+5 ) THEN
            MSG = 'Wrong number of P/L families for Tagged CO!'
            CALL ERROR_STOP( MSG, LOCATION )
         ENDIF
      ENDIF

      ! Offline aerosol -- turn off DO_SAVE_PL, since we use ND05,
      ! ND06, ND07, ND08, ND13 etc diagnostics instead of ND65
      IF ( ITS_AN_AEROSOL_SIM() ) THEN 
         DO_SAVE_PL    = .FALSE.
         DO_SAVE_O3    = .FALSE.
         ND65          = 0
         NFAM          = 0
         FAM_NAME(:)   = ''
         FAM_TYPE(:)   = ''
         FAM_NMEM(:)   = 0
         FAM_MEMB(:,:) = ''
         FAM_COEF(:,:) = 0d0
      ENDIF

      !=================================================================
      ! Print to screen
      !=================================================================
      WRITE( 6, '(/,a)' ) 'PROD & LOSS DIAGNOSTIC MENU'
      WRITE( 6, '(  a)' ) '---------------------------'      
      WRITE( 6, 100     ) 'Turn on prod & loss diag?   : ', DO_SAVE_PL
      WRITE( 6, 110     ) '# of levels for P/L diag    : ', ND65
      WRITE( 6, 100     ) 'Save P(Ox), L(Ox) for TagOx?: ', DO_SAVE_O3
      
      ! Loop over families
      DO F = 1, NFAM

         ! Write family name, type and # of members
         WRITE( 6, 120 ) FAM_NAME(F), FAM_TYPE(F)
         
         ! Write info about each constituent member
         DO M = 1, FAM_NMEM(F)
            WRITE( 6, 130 ) M, FAM_COEF(M,F), FAM_MEMB(M,F)
         ENDDO
      ENDDO

      ! FORMAT statements
 100  FORMAT( A, L5 )
 110  FORMAT( A, I5 )
 120  FORMAT( /, 'Family=', A5, '  Type=', a4 )
 130  FORMAT( I3, 1X, F4.1, 1X, A5 )

      !=================================================================
      ! Call setup routines from other F90 modules
      !=================================================================

      ! Pass variables to "diag_pl_mod.f" 
      CALL INIT_DIAG_PL( DO_SAVE_PL, DO_SAVE_O3, NFAM,     FAM_NAME,   
     &                   FAM_TYPE,   FAM_NMEM,   FAM_MEMB, FAM_COEF )

      ! Return to calling program
      END SUBROUTINE READ_PROD_LOSS_MENU

!------------------------------------------------------------------------------

      SUBROUTINE READ_UNIX_CMDS_MENU
!
!******************************************************************************
!  Subroutine READ_UNIX_CMDS_MENU reads the UNIX CMDS MENU section of the 
!  GEOS-CHEM input file. (bmy, 7/20/04)
!
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE CHARPAK_MOD, ONLY : STRSQUEEZE
      USE UNIX_CMDS_MOD

      ! Local variables
      LOGICAL            :: EOF
      INTEGER            :: N
      CHARACTER(LEN=255) :: SUBSTRS(MAXDIM)

      !=================================================================
      ! READ_UNIX_CMDS_MENU begins here!
      !=================================================================

      ! Background
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_unix_cmds_menu:1' )
      READ( SUBSTRS(1:N), '(a)' ) BACKGROUND

      ! Redirect
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_unix_cmds_menu:2' )
      READ( SUBSTRS(1:N), '(a)' ) REDIRECT

      ! Remove command
      REMOVE_CMD = READ_ONE_LINE( EOF,    'read_unix_cmds_menu:3' ) 
      REMOVE_CMD = REMOVE_CMD(FIRSTCOL:)
      CALL STRSQUEEZE( REMOVE_CMD )

      ! Separator
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_unix_cmds_menu:4' )
      READ( SUBSTRS(1:N), '(a)' ) SEPARATOR

      ! Wild Card
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_unix_cmds_menu:5' )
      READ( SUBSTRS(1:N), '(a)' ) WILD_CARD

      ! Unzip command
      UNZIP_CMD = READ_ONE_LINE( EOF,     'read_unix_cmds_menu:6' ) 
      UNZIP_CMD = UNZIP_CMD(FIRSTCOL:)
      CALL STRSQUEEZE( UNZIP_CMD )

      ! Zip suffix
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_unix_cmds_menu:7' )
      READ( SUBSTRS(1:N), '(a)' ) ZIP_SUFFIX

      ! Separator line
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_unix_cmds_menu:8' )

      ! Just hardwire the SPACE character
      SPACE = ' '

      !=================================================================
      ! Print to screen
      !=================================================================
      WRITE( 6, '(/,a)' ) 'UNIX CMDS MENU'
      WRITE( 6, '(  a)' ) '---------------'            
      WRITE( 6, 100     ) 'Unix BACKGROUND  command    : ', 
     &                    TRIM( BACKGROUND )
      WRITE( 6, 100     ) 'Unix REDIRECT    command    : ', 
     &                    TRIM( REDIRECT   )
      WRITE( 6, 100     ) 'Unix REMOVE      command    : ',
     &                    TRIM( REMOVE_CMD )
      WRITE( 6, 100     ) 'Unix SEPARATOR   command    : ',
     &                    TRIM( SEPARATOR  )
      WRITE( 6, 100     ) 'Unix WHITE SPACE command    : ',
     &                    TRIM( SPACE      )
      WRITE( 6, 100     ) 'Unix WILD CARD   command    : ',
     &                    TRIM( WILD_CARD  )
      WRITE( 6, 100     ) 'Unix UNZIP       command    : ',
     &                    TRIM( UNZIP_CMD  )
      
      ! FORMAT statements
 100  FORMAT( A, A )

      ! Return to calling program
      END SUBROUTINE READ_UNIX_CMDS_MENU

!------------------------------------------------------------------------------

      SUBROUTINE READ_NESTED_GRID_MENU
!
!******************************************************************************
!  Subroutine READ_NESTED_GRID_MENU reads the NESTED GRID MENU section of 
!  the GEOS-CHEM input file. (bmy, 7/20/04)
!
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE DIRECTORY_MOD, ONLY : TPBC_DIR
      USE LOGICAL_MOD,   ONLY : LWINDO
      USE TPCORE_BC_MOD, ONLY : INIT_TPCORE_BC
 
      ! Local variables
      INTEGER            :: I0W, J0W, I1, I2, J1, J2, N, TS
      CHARACTER(LEN=255) :: SUBSTRS(MAXDIM)

      !=================================================================
      ! READ_NESTED_GRID_MENU begins here!
      !=================================================================

      ! Save out TPCORE BC's at 4x5
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_nested_grid_menu:1' )
      READ( SUBSTRS(1:N), * ) LWINDO

      ! Directory where 4x5 TPCORE BC's are stored
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_nested_grid_menu:2' )
      READ( SUBSTRS(1:N), '(a)' ) TPBC_DIR

      ! Timestep for 4x5 TPCORE BC's
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_nested_grid_menu:3' )
      READ( SUBSTRS(1:N), * ) TS

      ! Lower left box
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 2, 'read_nested_grid_menu:4' )
      READ( SUBSTRS(1:N), * ) I1, J1

      ! Timestep for 4x5 TPCORE BC's
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 2, 'read_nested_grid_menu:5' )
      READ( SUBSTRS(1:N), * ) I2, J2

      ! 1x1 transport region offsets 
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 2, 'read_nested_grid_menu:6' )
      READ( SUBSTRS(1:N), * ) I0W, J0W     

      ! Separator line
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_nested_grid_menu:7' )

      !=================================================================
      ! Print to screen
      !=================================================================
      WRITE( 6, '(/,a)' ) 'NESTED GRID MENU'
      WRITE( 6, '(  a)' ) '----------------'
      WRITE( 6, 100     ) 'Save TPCORE BC''s at 4x5     : ',LWINDO
      WRITE( 6, 110     ) 'Dir w/ archived OH files    : ', 
     &                     TRIM( TPBC_DIR )
      WRITE( 6, 120     ) 'Timestep for 4x5 BC''s [min] : ', TS
      WRITE( 6, 120     ) 'Bot left  box of 4x5 BC area: ', I1,  J1
      WRITE( 6, 120     ) 'Top right box of 4x5 BC area: ', I2,  J2
      WRITE( 6, 120     ) '1 x 1 window offsets        : ', I0W, J0W

      ! FORMAT statements
 100  FORMAT( A, L5  )
 110  FORMAT( A, A   )
 120  FORMAT( A, 2I5 )
    
      !=================================================================
      ! Call setup routines from other F90 modules
      !=================================================================

      ! Pass values to "tpcore_bc_mod.f"
      CALL INIT_TPCORE_BC( TS, I0W, J0W, I1, J1, I2, J2 )

      ! Return to calling program
      END SUBROUTINE READ_NESTED_GRID_MENU

!------------------------------------------------------------------------------

      SUBROUTINE READ_BENCHMARK_MENU
!
!******************************************************************************
!  Subroutine READ_BENCHMARK_MENU reads the BENCHMARK MENU section of 
!  the GEOS-CHEM input file. (bmy, 7/20/04)
!
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE BENCHMARK_MOD, ONLY : INITIAL_FILE, FINAL_FILE
      USE LOGICAL_MOD,   ONLY : LSTDRUN
 
      ! Local variables
      INTEGER            :: N 
      CHARACTER(LEN=255) :: SUBSTRS(MAXDIM)

      !=================================================================
      ! READ_NESTED_GRID_MENU begins here!
      !=================================================================

      ! Save benchmark diagnostic output?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_benchmark_menu:1' )
      READ( SUBSTRS(1:N), * ) LSTDRUN

      ! Filename for initial tracer mass
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_benchmark_menu:2' )
      READ( SUBSTRS(1:N), '(a)' ) INITIAL_FILE

      ! Filename for initial tracer mass
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_benchmark_menu:3' )
      READ( SUBSTRS(1:N), '(a)' ) FINAL_FILE

      ! Separator line
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_benchmark_menu:4' )

      !=================================================================
      ! Print to screen
      !=================================================================
      WRITE( 6, '(/,a)' ) 'BENCHMARK MENU'
      WRITE( 6, '(  a)' ) '--------------'
      WRITE( 6, 100     ) 'Save benchmark diag output? : ', LSTDRUN
      WRITE( 6, 110     ) 'File for initial tracer mass: ',  
     &                     TRIM( INITIAL_FILE )
      WRITE( 6, 110     ) 'File for final tracer mass  : ',  
     &                     TRIM( FINAL_FILE )

      ! FORMAT statements
 100  FORMAT( A, L5  )
 110  FORMAT( A, A   )
    
      ! Return to calling program
      END SUBROUTINE READ_BENCHMARK_MENU

!------------------------------------------------------------------------------

      SUBROUTINE READ_ARCHIVED_OH_MENU
!
!******************************************************************************
!  Subroutine READ_ARCHIVED_OH_MENU reads the ARCHIVED OH MENU section of the 
!  GEOS-CHEM input file. (bmy, 7/20/04)
!
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE DIRECTORY_MOD, ONLY : OH_DIR
 
      ! Local variables
      INTEGER            :: N
      CHARACTER(LEN=255) :: SUBSTRS(MAXDIM)

      !=================================================================
      ! READ_ARCHIVED_OH_MENU begins here!
      !=================================================================

      ! Directory where offline OH files are stored
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_archived_oh_menu:1' )
      READ( SUBSTRS(1:N), '(a)' ) OH_DIR

      ! Separator line
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_archived_oh_menu:1' )

      !=================================================================
      ! Print to screen
      !=================================================================
      WRITE( 6, '(/,a)' ) 'ARCHIVED OH MENU'
      WRITE( 6, '(  a)' ) '-----------------'          
      WRITE( 6, 100     ) 'Dir w/ archived OH files    : ', TRIM(OH_DIR)

      ! FORMAT statements
 100  FORMAT( A, A )
    
      ! Return to calling program
      END SUBROUTINE READ_ARCHIVED_OH_MENU

!------------------------------------------------------------------------------

      SUBROUTINE READ_O3PL_MENU
!
!******************************************************************************
!  Subroutine READ_O3PL_MENU reads the O3 P/L MENU section of the 
!  GEOS-CHEM input file. (bmy, 7/20/04)
!
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE DIRECTORY_MOD, ONLY : O3PL_DIR
 
      ! Local variables
      INTEGER            :: N
      CHARACTER(LEN=255) :: SUBSTRS(MAXDIM)

      !=================================================================
      ! READ_O3PL_MENU begins here!
      !=================================================================

      ! Directory where archived P(O3) and L(O3) files are stored
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_archived_oh_menu:1' )
      READ( SUBSTRS(1:N), '(a)' ) O3PL_DIR

      ! Separator line
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_archived_oh_menu:1' )

      !=================================================================
      ! Print to screen
      !=================================================================
      WRITE( 6, '(/,a)' ) 'ARCHIVED O3PL MENU'
      WRITE( 6, '(  a)' ) '-------------------'          
      WRITE( 6, 100     ) 'Dir w/ archived O3 P/L files: ', 
     &                    TRIM( O3PL_DIR )

      ! FORMAT statements
 100  FORMAT( A, A )
    
      ! Return to calling program
      END SUBROUTINE READ_O3PL_MENU

!------------------------------------------------------------------------------

      SUBROUTINE VALIDATE_DIRECTORIES
!
!******************************************************************************
!  Subroutine VALIDATE_DIRECTORIES makes sure that each of the directories
!  that we have read from the GEOS-CHEM input file are valid.  Also, trailing
!  separator characters will be added. (bmy, 7/20/04)
!
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE DIRECTORY_MOD
      USE GRID_MOD,    ONLY : ITS_A_NESTED_GRID 
      USE LOGICAL_MOD, ONLY : LWINDO,      LUNZIP
      USE TIME_MOD,    ONLY : EXPAND_DATE, GET_NYMDb, GET_NYMDe

#     include "define.h" 

      ! Local variables
      INTEGER            :: NYMDb, NYMDe
      CHARACTER(LEN=255) :: DIR

      !=================================================================
      ! VALIDATE_DIRECTORIES begins here!
      !=================================================================

      ! Get starting & ending dates
      NYMDb = GET_NYMDb()
      NYMDe = GET_NYMDe()

      ! Check directories
      CALL CHECK_DIRECTORY( DATA_DIR )
      CALL CHECK_DIRECTORY( RUN_DIR  )
      CALL CHECK_DIRECTORY( OH_DIR   )
      CALL CHECK_DIRECTORY( O3PL_DIR )

      ! Only validate the TEMP_DIR if we are unzipping met fields
      IF ( LUNZIP ) CALL CHECK_DIRECTORY( TEMP_DIR )

      ! Only validate the TPCORE BC directory if we need it
      IF ( LWINDO .or. ITS_A_NESTED_GRID() ) THEN
         CALL CHECK_DIRECTORY( TPBC_DIR )
      ENDIF

#if   defined( GEOS_1 )

      ! Check GEOS-1 met field directory (starting date)
      DIR = GEOS_1_DIR
      CALL EXPAND_DATE( DIR, NYMDb, 000000 )
      DIR = TRIM( DATA_DIR ) // TRIM( DIR )
      CALL CHECK_DIRECTORY( DIR )

      ! Check GEOS-1 met field directory (ending date)
      DIR = GEOS_1_DIR
      CALL EXPAND_DATE( DIR, NYMDe, 000000 )
      DIR = TRIM( DATA_DIR ) // TRIM( DIR )
      CALL CHECK_DIRECTORY( DIR )

#elif defined( GEOS_STRAT )

      ! Check GEOS-STRAT met field directory (starting date)
      DIR = GEOS_S_DIR
      CALL EXPAND_DATE( DIR, NYMDb, 000000 )
      DIR = TRIM( DATA_DIR ) // TRIM( DIR )
      CALL CHECK_DIRECTORY( DIR )

      ! Check GEOS-STRAT met field directory (ending date)
      DIR = GEOS_S_DIR
      CALL EXPAND_DATE( DIR, NYMDe, 000000 )
      DIR = TRIM( DATA_DIR ) // TRIM( DIR )
      CALL CHECK_DIRECTORY( DIR )

#elif defined( GEOS_3 )

      ! Check GEOS-3 met field directory (starting date)
      DIR = GEOS_3_DIR
      CALL EXPAND_DATE( DIR, NYMDb, 000000 )
      DIR = TRIM( DATA_DIR ) // TRIM( DIR )
      CALL CHECK_DIRECTORY( DIR )

      ! Check GEOS-3 met field directory (ending date)
      DIR = GEOS_3_DIR
      CALL EXPAND_DATE( DIR, NYMDe, 000000 )
      DIR = TRIM( DATA_DIR ) // TRIM( DIR )
      CALL CHECK_DIRECTORY( DIR )

#elif defined( GEOS_4 )

      ! Check GEOS-4 met field directory (starting date)
      DIR = GEOS_4_DIR
      CALL EXPAND_DATE( DIR, NYMDb, 000000 )
      DIR = TRIM( DATA_DIR ) // TRIM( DIR )
      CALL CHECK_DIRECTORY( DIR )


      ! Check GEOS-4 met field directory (ending date)
      DIR = GEOS_4_DIR
      CALL EXPAND_DATE( DIR, NYMDe, 000000 )
      DIR = TRIM( DATA_DIR ) // TRIM( DIR )
      CALL CHECK_DIRECTORY( DIR )

#endif

      ! Return to calling program
      END SUBROUTINE VALIDATE_DIRECTORIES

!------------------------------------------------------------------------------

      SUBROUTINE CHECK_DIRECTORY( DIR )
!
!******************************************************************************
!  Subroutine CHECK_DIRECTORY makes sure that the given directory 
!  is valid.  Also a trailing slash character will be added if necessary. 
!  (bmy, 3/20/03, 7/20/04)
! 
!  Arguments as Input/Output:
!  ============================================================================
!  (1 ) DIR (CHARACTER) : Directory to be checked
!
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules 
      USE ERROR_MOD,     ONLY : ERROR_STOP
      USE UNIX_CMDS_MOD, ONLY : SEPARATOR

#     include "define.h"      ! C-preprocessor flags

      ! Arguments
      CHARACTER(LEN=*), INTENT(INOUT) :: DIR
      
      ! Local variables
      LOGICAL                         :: IT_EXISTS
      INTEGER                         :: C
      CHARACTER(LEN=255)              :: MSG
#if   defined( COMPAQ )
      INTEGER*4, EXTERNAL             :: ACCESS
#endif
      
      !=================================================================
      ! CHECK_DIRECTORY begins here!
      !=================================================================

      ! Locate the last non-white-space character of NEWDIR
      C = LEN_TRIM( DIR )

      ! Add the trailing directory separator if it is not present
      IF ( DIR(C:C) /= TRIM( SEPARATOR ) ) THEN 
         DIR(C+1:C+1) = TRIM( SEPARATOR )
      ENDIF
     
      !=================================================================
      ! Test if the directory actually exists
      !=================================================================

#if   defined( COMPAQ )

      ! Test whether directory exists for COMPAQ
      IT_EXISTS = ( ACCESS( TRIM( DIR ), ' ' ) == 0 )

#else

      ! Test whether directory exists w/ F90 INQUIRE function
      INQUIRE( FILE=TRIM( DIR ), EXIST=IT_EXISTS )

#endif     
 
      ! If the directory doesn't exist, write an error
      ! message, flush the std output buffer, and stop
      IF ( .not. IT_EXISTS ) THEN
         MSG = 'Invalid directory: ' // TRIM( DIR ) 
         CALL ERROR_STOP( MSG, 'CHECK_DIRECTORY ("input_mod.f")' )
      ENDIF

      ! Return to calling program
      END SUBROUTINE CHECK_DIRECTORY

!------------------------------------------------------------------------------

      SUBROUTINE CHECK_TIME_STEPS
!
!******************************************************************************
!  Subroutine CHECK_TIME_STEPS computes the smallest dynamic time step for the
!  model, based on which operation are turned on.  This is called from routine
!  READ_INPUT_FILE, after all of the timesteps and logical flags have been
!  read from "input.geos". (bmy, 7/20/04)
!
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE LOGICAL_MOD
      USE TIME_MOD,  ONLY : SET_TIMESTEPS
      
      ! Local variables
      INTEGER              :: I, J, K, L, TS_SMALLEST
      
      !=================================================================
      ! CHECK_TIME_STEPS begins here!
      !=================================================================

      ! NUNIT is time step in minutes for unit conversion
      TS_UNIT = -1

      ! Only do unit conversion if 
      IF ( LTRAN .or. LCONV .or. LTURB ) THEN
         TS_UNIT = MAX( TS_DYN, TS_CONV )
      ENDIF

      ! Compute NSMALLEST as the minimum of NDYN, NCONV, NSRCE, NCHEM
      I = TS_DYN
      J = TS_CONV
      K = TS_EMIS
      L = TS_CHEM

      IF ( .not. LTRAN                  ) I = 999999 
      IF ( .not. LCONV .and..not. LTURB ) J = 999999
      IF ( .not. LDRYD .and..not. LEMIS ) K = 999999
      IF ( .not. LCHEM                  ) L = 999999

      ! Get the smallest of all of the above
      TS_SMALLEST = MIN( I, J, K, L )

      ! If all of the operators above are turned off, 
      ! then set NSMALLEST to NDYN.
      IF ( TS_SMALLEST == 999999 ) THEN 
         TS_SMALLEST = TS_DYN
      ENDIF
            
      ! If NDYN is smaller than NSMALLEST, reset NSMALLEST
      ! to NDYN.  Also reset NTDT and NSTP accordingly.
      ! This is useful for runs where transport is turned off,
      ! but where chemistry is turned on. 
      IF ( TS_DYN < TS_SMALLEST ) THEN
         TS_DYN = TS_SMALLEST
      ENDIF

      ! Initialize timesteps in "time_mod.f"
      CALL SET_TIMESTEPS( CHEMISTRY=TS_CHEM, EMISSION=TS_EMIS, 
     &                    DYNAMICS=TS_DYN,   UNIT_CONV=TS_UNIT,
     &                    CONVECTION=TS_CONV )
      
      ! Return to MAIN program
      END SUBROUTINE CHECK_TIME_STEPS

!------------------------------------------------------------------------------

      SUBROUTINE IS_LAST_DAY_GOOD
!
!******************************************************************************
!  Suborutine IS_LAST_DAY_GOOD tests to see if there is output scheduled on 
!  the last day of the run.  (bmy, 1/11/05)
!
!  NOTES:
!  (1 ) Moved to "input_mod.f" from "main.f" (bmy, 1/11/05)
!******************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD,  ONLY : ERROR_STOP
      USE JULDAY_MOD, ONLY : JULDAY
      USE TIME_MOD,   ONLY : GET_NYMDe, ITS_A_LEAPYEAR, YMD_EXTRACT

#     include "CMN_SIZE"   ! Size parameters
#     include "CMN_DIAG"   ! NJDAY

      ! Local variables
      INTEGER             :: NYMDe, Y, M, D, LASTDAY
      REAL*8              :: JD, JD0

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
     &                    'IS_LAST_DAY_GOOD ("input_mod.f")' )
      ENDIF
     
      ! Return to calling program
      END SUBROUTINE IS_LAST_DAY_GOOD

!------------------------------------------------------------------------------

      SUBROUTINE INIT_INPUT
!
!******************************************************************************
!  Subroutine INIT_INPUT initializes all variables from "directory_mod.f" and
!  "logical_mod.f" for safety's sake. (bmy, 7/20/04, 11/5/04)
!
!  NOTES:
!  (1 ) Now also initialize LNEI99 from "logical_mod.f" (bmy, 11/5/04)
!  (2 ) Now also initialize LAVHRRLAI from "logical_mod.f" (bmy, 12/20/04)
!******************************************************************************
!
      ! References to F90 modules
      USE DIRECTORY_MOD
      USE LOGICAL_MOD

      !=================================================================
      ! INIT_INPUT begins here!
      !=================================================================

      ! Initialize directories
      DATA_DIR   = ''
      GEOS_1_DIR = ''
      GEOS_S_DIR = ''
      GEOS_3_DIR = ''
      GEOS_4_DIR = ''
      TEMP_DIR   = ''
      RUN_DIR    = ''
      OH_DIR     = ''
      O3PL_DIR   = ''
      TPBC_DIR   = ''

      ! Initialize logicals
      LATEQ      = .FALSE.
      LAVHRRLAI  = .FALSE.
      LCARB      = .FALSE.
      LDEAD      = .FALSE.
      LDUST      = .FALSE.
      LSULF      = .FALSE.
      LSOA       = .FALSE.
      LSSALT     = .FALSE.
      LCHEM      = .FALSE.
      LEMBED     = .FALSE.
      LCONV      = .FALSE.
      LDBUG      = .FALSE.
      LDIAG      = .FALSE.
      LPRT       = .FALSE.
      LSTDRUN    = .FALSE.
      LDRYD      = .FALSE. 
      LAIRNOX    = .FALSE.
      LANTHRO    = .FALSE.
      LBIONOX    = .FALSE.   
      LBIOMASS   = .FALSE.
      LBIOFUEL   = .FALSE.
      LBIOGENIC  = .FALSE.
      LBBSEA     = .FALSE.
      LEMIS      = .FALSE.
      LFFNOX     = .FALSE.               
      LFOSSIL    = .FALSE.  
      LLIGHTNOX  = .FALSE.
      LMONOT     = .FALSE.
      LNEI99     = .FALSE.
      LSHIPSO2   = .FALSE.
      LSOILNOX   = .FALSE.
      LTOMSAI    = .FALSE.
      LWOODCO    = .FALSE. 
      LFILL      = .FALSE.
      LMFCT      = .FALSE.
      LTRAN      = .FALSE.
      LTPFV      = .FALSE.
      LUPBD      = .FALSE.
      LWINDO     = .FALSE.
      LUNZIP     = .FALSE.
      LWAIT      = .FALSE.
      LTURB      = .FALSE.
      LSVGLB     = .FALSE.
      LSPLIT     = .FALSE.
      LWETD      = .FALSE.

      ! Initialize counters
      CT1        = 0
      CT2        = 0
      CT3        = 0

      ! Return to calling program
      END SUBROUTINE INIT_INPUT

!------------------------------------------------------------------------------

      ! End of module
      END MODULE INPUT_MOD
