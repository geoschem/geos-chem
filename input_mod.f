      MODULE INPUT_MOD
!
!******************************************************************************
!
!  LEMIS    : If = T, switch on EMISSIONS       subroutines
!  LDRYD    : If = T, switch on DRY DEPOSITION  subroutines
!  LCHEM    : If = T, switch on CHEMISTRY       subroutines
!  LTRAN    : If = T, switch on TRANSPORT       subroutines
!  LKZZ     : If = T, switch on VERTICAL DIFFUSION subroutines
!  LTURB    : If = T, switch on DRY CONVECTION  subroutines
!  LCONV    : If = T, switch on WET CONVECTION  subroutines
!  LWETD    : If = T, switch on WET DEPOSITION  subroutines
!  LDBUG    : If = T, write info to 'debug1' and 'debug2' files
!  LDIAG    : If = T, switch on NDxx DIAGNOSTIC subroutines
!  LRERD    : If = T, re-read data for the same day
!  LMFCT    : Argument for TPCORE (transport)
!  LFILL    : Argument for TPCORE (transport)
!  LWAIT    : If = T, wait for met fields to be unzipped in the foreground
!  LTOMSAI  : If = T, will use TOMS AI data for B.B. interannual variability
!  LBBSEA   : If = T, will use seasonal biomass burning emissions
!  LSTDRUN  : If = T, then save quantities for STDRUN model evaluation
!  LFORCE   : If = T, will force initialization to happen from punch file
!  LSPLIT   : If = T, will split CO emissions into separate regions
!  LMONOT   : If = T, will use monoterpenes emissions for acetone source
!  BACKGROUND  : String variable for background operator  (' &'    in Unix) 
!  REDIRECT    : String variable for redirection operator (' >'    in Unix)
!  REMOVE_CMD  : String variable for remove command       ('rm'    in Unix)
!  SEPARATOR   : String variable for dir path separator   ('/'     in Unix)
!  SPACE       : String variable for blank spaces         (' '     in Unix)    
!  STAR        : String variable for wild card operator   ('*'     in Unix)
!  UNZIP_CMD   : String variable for unzip command        ('gzcat' in Unix) 
!
!  A3_SUFFIX   : Suffix for DAO A-3  (Average 3h      ) met fields file
!  A6_SUFFIX   : Suffix for DAO A-6  (Average 6h      ) met fields file
!  I6_SUFFIX   : Suffix for DAO I-6  (Instantaneous 6h) met fields file
!  PH_SUFFIX   : Suffix for DAO PHIS (geopotential hts) met fields file
!  KZZ_SUFFIX  : Suffix for DAO KZZ  (Average 3h      ) met fields file
!  GRID_SUFFIX : Suffix for grid resolution
!  ZIP_SUFFIX  : Suffix for denoting compressed files
!
!  DATA_DIR    : Main DAO met field directory
!  GEOS_1_DIR  : Subdirectory of DATA_DIR where GEOS-1     data are stored
!  GEOS_S_DIR  : Subdirectory of DATA_DIR where GEOS-STRAT data are stored
!  TEMP_DIR    : Directory for temporary storage of met field files
!******************************************************************************
!
      IMPLICIT NONE

      !=================================================================
      ! MODULE PRIVATE DECLARATIONS -- keep certain internal variables 
      ! and routines from being seen outside "input_mod.f"
      !=================================================================
      
      ! PRIVATE module routines
      PRIVATE :: READ_ONE_LINE,        READ_SIMULATION_MENU
      PRIVATE :: READ_GRID_SIZE_MENU,  READ_MET_FIELDS_MENU
      PRIVATE :: READ_TRACER_MENU,     READ_EMISSIONS_MENU
      PRIVATE :: READ_CHEMISTRY_MENU,  READ_TRANSPORT_MENU
      PRIVATE :: READ_CONVECTION_MENU, READ_DEPOSITION_MENU
      PRIVATE :: READ_OUTPUT_MENU,     READ_DIAGNOSTIC_MENU
      PRIVATE :: READ_ND49_MENU,       READ_ND50_MENU
      PRIVATE :: READ_ND51_MENU,       READ_PRODLOSS_MENU
      PRIVATE :: READ_UNIX_CMDS_MENU,  SPLIT_ONE_LINE

      !=================================================================
      ! MODULE VARIABLES 
      !=================================================================
      LOGICAL                       :: DEBUG    = .TRUE.
      LOGICAL                       :: VERBOSE  = .TRUE.
      INTEGER, PARAMETER            :: FIRSTCOL = 26
      INTEGER, PARAMETER            :: MAXDIM   = 255
      CHARACTER(LEN=255)            :: FILENAME = 'input.geos.new'

      ! Top Title
      CHARACTER(LEN=255)            :: TOPTITLE

      ! Simulation Menu
      LOGICAL                       :: LSVGLB
      INTEGER                       :: NYMDb, NHMSb 
      INTEGER                       :: NYMDe, NHMSe
      CHARACTER(LEN=255)            :: DATA_DIR
      CHARACTER(LEN=255)            :: IN_RESTART_FILE
      CHARACTER(LEN=255)            :: OUT_RESTART_FILE
      CHARACTER(LEN=255)            :: GEOS_1_DIR,  GEOS_S_DIR
      CHARACTER(LEN=255)            :: GEOS_3_DIR,  GEOS_4_DIR
      CHARACTER(LEN=255)            :: RUN_DIR,     TEMP_DIR

      ! Met fields menu
      LOGICAL                       :: LWAIT,     LRERD
      CHARACTER(LEN=255)            :: A3_FILE,   A6_FILE,   I6_FILE   
      CHARACTER(LEN=255)            :: PHIS_FILE, GWET_FILE, KZZ_FILE

      ! Tracer menu     
      INTEGER                       :: NSRCX, NTRACE, N_EMISSION
      INTEGER                       :: N_BIOMASS, N_BIOFUEL
      INTEGER,          ALLOCATABLE :: ID_TRACER(:)
      INTEGER,          ALLOCATABLE :: ID_EMISSION(:)
      INTEGER,          ALLOCATABLE :: ID_BIOMASS(:)  
      INTEGER,          ALLOCATABLE :: ID_BIOFUEL(:)
      CHARACTER(LEN=9), ALLOCATABLE :: TRACER_NAME(:)
      REAL*8,           ALLOCATABLE :: TRACER_MOLWT(:)

      ! Emission menu
      INTEGER                       :: NSRCE,     FSCALYR
      LOGICAL                       :: LEMIS,     LFOSSIL,  LBIOFUEL 
      LOGICAL                       :: LBIOGENIC, LMONOT,   LBIOMASS
      LOGICAL                       :: LBBSEA,    LTOMSAI,  LAIRNOX
      LOGICAL                       :: LLIGHTNOX, LSOILNOX

      ! Chemistry menu
      LOGICAL                       :: LCHEM, LEMBED
      INTEGER                       :: NCHEM, IEBD1, JEBD1, IEBD2, JEBD2

      ! Transport menu
      LOGICAL                       :: LTRAN, LUPBD
      INTEGER                       :: NDYN,  IORD, JORD, KORD, J1, KS
      REAL*8                        :: UMAX

      ! Convection menu
      LOGICAL                       :: LCONV, LTURB, LKZZ
      INTEGER                       :: NCONV

      ! Deposition menu
      LOGICAL                       :: LDRYD, LWETD
      INTEGER                       :: NDEP

      ! Output menu
      INTEGER                       :: NJDAY(366)

      ! Diagnostic menu
      CHARACTER(LEN=255)            :: BPCH_FILE

      ! ND49 menu (temporary)
      INTEGER                       :: ND49
      CHARACTER(LEN=255)            :: ND49_OUTPUT_FILE
      REAL*8,           ALLOCATABLE :: ND49_TRACERS(:)
      INTEGER                       :: ND49_IMIN, ND49_IMAX
      INTEGER                       :: ND49_JMIN, ND49_JMAX
      INTEGER                       :: ND49_LMIN, ND49_LMAX
      INTEGER                       :: ND49_DMIN, ND49_DMAX
      
      ! ND50 menu (temporary)
      INTEGER                       :: ND50
      CHARACTER(LEN=255)            :: ND50_OUTPUT_FILE
      REAL*8,           ALLOCATABLE :: ND50_TRACERS(:)
      INTEGER                       :: ND50_IMIN, ND50_IMAX
      INTEGER                       :: ND50_JMIN, ND50_JMAX
      INTEGER                       :: ND50_LMIN, ND50_LMAX
      INTEGER                       :: ND50_DMIN, ND50_DMAX

      ! ND41 menu (temporary)
      INTEGER                       :: ND51
      CHARACTER(LEN=255)            :: ND51_OUTPUT_FILE
      REAL*8,           ALLOCATABLE :: ND51_TRACERS(:)
      INTEGER                       :: ND51_IMIN, ND51_IMAX
      INTEGER                       :: ND51_JMIN, ND51_JMAX
      INTEGER                       :: ND51_LMIN, ND51_LMAX
      INTEGER                       :: ND51_DMIN, ND51_DMAX

      ! Unix cmds menu
      CHARACTER(LEN=1)              :: SPACE       
      CHARACTER(LEN=255)            :: BACKGROUND,  REDIRECT 
      CHARACTER(LEN=255)            :: REMOVE_CMD,  SEPARATOR   
      CHARACTER(LEN=255)            :: UNZIP_CMD,   WILD_CARD 
      CHARACTER(LEN=255)            :: GZIP_SUFFIX

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement 
      !=================================================================
      CONTAINS

!------------------------------------------------------------------------------

      SUBROUTINE READ_INPUT_FILE
!
!******************************************************************************
!  Subroutine READ_INPUT_FILE is the driver program for reading the GEOS_CHEM
!  input file from disk. (bmy, 7/11/02)
!
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE FILE_MOD, ONLY : IU_GEOS, IOERROR
      
      ! Local variables
      LOGICAL            :: EOF
      INTEGER            :: IOS
      CHARACTER(LEN=255) :: LINE

      !=================================================================
      ! READ_INPUT_FILE begins here!
      !=================================================================  

      ! Echo output
      WRITE( 6, '(a)' ) REPEAT( '=', 79 )
      WRITE( 6, 100   ) TRIM( FILENAME )
 100  FORMAT( 'READ_INPUT_FILE: Reading ', a, / )

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
         
         !=============================================================
         ! Call individual subroutines to read sections of the file
         ! NOTE: The order in which the menus are listed in the file
         ! does not matter!  You can change them around!
         !=============================================================
         IF      ( INDEX( LINE, 'SIMULATION MENU' ) > 0 ) THEN
            CALL READ_SIMULATION_MENU             
                                                  
         ELSE IF ( INDEX( LINE, 'GRID SIZE MENU'  ) > 0 ) THEN
            CALL READ_GRID_SIZE_MENU              
                                                  
         ELSE IF ( INDEX( LINE, 'MET FIELDS MENU' ) > 0 ) THEN
            CALL READ_MET_FIELDS_MENU             
                                                  
         ELSE IF ( INDEX( LINE, 'TRACER MENU'     ) > 0 ) THEN
            CALL READ_TRACER_MENU                 
                                                  
         ELSE IF ( INDEX( LINE, 'EMISSIONS MENU'  ) > 0 ) THEN
            CALL READ_EMISSIONS_MENU              
                                                  
         ELSE IF ( INDEX( LINE, 'CHEMISTRY MENU'  ) > 0 ) THEN
            CALL READ_CHEMISTRY_MENU              
                                                  
         ELSE IF ( INDEX( LINE, 'TRANSPORT MENU'  ) > 0 ) THEN
            CALL READ_TRANSPORT_MENU              
                                                  
         ELSE IF ( INDEX( LINE, 'CONVECTION MENU' ) > 0 ) THEN
            CALL READ_CONVECTION_MENU             
                                                  
         ELSE IF ( INDEX( LINE, 'DEPOSITION MENU' ) > 0 ) THEN
            CALL READ_DEPOSITION_MENU             
                                                  
         ELSE IF ( INDEX( LINE, 'OUTPUT MENU'     ) > 0 ) THEN
            CALL READ_OUTPUT_MENU                 
                                                  
         ELSE IF ( INDEX( LINE, 'DIAGNOSTIC MENU' ) > 0 ) THEN
            CALL READ_DIAGNOSTIC_MENU             
                                                  
         ELSE IF ( INDEX( LINE, 'ND49 MENU'       ) > 0 ) THEN
            CALL READ_ND49_MENU                   
                                                  
         ELSE IF ( INDEX( LINE, 'ND50 MENU'       ) > 0 ) THEN
            CALL READ_ND50_MENU                   
                                                  
         ELSE IF ( INDEX( LINE, 'ND51 MENU'       ) > 0 ) THEN
            CALL READ_ND51_MENU                   
                                                  
         ELSE IF ( INDEX( LINE, 'PRODLOSS MENU'   ) > 0 ) THEN
            CALL READ_PRODLOSS_MENU               
                                                  
         ELSE IF ( INDEX( LINE, 'UNIX CMDS MENU'  ) > 0 ) THEN 
            CALL READ_UNIX_CMDS_MENU              

         ELSE IF ( INDEX( LINE, 'OFFLINE OH MENU' ) > 0 ) THEN 
            CALL READ_OFFLINE_OH_MENU
                                                  
         ELSE IF ( INDEX( LINE, 'END OF FILE'     ) > 0 ) THEN 
            EXIT

         ELSE
            ! Nothing
            
         ENDIF  
      ENDDO

      ! Close input file
      CLOSE( IU_GEOS )

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
!  be done. (bmy, 7/11/02)
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
      ! READ_LINE begins here!
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
!  READ_ONE_LINE), and separates it into substrings.  
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

      !=================================================================
      ! Subroutine READ_SIMULATION_MENU reads the SIMULATION MENU
      ! section of the GEOS-CHEM input file (bmy, 7/11/02)
      !
      ! NOTES:
      !=================================================================

      ! Local variables
      INTEGER            :: N
      CHARACTER(LEN=255) :: SUBSTRS(MAXDIM)

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
      READ( SUBSTRS(1:N), '(a)' ) IN_RESTART_FILE

      ! Make new restart file?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_simulation_menu:5' )
      READ( SUBSTRS(1:N), '(a)' ) LSVGLB

      ! Output restart file(s)
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_simulation_menu:6' )
      READ( SUBSTRS(1:N), '(a)' ) OUT_RESTART_FILE

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

      ! Separator line
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_simulation_menu:13' )

      !=================================================================
      ! Print to screen
      !=================================================================
      IF ( DEBUG ) THEN
         PRINT*, NYMDb
         PRINT*, NHMSb
         PRINT*, NYMDe
         PRINT*, NHMSe
         PRINT*, TRIM( RUN_DIR )
         PRINT*, TRIM( IN_RESTART_FILE )
         PRINT*, LSVGLB
         PRINT*, TRIM( OUT_RESTART_FILE )
         PRINT*, TRIM( GEOS_1_DIR )
         PRINT*, TRIM( GEOS_S_DIR )
         PRINT*, TRIM( GEOS_3_DIR )
         PRINT*, TRIM( GEOS_4_DIR )
         PRINT*, TRIM( TEMP_DIR )
      ENDIF

      ! Return to calling program
      END SUBROUTINE READ_SIMULATION_MENU

!------------------------------------------------------------------------------

      SUBROUTINE READ_GRID_SIZE_MENU

      !=================================================================
      ! Subroutine READ_GRID_SIZE_MENU reads the GRID SIZE MENU
      ! section of the GEOS-CHEM input file (bmy, 7/11/02)
      !
      ! NOTES:
      !=================================================================

      ! References to F90 modules
      USE GRID_MOD, ONLY : IM, JM, LM, I0, J0, COMPUTE_GRID

      ! Local variables
      INTEGER            :: N, AS
      CHARACTER(LEN=255) :: SUBSTRS(MAXDIM)

      !=================================================================
      ! READ_GRID_SIZE_MENU begins here!
      !=================================================================

      ! IM, JM, LM
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 3, 'read_grid_size_menu:1' )
      READ( SUBSTRS(1:N), * ) IM, JM, LM

      ! I0, J0
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 2, 'read_grid_size_menu:2' )
      READ( SUBSTRS(1:N), * ) I0, J0

      ! Separator line
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_grid_size_menu:3' )

      ! Compute lat/lon/surface area variables
      CALL COMPUTE_GRID

      !=================================================================
      ! Print to screen
      !=================================================================
      IF ( DEBUG ) THEN
         PRINT*, IM, JM, LM
         PRINT*, I0, J0
      ENDIF

      ! Return to calling program
      END SUBROUTINE READ_GRID_SIZE_MENU

!------------------------------------------------------------------------------

      SUBROUTINE READ_MET_FIELDS_MENU

      !=================================================================
      ! Subroutine READ_MET_FIELDS_MENU reads the MET FIELDS MENU
      ! section of the GEOS-CHEM input file (bmy, 7/11/02)
      !
      ! NOTES:
      !=================================================================

      ! Local variables
      INTEGER            :: N
      CHARACTER(LEN=255) :: SUBSTRS(MAXDIM)

      !=================================================================
      ! READ_MET_FIELDS_MENU begins here!
      !=================================================================

      ! Wait for met fields?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_met_fields_menu:1' )
      READ( SUBSTRS(1:N), * ) LWAIT

      ! Reread one day of met data?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_met_fields_menu:2' )
      READ( SUBSTRS(1:N), * ) LRERD

      ! A-3 file
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_met_fields_menu:3' )
      READ( SUBSTRS(1:N), '(a)' ) A3_FILE

      ! A-6 file
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_met_fields_menu:4' )
      READ( SUBSTRS(1:N), '(a)' ) A6_FILE

      ! I-6 file
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_met_fields_menu:5' )
      READ( SUBSTRS(1:N), '(a)' ) I6_FILE

      ! PHIS file
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_met_fields_menu:6' )
      READ( SUBSTRS(1:N), '(a)' ) PHIS_FILE

      ! GWET file (GEOS-3 only)
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_met_fields_menu:7' )
      READ( SUBSTRS(1:N), '(a)' ) GWET_FILE

      ! KZZ file
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_met_fields_menu:8' )
      READ( SUBSTRS(1:N), '(a)' ) KZZ_FILE

      ! Divider
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_met_fields_menu:9' )

      ! Met fields menu
      IF ( DEBUG ) THEN
         PRINT*, LWAIT
         PRINT*, LRERD
         PRINT*, TRIM( A3_FILE   )
         PRINT*, TRIM( A6_FILE   )
         PRINT*, TRIM( I6_FILE   )
         PRINT*, TRIM( PHIS_FILE )
         PRINT*, TRIM( GWET_FILE )
         PRINT*, TRIM( KZZ_FILE  )
      ENDIF

      ! Return to calling program
      END SUBROUTINE READ_MET_FIELDS_MENU

!------------------------------------------------------------------------------

      SUBROUTINE READ_TRACER_MENU

      !=================================================================
      ! Subroutine READ_TRACER_MENU reads the TRACER MENU
      ! section of the GEOS-CHEM input file (bmy, 7/11/02)
      !
      ! NOTES:
      !=================================================================

      ! Local variables
      LOGICAL            :: EOF
      INTEGER            :: N, AS, I, I1, I2, I3, I4
      CHARACTER(LEN=255) :: SUBSTRS(MAXDIM), LINE, C1
      REAL*8             :: R1

      !=================================================================
      ! READ_TRACER_MENU begins here!
      !=================================================================

      ! NSRCX
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_tracer_menu:1' )
      READ( SUBSTRS(1:N), * ) NSRCX

      ! NTRACE
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_tracer_menu:2' )
      READ( SUBSTRS(1:N), * ) NTRACE
      
      ! Separator line
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 6, 'read_tracer_menu:3' )

      !=================================================================
      ! Allocate arrays 
      !=================================================================
      ALLOCATE( ID_TRACER( NTRACE ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'ID_TRACER' )
      ID_TRACER = 0

      ALLOCATE( ID_EMISSION( NTRACE ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'ID_EMISSION' )
      ID_EMISSION = 0

      ALLOCATE( ID_BIOMASS( NTRACE ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'ID_BIOMASS' )
      ID_BIOMASS = 0

      ALLOCATE( ID_BIOFUEL( NTRACE ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'ID_BIOFUEL' )
      ID_BIOFUEL = 0

      ALLOCATE( TRACER_NAME( NTRACE ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'TRACER_NAME' )
      TRACER_NAME = ''

      ALLOCATE( TRACER_MOLWT( NTRACE ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'TRACER_MOLWT' )
      TRACER_MOLWT = 0d0

      ! Also initialize counter variables
      N_EMISSION = 0
      N_BIOFUEL  = 0
      N_BIOMASS  = 0
     
      !=================================================================
      ! Read tracer ID, emission ID, biomass ID, biofuel ID, & mol wt
      !=================================================================
      DO N = 1, NTRACE

         ! Read line, trap premature EOF
         LINE = READ_ONE_LINE( EOF, 'read_tracer_menu:1' ) 

         ! Break LINE up into individual variables
         READ( LINE(FIRSTCOL:), '(4i4,1x,a10,f6.1)' ) I1,I2,I3,I4,C1,R1

         ! Copy into arrays
         ID_TRACER(N)    = I1
         ID_EMISSION(N)  = I2
         ID_BIOMASS(N)   = I3
         ID_BIOFUEL(N)   = I4
         TRACER_NAME(N)  = TRIM( C1 )
         TRACER_MOLWT(N) = R1

         ! Keep track of the # of biofuel, biomass, and emission tracers
         IF ( ID_BIOFUEL(N)  > 0 ) N_BIOFUEL  = N_BIOFUEL  + 1
         IF ( ID_BIOMASS(N)  > 0 ) N_BIOMASS  = N_BIOMASS  + 1
         IF ( ID_EMISSION(N) > 0 ) N_EMISSION = N_EMISSION + 1

      ENDDO

      ! Separator line
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_tracer_menu:5' )

      !=================================================================
      ! Print values
      !=================================================================
      IF ( DEBUG ) THEN
         DO N = 1, NTRACE
            WRITE( 6, '(4i4,1x,a10,f6.1)' )
     &           ID_TRACER(N),  ID_EMISSION(N), ID_BIOMASS(N), 
     &           ID_BIOFUEL(N), TRACER_NAME(N), TRACER_MOLWT(N)   
         ENDDO

         PRINT*, 'N_EMISSION: ', N_EMISSION
         PRINT*, 'N_BIOFUEL : ', N_BIOFUEL
         PRINT*, 'N_BIOMASS : ', N_BIOMASS
      ENDIF
     
      ! Return to calling program
      END SUBROUTINE READ_TRACER_MENU

!------------------------------------------------------------------------------

      SUBROUTINE READ_EMISSIONS_MENU

      !=================================================================
      ! Subroutine READ_EMISSIONS_MENU reads the EMISSIONS MENU
      ! section of the GEOS-CHEM input file (bmy, 7/11/02)
      !
      ! NOTES:
      !=================================================================

      ! Local variables
      INTEGER            :: N
      CHARACTER(LEN=255) :: SUBSTRS(MAXDIM)

      !=================================================================
      ! READ_EMISSIONS_MENU begins here!
      !=================================================================

      ! Turn on emissions?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_emissions_menu:1' )
      READ( SUBSTRS(1:N), * ) LEMIS

      ! Emissions timestep
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_emissions_menu:2' )
      READ( SUBSTRS(1:N), * ) NSRCE

      ! Include anthropogenic emissions?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_emissions_menu:3' )
      READ( SUBSTRS(1:N), * ) LFOSSIL

      ! Scale 1985 to year ???
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

      ! Use Aircraft NOx
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_emissions_menu:12' )
      READ( SUBSTRS(1:N), * ) LAIRNOX

      ! Use lightning NOx
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_emissions_menu:13' )
      READ( SUBSTRS(1:N), * ) LLIGHTNOX

      ! Use soil NOx
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_emissions_menu:14' )
      READ( SUBSTRS(1:N), * ) LSOILNOX

      ! Separator line
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_emissions_menu:15' )
      
      !### Debug
      IF ( DEBUG ) THEN
         PRINT*, LEMIS
         PRINT*, NSRCE
         PRINT*, LFOSSIL
         PRINT*, FSCALYR
         PRINT*, LBIOFUEL
         PRINT*, LBIOMASS
         PRINT*, LMONOT
         PRINT*, LBIOMASS
         PRINT*, LBBSEA
         PRINT*, LTOMSAI
         PRINT*, LAIRNOX
         PRINT*, LAIRNOX
         PRINT*, LLIGHTNOX
         PRINT*, LLIGHTNOX
      ENDIF

      ! Return to calling program
      END SUBROUTINE READ_EMISSIONS_MENU

!------------------------------------------------------------------------------

      SUBROUTINE READ_CHEMISTRY_MENU

      !=================================================================
      ! Subroutine READ_CHEMISTRY_MENU reads the CHEMISTRY MENU
      ! section of the GEOS-CHEM input file (bmy, 7/11/02)
      !
      ! NOTES:
      !=================================================================

      ! Local variables
      INTEGER            :: N
      CHARACTER(LEN=255) :: SUBSTRS(MAXDIM), I1, I2, I3, I4

      !=================================================================
      ! READ_CHEMISTRY_MENU begins here!
      !=================================================================

      ! Turn on chemistry?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_chemistry_menu:1' )
      READ( SUBSTRS(1:N), * ) LCHEM

      ! Chemistry timestep
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_chemistry_menu:2' )
      READ( SUBSTRS(1:N), * ) NCHEM

      ! Use embedded chemistry?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_chemistry_menu:3' )
      READ( SUBSTRS(1:N), * ) LEMBED

      ! Embedded chemistry limits?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 4, 'read_chemistry_menu:4' )
      READ( SUBSTRS(1:N), * ) IEBD1, JEBD1, IEBD2, JEBD2

      ! Separator line
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_chemistry_menu:5' )

      !### Debug
      IF ( DEBUG ) THEN
         PRINT*, LCHEM
         PRINT*, NCHEM
         PRINT*, LEMBED
         PRINT*, IEBD1, JEBD1, IEBD2, JEBD2
      ENDIF

      ! Return to calling program
      END SUBROUTINE READ_CHEMISTRY_MENU

!------------------------------------------------------------------------------

      SUBROUTINE READ_TRANSPORT_MENU

      !=================================================================
      ! Subroutine READ_TRANSPORT_MENU reads the CHEMISTRY MENU
      ! section of the GEOS-CHEM input file (bmy, 7/11/02)
      !
      ! NOTES:
      !=================================================================

      ! Local variables
      INTEGER            :: N
      CHARACTER(LEN=255) :: SUBSTRS(MAXDIM)

      !=================================================================
      ! READ_TRANSPORT_MENU begins here!
      !=================================================================

      ! Turn on transport?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_transport_menu:1' )
      READ( SUBSTRS(1:N), * ) LTRAN

      ! Transport timestep
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_transport_menu:2' )
      READ( SUBSTRS(1:N), * ) NDYN

      ! Include strat O3/NOy
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_transport_menu:3' )
      READ( SUBSTRS(1:N), * ) LUPBD

      ! IORD, JORD, KORD
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 3, 'read_transport_menu:4' )
      READ( SUBSTRS(1:N), * ) IORD, JORD, KORD

      ! UMAX, J1, KS
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 3, 'read_transport_menu:4' )
      READ( SUBSTRS(1:N), * ) UMAX, J1, KS

      ! Separator line
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_transport_menu:5' )

      !### Debug
      IF ( DEBUG ) THEN
         PRINT*, LTRAN
         PRINT*, NDYN
         PRINT*, LUPBD
         PRINT*, IORD, JORD, KORD
         PRINT*, UMAX, J1, KS
      ENDIF
      
      ! Return to calling program
      END SUBROUTINE READ_TRANSPORT_MENU

!------------------------------------------------------------------------------

      SUBROUTINE READ_CONVECTION_MENU

      !=================================================================
      ! Subroutine READ_CONVECTION_MENU reads the CHEMISTRY MENU
      ! section of the GEOS-CHEM input file (bmy, 7/11/02)
      !
      ! NOTES:
      !=================================================================

      ! Local variables
      INTEGER            :: N
      CHARACTER(LEN=255) :: SUBSTRS(MAXDIM)

      !=================================================================
      ! READ_CONVECTION_MENU begins here!
      !=================================================================

      ! Turn on convection?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_convection_menu:1' )
      READ( SUBSTRS(1:N), * ) LCONV

      ! Turn on BL mixing
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_convection_menu:2' )
      READ( SUBSTRS(1:N), * ) LTURB

      ! Turn on KZZ mixing?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_convection_menu:3' )
      READ( SUBSTRS(1:N), * ) LKZZ

      ! Convection timestep
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_convection_menu:4' )
      READ( SUBSTRS(1:N), * ) NCONV

      ! Separator line
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_convection_menu:5' )

      !### Debug
      IF ( DEBUG ) THEN
         PRINT*, LCONV
         PRINT*, LTURB
         PRINT*, LKZZ
         PRINT*, NCONV
      ENDIF

      ! Return to calling program 
      END SUBROUTINE READ_CONVECTION_MENU

!------------------------------------------------------------------------------

      SUBROUTINE READ_DEPOSITION_MENU

      !=================================================================
      ! Subroutine READ_DEPOSITION_MENU reads the CHEMISTRY MENU
      ! section of the GEOS-CHEM input file (bmy, 7/11/02)
      !
      ! NOTES:
      !=================================================================

      ! Local variables
      INTEGER            :: N
      CHARACTER(LEN=255) :: SUBSTRS(MAXDIM)

      !=================================================================
      ! READ_DEPOSITION_MENU begins here!
      !=================================================================

      ! Turn on drydep?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_deposition_menu:1' )
      READ( SUBSTRS(1:N), * ) LDRYD

      ! Turn on wetdep?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_deposition_menu:2' )
      READ( SUBSTRS(1:N), * ) LWETD

      ! Deposition timestep 
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_deposition_menu:3' )
      READ( SUBSTRS(1:N), * ) NDEP

      ! Separator line
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_deposition_menu:4' )

      !### Debug
      IF ( DEBUG ) THEN
         PRINT*, LDRYD
         PRINT*, LWETD
         PRINT*, NDEP
      ENDIF

      ! Return to calling program
      END SUBROUTINE READ_DEPOSITION_MENU

!------------------------------------------------------------------------------

      SUBROUTINE READ_OUTPUT_MENU

      !=================================================================
      ! Subroutine READ_OUTPUT_MENU reads the OUTPUT MENU
      ! section of the GEOS-CHEM input file (bmy, 7/11/02)
      !
      ! NOTES:
      !=================================================================

      ! References to F90 modules
      USE DIAG_MOD, ONLY : NJDAY
      USE FILE_MOD, ONLY : IU_GEOS, IOERROR

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

      !### Debug
      IF ( DEBUG ) THEN 
         WRITE( 6, 100 ) NJDAY
      ENDIF

      ! Return to calling program
      END SUBROUTINE READ_OUTPUT_MENU

!------------------------------------------------------------------------------

      SUBROUTINE SET_TINDEX( N_DIAG, SUBSTRS, N, NMAX  )
      
      !=================================================================
      ! Subroutine SET_TINDEX sets the TINDEX and TMAX arrays, which
      ! determine how many tracers to print to the punch file.
      ! (bmy, 7/30/02)
      !=================================================================

      ! References to F90 modules
      USE DIAG_MOD, ONLY : TMAX, TINDEX

      ! Arguments
      INTEGER,            INTENT(IN) :: N_DIAG, N, NMAX
      CHARACTER(LEN=255), INTENT(IN) :: SUBSTRS(N)

      ! Local variables
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
      IF ( TRIM( SUBSTRS(1) ) == 'all' ) THEN 

         ! TMAX is the max # of tracers to print out
         TMAX(N_DIAG) = NMAX 

         ! Set TINDEX = 1..NTRACE
         DO M = 1, TMAX(N_DIAG)
            TINDEX(N_DIAG,M) = M
         ENDDO

      ELSE 

         ! Otherwise, set TMAX, TINDEX to the # of tracers
         ! listed in "input.ctm" -- need some error checks too
         TMAX(N_DIAG) = N-1
         READ( SUBSTRS(1:N), * ) ( TINDEX(N_DIAG,M), M=1,N )

      ENDIF

      ! Return to calling program
      END SUBROUTINE SET_TINDEX

!------------------------------------------------------------------------------

      SUBROUTINE READ_DIAGNOSTIC_MENU

      !=================================================================
      ! Subroutine READ_DIAGNOSTIC_MENU reads the DIAGNOSTIC MENU
      ! section of the GEOS-CHEM input file (bmy, 7/30/02)
      !
      ! NOTES:
      !=================================================================

      ! References to F90 modules
      USE DIAG_MOD

#     include "CMN_SIZE" ! Size parameters

      ! Local variables
      INTEGER            :: N, M
      CHARACTER(LEN=255) :: SUBSTRS(MAXDIM)

      !=================================================================
      ! READ_DIAGNOSTIC_MENU begins here!
      !=================================================================

      ! Binary punch file name
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1,  'read_diagnostic_menu:1' )
      READ( SUBSTRS(1:N), '(a)' ) BPCH_FILE

      ! Separator line
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:2' )

      ! ND01
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:3' )
      READ( SUBSTRS(1), * ) ND01
      CALL SET_TINDEX( 01, SUBSTRS(2:N), N-1, NTRACE )

      ! ND02
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:4' )
      READ( SUBSTRS(1), * ) ND02
      CALL SET_TINDEX( 02, SUBSTRS(2:N), N-1, NTRACE )

      ! ND03
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:5' )
      READ( SUBSTRS(1), * ) ND03
      CALL SET_TINDEX( 03, SUBSTRS(2:N), N-1, NTRACE )

      ! ND04
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:6' )
      READ( SUBSTRS(1), * ) ND04
      CALL SET_TINDEX( 04, SUBSTRS(2:N), N-1, NTRACE )

      ! ND05
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:7' )
      READ( SUBSTRS(1), * ) ND05
      CALL SET_TINDEX( 05, SUBSTRS(2:N), N-1, NTRACE )

      ! ND06
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:8' )
      READ( SUBSTRS(1), * ) ND06
      CALL SET_TINDEX( 06, SUBSTRS(2:N), N-1, NTRACE )

      ! ND07
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:9' )
      READ( SUBSTRS(1), * ) ND07
      CALL SET_TINDEX( 07, SUBSTRS(2:N), N-1, NTRACE )

      ! ND08
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:10' )
      READ( SUBSTRS(1), * ) ND08
      CALL SET_TINDEX( 08, SUBSTRS(2:N), N-1, NTRACE )

      ! ND09
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:11' )
      READ( SUBSTRS(1), * ) ND09
      CALL SET_TINDEX( 09, SUBSTRS(2:N), N-1, NTRACE )

      ! ND10
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:12' )
      READ( SUBSTRS(1), * ) ND10
      CALL SET_TINDEX( 10, SUBSTRS(2:N), N, NTRACE )

      ! ND11
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:13' )
      READ( SUBSTRS(1), * ) ND11
      CALL SET_TINDEX( 11, SUBSTRS(2:N), N-1, 7 )

      ! ND12
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:14' )
      READ( SUBSTRS(1), * ) ND12
      CALL SET_TINDEX( 12, SUBSTRS(2:N), N-1, 1 )

      ! ND13
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:15' )
      READ( SUBSTRS(1), * ) ND13
      CALL SET_TINDEX( 13, SUBSTRS(2:N), N-1, 7 )

      ! ND14
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:16' )
      READ( SUBSTRS(1), * ) ND14
      CALL SET_TINDEX( 14, SUBSTRS(2:N), N-1, NTRACE )

      ! ND15
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:17' )
      READ( SUBSTRS(1), * ) ND15
      CALL SET_TINDEX( 15, SUBSTRS(2:N), N-1, NTRACE )

      ! ND16
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:18' )
      READ( SUBSTRS(1), * ) ND16
      CALL SET_TINDEX( 16, SUBSTRS(2:N), N-1, NTRACE )

      ! ND17 -- revisit this
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:19' )
      READ( SUBSTRS(1), * ) ND17
      CALL SET_TINDEX( 17, SUBSTRS(2:N), N-1, 4 )

      ! ND18 -- revisit this
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:20' )
      READ( SUBSTRS(1), * ) ND18
      CALL SET_TINDEX( 18, SUBSTRS(2:N), N-1, 4 )

      ! ND19
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:21' )
      READ( SUBSTRS(1), * ) ND19
      CALL SET_TINDEX( 19, SUBSTRS(2:N), N, NTRACE )

      ! ND20
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:22' )
      READ( SUBSTRS(1), * ) ND20
      CALL SET_TINDEX( 20, SUBSTRS(2:N), N-1, 2 )      

      ! ND21
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:23' )
      READ( SUBSTRS(1), * ) ND21
      CALL SET_TINDEX( 21, SUBSTRS(2:N), N-1, 10 )

      ! ND22
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:24' )
      READ( SUBSTRS(1), * ) ND22
      CALL SET_TINDEX( 22, SUBSTRS(2:N), N-1, 6 )

      CALL SPLIT_ONE_LINE( SUBSTRS, N, 2,  'read_diagnostic_menu:25' ) 
      READ( SUBSTRS(1:N), * ) HR1_JV, HR2_JV

      ! ND23
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:26' )
      READ( SUBSTRS(1), * ) ND23
      CALL SET_TINDEX( 23, SUBSTRS(2:N), N-1, 1 )

      ! ND24
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:27' )
      READ( SUBSTRS(1), * ) ND24
      CALL SET_TINDEX( 24, SUBSTRS(2:N), N-1, NTRACE )

      ! ND25
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:28' )
      READ( SUBSTRS(1), * ) ND25
      CALL SET_TINDEX( 25, SUBSTRS(2:N), N-1, NTRACE )

      ! ND26
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:29' )
      READ( SUBSTRS(1), * ) ND26
      CALL SET_TINDEX( 26, SUBSTRS(2:N), N-1, NTRACE )

      ! ND27
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:30' )
      READ( SUBSTRS(1), * ) ND28
      CALL SET_TINDEX( 27, SUBSTRS(2:N), N-1, 3 )

      ! ND28
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:31' )
      READ( SUBSTRS(1), * ) ND28
      CALL SET_TINDEX( 28, SUBSTRS(2:N), N-1, N_BIOMASS )

      ! ND29
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:32' )
      READ( SUBSTRS(1), * ) ND29
      CALL SET_TINDEX( 29, SUBSTRS(2:N), N-1, 5 )

      ! ND30
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:33' )
      READ( SUBSTRS(1), * ) ND30
      CALL SET_TINDEX( 30, SUBSTRS(2:N), N-1, 1 )

      ! ND31
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:34' )
      READ( SUBSTRS(1), * ) ND31
      CALL SET_TINDEX( 31, SUBSTRS(2:N), N-1, 1 )

      ! ND32
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:34' )
      READ( SUBSTRS(1), * ) ND32
      CALL SET_TINDEX( 32, SUBSTRS(2:N), N-1, 1 )

      ! ND33
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:35' )
      READ( SUBSTRS(1), * ) ND33
      CALL SET_TINDEX( 33, SUBSTRS(2:N), N-1, NTRACE )

      ! ND34
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:36' )
      READ( SUBSTRS(1), * ) ND34
      CALL SET_TINDEX( 34, SUBSTRS(2:N), N-1, N_BIOFUEL )

      ! ND35
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:37' )
      READ( SUBSTRS(1), * ) ND35
      CALL SET_TINDEX( 35, SUBSTRS(2:N), N-1, NTRACE )

      ! ND36
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:38' )
      READ( SUBSTRS(1), * ) ND36
      CALL SET_TINDEX( 36, SUBSTRS(2:N), N-1, N_EMISSION )

      ! ND37
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:39' )
      READ( SUBSTRS(1), * ) ND37
      CALL SET_TINDEX( 37, SUBSTRS(2:N), N-1, 4 )

      ! ND38
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:40' )
      READ( SUBSTRS(1), * ) ND38
      CALL SET_TINDEX( 38, SUBSTRS(2:N), N-1, 4 )

      ! ND39
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:41' )
      READ( SUBSTRS(1), * ) ND39
      CALL SET_TINDEX( 39, SUBSTRS(2:N), N-1, 4 )

      ! ND40
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:42' )
      READ( SUBSTRS(1), * ) ND40
      CALL SET_TINDEX( 40, SUBSTRS(2:N), N-1, 1 )

      ! ND41
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:43' )
      READ( SUBSTRS(1), * ) ND41
      CALL SET_TINDEX( 41, SUBSTRS(2:N), N-1, 1 )

      ! ND42
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:44' )
      READ( SUBSTRS(1), * ) ND42
      CALL SET_TINDEX( 42, SUBSTRS(2:N), N-1, 0 )

      ! ND43
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:45' )
      READ( SUBSTRS(1), * ) ND43
      CALL SET_TINDEX( 43, SUBSTRS(2:N), N-1, 4 )

      CALL SPLIT_ONE_LINE( SUBSTRS, N, 2,  'read_diagnostic_menu:46' )
      READ( SUBSTRS(1:N), * ) HR1_OH, HR2_OH

      CALL SPLIT_ONE_LINE( SUBSTRS, N, 2,  'read_diagnostic_menu:47' )
      READ( SUBSTRS(1:N), * ) HR1_NO, HR2_NO

      PRINT*, 'HR1_OH, HR2_OH: ', HR1_OH, HR2_OH
      PRINT*, 'HR1_NO, HR2_NO: ', HR1_NO, HR2_NO

      ! ND44
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:48' )
      READ( SUBSTRS(1), * ) ND44
      CALL SET_TINDEX( 44, SUBSTRS(2:N), N-1, 46 )

      ! ND45
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:49' )
      READ( SUBSTRS(1), * ) ND45
      CALL SET_TINDEX( 45, SUBSTRS(2:N), N-1, NTRACE )

      ! ND46
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:50' )
      READ( SUBSTRS(1), * ) ND46
      CALL SET_TINDEX( 46, SUBSTRS(2:N), N-1, 4 )

      ! ND47
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:51' )
      READ( SUBSTRS(1), * ) ND47
      CALL SET_TINDEX( 47, SUBSTRS(2:N), N-1, NTRACE )

      ! ND52
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:52' )
      READ( SUBSTRS(1), * ) ND52
      CALL SET_TINDEX( 52, SUBSTRS(2:N), N-1, NTRACE )

      ! ND53
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:53' )
      READ( SUBSTRS(1), * ) ND53
      CALL SET_TINDEX( 53, SUBSTRS(2:N), N-1, NTRACE )

      ! ND54
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:54' )
      READ( SUBSTRS(1), * ) ND54
      CALL SET_TINDEX( 54, SUBSTRS(2:N), N-1, NTRACE )

      ! ND55
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:55' )
      READ( SUBSTRS(1), * ) ND55
      CALL SET_TINDEX( 55, SUBSTRS(2:N), N-1, 3 )

      ! ND56
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:56' )
      READ( SUBSTRS(1), * ) ND56
      CALL SET_TINDEX( 56, SUBSTRS(2:N), N-1, NTRACE )

      ! ND57
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:57' )
      READ( SUBSTRS(1), * ) ND57
      CALL SET_TINDEX( 57, SUBSTRS(2:N), N-1, NTRACE )

      ! ND58
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:58' )
      READ( SUBSTRS(1), * ) ND58
      CALL SET_TINDEX( 58, SUBSTRS(2:N), N-1, NTRACE )

      ! ND59
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:59' )
      READ( SUBSTRS(1), * ) ND59
      CALL SET_TINDEX( 59, SUBSTRS(2:N), N-1, NTRACE )

      ! ND60
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:60' )
      READ( SUBSTRS(1), * ) ND60
      CALL SET_TINDEX( 60, SUBSTRS(2:N), N-1, NTRACE )

      ! ND61
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:61' )
      READ( SUBSTRS(1), * ) ND61
      CALL SET_TINDEX( 61, SUBSTRS(2:N), N-1, NTRACE )

      ! ND62
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:62' )
      READ( SUBSTRS(1), * ) ND62
      CALL SET_TINDEX( 62, SUBSTRS(2:N), N-1, NTRACE )

      ! ND63
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:63' )
      READ( SUBSTRS(1), * ) ND63
      CALL SET_TINDEX( 63, SUBSTRS(2:N), N-1, NTRACE )

      ! ND64
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:64' )
      READ( SUBSTRS(1), * ) ND64
      CALL SET_TINDEX( 64, SUBSTRS(2:N), N-1, NTRACE )

      ! ND66
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:66' )
      READ( SUBSTRS(1), * ) ND66
      CALL SET_TINDEX( 66, SUBSTRS(2:N), N-1, 6 )

      ! ND67
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:67' )
      READ( SUBSTRS(1), * ) ND67
      CALL SET_TINDEX( 67, SUBSTRS(2:N), N-1, 18 )

      ! ND68
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:68' )
      READ( SUBSTRS(1), * ) ND68
      CALL SET_TINDEX( 68, SUBSTRS(2:N), N-1, 4 )

      ! ND69
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:69' )
      READ( SUBSTRS(1), * ) ND69
      CALL SET_TINDEX( 69, SUBSTRS(2:N), N-1, 1 )

      ! ND70
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:70' )
      READ( SUBSTRS(1), * ) ND70
      CALL SET_TINDEX( 70, SUBSTRS(2:N), N-1, 0 )

      ! Separator line
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1,  'read_diagnostic_menu:71' )

      ! print TINDEX
      IF ( DEBUG ) THEN
         DO N = 1, 70
            PRINT*, N, ' : ', TINDEX(N,1:24)
         ENDDO
      ENDIF

      ! Return to calling program
      END SUBROUTINE READ_DIAGNOSTIC_MENU

!------------------------------------------------------------------------------

      SUBROUTINE READ_ND49_MENU

      !=================================================================
      ! Subroutine READ_ND49_MENU reads the UNIX COMMANDS MENU
      ! section of the GEOS-CHEM input file (bmy, 7/11/02)
      !
      ! NOTES:
      !=================================================================

      ! Local variables
      LOGICAL            :: LND49
      INTEGER            :: N, I, AS
      CHARACTER(LEN=255) :: SUBSTRS(MAXDIM)

      !=================================================================
      ! READ_ND49_MENU begins here!
      !=================================================================
      
      ! Initialize
      ND49 = 0

      ! Turn on ND49 diagnostic
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1,  'read_nd49_menu:1' )
      READ( SUBSTRS(1:N), * ) LND49
      IF ( LND49 ) ND49 = 1

      ! Instantaneous 3-D timeseries file
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1,  'read_nd49_menu:2' )
      READ( SUBSTRS(1:N), '(a)' ) ND49_OUTPUT_FILE

      ! Tracers to include
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_nd49_menu:3' )
 
      ALLOCATE( ND49_TRACERS( N ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'ND49_TRACERS' )
      ND49_TRACERS = 0

      READ( SUBSTRS(1:N), * ) ( ND49_TRACERS(I), I=1,N )

      ! IMIN, IMAX
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 2,  'read_nd49_menu:4' )
      READ( SUBSTRS(1:N), * ) ND49_IMIN, ND49_IMAX

      ! JMIN, JMAX
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 2,  'read_nd49_menu:5' )
      READ( SUBSTRS(1:N), * ) ND49_JMIN, ND49_JMAX

      ! LMIN, LMAX
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 2,  'read_nd49_menu:6' )
      READ( SUBSTRS(1:N), * ) ND49_LMIN, ND49_LMAX

      ! DMIN, DMAX
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 2,  'read_nd49_menu:7' )
      READ( SUBSTRS(1:N), * ) ND49_DMIN, ND49_DMAX

      ! Separator line
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1,  'read_nd49_menu:8' )

      !### Debug
      IF ( DEBUG ) THEN
         PRINT*, ND49
         PRINT*, TRIM( ND49_OUTPUT_FILE )
         PRINT*, ND49_TRACERS
         PRINT*, ND49_IMIN, ND49_IMAX
         PRINT*, ND49_JMIN, ND49_JMAX
         PRINT*, ND49_LMIN, ND49_LMAX
         PRINT*, ND49_DMIN, ND49_DMAX
      ENDIF

      ! Return to calling program
      END SUBROUTINE READ_ND49_MENU

!------------------------------------------------------------------------------

      SUBROUTINE READ_ND50_MENU

      ! Local variables
      LOGICAL            :: LND50
      INTEGER            :: N, I, AS
      CHARACTER(LEN=255) :: SUBSTRS(MAXDIM)

      !=================================================================
      ! READ_ND50_MENU begins here!
      !=================================================================
      
      ! Initialize
      ND50 = 0

      ! Turn on ND49 diagnostic
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1,  'read_nd50_menu:1' )
      READ( SUBSTRS(1:N), * ) LND50
      IF ( LND50 ) ND50 = 1

      ! Instantaneous 3-D timeseries file
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1,  'read_nd50_menu:2' )
      READ( SUBSTRS(1:N), '(a)' ) ND50_OUTPUT_FILE

      ! Tracers to include
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_nd50_menu:3' )
 
      ALLOCATE( ND50_TRACERS( N ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'ND50_TRACERS' )
      ND50_TRACERS = 0

      READ( SUBSTRS(1:N), * ) ( ND50_TRACERS(I), I=1,N )

      ! IMIN, IMAX
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 2,  'read_nd50_menu:4' )
      READ( SUBSTRS(1:N), * ) ND50_IMIN, ND50_IMAX

      ! JMIN, JMAX
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 2,  'read_nd50_menu:5' )
      READ( SUBSTRS(1:N), * ) ND50_JMIN, ND50_JMAX

      ! LMIN, LMAX
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 2,  'read_nd50_menu:6' )
      READ( SUBSTRS(1:N), * ) ND50_LMIN, ND50_LMAX

      ! DMIN, DMAX
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 2,  'read_nd50_menu:7' )
      READ( SUBSTRS(1:N), * ) ND50_DMIN, ND50_DMAX

      ! Separator line
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1,  'read_nd50_menu:8' )

      !### Debug
      IF ( DEBUG ) THEN
         PRINT*, ND50
         PRINT*, TRIM( ND50_OUTPUT_FILE )
         PRINT*, ND50_TRACERS
         PRINT*, ND50_IMIN, ND50_IMAX
         PRINT*, ND50_JMIN, ND50_JMAX
         PRINT*, ND50_LMIN, ND50_LMAX
         PRINT*, ND50_DMIN, ND50_DMAX
      ENDIF

      ! Return to calling program
      END SUBROUTINE READ_ND50_MENU

!------------------------------------------------------------------------------

      SUBROUTINE READ_ND51_MENU

      ! Local variables
      LOGICAL            :: LND51
      INTEGER            :: N, I, AS
      CHARACTER(LEN=255) :: SUBSTRS(MAXDIM)

      !=================================================================
      ! READ_ND51_MENU begins here!
      !=================================================================
      
      ! Initialize
      ND51 = 0

      ! Turn on ND49 diagnostic
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1,  'read_nd51_menu:1' )
      READ( SUBSTRS(1:N), * ) LND51
      IF ( LND51 ) ND51 = 1

      ! Instantaneous 3-D timeseries file
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1,  'read_nd51_menu:2' )
      READ( SUBSTRS(1:N), '(a)' ) ND51_OUTPUT_FILE

      ! Tracers to include
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_nd51_menu:3' )
 
      ALLOCATE( ND51_TRACERS( N ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'ND51_TRACERS' )
      ND51_TRACERS = 0

      READ( SUBSTRS(1:N), * ) ( ND51_TRACERS(I), I=1,N )

      ! IMIN, IMAX
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 2,  'read_nd51_menu:4' )
      READ( SUBSTRS(1:N), * ) ND51_IMIN, ND51_IMAX

      ! JMIN, JMAX
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 2,  'read_nd51_menu:5' )
      READ( SUBSTRS(1:N), * ) ND51_JMIN, ND51_JMAX

      ! LMIN, LMAX
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 2,  'read_nd51_menu:6' )
      READ( SUBSTRS(1:N), * ) ND51_LMIN, ND51_LMAX

      ! DMIN, DMAX
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 2,  'read_nd51_menu:7' )
      READ( SUBSTRS(1:N), * ) ND51_DMIN, ND51_DMAX

      ! Separator line
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1,  'read_nd51_menu:8' )

      !### Debug
      IF ( DEBUG ) THEN
         PRINT*, ND51
         PRINT*, TRIM( ND51_OUTPUT_FILE )
         PRINT*, ND51_TRACERS
         PRINT*, ND51_IMIN, ND51_IMAX
         PRINT*, ND51_JMIN, ND51_JMAX
         PRINT*, ND51_LMIN, ND51_LMAX
         PRINT*, ND51_DMIN, ND51_DMAX
      ENDIF

      ! Return to calling program
      END SUBROUTINE READ_ND51_MENU

!------------------------------------------------------------------------------

      SUBROUTINE READ_PRODLOSS_MENU

      ! Local variables
      LOGICAL            :: EOF
      INTEGER            :: N
      CHARACTER(LEN=255) :: LINE

      ! Skip for now
      DO N = 1, 38
         LINE = READ_ONE_LINE( EOF ) 
         IF ( EOF ) STOP
      ENDDO

      ! Return to calling program
      END SUBROUTINE READ_PRODLOSS_MENU

!------------------------------------------------------------------------------

      SUBROUTINE READ_UNIX_CMDS_MENU

      !=================================================================
      ! Subroutine READ_UNIX_CMDS_MENU reads the UNIX COMMANDS MENU
      ! section of the GEOS-CHEM input file (bmy, 7/11/02)
      !
      ! NOTES:
      !=================================================================

      ! References to F90 modules
      USE CHARPAK_MOD, ONLY : STRSQUEEZE

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

      ! Space
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_unix_cmds_menu:5' )
      READ( SUBSTRS(1:N), '(a)' ) SPACE

      ! Wild Card
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_unix_cmds_menu:6' )
      READ( SUBSTRS(1:N), '(a)' ) WILD_CARD

      ! Unzip command
      UNZIP_CMD = READ_ONE_LINE( EOF,     'read_unix_cmds_menu:7' ) 
      UNZIP_CMD = UNZIP_CMD(FIRSTCOL:)
      CALL STRSQUEEZE( UNZIP_CMD )

      ! Gzip suffix
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_unix_cmds_menu:8' )
      READ( SUBSTRS(1:N), '(a)' ) GZIP_SUFFIX

      ! Separator line
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_unix_cmds_menu:9' )

      !### Debug
      IF ( DEBUG ) THEN
         PRINT*, TRIM( BACKGROUND  )     
         PRINT*, TRIM( REDIRECT    )    
         PRINT*, TRIM( REMOVE_CMD  )    
         PRINT*, TRIM( SEPARATOR   )    
         PRINT*, SPACE            
         PRINT*, TRIM( WILD_CARD   )        
         PRINT*, TRIM( UNZIP_CMD   )
         PRINT*, TRIM( GZIP_SUFFIX )
      ENDIF

      ! Return to calling program
      END SUBROUTINE READ_UNIX_CMDS_MENU

!------------------------------------------------------------------------------
      
      SUBROUTINE CLEANUP_INPUT

      !=================================================================
      ! Subroutine CLEANUP_INPUT deallocates all allocatable arrays
      !=================================================================
      IF ( ALLOCATED( ID_TRACER    ) ) DEALLOCATE( ID_TRACER    )
      IF ( ALLOCATED( ID_EMISSION  ) ) DEALLOCATE( ID_EMISSION  ) 
      IF ( ALLOCATED( ID_BIOMASS   ) ) DEALLOCATE( ID_BIOMASS   ) 
      IF ( ALLOCATED( ID_BIOFUEL   ) ) DEALLOCATE( ID_BIOFUEL   ) 
      IF ( ALLOCATED( TRACER_NAME  ) ) DEALLOCATE( TRACER_NAME  )
      IF ( ALLOCATED( TRACER_MOLWT ) ) DEALLOCATE( TRACER_MOLWT )
      IF ( ALLOCATED( ND49_TRACERS ) ) DEALLOCATE( ND49_TRACERS )

      ! Return to calling program
      END SUBROUTINE CLEANUP_INPUT

!------------------------------------------------------------------------------

      END MODULE INPUT_MOD
