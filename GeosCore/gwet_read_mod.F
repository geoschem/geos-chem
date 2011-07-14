! $Id: gwet_read_mod.f,v 1.1 2009/09/16 14:06:25 bmy Exp $
      MODULE GWET_READ_MOD
!
!******************************************************************************
!  Module GWET_READ_MOD contains routines that unzip, open, and 
!  read the GEOS-CHEM GWET (avg 3-hour) met fields from disk. 
!  (tdf, bmy, 3/30/04, 8/4/06)
!
!  NOTE: GWET fields are included in GEOS-4 met data, so we only need this
!        module to read GWET data from the GEOS-3 data sets. (bmy, 3/30/04)
! 
!  Module Routines:
!  =========================================================================
!  (1 ) UNZIP_GWET_FIELDS : Unzips & copies met field files to a temp dir
!  (2 ) DO_OPEN_GWET      : Returns TRUE if it's time to read GWET fields
!  (3 ) OPEN_GWET_FIELDS  : Opens met field files residing in the temp dir
!  (4 ) GET_GWET_FIELDS   : Wrapper for routine READ_GWET
!  (5 ) CHECK_TIME        : Tests if GWET met field times equal current time
!  (6 ) READ_GWET         : Reads GWET fields from disk
!  (7 ) GWET_CHECK        : Checks if we have found all of the GWET fields
! 
!  GEOS-CHEM modules referenced by gwet_read_mod.f
!  ============================================================================
!  (1 ) bpch2_mod.f     : Module containing routines for binary punch file I/O
!  (2 ) dao_mod.f       : Module containing arrays for DAO met fields
!  (3 ) diag_mod.f      : Module containing GEOS-CHEM diagnostic arrays
!  (4 ) directory_mod.f : Module containing GEOS-CHEM data & met field dirs
!  (5 ) error_mod.f     : Module containing NaN and other error check routines
!  (6 ) logical_mod.f   : Module containing GEOS-CHEM logical switches 
!  (7 ) file_mod.f      : Module containing file unit #'s and error checks
!  (8 ) time_mod.f      : Module containing routines for computing time & date
!  (9 ) transfer_mod.f  : Module containing routines to cast & resize arrays
!  (10) unix_cmds_mod.f : Module containing Unix commands for unzipping files
!
!  NOTES:
!  (1 ) Adapted from "a3_read_mod.f" (tdf, rjp, 6/30/04)
!  (2 ) Now references "directory_mod.f", "logical_mod.f", and 
!        "unix_cmds_mod.f" (bmy, 7/20/04)
!  (3 ) Now references FILE_EXISTS from "file_mod.f" (bmy, 3/23/05)
!  (4 ) Updated comments (bmy, 10/24/05)
!  (5 ) Remove support for GEOS-1 and GEOS-STRAT met fields (bmy, 8/4/06)
!******************************************************************************
!
      IMPLICIT NONE

      !=================================================================
      ! MODULE PRIVATE DECLARATIONS -- keep certain internal variables 
      ! and routines from being seen outside "gwet_read_mod.f"
      !=================================================================

      ! PRIVATE module routines
      PRIVATE :: GWET_CHECK, CHECK_TIME, DO_OPEN_GWET, READ_GWET

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement 
      !=================================================================
      CONTAINS

!------------------------------------------------------------------------------

      SUBROUTINE UNZIP_GWET_FIELDS( OPTION, NYMD )
!
!******************************************************************************
!  Subroutine UNZIP_GWET_FIELDS invokes a FORTRAN system call to uncompress
!  GEOS-CHEM GWET met field files and store the uncompressed data in a 
!  temporary directory, where GEOS-CHEM can read them.  The original data 
!  files are not disturbed.  (tdf, bmy, 3/30/04, 8/4/06)
!
!  Arguments as input:
!  ============================================================================
!  (1 ) OPTION (CHAR*(*)) : Option
!  (2 ) NYMD   (INTEGER ) : YYYYMMDD of GWET file to be unzipped (optional)
!
!  NOTES:
!  (1 ) Adapted from "a3_read_mod.f" (tdf, bmy, 3/30/04)
!  (2 ) Now reference "directory_mod.f" and "unix_cmds_mod.f".  Also remove 
!        reference to CMN_SETUP. (bmy, 7/20/04)
!  (3 ) Now reference "directory_mod.f" and "unix_cmds_mod.f". Now prevent 
!        EXPAND_DATE from overwriting directory paths with Y/M/D tokens in 
!        them (bmy, 7/20/04)
!  (4 ) Remove support for GEOS-1 and GEOS-STRAT met fields (bmy, 8/4/06)
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD,     ONLY : GET_RES_EXT
      USE DIRECTORY_MOD, ONLY : DATA_DIR,   GEOS_3_DIR, TEMP_DIR
      USE ERROR_MOD,     ONLY : ERROR_STOP
      USE TIME_MOD,      ONLY : EXPAND_DATE
      USE UNIX_CMDS_MOD, ONLY : BACKGROUND, REDIRECT,   REMOVE_CMD 
      USE UNIX_CMDS_MOD, ONLY : UNZIP_CMD,  WILD_CARD,  ZIP_SUFFIX

#     include "CMN_SIZE"            ! Size parameters

      ! Arguments
      CHARACTER(LEN=*),  INTENT(IN) :: OPTION
      INTEGER, OPTIONAL, INTENT(IN) :: NYMD

      ! Local variables
      CHARACTER(LEN=255)            :: GEOS_DIR,     GWET_STR
      CHARACTER(LEN=255)            :: GWET_FILE_GZ, GWET_FILE
      CHARACTER(LEN=255)            :: UNZIP_BG,     UNZIP_FG
      CHARACTER(LEN=255)            :: REMOVE_ALL,   REMOVE_DATE

      !=================================================================
      ! UNZIP_GWET_FIELDS begins here!
      !=================================================================
      IF ( PRESENT( NYMD ) ) THEN
 
         ! String w/ date & resolution
         GEOS_DIR = TRIM( GEOS_3_DIR )
         GWET_STR = 'YYYYMMDD.gwet.' // GET_RES_EXT() 

         ! Replace date tokens
         CALL EXPAND_DATE( GEOS_DIR, NYMD, 000000 )
         CALL EXPAND_DATE( GWET_STR, NYMD, 000000 )

         ! Location of zipped GWET file in data dir
         GWET_FILE_GZ = TRIM( DATA_DIR   ) // TRIM( GEOS_DIR   ) // 
     &                  TRIM( GWET_STR   ) // TRIM( ZIP_SUFFIX )

         ! Location of unzipped GWET file in temp dir
         GWET_FILE    = TRIM( TEMP_DIR   ) // TRIM( GWET_STR   )
         
         ! Remove GWET files for this date from temp dir
         REMOVE_DATE  = TRIM( REMOVE_CMD ) // ' '                // 
     &                  TRIM( TEMP_DIR   ) // TRIM( GWET_STR   )

         !==============================================================
         ! Define the foreground and background UNZIP commands
         !==============================================================

         ! Foreground unzip
         UNZIP_FG = TRIM( UNZIP_CMD ) // ' ' // TRIM( GWET_FILE_GZ ) // 
     &              TRIM( REDIRECT  ) // ' ' // TRIM( GWET_FILE    )  

         ! Background unzip
         UNZIP_BG  = TRIM( UNZIP_FG ) // TRIM( BACKGROUND )
      ENDIF

      !=================================================================
      ! Define command to remove all GWET files from the TEMP dir
      !=================================================================
      REMOVE_ALL = TRIM( REMOVE_CMD ) // ' '    // TRIM( TEMP_DIR   ) // 
     &             TRIM( WILD_CARD  ) // '.gwet.' // TRIM( WILD_CARD  ) 

      !=================================================================
      ! Perform an F90 system call to do the desired operation
      !=================================================================
      SELECT CASE ( TRIM( OPTION ) )
         
         ! Unzip GWET fields in the Unix foreground
         CASE ( 'unzip foreground' )
            WRITE( 6, 100 ) TRIM( GWET_FILE_GZ )
            CALL SYSTEM( TRIM( UNZIP_FG ) )

         ! Unzip A-3 fields in the Unix background
         CASE ( 'unzip background' )
            WRITE( 6, 100 ) TRIM( GWET_FILE_GZ )
            CALL SYSTEM( TRIM( UNZIP_BG ) )

         ! Remove A-3 field for this date in temp dir
         CASE ( 'remove date' )
            WRITE( 6, 110 ) TRIM( GWET_FILE )
            CALL SYSTEM( TRIM( REMOVE_DATE ) )
            
         ! Remove all A-3 fields in temp dir
         CASE ( 'remove all' )
            WRITE( 6, 120 ) TRIM( REMOVE_ALL )
            CALL SYSTEM( TRIM( REMOVE_ALL ) )

         ! Error -- bad option!
         CASE DEFAULT
            CALL ERROR_STOP( 'Invalid value for OPTION!', 
     &                       'UNZIP_GWET_FIELDS (gwet_read_mod.f)' )
            
      END SELECT

      ! FORMAT strings
 100  FORMAT( '     - Unzipping: ', a )
 110  FORMAT( '     - Removing: ', a )
 120  FORMAT( '     - About to execute command: ', a )

      ! Return to calling program
      END SUBROUTINE UNZIP_GWET_FIELDS

!------------------------------------------------------------------------------

      FUNCTION DO_OPEN_GWET( NYMD, NHMS ) RESULT( DO_OPEN )
!
!******************************************************************************
!  Function DO_OPEN_GWET returns TRUE if is time to open the GWET met field 
!  file or FALSE otherwise.  This prevents us from opening a file which has 
!  already been opened. (tdf, bmy, 3/30/04)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) NYMD (INTEGER) : YYYYMMDD 
!  (2 ) NHMS (INTEGER) :  and HHMMSS to be tested for GWET file open
!
!  NOTES:
!******************************************************************************
!
      ! Arguments
      INTEGER, INTENT(IN) :: NYMD, NHMS 

      ! Local variables
      LOGICAL             :: DO_OPEN
      LOGICAL, SAVE       :: FIRST    = .TRUE.
      INTEGER, SAVE       :: LASTNYMD = -1
      INTEGER, SAVE       :: LASTNHMS = -1
      
      !=================================================================
      ! DO_OPEN_GWET begins here!
      !=================================================================

      ! Initialize
      DO_OPEN = .FALSE.

      ! Return if we have already opened the file
      IF ( NYMD == LASTNYMD .and. NHMS == LASTNHMS ) THEN
         DO_OPEN = .FALSE. 
         GOTO 999
      ENDIF

      ! Open GWET file if it's 00:00 GMT, or on the first call
      IF ( NHMS == 000000 .or. FIRST ) THEN
         DO_OPEN = .TRUE. 
         GOTO 999
      ENDIF

      !=================================================================
      ! Reset quantities for next call
      !=================================================================
 999  CONTINUE
      LASTNYMD = NYMD
      LASTNHMS = NHMS
      FIRST    = .FALSE.

      ! Return to calling program
      END FUNCTION DO_OPEN_GWET

!------------------------------------------------------------------------------

      SUBROUTINE OPEN_GWET_FIELDS( NYMD, NHMS )
!
!******************************************************************************
!  Subroutine OPEN_gwet_FIELDS opens the A-3 met fields file for date NYMD 
!  and time NHMS. (tdf, bmy, 3/30/04, 8/4/06)
!  
!  Arguments as Input:
!  ===========================================================================
!  (1 ) NYMD (INTEGER) : YYYYMMDD
!  (2 ) NHMS (INTEGER) :  and HHMMSS timestamps for A-3 file
!
!  NOTES:
!  (1 ) Adapted from "a3_read_mod.f" (tdf, bmy, 3/30/04)
!  (2 ) Now references "directory_mod.f" instead of CMN_SETUP.  Also now
!        references LUNZIP from "logical_mod.f".  Also now prevents EXPAND_DATE
!        from overwriting Y/M/D tokens in directory paths. (bmy, 7/20/04)
!  (3 ) Now use FILE_EXISTS from "file_mod.f" to determine if file unit 
!        IU_GWET refers to a valid file on disk (bmy, 3/23/05)
!  (4 ) Remove support for GEOS-1 and GEOS-STRAT met fields (bmy, 8/4/06)
!******************************************************************************
!      
      ! References to F90 modules
      USE BPCH2_MOD,     ONLY : GET_RES_EXT
      USE DIRECTORY_MOD, ONLY : DATA_DIR,   GEOS_3_DIR, TEMP_DIR
      USE ERROR_MOD,     ONLY : ERROR_STOP
      USE LOGICAL_MOD,   ONLY : LUNZIP
      USE FILE_MOD,      ONLY : IU_GWET,    IOERROR,    FILE_EXISTS
      USE TIME_MOD,      ONLY : EXPAND_DATE

#     include "CMN_SIZE"      ! Size parameters

      ! Arguments
      INTEGER, INTENT(IN)    :: NYMD, NHMS

      ! Local variables
      LOGICAL                :: DO_OPEN
      LOGICAL                :: IT_EXISTS
      INTEGER                :: IOS
      CHARACTER(LEN=8)       :: IDENT
      CHARACTER(LEN=255)     :: GEOS_DIR
      CHARACTER(LEN=255)     :: GWET_FILE
      CHARACTER(LEN=255)     :: PATH

      !=================================================================
      ! OPEN_GWET_FIELDS begins here!
      !=================================================================

      ! Open A-3 fields at the proper time, or on the first call
      IF ( DO_OPEN_GWET( NYMD, NHMS ) ) THEN

         ! String w/ date and resolution
         GEOS_DIR  = TRIM( GEOS_3_DIR )
         GWET_FILE = 'YYYYMMDD.gwet.' // GET_RES_EXT()

         ! Replace date tokens
         CALL EXPAND_DATE( GEOS_DIR,  NYMD, NHMS )
         CALL EXPAND_DATE( GWET_FILE, NYMD, NHMS )

         ! If unzipping, open GEOS-1 file in TEMP dir
         ! If not unzipping, open GEOS-1 file in DATA dir
         IF ( LUNZIP ) THEN
            PATH = TRIM( TEMP_DIR ) // TRIM( GWET_FILE )
         ELSE
            PATH = TRIM( DATA_DIR ) // 
     &             TRIM( GEOS_DIR ) // TRIM( GWET_FILE )
         ENDIF

         ! Close previously opened GWET file
         CLOSE( IU_GWET )
         
         ! Make sure file unit is valid before we open the file
         IF ( .not. FILE_EXISTS( IU_GWET ) ) THEN
            CALL ERROR_STOP( 'Could not find file!', 
     &                       'OPEN_GWET_FIELDS (gwet_read_mod.f)' )
         ENDIF

         ! Open the file
         OPEN( UNIT   = IU_GWET,       FILE   = TRIM( PATH ),
     &         STATUS = 'OLD',         ACCESS = 'SEQUENTIAL',  
     &         FORM   = 'UNFORMATTED', IOSTAT = IOS )
               
         IF ( IOS /= 0 ) THEN
            CALL IOERROR( IOS, IU_GWET, 'open_gwet_fields:1' )
         ENDIF

         ! Echo info
         WRITE( 6, 100 ) TRIM( PATH )
 100     FORMAT( '     - Opening: ', a )         

      ENDIF

      ! Return to calling program
      END SUBROUTINE OPEN_GWET_FIELDS

!------------------------------------------------------------------------------

      SUBROUTINE GET_GWET_FIELDS( NYMD, NHMS )
!
!******************************************************************************
!  Subroutine GET_GWET_FIELDS is a wrapper for routine READ_GWET.  
!  GET_GWET_FIELDS calls READ_GWET properly for reading GEOS-1, 
!  GEOS-STRAT, GEOS-3, or GEOS-4 met data sets. (tdf, bmy, 3/30/04)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) NYMD (INTEGER) : YYYYMMDD
!  (2 ) NHMS (INTEGER) :  and HHMMSS of A-3 fields to be read from disk
!
!  NOTES:
!  (1 ) Adapted from "a3_read_mod.f" (bmy, 3/30/04)
!******************************************************************************
!
      ! References to F90 modules
      USE DAO_MOD, ONLY : GWETTOP

#     include "CMN_SIZE"  ! Size parameters

      ! Arguments
      INTEGER, INTENT(IN) :: NYMD, NHMS 

      ! Local variables
      INTEGER, SAVE       :: LASTNYMD = -1, LASTNHMS = -1

      !=================================================================
      ! GET_GWET_FIELDS begins here!
      !=================================================================

      ! Skip over previously-read GWET fields
      IF ( NYMD == LASTNYMD .and. NHMS == LASTNHMS ) THEN
         WRITE( 6, 100 ) NYMD, NHMS
 100     FORMAT( '     - GWET met fields for NYMD, NHMS = ', 
     &           i8.8, 1x, i6.6, ' have been read already' ) 
         RETURN
      ENDIF

      !=================================================================
      ! Read the GWET data!
      !=================================================================
      CALL READ_GWET( NYMD=NYMD, NHMS=NHMS, GWET=GWETTOP )

      ! Save NYMD, NHMS for next call
      LASTNYMD = NYMD
      LASTNHMS = NHMS

      ! Return to MAIN program
      END SUBROUTINE GET_GWET_FIELDS

!------------------------------------------------------------------------------

      FUNCTION CHECK_TIME( XYMD, XHMS, NYMD, NHMS ) RESULT( ITS_TIME )
!
!******************************************************************************
!  Function CHECK_TIME checks to see if the timestamp of the GWET field just
!  read from disk matches the current time.  If so, then it's time to return
!  the GWET field to the calling program. (tdf, bmy, 3/30/04, 8/4/06)
!  
!  Arguments as Input:
!  ============================================================================
!  (1 ) XYMD (INTEGER) : YYYYMMDD timestamp for A-3 field in file
!  (2 ) XHMS (INTEGER) : HHMMSS   timestamp for A-3 field in file
!  (3 ) NYMD (INTEGER) : YYYYMMDD  at which A-3 field is to be read
!  (4 ) NHMS (INTEGER) : HHMMSS    at which A-3 field is to be read
!
!  NOTES:
!  (1 ) Adapted from "a3_read_mod.f" (bmy, 3/30/04)
!  (2 ) Remove support for GEOS-1 and GEOS-STRAT met fields (bmy, 8/4/06)
!******************************************************************************
!
#     include "CMN_SIZE"

      ! Arguments 
      INTEGER, INTENT(IN) :: XYMD, XHMS, NYMD, NHMS
      
      ! Function value
      LOGICAL             :: ITS_TIME

      !=================================================================
      ! CHECK_TIME begins here!
      !=================================================================
      IF ( XYMD == NYMD .AND. XHMS == NHMS ) THEN
         ITS_TIME = .TRUE.
      ELSE
         ITS_TIME = .FALSE.
      ENDIF

      ! Return to calling program
      END FUNCTION CHECK_TIME

!-----------------------------------------------------------------------------

      SUBROUTINE READ_GWET( NYMD, NHMS, GWET )
!
!******************************************************************************
!  Subroutine READ_GWET reads GEOS GWET (3-hr avg) fields from disk.
!  (tdf, bmy, 3/30/04, 8/4/06)
! 
!  Arguments as input:
!  ============================================================================
!  (1 ) NYMD : YYYYMMDD
!  (2 ) NHMS :  and HHMMSS of A-3 met fields to be accessed 
!
!  Arguments as Output:
!  ============================================================================
!  (1 ) GWET : (2-D) GMAO topsoil wetness [unitless]
!
!  NOTES:
!  (1 ) Adapted from "a3_read_mod.f" (tdf, bmy, 3/30/04)
!  (2 ) Remove support for GEOS-1 and GEOS-STRAT met fields (bmy, 8/4/06)
!******************************************************************************
!
      ! References to F90 modules
      USE DIAG_MOD,     ONLY : AD67
      USE FILE_MOD,     ONLY : IOERROR,     IU_GWET
      USE TIME_MOD,     ONLY : TIMESTAMP_STRING
      USE TRANSFER_MOD, ONLY : TRANSFER_2D, TRANSFER_TO_1D

#     include "CMN_SIZE"             ! Size parameters
#     include "CMN_DIAG"             ! ND67

      ! Arguments
      INTEGER, INTENT(IN)            :: NYMD, NHMS
      REAL*8,  INTENT(OUT), OPTIONAL :: GWET(IIPAR,JJPAR) 

      ! Local Variables
      INTEGER, PARAMETER             :: N_GWET=1
      INTEGER                        :: I, IJLOOP, IOS, J, NFOUND 
      REAL*4                         :: Q2(IIPAR,JJPAR)
      CHARACTER(LEN=8)               :: NAME
      CHARACTER(LEN=16)              :: STAMP
      INTEGER                        :: XYMD, XHMS

      !=================================================================
      ! READ_gwet begins here!      
      !=================================================================

      ! Zero the number of A-3 fields that we have found
      NFOUND = 0

      !=================================================================
      ! Read the GWET fields from disk
      !=================================================================
      DO

         ! Read the GWET field name
         READ( IU_GWET, IOSTAT=IOS ) NAME

         ! End of file test -- make sure we have found all fields
         IF ( IOS < 0 ) THEN
            CALL GWET_CHECK( NFOUND, N_GWET )
            EXIT
         ENDIF

         ! IOS > 0: True I/O error; stop w/ err msg
         IF ( IOS > 0 ) CALL IOERROR( IOS, IU_GWET, 'read_gwet:1' )

         ! CASE statement for GWET fields
         SELECT CASE ( TRIM( NAME ) )

            !--------------------------------
            ! GWET: Ground wetness (0-1)   
            !--------------------------------
            CASE ( 'GWET', 'GWETTOP' )
               READ( IU_GWET, IOSTAT=IOS ) XYMD, XHMS, Q2
               IF ( IOS /= 0 ) 
     &              CALL IOERROR( IOS, IU_GWET, 'read_gwet:14' )
             
               IF ( CHECK_TIME( XYMD, XHMS, NYMD, NHMS ) ) THEN
                  IF ( PRESENT( GWET ) ) CALL TRANSFER_2D( Q2, GWET )
                  NFOUND = NFOUND + 1
               ENDIF

         END SELECT
               
         !==============================================================
         ! If we have found all the fields for this time, then exit 
         ! the loop.  Otherwise, go on to the next iteration.
         !==============================================================
         IF ( CHECK_TIME( XYMD, XHMS, NYMD, NHMS ) .and. 
     &        NFOUND == N_GWET ) THEN 
            STAMP = TIMESTAMP_STRING( NYMD, NHMS )
            WRITE( 6, 210 ) NFOUND, STAMP
 210        FORMAT( '     - Found all ', i3, ' GWET met fields for ', a)
            EXIT
         ENDIF
      ENDDO

      !=================================================================
      ! ND67 diagnostic: A-3 surface fields:
      !
      ! (22) GWET : Top soil wetness [unitless]
      !=================================================================
      IF ( ND67 > 0 ) THEN
         IF ( PRESENT( GWET ) ) AD67(:,:,22) = AD67(:,:,22) + GWET
      ENDIF
         
      ! Return to calling program
      END SUBROUTINE READ_GWET

!------------------------------------------------------------------------------

      SUBROUTINE GWET_CHECK( NFOUND, N_GWET )
!
!******************************************************************************
!  Subroutine GWET_CHECK prints an error message if not all of the GWET met 
!  fields are found.  The run is also terminated. (tdf, bmy, 3/30/04)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) NFOUND (INTEGER) : # of GWET met fields read from disk
!  (2 ) N_GWET (INTEGER) : # of GWET met fields expected to be read from disk
!
!  NOTES
!  (1 ) Adapted from "a3_read_mod.f" (bmy, 3/30/04)
!******************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD, ONLY : GEOS_CHEM_STOP

      ! Arguments
      INTEGER, INTENT(IN) :: NFOUND, N_GWET

      !=================================================================
      ! GWET_CHECK begins here!
      !=================================================================
      IF ( NFOUND /= N_GWET ) THEN
         WRITE( 6, '(a)' ) REPEAT( '=', 79 )
         WRITE( 6, '(a)' ) 'ERROR -- not enough A-3 fields found!'      

         WRITE( 6, 120   ) N_GWET, NFOUND
 120     FORMAT( 'There are ', i2, ' fields but only ', i2 ,
     &           ' were found!' )

         WRITE( 6, '(a)' ) '### STOP in GWET_CHECK (dao_read_mod.f)'
         WRITE( 6, '(a)' ) REPEAT( '=', 79 )

         ! Deallocate arrays and stop (bmy, 10/15/02)
         CALL GEOS_CHEM_STOP
      ENDIF

      ! Return to calling program
      END SUBROUTINE GWET_CHECK

!------------------------------------------------------------------------------

      ! End of module
      END MODULE GWET_READ_MOD
