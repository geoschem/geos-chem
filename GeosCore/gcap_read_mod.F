! $Id: gcap_read_mod.f,v 1.1 2009/09/16 14:06:26 bmy Exp $
      MODULE GCAP_READ_MOD
!
!******************************************************************************
!  Module PHIS_READ_MOD contains subroutines that unzip, open, and read the 
!  GCAP PHIS and LWI_GISS fields from disk. (bmy, swu, 2/1/06)
! 
!  Module Routines:
!  ============================================================================
!  (1 ) UNZIP_GCAP_FIELDS : Unzips & copies met field files to a temp dir
!  (2 ) OPEN_GCAP_FIELDS  : Opens met field files residing in the temp dir
!  (3 ) GET_GCAP_FIELDS   : Wrapper for routine READ_I6
!  (4 ) CHECK_TIME        : Tests if met field timestamps equal current time
!  (5 ) READ_GCAP         : Reads PHIS fields from disk
!  (6 ) GCAP_CHECK        : Checks if we have found all the PHIS field
! 
!  GEOS-CHEM modules referenced by phis_read_mod.f
!  ============================================================================
!  (1 ) bpch2_mod.f       : Module w/ routines for binary punch file I/O
!  (2 ) dao_mod.f         : Module w/ arrays for DAO met fields
!  (3 ) diag_mod.f        : Module w/ GEOS-CHEM diagnostic arrays
!  (4 ) directory_mod.f   : Module w/ GEOS-CHEM data & met field dirs
!  (5 ) error_mod.f       : Module w/ NaN and other error check routines
!  (6 ) logical_mod.f     : Module w/ GEOS-CHEM logical switches 
!  (7 ) file_mod.f        : Module w/ file unit #'s and error checks
!  (8 ) time_mod.f        : Module w/ routines for computing time & date
!  (9 ) transfer_mod.f    : Module w/ routines to cast & resize arrays
!  (10) unix_cmds_mod.f   : Module w/ Unix commands for unzipping
!
!  NOTES:
!  (1 ) Adapted from the obsolete "phis_read_mod.f" (bmy, 2/1/06)
!******************************************************************************
!
!%%% NOTE: The GEOS-5 GCM won't call this routine to read data, so we
!%%% don't have to add the reference to findFreeLUN (bmy, 8/3/12)
!#if defined( DEVEL )
!      USE inquireMod, ONLY : findFreeLUN
!#endif

      IMPLICIT NONE
#     include "define.h"

      !=================================================================
      ! MODULE PRIVATE DECLARATIONS -- keep certain internal variables 
      ! and routines from being seen outside "gcap_read_mod.f"
      !=================================================================

      ! Make everything PRIVATE ...
      PRIVATE 

      ! ... except these routines
      PUBLIC :: GET_GCAP_FIELDS
      PUBLIC :: OPEN_GCAP_FIELDS
      PUBLIC :: UNZIP_GCAP_FIELDS

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement 
      !=================================================================
      CONTAINS

!------------------------------------------------------------------------------

      SUBROUTINE UNZIP_GCAP_FIELDS( OPTION )
!
!******************************************************************************
!  Subroutine UNZIP_GCAP_FIELDS invokes a FORTRAN system call to uncompress
!  GCAP PHIS met field files and store the uncompressed data in a 
!  temporary directory, where GEOS-CHEM can read them.  The original data 
!  files are not disturbed.  (bmy, bdf, 6/15/98, 5/25/05)
!
!  Arguments as input:
!  ===========================================================================
!  (1 ) OPTION (CHAR*(*)) : Option
!
!  NOTES:
!  (1 ) Adapted from UNZIP_MET_FIELDS of "dao_read_mod.f" (bmy, 6/16/03)
!  (2 ) Directory information YYYY/MM or YYYYMM is now contained w/in 
!        GEOS_1_DIR, GEOS_S_DIR, GEOS_3_DIR, GEOS_4_DIR (bmy, 12/11/03)
!  (3 ) Now reference "directory_mod.f" and "unix_cmds_mod.f". Now prevent 
!        EXPAND_DATE from overwriting directory paths with Y/M/D tokens in 
!        them (bmy, 7/20/04)
!  (4 ) Now modified for GEOS-5 and GCAP met fields (swu, bmy, 5/25/05)
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD,    ONLY : GET_RES_EXT
      USE DIRECTORY_MOD
      USE ERROR_MOD,    ONLY : ERROR_STOP
      USE TIME_MOD,     ONLY : EXPAND_DATE
      USE UNIX_CMDS_MOD

      USE CMN_SIZE_MOD     ! Size parameters

      ! Arguments
      CHARACTER(LEN=*),  INTENT(IN) :: OPTION

      ! Local variables
      INTEGER                       :: NYMD
      CHARACTER(LEN=255)            :: GEOS_DIR,     PHIS_STR
      CHARACTER(LEN=255)            :: PHIS_FILE_GZ, PHIS_FILE
      CHARACTER(LEN=255)            :: UNZIP_BG,     UNZIP_FG
      CHARACTER(LEN=255)            :: REMOVE_ALL,   REMOVE_DATE

      !=================================================================
      ! UNZIP_GCAP_FIELD begins here!
      !=================================================================

      ! Date for PHIS field
      NYMD     = 20000101

      ! Strings for directory & filename
      GEOS_DIR = TRIM( GCAP_DIR )
      PHIS_STR = 'YYYYMMDD.phis.' // GET_RES_EXT() 

      ! Replace date tokens
      CALL EXPAND_DATE( GEOS_DIR, NYMD, 000000 )
      CALL EXPAND_DATE( PHIS_STR, NYMD, 000000 )

      ! Location of zipped A-3 file in data dir
      PHIS_FILE_GZ = TRIM( DATA_DIR   ) // TRIM( GEOS_DIR   ) //
     &               TRIM( PHIS_STR   ) // TRIM( ZIP_SUFFIX )

      ! Location of unzipped A-3 file in temp dir
      PHIS_FILE    = TRIM( TEMP_DIR   ) // TRIM( PHIS_STR   )
         
      ! Remove A-3 files for this date from temp dir 
      REMOVE_DATE  = TRIM( REMOVE_CMD ) // ' '                // 
     &               TRIM( TEMP_DIR   ) // TRIM( PHIS_STR   )

      !=================================================================
      ! Define the foreground and background UNZIP commands
      !=================================================================

      ! Foreground unzip
      UNZIP_FG = TRIM( UNZIP_CMD ) // ' ' // TRIM( PHIS_FILE_GZ ) // 
     &           TRIM( REDIRECT  ) // ' ' // TRIM( PHIS_FILE    )  

      ! Background unzip
      UNZIP_BG  = TRIM( UNZIP_FG ) // TRIM( BACKGROUND )

      !=================================================================
      ! Define command to remove all PHIS files from the TEMP dir
      !=================================================================
      REMOVE_ALL = TRIM( REMOVE_CMD ) // ' '    // TRIM( TEMP_DIR  ) // 
     &             TRIM( WILD_CARD  ) //'.phis.'// TRIM( WILD_CARD ) 

      !=================================================================
      ! Perform an F90 system call to do the desired operation
      !=================================================================
      SELECT CASE ( TRIM( OPTION ) )
         
         ! Unzip PHIS field in the Unix foreground
         CASE ( 'unzip foreground' )
            WRITE( 6, 100 ) TRIM( PHIS_FILE_GZ )
            CALL SYSTEM( TRIM( UNZIP_FG ) )

         ! Unzip PHIS field in the Unix background
         CASE ( 'unzip background' )
            WRITE( 6, 100 ) TRIM( PHIS_FILE_GZ )
            CALL SYSTEM( TRIM( UNZIP_BG ) )

         ! Remove A-3 field for this date in temp dir
         CASE ( 'remove date' )
            WRITE( 6, 110 ) TRIM( PHIS_FILE )
            CALL SYSTEM( TRIM( REMOVE_DATE ) )
            
         ! Remove all A-3 fields in temp dir
         CASE ( 'remove all' )
            WRITE( 6, 120 ) TRIM( REMOVE_ALL )
            CALL SYSTEM( TRIM( REMOVE_ALL ) )

         ! Error -- bad option!
         CASE DEFAULT
            CALL ERROR_STOP( 'Invalid value for OPTION!', 
     &                       'UNZIP_PHIS_FIELDS (phis_read_mod.f)' )
            
      END SELECT

      ! FORMAT strings
 100  FORMAT( '     - Unzipping: ', a )
 110  FORMAT( '     - Removing: ', a )
 120  FORMAT( '     - About to execute command: ', a )

      ! Return to calling program
      END SUBROUTINE UNZIP_GCAP_FIELDS

!------------------------------------------------------------------------------

      SUBROUTINE OPEN_GCAP_FIELDS 
!
!******************************************************************************
!  Subroutine OPEN_GCAP_FIELDS opens the PHIS and LWI met fields file for 
!  date 2000/01/01. (swu, bmy, 2/1/06)
!  
!  NOTES:
!  (1 ) Adapted from OPEN_MET_FIELDS of "dao_read_mod.f" (bmy, 6/13/03)
!  (2 ) Now opens either zipped or unzipped files (bmy, 12/11/03)
!  (3 ) Now skips past the GEOS-4 ident string (bmy, 12/12/04)
!  (4 ) Now references "directory_mod.f" instead of CMN_SETUP.  Also now
!        references LUNZIP from "logical_mod.f".  Also now prevents EXPAND_DATE
!        from overwriting Y/M/D tokens in directory paths. (bmy, 7/20/04)
!  (5 ) Now use FILE_EXISTS from "file_mod.f" to determine if file unit IU_PH 
!        refers to a valid file on disk (bmy, 3/23/05)
!  (6 ) Now modified for GEOS-5 and GCAP met fields (swu, bmy, 5/25/05)
!******************************************************************************
!      
      ! References to F90 modules
      USE BPCH2_MOD,    ONLY : GET_RES_EXT
      USE DIRECTORY_MOD
      USE ERROR_MOD,    ONLY : ERROR_STOP
      USE LOGICAL_MOD,  ONLY : LUNZIP
      USE FILE_MOD,     ONLY : IOERROR, FILE_EXISTS
      USE TIME_MOD,     ONLY : EXPAND_DATE
      USE CMN_SIZE_MOD       ! Size parameters
!#if   !defined( DEVEL )
      USE FILE_MOD,     ONLY : IU_PH
!#endif

      ! Local variables
      LOGICAL, SAVE         :: FIRST = .TRUE.
      LOGICAL               :: IT_EXISTS, IS_OPEN
!#if   defined( DEVEL )
!      INTEGER               :: IU_PH
!#endif
      INTEGER               :: NYMD, NHMS
      INTEGER               :: IOS,  IUNIT
      CHARACTER(LEN=8)      :: IDENT
      CHARACTER(LEN=255)    :: GEOS_DIR
      CHARACTER(LEN=255)    :: PHIS_FILE
      CHARACTER(LEN=255)    :: PATH

      !=================================================================
      ! OPEN_PHIS_FIELDS begins here!
      !=================================================================
      
      ! Define date and hour
      NYMD = 20000101
      NHMS = 000000

      ! Open the A-3 file 0 GMT of each day, or on the first call
      IF ( NHMS == 000000 .or. FIRST ) THEN

         ! Strings for directory & filename
         GEOS_DIR  = TRIM( GCAP_DIR )
         PHIS_FILE = 'YYYYMMDD.phis.' // GET_RES_EXT() 

         ! Replace date tokens
         CALL EXPAND_DATE( GEOS_DIR,  NYMD, NHMS )
         CALL EXPAND_DATE( PHIS_FILE, NYMD, NHMS )

         ! If unzipping, open GEOS-1 file in TEMP dir
         ! If not unzipping, open GEOS-1 file in DATA dir
         IF ( LUNZIP ) THEN
            PATH = TRIM( TEMP_DIR ) // TRIM( PHIS_FILE )
         ELSE
            PATH = TRIM( DATA_DIR ) // 
     &             TRIM( GEOS_DIR ) // TRIM( PHIS_FILE )
         ENDIF

         ! Close previously opened A-3 file
	 INQUIRE(IU_PH, OPENED=IS_OPEN)
         IF(IS_OPEN) CLOSE( IU_PH )

         ! Make sure the file unit is valid before we open the file
         IF ( .not. FILE_EXISTS( IU_PH ) ) THEN
            CALL ERROR_STOP( 'Could not find file!', 
     &                       'OPEN_PHIS_FIELD (phis_read_mod.f)' )
         ENDIF

         ! Open the file
!#if defined( DEVEL )
!	 IU_PH = findFreeLUN()
!#endif
         OPEN( UNIT   = IU_PH,         FILE   = TRIM( PATH ),
     &         STATUS = 'OLD',         ACCESS = 'SEQUENTIAL',  
     &         FORM   = 'UNFORMATTED', IOSTAT = IOS )
               
         IF ( IOS /= 0 ) THEN
            CALL IOERROR( IOS, IU_PH, 'open_phis_fields:1' )
         ENDIF

         ! Echo info
         WRITE( 6, 100 ) TRIM( PATH )
 100     FORMAT( '     - Opening: ', a )
         
         ! Set the proper first-time-flag false
         FIRST = .FALSE.

         ! Skip past the GEOS-4 ident string
         READ( IU_PH, IOSTAT=IOS ) IDENT

         IF ( IOS /= 0 ) THEN
            CALL IOERROR( IOS, IU_PH, 'open_phis_fields:2' )
         ENDIF
      ENDIF

      ! Return to calling program
      END SUBROUTINE OPEN_GCAP_FIELDS

!------------------------------------------------------------------------------

      SUBROUTINE GET_GCAP_FIELDS
!
!******************************************************************************
!  Subroutine GET_GCAP_FIELDS is a wrapper for routine READ_PHIS.  This routine
!  calls READ_PHIS properly for reading PHIS fields from GEOS-1, GEOS-STRAT, 
!  GEOS-3, or GEOS-4 met data sets at the START of a GEOS-CHEM run. 
!  (bmy, swu, 2/1/06)
!
!  NOTES:
!  (1 ) Now also read LWI_GISS for GCAP met fields (swu, bmy, 5/25/05)
!******************************************************************************
! 
      ! References to F90 modules
      USE DAO_MOD, ONLY : PHIS, LWI_GISS

      ! Local variables
      INTEGER          :: NYMD, NHMS 

      !=================================================================
      ! GET_PHIS_FIELD begins here!
      !=================================================================
      
      ! Date and time
      NYMD = 20000101
      NHMS = 000000

      ! For GCAP met fields: read PHIS and LWI_GISS
      CALL READ_GCAP( NYMD=NYMD, NHMS=NHMS, PHIS=PHIS, LWI=LWI_GISS )

      ! Return to calling program
      END SUBROUTINE GET_GCAP_FIELDS

!---------------------------------------------------------------------------

      FUNCTION CHECK_TIME( XYMD, XHMS, NYMD, NHMS ) RESULT( ITS_TIME )
!
!******************************************************************************
!  Function CHECK_TIME checks to see if the timestamp of the A-3 field just
!  read from disk matches the current time.  If so, then it's time to return
!  the A-3 field to the calling program. (bmy, 6/16/03)
!  
!  Arguments as Input:
!  ============================================================================
!  (1 ) XYMD (REAL*4 or INTEGER) : (YY)YYMMDD timestamp for A-3 field in file
!  (2 ) XHMS (REAL*4 or INTEGER) : HHMMSS     timestamp for A-3 field in file
!  (3 ) NYMD (INTEGER          ) : YYYYMMDD   at which A-3 field is to be read
!  (4 ) NHMS (INTEGER          ) : HHMMSS     at which A-3 field is to be read
!
!  NOTES:
!******************************************************************************
!
      USE CMN_SIZE_MOD

      ! Arguments 
      INTEGER, INTENT(IN) :: XYMD, XHMS, NYMD, NHMS
      
      ! Function value
      LOGICAL             :: ITS_TIME

      !=================================================================
      ! GEOS-3, GEOS-4: XYMD and XHMS are integers
      !=================================================================
      IF ( XYMD == NYMD .AND. XHMS == NHMS ) THEN
         ITS_TIME = .TRUE.
      ELSE
         ITS_TIME = .FALSE.
      ENDIF

      ! Return to calling program
      END FUNCTION CHECK_TIME

!------------------------------------------------------------------------------

      SUBROUTINE READ_GCAP( NYMD, NHMS, PHIS, LWI )
!
!******************************************************************************
!  Subroutine READ_PHIS reads DAO PHIS (surface geopotential heights) field 
!  from disk.  PHIS is an I-6 field, but is time-independent.  Thus READ_PHIS
!  only needs to be called once at the beginning of the model run.
!  (bmy, swu, 2/1/06)
!
!  Arguments as input:
!  ============================================================================
!  (1 ) NYMD   : YYMMDD
!  (2 ) NHMS   :  and HHMMSS of PHIS field to be read from disk
!
!  Arguments as output:
!  ============================================================================
!  (3 ) PHIS   : DAO field for surface geopotential height (= g0 * m)
!                in units of m^2 / s^2, where g0 = 9.8 m / s^2.
!
!  NOTES:
!  (1 ) Adapted from READ_PHIS from "dao_read_mod.f" (bmy, 6/16/03)
!  (2 ) Now use function TIMESTAMP_STRING from "time_mod.f" for formatted 
!        date/time output. (bmy, 10/28/03)
!  (3 ) Now also read LWI_GISS for GCAP met fields.  Added optional variable
!        LWI to the arg list. (swu, bmy, 5/25/05)
!******************************************************************************
!
      ! References to F90 modules
      USE DIAG_MOD,     ONLY : AD67
      USE FILE_MOD,     ONLY : IOERROR, IU_PH
      USE TIME_MOD,     ONLY : TIMESTAMP_STRING
      USE TRANSFER_MOD, ONLY : TRANSFER_2D

      USE CMN_SIZE_MOD   ! Size parameters
      USE CMN_DIAG_MOD   ! ND67
      USE CMN_GCTM_MOD   ! g0

      ! Arguments
      INTEGER, INTENT(IN)            :: NYMD, NHMS
      REAL*8,  INTENT(OUT)           :: PHIS(IIPAR,JJPAR) 
      REAL*8,  INTENT(OUT), OPTIONAL :: LWI(IIPAR,JJPAR) 

      ! Local Variables
      INTEGER                        :: NFOUND, IOS
      INTEGER                        :: XYMD,   XHMS 
      REAL*4                         :: Q2(IIPAR,JJPAR)
      CHARACTER(LEN=8)               :: NAME
      CHARACTER(LEN=16)              :: STAMP

      ! Number of fields in the file
      INTEGER, PARAMETER             :: N_PHIS = 2

      !=================================================================
      ! READ_PHIS begins here!
      !=================================================================

      ! Zero number of PHIS fields we have found
      NFOUND = 0

      !=================================================================
      ! Read PHIS field from disk
      !=================================================================      
      DO

         ! PHIS field name
         READ( IU_PH, IOSTAT=IOS ) NAME

         ! IOS < 0: EOF, but make sure we have found everything
         IF ( IOS < 0 ) THEN
            CALL GCAP_CHECK( NFOUND, N_PHIS )
            EXIT
         ENDIF

         ! IOS > 0: True I/O error
         IF ( IOS > 0 ) CALL IOERROR( IOS, IU_PH, 'read_phis:1' )

         ! CASE statement for met fields
         SELECT CASE ( TRIM( NAME ) )

            !---------------------------------
            ! PHIS: geopotential heights
            !---------------------------------
            CASE ( 'PHIS' ) 
               READ( IU_PH, IOSTAT=IOS ) XYMD, XHMS, Q2
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_PH, 'read_phis:2' )

               IF ( CHECK_TIME( XYMD, XHMS, NYMD, NHMS ) ) THEN
                  CALL TRANSFER_2D( Q2, PHIS )
                  NFOUND = NFOUND + 1
               ENDIF

            !---------------------------------
            ! LWI_GISS: GCAP land/water flags
            !---------------------------------
            CASE ( 'LWI', 'LWI_GISS' ) 
               READ( IU_PH, IOSTAT=IOS ) XYMD, XHMS, Q2
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_PH, 'read_phis:3' )
               
               IF ( CHECK_TIME( XYMD, XHMS, NYMD, NHMS ) ) THEN
                  IF ( PRESENT( LWI ) ) CALL TRANSFER_2D( Q2, LWI )
                  NFOUND = NFOUND + 1
               ENDIF

            !--------------------------------
            ! Field not found
            !--------------------------------
            CASE DEFAULT
               WRITE( 6, '(a)' ) 'Searching for next field!'
         END SELECT

         !==============================================================
         ! If we have found all the fields for this time, then exit 
         ! the loop.  Otherwise, go on to the next iteration.
         !==============================================================
         IF ( CHECK_TIME( XYMD, XHMS, NYMD, NHMS ) .AND. 
     &        NFOUND == N_PHIS ) THEN
            STAMP = TIMESTAMP_STRING( NYMD, NHMS )
            WRITE( 6, 200 ) STAMP
 200        FORMAT( '     - Found PHIS met fields for ', a )
            EXIT
         ENDIF                  
      ENDDO

      !=================================================================
      ! Divide PHIS by 9.8 m / s^2 to obtain surface heights in meters. 
      !
      ! ND67 diagnostic: Accumulating DAO surface fields:
      ! Field #15 in the ND67 diagnostic is the geopotential heights
      !=================================================================
      PHIS = PHIS / g0

      IF ( ND67 > 0 ) THEN
         AD67(:,:,15) = AD67(:,:,15) + PHIS
      ENDIF  

      ! Since we only read PHIS at the start of the run,
      ! close the file unit (bmy, 6/16/03)
      CLOSE( IU_PH )

      ! Return to calling program      
      END SUBROUTINE READ_GCAP

!------------------------------------------------------------------------------

      SUBROUTINE GCAP_CHECK( NFOUND, N_PHIS )
!
!******************************************************************************
!  Subroutine PHIS_CHECK prints an error message if not all of the A-3 met 
!  fields are found.  The run is also terminated. (bmy, 10/27/00, 6/16/03)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) NFOUND (INTEGER) : # of met fields read from disk
!  (2 ) N_PHIS (INTEGER) : # of met fields expected to be read from disk
!
!  NOTES
!  (1 ) Adapted from DAO_CHECK from "dao_read_mod.f" (bmy, 6/16/03)
!******************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD, ONLY : GEOS_CHEM_STOP

      ! Arguments
      INTEGER, INTENT(IN) :: NFOUND, N_PHIS

      !=================================================================
      ! PHIS_CHECK begins here!
      !=================================================================
      IF ( NFOUND /= N_PHIS ) THEN
         WRITE( 6, '(a)' ) REPEAT( '=', 79 )
         WRITE( 6, '(a)' ) 'ERROR -- not enough PHIS fields found!'      

         WRITE( 6, 120   ) N_PHIS, NFOUND
 120     FORMAT( 'There are ', i2, ' fields but only ', i2 ,
     &           ' were found!' )

         WRITE( 6, '(a)' ) '### STOP in PHIS_CHECK (dao_read_mod.f)'
         WRITE( 6, '(a)' ) REPEAT( '=', 79 )

         ! Deallocate arrays and stop (bmy, 10/15/02)
         CALL GEOS_CHEM_STOP
      ENDIF

      ! Return to calling program
      END SUBROUTINE GCAP_CHECK

!------------------------------------------------------------------------------

      ! End of module
      END MODULE GCAP_READ_MOD
