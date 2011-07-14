! $ Id: xtra_read_mod.f v2.2 2005/4/20 21:17:00 tmf Exp $
      MODULE XTRA_READ_MOD
!
!******************************************************************************
!  Module XTRA_READ_MOD contains routines that unzip, open, and read the
!  GEOS-CHEM XTRA (avg 3-hour) met fields from disk. (dsa, tmf, bmy, 10/20/05)
! 
!  Module Routines:
!  =========================================================================
!  (1 ) UNZIP_XTRA_FIELDS : Unzips & copies met field files to a temp dir
!  (2 ) DO_OPEN_XTRA      : Returns TRUE if it's time to read XTRA fields
!  (3 ) OPEN_XTRA_FIELDS  : Opens met field files residing in the temp dir
!  (4 ) GET_XTRA_FIELDS   : Wrapper for routine READ_XTRA
!  (5 ) GET_N_XTRA        : Returns # of XTRA fields for each DAO data set 
!  (6 ) CHECK_TIME        : Tests if XTRA met field timestamps = current time
!  (7 ) READ_XTRA         : Reads XTRA fields from disk
!  (8 ) XTRA_CHECK        : Checks if we have found all of the XTRA fields
! 
!  GEOS-CHEM modules referenced by xtra_read_mod.f
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
!  (10) unix_cmds_mod.f   : Module w/ Unix commands for unzipping etc.
!
!  NOTES:
!******************************************************************************
!
      IMPLICIT NONE

      !=================================================================
      ! MODULE PRIVATE DECLARATIONS -- keep certain internal variables 
      ! and routines from being seen outside "xtra_read_mod.f"
      !=================================================================

      ! Make everything PRIVATE ...
      PRIVATE

      ! ... except these routines
      PUBLIC  :: GET_XTRA_FIELDS
      PUBLIC  :: OPEN_XTRA_FIELDS
      PUBLIC  :: UNZIP_XTRA_FIELDS

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement 
      !=================================================================
      CONTAINS

!------------------------------------------------------------------------------

      SUBROUTINE UNZIP_XTRA_FIELDS( OPTION, NYMD )
!
!*****************************************************************************
!  Subroutine UNZIP_XTRA_FIELDS invokes a FORTRAN system call to uncompress
!  GEOS-CHEM GEOS-3 XTRA met field files and store the uncompressed data in a 
!  temporary directory, where GEOS-CHEM can read them.  The original data 
!  files are not disturbed.  (dsa, tmf, bmy, 10/20/05)
!
!  Arguments as input:
!  ===========================================================================
!  (1 ) OPTION (CHAR*(*)) : Option
!  (2 ) NYMD   (INTEGER ) : YYYYMMDD of XTRA file to be unzipped (optional)
!
!  NOTES:
!*****************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD,     ONLY : GET_RES_EXT
      USE DIRECTORY_MOD, ONLY : DATA_DIR,   GEOS_3_DIR, TEMP_DIR
      USE ERROR_MOD,     ONLY : ERROR_STOP
      USE TIME_MOD,      ONLY : EXPAND_DATE
      USE UNIX_CMDS_MOD, ONLY : BACKGROUND, REDIRECT,   REMOVE_CMD
      USE UNIX_CMDS_MOD, ONLY : UNZIP_CMD,  WILD_CARD,  ZIP_SUFFIX 

#     include "CMN_SIZE"

      ! Arguments
      CHARACTER(LEN=*),  INTENT(IN) :: OPTION
      INTEGER, OPTIONAL, INTENT(IN) :: NYMD

      ! Local variables
      CHARACTER(LEN=255)            :: XTRA_STR,     GEOS_DIR
      CHARACTER(LEN=255)            :: XTRA_FILE_GZ, XTRA_FILE
      CHARACTER(LEN=255)            :: UNZIP_BG,   UNZIP_FG
      CHARACTER(LEN=255)            :: REMOVE_ALL, REMOVE_DATE

      !=================================================================
      ! UNZIP_XTRA_FIELDS begins here!
      !=================================================================
      IF ( PRESENT( NYMD ) ) THEN
      
         ! Strings for directory & filename
         GEOS_DIR = TRIM( GEOS_3_DIR )
         XTRA_STR = 'YYYYMMDD.xtra.' // GET_RES_EXT()

         ! Replace date tokens
         CALL EXPAND_DATE( GEOS_DIR, NYMD, 000000 )
         CALL EXPAND_DATE( XTRA_STR, NYMD, 000000 )

         ! Location of zipped XTRA file in data dir 
         XTRA_FILE_GZ  = TRIM( DATA_DIR   ) // TRIM( GEOS_DIR   ) // 
     &                   TRIM( XTRA_STR   ) // TRIM( ZIP_SUFFIX )

         ! Location of unzipped XTRA file in temp dir 
         XTRA_FILE     = TRIM( TEMP_DIR   ) // TRIM( XTRA_STR   )     
         
         ! Remove XTRA files for this date from temp dir
         REMOVE_DATE   = TRIM( REMOVE_CMD ) // ' '                // 
     &                   TRIM( TEMP_DIR   ) // TRIM( XTRA_STR   )   

         !==============================================================
         ! Define the foreground and background UNZIP commands
         !==============================================================

         ! Foreground unzip
         UNZIP_FG = TRIM( UNZIP_CMD ) // ' ' // TRIM( XTRA_FILE_GZ ) // 
     &              TRIM( REDIRECT  ) // ' ' // TRIM( XTRA_FILE    )  

         ! Background unzip
         UNZIP_BG  = TRIM( UNZIP_FG ) // TRIM( BACKGROUND )
      ENDIF

      !=================================================================
      ! Define command to remove all XTRA files from the TEMP dir
      !=================================================================
      REMOVE_ALL = TRIM( REMOVE_CMD ) // ' '  // TRIM( TEMP_DIR   ) // 
     &             TRIM( WILD_CARD  ) // '.xtra.' // TRIM( WILD_CARD  ) 

      !=================================================================
      ! Perform an F90 system call to do the desired operation
      !=================================================================
      SELECT CASE ( TRIM( OPTION ) )
         
         ! Unzip XTRA fields in the Unix foreground
         CASE ( 'unzip foreground' )
            WRITE( 6, 100 ) TRIM( XTRA_FILE_GZ )
            CALL SYSTEM( TRIM( UNZIP_FG ) )

         ! Unzip XTRA fields in the Unix background
         CASE ( 'unzip background' )
            WRITE( 6, 100 ) TRIM( XTRA_FILE_GZ )
            CALL SYSTEM( TRIM( UNZIP_BG ) )

         ! Remove XTRA field for this date in temp dir
         CASE ( 'remove date' )
            WRITE( 6, 110 ) TRIM( XTRA_FILE )
            CALL SYSTEM( TRIM( REMOVE_DATE ) )
            
         ! Remove all XTRA fields in temp dir
         CASE ( 'remove all' )
            WRITE( 6, 120 ) TRIM( REMOVE_ALL )
            CALL SYSTEM( TRIM( REMOVE_ALL ) )

         ! Error -- bad option!
         CASE DEFAULT
            CALL ERROR_STOP( 'Invalid value for OPTION!', 
     &                       'UNZIP_XTRA_FIELDS (xtra_read_mod.f)' )
            
      END SELECT

      ! FORMAT strings
 100  FORMAT( '     - Unzipping: ', a )
 110  FORMAT( '     - Removing: ', a )
 120  FORMAT( '     - About to execute command: ', a )

      ! Return to calling program
      END SUBROUTINE UNZIP_XTRA_FIELDS

!------------------------------------------------------------------------------

      FUNCTION DO_OPEN_XTRA( NYMD, NHMS ) RESULT( DO_OPEN )
!
!******************************************************************************
!  Function DO_OPEN_XTRA returns TRUE if is time to open the XTRA met field 
!  file or FALSE otherwise.  This prevents us from opening a file which has 
!  already been opened. (dsa, tmf, bmy, 10/20/05)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) NYMD (INTEGER) : YYYYMMDD 
!  (2 ) NHMS (INTEGER) :  and HHMMSS to be tested for A-3 file open
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
      ! DO_OPEN_XTRA begins here!
      !=================================================================

      ! Initialize
      DO_OPEN = .FALSE.

      ! Return if we have already opened the file
      IF ( NYMD == LASTNYMD .and. NHMS == LASTNHMS ) THEN
         DO_OPEN = .FALSE. 
         GOTO 999
      ENDIF

      ! Open XTRA file if it's 00:00 GMT, or on the first call
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
      END FUNCTION DO_OPEN_XTRA

!------------------------------------------------------------------------------

      SUBROUTINE OPEN_XTRA_FIELDS( NYMD, NHMS )
!
!******************************************************************************
!  Subroutine OPEN_XTRA_FIELDS opens the XTRA met fields file for date NYMD 
!  and time NHMS. (dsa, tmf, bmy, 10/20/05)
!  
!  Arguments as input:
!  ===========================================================================
!  (1 ) NYMD (INTEGER) : YYYYMMDD
!  (2 ) NHMS (INTEGER) :  and HHMMSS timestamps for XTRA file
!
!  NOTES:
!******************************************************************************
!      
      ! References to F90 modules
      USE BPCH2_MOD,     ONLY : GET_RES_EXT
      USE DIRECTORY_MOD, ONLY : DATA_DIR, GEOS_3_DIR, TEMP_DIR
      USE ERROR_MOD,     ONLY : ERROR_STOP
      USE LOGICAL_MOD,   ONLY : LUNZIP
      USE FILE_MOD,      ONLY : IU_XT, IOERROR, FILE_EXISTS
      USE TIME_MOD,      ONLY : EXPAND_DATE

#     include "CMN_SIZE"      ! Size parameters

      ! Arguments
      INTEGER, INTENT(IN)    :: NYMD, NHMS

      ! Local variables
      LOGICAL                :: DO_OPEN
      LOGICAL                :: IT_EXISTS
      INTEGER                :: IOS
      CHARACTER(LEN=8)       :: IDENT
      CHARACTER(LEN=255)     :: XTRA_FILE
      CHARACTER(LEN=255)     :: GEOS_DIR
      CHARACTER(LEN=255)     :: PATH

      !=================================================================
      ! OPEN_XTRA_FIELDS begins here!
      !=================================================================

      ! Open XTRA fields at the proper time, or on the first call
      IF ( DO_OPEN_XTRA( NYMD, NHMS ) ) THEN

         ! Strings for directory & filename
         GEOS_DIR = TRIM( GEOS_3_DIR )
         XTRA_FILE  = 'YYYYMMDD.xtra.' // GET_RES_EXT()

         ! Replace date tokens
         CALL EXPAND_DATE( XTRA_FILE,  NYMD, NHMS )
         CALL EXPAND_DATE( GEOS_DIR, NYMD, NHMS )

         ! If unzipping, open GEOS-4 file in TEMP dir
         ! If not unzipping, open GEOS-4 file in DATA dir
         IF ( LUNZIP ) THEN
            PATH = TRIM( TEMP_DIR ) // TRIM( XTRA_FILE )
         ELSE
            PATH = TRIM( DATA_DIR ) // 
     &             TRIM( GEOS_DIR ) // TRIM( XTRA_FILE )
         ENDIF

         ! Close previously opened XTRA file
         CLOSE( IU_XT )

         ! Make sure the file unit is valid before we open the file
         IF ( .not. FILE_EXISTS( IU_XT ) ) THEN
            CALL ERROR_STOP( 'Could not find file!', 
     &                       'OPEN_XTRA_FIELDS (xtra_read_mod.f)' )
         ENDIF

         ! Open the file
         OPEN( UNIT   = IU_XT,         FILE   = TRIM( PATH ),
     &         STATUS = 'OLD',         ACCESS = 'SEQUENTIAL',  
     &         FORM   = 'UNFORMATTED', IOSTAT = IOS )
               
         IF ( IOS /= 0 ) THEN
            CALL IOERROR( IOS, IU_XT, 'open_xtra_fields:1' )
         ENDIF

         ! Echo info
         WRITE( 6, 100 ) TRIM( PATH )
 100     FORMAT( '     - Opening: ', a )
         
      ENDIF

      ! Return to calling program
      END SUBROUTINE OPEN_XTRA_FIELDS

!------------------------------------------------------------------------------

      SUBROUTINE GET_XTRA_FIELDS( NYMD, NHMS )
!
!******************************************************************************
!  Subroutine GET_XTRA_FIELDS is a wrapper for routine READ_XTRA.  
!  GET_XTRA_FIELDS calls READ_XTRA properly for reading the GEOS-3 met data
!  set. (dsa, tmf, bmy, 10/20/05)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) NYMD (INTEGER) : YYYYMMDD
!  (2 ) NHMS (INTEGER) :  and HHMMSS of XTRA fields to be read from disk
!
!  NOTES:
!  (1 ) Now extract only PARDR, PARDF for MEGAN biogenics inventory and SNOW
!        for dust emissions from GEOS3. (tmf, 6/23/05)
!******************************************************************************
!
      ! References to F90 modules
      USE DAO_MOD,    ONLY : PARDR, PARDF, SNOW
      USE FILE_MOD,   ONLY : IU_XT

#     include "CMN_SIZE"   ! Size parameters

      ! Arguments
      INTEGER, INTENT(IN) :: NYMD, NHMS 

      ! Local variables
      INTEGER, SAVE       :: LASTNYMD = -1, LASTNHMS = -1

      !=================================================================
      ! GET_XTRA_FIELDS begins here!
      !=================================================================

      ! Skip over previously-read XTRA fields
      IF ( NYMD == LASTNYMD .and. NHMS == LASTNHMS ) THEN
         WRITE( 6, 100 ) NYMD, NHMS
 100     FORMAT( '     - XTRA met fields for NYMD, NHMS = ', 
     &           i8.8, 1x, i6.6, ' have been read already' ) 
         RETURN
      ENDIF

      ! Read PARDR, PARDF fields
      CALL READ_XTRA( NYMD=NYMD,   NHMS=NHMS,
     &                PARDR=PARDR, PARDF=PARDF, SNOW=SNOW  )

      ! Save NYMD, NHMS for next call
      LASTNYMD = NYMD
      LASTNHMS = NHMS

      ! Return to MAIN program
      END SUBROUTINE GET_XTRA_FIELDS

!---------------------------------------------------------------------------

      FUNCTION CHECK_TIME( XYMD, XHMS, NYMD, NHMS ) RESULT( ITS_TIME )
!
!******************************************************************************
!  Function CHECK_TIME checks to see if the timestamp of the XTRA field just
!  read from disk matches the current time.  If so, then it's time to return
!  the XTRA field to the calling program. (dsa, tmf, bmy, 10/20/05)
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
#     include "CMN_SIZE"

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

!-----------------------------------------------------------------------------

      SUBROUTINE READ_XTRA( NYMD,  NHMS, 
     &                      PARDR, PARDF,  TSKIN, LAI, 
     &                      EVAP,  RADLWG, SNOW        )
!
!******************************************************************************
!  Subroutine READ_XTRA reads GEOS-3 XTRA (3-hr avg) fields from disk.
!  (dsa, tmf, bmy, 10/20/05)
! 
!  Arguments as input:
!  ============================================================================
!  (1 ) NYMD    : YYYYMMDD
!  (2 ) NHMS    :  and HHMMSS of XTRA met fields to be accessed 
!
!  XTRA Met Fields as Output:
!  ============================================================================
!  (1 ) PARDF   : (2-D) GMAO Photosyn active diffuse radiation   [W/m2]
!  (2 ) PARDR   : (2-D) GMAO Photosyn active direct radiation    [W/m2]
!  (3 ) TSKIN   : (2-D) GMAO Surface ground/sea surface temp     [K]
!  (4 ) LAI     : (2-D) GMAO Leaf area indices                   [unitless]
!  (5 ) EVAP    : (2-D) GMAO Evaporation                         [mm/day]
!  (6 ) RADLWG  : (2-D) GMAO Net upward LW rad at the ground     [W/m2]
!  (7 ) SNOW    : (2-D) GMAO Snow cover (H2O equivalent)         [mm H2O]
!
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE DIAG_MOD,              ONLY : AD67
      USE FILE_MOD,              ONLY : IOERROR,     IU_XT
      USE TIME_MOD,              ONLY : SET_CT_XTRA, TIMESTAMP_STRING
      USE TRANSFER_MOD,          ONLY : TRANSFER_2D, TRANSFER_TO_1D

#     include "CMN_SIZE"              ! Size parameters
#     include "CMN_DIAG"              ! ND67

      ! Arguments
      INTEGER, INTENT(IN)            :: NYMD, NHMS
      REAL*8,  INTENT(OUT), OPTIONAL :: PARDR (IIPAR,JJPAR) 
      REAL*8,  INTENT(OUT), OPTIONAL :: PARDF (IIPAR,JJPAR) 
      REAL*8,  INTENT(OUT), OPTIONAL :: TSKIN (IIPAR,JJPAR) 
      REAL*8,  INTENT(OUT), OPTIONAL :: LAI   (IIPAR,JJPAR)
      REAL*8,  INTENT(OUT), OPTIONAL :: EVAP  (IIPAR,JJPAR)
      REAL*8,  INTENT(OUT), OPTIONAL :: RADLWG(IIPAR,JJPAR)
      REAL*8,  INTENT(OUT), OPTIONAL :: SNOW  (IIPAR,JJPAR)

      ! Local Variables
      INTEGER                        :: I, IJLOOP, IOS, J 
      INTEGER                        :: N_XTRA, NFOUND 
      INTEGER                        :: XYMD, XHMS
      REAL*4                         :: Q2(IIPAR,JJPAR)
      CHARACTER(LEN=8)               :: NAME
      CHARACTER(LEN=16)              :: STAMP
    
      !=================================================================
      ! READ_XTRA begins here!      
      !=================================================================

      ! Get the number of XTRA fields stored in this data set
      N_XTRA = 7

      ! Zero the number of A-3 fields that we have found
      NFOUND = 0

      !=================================================================
      ! Read the XTRA fields from disk
      !=================================================================
      DO

         ! Read the XTRA field name
         READ( IU_XT, IOSTAT=IOS ) NAME

         ! End of file test -- make sure we have found all fields
         IF ( IOS < 0 ) THEN
            CALL XTRA_CHECK( NFOUND, N_XTRA )
            EXIT
         ENDIF

         ! IOS > 0: True I/O error; stop w/ err msg
         IF ( IOS > 0 ) CALL IOERROR( IOS, IU_XT, 'read_xtra:1' )

         ! CASE statement for XTRA fields
         SELECT CASE ( TRIM( NAME ) )

            !--------------------------------
            ! PARDR: Photosyn active direct radiation
            !--------------------------------
            CASE ( 'PARDR' ) 
               READ( IU_XT, IOSTAT=IOS ) XYMD, XHMS, Q2
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_XT, 'read_xtra:2' )

               IF ( CHECK_TIME( XYMD, XHMS, NYMD, NHMS ) ) THEN
                  IF ( PRESENT( PARDR ) ) CALL TRANSFER_2D( Q2, PARDR )
                  NFOUND = NFOUND + 1
               ENDIF

            !--------------------------------
            ! PARDF: Photosyn active diffuse radiation
            !--------------------------------
            CASE ( 'PARDF' ) 
               READ( IU_XT, IOSTAT=IOS ) XYMD, XHMS, Q2
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_XT, 'read_xtra:3' )

               IF ( CHECK_TIME( XYMD, XHMS, NYMD, NHMS ) ) THEN
                  IF ( PRESENT( PARDF ) ) CALL TRANSFER_2D( Q2, PARDF )
                  NFOUND = NFOUND + 1
               ENDIF

            !--------------------------------
            ! TSKIN, TGROUND: Surface ground/sea surface temp
            !--------------------------------
            CASE ( 'TSKIN', 'TGROUND' ) 
               READ( IU_XT, IOSTAT=IOS ) XYMD, XHMS, Q2
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_XT, 'read_xtra:4' )

               IF ( CHECK_TIME( XYMD, XHMS, NYMD, NHMS ) ) THEN
                  IF ( PRESENT( TSKIN ) ) CALL TRANSFER_2D( Q2, TSKIN )
                  NFOUND = NFOUND + 1
               ENDIF

            !--------------------------------
            ! EVAP: Evaporation
            !--------------------------------
            CASE ( 'EVAP' ) 
               READ( IU_XT, IOSTAT=IOS ) XYMD, XHMS, Q2
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_XT, 'read_xtra:5' )

               IF ( CHECK_TIME( XYMD, XHMS, NYMD, NHMS ) ) THEN
                  IF ( PRESENT( EVAP ) ) CALL TRANSFER_2D( Q2, EVAP )
                  NFOUND = NFOUND + 1
               ENDIF

            !--------------------------------
            ! LAI: Leaf area indices
            !--------------------------------
            CASE ( 'LAI' ) 
               READ( IU_XT, IOSTAT=IOS ) XYMD, XHMS, Q2
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_XT, 'read_xtra:6' )

               IF ( CHECK_TIME( XYMD, XHMS, NYMD, NHMS ) ) THEN
                  IF ( PRESENT( LAI ) ) CALL TRANSFER_2D( Q2, LAI )
                  NFOUND = NFOUND + 1
               ENDIF

            !--------------------------------
            ! SNOW: Snow cover (H2O equivalent)
            !--------------------------------
            CASE ( 'SNOW' ) 
               READ( IU_XT, IOSTAT=IOS ) XYMD, XHMS, Q2
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_XT, 'read_xtra:7' )

               IF ( CHECK_TIME( XYMD, XHMS, NYMD, NHMS ) ) THEN
                  IF ( PRESENT( SNOW ) ) CALL TRANSFER_2D( Q2, SNOW )
                  NFOUND = NFOUND + 1
               ENDIF

            !--------------------------------
            ! RADLWG: Net upward LW rad at the ground
            !--------------------------------
            CASE ( 'RADLWG' ) 
               READ( IU_XT, IOSTAT=IOS ) XYMD, XHMS, Q2
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_XT, 'read_xtra:8' )

               IF ( CHECK_TIME( XYMD, XHMS, NYMD, NHMS ) ) THEN
                  IF ( PRESENT( RADLWG ) ) 
     &               CALL TRANSFER_2D( Q2, RADLWG )
                  NFOUND = NFOUND + 1
               ENDIF

         END SELECT
               
         !==============================================================
         ! If we have found all the fields for this time, then exit 
         ! the loop.  Otherwise, go on to the next iteration.
         !==============================================================
         IF ( CHECK_TIME( XYMD, XHMS, NYMD, NHMS ) .and. 
     &        NFOUND == N_XTRA ) THEN 
            STAMP = TIMESTAMP_STRING( NYMD, NHMS )
            WRITE( 6, 210 ) NFOUND, STAMP
 210        FORMAT( '     - Found all ', i3, 
     &                    ' XTRA met fields for ', a )
            EXIT
         ENDIF
      ENDDO

      !=================================================================
      ! ND67 diagnostic: A-3 surface fields:
      !
      ! (19) TSKIN  : Ground/sea surface temp.           [hPa]  
      ! (20) PARDF  : Photosyn active diffuse radiation  [W/m2]
      ! (21) PARDR  : Photosyn active direct  radiation  [W/m2]
      !=================================================================
      IF ( ND67 > 0 ) THEN
         IF ( PRESENT( TSKIN   ) ) AD67(:,:,19) = AD67(:,:,19) + TSKIN
         IF ( PRESENT( PARDF   ) ) AD67(:,:,20) = AD67(:,:,20) + PARDF
         IF ( PRESENT( PARDR   ) ) AD67(:,:,21) = AD67(:,:,21) + PARDR
      ENDIF
         
      ! Increment # of times READ_XTRA is called
      CALL SET_CT_XTRA( INCREMENT=.TRUE. )

      ! Return to calling program
      END SUBROUTINE READ_XTRA

!------------------------------------------------------------------------------

      SUBROUTINE XTRA_CHECK( NFOUND, N_XTRA )
!
!******************************************************************************
!  Subroutine XTRA_CHECK prints an error message if not all of the XTRA met 
!  fields are found.  The run is also terminated. (bmy, 10/27/00, 6/23/03)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) NFOUND (INTEGER) : # of XTRA met fields read from disk
!  (2 ) N_XTRA   (INTEGER) : # of XTRA met fields expected to be read from disk
!
!  NOTES
!  (1 ) Adapted from DAO_CHECK from "dao_read_mod.f" (bmy, 6/23/03)
!******************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD, ONLY : GEOS_CHEM_STOP

      ! Arguments
      INTEGER, INTENT(IN) :: NFOUND, N_XTRA

      !=================================================================
      ! XTRA_CHECK begins here!
      !=================================================================
      IF ( NFOUND /= N_XTRA ) THEN
         WRITE( 6, '(a)' ) REPEAT( '=', 79 )
         WRITE( 6, '(a)' ) 'ERROR -- not enough XTRA fields found!'      

         WRITE( 6, 120   ) N_XTRA, NFOUND
 120     FORMAT( 'There are ', i2, ' fields but only ', i2 ,
     &           ' were found!' )

         WRITE( 6, '(a)' ) '### STOP in XTRA_CHECK (xtra_read_mod.f)'
         WRITE( 6, '(a)' ) REPEAT( '=', 79 )

         ! Deallocate arrays and stop (bmy, 10/15/02)
         CALL GEOS_CHEM_STOP
      ENDIF

      ! Return to calling program
      END SUBROUTINE XTRA_CHECK

!------------------------------------------------------------------------------

      END MODULE XTRA_READ_MOD
