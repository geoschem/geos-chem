! $Id: a6_read_mod.f,v 1.3 2003/12/11 21:54:08 bmy Exp $
      MODULE A6_READ_MOD
!
!******************************************************************************
!  Module A6_READ_MOD contains subroutines that unzip, open, and read
!  GEOS-CHEM A-3 (avg 3-hour) met fields from disk. (bmy, 6/19/03, 12/11/03)
! 
!  Module Routines:
!  ============================================================================
!  (1 ) UNZIP_A6_FIELDS : Unzips & copies met field files to a temp dir
!  (2 ) DO_OPEN_A6      : Returns TRUE if it's time to open A-6 fields
!  (3 ) OPEN_A6_FIELDS  : Opens met field files residing in the temp dir
!  (4 ) GET_A6_FIELDS   : Wrapper for routine READ_A6
!  (5 ) MAKE_CLDFRC     : Computes CFRAC from CLMO and CLRO for GEOS-STRAT
!  (6 ) GET_N_A6        : Returns # of A-6 fields for each DAO data set
!  (7 ) CHECK_TIME      : Tests if A-6 et field timestamps equal current time
!  (8 ) READ_A6         : Reads A-6 fields from disk
!  (9 ) A6_CHECK        : Checks if we have found all of the A-6 fields
! 
!  GEOS-CHEM modules referenced by a6_read_mod.f
!  ============================================================================
!  (1 ) bpch2_mod.f     : Module containing routines for binary punch file I/O
!  (2 ) dao_mod.f       : Module containing arrays for DAO met fields
!  (3 ) diag_mod.f      : Module containing GEOS-CHEM diagnostic arrays
!  (4 ) error_mod.f     : Module containing NaN and other error check routines
!  (5 ) file_mod.f      : Module containing file unit #'s and error checks
!  (6 ) time_mod.f      : Module containing routines for computing time & date
!  (7 ) transfer_mod.f  : Module containing routines to cast & resize arrays
!
!  NOTES:
!  (1 ) Adapted from "dao_read_mod.f" (bmy, 6/19/03)
!  (2 ) Now use TIMESTAMP_STRING for formatted output (bmy, 10/28/03)
!  (3 ) CLDFRC is now a 2-D array in MAKE_CLDFRC< GET_A6_FIELDS.  Also now
!        read from either zipped or unzipped files. (bmy, 12/9/03)
!******************************************************************************
!
      IMPLICIT NONE

      !=================================================================
      ! MODULE PRIVATE DECLARATIONS -- keep certain internal variables 
      ! and routines from being seen outside "a6_read_mod.f"
      !=================================================================

      ! PRIVATE module routines
      PRIVATE :: A6_CHECK, CHECK_TIME,  DO_OPEN_A6
      PRIVATE :: GET_N_A6, MAKE_CLDFRC, READ_A6

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement 
      !=================================================================
      CONTAINS

!------------------------------------------------------------------------------

      SUBROUTINE UNZIP_A6_FIELDS( OPTION, NYMD )
!
!*****************************************************************************
!  Subroutine UNZIP_A6_FIELDS invokes a FORTRAN system call to uncompress
!  GEOS-CHEM A-6 met field files and store the uncompressed data in a 
!  temporary directory, where GEOS-CHEM can read them.  The original data 
!  files are not disturbed.  (bmy, bdf, 6/15/98, 12/11/03)
!
!  Arguments as input:
!  ===========================================================================
!  (1 ) OPTION (CHAR*(*)) : Option
!  (2 ) NYMD   (INTEGER ) : YYYYMMDD of A-6 file to be unzipped (optional)
!
!  NOTES:
!  (1 ) Adapted from UNZIP_MET_FIELDS of "dao_read_mod.f" (bmy, 6/19/03)
!  (2 ) Directory information YYYY/MM or YYYYMM is now contained w/in 
!        GEOS_1_DIR, GEOS_S_DIR, GEOS_3_DIR, GEOS_4_DIR (bmy, 12/11/03)
!*****************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD, ONLY : GET_RES_EXT
      USE ERROR_MOD, ONLY : ERROR_STOP
      USE TIME_MOD,  ONLY : EXPAND_DATE

#     include "CMN_SIZE"
#     include "CMN_SETUP"

      ! Arguments
      CHARACTER(LEN=*),  INTENT(IN) :: OPTION
      INTEGER, OPTIONAL, INTENT(IN) :: NYMD

      ! Local variables
      CHARACTER(LEN=255)            :: A6_FILE_GZ, A6_FILE
      CHARACTER(LEN=255)            :: UNZIP_BG,   UNZIP_FG
      CHARACTER(LEN=255)            :: REMOVE_ALL, REMOVE_DATE

      !=================================================================
      ! UNZIP_A6_FIELDS begins here!
      !=================================================================
      IF ( PRESENT( NYMD ) ) THEN
      
#if   defined( GEOS_1 )

         ! Location of zipped A-3 file in data dir (GEOS-1)
         A6_FILE_GZ  = TRIM( DATA_DIR )      // TRIM( GEOS_1_DIR ) // 
!-----------------------------------------------------------------------------
! Prior to 12/11/03:
!     &                 'YYMM/YYMMDD.a6.'     // GET_RES_EXT()      // 
!-----------------------------------------------------------------------------
     &                 'YYMMDD.a6.'          // GET_RES_EXT()      // 
     &                 TRIM( ZIP_SUFFIX )

         ! Location of unzipped A-3 file in temp dir (GEOS-1)
         A6_FILE     = TRIM( TEMP_DIR )      // 'YYMMDD.a6.'       // 
     &                 GET_RES_EXT()
         
         ! Remove A-3 files for this date from temp dir (GEOS-1)
         REMOVE_DATE = TRIM( REMOVE_CMD )    // ' '                // 
     &                 TRIM( TEMP_DIR   )    // 'YYMMDD.a6.'       // 
     &                 GET_RES_EXT()  

#elif defined( GEOS_STRAT )

         ! Location of zipped A-3 file in data dir (GEOS-STRAT)
         A6_FILE_GZ  = TRIM( DATA_DIR )      // TRIM( GEOS_S_DIR ) // 
!-----------------------------------------------------------------------------
! Prior to 12/11/03:
!     &                 'YYMM/YYMMDD.a6.'     // GET_RES_EXT()      // 
!-----------------------------------------------------------------------------
     &                 'YYMMDD.a6.'          // GET_RES_EXT()      // 
     &                 TRIM( ZIP_SUFFIX )

         ! Location of unzipped A-3 file in temp dir (GEOS-STRAT)
         A6_FILE     = TRIM( TEMP_DIR )      // 'YYMMDD.a6.'       // 
     &                 GET_RES_EXT()

         ! Remove A-3 files for this date from temp dir (GEOS-STRAT)
         REMOVE_DATE = TRIM( REMOVE_CMD )    // ' '                // 
     &                 TRIM( TEMP_DIR   )    // 'YYMMDD.a6.'       // 
     &                 GET_RES_EXT()  

#elif defined( GEOS_3 )

         ! Location of zipped A-3 file in data dir (GEOS-3)
         A6_FILE_GZ  = TRIM( DATA_DIR )      // TRIM( GEOS_3_DIR ) // 
!-----------------------------------------------------------------------------
! Prior to 12/11/03:
!     &                 'YYYYMM/YYYYMMDD.a6.' // GET_RES_EXT()      // 
!-----------------------------------------------------------------------------
     &                 'YYYYMMDD.a6.'        // GET_RES_EXT()      // 
     &                 TRIM( ZIP_SUFFIX )

         ! Location of unzipped A-3 file in temp dir (GEOS-3)
         A6_FILE     = TRIM( TEMP_DIR )      // 'YYYYMMDD.a6.'     // 
     &                 GET_RES_EXT()

         ! Remove A-3 files for this date from temp dir (GEOS-3) 
         REMOVE_DATE = TRIM( REMOVE_CMD )    // ' '                // 
     &                 TRIM( TEMP_DIR   )    // 'YYYYMMDD.a6.'     // 
     &                 GET_RES_EXT()  

#elif defined( GEOS_4 )

         ! Location of zipped A-3 file in data dir (GEOS-4)
         A6_FILE_GZ  = TRIM( DATA_DIR )      // TRIM( GEOS_4_DIR ) // 
!-----------------------------------------------------------------------------
!     &                 'YYYYMM/YYYYMMDD.a6.' // GET_RES_EXT()      // 
!-----------------------------------------------------------------------------
     &                 'YYYYMMDD.a6.'        // GET_RES_EXT()      // 
     &                 TRIM( ZIP_SUFFIX )

         ! Location of unzipped A-3 file in temp dir (GEOS-4)
         A6_FILE     = TRIM( TEMP_DIR )      // 'YYYYMMDD.a6.'     // 
     &                 GET_RES_EXT()

         ! Remove A-3 files for this date from temp dir (GEOS-4)
         REMOVE_DATE = TRIM( REMOVE_CMD )    // ' '                // 
     &                 TRIM( TEMP_DIR   )    // 'YYYYMMDD.a6.'     // 
     &                 GET_RES_EXT()  

#endif

         !==============================================================
         ! Replace tokens in character variables w/ year & month values
         !==============================================================
         CALL EXPAND_DATE( A6_FILE_GZ,  NYMD, 000000 )
         CALL EXPAND_DATE( A6_FILE,     NYMD, 000000 )
         CALL EXPAND_DATE( REMOVE_DATE, NYMD, 000000 )

         !==============================================================
         ! Define the foreground and background UNZIP commands
         !==============================================================

         ! Foreground unzip
         UNZIP_FG = TRIM( UNZIP_CMD ) // ' ' // TRIM( A6_FILE_GZ ) // 
     &              TRIM( REDIRECT  ) // ' ' // TRIM( A6_FILE    )  

         ! Background unzip
         UNZIP_BG  = TRIM( UNZIP_FG ) // TRIM( BACKGROUND )
      ENDIF

      !=================================================================
      ! Define command to remove all A-6 files from the TEMP dir
      !=================================================================
      REMOVE_ALL = TRIM( REMOVE_CMD ) // ' '    // TRIM( TEMP_DIR  ) // 
     &             TRIM( WILD_CARD  ) // '.a6.' // TRIM( WILD_CARD ) 

      !=================================================================
      ! Perform an F90 system call to do the desired operation
      !=================================================================
      SELECT CASE ( TRIM( OPTION ) )
         
         ! Unzip A-3 fields in the Unix foreground
         CASE ( 'unzip foreground' )
            WRITE( 6, 100 ) TRIM( A6_FILE_GZ )
            CALL SYSTEM( TRIM( UNZIP_FG ) )

         ! Unzip A-3 fields in the Unix background
         CASE ( 'unzip background' )
            WRITE( 6, 100 ) TRIM( A6_FILE_GZ )
            CALL SYSTEM( TRIM( UNZIP_BG ) )

         ! Remove A-3 field for this date in temp dir
         CASE ( 'remove date' )
            WRITE( 6, 110 ) TRIM( A6_FILE )
            CALL SYSTEM( TRIM( REMOVE_DATE ) )
            
         ! Remove all A-3 fields in temp dir
         CASE ( 'remove all' )
            WRITE( 6, 120 ) TRIM( REMOVE_ALL )
            CALL SYSTEM( TRIM( REMOVE_ALL ) )

         ! Error -- bad option!
         CASE DEFAULT
            CALL ERROR_STOP( 'Invalid value for OPTION!', 
     &                       'UNZIP_A6_FIELDS (a6_read_mod.f)' )
            
      END SELECT

      ! FORMAT strings
 100  FORMAT( '     - Unzipping: ', a )
 110  FORMAT( '     - Removing: ', a )
 120  FORMAT( '     - About to execute command: ', a )

      ! Return to calling program
      END SUBROUTINE UNZIP_A6_FIELDS

!------------------------------------------------------------------------------

      FUNCTION DO_OPEN_A6( NYMD, NHMS ) RESULT( DO_OPEN )
!
!******************************************************************************
!  Function DO_OPEN_A6 returns TRUE if is time to open the A-6 met field file
!  or FALSE otherwise.  This prevents us from opening a file which has already
!  been opened. (bmy, 6/19/03)
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
      ! DO_OPEN_A6 begins here!
      !=================================================================

      ! Initialize
      DO_OPEN = .FALSE.

      ! Return if we have already opened the file
      IF ( NYMD == LASTNYMD .and. NHMS == LASTNHMS ) THEN
         DO_OPEN = .FALSE. 
         GOTO 999
      ENDIF

#if   defined( GEOS_4 )

      ! Open file if it's 01:30 GMT or first call (GEOS-4 only)
      IF ( NHMS == 030000 .or. FIRST ) THEN
         DO_OPEN = .TRUE. 
         GOTO 999
      ENDIF

#else

      ! Open file if it's 00:00 GMT or first call (GEOS-1, GEOS-S, GEOS-3)
      IF ( NHMS == 000000 .or. FIRST ) THEN
         DO_OPEN = .TRUE. 
         GOTO 999
      ENDIF

#endif

      !=================================================================
      ! Reset quantities for next call
      !=================================================================
 999  CONTINUE
      LASTNYMD = NYMD
      LASTNHMS = NHMS
      FIRST    = .FALSE.
      
      ! Return to calling program
      END FUNCTION DO_OPEN_A6

!------------------------------------------------------------------------------

      SUBROUTINE OPEN_A6_FIELDS( NYMD, NHMS )
!
!******************************************************************************
!  Subroutine OPEN_A6_FIELDS opens the A-6 met fields file for date NYMD and 
!  time NHMS. (bmy, bdf, 6/15/98, 12/11/03)
!  
!  Arguments as input:
!  ===========================================================================
!  (1 ) NYMD (INTEGER)   : Current value of YYYYMMDD
!  (2 ) NHMS (INTEGER)   : Current value of HHMMSS
!
!  NOTES:
!  (1 ) Adapted from OPEN_MET_FIELDS of "dao_read_mod.f" (bmy, 6/19/03)
!  (2 ) Now opens either zipped or unzipped files (bmy, 12/11/03)
!******************************************************************************
!      
      ! References to F90 modules
      USE BPCH2_MOD, ONLY : GET_RES_EXT
      USE ERROR_MOD, ONLY : ERROR_STOP
      USE FILE_MOD,  ONLY : IU_A6, IOERROR
      USE TIME_MOD,  ONLY : EXPAND_DATE

#     include "CMN_SIZE"           ! Size parameters
#     include "CMN_SETUP"          ! GEOS directories

      ! Arguments
      INTEGER,          INTENT(IN) :: NYMD, NHMS

      ! Local variables
      LOGICAL, SAVE                :: FIRST = .TRUE.
      LOGICAL                      :: IT_EXISTS
      INTEGER                      :: IOS, IUNIT
      CHARACTER(LEN=255)           :: INPUT_DIR
      CHARACTER(LEN=255)           :: PATH

      !=================================================================
      ! OPEN_A6_FIELDS begins here!
      !=================================================================

      ! Open A-6 file at the proper time, or on the first call
      IF ( DO_OPEN_A6( NYMD, NHMS ) ) THEN

!-----------------------------------------------------------------------------
! Prior to 12/11/03:
! Now open either zipped or unzipped data (bmy, 12/11/03)
!#if   defined( GEOS_1 ) || defined( GEOS_STRAT )
!
!         ! Location of A-3 file in temp dir (GEOS-1, GEOS-S)
!         PATH = TRIM( TEMP_DIR )  // 'YYMMDD.a6.' // GET_RES_EXT()
!
!#else
!
!         ! Location of A-3 file in temp dir (GEOS-3, GEOS_4)
!         PATH = TRIM( TEMP_DIR )  // 'YYYYMMDD.a6.' // GET_RES_EXT()
!
!#endif
!-----------------------------------------------------------------------------

#if   defined( GEOS_1 ) 

         ! If unzipping, open GEOS-1 file in TEMP dir
         ! If not unzipping, open GEOS-1 file in DATA dir
         IF ( LUNZIP ) THEN
            PATH = TRIM( TEMP_DIR ) // 'YYMMDD.a6.' // GET_RES_EXT()
         ELSE
            PATH = TRIM( DATA_DIR ) // TRIM( GEOS_1_DIR ) // 
     &             'YYMMDD.a6.'     // GET_RES_EXT() 
         ENDIF

#elif defined( GEOS_STRAT )

         ! If unzipping, open GEOS-STRAT file in TEMP dir
         ! If not unzipping, open GEOS-STRAT file in DATA dir
         IF ( LUNZIP ) THEN
            PATH = TRIM( TEMP_DIR ) // 'YYMMDD.a6.' // GET_RES_EXT()
         ELSE
            PATH = TRIM( DATA_DIR ) // TRIM( GEOS_S_DIR ) // 
     &             'YYMMDD.a6.'     // GET_RES_EXT() 
         ENDIF

#elif defined( GEOS_3 )

         ! If unzipping, open GEOS-3 file in TEMP dir
         ! If not unzipping, open GEOS-3 file in DATA dir
         IF ( LUNZIP ) THEN
            PATH = TRIM( TEMP_DIR ) // 'YYYYMMDD.a6.' // GET_RES_EXT()
         ELSE
            PATH = TRIM( DATA_DIR ) // TRIM( GEOS_3_DIR ) // 
     &             'YYYYMMDD.a6.'   // GET_RES_EXT() 
         ENDIF

#elif defined( GEOS_4 )

         ! If unzipping, open GEOS-4 file in TEMP dir
         ! If not unzipping, open GEOS-4 file in DATA dir
         IF ( LUNZIP ) THEN
            PATH = TRIM( TEMP_DIR ) // 'YYYYMMDD.a6.' // GET_RES_EXT()
         ELSE
            PATH = TRIM( DATA_DIR ) // TRIM( GEOS_4_DIR ) // 
     &             'YYYYMMDD.a6.'   // GET_RES_EXT() 
         ENDIF

#endif

         ! Replace YYYYMMDD in PATH w/ actual date
         CALL EXPAND_DATE( PATH, NYMD, 000000 )

         ! Close previously opened A-3 file
         CLOSE( IU_A6 )

         ! Make sure the file exists before we open it!
         ! Maybe make this a function in ERROR_MOD (bmy, 6/19/03)
         INQUIRE( IU_A6, EXIST=IT_EXISTS )
            
         IF ( .not. IT_EXISTS ) THEN
            CALL ERROR_STOP( 'Could not find file!', 
     &                       'OPEN_A6_FIELDS (a6_read_mod.f)' )
         ENDIF

         ! Open the file
         OPEN( UNIT   = IU_A6,         FILE   = TRIM( PATH ),
     &         STATUS = 'OLD',         ACCESS = 'SEQUENTIAL',  
     &         FORM   = 'UNFORMATTED', IOSTAT = IOS )
               
         IF ( IOS /= 0 ) THEN
            CALL IOERROR( IOS, IU_A6, 'open_a6_fields:1' )
         ENDIF

         WRITE( 6, '( ''     - Opening: '', a )' ) TRIM( PATH )
      ENDIF

      ! Return to calling program
      END SUBROUTINE OPEN_A6_FIELDS

!------------------------------------------------------------------------------

      SUBROUTINE GET_A6_FIELDS( NYMD, NHMS )
!
!******************************************************************************
!  Subroutine GET_A6_FIELDS is a wrapper for routine READ_A6.  GET_A6_FIELDS
!  calls READ_A6 properly for reading A-6 fields from GEOS-1, GEOS-STRAT, 
!  GEOS-3, or GEOS-4 met data sets. (bmy, 6/19/03)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) NYMD (INTEGER) : YYYYMMDD
!  (2 ) NHMS (INTEGER) :  and HHMMSS of A-6 fields to be read from disk
!
!  NOTES:
!  (1 ) CFRAC has been removed from CMN_DEP.  Now use CLDFRC(I,J) from
!        "dao_mod.f" (bmy, 12/9/03)
!******************************************************************************
!
      ! References to F90 modules
      USE DAO_MOD, ONLY : CLDF,    CLDFRC, CLDMAS, CLMOSW, CLROSW, 
     &                    CLDTOPS, DTRAIN, HKBETA, HKETA,  MOISTQ, 
     &                    OPTDEP,  SPHU,   T,      UWND,   VWND,   
     &                    ZMEU,    ZMMD,   ZMMU  

#     include "CMN_SIZE"  ! Size parameters
!-----------------------------------------------------
! Prior to 12/9/03:
!#     include "CMN_DEP"   ! CFRAC
!-----------------------------------------------------

      ! Arguments
      INTEGER, INTENT(IN) :: NYMD, NHMS 

      ! Local variables
      INTEGER, SAVE       :: LASTNYMD = -1, LASTNHMS = -1

      !=================================================================
      ! GET_A6_FIELDS begins here!
      !=================================================================

      ! Skip over previously-read A-6 fields
      IF ( NYMD == LASTNYMD .and. NHMS == LASTNHMS ) THEN
         WRITE( 6, 100 ) NYMD, NHMS
 100     FORMAT( '     - A-6 met fields for NYMD, NHMS = ', 
     &           i8.8, 1x, i6.6, ' have been read already' ) 
         RETURN
      ENDIF

#if   defined( GEOS_1 )

      !=================================================================
      ! GEOS-1 and GEOS-STRAT: 
      ! read CLDF, CLDMAS, CLDTOPS, CLMO, CLRO, DTRAIN, MOISTQ
      !=================================================================
      CALL READ_A6( NYMD=NYMD,     NHMS=NHMS,       CLDF=CLDF,             
     &              CLDMAS=CLDMAS, CLDTOPS=CLDTOPS, CLMOLW=CLMOSW, 
     &              CLROLW=CLROSW, DTRAIN=DTRAIN,   MOISTQ=MOISTQ ) 

#elif defined( GEOS_STRAT )

      !=================================================================
      ! GEOS-S: read CLDF, CLDMAS, CLDTOPS, CLMO, CLRO, DTRAIN, MOISTQ
      !=================================================================
      CALL READ_A6( NYMD=NYMD,     NHMS=NHMS,       CLDF=CLDF,             
     &              CLDMAS=CLDMAS, CLDTOPS=CLDTOPS, CLMOLW=CLMOSW, 
     &              CLROLW=CLROSW, DTRAIN=DTRAIN,   MOISTQ=MOISTQ ) 

      ! Construct the 2-D CFRAC field from CLMO and CLRO
      ! since this field is missing from GEOS-STRAT data
      !-----------------------------------------------------
      ! Prior to 12/9/03:
      !CALL MAKE_CLDFRC( CLMOSW, CLROSW, CFRAC )
      !-----------------------------------------------------
      CALL MAKE_CLDFRC( CLMOSW, CLROSW, CLDFRC )

#elif defined( GEOS_3 ) 

      !=================================================================      
      ! GEOS-3: read CLDF, CLDMAS, CLDTOPS, DTRAIN, MOISTQ, OPTDEP
      !=================================================================
      CALL READ_A6( NYMD=NYMD,     NHMS=NHMS,       
     &              CLDF=CLDF,     CLDMAS=CLDMAS, CLDTOPS=CLDTOPS, 
     &              DTRAIN=DTRAIN, MOISTQ=MOISTQ, OPTDEPTH=OPTDEP )

#elif defined( GEOS_4 ) 

      !=================================================================      
      ! GEOS-4: read CLDF, HKBETA, HKETA, MOISTQ, OPTDEP, SPHU
      !              TMPU, UWND,   VWND,  ZMEU,   ZMMD,   ZMMU
      !=================================================================
      CALL READ_A6( NYMD=NYMD,     NHMS=NHMS,     
     &              CLDF=CLDF,     HKBETA=HKBETA,   HKETA=HKETA,   
     &              MOISTQ=MOISTQ, OPTDEPTH=OPTDEP, Q=SPHU,        
     &              T=T,           U=UWND,          V=VWND,        
     &              ZMEU=ZMEU,     ZMMD=ZMMD,       ZMMU=ZMMU ) 

#endif

      ! Save NYMD and NHMS for next call
      LASTNYMD = NYMD
      LASTNHMS = NHMS

      ! Return to calling program
      END SUBROUTINE GET_A6_FIELDS

!------------------------------------------------------------------------------

      SUBROUTINE MAKE_CLDFRC( CLMOSW, CLROSW, CLDFRC )
!
!******************************************************************************
!  Subroutine MAKE_CLDFRC constructs the DAO CLDFRC field from the 
!  GEOS-STRAT CLMOSW and CLROSW fields. (bmy, 3/17/99, 12/9/03) 
!
!  Arguments as Input:
!  ===========================================================================
!  (1) CLMOSW (REAL*8) : GMAO Maximum Overlap Cloud Fraction Field
!  (2) CLROSW (REAL*8) : GMAO Random  Overlap Cloud Fraction Field
!
!  Arguments as Output:
!  ===========================================================================
!  (3) CLDFRC (REAL*8) : GMAO Column Cloud Fraction
!
!  NOTES:
!  (1 ) CLDFRC is not archived for GEOS-STRAT data, so we must compute it
!        from the CLMO and CLRO cloud fraction fields. (bmy, 6/26/00)
!  (2 ) CLDFRC is dimensioned (MAXIJ = IIPAR*JJPAR) for compatibility with
!        the Harvard dry deposition subroutines (bmy, 6/26/00)
!  (3 ) Save CLDFRC to ND67 diagnostic (bmy, 6/28/00)
!  (4 ) Updated comments (bmy, 4/4/01)
!  (5 ) Replaced all instances of IM with IIPAR and JM with JJPAR, in order
!        to prevent namespace confusion for the new TPCORE (bmy, 6/25/02)
!  (6 ) Moved from "dao_mod.f" to "a6_read_mod.f" (bmy, 6/19/03)
!  (7 ) CLDFRC is now a 2-D array (bmy, 12/9/03)
!******************************************************************************
!
      ! Reference to F90 modules
      USE DIAG_MOD, ONLY : AD67

#     include "CMN_SIZE"
#     include "CMN_DIAG"

      ! Arguments
      REAL*8, INTENT(IN)  :: CLMOSW(LLPAR,IIPAR,JJPAR)
      REAL*8, INTENT(IN)  :: CLROSW(LLPAR,IIPAR,JJPAR)
      !-------------------------------------------
      ! Prior to 12/9/03:
      !REAL*8, INTENT(OUT) :: CLDFRC(MAXIJ)
      !-------------------------------------------
      REAL*8, INTENT(OUT) :: CLDFRC(IIPAR,JJPAR)

      ! Local variables
      !-------------------------------------------
      ! Prior to 12/9/03
      !INTEGER             :: I, IJLOOP, J, L
      !-------------------------------------------
      INTEGER             :: I, J, L
      REAL*8              :: C1, C2

      !=================================================================
      ! MAKE_CLDFRC begins here!!!
      !
      ! Compute CLDFRC value for each location (I,J)
      !
      ! NOTES:
      ! (1) CLDFRC = the fractional cloud cover as seen from the 
      !     surface looking up, that is, along a vertical line of sight
      !     extending from the surface through the top of the atmosphere. 
      !
      !     The maximum overlap clear sky probability computed along 
      !     a line of sight from a grid box (I,J,L=1) at the surface to
      !     a grid box (I,J,L=LLPAR) at the top of the atmosphere is:
      ! 
      !         C1 = 1 - (maximum of CLRO(L,I,J)), L = 1, LLPAR
      !
      !     The random overlap clear sky probability computed along 
      !     a line of sight from a grid box (I,J,L=1) at the surface to 
      !     a grid box (I,J,L=LLPAR) at the top of the atmosphere is:
      !         
      !         C2 = product of (1 - CLRO(L,I,J)), L = 1, LLPAR
      !
      !     Thus the fractional cloud cover as seen from grid box (I,J,L=1)
      !     at the surface looking up is:
      !
      !         CLDFRC( @ L, I, J ) = 1.0 - (C1 * C2)
      !
      !     CLDFRC is used by the Harvard CTM dry deposition routines.
      !
      ! (2) In GEOS-1 and GEOS-STRAT, CLMOLW=CLMOSW and CLROLW=CLROSW.
      !=================================================================
      !--------------------------
      ! Prior to 12/9/03:
      !IJLOOP = 0
      !--------------------------
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         !-----------------------
         ! Prior to 12/9/03:
         !IJLOOP = IJLOOP + 1
         !-----------------------

         C1 = CLMOSW(1,I,J)       
         C2 = 1.0d0 - CLROSW(1,I,J) 

         DO L = 2, LLPAR
            IF ( CLMOSW(L,I,J) > CLMOSW(L-1,I,J) ) THEN
               C1 = CLMOSW(L,I,J)
            ENDIF

            C2 = C2 * ( 1.0d0 - CLROSW(L,I,J) )
         ENDDO

         C1             = 1.0d0 - C1
         !--------------------------------------------------------
         ! Prior to 12/9/03:
         !CLDFRC(IJLOOP) = 1.0d0 - (C1 * C2)
         !--------------------------------------------------------
         CLDFRC(I,J) = 1.0d0 - (C1 * C2)
      ENDDO
      ENDDO
      
      !=================================================================
      ! ND67 diagnostic -- save CLDFRC as tracer #10
      !=================================================================
      IF ( ND67 > 0 ) THEN
         !-------------------------------
         ! Prior to 12/9/03:
         !IJLOOP = 0
         !-------------------------------
         DO J = 1, JJPAR
         DO I = 1, IIPAR
            !---------------------------------------------------
            ! Prior to 12/9/03:
            !IJLOOP       = IJLOOP + 1
            !AD67(I,J,10) = AD67(I,J,10) + CLDFRC(IJLOOP)
            !---------------------------------------------------
            AD67(I,J,10) = AD67(I,J,10) + CLDFRC(I,J)
         ENDDO
         ENDDO
      ENDIF

      ! Return to calling program
      END SUBROUTINE MAKE_CLDFRC

!------------------------------------------------------------------------------

      FUNCTION GET_N_A6() RESULT( N_A6 )
!
!******************************************************************************
!  Function GET_N_A6 returns the number of A-6 fields per met data set
!  (GEOS-1, GEOS-STRAT, GEOS-3, GEOS-4). (bmy, 6/19/03) 
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) NYMD (INTEGER) : YYYYMMDD for which to read in A-6 fields
!
!  NOTES:
!******************************************************************************
!
#     include "CMN_SIZE" 

      ! Function value
      INTEGER :: N_A6

      !=================================================================
      ! GET_N_A6 begins here!
      !=================================================================
#if   defined( GEOS_1 ) || defined( GEOS_STRAT )

      ! GEOS-1 and GEOS-STRAT have 5 A-6 fields
      N_A6 = 5

#elif defined( GEOS_3 )

      ! GEOS-3 has 6 A-6 fields
      N_A6 = 6

#elif defined( GEOS_4 )

      ! GEOS-4 has 12 A-6 fields
      N_A6 = 12

#endif

      ! Return to calling program
      END FUNCTION GET_N_A6

!---------------------------------------------------------------------------

      FUNCTION CHECK_TIME( XYMD, XHMS, NYMD, NHMS ) RESULT( ITS_TIME )
!
!******************************************************************************
!  Function CHECK_TIME checks to see if the timestamp of the A-3 field just
!  read from disk matches the current time.  If so, then it's time to return
!  the A-3 field to the calling program. (bmy, 6/19/03)
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

#if   defined( GEOS_1 ) || defined( GEOS_STRAT )

      ! Arguments
      REAL*4,  INTENT(IN) :: XYMD, XHMS 
      INTEGER, INTENT(IN) :: NYMD, NHMS

      ! Function value
      LOGICAL             :: ITS_TIME

      !=================================================================
      ! GEOS-1 and GEOS-STRAT: XYMD and XHMS are REAL*4
      !=================================================================
      IF ( INT(XYMD) == NYMD-19000000 .AND. INT(XHMS) == NHMS ) THEN
         ITS_TIME = .TRUE.
      ELSE
         ITS_TIME = .FALSE.
      ENDIF

#else

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

#endif

      ! Return to calling program
      END FUNCTION CHECK_TIME

!-----------------------------------------------------------------------------

      SUBROUTINE READ_A6( NYMD,   NHMS,   
     &                    CLDF,   CLDMAS,   CLDTOPS, CLMOLW, 
     &                    CLROLW, DTRAIN,   HKBETA,  HKETA,
     &                    MOISTQ, OPTDEPTH, Q,       T,       
     &                    U,      V,        ZMEU,    ZMMD,
     &                    ZMMU )
!
!*****************************************************************************
!  Subroutine READ_A6 reads A-6 (avg 6-hr) met fields from disk. 
!  (bmy, 6/5/98, 10/28/03)
! 
!  Arguments as input:
!  ===========================================================================
!  (1 ) NYMD    : YYYYMMDD
!  (2 ) NHMS    :  and HHMMSS of A-6 met fields to be accessed
!
!  A-6 Met Fields as Output (Optional Arguments):
!  ============================================================================
!  (3 ) CLDF    : (3-D) Total cloud fractions                [unitless]
!  (4 ) CLDMAS  : (3-D) Cloud mass flux field                [kg/m2/600s]
!  (5 ) CLDTOPS : (2-D) CTM Level in which cloud top occurs  [unitless]
!  (6 ) CLMOLW  : (3-D) GEOS-1 LW max-overlap cloud frac     [unitless]
!  (7 ) CLROLW  : (3-D) GEOS-1 LW random-overlap cloud frac  [unitless]
!  (8 ) DTRAIN  : (3-D) Detrainment field                    [kg/m2/600s]
!  (9 ) HKBETA  : (3-D) Hack overshoot parameter             [unitless]
!  (10) HKETA   : (3-D) Hack convective mass flux            [kg/m2/s]
!  (11) MOISTQ  : (3-D) DAO water vapor tendency d           [g H2O/kg air/day]
!  (12) OPTDEPTH: (3-D) GEOS-2 grid box optical depth        [unitless]
!  (13) Q       : (3-D) Specific humidity                    [g H2O/kg air]
!  (14) T       : (3-D) Temperature                          [K]
!  (15) U       : (3-D) Zonal winds                          [m/s]
!  (16) V       : (3-D) Meridional winds                     [m/s]
!  (17) ZMEU    : (3-D) Zhang/McFarlane updraft entrainment  [Pa/s]
!  (18) ZMMD    : (3-D) Zhang/McFarlane downdraft mass flux  [Pa/s]
!  (19) ZMMU    : (3-D) Zhang/McFarlane updraft mass flux    [Pa/s]
!
!  NOTES:
!  (1 ) Adapted from READ_A6 of "dao_read_mod.f" (bmy, 6/19/03)
!  (2 ) Now use function TIMESTAMP_STRING from "time_mod.f" for formatted 
!        date/time output. (bmy, 10/28/03)
!******************************************************************************
!
      ! References to F90 modules
      USE DIAG_MOD,     ONLY : AD66,        AD67
      USE FILE_MOD,     ONLY : IOERROR,     IU_A6
      USE TIME_MOD,     ONLY : SET_CT_A6,   TIMESTAMP_STRING
      USE TRANSFER_MOD, ONLY : TRANSFER_A6, TRANSFER_3D

#     include "CMN_SIZE"             ! Size parameters
#     include "CMN_DIAG"             ! ND66, ND67

      ! Arguments
      INTEGER, INTENT(IN)            :: NYMD, NHMS
      INTEGER, INTENT(OUT), OPTIONAL :: CLDTOPS(IIPAR,JJPAR)
      REAL*8,  INTENT(OUT), OPTIONAL :: CLDF(LLPAR,IIPAR,JJPAR)
      REAL*8,  INTENT(OUT), OPTIONAL :: CLDMAS(IIPAR,JJPAR,LLPAR)
      REAL*8,  INTENT(OUT), OPTIONAL :: CLMOLW(LLPAR,IIPAR,JJPAR)
      REAL*8,  INTENT(OUT), OPTIONAL :: CLROLW(LLPAR,IIPAR,JJPAR)
      REAL*8,  INTENT(OUT), OPTIONAL :: DTRAIN(IIPAR,JJPAR,LLPAR)
      REAL*8,  INTENT(OUT), OPTIONAL :: HKBETA(IIPAR,JJPAR,LLPAR)
      REAL*8,  INTENT(OUT), OPTIONAL :: HKETA(IIPAR,JJPAR,LLPAR)
      REAL*8,  INTENT(OUT), OPTIONAL :: MOISTQ(LLPAR,IIPAR,JJPAR) 
      REAL*8,  INTENT(OUT), OPTIONAL :: OPTDEPTH(LLPAR,IIPAR,JJPAR)
      REAL*8,  INTENT(OUT), OPTIONAL :: Q(IIPAR,JJPAR,LLPAR)
      REAL*8,  INTENT(OUT), OPTIONAL :: T(IIPAR,JJPAR,LLPAR)
      REAL*8,  INTENT(OUT), OPTIONAL :: U(IIPAR,JJPAR,LLPAR)
      REAL*8,  INTENT(OUT), OPTIONAL :: V(IIPAR,JJPAR,LLPAR)
      REAL*8,  INTENT(OUT), OPTIONAL :: ZMEU(IIPAR,JJPAR,LLPAR)
      REAL*8,  INTENT(OUT), OPTIONAL :: ZMMD(IIPAR,JJPAR,LLPAR)
      REAL*8,  INTENT(OUT), OPTIONAL :: ZMMU(IIPAR,JJPAR,LLPAR)

      ! Local variables
      INTEGER                        :: I, IJLOOP, J, K, L
      INTEGER                        :: IOS, NFOUND, N_A6
      REAL*4                         :: Q3(IGLOB,JGLOB,LGLOB)       
      REAL*8                         :: C1, C2
      REAL*8                         :: TAUCLD(LLPAR,IIPAR,JJPAR)
      REAL*8                         :: CLDTOT(LLPAR,IIPAR,JJPAR)
      CHARACTER(LEN=8)               :: NAME
      CHARACTER(LEN=16)              :: STAMP

      ! XYMD, XHMS must be REAL*4 for GEOS-1 and GEOS-STRAT
      ! but INTEGER for GEOS-3 and GEOS-4 (bmy, 6/19/03)
#if   defined( GEOS_1 ) || defined( GEOS_STRAT )
      REAL*4                         :: XYMD, XHMS
#else
      INTEGER                        :: XYMD, XHMS
#endif

      !=================================================================
      ! READ_A6 begins here!      
      !=================================================================

      ! Get number of A-6 fields
      N_A6   = GET_N_A6()

      ! Zero number of fields that we have found
      NFOUND = 0

      !=================================================================
      ! Read the A-6 fields from disk
      !=================================================================
      DO

         ! A-6 field name
         READ( IU_A6, IOSTAT=IOS ) NAME

         ! IOS < 0: End-of-file; make sure we've found 
         ! all the A-6 fields before exiting this loop
         IF ( IOS < 0 ) THEN
            CALL A6_CHECK( NFOUND, N_A6 )
            EXIT
         ENDIF

         ! IOS > 0: True I/O Error, stop w/ error msg 
         IF ( IOS > 0 ) CALL IOERROR( IOS, IU_A6, 'read_a6:1' )

         ! CASE statement for A-6 fields
         SELECT CASE ( TRIM( NAME ) )

            !--------------------------------
            ! CLDMAS: cloud mass flux
            ! (GEOS-1, GEOS-STRAT, GEOS-3)
            !--------------------------------
            CASE ( 'CLDMAS' ) 
               READ( IU_A6, IOSTAT=IOS ) XYMD, XHMS, Q3
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A6, 'read_a6:2' )
             
               IF ( CHECK_TIME( XYMD, XHMS, NYMD, NHMS ) ) THEN
                  IF ( PRESENT( CLDMAS ) ) CALL TRANSFER_3D( Q3,CLDMAS )
                  NFOUND = NFOUND + 1
               ENDIF

            !--------------------------------
            ! CLDTOT: 3-D total cloud frac
            ! (GEOS-3, GEOS-4 only)
            !--------------------------------
            CASE ( 'CLDTOT' ) 
               READ( IU_A6, IOSTAT=IOS ) XYMD, XHMS, Q3
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A6, 'read_a6:3' )

               IF ( CHECK_TIME( XYMD, XHMS, NYMD, NHMS ) ) THEN
                  CALL TRANSFER_A6( Q3, CLDTOT )
                  NFOUND = NFOUND +1 
               ENDIF

            !--------------------------------
            ! CLMOLW: max overlap cloud frac
            ! (GEOS-1 and GEOS-STRAT only)
            !--------------------------------
            CASE ( 'CLMOLW' ) 
               READ( IU_A6, IOSTAT=IOS ) XYMD, XHMS, Q3
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A6, 'read_a6:4' )

               IF ( CHECK_TIME( XYMD, XHMS, NYMD, NHMS ) ) THEN
                  IF ( PRESENT( CLMOLW ) ) CALL TRANSFER_A6( Q3,CLMOLW )
                  NFOUND = NFOUND + 1 
               ENDIF

            !--------------------------------
            ! CLROLW: max overlap cloud frac
            ! (GEOS-1 and GEOS-STRAT only)
            !--------------------------------
            CASE ( 'CLROLW' ) 
               READ( IU_A6, IOSTAT=IOS ) XYMD, XHMS, Q3
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A6, 'read_a6:5' )

               IF ( CHECK_TIME( XYMD, XHMS, NYMD, NHMS ) ) THEN
                  IF ( PRESENT( CLROLW ) ) CALL TRANSFER_A6( Q3,CLROLW )
                  NFOUND = NFOUND + 1 
               ENDIF
               
            !--------------------------------
            ! DTRAIN: cloud detrainment
            ! (GEOS-1, GEOS-STRAT, GEOS-3)
            !--------------------------------
            CASE ( 'DTRAIN' ) 
               READ( IU_A6, IOSTAT=IOS ) XYMD, XHMS, Q3
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A6, 'read_a6:6' )
 
               IF ( CHECK_TIME( XYMD, XHMS, NYMD, NHMS ) ) THEN
                  IF ( PRESENT( DTRAIN ) ) CALL TRANSFER_3D( Q3,DTRAIN )
                  NFOUND = NFOUND + 1
               ENDIF

            !--------------------------------
            ! HKBETA: Hack overshoot param. 
            ! (GEOS-4 only)
            !--------------------------------
            CASE ( 'HKBETA' ) 
               READ( IU_A6, IOSTAT=IOS ) XYMD, XHMS, Q3
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A6, 'read_a6:7' )
 
               IF ( CHECK_TIME( XYMD, XHMS, NYMD, NHMS ) ) THEN
                  IF ( PRESENT( HKBETA ) ) CALL TRANSFER_3D( Q3,HKBETA )
                  NFOUND = NFOUND + 1
               ENDIF

            !--------------------------------
            ! HKETA: Hack conv mass flux 
            ! (GEOS-4 only)
            !--------------------------------
            CASE ( 'HKETA' ) 
               READ( IU_A6, IOSTAT=IOS ) XYMD, XHMS, Q3
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A6, 'read_a6:8' )
 
               IF ( CHECK_TIME( XYMD, XHMS, NYMD, NHMS ) ) THEN
                  IF ( PRESENT( HKETA ) ) CALL TRANSFER_3D( Q3, HKETA )
                  NFOUND = NFOUND + 1
               ENDIF

            !--------------------------------
            ! MOISTQ: water vapor tendency 
            ! (all GEOS versions)
            !--------------------------------
            CASE ( 'MOISTQ' ) 
               READ( IU_A6, IOSTAT=IOS ) XYMD, XHMS, Q3
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A6, 'read_a6:7' )

               IF ( CHECK_TIME( XYMD, XHMS, NYMD, NHMS ) ) THEN
                  IF ( PRESENT( MOISTQ ) ) CALL TRANSFER_A6( Q3,MOISTQ )
                  NFOUND = NFOUND + 1 
               ENDIF

            !--------------------------------
            ! OPTDEPTH: grid box optical depth
            ! (GEOS-3 and GEOS-4 only)
            !--------------------------------
            CASE ( 'OPTDEPTH' ) 
               READ( IU_A6, IOSTAT=IOS ) XYMD, XHMS, Q3
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A6, 'read_a6:8' )

               IF ( CHECK_TIME( XYMD, XHMS, NYMD, NHMS ) ) THEN
                  IF ( PRESENT( OPTDEPTH ) ) THEN
                     CALL TRANSFER_A6( Q3, OPTDEPTH )
                  ENDIF
                  NFOUND = NFOUND + 1 
               ENDIF

            !--------------------------------
            ! Q: 6-h avg specific humidity
            ! (GEOS-4 only)
            !--------------------------------
            CASE ( 'Q' ) 
               READ( IU_A6, IOSTAT=IOS ) XYMD, XHMS, Q3
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A6, 'read_a6:9' )

               IF ( CHECK_TIME( XYMD, XHMS, NYMD, NHMS ) ) THEN
                  IF ( PRESENT( Q ) ) CALL TRANSFER_3D( Q3, Q )
                  NFOUND = NFOUND + 1 
               ENDIF

            !--------------------------------
            ! T: 6-h avg temperature
            ! (GEOS-4 only)
            !--------------------------------
            CASE ( 'T' ) 
               READ( IU_A6, IOSTAT=IOS ) XYMD, XHMS, Q3
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A6, 'read_a6:10' )

               IF ( CHECK_TIME( XYMD, XHMS, NYMD, NHMS ) ) THEN
                  IF ( PRESENT( T ) ) CALL TRANSFER_3D( Q3, T )
                  NFOUND = NFOUND + 1 
               ENDIF

            !--------------------------------
            ! U: 6-h avg zonal wind
            ! (GEOS-4 only)
            !--------------------------------
            CASE ( 'U' ) 
               READ( IU_A6, IOSTAT=IOS ) XYMD, XHMS, Q3
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A6, 'read_a6:11' )

               IF ( CHECK_TIME( XYMD, XHMS, NYMD, NHMS ) ) THEN
                  IF ( PRESENT( U ) ) CALL TRANSFER_3D( Q3, U )
                  NFOUND = NFOUND + 1 
               ENDIF

            !--------------------------------
            ! V: 6-h avg meridional wind
            ! (GEOS-4 only)
            !--------------------------------
            CASE ( 'V' ) 
               READ( IU_A6, IOSTAT=IOS ) XYMD, XHMS, Q3
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A6, 'read_a6:12' )

               IF ( CHECK_TIME( XYMD, XHMS, NYMD, NHMS ) ) THEN
                  IF ( PRESENT( V ) ) CALL TRANSFER_3D( Q3, V )
                  NFOUND = NFOUND + 1 
               ENDIF

            !--------------------------------
            ! ZMEU: Z&M updraft entrainment
            ! (GEOS-4 only)
            !--------------------------------
            CASE ( 'ZMEU' ) 
               READ( IU_A6, IOSTAT=IOS ) XYMD, XHMS, Q3
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A6, 'read_a6:13' )

               IF ( CHECK_TIME( XYMD, XHMS, NYMD, NHMS ) ) THEN
                  IF ( PRESENT( ZMEU ) ) CALL TRANSFER_3D( Q3, ZMEU )
                  NFOUND = NFOUND + 1 
               ENDIF

            !--------------------------------
            ! ZMMD: Z&M downdraft mass flux
            ! (GEOS-4 only)
            !--------------------------------
            CASE ( 'ZMMD' ) 
               READ( IU_A6, IOSTAT=IOS ) XYMD, XHMS, Q3
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A6, 'read_a6:14' )

               IF ( CHECK_TIME( XYMD, XHMS, NYMD, NHMS ) ) THEN
                  IF ( PRESENT( ZMMD ) ) CALL TRANSFER_3D( Q3, ZMMD )
                  NFOUND = NFOUND + 1 
               ENDIF

            !--------------------------------
            ! ZMMU: Z&M updraft mass flux
            ! (GEOS-4 only)
            !--------------------------------
            CASE ( 'ZMMU' ) 
               READ( IU_A6, IOSTAT=IOS ) XYMD, XHMS, Q3
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A6, 'read_a6:15' )

               IF ( CHECK_TIME( XYMD, XHMS, NYMD, NHMS ) ) THEN
                  IF ( PRESENT( ZMMU ) ) CALL TRANSFER_3D( Q3, ZMMU )
                  NFOUND = NFOUND + 1 
               ENDIF

            !--------------------------------
            ! TAUCLD: in-cloud optical depth 
            ! Just skip over this
            !--------------------------------
            CASE ( 'TAUCLD' ) 
               READ( IU_A6, IOSTAT=IOS ) XYMD, XHMS, Q3
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A6, 'read_a6:16' )

               IF ( CHECK_TIME( XYMD, XHMS, NYMD, NHMS ) ) THEN
                  NFOUND = NFOUND + 1 
               ENDIF

            !--------------------------------
            ! KH: Just skip over this
            !--------------------------------
            CASE ( 'KH' ) 
               READ( IU_A6, IOSTAT=IOS ) XYMD, XHMS, Q3
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A6, 'read_a6:17' )

               IF ( CHECK_TIME( XYMD, XHMS, NYMD, NHMS ) ) THEN
                  NFOUND = NFOUND + 1 
               ENDIF

            ! Field not found
            CASE DEFAULT
               WRITE ( 6, '(a)' ) 'Searching for next A-6 field!'
               
         END SELECT

         !==============================================================
         ! If we have found all the fields for this time, then exit 
         ! the loop.  Otherwise, go on to the next iteration.
         !==============================================================
         IF ( CHECK_TIME( XYMD, XHMS, NYMD, NHMS ) .AND. 
     &        NFOUND == N_A6 ) THEN
!------------------------------------------------------------------------
! Prior to 10/28/03:
! Now use TIMESTAMP_STRING for formatted output (bmy, 10/28/03)
!            WRITE( 6, 210 ) NFOUND, NYMD, NHMS
! 210        FORMAT( '     - Found all ', i3, 
!     &              ' A-6 met fields for NYMD, NHMS = ', 
!     &              i8.8, 1x, i6.6 )
!------------------------------------------------------------------------
            STAMP = TIMESTAMP_STRING( NYMD, NHMS )
            WRITE( 6, 210 ) NFOUND, STAMP
 210        FORMAT( '     - Found all ', i3, ' A-6 met fields for ', a )
            EXIT
         ENDIF
      ENDDO

      !=================================================================
      ! Due to an error in the DAO archiving process, the CLDMAS and 
      ! DTRAIN fields have units of [kg/m2/600s].  Divide here by 600 
      ! to convert CLDMAS and DTRAIN into units of [kg/m2/s].
      !=================================================================
      IF ( PRESENT( CLDMAS ) ) CLDMAS = CLDMAS / 600d0
      IF ( PRESENT( DTRAIN ) ) DTRAIN = DTRAIN / 600d0

      !=================================================================
      ! CLDTOPS(I,J) = level of convective cloud top at (I,J).
      ! GEOS-CHEM cloud top at (I,J) is at top of first level where 
      ! CLDMAS goes from being nonzero to zero.  In other words:
      !
      ! If at vertical level L, CLDMAS(L,  :,:) == 0 and 
      !                         CLDMAS(L-1,:,:) /= 0 
      !
      ! then CLDTOPS(I,J) = L.
      !=================================================================
      IF ( PRESENT( CLDTOPS ) .and. PRESENT( CLDMAS ) ) THEN
         DO J = 1, JJPAR
            DO I = 1, IIPAR
               K = 1
               DO L = 1, LLPAR
                  IF ( CLDMAS(I,J,L) > 0d0 ) THEN
                     K = K + 1
                  ENDIF
               ENDDO
               CLDTOPS(I,J) = K
            ENDDO
         ENDDO
      ENDIF

      !=================================================================
      ! CLDF(IIPAR,JJPAR,LLPAR), is the DAO total cloud fraction.
      ! 
      ! For GEOS-1 and GEOS-STRAT, CLDF is computed from the CLMO and
      ! CLRO cloud fraction fields as follows:
      !     
      ! The total clear sky probability at grid box (I,J,L) is:  
      !                                                                       
      !       ( 1 - CLMO(I,J,L) ) * ( 1 - CLRO(I,J,L) ),            
      !                                                             
      ! thus the total cloudy-sky probability (cloud fraction) is:
      !                                                                       
      !       1 - ( 1 - CLMO(I,J,L) ) * ( 1 - CLRO(I,J,L) ).   
      !
      ! For GEOS-3 and GEOS-4, CLDF is read directly from disk 
      ! as the CLDTOT met field from the loop above.
      !=================================================================
#if   defined( GEOS_1 ) || defined( GEOS_STRAT )
      IF ( PRESENT( CLDF ) ) THEN
         CLDF = 1d0 - ( ( 1d0 - CLMOLW ) * ( 1d0 - CLROLW ) )
      ENDIF

#else
      IF ( PRESENT( CLDF ) ) THEN
         CLDF = CLDTOT
      ENDIF

#endif

      !=================================================================
      ! For 1998 GEOS-3 fields only, create OPTDEPTH = TAUCLD * CLDTOT
      !
      ! The 1998 fields only store TAUCLD, which is the in-cloud 
      ! optical depth.  The actual grid box optical depth is 
      ! TAUCLD * CLDTOT, which is what FAST-J needs. (bmy, 10/11/01)
      !=================================================================
#if   defined( GEOS_3 )
      IF ( PRESENT( OPTDEPTH ) .and. ( NYMD / 10000 ) == 1998 ) THEN
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
         DO J = 1, JJPAR
         DO I = 1, IIPAR
         DO L = 1, LLPAR
            OPTDEPTH(L,I,J) = TAUCLD(L,I,J) * CLDTOT(L,I,J) 
         ENDDO
         ENDDO
         ENDDO
!$OMP END PARALLEL DO
      ENDIF
#endif

      !=================================================================
      ! MOISTQ < 0 denotes precipitation.  Convert negative values to
      ! positives, and then divide by 8.64d7 to convert to units of
      ! [kg H2O/kg air/s].  (bmy, 4/5/99)
      !=================================================================
      IF ( PRESENT( MOISTQ ) ) MOISTQ = -MOISTQ / 8.64d7

      !=================================================================
      ! ND66 diagnostic: A-6 fields
      !
      ! (1 ) UWND   : 6-h average U-winds             [m/s]
      ! (2 ) VWND   : 6=h average V-winds             [m/s]
      ! (3 ) TMPU   : 6-h average Temperature         [K]
      ! (4 ) SPHU   : 6-h average Specific humidity   [g H20/kg air]   
      ! (5 ) CLDMAS : Convective Mass Flux            [kg/m2/s] 
      !=================================================================
      IF ( ND66 > 0 ) THEN
         IF ( PRESENT( U ) ) THEN 
            AD66(:,:,1:LD66,1) = AD66(:,:,1:LD66,1) + U(:,:,1:LD66)
         ENDIF  

         IF ( PRESENT( V ) ) THEN 
            AD66(:,:,1:LD66,2) = AD66(:,:,1:LD66,2) + V(:,:,1:LD66)
         ENDIF  

         IF ( PRESENT( T ) ) THEN 
            AD66(:,:,1:LD66,3) = AD66(:,:,1:LD66,3) + T(:,:,1:LD66)
         ENDIF  

         IF ( PRESENT( Q ) ) THEN 
            AD66(:,:,1:LD66,4) = AD66(:,:,1:LD66,4) + Q(:,:,1:LD66)
         ENDIF  
         
         IF ( PRESENT( CLDMAS ) ) THEN 
            AD66(:,:,1:LD66,5) = AD66(:,:,1:LD66,5) + CLDMAS(:,:,1:LD66)
         ENDIF  
      ENDIF

      !=================================================================
      ! ND67 diagnostic: Accumulating DAO surface fields
      ! Field # 16 is the cloud top heights
      !=================================================================
      IF ( ND67 > 0 ) THEN 
         IF ( PRESENT( CLDTOPS ) ) AD67(:,:,16) = AD67(:,:,16) + CLDTOPS
      ENDIF  

      !=================================================================
      ! Update A-6 fields diagnostic counter
      !=================================================================
      CALL SET_CT_A6( INCREMENT=.TRUE. )

      ! Return to calling program
      END SUBROUTINE READ_A6

!------------------------------------------------------------------------------

      SUBROUTINE A6_CHECK( NFOUND, N_A6 )
!
!******************************************************************************
!  Subroutine A6_CHECK prints an error message if not all of the A-6 met 
!  fields are found.  The run is also terminated. (bmy, 10/27/00, 6/19/03)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) NFOUND (INTEGER) : # of A-6 met fields read from disk
!  (2 ) N_A6   (INTEGER) : # of A-6 met fields expected to be read from disk
!
!  NOTES
!  (1 ) Adapted from DAO_CHECK from "dao_read_mod.f" (bmy, 6/19/03)
!******************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD, ONLY : GEOS_CHEM_STOP

      ! Arguments
      INTEGER, INTENT(IN) :: NFOUND, N_A6

      !=================================================================
      ! A6_CHECK begins here!
      !=================================================================
      IF ( NFOUND /= N_A6 ) THEN
         WRITE( 6, '(a)' ) REPEAT( '=', 79 )
         WRITE( 6, '(a)' ) 'ERROR -- not enough A-6 fields found!'      

         WRITE( 6, 120   ) N_A6, NFOUND
 120     FORMAT( 'There are ', i2, ' fields but only ', i2 ,
     &           ' were found!' )

         WRITE( 6, '(a)' ) '### STOP in A6_CHECK (dao_read_mod.f)'
         WRITE( 6, '(a)' ) REPEAT( '=', 79 )

         ! Deallocate arrays and stop (bmy, 10/15/02)
         CALL GEOS_CHEM_STOP
      ENDIF

      ! Return to calling program
      END SUBROUTINE A6_CHECK

!------------------------------------------------------------------------------

      END MODULE A6_READ_MOD
