! $Id: a3_read_mod.f,v 1.2 2003/10/30 16:17:16 bmy Exp $
      MODULE A3_READ_MOD
!
!******************************************************************************
!  Module A3_READ_MOD contains routines that unzip, open, and read the
!  GEOS-CHEM A-3 (avg 3-hour) met fields from disk. (bmy, 6/23/03)
! 
!  Module Routines:
!  =========================================================================
!  (1 ) UNZIP_A3_FIELDS : Unzips & copies met field files to a temp dir
!  (2 ) DO_OPEN_A3      : Returns TRUE if it's time to read A-3 fields
!  (3 ) OPEN_A3_FIELDS  : Opens met field files residing in the temp dir
!  (4 ) GET_A3_FIELDS   : Wrapper for routine READ_A3
!  (5 ) GET_N_A3        : Returns # of A-3 fields for each DAO data set 
!  (6 ) CHECK_TIME      : Tests if A-3 met field timestamps equal current time
!  (7 ) READ_A3         : Reads A-3 fields from disk
!  (8 ) ARCHIVE_ND67_1D : Archives 1-D data arrays into the ND67 diagnostic
!  (9 ) A3_CHECK        : Checks if we have found all of the A-3 fields
! 
!  GEOS-CHEM modules referenced by a3_read_mod.f
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
!  (1 ) Adapted from "dao_read_mod.f" (bmy, 6/23/03)
!******************************************************************************
!
      IMPLICIT NONE

      !=================================================================
      ! MODULE PRIVATE DECLARATIONS -- keep certain internal variables 
      ! and routines from being seen outside "a3_read_mod.f"
      !=================================================================

      ! PRIVATE module routines
      PRIVATE :: A3_CHECK, CHECK_TIME, DO_OPEN_A3, GET_N_A3, READ_A3

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement 
      !=================================================================
      CONTAINS

!------------------------------------------------------------------------------

      SUBROUTINE UNZIP_A3_FIELDS( OPTION, NYMD )
!
!*****************************************************************************
!  Subroutine UNZIP_A3_FIELDS invokes a FORTRAN system call to uncompress
!  GEOS-CHEM A-3 met field files and store the uncompressed data in a 
!  temporary directory, where GEOS-CHEM can read them.  The original data 
!  files are not disturbed.  (bmy, bdf, 6/15/98, 6/23/03)
!
!  Arguments as input:
!  ===========================================================================
!  (1 ) OPTION (CHAR*(*)) : Option
!  (2 ) NYMD   (INTEGER ) : YYYYMMDD of A-3 file to be unzipped (optional)
!
!  NOTES:
!  (1 ) Adapted from UNZIP_MET_FIELDS of "dao_read_mod.f" (bmy, 6/23/03)
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
      CHARACTER(LEN=255)            :: A3_FILE_GZ, A3_FILE
      CHARACTER(LEN=255)            :: UNZIP_BG,   UNZIP_FG
      CHARACTER(LEN=255)            :: REMOVE_ALL, REMOVE_DATE

      !=================================================================
      ! UNZIP_A3_FIELDS begins here!
      !=================================================================
      IF ( PRESENT( NYMD ) ) THEN
      
#if   defined( GEOS_1 )

         ! Location of zipped A-3 file in data dir (GEOS-1)
         A3_FILE_GZ  = TRIM( DATA_DIR )      // TRIM( GEOS_1_DIR ) // 
     &                 'YYMM/YYMMDD.a3.'     // GET_RES_EXT()      // 
     &                 TRIM( ZIP_SUFFIX )

         ! Location of unzipped A-3 file in temp dir (GEOS-1)
         A3_FILE     = TRIM( TEMP_DIR )      // 'YYMMDD.a3.'       // 
     &                 GET_RES_EXT()
         
         ! Remove A-3 files for this date from temp dir (GEOS-1)
         REMOVE_DATE = TRIM( REMOVE_CMD )    // ' '                // 
     &                 TRIM( TEMP_DIR   )    // 'YYMMDD.a3.'       // 
     &                 GET_RES_EXT()  

#elif defined( GEOS_STRAT )

         ! Location of zipped A-3 file in data dir (GEOS-STRAT)
         A3_FILE_GZ  = TRIM( DATA_DIR )      // TRIM( GEOS_S_DIR ) // 
     &                 'YYMM/YYMMDD.a3.'     // GET_RES_EXT()      // 
     &                 TRIM( ZIP_SUFFIX )

         ! Location of unzipped A-3 file in temp dir (GEOS-STRAT)
         A3_FILE     = TRIM( TEMP_DIR )      // 'YYMMDD.a3.'       // 
     &                 GET_RES_EXT()

         ! Remove A-3 files for this date from temp dir (GEOS-STRAT)
         REMOVE_DATE = TRIM( REMOVE_CMD )    // ' '                // 
     &                 TRIM( TEMP_DIR   )    // 'YYMMDD.a3.'       // 
     &                 GET_RES_EXT()  

#elif defined( GEOS_3 )

         ! Location of zipped A-3 file in data dir (GEOS-3)
         A3_FILE_GZ  = TRIM( DATA_DIR )      // TRIM( GEOS_3_DIR ) // 
     &                 'YYYYMM/YYYYMMDD.a3.' // GET_RES_EXT()      // 
     &                 TRIM( ZIP_SUFFIX )

         ! Location of unzipped A-3 file in temp dir (GEOS-3)
         A3_FILE     = TRIM( TEMP_DIR )      // 'YYYYMMDD.a3.'     // 
     &                 GET_RES_EXT() 

         ! Remove A-3 files for this date from temp dir (GEOS-3)
         REMOVE_DATE = TRIM( REMOVE_CMD )    // ' '                // 
     &                 TRIM( TEMP_DIR   )    // 'YYYYMMDD.a3.'     // 
     &                 GET_RES_EXT()  

#elif defined( GEOS_4 )

         ! Location of zipped A-3 file in data dir (GEOS-4)
         A3_FILE_GZ  = TRIM( DATA_DIR )      // TRIM( GEOS_4_DIR ) // 
     &                 'YYYYMM/YYYYMMDD.a3.' // GET_RES_EXT()      // 
     &                 TRIM( ZIP_SUFFIX )

         ! Location of unzipped A-3 file in temp dir (GEOS-4)
         A3_FILE     = TRIM( TEMP_DIR )      // 'YYYYMMDD.a3.'     // 
     &                 GET_RES_EXT()

         ! Remove A-3 files for this date from temp dir (GEOS-4)
         REMOVE_DATE = TRIM( REMOVE_CMD )    // ' '                // 
     &                 TRIM( TEMP_DIR   )    // 'YYYYMMDD.a3.'     // 
     &                 GET_RES_EXT()  

#endif

         !==============================================================
         ! Replace tokens in character variables w/ year & month values
         !==============================================================
         CALL EXPAND_DATE( A3_FILE_GZ,  NYMD, 000000 )
         CALL EXPAND_DATE( A3_FILE,     NYMD, 000000 )
         CALL EXPAND_DATE( REMOVE_DATE, NYMD, 000000 )

         !==============================================================
         ! Define the foreground and background UNZIP commands
         !==============================================================

         ! Foreground unzip
         UNZIP_FG = TRIM( UNZIP_CMD ) // ' ' // TRIM( A3_FILE_GZ ) // 
     &              TRIM( REDIRECT  ) // ' ' // TRIM( A3_FILE    )  

         ! Background unzip
         UNZIP_BG  = TRIM( UNZIP_FG ) // TRIM( BACKGROUND )
      ENDIF

      !=================================================================
      ! Define command to remove all A-3 files from the TEMP dir
      !=================================================================
      REMOVE_ALL = TRIM( REMOVE_CMD ) // ' '    // TRIM( TEMP_DIR   ) // 
     &             TRIM( WILD_CARD  ) // '.a3.' // TRIM( WILD_CARD  ) 

      !=================================================================
      ! Perform an F90 system call to do the desired operation
      !=================================================================
      SELECT CASE ( TRIM( OPTION ) )
         
         ! Unzip A-3 fields in the Unix foreground
         CASE ( 'unzip foreground' )
            WRITE( 6, 100 ) TRIM( A3_FILE_GZ )
            CALL SYSTEM( TRIM( UNZIP_FG ) )

         ! Unzip A-3 fields in the Unix background
         CASE ( 'unzip background' )
            WRITE( 6, 100 ) TRIM( A3_FILE_GZ )
            CALL SYSTEM( TRIM( UNZIP_BG ) )

         ! Remove A-3 field for this date in temp dir
         CASE ( 'remove date' )
            WRITE( 6, 110 ) TRIM( A3_FILE )
            CALL SYSTEM( TRIM( REMOVE_DATE ) )
            
         ! Remove all A-3 fields in temp dir
         CASE ( 'remove all' )
            WRITE( 6, 120 ) TRIM( REMOVE_ALL )
            CALL SYSTEM( TRIM( REMOVE_ALL ) )

         ! Error -- bad option!
         CASE DEFAULT
            CALL ERROR_STOP( 'Invalid value for OPTION!', 
     &                       'UNZIP_A3_FIELDS (a3_read_mod.f)' )
            
      END SELECT

      ! FORMAT strings
 100  FORMAT( '     - Unzipping: ', a )
 110  FORMAT( '     - Removing: ', a )
 120  FORMAT( '     - About to execute command: ', a )

      ! Return to calling program
      END SUBROUTINE UNZIP_A3_FIELDS

!------------------------------------------------------------------------------

      FUNCTION DO_OPEN_A3( NYMD, NHMS ) RESULT( DO_OPEN )
!
!******************************************************************************
!  Function DO_OPEN_A3 returns TRUE if is time to open the A-3 met field file
!  or FALSE otherwise.  This prevents us from opening a file which has already
!  been opened. (bmy, 6/23/03)
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
      ! DO_OPEN_A3 begins here!
      !=================================================================

      ! Initialize
      DO_OPEN = .FALSE.

      ! Return if we have already opened the file
      IF ( NYMD == LASTNYMD .and. NHMS == LASTNHMS ) THEN
         DO_OPEN = .FALSE. 
         GOTO 999
      ENDIF

#if   defined( GEOS_4 )

      ! Open A-3 file if it's 01:30 GMT,  or on the first call
      IF ( NHMS == 013000 .or. FIRST ) THEN
         DO_OPEN = .TRUE. 
         GOTO 999
      ENDIF

#else

      ! Open A-3 file if it's 00:00 GMT, or on the first call
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
      END FUNCTION DO_OPEN_A3

!------------------------------------------------------------------------------

      SUBROUTINE OPEN_A3_FIELDS( NYMD, NHMS )
!
!******************************************************************************
!  Subroutine OPEN_A3_FIELDS opens the A-3 met fields file for date NYMD and 
!  time NHMS. (bmy, bdf, 6/15/98, 6/23/03)
!  
!  Arguments as input:
!  ===========================================================================
!  (1 ) NYMD (INTEGER) : YYYYMMDD
!  (2 ) NHMS (INTEGER) :  and HHMMSS timestamps for A-3 file
!
!  NOTES:
!  (1 ) Adapted from OPEN_MET_FIELDS of "dao_read_mod.f" (bmy, 6/13/03)
!******************************************************************************
!      
      ! References to F90 modules
      USE BPCH2_MOD, ONLY : GET_RES_EXT
      USE ERROR_MOD, ONLY : ERROR_STOP
      USE FILE_MOD,  ONLY : IU_A3, IOERROR
      USE TIME_MOD,  ONLY : EXPAND_DATE

#     include "CMN_SIZE"  ! Size parameters
#     include "CMN_SETUP" ! GEOS directories

      ! Arguments
      INTEGER, INTENT(IN) :: NYMD, NHMS

      ! Local variables
      LOGICAL             :: DO_OPEN
      LOGICAL             :: IT_EXISTS
      INTEGER             :: IOS
      CHARACTER(LEN=255)  :: PATH

      !=================================================================
      ! OPEN_A3_FIELDS begins here!
      !=================================================================

      ! Open A-3 fields at the proper time, or on the first call
      IF ( DO_OPEN_A3( NYMD, NHMS ) ) THEN

#if   defined( GEOS_1 ) || defined( GEOS_STRAT )

         ! Location of A-3 file in temp dir (GEOS-1, GEOS-S)
         PATH = TRIM( TEMP_DIR )  // 'YYMMDD.a3.' // GET_RES_EXT()

#else

         ! Location of A-3 file in temp dir (GEOS-3, GEOS-4)
         PATH = TRIM( TEMP_DIR )  // 'YYYYMMDD.a3.' // GET_RES_EXT()

#endif

         ! Replace YYYYMMDD in PATH w/ actual date
         CALL EXPAND_DATE( PATH, NYMD, 000000 )

         ! Close previously opened A-3 file
         CLOSE( IU_A3 )

         ! Make sure the file exists before we open it!
         ! Maybe make this a function in ERROR_MOD (bmy, 6/23/03)
         INQUIRE( IU_A3, EXIST=IT_EXISTS )
            
         IF ( .not. IT_EXISTS ) THEN
            CALL ERROR_STOP( 'Could not find file!', 
     &                       'OPEN_A3_FIELDS (a3_read_mod.f)' )
         ENDIF

         ! Open the file
         OPEN( UNIT   = IU_A3,         FILE   = TRIM( PATH ),
     &         STATUS = 'OLD',         ACCESS = 'SEQUENTIAL',  
     &         FORM   = 'UNFORMATTED', IOSTAT = IOS )
               
         IF ( IOS /= 0 ) THEN
            CALL IOERROR( IOS, IU_A3, 'open_a3_fields:1' )
         ENDIF

         WRITE( 6, '( ''     - Opening: '', a )' ) TRIM( PATH )
      ENDIF

      ! Return to calling program
      END SUBROUTINE OPEN_A3_FIELDS

!------------------------------------------------------------------------------

      SUBROUTINE GET_A3_FIELDS( NYMD, NHMS )
!
!******************************************************************************
!  Subroutine GET_A3_FIELDS is a wrapper for routine READ_A3.  GET_A3_FIELDS
!  calls READ_A3 properly for reading GEOS-1, GEOS-STRAT, GEOS-3, or GEOS-4
!  met data sets. (bmy, 6/23/03)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) NYMD (INTEGER) : YYYYMMDD
!  (2 ) NHMS (INTEGER) :  and HHMMSS of A-3 fields to be read from disk
!
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE DAO_MOD, ONLY : ALBD, GWETTOP, HFLUX,  PARDF, PARDR,
     &                    PBL,  PREACC,  PRECON, TS,    TSKIN,   
     &                    U10M, USTAR,   V10M,   Z0

#     include "CMN_SIZE"  ! Size parameters
#     include "CMN_DEP"   ! CFRAC, RADIAT

      ! Arguments
      INTEGER, INTENT(IN) :: NYMD, NHMS 

      ! Local variables
      INTEGER, SAVE       :: LASTNYMD = -1, LASTNHMS = -1

      !=================================================================
      ! GET_A3_FIELDS begins here!
      !=================================================================

      ! Skip over previously-read A-3 fields
      IF ( NYMD == LASTNYMD .and. NHMS == LASTNHMS ) THEN
         WRITE( 6, 100 ) NYMD, NHMS
 100     FORMAT( '     - A-3 met fields for NYMD, NHMS = ', 
     &           i8.8, 1x, i6.6, ' have been read already' ) 
         RETURN
      ENDIF

#if   defined( GEOS_1 ) || defined( GEOS_3 ) 

      !=================================================================
      ! For GEOS-1 and GEOS-3, read the following fields:
      !    HFLUX, RADSWG, PREACC, PRECON, USTAR, 
      !    Z0,    PBL,    CLDFRC, U10M,   V10M
      !=================================================================
      CALL READ_A3( NYMD=NYMD,     NHMS=NHMS,     CLDFRC=CFRAC,  
     &              HFLUX=HFLUX,   PBL=PBL,       PREACC=PREACC, 
     &              PRECON=PRECON, RADSWG=RADIAT, TS=TS,
     &              U10M=U10M,     USTAR=USTAR,   V10M=V10M, 
     &              Z0=Z0 )              

#elif defined( GEOS_STRAT )

      !=================================================================
      ! For GEOS-STRAT, always read the following fields:
      !    HFLUX, RADSWG, PREACC, PRECON, USTAR, Z0, PBL   
      !
      ! Prior to 13 Feb 1997, U10M and V10M were not archived in the 
      ! GEOS-STRAT data.  If U10M and V10M are available, read them. 
      !=================================================================
      IF ( NYMD < 19970213 ) THEN

         ! Prior to 13 Feb 1997, read 7 A-3 fields
         CALL READ_A3( NYMD=NYMD,     NHMS=NHMS,     HFLUX=HFLUX,   
     &                 PBL=PBL,       PREACC=PREACC, PRECON=PRECON, 
     &                 RADSWG=RADIAT, TS=TS,         USTAR=USTAR,
     &                 Z0=Z0 )  

      ELSE
         
         ! On or after 13 Feb 1997, read 9 A-3 fields
         CALL READ_A3( NYMD=NYMD,     NHMS=NHMS,     HFLUX=HFLUX,   
     &                 PBL=PBL,       PREACC=PREACC, PRECON=PRECON, 
     &                 RADSWG=RADIAT, TS=TS,         U10M=U10M,
     &                 USTAR=USTAR,   V10M=V10M,     Z0=Z0 )         
 
      ENDIF

#elif defined( GEOS_4 )

      !================================================================
      ! For GEOS-4/fvDAS read the following fields
      !================================================================
      CALL READ_A3( NYMD=NYMD,       NHMS=NHMS,     
     &              ALBEDO=ALBD,     CLDFRC=CFRAC,  HFLUX=HFLUX,   
     &              GWETTOP=GWETTOP, PARDF=PARDF,   PARDR=PARDR,
     &              PBL=PBL,         PREACC=PREACC, PRECON=PRECON,   
     &              RADSWG=RADIAT,   TS=TS,         TSKIN=TSKIN,     
     &              U10M=U10M,       USTAR=USTAR,   V10M=V10M,       
     &              Z0=Z0 )

#endif

      ! Save NYMD, NHMS for next call
      LASTNYMD = NYMD
      LASTNHMS = NHMS

      ! Return to MAIN program
      END SUBROUTINE GET_A3_FIELDS

!------------------------------------------------------------------------------

      FUNCTION GET_N_A3( NYMD ) RESULT( N_A3 )
!
!******************************************************************************
!  Function GET_N_A3 returns the number of A-3 fields per met data set
!  (GEOS-1, GEOS-STRAT, GEOS-3, GEOS-4). (bmy, 6/23/03) 
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) NYMD (INTEGER) : YYYYMMDD for which to read in A-3 fields
!
!  NOTES:
!******************************************************************************
!
#     include "CMN_SIZE" 

      ! Arguments
      INTEGER, INTENT(IN) :: NYMD

      ! Function value
      INTEGER             :: N_A3

      !=================================================================
      ! GET_N_A3 begins here!
      !=================================================================
#if   defined( GEOS_1 ) 

      ! GEOS-1 always has 12 A-3 fields per file. 
      N_A3 = 12

#elif defined( GEOS_STRAT )

      ! GEOS-STRAT before   13 Feb 1997 has 9  A-3 fields per file
      ! GEOS-STRAT on/after 13 Feb 1997 has 11 A-3 fields per file
      IF ( NYMD < 19970213 ) THEN
         N_A3 = 9
      ELSE
         N_A3 = 11
      ENDIF

#elif defined( GEOS_3 )

      ! GEOS-3 for 1998, 1999 has 14 A-3 fields per file
      ! GEOS-3 for 2000+      has 11 A-3 fields per file
      SELECT CASE ( NYMD / 10000 )
         CASE ( 1998, 1999 )
            N_A3 = 14

         CASE DEFAULT
            N_A3 = 11
      
      END SELECT

#elif defined( GEOS_4 )

      ! GEOS-4/fvDAS has 16 fields
      N_A3 = 16

#endif

      ! Return to calling program
      END FUNCTION GET_N_A3

!---------------------------------------------------------------------------

      FUNCTION CHECK_TIME( XYMD, XHMS, NYMD, NHMS ) RESULT( ITS_TIME )
!
!******************************************************************************
!  Function CHECK_TIME checks to see if the timestamp of the A-3 field just
!  read from disk matches the current time.  If so, then it's time to return
!  the A-3 field to the calling program. (bmy, 6/23/03)
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

      SUBROUTINE READ_A3( NYMD,   NHMS, 
     &                    ALBEDO, CLDFRC, GWETTOP, HFLUX,  PARDF,
     &                    PARDR,  PBL,    PREACC,  PRECON, RADSWG,  
     &                    RADSWT, TS,     TSKIN,   U10M,   USTAR,   
     &                    V10M,   Z0 )
!
!******************************************************************************
!  Subroutine READ_A3 reads GEOS A-3 (3-hr avg) fields from disk.
!  (bmy, 5/8/98, 10/28/03)
! 
!  Arguments as input:
!  ============================================================================
!  (1 ) NYMD    : YYYYMMDD
!  (2 ) NHMS    :  and HHMMSS of A-3 met fields to be accessed 
!
!  A-3 Met Fields as Output:
!  ============================================================================
!  (1 ) ALBEDO  : (2-D) GMAO surface albedo at 10 m            [unitless]
!  (2 ) CLDFRC  : (2-D) GMAO column cloud fraction @ ground    [unitless]
!  (3 ) GWETTOP : (2-D) GMAO topsoil wetness                   [unitless]
!  (4 ) HFLUX   : (2-D) GMAO sensible heat flux                [W/m2] 
!  (5 ) PARDF   : (2-D) GMAO photosyn active diffuse radiation [W/m2]
!  (6 ) PARDR   : (2-D) GMAO photosyn active direct radiation  [W/m2]
!  (7 ) PBL     : (2-D) GMAO planetary boundary layer depth    [mb] 
!  (8 ) PREACC  : (2-D) GMAO accumulated precip @ ground       [mm H2O/day]
!  (9 ) PRECON  : (2-D) GMAO convective  precip @ ground       [mm H2O/day]
!  (10) RADSWG  : (2-D) GMAO solar insolation @ ground         [W/m2]
!  (11) RADSWT  : (2-D) GMAO solar insolation @ atm. top       [W/m2]
!  (12) TS      : (2-D) GMAO surface air temperature           [K]
!  (13) TSKIN   : (2-D) GMAO radiance temperature              [K]
!  (14) USTAR   : (2-D) GMAO friction velocity                 [m/s]
!  (13) U10M    : (2-D) GMAO U-wind at 10 m                    [m/s]
!  (14) V10M    : (2-D) GMAO V-wind at 10 m                    [m/s]
!  (15) Z0      : (2-D) GMAO roughness height                  [m] 
!
!  NOTES:
!  (1 ) Now use function TIMESTAMP_STRING from "time_mod.f" for formatted 
!        date/time output. (bmy, 10/28/03)
!******************************************************************************
!
      ! References to F90 modules
      USE DIAG_MOD,     ONLY : AD67
      USE FILE_MOD,     ONLY : IOERROR,     IU_A3
      USE TIME_MOD,     ONLY : SET_CT_A3,   TIMESTAMP_STRING
      USE TRANSFER_MOD, ONLY : TRANSFER_2D, TRANSFER_TO_1D

#     include "CMN_SIZE"             ! Size parameters
#     include "CMN_DIAG"             ! ND67

      ! Arguments
      INTEGER, INTENT(IN)            :: NYMD, NHMS
      REAL*8,  INTENT(OUT), OPTIONAL :: ALBEDO(IIPAR,JJPAR) 
      REAL*8,  INTENT(OUT), OPTIONAL :: CLDFRC(MAXIJ)
      REAL*8,  INTENT(OUT), OPTIONAL :: GWETTOP(IIPAR,JJPAR) 
      REAL*8,  INTENT(OUT), OPTIONAL :: HFLUX(IIPAR,JJPAR)
      REAL*8,  INTENT(OUT), OPTIONAL :: PARDF(IIPAR,JJPAR) 
      REAL*8,  INTENT(OUT), OPTIONAL :: PARDR(IIPAR,JJPAR) 
      REAL*8,  INTENT(OUT), OPTIONAL :: PBL(IIPAR,JJPAR) 
      REAL*8,  INTENT(OUT), OPTIONAL :: PREACC(IIPAR,JJPAR) 
      REAL*8,  INTENT(OUT), OPTIONAL :: PRECON(IIPAR,JJPAR)
      REAL*8,  INTENT(OUT), OPTIONAL :: RADSWG(MAXIJ)
      REAL*8,  INTENT(OUT), OPTIONAL :: RADSWT(IIPAR,JJPAR)
      REAL*8,  INTENT(OUT), OPTIONAL :: TS(IIPAR,JJPAR) 
      REAL*8,  INTENT(OUT), OPTIONAL :: TSKIN(IIPAR,JJPAR) 
      REAL*8,  INTENT(OUT), OPTIONAL :: U10M(IIPAR,JJPAR) 
      REAL*8,  INTENT(OUT), OPTIONAL :: USTAR(MAXIJ) 
      REAL*8,  INTENT(OUT), OPTIONAL :: V10M(IIPAR,JJPAR)
      REAL*8,  INTENT(OUT), OPTIONAL :: Z0(MAXIJ)

      ! Local Variables
      INTEGER                        :: I, IJLOOP, IOS, J, N_A3, NFOUND 
      REAL*4                         :: Q2(IGLOB,JGLOB)
      CHARACTER(LEN=8)               :: NAME
      CHARACTER(LEN=16)              :: STAMP

      ! XYMD, XHMS must be REAL*4 for GEOS-1, GEOS-STRAT
      ! but INTEGER for GEOS-3 and GEOS-4 (bmy, 6/23/03)
#if   defined( GEOS_1 ) || defined( GEOS_STRAT )
      REAL*4                         :: XYMD, XHMS 
#else
      INTEGER                        :: XYMD, XHMS
#endif

      !=================================================================
      ! READ_A3 begins here!      
      !=================================================================

      ! Get the number of A-3 fields stored in this data set
      N_A3 = GET_N_A3( NYMD )

      ! Zero the number of A-3 fields that we have found
      NFOUND = 0

      !=================================================================
      ! Read the A-3 fields from disk
      !=================================================================
      DO

         ! Read the A-3 field name
         READ( IU_A3, IOSTAT=IOS ) NAME

         ! End of file test -- make sure we have found all fields
         IF ( IOS < 0 ) THEN
            CALL A3_CHECK( NFOUND, N_A3 )
            EXIT
         ENDIF

         ! IOS > 0: True I/O error; stop w/ err msg
         IF ( IOS > 0 ) CALL IOERROR( IOS, IU_A3, 'read_a3:1' )

         ! CASE statement for A-3 fields
         SELECT CASE ( TRIM( NAME ) )

            !--------------------------------
            ! ALBEDO: surface albedo
            !--------------------------------
            CASE ( 'ALBEDO' ) 
               READ( IU_A3, IOSTAT=IOS ) XYMD, XHMS, Q2
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A3, 'read_a3:2' )

               IF ( CHECK_TIME( XYMD, XHMS, NYMD, NHMS ) ) THEN
                  IF ( PRESENT( ALBEDO ) ) CALL TRANSFER_2D( Q2,ALBEDO )
                  NFOUND = NFOUND + 1
               ENDIF

            !---------------------------------
            ! CLDFRC: column cloud fraction
            ! Store in 1-D array
            !---------------------------------
            CASE ( 'CLDFRC' )
               READ ( IU_A3, IOSTAT=IOS ) XYMD, XHMS, Q2
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A3, 'read_a3:3' )
             
               IF ( CHECK_TIME( XYMD, XHMS, NYMD, NHMS ) ) THEN
                  IF ( PRESENT( CLDFRC ) ) THEN
                     CALL TRANSFER_TO_1D( Q2, CLDFRC )
                  ENDIF
                  NFOUND = NFOUND + 1
               ENDIF

            !--------------------------------
            ! GWETTOP: topsoil wetness 
            !--------------------------------
            CASE ( 'GWETTOP' ) 
               READ( IU_A3, IOSTAT=IOS ) XYMD, XHMS, Q2
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A3, 'read_a3:4' )

               IF ( CHECK_TIME( XYMD, XHMS, NYMD, NHMS ) ) THEN
                  IF ( PRESENT( GWETTOP ) ) CALL TRANSFER_2D(Q2,GWETTOP)
                  NFOUND = NFOUND + 1
               ENDIF

            !--------------------------------
            ! HFLUX:  sensible heat flux 
            !--------------------------------
            CASE ( 'HFLUX' ) 
               READ( IU_A3, IOSTAT=IOS ) XYMD, XHMS, Q2
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A3, 'read_a3:5' )

               IF ( CHECK_TIME( XYMD, XHMS, NYMD, NHMS ) ) THEN
                  IF ( PRESENT( HFLUX ) ) CALL TRANSFER_2D( Q2, HFLUX )
                  NFOUND = NFOUND + 1
               ENDIF

            !--------------------------------
            ! PARDF: photosyn active diff rad
            !--------------------------------
            CASE ( 'PARDF' )
               READ( IU_A3, IOSTAT=IOS ) XYMD, XHMS, Q2
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A3, 'read_a3:6' )
             
               IF ( CHECK_TIME( XYMD, XHMS, NYMD, NHMS ) ) THEN
                  IF ( PRESENT( PARDF ) ) CALL TRANSFER_2D( Q2, PARDF )
                  NFOUND = NFOUND + 1
               ENDIF  

            !--------------------------------
            ! PARDR: photosyn active dir rad
            !--------------------------------
            CASE ( 'PARDR' )
               READ( IU_A3, IOSTAT=IOS ) XYMD, XHMS, Q2
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A3, 'read_a3:7' )
             
               IF ( CHECK_TIME( XYMD, XHMS, NYMD, NHMS ) ) THEN
                  IF ( PRESENT( PARDR ) ) CALL TRANSFER_2D( Q2, PARDR )
                  NFOUND = NFOUND + 1
               ENDIF 

            !--------------------------------
            ! PBL: boundary layer depth
            !--------------------------------
            CASE ( 'PBL', 'PBLH' )
               READ( IU_A3, IOSTAT=IOS ) XYMD, XHMS, Q2
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A3, 'read_a3:8' )
             
               IF ( CHECK_TIME( XYMD, XHMS, NYMD, NHMS ) ) THEN
                  IF ( PRESENT( PBL ) ) CALL TRANSFER_2D( Q2, PBL )
                  NFOUND = NFOUND + 1
               ENDIF       

            !--------------------------------
            ! PREACC: total precip at ground
            !--------------------------------
            CASE ( 'PREACC' )
               READ( IU_A3, IOSTAT=IOS ) XYMD, XHMS, Q2
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A3, 'read_a3:9' )
             
               IF ( CHECK_TIME( XYMD, XHMS, NYMD, NHMS ) ) THEN
                  IF ( PRESENT( PREACC ) ) CALL TRANSFER_2D( Q2,PREACC )
                  NFOUND = NFOUND + 1
               ENDIF

            !--------------------------------
            ! PRECON: conv precip at ground
            !--------------------------------
            CASE ( 'PRECON' )
               READ( IU_A3, IOSTAT=IOS ) XYMD, XHMS, Q2
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A3, 'read_a3:10' )
              
               IF ( CHECK_TIME( XYMD, XHMS, NYMD, NHMS ) ) THEN
                  IF ( PRESENT( PRECON ) ) CALL TRANSFER_2D( Q2,PRECON )
                  NFOUND = NFOUND + 1
               ENDIF

            !--------------------------------
            ! RADSWG: solar rad at ground 
            ! Save to 1-D array 
            !--------------------------------
            CASE ( 'RADSWG' ) 
               READ( IU_A3, IOSTAT=IOS ) XYMD, XHMS, Q2
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A3, 'read_a3:11' )
             
               IF ( CHECK_TIME( XYMD, XHMS, NYMD, NHMS ) ) THEN
                  IF ( PRESENT( RADSWG ) ) THEN
                     CALL TRANSFER_TO_1D( Q2, RADSWG )
                  ENDIF
                  NFOUND = NFOUND + 1
               ENDIF

            !--------------------------------
            ! RADSWG -- solar rad at atm top
            !--------------------------------
            CASE ( 'RADSWT' )
               READ( IU_A3, IOSTAT=IOS ) XYMD, XHMS, Q2
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A3, 'read_a3:12' )
             
               IF ( CHECK_TIME( XYMD, XHMS, NYMD, NHMS ) ) THEN
                  IF ( PRESENT( RADSWT ) ) CALL TRANSFER_2D( Q2,RADSWT )
                  NFOUND = NFOUND + 1
               ENDIF

            !--------------------------------            
            ! TS: surface air temperature
            !--------------------------------
            CASE ( 'TS', 'TGROUND', 'T2M' )
               READ( IU_A3, IOSTAT=IOS ) XYMD, XHMS, Q2
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A3, 'read_a3:13' )
             
               IF ( CHECK_TIME( XYMD, XHMS, NYMD, NHMS ) ) THEN
                  IF ( PRESENT( TS ) ) CALL TRANSFER_2D( Q2, TS )
                  NFOUND = NFOUND + 1
               ENDIF

            !--------------------------------            
            ! TSKIN: radiance temperature
            !--------------------------------
            CASE ( 'TSKIN' )
               READ( IU_A3, IOSTAT=IOS ) XYMD, XHMS, Q2
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A3, 'read_a3:14' )
             
               IF ( CHECK_TIME( XYMD, XHMS, NYMD, NHMS ) ) THEN
                  IF ( PRESENT( TSKIN ) ) CALL TRANSFER_2D( Q2, TSKIN )
                  NFOUND = NFOUND + 1
               ENDIF

            !--------------------------------            
            ! U10M -- U-wind at 10 m
            !--------------------------------            
            CASE ( 'U10M' )
               READ( IU_A3, IOSTAT=IOS ) XYMD, XHMS, Q2
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A3, 'read_a3:15' )
             
               IF ( CHECK_TIME( XYMD, XHMS, NYMD, NHMS ) ) THEN
                  IF ( PRESENT( U10M ) ) CALL TRANSFER_2D( Q2, U10M )
                  NFOUND = NFOUND + 1
               ENDIF

            !--------------------------------            
            ! USTAR -- friction velocity
            ! Save into 1-D array
            !--------------------------------            
            CASE ( 'USTAR' )
               READ( IU_A3, IOSTAT=IOS ) XYMD, XHMS, Q2
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A3, 'read_a3:16' )
             
               IF ( CHECK_TIME( XYMD, XHMS, NYMD, NHMS ) ) THEN
                  IF ( PRESENT( USTAR ) ) THEN
                     CALL TRANSFER_TO_1D( Q2, USTAR )
                  ENDIF
                  NFOUND = NFOUND + 1
               ENDIF
            
            !--------------------------------            
            ! V10M -- V-wind at 10 m
            !--------------------------------            
            CASE ( 'V10M' )
               READ( IU_A3, IOSTAT=IOS ) XYMD, XHMS, Q2
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A3, 'read_a3:17' )
             
               IF ( CHECK_TIME( XYMD, XHMS, NYMD, NHMS ) ) THEN
                  IF ( PRESENT( V10M ) ) CALL TRANSFER_2D( Q2, V10M )
                  NFOUND = NFOUND + 1
               ENDIF

            !--------------------------------            
            ! ZO: roughness heights
            ! Save into 1-D array
            !--------------------------------            
            CASE ( 'Z0', 'Z0M' )
               READ( IU_A3, IOSTAT=IOS ) XYMD, XHMS, Q2
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A3, 'read_a3:18' )
             
               IF ( CHECK_TIME( XYMD, XHMS, NYMD, NHMS ) ) THEN
                  IF ( PRESENT( Z0 ) ) THEN
                     CALL TRANSFER_TO_1D( Q2, Z0 )
                  ENDIF
                  NFOUND = NFOUND + 1
               ENDIF
  
            !--------------------------------
            ! TPW: just skip over this
            !--------------------------------
            CASE ( 'TPW' )
               READ( IU_A3, IOSTAT=IOS ) XYMD, XHMS, Q2
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A3, 'read_a3:19' )

               IF ( CHECK_TIME( XYMD, XHMS, NYMD, NHMS ) ) THEN
                  NFOUND = NFOUND + 1
               ENDIF

            !--------------------------------
            ! CLDTMP: just skip over this
            !--------------------------------
            CASE ( 'CLDTMP' )
               READ( IU_A3, IOSTAT=IOS ) XYMD, XHMS, Q2
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A3, 'read_a3:20' )
             
               IF ( CHECK_TIME( XYMD, XHMS, NYMD, NHMS ) ) THEN
                  NFOUND = NFOUND + 1
               ENDIF

         END SELECT
               
         !==============================================================
         ! If we have found all the fields for this time, then exit 
         ! the loop.  Otherwise, go on to the next iteration.
         !==============================================================
         IF ( CHECK_TIME( XYMD, XHMS, NYMD, NHMS ) .and. 
     &        NFOUND == N_A3 ) THEN 
!----------------------------------------------------------------------------
! Prior to 10/28/03:
! Now use TIMESTAMP_STRING for formatted output (bmy, 10/28/03)
!            WRITE( 6, 210 ) NFOUND, NYMD, NHMS
! 210        FORMAT( '     - Found all ', i3, 
!     &              ' A-3 met fields for NYMD, NHMS = ', 
!     &              i8.8, 1x, i6.6 )
!----------------------------------------------------------------------------
            STAMP = TIMESTAMP_STRING( NYMD, NHMS )
            WRITE( 6, 210 ) NFOUND, STAMP
 210        FORMAT( '     - Found all ', i3, ' A-3 met fields for ', a )
            EXIT
         ENDIF
      ENDDO

      !=================================================================
      ! ND67 diagnostic: A-3 surface fields:
      !
      ! (1 ) HFLUX  : sensible heat flux from surface    [W/m2]  
      ! (2 ) RADSWG : solar radiation @ ground           [W/m2]  
      ! (3 ) PREACC : total precip. @ ground             [mm/day] 
      ! (4 ) PRECON : conv.  precip. @ ground            [mm/day] 
      ! (5 ) TS     : surface air temperature            [K]    
      ! (6 ) RADSWT : solar radiation @ atm. top         [W/m2]  
      ! (7 ) USTAR  : friction velocity                  [m/s]   
      ! (8 ) Z0     : surface roughness height           [m]    
      ! (9 ) PBL    : planetary boundary layer depth     [hPa]   
      ! (10) CLDFRC : column cloud fraction              [unitless]  
      ! (11) U10M   : U-winds @ 10 meters altitude       [m/s]   
      ! (12) V10M   : V-winds @ 10 meters altitude       [m/s]   
      ! (13) PS-PBL : Boundary Layer Top Pressure        [hPa]
      ! (14) ALBD   : Surface Albedo                     [unitless]
      ! (15) PHIS   : Geopotential Heights               [m]      
      ! (16) CLTOP  : Cloud Top Height                   [levels]
      ! (17) TROPP  : Tropopause pressure                [hPa] 
      ! (18) SLP    : Sea Level pressure                 [hPa]
      ! (19) TSKIN  : Ground/sea surface temp.           [hPa]  
      ! (20) PARDF  : Photosyn active diffuse radiation  [W/m2]
      ! (21) PARDR  : Photosyn active direct  radiation  [W/m2]
      ! (22) GWET   : Top soil wetness                   [unitless]
      !=================================================================
      IF ( ND67 > 0 ) THEN

         ! Archive 2-D arrays
         IF ( PRESENT( HFLUX   ) ) AD67(:,:,1 ) = AD67(:,:,1 ) + HFLUX
         IF ( PRESENT( PREACC  ) ) AD67(:,:,3 ) = AD67(:,:,3 ) + PREACC
         IF ( PRESENT( PRECON  ) ) AD67(:,:,4 ) = AD67(:,:,4 ) + PRECON
         IF ( PRESENT( TS      ) ) AD67(:,:,5 ) = AD67(:,:,5 ) + TS
         IF ( PRESENT( RADSWT  ) ) AD67(:,:,6 ) = AD67(:,:,6 ) + RADSWT
         IF ( PRESENT( PBL     ) ) AD67(:,:,9 ) = AD67(:,:,9 ) + PBL
         IF ( PRESENT( U10M    ) ) AD67(:,:,11) = AD67(:,:,11) + U10M
         IF ( PRESENT( V10M    ) ) AD67(:,:,12) = AD67(:,:,12) + V10M
         IF ( PRESENT( ALBEDO  ) ) AD67(:,:,14) = AD67(:,:,14) + ALBEDO
         IF ( PRESENT( TSKIN   ) ) AD67(:,:,19) = AD67(:,:,19) + TSKIN
         IF ( PRESENT( PARDF   ) ) AD67(:,:,20) = AD67(:,:,20) + PARDF
         IF ( PRESENT( PARDR   ) ) AD67(:,:,21) = AD67(:,:,21) + PARDR
         IF ( PRESENT( GWETTOP ) ) AD67(:,:,22) = AD67(:,:,22) + GWETTOP

         ! Archive 1-D arrays
         IF ( PRESENT( RADSWG  ) ) CALL ARCHIVE_ND67_1D( RADSWG,  2 ) 
         IF ( PRESENT( USTAR   ) ) CALL ARCHIVE_ND67_1D( USTAR,   7 ) 
         IF ( PRESENT( Z0      ) ) CALL ARCHIVE_ND67_1D( Z0,      8 ) 
         IF ( PRESENT( CLDFRC  ) ) CALL ARCHIVE_ND67_1D( CLDFRC, 10 ) 

      ENDIF
         
      ! Increment KDA3FLDS -- this is the # of times READ_A3 is called
      CALL SET_CT_A3( INCREMENT=.TRUE. )

      ! Return to calling program
      END SUBROUTINE READ_A3

!------------------------------------------------------------------------------

      SUBROUTINE ARCHIVE_ND67_1D( FIELD, N )
!
!******************************************************************************
!  Subroutine ARCHIVE_ND67_1D saves 1-D arrays for the ND67 diagnostic.
!  (bmy, 6/23/03)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) FIELD (REAL*8 ) : Met field array to be archived into ND67
!  (2 ) N     (INTEGER) : Index of AD67 array under which to archive data
!
!  NOTES
!******************************************************************************
!
      ! References to F90 modules
      USE DIAG_MOD, ONLY : AD67

#     include "CMN_SIZE"  ! Size parameters

      ! Arguments
      INTEGER, INTENT(IN) :: N
      REAL*8,  INTENT(IN) :: FIELD(MAXIJ)

      ! Local variables
      INTEGER             :: I, IJLOOP, J

      !=================================================================
      ! ARCHIVE_1D begins here
      !=================================================================

      ! Zero 1-D counter
      IJLOOP = 0

      ! Archive for diagnostic
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         IJLOOP      = IJLOOP + 1
         AD67(I,J,N) = AD67(I,J,N) + FIELD(IJLOOP)
      ENDDO
      ENDDO

      ! Return to calling program
      END SUBROUTINE ARCHIVE_ND67_1D

!------------------------------------------------------------------------------

      SUBROUTINE A3_CHECK( NFOUND, N_A3 )
!
!******************************************************************************
!  Subroutine A3_CHECK prints an error message if not all of the A-3 met 
!  fields are found.  The run is also terminated. (bmy, 10/27/00, 6/23/03)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) NFOUND (INTEGER) : # of A-3 met fields read from disk
!  (2 ) N_A3   (INTEGER) : # of A-3 met fields expected to be read from disk
!
!  NOTES
!  (1 ) Adapted from DAO_CHECK from "dao_read_mod.f" (bmy, 6/23/03)
!******************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD, ONLY : GEOS_CHEM_STOP

      ! Arguments
      INTEGER, INTENT(IN) :: NFOUND, N_A3

      !=================================================================
      ! A3_CHECK begins here!
      !=================================================================
      IF ( NFOUND /= N_A3 ) THEN
         WRITE( 6, '(a)' ) REPEAT( '=', 79 )
         WRITE( 6, '(a)' ) 'ERROR -- not enough A-3 fields found!'      

         WRITE( 6, 120   ) N_A3, NFOUND
 120     FORMAT( 'There are ', i2, ' fields but only ', i2 ,
     &           ' were found!' )

         WRITE( 6, '(a)' ) '### STOP in A3_CHECK (dao_read_mod.f)'
         WRITE( 6, '(a)' ) REPEAT( '=', 79 )

         ! Deallocate arrays and stop (bmy, 10/15/02)
         CALL GEOS_CHEM_STOP
      ENDIF

      ! Return to calling program
      END SUBROUTINE A3_CHECK

!------------------------------------------------------------------------------

      END MODULE A3_READ_MOD
