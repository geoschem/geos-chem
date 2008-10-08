! $Id: a3_read_mod.f,v 1.21 2008/10/08 18:30:33 bmy Exp $
      MODULE A3_READ_MOD
!
!******************************************************************************
!  Module A3_READ_MOD contains routines that unzip, open, and read the
!  GEOS-Chem A-3 (avg 3-hour) met fields from disk. (bmy, 6/23/03, 10/7/08)
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
!  (4 ) directory_mod.f : Module containing GEOS-CHEM data & met field dirs
!  (5 ) error_mod.f     : Module containing NaN and other error check routines
!  (6 ) logical_mod.f   : Module containing GEOS-CHEM logical switches 
!  (7 ) file_mod.f      : Module containing file unit #'s and error checks
!  (8 ) time_mod.f      : Module containing routines for computing time & date
!  (9 ) transfer_mod.f  : Module containing routines to cast & resize arrays
!  (10) unix_cmds_mod.f : Module containing Unix commands for unzipping etc.
!
!  NOTES:
!  (1 ) Adapted from "dao_read_mod.f" (bmy, 6/23/03)
!  (2 ) Now can read from either zipped or unzipped files. (bmy, 12/11/03)
!  (3 ) Now skips past the GEOS-4 met field ident string (bmy, 12/12/03)
!  (4 ) Now references "unix_cmds_mod.f", "directory_mod.f", and 
!        "logical_mod.f" (bmy, 7/20/04)
!  (5 ) Now references FILE_EXISTS from "file_mod.f" (bmy, 3/23/05)
!  (6 ) Now modified for GEOS-5 and GCAP met fields (bmy, 5/25/05)
!  (7 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (8 ) Fixed typos for GCAP fields and ND67 diagnostics (bmy, 2/9/06)
!  (9 ) Remove support for GEOS-1 and GEOS-STRAT met fields (bmy, 8/4/06)
!  (10) Now read PARDF, PARDR for GCAP met fields (swu, bmy, 10/4/06)
!  (11) Extra modifications for GEOS-5 met fields (bmy, 1/17/07)
!  (12) Now get the # of A-3 fields from the file ident string (bmy, 10/7/08)
!******************************************************************************
!
      IMPLICIT NONE

      !=================================================================
      ! MODULE PRIVATE DECLARATIONS -- keep certain internal variables 
      ! and routines from being seen outside "a3_read_mod.f"
      !=================================================================

      ! Make everything PRIVATE ...
      PRIVATE

      ! ... except these routines
      PUBLIC :: ARCHIVE_ND67_1D 
      PUBLIC :: GET_A3_FIELDS
      PUBLIC :: OPEN_A3_FIELDS
      PUBLIC :: UNZIP_A3_FIELDS

      !=================================================================
      ! MODULE VARIABLES
      !=================================================================

      ! Number of A3 fields in the file
      INTEGER :: N_A3_FIELDS

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
!  files are not disturbed.  (bmy, bdf, 6/15/98, 8/4/06)
!
!  Arguments as input:
!  ===========================================================================
!  (1 ) OPTION (CHAR*(*)) : Option
!  (2 ) NYMD   (INTEGER ) : YYYYMMDD of A-3 file to be unzipped (optional)
!
!  NOTES:
!  (1 ) Adapted from UNZIP_MET_FIELDS of "dao_read_mod.f" (bmy, 6/23/03)
!  (2 ) Directory information YYYY/MM or YYYYMM is now contained w/in 
!        GEOS_1_DIR, GEOS_S_DIR, GEOS_3_DIR, GEOS_4_DIR (bmy, 12/11/03)
!  (3 ) Now reference "directory_mod.f" and "unix_cmds_mod.f". Now prevent 
!        EXPAND_DATE from overwriting directory paths with Y/M/D tokens in 
!        them (bmy, 7/20/04)
!  (4 ) Now modified for GEOS-5 and GCAP met fields (swu, bmy, 5/25/05)
!  (5 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (6 ) Remove support for GEOS-1 and GEOS-STRAT met fields (bmy, 8/4/06)
!*****************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD,     ONLY : GET_RES_EXT
      USE DIRECTORY_MOD, ONLY : DATA_DIR,   GCAP_DIR,   GEOS_3_DIR 
      USE DIRECTORY_MOD, ONLY : GEOS_4_DIR, GEOS_5_DIR, TEMP_DIR 
      USE ERROR_MOD,     ONLY : ERROR_STOP
      USE TIME_MOD,      ONLY : EXPAND_DATE
      USE UNIX_CMDS_MOD, ONLY : BACKGROUND, REDIRECT,   REMOVE_CMD 
      USE UNIX_CMDS_MOD, ONLY : UNZIP_CMD,  WILD_CARD,  ZIP_SUFFIX

#     include "CMN_SIZE"

      ! Arguments
      CHARACTER(LEN=*),  INTENT(IN) :: OPTION
      INTEGER, OPTIONAL, INTENT(IN) :: NYMD

      ! Local variables
      CHARACTER(LEN=255)            :: A3_STR,     GEOS_DIR
      CHARACTER(LEN=255)            :: A3_FILE_GZ, A3_FILE
      CHARACTER(LEN=255)            :: UNZIP_BG,   UNZIP_FG
      CHARACTER(LEN=255)            :: REMOVE_ALL, REMOVE_DATE

      !=================================================================
      ! UNZIP_A3_FIELDS begins here!
      !=================================================================
      IF ( PRESENT( NYMD ) ) THEN

#if   defined( GEOS_3 )

         ! Strings for directory & filename
         GEOS_DIR = TRIM( GEOS_3_DIR )
         A3_STR   = 'YYYYMMDD.a3.' // GET_RES_EXT()

#elif defined( GEOS_4 )

         ! Strings for directory & filename
         GEOS_DIR = TRIM( GEOS_4_DIR )
         A3_STR   = 'YYYYMMDD.a3.' // GET_RES_EXT()

#elif defined( GEOS_5 )

         ! Strings for directory & filename
         GEOS_DIR = TRIM( GEOS_5_DIR )
         A3_STR   = 'YYYYMMDD.a3.' // GET_RES_EXT()

#elif defined( GCAP )

         ! Strings for directory & filename
         GEOS_DIR = TRIM( GCAP_DIR )
         A3_STR   = 'YYYYMMDD.a3.' // GET_RES_EXT()

#endif

         ! Replace date tokens
         CALL EXPAND_DATE( GEOS_DIR, NYMD, 000000 )
         CALL EXPAND_DATE( A3_STR,   NYMD, 000000 )

         ! Location of zipped A-3 file in data dir 
         A3_FILE_GZ  = TRIM( DATA_DIR   ) // TRIM( GEOS_DIR   ) // 
     &                 TRIM( A3_STR     ) // TRIM( ZIP_SUFFIX )

         ! Location of unzipped A-3 file in temp dir 
         A3_FILE     = TRIM( TEMP_DIR   ) // TRIM( A3_STR     )     
         
         ! Remove A-3 files for this date from temp dir
         REMOVE_DATE = TRIM( REMOVE_CMD ) // ' '                // 
     &                 TRIM( TEMP_DIR   ) // TRIM( A3_STR     )   

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
!  been opened. (bmy, 6/23/03, 5/25/05)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) NYMD (INTEGER) : YYYYMMDD 
!  (2 ) NHMS (INTEGER) :  and HHMMSS to be tested for A-3 file open
!
!  NOTES:
!  (1 ) Now modified for GEOS-5 and GCAP met fields (swu, bmy, 5/25/05)
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

#if   defined( GEOS_3 ) 

      ! Open A-3 file if it's 00:00 GMT, or on the first call
      IF ( NHMS == 000000 .or. FIRST ) THEN
         DO_OPEN = .TRUE. 
         GOTO 999
      ENDIF

#else

      ! Open A-3 file if it's 01:30 GMT,  or on the first call
      IF ( NHMS == 013000 .or. FIRST ) THEN
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
!  time NHMS. (bmy, bdf, 6/15/98, 10/7/08)
!  
!  Arguments as input:
!  ===========================================================================
!  (1 ) NYMD (INTEGER) : YYYYMMDD
!  (2 ) NHMS (INTEGER) :  and HHMMSS timestamps for A-3 file
!
!  NOTES:
!  (1 ) Adapted from OPEN_MET_FIELDS of "dao_read_mod.f" (bmy, 6/13/03)
!  (2 ) Now opens either zipped or unzipped files (bmy, 12/11/03)
!  (3 ) Now skips past the GEOS-4 ident string (bmy, 12/12/03)
!  (4 ) Now references "directory_mod.f" instead of CMN_SETUP.  Also now
!        references LUNZIP from "logical_mod.f".  Also now prevents EXPAND_DATE
!        from overwriting Y/M/D tokens in directory paths. (bmy, 7/20/04)
!  (5 ) Now use FILE_EXISTS from "file_mod.f" to determine if file unit IU_A3
!        refers to a valid file on disk (bmy, 3/23/05)
!  (6 ) Now modified for GEOS-5 and GCAP met fields (swu, bmy, 5/25/05)  
!  (7 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (8 ) Remove support for GEOS-1 and GEOS-STRAT met fields (bmy, 8/4/06)
!  (9 ) Now get the # of A-3 fields from the file ident string (bmy, 10/7/08)
!******************************************************************************
!      
      ! References to F90 modules
      USE BPCH2_MOD,     ONLY : GET_RES_EXT
      USE DIRECTORY_MOD, ONLY : DATA_DIR,   GCAP_DIR,   GEOS_3_DIR 
      USE DIRECTORY_MOD, ONLY : GEOS_4_DIR, GEOS_5_DIR, TEMP_DIR 
      USE LOGICAL_MOD,   ONLY : LUNZIP
      USE ERROR_MOD,     ONLY : ERROR_STOP
      USE FILE_MOD,      ONLY : IU_A3,      IOERROR,    FILE_EXISTS
      USE TIME_MOD,      ONLY : EXPAND_DATE

#     include "CMN_SIZE"      ! Size parameters

      ! Arguments
      INTEGER, INTENT(IN)    :: NYMD, NHMS

      ! Local variables
      LOGICAL                :: IT_EXISTS
      INTEGER                :: IOS
      CHARACTER(LEN=2)       :: DUM
      CHARACTER(LEN=8)       :: IDENT
      CHARACTER(LEN=255)     :: A3_FILE
      CHARACTER(LEN=255)     :: GEOS_DIR
      CHARACTER(LEN=255)     :: PATH

      !=================================================================
      ! OPEN_A3_FIELDS begins here!
      !=================================================================
 
      ! Open A-3 fields at the proper time, or on the first call
      IF ( DO_OPEN_A3( NYMD, NHMS ) ) THEN

#if   defined( GEOS_3 )

         ! Strings for directory & filename
         GEOS_DIR = TRIM( GEOS_3_DIR )
         A3_FILE  = 'YYYYMMDD.a3.' // GET_RES_EXT()

#elif defined( GEOS_4 )

         ! Strings for directory & filename
         GEOS_DIR = TRIM( GEOS_4_DIR )
         A3_FILE  = 'YYYYMMDD.a3.' // GET_RES_EXT()

#elif defined( GEOS_5 )

         ! Strings for directory & filename
         GEOS_DIR = TRIM( GEOS_5_DIR )
         A3_FILE  = 'YYYYMMDD.a3.' // GET_RES_EXT()

#elif defined( GCAP )

         ! Strings for directory & filename
         GEOS_DIR = TRIM( GCAP_DIR )
         A3_FILE  = 'YYYYMMDD.a3.' // GET_RES_EXT()

#endif

         ! Replace date tokens
         CALL EXPAND_DATE( A3_FILE,  NYMD, NHMS )
         CALL EXPAND_DATE( GEOS_DIR, NYMD, NHMS )

         ! If unzipping, open GEOS-4 file in TEMP dir
         ! If not unzipping, open GEOS-4 file in DATA dir
         IF ( LUNZIP ) THEN
            PATH = TRIM( TEMP_DIR ) // TRIM( A3_FILE )
         ELSE
            PATH = TRIM( DATA_DIR ) // 
     &             TRIM( GEOS_DIR ) // TRIM( A3_FILE )
         ENDIF

         ! Close previously opened A-3 file
         CLOSE( IU_A3 )

         ! Make sure the file unit is valid before we open the file
         IF ( .not. FILE_EXISTS( IU_A3 ) ) THEN
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

         ! Echo info
         WRITE( 6, 100 ) TRIM( PATH )
 100     FORMAT( '     - Opening: ', a )
         
#if   defined( GEOS_4 ) || defined( GEOS_5 ) || defined( GCAP )

         ! Skip past the ident string
         READ( IU_A3, IOSTAT=IOS ) IDENT

         IF ( IOS /= 0 ) THEN
            CALL IOERROR( IOS, IU_A3, 'open_a3_fields:2' )
         ENDIF

         ! The last 2 digits of the ident string
         ! is the # of fields contained in the file
         READ( IDENT(7:8), '(i2.2)' ) N_A3_FIELDS        
        
#endif

      ENDIF

      ! Return to calling program
      END SUBROUTINE OPEN_A3_FIELDS

!------------------------------------------------------------------------------

      SUBROUTINE GET_A3_FIELDS( NYMD, NHMS )
!
!******************************************************************************
!  Subroutine GET_A3_FIELDS is a wrapper for routine READ_A3.  GET_A3_FIELDS
!  calls READ_A3 properly for reading GEOS-3, GEOS-4, GEOS-5, or GCAP
!  met data sets. (bmy, 6/23/03, 1/17/07)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) NYMD (INTEGER) : YYYYMMDD
!  (2 ) NHMS (INTEGER) :  and HHMMSS of A-3 fields to be read from disk
!
!  NOTES:
!  (1 ) Now save RADSWG to the RADSWG array (instead of RADIAT).  Now save
!        CLDFRC to the CLDFRC array (instead of CFRAC).  Now get RADLWG, 
!        SNOW arrays.  Also updated comments. (bmy, 12/9/03)
!  (2 ) Now modified for GEOS-5 and GCAP met fields (swu, bmy, 5/25/05)
!  (3 ) Bug fix: replace RADSWG in call to READ_A3 for GCAP met fields.
!        (bmy, 2/9/06)
!  (4 ) Remove support for GEOS-1 and GEOS-STRAT met fields (bmy, 8/4/06)
!  (5 ) Now read PARDF, PARDR for GCAP met fields (swu, bmy, 10/4/06)
!  (6 ) Now read SNOW and GETWETTOP for GCAP met fields (swu, phs, 11/15/06)
!  (7 ) Now read extra fields for GEOS-5 (bmy, 1/17/07)
!******************************************************************************
!
      ! References to F90 modules
      USE DAO_MOD, ONLY : ALBD,    CLDFRC, EVAP,   GRN,      GWETROOT 
      USE DAO_MOD, ONLY : GWETTOP, HFLUX,  LAI,    MOLENGTH, OICE    
      USE DAO_MOD, ONLY : PARDF,   PARDR,  PBL,    PREACC,   PRECON
      USE DAO_MOD, ONLY : PRECSNO, RADLWG, RADSWG, SNICE,    SNOMAS
      USE DAO_MOD, ONLY : SNODP,   SNOW,   TROPP,  TS,       TSKIN
      USE DAO_MOD, ONLY : U10M,    USTAR,  V10M,   Z0

#     include "CMN_SIZE"  ! Size parameters

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

#if   defined( GEOS_3 ) 

      !=================================================================
      ! For GEOS-3, read the following fields:
      !    HFLUX, RADSWG, PREACC, PRECON, USTAR, 
      !    Z0,    PBL,    CLDFRC, U10M,   V10M
      !=================================================================
      CALL READ_A3( NYMD=NYMD,     NHMS=NHMS,     CLDFRC=CLDFRC,  
     &              HFLUX=HFLUX,   PBL=PBL,       PREACC=PREACC, 
     &              PRECON=PRECON, RADSWG=RADSWG, TS=TS,
     &              U10M=U10M,     USTAR=USTAR,   V10M=V10M, 
     &              Z0=Z0 )              

#elif defined( GEOS_4 )

      !================================================================
      ! For GEOS-4, read the following fields:
      !
      !    ALBEDO, CLDFRC, HFLUX,  GWETTOP, PARDF,  PARDR, 
      !    PBLH,   PREACC, PRECON, RADLWG,  RADSWG, SNOW,  
      !    T2M,    TSKIN,  U10M,   USTAR,   V10M,   Z0 
      !
      ! NOTES: 
      ! (1 ) ALBEDO is an A-3 field in GEOS-4.
      ! (2 ) T2M is used as a proxy for TS in GEOS-4.
      !================================================================
      CALL READ_A3( NYMD=NYMD,       NHMS=NHMS,     
     &              ALBEDO=ALBD,     CLDFRC=CLDFRC, HFLUX=HFLUX,   
     &              GWETTOP=GWETTOP, PARDF=PARDF,   PARDR=PARDR,  
     &              PBL=PBL,         PREACC=PREACC, PRECON=PRECON, 
     &              RADLWG=RADLWG,   RADSWG=RADSWG, SNOW=SNOW,  
     &              TS=TS,           TSKIN=TSKIN,   U10M=U10M,       
     &              USTAR=USTAR ,    V10M=V10M,     Z0=Z0 )

#elif defined( GEOS_5 )

      !=================================================================
      ! For GEOS-5, read the following fields:
      !
      !    ALBEDO, CLDFRC, EVAP,   GRN,    GWETROOT, GWETTOP, HFLUX, 
      !    LAI,    PARDF,  PARDR,  PBLH,   PREACC,   PRECON,  PRECSNO,
      !    RADLWG, RADSWG, SNODP,  SNOMAS, T2M,      TROPP,   TSKIN,   
      !    U10M,   USTAR,  V10M,   Z0 
      !
      !=================================================================
      CALL READ_A3( NYMD=NYMD,         NHMS=NHMS,       ALBEDO=ALBD,   
     &              CLDFRC=CLDFRC,     EVAP=EVAP,       GRN=GRN,       
     &              GWETROOT=GWETROOT, GWETTOP=GWETTOP, HFLUX=HFLUX,   
     &              LAI=LAI,           PARDF=PARDF,     PARDR=PARDR,       
     &              PBL=PBL,           PREACC=PREACC,   PRECON=PRECON, 
     &              PRECSNO=PRECSNO,   RADLWG=RADLWG,   RADSWG=RADSWG, 
     &              SNODP=SNODP,       SNOMAS=SNOMAS,   TROPP=TROPP,   
     &              TS=TS,             TSKIN=TSKIN,     U10M=U10M,     
     &              USTAR=USTAR,       V10M=V10M,       Z0=Z0 )

#elif defined( GCAP )

      !================================================================
      ! For GCAP, read the following fields:
      !
      !    ALBEDO, GETWETTOP, MOLENGTH, OICE,  PARDF, PARDR, PBL,
      !    PREACC, PRECON,    RADSWG,   SNICE, SNOW,  TS,    U10M,  
      !    USTAR,  V10M
      !
      ! NOTES: 
      !================================================================
      CALL READ_A3( NYMD=NYMD,     NHMS=NHMS,     
     &              ALBEDO=ALBD,   GWETTOP=GWETTOP,   MOLENGTH=MOLENGTH,      
     &              OICE=OICE,     PARDF=PARDF,       PARDR=PARDR,
     &              PBL=PBL,       PREACC=PREACC,     PRECON=PRECON,
     &              RADSWG=RADSWG, SNICE=SNICE,       SNOW=SNOW, 
     &              TS=TS,         U10M=U10M,         USTAR=USTAR,
     &              V10M=V10M  )

      ! Convert [m] to [mm]
      SNOW = SNOW*1000.D0

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
!  (bmy, 6/23/03, 1/17/07) 
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) NYMD (INTEGER) : YYYYMMDD for which to read in A-3 fields
!
!  NOTES:
!  (1 ) GEOS-4/fvDAS now has 19 A-3 fields; we added LAI, RADLWG, SNOW.
!        (bmy, 12/9/03)
!  (2 ) Now modified for GEOS-5 and GCAP met fields (bmy, 5/25/05)
!  (3 ) Remove support for GEOS-1 and GEOS-STRAT met fields (bmy, 8/4/06)
!  (4 ) Increase # of fields for GCAP from 12 to 16 (swu, bmy, 10/4/06)
!  (5 ) Increase # of fields for GEOS-5 to 25 (bmy, 1/17/07)
!******************************************************************************
!
#     include "CMN_SIZE"   ! Size parameters

      ! Arguments
      INTEGER, INTENT(IN) :: NYMD

      ! Function value
      INTEGER             :: N_A3

      !=================================================================
      ! GET_N_A3 begins here!
      !=================================================================
#if   defined( GEOS_3 )

      ! GEOS-3 for 1998, 1999 has 14 A-3 fields per file
      ! GEOS-3 for 2000+      has 11 A-3 fields per file
      SELECT CASE ( NYMD / 10000 )
         CASE ( 1998, 1999 )
            N_A3 = 14

         CASE DEFAULT
            N_A3 = 11
      
      END SELECT

#elif defined( GEOS_4 )

      ! GEOS-4 has 19 A-3 fields
      N_A3 = 19

#elif defined( GEOS_5 )

      ! GEOS-5 has 25 A-3 fields
      N_A3 = 25

#elif defined( GCAP )
      
      ! GCAP has 16 fields
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
!  the A-3 field to the calling program. (bmy, 6/23/03, 8/4/06)
!  
!  Arguments as Input:
!  ============================================================================
!  (1 ) XYMD (INTEGER) : YYYYMMDD timestamp for A-3 field in file
!  (2 ) XHMS (INTEGER) : HHMMSS   timestamp for A-3 field in file
!  (3 ) NYMD (INTEGER) : YYYYMMDD at which A-3 field is to be read
!  (4 ) NHMS (INTEGER) : HHMMSS   at which A-3 field is to be read
!
!  NOTES:
!  (1 ) Remove support for GEOS-1 and GEOS-STRAT met fields (bmy, 8/4/06)
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

      SUBROUTINE READ_A3( NYMD,    NHMS, 
     &                    ALBEDO,  CLDFRC, EVAP,   GRN,      GWETROOT,
     &                    GWETTOP, HFLUX,  LAI,    MOLENGTH, OICE,   
     &                    PARDF,   PARDR,  PBL,    PREACC,   PRECON,  
     &                    PRECSNO, RADLWG, RADSWG, RADSWT,   SNICE,    
     &                    SNODP,   SNOMAS, SNOW,   TROPP,    TS,       
     &                    TSKIN,   U10M,   USTAR,  V10M,     Z0 )
!
!******************************************************************************
!  Subroutine READ_A3 reads GEOS A-3 (3-hr avg) fields from disk.
!  (bmy, 5/8/98, 10/7/08)
! 
!  Arguments as input:
!  ============================================================================
!  (1 ) NYMD     : YYYYMMDD
!  (2 ) NHMS     :  and HHMMSS of A-3 met fields to be accessed 
!
!  A-3 Met Fields as Output:
!  ============================================================================
!  (1 ) ALBEDO   : (2-D) GMAO surface albedo at 10 m            [unitless]
!  (2 ) CLDFRC   : (2-D) GMAO column cloud fraction @ ground    [unitless]
!  (3 ) EVAP     : (2-D) GMAO evapotranspiration flux 
!  (4 ) GRN      : (2-D) GMAO greenness index
!  (5 ) GWETROOT : (2-D) GMAO root soil wetness                 [unitless]
!  (6 ) GWETTOP  : (2-D) GMAO topsoil wetness                   [unitless]
!  (7 ) HFLUX    : (2-D) GMAO sensible heat flux                [W/m2] 
!  (8 ) LAI      : (2-D) GMAO leaf area index                   [m2/m2]
!  (9 ) MOLENGTH : (2-D) GCAP Monin-Obhukov length              [m]
!  (10) OICE     : (2-D) GCAP fraction of ocean ice             [unitless]
!  (11) PARDF    : (2-D) GMAO photosyn active diffuse radiation [W/m2]
!  (12) PARDR    : (2-D) GMAO photosyn active direct radiation  [W/m2]
!  (13) PBL      : (2-D) GMAO planetary boundary layer depth    [mb] 
!  (14) PREACC   : (2-D) GMAO accumulated precip @ ground       [mm H2O/day]
!  (15) PRECON   : (2-D) GMAO convective  precip @ ground       [mm H2O/day]
!  (16) PRECSNO  : (2-D) GMAO "snow" precip @ ground            
!  (17) RADLWG   : (2-D) GMAO upward LW flux @ ground           [W/m2]
!  (18) RADSWG   : (2-D) GMAO downward SW flux @ ground         [W/m2]
!  (19) RADSWT   : (2-D) GMAO downward SW flux @ atm top        [W/m2]
!  (20) SNICE    : (2-D) GCAP fraction of snow/ice              [unitless]
!  (21) SNODP    : (2-D) GMAO GEOS-5 geometric snow depth       [m]
!  (22) SNOMAS   : (2-D) GMAO GEOS-5 H2O equiv snow depth       [m]
!  (23) SNOW     : (2-D) GMAO snow depth (H2O equivalent)       [mm H2O]
!  (24) TROPP    : (2-D) GMAO tropopause pressure               [hPa]
!  (25) TS       : (2-D) GMAO surface air temperature           [K]
!  (26) TSKIN    : (2-D) GMAO radiance temperature              [K]
!  (27) USTAR    : (2-D) GMAO friction velocity                 [m/s]
!  (28) U10M     : (2-D) GMAO U-wind at 10 m                    [m/s]
!  (29) V10M     : (2-D) GMAO V-wind at 10 m                    [m/s]
!  (30) Z0       : (2-D) GMAO roughness height                  [m] 
!
!  NOTES:
!  (1 ) Now use function TIMESTAMP_STRING from "time_mod.f" for formatted 
!        date/time output. (bmy, 10/28/03)
!  (2 ) RADSWG, CLDFRC, USTAR, and Z0. are now 2-D arrays.  Also added RADLWG 
!        and SNOW arrays via the arg list.  Now skip over LAI. (bmy, 12/9/03)
!  (3 ) Now modified for GEOS-5 and GCAP met fields.  Added GCAP MOLENGTH,
!        SNICE, OICE optional arguments. (swu, bmy, 5/25/05)
!  (4 ) Fixed typo in the ND67 diagnostic for RADSWG (swu, bmy, 2/9/06)
!  (5 ) Remove support for GEOS-1 and GEOS-STRAT met fields (bmy, 8/4/06)
!  (6 ) Add "PARDIF", "PARDIR" to case statement for GCAP (swu, bmy, 10/4/06)
!  (7 ) Add EVAP, GRN, GWETROOT, LAI, PRECSNO, SNODP, SNOMAS, and TROPP as 
!        optional arguments.  Also update the CASE statement accordingly for 
!        GEOS-5 fields.  Convert GEOS-5 PRECTOT and PRECCON fields from
!        [kg/m2/s] to [mm/day] for backwards compatibility. (bmy, 1/17/07)
!  (8 ) Now get the # of A-3 fields from the file ident string (bmy, 10/7/08)
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
      REAL*8,  INTENT(OUT), OPTIONAL :: CLDFRC(IIPAR,JJPAR)
      REAL*8,  INTENT(OUT), OPTIONAL :: EVAP(IIPAR,JJPAR) 
      REAL*8,  INTENT(OUT), OPTIONAL :: GRN(IIPAR,JJPAR)
      REAL*8,  INTENT(OUT), OPTIONAL :: GWETROOT(IIPAR,JJPAR) 
      REAL*8,  INTENT(OUT), OPTIONAL :: GWETTOP(IIPAR,JJPAR) 
      REAL*8,  INTENT(OUT), OPTIONAL :: HFLUX(IIPAR,JJPAR)
      REAL*8,  INTENT(OUT), OPTIONAL :: LAI(IIPAR,JJPAR)
      REAL*8,  INTENT(OUT), OPTIONAL :: MOLENGTH(IIPAR,JJPAR)
      REAL*8,  INTENT(OUT), OPTIONAL :: OICE(IIPAR,JJPAR)
      REAL*8,  INTENT(OUT), OPTIONAL :: PARDF(IIPAR,JJPAR) 
      REAL*8,  INTENT(OUT), OPTIONAL :: PARDR(IIPAR,JJPAR) 
      REAL*8,  INTENT(OUT), OPTIONAL :: PBL(IIPAR,JJPAR) 
      REAL*8,  INTENT(OUT), OPTIONAL :: PREACC(IIPAR,JJPAR) 
      REAL*8,  INTENT(OUT), OPTIONAL :: PRECON(IIPAR,JJPAR)
      REAL*8,  INTENT(OUT), OPTIONAL :: PRECSNO(IIPAR,JJPAR)
      REAL*8,  INTENT(OUT), OPTIONAL :: RADLWG(IIPAR,JJPAR)
      REAL*8,  INTENT(OUT), OPTIONAL :: RADSWG(IIPAR,JJPAR)
      REAL*8,  INTENT(OUT), OPTIONAL :: RADSWT(IIPAR,JJPAR)
      REAL*8,  INTENT(OUT), OPTIONAL :: SNICE(IIPAR,JJPAR)
      REAL*8,  INTENT(OUT), OPTIONAL :: SNODP(IIPAR,JJPAR)
      REAL*8,  INTENT(OUT), OPTIONAL :: SNOMAS(IIPAR,JJPAR)
      REAL*8,  INTENT(OUT), OPTIONAL :: SNOW(IIPAR,JJPAR)
      REAL*8,  INTENT(OUT), OPTIONAL :: TROPP(IIPAR,JJPAR) 
      REAL*8,  INTENT(OUT), OPTIONAL :: TS(IIPAR,JJPAR) 
      REAL*8,  INTENT(OUT), OPTIONAL :: TSKIN(IIPAR,JJPAR) 
      REAL*8,  INTENT(OUT), OPTIONAL :: U10M(IIPAR,JJPAR) 
      REAL*8,  INTENT(OUT), OPTIONAL :: USTAR(IIPAR,JJPAR) 
      REAL*8,  INTENT(OUT), OPTIONAL :: V10M(IIPAR,JJPAR)
      REAL*8,  INTENT(OUT), OPTIONAL :: Z0(IIPAR,JJPAR)

      ! Local Variables
      INTEGER                        :: I, IJLOOP, IOS, J, N_A3, NFOUND 
      REAL*4                         :: Q2(IGLOB,JGLOB)
      CHARACTER(LEN=8)               :: NAME
      CHARACTER(LEN=16)              :: STAMP
      INTEGER                        :: XYMD, XHMS

      !=================================================================
      ! READ_A3 begins here!      
      !=================================================================

      ! Get the number of A-3 fields stored in this data set
#if   defined( GEOS_5 ) 
      N_A3 = N_A3_FIELDS
#else
      N_A3 = GET_N_A3( NYMD )
#endif

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
            ! (aka CLDTOT in GEOS-5)
            !---------------------------------
            CASE ( 'CLDFRC', 'CLDTOT' )
               READ ( IU_A3, IOSTAT=IOS ) XYMD, XHMS, Q2
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A3, 'read_a3:3' )
             
               IF ( CHECK_TIME( XYMD, XHMS, NYMD, NHMS ) ) THEN
                  IF ( PRESENT( CLDFRC ) ) CALL TRANSFER_2D( Q2,CLDFRC )
                  NFOUND = NFOUND + 1
               ENDIF

            !---------------------------------
            ! EVAP: evapotranspiration flux
            !---------------------------------
            CASE ( 'EVAP' )
               READ ( IU_A3, IOSTAT=IOS ) XYMD, XHMS, Q2
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A3, 'read_a3:4' )
             
               IF ( CHECK_TIME( XYMD, XHMS, NYMD, NHMS ) ) THEN
                  IF ( PRESENT( EVAP ) ) CALL TRANSFER_2D( Q2, EVAP )
                  NFOUND = NFOUND + 1
               ENDIF

            !---------------------------------
            ! GRN: evapotranspiration flux
            !---------------------------------
            CASE ( 'GRN' )
               READ ( IU_A3, IOSTAT=IOS ) XYMD, XHMS, Q2
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A3, 'read_a3:5' )
             
               IF ( CHECK_TIME( XYMD, XHMS, NYMD, NHMS ) ) THEN
                  IF ( PRESENT( GRN ) ) CALL TRANSFER_2D( Q2, GRN )
                  NFOUND = NFOUND + 1
               ENDIF

            !--------------------------------
            ! GWETROOT: root soil wetness
            !--------------------------------
            CASE ( 'GWETROOT' ) 
               READ( IU_A3, IOSTAT=IOS ) XYMD, XHMS, Q2
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A3, 'read_a3:6' )

               IF ( CHECK_TIME( XYMD, XHMS, NYMD, NHMS ) ) THEN
                  IF ( PRESENT( GWETROOT ) ) THEN
                     CALL TRANSFER_2D( Q2,GWETROOT )
                  ENDIF
                  NFOUND = NFOUND + 1
               ENDIF

            !--------------------------------
            ! GWETTOP: topsoil wetness 
            !--------------------------------
            CASE ( 'GWETTOP', 'SOIL' ) 
               READ( IU_A3, IOSTAT=IOS ) XYMD, XHMS, Q2
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A3, 'read_a3:7' )

               IF ( CHECK_TIME( XYMD, XHMS, NYMD, NHMS ) ) THEN
                  IF ( PRESENT( GWETTOP ) ) CALL TRANSFER_2D(Q2,GWETTOP)
                  NFOUND = NFOUND + 1
               ENDIF

            !--------------------------------
            ! HFLUX:  sensible heat flux 
            !--------------------------------
            CASE ( 'HFLUX' ) 
               READ( IU_A3, IOSTAT=IOS ) XYMD, XHMS, Q2
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A3, 'read_a3:8' )

               IF ( CHECK_TIME( XYMD, XHMS, NYMD, NHMS ) ) THEN
                  IF ( PRESENT( HFLUX ) ) CALL TRANSFER_2D( Q2, HFLUX )
                  NFOUND = NFOUND + 1
               ENDIF

            !--------------------------------
            ! LAI: GMAO leaf area index
            !--------------------------------
            CASE ( 'LAI' ) 
               READ( IU_A3, IOSTAT=IOS ) XYMD, XHMS, Q2
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A3, 'read_a3:9' )

               IF ( CHECK_TIME( XYMD, XHMS, NYMD, NHMS ) ) THEN
                  IF ( PRESENT( LAI ) ) CALL TRANSFER_2D( Q2, LAI )
                  NFOUND = NFOUND + 1
               ENDIF

            !--------------------------------
            ! MOLENGTH: GCAP M-O length
            !--------------------------------
            CASE ( 'MOLENGTH' ) 
               READ( IU_A3, IOSTAT=IOS ) XYMD, XHMS, Q2
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A3, 'read_a3:10' )

               IF ( CHECK_TIME( XYMD, XHMS, NYMD, NHMS ) ) THEN
                  IF ( PRESENT( MOLENGTH ) ) THEN
                     CALL TRANSFER_2D( Q2, MOLENGTH )
                  ENDIF
                  NFOUND = NFOUND + 1
               ENDIF

            !--------------------------------
            ! OICE: GCAP frac of ocean ice
            !--------------------------------
            CASE ( 'OICE' ) 
               READ( IU_A3, IOSTAT=IOS ) XYMD, XHMS, Q2
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A3, 'read_a3:11' )

               IF ( CHECK_TIME( XYMD, XHMS, NYMD, NHMS ) ) THEN
                  IF ( PRESENT( OICE ) ) CALL TRANSFER_2D( Q2, OICE )
                  NFOUND = NFOUND + 1
               ENDIF

            !--------------------------------
            ! PARDF: photosyn active diff rad
            !--------------------------------
            CASE ( 'PARDF', 'PARDIF' )
               READ( IU_A3, IOSTAT=IOS ) XYMD, XHMS, Q2
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A3, 'read_a3:12' )
             
               IF ( CHECK_TIME( XYMD, XHMS, NYMD, NHMS ) ) THEN
                  IF ( PRESENT( PARDF ) ) CALL TRANSFER_2D( Q2, PARDF )
                  NFOUND = NFOUND + 1
               ENDIF  

            !--------------------------------
            ! PARDR: photosyn active dir rad
            !--------------------------------
            CASE ( 'PARDR', 'PARDIR' )
               READ( IU_A3, IOSTAT=IOS ) XYMD, XHMS, Q2
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A3, 'read_a3:13' )
             
               IF ( CHECK_TIME( XYMD, XHMS, NYMD, NHMS ) ) THEN
                  IF ( PRESENT( PARDR ) ) CALL TRANSFER_2D( Q2, PARDR )
                  NFOUND = NFOUND + 1
               ENDIF 

            !--------------------------------
            ! PBL: boundary layer depth
            !--------------------------------
            CASE ( 'PBL', 'PBLH' )
               READ( IU_A3, IOSTAT=IOS ) XYMD, XHMS, Q2
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A3, 'read_a3:14' )
             
               IF ( CHECK_TIME( XYMD, XHMS, NYMD, NHMS ) ) THEN
                  IF ( PRESENT( PBL ) ) CALL TRANSFER_2D( Q2, PBL )
                  NFOUND = NFOUND + 1
               ENDIF       

            !--------------------------------
            ! PREACC: total precip at ground
            ! (aka PRECTOT in GEOS-5)
            !--------------------------------
            CASE ( 'PREACC', 'PRECTOT' )
               READ( IU_A3, IOSTAT=IOS ) XYMD, XHMS, Q2
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A3, 'read_a3:15' )
             
               IF ( CHECK_TIME( XYMD, XHMS, NYMD, NHMS ) ) THEN
                  IF ( PRESENT( PREACC ) ) CALL TRANSFER_2D( Q2,PREACC )
                  NFOUND = NFOUND + 1
               ENDIF

            !--------------------------------
            ! PRECON: conv precip at ground
            ! (aka PRECCON in GEOS-5)
            !--------------------------------
            CASE ( 'PRECON', 'PRECCON' )
               READ( IU_A3, IOSTAT=IOS ) XYMD, XHMS, Q2
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A3, 'read_a3:16' )
              
               IF ( CHECK_TIME( XYMD, XHMS, NYMD, NHMS ) ) THEN
                  IF ( PRESENT( PRECON ) ) CALL TRANSFER_2D( Q2,PRECON )
                  NFOUND = NFOUND + 1
               ENDIF

            !--------------------------------
            ! PRECSNO: snow precip at ground
            !--------------------------------
            CASE ( 'PRECSNO' )
               READ( IU_A3, IOSTAT=IOS ) XYMD, XHMS, Q2
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A3, 'read_a3:17' )
              
               IF ( CHECK_TIME( XYMD, XHMS, NYMD, NHMS ) ) THEN
                  IF ( PRESENT( PRECSNO ) ) CALL TRANSFER_2D(Q2,PRECSNO)
                  NFOUND = NFOUND + 1
               ENDIF

            !--------------------------------
            ! RADLWG: solar rad at ground 
            ! (aka LWGNET in GEOS-5)
            !--------------------------------
            CASE ( 'RADLWG', 'LWGNET' ) 
               READ( IU_A3, IOSTAT=IOS ) XYMD, XHMS, Q2
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A3, 'read_a3:18' )
             
               IF ( CHECK_TIME( XYMD, XHMS, NYMD, NHMS ) ) THEN
                  IF ( PRESENT( RADLWG ) ) CALL TRANSFER_2D( Q2,RADLWG )
                  NFOUND = NFOUND + 1
               ENDIF

            !--------------------------------
            ! RADSWG: solar rad at ground 
            ! (aka SWGNET in GEOS-5)
            !--------------------------------
            CASE ( 'RADSWG', 'SWGNET' ) 
               READ( IU_A3, IOSTAT=IOS ) XYMD, XHMS, Q2
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A3, 'read_a3:19' )
             
               IF ( CHECK_TIME( XYMD, XHMS, NYMD, NHMS ) ) THEN
                  IF ( PRESENT( RADSWG ) ) CALL TRANSFER_2D( Q2,RADSWG )
                  NFOUND = NFOUND + 1
               ENDIF

            !--------------------------------
            ! RADSWT: solar rad at atm top
            !--------------------------------
            CASE ( 'RADSWT' )
               READ( IU_A3, IOSTAT=IOS ) XYMD, XHMS, Q2
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A3, 'read_a3:20' )
             
               IF ( CHECK_TIME( XYMD, XHMS, NYMD, NHMS ) ) THEN
                  IF ( PRESENT( RADSWT ) ) CALL TRANSFER_2D( Q2,RADSWT )
                  NFOUND = NFOUND + 1
               ENDIF

            !--------------------------------
            ! SNICE: GCAP frac of snow/ice
            !--------------------------------
            CASE ( 'SNICE' ) 
               READ( IU_A3, IOSTAT=IOS ) XYMD, XHMS, Q2
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A3, 'read_a3:21' )

               IF ( CHECK_TIME( XYMD, XHMS, NYMD, NHMS ) ) THEN
                  IF ( PRESENT( SNICE ) ) CALL TRANSFER_2D( Q2, SNICE )
                  NFOUND = NFOUND + 1
               ENDIF

            !--------------------------------
            ! SNOMAS: GEOS-5 snow depth 
            ! (H2O equivalent)
            !--------------------------------
            CASE ( 'SNOMAS' )
               READ( IU_A3, IOSTAT=IOS ) XYMD, XHMS, Q2
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A3, 'read_a3:22' )
             
               IF ( CHECK_TIME( XYMD, XHMS, NYMD, NHMS ) ) THEN
                  IF ( PRESENT( SNOMAS ) ) CALL TRANSFER_2D( Q2,SNOMAS )
                  NFOUND = NFOUND + 1
               ENDIF

            !--------------------------------
            ! SNODP: GEOS-5 snow depth 
            ! (geometric & snowpack)
            !--------------------------------
            CASE ( 'SNODP' )
               READ( IU_A3, IOSTAT=IOS ) XYMD, XHMS, Q2
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A3, 'read_a3:23' )
             
               IF ( CHECK_TIME( XYMD, XHMS, NYMD, NHMS ) ) THEN
                  IF ( PRESENT( SNODP ) ) CALL TRANSFER_2D( Q2, SNODP )
                  NFOUND = NFOUND + 1
               ENDIF

            !--------------------------------
            ! SNOW: snow depth (H2O equiv.)
            !--------------------------------
            CASE ( 'SNOW', 'SNOWD' )
               READ( IU_A3, IOSTAT=IOS ) XYMD, XHMS, Q2
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A3, 'read_a3:24' )
             
               IF ( CHECK_TIME( XYMD, XHMS, NYMD, NHMS ) ) THEN
                  IF ( PRESENT( SNOW ) ) CALL TRANSFER_2D( Q2, SNOW )
                  NFOUND = NFOUND + 1
               ENDIF

            !--------------------------------            
            ! TROPP: tropopause pressure
            !--------------------------------
            CASE ( 'TROPP' )
               READ( IU_A3, IOSTAT=IOS ) XYMD, XHMS, Q2
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A3, 'read_a3:25' )
             
               IF ( CHECK_TIME( XYMD, XHMS, NYMD, NHMS ) ) THEN
                  IF ( PRESENT( TROPP ) ) CALL TRANSFER_2D( Q2, TROPP )
                  NFOUND = NFOUND + 1
               ENDIF

            !--------------------------------            
            ! TS: surface air temperature
            !--------------------------------
            CASE ( 'TS', 'TGROUND', 'T2M' )
               READ( IU_A3, IOSTAT=IOS ) XYMD, XHMS, Q2
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A3, 'read_a3:26' )
             
               IF ( CHECK_TIME( XYMD, XHMS, NYMD, NHMS ) ) THEN
                  IF ( PRESENT( TS ) ) CALL TRANSFER_2D( Q2, TS )
                  NFOUND = NFOUND + 1
               ENDIF

            !--------------------------------            
            ! TSKIN: radiance temperature
            !--------------------------------
            CASE ( 'TSKIN' )
               READ( IU_A3, IOSTAT=IOS ) XYMD, XHMS, Q2
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A3, 'read_a3:27' )
             
               IF ( CHECK_TIME( XYMD, XHMS, NYMD, NHMS ) ) THEN
                  IF ( PRESENT( TSKIN ) ) CALL TRANSFER_2D( Q2, TSKIN )
                  NFOUND = NFOUND + 1
               ENDIF

            !--------------------------------            
            ! U10M: U-wind at 10 m
            !--------------------------------            
            CASE ( 'U10M', 'USS' )
               READ( IU_A3, IOSTAT=IOS ) XYMD, XHMS, Q2
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A3, 'read_a3:28' )
             
               IF ( CHECK_TIME( XYMD, XHMS, NYMD, NHMS ) ) THEN
                  IF ( PRESENT( U10M ) ) CALL TRANSFER_2D( Q2, U10M )
                  NFOUND = NFOUND + 1
               ENDIF

            !--------------------------------            
            ! USTAR: friction velocity
            !--------------------------------            
            CASE ( 'USTAR' )
               READ( IU_A3, IOSTAT=IOS ) XYMD, XHMS, Q2
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A3, 'read_a3:29' )
             
               IF ( CHECK_TIME( XYMD, XHMS, NYMD, NHMS ) ) THEN
                  IF ( PRESENT( USTAR ) ) CALL TRANSFER_2D( Q2, USTAR )
                  NFOUND = NFOUND + 1
               ENDIF
            
            !--------------------------------            
            ! V10M: V-wind at 10 m
            !--------------------------------            
            CASE ( 'V10M', 'VSS' )
               READ( IU_A3, IOSTAT=IOS ) XYMD, XHMS, Q2
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A3, 'read_a3:30' )
             
               IF ( CHECK_TIME( XYMD, XHMS, NYMD, NHMS ) ) THEN
                  IF ( PRESENT( V10M ) ) CALL TRANSFER_2D( Q2, V10M )
                  NFOUND = NFOUND + 1
               ENDIF

            !--------------------------------            
            ! Z0: roughness heights
            !--------------------------------            
            CASE ( 'Z0', 'Z0M' )
               READ( IU_A3, IOSTAT=IOS ) XYMD, XHMS, Q2
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A3, 'read_a3:31' )
             
               IF ( CHECK_TIME( XYMD, XHMS, NYMD, NHMS ) ) THEN
                  IF ( PRESENT( Z0 ) ) CALL TRANSFER_2D( Q2, Z0 ) 
                  NFOUND = NFOUND + 1
               ENDIF
  
            !--------------------------------
            ! TPW: just skip over this
            !--------------------------------
            CASE ( 'TPW' )
               READ( IU_A3, IOSTAT=IOS ) XYMD, XHMS, Q2
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A3, 'read_a3:32' )

               IF ( CHECK_TIME( XYMD, XHMS, NYMD, NHMS ) ) THEN
                  NFOUND = NFOUND + 1
               ENDIF

            !--------------------------------
            ! CLDTMP: just skip over this
            !--------------------------------
            CASE ( 'CLDTMP' )
               READ( IU_A3, IOSTAT=IOS ) XYMD, XHMS, Q2
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A3, 'read_a3:33' )
             
               IF ( CHECK_TIME( XYMD, XHMS, NYMD, NHMS ) ) THEN
                  NFOUND = NFOUND + 1
               ENDIF

            !--------------------------------
            ! Extra GEOS-5 fields
            ! Skip over for now, add later
            ! (bmy, 10/7/08)
            !--------------------------------
            CASE ( 'EFLUX',   'FRLAKE',  'FRLAND', 'FRLANDIC',
     &             'FROCEAN', 'PRECANV', 'LWTUP'  )
               READ( IU_A3, IOSTAT=IOS ) XYMD, XHMS, Q2
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A3, 'read_a3:33' )
             
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
            STAMP = TIMESTAMP_STRING( NYMD, NHMS )
            WRITE( 6, 210 ) NFOUND, STAMP
 210        FORMAT( '     - Found all ', i3, ' A-3 met fields for ', a )
            EXIT
         ENDIF
      ENDDO

#if   defined( GEOS_5 )

      !=================================================================
      ! SPECIAL HANDLING FOR GEOS-5
      ! ---------------------------
      ! In GEOS-5, the PRECTOT (cf PREACC) and PRECCON (cf PRECON) met
      ! fields have units of [kg/m2/s].  In all other GEOS versions, 
      ! PREACC and PRECON have units of [mm/day].  
      !
      ! Therefore, for backwards compatibility with existing code, 
      ! apply the following unit conversion to the GEOS-5 PRECTOT and
      ! PRECCON fields:
      !
      !
      !     kg  |    m3    | 86400 s | 1000 mm
      !   ------+----------+---------+--------- = 86400 
      !    m2 s |  1000 kg |  day    |   m
      !              ^
      !              |
      !       1 / density of water 
      !=================================================================
      
      ! Convert GEOS-5 PRECTOT from [kg/m2/s] to [mm/day]
      IF ( PRESENT( PREACC ) ) PREACC = PREACC * 86400d0

      ! Convert GEOS-5 PRECCON from [kg/m2/s] to [mm/day]
      IF ( PRESENT( PRECON ) ) PRECON = PRECON * 86400d0

#endif

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
         IF ( PRESENT( HFLUX   ) ) AD67(:,:,1 ) = AD67(:,:,1 ) + HFLUX
         IF ( PRESENT( RADSWG  ) ) AD67(:,:,2 ) = AD67(:,:,2 ) + RADSWG
         IF ( PRESENT( PREACC  ) ) AD67(:,:,3 ) = AD67(:,:,3 ) + PREACC
         IF ( PRESENT( PRECON  ) ) AD67(:,:,4 ) = AD67(:,:,4 ) + PRECON
         IF ( PRESENT( TS      ) ) AD67(:,:,5 ) = AD67(:,:,5 ) + TS
         IF ( PRESENT( RADSWT  ) ) AD67(:,:,6 ) = AD67(:,:,6 ) + RADSWT
         IF ( PRESENT( USTAR   ) ) AD67(:,:,7 ) = AD67(:,:,7 ) + USTAR
         IF ( PRESENT( Z0      ) ) AD67(:,:,8 ) = AD67(:,:,8 ) + Z0
         IF ( PRESENT( PBL     ) ) AD67(:,:,9 ) = AD67(:,:,9 ) + PBL
         IF ( PRESENT( CLDFRC  ) ) AD67(:,:,10) = AD67(:,:,10) + CLDFRC
         IF ( PRESENT( U10M    ) ) AD67(:,:,11) = AD67(:,:,11) + U10M
         IF ( PRESENT( V10M    ) ) AD67(:,:,12) = AD67(:,:,12) + V10M
         IF ( PRESENT( ALBEDO  ) ) AD67(:,:,14) = AD67(:,:,14) + ALBEDO
         IF ( PRESENT( TROPP   ) ) AD67(:,:,17) = AD67(:,:,17) + TROPP
         IF ( PRESENT( TSKIN   ) ) AD67(:,:,19) = AD67(:,:,19) + TSKIN
         IF ( PRESENT( PARDF   ) ) AD67(:,:,20) = AD67(:,:,20) + PARDF
         IF ( PRESENT( PARDR   ) ) AD67(:,:,21) = AD67(:,:,21) + PARDR
         IF ( PRESENT( GWETTOP ) ) AD67(:,:,22) = AD67(:,:,22) + GWETTOP
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

      ! End of module
      END MODULE A3_READ_MOD
