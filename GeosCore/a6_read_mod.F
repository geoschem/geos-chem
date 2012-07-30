! $Id: a6_read_mod.f,v 1.3 2010/02/02 16:57:55 bmy Exp $
      MODULE A6_READ_MOD
!
!******************************************************************************
!  Module A6_READ_MOD contains subroutines that unzip, open, and read
!  GEOS-CHEM A-6 (avg 6-hour) met fields from disk. (bmy, 6/19/03, 10/15/09)
! 
!  Module Routines:
!  ============================================================================
!  (1 ) UNZIP_A6_FIELDS : Unzips & copies met field files to a temp dir
!  (2 ) DO_OPEN_A6      : Returns TRUE if it's time to open A-6 fields
!  (3 ) OPEN_A6_FIELDS  : Opens met field files residing in the temp dir
!  (4 ) GET_A6_FIELDS   : Wrapper for routine READ_A6
!  (5 ) MAKE_GCAP_CLDFRC: Computes CLDFRC from 3-D CLDF field for GCAP
!  (6 ) GET_N_A6        : Returns # of A-6 fields for each DAO data set
!  (7 ) CHECK_TIME      : Tests if A-6 et field timestamps equal current time
!  (8 ) READ_A6         : Reads A-6 fields from disk
!  (9 ) A6_CHECK        : Checks if we have found all of the A-6 fields
! 
!  GEOS-CHEM modules referenced by a6_read_mod.f
!  ============================================================================
!  (1 ) bpch2_mod.f     : Module w/ routines for binary punch file I/O
!  (2 ) dao_mod.f       : Module w/ arrays for DAO met fields
!  (3 ) diag_mod.f      : Module w/ GEOS-CHEM diagnostic arrays
!  (4 ) directory_mod.f : Module w/ GEOS-CHEM data & met field dirs
!  (5 ) error_mod.f     : Module w/ NaN and other error check routines
!  (6 ) logical_mod.f   : Module w/ GEOS-CHEM logical switches 
!  (7 ) file_mod.f      : Module w/ file unit #'s and error checks
!  (8 ) pressure_mod.f  : Module w/ routines to compute P(I,J,L)
!  (9 ) time_mod.f      : Module w/ routines for computing time & date
!  (10) transfer_mod.f  : Module w/ routines to cast & resize arrays
!  (11) unix_cmds_mod.f : Module w/ Unix commands for unzipping etc.
!
!  NOTES:
!  (1 ) Adapted from "dao_read_mod.f" (bmy, 6/19/03)
!  (2 ) Now use TIMESTAMP_STRING for formatted output (bmy, 10/28/03)
!  (3 ) CLDFRC is now a 2-D array in MAKE_CLDFRC< GET_A6_FIELDS.  Also now
!        read from either zipped or unzipped files. (bmy, 12/9/03)
!  (4 ) Now skips past the GEOS-4 ident string (bmy, 12/12/03)
!  (5 ) Bug fix: need to determine CLDTOPS for GEOS-4.  (bmy, 3/4/04)
!  (6 ) Now modified for GEOS-4 "a_llk_03" and "a_llk_04" data (bmy, 3/4/04)
!  (7 ) Now references "unix_cmds_mod.f", "directory_mod.f" and
!        "logical_mod.f" (bmy, 7/20/04)
!  (8 ) Now references FILE_EXISTS from "file_mod.f" (bmy, 3/23/05)
!  (9 ) Now modified for GEOS-5 and GCAP met fields.  Added MAKE_GCAP_CLDFRC
!        routine. (swu, bmy, 5/25/05)
!  (10) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (11) Bug fix in ND66 diagnostic for ZMMU (bmy, 2/1/06)
!  (12) Remove support for GEOS-1 and GEOS-STRAT met fields (bmy, 8/4/06)
!  (13) Now set negative Q (i.e. SPHU)to a small positive # (bmy, 9/8/06)
!  (14) Now read extra fields for GEOS-5.  Bug fix: we must convert RH from 
!        unitless to % to be compatible w/ present drydep etc. algorithms. 
!        (phs, bmy, 3/28/08)
!  (15) Now get the # of A-6 fields from the file ident string (bmy, 10/7/08)
!  (16) Remove references to IN_CLOUD_OD (bmy, 10/15/09)
!  28 Feb 2012 - R. Yantosca - Removed support for GEOS-3
!******************************************************************************
!
#if defined( DEVEL )
      USE inquireMod, ONLY : findFreeLUN
#endif
      IMPLICIT NONE
#     include "define.h"

      !=================================================================
      ! MODULE PRIVATE DECLARATIONS -- keep certain internal variables 
      ! and routines from being seen outside "a6_read_mod.f"
      !=================================================================

      ! Make everything PRIVATE ...
      PRIVATE

      ! ... except these routines
      PUBLIC :: GET_A6_FIELDS   
      PUBLIC :: OPEN_A6_FIELDS  
      PUBLIC :: UNZIP_A6_FIELDS

      !=================================================================
      ! MODULE VARIABLES
      !=================================================================

      ! Number of A6 fields in the file
      INTEGER :: N_A6_FIELDS

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement 
      !=================================================================
      CONTAINS

!------------------------------------------------------------------------------

      SUBROUTINE UNZIP_A6_FIELDS( OPTION, NYMD )
!
!******************************************************************************
!  Subroutine UNZIP_A6_FIELDS invokes a FORTRAN system call to uncompress
!  GEOS-CHEM A-6 met field files and store the uncompressed data in a 
!  temporary directory, where GEOS-CHEM can read them.  The original data 
!  files are not disturbed.  (bmy, bdf, 6/15/98, 8/4/06)
!
!  Arguments as input:
!  ============================================================================
!  (1 ) OPTION (CHAR*(*)) : Option
!  (2 ) NYMD   (INTEGER ) : YYYYMMDD of A-6 file to be unzipped (optional)
!
!  NOTES:
!  (1 ) Adapted from UNZIP_MET_FIELDS of "dao_read_mod.f" (bmy, 6/19/03)
!  (2 ) Directory information YYYY/MM or YYYYMM is now contained w/in 
!        GEOS_1_DIR, GEOS_S_DIR, GEOS_3_DIR, GEOS_4_DIR (bmy, 12/11/03)
!  (3 ) Now reference "directory_mod.f" and "unix_cmds_mod.f". Now prevent 
!        EXPAND_DATE from overwriting directory paths with Y/M/D tokens in 
!        them (bmy, 7/20/04)
!  (4 ) Removed code for GEOS-4 a_llk_03 data.  Also modified for GEOS-5
!        and GCAP met fields. (bmy, 5/25/05)
!  (5 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (6 ) Remove support for GEOS-1 and GEOS-STRAT met fields (bmy, 8/4/06)
!  28 Feb 2012 - R. Yantosca - Removed support for GEOS-3
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD,     ONLY : GET_RES_EXT
      USE DIRECTORY_MOD, ONLY : DATA_DIR,   GCAP_DIR
      USE DIRECTORY_MOD, ONLY : GEOS_4_DIR, GEOS_5_DIR, TEMP_DIR 
      USE ERROR_MOD,     ONLY : ERROR_STOP
      USE TIME_MOD,      ONLY : EXPAND_DATE
      USE UNIX_CMDS_MOD, ONLY : BACKGROUND, REDIRECT,   REMOVE_CMD 
      USE UNIX_CMDS_MOD, ONLY : UNZIP_CMD,  WILD_CARD,  ZIP_SUFFIX

      USE CMN_SIZE_MOD

      ! Arguments
      CHARACTER(LEN=*),  INTENT(IN) :: OPTION
      INTEGER, OPTIONAL, INTENT(IN) :: NYMD

      ! Local variables
      CHARACTER(LEN=255)            :: GEOS_DIR,   A6_STR
      CHARACTER(LEN=255)            :: A6_FILE_GZ, A6_FILE
      CHARACTER(LEN=255)            :: UNZIP_BG,   UNZIP_FG
      CHARACTER(LEN=255)            :: REMOVE_ALL, REMOVE_DATE

      !=================================================================
      ! UNZIP_A6_FIELDS begins here!
      !=================================================================
      IF ( PRESENT( NYMD ) ) THEN

#if   defined( GEOS_4 )

         ! Strings for directory & filename
         GEOS_DIR = TRIM( GEOS_4_DIR )
         A6_STR   = 'YYYYMMDD.a6.' // GET_RES_EXT() 

#elif defined( GEOS_5 )

         ! Strings for directory & filename
         GEOS_DIR = TRIM( GEOS_5_DIR )
         A6_STR   = 'YYYYMMDD.a6.' // GET_RES_EXT() 

#elif defined( GCAP )

         ! Strings for directory & filename
         GEOS_DIR = TRIM( GCAP_DIR )
         A6_STR   = 'YYYYMMDD.a6.' // GET_RES_EXT() 

#endif

         ! Replace date tokens
         CALL EXPAND_DATE( GEOS_DIR, NYMD, 000000 )
         CALL EXPAND_DATE( A6_STR,   NYMD, 000000 )

         ! Location of zipped A-3 file in data dir
         A6_FILE_GZ  = TRIM( DATA_DIR  ) // TRIM( GEOS_DIR   ) // 
     &                 TRIM( A6_STR    ) // TRIM( ZIP_SUFFIX )

         ! Location of unzipped A-3 file in temp dir
         A6_FILE     = TRIM( TEMP_DIR  ) // TRIM( A6_STR     )
         
         ! Remove A-3 files for this date from temp dir 
         REMOVE_DATE = TRIM( REMOVE_CMD ) // ' '               // 
     &                 TRIM( TEMP_DIR   ) // TRIM( A6_STR    ) 

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
!  been opened. (bmy, 6/19/03, 5/25/05)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) NYMD (INTEGER) : YYYYMMDD 
!  (2 ) NHMS (INTEGER) :  and HHMMSS to be tested for A-3 file open
!
!  NOTES:
!  (1 ) Now modified for GEOS-4 "a_llk_03" or "a_llk_04" data (bmy, 3/22/04)
!  (2 ) Remove code for obsolete GEOS-4 a_llk_03 data.  Also modified for
!        GEOS-5 and GCAP met fields. (swu, bmy, 5/25/05)
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

#if   defined( GCAP )

      ! Open file if it's 03 GMT or first call (GCAP only) 
      IF ( NHMS == 030000 .or. FIRST ) THEN
         DO_OPEN = .TRUE. 
         GOTO 999
      ENDIF

#else

      ! Open file if it's 00:00 GMT or first call (all GEOS data)
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
!  time NHMS. (bmy, bdf, 6/15/98, 10/15/09)
!  
!  Arguments as input:
!  ===========================================================================
!  (1 ) NYMD (INTEGER)   : Current value of YYYYMMDD
!  (2 ) NHMS (INTEGER)   : Current value of HHMMSS
!
!  NOTES:
!  (1 ) Adapted from OPEN_MET_FIELDS of "dao_read_mod.f" (bmy, 6/19/03)
!  (2 ) Now opens either zipped or unzipped files (bmy, 12/11/03)
!  (3 ) Now skips past the GEOS-4 ident string (bmy, 12/12/03)
!  (4 ) Now references "directory_mod.f" instead of CMN_SETUP.  Also now
!        references LUNZIP from "logical_mod.f".  Also now prevents EXPAND_DATE
!        from overwriting Y/M/D tokens in directory paths. (bmy, 7/20/04)
!  (5 ) Now use FILE_EXISTS from "file_mod.f" to determine if file unit IU_A6 
!        refers to a valid file on disk (bmy, 3/23/05)
!  (6 ) Now modified for GEOS-5 and GCAP met fields (swu, bmy, 5/25/05)
!  (7 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (8 ) Remove support for GEOS-1 and GEOS-STRAT met fields (bmy, 8/4/06)
!  (9 ) Now get the # of A-3 fields from the file ident string (bmy, 10/7/08)
!  (10) Set N_A6_FIELDS=21 for GEOS-5 and IN_CLOUD_OD (jmao, bmy, 2/12/09)
!  (11) Remove references to IN_CLOUD_OD (bmy, 10/15/09)
!  28 Feb 2012 - R. Yantosca - Removed support for GEOS-3
!******************************************************************************
!      
      ! References to F90 modules
      USE BPCH2_MOD,     ONLY : GET_RES_EXT
      USE DIRECTORY_MOD, ONLY : DATA_DIR,   GCAP_DIR
      USE DIRECTORY_MOD, ONLY : GEOS_4_DIR, GEOS_5_DIR, TEMP_DIR 
      USE ERROR_MOD,     ONLY : ERROR_STOP
      USE LOGICAL_MOD,   ONLY : LUNZIP
      USE FILE_MOD,      ONLY : IU_A6, IOERROR, FILE_EXISTS
      USE TIME_MOD,      ONLY : EXPAND_DATE

      USE CMN_SIZE_MOD        ! Size parameters

      ! Arguments
      INTEGER, INTENT(IN)    :: NYMD, NHMS

      ! Local variables
      LOGICAL, SAVE          :: FIRST = .TRUE.
      LOGICAL                :: IT_EXISTS, IS_OPEN
      INTEGER                :: IOS, IUNIT
      CHARACTER(LEN=8)       :: IDENT
      CHARACTER(LEN=255)     :: A6_FILE
      CHARACTER(LEN=255)     :: GEOS_DIR
      CHARACTER(LEN=255)     :: PATH

      !=================================================================
      ! OPEN_A6_FIELDS begins here!
      !=================================================================

      ! Open A-6 file at the proper time, or on the first call
      IF ( DO_OPEN_A6( NYMD, NHMS ) ) THEN

#if   defined( GEOS_4 )

         ! Strings for directory & filename
         GEOS_DIR = TRIM( GEOS_4_DIR )
         A6_FILE  = 'YYYYMMDD.a6.' // GET_RES_EXT()

#elif defined( GEOS_5 )

         ! Strings for directory & filename
         GEOS_DIR = TRIM( GEOS_5_DIR )
         A6_FILE  = 'YYYYMMDD.a6.' // GET_RES_EXT()

#elif defined( GCAP )

         ! Strings for directory & filename
         GEOS_DIR = TRIM( GCAP_DIR )
         A6_FILE  = 'YYYYMMDD.a6.' // GET_RES_EXT()

#endif

         ! Replace date tokens
         CALL EXPAND_DATE( GEOS_DIR, NYMD, NHMS )
         CALL EXPAND_DATE( A6_FILE,  NYMD, NHMS )

         ! If unzipping, open GEOS-1 file in TEMP dir
         ! If not unzipping, open GEOS-1 file in DATA dir
         IF ( LUNZIP ) THEN
            PATH = TRIM( TEMP_DIR ) // TRIM( A6_FILE )
         ELSE
            PATH = TRIM( DATA_DIR ) // 
     &             TRIM( GEOS_DIR ) // TRIM( A6_FILE )
         ENDIF

         ! Close previously opened A-3 file
	 INQUIRE(IU_A6, OPENED=IS_OPEN)
         IF(IS_OPEN) CLOSE( IU_A6 )

         ! Make sure the file unit is valid before we open the file
         IF ( .not. FILE_EXISTS( IU_A6 ) ) THEN
            CALL ERROR_STOP( 'Could not find file!', 
     &                       'OPEN_A6_FIELDS (a6_read_mod.f)' )
         ENDIF

         ! Open the file
#if defined( DEVEL )
         IU_A6 = findFreeLUN()
#endif
         OPEN( UNIT   = IU_A6,         FILE   = TRIM( PATH ),
     &         STATUS = 'OLD',         ACCESS = 'SEQUENTIAL',  
     &         FORM   = 'UNFORMATTED', IOSTAT = IOS )
               
         IF ( IOS /= 0 ) THEN
            CALL IOERROR( IOS, IU_A6, 'open_a6_fields:1' )
         ENDIF

         ! Echo info
         WRITE( 6, 100 ) TRIM( PATH )
 100     FORMAT( '     - Opening: ', a ) 

         ! Skip past the ident string
         READ( IU_A6, IOSTAT=IOS ) IDENT

         IF ( IOS /= 0 ) THEN
            CALL IOERROR( IOS, IU_A6, 'open_a6_fields:2' )
         ENDIF

         ! The last 2 digits of the ident string
         ! is the # of fields contained in the file
         READ( IDENT(7:8), '(i2.2)' ) N_A6_FIELDS

#if   defined( GEOS_5 )
         !%%% KLUDGE: set N_A6_FIELDS=21 when using the reprocessed
         !%%% GEOS-5 met.   This accounts for CMFMC (which doesn't seem
         !%%% to get counted) as well as for MOISTQ, which is an extra
         !%%% derived field.  (jmao, bmy, 2/12/09)
         N_A6_FIELDS = 21
#endif

      ENDIF

      ! Return to calling program
      END SUBROUTINE OPEN_A6_FIELDS

!------------------------------------------------------------------------------

#if defined( DEVEL )
      SUBROUTINE GET_A6_FIELDS( NYMD, NHMS, LOCAL_MET )
#else
      SUBROUTINE GET_A6_FIELDS( NYMD, NHMS )
#endif
!
!******************************************************************************
!  Subroutine GET_A6_FIELDS is a wrapper for routine READ_A6.  GET_A6_FIELDS
!  calls READ_A6 properly for reading A-6 fields from GEOS-1, GEOS-STRAT, 
!  GEOS-3, GEO b  S-4, GEOS-5, or GCAP met data sets. (bmy, 6/19/03, 10/30/07)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) NYMD (INTEGER) : YYYYMMDD
!  (2 ) NHMS (INTEGER) :  and HHMMSS of A-6 fields to be read from disk
!
!  NOTES:
!  (1 ) CFRAC has been removed from CMN_DEP.  Now use CLDFRC(I,J) from
!        "dao_mod.f" (bmy, 12/9/03)
!  (2 ) Now pass CLDTOPS to READ_A6 for GEOS-4 (bmy, 3/4/04)
!  (3 ) Now modified for GEOS-5 and GCAP met fields (swu, bmy, 5/25/05)
!  (4 ) Remove support for GEOS-1 and GEOS-STRAT met fields (bmy, 8/4/06)
!  (5 ) Now read CMFMC, DQIDTMST, DQLDTMST, DQRCON, DQRLSC, DQVDTMST, MFXC,
!        MFYC, MFZ, PV, QI, QL, RH, TAUCLI, TAUCLW for GEOS-5 
!        (bmy, 10/30/07)
!  28 Feb 2012 - R. Yantosca - Removed support for GEOS-3
!******************************************************************************
!
      ! References to F90 modules
      USE DAO_MOD,      ONLY : CLDF,    CLDFRC,   CLDMAS,   CLDTOPS 
      USE DAO_MOD,      ONLY : CMFMC,   DETRAINE, DETRAINN, DNDE    
      USE DAO_MOD,      ONLY : DNDN,    DQIDTMST, DQLDTMST, DQRCON
      USE DAO_MOD,      ONLY : DQRLSC,  DQVDTMST, DTRAIN,   ENTRAIN
      USE DAO_MOD,      ONLY : HKBETA,  HKETA,    MFXC,     MFYC
      USE DAO_MOD,      ONLY : MFZ,     MOISTQ,   OPTDEP,   PV
      USE DAO_MOD,      ONLY : QI,      QL,       RH,       SPHU 
      USE DAO_MOD,      ONLY : T,       TAUCLI,   TAUCLW,   UPDE
      USE DAO_MOD,      ONLY : UPDN,    UWND,     VWND,     ZMEU
      USE DAO_MOD,      ONLY : ZMMD,    ZMMU

      USE CMN_SIZE_MOD  ! Size parameters
#if defined( DEVEL )
      USE GC_TYPE_MOD,  ONLY : GC_MET_LOCAL
#endif

      ! Arguments
      INTEGER, INTENT(IN) :: NYMD, NHMS 
#if defined( DEVEL )
      TYPE(GC_MET_LOCAL), INTENT(INOUT) :: LOCAL_MET
#endif

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

#if   defined( GEOS_4 )

      !=================================================================      
      ! GEOS-4: get CLDF, CLDTOPS, HKBETA, HKETA, MOISTQ, OPTDEP, SPHU
      !             TMPU, UWND,    VWND,  ZMEU,   ZMMD,   ZMMU
      !=================================================================
      CALL READ_A6( NYMD=NYMD,     NHMS=NHMS,       CLDTOPS=CLDTOPS,
     &              CLDF=CLDF,     HKBETA=HKBETA,   HKETA=HKETA,   
     &              MOISTQ=MOISTQ, OPTDEPTH=OPTDEP, Q=SPHU,        
     &              T=T,           U=UWND,          V=VWND,        
     &              ZMEU=ZMEU,     ZMMD=ZMMD,       ZMMU=ZMMU ) 

#elif defined( GEOS_5 )

      !=================================================================      
      ! GEOS-5: get CLDF,   CLDTOPS, CMFMC,    DQIDTMST, DQLDTMST,
      !             DQRCON, DQRLSC,  DQVDTMST, MFXC,     MFYC,
      !             MFZ,    MOISTQ,  OPTDEPTH, PLE,      PV,
      !             RH,     QV,      T,        TAUCLI,   TAUCLW,
      !             U,      V    fields
      !=================================================================
      CALL READ_A6( NYMD=NYMD,          NHMS=NHMS,       
     &              CLDF=CLDF,          CLDTOPS=CLDTOPS, 
     &              CMFMC=CMFMC,        DQIDTMST=DQIDTMST, 
     &              DQLDTMST=DQLDTMST,  DQRCON=DQRCON, 
     &              DQRLSC=DQRLSC,      DQVDTMST=DQVDTMST,
     &              DTRAIN=DTRAIN,     
     &              MOISTQ=MOISTQ,      OPTDEPTH=OPTDEP, 
     &              PV=PV,              RH=RH,             
     &              Q=SPHU,             QL=QL,             
     &              QI=QI,              T=T,               
     &              TAUCLI=TAUCLI,      TAUCLW=TAUCLW,     
     &              U=UWND,             V=VWND         ) 


#elif defined( GCAP ) 

      !=================================================================
      ! GCAP: read CLDF,   DETRAINE, DETRAIN, DNDE, DNDN, ENTRAIN,   
      !            MOISTQ, OPTDEPTH, SPHU,    T=T,  UWND, UPDE,
      !            UPDN,   VWND, and compute CLDTOPS & CLDFRC
      !=================================================================
      CALL READ_A6( NYMD=NYMD,          NHMS=NHMS,         
     &              CLDF=CLDF,          CLDTOPS=CLDTOPS, 
     &              DETRAINE=DETRAINE,  DETRAINN=DETRAINN, 
     &              DNDE=DNDE,          DNDN=DNDN,         
     &              ENTRAIN=ENTRAIN,    MOISTQ=MOISTQ,
     &              OPTDEPTH=OPTDEP,    Q=SPHU,            
     &              T=T,                U=UWND,            
     &              UPDE=UPDE,          UPDN=UPDN,  
     &              V=VWND )
          
      ! Create 2-D CLDFRC field from 3-D CLDF field
      CALL MAKE_GCAP_CLDFRC( CLDF, CLDFRC )
      
#endif

      ! Save NYMD and NHMS for next call
      LASTNYMD = NYMD
      LASTNHMS = NHMS

#if defined( DEVEL )
      LOCAL_MET%CLDF     = CLDF
      LOCAL_MET%CMFMC    = CMFMC
      LOCAL_MET%DTRAIN   = DTRAIN
      LOCAL_MET%MOISTQ   = MOISTQ
      LOCAL_MET%OPTD     = OPTDEP
      LOCAL_MET%RH       = RH
      LOCAL_MET%SPHU     = SPHU
      LOCAL_MET%T        = T
      LOCAL_MET%DQIDTMST = DQIDTMST
      LOCAL_MET%DQLDTMST = DQLDTMST
      LOCAL_MET%TAUCLI   = TAUCLI
      LOCAL_MET%TAUCLW   = TAUCLW
      LOCAL_MET%DQVDTMST = DQVDTMST
#endif

      ! Return to calling program
      END SUBROUTINE GET_A6_FIELDS

!------------------------------------------------------------------------------

      SUBROUTINE MAKE_GCAP_CLDFRC( CLDF, CLDFRC )
!
!******************************************************************************
!  Subroutine MAKE_CLDFRC constructs the GCAP CLDFRC field from the 3-D
!  cloud fraction field. (swu, bmy, 5/25/05)
!
!  Arguments as Input:
!  ===========================================================================
!  (1 ) CLDF   (REAL*8) : GCAP 3-D cloud fraction field [unitless]
!
!  Arguments as Output:
!  ===========================================================================
!  (2 ) CLDFRC (REAL*8) : GCAP column cloud fraction field [unitless]
!
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE DIAG_MOD, ONLY : AD67

      USE CMN_SIZE_MOD   ! Size parameters
      USE CMN_DIAG_MOD   ! ND67

      ! Arguments
      REAL*8, INTENT(IN)  :: CLDF(LLPAR,IIPAR,JJPAR)
      REAL*8, INTENT(OUT) :: CLDFRC(IIPAR,JJPAR)

      ! Local variables
      LOGICAL             :: IS_ND67 
      INTEGER             :: I, J
      
      !=================================================================
      ! MAKE_GCAP_CLDFRC begins here!
      !=================================================================

      ! Is the ND67 diagnostic turned on?
      IS_ND67 = ( ND67 > 0 )

      ! Make 2-D cloud fraction
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J )
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Take max value
         CLDFRC(I,J) = MAXVAL( CLDF(:,I,J) )

         ! Store in ND67 diagnostic if necessary
         IF ( IS_ND67 ) AD67(I,J,10) = AD67(I,J,10) + CLDFRC(I,J)
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Return to calling program
      END SUBROUTINE MAKE_GCAP_CLDFRC

!------------------------------------------------------------------------------

      FUNCTION GET_N_A6() RESULT( N_A6 )
!
!******************************************************************************
!  Function GET_N_A6 returns the number of A-6 fields per met data set
!  (GEOS-3, GEOS-4, GEOS-5, or GCAP). (bmy, 6/19/03, 5/15/07) 
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) NYMD (INTEGER) : YYYYMMDD for which to read in A-6 fields
!
!  NOTES:
!  (1 ) Now modified for GCAP and GEOS-5 met fields (swu, bmy, 5/25/05)
!  (2 ) Remove support for GEOS-1 and GEOS-STRAT met fields (bmy, 8/4/06)
!  (3 ) Increase number of A-6 fields for GEOS-5 to 21 (bmy, 5/15/07)
!  28 Feb 2012 - R. Yantosca - Removed support for GEOS-3
!******************************************************************************
!
      USE CMN_SIZE_MOD 

      ! Function value
      INTEGER :: N_A6

      !=================================================================
      ! GET_N_A6 begins here!
      !=================================================================
#if   defined( GEOS_4 )

      ! GEOS-4 has 12 A-6 fields
      N_A6 = 12

#elif defined( GEOS_5 )

      ! GEOS-5 has 19 A-6 fields
      N_A6 = 21

#elif defined( GCAP )
      
      ! GCAP has 14 A-6 fields
      N_A6 = 14

#endif

      ! Return to calling program
      END FUNCTION GET_N_A6

!------------------------------------------------------------------------------

      FUNCTION CHECK_TIME( XYMD, XHMS, NYMD, NHMS ) RESULT( ITS_TIME )
!
!******************************************************************************
!  Function CHECK_TIME checks to see if the timestamp of the A-3 field just
!  read from disk matches the current time.  If so, then it's time to return
!  the A-3 field to the calling program. (bmy, 6/19/03, 8/4/06)
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
      USE CMN_SIZE_MOD

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

      SUBROUTINE READ_A6( NYMD,      NHMS,   
     &                    CLDF,      CLDMAS,    CLDTOPS,   CMFMC, 
     &                    DETRAINE,  DETRAINN,  DNDE,      DNDN,     
     &                    DQIDTMST,  DQLDTMST,  DQRCON,    DQRLSC,      
     &                    DQVDTMST,  DTRAIN,    ENTRAIN,   HKBETA,    
     &                    HKETA,     MFXC,      MFYC,      MFZ,
     &                    MOISTQ,    OPTDEPTH,  PLE,       PV,
     &                    Q,         QI,        QL,        RH,        
     &                    T,         TAUCLI,    TAUCLW,    U,         
     &                    UPDE,      UPDN,      V,         ZMEU,      
     &                    ZMMD,      ZMMU )
!
!******************************************************************************
!  Subroutine READ_A6 reads A-6 (avg 6-hr) met fields from disk. 
!  (bmy, 6/5/98, 10/15/09)
! 
!  Arguments as input:
!  ===========================================================================
!  (1 ) NYMD     : YYYYMMDD
!  (2 ) NHMS     :  and HHMMSS of A-6 met fields to be accessed
!
!  A-6 Met Fields as Output (Optional Arguments):
!  ============================================================================
!  (3 ) CLDF     : (3-D) Total cloud fractions                 [unitless]
!  (4 ) CLDMAS   : (3-D) Cloud mass flux field                 [kg/m2/600s]
!  (5 ) CLDTOPS  : (2-D) CTM Level in which cloud top occurs   [unitless]
!  (6 ) CMFMC    : (3-D) GEOS-5 cloud mass flux                [kg/m2/s]
!  (7 ) DETRAINE : (3-D) GCAP detrainment (entraining plume)   [kg/m2/s]
!  (8 ) DETRAINN : (3-D) GCAP detrainment (non-entr'n plume)
!  (9 ) DNDE     : (3-D) GCAP downdraft   (entraining plume)
!  (10) DNDN     : (3-D) GCAP downdraft   (non-entr'n plume)
!  (11) DQIDTMST : (3-D) GEOS-5 ice tendency in moist proc     [kg/kg/s]
!  (12) DQLDTMST : (3-D) GEOS-5 liquid tendency in moist proc  [kg/kg/s] 
!  (13) DQRCON   : (3-D) GEOS-5 precip formation rate / conv
!  (14) DQRLSC   : (3-D) GEOS-5 precip formation rate / lg scl 
!  (15) DQVDTMST : (3-D) GEOS-5 vapor tendency in moist proc   [kg/kg/s] 
!  (16) DTRAIN   : (3-D) Detrainment field                     [kg/m2/s]
!  (17) ENTRAIN  : (3-D) GCAP entrainment 
!  (18) HKBETA   : (3-D) Hack overshoot parameter              [unitless]
!  (19) HKETA    : (3-D) Hack convective mass flux             [kg/m2/s]
!  (20) MFXC     : (3-D) GEOS-5 E-W mass flux                  [Pa*m2/s]
!  (21) MFYC     : (3-D) GEOS-5 N-S mass flux                  [Pa*m2/s]
!  (22) MFZ      : (3-D) GEOS-5 up/down mass flux              [kg/m2/s]
!  (23) MOISTQ   : (3-D) DAO water vapor tendency d            [g/kg/day]
!  (24) OPTDEPTH : (3-D) GEOS grid box optical depth           [unitless]
!  (25) PLE      : (3-D) GEOS-5 pressure edges                 [hPa]
!  (26) PV       : (3-D) GEOS-5 potential vorticity            [kg*m2/kg/s]
!  (27) Q        : (3-D) Specific humidity                     [g H2O/kg air]
!  (28) T        : (3-D) Temperature                           [K]
!  (29) TAUCLI   : (3-D) GEOS ice path optical depth           [unitless]
!  (30) TAUCLW   : (3-D) GEOS water path optical depth         [unitless]
!  (31) U        : (3-D) Zonal winds                           [m/s]
!  (32) UPDE     : (3-D) GCAP updraft (entraining plume)
!  (33) UPDN     : (3-D) GCAP updraft (non-entr'n plume)
!  (34) V        : (3-D) Meridional winds                      [m/s]
!  (35) ZMEU     : (3-D) Zhang/McFarlane updraft entrainment   [Pa/s]
!  (36) ZMMD     : (3-D) Zhang/McFarlane downdraft mass flux   [Pa/s]
!  (37) ZMMU     : (3-D) Zhang/McFarlane updraft mass flux     [Pa/s]
!
!  NOTES:
!  (1 ) Adapted from READ_A6 of "dao_read_mod.f" (bmy, 6/19/03)
!  (2 ) Now use function TIMESTAMP_STRING from "time_mod.f" for formatted 
!        date/time output. (bmy, 10/28/03)
!  (3 ) Now compute CLDTOPS using ZMMU for GEOS-4 (bmy, 3/4/04)
!  (4 ) Now modified for GEOS-5 and GCAP fields.  Added DETRAINE, 
!        DETRAINN, DNDE, DNDN, ENTRAIN, UPDE, UPDN as optional arguments.
!        Now references "CMN_DIAG". (swu, bmy, 5/25/05)
!  (5 ) Bug fix in ND66 diagnostic for GEOS-4 (bmy, 2/1/06)
!  (6 ) Remove support for GEOS-1 and GEOS-STRAT met fields (bmy, 8/4/06)
!  (7 ) Now set negative SPHU to a small positive # (1d-32) instead of zero, 
!        so as not to blow up logarithms (bmy, 9/8/06)
!  (8 ) Add CMFMC, DQIDTMST, DQLDTMST, DQRCON, DQRLSC, DQVDTMST, MFXC, MFYC, 
!        MFZ, PLE, PV, RH, TAUCLI, and TAUCLW as optional arguments.  Also 
!        update the CASE statement accordingly for GEOS-5 met fields. 
!        Now reference TRANSFER_3D_Lp1 from "transfer_mod.f".  Now convert
!        GEOS-5 specific humidity from [kg/kg] to [g/kg] for compatibility
!        with existing routines.  Also recognize EPV, which is an alternate 
!        name for PV.  Bug fix: convert GEOS-5 RH from unitless to %.
!        (phs, bmy, 3/28/08)
!  (8 ) Now get the # of A-6 fields from the file ident string (bmy, 10/7/08)
!  (9 ) Remove references to IN_CLOUD_OD (bmy, 10/15/09)
!  28 Feb 2012 - R. Yantosca - Removed support for GEOS-3
!******************************************************************************
!
      ! References to F90 modules
      USE DIAG_MOD,     ONLY : AD66,        AD67
      USE FILE_MOD,     ONLY : IOERROR,     IU_A6
      USE TIME_MOD,     ONLY : SET_CT_A6,   TIMESTAMP_STRING
      USE TRANSFER_MOD, ONLY : TRANSFER_A6, TRANSFER_3D_Lp1
      USE TRANSFER_MOD, ONLY : TRANSFER_3D, TRANSFER_G5_PLE

      USE CMN_SIZE_MOD               ! Size parameters
      USE CMN_DIAG_MOD               ! ND66, ND67
      USE CMN_GCTM_MOD               ! g0

      ! Arguments
      INTEGER, INTENT(IN)            :: NYMD, NHMS
      INTEGER, INTENT(OUT), OPTIONAL :: CLDTOPS(IIPAR,JJPAR)
      REAL*8,  INTENT(OUT), OPTIONAL :: CLDF(LLPAR,IIPAR,JJPAR)
      REAL*8,  INTENT(OUT), OPTIONAL :: CLDMAS(IIPAR,JJPAR,LLPAR)
      REAL*8,  INTENT(OUT), OPTIONAL :: CMFMC(IIPAR,JJPAR,LLPAR+1)
      REAL*8,  INTENT(OUT), OPTIONAL :: DETRAINE(IIPAR,JJPAR,LLPAR)
      REAL*8,  INTENT(OUT), OPTIONAL :: DETRAINN(IIPAR,JJPAR,LLPAR)
      REAL*8,  INTENT(OUT), OPTIONAL :: DNDE(IIPAR,JJPAR,LLPAR)
      REAL*8,  INTENT(OUT), OPTIONAL :: DNDN(IIPAR,JJPAR,LLPAR)
      REAL*8,  INTENT(OUT), OPTIONAL :: DQIDTMST(IIPAR,JJPAR,LLPAR)
      REAL*8,  INTENT(OUT), OPTIONAL :: DQLDTMST(IIPAR,JJPAR,LLPAR)
      REAL*8,  INTENT(OUT), OPTIONAL :: DQRCON(IIPAR,JJPAR,LLPAR)
      REAL*8,  INTENT(OUT), OPTIONAL :: DQRLSC(IIPAR,JJPAR,LLPAR)
      REAL*8,  INTENT(OUT), OPTIONAL :: DQVDTMST(IIPAR,JJPAR,LLPAR)
      REAL*8,  INTENT(OUT), OPTIONAL :: DTRAIN(IIPAR,JJPAR,LLPAR)
      REAL*8,  INTENT(OUT), OPTIONAL :: ENTRAIN(IIPAR,JJPAR,LLPAR)
      REAL*8,  INTENT(OUT), OPTIONAL :: HKBETA(IIPAR,JJPAR,LLPAR)
      REAL*8,  INTENT(OUT), OPTIONAL :: HKETA(IIPAR,JJPAR,LLPAR)
      REAL*8,  INTENT(OUT), OPTIONAL :: MFXC(IIPAR,JJPAR,LLPAR)
      REAL*8,  INTENT(OUT), OPTIONAL :: MFYC(IIPAR,JJPAR,LLPAR)
      REAL*8,  INTENT(OUT), OPTIONAL :: MFZ(IIPAR,JJPAR,LLPAR+1)
      REAL*8,  INTENT(OUT), OPTIONAL :: MOISTQ(LLPAR,IIPAR,JJPAR) 
      REAL*8,  INTENT(OUT), OPTIONAL :: OPTDEPTH(LLPAR,IIPAR,JJPAR)
      REAL*8,  INTENT(OUT), OPTIONAL :: PLE(IIPAR,JJPAR,LLPAR+1)
      REAL*8,  INTENT(OUT), OPTIONAL :: PV(IIPAR,JJPAR,LLPAR)
      REAL*8,  INTENT(OUT), OPTIONAL :: Q(IIPAR,JJPAR,LLPAR)
      REAL*8,  INTENT(OUT), OPTIONAL :: QI(IIPAR,JJPAR,LLPAR)
      REAL*8,  INTENT(OUT), OPTIONAL :: QL(IIPAR,JJPAR,LLPAR)
      REAL*8,  INTENT(OUT), OPTIONAL :: RH(IIPAR,JJPAR,LLPAR)
      REAL*8,  INTENT(OUT), OPTIONAL :: T(IIPAR,JJPAR,LLPAR)
      REAL*8,  INTENT(OUT), OPTIONAL :: TAUCLI(IIPAR,JJPAR,LLPAR)
      REAL*8,  INTENT(OUT), OPTIONAL :: TAUCLW(IIPAR,JJPAR,LLPAR)
      REAL*8,  INTENT(OUT), OPTIONAL :: U(IIPAR,JJPAR,LLPAR)
      REAL*8,  INTENT(OUT), OPTIONAL :: UPDE(IIPAR,JJPAR,LLPAR)
      REAL*8,  INTENT(OUT), OPTIONAL :: UPDN(IIPAR,JJPAR,LLPAR)
      REAL*8,  INTENT(OUT), OPTIONAL :: V(IIPAR,JJPAR,LLPAR)
      REAL*8,  INTENT(OUT), OPTIONAL :: ZMEU(IIPAR,JJPAR,LLPAR)
      REAL*8,  INTENT(OUT), OPTIONAL :: ZMMD(IIPAR,JJPAR,LLPAR)
      REAL*8,  INTENT(OUT), OPTIONAL :: ZMMU(IIPAR,JJPAR,LLPAR)

      ! Local variables
      INTEGER                        :: I, IJLOOP, J, K, L
      INTEGER                        :: IOS, NFOUND, N_A6
      REAL*4                         :: D(IIPAR,JJPAR,LGLOB)       
      REAL*4                         :: D1(IIPAR,JJPAR,LGLOB+1)       
      REAL*8                         :: C1, C2
      REAL*8                         :: TAUCLD(LLPAR,IIPAR,JJPAR)
      REAL*8                         :: CLDTOT(LLPAR,IIPAR,JJPAR)
      CHARACTER(LEN=8)               :: NAME
      CHARACTER(LEN=16)              :: STAMP
      INTEGER                        :: XYMD, XHMS

      !=================================================================
      ! READ_A6 begins here!      
      !=================================================================

      ! Get number of A-6 fields
#if   defined( GEOS_5 )
      N_A6 = N_A6_FIELDS
#else
      N_A6 = GET_N_A6()
#endif

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

            !--------------------------------------
            ! CLDMAS: GEOS-3 cloud mass flux
            !--------------------------------------
            CASE ( 'CLDMAS' ) 
               READ( IU_A6, IOSTAT=IOS ) XYMD, XHMS, D
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A6, 'read_a6:2' )
             
               IF ( CHECK_TIME( XYMD, XHMS, NYMD, NHMS ) ) THEN
                  IF ( PRESENT( CLDMAS ) ) CALL TRANSFER_3D( D, CLDMAS )
                  NFOUND = NFOUND + 1
               ENDIF

            !--------------------------------------
            ! CLDF  : GEOS-3 3-D total cloud frac
            ! CLDTOT: GEOS-4 3-D total cloud frac
            ! CLOUD : GEOS-5 3-D total cloud frac
            !--------------------------------------
            CASE ( 'CLDTOT', 'CLDF', 'CLOUD' ) 
               READ( IU_A6, IOSTAT=IOS ) XYMD, XHMS, D
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A6, 'read_a6:3' )

               IF ( CHECK_TIME( XYMD, XHMS, NYMD, NHMS ) ) THEN
                  CALL TRANSFER_A6( D, CLDTOT )
                  NFOUND = NFOUND +1 
               ENDIF

            !------------------------------------
            ! CMFMC: GEOS-5 cloud mass flux
            !------------------------------------
            CASE ( 'CMFMC' ) 
               READ( IU_A6, IOSTAT=IOS ) XYMD, XHMS, D1
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A6, 'read_a6:4' )

               IF ( CHECK_TIME( XYMD, XHMS, NYMD, NHMS ) ) THEN
                  IF ( PRESENT( CMFMC ) ) THEN
                     CALL TRANSFER_3D_Lp1( D1, CMFMC )
                  ENDIF
                  NFOUND = NFOUND + 1 
               ENDIF

            !--------------------------------------
            ! DETRAINE: GCAP Detrainment (ent pl)
            !--------------------------------------
            CASE ( 'DETRAINE' ) 
               READ( IU_A6, IOSTAT=IOS ) XYMD, XHMS, D
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A6, 'read_a6:5' )
 
               IF ( CHECK_TIME( XYMD, XHMS, NYMD, NHMS ) ) THEN
                  IF ( PRESENT( DETRAINE ) ) THEN
                     CALL TRANSFER_3D( D, DETRAINE )
                  ENDIF
                  NFOUND = NFOUND + 1
               ENDIF

            !--------------------------------------
            ! DETRAINN: GCAP Detrainment (non-ent)
            !--------------------------------------
            CASE ( 'DETRAINN' ) 
               READ( IU_A6, IOSTAT=IOS ) XYMD, XHMS, D
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A6, 'read_a6:6' )
 
               IF ( CHECK_TIME( XYMD, XHMS, NYMD, NHMS ) ) THEN
                  IF ( PRESENT( DETRAINN ) ) THEN
                     CALL TRANSFER_3D( D, DETRAINN )
                  ENDIF
                  NFOUND = NFOUND + 1
               ENDIF

            !--------------------------------------
            ! DNDE: GCAP Downdraft (ent plume)
            !--------------------------------------
            CASE ( 'DNDE' ) 
               READ( IU_A6, IOSTAT=IOS ) XYMD, XHMS, D
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A6, 'read_a6:7' )
 
               IF ( CHECK_TIME( XYMD, XHMS, NYMD, NHMS ) ) THEN
                  IF ( PRESENT( DNDE ) ) CALL TRANSFER_3D( D, DNDE )
                  NFOUND = NFOUND + 1
               ENDIF
               
            !--------------------------------------
            ! DNDN: GCAP Downdraft (non-ent plume)
            !--------------------------------------
            CASE ( 'DNDN' ) 
               READ( IU_A6, IOSTAT=IOS ) XYMD, XHMS, D
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A6, 'read_a6:8' )
 
               IF ( CHECK_TIME( XYMD, XHMS, NYMD, NHMS ) ) THEN
                  IF ( PRESENT( DNDN ) ) CALL TRANSFER_3D( D, DNDN )
                  NFOUND = NFOUND + 1
               ENDIF

            !--------------------------------------
            ! DQIDTMST: GEOS-5 ice tend in moist p 
            !--------------------------------------
            CASE ( 'DQIDTMST' ) 
               READ( IU_A6, IOSTAT=IOS ) XYMD, XHMS, D
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A6, 'read_a6:9' )
 
               IF ( CHECK_TIME( XYMD, XHMS, NYMD, NHMS ) ) THEN
                  IF ( PRESENT( DQIDTMST ) ) THEN 
                     CALL TRANSFER_3D( D, DQIDTMST )
                  ENDIF
                  NFOUND = NFOUND + 1
               ENDIF

            !--------------------------------------
            ! DQLDTMST: GEOS-5 liq tend in moist p
            !--------------------------------------
            CASE ( 'DQLDTMST' ) 
               READ( IU_A6, IOSTAT=IOS ) XYMD, XHMS, D
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A6, 'read_a6:10' )
 
               IF ( CHECK_TIME( XYMD, XHMS, NYMD, NHMS ) ) THEN
                  IF ( PRESENT( DQLDTMST ) ) THEN
                     CALL TRANSFER_3D( D, DQLDTMST )
                  ENDIF
                  NFOUND = NFOUND + 1
               ENDIF

            !--------------------------------------
            ! DQRCON: GEOS-5 conv rain prod rate
            !--------------------------------------
            CASE ( 'DQRCON' ) 
               READ( IU_A6, IOSTAT=IOS ) XYMD, XHMS, D
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A6, 'read_a6:11' )
 
               IF ( CHECK_TIME( XYMD, XHMS, NYMD, NHMS ) ) THEN
                  IF ( PRESENT( DQRCON ) ) THEN 
                     CALL TRANSFER_3D( D, DQRCON )
                  ENDIF
                  NFOUND = NFOUND + 1
               ENDIF

            !--------------------------------------
            ! DQRLSC: GEOS-5 lg scl rain prod rate
            !--------------------------------------
            CASE ( 'DQRLSC' ) 
               READ( IU_A6, IOSTAT=IOS ) XYMD, XHMS, D
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A6, 'read_a6:12' )
 
               IF ( CHECK_TIME( XYMD, XHMS, NYMD, NHMS ) ) THEN
                  IF ( PRESENT( DQRLSC ) ) THEN 
                     CALL TRANSFER_3D( D, DQRLSC )
                  ENDIF
                  NFOUND = NFOUND + 1
               ENDIF

            !--------------------------------------
            ! DQVDTMST: GEOS-5 vap tend in moist p
            !--------------------------------------
            CASE ( 'DQVDTMST' ) 
               READ( IU_A6, IOSTAT=IOS ) XYMD, XHMS, D
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A6, 'read_a6:13' )
 
               IF ( CHECK_TIME( XYMD, XHMS, NYMD, NHMS ) ) THEN
                  IF ( PRESENT( DQVDTMST ) ) THEN 
                     CALL TRANSFER_3D( D, DQVDTMST )
                  ENDIF
                  NFOUND = NFOUND + 1
               ENDIF

            !--------------------------------------
            ! DTRAIN: GEOS-3 & GEOS-5 detrainment
            !--------------------------------------
            CASE ( 'DTRAIN' ) 
               READ( IU_A6, IOSTAT=IOS ) XYMD, XHMS, D
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A6, 'read_a6:14' )
 
               IF ( CHECK_TIME( XYMD, XHMS, NYMD, NHMS ) ) THEN
                  IF ( PRESENT( DTRAIN ) ) CALL TRANSFER_3D( D, DTRAIN )
                  NFOUND = NFOUND + 1
               ENDIF

            !--------------------------------------
            ! ENTRAIN: GCAP Entrainment
            !--------------------------------------
            CASE ( 'ENTRAIN' ) 
               READ( IU_A6, IOSTAT=IOS ) XYMD, XHMS, D
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A6, 'read_a6:15' )
 
               IF ( CHECK_TIME( XYMD, XHMS, NYMD, NHMS ) ) THEN
                  IF ( PRESENT( ENTRAIN ) ) THEN
                     CALL TRANSFER_3D( D, ENTRAIN )
                  ENDIF
                  NFOUND = NFOUND + 1
               ENDIF

            !--------------------------------------
            ! HKBETA: GEOS-4 Hack overshoot param. 
            !--------------------------------------
            CASE ( 'HKBETA' ) 
               READ( IU_A6, IOSTAT=IOS ) XYMD, XHMS, D
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A6, 'read_a6:16' )
 
               IF ( CHECK_TIME( XYMD, XHMS, NYMD, NHMS ) ) THEN
                  IF ( PRESENT( HKBETA ) ) CALL TRANSFER_3D( D, HKBETA )
                  NFOUND = NFOUND + 1
               ENDIF

            !--------------------------------------
            ! HKETA: GEOS-4 Hack conv mass flux 
            !--------------------------------------
            CASE ( 'HKETA' ) 
               READ( IU_A6, IOSTAT=IOS ) XYMD, XHMS, D
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A6, 'read_a6:17' )
 
               IF ( CHECK_TIME( XYMD, XHMS, NYMD, NHMS ) ) THEN
                  IF ( PRESENT( HKETA ) ) CALL TRANSFER_3D( D, HKETA )
                  NFOUND = NFOUND + 1
               ENDIF

            !--------------------------------------
            ! MFXC: GEOS-5 E-W mass flux (C-grid)
            !--------------------------------------
            CASE ( 'MFXC' ) 
               READ( IU_A6, IOSTAT=IOS ) XYMD, XHMS, D
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A6, 'read_a6:18' )

               IF ( CHECK_TIME( XYMD, XHMS, NYMD, NHMS ) ) THEN
                  IF ( PRESENT( MFXC ) ) CALL TRANSFER_3D( D, MFXC )
                  NFOUND = NFOUND + 1 
               ENDIF

            !--------------------------------------
            ! MFYC: GEOS-5 N-S mass flux (C-grid)
            !--------------------------------------
            CASE ( 'MFYC' ) 
               READ( IU_A6, IOSTAT=IOS ) XYMD, XHMS, D
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A6, 'read_a6:19' )

               IF ( CHECK_TIME( XYMD, XHMS, NYMD, NHMS ) ) THEN
                  IF ( PRESENT( MFYC ) ) CALL TRANSFER_3D( D, MFYC )
                  NFOUND = NFOUND + 1 
               ENDIF

            !--------------------------------------
            ! MFZ: GEOS-5 vert mass flux (C-grid)
            !--------------------------------------
            CASE ( 'MFZ' ) 
               READ( IU_A6, IOSTAT=IOS ) XYMD, XHMS, D1
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A6, 'read_a6:20' )

               IF ( CHECK_TIME( XYMD, XHMS, NYMD, NHMS ) ) THEN
                  IF ( PRESENT( MFZ ) ) CALL TRANSFER_3D_Lp1( D1, MFZ )
                  NFOUND = NFOUND + 1 
               ENDIF

            !--------------------------------------
            ! MOISTQ: tendency of SPHU
            !--------------------------------------
            CASE ( 'MOISTQ' ) 
               READ( IU_A6, IOSTAT=IOS ) XYMD, XHMS, D
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A6, 'read_a6:21' )

               IF ( CHECK_TIME( XYMD, XHMS, NYMD, NHMS ) ) THEN
                  IF ( PRESENT( MOISTQ ) ) CALL TRANSFER_A6( D, MOISTQ )
                  NFOUND = NFOUND + 1 
               ENDIF

            !--------------------------------------
            ! OPTDEPTH: grid box optical depth
            !--------------------------------------
            CASE ( 'OPTDEPTH' ) 
               READ( IU_A6, IOSTAT=IOS ) XYMD, XHMS, D
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A6, 'read_a6:22' )

               IF ( CHECK_TIME( XYMD, XHMS, NYMD, NHMS ) ) THEN
                  IF ( PRESENT( OPTDEPTH ) ) THEN
                     CALL TRANSFER_A6( D, OPTDEPTH )
                  ENDIF
                  NFOUND = NFOUND + 1 
               ENDIF

            !--------------------------------------
            ! PLE: GEOS-5 pressure edges
            !--------------------------------------
            CASE ( 'PLE' ) 
               READ( IU_A6, IOSTAT=IOS ) XYMD, XHMS, D1
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A6, 'read_a6:23' )

               IF ( CHECK_TIME( XYMD, XHMS, NYMD, NHMS ) ) THEN
                  IF ( PRESENT( PLE ) ) CALL TRANSFER_G5_PLE( D1, PLE )
                  NFOUND = NFOUND + 1 
               ENDIF

            !--------------------------------------
            ! PV: GEOS-5 Ertel potential vorticity
            !--------------------------------------
            CASE ( 'PV', 'EPV' ) 
               READ( IU_A6, IOSTAT=IOS ) XYMD, XHMS, D
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A6, 'read_a6:24' )

               IF ( CHECK_TIME( XYMD, XHMS, NYMD, NHMS ) ) THEN
                  IF ( PRESENT( PV ) ) CALL TRANSFER_3D( D, PV )
                  NFOUND = NFOUND + 1 
               ENDIF

            !--------------------------------------
            ! Q:  GEOS-4 specific humidity [g/kg]
            ! QV: GEOS-5 specific humidity [kg/kg]
            !--------------------------------------
            CASE ( 'Q', 'QV' ) 
               READ( IU_A6, IOSTAT=IOS ) XYMD, XHMS, D
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A6, 'read_a6:25' )

               IF ( CHECK_TIME( XYMD, XHMS, NYMD, NHMS ) ) THEN
                  IF ( PRESENT( Q ) ) CALL TRANSFER_3D( D, Q )
                  NFOUND = NFOUND + 1 

                  ! NOTE: Now set negative Q to a small positive # 
                  ! instead of zero, so as not to blow up logarithms
                  ! (bmy, 9/8/06)
                  WHERE ( Q < 0d0 ) Q = 1d-32
               ENDIF

            !--------------------------------------
            ! QI: GEOS-5 ice mixing ratio [kg/kg]
            !--------------------------------------
            CASE ( 'QI' ) 
               READ( IU_A6, IOSTAT=IOS ) XYMD, XHMS, D
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A6, 'read_a6:26' )
            
               IF ( CHECK_TIME( XYMD, XHMS, NYMD, NHMS ) ) THEN
                  IF ( PRESENT( QI ) ) CALL TRANSFER_3D( D, QI )
                  NFOUND = NFOUND + 1 
               ENDIF
            
            !--------------------------------------
            ! QL: GEOS-5 water mix ratio [kg/kg]
            !--------------------------------------
            CASE ( 'QL' ) 
               READ( IU_A6, IOSTAT=IOS ) XYMD, XHMS, D
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A6, 'read_a6:27' )
            
               IF ( CHECK_TIME( XYMD, XHMS, NYMD, NHMS ) ) THEN
                  IF ( PRESENT( QL ) ) CALL TRANSFER_3D( D, QL )
                  NFOUND = NFOUND + 1 
               ENDIF

            !--------------------------------------
            ! RH: GEOS-5 relative humidity [%]
            !--------------------------------------
            CASE ( 'RH' ) 
               READ( IU_A6, IOSTAT=IOS ) XYMD, XHMS, D
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A6, 'read_a6:28' )

               IF ( CHECK_TIME( XYMD, XHMS, NYMD, NHMS ) ) THEN
                  IF ( PRESENT( RH ) ) CALL TRANSFER_3D( D, RH )
                  NFOUND = NFOUND + 1 
               ENDIF

            !--------------------------------------
            ! T: 3-D temperature
            !--------------------------------------
            CASE ( 'T' ) 
               READ( IU_A6, IOSTAT=IOS ) XYMD, XHMS, D
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A6, 'read_a6:29' )

               IF ( CHECK_TIME( XYMD, XHMS, NYMD, NHMS ) ) THEN
                  IF ( PRESENT( T ) ) CALL TRANSFER_3D( D, T )
                  NFOUND = NFOUND + 1 
               ENDIF

            !--------------------------------------
            ! TAUCLI: GEOS-5 ice path opt depth
            !--------------------------------------
            CASE ( 'TAUCLI' ) 
               READ( IU_A6, IOSTAT=IOS ) XYMD, XHMS, D
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A6, 'read_a6:30' )

               IF ( CHECK_TIME( XYMD, XHMS, NYMD, NHMS ) ) THEN
                  IF ( PRESENT( TAUCLI ) ) CALL TRANSFER_3D( D, TAUCLI )
                  NFOUND = NFOUND + 1 
               ENDIF

            !--------------------------------------
            ! TAUCLW: GEOS-5 water path opt depth
            !--------------------------------------
            CASE ( 'TAUCLW' ) 
               READ( IU_A6, IOSTAT=IOS ) XYMD, XHMS, D
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A6, 'read_a6:31' )

               IF ( CHECK_TIME( XYMD, XHMS, NYMD, NHMS ) ) THEN
                  IF ( PRESENT( TAUCLW ) ) CALL TRANSFER_3D( D, TAUCLW )
                  NFOUND = NFOUND + 1 
               ENDIF

            !--------------------------------------
            ! U: GEOS-4 & GEOS-5 zonal wind
            !--------------------------------------
            CASE ( 'U' ) 
               READ( IU_A6, IOSTAT=IOS ) XYMD, XHMS, D
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A6, 'read_a6:32' )

               IF ( CHECK_TIME( XYMD, XHMS, NYMD, NHMS ) ) THEN
                  IF ( PRESENT( U ) ) CALL TRANSFER_3D( D, U )
                  NFOUND = NFOUND + 1 
               ENDIF

            !--------------------------------------
            ! UPDE: GCAP Downdraft (ent plume)
            !--------------------------------------
            CASE ( 'UPDE' ) 
               READ( IU_A6, IOSTAT=IOS ) XYMD, XHMS, D
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A6, 'read_a6:33' )
 
               IF ( CHECK_TIME( XYMD, XHMS, NYMD, NHMS ) ) THEN
                  IF ( PRESENT( UPDE ) ) CALL TRANSFER_3D( D, UPDE )
                  NFOUND = NFOUND + 1
               ENDIF
               
            !--------------------------------------
            ! UPDN: Downdraft (non-ent plume)
            !--------------------------------------
            CASE ( 'UPDN' ) 
               READ( IU_A6, IOSTAT=IOS ) XYMD, XHMS, D
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A6, 'read_a6:34' )
 
               IF ( CHECK_TIME( XYMD, XHMS, NYMD, NHMS ) ) THEN
                  IF ( PRESENT( UPDN ) ) CALL TRANSFER_3D( D, UPDN )
                  NFOUND = NFOUND + 1
               ENDIF

            !--------------------------------------
            ! V: GEOS-4 & GEOS-5 meridional wind
            !--------------------------------------
            CASE ( 'V' ) 
               READ( IU_A6, IOSTAT=IOS ) XYMD, XHMS, D
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A6, 'read_a6:35' )

               IF ( CHECK_TIME( XYMD, XHMS, NYMD, NHMS ) ) THEN
                  IF ( PRESENT( V ) ) CALL TRANSFER_3D( D, V )
                  NFOUND = NFOUND + 1 
               ENDIF

            !--------------------------------------
            ! ZMEU: GEOS-4 Z&M updraft entrainment
            !--------------------------------------
            CASE ( 'ZMEU' ) 
               READ( IU_A6, IOSTAT=IOS ) XYMD, XHMS, D
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A6, 'read_a6:36' )

               IF ( CHECK_TIME( XYMD, XHMS, NYMD, NHMS ) ) THEN
                  IF ( PRESENT( ZMEU ) ) CALL TRANSFER_3D( D, ZMEU )
                  NFOUND = NFOUND + 1 
               ENDIF

            !--------------------------------------
            ! ZMMD: GEOS-4 Z&M downdraft mass flux
            !--------------------------------------
            CASE ( 'ZMMD' ) 
               READ( IU_A6, IOSTAT=IOS ) XYMD, XHMS, D
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A6, 'read_a6:37' )

               IF ( CHECK_TIME( XYMD, XHMS, NYMD, NHMS ) ) THEN
                  IF ( PRESENT( ZMMD ) ) CALL TRANSFER_3D( D, ZMMD )
                  NFOUND = NFOUND + 1 
               ENDIF

            !--------------------------------------
            ! ZMMU: GEOS-4 Z&M updraft mass flux
            !--------------------------------------
            CASE ( 'ZMMU' ) 
               READ( IU_A6, IOSTAT=IOS ) XYMD, XHMS, D
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A6, 'read_a6:38' )

               IF ( CHECK_TIME( XYMD, XHMS, NYMD, NHMS ) ) THEN
                  IF ( PRESENT( ZMMU ) ) CALL TRANSFER_3D( D, ZMMU )
                  NFOUND = NFOUND + 1 
               ENDIF

            !--------------------------------------
            ! TAUCLD: in-cloud optical depth 
            ! Just skip over this
            !--------------------------------------
            CASE ( 'TAUCLD' ) 
               READ( IU_A6, IOSTAT=IOS ) XYMD, XHMS, D
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A6, 'read_a6:39' )

               IF ( CHECK_TIME( XYMD, XHMS, NYMD, NHMS ) ) THEN
                  NFOUND = NFOUND + 1 
               ENDIF

            !--------------------------------------
            ! KH: Just skip over this
            !--------------------------------------
            CASE ( 'KH' ) 
               READ( IU_A6, IOSTAT=IOS ) XYMD, XHMS, D
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A6, 'read_a6:40' )

               IF ( CHECK_TIME( XYMD, XHMS, NYMD, NHMS ) ) THEN
                  NFOUND = NFOUND + 1 
               ENDIF

            !--------------------------------------
            ! Extra GEOS-5 fields
            ! Skip over these now; add later
            !--------------------------------------
            CASE ( 'OMEGA' ) 
               READ( IU_A6, IOSTAT=IOS ) XYMD, XHMS, D
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A6, 'read_a6:40' )

               IF ( CHECK_TIME( XYMD, XHMS, NYMD, NHMS ) ) THEN
                  NFOUND = NFOUND + 1 
               ENDIF

            ! Field not found -- skip over
            CASE DEFAULT
               WRITE ( 6, '(a)' ) 'Searching for next A-6 field!'
               READ( IU_A6, IOSTAT=IOS ) XYMD, XHMS, D
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A6, 'read_a6:41' )

         END SELECT

         !==============================================================
         ! If we have found all the fields for this time, then exit 
         ! the loop.  Otherwise, go on to the next iteration.
         !==============================================================
         IF ( CHECK_TIME( XYMD, XHMS, NYMD, NHMS ) .AND. 
     &        NFOUND == N_A6 ) THEN
            STAMP = TIMESTAMP_STRING( NYMD, NHMS )
            WRITE( 6, 210 ) NFOUND, STAMP
 210        FORMAT( '     - Found all ', i3, ' A-6 met fields for ', a )
            EXIT
         ENDIF
      ENDDO

      !=================================================================
      ! CLDTOPS(I,J) = level of convective cloud top at (I,J).
      ! GEOS-CHEM cloud top at (I,J) is at top of first level where 
      ! cloud mass flux goes from being nonzero to zero.  
      !
      ! For GEOS-3 : mass flux is "CLDMAS" field
      ! For GEOS-4 : mass flux is "ZMMU"   field
      ! For GEOS-5 : mass flux is "CMFMC"  field
      ! For GCAP   : mass flux is "UPDN"   field
      !=================================================================
#if   defined( GCAP )

      !------------------------------
      ! Special handling for GCAP
      !------------------------------

      ! CLDTOPS is highest location of ZMMU in the column (I,J)
      IF ( PRESENT( CLDTOPS ) .and. PRESENT( UPDN ) ) THEN
         DO J = 1, JJPAR
         DO I = 1, IIPAR
            CLDTOPS(I,J) = 1
            DO L = LLPAR, 1, -1
               IF ( UPDN(I,J,L) > 0d0 ) THEN
                  CLDTOPS(I,J) = L + 1
                  EXIT
               ENDIF
            ENDDO         
         ENDDO
         ENDDO
      ENDIF     

#elif defined( GEOS_4 )

      !------------------------------
      ! Special handling for GEOS-4
      !------------------------------

      ! CLDTOPS is highest location of ZMMU in the column (I,J)
      IF ( PRESENT( CLDTOPS ) .and. PRESENT( ZMMU ) ) THEN
         DO J = 1, JJPAR
         DO I = 1, IIPAR
            CLDTOPS(I,J) = 1
            DO L = LLPAR, 1, -1
               IF ( ZMMU(I,J,L) > 0d0 ) THEN
                  CLDTOPS(I,J) = L + 1
                  EXIT
               ENDIF
            ENDDO         
         ENDDO
         ENDDO
      ENDIF

#elif defined( GEOS_5 )

      !------------------------------
      ! Special handling for GEOS-5
      !------------------------------

      ! Convert RH from unitless to percent (phs, bmy, 3/28/08)
      ! %%% NOTE: GEOS-5 file spec says units of RH are % but that's wrong!
      ! Temporary fix: force RH to be positive (phs, 5/1/08)
      IF ( PRESENT( RH ) ) THEN
         RH = RH * 100d0
         RH = MAX(RH, 0D0)
      ENDIF

      ! Convert GEOS-5 specific humidity from [kg/kg] to [g/kg]
      IF ( PRESENT( Q ) ) Q = Q * 1000d0

      ! CLDTOPS highest location of CMFMC in the column (I,J)
      IF ( PRESENT( CLDTOPS ) .and. PRESENT( CMFMC ) ) THEN
         DO J = 1, JJPAR
         DO I = 1, IIPAR
            CLDTOPS(I,J) = 1
            DO L = LLPAR, 1, -1
               IF ( CMFMC(I,J,L) > 0d0 ) THEN
                  CLDTOPS(I,J) = L + 1
                  EXIT
               ENDIF
            ENDDO
         ENDDO
         ENDDO
      ENDIF

#endif

      ! CLDF is read directly from disk CLDTOT met field
      IF ( PRESENT( CLDF ) ) THEN
         CLDF = CLDTOT
      ENDIF

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
      ! (6 ) DTRAIN : Detrainment mass flux           [kg/m2/s]
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
         
         ! GEOS-3 cloud mass flux
         IF ( PRESENT( CLDMAS ) ) THEN 
            AD66(:,:,1:LD66,5) = AD66(:,:,1:LD66,5) + CLDMAS(:,:,1:LD66)
         ENDIF  
     
         ! GEOS-4 cloud mass flux
         IF ( PRESENT( ZMMU ) ) THEN
            AD66(:,:,1:LD66,5) = AD66(:,:,1:LD66,5) + ZMMU(:,:,1:LD66)
         ENDIF
      
         ! GEOS-5 cloud mass flux
         IF ( PRESENT( CMFMC ) ) THEN
            AD66(:,:,1:LD66,5) = AD66(:,:,1:LD66,5) + CMFMC(:,:,1:LD66)
         ENDIF

         ! GCAP cloud mass flux 
         IF ( PRESENT( UPDE ) ) THEN
            AD66(:,:,1:LD66,5) = AD66(:,:,1:LD66,5) +UPDE(:,:,1:LD66)/g0
         ENDIF

         ! GCAP cloud mass flux 
         IF ( PRESENT( UPDN ) ) THEN
            AD66(:,:,1:LD66,5) = AD66(:,:,1:LD66,5) +UPDN(:,:,1:LD66)/g0
         ENDIF

         ! GEOS-3 & GEOS-5 detrainment
         IF ( PRESENT( DTRAIN ) ) THEN 
            AD66(:,:,1:LD66,6) = AD66(:,:,1:LD66,6) + DTRAIN(:,:,1:LD66)
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
      
      ! End of module
      END MODULE A6_READ_MOD
