!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: a3_read_mod
!
! !DESCRIPTION: Module A3\_READ\_MOD contains routines that unzip, open, and 
!  read the GEOS-Chem A3 (avg 3-hour) met fields from disk. 
!\\
!\\
! !INTERFACE:
!
      MODULE A3_READ_MOD
!
! !USES:
!
      USE inquireMod, ONLY : findFreeLUN

      IMPLICIT NONE
      PRIVATE
!
! !PUBLIC DATA MEMBERS:
!
      PUBLIC  :: ARCHIVE_ND67_1D 
      PUBLIC  :: GET_A3_FIELDS
      PUBLIC  :: OPEN_A3_FIELDS
      PUBLIC  :: UNZIP_A3_FIELDS
!
! !PRIVATE DATA MEMBERS:
!
      PRIVATE :: A3_CHECK        
      PRIVATE :: CHECK_TIME      
      PRIVATE :: DO_OPEN_A3      
      PRIVATE :: GET_N_A3        
      PRIVATE :: READ_A3         
!
! !REMARKS:
!  This module reads GEOS-4, GEOS-5, and GCAP met fields
!  MERRA met fields are read in routines merra_*_mod.F
!  GEOS-FP met fields are read in geosfp_read_mod.F
!
! !REVISION HISTORY:
!  23 Jun 2003 - R. Yantosca - Initial version
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
!  (13) Remove references to IN_CLOUD_OD (bmy, 10/15/09)
!  21 Sep 2010 - R. Yantosca - Added ProTeX headers
!  28 Feb 2012 - R. Yantosca - Removed support for GEOS-3
!  03 Aug 2012 - R. Yantosca - Move calls to findFreeLUN out of DEVEL block
!  03 Aug 2012 - R. Yantosca - Now make IU_A3 a private module variable
!  15 Nov 2012 - R. Yantosca - Now replace dao_mod.F arrays with State_Met
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE TYPES:
!
      ! Number of A3 fields in the file
      INTEGER :: N_A3_FIELDS

      ! File unit for A-3 met fields
      INTEGER :: IU_A3

      CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: unzip_a3_fields 
!
! !DESCRIPTION: Subroutine UNZIP\_A3\_FIELDS invokes a FORTRAN system call 
!  to uncompress GEOS-Chem A3 met field files and store the uncompressed data 
!  in a temporary directory, where GEOS-CHEM can read them.  The original data 
!  files are not disturbed.  
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE UNZIP_A3_FIELDS( Input_Opt, OPTION, NYMD )
!
! !USES:
!
      USE BPCH2_MOD,          ONLY : GET_RES_EXT
      USE CMN_SIZE_MOD
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE ERROR_MOD,          ONLY : ERROR_STOP
      USE TIME_MOD,           ONLY : EXPAND_DATE
!
! !INPUT PARAMETERS: 
!
      CHARACTER(LEN=*),  INTENT(IN) :: OPTION    ! Unzip option
      INTEGER, OPTIONAL, INTENT(IN) :: NYMD      ! YYYY/MM/DD of file to unzip
      TYPE(OptInput),    INTENT(IN) :: Input_Opt ! Input Options object
!     
! !REVISION HISTORY: 
!  15 Jun 1998 - R. Yantosca - Initial version
!  (1 ) Adapted from UNZIP_MET_FIELDS of "dao_read_mod.f" (bmy, 6/23/03)
!  (2 ) Directory information YYYY/MM or YYYYMM is now contained w/in 
!        GEOS_1_DIR, GEOS_S_DIR, GEOS_3_DIR, GEOS_4_DIR (bmy, 12/11/03)
!  (3 ) Now reference "directory_mod.f" and "unix_cmds_mod.f". Now prevent 
!        EXPAND_DATE from overwriting directory paths with Y/M/D tokens in 
!        them (bmy, 7/20/04)
!  (4 ) Now modified for GEOS-5 and GCAP met fields (swu, bmy, 5/25/05)
!  (5 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (6 ) Remove support for GEOS-1 and GEOS-STRAT met fields (bmy, 8/4/06)
!  21 Sep 2010 - R. Yantosca - Added ProTeX headers
!  28 Feb 2012 - R. Yantosca - Removed support for GEOS-3
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      CHARACTER(LEN=255) :: A3_STR,     GEOS_DIR
      CHARACTER(LEN=255) :: A3_FILE_GZ, A3_FILE
      CHARACTER(LEN=255) :: UNZIP_BG,   UNZIP_FG
      CHARACTER(LEN=255) :: REMOVE_ALL, REMOVE_DATE

      !=================================================================
      ! UNZIP_A3_FIELDS begins here!
      !=================================================================
      IF ( PRESENT( NYMD ) ) THEN

#if   defined( GEOS_4 )

         ! Strings for directory & filename
         GEOS_DIR = TRIM( Input_Opt%GEOS_4_DIR )
         A3_STR   = 'YYYYMMDD.a3.' // GET_RES_EXT()

#elif defined( GEOS_5 )

         ! Strings for directory & filename
         GEOS_DIR = TRIM( Input_Opt%GEOS_5_DIR )
         A3_STR   = 'YYYYMMDD.a3.' // GET_RES_EXT()

#elif defined( GCAP )

         ! Strings for directory & filename
         GEOS_DIR = TRIM( Input_Opt%GCAP_DIR )
         A3_STR   = 'YYYYMMDD.a3.' // GET_RES_EXT()

#endif

         ! Replace date tokens
         CALL EXPAND_DATE( GEOS_DIR, NYMD, 000000 )
         CALL EXPAND_DATE( A3_STR,   NYMD, 000000 )

         ! Location of zipped A-3 file in data dir 
         A3_FILE_GZ  = TRIM( Input_Opt%DATA_DIR   ) // 
     &                 TRIM( GEOS_DIR             ) // 
     &                 TRIM( A3_STR               ) // 
     &                 TRIM( Input_Opt%ZIP_SUFFIX )

         ! Location of unzipped A-3 file in temp dir 
         A3_FILE     = TRIM( Input_Opt%TEMP_DIR   ) // 
     &                 TRIM( A3_STR               )

         ! Remove A-3 files for this date from temp dir
         REMOVE_DATE = TRIM( Input_Opt%REMOVE_CMD ) // ' '    // 
     &                 TRIM( Input_Opt%TEMP_DIR   ) // 
     &                 TRIM( A3_STR               )   

         !==============================================================
         ! Define the foreground and background UNZIP commands
         !==============================================================

         ! Foreground unzip
         UNZIP_FG = TRIM( Input_Opt%UNZIP_CMD     ) // ' '    // 
     &              TRIM( A3_FILE_GZ              ) // 
     &              TRIM( Input_Opt%REDIRECT      ) // ' '    //
     &              TRIM( A3_FILE                 )  

         ! Background unzip
         UNZIP_BG  = TRIM( UNZIP_FG               ) // 
     &               TRIM( Input_Opt%BACKGROUND   )
      ENDIF

      !=================================================================
      ! Define command to remove all A-3 files from the TEMP dir
      !=================================================================
      REMOVE_ALL = TRIM( Input_Opt%REMOVE_CMD     ) // ' '    // 
     &             TRIM( Input_Opt%TEMP_DIR       ) // 
     &             TRIM( Input_Opt%WILD_CARD      ) // '.a3.' // 
     &             TRIM( Input_Opt%WILD_CARD      ) 

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

      END SUBROUTINE UNZIP_A3_FIELDS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: do_open_a3
!
! !DESCRIPTION: Function DO\_OPEN\_A3 returns TRUE if is time to open the A3 
!  met field file or FALSE otherwise.  This prevents us from opening a file 
!  which has already been opened. 
!\\
!\\
! !INTERFACE:
!
      FUNCTION DO_OPEN_A3( NYMD, NHMS, RESET ) RESULT( DO_OPEN )
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN)            :: NYMD    ! YYYY/MM/DD date
      INTEGER, INTENT(IN)            :: NHMS    ! hh:mm:ss time
      LOGICAL, INTENT(IN),  OPTIONAL :: RESET   ! Reset first-time flag?
! 
! !REVISION HISTORY: 
!  23 Jun 2003 - R. Yantosca - Initial version
!  (1 ) Now modified for GEOS-5 and GCAP met fields (swu, bmy, 5/25/05)
!  21 Sep 2010 - R. Yantosca - Added ProTeX headers
!  21 Sep 2010 - R. Yantosca - Bug fix: If we are using MEGAN (which reads many
!                              days of A3 data @ start of run), then reset the
!                              first-time flag.  This will prevent an error if
!                              if the start time is not 00 GMT.
!  28 Feb 2012 - R. Yantosca - Removed support for GEOS-3
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL       :: DO_OPEN
      LOGICAL, SAVE :: FIRST    = .TRUE.
      INTEGER, SAVE :: LASTNYMD = -1
      INTEGER, SAVE :: LASTNHMS = -1
      
      !=================================================================
      ! DO_OPEN_A3 begins here!
      !=================================================================

      ! Reset the FIRST flag if 
      IF ( PRESENT( RESET ) ) THEN
         IF ( RESET ) FIRST = .TRUE.
      ENDIF

      ! Initialize
      DO_OPEN = .FALSE.

      ! Return if we have already opened the file
      IF ( NYMD == LASTNYMD .and. NHMS == LASTNHMS ) THEN
         DO_OPEN = .FALSE. 
         GOTO 999
      ENDIF

      ! Open A-3 file if it's 01:30 GMT,  or on the first call
      IF ( NHMS == 013000 .or. FIRST ) THEN
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

      END FUNCTION DO_OPEN_A3
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: open_a3_fields
!
! !DESCRIPTION: Subroutine OPEN\_A3\_FIELDS opens the A3 met fields file 
!  for date NYMD and time NHMS.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE OPEN_A3_FIELDS( NYMD, NHMS, Input_Opt, RESET, IUNIT )
!
! !USES:
!
      USE BPCH2_MOD,          ONLY : GET_RES_EXT
      USE CMN_SIZE_MOD                          ! Size parameters
      USE ERROR_MOD,          ONLY : ERROR_STOP
      USE FILE_MOD,           ONLY : IOERROR,    FILE_EXISTS
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE TIME_MOD,           ONLY : EXPAND_DATE
!
! !INPUT PARAMETERS: 
!
      TYPE(OptInput), INTENT(IN)            :: Input_Opt  ! Input_Options
      INTEGER,        INTENT(IN)            :: NYMD       ! YYYY/MM/dd and 
      INTEGER,        INTENT(IN)            :: NHMS       !  hh:mm:ss of data
      LOGICAL,        INTENT(IN),  OPTIONAL :: RESET      ! Reset first flag?
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT), OPTIONAL :: IUNIT      ! Returns IU_A1 
! 
! !REVISION HISTORY: 
!  15 Jun 1998 - R. Yantosca - Initial version
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
!  21 Sep 2010 - R. Yantosca - Added ProTeX headers
!  21 Sep 2010 - R. Yantosca - Now pass RESET flag to DO_OPEN_A3
!  28 Feb 2012 - R. Yantosca - Removed support for GEOS-3
!  03 Aug 2012 - R. Yantosca - Move calls to findFreeLUN out of DEVEL block
!  03 Aug 2012 - R. Yantosca - Now use findFreeLUN to define IU_A3 locally
!  06 Aug 2012 - R. Yantosca - Add IUNIT to pass IU_A3 to calling routine
!  07 Aug 2012 - R. Yantosca - Now print LUN used to open file
!  11 Apr 2013 - R. Yantosca - Now pass Input_Opt object
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL            :: DO_RESET
      LOGICAL            :: IT_EXISTS, IS_OPEN
      INTEGER            :: IOS
      CHARACTER(LEN=2)   :: DUM
      CHARACTER(LEN=8)   :: IDENT
      CHARACTER(LEN=255) :: A3_FILE
      CHARACTER(LEN=255) :: GEOS_DIR
      CHARACTER(LEN=255) :: PATH

      !=================================================================
      ! OPEN_A3_FIELDS begins here!
      !=================================================================

      IF ( PRESENT( RESET ) ) THEN
         DO_RESET = RESET
      ELSE
         DO_RESET = .FALSE. 
      ENDIF

      ! Open A-3 fields at the proper time, or on the first call
      IF ( DO_OPEN_A3( NYMD, NHMS, DO_RESET ) ) THEN

#if   defined( GEOS_4 )

         ! Strings for directory & filename
         GEOS_DIR = TRIM( Input_Opt%GEOS_4_DIR )
         A3_FILE  = 'YYYYMMDD.a3.' // GET_RES_EXT()

#elif defined( GEOS_5 )

         ! Strings for directory & filename
         GEOS_DIR = TRIM( Input_Opt%GEOS_5_DIR )
         A3_FILE  = 'YYYYMMDD.a3.' // GET_RES_EXT()

#elif defined( GCAP )

         ! Strings for directory & filename
         GEOS_DIR = TRIM( Input_Opt%GCAP_DIR )
         A3_FILE  = 'YYYYMMDD.a3.' // GET_RES_EXT()

#endif

         ! Replace date tokens
         CALL EXPAND_DATE( A3_FILE,  NYMD, NHMS )
         CALL EXPAND_DATE( GEOS_DIR, NYMD, NHMS )

         ! If unzipping, open GEOS-4 file in TEMP dir
         ! If not unzipping, open GEOS-4 file in DATA dir
         IF ( Input_Opt%LUNZIP ) THEN
            PATH = TRIM( Input_Opt%TEMP_DIR ) // 
     &             TRIM( A3_FILE            )
         ELSE
            PATH = TRIM( Input_Opt%DATA_DIR ) // 
     &             TRIM( GEOS_DIR           ) // 
     &             TRIM( A3_FILE            )
         ENDIF

         ! Close previously opened A-3 file
	 INQUIRE( IU_A3, OPENED=IS_OPEN )
         IF( IS_OPEN ) CLOSE( IU_A3 )

         ! Make sure the file unit is valid before we open the file
         IF ( .not. FILE_EXISTS( IU_A3 ) ) THEN
            CALL ERROR_STOP( 'Could not find file!', 
     &                       'OPEN_A3_FIELDS (a3_read_mod.f)' )
         ENDIF

         ! Find a free file LUN
         IU_A3 = findFreeLUN( IU_A3 )

         ! Pass the value of IU_A3 to back to the calling routine
         IF ( PRESENT( IUNIT ) ) THEN
            IUNIT = IU_A3
         ENDIF

         ! Open the A3 file
         OPEN( UNIT   = IU_A3,         FILE   = TRIM( PATH ),
     &         STATUS = 'OLD',         ACCESS = 'SEQUENTIAL',  
     &         FORM   = 'UNFORMATTED', IOSTAT = IOS )
               
         IF ( IOS /= 0 ) THEN
            CALL IOERROR( IOS, IU_A3, 'open_a3_fields:1' )
         ENDIF

         ! Echo info
         WRITE( 6, 100 ) TRIM( PATH ), IU_A3
 100     FORMAT( '     - Opening: ', a,  ' on unit ', i4 )
         
         ! Skip past the ident string
         READ( IU_A3, IOSTAT=IOS ) IDENT

         IF ( IOS /= 0 ) THEN
            CALL IOERROR( IOS, IU_A3, 'open_a3_fields:2' )
         ENDIF

         ! The last 2 digits of the ident string
         ! is the # of fields contained in the file
         READ( IDENT(7:8), '(i2.2)' ) N_A3_FIELDS        
        
      ENDIF

      END SUBROUTINE OPEN_A3_FIELDS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_a3_fields
!
! !DESCRIPTION: Subroutine GET\_A3\_FIELDS is a wrapper for routine READ\_A3.  
!  GET\_A3\_FIELDS calls READ\_A3 properly for reading GEOS-3, GEOS-4, GEOS-5, 
!  or GCAP met data sets. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE GET_A3_FIELDS( NYMD, NHMS, State_Met )
!
! !USES:
!
      USE CMN_SIZE_MOD
      USE GIGC_State_Met_Mod, ONLY : MetState
#if   defined( APM )
      USE DAO_MOD, ONLY : PBL0a, PBL0b
#endif
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN) :: NYMD   ! YYYY/MM/DD
      INTEGER, INTENT(IN) :: NHMS   !  and hh:mm:ss of desired data
!
! !INPUT/OUTPUT PARAMETERS:
!
      TYPE(MetState), INTENT(INOUT) :: State_Met   ! Meteorology State object
! 
! !REVISION HISTORY: 
!  23 Jun 2003 - R. Yantosca - Initial version
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
!  (8 ) Now read EFLUX field for non-local PBL scheme (only GEOS5).
!        (ccc, 5/14/09)
!  (9 ) Now read FRLAND, FROCEAN, FRLANDIC and FRLAKE for methane
!        (kjw, 8/18/09)
!  21 Sep 2010 - R. Yantosca - Added ProTeX headers
!  28 Feb 2012 - R. Yantosca - Removed support for GEOS-3
!  09 Nov 2012 - M. Payer    - Copy all met fields to the State_Met derived type
!                              object
!  15 Nov 2012 - R. Yantosca - Now replace dao_mod.F arrays with State_Met
!  02 Jul 2013 - R. Yantosca - For GCAP, PBL is now State_Met%PBLH (similarly
!                              for State_Met%PRECCON and State_Met%PRECTOT)
!  13 Aug 2013 - R. Yantosca - For GEOS-4, PBL is now State_Met%PBLH (similarly
!                              for State_Met%PRECCON and State_Met%PRECTOT)
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER, SAVE :: LASTNYMD = -1, LASTNHMS = -1
#if   defined( APM )
      LOGICAL, SAVE :: FIRST=.TRUE.
#endif

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

#if   defined( GEOS_4 )

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
      CALL READ_A3( NYMD    = NYMD,       
     &              NHMS    = NHMS,     
     &              ALBEDO  = State_Met%ALBD,     
     &              CLDFRC  = State_Met%CLDFRC, 
     &              HFLUX   = State_Met%HFLUX,   
     &              GWETTOP = State_Met%GWETTOP, 
     &              PARDF   = State_Met%PARDF,   
     &              PARDR   = State_Met%PARDR,  
     &              PBL     = State_Met%PBLH,         
     &              PREACC  = State_Met%PRECTOT, 
     &              PRECON  = State_Met%PRECCON, 
     &              RADLWG  = State_Met%RADLWG,   
     &              RADSWG  = State_Met%RADSWG, 
     &              SNOW    = State_Met%SNOW,  
     &              TS      = State_Met%TS,           
     &              TSKIN   = State_Met%TSKIN,   
     &              U10M    = State_Met%U10M,       
     &              USTAR   = State_Met%USTAR,    
     &              V10M    = State_Met%V10M,     
     &              Z0      = State_Met%Z0 )

#elif defined( GEOS_5 )

      !=================================================================
      ! For GEOS-5, read the following fields:
      !
      !    ALBEDO, CLDFRC, EVAP,   GRN,    GWETROOT, GWETTOP, HFLUX, 
      !    LAI,    PARDF,  PARDR,  PBLH,   PREACC,   PRECON,  PRECSNO,
      !    RADLWG, RADSWG, SNODP,  SNOMAS, T2M,      TROPP,   TSKIN,   
      !    U10M,   USTAR,  V10M,   Z0,
      !    FRLAND, FRLANDIC,      FROCEAN,          FRLAKE
      !
      !=================================================================
      CALL READ_A3( NYMD     = NYMD,         
     &              NHMS     = NHMS,       
     &              ALBEDO   = State_Met%ALBD,   
     &              CLDFRC   = State_Met%CLDFRC,     
     &              EVAP     = State_Met%EVAP,       
     &              GRN      = State_Met%GRN,       
     &              GWETROOT = State_Met%GWETROOT, 
     &              GWETTOP  = State_Met%GWETTOP, 
     &              HFLUX    = State_Met%HFLUX,   
     &              LAI      = State_Met%LAI,           
     &              PARDF    = State_Met%PARDF,     
     &              PARDR    = State_Met%PARDR,       
     &              PBL      = State_Met%PBLH,           
     &              PREACC   = State_Met%PRECTOT,   
     &              PRECON   = State_Met%PRECCON, 
     &              PRECSNO  = State_Met%PRECSNO,   
     &              RADLWG   = State_Met%RADLWG,   
     &              RADSWG   = State_Met%RADSWG, 
     &              SNODP    = State_Met%SNODP,       
     &              SNOMAS   = State_Met%SNOMAS,   
     &              TROPP    = State_Met%TROPP,   
     &              TS       = State_Met%TS,             
     &              TSKIN    = State_Met%TSKIN,     
     &              U10M     = State_Met%U10M,     
     &              USTAR    = State_Met%USTAR,       
     &              V10M     = State_Met%V10M,       
     &              Z0       = State_Met%Z0, 
     &              EFLUX    = State_Met%EFLUX,       
     &              FRLAND   = State_Met%FRLAND,   
     &              FROCEAN  = State_Met%FROCEAN, 
     &              FRLAKE   = State_Met%FRLAKE,     
     &              FRLANDIC = State_Met%FRLANDIC )

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
      CALL READ_A3( NYMD     = NYMD,     
     &              NHMS     = NHMS,     
     &              ALBEDO   = State_Met%ALBD,   
     &              GWETTOP  = State_Met%GWETTOP,   
     &              MOLENGTH = State_Met%MOLENGTH,      
     &              OICE     = State_Met%OICE,     
     &              PARDF    = State_Met%PARDF,       
     &              PARDR    = State_Met%PARDR,
     &              PBL      = State_Met%PBLH,
     &              PREACC   = State_Met%PRECTOT,
     &              PRECON   = State_Met%PRECCON,
     &              RADSWG   = State_Met%RADSWG, 
     &              SNICE    = State_Met%SNICE,
     &              SNOW     = State_Met%SNOW, 
     &              TS       = State_Met%TS,
     &              U10M     = State_Met%U10M,         
     &              USTAR    = State_Met%USTAR,
     &              V10M     = State_Met%V10M  )

      ! Convert [m] to [mm]
      State_Met%SNOW = State_Met%SNOW*1000.D0

#endif

      ! Save NYMD, NHMS for next call
      LASTNYMD = NYMD
      LASTNHMS = NHMS

#if   defined( APM )
      !APM_GanLuo+ For PBL interpolation
      IF ( FIRST ) THEN
         PBL0b = State_Met%PBLH
         FIRST = .FALSE.
      ENDIF

      PBL0a = PBL0b
      PBL0b = State_Met%PBLH
#endif

      END SUBROUTINE GET_A3_FIELDS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_n_a3
!
! !DESCRIPTION: Function GET\_N\_A3 returns the number of A3 fields per 
!  met data set.
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_N_A3( NYMD ) RESULT( N_A3 )
!
! !USES:
!
      USE CMN_SIZE_MOD              ! Size parameters
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN) :: NYMD   ! YYYY/MM/DD date
!
! !RETURN VALUE:
!
      INTEGER             :: N_A3   ! Number of A3 fields in file
!
! !REMARKS:
! 
! 
! !REVISION HISTORY: 
!  23 Jun 2003 - R. Yantosca - Initial version
!  (1 ) GEOS-4/fvDAS now has 19 A-3 fields; we added LAI, RADLWG, SNOW.
!        (bmy, 12/9/03)
!  (2 ) Now modified for GEOS-5 and GCAP met fields (bmy, 5/25/05)
!  (3 ) Remove support for GEOS-1 and GEOS-STRAT met fields (bmy, 8/4/06)
!  (4 ) Increase # of fields for GCAP from 12 to 16 (swu, bmy, 10/4/06)
!  (5 ) Increase # of fields for GEOS-5 to 25 (bmy, 1/17/07)
!  (6 ) Increase # of fields for GEOS-5 to 26 (EFLUX) (ccc, 5/21/09)
!  21 Sep 2010 - R. Yantosca - Added ProTeX headers
!  28 Feb 2012 - R. Yantosca - Removed support for GEOS-3
!EOP
!------------------------------------------------------------------------------
!BOC
      !=================================================================
      ! GET_N_A3 begins here!
      !=================================================================
#if   defined( GEOS_4 )

      ! GEOS-4 has 19 A-3 fields
      N_A3 = 19

#elif defined( GEOS_5 )

      ! GEOS-5 has 33 A-3 fields (only 30 used in GEOS-Chem)
      N_A3 = 33

#elif defined( GCAP )
      
      ! GCAP has 16 fields
      N_A3 = 16

#endif

      END FUNCTION GET_N_A3
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: check_time
!
! !DESCRIPTION: Function CHECK\_TIME checks to see if the timestamp of the A3 
!  field just read from disk matches the current time.  If so, then it's 
!  time to return the A3 field to the calling program. 
!\\
!\\
! !INTERFACE:
!
      FUNCTION CHECK_TIME( XYMD, XHMS, NYMD, NHMS ) RESULT( ITS_TIME )
!
! !USES:
!
      USE CMN_SIZE_MOD
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN) :: XYMD       ! YYYY/MM/DD and hh:mm:ss
      INTEGER, INTENT(IN) :: XHMS       !  timestamp of A3 data in file
      INTEGER, INTENT(IN) :: NYMD       ! YYYY/MM/DD and hh:mm:ss
      INTEGER, INTENT(IN) :: NHMS       !  timestamp for desired data
! 
! !RETURN VALUE:
!
      LOGICAL             :: ITS_TIME   ! =T if XYMD & XHMS match NYMD & NHMS
!
! !REVISION HISTORY: 
!  23 Jun 2003 - R. Yantosca - Initial version
!  (1 ) Remove support for GEOS-1 and GEOS-STRAT met fields (bmy, 8/4/06)
!  21 Sep 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
      !=================================================================
      ! CHECK_TIME begins here!
      !=================================================================
      IF ( XYMD == NYMD .AND. XHMS == NHMS ) THEN
         ITS_TIME = .TRUE.
      ELSE
         ITS_TIME = .FALSE.
      ENDIF

      END FUNCTION CHECK_TIME
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_a3
!
! !DESCRIPTION: Subroutine READ\_A3 reads GEOS A-3 (3-hr avg) fields from disk.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE READ_A3( NYMD,    NHMS, 
     &                    ALBEDO,  CLDFRC, EVAP,   GRN,      GWETROOT,
     &                    GWETTOP, HFLUX,  LAI,    MOLENGTH, OICE,   
     &                    PARDF,   PARDR,  PBL,    PREACC,   PRECON,  
     &                    PRECSNO, RADLWG, RADSWG, RADSWT,   SNICE,    
     &                    SNODP,   SNOMAS, SNOW,   TROPP,    TS,       
     &                    TSKIN,   U10M,   USTAR,  V10M,     Z0,
     &                    EFLUX,   FRLAND, FRLAKE, FROCEAN,  FRLANDIC )
!
! !USES:
!
      USE DIAG_MOD,     ONLY : AD67
      USE FILE_MOD,     ONLY : IOERROR
      USE TIME_MOD,     ONLY : SET_CT_A3
      USE TIME_MOD,     ONLY : TIMESTAMP_STRING
      USE TRANSFER_MOD, ONLY : TRANSFER_2D
      USE TRANSFER_MOD, ONLY : TRANSFER_TO_1D

      USE CMN_SIZE_MOD                         ! Size parameters
      USE CMN_DIAG_MOD                         ! ND67
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN)            :: NYMD   ! YYYYMMDD
      INTEGER, INTENT(IN)            :: NHMS   !  and hhmmss of desired data
!
! !OUTPUT PARAMETERS:
!
      REAL*8,  INTENT(OUT), OPTIONAL :: ALBEDO  (IIPAR,JJPAR) 
      REAL*8,  INTENT(OUT), OPTIONAL :: CLDFRC  (IIPAR,JJPAR)
      REAL*8,  INTENT(OUT), OPTIONAL :: EVAP    (IIPAR,JJPAR) 
      REAL*8,  INTENT(OUT), OPTIONAL :: GRN     (IIPAR,JJPAR)
      REAL*8,  INTENT(OUT), OPTIONAL :: GWETROOT(IIPAR,JJPAR) 
      REAL*8,  INTENT(OUT), OPTIONAL :: GWETTOP (IIPAR,JJPAR) 
      REAL*8,  INTENT(OUT), OPTIONAL :: HFLUX   (IIPAR,JJPAR)
      REAL*8,  INTENT(OUT), OPTIONAL :: LAI     (IIPAR,JJPAR)
      REAL*8,  INTENT(OUT), OPTIONAL :: MOLENGTH(IIPAR,JJPAR)
      REAL*8,  INTENT(OUT), OPTIONAL :: OICE    (IIPAR,JJPAR)
      REAL*8,  INTENT(OUT), OPTIONAL :: PARDF   (IIPAR,JJPAR) 
      REAL*8,  INTENT(OUT), OPTIONAL :: PARDR   (IIPAR,JJPAR) 
      REAL*8,  INTENT(OUT), OPTIONAL :: PBL     (IIPAR,JJPAR) 
      REAL*8,  INTENT(OUT), OPTIONAL :: PREACC  (IIPAR,JJPAR) 
      REAL*8,  INTENT(OUT), OPTIONAL :: PRECON  (IIPAR,JJPAR)
      REAL*8,  INTENT(OUT), OPTIONAL :: PRECSNO (IIPAR,JJPAR)
      REAL*8,  INTENT(OUT), OPTIONAL :: RADLWG  (IIPAR,JJPAR)
      REAL*8,  INTENT(OUT), OPTIONAL :: RADSWG  (IIPAR,JJPAR)
      REAL*8,  INTENT(OUT), OPTIONAL :: RADSWT  (IIPAR,JJPAR)
      REAL*8,  INTENT(OUT), OPTIONAL :: SNICE   (IIPAR,JJPAR)
      REAL*8,  INTENT(OUT), OPTIONAL :: SNODP   (IIPAR,JJPAR)
      REAL*8,  INTENT(OUT), OPTIONAL :: SNOMAS  (IIPAR,JJPAR)
      REAL*8,  INTENT(OUT), OPTIONAL :: SNOW    (IIPAR,JJPAR)
      REAL*8,  INTENT(OUT), OPTIONAL :: TROPP   (IIPAR,JJPAR) 
      REAL*8,  INTENT(OUT), OPTIONAL :: TS      (IIPAR,JJPAR) 
      REAL*8,  INTENT(OUT), OPTIONAL :: TSKIN   (IIPAR,JJPAR) 
      REAL*8,  INTENT(OUT), OPTIONAL :: U10M    (IIPAR,JJPAR) 
      REAL*8,  INTENT(OUT), OPTIONAL :: USTAR   (IIPAR,JJPAR) 
      REAL*8,  INTENT(OUT), OPTIONAL :: V10M    (IIPAR,JJPAR)
      REAL*8,  INTENT(OUT), OPTIONAL :: Z0      (IIPAR,JJPAR)
      REAL*8,  INTENT(OUT), OPTIONAL :: EFLUX   (IIPAR,JJPAR)
      REAL*8,  INTENT(OUT), OPTIONAL :: FRLAND  (IIPAR,JJPAR)
      REAL*8,  INTENT(OUT), OPTIONAL :: FRLAKE  (IIPAR,JJPAR)
      REAL*8,  INTENT(OUT), OPTIONAL :: FROCEAN (IIPAR,JJPAR)
      REAL*8,  INTENT(OUT), OPTIONAL :: FRLANDIC(IIPAR,JJPAR)
!
! !REMARKS:
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
!  (31) EFLUX    : (2-D) GMAO latent heat flux                  [W/m2]
!  (32) FRLAND   : (2-D) GMAO fraction of land                  [unitless]
!  (33) FROCEAN  : (2-D) GMAO fraction of ocean                 [unitless]
!  (34) FRLANDIC : (2-D) GMAO fraction of land ice              [unitless]
!  (35) FRLAKE   : (2-D) GMAO fraction of lake water            [unitless]
! 
! !REVISION HISTORY: 
!  08 May 1998 - R. Yantosca - Initial version
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
!  (9 ) Now read EFLUX for non-local PBL scheme for GEOS5 (ccc, 5/14/09)
!  (10) Now read FRLAND, FROCEAN, FRLANDIC, FRLAKE for methane (kjw, 8/18/09)
!  (11) Remove reference to IN_CLOUD_OD (bmy, 10/15/09)
!  21 Sep 2010 - R. Yantosca - Added ProTeX headers
!  08 Jun 2012 - S. Philip   - Correction for GEOS_5 boundary layer height
!  03 Aug 2012 - R. Yantosca - Now use locally-defined IU_A3 file LUN
!  07 Aug 2012 - R. Yantosca - Now print LUN used to open file
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER           :: I, IJLOOP, IOS, J, N_A3, NFOUND 
      REAL*4            :: Q2(IIPAR,JJPAR)
      CHARACTER(LEN=8)  :: NAME
      CHARACTER(LEN=16) :: STAMP
      INTEGER           :: XYMD, XHMS

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
            ! Get EFLUX (Lin, 03/31/09)
            ! EFLUX: Latent heat flux
            !--------------------------------            
            CASE ( 'EFLUX' )
               READ( IU_A3, IOSTAT=IOS ) XYMD, XHMS, Q2
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A3, 'read_a3:32' )
             
               IF ( CHECK_TIME( XYMD, XHMS, NYMD, NHMS ) ) THEN
                  IF ( PRESENT( EFLUX ) ) CALL TRANSFER_2D( Q2, EFLUX )
                  NFOUND = NFOUND + 1
               ENDIF

            !--------------------------------
            ! TPW: just skip over this
            !--------------------------------
            CASE ( 'TPW' )
               READ( IU_A3, IOSTAT=IOS ) XYMD, XHMS, Q2
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A3, 'read_a3:33' )

               IF ( CHECK_TIME( XYMD, XHMS, NYMD, NHMS ) ) THEN
                  NFOUND = NFOUND + 1
               ENDIF

            !--------------------------------
            ! CLDTMP: just skip over this
            !--------------------------------
            CASE ( 'CLDTMP' )
               READ( IU_A3, IOSTAT=IOS ) XYMD, XHMS, Q2
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A3, 'read_a3:34' )
             
               IF ( CHECK_TIME( XYMD, XHMS, NYMD, NHMS ) ) THEN
                  NFOUND = NFOUND + 1
               ENDIF

            !--------------------------------            
            ! FRLAND: land fraction
            !--------------------------------            
            CASE ( 'FRLAND' )
               READ( IU_A3, IOSTAT=IOS ) XYMD, XHMS, Q2
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A3, 'read_a3:35')
             
               IF ( CHECK_TIME( XYMD, XHMS, NYMD, NHMS ) ) THEN
                  IF ( PRESENT( FRLAND ) ) CALL TRANSFER_2D( Q2, FRLAND) 
                  NFOUND = NFOUND + 1
               ENDIF         

            !--------------------------------            
            ! FRLAKE: lake fraction
            !--------------------------------            
            CASE ( 'FRLAKE' )
               READ( IU_A3, IOSTAT=IOS ) XYMD, XHMS, Q2
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A3, 'read_a3:36')
             
               IF ( CHECK_TIME( XYMD, XHMS, NYMD, NHMS ) ) THEN
                  IF ( PRESENT( FRLAKE ) ) CALL TRANSFER_2D( Q2, FRLAKE) 
                  NFOUND = NFOUND + 1
               ENDIF 

            !--------------------------------            
            ! FROCEAN: ocean fraction
            !--------------------------------            
            CASE ( 'FROCEAN' )
               READ( IU_A3, IOSTAT=IOS ) XYMD, XHMS, Q2
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A3, 'read_a3:37')
             
               IF ( CHECK_TIME( XYMD, XHMS, NYMD, NHMS ) ) THEN
                  IF ( PRESENT( FROCEAN ) ) THEN
                     CALL TRANSFER_2D( Q2, FROCEAN ) 
                     NFOUND = NFOUND + 1
                  ENDIF
               ENDIF 

            !--------------------------------            
            ! FRLANDIC: land ice fraction
            !--------------------------------            
            CASE ( 'FRLANDIC' )
               READ( IU_A3, IOSTAT=IOS ) XYMD, XHMS, Q2
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A3, 'read_a3:38')
             
               IF ( CHECK_TIME( XYMD, XHMS, NYMD, NHMS ) ) THEN
                  IF ( PRESENT( FRLANDIC ) ) THEN
                     CALL TRANSFER_2D( Q2, FRLANDIC ) 
                     NFOUND = NFOUND + 1
                  ENDIF
               ENDIF 

            !--------------------------------
            ! Extra GEOS-5 fields
            ! Skip over for now, add later
            ! (bmy, 10/7/08)
            !--------------------------------
            CASE ( 'PRECANV', 'LWTUP',  'QV2M'  )
               READ( IU_A3, IOSTAT=IOS ) XYMD, XHMS, Q2
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A3, 'read_a3:39' )
             
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
            WRITE( 6, 210 ) NFOUND, STAMP, IU_A3
 210        FORMAT( '     - Found all ', i3, ' A-3 met fields for ', 
     &                   a, ' on unit ', i4 )
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

      ! Correction for GEOS_5 PBLH:
      ! PBLH must be greater than a minimum mechanical mixing depth,
      ! defined as 700*friction velocity
      ! (Koracin and Berkowicz, 1988; Lin and McElroy, 2010)
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         PBL(I,J) = max( PBL(I,J), 700d0*USTAR(I,J) )
      ENDDO
      ENDDO

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
      ! (23) EFLUX  : latent heat flux from surface      [W/m2]
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
         IF ( PRESENT( EFLUX   ) ) AD67(:,:,23) = AD67(:,:,23) + EFLUX
      ENDIF
         
      ! Increment KDA3FLDS -- this is the # of times READ_A3 is called
      CALL SET_CT_A3( INCREMENT=.TRUE. )

      END SUBROUTINE READ_A3
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: archive_nd67_1d
!
! !DESCRIPTION: Subroutine ARCHIVE\_ND67\_1D saves 1-D arrays for the ND67 
!  diagnostic.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE ARCHIVE_ND67_1D( FIELD, N )
!
! !USES:
!
      USE DIAG_MOD, ONLY : AD67             ! ND67 diagnostic array

      USE CMN_SIZE_MOD                      ! Size parameters
!
! !INPUT PARAMETERS: 
!
      REAL*8,  INTENT(IN) :: FIELD(MAXIJ)   ! Array to be archived in ND67
      INTEGER, INTENT(IN) :: N              ! ND67 index in which to store data
! 
! !REVISION HISTORY: 
!  23 Jun 2003 - R. Yantosca - Initial version
!  21 Sep 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER :: I, IJLOOP, J

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

      END SUBROUTINE ARCHIVE_ND67_1D
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: a3_check
!
! !DESCRIPTION: Subroutine A3\_CHECK prints an error message if not all 
!  of the A-3 met fields are found.  The run is also terminated. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE A3_CHECK( NFOUND, N_A3 )
!
! !USES:
!
      USE ERROR_MOD, ONLY : GEOS_CHEM_STOP
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN) :: NFOUND   ! Number of A3 fields found in file
      INTEGER, INTENT(IN) :: N_A3     ! Expected number of A3 fields
! 
! !REVISION HISTORY: 
!  27 Oct 2000 - R. Yantosca - Initial version
!  (1 ) Adapted from DAO_CHECK from "dao_read_mod.f" (bmy, 6/23/03)
!  21 Sep 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
      !=================================================================
      ! A3_CHECK begins here!
      !=================================================================
      IF ( NFOUND /= N_A3 ) THEN
         WRITE( 6, '(a)' ) REPEAT( '=', 79 )
         WRITE( 6, '(a)' ) 'ERROR -- not enough A-3 fields found!'      

         WRITE( 6, 120   ) N_A3, NFOUND
 120     FORMAT( 'There are ', i2, ' fields but only ', i2 ,
     &           ' were found!' )

         WRITE( 6, '(a)' ) '### STOP in A3_CHECK (a3_read_mod.f)'
         WRITE( 6, '(a)' ) REPEAT( '=', 79 )

         ! Deallocate arrays and stop (bmy, 10/15/02)
         CALL GEOS_CHEM_STOP
      ENDIF

      END SUBROUTINE A3_CHECK
!EOC
      END MODULE A3_READ_MOD
