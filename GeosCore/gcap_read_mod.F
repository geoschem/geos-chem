!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: gcap_read_mod.f
!
! !DESCRIPTION: Module GCAP\_READ\_MOD contains file unit numbers, as well 
!  as file I/O routines for GEOS-Chem.  FILE\_MOD keeps all of the I/O unit 
!  numbers in a single location for convenient access.
!\\
!\\
! !INTERFACE: 
!
      MODULE GCAP_READ_MOD
!
! !USES:
! 
      USE inquireMod, ONLY : findFreeLUN

      IMPLICIT NONE
      PRIVATE 
!
! !PRIVATE MEMBER FUNCTIONS:
!
      PRIVATE :: CHECK_TIME
      PRIVATE :: READ_GCAP
      PRIVATE :: GCAP_CHECK
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC :: GET_GCAP_FIELDS
      PUBLIC :: OPEN_GCAP_FIELDS
      PUBLIC :: UNZIP_GCAP_FIELDS
!
! !REVISION HISTORY:
!  (1 ) Adapted from the obsolete "phis_read_mod.f" (bmy, 2/1/06)
!  03 Aug 2012 - R. Yantosca - Move calls to findFreeLUN out of DEVEL block
!  06 Aug 2012 - R. Yantosca - Now make IU_PH a local variable
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE TYPES:
!
      ! LUN for the GCAP "phis" file
      INTEGER :: IU_PH   

      CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: unzip_gcap_fields
!
! !DESCRIPTION: Subroutine UNZIP\_GCAP\_FIELDS invokes a FORTRAN system call
!  to uncompress GCAP PHIS met field files and store the uncompressed data in a 
!  temporary directory, where GEOS-CHEM can read them.  The original data 
!  files are not disturbed. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE UNZIP_GCAP_FIELDS( Input_Opt, OPTION )
!
! !USES:
!
      USE CMN_SIZE_MOD
      USE BPCH2_MOD,          ONLY : GET_RES_EXT
      USE ERROR_MOD,          ONLY : ERROR_STOP
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE TIME_MOD,           ONLY : EXPAND_DATE
!
! !INPUT PARAMETERS: 
!
      TYPE(OptInput),    INTENT(IN) :: Input_Opt   ! Input Options object
      CHARACTER(LEN=*),  INTENT(IN) :: OPTION      ! Unzipping option
! 
! !REVISION HISTORY: 
!  15 Jun 1998 - R. Yantosca - Initial version
!  (1 ) Adapted from UNZIP_MET_FIELDS of "dao_read_mod.f" (bmy, 6/16/03)
!  (2 ) Directory information YYYY/MM or YYYYMM is now contained w/in 
!        GEOS_1_DIR, GEOS_S_DIR, GEOS_3_DIR, GEOS_4_DIR (bmy, 12/11/03)
!  (3 ) Now reference "directory_mod.f" and "unix_cmds_mod.f". Now prevent 
!        EXPAND_DATE from overwriting directory paths with Y/M/D tokens in 
!        them (bmy, 7/20/04)
!  (4 ) Now modified for GEOS-5 and GCAP met fields (swu, bmy, 5/25/05)
!  06 Aug 2012 - R. Yantosca - Added ProTeX headers
!  11 Apr 2013 - R. Yantosca - Now pass fields via Input_Opt
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER            :: NYMD
      CHARACTER(LEN=255) :: GEOS_DIR,     PHIS_STR
      CHARACTER(LEN=255) :: PHIS_FILE_GZ, PHIS_FILE
      CHARACTER(LEN=255) :: UNZIP_BG,     UNZIP_FG
      CHARACTER(LEN=255) :: REMOVE_ALL,   REMOVE_DATE

      !=================================================================
      ! UNZIP_GCAP_FIELD begins here!
      !=================================================================

      ! Date for PHIS field
      NYMD     = 20000101

      ! Strings for directory & filename
      GEOS_DIR = TRIM( Input_Opt%GCAP_DIR )
      PHIS_STR = 'YYYYMMDD.phis.' // GET_RES_EXT() 

      ! Replace date tokens
      CALL EXPAND_DATE( GEOS_DIR, NYMD, 000000 )
      CALL EXPAND_DATE( PHIS_STR, NYMD, 000000 )

      ! Location of zipped A-3 file in data dir
      PHIS_FILE_GZ = TRIM( Input_Opt%DATA_DIR   ) // 
     &               TRIM( GEOS_DIR             ) //
     &               TRIM( PHIS_STR             ) // 
     &               TRIM( Input_Opt%ZIP_SUFFIX )

      ! Location of unzipped A-3 file in temp dir
      PHIS_FILE    = TRIM( Input_Opt%TEMP_DIR   ) // 
     &               TRIM( PHIS_STR             )
         
      ! Remove A-3 files for this date from temp dir 
      REMOVE_DATE  = TRIM( Input_Opt%REMOVE_CMD ) // ' '    // 
     &               TRIM( Input_Opt%TEMP_DIR   ) // 
     &               TRIM( PHIS_STR             )

      !=================================================================
      ! Define the foreground and background UNZIP commands
      !=================================================================

      ! Foreground unzip
      UNZIP_FG = TRIM( Input_Opt%UNZIP_CMD      ) // ' '    // 
     &           TRIM( PHIS_FILE_GZ             ) // 
     &           TRIM( Input_Opt%REDIRECT       ) // ' '    // 
     &           TRIM( PHIS_FILE                )  

      ! Background unzip
      UNZIP_BG  = TRIM( UNZIP_FG                ) // 
     &            TRIM( Input_Opt%BACKGROUND    )

      !=================================================================
      ! Define command to remove all PHIS files from the TEMP dir
      !=================================================================
      REMOVE_ALL = TRIM( Input_Opt%REMOVE_CMD   ) // ' '    // 
     &             TRIM( Input_Opt%TEMP_DIR     ) // 
     &             TRIM( Input_Opt%WILD_CARD    ) //'.phis.'// 
     &             TRIM( Input_Opt%WILD_CARD    ) 

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

      END SUBROUTINE UNZIP_GCAP_FIELDS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: open_gcap_fields
!
! !DESCRIPTION: Subroutine OPEN\_GCAP\_FIELDS opens the PHIS and LWI met 
!  fields file.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE OPEN_GCAP_FIELDS( Input_Opt )
!
! !USES:
!
      USE BPCH2_MOD,          ONLY : GET_RES_EXT
      USE CMN_SIZE_MOD
      USE ERROR_MOD,          ONLY : ERROR_STOP
      USE FILE_MOD,           ONLY : IOERROR, FILE_EXISTS
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE TIME_MOD,           ONLY : EXPAND_DATE
!
! !INPUT PARAMETERS:
!
      TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
!
! !REVISION HISTORY: 
!  01 Feb 2006 - S. Wu        - Initial version
!  (1 ) Adapted from OPEN_MET_FIELDS of "dao_read_mod.f" (bmy, 6/13/03)
!  (2 ) Now opens either zipped or unzipped files (bmy, 12/11/03)
!  (3 ) Now skips past the GEOS-4 ident string (bmy, 12/12/04)
!  (4 ) Now references "directory_mod.f" instead of CMN_SETUP.  Also now
!        references LUNZIP from "logical_mod.f".  Also now prevents EXPAND_DATE
!        from overwriting Y/M/D tokens in directory paths. (bmy, 7/20/04)
!  (5 ) Now use FILE_EXISTS from "file_mod.f" to determine if file unit IU_PH 
!        refers to a valid file on disk (bmy, 3/23/05)
!  (6 ) Now modified for GEOS-5 and GCAP met fields (swu, bmy, 5/25/05)
!  06 Aug 2012 - R. Yantosca - Move calls to findFreeLUN out of DEVEL block
!  06 Aug 2012 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL, SAVE      :: FIRST = .TRUE.
      LOGICAL            :: IT_EXISTS, IS_OPEN
      INTEGER            :: NYMD, NHMS
      INTEGER            :: IOS,  IUNIT
      CHARACTER(LEN=8)   :: IDENT
      CHARACTER(LEN=255) :: GEOS_DIR
      CHARACTER(LEN=255) :: PHIS_FILE
      CHARACTER(LEN=255) :: PATH

      !=================================================================
      ! OPEN_PHIS_FIELDS begins here!
      !=================================================================
      
      ! Define date and hour
      NYMD = 20000101
      NHMS = 000000

      ! Open the A-3 file 0 GMT of each day, or on the first call
      IF ( NHMS == 000000 .or. FIRST ) THEN

         ! Strings for directory & filename
         GEOS_DIR  = TRIM( Input_Opt%GCAP_DIR )
         PHIS_FILE = 'YYYYMMDD.phis.' // GET_RES_EXT() 

         ! Replace date tokens
         CALL EXPAND_DATE( GEOS_DIR,  NYMD, NHMS )
         CALL EXPAND_DATE( PHIS_FILE, NYMD, NHMS )

         ! If unzipping, open GEOS-1 file in TEMP dir
         ! If not unzipping, open GEOS-1 file in DATA dir
         IF ( Input_Opt%LUNZIP ) THEN
            PATH = TRIM( Input_Opt%TEMP_DIR ) // 
     &             TRIM( PHIS_FILE          )
         ELSE
            PATH = TRIM( Input_Opt%DATA_DIR ) // 
     &             TRIM( GEOS_DIR           ) // 
     &             TRIM( PHIS_FILE          )
         ENDIF

         ! Close previously opened A-3 file
	 INQUIRE( IU_PH, OPENED=IS_OPEN )
         IF ( IS_OPEN ) CLOSE( IU_PH )

         ! Make sure the file unit is valid before we open the file
         IF ( .not. FILE_EXISTS( IU_PH ) ) THEN
            CALL ERROR_STOP( 'Could not find file!', 
     &                       'OPEN_PHIS_FIELD (phis_read_mod.f)' )
         ENDIF

         ! Open the file
	 IU_PH = findFreeLUN()

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

      END SUBROUTINE OPEN_GCAP_FIELDS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_gcap_fields
!
! !DESCRIPTION: Subroutine GET\_GCAP\_FIELDS calls READ\_GCAP to read GCAP
!  fields from disk at the start of a GEOS-Chem run. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE GET_GCAP_FIELDS( Input_Opt, State_Met )
!
! !USES:
!
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GIGC_State_Met_Mod, ONLY : MetState
!
! !INPUT PARAMETERS:
!
      TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
!
! !INPUT/OUTPUT PARAMETERS:
!
      TYPE(MetState), INTENT(INOUT) :: State_Met   ! Meteorology State object
! 
! !REVISION HISTORY: 
!  01 Feb 2006 - S. Wu       - Initial version
!  (1 ) Now also read LWI_GISS for GCAP met fields (swu, bmy, 5/25/05)
!  06 Aug 2012 - R. Yantosca - Added ProTeX headers
!  09 Nov 2012 - M. Payer    - Copy all met fields to the State_Met derived type
!                              object
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER :: NYMD, NHMS 

      !=================================================================
      ! GET_PHIS_FIELD begins here!
      !=================================================================
      
      ! Date and time
      NYMD = 20000101
      NHMS = 000000

      ! For GCAP met fields: read PHIS and LWI_GISS
      CALL READ_GCAP( NYMD = NYMD, 
     &                NHMS = NHMS, 
     &                PHIS = State_Met%PHIS, 
     &                LWI  = State_Met%LWI_GISS )

      ! Return to calling program
      END SUBROUTINE GET_GCAP_FIELDS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: check_time
!
! !DESCRIPTION: Function CHECK\_TIME checks to see if the timestamp of the 
!  GCAP field just read from disk matches the current time.  If so, then it's 
!  time to return the GCAP field to the calling program. (bmy, 6/16/03)
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
      INTEGER, INTENT(IN) :: XYMD       ! Date stamp in file (YYYYMMDD)
      INTEGER, INTENT(IN) :: XHMS       ! Time stamp in file (hhmmss)
      INTEGER, INTENT(IN) :: NYMD       ! Current model date (YYYYMMDD)
      INTEGER, INTENT(IN) :: NHMS       ! Current model time (hhmmss) 
!
! !RETURN VALUE:
!
      ! Function value
      LOGICAL             :: ITS_TIME   ! = T is time to return fields
! 
! !REVISION HISTORY: 
!  16 Jun 2003 - R. Yantosca - Initial version
!  06 Aug 2012 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
      !=================================================================
      ! GEOS-3, GEOS-4: XYMD and XHMS are integers
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
! !IROUTINE: read_gcap
!
! !DESCRIPTION: Subroutine READ\_GCAP reads the PHIS (surface geopotential 
!  heights) field from disk.  PHIS is an I-6 field, but is time-independent.  
!  Thus READ\_GCAP only needs to be called once at the beginning of the model 
!  run.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE READ_GCAP( NYMD, NHMS, PHIS, LWI )
!
! !USES:
!
      USE DIAG_MOD,     ONLY : AD67
      USE FILE_MOD,     ONLY : IOERROR
      USE TIME_MOD,     ONLY : TIMESTAMP_STRING
      USE TRANSFER_MOD, ONLY : TRANSFER_2D

      USE CMN_SIZE_MOD                                     ! Size parameters
      USE CMN_DIAG_MOD                                     ! ND67
      USE CMN_GCTM_MOD                                     ! g0
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN)            :: NYMD               ! YYYYMMDD date
      INTEGER, INTENT(IN)            :: NHMS               ! hhmmss time
!
! !OUTPUT PARAMETERS:
!
      REAL*8,  INTENT(OUT)           :: PHIS(IIPAR,JJPAR)  ! PHIS [m2/s2]
      REAL*8,  INTENT(OUT), OPTIONAL :: LWI(IIPAR,JJPAR)   ! LWI flags
! 
! !REVISION HISTORY: 
!  01 Feb 2006 - S. Wu       - Initial version
!  (1 ) Adapted from READ_PHIS from "dao_read_mod.f" (bmy, 6/16/03)
!  (2 ) Now use function TIMESTAMP_STRING from "time_mod.f" for formatted 
!        date/time output. (bmy, 10/28/03)
!  (3 ) Now also read LWI_GISS for GCAP met fields.  Added optional variable
!        LWI to the arg list. (swu, bmy, 5/25/05)
!  06 Aug 2012 - R. Yantosca - Added ProTeX headers
!  07 Aug 2012 - R. Yantosca - Now print LUN used to open file
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER            :: NFOUND, IOS
      INTEGER            :: XYMD,   XHMS 
      REAL*4             :: Q2(IIPAR,JJPAR)
      CHARACTER(LEN=8)   :: NAME
      CHARACTER(LEN=16)  :: STAMP

      ! Number of fields in the file
      INTEGER, PARAMETER :: N_PHIS = 2

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
            WRITE( 6, 200 ) STAMP, IU_PH 
 200        FORMAT( '     - Found GCAP met fields for ', 
     &                    a, ' on unit ', i4 )
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

      END SUBROUTINE READ_GCAP
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: gcap_check
!
! !DESCRIPTION: Subroutine GCAP\_CHECK prints an error message if not all of 
!  the GCAP met fields are found.  The run is also terminated. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE GCAP_CHECK( NFOUND, N_PHIS )
!
! !USES:
!
      USE ERROR_MOD, ONLY : GEOS_CHEM_STOP

!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN) :: NFOUND   ! # of fields read from disk
      INTEGER, INTENT(IN) :: N_PHIS   ! # of of fields expected to be found
!
! !REMARKS:
! 
! 
! !REVISION HISTORY: 
!  15 Jun 1998 - R. Yantosca - Initial version
!  (1 ) Adapted from DAO_CHECK from "dao_read_mod.f" (bmy, 6/16/03)
!  06 Aug 2012 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
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

      END SUBROUTINE GCAP_CHECK
!EOC
      END MODULE GCAP_READ_MOD
