! $Id: bpch2_mod.f,v 1.3 2003/12/05 21:13:57 bmy Exp $
      MODULE BPCH2_MOD
!
!******************************************************************************
!  Module BPCH2_MOD contains the routines used to read data from and write
!  data to binary punch (BPCH) file format (v. 2.0). (bmy, 6/28/00, 11/3/03)
!
!  Module Routines:
!  ============================================================================
!  (1 ) OPEN_BPCH2_FOR_READ : Opens binary punch file for input
!  (2 ) OPEN_BPCH_FOR_WRITE : Opens binary punch file for output
!  (3 ) BPCH2_HDR           : writes "top-of-file" header to BPCH file
!  (4 ) BPCH2               : writes a data block to a binary punch file
!  (5 ) READ_BPCH2          : reads a data block from a binary punch file
!  (6 ) GET_MODELNAME       : returns MODELNAME for the given met field
!  (7 ) GET_NAME_EXT        : returns file extension string for model name
!  (8 ) GET_RES_EXT         : returns file extension string for model res.
!  (9 ) GET_TAU0_6A         : computes TAU0 from MONTH, DAY, YEAR (, H, M, SEC)
!
!  Module Interfaces
!  ============================================================================
!  (1 ) GET_TAU0            : Overloads GET_TAU0_6A
!
!  GEOS-CHEM modules referenced by bpch2_mod.f
!  ============================================================================
!  (1 ) error_mod.f         : Contains NaN and other error-check routines
!  (2 ) file_mod.f          : Contains file unit numbers and error checks
!  (3 ) julday_mod.f        : Contains astronomical Julian date routines
!
!  NOTES:
!  (1 ) Added routine GET_TAU0 (bmy, 7/20/00)
!  (2 ) Added years 1985-2001 for routine GET_TAU0 (bmy, 8/1/00)
!  (3 ) Use IOS /= 0 criterion to also check for EOF (bmy, 9/12/00)
!  (4 ) Removed obsolete code in "read_bpch2.f" (bmy, 12/18/00)
!  (5 ) Correct error for 1991 TAU values in GET_TAU0 (bnd, bmy, 1/4/01)
!  (6 ) BPCH2_MOD is now independent of any GEOS-CHEM size parameters.
!        (bmy, 4/18/01)
!  (7 ) Now have 2 versions of "GET_TAU0" overloaded by an interface.  The
!        original version takes 2 arguments (MONTH, YEAR).  The new version
!        takes 3 arguments (MONTH, DAY, YEAR). (bmy, 8/22/01)
!  (8 ) Updated comments (bmy, 9/4/01)
!  (9 ) Renamed GET_TAU0_3A to GET_TAU0_6A, and updated the GET_TAU0 
!        interface.  Also updated comments (bmy, 9/26/01)
!  (10) Now use special model name for GEOS-3 w/ 30 layers (bmy, 10/9/01)
!  (11) Minor bug fix in GET_TAU0_2A.  Also deleted obsolete code from 9/01.
!        (bmy, 11/15/01)
!  (12) Moved routines JULDAY, MINT, CALDATE to "julian_mod.f".  Now 
!        references routine JULDAY from "julday_mod.f".  Also added code
!        for GEOS-4/fvDAS model type. (bmy, 11/20/01)
!  (23) Now divide module header into MODULE PRIVATE, MODULE VARIABLES, and
!        MODULE ROUTINES sections.  Also add MODULE INTERFACES section,
!        since we have an interface here. (bmy, 5/28/02)
!  (24) Added OPEN_BPCH2_FOR_READ and OPEN_BPCH2_FOR_WRITE.  Also now 
!        reference IU_FILE and IOERROR from "file_mod.f". (bmy, 7/30/02)
!  (25) Now references "error_mod.f".  Also obsoleted routine GET_TAU0_2A.
!        (bmy, 10/15/02)
!  (26) Made modification in READ_BPCH2 for 1x1 nested grids (bmy, 3/11/03)
!  (27) Modifications for GEOS-4, 30-layer grid (bmy, 11/3/03)
!******************************************************************************
!
      IMPLICIT NONE

      !=================================================================
      ! MODULE PRIVATE DECLARATIONS -- keep certain internal variables 
      ! and routines from being seen outside "bpch2_mod.f"
      !=================================================================

      ! PRIVATE module routines
      PRIVATE :: GET_TAU0_6A

      !=================================================================
      ! MODULE INTERFACES -- "bind" two or more routines with different
      ! argument types or # of arguments under one unique name
      !================================================================= 
      INTERFACE GET_TAU0
         MODULE PROCEDURE GET_TAU0_6A
      END INTERFACE

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement 
      !=================================================================
      CONTAINS

!------------------------------------------------------------------------------

      SUBROUTINE OPEN_BPCH2_FOR_READ( IUNIT, FILENAME, TITLE )
!
!******************************************************************************
!  Subroutine OPEN_BPCH2_FOR_READ opens a binary punch file (version 2.0 
!  format) for reading only.  Also reads FTI and TITLE strings. 
!  (bmy, 7/30/02, 10/15/02)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) IUNIT    (INTEGER )  : Logical unit number of the file to be opened
!  (2 ) FILENAME (CHARACTER) : Name of the file to be opened
!
!  Arguments as Output:
!  ============================================================================
!  (3 ) TITLE    (CHARACTER) : OPTIONAL: returns TITLE to calling program
!
!  NOTES:
!  (1 ) Now references ERROR_STOP from "error_mod.f" (bmy, 10/15/02)
!******************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD, ONLY : ERROR_STOP
      USE FILE_MOD,  ONLY : IOERROR

      ! Arguments
      INTEGER,           INTENT(IN)            :: IUNIT
      CHARACTER(LEN=*),  INTENT(IN)            :: FILENAME
      CHARACTER(LEN=80), INTENT(OUT), OPTIONAL :: TITLE

      ! Local variables
      INTEGER                                  :: IOS
      CHARACTER(LEN=40)                        :: FTI
      CHARACTER(LEN=80)                        :: TMP_TITLE

      !=================================================================
      ! OPEN_BPCH2_FOR_READ begins here!
      !=================================================================

      ! Open file for input -- readonly
      OPEN( IUNIT,      FILE=TRIM( FILENAME ), STATUS='OLD',
     &      IOSTAT=IOS, FORM='UNFORMATTED',    ACCESS='SEQUENTIAL' )

      ! Error check
      IF ( IOS /= 0 ) CALL IOERROR( IOS, IUNIT, 'open_bpch2_for_read:1')

      ! Read file type identifier
      READ( IUNIT, IOSTAT=IOS ) FTI

      ! Error check
      IF ( IOS /= 0 ) CALL IOERROR( IOS, IUNIT, 'open_bpch2_for_read:2')

      ! Stop if this is not a binary punch file
      IF ( TRIM( FTI ) /= 'CTM bin 02' ) THEN
         CALL ERROR_STOP( 'Invalid file format!', 
     &                    'OPEN_BPCH2_FOR_READ (bpch2_mod.f)')
      ENDIF

      ! Read top title
      READ( IUNIT, IOSTAT=IOS ) TMP_TITLE

      ! Error check
      IF ( IOS /= 0 ) CALL IOERROR( IOS, IUNIT, 'open_bpch2_for_read:3')

      ! Copy value of TMP_TITLE to TITLE for return 
      IF ( PRESENT( TITLE ) ) TITLE = TMP_TITLE

      ! Return to calling program
      END SUBROUTINE OPEN_BPCH2_FOR_READ

!------------------------------------------------------------------------------

      SUBROUTINE OPEN_BPCH2_FOR_WRITE( IUNIT, FILENAME, TITLE )
!
!******************************************************************************
!  Subroutine OPEN_BPCH2_FOR_WRITE opens a binary punch file (version 2.0)
!  for writing. (bmy, 7/30/02)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) IUNIT    (INTEGER )  : Logical unit number of the file to be opened
!  (2 ) FILENAME (CHARACTER) : Name of the file to be opened
!  (3 ) TITLE    (CHARACTER) : Optional: title for top of file
!
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE FILE_MOD, ONLY : IOERROR

      ! Arguments
      INTEGER,           INTENT(IN)           :: IUNIT
      CHARACTER(LEN=*),  INTENT(IN)           :: FILENAME
      CHARACTER(LEN=80), INTENT(IN), OPTIONAL :: TITLE

      ! Local variables
      INTEGER                                 :: IOS
      CHARACTER(LEN=80)                       :: TMP_TITLE

      !=================================================================
      ! OPEN_BPCH2_FOR_WRITE begins here!
      !=================================================================

      ! If TITLE is not passed, create a default title string
      IF ( PRESENT( TITLE ) ) THEN
         TMP_TITLE = TITLE
      ELSE
         TMP_TITLE = 'GEOS-CHEM binary punch file v. 2.0'
      ENDIF

      ! Open file for output
      OPEN( IUNIT,      FILE=TRIM( FILENAME ), STATUS='UNKNOWN',
     &      IOSTAT=IOS, FORM='UNFORMATTED',    ACCESS='SEQUENTIAL' )

      ! Error check
      IF ( IOS /= 0 ) CALL IOERROR( IOS, IUNIT,'open_bpch2_for_write:1')

      ! Write the top-of-file title to disk
      CALL BPCH2_HDR( IUNIT, TMP_TITLE )

      ! Return to calling program
      END SUBROUTINE OPEN_BPCH2_FOR_WRITE

!------------------------------------------------------------------------------

      SUBROUTINE BPCH2_HDR ( IUNIT, TITLE )
!
!******************************************************************************
!  Subroutine BPCH2_HDR writes a header at the top of the binary
!  punch file, version 2.0 (bmy, 5/27/99, 7/30/02).
!
!  Arguments as input:
!  ============================================================================
!  (1) IUNIT : INTEGER - logical unit number of binary punch file
!  (2) TITLE : CHAR*80 - description of data contained in binary punch file
!
!  NOTES:
!  (1 ) Added this routine to "bpch_mod.f" (bmy, 6/28/00)
!  (2 ) Use IOS /= 0 criterion to also check for EOF condition (bmy, 9/12/00)
!  (3 ) Now reference IOERROR from "file_mod.f". (bmy, 6/26/02)
!******************************************************************************
!
      ! References to F90 modules
      USE FILE_MOD, ONLY : IOERROR

      ! Arguments
      INTEGER,           INTENT(IN) :: IUNIT
      CHARACTER(LEN=80), INTENT(IN) :: TITLE

      ! Local variable
      INTEGER                       :: IOS
      CHARACTER(LEN=40)             :: FTI = 'CTM bin 02'

      !=================================================================
      ! BPCH2_HDR begins here!
      !
      ! Write header information to binary punch file 
      ! Also be sure to trap I/O Error conditions
      !=================================================================
      WRITE ( IUNIT, IOSTAT=IOS ) FTI
      IF ( IOS /= 0 ) CALL IOERROR( IOS, IUNIT, 'bpch2_hdr:1' )

      WRITE ( IUNIT, IOSTAT=IOS ) TITLE
      IF ( IOS /= 0 ) CALL IOERROR( IOS, IUNIT, 'bpch2_hdr:2' )

      ! Return to calling program    
      END SUBROUTINE BPCH2_HDR

!------------------------------------------------------------------------------

      SUBROUTINE BPCH2( IUNIT,     MODELNAME, LONRES,   LATRES,
     &                  HALFPOLAR, CENTER180, CATEGORY, NTRACER,    
     &                  UNIT,      TAU0,      TAU1,     RESERVED,   
     &                  NI,        NJ,        NL,       IFIRST,     
     &                  JFIRST,    LFIRST,    ARRAY )
!
!******************************************************************************
!  Subroutine BPCH2 writes binary punch file (version 2.0) to disk.
!  Information about the model grid is also stored with each data block.
!  (bmy, 5/27/99, 7/30/02)
!
!  Arguments as input:
!  ============================================================================
!  (1    ) IUNIT      : INTEGER  - logical unit number of the file 
!  (2    ) MODELNAME  : CHAR*40  - Name of model used to create output
!  (3    ) LONRES     : REAL*4   - Longitude resolution of grid, in degrees
!  (4    ) LATRES     : REAL*4   - Latitude resolution of grid, in degrees
!  (4    ) HALFPOLAR  : INTEGER  - flag, =1 if model has half-polar boxes
!  (5    ) CENTER180  : INTEGER  - flag, =1 if model has lon center on 180 deg
!  (6    ) CATEGORY   : CHAR*40  - diagnostic category name
!  (7    ) NTRACER    : INTEGER  - number of tracer
!  (8    ) UNIT       : CHAR*40  - units of data
!  (9    ) TAU0       : REAL*8   - TAU at start of diagnostic interval
!  (10   ) TAU1       : REAL*8   - TAU at end   of diagnostic interval
!  (11   ) RESERVED   : CHAR*40  - Reserved for future use
!  (12-14) NI,NJ,NL   : INTEGER  - dimensions of ARRAY
!  (15   ) IFIRST     : INTEGER  - I-index of the first grid box
!  (16   ) JFIRST     : INTEGER  - J-index of the first grid box
!  (17   ) LFIRST     : INTEGER  - L-index of the first grid box
!  (18   ) ARRAY      : REAL*4   - data block to be written to the file
!
!  NOTES:
!  (1 ) Added indices to IOERROR calls (e.g. "bpch2:1", "bpch2:2", etc.) 
!        (bmy, 10/4/99)
!  (2 ) Added this routine to "bpch_mod.f" (bmy, 6/28/00)
!  (3 ) Use IOS /= 0 criterion to also check for EOF condition (bmy, 9/12/00)
!  (4 ) Now reference IOERROR from "file_mod.f". (bmy, 6/26/02)
!******************************************************************************
!  
      ! References to F90 modules
      USE FILE_MOD, ONLY : IOERROR

      ! Arguments
      INTEGER,           INTENT(IN) :: IUNIT
      INTEGER,           INTENT(IN) :: NTRACER 
      INTEGER,           INTENT(IN) :: NI, NJ, NL 
      INTEGER,           INTENT(IN) :: IFIRST, JFIRST, LFIRST
      INTEGER,           INTENT(IN) :: HALFPOLAR, CENTER180
      REAL*4,            INTENT(IN) :: ARRAY( NI, NJ, NL )
      REAL*4,            INTENT(IN) :: LONRES, LATRES
      REAL*8,            INTENT(IN) :: TAU0,   TAU1
      CHARACTER(LEN=20), INTENT(IN) :: MODELNAME
      CHARACTER(LEN=40), INTENT(IN) :: CATEGORY
      CHARACTER(LEN=40), INTENT(IN) :: RESERVED
      CHARACTER(LEN=40), INTENT(IN) :: UNIT

      ! Local variables
      INTEGER                       :: I, J, L, NSKIP, IOS

      ! For computing NSKIP
      INTEGER, PARAMETER            :: BYTES_PER_NUMBER = 4
      INTEGER, PARAMETER            :: END_OF_RECORD    = 8

      !=================================================================
      ! BPCH2 begins here!!  
      !
      ! Compute the number of bytes to skip between the end of one 
      ! data block and the beginning of the next data header line
      !=================================================================
      NSKIP = ( BYTES_PER_NUMBER * ( NI * NJ * NL ) ) + END_OF_RECORD

      !=================================================================
      ! Write data block to binary punch file
      ! Check for I/O errors
      !=================================================================
      WRITE( IUNIT, IOSTAT=IOS ) 
     &     MODELNAME, LONRES, LATRES, HALFPOLAR, CENTER180

      IF ( IOS /= 0 ) CALL IOERROR( IOS, IUNIT, 'bpch2:1' )

      WRITE( IUNIT, IOSTAT = IOS ) 
     &     CATEGORY, NTRACER,  UNIT, TAU0,   TAU1,   RESERVED,
     &     NI,       NJ,       NL,   IFIRST, JFIRST, LFIRST,
     &     NSKIP

      IF ( IOS /= 0 ) CALL IOERROR( IOS, IUNIT, 'bpch2:2' )

      WRITE( IUNIT, IOSTAT=IOS ) 
     &     ( ( ( ARRAY(I,J,L), I=1,NI ), J=1,NJ ), L=1,NL )

      IF ( IOS /= 0 ) CALL IOERROR( IOS, IUNIT, 'bpch2:3' )

      !=================================================================
      ! Return to calling program      
      !=================================================================
      END SUBROUTINE BPCH2

!------------------------------------------------------------------------------

      SUBROUTINE READ_BPCH2( FILENAME, CATEGORY_IN, TRACER_IN, 
     &                       TAU0_IN,  IX,          JX,          
     &                       LX,       ARRAY,       QUIET ) 
!
!******************************************************************************
!  Subroutine READ_BPCH2 reads a binary punch file (v. 2.0) and extracts
!  a data block that matches the given category, tracer, and tau value.
!  (bmy, 12/10/99, 3/14/03)
!
!  Arguments as Input:
!  ============================================================================
!  (1  ) FILENAME    : (CHARACTER) String for input file name
!  (2  ) CATEGORY_IN : (CHARACTER) Category name for the desired data block
!  (3  ) TRACER_IN   : (INTEGER  ) Tracer number for which to extract data
!  (4  ) TAU0_IN     : (REAL*8   ) TAU value for which to extract data
!  (5-7) IX, JX, LX  : (INTEGER  ) Dimensions of ARRAY (see below) 
!  (9  ) QUIET       : (LOGICAL  ) Optional flag for suppressing printing
!
!  Arguments as Output:
!  ============================================================================
!  (8  ) ARRAY       : (REAL*4   ) Array to hold extracted data values
!
!  NOTES:
!  (1 ) Assumes that we are reading in a global-size data block.
!  (2 ) Trap all I/O errors with subroutine IOERROR.F.
!  (3 ) Now stop with an error message if no matches are found. (bmy, 3/9/00)
!  (4 ) Added this routine to "bpch_mod.f" (bmy, 6/28/00)
!  (5 ) Use IOS /= 0 criterion to also check for EOF condition (bmy, 9/12/00)
!  (6 ) TEMPARRAY now dimensioned to be of global size (bmy, 10/12/00) 
!  (7 ) Removed obsolete code from 10/12/00 (bmy, 12/18/00)
!  (8 ) Now make TEMPARRAY independent of CMN_SIZE parameters (bmy, 4/17/01)
!  (9 ) Removed old commented-out code (bmy, 4/20/01)
!  (10) Now reference IU_FILE and IOERROR from "file_mod.f".  Now call 
!        OPEN_BPCH2_FOR_READ to open the binary punch file.  Now use IU_FILE
!        as the unit number instead of a locally-defined IUNIT. (bmy, 7/30/02)
!  (11) Now references ERROR_STOP from "error_mod.f" (bmy, 10/15/02)
!  (12) Now set IFIRST=1, JFIRST=1 for 1x1 nested grids.  Now needs to
!        reference "define.h".  Added OPTIONAL QUIET flag. (bmy, 3/14/03) 
!******************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD, ONLY : ERROR_STOP
      USE FILE_MOD,  ONLY : IU_FILE, IOERROR

#     include "define.h" 

      ! Arguments
      LOGICAL, OPTIONAL, INTENT(IN)  :: QUIET
      INTEGER,           INTENT(IN)  :: IX, JX, LX, TRACER_IN
      CHARACTER(LEN=*),  INTENT(IN)  :: FILENAME, CATEGORY_IN 
      REAL*8,            INTENT(IN)  :: TAU0_IN
      REAL*4,            INTENT(OUT) :: ARRAY(IX, JX, LX)      

      ! Local variables
      LOGICAL            :: FOUND, TMP_QUIET
      INTEGER            :: I,  J,  L,  N,  IOS, M
      INTEGER            :: I1, I2, J1, J2, L1,  L2
      CHARACTER(LEN=255) :: MSG
      
      ! Make TEMPARRAY big enough to for a 1x1 grid (bmy, 4/17/01)
      REAL*4             :: TEMPARRAY(360,181,70)

      ! For binary punch file, version 2.0
      INTEGER            :: NTRACER,   NSKIP
      INTEGER            :: HALFPOLAR, CENTER180
      INTEGER            :: NI,        NJ,        NL
      INTEGER            :: IFIRST,    JFIRST,    LFIRST
      REAL*4             :: LONRES,    LATRES
      REAL*8             :: ZTAU0,     ZTAU1
      CHARACTER(LEN=20)  :: MODELNAME
      CHARACTER(LEN=40)  :: CATEGORY
      CHARACTER(LEN=40)  :: UNIT     
      CHARACTER(LEN=40)  :: RESERVED

      !=================================================================
      ! READ_BPCH2 begins here!
      !  
      ! Initialize some variables
      !=================================================================
      FOUND            = .FALSE.
      ARRAY(:,:,:)     = 0e0
      TEMPARRAY(:,:,:) = 0e0

      ! Define a temporary variable for QUIET
      IF ( PRESENT( QUIET ) ) THEN
         TMP_QUIET = QUIET
      ELSE
         TMP_QUIET = .FALSE.
      ENDIF

      !=================================================================
      ! Open binary punch file and read top-of-file header.
      ! Do some error checking to make sure the file is the right format.
      !=================================================================
      CALL OPEN_BPCH2_FOR_READ( IU_FILE, FILENAME )

      !=================================================================
      ! Read data from the binary punch file 
      !
      ! NOTE: IOS < 0 is end-of-file, IOS > 0 is error condition
      !=================================================================
      DO
         READ( IU_FILE, IOSTAT=IOS ) 
     &        MODELNAME, LONRES, LATRES, HALFPOLAR, CENTER180
         
         IF ( IOS < 0 ) EXIT
         IF ( IOS > 0 ) CALL IOERROR( IOS, IU_FILE, 'read_bpch2:4' )

         READ( IU_FILE, IOSTAT=IOS ) 
     &        CATEGORY, NTRACER,  UNIT, ZTAU0,  ZTAU1,  RESERVED,
     &        NI,       NJ,       NL,   IFIRST, JFIRST, LFIRST,
     &        NSKIP

         IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'read_bpch2:5' )

         READ( IU_FILE, IOSTAT=IOS ) 
     &        ( ( ( TEMPARRAY(I,J,L), I=1,NI ), J=1,NJ ), L=1,NL )

         IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'read_bpch2:6' )
         
         ! Test for a match
         IF ( TRIM( CATEGORY_IN ) == TRIM( CATEGORY ) .and. 
     &        TRACER_IN           == NTRACER          .and.
     &        TAU0_IN             == ZTAU0 ) THEN
            FOUND = .TRUE.
            EXIT
         ENDIF

      ENDDO

      !=================================================================
      ! We have found a match!  Copy TEMPARRAY to ARRAY, taking into 
      ! account the starting positions (IFIRST, JFIRST, LFIRST) of 
      ! the data block.
      !=================================================================
      IF ( FOUND ) THEN 

#if   defined( GRID1x1 ) 
         ! This is a kludge to overwrite the IFIRST, JFIRST, LFIRST For
         ! the 1x1 nested grid.  1x1 met fields & other data are already
         ! cut down to size to save space. (bmy, 3/11/03)
         I1 = 1
         J1 = 1
         L1 = LFIRST
#else
         ! Otherwise IFIRST, JFIRST, FIRST from the file (bmy, 3/11/03)
         I1 = IFIRST
         J1 = JFIRST
         L1 = LFIRST
#endif     
 
         I2 = NI + I1 - 1
         J2 = NJ + J1 - 1
         L2 = NL + L1 - 1
                  
         ARRAY( I1:I2, J1:J2, L1:L2 ) = TEMPARRAY( 1:NI, 1:NJ, 1:NL )

         ! Flag to decide whether or not we will echo info (bmy, 3/14/03)
         IF ( .not. TMP_QUIET ) THEN 
            WRITE( 6, 100 ) ZTAU0, NTRACER
 100        FORMAT( 'READ_BPCH2: Found data for TAU = ', f10.2, 
     &              ' and tracer # ', i6 )
         ENDIF

      ELSE
         MSG = 'No matches found for file ' // TRIM( FILENAME ) // '!'
         CALL ERROR_STOP( MSG, 'READ_BPCH2 (bpch2_mod.f)!' )
      ENDIF

      !=================================================================
      ! Close file and quit
      !=================================================================
      CLOSE( IU_FILE )

      ! Return to calling program
      END SUBROUTINE READ_BPCH2

!------------------------------------------------------------------------------

      FUNCTION GET_MODELNAME() RESULT( MODELNAME )
!
!******************************************************************************
!  Function GET_MODELNAME returns the proper value of MODELNAME for GEOS-1,
!  GEOS-STRAT, GEOS-2, or GEOS-3 data.  MODELNAME is written to the binary
!  punch file and is used by the GAMAP package. (bmy, 6/22/00, 11/20/01)
!
!  NOTES:
!  (1 ) Now use special model name for GEOS-3 w/ 30 layers (bmy, 10/9/01)
!  (2 ) Added modelname for GEOS-4/fvDAS model type (bmy, 11/20/01)
!  (2 ) Added "GEOS4_30L" for reduced GEOS-4 grid.  Also now use C-preprocessor
!        switch "GRID30LEV" instead of IF statements. (bmy, 11/3/03)
!******************************************************************************
!
#     include "CMN_SIZE"

      ! MODELNAME holds the return value for the function
      CHARACTER(LEN=20) :: MODELNAME

      !=================================================================
      ! GET_MODELNAME begins here!
      !=================================================================

#if   defined( GEOS_1 ) 
      MODELNAME = 'GEOS1'
     
#elif defined( GEOS_STRAT ) 
      MODELNAME = 'GEOS_STRAT'

!------------------------------------------------------------
! Prior to 11/3/03:
! Remove GEOS_2 model name -- it's obsolete! (bmy, 11/3/03)
!#elif defined( GEOS_2 ) 
!      MODELNAME = 'GEOS2'
!------------------------------------------------------------

#elif defined( GEOS_3 )
      
      !----------------------------------------------------------------------
      ! Prior to 11/3/03:
      !! Write a special model name to the punch file for GAMAP
      !! if we are using regridded vertical resolution (bmy, 10/9/01)
      ! Now use GRID30LEV Cpp switch -- eliminate IF statement (bmy, 11/3/03)
      !IF ( LLPAR == 30 ) THEN 
      !   MODELNAME = 'GEOS3_30L'
      !ELSE
      !   MODELNAME = 'GEOS3'
      !ENDIF
      !-----------------------------------------------------------------------
#if   defined( GRID30LEV )
      MODELNAME = 'GEOS3_30L'   ! Special modelname for GAMAP (GEOS-3 30L)
#else
      MODELNAME = 'GEOS3'       ! Normal modelname (GEOS-3 48L)
#endif


#elif defined( GEOS_4 )
      !--------------------------------------------------------------
      ! Prior to 11/3/03:
      ! Need to save the bpch file under the name "GEOS4_30L" for
      ! the reduced 30-level grid (bmy, 11/3/03)
      !MODELNAME = 'GEOS4'
      !--------------------------------------------------------------
#if   defined( GRID30LEV )
      MODELNAME = 'GEOS4_30L'   ! Special modelname for GAMAP (GEOS-4 30L)
#else 
      MODELNAME = 'GEOS4'       ! Original modelname (GEOS-4 55L)
#endif

#endif

      ! Return to calling program
      END FUNCTION GET_MODELNAME

!------------------------------------------------------------------------------

      FUNCTION GET_NAME_EXT() RESULT( NAME_EXT )
!
!******************************************************************************
!  Function GET_NAME_EXT returns the proper filename extension for CTM
!  model name (i.e. "geos1", "geoss", "geos3", or "geos4").  
!  (bmy, 6/28/00, 11/3/03)
!  
!  NOTES:
!  (1 ) Added name string for GEOS-4/fvDAS model type (bmy, 11/20/01)
!  (2 ) Remove obsolete "geos2" model name strning (bmy, 11/3/03)
!******************************************************************************
!
#     include "define.h"

      ! EXTENSION holds the return value for the function
      CHARACTER(LEN=5) :: NAME_EXT

#if   defined( GEOS_1 ) 
      NAME_EXT = 'geos1'
     
#elif defined( GEOS_STRAT ) 
      NAME_EXT = 'geoss'

!-----------------------------------------------------------------
! Prior to 11/3/03:
! Remove obsolete "geos2" model name strning (bmy, 11/3/03)
!#elif defined( GEOS_2 ) 
!      NAME_EXT = 'geos2'
!-----------------------------------------------------------------

#elif defined( GEOS_3 )
      NAME_EXT = 'geos3'

#elif defined( GEOS_4 )
      NAME_EXT = 'geos4'

#endif

      ! Return to calling program
      END FUNCTION GET_NAME_EXT

!------------------------------------------------------------------------------

      FUNCTION GET_RES_EXT() RESULT( RES_EXT )
!
!******************************************************************************
!  Function GET_RES_EXT returns the proper filename extension for
!  CTM grid resolution (i.e. "1x1", "2x25", "4x5").  (bmy, 6/28/00)
!******************************************************************************
!
#     include "define.h"

#if   defined( GRID4x5 )
      CHARACTER(LEN=3) :: RES_EXT
      RES_EXT = '4x5'
     
#elif defined( GRID2x25 ) 
      CHARACTER(LEN=4) :: RES_EXT
      RES_EXT = '2x25'

#elif defined( GRID1x1 ) 
      CHARACTER(LEN=3) :: RES_EXT
      RES_EXT = '1x1'

#endif

      END FUNCTION GET_RES_EXT

!------------------------------------------------------------------------------

      FUNCTION GET_TAU0_6A( MONTH, DAY, YEAR, HOUR, MIN, SEC ) 
     &         RESULT( THIS_TAU0 )
!
!******************************************************************************
!  Function GET_TAU0_6A returns the corresponding TAU0 value for the first 
!  day of a given MONTH of a given YEAR.  This is necessary to index monthly 
!  mean binary punch files, which are used as input to GEOS-CHEM.
!  (bmy, 9/26/01) 
!
!  This function takes 3 mandatory arguments (MONTH, DAY, YEAR) and 3 
!  optional arguments (HOUR, MIN, SEC).  It is intended to replace the current 
!  2-argument version of GET_TAU0.  The advantage being that GET_TAU0_3A can 
!  compute a TAU0 for any date and time in the GEOS-CHEM epoch, rather than 
!  just the first day of each month.  Overload this w/ an interface so that 
!  the user can also choose the version of GET_TAU0 w/ 2 arguments 
!  (MONTH, YEAR), which is the prior version.
!
!  Arguments as Input:
!  ===========================================================================
!  (1 ) MONTH (INTEGER) : Month of year (1-12)
!  (2 ) DAY   (INTEGER) : Day of month (1-31)
!  (3 ) YEAR  (INTEGER) : 4-digit year number (e.g. 1985,2001)
!  (4 ) HOUR  (INTEGER) : OPTIONAL: Hour of day (0-24)
!  (5 ) MIN   (INTEGER) : OPTIONAL: Minute of hour (0-59)
!  (6 ) SEC   (INTEGER) : OPTIONAL: Seconds of minute (0-59)
!
!  NOTES: 
!  (1 ) 1985 is the first year of the GEOS epoch.
!  (2 ) Add TAU0 values for years 1985-2001 (bmy, 8/1/00)
!  (3 ) Correct error for 1991 TAU values.  Also added 2002 and 2003.
!        (bnd, bmy, 1/4/01)
!  (4 ) Updated comments  (bmy, 9/26/01)
!  (5 ) Now references JULDAY from "julday_mod.f" (bmy, 11/20/01)
!  (6 ) Now references ERROR_STOP from "error_mod.f"  (bmy, 10/15/02)
!******************************************************************************
!
      ! Reference to F90 modules
      USE ERROR_MOD,  ONLY : ERROR_STOP
      USE JULDAY_MOD, ONLY : JULDAY

      ! Arguments
      INTEGER, INTENT(IN)           :: MONTH
      INTEGER, INTENT(IN)           :: DAY
      INTEGER, INTENT(IN)           :: YEAR
      INTEGER, INTENT(IN), OPTIONAL :: HOUR
      INTEGER, INTENT(IN), OPTIONAL :: MIN
      INTEGER, INTENT(IN), OPTIONAL :: SEC

      ! Local variables
      INTEGER                       :: TMP_HOUR, TMP_MIN, TMP_SEC
      REAL*8                        :: DAYS

      ! Return value
      REAL*8                        :: THIS_TAU0
      
      !=================================================================
      ! GET_TAU0_6A begins here!
      !=================================================================

      ! Error checking 
      IF ( MONTH < 1 .or. MONTH > 12 ) THEN
         CALL ERROR_STOP ( 'Invalid MONTH selection!', 'GET_TAU0' )
      ENDIF

      ! Error checking 
      IF ( DAY < 1 .or. DAY > 31 ) THEN
         CALL ERROR_STOP ( 'Invalid DAY selection!', 'GET_TAU0' )
      ENDIF

      ! If HOUR isn't passed, default to 0
      IF ( PRESENT( HOUR ) ) THEN
         TMP_HOUR = HOUR
      ELSE
         TMP_HOUR = 0
      ENDIF 

      ! If MIN isn't passed, default to 0
      IF ( PRESENT( MIN ) ) THEN
         TMP_MIN = MIN
      ELSE
         TMP_MIN = 0 
      ENDIF 

      ! If SEC isn't passed, default to 0
      IF ( PRESENT( SEC ) ) THEN
         TMP_SEC = SEC
      ELSE
         TMP_SEC = 0 
      ENDIF 

      ! Number of days since midnight on 1/1/1985
      THIS_TAU0 = JULDAY( YEAR, MONTH, DBLE( DAY ) ) - 2446066.5d0

      ! Multiply by 24 to get hours since 1/1/1985
      ! Also add in the hours elapsed since midnight on this date
      THIS_TAU0 = ( THIS_TAU0 * 24d0 ) + ( TMP_HOUR         ) + 
     &            ( TMP_MIN   / 60d0 ) + ( TMP_SEC / 3600d0 )

      ! Return to calling program
      END FUNCTION GET_TAU0_6A

!------------------------------------------------------------------------------

      END MODULE BPCH2_MOD
