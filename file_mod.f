! $Id: file_mod.f,v 1.5 2004/05/03 14:46:17 bmy Exp $
      MODULE FILE_MOD
!
!******************************************************************************
!  Module FILE_MOD contains file unit numbers, as well as file I/O routines
!  for GEOS-CHEM.  FILE_MOD keeps all of the I/O unit numbers in a single
!  location for convenient access. (bmy, 7/1/02, 4/23/04)
!
!  Module Variables:
!  ============================================================================
!  (1 ) IU_RST      : Unit # for file "gctm.trc.YYYYMMDD"
!  (2 ) IU_INPUT    : Unit # for file "input.ctm"       
!  (3 ) IU_CHEMDAT  : Unit # for file "chem.dat"        
!  (4 ) IU_FASTJ    : Unit # for file "ratj.d", "jv_atms.dat", "jv_spec.dat"
!  (5 ) IU_INPTR    : Unit # for file "inptr.ctm"       
!  (6 ) IU_GEOS     : Unit # for file "input.geos"      
!  (7 ) IU_TS       : Unit # for file "ctm.ts"          
!  (8 ) IU_BPCH     : Unit # for file "ctm.bpch"        
!  (9 ) IU_ND20     : Unit # for file "rate.YYYYMMDD"   
!  (10) IU_ND49     : Unit # for file "tsYYYYMMDD.bpch" 
!  (11) IU_ND50     : Unit # for file "ts24h.bpch"
!  (12) IU_ND51     : Unit # for file "ts10_12am.bpch" or "ts1_4pm.bpch"
!  (13) IU_PLANE    : Unit # for plane flight diagnostic output file
!  (14) IU_FILE     : Unit # for files opened & closed in same routine
!  (15) IU_PH       : Unit # for GEOS-CHEM PHIS met field file
!  (16) IU_I6       : Unit # for GEOS-CHEM I-6  met field file
!  (17) IU_A6       : Unit # for GEOS-CHEM A-6  met field file
!  (18) IU_A3       : Unit # for GEOS-CHEM A-3  met field file
!  (19) IU_KZZ      : Unit # for GEOS-CHEM KZZ  met field file
!  (20) IU_GWET     : Unit # for GEOS-CHEM GWET met field file
!  (21) IU_SMV2LOG  : Unit # for "smv2.log" file -- SMVGEAR II rxns & species
!  (22) IU_DEBUG    : Unit # left for debugging purposes
!
!  Module Routines
!  ============================================================================
!  (1 ) IOERROR     : Stops w/ error msg output if I/O errors are detected
!  (2 ) CLOSE_FILES : Closes all files at the end of a GEOS-CHEM run
!
!  GEOS-CHEM modules referenced by file_mod.f
!  ============================================================================
!  (1 ) error_mod.f : Module containing NaN and other error check routines
!
!  NOTES:
!  (1 ) Moved "ioerror.f" into this module. (bmy, 7/1/02)
!  (2 ) Now references "error_mod.f" (bmy, 10/15/02)
!  (3 ) Renamed cpp switch from DEC_COMPAQ to COMPAQ.  Also added code to 
!        trap I/O errors on SUN/Sparc platform. (bmy, 3/23/03)
!  (4 ) Now added IU_BC for nested boundary conditions as unit 18 
!        (bmy, 3/27/03)
!  (5 ) Renamed IU_CTMCHEM to IU_SMV2LOG (bmy, 4/21/03)
!  (6 ) Now print out I/O errors for IBM and INTEL_FC compilers (bmy, 11/6/03)
!  (7 ) Changed the name of some cpp switches in "define.h" (bmy, 12/2/03)
!******************************************************************************
!
      IMPLICIT NONE

      !=================================================================
      ! MODULE VARIABLES
      !=================================================================    
      INTEGER, PARAMETER :: IU_RST     = 1
      INTEGER, PARAMETER :: IU_INPUT   = 5
      INTEGER, PARAMETER :: IU_CHEMDAT = 7
      INTEGER, PARAMETER :: IU_FASTJ   = 8
      INTEGER, PARAMETER :: IU_INPTR   = 9
      INTEGER, PARAMETER :: IU_GEOS    = 10
      INTEGER, PARAMETER :: IU_TS      = 11      
      INTEGER, PARAMETER :: IU_BPCH    = 12
      INTEGER, PARAMETER :: IU_ND20    = 13
      INTEGER, PARAMETER :: IU_ND49    = 14
      INTEGER, PARAMETER :: IU_ND50    = 15
      INTEGER, PARAMETER :: IU_ND51    = 16
      INTEGER, PARAMETER :: IU_ND52    = 20  ! For ICARTT (bmy, 4/21/04)
      INTEGER, PARAMETER :: IU_PLANE   = 17
      INTEGER, PARAMETER :: IU_BC      = 18
      INTEGER, PARAMETER :: IU_FILE    = 65
      INTEGER, PARAMETER :: IU_PH      = 70
      INTEGER, PARAMETER :: IU_I6      = 71
      INTEGER, PARAMETER :: IU_A6      = 72
      INTEGER, PARAMETER :: IU_A3      = 73
      INTEGER, PARAMETER :: IU_KZZ     = 74
      INTEGER, PARAMETER :: IU_GWET    = 75
      INTEGER, PARAMETER :: IU_SMV2LOG = 93
      INTEGER, PARAMETER :: IU_DEBUG   = 98

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement 
      !=================================================================
      CONTAINS

!------------------------------------------------------------------------------

      SUBROUTINE IOERROR( ERROR_NUM, UNIT, LOCATION )
!
!******************************************************************************
!  Subroutine IOERRROR prints out I/O error messages.  The error number, 
!  file unit, location, and a brief description will be printed, and 
!  program execution will be halted. (bmy, 5/28/99, 4/23/04)
!
!  Arguments as input:
!  ===========================================================================
!  (1 ) ERROR_NUM : I/O error number (output from the IOSTAT flag)
!  (2 ) UNIT      : Unit # of the file where the I/O error occurred
!  (3 ) LOCATION  : Name of the routine in which the error occurred
!
!  NOTES:
!  (1 ) Now flush the standard output buffer before stopping.  
!        Also updated comments. (bmy, 2/7/00)
!  (2 ) Changed ROUTINE_NAME to LOCATION.  Now also use C-library routines
!        gerror and strerror() to get the error string corresponding to 
!        ERROR_NUM.  For SGI platform, also print the command string that
!        will call the SGI "explain" command, which will yield additional
!        information about the error.  Updated comments, cosmetic changes.  
!        Now also reference "define.h". (bmy, 3/21/02)
!  (3 ) Moved into "file_mod.f".  Now reference GEOS_CHEM_STOP from module
!        "error_mod.f".  Updated comments, cosmetic changes. (bmy, 10/15/02)
!  (4 ) Renamed cpp switch from DEC_COMPAQ to COMPAQ.  Also added code to 
!        display I/O errors on SUN platform. (bmy, 3/23/03)
!  (5 ) Now call GERROR for IBM and INTEL_FC compilers (bmy, 11/6/03)
!  (6 ) Renamed SGI to SGI_MIPS, LINUX to LINUX_PGI, INTEL_FC to INTEL_IFC, 
!        and added LINUX_EFC. (bmy, 12/2/03)
!  (7 ) Now don't flush the buffer for LINUX_EFC (bmy, 4/23/04)
!******************************************************************************
!  
      ! References to F90 modules
      USE ERROR_MOD, ONLY : GEOS_CHEM_STOP

      IMPLICIT NONE

#     include "define.h"           ! C-preprocessor switches
    
      ! Arguments
      INTEGER,          INTENT(IN) :: ERROR_NUM, UNIT
      CHARACTER(LEN=*), INTENT(IN) :: LOCATION 

      ! Local variables
      CHARACTER(LEN=10)            :: ERROR_NUMSTR 
      CHARACTER(LEN=255)           :: ERROR_MSG
      CHARACTER(LEN=255)           :: EXPLAIN_CMD

      ! External functions
      CHARACTER(LEN=255), EXTERNAL :: GERROR

      !=================================================================
      ! IOERROR begins here!
      !=================================================================

      ! Fancy output
      WRITE( 6, '(a)' ) REPEAT( '=', 79 )

      ! Write error number, unit, location
      WRITE( 6, 110 ) ERROR_NUM, UNIT, TRIM( LOCATION )
 110  FORMAT( 'GEOS-CHEM I/O ERROR ', i5, ' in file unit ', i5, /, 
     &        'Encountered at routine:location ', a )

#if   defined( SGI_MIPS )
      
      !=================================================================
      ! For SGI: print error msg and construct explain command string
      !=================================================================
      IF ( ERROR_NUM == 2 ) THEN   

         ! Error 2 is "file not found", so handle that separately.
         ! You can't use the explain command w/ error 2.
         WRITE( 6, '(/,a)' ) 'Error: No such file or directory'
         
      ELSE

         ! Call SGI strerror routine to convert ERROR_NUM to ERROR_MSG
         ERROR_MSG = GERROR()

         ! Print error message to std output
         WRITE( 6, 120 ) TRIM( ERROR_MSG )
 120     FORMAT( /, 'Error: ', a )
      
         ! Convert ERROR_NUM to string format
         WRITE( ERROR_NUMSTR, '(i10)' ) ERROR_NUM

         ! Construct argument for SGI explain command
         IF ( ERROR_NUM >= 1000 .and. ERROR_NUM < 4000 ) THEN 
            EXPLAIN_CMD = 'explain cf90-' // 
     &                     TRIM( ADJUSTL( ERROR_NUMSTR ))

         ELSE IF ( ERROR_NUM >= 4000 ) THEN
            EXPLAIN_CMD = 'explain lib-'  // 
     &                    TRIM( ADJUSTL( ERROR_NUMSTR ))
         ENDIF
      
         ! Print command string for the SGI explain command
         WRITE( 6, 130 ) TRIM( EXPLAIN_CMD )
 130     FORMAT( /, 'Type "', a, '" at the Unix prompt for an ',
     &              'explanation of the error.' )
      ENDIF

#elif defined( COMPAQ )

      !=================================================================
      ! For COMPAQ/Alpha: call gerror() to get the I/O error msg
      !=================================================================

      ! GERROR returns ERROR_MSG corresponding to ERROR_NUM 
      ERROR_MSG = GERROR()
 
      ! Print error message to std output
      WRITE( 6, 120 ) TRIM( ERROR_MSG )
 120  FORMAT( /, 'Error: ', a )

#elif defined( LINUX_PGI ) || defined( LINUX_IFC ) || defined( LINUX_EFC )

      !=================================================================
      ! For LINUX platform: call gerror() to get the I/O error msg
      !=================================================================

      ! GERROR returns ERROR_MSG corresponding to ERROR_NUM 
      ERROR_MSG = GERROR()
 
      ! Print error message to std output
      WRITE( 6, 120 ) TRIM( ERROR_MSG )
 120  FORMAT( /, 'Error: ', a )

#elif defined( SPARC ) 

      !=================================================================
      ! For SUN/Sparc platform: call gerror() to get the I/O error msg
      !=================================================================

      ! GERROR returns ERROR_MSG corresponding to ERROR_NUM 
      ERROR_MSG = GERROR()
 
      ! Print error message to std output
      WRITE( 6, 120 ) TRIM( ERROR_MSG )
 120  FORMAT( /, 'Error: ', a )

#elif defined( IBM_AIX ) 

      !=================================================================
      ! For IBM/AIX platform: call gerror() to get the I/O error msg
      !=================================================================

      ! GERROR returns ERROR_MSG corresponding to ERROR_NUM 
      ERROR_MSG = GERROR()
 
      ! Print error message to std output
      WRITE( 6, 120 ) TRIM( ERROR_MSG )
 120  FORMAT( /, 'Error: ', a )

#elif defined( INTEL_FC ) 

      !=================================================================
      ! For INTEL F90 Compiler: call gerror() to get the I/O error msg
      !=================================================================

      ! GERROR returns ERROR_MSG corresponding to ERROR_NUM 
      ERROR_MSG = GERROR()
 
      ! Print error message to std output
      WRITE( 6, 120 ) TRIM( ERROR_MSG )
 120  FORMAT( /, 'Error: ', a )

#endif
      
      ! Fancy output
      WRITE( 6, '(a)' ) REPEAT( '=', 79 )

#if   !defined( LINUX_EFC )
      CALL FLUSH( 6 )
#endif

      ! Deallocate arrays and stop safely 
      CALL GEOS_CHEM_STOP
      
      ! End of program
      END SUBROUTINE IOERROR

!------------------------------------------------------------------------------

      SUBROUTINE CLOSE_FILES
!
!*****************************************************************************
!  Subroutine CLOSE_FILES closes files used by GEOS-CHEM.  This should be
!  called only from the end of the "main.f" program. (bmy, 3/4/98, 3/27/03)
!
!  NOTES:
!  (1 ) Moved into "file_mod.f" (bmy, 6/27/02)
!  (2 ) Also close IU_BC (bmy, 3/27/03)
!*****************************************************************************
!     
      !=================================================================
      ! CLOSE_FILES begins here!
      !=================================================================
     
      ! Close all file units, regardless of whether they were
      ! ever opened or not.  This is for safety's sake. 
      CLOSE( IU_RST     )
      CLOSE( IU_INPUT   )
      CLOSE( IU_CHEMDAT )
      CLOSE( IU_FASTJ   )
      CLOSE( IU_INPTR   )
      CLOSE( IU_GEOS    )
      CLOSE( IU_TS      )     
      CLOSE( IU_BPCH    )
      CLOSE( IU_ND20    )
      CLOSE( IU_ND49    )
      CLOSE( IU_ND50    )
      CLOSE( IU_ND51    )
      CLOSE( IU_PLANE   )
      CLOSE( IU_BC      )
      CLOSE( IU_FILE    )
      CLOSE( IU_PH      )
      CLOSE( IU_I6      )
      CLOSE( IU_A6      )
      CLOSE( IU_A3      )
      CLOSE( IU_KZZ     )
      CLOSE( IU_GWET    )
      CLOSE( IU_SMV2LOG )
      CLOSE( IU_DEBUG   )

      ! For ICARTT (bmy, 4/21/04)
      CLOSE( IU_ND52    )

      ! Return to calling program
      END SUBROUTINE CLOSE_FILES

!------------------------------------------------------------------------------

      END MODULE FILE_MOD
