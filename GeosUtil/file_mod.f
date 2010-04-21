! $Id: file_mod.f,v 1.3 2010/03/15 19:33:19 ccarouge Exp $
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: file_mod.f
!
! !DESCRIPTION: Module FILE\_MOD contains file unit numbers, as well as file 
!  I/O routines for GEOS-Chem.  FILE\_MOD keeps all of the I/O unit numbers 
!  in a single location for convenient access.
!\\
!\\
! !INTERFACE: 
!
      MODULE FILE_MOD
! 
! !USES:
!
      IMPLICIT NONE
      PRIVATE
!
! !DEFINED PARAMETERS:
!
      ! Logical file unit numbers for ...
      INTEGER, PUBLIC, PARAMETER :: IU_RST     = 1   ! Tracer restart file
      INTEGER, PUBLIC, PARAMETER :: IU_CHEMDAT = 7   ! "chem.dat" 
      INTEGER, PUBLIC, PARAMETER :: IU_FASTJ   = 8   ! FAST-J input files
      INTEGER, PUBLIC, PARAMETER :: IU_GEOS    = 10  ! "input.geos" 
      INTEGER, PUBLIC, PARAMETER :: IU_BPCH    = 11  ! "ctm.bpch" 
      INTEGER, PUBLIC, PARAMETER :: IU_ND20    = 12  ! "rate.YYYYMMDD"   
      INTEGER, PUBLIC, PARAMETER :: IU_ND48    = 13  ! ND48 output
      INTEGER, PUBLIC, PARAMETER :: IU_ND49    = 14  ! "tsYYYYMMDD.bpch" 
      INTEGER, PUBLIC, PARAMETER :: IU_ND50    = 15  ! "ts24h.bpch"
      INTEGER, PUBLIC, PARAMETER :: IU_ND51    = 16  ! "ts10_12am.bpch" etc.
      INTEGER, PUBLIC, PARAMETER :: IU_ND51b   = 23  ! for ND51b diagnostic
      INTEGER, PUBLIC, PARAMETER :: IU_ND52    = 17  ! ND52 output (NRT only)
      INTEGER, PUBLIC, PARAMETER :: IU_PLANE   = 18  ! "plane.log"
      INTEGER, PUBLIC, PARAMETER :: IU_BC      = 19  ! TPCORE BC files
      INTEGER, PUBLIC, PARAMETER :: IU_BC_NA   = 20  ! TPCORE BC files: NA grid
      INTEGER, PUBLIC, PARAMETER :: IU_BC_EU   = 21  ! TPCORE BC files: EU grid
      INTEGER, PUBLIC, PARAMETER :: IU_BC_CH   = 22  ! TPCORE BC files: CH grid
      INTEGER, PUBLIC, PARAMETER :: IU_FILE    = 65  ! Generic file
      INTEGER, PUBLIC, PARAMETER :: IU_TP      = 69  ! "YYYYMMDD.tropp.*"
      INTEGER, PUBLIC, PARAMETER :: IU_PH      = 70  ! "YYYYMMDD.phis.*"
      INTEGER, PUBLIC, PARAMETER :: IU_I6      = 71  ! "YYYYMMDD.i6.*"
      INTEGER, PUBLIC, PARAMETER :: IU_A6      = 72  ! "YYYYMMDD.a6.*"
      INTEGER, PUBLIC, PARAMETER :: IU_A3      = 73  ! "YYYYMMDD.a3.*"
      INTEGER, PUBLIC, PARAMETER :: IU_KZZ     = 74  ! %%% NOW OBSOLETE %%%
      INTEGER, PUBLIC, PARAMETER :: IU_GWET    = 75  ! "YYYYMMDD.gwet.*"
      INTEGER, PUBLIC, PARAMETER :: IU_XT      = 76  ! "YYYYMMDD.xtra.*"
      INTEGER, PUBLIC, PARAMETER :: IU_SMV2LOG = 93  ! "smv2.log"
      INTEGER, PUBLIC, PARAMETER :: IU_DEBUG   = 98  ! Reserved for debugging
      ! For soaprod
      ! (dkh, 03/26/07)  
      INTEGER, PUBLIC, PARAMETER :: IU_OAP     = 99
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC  :: CLOSE_FILES
      PUBLIC  :: FILE_EXISTS
      PUBLIC  :: IOERROR

      INTERFACE FILE_EXISTS
         MODULE PROCEDURE FILE_EX_C
         MODULE PROCEDURE FILE_EX_I
      END INTERFACE
!
! !PRIVATE MEMBER FUNCTIONS:
!
      PRIVATE :: FILE_EX_C
      PRIVATE :: FILE_EX_I
!
! !REVISION HISTORY:
!  (1 ) Moved "ioerror.f" into this module. (bmy, 7/1/02)
!  (2 ) Now references "error_mod.f" (bmy, 10/15/02)
!  (3 ) Renamed cpp switch from DEC_COMPAQ to COMPAQ.  Also added code to 
!        trap I/O errors on SUN/Sparc platform. (bmy, 3/23/03)
!  (4 ) Now added IU_BC for nested boundary conditions as unit 18 
!        (bmy, 3/27/03)
!  (5 ) Renamed IU_CTMCHEM to IU_SMV2LOG (bmy, 4/21/03)
!  (6 ) Now print out I/O errors for IBM and INTEL_FC compilers (bmy, 11/6/03)
!  (7 ) Changed the name of some cpp switches in "define.h" (bmy, 12/2/03)
!  (8 ) Renumbered the order of the files.  Also removed IU_INPTR and 
!        IU_INPUT since they are now obsolete. (bmy, 7/20/04)
!  (9 ) Added overloaded routines FILE_EX_C and FILE_EX_I (bmy, 3/23/05)
!  (10) Added LINUX_IFORT switch for Intel v8 & v9 compilers (bmy, 10/18/05)
!  (11) Added IU_XT for GEOS3 XTRA met fields files for MEGAN (tmf, 10/20/05)
!  (12) Extra modification for Intel v9 compiler (bmy, 11/2/05)
!  (13) Now print IFORT error messages (bmy, 11/30/05)
!  (14) Remove support for LINUX_IFC & LINUX_EFC compilers (bmy, 8/4/06)
!  (15) Remove support for SGI & COMPAQ compilers (bmy, 7/8/09)
!  20 Nov 2009 - R. Yantosca - Added ProTeX headers
!  18 Dec 2009 - Aaron van D - Added file units IU_BC_NA, IU_BC_EU, IU_BC_CH
!  15 Mar 2010 - D. Henze    - Add IU_OAP for SOA restart file.  
!EOP
!------------------------------------------------------------------------------
!BOC
      CONTAINS
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ioerror
!
! !DESCRIPTION: Subroutine IOERROR prints out I/O error messages.  
!  The error number, file unit, location, and a brief description will 
!  be printed, and program execution will be halted. (bmy, 5/28/99, 7/4/09)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE IOERROR( ERROR_NUM, UNIT, LOCATION )
!
! !USES:
!
      USE ERROR_MOD, ONLY : GEOS_CHEM_STOP

#     include "define.h"                         ! C-preprocessor switches
!
! !INPUT PARAMETERS: 
!
      INTEGER,          INTENT(IN) :: ERROR_NUM  ! I/O error from IOSTAT
      INTEGER,          INTENT(IN) :: UNIT       ! Logical unit # for file
      CHARACTER(LEN=*), INTENT(IN) :: LOCATION   ! Descriptive message
!
! !REVISION HISTORY:
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
!  (8 ) Modifications for Linux/IFORT Intel v9 compiler (bmy, 11/2/05)
!  (9 ) Now call IFORT_ERRMSG to get the IFORT error messages (bmy, 11/30/05)
!  (10) Remove support for LINUX_IFC & LINUX_EFC compilers (bmy, 8/4/06)
!  (10) Remove support for SGI & COMPAQ compilers.  Add IBM_XLF switch.
!        (bmy, 7/8/09)
!  20 Nov 2009 - R. Yantosca - Removed commented-out code for SGI, COMPAQ
!  20 Nov 2009 - R. Yantosca - Added ProTeX header
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
! 
      CHARACTER(LEN=10)            :: ERROR_NUMSTR 
      CHARACTER(LEN=255)           :: ERROR_MSG
      CHARACTER(LEN=255)           :: EXPLAIN_CMD

      ! External functions
      CHARACTER(LEN=255), EXTERNAL :: GERROR, IFORT_ERRMSG

      !=================================================================
      ! IOERROR begins here!
      !=================================================================

      ! Fancy output
      WRITE( 6, '(a)' ) REPEAT( '=', 79 )

      ! Write error number, unit, location
      WRITE( 6, 110 ) ERROR_NUM, UNIT, TRIM( LOCATION )
 110  FORMAT( 'GEOS-CHEM I/O ERROR ', i5, ' in file unit ', i5, /, 
     &        'Encountered at routine:location ', a )

#if   defined( LINUX_PGI )

      !=================================================================
      ! For LINUX platform w/ PGI compiler
      ! Call gerror() to get the I/O error msg
      !=================================================================

      ! GERROR returns ERROR_MSG corresponding to ERROR_NUM 
      ERROR_MSG = GERROR()
 
      ! Print error message to std output
      WRITE( 6, 120 ) TRIM( ERROR_MSG )
 120  FORMAT( /, 'Error: ', a )

#elif defined( LINUX_IFORT ) 

      !=================================================================
      ! For LINUX platform w/ IFORT v8/v9 compiler:
      ! Call IFORT_ERRMSG to get the error number and message
      !=================================================================

      ! Get an error msg corresponding to this error number
      ERROR_MSG = IFORT_ERRMSG( ERROR_NUM )

      ! Print error message to std output
      WRITE( 6, 120 ) ERROR_NUM, TRIM( ERROR_MSG )
 120  FORMAT( /, 'Error ', i4, ': ', a )

#elif defined( SPARC ) 

      !=================================================================
      ! For SUN/Sparc platform: call gerror() to get the I/O error msg
      !=================================================================

      ! GERROR returns ERROR_MSG corresponding to ERROR_NUM 
      ERROR_MSG = GERROR()
 
      ! Print error message to std output
      WRITE( 6, 120 ) TRIM( ERROR_MSG )
 120  FORMAT( /, 'Error: ', a )

#elif defined( IBM_AIX ) || defined( IBM_XLF )

      !=================================================================
      ! For IBM/AIX platform: call gerror() to get the I/O error msg
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
      
      END SUBROUTINE IOERROR
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: file_ex_c
!
! !DESCRIPTION: Function FILE\_EX\_C returns TRUE if FILENAME exists or FALSE 
!  otherwise.  This is handled in a platform-independent way.  The argument 
!  is of CHARACTER type.
!\\
!\\
! !INTERFACE:
!
      FUNCTION FILE_EX_C( FILENAME ) RESULT( IT_EXISTS )
!
! !USES:
!
#     include "define.h"
!
! !INPUT PARAMETERS: 
!
      CHARACTER(LEN=*), INTENT(IN) :: FILENAME   ! Name of file or dir to test
!
! !RETURN VALUE:
!
      LOGICAL                      :: IT_EXISTS  ! =T if the file/dir exists
!
! !REMARKS:
!  This routine is overloaded by public interface FILE_EXISTS.
!
! !REVISION HISTORY:
!  23 Mar 2005 - R. Yantosca - Initial version
!  20 Nov 2009 - R. Yantosca - Updated for LINUX/IFORT Intel v9 compiler
!  20 Nov 2009 - R. Yantosca - Added ProTeX header
!EOP
!------------------------------------------------------------------------------
!BOC

#if   defined( COMPAQ )

      !------------------
      ! COMPAQ compiler
      !------------------

      ! Reference external library function
      INTEGER*4, EXTERNAL :: ACCESS

      ! Test whether directory exists for COMPAQ
      IT_EXISTS = ( ACCESS( TRIM( FILENAME ), ' ' ) == 0 )

#else

      !------------------
      ! Other compilers
      !------------------

      ! Test whether directory exists w/ F90 INQUIRE function
      INQUIRE( FILE=TRIM( FILENAME ), EXIST=IT_EXISTS )

#if   defined( LINUX_IFORT ) 
      
      ! Intel IFORT v9 compiler requires use of the DIRECTORY keyword to 
      ! INQUIRE for checking existence of directories.  (bmy, 11/2/05)
      IF ( .not. IT_EXISTS ) THEN
         INQUIRE( DIRECTORY=TRIM( FILENAME ), EXIST=IT_EXISTS )
      ENDIF

#endif

#endif 

      END FUNCTION FILE_EX_C
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: file_ex_i
!
! !DESCRIPTION: Function FILE\_EX\_I returns TRUE if FILENAME exists or FALSE 
!  otherwise.  This is handled in a platform-independent way.  The argument 
!  is of INTEGER type.
!\\
!\\
! !INTERFACE:
!
      FUNCTION FILE_EX_I( IUNIT ) RESULT( IT_EXISTS )
!
! !USES:
!
#     include "define.h"
!
! !INPUT PARAMETERS: 
!
      ! Arguments
      INTEGER, INTENT(IN) :: IUNIT      ! LUN of file to be tested
!
! !RETURN VALUE:
!
      LOGICAL             :: IT_EXISTS  ! =T if the file/dir exists
!
! !REMARKS:
!  This routine is overloaded by public interface FILE_EXISTS.
!
! !REVISION HISTORY:
!  23 Mar 2005 - R. Yantosca - Initial version
!  20 Nov 2009 - R. Yantosca - Added ProTeX header
!EOP
!------------------------------------------------------------------------------
!BOC
      ! Test whether file unit exists w/ F90 INQUIRE function
      INQUIRE( IUNIT, EXIST=IT_EXISTS )

      END FUNCTION FILE_EX_I
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: close_files
!
! !DESCRIPTION: Subroutine CLOSE\_FILES closes files used by GEOS-Chem.  This 
!  should be called only from the end of the "main.f" program. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CLOSE_FILES
!
! !REVISION HISTORY:
!  04 Mar 1998 - R. Yantosca - Initial version
!  27 Jun 2002 - R. Yantosca - Moved into "file_mod.f"  
!  27 Mar 2003 - R. Yantosca - Also close IU_BC
!  20 Jul 2004 - R. Yantosca - Removed obsolete IU_INPUT and IU_INPTR.
!  20 Jul 2004 - R. Yantosca - Also renamed IU_TS to IU_ND48.
!  20 Oct 2005 - R. Yantosca - Also close IU_XT.
!  20 Nov 2009 - R. Yantosca - Added ProTeX header
!  18 Dec 2009 - Aaron van D - Now close files IU_BC_NA, IU_BC_EU, IU_BC_CH   
!EOP
!------------------------------------------------------------------------------
!BOC
      CLOSE( IU_RST     )
      CLOSE( IU_CHEMDAT )
      CLOSE( IU_FASTJ   )
      CLOSE( IU_GEOS    )
      CLOSE( IU_BPCH    )
      CLOSE( IU_ND20    )
      CLOSE( IU_ND48    )
      CLOSE( IU_ND49    )
      CLOSE( IU_ND50    )
      CLOSE( IU_ND51    )
      CLOSE( IU_ND51b   )
      CLOSE( IU_ND52    )
      CLOSE( IU_PLANE   )
      CLOSE( IU_BC      )
      CLOSE( IU_BC_NA   )
      CLOSE( IU_BC_EU   )
      CLOSE( IU_BC_CH   )
      CLOSE( IU_FILE    )
      CLOSE( IU_PH      )
      CLOSE( IU_TP      )
      CLOSE( IU_I6      )
      CLOSE( IU_A6      )
      CLOSE( IU_A3      )
      CLOSE( IU_KZZ     )
      CLOSE( IU_GWET    )
      CLOSE( IU_XT      )
      CLOSE( IU_SMV2LOG )
      CLOSE( IU_DEBUG   )

      END SUBROUTINE CLOSE_FILES
!EOC      
      ! End of module
      END MODULE FILE_MOD
