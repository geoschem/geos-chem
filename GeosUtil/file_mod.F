!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: file_mod
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
      !----------------------------------------------------------------
      ! In the GEOS-5 GCM, the unit numbers cannot be PARAMETERs.
      ! Instead,  use INQUIREs to find open LUNs at the point of 
      ! request.  References to most IU_* variables have now been
      ! made local.  IU_BPCH is the only LUN that needs to be seen
      ! across several variables.
      !----------------------------------------------------------------

      ! Logical file unit numbers for ...
      INTEGER, PUBLIC :: IU_BPCH      ! "ctm.bpch" 
      INTEGER, PUBLIC :: IU_FILE
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
!  19 Aug 2010 - R. Yantosca - Added IU_CN and IU_A1 parameters for MERRA
!  19 Aug 2010 - R. Yantosca - Remove IU_KZZ
!  29 May 2010 - S. Kim      - Add IU_BC_SE for the SEAC4RS grid
!  06 Aug 2012 - R. Yantosca - Remove several IU_* variables, as these have
!                              now been moved to various other modules
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!EOP
!------------------------------------------------------------------------------
!BOC
      CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
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
!
! !INPUT PARAMETERS: 
!
      INTEGER,          INTENT(IN) :: ERROR_NUM  ! I/O error from IOSTAT
      INTEGER,          INTENT(IN) :: UNIT       ! Logical unit # for file
      CHARACTER(LEN=*), INTENT(IN) :: LOCATION   ! Descriptive message
!
! !REMARKS:
!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!  %%%  NOTE: "LINUX_IFORT" and "LINUX_PGI" ARE CURRENTLY THE ONLY       %%%
!  %%%  SUPPORTED COMPILER OPTIONS.  MOST SYSTEMS NOW USE A UNIX VERSION %%%
!  %%%  BASED ON LINUX, SO OTHER COMPILERS (IBM/AIX, SUN/SPARC, etc.)    %%%
!  %%%  ARE GENERALLY NOT USED ANYMORE.  LEAVE THE OLDER CODE HERE JUST  %%%
!  %%%  IN CASE WE NEED TO REVERT TO IT AGAIN IN THE FUTURE.             %%%
!  %%%     -- Bob Yantosca, 20 Aug 2013                                  %%%
!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
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
!                  GEOS-Chem Global Chemical Transport Model                  !
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
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
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
!                  GEOS-Chem Global Chemical Transport Model                  !
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
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!EOP
!------------------------------------------------------------------------------
!BOC
      ! Test whether file unit exists w/ F90 INQUIRE function
      INQUIRE( IUNIT, EXIST=IT_EXISTS )

      END FUNCTION FILE_EX_I
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
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
!  19 Aug 2010 - R. Yantosca - Remove IU_KZZ  
!  19 Aug 2010 - R. Yantosca - Now close IU_A1
!  29 May 2010 - S. Kim      - Now close IU_BC_SE
!EOP
!------------------------------------------------------------------------------
!BOC
      CLOSE( IU_BPCH    )
      CLOSE( IU_FILE    )

      END SUBROUTINE CLOSE_FILES
!EOC      
      ! End of module
      END MODULE FILE_MOD
