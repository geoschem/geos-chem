!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: file_mod.F90
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
  END INTERFACE FILE_EXISTS
!
! !PRIVATE MEMBER FUNCTIONS:
!
  PRIVATE :: FILE_EX_C
  PRIVATE :: FILE_EX_I
!
! !REVISION HISTORY:
! See https://github.com/geoschem/geos-chem for complete history
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
! !IROUTINE: IoError
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
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
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
110 FORMAT( 'GEOS-CHEM I/O ERROR ', i5, ' in file unit ', i5, /, &
            'Encountered at routine:location ', a )

#ifdef LINUX_IFORT
    !=================================================================
    ! For LINUX platform w/ IFORT v8/v9 compiler:
    ! Call IFORT_ERRMSG to get the error number and message
    !=================================================================

    ! Get an error msg corresponding to this error number
    ERROR_MSG = IFORT_ERRMSG( ERROR_NUM )

    ! Print error message to std output
    WRITE( 6, 120 ) ERROR_NUM, TRIM( ERROR_MSG )
120 FORMAT( /, 'Error ', i4, ': ', a )
#endif

    ! Fancy output
    WRITE( 6, '(a)' ) REPEAT( '=', 79 )

    ! Deallocate arrays and stop safely
    CALL GEOS_CHEM_STOP

  END SUBROUTINE IOERROR
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: File_Ex_C
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
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC

    ! Test whether directory exists w/ F90 INQUIRE function
    INQUIRE( FILE=TRIM( FILENAME ), EXIST=IT_EXISTS )

#ifdef LINUX_IFORT
    ! Intel IFORT v9 compiler requires use of the DIRECTORY keyword to 
    ! INQUIRE for checking existence of directories.  (bmy, 11/2/05)
    IF ( .not. IT_EXISTS ) THEN
       INQUIRE( DIRECTORY=TRIM( FILENAME ), EXIST=IT_EXISTS )
    ENDIF
#endif

  END FUNCTION FILE_EX_C
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: File_Ex_I
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
!  See https://github.com/geoschem/geos-chem for complete history
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
! !IROUTINE: Close_Files
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
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
    CLOSE( IU_BPCH    )
    CLOSE( IU_FILE    )

  END SUBROUTINE CLOSE_FILES
!EOC
END MODULE FILE_MOD
