! $Id: rdlight.f,v 1.2 2004/09/21 18:04:17 bmy Exp $
      SUBROUTINE RDLIGHT
!
!******************************************************************************
!  Subroutine RDLIGHT reads the polynomial coefficients for isoprene
!  emissions from disk. (yhw, bmy, 7/6/01, 7/20/04)
!
!  NOTES:
!  (1 ) Now use F90 syntax.  Now reads the file "light.table" directly
!        from DATA_DIR so that symbolic links are unnecessary.  Also use 
!        IOERROR to trap I/O errors.  Updated comments and made cosmetic 
!        changes (bmy, 7/6/01)
!  (2 ) Deleted obsolete code from ages ago.  Also print full pathname
!        of the "light.table" file. (bmy, 9/4/01)
!  (3 ) Now read file "light.table" from the DATA_DIR/biogenic_200203/ 
!        directory.  Added FILENAME variable. (bmy, 3/29/02)
!  (4 ) Deleted obsolete code from March 2002.  Now reference IU_FILE and
!        IOERROR from "file_mod.f".  Now use IU_FILE instead of IUNIT as
!        the file unit number. (bmy, 6/27/02)
!  (5 ) Now references DATA_DIR from "directory_mod.f" (bmy, 7/20/04)
!******************************************************************************
!
      ! References to F90 modules
      USE DIRECTORY_MOD, ONLY : DATA_DIR 
      USE FILE_MOD,      ONLY : IU_FILE, IOERROR

      IMPLICIT NONE

#     include "CMN_SIZE"  ! Size parameters
#     include "CMN_ISOP"  ! SOPCOEFF
!--------------------------------------------
! Prior to 7/20/04:
!#     include "CMN_SETUP" ! DATA_DIR 
!--------------------------------------------

      INTEGER            :: I, IOS
      CHARACTER(LEN=80)  :: DUM
      CHARACTER(LEN=255) :: FILENAME

      !=================================================================
      ! RDLIGHT begins here!
      !=================================================================
      
      ! File containing polynomial data
      FILENAME = TRIM( DATA_DIR ) // 'biogenic_200203/light.table' 

      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( '     - RDLIGHT: Reading ', a )

      ! Open the "light.table" file in DATA_DIR/biogenic_200203/
      OPEN( IU_FILE, FILE=TRIM( FILENAME ), STATUS='OLD', IOSTAT=IOS ) 
      IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'rdlight:1' )

      ! Read header line
      READ( IU_FILE, '(a80)', IOSTAT=IOS ) DUM
      IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'rdlight:2' )
      
      ! Read data
      READ( IU_FILE,'(8(1PE10.2))',IOSTAT=IOS ) (SOPCOEFF(I), I=1,NPOLY)
      IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'rdlight:3' )

      ! Close file
      CLOSE( IU_FILE )

      ! Return to calling program
      END SUBROUTINE RDLIGHT
