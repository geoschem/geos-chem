!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: rdlight
!
! !DESCRIPTION: Subroutine RDLIGHT reads the polynomial coefficients for 
!  isoprene emissions from disk.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE RDLIGHT
!
! !USES:
!
      USE DIRECTORY_MOD, ONLY : DATA_DIR 
      USE FILE_MOD,      ONLY : IU_FILE, IOERROR

      USE CMN_SIZE_MOD  ! Size parameters
      USE CMN_ISOP_MOD  ! SOPCOEFF

      IMPLICIT NONE
#     include "define.h"
#     include "define.h"

! 
! !REVISION HISTORY: 
!  06 Jul 2001 - Y. H. Wang, R. Yantosca - Initial version
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
!  02 Dec 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
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

      END SUBROUTINE RDLIGHT
!EOC
