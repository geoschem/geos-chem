! $Id: rdmonot.f,v 1.2 2004/09/21 18:04:17 bmy Exp $
      SUBROUTINE RDMONOT( GMONOT )
!
!******************************************************************************
!  Subroutine RDMONOT reads baseline monoterpene emission values from 
!  Guenther et al. (1995), as a function of Olson landtype area.
!  (bdf, bmy, 7/6/01, 7/20/04) 
!
!  Arguments as Output:
!  ============================================================================
!  (1 ) GMONOT: Monoterpene emissions for each landtype [atoms C/cm2 leaf/s]
!
!  NOTES:
!  (1 ) Now read updated file "monotemis.v4-13.table" (bdf, bmy, 6/6/01)
!  (2 ) Now reference DATA_DIR from "CMN_SETUP. (bmy, 6/6/01)
!  (3 ) Now use IOERROR to trap I/O errors (bmy, 6/6/01)
!  (4 ) IUNIT=65 is now a parameter (bmy, 7/6/01)
!  (5 ) Now read file "monotemis.v4-13.table" from the 
!        DATA_DIR/biogenic_200203 directory (bmy, 3/29/02)
!  (6 ) Removed obsolete code from March 2002.  Now reference IU_FILE and 
!        IOERROR from "file_mod.f".  Now use IU_FILE as the file unit number
!        instead of IUNIT. (bmy, 6/27/02)
!  (7 ) Now references DATA_DIR from "directory_mod.f" (bmy, 7/20/04)
!******************************************************************************
!
      ! References to F90 modules
      USE DIRECTORY_MOD, ONLY : DATA_DIR
      USE FILE_MOD,      ONLY : IU_FILE, IOERROR

      IMPLICIT NONE

#     include "CMN_SIZE"  ! Size parameters
!----------------------------------------------
! Prior to 7/20/04:
!#     include "CMN_SETUP" ! DATA_DIR 
!----------------------------------------------

      ! Arguments
      REAL*8, INTENT(OUT) :: GMONOT(NVEGTYPE)

      ! Local variables
      INTEGER             :: N, T, IOS
      CHARACTER(LEN=255)  :: FILENAME

      !=================================================================
      ! RDMONOT begins here!
      !=================================================================

      ! Monoterpene file name
      FILENAME = TRIM( DATA_DIR ) // 
     &           'biogenic_200203/monotemis.v4-13.table'

      ! Echo output
      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( '     - RDMONOT: Reading ', a )

      ! Open file
      OPEN( IU_FILE, FILE=TRIM( FILENAME ), STATUS='OLD', 
     &               FORM='FORMATTED',      IOSTAT=IOS )
      IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'rdmonot:1' )

      ! Read header line
      READ( IU_FILE, * ) 
      IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'rdmonot:2' )

      ! Loop over vegetation types and read emissions [atoms C/cm2 leaf/s]
      DO N = 1, NVEGTYPE
         READ( IU_FILE, *, IOSTAT=IOS ) T, GMONOT(N)
         IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'rdmonot:3' )
      ENDDO

      ! Close IU_FILE
      CLOSE( IU_FILE )

      ! Return to calling program
      END SUBROUTINE RDMONOT
