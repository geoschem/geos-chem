! $Id: rdisopt.f,v 1.2 2004/09/21 18:04:16 bmy Exp $
      SUBROUTINE RDISOPT( CONVERT )
!
!******************************************************************************
!  Subroutine RDISOPT reads in the baseline emissions for Isoprene, as
!  a function of Olson land type. (yhw, bmy, 7/6/01, 7/20/04)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) CONVERT (REAL*8) : Base emissions for Isoprene by Olson land type
!                          [atoms C/cm2 leaf/s]
!
!  NOTES:
!  (1 ) Now use F90 syntax.  Use IOERROR to trap I/O errors.  Now read the
!        "isopemis.table" file directly from DATA_DIR.  Updated comments
!        and made cosmetic changes.  CMN_ISOP is not needed. (bmy, 7/6/01)
!  (2 ) Deleted obsolete code from ages past (bmy, 9/4/01)
!  (3 ) Now read the "isopemis.table" file from the DATA_DIR/biogenic_200203/
!        directory (bmy, 3/29/02)
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
!-----------------------------------------------
! Prior to 7/20/04:
!#     include "CMN_SETUP" ! DATA_DIR
!-----------------------------------------------

      ! Arguments
      REAL*8, INTENT(OUT) :: CONVERT(NVEGTYPE)

      ! Local variables
      INTEGER             :: I, J, IOS
      !INTEGER, PARAMETER  :: IUNIT=65
      CHARACTER(LEN=80)   :: DUM
      CHARACTER(LEN=255)  :: FILENAME

      !=================================================================
      ! RDISOPT begins here!
      !=================================================================

      ! Define the file name
      FILENAME = TRIM( DATA_DIR ) // 'biogenic_200203/isopemis.table'

      ! Echo info to stdout
      WRITE( 6, 10 ) TRIM( FILENAME )
 10   FORMAT( '     - RDISOPT: Reading ', a )

      ! Open file
      OPEN( IU_FILE, FILE=TRIM( FILENAME ), STATUS='OLD',
     &               FORM='FORMATTED',      IOSTAT=IOS )
      IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'rdisopt:1' )

      ! Read header line
      READ( IU_FILE, '(a80)', IOSTAT=IOS ) DUM
      IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'rdisopt:2' )

      ! Read base isoprene emissons by landtype [atoms C/cm2 leaf/s]
      DO I = 1, NVEGTYPE
         READ( IU_FILE, *, IOSTAT=IOS ) J, CONVERT(I)
         IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'rdisopt:3' )
      ENDDO

      ! Close file
      CLOSE( IU_FILE ) 

      ! Return to calling program
      END SUBROUTINE RDISOPT
