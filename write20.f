! $Id: write20.f,v 1.1 2003/06/30 20:26:02 bmy Exp $
      SUBROUTINE WRITE20( NYMD, YTAU1, YTAU2, PL24H )
!
!******************************************************************************
!  Subroutine WRITE20 saves production and loss rates to disk, where they 
!  will be later read by subroutine CHEMO3. (bey, bmy, 6/9/99, 6/27/02)
!
!  NOTE: Need to change this to binary punch format (bmy, 11/15/02)
!
!  Arguments as input:
!  ==========================================================================
!  (1) NYMD  : Current value of YYMMDD
!  (2) YTAU1 : TAU value (elapsed hours) at start of diagnostic interval 
!  (3) YTAU2 : TAU value (elapsed hours) at start of diagnostic interval 
!  (4) PL24H : Array containing prod/loss rates of O3 (molec/cm3/s)
!
!  NOTES:
!  (1 ) Now use function NYMD_STRING from "time_mod.f" to generate a
!        Y2K compliant string for all data sets. (bmy, 6/22/00)
!  (2 ) Now use IOS /= 0 to trap both I/O errors and EOF. (bmy, 9/13/00) 
!  (3 ) Now reference IU_FILE and IOERROR from "file_mod.f".  Now use IU_FILE
!        instead of IUT as the file unit #.  Updated comments, cosmetic
!        changes. (bmy, 6/27/02)
!  (4 ) Remove reference to obsolete IJLN (bmy, 11/15/02)
!******************************************************************************
!
      ! References to F90 modules
      USE FILE_MOD, ONLY : IU_FILE, IOERROR
      USE TIME_MOD, ONLY : NYMD_STRING

      IMPLICIT NONE

#     include "CMN_SIZE"
#     include "CMN" 
#     include "CMN_SETUP"  

      ! Arguments
      INTEGER, INTENT(IN) :: NYMD
      REAL*8,  INTENT(IN) :: YTAU1, YTAU2 
      REAL*8,  INTENT(IN) :: PL24H(IIPAR, JJPAR, LLTROP, 2)

      ! Local variables
      INTEGER             :: I, I1, I2, J, J1, J2, L, N, IOS

      CHARACTER(LEN=60)   :: FILENAME

      !=================================================================
      ! WRITE20 begins here!
      !=================================================================

      ! Get output ranges for I and J from the "input.geos" file
      I1 = 1
      I2 = IIPAR
      J1 = 1
      J2 = JJPAR
      
      !=================================================================
      ! Open the output file 
      !=================================================================
      FILENAME = '/scratch/bmy/nd20/rate.' // NYMD_STRING( NYMD )

      WRITE(6,222) FILENAME
 222  FORMAT('---- DIAG20 -- open file  : ',a30)

      ! Open the file
      OPEN( IU_FILE, FILE=TRIM( FILENAME ), STATUS='UNKNOWN', 
     &               FORM='UNFORMATTED',    IOSTAT=IOS )

      ! Error check
      IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'write20:1' )
                  
      ! Write TAU values at beginning and end of interval to disk
      WRITE( IU_FILE, IOSTAT=IOS ) YTAU1, YTAU2

      ! Error check
      IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'write20:2' )

      !=================================================================
      ! Write production rates to disk
      !=================================================================

      ! Write tracer # 
      N = 1
      WRITE( IU_FILE, IOSTAT=IOS ) N

      ! Error check
      IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'write20:3' )
      
      ! Write production rates (N=1) to disk
      WRITE( IU_FILE, IOSTAT=IOS ) 
     &   ( ( ( PL24H(I,J,L,N), I=I1,I2 ), J=J1,J2 ), L=1,LLTROP )
      IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'write20:4' )

      !=================================================================
      ! Write loss rates to disk
      !=================================================================

      ! Write tracer #
      N = 2
      WRITE( IU_FILE, IOSTAT=IOS ) N

      ! Error check
      IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'write20:5' )
      
      ! Write loss rates (N=2)
      WRITE( IU_FILE, IOSTAT=IOS  )  
     &   ( ( ( PL24H(I,J,L,N), I=I1,I2 ), J=J1,J2 ), L=1,LLTROP )
      
      ! Error check
      IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'write20:6' )

      ! Close file
      CLOSE( IU_FILE )

      ! Return to calling program
      END SUBROUTINE WRITE20



