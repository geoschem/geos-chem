! $Id: numo3.f,v 1.1 2003/06/30 20:26:08 bmy Exp $
      FUNCTION NUMO3( ISEASON, LATI )
!
!******************************************************************************
!  Function NUMO3 returns the appropriate index in the STDO3 array for a
!  given latitude and season.  NUMO3 is only needed for SLOW-J photolysis.
!  (lwh, jyl, gmg, djj, 1990's; bmy, 4/4/03)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) ISEASON (INTEGER) : Season index (1-4)
!  (2 ) LATI    (INTEGER) : Grid box latitude [degrees]
!
!  NOTES:
!  (1 ) Updated comments, cosmetic changes (bmy, 4/4/03)
!******************************************************************************
!
      IMPLICIT NONE

      ! Arguments
      INTEGER, INTENT(IN) :: ISEASON, LATI
      
      ! Function value
      INTEGER             :: NUMO3

      !=================================================================
      ! NUMO3 begins here!
      !=================================================================

      ! Return the appropriate index in the STDO3 array
      ! Adjust for latitude (9-30-45)
      NUMO3 = 5
      IF ( LATI < 20 ) NUMO3 = 3
      IF ( LATI > 37 ) NUMO3 = 7

      ! Adjust for season (except for 30N, where we only have one profile)
      IF ( NUMO3 /= 5 ) THEN
         IF ( ISEASON == 2 ) NUMO3 = NUMO3 - 1
         IF ( ISEASON == 4 ) NUMO3 = NUMO3 + 1
      ENDIF

      ! Return to calling program
      END FUNCTION NUMO3
