! $Id: numair.f,v 1.1 2003/06/30 20:26:01 bmy Exp $
      FUNCTION NUMAIR( ISEASON, LATI )
!
!******************************************************************************
!  Function NUMAIR returns the appropriate index in the STDAIR array for
!  a given latitude and season.  NUMAIR is only needed for SLOW-J photolysis.
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
      INTEGER             :: NUMAIR

      !=================================================================
      ! NUMAIR begins here!
      !=================================================================

      ! Return the appropriate index in the STDAIR array
      ! adjust for latitude (0-15-30-45-60-75)
      NUMAIR = 9
      IF ( LATI > 7  ) NUMAIR = 13
      IF ( LATI > 22 ) NUMAIR = 17
      IF ( LATI > 37 ) NUMAIR = 21
      IF ( LATI > 52 ) NUMAIR = 25
      IF ( LATI > 67 ) NUMAIR = 29

      ! adjust for season
      IF ( ISEASON < 4 ) NUMAIR = NUMAIR + ISEASON

      ! Return to calling program
      END FUNCTION NUMAIR
