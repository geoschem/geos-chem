! $Id: findmon.f,v 1.1 2009/09/16 14:06:29 bmy Exp $
      SUBROUTINE FINDMON( JDAY, INMONTH, MM, STARTDAY )
!
!******************************************************************************
!  Function FINDMON finds which month JDAY (day of this year) is in.  
!  FINDMON is called by the Leaf Area Index routine rdlai.f.
!  (yhw, gmg, djj, 1994; bmy, 4/4/03)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) JDAY     (INTEGER) : Current day of year (0-365 or 0-366, leap years)
!  (2 ) INMONTH  (INTEGER) : Current month (1-12)
!  (4 ) STARTDAY (INTEGER) : Array of starting days for LAI monthly data
!
!  Arguments as Output:
!  ============================================================================
!  (3 ) MM       (INTEGER) : Output month number (1-12)
!
!  NOTES:
!  (1 ) Updated comments, cosmetic changes (bmy, 4/4/03)
!******************************************************************************
!
      IMPLICIT NONE

      ! Arguments
      INTEGER, INTENT(IN)  :: JDAY, INMONTH, STARTDAY(13)
      INTEGER, INTENT(OUT) :: MM

      !=================================================================
      ! FINDMON begins here!
      !=================================================================
      IF ( JDAY < STARTDAY(1) ) THEN
         MM = 12
      ELSE IF ( JDAY < STARTDAY(INMONTH) ) THEN
         MM = INMONTH-1
      ELSE
         MM = INMONTH
      ENDIF

      ! Return to calling program
      END SUBROUTINE FINDMON
