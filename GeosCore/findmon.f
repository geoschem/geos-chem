! $Id: findmon.f,v 1.2 2009/11/30 19:57:57 ccarouge Exp $
      SUBROUTINE FINDMON( JDAY, INMONTH, INYEAR, MM, YYYY, STARTDAY )
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
!  (3 ) INYEAR   (INTEGER) : Current simulation year
!  (6 ) STARTDAY (INTEGER) : Array of starting days for LAI monthly data
!
!  Arguments as Output:
!  ============================================================================
!  (4 ) MM       (INTEGER) : Output month number (1-12)
!  (5 ) YYYY     (INTEGER) : Output correct year for the LAI data 
!
!  NOTES:
!  (1 ) Updated comments, cosmetic changes (bmy, 4/4/03)
!  (2 ) Add the current simulation year as input & the current LAI as output.
!       This is necessary for reading in MODIS LAI (mpb,2009).
!******************************************************************************
!
      IMPLICIT NONE

      ! Arguments
      INTEGER, INTENT(IN)  :: JDAY, INMONTH, STARTDAY(13)
      INTEGER, INTENT(IN)  :: INYEAR ! (mpb,2008)
      INTEGER, INTENT(OUT) :: MM
      INTEGER, INTENT(OUT) :: YYYY   ! (mpb,2008)

      !=================================================================
      ! FINDMON begins here!
      !=================================================================
      IF ( JDAY < STARTDAY(1) ) THEN
         MM = 12
         YYYY = INYEAR - 1      ! (mpb,2008)
      ELSE IF ( JDAY < STARTDAY(INMONTH) ) THEN
         MM = INMONTH-1
         YYYY = INYEAR           ! (mpb,2008)
      ELSE
         MM = INMONTH
         YYYY = INYEAR           ! (mpb,2008)
      ENDIF

      ! Return to calling program
      END SUBROUTINE FINDMON
