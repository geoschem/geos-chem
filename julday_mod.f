! $Id: julday_mod.f,v 1.1 2003/06/30 20:26:07 bmy Exp $
      MODULE JULDAY_MOD
!
!******************************************************************************
!  Module JULDAY_MOD contains routines used to convert from month/day/year
!  to Astronomical Julian Date and back again. (bmy, 11/26/01, 6/26/02)
!
!  Module Routines:
!  ============================================================================
!  (1 ) JULDAY  : Given month/day/year, computes astronomical Julian date
!  (2 ) MINT    : Modified integer function, used by JULDAY
!  (3 ) CALDATE : Given astronomical julian date, computes YYYYMMDD, HHMMSS
!
!  NOTES:
!  (1 ) Moved JULDAY, MINT, CALDATE here from "bpch2_mod.f" (bmy, 11/20/01)
!  (2 ) Bug fix: now compute NHMS correctly.  Also use REAL*4 variables to
!        avoid roundoff errors. (bmy, 11/26/01)
!  (3 ) Updated comments (bmy, 5/28/02)
!  (4 ) Renamed arguments for clarity (bmy, 6/26/02)
!******************************************************************************
!
      IMPLICIT NONE

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement 
      !=================================================================
      CONTAINS

!------------------------------------------------------------------------------
     
      FUNCTION JULDAY( YYYY, MM, DD ) RESULT( JULIANDAY )
!
!******************************************************************************
!  Function JULDAY returns the astronomical Julian day (bmy, 11/26/01, 6/26/02)
!
!  Algorithm taken from "Practical Astronomy With Your Calculator",
!  Third Edition, by Peter Duffett-Smith, Cambridge UP, 1992.
! 
!  Arguments as Input:
!  ===========================================================================
!  (1 ) YYYY  (INTEGER) : Current year 
!  (2 ) MM    (INTEGER) : Current month
!  (3 ) DD    (REAL*8 ) : Current day (can be fractional, e.g. 17.25)
!
!  NOTES:
!  (1 ) JULDAY requires the external function MINT.F.
!  (2 ) JULDAY will compute the correct Julian day for any 
!        BC or AD date.
!  (3 ) For BC dates, subtract 1 from the year and append a minus 
!        sign.  For example, 1 BC is 0, 2 BC is -1, etc.  This is 
!        necessary for the algorithm.  
!  (4 ) Changed YEAR to YYYY, MONTH to MM, and DAY to DD for documentation
!        purposes. (bmy, 6/26/02)
!******************************************************************************
!   
      ! Arguments
      INTEGER, INTENT(IN) :: YYYY, MM
      REAL*8              :: DD,   JULIANDAY
   
      ! Local variables
      INTEGER             :: YEAR1, MONTH1
      REAL*8              :: X1, A, B, C, D
      LOGICAL             :: ISGREGORIAN
   
      !==================================================================
      ! JULDAY begins here!
      !
      ! Follow algorithm from Peter Duffett-Smith (1992)
      !==================================================================
   
      ! Compute YEAR and MONTH1
      IF ( ( MM == 1 ) .OR. ( MM == 2 ) ) THEN
         YEAR1  = YYYY  - 1
         MONTH1 = MM    + 12 
      ELSE
         YEAR1  = YYYY
         MONTH1 = MM
      ENDIF
   
      ! Compute the "A" term. 
      X1 = DBLE( YYYY ) / 100.0d0
      A  = MINT( X1 )
   
      ! The Gregorian calendar begins on 10 October 1582
      ! Any dates prior to this will be in the Julian calendar
      IF ( YYYY > 1582 ) THEN
         ISGREGORIAN = .TRUE.
      ELSE
         IF ( ( YYYY   == 1582 )  .AND. 
     &        ( MONTH1 >= 10   )  .AND. 
     &        ( DD     >= 15.0 ) ) THEN 
            ISGREGORIAN = .TRUE.
         ELSE
            ISGREGORIAN = .FALSE.
         ENDIF
      ENDIF
            
      ! Compute the "B" term according to Gregorian or Julian calendar
      IF ( ISGREGORIAN ) THEN
         B = 2.0d0 - A + MINT( A / 4.0d0 )
      ELSE
         B = 0.0d0
      ENDIF
   
      ! Compute the "C" term for BC dates (YEAR1 <= 0 ) 
      ! or AD dates (YEAR1 > 0)
      IF ( YEAR1 < 0 ) THEN
         X1 = ( 365.25d0 * YEAR1 ) - 0.75d0
         C  = MINT( X1 )
      ELSE
         X1 = 365.25d0 * YEAR1
         C  = MINT( X1 ) 
      ENDIF
   
      ! Compute the "D" term    
      X1 = 30.6001d0 * DBLE( MONTH1 + 1 )
      D  = MINT( X1 )
   
      ! Add the terms to get the Julian Day number 
      JULIANDAY = B + C + D + DD + 1720994.5d0
   
      ! Return to calling program
      END FUNCTION JULDAY
   
!------------------------------------------------------------------------------
   
      FUNCTION MINT( X ) RESULT ( VALUE )
!
!******************************************************************************
!  Function MINT (bmy, 11/20/01) defines the MODIFIED INTEGER FUNCTION:
! 
!  MINT = -INT( ABS( X ) ), X <  0
!  MINT =  INT( ABS( X ) ), X >= 0
! 
!  Arguments as Input:
!  ============================================================================
!  (1) X : (REAL*8) Argument for the function MINT
!
!  NOTES:
!  (1) MINT is primarily intended for use with routine JULDAY.
!******************************************************************************
!
      ! Arguments
      REAL*8, INTENT(IN) :: X
         
      ! Return value
      REAL*8             :: value
   
      !=================================================================
      ! MINT begins here!
      !=================================================================
      IF ( X < 0d0 ) THEN 
         VALUE = -INT( ABS( X ) )        
      ELSE
         VALUE =  INT( ABS( X ) )        
      ENDIF
   
      ! Return to calling program
      END FUNCTION MINT
   
!------------------------------------------------------------------------------
   
      SUBROUTINE CALDATE( JULIANDAY, YYYYMMDD, HHMMSS )
!
!******************************************************************************
!  Subroutine CALDATE converts an astronomical Julian day to a YYYYMMDD
!  (year-month-day) and HHMMSS (hour-min-sec) format. (bmy, 11/26/01, 6/26/02)
!   
!  Algorithm taken from "Practical Astronomy With Your Calculator",
!  Third Edition, by Peter Duffett-Smith, Cambridge UP, 1992.
!   
!  Arguments as Input:
!  ============================================================================
!  (1 ) JULIANDAY : REAL*8  : Astronomical julian day
!   
!  Arguments as Output:
!  ============================================================================
!  (1 ) YYYYMMDD (INTEGER) : year-month-day corresponding to JULIANDAY
!  (2 ) HHMMSS   (INTEGER) : hour-min-sec   corresponding to JULIANDAY
!
!  NOTES:
!  (1 ) Now compute HHMMSS correctly.  Also use REAL*4 variables HH, MM, SS
!        to avoid roundoff errors. (bmy, 11/21/01)
!  (2 ) Renamed NYMD to YYYYMMDD and NHMS to HHMMSS for documentation
!        purposes (bmy, 6/26/02)
!******************************************************************************
!
      ! Arguments
      REAL*8,  INTENT(IN)  :: JULIANDAY
      INTEGER, INTENT(OUT) :: YYYYMMDD, HHMMSS
    
      ! Local variables
      REAL*4               :: HH, MM, SS
      REAL*8               :: A, B, C, D, DAY, E, F 
      REAL*8               :: FDAY, G, I, J, JD, M, Y

      !=================================================================
      ! CALDATE begins here!
      ! See "Practical astronomy with your calculator", Peter Duffett-
      ! Smith 1992, for an explanation of the following algorithm.
      !=================================================================
      JD = JULIANDAY + 0.5d0
      I  = INT( JD )
      F  = JD - INT( I )
    
      IF ( I > 2299160d0 ) THEN
         A = INT( ( I - 1867216.25d0 ) / 36524.25d0 )
         B = I + 1 + A - INT( A / 4 )
      ELSE
         B = I
      ENDIF
    
      C = B + 1524d0
      
      D = INT( ( C - 122.1d0 ) / 365.25d0 )
    
      E = INT( 365.25d0 * D )
    
      G = INT( ( C - E ) / 30.6001d0 )
   
      ! DAY is the day number
      DAY  = C - E + F - INT( 30.6001d0 * G ) 
    
      ! FDAY is the fractional day number
      FDAY = DAY - INT( DAY )
      
      ! M is the month number
      IF ( G < 13.5d0 ) THEN
         M = G - 1d0
      ELSE
         M = G - 13d0
      ENDIF
    
      ! Y is the year number
      IF ( M > 2.5d0 ) THEN
         Y = D - 4716d0
      ELSE
         Y = D - 4715d0
      ENDIF
    
      ! Year-month-day value
      YYYYMMDD = ( INT( Y ) * 10000 ) + ( INT( M ) * 100 ) + INT( DAY )
      
      ! Hour-minute-second value
      ! NOTE: HH, MM, SS are REAL*4 to avoid numerical roundoff errors
      HH     = FDAY * 24d0 
      MM     = ( HH - INT( HH ) ) * 60d0
      SS     = ( MM - INT( MM ) ) * 60d0   
      HHMMSS = ( INT( HH ) * 10000 ) + ( INT( MM ) * 100 ) + INT( SS )
      
      ! Return to calling program
      END SUBROUTINE CALDATE

!------------------------------------------------------------------------------

      END MODULE JULDAY_MOD
