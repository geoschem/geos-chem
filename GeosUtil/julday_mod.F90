!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: julday_mod.F90
!
! !DESCRIPTION: Module JULDAY\_MOD contains routines used to convert from
!  month/day/year to Astronomical Julian Date and back again.
!\\
!\\
! !INTERFACE:
!
MODULE JULDAY_MOD
!
! !USES:
!
  USE PRECISION_MOD    ! For GEOS-Chem Precision (fp)

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC  :: JULDAY
  PUBLIC  :: CALDATE
!
! !PRIVATE MEMBER FUNCTIONS:
!
  PRIVATE :: MINT
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: JulDay
!
! !DESCRIPTION: Function JULDAY returns the astronomical Julian day.
!\\
!\\
! !INTERFACE:
!
  FUNCTION JULDAY( YYYY, MM, DD ) RESULT( JULIANDAY )
!
! !INPUT PARAMETERS:
!
    INTEGER, INTENT(IN) :: YYYY        ! Year (must be in 4-digit format!)
    INTEGER, INTENT(IN) :: MM          ! Month (1-12)
    REAL*8,  INTENT(IN) :: DD          ! Day of month (may be fractional!)
!
! !RETURN VALUE:
!
    REAL*8              :: JULIANDAY   ! Astronomical Julian Date
!
! !REMARKS:
!  (1) Algorithm taken from "Practical Astronomy With Your Calculator",
!       Third Edition, by Peter Duffett-Smith, Cambridge UP, 1992.
!  (2) Requires the external function MINT.F.
!  (3) JulDay will compute the correct Julian day for any BC or AD date.
!  (4) For BC dates, subtract 1 from the year and append a minus sign.
!       For example, 1 BC is 0, 2 BC is -1, etc.  This is necessary for
!       the algorithm.
!
! !REVISION HISTORY:
!  26 Nov 2001 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
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
    X1 = DBLE( YEAR1 ) / 100.0d0
    A  = MINT( X1 )

    ! The Gregorian calendar begins on 10 October 1582
    ! Any dates prior to this will be in the Julian calendar
    IF ( YYYY > 1582 ) THEN
       ISGREGORIAN = .TRUE.
    ELSE
       IF ( ( YYYY   == 1582 )  .AND. &
            ( MONTH1 >= 10   )  .AND. &
            ( DD     >= 15.0 ) ) THEN
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

  END FUNCTION JULDAY
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Mint
!
! !DESCRIPTION: Function MINT is the modified integer function.
!\\
!\\
! !INTERFACE:
!
  FUNCTION MINT( X ) RESULT ( VALUE )
!
! !INPUT PARAMETERS:
!
    REAL*8, INTENT(IN) :: X
!
! !RETURN VALUE:
!
    REAL*8             :: VALUE
!
! !REMARKS:
!  The modified integer function is defined as follows:
!
!            { -INT( ABS( X ) )   for X < 0
!     MINT = {
!            {  INT( ABS( X ) )   for X >= 0
!
! !REVISION HISTORY:
!  20 Nov 2001 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
    IF ( X < 0d0 ) THEN
       VALUE = -INT( ABS( X ) )
    ELSE
       VALUE =  INT( ABS( X ) )
    ENDIF

  END FUNCTION MINT
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: CalDate
!
! !DESCRIPTION: Subroutine CALDATE converts an astronomical Julian day to
!  the YYYYMMDD and HHMMSS format.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE CALDATE( JULIANDAY, YYYYMMDD, HHMMSS )
!
! !INPUT PARAMETERS:
!
    REAL*8,  INTENT(IN)  :: JULIANDAY  ! Astronomical Julian Date
!
! !OUTPUT PARAMETERS:
!
    INTEGER, INTENT(OUT) :: YYYYMMDD   ! Date in YYYY/MM/DD format
    INTEGER, INTENT(OUT) :: HHMMSS     ! Time in hh:mm:ss format
!
! !REMARKS:
!   Algorithm taken from "Practical Astronomy With Your Calculator",
!   Third Edition, by Peter Duffett-Smith, Cambridge UP, 1992.
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
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
    !------------------------------------------------------------------
    ! NOTE: Some times (like 40min = 0.6666 hrs) will cause a roundoff
    ! error that will make the minutes eg. 39.9999 instead of 40.
    ! For now put in a kludge to rectify this situation.
    IF ( INT(SS) == 59 ) THEN
       SS = 0.0e0
       MM = NINT( MM )
    ENDIF
    !---------------------------------------------------------------
    HHMMSS = ( INT( HH ) * 10000 ) + ( INT( MM ) * 100 ) + INT( SS )

  END SUBROUTINE CALDATE
!EOC
END MODULE JULDAY_MOD
