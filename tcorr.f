! $Id: tcorr.f,v 1.2 2003/07/08 15:29:45 bmy Exp $
      FUNCTION TCORR( TEMP )
!
!******************************************************************************
!  Function TCORR applies the temperature correction for isoprene emissions, 
!  according to Guenther et al.(92) (yhw, 11/15/93; bmy, 4/4/03)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) TEMP (REAL*8) : Temperature [K]
!
!  References:
!  ============================================================================
!  Guenther et al, 1992, ... 
!
!  NOTES:
!  (1 ) Removed DATA statements, replaced w/ F90 syntax.  Updated comments
!        and made cosmetic changes (bmy, 4/4/03)
!******************************************************************************
!
      IMPLICIT NONE

      ! Arguments
      REAL*8, INTENT(IN) :: TEMP

      ! Local variables
      REAL*8, PARAMETER  :: R   = 8.314
      REAL*8, PARAMETER  :: CT1 = 95000.
      REAL*8, PARAMETER  :: CT2 = 230000.
      REAL*8, PARAMETER  :: T1  = 303.
      REAL*8, PARAMETER  :: T3  = 314.
      
      ! Function value
      REAL*8             :: TCORR
      
      !=================================================================
      ! TCORR begins here!
      !=================================================================
      TCORR =
     &     EXP( CT1/(R*T1*TEMP) * (TEMP-T1) ) /
     &     ( 1 + EXP( CT2/(R*T1*TEMP) * (TEMP-T3) ) )

      ! Return to calling program
      END FUNCTION TCORR
