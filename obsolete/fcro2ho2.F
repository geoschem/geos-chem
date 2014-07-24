!fgap
!based on saunder 2003 k14
      REAL*8 FUNCTION FCRO2HO2( XCARBN )

      IMPLICIT NONE

      ! Arguments
      REAL*8, INTENT(IN) :: XCARBN

      FCRO2HO2 = 1D0-EXP(-0.245D0*XCARBN)
     
      ! Return to calling program
      END FUNCTION FCRO2HO2
