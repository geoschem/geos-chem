! $Id: arsl1k.f,v 1.1 2003/06/30 20:26:07 bmy Exp $
      REAL*8 FUNCTION ARSL1K( AREA, RADIUS, DENAIR, STKCF, STK, SQM )
!
!******************************************************************************
!  Subroutine ARSL1K calculates the 1st-order loss rate of species on 
!  wet aerosol surface. (lwh, jyl, gmg, djj, 7/1/94; bmy, 4/1/03)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) AREA   (REAL*8) : sfc area of wet aerosols/volume of air [cm2/cm3]
!  (2 ) RADIUS (REAL*8) : radius of wet aerosol [cm], order of 0.01-10 um;
!                           note that radius here is Rd, not Ro
!  (3 ) DENAIR (REAL*8) : Density of air [#/cm3]
!  (4 ) STKCF  (REAL*8) : Sticking coefficient [unitless], order of 0.1
!  (5 ) STK    (REAL*8) : Square root of temperature [K]
!  (6 ) SQM    (REAL*8) : Square root of molecular weight [g/mole]
!
!  References:
!  ============================================================================
!  The 1st-order loss rate on wet aerosol (Dentener's Thesis, p. 14)
!  is computed as:
!
!      ARSL1K [1/s] = area / [ radius/dfkg + 4./(stkcf * xmms) ]        
! 
!  where XMMS = Mean molecular speed [cm/s] = sqrt(8R*TK/pi/M) for Maxwell 
!        DFKG = Gas phase diffusion coeff [cm2/s] (order of 0.1)
! 
!  NOTES:
!  (1 ) Updated comments, cosmetic changes (bmy, 4/4/03)
!******************************************************************************
!
      IMPLICIT NONE
      
      ! Arguments
      REAL*8, INTENT(IN) :: STKCF, AREA, RADIUS, STK, SQM, DENAIR

      ! Local variables
      REAL*8             :: DFKG

      !=================================================================
      ! ARSL1K begins here!
      !=================================================================
      IF ( AREA < 0.0 .OR. RADIUS < 0.0 ) THEN

         ! Use default value if AREA or RADIUS is negative
         ARSL1K = 1.D-3

      ELSE

         ! DFKG = Gas phase diffusion coeff [cm2/s] (order of 0.1)
         DFKG  = 9.45D17/DENAIR * STK * SQRT(3.472D-2 + 1.D0/(SQM*SQM))

         ! Compute ARSL1K according to the formula listed above
         ARSL1K = AREA / ( RADIUS/DFKG + 2.749064E-4*SQM/(STKCF*STK) )

      ENDIF

      ! Return to calling program
      END FUNCTION ARSL1K
