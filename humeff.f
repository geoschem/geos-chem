! $Id: humeff.f,v 1.1 2003/06/30 20:26:04 bmy Exp $
      REAL*8 FUNCTION HUMEFF( RH, RU, MWARSL )
!
!******************************************************************************
!  Function HUMEFF calculates the fractional increase in the radius of 
!  the wet aerosol compared to dry aerosol (modified bmy 2/27/02)
!
!  Assume we start with NH4HSO4 (which dissociates to NH4+ and HSO4-).
!  Assume the mole fraction of H2O in the bulk wet aerosol equals R.H.
!  Note that mole fraction of NH4HSO4 = (1-RH)/2, b/c of the dissociation.
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) RH     (REAL*8)  : Relative humidity          [fraction]
!  (2 ) RU     (REAL*8)  : Aerosol density            [g/cm3]
!  (3 ) MWARSL (INTEGER) : Dry aerosol molecular mass [g/mole]
!
!  Function Value:
!  ===========================================================================
!  (4 ) HUMEFF (REAL*8)  : Fractional inrease in radius [dr/r]
!
!  NOTES:
!  (1 ) Now cap at 90% relative humidity.  Updated comments, cosmetic 
!        changes. (bmy, rvm, 2/27/02)
!******************************************************************************
!
      IMPLICIT NONE

      ! Arguments
      REAL*8  :: RH,RU
      INTEGER :: MWARSL

      ! Local variables
      REAL*8  :: MWH2O
      REAL*8  :: VH2O,VDRY,VOLRAT

      !=================================================================
      ! HUMEFF begins here!
      !=================================================================

      ! molecular mass of H2O (g/mole)
      MWH2O = 18.0

      ! Cap at 90% RH (rvm, bmy, 2/27/02)
      IF ( RH > 0.90 ) RH = 0.90
	
      ! First find volume ratio (V+dV)/V, which equals 1 + V(H2O)/V(NH4HSO4).
      ! Volume fraction = mole fraction * MW / density
      VH2O   = RH *           MWH2O
      VDRY   = (1.0-RH)*0.5 * MWARSL / RU
      VOLRAT = (1. + VH2O/VDRY)

      ! Then convert to radius ratio, and subtract 1 to get dr/r.
      HUMEFF = VOLRAT**0.333333333 - 1.

      ! Return to calling program
      END FUNCTION HUMEFF
