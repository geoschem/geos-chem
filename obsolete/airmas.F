! $Id: airmas.f,v 1.1 2009/09/16 14:06:42 bmy Exp $
      FUNCTION AIRMAS( G, H )
!
!******************************************************************************
!  Function AIRMAS corrects the optical path through the layer for the 
!  curvature of the earth.  This correction is very small and according to 
!  RJS fairly laborious to derive --- Do it someday.  Neglecting this 
!  correction, AIRMAS = 1/GMU. (lwh, jyl, gmg, djj, 1990's; bmy, 4/4/03)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) G (REAL*8) : Cosine of solar zenith angle [unitless]
!  (2 ) H (REAL*8) : Scale height of atmosphere [m]
!
!  NOTES:
!  (1 ) Updated comments, cosmetic changes (bmy, 4/4/03)
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!******************************************************************************
! 
      IMPLICIT NONE

      ! Arguments
      REAL*8, INTENT(IN) :: G, H

      ! Function value
      REAL*8             :: AIRMAS

      !=================================================================
      ! AIRMAS begins here!
      !=================================================================
      AIRMAS = (1.0D0+H)/SQRT(G*G + 2.0D0*H*(1.0D0 -
     &        0.6817D0*EXP(-57.3D0*ABS(G)/SQRT(1.0D0+5500.D0*H))/
     &                                        (1.0D0+0.625D0*H)))

      ! Return to calling program
      END FUNCTION AIRMAS
