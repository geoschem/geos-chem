! $Id: NOABS.f,v 1.1 2009/09/16 14:06:46 bmy Exp $
      SUBROUTINE NOABS(XLO3,XLO2,XLRAY,BCAER,RFLECT)
!-----------------------------------------------------------------------
!  Zero out absorption terms to check scattering code. Leave a little
!  Rayleigh to provide a minimal optical depth, and set surface albedo
!  to unity.
!-----------------------------------------------------------------------
      IMPLICIT NONE
      real*8 XLO3,XLO2,XLRAY,BCAER,RFLECT
      XLO3=0.d0
      XLO2=0.d0
      XLRAY=XLRAY*1.d-10
      BCAER=0.d0      
      RFLECT=1.d0
      RETURN
      END
