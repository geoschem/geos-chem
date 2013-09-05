! $Id: XSECO3.f,v 1.1 2009/09/16 14:06:45 bmy Exp $
      FUNCTION XSECO3(K,TTT)
!-----------------------------------------------------------------------
c  Cross-sections for O3 for all processes interpolated across 3 temps
!-----------------------------------------------------------------------

      USE JV_CMN_MOD

      IMPLICIT NONE

      integer k
      real*8 ttt, flint, xseco3
      XSECO3  = 
     F  FLINT(TTT,TQQ(1,2),TQQ(2,2),TQQ(3,2),QO3(K,1),QO3(K,2),QO3(K,3))
      return
      end
