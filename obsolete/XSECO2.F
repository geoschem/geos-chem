! $Id: XSECO2.f,v 1.1 2009/09/16 14:06:45 bmy Exp $
      FUNCTION XSECO2(K,TTT)
!-----------------------------------------------------------------------
c  Cross-sections for O2 interpolated across 3 temps; No S_R Bands yet!
!-----------------------------------------------------------------------
      
      USE JV_CMN_MOD

      IMPLICIT NONE

      integer k
      real*8 ttt, flint, xseco2
      XSECO2 =
     F  FLINT(TTT,TQQ(1,1),TQQ(2,1),TQQ(3,1),QO2(K,1),QO2(K,2),QO2(K,3))
      return
      end
