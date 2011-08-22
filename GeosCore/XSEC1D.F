C $Id: XSEC1D.f,v 1.1 2009/09/16 14:06:45 bmy Exp $
      FUNCTION XSEC1D(K,TTT)
C-----------------------------------------------------------------------
c  Quantum yields for O3 --> O2 + O(1D) interpolated across 3 temps
C-----------------------------------------------------------------------
#     include "define.h"

      USE CMN_FJ_MOD, ONLY : IPAR, JPAR, LPAR, JPMAX, JPPJ
      USE JV_CMN_MOD

      IMPLICIT NONE
#     include "define.h"

      integer k
      real*8 ttt, flint, xsec1d
      XSEC1D =
     F  FLINT(TTT,TQQ(1,3),TQQ(2,3),TQQ(3,3),Q1D(K,1),Q1D(K,2),Q1D(K,3))
      return
      end
