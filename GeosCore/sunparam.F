C $Id: sunparam.f,v 1.1 2009/09/16 14:06:02 bmy Exp $
      SUBROUTINE SUNPARAM(X)

      IMPLICIT NONE
#     include "define.h"
#     include "define.h"

C===============================================
C the sequence is lai,suncos,cloud fraction
C===============================================
C  NN = number of variables (lai,suncos,cloud fraction)
      INTEGER NN
      PARAMETER(NN=3)
C  ND = scaling factor for each variable
      INTEGER ND(NN),I
      DATA ND /55,20,11/
C  X0 = maximum for each variable
      REAL*8 X(NN),X0(NN),XLOW
      DATA X0 /11.,1.,1./

      DO I=1,NN
        X(I)=MIN(X(I),X0(I))
C XLOW = minimum for each variable
        IF (I.NE.3) THEN
          XLOW=X0(I)/REAL(ND(I))
        ELSE
          XLOW= 0.
        END IF
        X(I)=MAX(X(I),XLOW)
        X(I)=X(I)/X0(I)
      END DO

      RETURN
      END
