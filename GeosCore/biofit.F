C $Id: biofit.f,v 1.1 2009/09/16 14:06:40 bmy Exp $
      REAL*8 FUNCTION BIOFIT(COEFF1,XLAI1,SUNCOS1,CFRAC1)

      USE CMN_SIZE_MOD
      USE CMN_DEP_MOD

      IMPLICIT NONE
#     include "define.h"
C===============================================
C Calculate the light correction
C===============================================
C* BIOFIT and SUNPARAM were written by Y.H. Wang.  See comment
C* in subroutine DEPVEL on what these subroutines do.
C*************************************************************
      INTEGER KK
      PARAMETER (KK=4)
      REAL*8 COEFF1(NPOLY),TERM(KK),REALTERM(NPOLY)
      REAL*8 XLAI1,SUNCOS1,CFRAC1
      INTEGER K,K1,K2,K3

      TERM(1)=1.
      TERM(2)=XLAI1
      TERM(3)=SUNCOS1
      TERM(4)=CFRAC1
      CALL SUNPARAM(TERM(2))
      K=0
      DO K3=1,KK
        DO K2=K3,KK
          DO K1=K2,KK
            K=K+1
            REALTERM(K)=TERM(K1)*TERM(K2)*TERM(K3)
          END DO
        END DO
      END DO
      BIOFIT=0
      DO K=1,NPOLY
        BIOFIT=BIOFIT+COEFF1(K)*REALTERM(K)
      END DO
      IF (BIOFIT.LT.0.1) BIOFIT=0.1

      RETURN
      END
