C $Id: VTRI0.f,v 1.1 2003/06/30 20:26:02 bmy Exp $
      SUBROUTINE VTRI0 ( A,B,C,F,Y,K,irun)
!**** Used for DIFFUSE - KZZ code
!**** implicit none added (bdf 2/24/99)
      IMPLICIT NONE

      INTEGER irun, K, I, L, LM1
      REAL*8 A(irun,K),B(irun,K),C(irun,K),Y(irun,K+1)
      REAL*8 F(irun,K)
C

      DO 9000 I = 1,irun
         A(I,1) = 1. / A(I,1)
 9000 CONTINUE
C
      DO 100 L = 2,K
         LM1 = L - 1
         DO 9002 I = 1,irun
            C(I,L) = C(I,L) * A(I,LM1)
            A(I,L) = 1. / ( A(I,L) - B(I,LM1) * C(I,L) )
            F(I,L) = F(I,L) + F(I,LM1) * C(I,L)
 9002    CONTINUE
 100  CONTINUE
C
      DO 200 L = K,1,-1
         DO 9004 I = 1,irun
            Y(I,L) = (F(I,L) + B(I,L) * Y(I,L+1)) * A(I,L)
 9004    CONTINUE
 200  CONTINUE
C
      RETURN
      END
