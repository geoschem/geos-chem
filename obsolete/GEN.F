! $Id: GEN.f,v 1.1 2009/09/16 14:06:48 bmy Exp $
      SUBROUTINE GEN(ID)
!-----------------------------------------------------------------------
!  Generates coefficient matrices for the block tri-diagonal system:
!               A(I)*X(I-1) + B(I)*X(I) + C(I)*X(I+1) = H(I)
!-----------------------------------------------------------------------
      USE JV_MIE_MOD

      IMPLICIT NONE


      integer id, id0, id1, im, i, j, k, mstart
      real*8  sum0, sum1, sum2, sum3
      real*8  deltau, d1, d2, surfac
!---------------------------------------------
      IF(ID.EQ.1 .OR. ID.EQ.ND) THEN
!---------calculate generic 2nd-order terms for boundaries
       ID0 = ID
       ID1 = ID+1
       IF(ID.GE.ND) ID1 = ID-1
       DO 10 I=1,N
          SUM0 = 0.0d0
          SUM1 = 0.0d0
          SUM2 = 0.0d0
          SUM3 = 0.0d0
        DO IM=M,MFIT,2
          SUM0 = SUM0 + POMEGA(IM,ID0)*PM(I,IM)*PM0(IM)
          SUM2 = SUM2 + POMEGA(IM,ID1)*PM(I,IM)*PM0(IM)
        ENDDO
        DO IM=M+1,MFIT,2
          SUM1 = SUM1 + POMEGA(IM,ID0)*PM(I,IM)*PM0(IM)
          SUM3 = SUM3 + POMEGA(IM,ID1)*PM(I,IM)*PM0(IM)
        ENDDO
         H(I) = 0.5d0*(SUM0*FZ(ID0) + SUM2*FZ(ID1))
         A(I) = 0.5d0*(SUM1*FZ(ID0) + SUM3*FZ(ID1))
        DO J=1,I
          SUM0 = 0.0d0
          SUM1 = 0.0d0
          SUM2 = 0.0d0
          SUM3 = 0.0d0
         DO IM=M,MFIT,2
          SUM0 = SUM0 + POMEGA(IM,ID0)*PM(I,IM)*PM(J,IM)
          SUM2 = SUM2 + POMEGA(IM,ID1)*PM(I,IM)*PM(J,IM)
         ENDDO
         DO IM=M+1,MFIT,2
          SUM1 = SUM1 + POMEGA(IM,ID0)*PM(I,IM)*PM(J,IM)
          SUM3 = SUM3 + POMEGA(IM,ID1)*PM(I,IM)*PM(J,IM)
         ENDDO
         S(I,J) = - SUM2*WT(J)
         S(J,I) = - SUM2*WT(I)
         W(I,J) = - SUM1*WT(J)
         W(J,I) = - SUM1*WT(I)
         U1(I,J) = - SUM3*WT(J)
         U1(J,I) = - SUM3*WT(I)
          SUM0 = 0.5d0*(SUM0 + SUM2)
         B(I,J) = - SUM0*WT(J)
         B(J,I) = - SUM0*WT(I)
        ENDDO
         S(I,I) = S(I,I) + 1.0d0
         W(I,I) = W(I,I) + 1.0d0
         U1(I,I) = U1(I,I) + 1.0d0
         B(I,I) = B(I,I) + 1.0d0
   10  CONTINUE
       DO I=1,N
         SUM0 = 0.0d0
        DO J=1,N
         SUM0 = SUM0 + S(I,J)*A(J)/EMU(J)
        ENDDO
        C1(I) = SUM0
       ENDDO
       DO I=1,N
        DO J=1,N
          SUM0 = 0.0d0
          SUM2 = 0.0d0
         DO K=1,N
          SUM0 = SUM0 + S(J,K)*W(K,I)/EMU(K)
          SUM2 = SUM2 + S(J,K)*U1(K,I)/EMU(K)
         ENDDO
         A(J) = SUM0
         V1(J) = SUM2
        ENDDO
        DO J=1,N
         W(J,I) = A(J)
         U1(J,I) = V1(J)
        ENDDO
       ENDDO
       IF (ID.EQ.1) THEN
!-------------upper boundary, 2nd-order, C-matrix is full (CC)
        DELTAU = ZTAU(2) - ZTAU(1)
        D2 = 0.25d0*DELTAU
        DO I=1,N
          D1 = EMU(I)/DELTAU
          DO J=1,N
           B(I,J) = B(I,J) + D2*W(I,J)
           CC(I,J) = D2*U1(I,J)
          ENDDO
          B(I,I) = B(I,I) + D1
          CC(I,I) = CC(I,I) - D1
!         H(I) = H(I) + 2.0d0*D2*C1(I) + D1*SISOTP
          H(I) = H(I) + 2.0d0*D2*C1(I)
          A(I) = 0.0d0
        ENDDO
       ELSE
!-------------lower boundary, 2nd-order, A-matrix is full (AA)
        DELTAU = ZTAU(ND) - ZTAU(ND-1)
        D2 = 0.25d0*DELTAU
        SURFAC = 4.0d0*ZREFL/(1.0d0 + ZREFL)
        DO I=1,N
          D1 = EMU(I)/DELTAU
          H(I) = H(I) - 2.0d0*D2*C1(I)
           SUM0 = 0.0d0
          DO J=1,N
           SUM0 = SUM0 + W(I,J)
          ENDDO
           SUM0 = D1 + D2*SUM0
           SUM1 = SURFAC*SUM0
          DO J=1,N
           B(I,J) = B(I,J) + D2*W(I,J) - SUM1*EMU(J)*WT(J)
          ENDDO
          B(I,I) = B(I,I) + D1
          H(I) = H(I) + SUM0*ZFLUX
          DO J=1,N
           AA(I,J) = - D2*U1(I,J)
          ENDDO
           AA(I,I) = AA(I,I) + D1
           C1(I) = 0.0d0
        ENDDO
       ENDIF
!------------intermediate points:  can be even or odd, A & C diagonal
      ELSE
        DELTAU = ZTAU(ID+1) - ZTAU(ID-1)
        MSTART = M + MOD(ID+1,2)
        DO I=1,N
          A(I) = EMU(I)/DELTAU
          C1(I) = -A(I)
           SUM0 = 0.0d0
          DO IM=MSTART,MFIT,2
           SUM0 = SUM0 + POMEGA(IM,ID)*PM(I,IM)*PM0(IM)
          ENDDO
          H(I) = SUM0*FZ(ID)
          DO J=1,I
            SUM0 = 0.0d0
           DO IM=MSTART,MFIT,2
            SUM0 = SUM0 + POMEGA(IM,ID)*PM(I,IM)*PM(J,IM)
           ENDDO
            B(I,J) =  - SUM0*WT(J)
            B(J,I) =  - SUM0*WT(I)
          ENDDO
          B(I,I) = B(I,I) + 1.0d0
        ENDDO
      ENDIF
      RETURN
      END
