C $Id: polint.f,v 1.1 2003/06/30 20:26:06 bmy Exp $
C*************************************************************
      SUBROUTINE POLINT(XA,YA,N,X,Y,DY)
C*************************************************************
      IMPLICIT NONE
C****************************************
C Called by SR getrefl to interpolate reflectivities
C   from Srefl.table.
C Written by cms (?/?) and adapted by bnd (4/98).
C****************************************
C
      INTEGER N,NMAX,I,M,Ns
      PARAMETER (NMAX=5)
      REAL*8 DY,X,Y,XA(N),YA(N)
      REAL*8 DEN,DIF,DIFT,HO,HP,W,C(NMAX),D(NMAX)

      NS=1
      DIF=ABS(X-XA(1))
      DO 11 I=1,N
        DIFT=ABS(X-XA(I))
        IF (DIFT.LT.DIF) THEN
          NS=I
          DIF=DIFT
        ENDIF
        C(I)=YA(I)
        D(I)=YA(I)
11    CONTINUE
      Y=YA(NS)
      NS=NS-1
      DO 13 M=1,N-1
        DO 12 I=1,N-M
          HO=XA(I)-X
          HP=XA(I+M)-X
          W=C(I+1)-D(I)
          DEN=HO-HP
          IF(DEN.EQ.0.)PAUSE 'FAILURE IN POLINT'
          DEN=W/DEN
          D(I)=HP*DEN
          C(I)=HO*DEN
12      CONTINUE
        IF (2*NS.LT.N-M)THEN
          DY=C(NS+1)
        ELSE
          DY=D(NS)
          NS=NS-1
        ENDIF
        Y=Y+DY
13    CONTINUE
C
      RETURN
      END

