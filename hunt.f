C $Id: hunt.f,v 1.1 2003/06/30 20:26:04 bmy Exp $
      SUBROUTINE HUNT(XA,NDIM,XPT,ILO)

      IMPLICIT NONE
C**********************************************************************
C                                                                     *
C  HARVARD TROPOSPHERIC CHEMISTRY MODULE FOR 3-D APPLICATIONS         *
C  by Larry Horowitz, Jinyou Liang, Gerry Gardner, Prof. Daniel Jacob *
C                           (Release V2.0)                            *
C**********************************************************************
C
C       For monotonically increasing OR decreasing array XA, dimension NDIM,
C       returns value ILO such that XPT is between XA(I) and XA(I+1).
C       ILO serves as initial guess.
C       If XPT is out of range, I=0 or I=NDIM is returned.
C       Page 91, Numerical Recipies.
C
      INTEGER NDIM,ILO,IHI,INC,IM
      REAL*8 XPT
      REAL*8 XA(NDIM)

      LOGICAL LASCND

      
      LASCND=XA(NDIM).GT.XA(1)
      IF(ILO.LE.0.OR.ILO.GT.NDIM) THEN
         ILO=0
         IHI=NDIM+1
         GOTO 3
      ENDIF
      INC=1
      IF(XPT.GE.XA(ILO).EQV.LASCND) THEN
 1       IHI=ILO+INC
         IF(IHI.GT.NDIM) THEN
            IHI=NDIM+1
         ELSE IF(XPT.GE.XA(IHI).EQV.LASCND) THEN
            ILO=IHI
            INC=INC+INC
            GOTO 1
         ENDIF
      ELSE
         IHI=ILO
 2       ILO=IHI-INC
         IF(ILO.LT.1)THEN
            ILO=0
         ELSE IF(XPT.LT.XA(ILO).EQV.LASCND) THEN
            IHI=ILO
            INC=INC+1
            GOTO 2
         ENDIF
      ENDIF
      
 3    DO WHILE(IHI-ILO.GT.1)
         IM=(IHI+ILO)/2
         IF(XPT.GT.XA(IM).EQV.LASCND) THEN
            ILO=IM
         ELSE
            IHI=IM
         ENDIF
      ENDDO
      
      IF(ILO.EQ.0) THEN
         IF(ABS(XPT-XA(1)).LT.1.E-3) ILO=1
      ENDIF
      IF(ILO.EQ.NDIM) THEN
         IF(ABS(XPT-XA(NDIM)).LT.1.E-3) ILO=NDIM-1
      ENDIF
      
      RETURN
      END

