C $Id: FLINT.f,v 1.1 2003/06/30 20:26:05 bmy Exp $
      REAL*8 FUNCTION FLINT (TINT,T1,T2,T3,F1,F2,F3)
C-----------------------------------------------------------------------
c  Three-point linear interpolation function
C-----------------------------------------------------------------------
      real*8 TINT,T1,T2,T3,F1,F2,F3
      IF (TINT .LE. T2)  THEN
        IF (TINT .LE. T1)  THEN
          FLINT  = F1
        ELSE
          FLINT = F1 + (F2 - F1)*(TINT -T1)/(T2 -T1)
        ENDIF
      ELSE
        IF (TINT .GE. T3)  THEN
          FLINT  = F3
        ELSE
          FLINT = F2 + (F3 - F2)*(TINT -T2)/(T3 -T2)
        ENDIF
      ENDIF
      return
      end
