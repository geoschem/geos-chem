! $Id: getifsun.f,v 1.1 2003/06/30 20:26:08 bmy Exp $
      INTEGER FUNCTION GETIFSUN(SUNCOS)

      ! References to F90 modules (bmy, 10/19/00)
      USE COMODE_MOD, ONLY : IXSAVE, IYSAVE, JLOP

      IMPLICIT NONE
#     include "CMN_SIZE"
#     include "comode.h"

      INTEGER   I,J,K,JLOOP,IJWINDOW,IX,IY
      REAL*8    SUNCOS(MAXIJ)
C
*** see if photolysis should be considered.
c  Get the right index for SUNCOS, which is calculated
c  outside of chemistry module.
C  (This works for LEMBED= .TRUE. or .FALSE.)

      K       = 0
      DO 240 J = 1, NLAT
         DO 230 I = 1, NLONG
            JLOOP = JLOP(I,J,1)
            IF (JLOOP.EQ.0) GOTO 230
            IX=IXSAVE(JLOOP)
            IY=IYSAVE(JLOOP)
            IJWINDOW         = (IY-1)*IIPAR + IX
            IF(SUNCOS(IJWINDOW).GT.0.D0) K = 1
 230     CONTINUE
 240  CONTINUE
      GETIFSUN   = 2 - K

      RETURN
      END

