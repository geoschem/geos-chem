! $Id: getifsun.f,v 1.2 2010/02/26 18:19:59 bmy Exp $
      INTEGER FUNCTION GETIFSUN(SUNCOS)

      ! References to F90 modules (bmy, 10/19/00)
      USE COMODE_MOD, ONLY : IXSAVE, IYSAVE, JLOP

      USE CMN_SIZE_MOD
      USE COMODE_LOOP_MOD

      IMPLICIT NONE
#     include "define.h"

      INTEGER   I,J,K,JLOOP,IJWINDOW,IX,IY
      REAL*8    SUNCOS(MAXIJ)
!
*** see if photolysis should be considered.
c  Get the right index for SUNCOS, which is calculated
c  outside of chemistry module.

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

