C $Id: o3comp.f,v 1.1 2003/06/30 20:26:09 bmy Exp $
      SUBROUTINE O3COMP(I,J,L)

      IMPLICIT NONE
C New version of O3COMP: Called from DRYDEP; finds fractions
C of O3 and NO2 for rural box and all 4 plumes boxes.

#     include "CMN_SIZE"
#     include "CMN"
#     include "CMN_O3"
#     include "CMN_DEP"
      
!**** Remove unnecessary variables (bmy, 3/16/99)
      INTEGER  I,J,L,IPL

      IF (.NOT.LTROP) THEN
         FO3(5) = 1.0
         FNO2(5)= 1.0
         FNO(5) = 1.0
         GOTO 110
      END IF

      FO3(5) = FRACO3(I,J,L)
      FNO(5) = FRACNO(I,J,L)

C  We're not using FNO2() for now
      DO IPL=1, 5
         FNO2(IPL) = 0.0
      ENDDO

 110  RETURN
      END
