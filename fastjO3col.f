C $Id: fastjO3col.f,v 1.1 2003/06/30 20:26:08 bmy Exp $
      SUBROUTINE fastjO3col(o3up,J,I)
C-----------------------------------------------------------------------
C SR to integrate O3 profile from fastj to calculate O3 column above
C each box
C-----------------------------------------------------------------------
c     DO3    Ozone number density at each pressure level (cm-3)
C     o3up   Ozone column above each box (DU)
C-----------------------------------------------------------------------
      IMPLICIT NONE
#     include "cmn_fj.h"
#     include "jv_cmn.h"

      INTEGER I,J,L,LL
      REAL*8 o3up(IIPAR,JJPAR,NB-1)
C
C o3up = O3 column (DU) above.
C
C Integrate overhead O3 column for each box.
C
      DO L=1,NB-1
C
         DO LL=1,NB
          
            IF(LL.GT.L) THEN
C
               o3up(I,J,L)=o3up(I,J,L)+DO3(LL)
C     
            ENDIF
C     
         ENDDO
C
         o3up(I,J,L)=o3up(I,J,L)/2.69E16
C
      ENDDO
C
      RETURN
      END
