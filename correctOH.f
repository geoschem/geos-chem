C**********************************************************************
      SUBROUTINE CORRECTOH(BOH,COLON,COLAT,COALT)
C**********************************************************************
      IMPLICIT NONE
#     include "CMN_OH"
      INTEGER N,I,J,K
      INTEGER COLON,COLAT,COALT
      REAL*8 BOH(MLONBX,MLATBX,MVRTBX)
C**********************************************************************
C Created by Bryan Duncan.     
C**********************************************************************
c This SR corrects OH to remove the problem with high
c background aerosol (tau=0.1) in the OH parameterization.
c This file contains ratios of OH calculated assuming a uniform
c tau=0.01 over the globe and assuming tau=0.1. OH calculated by
c CMS, June, 2001. Ratio created by bnd, June, 2001.
c*************************************************************************

c      print*,'COLON,COLAT,COALT=',COLON,COLAT,COALT
c      print*,'BOH(COLON,COLAT,COALT)=',BOH(COLON,COLAT,COALT)
c      print*,'OH_SEASON=',OH_SEASON
c      print*,'OH_LAT=',OH_LAT
c      print*,'OH_PRESS=',OH_PRESS

      N=OH_SEASON
C
      DO J=1,NCMSLATS2
C
       DO K=NCMSALTS2,1,-1
C
C         print*,'a',j,k,OH_LAT,CMSLATS2(J)
         IF(OH_LAT.GE.CMSLATS2(J)) THEN         
C
C         print*,'b',j,k,OH_PRESS,CMSALTS2(K)
           IF(CMSALTS2(K).GE.OH_PRESS) THEN
C         print*,'c',j,k,n,correction(N,J,K)
            BOH(COLON,COLAT,COALT)=BOH(COLON,COLAT,COALT)*
     *           correction(N,J,K)
              GOTO 2
           ENDIF
C
           IF(K.EQ.1) THEN
C         print*,'d',j,k,n,correction(N,J,K)
            BOH(COLON,COLAT,COALT)=BOH(COLON,COLAT,COALT)*
     *           correction(N,J,K)
              GOTO 2
           ENDIF
C
         ENDIF
C
       ENDDO
C
      ENDDO
C
C Error Check.
      PRINT*,'STOPPED IN SR correctOH!'
      PRINT*,'Point lies nowhere!'
      STOP
C
 2    CONTINUE
C         print*,'e',BOH(COLON,COLAT,COALT)
      IF(BOH(COLON,COLAT,COALT).LT.0.) THEN
         print*,COLON,COLAT,COALT,BOH(COLON,COLAT,COALT)
         print*,'went 0!'
         stop
      ENDIF
C
      RETURN
      END
