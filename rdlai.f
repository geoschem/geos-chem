! $Id: rdlai.f,v 1.1 2003/06/30 20:26:08 bmy Exp $
      SUBROUTINE RDLAI( JDAY, MONTH )

C**********************************************************************
C                                                                     *
C  HARVARD ATMOSPHERIC CHEMISTRY MODELING GROUP                       *
C  MODULE FOR SOIL NOx EMISSIONS                                      *
C  by Yuhang Wang, Gerry Gardner and Prof. Daniel Jacob               *
C  (Release V2.1)                                                     *
C                                                                     *
C  Contact person: Bob Yantosca (bmy@io.harvard.edu)                  *
C                                                                     *
C**********************************************************************
C Be sure to force double precision with the DBLE function            *
C and the "D" exponent, wherever necessary. (bmy, 10/6/99)            *
C**********************************************************************
C Replace IMX with IGLOB and JMX with JGLOB (bmy, 6/25/02)            *
C**********************************************************************

      ! References to F90 modules (bmy, 2/11/03)
      USE GRID_MOD, ONLY : GET_XOFFSET, GET_YOFFSET

      IMPLICIT NONE

C**********************************************************************
C update daily the LAIs (Leaf Area Index)                             *
C**********************************************************************
C IREG      = Number of landtypes in grid square (I,J)                *
C XLAI      = Leaf Area Index of land type element K  (I,J,K)         *
C             current month                                           *
C XLAI2     = Leaf Area Index of land type element K  (I,J,K)         *
C             following month                                         *
C XYLAI     = Leaf Area Index of land type element K  (IJLOOP,K)      *
C**********************************************************************

#     include "CMN_SIZE"
#     include "CMN_DEP"
#     include "CMN_VEL"

      INTEGER STARTDAY(13),ISAVE
      DATA STARTDAY /15,45,74,105,135,166,196,227,258,288,319,349,380/
      DATA ISAVE /0/
      SAVE ISAVE

      INTEGER IMUL
      INTEGER I,J,K,IJLOOP,MM,ITD
      INTEGER JDAY,MONTH,IREF,JREF

      ! Need to add I0, J0 as local variables (bmy, 2/11/03)
      INTEGER I0, J0
      
      ! Get nested-grid offsets (bmy, 2/11/03)
      I0 = GET_XOFFSET()
      J0 = GET_YOFFSET()

      IF (ISAVE.EQ.0) THEN
         ISAVE=1
         CALL FINDMON(JDAY,MONTH,MM,STARTDAY)
         IF (JDAY.LT.STARTDAY(1)) THEN
            IMUL=365-STARTDAY(12)+JDAY
            ITD = 31
         ELSE
            IMUL=JDAY-STARTDAY(MM)
            ITD = STARTDAY(MM+1) - STARTDAY(MM)
         END IF
         CALL READLAI(MM)
         DO J=1,JGLOB
         DO I=1,IGLOB
         DO K=1,IREG(I,J)
            XLAI2(I,J,K) = (XLAI2(I,J,K)-XLAI(I,J,K))/(DBLE(ITD))
            XLAI(I,J,K)=XLAI(I,J,K)+ XLAI2(I,J,K) * DBLE(IMUL)
         END DO
         END DO
         END DO
      ELSE
         CALL FINDMON(JDAY,MONTH,MM,STARTDAY)
         IF (JDAY.EQ.STARTDAY(MM)) THEN
            ITD = STARTDAY(MM+1) - STARTDAY(MM)
            CALL READLAI(MM)
            DO J=1,JGLOB
            DO I=1,IGLOB
            DO K=1,IREG(I,J)
               XLAI2(I,J,K) = (XLAI2(I,J,K)-XLAI(I,J,K))/(DBLE(ITD))
            END DO
            END DO
            END DO
         ELSE
            DO J=1,JGLOB
            DO I=1,IGLOB
            DO K=1,IREG(I,J)
               XLAI(I,J,K)=XLAI(I,J,K)+ XLAI2(I,J,K)
            END DO
            END DO
            END DO
         END IF
      END IF

      IJLOOP = 0
      DO J = 1, JJPAR
         JREF = J + J0
         DO I = 1, IIPAR
            IJLOOP = IJLOOP + 1
            DO K=1,IJREG(IJLOOP)
               IREF = I + I0
               XYLAI(IJLOOP,K)=XLAI(IREF,JREF,K)
            END DO
         END DO
      END DO

      ! Return to calling program
      END SUBROUTINE RDLAI
