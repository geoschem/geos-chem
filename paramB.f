C $Id: paramB.f,v 1.1 2003/06/30 20:26:01 bmy Exp $
C**********************************************************************
      SUBROUTINE PARAMB(BLAT,BSEASON,COCOUNT,NSEAS,NLATS,ALATS)
C**********************************************************************
      IMPLICIT NONE
      INTEGER K,M,J,BSEASON,COCOUNT,NSEAS,NLATS
      REAL*8 BLAT,ALATS(NLATS)
C
C Information used in Error Check.
C
      INTEGER NERROR,CSEASON1(18),CSEASON2(21),CSEASON3(21),
     *        CSEASON4(21),NSTOP
      REAL*8 DLATS1(6),DLATS2(12),DLATS3(12),
     * DLATS4(12),DLATS5(12),DLATS6(12),DLATS7(9),DLATS8(6)
      DATA CSEASON1/1,3,7,11,15,19,28,30,34,38,42,46,55,57,61,65,
     *              69,73/
      DATA CSEASON2/4,8,12,16,20,23,26,31,35,39,43,47,50,53,58,62,
     *              66,70,74,77,80/
      DATA CSEASON3/5,9,13,17,21,24,27,32,36,40,44,48,51,54,
     *              59,63,67,71,75,78,81/
      DATA CSEASON4/2,6,10,14,18,22,25,29,33,37,41,45,49,52,56,
     *              60,64,48,72,76,79/
      DATA DLATS1/1,2,28,29,55,56/
      DATA DLATS2/3,30,57,4,31,58,5,32,50,6,33,60/
      DATA DLATS3/7,34,61,8,35,62,9,36,63,10,37,64/
      DATA DLATS4/11,38,65,12,39,66,13,40,67,14,41,68/
      DATA DLATS5/15,42,69,16,43,70,17,44,71,18,45,72/
      DATA DLATS6/19,46,73,20,47,74,21,48,75,22,49,76/
      DATA DLATS7/23,50,77,24,51,78,25,52,79/
      DATA DLATS8/27,54,81,26,53,80/
C**********************************************************************
C Created by Bryan Duncan.
C**********************************************************************
C SR PARAMB is called by SR REMOTE.  Based on season and latitude, this
C subroutine choses the appropriate identification number of the 
C parameterization required.  
C
C  ----------------------------------------------------------
C  | Latitude   |   DJF    |   MAM    |   JJA    |   SON    | 
C  ----------------------------------------------------------
C  |   > 80 N   |   AVG    |   AVG    |   AVG    |   AVG    |
C  |  60-80 N   |   AVG    | 26 43 80 | 27 54 81 |   AVG    |
C  |  40-60 N   |   AVG    | 23 50 77 | 24 51 78 | 25 52 79 |  
C  |  30-40 N   | 19 46 73 | 20 47 74 | 21 48 75 | 22 49 76 |
C  |   0-30 N   | 15 42 69 | 16 43 70 | 17 44 71 | 18 45 72 |  
C  |   0-30 S   | 11 38 65 | 12 39 66 | 13 40 67 | 14 41 68 |  
C  |  30-40 S   |  7 34 61 |  8 35 62 |  9 36 63 | 10 37 64 |  
C  |  40-60 S   |  3 30 57 |  4 31 58 |  5 32 59 |  6 33 60 |  
C  |  60-80 S   |  1 28 55 |   AVG    |   AVG    |  2 29 56 |  
C  |   > 80 S   |   AVG    |   AVG    |   AVG    |    AVG   |  
C  ----------------------------------------------------------
C
C The table above illustrates the spatial and temporal coverage of the
C parameterizations found in the "REMOTE" subdomain.  In the boxes labeled 
C "AVG", the OH is sufficiently low so that the climatological mean OH for
C that season and latitude is used (Spivakovsky et al., 1999).  (For these 
C boxes, the parameterized OH is not calculated.)  For all other boxes,
C the parameterized OH is calculated.  The three numbers in the box correspond
C to the identification numbers of the parameterization for that box.
C The first number is for the surface layer, the second for the middle 
C tropospheric layer and the third for the upper tropospheric layer.
C For instance, the parameterization the code would use to predict OH for
C April, 15N, middle troposphere has the identification number 43.
C**********************************************************************
C List of Variables & Arrays
C
C COCOUNT = identification number of parameterization. 
C
C**********************************************************************
C Loop over latitude bands and then seasons to find appropriate
C parameterization.
C
      DO K=1,NLATS
C
C***************
C
        IF(BLAT.GE.ALATS(K).AND.K.NE.NLATS) THEN
          GOTO 298
        ENDIF
C
           NSTOP=0
           NERROR=1
C
        IF(BLAT.LT.ALATS(1)) THEN
          NERROR=0
          IF(BSEASON.EQ.2) NSTOP=1
          IF(BSEASON.EQ.3) NSTOP=1
        ENDIF 
        IF(BLAT.LT.ALATS(7).AND.BLAT.GE.ALATS(1)) THEN
         COCOUNT=(K-1)*NSEAS-2
        ENDIF 
        IF(BLAT.GE.ALATS(7)) THEN
          IF(BSEASON.EQ.1) NSTOP=1
          IF(BSEASON.EQ.4) NSTOP=1
         COCOUNT=(K)*NSEAS-2-1
        ENDIF 
C
C**********************************************************************
C
C Error Check.
C
        IF(NERROR.EQ.1.AND.COCOUNT.EQ.0) THEN
          NSTOP=1
        ENDIF
C
C NSTOP is an error check to see if a parameterization which
C       does not exist has been called.  This error indicates
C       that there is a problem in SR SKIPS.

        IF(NSTOP.EQ.1) THEN
          PRINT*,'***********************************'
          PRINT*,'Stopped in SR PARAMB/SR REMOTE!'
          PRINT*,'A parameterization which does not exist has been'
          PRINT*,'called.'
          PRINT*,'Season=',BSEASON
          PRINT*,'Latitude=',BLAT
          PRINT*,'***********************************'
          STOP
        ENDIF
C
C End Error Check.
C
C**********************************************************************
C Loop over seasons.
C
        DO M=1,NSEAS
C
C***********
C
        IF(BLAT.LT.ALATS(7).AND.BLAT.GE.ALATS(6).AND.
     *       M.EQ.1) GOTO 299 
        IF(BLAT.LT.ALATS(1).AND.(M.EQ.2.OR.M.EQ.3)) GOTO 299 
        IF(BLAT.GE.ALATS(7).AND.M.EQ.1) GOTO 299
        IF(BLAT.GE.ALATS(7).AND.M.EQ.4) GOTO 299

          IF(BSEASON.GT.M) THEN
            COCOUNT=COCOUNT+1
            GOTO 299
          ENDIF
C
                 COCOUNT=COCOUNT+1
C
                 GOTO 300
C
 299        CONTINUE
C***********
        ENDDO ! End Loop over Seasons.
C***********
 298        CONTINUE
C***************
      ENDDO   ! End Loop over Latitude Bands.
C***************
C
 300  CONTINUE
C
C**********************************************************************
C
C Error Check.
C
      IF(BSEASON.EQ.1) THEN
        DO J=1,18
         IF(COCOUNT.EQ.CSEASON1(J)) GOTO 301
        ENDDO
        PRINT*,'Season is DJF:  Did not call a DJF parameterization!'
      ENDIF
      IF(BSEASON.EQ.2) THEN
        DO J=1,21
         IF(COCOUNT.EQ.CSEASON2(J)) GOTO 301
        ENDDO
        PRINT*,'Season is MAM:  Did not call a MAM parameterization!'
      ENDIF
      IF(BSEASON.EQ.3) THEN
        DO J=1,21
         IF(COCOUNT.EQ.CSEASON3(J)) GOTO 301
        ENDDO
        PRINT*,'Season is JJA:  Did not call a JJA parameterization!'
      ENDIF
      IF(BSEASON.EQ.4) THEN
        DO J=1,21
         IF(COCOUNT.EQ.CSEASON4(J)) GOTO 301
        ENDDO
        PRINT*,'Season is SON:  Did not call a SON parameterization!'
      ENDIF

        PRINT*,'Stopped in SR PARAMB.'
        PRINT*,'COCOUNT = ',COCOUNT
        STOP
C
 301  CONTINUE
C
      IF(BLAT.LT.ALATS(1)) THEN
        DO J=1,9
         IF(COCOUNT.EQ.DLATS1(J)) GOTO 302 
        ENDDO
      PRINT*,'Latitude is not lt 60S: Did not call 60-80S param!'
      ENDIF

      IF(BLAT.LT.ALATS(2).AND.BLAT.GE.ALATS(1)) THEN
        DO J=1,12
         IF(COCOUNT.EQ.DLATS2(J)) GOTO 302 
        ENDDO
      PRINT*,'Latitude is not 40-60S: Did not call 40-60S param!'
      ENDIF

      IF(BLAT.LT.ALATS(3).AND.BLAT.GE.ALATS(2)) THEN
        DO J=1,12
         IF(COCOUNT.EQ.DLATS3(J)) GOTO 302 
        ENDDO
      PRINT*,'Latitude is not 30-40S: Did not call 30-40S param!'
      ENDIF

      IF(BLAT.LT.ALATS(4).AND.BLAT.GE.ALATS(3)) THEN
        DO J=1,12
         IF(COCOUNT.EQ.DLATS4(J)) GOTO 302 
        ENDDO
      PRINT*,'Latitude is not 0-30S: Did not call 0-30S param!'
      ENDIF

      IF(BLAT.LT.ALATS(5).AND.BLAT.GE.ALATS(4)) THEN
        DO J=1,12
         IF(COCOUNT.EQ.DLATS5(J)) GOTO 302 
        ENDDO
      PRINT*,'Latitude is not 0-30N: Did not call 0-30N param!'
      ENDIF

      IF(BLAT.LT.ALATS(6).AND.BLAT.GE.ALATS(5)) THEN
        DO J=1,12
         IF(COCOUNT.EQ.DLATS6(J)) GOTO 302 
        ENDDO
      PRINT*,'Latitude is not 30-40N: Did not call 30-40N param!'
      ENDIF

      IF(BLAT.LT.ALATS(7).AND.BLAT.GE.ALATS(6)) THEN
        DO J=1,9
         IF(COCOUNT.EQ.DLATS7(J)) GOTO 302 
        ENDDO
      PRINT*,'Latitude is not 40-60N: Did not call 40-60N param!'
      ENDIF

      IF(BLAT.GE.ALATS(7)) THEN
        DO J=1,6
         IF(COCOUNT.EQ.DLATS8(J)) GOTO 302 
        ENDDO
      PRINT*,'Latitude is not gt 60N: Did not call 60-80N param!'
      ENDIF

        PRINT*,'Stopped in SR PARAMB.'
        PRINT*,'COCOUNT = ',COCOUNT
        STOP
C
C End Error Check.
C
C**********************************************************************
C
 302  RETURN
      END
