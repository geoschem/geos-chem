C $Id: paramC.f,v 1.1 2003/06/30 20:26:01 bmy Exp $
C**********************************************************************
      SUBROUTINE PARAMC(BLAT,BSEASON,COCOUNT,NSEAS,NLATS,ALATS)
C**********************************************************************
      IMPLICIT NONE
      INTEGER K,M,J,BSEASON,COCOUNT,NSEAS,NLATS
      REAL*8 BLAT,ALATS(NLATS)
C
C Information used in Error Check.
C
      INTEGER NSTOP,NERROR,CSEASON1(5),CSEASON2(5),CSEASON3(6),
     *        CSEASON4(6)
      REAL*8 DLATS2(2),DLATS3(4),DLATS4(4),DLATS5(4),DLATS6(4),
     *        DLATS7(3),DLATS8(1)
      DATA CSEASON1/1,3,7,11,15/
      DATA CSEASON2/4,8,12,16,19/
      DATA CSEASON3/5,9,13,17,20,22/
      DATA CSEASON4/2,6,10,14,18,21/
      DATA DLATS2/1,2/
      DATA DLATS3/3,4,5,6/
      DATA DLATS4/7,8,9,10/
      DATA DLATS5/11,12,13,14/
      DATA DLATS6/15,16,17,18/
      DATA DLATS7/19,20,21/
      DATA DLATS8/22/
C**********************************************************************
C Created by Bryan Duncan.
C**********************************************************************
C SR PARAMC is called by SR CLEANISOP.  Based on season and latitude, this
C subroutine choses the appropriate identification number of the
C parameterization required.
C
C  ----------------------------------------------------------
C  | Latitude   |   DJF    |   MAM    |   JJA    |   SON    |
C  ----------------------------------------------------------
C  |   > 80 N   |   AVG    |   AVG    |    *     |   AVG    |
C  |  60-80 N   |   AVG    |    *     | 22 44 66 |   AVG    |
C  |  40-60 N   |   AVG    | 19 41 63 | 20 42 64 | 21 43 65 |
C  |  30-40 N   | 15 37 59 | 16 38 60 | 17 39 61 | 18 40 62 |
C  |   0-30 N   | 11 33 55 | 12 34 56 | 13 35 57 | 14 36 58 |
C  |   0-30 S   |  7 29 51 |  8 30 52 |  9 31 53 | 10 32 54 |
C  |  30-40 S   |  3 25 47 |  4 26 48 |  5 27 49 |  6 28 50 |
C  |  40-60 S   |  1 23 45 |    *     |    *     |  2 24 46 |
C  |  60-80 S   |     *    |   AVG    |   AVG    |     *    |
C  |   > 80 S   |     *    |   AVG    |   AVG    |    AVG   |
C  ----------------------------------------------------------
C
C The table above illustrates the spatial and temporal coverage of the
C parameterizations found in the "CLEANISOP" subdomain.  In the boxes labeled
C "AVG", the OH is sufficiently low so that the climatological mean OH for
C that season and latitude is used (Spivakovsky et al., 1999).  (For these
C boxes, the parameterized OH is not calculated.)  For all other boxes, the
C parameterized OH is calculated. The 3 numbers in the box correspond
C to the identification numbers of the parameterization for that box.
C The first number is for the surface layer, the second for the middle
C tropospheric layer and the third for the upper tropospheric layer.
C For instance, the parameterization the code would use to predict OH for
C April, 15N, middle troposphere has the identification number 34.
C The boxes designated by an asterisk are not covered by this domain.
C**********************************************************************
C List of Variables & Arrays
C
C COCOUNT = identification number of parameterization.
C
C**********************************************************************
C Loop over latitudes.
C
      DO K=2,NLATS
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
        IF(BLAT.LT.ALATS(1)) NSTOP=1
        IF(BLAT.LT.ALATS(2).AND.BLAT.GE.ALATS(1)) THEN
           NERROR=0
          IF(BSEASON.EQ.2) NSTOP=1
          IF(BSEASON.EQ.3) NSTOP=1
        ENDIF
C
        IF(BLAT.LT.ALATS(3).AND.BLAT.GE.ALATS(2)) THEN
          COCOUNT=2
        ENDIF
C
        IF(BLAT.LT.ALATS(7).AND.BLAT.GE.ALATS(3)) THEN
          COCOUNT=(K-1)*NSEAS-6
          IF(BLAT.GE.ALATS(6)) THEN
c              COCOUNT=COCOUNT-1
              IF(BSEASON.EQ.1) NSTOP=1
          ENDIF
        ENDIF
        IF(BLAT.GE.ALATS(7)) THEN
          IF(BSEASON.NE.3) NSTOP=1
          COCOUNT=(K-1)*NSEAS-2-1
        ENDIF
C
C Error Check.
C
        IF(NERROR.EQ.1.AND.COCOUNT.EQ.0) THEN
          PRINT*,'Latitude band not found!'
          NSTOP=1
        ENDIF
C
C NSTOP is an error check to see if a parameterization which
C       does not exist has been called.  This error indicates
C       that there is a problem in SR SKIPS.

        IF(NSTOP.EQ.1) THEN
          PRINT*,'STOPPED IN SR PARAMC'
          STOP
        ENDIF
C
C***********
C Loop over seasons.
C
        DO M=1,NSEAS
C
C***********
        IF(BLAT.LT.ALATS(2).AND.BLAT.GE.ALATS(1).AND.
     *       (M.EQ.2.OR.M.EQ.3)) GOTO 299
        IF(BLAT.LT.ALATS(7).AND.BLAT.GE.ALATS(6).AND.
     *       M.EQ.1) GOTO 299
        IF(BLAT.GE.ALATS(7).AND.M.NE.3) GOTO 299
C
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
 298     CONTINUE
C***************
      ENDDO   ! End Loop over Latitude Bands.
C***************
C
 300        CONTINUE
C
C Error Check.

      IF(BSEASON.EQ.1) THEN
        DO J=1,5
         IF(COCOUNT.EQ.CSEASON1(J)) GOTO 301
         IF(COCOUNT.EQ.CSEASON1(J)+22) GOTO 301
         IF(COCOUNT.EQ.CSEASON1(J)+44) GOTO 301
        ENDDO
        PRINT*,'Season is DJF:  Did not call a DJF param!'
      ENDIF
      IF(BSEASON.EQ.2) THEN
        DO J=1,5
         IF(COCOUNT.EQ.CSEASON2(J)) GOTO 301
         IF(COCOUNT.EQ.CSEASON2(J)+22) GOTO 301
         IF(COCOUNT.EQ.CSEASON2(J)+44) GOTO 301
        ENDDO
        PRINT*,'Season is MAM:  Did not call a MAM param!'
      ENDIF
      IF(BSEASON.EQ.3) THEN
        DO J=1,6
         IF(COCOUNT.EQ.CSEASON3(J)) GOTO 301
         IF(COCOUNT.EQ.CSEASON3(J)+22) GOTO 301
         IF(COCOUNT.EQ.CSEASON3(J)+44) GOTO 301
        ENDDO
        PRINT*,'Season is JJA:  Did not call a JJA param!'
      ENDIF
      IF(BSEASON.EQ.4) THEN
        DO J=1,6
         IF(COCOUNT.EQ.CSEASON4(J)) GOTO 301
         IF(COCOUNT.EQ.CSEASON4(J)+22) GOTO 301
         IF(COCOUNT.EQ.CSEASON4(J)+44) GOTO 301
        ENDDO
        PRINT*,'Season is SON:  Did not call a SON param!'
      ENDIF

        PRINT*,'Stopped in SR paramC.'
        PRINT*,'COCOUNT = ',COCOUNT
        STOP

 301  CONTINUE
C
      IF(BLAT.LT.ALATS(2).AND.BLAT.GE.ALATS(1)) THEN
        DO J=1,2
         IF(COCOUNT.EQ.DLATS2(J)) GOTO 302
        ENDDO
        PRINT*,'Latitude is not 40-60S: Did not call a 40-60S param!'
      ENDIF

      IF(BLAT.LT.ALATS(3).AND.BLAT.GE.ALATS(2)) THEN
        DO J=1,4
         IF(COCOUNT.EQ.DLATS3(J)) GOTO 302
        ENDDO
        PRINT*,'Latitude is not 30-40S: Did not call a 30-40S param!'
      ENDIF

      IF(BLAT.LT.ALATS(4).AND.BLAT.GE.ALATS(3)) THEN
        DO J=1,4
         IF(COCOUNT.EQ.DLATS4(J)) GOTO 302
        ENDDO
        PRINT*,'Latitude is not 0-30S: Did not call a 0-30S param!'
      ENDIF

      IF(BLAT.LT.ALATS(5).AND.BLAT.GE.ALATS(4)) THEN
        DO J=1,4
         IF(COCOUNT.EQ.DLATS5(J)) GOTO 302
        ENDDO
        PRINT*,'Latitude is not 0-30N: Did not call a 0-30N param!'
      ENDIF

      IF(BLAT.LT.ALATS(6).AND.BLAT.GE.ALATS(5)) THEN
        DO J=1,4
         IF(COCOUNT.EQ.DLATS6(J)) GOTO 302
        ENDDO
        PRINT*,'Latitude is not 30-40N: Did not call a 30-40N param!'
      ENDIF

      IF(BLAT.LT.ALATS(7).AND.BLAT.GE.ALATS(6)) THEN
        DO J=1,3
         IF(COCOUNT.EQ.DLATS7(J)) GOTO 302
        ENDDO
        PRINT*,'Latitude is not 40-60N: Did not call a 40-60N param!'
      ENDIF

      IF(BLAT.GE.ALATS(7)) THEN
        DO J=1,1
         IF(COCOUNT.EQ.DLATS8(J)) GOTO 302
        ENDDO
        PRINT*,'Latitude is not gt 60N: Did not call a 60-80N param!'
      ENDIF

        PRINT*,'Stopped in SR paramC.'
        PRINT*,'COCOUNT = ',COCOUNT
        STOP
C
 302  RETURN
      END
