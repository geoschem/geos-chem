C $Id: paramE.f,v 1.1 2003/06/30 20:26:01 bmy Exp $
C**********************************************************************
      SUBROUTINE PARAME(BLAT,BSEASON,COCOUNT,NSEAS,NLATS,ALATS)
C**********************************************************************
      IMPLICIT NONE
      INTEGER K,M,BSEASON,COCOUNT
      INTEGER NSEAS,NLATS,NSTOP
      REAL*8 BLAT,ALATS(NLATS)
C**********************************************************************
C Created by Bryan Duncan.
C**********************************************************************
C SR PARAME is called by SR HIPOLLUTEDISOP.  Based on season and latitude,
C this subroutine choses the appropriate identification number of the
C parameterization required for surface only.
C
C  ----------------------------------------------------------
C  | Latitude   |   DJF    |   MAM    |   JJA    |   SON    |
C  ----------------------------------------------------------
C  |   > 80 N   |   AVG    |   AVG    |    *     |   AVG    |
C  |  60-80 N   |   AVG    |    *     |    *     |   AVG    |
C  |  40-60 N   |   AVG    |    *     |   12     |   13     |
C  |  30-40 N   |     *    |    9     |   10     |   11     |
C  |   0-30 N   |     5    |    6     |    7     |    8     |
C  |   0-30 S   |     1    |    2     |    3     |    4     |
C  |  30-40 S   |     *    |    *     |    *     |    *     |
C  |  40-60 S   |     *    |    *     |    *     |    *     |
C  |  60-80 S   |     *    |   AVG    |   AVG    |    *     |
C  |   > 80 S   |     *    |   AVG    |   AVG    |   AVG    |
C  ----------------------------------------------------------
C
C The table above illustrates the spatial and temporal coverage of the
C parameterizations found in the "HIPOLLUTEDISOP" subdomain. In boxes labeled
C "AVG", the OH is sufficiently low so that the climatological mean OH for
C that season and latitude is used (Spivakovsky et al., 1999).  (For these
C boxes, the parameterized OH is not calculated.)  For all other boxes,
C the parameterized OH is calculated.  The number in the box corresponds
C to the identification number of the parameterization for that box.
C For instance, the parameterization the code would use to predict OH for
C 35N in April has the identification number 3.
C The boxes designated by an asterisk are not covered by this domain.
C**********************************************************************
C List of Variables & Arrays
C
C COCOUNT = identification number of parameterization.
C
C**********************************************************************
C
          NSTOP=0
C
        IF(BLAT.LE.ALATS(3)) NSTOP=1 
        IF(BLAT.GE.ALATS(7)) NSTOP=1
C
        IF(BLAT.LT.ALATS(7).AND.BLAT.GE.ALATS(6)) THEN
          IF(BSEASON.LE.2) NSTOP=1
          IF(BSEASON.EQ.3) COCOUNT=12
          IF(BSEASON.EQ.4) COCOUNT=13
        ENDIF
        IF(BLAT.LT.ALATS(6).AND.BLAT.GE.ALATS(5)) THEN
          IF(BSEASON.EQ.1) NSTOP=1
          IF(BSEASON.EQ.2) COCOUNT=9
          IF(BSEASON.EQ.3) COCOUNT=10
          IF(BSEASON.EQ.4) COCOUNT=11
        ENDIF
        IF(BLAT.LT.ALATS(5).AND.BLAT.GE.ALATS(4)) THEN
          IF(BSEASON.EQ.1) COCOUNT=5
          IF(BSEASON.EQ.2) COCOUNT=6
          IF(BSEASON.EQ.3) COCOUNT=7
          IF(BSEASON.EQ.4) COCOUNT=8
        ENDIF
        IF(BLAT.LT.ALATS(4).AND.BLAT.GE.ALATS(3)) THEN
          IF(BSEASON.EQ.1) COCOUNT=1
          IF(BSEASON.EQ.2) COCOUNT=2
          IF(BSEASON.EQ.3) COCOUNT=3
          IF(BSEASON.EQ.4) COCOUNT=4
        ENDIF
C
C**********************************************************************
C Error Check.
C
C NSTOP is an error check to see if a parameterization which
C       does not exist has been called.  This error indicates
C       that there is a problem in SR SKIPS.
C
        IF(NSTOP.EQ.1) THEN
          PRINT*,'STOPPED IN SR PARAME'
          STOP
        ENDIF
C
      IF(BSEASON.EQ.1)THEN
         IF(COCOUNT.NE.1.OR.COCOUNT.NE.5) GOTO 301
        PRINT*,'Season is DJF:  Did not call a DJF param!'
      ENDIF
      IF(BSEASON.EQ.2) THEN
         IF(COCOUNT.EQ.2.OR.COCOUNT.EQ.6.OR.COCOUNT.EQ.9)
     *     GOTO 301
        PRINT*,'Season is MAM:  Did not call a MAM param!'
      ENDIF
C
      IF(BSEASON.EQ.3) THEN
         IF(COCOUNT.EQ.3.OR.COCOUNT.EQ.7.OR.COCOUNT.EQ.10.
     *      OR.COCOUNT.EQ.12) GOTO 301
        PRINT*,'Season is JJA:  Did not call a JJA param!'
      ENDIF
C
      IF(BSEASON.EQ.4) THEN
         IF(COCOUNT.EQ.4.OR.COCOUNT.EQ.8.OR.COCOUNT.EQ.11.
     *      OR.COCOUNT.EQ.13) GOTO 301
        PRINT*,'Season is SON:  Did not call a SON param!'
      ENDIF
C
        PRINT*,'Stopped in SR paramE.'
        PRINT*,'COCOUNT = ',COCOUNT
        STOP
C
 301  CONTINUE
C
      IF(BLAT.LT.ALATS(7).AND.BLAT.GE.ALATS(6)) THEN
       IF(COCOUNT.GE.12) GOTO 302
      ENDIF
      IF(BLAT.LT.ALATS(6).AND.BLAT.GE.ALATS(5)) THEN
       IF(COCOUNT.GE.9.AND.COCOUNT.LE.11) GOTO 302
      ENDIF
      IF(BLAT.LT.ALATS(5).AND.BLAT.GE.ALATS(4)) THEN
       IF(COCOUNT.GE.5.AND.COCOUNT.LE.8) GOTO 302
      ENDIF
      IF(BLAT.LT.ALATS(4).AND.BLAT.GE.ALATS(3)) THEN
       IF(COCOUNT.LE.4) GOTO 302
      ENDIF
C
        PRINT*,'Stopped in SR paramE.'
        PRINT*,'COCOUNT = ',COCOUNT
        STOP
C
C**********************************************************************
C
 302  RETURN
      END
