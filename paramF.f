C $Id: paramF.f,v 1.1 2003/06/30 20:26:01 bmy Exp $
C**********************************************************************
      SUBROUTINE PARAMF(BLAT,BSEASON,COCOUNT,NSEAS,NLATS,ALATS)
C**********************************************************************
      IMPLICIT NONE
      INTEGER K,M,BSEASON,COCOUNT
      INTEGER NSEAS,NLATS,NSTOP
      REAL*8 BLAT,ALATS(NLATS)
C**********************************************************************
C Created by Bryan Duncan.
C**********************************************************************
C SR PARAMF is called by SR MODPOLLUTEDISOP.  Based on season and latitude,
C this subroutine choses the appropriate identification number of the
C parameterization required for middle troposphere only.
C
C  ----------------------------------------------------------
C  | Latitude   |   DJF    |   MAM    |   JJA    |   SON    |
C  ----------------------------------------------------------
C  |   > 80 N   |   AVG    |   AVG    |    *     |   AVG    |
C  |  60-80 N   |   AVG    |    *     |    *     |   AVG    |
C  |  40-60 N   |   AVG    |    5     |    6     |    7     |
C  |  30-40 N   |    1     |    2     |    3     |    4     |
C  |   0-30 N   |    *     |    *     |    *     |    *     |
C  |   0-30 S   |    *     |    *     |    *     |    *     |
C  |  30-40 S   |    *     |    *     |    *     |    *     |
C  |  40-60 S   |    *     |    *     |    *     |    *     |
C  |  60-80 S   |    *     |   AVG    |   AVG    |    *     |
C  |   > 80 S   |    *     |   AVG    |   AVG    |   AVG    |
C  ----------------------------------------------------------
C
C The table above illustrates the spatial and temporal coverage of the
C parameterizations found in the "MODPOLLUTEDISOP" subdomain. In boxes labeled
C "AVG", the OH is sufficiently low so that the climatological mean OH for
C that season and latitude is used (Spivakovsky et al., 1999).  (For these
C boxes, the parameterized OH is not calculated.)  For all other boxes,
C the parameterized OH is calculated.  The number in the box corresponds
C to the identification number of the parameterization for that box.
C For instance, the parameterization the code would use to predict OH for
C 35N in April has the identification number 2.
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
        IF(BLAT.LE.ALATS(5)) NSTOP=1 
        IF(BLAT.GE.ALATS(7)) NSTOP=1
C
        IF(BLAT.LE.ALATS(6)) THEN
          IF(BSEASON.EQ.1) COCOUNT=1
          IF(BSEASON.EQ.2) COCOUNT=2
          IF(BSEASON.EQ.3) COCOUNT=3
          IF(BSEASON.EQ.4) COCOUNT=4
        ELSE
          IF(BSEASON.EQ.1) NSTOP=1
          IF(BSEASON.EQ.2) COCOUNT=5
          IF(BSEASON.EQ.3) COCOUNT=6
          IF(BSEASON.EQ.4) COCOUNT=7
        ENDIF
C
C**********************************************************************
C Error Check.
C
C NSTOP is an error check to see if a parameterization which
C       does not exist has been called.  This error indicates
C       that there is a problem in SR SKIPS.

        IF(NSTOP.EQ.1) THEN
          PRINT*,'STOPPED IN SR PARAMF'
          STOP
        ENDIF
C
      IF(BSEASON.EQ.1) THEN
         IF(COCOUNT.EQ.1) GOTO 301
        PRINT*,'Season is DJF:  Did not call a DJF param!'
      ENDIF
C
      IF(BSEASON.EQ.2) THEN
         IF(COCOUNT.EQ.2.OR.COCOUNT.EQ.5) GOTO 301
        PRINT*,'Season is MAM:  Did not call a MAM param!'
      ENDIF
C
      IF(BSEASON.EQ.3) THEN
         IF(COCOUNT.EQ.3.OR.COCOUNT.EQ.6) GOTO 301
        PRINT*,'Season is JJA:  Did not call a JJA param!'
      ENDIF
C
      IF(BSEASON.EQ.4) THEN
         IF(COCOUNT.EQ.4.OR.COCOUNT.EQ.7) GOTO 301
        PRINT*,'Season is SON:  Did not call a SON param!'
      ENDIF
C
        PRINT*,'Stopped in SR paramF - 1.'
        PRINT*,'COCOUNT = ',COCOUNT
        STOP
C
 301  CONTINUE
C
      IF(BLAT.LE.ALATS(6)) THEN
       IF(COCOUNT.LE.4) GOTO 302
      ELSE
       IF(COCOUNT.GE.5) GOTO 302
      ENDIF
C
        PRINT*,'Stopped in SR paramF - 2.'
        PRINT*,'COCOUNT = ',COCOUNT
        STOP
C
C**********************************************************************
C
 302  RETURN
      END
