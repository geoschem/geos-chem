C $Id: modpolluted.f,v 1.1 2003/06/30 20:26:02 bmy Exp $
C**********************************************************************
      SUBROUTINE MODPOLLUTED(NDFUNCS)
C**********************************************************************
      IMPLICIT NONE
#     include "CMN_OH"
      INTEGER NDFUNCS,NALL,K
C**********************************************************************
C Select proper identification number of parameterization needed for
C   box characterized by unpolluted conditions with isoprene using
C   SR PARAMD-F. Store this value in COCOUNT. Then select the corresponding
C   polynomials and store in NDFUNCS.
C**********************************************************************
C   Created by Bryan Duncan.
C**********************************************************************
C D) Moderately Polluted
C***************************************************
C  D1) Surface
C***************************************************
C
      IF(OH_PRESS.GE.PRESSES(1)) THEN
C
C Determine the identification number of the appropriate parameterization.
C
        CALL PARAMD(OH_LAT,OH_SEASON,COCOUNT,NSEAS,NLATS,ALATS)
C
C Error Check.
C
         IF(COCOUNT.EQ.0) THEN
           PRINT*,'No parameterization was chosen!'
           PRINT*,'COCOUNT=0 in SR MODPOLLUTED'
           STOP
         ENDIF
C
C Determine the polynomials associated with selected identification
C   number.
C
            NALL=3*NPARAMA+NPARAMB+3*NPARAMC
         IF(COCOUNT.GT.1) THEN
          DO K=NALL+1,NALL+COCOUNT-1
           NDFUNCS=NDFUNCS+NDONE(K)
          ENDDO
         ENDIF
C
          COCOUNT=COCOUNT+NALL
C
C Error Check.
C
            NALL=NALL+NPARAMD
         IF(COCOUNT.GT.NALL) THEN
           PRINT*,'COCOUNT=',COCOUNT
           PRINT*,'NALL   =',NALL
           PRINT*,'COCOUNT > NALL in SR MODPOLLUTED'
           STOP
         ENDIF
C
            GOTO 321
C
      ENDIF
C
C***************************************************
C  END D1
C***************************************************
C  F1) Middle Troposphere
C***************************************************
C
      IF(OH_PRESS.LT.PRESSES(1).AND.OH_PRESS.GE.PRESSES(2)) THEN
C
C Determine the identification number of the appropriate parameterization.
C
        CALL PARAMF(OH_LAT,OH_SEASON,COCOUNT,NSEAS,NLATS,ALATS)
C
C Error Check.
C
         IF(COCOUNT.EQ.0) THEN
           PRINT*,'COCOUNT=0 in SR MODPOLLUTED'
           STOP
         ENDIF
C
C Determine the polynomials associated with selected identification
C   number.
C
            NALL=3*NPARAMA+NPARAMB+3*NPARAMC+NPARAMD+NPARAME
         IF(COCOUNT.GT.1) THEN
          DO K=NALL+1,NALL+COCOUNT-1
           NDFUNCS=NDFUNCS+NDONE(K)
          ENDDO
         ENDIF
C
          COCOUNT=COCOUNT+NALL
C
C Error Check.
C
            NALL=NALL+NPARAMF
         IF(COCOUNT.GT.NALL) THEN
           PRINT*,'COCOUNT=',COCOUNT
           PRINT*,'NALL   =',NALL
           PRINT*,'COCOUNT > NALL in SR MODPOLLUTED'
           STOP
         ENDIF
C
            GOTO 321
C
      ENDIF
C
C***************************************************
C  END F1
C***************************************************
C Error check.
C
      PRINT*,'********************************'
      PRINT*,'Stopped in SR MODPOLLUTED!'
      PRINT*,'No parameterization was chosen!'
      PRINT*,'********************************'
      STOP
C
C***************************************************
 321  RETURN
      END
