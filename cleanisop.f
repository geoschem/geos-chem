C $Id: cleanisop.f,v 1.1 2003/06/30 20:26:03 bmy Exp $
C**********************************************************************
      SUBROUTINE CLEANISOP(NDFUNCS)
C**********************************************************************
      IMPLICIT NONE
#     include "CMN_OH"
      INTEGER NDFUNCS,NALL,K
C**********************************************************************
C Select proper identification number of parameterization needed for
C   box characterized by unpolluted conditions with isoprene using
C   SR PARAMC. Store this value in COCOUNT. Then select the corresponding
C   polynomials and store in NDFUNCS.
C**********************************************************************
C   Created by Bryan Duncan.
C**********************************************************************
C C) Clean
C***************************************************
C  C1) Surface
C***************************************************
C
      IF(OH_PRESS.GE.PRESSES(1)) THEN
C
C Determine the identification number of the appropriate parameterization.
C
        CALL PARAMC(OH_LAT,OH_SEASON,COCOUNT,NSEAS,NLATS,ALATS)
C
C Error Check.
C
         IF(COCOUNT.EQ.0) THEN
           PRINT*,'COCOUNT=0 in SR CLEANISOP'
           STOP
         ENDIF
C
C Determine the polynomials associated with selected identification
C   number.
C
            NALL=3*NPARAMA+NPARAMB
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
            NALL=NALL+NPARAMC
         IF(COCOUNT.GT.NALL) THEN
           PRINT*,'COCOUNT=',COCOUNT
           PRINT*,'NALL   =',NALL
           PRINT*,'COCOUNT > NALL in SR CLEANISOP'
           STOP
         ENDIF
C
            GOTO 321
C
      ENDIF
C
C***************************************************
C  END C1
C***************************************************
C Sum of all NDONE functions in C).
C
            NDFUNCS=NDFUNCS+NDFUNCSC1
C
C***************************************************
C  C2) Middle Troposphere
C***************************************************
C
      IF(OH_PRESS.LT.PRESSES(1).AND.OH_PRESS.GE.PRESSES(2)) THEN
C
C Determine the identification number of the appropriate parameterization.
C
        CALL PARAMC(OH_LAT,OH_SEASON,COCOUNT,NSEAS,NLATS,ALATS)
C
C Error Check.
C
         IF(COCOUNT.EQ.0) THEN
           PRINT*,'COCOUNT=0 in SR CLEANISOP'
           STOP
         ENDIF
C
C Determine the polynomials associated with selected identification
C   number.
C
            NALL=3*NPARAMA+NPARAMB+NPARAMC
         IF(COCOUNT.GT.1) THEN
          DO K=NALL+1,NALL+COCOUNT-1
           NDFUNCS=NDFUNCS+NDONE(K)
          ENDDO
         ENDIF

          COCOUNT=COCOUNT+NALL

C Error Check.
C
            NALL=NALL+NPARAMC
         IF(COCOUNT.GT.NALL) THEN
           PRINT*,'COCOUNT=',COCOUNT
           PRINT*,'NALL   =',NALL
           PRINT*,'COCOUNT > NALL in SR CLEANISOP'
           STOP
         ENDIF
C
            GOTO 321
C
      ENDIF
C
C***************************************************
C  END C2
C***************************************************
C
            NDFUNCS=NDFUNCS+NDFUNCSC2
C
C***************************************************
C  C3) Upper Troposphere
C***************************************************
C
      IF(OH_PRESS.LT.PRESSES(2)) THEN
C
C Determine the identification number of the appropriate parameterization.
C
        CALL PARAMC(OH_LAT,OH_SEASON,COCOUNT,NSEAS,NLATS,ALATS)
C
C Error Check.
C
         IF(COCOUNT.EQ.0) THEN
           PRINT*,'COCOUNT=0 in SR CLEANISOP'
           STOP
         ENDIF
C
C Determine the polynomials associated with selected identification
C   number.
C
            NALL=3*NPARAMA+NPARAMB+2*NPARAMC
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
            NALL=NALL+NPARAMC
         IF(COCOUNT.GT.NALL) THEN
           PRINT*,'COCOUNT=',COCOUNT
           PRINT*,'NALL   =',NALL
           PRINT*,'COCOUNT > NALL in SR CLEANISOP'
           STOP
         ENDIF
C
            GOTO 321
C
      ENDIF
C
C***************************************************
C  END C3
C***************************************************
C
            NDFUNCS=NDFUNCS+NDFUNCSC3
C
C***************************************************
C Error check.
C
      PRINT*,'********************************'
      PRINT*,'Stopped in SR CLEANISOP!'
      PRINT*,'No parameterization was chosen!'
      PRINT*,'********************************'
      STOP
C
C***************************************************
 321  RETURN
      END
