C $Id: clean.f,v 1.1 2003/06/30 20:26:04 bmy Exp $
C**********************************************************************
      SUBROUTINE CLEAN(NDFUNCS)
C**********************************************************************
      IMPLICIT NONE
#     include "CMN_OH"
      INTEGER NDFUNCS,NALL,K
C**********************************************************************
C Select proper identification number of parameterization needed for 
C   box characterized by unpolluted conditions and no isoprene using
C   SR PARAMA. Store this value in COCOUNT. Then select the corresponding 
C   polynomials and store in NDFUNCS.
C**********************************************************************
C   Created by Bryan Duncan.
C**********************************************************************
C A) Clean
C***************************************************
C  A1) Surface
C***************************************************
C
      IF(OH_PRESS.GE.PRESSES(1)) THEN
C
C Determine the identification number of the appropriate parameterization.
C
        CALL PARAMA(OH_LAT,OH_SEASON,COCOUNT,NSEAS,NLATS,ALATS)
C
C Error Check.
C
         IF(COCOUNT.EQ.0) THEN
           PRINT*,'COCOUNT=0!'
           PRINT*,'No parameterization was chosen!'
           PRINT*,'Stopped in SR CLEAN.'
           STOP 
         ENDIF
C
C Determine the polynomials associated with selected identification
C   number.
C
            NALL=NPARAMB
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
            NALL=NALL+NPARAMA
         IF(COCOUNT.GT.NALL) THEN
           PRINT*,'COCOUNT=',COCOUNT
           PRINT*,'NALL   =',NALL
           PRINT*,'COCOUNT > NALL'
           PRINT*,'Stopped in SR CLEAN.'
           STOP 
         ENDIF
C
            GOTO 321
C
      ENDIF
C
C***************************************************
C  END A1
C***************************************************
C Sum of all NDONE functions in A).
C
            NDFUNCS=NDFUNCS+NDFUNCSA1
C
C***************************************************
C  A2) Middle Troposphere
C***************************************************
C
      IF(OH_PRESS.LT.PRESSES(1).AND.OH_PRESS.GE.PRESSES(2)) THEN
C
C Determine the identification number of the appropriate parameterization.
C
        CALL PARAMA(OH_LAT,OH_SEASON,COCOUNT,NSEAS,NLATS,ALATS)
C
C Error Check.
C
         IF(COCOUNT.EQ.0) THEN
           PRINT*,'COCOUNT=0!'
           PRINT*,'No parameterization was chosen!'
           PRINT*,'Stopped in SR CLEAN.'
           STOP 
         ENDIF
C
C Determine the polynomials associated with selected identification
C   number.
C
            NALL=NPARAMB+NPARAMA
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
            NALL=NALL+NPARAMA
         IF(COCOUNT.GT.NALL) THEN
           PRINT*,'COCOUNT=',COCOUNT
           PRINT*,'NALL   =',NALL
           PRINT*,'COCOUNT > NALL in SR CLEAN'
           STOP
         ENDIF
C
            GOTO 321
C
      ENDIF
C
C***************************************************
C  END A2
C***************************************************
C
            NDFUNCS=NDFUNCS+NDFUNCSA2
C
C***************************************************
C  A3) Upper Troposphere
C***************************************************
C
      IF(OH_PRESS.LT.PRESSES(2)) THEN
C
C Determine the identification number of the appropriate parameterization.
C
        CALL PARAMA(OH_LAT,OH_SEASON,COCOUNT,NSEAS,NLATS,ALATS)
C
C Error Check.
C
         IF(COCOUNT.EQ.0) THEN
           PRINT*,'COCOUNT=0!'
           PRINT*,'No parameterization was chosen!'
           PRINT*,'Stopped in SR CLEAN.'
           STOP 
         ENDIF
C
C Determine the polynomials associated with selected identification
C   number.
C
            NALL=NPARAMB+2*NPARAMA
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
           NALL=NALL+NPARAMA
         IF(COCOUNT.GT.NALL) THEN
           PRINT*,'COCOUNT=',COCOUNT
           PRINT*,'NALL   =',NALL
           PRINT*,'COCOUNT > NALL in SR CLEAN'
           STOP
         ENDIF
C
            GOTO 321
C
      ENDIF
C
C***************************************************
C  END A3
C***************************************************
C Error check.
C    
      PRINT*,'********************************'
      PRINT*,'Stopped in SR CLEAN!'
      PRINT*,'No parameterization was chosen!'
      PRINT*,'********************************'
      STOP
C
C***************************************************
 321  RETURN
      END
