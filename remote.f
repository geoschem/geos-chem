C $Id: remote.f,v 1.1 2003/06/30 20:26:10 bmy Exp $
C**********************************************************************
      SUBROUTINE REMOTE(NDFUNCS)
C**********************************************************************
      IMPLICIT NONE
#     include "CMN_OH"
      INTEGER NDFUNCS,NALL,K
C**********************************************************************
C Select proper identification number of parameterization needed for 
C   box characterized by remote conditions and no isoprene using
C   SR PARAMB. Store this value in COCOUNT. Then select the corresponding 
C   polynomials and store in NDFUNCS.
C**********************************************************************
C   Created by Bryan Duncan.
C**********************************************************************
C B) Remote
C***************************************************
C  B1) Surface
C***************************************************
C
      IF(OH_PRESS.GE.PRESSES(1)) THEN
C
C Determine the identification number of the appropriate parameterization.
C
        CALL PARAMB(OH_LAT,OH_SEASON,COCOUNT,NSEAS,NLATS,ALATS)
C
C Error Check.
C
         IF(COCOUNT.EQ.0) THEN
           PRINT*,'COCOUNT=0!'
           PRINT*,'No parameterization was chosen!'
           PRINT*,'Stopped in SR REMOTE.'
           STOP 
         ENDIF
C
C Determine the polynomials associated with selected identification
C   number.
C
         IF(COCOUNT.GT.1) THEN
          DO K=1,COCOUNT-1
           NDFUNCS=NDFUNCS+NDONE(K)
          ENDDO
         ENDIF
C
C Error Check.
C
            NALL=NPARAMB
         IF(COCOUNT.GT.NALL) THEN
           PRINT*,'COCOUNT=',COCOUNT
           PRINT*,'NALL   =',NALL
           PRINT*,'COCOUNT > NALL in SR REMOTE'
           STOP
         ENDIF
C
            GOTO 321
C
      ENDIF
C
C***************************************************
C  END B1
C***************************************************
C Sum of all NDONE functions in B).
C
            NDFUNCS=NDFUNCS+NDFUNCSB1
C
C***************************************************
C Error check.
C    
      PRINT*,'********************************'
      PRINT*,'Stopped in SR REMOTE!'
      PRINT*,'No parameterization was chosen!'
      PRINT*,'********************************'
      STOP
C
C***************************************************
 321  RETURN
      END
