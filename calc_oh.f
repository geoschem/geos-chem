C $Id: calc_oh.f,v 1.1 2003/06/30 20:26:06 bmy Exp $
C **********************************************************************
      SUBROUTINE CALC_OH(NDFUNCS,OOB)
C **********************************************************************
      IMPLICIT NONE
#     include "CMN_OH"
      INTEGER MM,M,KDONE,ID,JKEEP,NOROWS,KKRW,KROW,JVAR,
     *       K,NDFUNCS,NOTDONE,OOB,OOA
      REAL*8 FNCVAL,SYMBOL,BMINSA,MIDWAY,C(MXCOL),X(MXVAR),ERRORMINMAX,
     *       DIFFERR
      CHARACTER*18 ERRORCHECK(5),VARNAMES(19)
      DATA ERRORCHECK/'DECLINATION ANGLE','NOx','PRESSURE',
     *                'LATITUDE','O3'/
      DATA VARNAMES/'DECLINATION ANGLE','O3 COLUMN','NOx','OZONE',
     *   'CO','CH4','H2O','ACETONE','PRESSURE','PROPENE','ETHANE',
     *   'PROPANE','ALK4','ALBEDO','LATITUDE','TEMPERATURE',
     *   'CLOUD ABOVE','CLOUD BELOW','ISOPRENE'/
C **********************************************************************
C Created by Bryan Duncan.
C **********************************************************************
C Calculate the parameterized OH.
C **********************************************************************
C
c           print*,'*************'
c           print*,'inside SR calc_oh'
c           print*,'NPARAM=',NPARAM
c           print*,'COCOUNT=',COCOUNT
c           print*,'NDONE=',NDONE(COCOUNT)
c           print*,'*************'

       IF(COCOUNT.GT.NPARAM) THEN
            PRINT*,'*************'
            PRINT*,'COCOUNT=',COCOUNT
            PRINT*,'NPARAM=',NPARAM
            PRINT*,'SR CALC_OH: COCOUNT.GT.NPARAM!'
            PRINT*,'*************'
            STOP
       ENDIF
C
C NOTDONE = error check.
C         = 0 then done function not found.
C         = 1 then done function found.
C
       NOTDONE=0

       IF(COCOUNT.LT.1) THEN
           print*,'*************'
           print*,'STOPPED IN SR CALC_OH!'
           print*,'The box lies nowhere!' 
           print*,'NPARAM=',NPARAM
           print*,'COCOUNT=',COCOUNT
           print*,'NDONE=',NDONE(COCOUNT)
           print*,'*************'
           STOP
       ENDIF
C
C **********************************************************************
C 1) Independent variables are rescaled to fit range [-0.9, 0.9].  The 
C    parameterization was constructed with the independent variables rescaled.
C **********************************************************************
C Loop through MOVAR independent variables.
C
       DO M=1,MOVAR(COCOUNT)
C        print*,RANGEM(COCOUNT,M,1),INDVAR(M),RANGEM(COCOUNT,M,2)
C
        BMINSA=RANGEM(COCOUNT,M,2)-RANGEM(COCOUNT,M,1)
C
C If point falls out of prescribed ranges, reset to just in range.
C
        IF(INDVAR(M).LE.RANGEM(COCOUNT,M,1)) THEN
C
         IF(INDVAR(M).LT.RANGEM(COCOUNT,M,1)) THEN 
C
C************************
C Error Check.
C
C The following four independent variables should not be out of bounds
C  of the prescribed ranges.  An error here is a good indication that
C  the wrong parameterization has been chosen by the code.
C   M=1 : Declination Angle (i.e., Season)
C   M=3 : NOx
C   M=9 : Pressure (i.e., Altitude)
C   M=15: Latitude
C  
         IF(M.EQ.1.OR.M.EQ.15.OR.M.EQ.3.OR.M.EQ.9)THEN
           PRINT*,'*******************'
           PRINT*,'RANGE OUT OF BOUNDS FOR:'
           IF(M.EQ.1)  PRINT*,ERRORCHECK(1)
           IF(M.EQ.3)  PRINT*,ERRORCHECK(2)
           IF(M.EQ.9)  PRINT*,ERRORCHECK(3)
           IF(M.EQ.15) PRINT*,ERRORCHECK(4)
           IF(M.EQ.4)  PRINT*,ERRORCHECK(5)
           PRINT*,''
           PRINT*,'Minimum Range = ',RANGEM(COCOUNT,M,1) 
           PRINT*,'Maximum Range = ',RANGEM(COCOUNT,M,2) 
           PRINT*,'Variable Value= ',INDVAR(M)
           PRINT*,''
           PRINT*,'This variable should not be out of bounds'
           PRINT*,'of the prescribed ranges.  An error here is a good'
           PRINT*,'indication that the wrong parameterization has been'
           PRINT*,'chosen by the code.'
           PRINT*,''
           PRINT*,'STOPPED IN SR CALC_OH!'
           PRINT*,'*******************'
           STOP
         ENDIF
C
C End Error Check.
C************************
C
         ENDIF
C
             DIFFERR=RANGEM(COCOUNT,M,1)-INDVAR(M)
C
              OOA=M 
           IF(OOB.EQ.3.AND.M.EQ.18) OOA=19
C
             ERRORMINMAX=INDVAR(M)
             INDVAR(M)=RANGEM(COCOUNT,M,1)+1D-2*BMINSA
C
         IF(ERRORON.EQ.0.AND.DIFFERR.GT.1.D-8) THEN
           PRINT*,'*******************'
           PRINT*,'WARNING! RANGE OUT OF BOUNDS FOR:  ',VARNAMES(OOA)
           PRINT*,'Variable reset to value inside prescribed range.' 
           PRINT*,'Minimum Range        = ',RANGEM(COCOUNT,M,1) 
           PRINT*,'Maximum Range        = ',RANGEM(COCOUNT,M,2) 
           PRINT*,'Variable Value - NEW = ',INDVAR(M)
           PRINT*,'Variable Value - OLD = ',ERRORMINMAX
           PRINT*,'Latitude = ', INDVAR(15)
           PRINT*,'Longitude= ', INDVAR(20)
           PRINT*,'Pressure = ', INDVAR(9)
           PRINT*,'Dec. >   = ', INDVAR(1)
           PRINT*,'NOx      = ', INDVAR(3)

c      if(INDVAR(20).GT.-135.and.INDVAR(20).lt.-90.and.INDVAR(15).
c     *  lt.60.and.INDVAR(15).gt.40.) print*,'bryan'

           PRINT*,'*******************'
         ENDIF
C
        ENDIF
C
        IF(INDVAR(M).GE.RANGEM(COCOUNT,M,2)) THEN
C
C************************
C
         IF(INDVAR(M).GT.RANGEM(COCOUNT,M,2)) THEN
C
C************************
C Error Check.
C
         IF(M.EQ.1.OR.M.EQ.15.OR.M.EQ.3.OR.M.EQ.9)THEN
           PRINT*,'*******************'
           PRINT*,'RANGE OUT OF BOUNDS FOR:'
           IF(M.EQ.1)  PRINT*,ERRORCHECK(1)
           IF(M.EQ.3)  PRINT*,ERRORCHECK(2)
           IF(M.EQ.9)  PRINT*,ERRORCHECK(3)
           IF(M.EQ.15) PRINT*,ERRORCHECK(4)
           IF(M.EQ.4)  PRINT*,ERRORCHECK(5)
           PRINT*,''
           PRINT*,'Minimum Range = ',RANGEM(COCOUNT,M,1) 
           PRINT*,'Maximum Range = ',RANGEM(COCOUNT,M,2)
           PRINT*,'Variable Value= ',INDVAR(M)
           PRINT*,''
           PRINT*,'This variable should not be out of bounds'
           PRINT*,'of the prescribed ranges.  An error here is a good'
           PRINT*,'indication that the wrong parameterization has been'
           PRINT*,'chosen by the code.'
           PRINT*,''
           PRINT*,'STOPPED IN SR CALC_OH!'
           PRINT*,'*******************'
           STOP
         ENDIF
C************************
C
         ENDIF
C
             DIFFERR=INDVAR(M)-RANGEM(COCOUNT,M,2)
C
              OOA=M 
           IF(OOB.EQ.3.AND.M.EQ.18) OOA=19
C
             ERRORMINMAX=INDVAR(M)
             INDVAR(M)=RANGEM(COCOUNT,M,2)-1D-2*BMINSA
C
         IF(ERRORON.EQ.0.AND.DIFFERR.GT.1.D-8) THEN
           PRINT*,'*******************'
           PRINT*,'WARNING! RANGE OUT OF BOUNDS FOR:  ',VARNAMES(OOA)
           PRINT*,'Variable reset to value just prescribed range.' 
           PRINT*,'Minimum Range        = ',RANGEM(COCOUNT,M,1) 
           PRINT*,'Maximum Range        = ',RANGEM(COCOUNT,M,2) 
           PRINT*,'Variable Value - NEW = ',INDVAR(M)
           PRINT*,'Variable Value - OLD = ',ERRORMINMAX
           PRINT*,'Latitude = ', INDVAR(15)
           PRINT*,'Longitude= ', INDVAR(20)
           PRINT*,'Pressure = ', INDVAR(9)
           PRINT*,'Dec. >   = ', INDVAR(1)
           PRINT*,'NOx      = ', INDVAR(3)
c      if(INDVAR(20).GT.-135.and.INDVAR(20).lt.-90.and.INDVAR(15).
c     *  lt.60.and.INDVAR(15).gt.40.) print*,'bryan'
           PRINT*,'*******************'
         ENDIF
C
        ENDIF
C
        MIDWAY=0.9
C
         X(M)=-MIDWAY+
     *    (((INDVAR(M)-RANGEM(COCOUNT,M,1))/BMINSA)*MIDWAY*2.)
C
      ENDDO
C
C ***********************************************************
C ***********************************************************
C
       DO 2 KDONE=1,NDONE(COCOUNT)
C
C ***********************************************************
C 2) See if point within range of done function.
C ***********************************************************
C
          ID=IDENTOLD(COCOUNT,KDONE)
C
C *** Iterate number of inequalites (rows) this element.
C
          JKEEP=0
         NOROWS=IROWST(COCOUNT,ID,0)
C
         DO 30 KROW=1,NOROWS
            KKRW=IROWST(COCOUNT,ID,KROW)
C
C *** calculate plane equation values.
C
            FNCVAL=0.-ELTODO(COCOUNT,KKRW,MOVAR(COCOUNT)+2)

       DO 20 JVAR=1,MOVAR(COCOUNT)
           FNCVAL=FNCVAL+(ELTODO(COCOUNT,KKRW,JVAR)*INDVAR(JVAR))
   20  CONTINUE
C
C *** IF points above the "plane" are sought, verify that this point lies
C *** on or above plane.
C
            SYMBOL=ELTODO(COCOUNT,KKRW,MOVAR(COCOUNT)+1)
            IF (SYMBOL.GT.0..AND.FNCVAL.LT.0.) JKEEP=1

C *** If points below the "plane" are sought, verify that this point lies
C *** below plane.
C
            IF (SYMBOL.LT.0..AND.FNCVAL.GE.0.) JKEEP=1
C            IF (SYMBOL.LT.0..AND.FNCVAL.GT.0.) JKEEP=1

            IF(JKEEP.EQ.1) GOTO 2

   30    CONTINUE
C
C ***********************************************************
C 3) Calculate parameterized OH.
C ***********************************************************
C
         IF(JKEEP.EQ.0) THEN 

          NOTDONE=1
C
         do 3 MM=1,NENDW(COCOUNT,KDONE)
	   C(MM)=0.
           C(MM)=COEFF(COCOUNT,KDONE,MM)
 3       CONTINUE
C
           NDFUNCS=NDFUNCS+KDONE
C
           !--------------------------------------------------------------
           ! Turn off FUNC_IFS since we right now we are not using
           ! the parameterized OH (bmy, 4/14/00)
           !CALL FUNC_IFS(NDFUNCS,X,C,PARAMOH)
           !--------------------------------------------------------------
C
c          GOTO 2
          GOTO 33
C
         ENDIF
C
C ***********************************************************
 2      CONTINUE
C ***********************************************************
 33     CONTINUE
C Error Check
      IF(NOTDONE.EQ.0) THEN
        PRINT*,'SR CALC_OH:  Done function not found!!!!'
        STOP
      ENDIF
C
      RETURN
      END
