C $Id: getoh.f,v 1.1 2003/06/30 20:26:03 bmy Exp $
C*****************************************************************************
      SUBROUTINE GETOH(COLON,COLAT,COALT,FIRSTDT,BOH)
C*****************************************************************************
      IMPLICIT NONE
#     include "CMN_OH"
      LOGICAL FIRSTDT
      INTEGER I,J,K,L,M
      INTEGER NDFUNCS,NALL,NLOOK
      INTEGER COLON,COLAT,COALT,DOskip,OOB
      REAL*8 BOH(MLONBX,MLATBX,MVRTBX)
C
C*****************************************************************************
C Created by Bryan Duncan.
C*****************************************************************************
C This SR determines which parameterization is needed for an individual 
C box and calculates the parameterized [OH] (SR CALC_OH).
C
C The parameterizations are defined by latitude, pressure, [NOx], [ISOP]
C and season.  (The corresponding model variables are OH_LAT, OH_PRESS,
C OH_NOX, OH_ISOP and OH_SEASON, resp.)

            OH_LAT=INDVARA(COLON,COLAT,COALT,15)
            OH_LON=INDVARA(COLON,COLAT,COALT,20)
            OH_PRESS=INDVARA(COLON,COLAT,COALT,9)
            OH_NOX=INDVARA(COLON,COLAT,COALT,3)
            OH_ISOP=INDVARC(COLON,COLAT,COALT,18)
            OH_O3=INDVARA(COLON,COLAT,COALT,4)
C
C*****************************************************************************
C List of Variables & Arrays
C
C BOH = array holding parameterized OH.
C
C Dimensions of troposphere.  User specified values in CMN_OH.
C   COALT = position of box in vertical in troposphere
C   COLAT = position of box from pole to pole (following one longitude band)
C   COLON = position of box circling globe (following one latitude band)
C
C NDFUNCS = number of polynomials.  NDFUNCS is greater than the number
C           of parameterizations (NPARAM). A parameterization of
C           a subdomain may be described by more than one polynomial
C           due to domain divisions.
C
            NDFUNCS=0
C
C FIRSTDT = 0 on first chemistry time step.
C         = 1 on any time chemistry step, but first. 
C
C PARAMOH = parameterized [OH].
C
            PARAMOH=100.
C
C COCOUNT = identification number of parameterization.
C
            COCOUNT=0
C
C OOB     = bookkeeping for error check in SR CALC_OH.
C
            OOB=0
C
C INDVAR  = independent variables needed to calculate parameterized [OH].
C
            INDVAR(1:MXVAR)=0.
C
C DOskip  = Counter to skip parameterization if necessary.
C         = 0 : use parameterization to get [OH].
C         = 1 : set to low value and skip parameterization.
C         = 2 : use climatological [OH] (see SR READAVGOH & SR INTERPOH).
C         = 3 : use [OH] = 1x10^5 if isoprene > 0.6 ppbv and skip parameterization
C         = 4 : use [OH] = 1x10^4 if isoprene > 3 ppbv and skip parameterization
C
C Variables specific to box in question:
C
C OH_NOX    = NOx concentration
C OH_PRESS  = pressure
C OH_LAT    = latitude
C OH_ISOP   = isoprene concentration
C OH_LON    = longitude
C OH_MONTH = month (January-December = 1-12)
C OH_SEASON = (Northern Hemispheric) season
C        1 --> winter (Dec, Jan, Feb)
C        2 --> spring (Mar, Apr, May)
C        3 --> summer (Jun, Jul, Aug)
C        4 --> autumn (Sep, Oct, Nov)
C
C*****************************************************************************
C 1)  Determine which parameterization the box needs.
C*****************************************************************************
C Error Check.
C If you'd like to see the specifics for the model box you're
C calculating OH, change value of NLOOK.
C
              NLOOK=1
           IF(NLOOK.EQ.0) THEN

            print*,''
            print*,'OH_LAT   =',OH_LAT
            print*,'OH_LON   =',OH_LON
            print*,'OH_PRESS =',OH_PRESS
            print*,'OH_NOX   =',OH_NOX
            print*,'OH_ISOP  =',OH_ISOP
            print*,'OH_O3    =',OH_O3
            print*,'OH_SEASON=',OH_SEASON

           ENDIF
C
C*****************************************************************************
C SR SKIP prevents calling OH parameterizations in regions where OH is very low,
C resets the values of isoprene and NOx to keep from calling nonexistent
C parameterizations, and performs other miscellaneous functions.
C
      CALL SKIP(DOskip,COLON,COLAT,COALT)
C
C If DOskip = 1, set [OH] to low value (i.e., skip parameterization).
C
      IF(DOskip.EQ.1) GOTO 320
C
C In subdomains where OH is very low or negligible due to low
C  sunlight (e.g., high latitudes in winter), concentrations of OH are
C  set to climatological mean values as a function of latitude,
C  altitude and season. No parameterized OH is calculated. 
C
      IF(OH_O3.LT.1.) DOskip=2
      IF(DOskip.EQ.2) THEN
         CALL INTERPOH
         GOTO 320
      ENDIF
C
C*************
C Error Check.
c Surface
C       IF(OH_PRESS.LT.PRESSES(1)) goto 320
C
c MT
c       IF(OH_PRESS.GE.PRESSES(1)) goto 320
c       IF(OH_PRESS.LT.PRESSES(2)) goto 320
C
c UT
c       IF(OH_PRESS.GE.PRESSES(2)) goto 320
C
C
C*****************************************************************************
C B) REMOTE:  Unpolluted Domain:  Low NOx, No Isoprene (<10 ppt)
C             Surface Layer, MT, & UT (NOx <50 ppt)
C*****************************************************************************
C Surface & MT
      IF(OH_ISOP.LE.10..AND.OH_O3.LE.45.AND.
     *  OH_NOX.LE.CNOXS(7).AND.OH_PRESS.GE.PRESSES(2))THEN
C     
            CALL REMOTE(NDFUNCS)
C
C    Assign appropriate independent variables to be used in
C     OH calculation.
C
            IF(OH_PRESS.GE.PRESSES(1)) THEN
               INDVAR(1:MXVAR)=INDVARA(COLON,COLAT,COALT,1:MXVAR)
               OOB=1
            ELSE
               INDVAR(1:MXVAR)=INDVARB(COLON,COLAT,COALT,1:MXVAR)
               OOB=2
            ENDIF
C
            GOTO 321
C
      ENDIF
C UT
      IF(OH_ISOP.LE.10..AND.OH_O3.LE.100..AND.
     *  OH_NOX.LE.CNOXS(6).AND.OH_PRESS.LT.PRESSES(2))THEN
C
            CALL REMOTE(NDFUNCS)
C
C    Assign appropriate independent variables to be used in
C     OH calculation.
C
               INDVAR(1:MXVAR)=INDVARB(COLON,COLAT,COALT,1:MXVAR)
               OOB=2
C
            GOTO 321
C
      ENDIF
C
C*****************************************************************************
C END B) REMOTE
C*****************************************************************************
C
      NDFUNCS=NDFUNCS+NDFUNCSB1+NDFUNCSB2+NDFUNCSB3
C
C*****************************************************************************
C A) CLEAN:  Unpolluted Domain:  Low NOx, No Isoprene (<10 ppt) 
C            Surface Layer (NOx <300 ppt), Middle & Up Tropospheres
C*****************************************************************************
C
      IF(OH_ISOP.LE.10.) THEN
         IF((OH_NOX.LE.CNOXS(1).AND.OH_PRESS.GE.PRESSES(2)).OR.
     *        (OH_NOX.LE.CNOXS(2).AND.OH_PRESS.LT.PRESSES(2))) THEN
C     
            CALL CLEAN(NDFUNCS)
C
C    Assign appropriate independent variables to be used in
C     OH calculation.
C
            IF(OH_PRESS.GE.PRESSES(1)) THEN
               INDVAR(1:MXVAR)=INDVARA(COLON,COLAT,COALT,1:MXVAR)
               OOB=1
            ELSE
               INDVAR(1:MXVAR)=INDVARB(COLON,COLAT,COALT,1:MXVAR)
               OOB=2
            ENDIF
C
            GOTO 321
C     
         ENDIF
C
      ENDIF
C
C*****************************************************************************
C END A) CLEAN
C*****************************************************************************
C
      NDFUNCS=NDFUNCS+NDFUNCSA1+NDFUNCSA2+NDFUNCSA3
C
C*****************************************************************************
C C) CLEAN with ISOPRENE:  Low NOx, Isoprene (>10 ppt) 
C            Surface Layer (NOx <300 ppt), Mid Troposphere 
C            (NOx <300 ppt), Up Troposphere (NOx <500ppt)
C*****************************************************************************
C
      IF(OH_ISOP.GT.10.) THEN
C
         IF((OH_NOX.LE.CNOXS(1).AND.OH_PRESS.GE.PRESSES(2)).OR.
     *        (OH_NOX.LE.CNOXS(2).AND.OH_PRESS.LT.PRESSES(2))) THEN
C
           IF(OH_PRESS.GE.PRESSES(1)) THEN
            IF(DOskip.EQ.3) THEN
               BOH(COLON,COLAT,COALT)=1.E5
               GOTO 322
            ENDIF
C
            IF(DOskip.EQ.4) THEN
               BOH(COLON,COLAT,COALT)=1.E4
               GOTO 322
            ENDIF
           ENDIF
C
            CALL CLEANISOP(NDFUNCS)
C
            IF(OH_PRESS.GE.PRESSES(1)) THEN
               INDVAR(1:MXVAR)=INDVARC(COLON,COLAT,COALT,1:MXVAR)
               OOB=3
            ELSE
               INDVAR(1:MXVAR)=INDVARD(COLON,COLAT,COALT,1:MXVAR)
               OOB=4
            ENDIF
C     
            GOTO 321
C     
         ENDIF
C     
      ENDIF
C
C*****************************************************************************
C END C) CLEAN with ISOPRENE
C*****************************************************************************
C
      NDFUNCS=NDFUNCS+NDFUNCSC1+NDFUNCSC2+NDFUNCSC3
C
C*****************************************************************************
C D) MIDNOx SUBDOMAIN with ISOPRENE: 
C            Mid NOx (>300 ppt <1000), Isoprene
C            Surface (see F for Middle Troposphere)
C*****************************************************************************
C
      IF(OH_NOX.GT.CNOXS(1).AND.
     *     OH_NOX.LE.CNOXS(3).AND.OH_PRESS.GE.PRESSES(1)) THEN
C
         CALL MODPOLLUTED(NDFUNCS)
C
         INDVAR(1:MXVAR)=INDVARC(COLON,COLAT,COALT,1:MXVAR)
         OOB=3
C
         GOTO 321
C
      ENDIF
C     
C*****************************************************************************
C END D) MIDNOx SUBDOMAIN with ISOPRENE - Surface
C*****************************************************************************
C
      NDFUNCS=NDFUNCS+NDFUNCSD1
C
C*****************************************************************************
C E) HINOx SUBDOMAIN with ISOPRENE
C               Hi NOx (>1000), Isoprene
C*****************************************************************************
C
      IF(OH_NOX.GT.CNOXS(3).AND.OH_PRESS.GE.PRESSES(1)) THEN
C
         CALL HIPOLLUTED(NDFUNCS)
C
         INDVAR(1:MXVAR)=INDVARC(COLON,COLAT,COALT,1:MXVAR)
         OOB=3
C
         GOTO 321
C     
      ENDIF
C
C*****************************************************************************
C END E) HINOx SUBDOMAIN with ISOPRENE
C*****************************************************************************
C
      NDFUNCS=NDFUNCS+NDFUNCSE1
C
C*****************************************************************************
C F) MIDNOx SUBDOMAIN with ISOPRENE
C            Mid NOx (>300 ppt <1000), Isoprene
C            Middle Troposphere 
C****************************************************************************
C     
      IF(OH_NOX.GT.CNOXS(1).AND.OH_NOX.LE.CNOXS(3).AND.
     *   OH_PRESS.LT.PRESSES(1).AND.OH_PRESS.GE.PRESSES(2)) THEN
C     
         CALL MODPOLLUTED(NDFUNCS)
C     
         INDVAR(1:MXVAR)=INDVARD(COLON,COLAT,COALT,1:MXVAR)
         OOB=4
C     
         GOTO 321
C
      ENDIF
C
C*****************************************************************************
C END F) MIDNOx SUBDOMAIN with ISOPRENE - MT
C*****************************************************************************
C
      NDFUNCS=NDFUNCS+NDFUNCSF1
C
C*****************************************************************************
C Error Check.
C
      PRINT*,'BOX LIES NOWHERE!'
      PRINT*,'STOPPED IN SR GETOH!'
      PRINT*,'*********************'
      PRINT*,'OH_LAT=',OH_LAT
      PRINT*,'OH_PRESS=',OH_PRESS
      PRINT*,'OH_NOX=',OH_NOX
      PRINT*,'OH_ISOP=',OH_ISOP
      PRINT*,'OH_SEASON=',OH_SEASON
      PRINT*,'*********************'
      STOP
C     
C*****************************************************************************
C*****************************************************************************
C 2)  Calculate the parameterized OH.
C*****************************************************************************
C*****************************************************************************
C
 321  CALL CALC_OH(NDFUNCS,OOB)

 320  BOH(COLON,COLAT,COALT)=PARAMOH
C
C***********************************************************
C SR correctOH: Correct aerosol problem from chem1d.
C***********************************************************
      IF(BOH(COLON,COLAT,COALT).GT.0.) THEN
         CALL correctOH(BOH,COLON,COLAT,COALT)
      ENDIF
C***********************************************************
C
 322   RETURN
       END
