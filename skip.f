C $Id: skip.f,v 1.1 2003/06/30 20:26:07 bmy Exp $
C**********************************************************************
      SUBROUTINE SKIP(DOskip,COLON,COLAT,COALT)
C**********************************************************************
      IMPLICIT NONE
#     include "CMN_OH"
      INTEGER DOskip,COLON,COLAT,COALT,NODOMOD,NODOHI
C**********************************************************************
C   Created by Bryan Duncan.
C**********************************************************************
C SR SKIP performs the following functions:
C 1) Low OH :  Use if statements to keep calling OH parameterization
C              when OH is very low.  About 1/5 of the boxes
C              are skipped here.
C 2) RESET  :  Reset the values of isoprene and NOx to keep from calling
C              "nonexistent" parameterizations.
C 3) Other Miscellaneous Functions
C
C**********************************************************************
C List of Variables & Arrays
C
C DOskip = Counter to skip parameterization. 
C        = 0 : use parameterization to get [OH]
C        = 1 : set to low value and skip parameterization
C        = 2 : use average [OH] (see SR readavgOH)
C        = 3 : use [OH] = 1x10^5 if isoprene > 2 ppbv and skip parameterization
C        = 4 : use [OH] = 1x10^4 if isoprene > 5 ppbv and skip parameterization
C
            DOskip=0
C
C Dimensions of troposphere.  User specified values in CMN_OH.
C   COALT = position of box in vertical in troposphere
C   COLAT = position of box from pole to pole (following one longitude band)
C   COLON = position of box circling globe (following one latitude band)
C
C PARAMOH = parameterized [OH].
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
C***************************************************
C 1) LOW OH
C***************************************************
C SH (80-90S)
      IF(OH_LAT.LE.-80.) THEN
         COCOUNT=0
         DOskip=2
         GOTO 320
      ENDIF
C NH (80-90N)
      IF(OH_LAT.GE.80.) THEN
         COCOUNT=0
         DOskip=2
         GOTO 320
      ENDIF
C WIN NH (40-90N)
      IF(OH_LAT.GE.ALATS(6).AND.OH_SEASON.EQ.1) THEN
         COCOUNT=0
         DOskip=2
         GOTO 320
      ENDIF
C SUM & SPR SH (60-90S)
      IF(OH_LAT.LE.ALATS(2).AND.(OH_SEASON.EQ.3.OR.
     *   OH_SEASON.EQ.2)) THEN
         COCOUNT=0
         DOskip=2
         GOTO 320
      ENDIF
C AUT NH (60-90N)
      IF(OH_LAT.GE.ALATS(7).AND.OH_SEASON.EQ.4) THEN
         COCOUNT=0
         DOskip=2
         GOTO 320
      ENDIF
C
C***************************************************
C 2) RESET
C***************************************************
C***************************************************
C ISOPRENE
C***************************************************
C    [ISOP] > 2 ppb set [OH] = 1e10^5 and skip parameterization
C    [ISOP] > 5 ppb set [OH] = 1e10^4 and skip parameterization
C Assumptions based on:
C Mickley, L. J., D. J. Jacob, and D. Rind,
C Uncertainty in preindustrial abundance of tropospheric ozone:
C Implications for radiative forcing calculations.
C In press Journal of Geophysical Research, 2000.
C
c       IF(OH_ISOP.GE.2000.) THEN
       IF(OH_ISOP.GE.600.) THEN
         DOskip=3
         GOTO 320
       ENDIF
c       IF(OH_ISOP.GE.5000.) THEN
       IF(OH_ISOP.GE.3000.) THEN
         DOskip=4
         GOTO 320
       ENDIF
C
C These statements ensure isoprene parameterizations not called
C    for domains with no isoprene parameterizations.
C    40-90N WIN, 60-90N SPR, 60-90N AUT, 60-90S All Seasons
C
      IF(OH_LAT.GE.ALATS(NLATS-1).AND.OH_SEASON.EQ.1) OH_ISOP=1.
      IF(OH_LAT.GE.ALATS(NLATS).AND.OH_SEASON.NE.3) OH_ISOP=1.
      IF(OH_LAT.LE.ALATS(2).AND.OH_SEASON.EQ.2) OH_ISOP=1.
      IF(OH_LAT.LE.ALATS(2).AND.OH_SEASON.EQ.3) OH_ISOP=1.
      IF(OH_LAT.LE.ALATS(1)) OH_ISOP=1.
C
C***************************************************
C These statements ensure the MODERATELY POLLUTED parameterizations
C    not called for domains with no parameterization.
C    Surface and Middle Troposphere.
C***************************************************
C
C For surface.
C
      IF( OH_NOX.GT.CNOXS(1).AND.
     *  OH_NOX.LE.CNOXS(3).AND.OH_PRESS.GE.PRESSES(1)) THEN
C
         NODOMOD=0
C
       IF(OH_LAT.LE.ALATS(1)) NODOMOD=1
       IF(OH_LAT.LE.ALATS(2).AND.OH_LAT.GT.ALATS(1)) THEN
          IF(OH_SEASON.EQ.3.OR.OH_SEASON.EQ.2) NODOMOD=1
       ENDIF
       IF(OH_LAT.GT.ALATS(6).AND.OH_SEASON.EQ.1) NODOMOD=1
       IF(OH_LAT.GT.ALATS(7).AND.OH_SEASON.NE.3) NODOMOD=1
C
          IF(NODOMOD.NE.0) THEN
                OH_NOX=CNOXS(1)-10.
                INDVARA(COLON,COLAT,COALT,3)=OH_NOX
                INDVARB(COLON,COLAT,COALT,3)=OH_NOX
                INDVARC(COLON,COLAT,COALT,3)=OH_NOX
                INDVARD(COLON,COLAT,COALT,3)=OH_NOX
          ENDIF
C
      ENDIF
C
C******
C
C For Middle Troposphere.
C
      IF(OH_NOX.GT.CNOXS(1).AND.OH_PRESS.LT.PRESSES(1).
     *    AND.OH_PRESS.GE.PRESSES(2).AND.OH_NOX.LE.CNOXS(3)) THEN
C
         NODOMOD=0
C
       IF(OH_LAT.LE.ALATS(5)) NODOMOD=1
       IF(OH_LAT.GT.ALATS(6).AND.OH_SEASON.EQ.1) NODOMOD=1
       IF(OH_LAT.GT.ALATS(7)) NODOMOD=1
C
          IF(NODOMOD.NE.0) THEN
                OH_NOX=CNOXS(1)-10.
                INDVARA(COLON,COLAT,COALT,3)=OH_NOX
                INDVARB(COLON,COLAT,COALT,3)=OH_NOX
                INDVARC(COLON,COLAT,COALT,3)=OH_NOX
                INDVARD(COLON,COLAT,COALT,3)=OH_NOX
          ENDIF
C
      ENDIF
C
C No moderatly polluted parameterization done for upper trop.
C
      IF(OH_NOX.GT.CNOXS(3).AND.OH_PRESS.LT.PRESSES(2)) THEN
         DOskip=1
         GOTO 320
      ENDIF
C
C***************************************************
C These statements ensure the E) HEAVILY POLLUTED parameterizations
C    not called for domains with no parameterization.
C    Very hi NOx associated with low OH, especially in wintertime.
C***************************************************
C
C Winter months:  high NOx=low OH.
C
      IF(OH_NOX.GT.CNOXS(3)) THEN
        IF(OH_LAT.GE.ALATS(7).AND.OH_SEASON.EQ.2) THEN
            PARAMOH=0.
            DOskip=1
         GOTO 320
        ENDIF
        IF(OH_LAT.GE.ALATS(5).AND.OH_SEASON.EQ.1) THEN
            PARAMOH=0.
            DOskip=1
         GOTO 320
        ENDIF
        IF(OH_LAT.LE.ALATS(3).AND.OH_SEASON.EQ.3) THEN
            PARAMOH=0.
            DOskip=1
         GOTO 320
        ENDIF
      ENDIF
C
C No heavily polluted parameterization done for middle and upper trop.
C
      IF(OH_NOX.GT.CNOXS(4).AND.OH_PRESS.LT.PRESSES(1)) THEN
         DOskip=1
         GOTO 320
      ENDIF
C
      IF(OH_NOX.GT.CNOXS(4)) THEN
C
C Assume hi NOx = low OH for areas without heavily polluted parameterizations.
C  
       IF(OH_LAT.LT.ALATS(3).OR.OH_LAT.GE.ALATS(7).OR.
     *    (OH_LAT.GE.ALATS(3).AND.OH_LAT.LT.ALATS(4).AND.
     *    OH_SEASON.EQ.1).OR.(OH_LAT.GE.ALATS(6).AND.
     *    OH_SEASON.LE.2).OR.(OH_LAT.GE.ALATS(5).AND.
     *    OH_SEASON.EQ.1)) THEN

         DOskip=1
         GOTO 320
C
       ELSE
C
C Reset NOx < 5 ppbv in regions with heavily polluted parameterizations.
C
          OH_NOX=CNOXS(4)-10.
          INDVARA(COLON,COLAT,COALT,3)=OH_NOX
          INDVARB(COLON,COLAT,COALT,3)=OH_NOX
          INDVARC(COLON,COLAT,COALT,3)=OH_NOX
          INDVARD(COLON,COLAT,COALT,3)=OH_NOX
C
       ENDIF
C
      ENDIF
C
C Make sure not to call nonexistent parameterization. 
C
      IF(OH_NOX.GT.CNOXS(3).AND.OH_NOX.LT.CNOXS(4)) THEN
C       
         NODOHI=0
C
       IF(OH_LAT.LE.ALATS(3))  NODOHI=1
       IF(OH_LAT.LE.ALATS(4).AND.OH_SEASON.EQ.1) NODOHI=1
       IF(OH_LAT.GT.ALATS(7))  NODOHI=1
       IF(OH_LAT.GT.ALATS(6).AND.OH_SEASON.EQ.2) NODOHI=1

       IF(NODOHI.NE.0) THEN
          OH_NOX=CNOXS(3)-10.
          INDVARA(COLON,COLAT,COALT,3)=OH_NOX
          INDVARB(COLON,COLAT,COALT,3)=OH_NOX
          INDVARC(COLON,COLAT,COALT,3)=OH_NOX
          INDVARD(COLON,COLAT,COALT,3)=OH_NOX
       ENDIF
C
      ENDIF
C
C***************************************************
C 4) Miscellaneous Functions.
C***************************************************
C PRESS > 100 mb:  Don't calculate OH over 100 mb.

       IF(OH_PRESS.LT.PRESSES(3)) THEN
C
C Error Check.
C
         PRINT*,'************************************************'
         PRINT*,'WARNING MESSAGE - SR SKIP!'
         PRINT*,'The pressure for this box is less than 100 mb.'
         PRINT*,'No parameterization is designed for pressures'
         PRINT*,'less than 100 mb. Box is skipped and [OH] set'
         PRINT*,'to low value. '
         PRINT*,'LON BOX=',COLON
         PRINT*,'LAT BOX=',COLAT
         PRINT*,'ALT BOX=',COALT
         PRINT*,'************************************************'
C
         DOskip=1
C
         GOTO 320
C
       ENDIF
C
C NOX < 1 ppt:  Reset NOx to 1 ppt.
C
       IF(OH_NOX.LT.1.) THEN
         OH_NOX=1.1
         INDVARA(COLON,COLAT,COALT,3)=OH_NOX
         INDVARB(COLON,COLAT,COALT,3)=OH_NOX
         INDVARC(COLON,COLAT,COALT,3)=OH_NOX
         INDVARD(COLON,COLAT,COALT,3)=OH_NOX
       ENDIF
C
C***************************************************
C
 320  RETURN
      END
