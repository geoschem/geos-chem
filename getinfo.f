C* $Id: getinfo.f,v 1.1 2003/06/30 20:26:05 bmy Exp $
C ****************************************************************************
      SUBROUTINE GETINFO(OH_MONTH2,ALBD)
C ****************************************************************************
C THIS SUBROUTINE IS TO BE MODIFIED BY THE USER!
C The user can interface it with a model by passing necessary
C information through the subroutine's argument list or
C through the user's own common blocks.
C ****************************************************************************
C ****************************************************************************
C In this subroutine, the user supplies the parameterization code
C with information it needs to calculate the 24-hour averaged [OH].
C ****************************************************************************
C Created by Bryan Duncan.
C ****************************************************************************
C The user is to put whatever common blocks and declaration 
C statements here which are necessary to supply this subroutine
C with the model input it requires.
C
      ! References to F90 modules
      USE DAO_MOD,   ONLY : AIRVOL
      USE ERROR_MOD, ONLY : GEOS_CHEM_STOP

      IMPLICIT NONE

#     include "CMN_SIZE"       ! Size parameters
#     include "CMN"            ! STT
#     include "CMN_CO"         ! CO arrays
#     include "CMN_CO_BUDGET"  ! FMOL_CO

      REAL*8   PI,DEC,A0,A1,A2,A3,B1,B2,B3,R,BLAT,BLON,METHVAR
      PARAMETER ( PI = 3.141592 )
      INTEGER  II,KK,NDOCO,OH_MONTH2,NCHECKIT
      REAL*8   AISOP,ddd,ALBD(IIPAR,JJPAR),sumo3up

      ! Weight of air (taken from "comode.h") 
      REAL*8, PARAMETER :: WTAIR = 28.966d0
C
C ****************************************************************************
C This common block and these declaration statements are not to be
C modified by the user.
C
#     include "CMN_OH"
      INTEGER I,J,K,L,NSTOP
C
C ****************************************************************************
C The polynomials that are used to calculate the parameterized
C OH concentrations are a function of the following independent
C variables:
C
C   Independent Variables                   Units
C   -----------------------------------------------------------
C   Solar Declination Angle                 degrees
C   Total Ozone Column Above                DU
C   Total Nitrogen Oxides                   pptv
C   Ozone                                   ppbv
C   Carbon Monoxide                         ppbv
C   Methane                                 ppbv
C   Water Vapor                             ppmv
C   Acetone                                 pptv
C   Pressure                                mb
C   Propene                                 pptv
C   Ethane                                  pptv
C   Propane                                 pptv
C   ALK4                                    pptv
C   Surface Albedo                          unitless
C   Latitude                                degrees
C   Temperature                             K
C   Cloude Albedo - Above                   unitless
C   Cloude Albedo - Below                   unitless
C   Isoprene (daylight)                     pptv
C   -----------------------------------------------------------
C
C All these independent variables need to be 24-hour averages!!
C
C All the independent variables are REAL*8!
C
C ****************************************************************************
C
C The independent variables are stored in the arrays INDVARA-D:
C
C INDVARA is used in surface layer parameterizations without isoprene.
C INDVARB is used in middle and upper tropospheric layer parameterizations
C         without isoprene.
C INDVARC is used in surface layer parameterizations with isoprene.
C INDVARD is used in middle and upper tropospheric layer parameterizations
C         with isoprene.
C
C --------------------------------------------------------------------------
C   MXVAR    INDVARA        INDVARB        INDVARC         INDVARD
C --------------------------------------------------------------------------
C     1       DEC            DEC            DEC             DEC
C     2       O3COL          O3COL          O3COL           O3COL
C     3       NOt            NOt            NOt             NOt
C     4       O3             O3             O3              O3
C     5       CO             CO             CO              CO
C     6       CH4            CH4            CH4             CH4
C     7       H2O            H2O            H2O             H2O
C     8       ACET           ACET           ACET            ACET
C     9       PRESS          PRESS          PRESS           PRESS
C    10       PRPE           PRPE           PRPE            PRPE
C    11       ETHA           ETHA           ETHA            ETHA
C    12       PROP           PROP           PROP            PROP
C    13       ALK4           ALK4           ALK4            ALK4
C    14       ALB            ALB            ALB             ALB
C    15       LAT            LAT            LAT             LAT
C    16       TEMP           TEMP           TEMP            TEMP
C    17       CLOUDA         CLOUDA         CLOUDA          CLOUDA
C    18         *            CLOUDB         ISOP            CLOUDB
C    19         *              *              *             ISOP
C --------------------------------------------------------------------------
C   *No independent variable value placed in this position.
C
C ****************************************************************************
C Error Check.  If the user does not totally fill in the INDVARA-D arrays,
C then the code will stop to tell the user.
C
      INDVARA(:,:,:,:)=-1000.D0
      INDVARB(:,:,:,:)=-1000.D0
      INDVARC(:,:,:,:)=-1000.D0
      INDVARD(:,:,:,:)=-1000.D0
C
      INDVARA(:,:,:,18:19)=0.D0
      INDVARB(:,:,:,19)=0.D0
      INDVARC(:,:,:,19)=0.D0
C
      print*,'Filling INDVAR in SR getinfo.f'
C
C ****************************************************************************
C
C (1) Solar declination angle (low precision formula)
C     JDAY = julian day
C
      A0 = 0.006918
      A1 = 0.399912
      A2 = 0.006758
      A3 = 0.002697
      B1 = 0.070257
      B2 = 0.000907
      B3 = 0.000148
 
      R  = 2.* PI * float(JDAY-1) / 365.
C
      DEC = A0 - A1*COS(  R) + B1*SIN(  R)
     &         - A2*COS(2*R) + B2*SIN(2*R)
     &         - A3*COS(3*R) + B3*SIN(3*R)
C
      DEC = DEC*180./PI
C
      IF(DEC.GT.24..OR.DEC.LT.-24.) THEN
         PRINT*,'Stopped in SR GETINFO!'
         PRINT*,'  Declination angle is too high!'
         PRINT*,'DEC = ',DEC
         CALL GEOS_CHEM_STOP
      ENDIF
C
      INDVARA(:,:,:,1)=DEC
C
C (2) O3 Column Above (DU)
C
c      INDVARA(:,:,:,2)=o3up(:,:,:)

c this is temporary- only using the first longitude
c to keep consistent with 2-D nature of fastj o3 col.
! Make sure DO-loops are ordered L-J-I (bmy, 4/16/00)
cbnd      DO L=1,MVRTBX
cbnd      DO J=1,MLATBX
cbnd      DO I=1,MLONBX
cbnd         INDVARA(I,J,L,2)=o3up(1,J,L)
cbnd      ENDDO
cbnd      ENDDO
cbnd      ENDDO

c to keep consistent with 2-D nature of fastj o3 col.
       INDVARA(:,:,:,2)=0.D0

      DO J=1,MLATBX
      DO L=1,MVRTBX
        sumo3up=0.
      DO I=1,MLONBX
        sumo3up=sumo3up+o3up(I,J,L)
      IF(I.EQ.MLONBX) THEN
        sumo3up=sumo3up/REAL(MLONBX)
        IF(sumo3up.LT.5..and.L.LE.19) THEN
         print*,'stopped in SR getinfo'
         print*,'o3 column problem!'
         print*,I,L,sumo3up
         CALL GEOS_CHEM_STOP
        ENDIF
        DO II=1,MLONBX
          INDVARA(II,J,L,2)=sumo3up
        ENDDO
      ENDIF
      ENDDO
      ENDDO
      ENDDO
C
C Rescale O3 column to match TOMS data.
C   O3Col(L) = O3Col(L) * TOMSO3Col / O3Col(1)
C
      IF(NCLIMATOLOGY3.EQ.1) THEN
         DO I=1,MLONBX
         DO J=1,MLATBX
         DO L=1,MVRTBX

      INDVARA(II,J,L,2)=INDVARA(I,J,L,2)*ctm4x5month(I,J)/
     *                              INDVARA(I,J,1,2)

C Error check.
      IF(INDVARA(II,J,L,2).LT.1..OR.INDVARA(II,J,L,2).GT.600.) THEN
       PRINT*,'STOPPED IN SR GETINFO.F!'
       PRINT*,'The total o3 column (DU) =',INDVARA(II,J,L,2)
       PRINT*,'This could mean there is an error.'
       CALL GEOS_CHEM_STOP
      ENDIF

         ENDDO
         ENDDO
         ENDDO
      ENDIF

C
C (3) NOt (pptv)
C        NOt = NO + NO2 + 2N2O5 + NO3 + HNO2 + HNO4
C
C      INDVARA(:,:,:,3)=BBIJ(:,:,:,1)*1.D12
C
      INDVARA(:,:,:,3)=(BBIJ(:,:,:,1)+BBIJ(:,:,:,10)*2.D0+
     *     BBIJ(:,:,:,11))*1.D12
C       
C (4) O3 (ppbv)
C
      INDVARA(:,:,:,4)=BBIJ(:,:,:,8)*1.E9

      ! Make sure DO loops are ordered L-J-I (bmy, 4/16/00)
      DO L=1,MVRTBX
      DO J=1,MLATBX
      DO I=1,MLONBX
         IF(INDVARA(I,J,L,4).LT.1.) THEN
            INDVARA(I,J,L,4)=INDVARA(I,J,L-1,4) 
            IF(INDVARA(I,J,L,4).LT.2.) THEN
               PRINT*,'Stopped in SR Getinfo.f'
               PRINT*,'O3 < 1 ppb!'
               CALL GEOS_CHEM_STOP
            ENDIF
         ENDIF
      ENDDO
      ENDDO
      ENDDO
C
C (5) CO (ppbv)
C
      !NDOCO = 1
      NDOCO = 0 
C
      IF(NDOCO.EQ.0) THEN
C Convert STT(CO) from kgs/box to ppb.
C
         STTTOPPB(:,:,:)=1.D0/AD(:,:,:)*WTAIR/FMOL_CO/1000.D0*1.D9
C
C Variable CO.
C
         INDVARA(:,:,:,5)=STT(:,:,:,1)*STTTOPPB(:,:,:)
C 
      ELSE
C
C Monthly mean CO.
C
         INDVARA(:,:,:,5)=BBIJ(:,:,:,9)*1.E9
C
      ENDIF
C
C (6) CH4 (ppbv)
C
      DO L=1,MVRTBX
      DO J=1,MLATBX
      DO I=1,MLONBX

         METHVAR=0.D0

         IF(J.GE.1.AND.J.LE.12)  METHVAR=1645.D0
         IF(J.GE.13.AND.J.LE.23) METHVAR=1655.D0
         IF(J.GE.24.AND.J.LE.34) METHVAR=1715.D0
         IF(J.GE.35.AND.J.LE.46) METHVAR=1770.D0
         IF(METHVAR.EQ.0.) THEN
            PRINT*,'STOPPED IN SR GETINFO! Methane = 0!'
            CALL GEOS_CHEM_STOP
         ENDIF

         METHVAR=1700.D0

         INDVARA(:,J,:,6)=METHVAR
C
      ENDDO
      ENDDO
      ENDDO
C
C (7) H2O (ppmv)
C
      INDVARA(:,:,:,7)=Wavg(:,:,:)
C
C (8) Acetone (pptv) 
C
      INDVARA(:,:,:,8)=BBIJ(:,:,:,4)*0.333D0*1.E12
C
C (9) Pressure (mb) - 24 hr avg
C
      WHERE(Pavg(:,:,:).GT.1020.D0) Pavg(:,:,:)=1019.9D0

      DO L=1,MVRTBX
         INDVARA(:,:,L,9)=Pavg(:,:,L)  
      ENDDO
C
C (10) Propene (pptv)
C
      INDVARA(:,:,:,10)=BBIJ(:,:,:,5)*0.333D0*1.E12
C
C (11) Ethane (pptv)
C
      INDVARA(:,:,:,11)=BBIJ(:,:,:,7)*0.5D0*1.E12
C
C (12) Propane (pptv)
C
      INDVARA(:,:,:,12)=BBIJ(:,:,:,6)*0.333D0*1.E12
C
C (13) ALK4 (pptv)
C     
      INDVARA(:,:,:,13)=BBIJ(:,:,:,2)*0.25D0*1.E12
C
C (14) Albedo (fraction of 1)
C
Cvis       WHERE(AVGF(:,:).LT.0.06D0) AVGF(:,:)=0.061D0
      WHERE(ALBD(:,:).LT.0.06D0) ALBD(:,:)=0.061D0
      
Cvis       WHERE(AVGF(:,:).GT.0.80D0) AVGF(:,:)=0.79D0
      WHERE(ALBD(:,:).GT.0.80D0) ALBD(:,:)=0.79D0

      DO L=1,MVRTBX
Cvis       INDVARA(:,:,L,14)=AVGF(:,:)
         INDVARA(:,:,L,14)=ALBD(:,:)
      ENDDO
C
C (15) Latitude (deg = -90 to 90)
C
      BLAT=-91.D0
      DO J=1,MLATBX
         IF(J.NE.MLATBX.AND.J.NE.1) THEN
            BLAT=BLAT+4.D0
         ELSE
            BLAT=BLAT+2.D0
         ENDIF
         IF(J.EQ.2) BLAT=BLAT-1.D0
         IF(J.EQ.MLATBX) BLAT=BLAT+1.D0
         INDVARA(:,J,:,15)=BLAT
      ENDDO
C
C (16) Temperature (K)
C
      DO L=1,MVRTBX
         INDVARA(:,:,L,16)=Tavg(:,:,L)  
      ENDDO
C
C (17) Cloud flux fraction above box
C
      WHERE(Ravga(:,:,:).LT.0.D0) Ravga(:,:,:)=0.D0

      WHERE(Ravga(:,:,:).GT.0.6D0) Ravga(:,:,:)=0.59999D0

      INDVARA(:,:,:,17)=Ravga(:,:,:)
C*****************************************************************************
C
C The first 17 independent variables in INDVARA are the
C same for INDVARB, INDVARC and INDVARD.  However, they
C are different for the 18-19 independent variables.
C
      INDVARB(:,:,:,1:17)=INDVARA(:,:,:,1:17)
      INDVARC(:,:,:,1:17)=INDVARA(:,:,:,1:17)
      INDVARD(:,:,:,1:17)=INDVARA(:,:,:,1:17)
C
C*****************************************************************************
C
C (18) Cloud flux fraction below box
C
      WHERE(Ravgb(:,:,:).LT.0.D0) Ravgb(:,:,:)=0.D0

      WHERE(Ravgb(:,:,:).GT.0.6D0) Ravgb(:,:,:)=0.59999D0

      INDVARB(:,:,:,18)=Ravgb(:,:,:)
      INDVARD(:,:,:,18)=Ravgb(:,:,:)

c*****************************************************************************
C Resetting cloud flux fractions to 0 for test!
c      INDVARA(:,:,:,17)=0.0001
c      INDVARB(:,:,:,17)=0.0001
c      INDVARC(:,:,:,17)=0.0001
c      INDVARD(:,:,:,17)=0.0001
c      INDVARB(:,:,:,18)=0.0001
c      INDVARD(:,:,:,18)=0.0001
c*****************************************************************************
C
C (19) Isoprene (pptv) - daylight average?
C
      DO L=1,MVRTBX
      DO J=1,MLATBX
      DO I=1,MLONBX
C
         AISOP=BBIJ(I,J,L,3)*0.2D0*1.E12
C
c are you reading in 24-Hr averages or daylight averages?
c see CO_fillfields.f


         IF(AISOP.LT.2.D0) AISOP=2.D0

c        INDVARC(I,J,L,18)=DLOG(AISOP)
c        INDVARD(I,J,L,19)=DLOG(AISOP)
         INDVARC(I,J,L,18)=AISOP
         INDVARD(I,J,L,19)=AISOP

      ENDDO
      ENDDO
      ENDDO
C
C ****************************************************************************
C Error Check. Print out input variables.
C
       NCHECKIT=1
C
      IF(NCHECKIT.EQ.0) THEN       
C
      DO I=1,MLONBX
      DO J=1,MLATBX
      DO L=1,MVRTBX
      DO KK=1,19
cc       if(INDVARD(I,J,L,3).le.350..and.
cc     *    INDVARD(I,J,L,19).gt.11.) then     
      WRITE(819,*) INDVARD(I,J,L,KK)
cc       else
cc         ddd=0. 
cc        write(819,*) ddd
cc       endif
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      CALL GEOS_CHEM_STOP
C
      ENDIF
C
C ***************************************************************************
C ***************************************************************************
C Other parameters needed.
C ***************************************************************************
C ***************************************************************************
C
C A) Longitude (deg = -180 - 180)
C   Longitude is not an independent variable used in the polynomials
C   to calculate the parameterized OH, but is used to locate 
C   regions of biomass burning.  There are separate parameterizations
C   for regions affected by biomass burning.  For convenience, however,
C   longitude is stored in the INDVAR* arrays.
C
      BLON=-185.D0
      DO I=1,MLONBX
         BLON=BLON+5.D0

         INDVARA(I,:,:,20)=BLON
         INDVARB(I,:,:,20)=BLON
         INDVARC(I,:,:,20)=BLON
         INDVARD(I,:,:,20)=BLON
      ENDDO
C
C B) Season
C     OH_MONTH = month (January-December = 1-12)
C     OH_SEASON = Northern Hemispheric season
C        1 --> winter (Dec, Jan, Feb)
C        2 --> spring (Mar, Apr, May)
C        3 --> summer (Jun, Jul, Aug)
C        4 --> autumn (Sep, Oct, Nov)
C
C From month, determine season.  The parameterizations are
C functions of season and not month.
C
      OH_MONTH=OH_MONTH2
      IF(OH_MONTH.LE.2)  OH_SEASON=1
      IF(OH_MONTH.EQ.12) OH_SEASON=1
      IF(OH_MONTH.LE.5.AND.OH_MONTH.GE.3)  OH_SEASON=2
      IF(OH_MONTH.LE.8.AND.OH_MONTH.GE.6)  OH_SEASON=3
      IF(OH_MONTH.LE.11.AND.OH_MONTH.GE.9) OH_SEASON=4
C     
C*****************************************************************************
C Error Check.
C     
      NSTOP=0
C
      DO I = 1,MXVAR
      DO J = 1,MVRTBX
      DO K = 1,MLATBX
      DO L = 1,MLONBX
C
         IF(INDVARA(L,K,J,I).LT.-200.) THEN
            PRINT*,'******************************************'
            PRINT*,'ERROR IN SR GETINFO!'
            PRINT*,'The array INDVARA is not completely filled.'
            PRINT*,'Check independent variable number = ',I
            PRINT*,'INDVARA(L,K,J,I)=',INDVARA(L,K,J,I)
            NSTOP=1
         ENDIF
C
         IF(INDVARB(L,K,J,I).LT.-200.) THEN
            PRINT*,'******************************************'
            PRINT*,'ERROR IN SR GETINFO!'
            PRINT*,'The array INDVARB is not completely filled.'
            PRINT*,'Check independent variable number = ',I
            PRINT*,'INDVARB(L,K,J,I)=',INDVARB(L,K,J,I)
            NSTOP=1
         ENDIF
C
         IF(INDVARC(L,K,J,I).LT.-200.) THEN
            PRINT*,'******************************************'
            PRINT*,'ERROR IN SR GETINFO!'
            PRINT*,'The array INDVARC is not completely filled.'
            PRINT*,'Check independent variable number = ',I
            PRINT*,'INDVARC(L,K,J,I)=',INDVARC(L,K,J,I)
            NSTOP=1
         ENDIF
C
         IF(INDVARD(L,K,J,I).LT.-200.) THEN
            PRINT*,'******************************************'
            PRINT*,'ERROR IN SR GETINFO!'
            PRINT*,'The array INDVARD is not completely filled.'
            PRINT*,'Check independent variable number = ',I
            PRINT*,'INDVARD(L,K,J,I)=',INDVARD(L,K,J,I)
            NSTOP=1
         ENDIF
C
         IF(NSTOP.EQ.1) THEN
            PRINT*,'L,K,J=',L,K,J
            PRINT*,'******************************************'
            STOP
         ENDIF
C
      ENDDO
      ENDDO
      ENDDO
      ENDDO
C
      IF(OH_SEASON.LT.1.OR.OH_SEASON.GT.4) THEN
         PRINT*,'******************************************'
         PRINT*,'ERROR IN SR GETINFO!'
         PRINT*,'OH_SEASON does not have a value between 1 & 4.'
         PRINT*,'OH_SEASON = ',OH_SEASON
         PRINT*,'******************************************'
         CALL GEOS_CHEM_STOP
      ENDIF
C
C*****************************************************************************
C
      RETURN
      END
