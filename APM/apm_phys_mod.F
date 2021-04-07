#ifdef APM
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: apm_phys_mod
!
! !DESCRIPTION: Module APM\_PHYS\_MOD solves APM microphysics.
!\\
!\\
! !INTERFACE:
!
      MODULE APM_PHYS_MOD
!
! !USES:
!
      IMPLICIT NONE
      PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!  
      PUBLIC  :: APM_PHYS
!
! !PRIVATE MEMBER FUNCTIONS:
!
      PRIVATE :: DACTS
      PRIVATE :: GETGAMMAN2O5

! !REVISION HISTORY: 
!  23 Aug 2008 - F. Yu       - Initial version
!  08 Nov 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
      CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: apm_phys
!
! !DESCRIPTION: Subroutine APM\_PHYS is the driver routine for the APM
!  microphysics package.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE APM_PHYS( II,        JJ,      LL,       NCOAG1, 
     &                     NCOAG2,    IACT10,  IACT20,   IACT30,
     &                     RACT1,     RACT2,   RACT3,    NTEMPOUT1,
     &                     DENAIR,    PRESS,   YSIGMA,   TK,      RH,
     &                     XQ, PLVSOG01, PLVSOG1, CACID,  PACID, CNH3,
     &                     DT,        MMSA,    MNIT,     MNH4,
     &                     MBCS,      MOCS,    MDSTS,    MSALTS,
     &                     XMBC, XMOC, SOAT,   CLVSOG,   MSULFLV,
     &                     MBCLV,     MOCLV,   MDSTLV,   MSALTLV,
     &                     GFTOT1,    GFTOT2,  DENWET1,  DENWET2,
     &                     YSPGF ,    XBCLIFE, XOCLIFE,
     &                     VZ,        YCLDLIQ, XCDN,     XCDNSP,
     &                     XM1D,      XN1D,    TEMPOUT1, XMDST, FCLOUD1,
!OPT+     &                     FCLOUD1 )
     &                     ZBEXT,ZW,ZG,ZBABS,XBEXT1k,
     &                     YBEXT,YW,YG,ACS,XDMA,AERAREA,AERDRYR,
     &                     GAMMAPM,AEROCOMOUT,IFSITEOUT,
     &                     ATOM4N)
!
! !USES:
!
      USE APM_INIT_MOD, ONLY : IFNUCL,IFAG,IFATHN,IFNH3,IFNUCLORG
      USE APM_INIT_MOD, ONLY : NSO4,NSEA,NDSTB,NTYP,NBCOC
      USE APM_INIT_MOD, ONLY : RDRY, VDRY, RSALT, VSALT, YGF
      USE APM_INIT_MOD, ONLY : RBCOC, VBCOC
      USE APM_INIT_MOD, ONLY : TOTNUMBC,TOTNUMOC, TOTAREABC,TOTAREAOC
      USE APM_INIT_MOD, ONLY : DACT1, DACT2, DACT3, DENSULF, XMLVSOG
      USE APM_INIT_MOD, ONLY : V1LVSOG,V1ACID
      USE APM_INIT_MOD, ONLY : RDST,DENDST,VDST
      USE APM_INIT_MOD, ONLY : ONEPI
      USE APM_TIMN_MOD, ONLY : YUJTIMN1
      USE APM_ATHN_MOD, ONLY : YUJATHN

      USE APM_COAG_MOD, ONLY : APM_COAG,APM_COAGSCAV
      USE APM_GROW_MOD, ONLY : APM_GROW,APM_MOVEBIN
      USE APM_INIT_MOD, ONLY : ISITES,JSITES

      USE APM_OPTI_MOD, ONLY : APM_OPT, APM_OPT_LW   ! OPT+
      USE APM_INIT_MOD, ONLY : CEMITBCOC2 !OPT+
      USE APM_INIT_MOD, ONLY : IFCOAT !OPT+
      USE APM_INIT_MOD, ONLY : IFCOATBC !OPT+
      USE APM_INIT_MOD, ONLY : IFOPT, IFCDNPDF !OPT+
      USE APM_INIT_MOD, ONLY : DENBC   ! in kg/m3
      USE APM_INIT_MOD, ONLY : DTTEST
      USE APM_INIT_MOD, ONLY : BC_LIFEAPM,OC_LIFEAPM
      USE APM_INIT_MOD, ONLY : D0CONV 
      USE APM_INIT_MOD, ONLY : WBAR
      USE APM_INIT_MOD, ONLY : IFAMINEUP, AGAMA
      USE APM_INIT_MOD, ONLY : IFAEROCOMOUT

      USE apm_mixactivate, only: activate
!
! !INPUT PARAMETERS: 
!
      INTEGER :: II
      INTEGER :: JJ
      INTEGER :: LL
      INTEGER :: NTEMPOUT1
!
! !INPUT/OUTPUT PARAMETERS: 
! 
      INTEGER :: NCOAG1,NCOAG2, NCOAG4, NCOAG5
      INTEGER :: IACT10       ! bin index for cloud act 
      INTEGER :: IACT20       ! bin index for cloud act 
      INTEGER :: IACT30       ! bin index for cloud act 
      INTEGER :: RACT1        ! bin index for cloud act 
      INTEGER :: RACT2        ! bin index for cloud act 
      INTEGER :: RACT3        ! bin index for cloud act 
      REAL*8  :: VZ,YCLDLIQ, XCDN,XCDNSP
      REAL*8  :: DENAIR
      REAL*8  :: PRESS
      REAL*8  :: TK
      REAL*8  :: RH
      REAL*8  :: XQ
      REAL*8  :: PLVSOG01
      REAL*8  :: PLVSOG1
      REAL*8  :: CACID
      REAL*8  :: PACID
      REAL*8  :: CNH3 
      REAL*8  :: DT
      REAL*8  :: MMSA
      REAL*8  :: MNIT
      REAL*8  :: MNH4
      REAL*8  :: MBCS         ! mass of sulfate attached to primary particles
      REAL*8  :: MOCS         ! mass of sulfate attached to primary particles
      REAL*8  :: MDSTS        ! mass of sulfate attached to primary particles
      REAL*8  :: MSALTS       ! mass of sulfate attached to primary particles
      REAL*8  :: SOAT
      REAL*8  :: CLVSOG
      REAL*8  :: MSULFLV
      REAL*8  :: MBCLV
      REAL*8  :: MOCLV
      REAL*8  :: MDSTLV
      REAL*8  :: MSALTLV
      REAL*8  :: GFTOT1
      REAL*8  :: GFTOT2
      REAL*8  :: DENWET1
      REAL*8  :: DENWET2
      REAL*8  :: XM1D(NSO4+NSEA),XMDST(NDSTB),XMBC(NBCOC),XMOC(NBCOC)
      REAL*8  :: XN1D(NSO4)
      REAL*8  :: TEMPOUT1(NTEMPOUT1)
      REAL*8  :: FCLOUD1(NSO4+4)

      REAL*8  :: YSPGF,  XBCLIFE, XOCLIFE
      REAL*8  :: YSIGMA
      REAL*8  :: ACS, XDMA, YJATHN
      REAL*8  :: AERAREA(NSO4+NSEA+NDSTB+4), AERDRYR(NSO4+NSEA+NDSTB+4)
      REAL*8  :: PM25(8),GAMMAPM(NTYP),AEROCOMOUT(27)

      REAL*8  :: CHOM,XJN,XJI,YRATIO,YJBHN

      INTEGER :: IFSITEOUT
!
! !REMARKS:
!  APM AEROSOL Types (N=1,NTYP)
!  N=1:  Sulfate or Secondary particles (SO4 plus other species) 
!         (density: 1.7 g/cc)
!  N=2:  Sea Salt (density: 2.2 g/cc)
!  N=3:  DUST  (density: 2.5 g/cc for D<1 um, 2.65 g/cm3 for D>~1 um)
!  N=4:  BC  (density: ? 1.0 or 1.8 g/cc)
!  N=5:  OC  (density: 1.8 g/cc)
! 
! !REVISION HISTORY: 
!  17 Mar 2010 - F. Yu       - Initial version  
!  08 Nov 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER :: N,M,ITYPE, NI,NJ,NA
      INTEGER :: IRH
      REAL*8  :: PMB,XS,CACID1
      REAL*8  :: XN(NSO4),XMA(NSO4),XVA(NSO4)
      REAL*8  :: XNSALT(NSEA),XMSALT(NSEA),XVSALT(NSEA)
      REAL*8  :: XNDST(NDSTB),XVDST(NDSTB)
      REAL*8  :: XVBC(NBCOC),XVOC(NBCOC)
      REAL*8  :: RBCWET(NBCOC),ROCWET(NBCOC)
      REAL*8  :: YAREABC(NBCOC),YAREAOC(NBCOC)
      REAL*8  :: RBCGF(NBCOC),ROCGF(NBCOC)
      REAL*8  :: TOTN
      REAL*8  :: YJ, DVJ,  DV
      REAL*8  :: VGAS,FREEP,TAREA,AREA,YKN,FCORR
      REAL*8  :: YCS(NTYP), YCSDUST, TCS, SURF(NTYP)
      REAL*8  :: CSSULF,AREASULF,XRCM,YJAVE,YRSTAR,CACID0
      REAL*8  :: YYJ,YYJAVE,YJATHNAVE,CACIDAVE
      REAL*8  :: MSO4,MSOA
      REAL*8  :: MSO4B,MNITB,MNH4B,MMSAB,MSULFT,TOTMP,FSO4B
      REAL*8  :: MSP(NTYP),MCORE(NTYP),ZK(NTYP)
      REAL*8  :: ZDACT(NTYP,4),ZCCN(NTYP,3),ZTN(NTYP),TCCN(3)
      REAL*8  :: ATOM4N(4)
      REAL*8  :: DIAM1GAS,DIAM2GAS
      REAL*8  :: TNBC, TNOC,TNBC1,TNBC2,TNOC1,TNOC2
      REAL*8  :: GF2,GF3,GF4,FRACSOA,GFWATERVOL,GFWATER
      REAL*8  :: GFTOT,DENWET,DENWATER
      REAL*8  :: DENSALTWET,Y1,DENBCWET,DENOCWET
      REAL*8  :: TN10nm,TNOTHER
      REAL*8  :: RRATIO, RRATIOWET
      REAL*8  :: AREABC,AREAOC
      REAL*8  :: RWET(NSO4),RSALTWET(NSEA),RWETCM(NSO4),RSALTWETCM(NSEA)
      REAL*8  :: GFGAS, CSSALT,AREASALT
      REAL*8  :: MCONDOTH,TCSOTHER, TOTCONDOTH
      REAL*8  :: YAREASALT(NSEA), MSALTTOT, YSALTS(NSEA), MSALTSOA(NSEA)
      REAL*8  :: RSALTGF(NSEA),TAREASALT
      REAL*8  :: YAREADST(NDSTB),CSDST,AREADST,RDSTWET(NDSTB)
      REAL*8  :: FBCOC1,FBCOC2
      REAL*8  :: EPC, COAGN0, DTCOAG1,DTCOAG2,DTCOAG4,DTCOAG5
      INTEGER :: NCOAG(5), NCOAGMAX, I3nm, IERF, I10nm
      INTEGER :: I50nm(5),I80nm(5),I120nm(5)
      REAL*8  :: BINSIZE(5),BINMASS(5)
      REAL*8  :: RLOSULF
      REAL*8  :: SUMXVA0,SUMXVA1,DLVSOG,DACID
      REAL*8  :: XMCONDIN,CACIDIN,PACIDIN,TCSOTHERIN
      REAL*8  :: DTLEFT, DTNGC, DNMAX
      INTEGER :: IFNGC
      INTEGER :: ISITE,JSITE
      REAL*8  :: ZN1,ZN2,ZN3
      REAL*8  :: TCSLV,CLVSOG1,CLVSOG0
      INTEGER :: ICOND

!      REAL*8    :: MBCFF,MBCBB,MOCFF,MOCBB
      REAL*8    :: TK1
!OPT+
      INTEGER   :: ITEST,MBIN,I,IWL

      INTEGER, PARAMETER :: MWLS=16   !RRTMG 14 SW band + two wls (0.5 and 0.55 um)
      INTEGER, PARAMETER :: MWLL=9    !RRTMG 9 LW band

!      INTEGER, PARAMETER :: KWL1=6   !550 nm
      INTEGER, PARAMETER :: KWL1=4   !500 nm
!      INTEGER, PARAMETER :: KWL1=3   !390 nm
      REAL*8    :: WLS(MWLS)
      REAL*8    :: XBEXT(MWLS),XBEXT1k(40,5),AER1k(12)
      REAL*8    :: XW(MWLS),XG(MWLS)
      REAL*8    :: YBEXT(NTYP,MWLS),YW(NTYP,MWLS),YG(NTYP,MWLS)
      REAL*8    :: ZBEXT(MWLS),ZW(MWLS),ZG(MWLS)
!LW  for dust only
      REAL*8    :: WLL(MWLL)
      REAL*8    :: XBEXT1(MWLL),XW1(MWLL),XG1(MWLL)
!      REAL*8    :: ZBEXT1(MWLL),ZW1(MWLL),ZG1(MWLL)
      REAL*8    :: ZBABS(MWLL)

      REAL*8    :: vtotal,yvdry,vh2o
      REAL*8    :: xn_in(NSO4)
      REAL*8    :: XDCORE(NSO4)
      REAL*8    :: ZFV(6),ZFV_SP(6)

      REAL*8    :: RRATIOWETBC,RRATIOWETOC
      REAL*8    :: XNBC(NBCOC),XNOC(NBCOC), VRATIO
      REAL*8    :: YBCS,MBCSOA,MBCTOT,YOCS,MOCSOA,MOCTOT,CSBC,CSOC
      REAL*8    :: TAREAO, ZTEMP1, ZTEMP2
      REAL*8    :: MBCS0,MOCS0,MBCLV0,MOCLV0
      REAL*8    :: FLUXBC, FLUXOC,DRDTBC,DRDTOC
!OPT+
      REAL*8    :: RGTOT,XMTOT,RGM,SIG2TOT,SIG,YM(NSO4)
      REAL*8    :: MWATERSP
      REAL*8    :: DUST500

! for max sat
      INTEGER   :: NSIZE_AER(NTYP), IS, IACT40
      REAL*8    :: RRATIOBC, RRATIOOC, maxsatout,XSAT
      REAL*8    :: NAEROSOL(NSO4,NTYP),AM(NSO4,NTYP)
      REAL*8    :: HYGRO_AER(NSO4,NTYP)
      REAL*8    :: YCDN(NTYP)
      REAL*8    :: ZCDN(NTYP)

! for BHN_IIN
      REAL*8  :: YRHF, YAMOLF,YNTOT,YISTAR,YNWSTAR

      REAL*8  :: XN0, XNH3
      REAL*8  :: XJBN,XJBI,XJTN,XJTI
      REAL*8  :: XJBH,XJBIM,XJTH,XJTIM
      REAL*8  :: RBH,RBIM,RTH,RTIM
      REAL*8 :: YJAVE_BH,YJAVE_BIM,YJAVE_TH,YJAVE_OM,YJAVE_PON,YJAVE_POI
      REAL*8 :: YJAVE_TIM
      REAL*8 :: YJAVE_BN,YJAVE_BI,YJAVE_TN,YJAVE_TI

      REAL*8  :: YJOM,YJPON,YJPOI
      REAL*8  :: TAREABC,TAREAOC

!
! for PDF
      INTEGER   :: IW
      INTEGER, PARAMETER :: MW = 40
      real*8 :: DW,WW,SIGMA
      real*8 :: W(MW),PW(MW)

      REAL*8    :: XRATIO

      REAL*8  :: DENAER(5)
!YuBCden      DATA (DENAER(NA),NA=1,5)/1.7,2.2,2.65,1.8,1.8/  ! density (g/cm3)
!      DATA (DENAER(NA),NA=1,5)/1.7,2.2,2.65,1.2,1.8/  ! density (g/cm3)

!      DATA (WLS(IWL),IWL=1,MWLS)/0.3,0.55,0.94,1.785,3.19/  !OPT+
!      DATA (WLS(IWL),IWL=1,MWLS)/0.34,0.38,0.443,0.469,0.5,0.554,
!     &     0.645,0.675, 0.865,0.94,1.02,1.24,1.64,1.785,2.13,3.19/
      DATA (WLS(IWL),IWL=1,MWLS)/0.23,0.30,0.39,0.50,0.53,0.55,
     &       0.70,1.01,1.27,1.46,1.78,2.05,2.33,2.79,3.46,8.02/
      DATA (WLL(IWL),IWL=1,MWLL)/4.3,5.,6.,8.1,9.6,11.6,15.8,24.0,35./
!OPT+

!Luo for AEROCOM
      real*8 :: XSATAEROCOM(13)
      DATA XSATAEROCOM /0.0005,0.00075,0.001,0.0015,0.002,0.003,0.004
     &                 ,0.005,0.006,0.007,0.008,0.009,0.01/

      ATOM4N=0.d0

!      CLVSOG0 = 4.d6
      CLVSOG0 = 5.d5

      DENWATER = 1.0  ! density of water (g/cm3)
      DENAER(1)= 1.7
      DENAER(2)= 2.2
      DENAER(3)= 2.65
      DENAER(4) = DENBC/1000.0  !update BC density based on value from input.apm
      DENAER(5)= 1.8

      IF(PACID.LT.0.0) THEN
         WRITE(6,*) "PACID < 0, set it to 1.E-2"
         WRITE(6,99)II,JJ,LL,CACID,PACID
         PACID = 1.d-2
      ENDIF

      TEMPOUT1=1.d-30

      XN0 = 1.d-6*PRESS/(1.3807d-23*TK)   !#/cm3
      XNH3 = CNH3*1.d-12*XN0  !ppt to #/cm3

      IF(IFNH3.EQ.0) XNH3=1.d-30   ! no NH3 effect

      MCORE = 1.d-20
      MSP   = 1.d-21
      ZCCN  = 1.d-20
      PM25  = 1.d-10

!     
      MBCS0 = MBCS
      MOCS0 = MOCS
      MBCLV0 = MBCLV
      MOCLV0 = MOCLV

! The growth/movebin/coag subrountine is based on volume, convert mass to volume
! Sulfate
      MSO4 = 0.d0
      DO N=1,NSO4
        XN(N)=XN1D(N)
        XMA(N)=XM1D(N)
        MSO4 = MSO4 + XM1D(N)   ! total bin sulfate mass
        IF(XMA(N).LT.1.D-40)XMA(N)=1.D-40
        XVA(N)=XMA(N)*1.d-3/DENSULF   !XVA in cm3/cm3  XMA in kg/m3
      ENDDO

! Move particles across bins after cloud chem, XN is the values
! recorded before calling DO_CHEMISTRY in main.f
      CALL APM_MOVEBIN(NSO4,XN,XVA)  
!
! Seasalt
!
      DO N=1,NSEA
        XMSALT(N) = XM1D(NSO4+N)
        MCORE(2)=MCORE(2) + XMSALT(N)
        IF(XMSALT(N).LT.1.D-40)XMSALT(N)=1.D-40
        XVSALT(N)=XMSALT(N)*1.d-3/DENAER(2)   !XV in cm3/cm3  XM in kg/m3
        XNSALT(N) = XVSALT(N)/(1.E6*VSALT(N)) ! XN in #/cm3, VSALT in m3

        IF(RSALT(N).LT.1.25d-6) THEN
         PM25(7)=PM25(7) + XMSALT(N)*1.d9
        ENDIF

        IF(RSALT(N).GE.1.35d-9)THEN
        IF(RSALT(N).LE. 2.4d-6)THEN
         IF(RSALT(N).LT.10.d-9)THEN
          ATOM4N(1)=ATOM4N(1)+XNSALT(N)
         ELSE IF(RSALT(N).LT.60.d-9)THEN
          ATOM4N(2)=ATOM4N(2)+XNSALT(N)
         ELSE IF(RSALT(N).LT.0.5d-6)THEN
          ATOM4N(3)=ATOM4N(3)+XNSALT(N)
         ELSE
          ATOM4N(4)=ATOM4N(4)+XNSALT(N)
         ENDIF
        ENDIF
        ENDIF
      ENDDO
      ZTN(2) = SUM(XNSALT)
!
! Dust
      DO N=1,NDSTB
        MCORE(3)=MCORE(3) + XMDST(N)
        XVDST(N)=XMDST(N)*1.d-3/DENDST(N)   !XV in cm3/cm3,DENDST in g/cm3
        XNDST(N) = XVDST(N)/(1.E6*VDST(N)) ! XN in #/cm3, VDST in m3
        RDSTWET(N)=RDST(N)*100.   ! wet size in cm, use RDST for now

        IF(RDST(N).LT.1.25d-6) THEN
         PM25(8)=PM25(8) + XMDST(N)*1.d9
        ENDIF

        IF(RDST(N).GE.1.35d-9)THEN
        IF(RDST(N).LE. 2.4d-6)THEN
         IF(RDST(N).LT.10.d-9)THEN
          ATOM4N(1)=ATOM4N(1)+XNDST(N)
         ELSE IF(RDST(N).LT.60.d-9)THEN
          ATOM4N(2)=ATOM4N(2)+XNDST(N)
         ELSE IF(RDST(N).LT.0.5d-6)THEN
          ATOM4N(3)=ATOM4N(3)+XNDST(N)
         ELSE
          ATOM4N(4)=ATOM4N(4)+XNDST(N)
         ENDIF
        ENDIF
        ENDIF
      ENDDO
      ZTN(3) = SUM(XNDST)
      DUST500 = SUM(XNDST(6:NDSTB))

! BCOC
      DO N=1,NBCOC
        MCORE(4)=MCORE(4) + XMBC(N)
        XVBC(N)=XMBC(N)*1.D-3/DENAER(4)  !XV in cm3/cm3,DEN in g/cm3
        XNBC(N) = XVBC(N)/(1.E6*VBCOC(N)) ! XN in #/cm3, VBC in m3

        IF(RBCOC(N).LT.1.25D-6) THEN
         PM25(5)=PM25(5) + XMBC(N)*1.d9
        ENDIF
      ENDDO
      ZTN(4) = SUM(XNBC)

      DO N=1,NBCOC
        MCORE(5)=MCORE(5) + XMOC(N)
        XVOC(N)=XMOC(N)*1.D-3/DENAER(5)   !XV in cm3/cm3,DEN in g/cm3
        XNOC(N) = XVOC(N)/(1.E6*VBCOC(N)) ! XN in #/cm3, VOC in m3

        IF(RBCOC(N).LT.1.25D-6) THEN
         PM25(6)=PM25(6) + XMOC(N)*2.1d9
        ENDIF
      ENDDO
      ZTN(5) = SUM(XNOC)

!******************************************************************************
! Sulfate particle dry size increase due to uptake of NIT, NH4, SOA via
! equilibrium/partition  (ratio same for all sizes for now)
! assume SO4, NIT, NH4, SOA have same density for now

! MSULFLV teated as a part of MSO4 for now but not involved in isoropia calculation
      MSO4B = MSO4-MSULFLV
      IF(MSO4B.LE.0.) THEN
        WRITE(6,*)"MSO4.LE.MSULFLV",II,JJ,LL,MSO4,MSULFLV
! Reduce MSULFLV for now
        MSULFLV = 0.99*MSO4
        MSO4B = MSO4-MSULFLV
      ENDIF
      MCORE(1) = MSO4B

      DO N=1,NTYP
         IF(MCORE(N).LE.0.) THEN
!           WRITE(6,*)"MCORE.LE.0", N, II, JJ,LL
           MCORE(N)=1.d-20
         ENDIF
      ENDDO

! NIT, NH4, and MSA on SP
! When call inorganic equilibrium, MSULFT is used. Need scale to get
! NIT, NH4, MSA associated with SP.
      MSULFT = MSO4B+MBCS+MOCS+MDSTS+MSALTS  !Total SulfateNew  (kg/m3)
      MSULFT = MAX(MSULFT,1.d-20)
      FSO4B = MSO4B/MSULFT   ! fraction of SO4 in SP
      MNITB = MNIT * FSO4B
      MNH4B = MNH4 * FSO4B
      MMSAB = MMSA * FSO4B

! SV-SOA, MV-SOA on SP
! MOC*2.1+MSULFLV+MBCLV+MOCLV+MDSTLV+MSALTLV was used to get SOAT

      TOTMP = MCORE(5)*2.1+MSULFLV+MBCLV+MOCLV+MDSTLV+MSALTLV 
      TOTMP = MAX(TOTMP,1.d-20)
      MSOA = SOAT*MSULFLV/TOTMP ! SOA partioned into SP (kg/m3)

      MSP(1) = MSO4+MMSAB+MNITB+MNH4B+MSOA  ! total mass of SP
      MSP(2) = MSALTS*(1.+(MMSA+MNIT+MNH4)/MSULFT) ! SP mass on sea salt
     &        +MSALTLV*(1.+SOAT/TOTMP)  
      MSP(3) = MDSTS*(1.+(MMSA+MNIT+MNH4)/MSULFT) ! SP mass on dust
     &        +MDSTLV*(1.+SOAT/TOTMP)  
      MSP(4) = MBCS*(1.+(MMSA+MNIT+MNH4)/MSULFT) ! SP mass on BC
     &        +MBCLV*(1.+SOAT/TOTMP)  
      MSP(5) = MOCS*(1.+(MMSA+MNIT+MNH4)/MSULFT) ! SP mass on OC
     &        +MOCLV*(1.+SOAT/TOTMP)  
!     &        +MCORE(5)*2.*SOAT/TOTMP
     &        +MCORE(5)*2.1*SOAT/TOTMP

      XRATIO  = (MSO4B+MBCS+MOCS)/MSULFT
      PM25(1) = (MSO4B + MBCS+MOCS)*1.d9
      PM25(2) = (MNIT * XRATIO )*1.d9
      PM25(3) = (MNH4 * XRATIO)*1.d9
      PM25(4) = (SOAT + MSULFLV+MBCLV+MOCLV)*1.d9

!
! determine cloud activation dry diameters at three S based on composition

      XSAT = 0.001 ! assume 0.1%, will be updated in the next call
      CALL DACTS(TK,MSO4,MSULFLV,MMSAB,MNITB,MNH4B,MSOA,
     &              MCORE,MSP,XSAT,ZK,ZDACT)

!      CALL GETGAMMAN2O5(TK, RH, MSO4, MSULFLV, MMSAB,
!     &                  MNITB, MNH4B, MSOA, MCORE,
!     &                  MSP, GAMMAPM )

      IF(MSO4.LE.0.) THEN
         write(6,99)II,JJ,LL,MSO4,MSULFLV,MMSA,MNIT,MNH4,MSOA
         GFGAS = 1.
      ELSE
! MSULFLV included in MSO4
         GFGAS=(MSP(1)/MSO4)**(1./3.) 
      ENDIF

! for max sat
      ITYPE = 1
      NSIZE_AER(ITYPE)=NSO4
      DO N= 1, NSIZE_AER(ITYPE)
         NAEROSOL(N,ITYPE)=XN(N)*1.0d6  ! #/m3
         AM(N,ITYPE) = GFGAS*RDRY(N)    ! m
         HYGRO_AER(N,ITYPE)=ZK(ITYPE)
      ENDDO
!
! log-normal fitting to SP 
      RGTOT = 0.
      XMTOT = 1.d-20
      DO N= 1, NSO4                ! 
        YM(N) = GFGAS**3.0*XVA(N)*DENSULF*1.d3  !XVA in cm3/cm3  YM in kg/m3
        RGTOT = RGTOT + YM(N)*DLOG10(GFGAS*RDRY(N))
        XMTOT = XMTOT + YM(N)
      ENDDO
      RGM =10.**(RGTOT/XMTOT)   ! mass weighted median R in m

      SIG2TOT = 0.
      DO N= 1,NSO4
       SIG2TOT=SIG2TOT+YM(N)*DLOG10(GFGAS*RDRY(N)/RGM)
     &                     *DLOG10(GFGAS*RDRY(N)/RGM)
      ENDDO
      SIG = 10.**(SQRT(SIG2TOT/XMTOT))  ! standard deviation


! Growth due to uptake of water, need update for LV-SOA later
      IRH = INT(RH+0.5)
!      IRH= MIN0(99,IRH)
!OPT+      IRH= MIN0(99,IRH)
!OPT+ GF calculate large uncertainty at high RH, also should avoid AOD at the
! presence of cloud. set RH to a max of 95% for now.
      IRH= MIN0(99,IRH)
      IRH= MAX0(1,IRH)
!      GF2 = YGF(IRH,2)  ! growth factor of SO4+NIT+NH4 component 
! Based on ISOROPIA
      GF2 = YSPGF

      GF3 = YGF(IRH,3)  ! growth factor of SOA component 
      FRACSOA = MSOA/MSP(1)
      GFWATERVOL = (1.-FRACSOA)*GF2**3.0 + FRACSOA*GF3**3.0
      GFWATER = GFWATERVOL**(1./3.)
!
! Total growth factor
      GFTOT = GFGAS*GFWATER
      IF(GFTOT<1.D0)THEN
         WRITE(*,*)'GFTOT < 1, set to 1',II,JJ,LL,GFTOT
         GFTOT = 1.D0
      ENDIF
      GFTOT1 = GFTOT

      DO N=1,NSO4
         RWET(N) = RDRY(N)*GFTOT  ! Consider uptake of NIT,NH4,SOA and H20
         RWETCM(N) = RWET(N)*100.
      ENDDO

! Average density of wet particles
      DENWET=DENSULF/GFWATERVOL+DENWATER*(1.-1./GFWATERVOL)
      DENWET1 = DENWET
! Then find corresponding act bins 
      I3nm=1
      I10nm=1
      IACT10=1
      IACT20=1
      IACT30=1
      RACT1=1
      RACT2=1
      RACT3=1

      DO N   = 2, NSO4
         DIAM1GAS = GFGAS*RDRY(N-1)*2.
         DIAM2GAS = GFGAS*RDRY(N)*2.
         IF(3.d-9.GE.DIAM1GAS.and.3.d-9.LT.DIAM2GAS) I3nm=N
         IF(1.d-8.GE.DIAM1GAS.and.1.d-8.LT.DIAM2GAS) I10nm=N
         IF(ZDACT(1,1).GE.DIAM1GAS.and.ZDACT(1,1).LT.DIAM2GAS) IACT10=N
         IF(ZDACT(1,2).GE.DIAM1GAS.and.ZDACT(1,2).LT.DIAM2GAS) IACT20=N
         IF(ZDACT(1,3).GE.DIAM1GAS.and.ZDACT(1,3).LT.DIAM2GAS) IACT30=N

        IF(GFGAS*RDRY(N).GE.1.35d-9)THEN
        IF(GFGAS*RDRY(N).LE. 2.4d-6)THEN
         IF(GFGAS*RDRY(N).LT.10.d-9)THEN
          ATOM4N(1)=ATOM4N(1)+XN(N)
         ELSE IF(GFGAS*RDRY(N).LT.60.d-9)THEN
          ATOM4N(2)=ATOM4N(2)+XN(N)
         ELSE IF(GFGAS*RDRY(N).LT.0.5d-6)THEN
          ATOM4N(3)=ATOM4N(3)+XN(N)
         ELSE
          ATOM4N(4)=ATOM4N(4)+XN(N)
         ENDIF
        ENDIF
        ENDIF
      ENDDO
      RACT1=IACT20
!
      ZCCN(1,1)=SUM(XN(IACT10:NSO4))
      ZCCN(1,2)=SUM(XN(IACT20:NSO4))
      ZCCN(1,3)=SUM(XN(IACT30:NSO4))
      ZTN(1) = SUM(XN(I3nm:NSO4))

!
! Wet size and density of seasalt
!
      GF4 = YGF(IRH,4)  ! growth factor of seasalt component 
      IF(GF4.LT.1.D0)THEN
         WRITE(*,*)'GF4 < 1',II,JJ,LL,GF4
         GF4 = 1.0
      ENDIF
      GFTOT2= GF4

      DO N=1,NSEA
         RSALTWET(N) = RSALT(N)*GF4
         RSALTWETCM(N) = RSALTWET(N)*100.
      ENDDO
      Y1 = 1./(GF4**3.)
      DENSALTWET = DENAER(2)*Y1 + DENWATER*(1. - Y1)
      DENWET2 = DENSALTWET

! Wetsize of coated BCOC -- rough calculation for now
      IF(MCORE(4).GT.0.D0)THEN
         RRATIO = (1. + MSP(4)/MCORE(4))**(1./3.)
         RRATIOWET=(1.+MSP(4)*GFWATER**3./MCORE(4))**(1./3.)
      ELSE
         RRATIO=1.D0
         RRATIOWET=1.D0
      ENDIF
      RRATIOBC = RRATIO
      RRATIOWETBC = RRATIOWET

      IF(MCORE(5).GT.0.D0)THEN
         RRATIO = (1. + MSP(5)/MCORE(5))**(1./3.)
!bug         RRATIOWET=(1.+MOCS*(GFGAS*GFWATER)**3./MCORE(5))**(1./3.)
         RRATIOWET=(1.+MSP(5)*GFWATER**3./MCORE(5))**(1./3.)
      ELSE
         RRATIO=1.D0
         RRATIOWET=1.D0
      ENDIF

! OPT+
      RRATIOOC = RRATIO
      RRATIOWETOC = RRATIOWET
! OPT+

      DO N=1,NBCOC
        RBCWET(N)=RBCOC(N)*100.*RRATIOWETBC   ! wet size in cm
        ROCWET(N)=RBCOC(N)*100.*RRATIOWETOC   ! wet size in cm
      ENDDO
      Y1 = 1./(RRATIOWETBC**3.)
      DENBCWET = DENAER(4)*Y1 + DENWET*(1. - Y1)
      Y1 = 1./(RRATIOWETOC**3.)
      DENOCWET = DENAER(5)*Y1 + DENWET*(1. - Y1)

!******************************************************************************
! Condensation sink, coagulation sink  of various aerosols
!
      PMB = PRESS/100.    ! convert pa to mb
      VGAS = SQRT(216.03*TK)*100.  ! cm/s
      FREEP = PRESS*300./(1.013E5*TK)* 7.5E-6   ! H2SO4 mean free path (cm)

      TAREA = 1.d-20
      TCS = 1.d-20
      YAREASALT = 1.d-20
      DO N = 1, NTYP
        YCS(N) = 1.d-21
         IF(N.EQ.1) THEN   ! sulfate - use size distr to calculate CS
           CSSULF = 1.d-21
           AREASULF = 1.d-21
           DO NI = 1, NSO4
              XRCM = RWET(NI)*100.   ! cm
              AREA = 4.*ONEPI*XRCM*XRCM*XN(NI)  !cm2/cm3

              AERAREA(NI)=AREA
              AERDRYR(NI)=XRCM

              AREASULF = AREASULF + AREA
              YKN = FREEP/XRCM
              FCORR = YKN/(0.75+YKN) 
              CSSULF = CSSULF + 0.25*VGAS*AREA*FCORR  ! s-1
           ENDDO
           TAREA = TAREA + AREASULF
           SURF(N) = AREASULF
           YCS(N) = CSSULF
         ELSEIF(N.EQ.2) THEN   ! sea salt- use size distr to calculate CS
           CSSALT = 1.d-21
           AREASALT = 1.d-21
           DO NI = 1, NSEA
              XRCM = RSALTWET(NI)*100.   ! cm
              AREA = 4.*ONEPI*XRCM*XRCM*XNSALT(NI)  !cm2/cm3

              AERAREA(NSO4+NI)=AREA
              AERDRYR(NSO4+NI)=XRCM

              AREASALT = AREASALT + AREA
              YAREASALT(NI) = AREA
              YKN = FREEP/XRCM
              FCORR = YKN/(0.75+YKN) 
              CSSALT = CSSALT + 0.25*VGAS*AREA*FCORR  ! s-1
           ENDDO
           TAREA = TAREA + AREASALT
           SURF(N) = AREASALT
           YCS(N) = CSSALT
         ELSEIF(N.EQ.3) THEN   ! dust
           CSDST = 1.d-21
           AREADST = 1.d-21
           DO NI = 1, NDSTB
              XRCM = RDSTWET(NI)   ! cm
              AREA = 4.*ONEPI*XRCM*XRCM*XNDST(NI)  !cm2/cm3

              AERAREA(NSO4+NSEA+NI)=AREA
              AERDRYR(NSO4+NSEA+NI)=XRCM

              AREADST = AREADST + AREA
              YAREADST(NI) = AREA
              YKN = FREEP/XRCM
              FCORR = YKN/(0.75+YKN)
              CSDST = CSDST + 0.25*VGAS*AREA*FCORR  ! s-1
           ENDDO
           TAREA = TAREA + AREADST
           SURF(N) = AREADST
           YCS(N) = CSDST
         ELSEIF(N.EQ.4) THEN   ! BC
           CSBC = 1.D-21
           AREABC = 1.D-21
           DO NI = 1, NBCOC
              XRCM = RBCWET(NI)   ! cm
              AREA = 4.*ONEPI*XRCM*XRCM*XNBC(NI)  !cm2/cm3
              AREABC = AREABC + AREA
              YAREABC(NI) = AREA
              YKN = FREEP/XRCM
              FCORR = YKN/(0.75+YKN)
              CSBC = CSBC + 0.25*VGAS*AREA*FCORR  ! s-1
           ENDDO
           TAREA = TAREA + AREABC
           SURF(N) = AREABC
           YCS(N) = CSBC
           IF(IFCOATBC.EQ.0) THEN ! don't consider coating on BC, set CS to zero so that condensation and scavenging ignored
            YCS(N) = 1.D-20
           ENDIF

           SURF(N) = AREA
         ELSEIF(N.EQ.5) THEN   ! POC
           CSOC = 1.D-21
           AREAOC = 1.D-21
           DO NI = 1, NBCOC
              XRCM = ROCWET(NI)   ! cm
              AREA = 4.*ONEPI*XRCM*XRCM*XNOC(NI)  !cm2/cm3

              AREAOC = AREAOC + AREA
              YAREAOC(NI) = AREA
              YKN = FREEP/XRCM
              FCORR = YKN/(0.75+YKN)
              CSOC = CSOC + 0.25*VGAS*AREA*FCORR  ! s-1
           ENDDO
           TAREA = TAREA + AREAOC
           SURF(N) = AREAOC
           YCS(N) = CSOC

           SURF(N) = AREA
         ENDIF
         TCS = TCS + YCS(N)   ! Total CS 
      ENDDO
      XS =  TAREA*1.E8   ! XS is total surface area in um2/cm3
      TCSOTHER = TCS - CSSULF  ! CS of particles other than sulfate

! CS for amine uptake
      ACS = 0.
      IF(IFAMINEUP.EQ.1) THEN
        ACS = YCS(1) + YCS(2)*MSP(2)/(MSP(2)+MCORE(2)) 
     &               + YCS(3)*MSP(3)/(MSP(3)+MCORE(3))
     &               + YCS(4)*MSP(4)/(MSP(4)+MCORE(4))
     &               + YCS(5)*MSP(5)/(MSP(5)+MCORE(5))
        ACS = ACS * AGAMA   !AGAMA amine accommodation coef
      ENDIF

      IF(ZTN(2).GT.1.) THEN
       TAREASALT =SUM(YAREASALT)
       DO N=1, NSEA
! distribute MSALTS according to surface area
         YSALTS(N) = MSALTS * YAREASALT(N)/TAREASALT  
         MSALTSOA(N) = MSALTLV * YAREASALT(N)/TAREASALT  
         MSALTTOT = XMSALT(N) 
     &              + YSALTS(N)*GFGAS**3.0   ! sulfate and associated NIT, CH4, SOA
     &              + MSALTSOA(N)           ! SV-SOA not considered for now
         RSALTGF(N) = RSALT(N)*(MSALTTOT/XMSALT(N))**(1./3.)
       ENDDO
      ELSE
       RSALTGF = RSALT
      ENDIF

      N=1
      IF((2.*RSALTGF(N)).GE.ZDACT(2,1))
     &  ZCCN(2,1) = ZCCN(2,1) + XNSALT(N)
      IF((2.*RSALTGF(N)).GE.ZDACT(2,2))
     &  ZCCN(2,2) = ZCCN(2,2) + XNSALT(N)
      IF((2.*RSALTGF(N)).GE.ZDACT(2,3))
     &  ZCCN(2,3) = ZCCN(2,3) + XNSALT(N)

      DO N=2,NSEA  ! Sea-salt CCN at 3 S
         IF(ZDACT(2,2).GE.(2.*RSALTGF(N-1)).and.
     &      ZDACT(2,2).LT.(2.*RSALTGF(N))) RACT2=N
         IF((2.*RSALTGF(N)).GE.ZDACT(2,1)) 
     &     ZCCN(2,1) = ZCCN(2,1) + XNSALT(N)
         IF((2.*RSALTGF(N)).GE.ZDACT(2,2)) 
     &     ZCCN(2,2) = ZCCN(2,2) + XNSALT(N)
         IF((2.*RSALTGF(N)).GE.ZDACT(2,3)) 
     &     ZCCN(2,3) = ZCCN(2,3) + XNSALT(N)
      ENDDO

! for max sat
      ITYPE = 2
      NSIZE_AER(ITYPE)=NSEA
      DO N= 1, NSIZE_AER(ITYPE)
         NAEROSOL(N,ITYPE)=XNSALT(N)*1.0d6  ! #/m3
         AM(N,ITYPE) = RSALTGF(N)    ! m
         HYGRO_AER(N,ITYPE)=ZK(ITYPE)
      ENDDO

! RDST below to be updated later
      N=1
      IF((2.*RDST(N)).GE.ZDACT(3,1))
     &  ZCCN(3,1) = ZCCN(3,1) + XNDST(N)
      IF((2.*RDST(N)).GE.ZDACT(3,2))  ! dust CCN at S=~0.4%
     &  ZCCN(3,2) = ZCCN(3,2) + XNDST(N)
      IF((2.*RDST(N)).GE.ZDACT(3,3))
     &  ZCCN(3,3) = ZCCN(3,3) + XNDST(N)

      DO N=2,NDSTB  
         IF(ZDACT(3,2).GE.(2.*RDST(N-1)).and.
     &      ZDACT(3,2).LT.(2.*RDST(N))) RACT3=N
         IF((2.*RDST(N)).GE.ZDACT(3,1))
     &     ZCCN(3,1) = ZCCN(3,1) + XNDST(N)
         IF((2.*RDST(N)).GE.ZDACT(3,2))  ! dust CCN at S=~0.4%
     &     ZCCN(3,2) = ZCCN(3,2) + XNDST(N)
         IF((2.*RDST(N)).GE.ZDACT(3,3))
     &     ZCCN(3,3) = ZCCN(3,3) + XNDST(N)
      ENDDO

! for max sat
      ITYPE = 3
      NSIZE_AER(ITYPE)=NDSTB
      DO N= 1, NSIZE_AER(ITYPE)
         NAEROSOL(N,ITYPE)=XNDST(N)*1.0d6  ! #/m3
         AM(N,ITYPE) = RDST(N)    ! m
         HYGRO_AER(N,ITYPE)=ZK(ITYPE)
      ENDDO

      ITYPE = 4
      IF(ZTN(ITYPE).GT.1.) THEN
       TAREABC =SUM(YAREABC)
       DO N=1, NBCOC
! distribute MBCS according to surface area
         YBCS = MBCS * YAREABC(N)/TAREABC
         MBCSOA = MBCLV * YAREABC(N)/TAREABC
         MBCTOT = XMBC(N)
     &          + YBCS*GFGAS**3.0 ! sulfate and associated NIT, CH4, SOA
     &          + MBCSOA           ! SV-SOA not considered for now
         RBCGF(N) = RBCOC(N)*(MBCTOT/XMBC(N))**(1./3.)
       ENDDO
      ELSE
       RBCGF = RBCOC
      ENDIF

      DO N=1,NBCOC  ! BC CCN at 3 S
         IF((2.*RBCGF(N)).GE.ZDACT(ITYPE,1))
     &     ZCCN(ITYPE,1) = ZCCN(ITYPE,1) + XNBC(N)
         IF((2.*RBCGF(N)).GE.ZDACT(ITYPE,2))
     &     ZCCN(ITYPE,2) = ZCCN(ITYPE,2) + XNBC(N)
         IF((2.*RBCGF(N)).GE.ZDACT(ITYPE,3))
     &     ZCCN(ITYPE,3) = ZCCN(ITYPE,3) + XNBC(N)
      ENDDO

! for max sat
      NSIZE_AER(ITYPE)=NBCOC
      DO N= 1, NSIZE_AER(ITYPE)
         NAEROSOL(N,ITYPE)=XNBC(N)*1.0d6  ! #/m3
         AM(N,ITYPE) = RBCGF(N)    ! m
         HYGRO_AER(N,ITYPE)=ZK(ITYPE)
      ENDDO

      ITYPE = 5
      IF(ZTN(ITYPE).GT.1.) THEN
       TAREAOC =SUM(YAREAOC)
       DO N=1, NBCOC
! distribute MOCS according to surface area
         YOCS = MOCS * YAREAOC(N)/TAREAOC
         MOCSOA = MOCLV * YAREAOC(N)/TAREAOC
         MOCTOT = XMOC(N)
     &          + YOCS*GFGAS**3.0 ! sulfate and associated NIT, CH4, SOA
     &          + MOCSOA          ! SV-SOA not considered for now -- Need to consider
         ROCGF(N) = RBCOC(N)*(MOCTOT/XMOC(N))**(1./3.)
       ENDDO
      ELSE
       ROCGF = RBCOC
      ENDIF

      DO N=1,NBCOC  ! OC CCN at 3 S
         IF((2.*ROCGF(N)).GE.ZDACT(ITYPE,1))
     &     ZCCN(ITYPE,1) = ZCCN(ITYPE,1) + XNOC(N)
         IF((2.*ROCGF(N)).GE.ZDACT(ITYPE,2))
     &     ZCCN(ITYPE,2) = ZCCN(ITYPE,2) + XNOC(N)
         IF((2.*ROCGF(N)).GE.ZDACT(ITYPE,3))
     &     ZCCN(ITYPE,3) = ZCCN(ITYPE,3) + XNOC(N)
      ENDDO

! for max sat
      NSIZE_AER(ITYPE)=NBCOC
      DO N= 1, NSIZE_AER(ITYPE)
         NAEROSOL(N,ITYPE)=XNOC(N)*1.0d6  ! #/m3
         AM(N,ITYPE) = ROCGF(N)    ! m
         HYGRO_AER(N,ITYPE)=ZK(ITYPE)
      ENDDO

      ISITE = ISITES(42)
      JSITE = JSITES(42)
!      IF(II.EQ.ISITE.and.JJ.EQ.JSITE.AND.LL.LE.1) THEN
!        WRITE(1001,100)MSO4,MSULFLV,MMSAB,MNITB,MNH4B,MSOA,FBCOC1,FBCOC2
!        DO N =1, NTYP
!         WRITE(1001,100)MCORE(N),MSP(N),YCS(N),ZK(N),
!     &       ZDACT(N,1),ZDACT(N,2),ZDACT(N,3),
!     &       ZTN(N),ZCCN(N,1),ZCCN(N,2),ZCCN(N,3)
!        ENDDO
!        flush(1001)
!      ENDIF

      TCCN(2) = SUM(ZCCN(1:5,2))

      IF(TCCN(2)>1.d-9)THEN

      DO N=1,NSO4
         IF(N.LT.IACT20) THEN     ! average act diameter for sulfate aquous gain
            FCLOUD1(N) = 0.d0
         ELSE
            FCLOUD1(N) = XN(N)/TCCN(2)  ! parition based on number conc for now
         ENDIF
      ENDDO
      FCLOUD1(NSO4+1) = ZCCN(4,2)/TCCN(2)  ! partition based on number conc for now
      FCLOUD1(NSO4+2) = ZCCN(5,2)/TCCN(2)
      FCLOUD1(NSO4+3) = ZCCN(3,2)/TCCN(2)
      FCLOUD1(NSO4+4) = ZCCN(2,2)/TCCN(2)

      ELSE

      FCLOUD1=0.D0
      DO N=26,NSO4
        FCLOUD1(N)=1/15.D0
      ENDDO

      ENDIF

!******************************************************************************

!******************************************************************************
!OPT+
      IF(IFOPT.GT.0) THEN
! calculation of optical properties
      ITEST = 0
      IF(II.EQ.ISITES(27).and.JJ.EQ.JSITES(27).AND.LL.LE.1) THEN
         ITEST = 1
      ENDIF

!      WLS(1)=0.55  ! wavelength in um
      YBEXT=1.d-20
      XBEXT1k=0.d0
      YW=1.d-0
      YG=0.d-20

      ZBABS = 1.d-20
!SP   
      ITYPE = 1
      IF(YCS(ITYPE).GT.1.d-5.and.MSO4B.GT.1.d-12) THEN
       YVDRY = (MSO4B+MMSAB+MNITB+MNH4B)/DENSULF   ! INORGANICS
     &  + (MSULFLV+MSOA)/DENAER(5)               ! ORGANICS
      
       VH2O = YVDRY*(GFWATERVOL-1.0)
       VTOTAL = YVDRY + VH2O
     
       ZFV(1) = ((MSO4B+MMSAB)/DENSULF)/VTOTAL   !CORE, SO4
       ZFV(2) = ZFV(1)                   !SO4 
       ZFV(3) = (MNH4B/DENSULF)/VTOTAL   !NH4
       ZFV(4) = (MNITB/DENSULF)/VTOTAL   !NO3
       ZFV(5) = ((MSULFLV+MSOA)/DENAER(5))/VTOTAL   !SOA
       ZFV(6) = VH2O/VTOTAL              !H2O 

       IF(ITEST.EQ.1) WRITE(1001,*)"ZFV=",(ZFV(I),I=1,6)

       ZFV_SP = ZFV

       MBIN = 20
       DO I = 1, MBIN
        XDCORE(I)=RDRY(NSO4-MBIN+I)*2.E6   ! RDRY in m, XDCORE in um
        XN_IN(I) = XN(NSO4-MBIN+I)  
        IF(IFCOAT.EQ.0) THEN ! When coated SP on PP not considered in OPT calculation, put them in SP
          IF(FSO4B.LT.0.1) THEN
             VRATIO = 1./0.1
          ELSE
             VRATIO = 1./FSO4B
          ENDIF
          XN_IN(I) = XN_IN(I)*VRATIO
        ENDIF
       ENDDO

       IF(IFOPT==1)THEN
       CALL APM_OPT(ITEST,ITYPE,MBIN,XDCORE(1:MBIN),XN_IN(1:MBIN),
     &             ZFV,MWLS,WLS,XBEXT,XW,XG)
       YBEXT(ITYPE,:)=XBEXT(:)
       YW(ITYPE,:)=XW(:)
       YG(ITYPE,:)=XG(:)
       IF(ITEST.EQ.1) THEN
          WRITE(1001,191)ITYPE,sum(XN_IN(1:MBIN)),
     &        YBEXT(ITYPE,1),YW(ITYPE,1),YG(ITYPE,1)
       ENDIF
       ENDIF

      ENDIF

! sea salt
      ITYPE = 2
      IF(YCS(ITYPE).GT.1.d-5) THEN
       yvdry = SUM(XVSALT)
       vh2o = yvdry*(GF4**3.0-1.0)
       vtotal = yvdry + vh2o
 
       zfv(1) = yvdry/vtotal
       zfv(2) = 0.  !coated SP not considered for now
       zfv(3) = 0. 
       zfv(4) = 0.
       zfv(5) = 0.
       zfv(6) = vh2o/vtotal
 
       MBIN = 18
       DO I = 1, MBIN
        XDCORE(I)=RSALT(NSEA-MBIN+I)*2.E6   ! RSALT in m, XDCORE in um
        XN_IN(I) = XNSALT(NSEA-MBIN+I)  
       ENDDO
 
       IF(IFOPT==1)THEN
       CALL APM_OPT(ITEST,ITYPE,MBIN,XDCORE(1:MBIN),XN_IN(1:MBIN),
     &             ZFV,MWLS,WLS,XBEXT,XW,XG)
       YBEXT(ITYPE,:)=XBEXT(:)
       YW(ITYPE,:)=XW(:)
       YG(ITYPE,:)=XG(:)
       ENDIF

      ENDIF

     
!! dust
      ITYPE = 3
      IF(YCS(ITYPE).GT.1.d-5) THEN
       yvdry = SUM(XVDST)
       vh2o = 0.0
       vtotal = yvdry + vh2o
 
       zfv(1) = yvdry/vtotal
       zfv(2) = 0.  !coated SP not considered for now
       zfv(3) = 0.
       zfv(4) = 0.
       zfv(5) = 0.
       zfv(6) = vh2o/vtotal
 


       MBIN = NDSTB - 2
       DO I = 1, MBIN
        XDCORE(I)=RDST(NDSTB-MBIN+I)*2.E6   ! RSALT in m, XDCORE in um
        XN_IN(I) = XNDST(NDSTB-MBIN+I)
       ENDDO
       
       IF(IFOPT==1)THEN
       CALL APM_OPT(ITEST,ITYPE,MBIN,XDCORE(1:MBIN),XN_IN(1:MBIN),
     &             ZFV,MWLS,WLS,XBEXT,XW,XG)
       YBEXT(ITYPE,:)=XBEXT(:)
       YW(ITYPE,:)=XW(:)
       YG(ITYPE,:)=XG(:)

!LW  for dust only
       CALL APM_OPT_LW(ITEST,ITYPE,MBIN,XDCORE(1:MBIN),
     &   XN_IN(1:MBIN),MWLL,WLL,XBEXT1,XW1,XG1)

!       ZBEXT1(:)=XBEXT1(:)
!       ZW1(:)=XW1(:)
!       ZG1(:)=XG1(:)
       DO IWL = 1, MWLL
         ZBABS(IWL) = XBEXT1(IWL) * (1.-XW1(IWL))  ! LW calculation only need absorption coefficient
       ENDDO
       ENDIF
       
      ENDIF
 
!
!! BC 
      ITYPE=4
      IF(YCS(ITYPE).GT.1.D-5) THEN
       VRATIO = RRATIOWETBC**3.
       IF(IFCOAT.EQ.0.) VRATIO = 1.
!-----------------------------------------------------------------------------
!BMY KLUDGE:  Add this here to avoid floating point exception! (bmy, 6/28/19)
       zfv_sp = zfv
!-----------------------------------------------------------------------------

       zfv(1) = 1./VRATIO
       zfv(2) = ZFV_SP(2)*(1.-1./VRATIO)  !assume coated SP has same composition as SP for now
       zfv(3) = ZFV_SP(3)*(1.-1./VRATIO)
       zfv(4) = ZFV_SP(4)*(1.-1./VRATIO)
       zfv(5) = ZFV_SP(5)*(1.-1./VRATIO)
       zfv(6) = ZFV_SP(6)*(1.-1./VRATIO)

       MBIN = NBCOC
       DO I = 1, MBIN
        XDCORE(I)=RBCOC(I)*2.E6   ! RDRY in m, XDCORE in um
        XN_IN(I) = XNBC(I)
       ENDDO

       IF(IFOPT==1)THEN
       CALL APM_OPT(ITEST,ITYPE,MBIN,XDCORE(1:MBIN),XN_IN(1:MBIN),
     &             ZFV,MWLS,WLS,XBEXT,XW,XG)
       YBEXT(ITYPE,:)=XBEXT(:)
       YW(ITYPE,:)=XW(:)
       YG(ITYPE,:)=XG(:)
     
       IF(ITEST.EQ.1) THEN
          WRITE(1001,191)ITYPE,sum(XN_IN(1:MBIN)),
     &        YBEXT(ITYPE,1),YW(ITYPE,1),YG(ITYPE,1)

       ENDIF
       ENDIF

      ENDIF
 191  FORMAT(I3,10(1PE10.3))

!
!! POC  
      ITYPE=5
      IF(YCS(ITYPE).GT.1.D-5) THEN
       VRATIO = RRATIOWETOC**3.

       IF(IFCOAT.EQ.0.) VRATIO = 1.
!-----------------------------------------------------------------------------
!BMY KLUDGE: Add this here to avoid floating point exception! (bmy, 6/28/19)
       zfv_sp = zfv
!-----------------------------------------------------------------------------

       zfv(1) = 1./VRATIO
       zfv(2) = ZFV_SP(2)*(1.-1./VRATIO)  !assume coated SP has same composition as SP for now
       zfv(3) = ZFV_SP(3)*(1.-1./VRATIO)
       zfv(4) = ZFV_SP(4)*(1.-1./VRATIO)
       zfv(5) = ZFV_SP(5)*(1.-1./VRATIO)
       zfv(6) = ZFV_SP(6)*(1.-1./VRATIO)

       MBIN = NBCOC
       DO I = 1, MBIN
        XDCORE(I)=RBCOC(I)*2.E6   ! RBCOC in m, XDCORE in um
        XN_IN(I) = XNOC(I)
       ENDDO

       IF(IFOPT==1)THEN
       CALL APM_OPT(ITEST,ITYPE,MBIN,XDCORE(1:MBIN),XN_IN(1:MBIN),
     &             ZFV,MWLS,WLS,XBEXT,XW,XG)
       YBEXT(ITYPE,:)=XBEXT(:)
       YW(ITYPE,:)=XW(:)
       YG(ITYPE,:)=XG(:)
       ENDIF

      ENDIF

      IF(IFCDNPDF==1)THEN
! CDN calculation based on PDF

      SIGMA = 0.5
      DW = 4.0/float(MW)
      DO IW = 1, MW
       W(IW)=float(IW-MW/2-1)*DW + 0.05
       PW(IW)=exp(-0.5*W(IW)*W(IW)/(SIGMA*SIGMA))
     &    /(SQRT(2.*3.1416)*SIGMA)
      ENDDO

      YCDN = 0.
      DO IW = 1, MW
!       WW = W(IW) + WBAR
       WW = W(IW) + VZ
       ZCDN = 1.d-10
       IF(WW.GT.0.) THEN
        call activate(maxsatout,WW,tk,press, 
     &          NTYP, NTYP, NSO4, nsize_aer,    
     &          naerosol, am, hygro_aer)

        XSAT = maxsatout 
        CALL DACTS(TK,MSO4,MSULFLV,MMSAB,MNITB,MNH4B,MSOA,
     &              MCORE,MSP,XSAT,ZK,ZDACT)


! CDN at calculated S with assumed updarft velocity
        IS = 4
        IACT40=1
        RACT1=1
        RACT2=1
        RACT3=1

        DO N = 2, NSO4
         DIAM1GAS = GFGAS*RDRY(N-1)*2.
         DIAM2GAS = GFGAS*RDRY(N)*2.
         IF(ZDACT(1,IS).GE.DIAM1GAS.and.ZDACT(1,IS).LT.DIAM2GAS) 
     &             IACT40=N
        ENDDO
        ZCDN(1) = SUM(XN(IACT40:NSO4))
        RACT1 = IACT40

        N=1
        IF((2.*RSALTGF(N)).GE.ZDACT(2,IS))
     &    ZCDN(2) = XNSALT(N)
        IF((2.*RDST(N)).GE.ZDACT(3,IS))
     &    ZCDN(3) = XNDST(N)

        DO N=2,NSEA
         IF(ZDACT(2,IS).GE.(2.*RSALTGF(N-1)).and.
     &      ZDACT(2,IS).LT.(2.*RSALTGF(N))) RACT2=N
         IF((2.*RSALTGF(N)).GE.ZDACT(2,IS))THEN
           ZCDN(2) = ZCDN(2) + XNSALT(N)
         ENDIF
        ENDDO

        DO N=2,NDSTB  
         IF(ZDACT(3,IS).GE.(2.*RDST(N-1)).and.
     &      ZDACT(3,IS).LT.(2.*RDST(N))) RACT3=N  
         IF((2.*RDST(N)).GE.ZDACT(3,IS))THEN
           ZCDN(3) = ZCDN(3) + XNDST(N)
         ENDIF
        ENDDO

      DO N=1,NBCOC
         IF((2.*RBCGF(N)).GE.ZDACT(4,IS))
     &     ZCDN(4) = ZCDN(4) + XNBC(N)
         IF((2.*ROCGF(N)).GE.ZDACT(5,IS))
     &     ZCDN(5) = ZCDN(5) + XNOC(N)
      ENDDO



       ENDIF
       DO ITYPE = 1, NTYP
        YCDN(ITYPE)=YCDN(ITYPE) + PW(IW)*DW* ZCDN(ITYPE)
       ENDDO
      ENDDO

      XCDN = sum(YCDN)
      XCDNSP = YCDN(1)

      ELSE
        XCDN = 50.D6
        XCDNSP = 50.D6
      ENDIF

!      IF(ITEST.EQ.1) THEN
!        WRITE(1002,*)TK,PRESS,maxsatout
!        DO ITYPE = 1, NTYP
!          WRITE(1002,192)ITYPE,NSIZE_AER(ITYPE)
!          WRITE(1002,193)(NAEROSOL(N,ITYPE),N=1,NSIZE_AER(ITYPE))
!          WRITE(1002,193)(AM(N,ITYPE),N=1,NSIZE_AER(ITYPE))
!          WRITE(1002,193)(HYGRO_AER(N,ITYPE),N=1,NSIZE_AER(ITYPE))
!        ENDDO
!      ENDIF
! 192  FORMAT(I3,I3)
! 193  FORMAT(50(1PE9.2))

!
!! Calculate total optical properties
      ZBEXT = 1.d-20
      DO IWL=1, MWLS
         ZTEMP1 = 1.d-20
         ZTEMP2 = 1.d-20
         DO ITYPE=1,NTYP
          ZBEXT(IWL)=ZBEXT(IWL) + YBEXT(ITYPE,IWL)

          ZTEMP1 = ZTEMP1 + YBEXT(ITYPE,IWL)*YW(ITYPE,IWL)   ! total scattering
          ZTEMP2 = ZTEMP2 + YBEXT(ITYPE,IWL)*YW(ITYPE,IWL)*YG(ITYPE,IWL)   
         ENDDO
         ZW(IWL)=ZTEMP1/ZBEXT(IWL)
         ZG(IWL)=ZTEMP2/ZTEMP1

      ENDDO
!
!******************************************************************************
!OPT+
      ENDIF  !IFOPT
!******************************************************************************
! Reduced time step for nucleation/growth

      IFNGC = 0   ! timestep not reduced
      DNMAX = 500.  ! max change of XN1 due to nucl

      CACID0 = CACID

      CACIDAVE = 1.D-30
      YJAVE = 1.D-30
      YYJAVE = 1.D-30
      YJATHNAVE = 1.D-30
      YJATHN = 1.D-30
      YJAVE_BH = 1.D-30
      YJAVE_BIM = 1.D-30
      YJAVE_TH = 1.D-30
      YJAVE_TIM = 1.D-30
      YJAVE_OM = 1.D-30
      YJAVE_PON = 1.D-30
      YJAVE_POI = 1.D-30

      YJAVE_BN = 1.D-30
      YJAVE_BI = 1.D-30
      YJAVE_TN = 1.D-30
      YJAVE_TI = 1.D-30

      DACID = 0.
      DLVSOG = 0.

      DTLEFT = DT

      IF(PLVSOG1.GT.1.0) THEN
        CHOM = CLVSOG*PLVSOG01/PLVSOG1
      ELSE
        CHOM = CLVSOG
      ENDIF

      DTNGC = DT
!      IF(CACID.GT.5.d5)THEN
      IF(CACID.GT.5.d5.or.(IFNUCL.GE.13.and.CHOM.GT.5.d5))THEN
 80      CONTINUE
! Nucl   
         IF(IFNUCL.EQ.1) THEN   ! IMN
            CALL YUJTIMN1(CACID,RH,TK,XQ,XS,XNH3,XJBH,XJBIM,XJTH,XJTIM,
     &                    RBH,RBIM,RTH,RTIM)
            YJ = XJTIM


          YYJ = YJ
          YJ = YYJ
          
         ELSEIF(IFNUCL.EQ.2) THEN   ! KBHN (XQ=1.E-20, XNH3=0)
            XQ=1.E-20
            XNH3=1.E-20
            CALL YUJTIMN1(CACID,RH,TK,XQ,XS,XNH3,XJBH,XJBIM,XJTH,XJTIM,
     &                    RBH,RBIM,RTH,RTIM)
            YYJ = YJ
         ELSEIF(IFNUCL.EQ.3) THEN  !empirical activation nucl
            YJ = 3.5E-7 * CACID   !
         ELSEIF(IFNUCL.EQ.4) THEN  !empirical kinetic nucl
            YJ = 5.5E-14 * CACID * CACID   !
         ELSE
            write(6,*)"STOP at apm_phys_mod.f: Check IFNUCL value"
            stop
         ENDIF

         IF(DT.GT.(DNMAX/YJ)) THEN
            DTNGC = DNMAX/YJ    ! small timestep for large J
            DTNGC = MAX(DTNGC,1.8d2)  !set minimum DTNGC
         ELSE
            DTNGC = DT
         ENDIF

         IF(DTLEFT.GT.(1.5*DTNGC)) THEN
           DTLEFT = DTLEFT - DTNGC
         ELSE
           DTNGC = DTLEFT
           DTLEFT = 0.
         ENDIF
         IF(DTNGC.LT.DT)  IFNGC = 1    ! cut time step

         ZN1 =SUM(XN)

         YJAVE = YJAVE + YJ*DTNGC/DT

         YJATHNAVE = YJATHNAVE + YJATHN*DTNGC/DT

         YJAVE_BH = YJAVE_BH + XJBH*DTNGC/DT
         YJAVE_BIM = YJAVE_BIM + XJBIM*DTNGC/DT
         YJAVE_TH = YJAVE_TH + XJTH*DTNGC/DT
         YJAVE_TIM = YJAVE_TIM + XJTIM*DTNGC/DT

         !YJAVE_OM = YJAVE_OM + YJOM*DTNGC/DT
         !YJAVE_PON = YJAVE_PON + YJPON*DTNGC/DT
         !YJAVE_POI = YJAVE_POI + YJPOI*DTNGC/DT

         !YJAVE_BN = YJAVE_BN + XJBN*DTNGC/DT
         !YJAVE_BI = YJAVE_BI + XJBI*DTNGC/DT
         !YJAVE_TN = YJAVE_TN + XJTN*DTNGC/DT
         !YJAVE_TI = YJAVE_TI + XJTI*DTNGC/DT

         XN(1) = XN(1) + YJ*DTNGC
         DVJ = YJ * DTNGC * VDRY(1)*1.E6   ! cm3/cm3
         XVA(1) = XVA(1) + DVJ   ! cm3/cm3
         CACID = CACID - DVJ/V1ACID

! Calculate condensational growth of H2SO4
         IF(CACID.GT.5.d5) THEN
          ICOND = 1
          XMCONDIN = 96.0   !g/mol
          CACIDIN = CACID
          PACIDIN = PACID
          TCSOTHERIN = TCSOTHER
 
          SUMXVA0 = SUM(XVA)
 
          CALL APM_GROW(ICOND,NSO4,TK,PRESS,CACIDIN,PACIDIN, 
     &      TCSOTHERIN,DTNGC,XN,XVA,RWETCM,TOTCONDOTH,XMCONDIN)   ! XN in #/cm3, XVA in cm3/cm3
                                                 ! TOTCONDOTH in # H2SO4/ cm3

          SUMXVA1 = SUM(XVA)
          SUMXVA1 = MAX(SUM(XVA),1.D-30)
          DACID = SUMXVA1 - SUMXVA0   ! amount condensed (cm3/cm3)
          CACID = CACIDIN
! SULFATE MASS condensing on particles other than sulfate
          IF(TCSOTHER.GT.1.D-5) THEN
           MCONDOTH = TOTCONDOTH*V1ACID*DENSULF*1.d3  ! kg/m3
           MSALTS = MSALTS + MCONDOTH * YCS(2)/TCSOTHER
           MDSTS = MDSTS + MCONDOTH * YCS(3)/TCSOTHER
           MBCS   = MBCS +   MCONDOTH * YCS(4)/TCSOTHER
           MOCS   = MOCS +   MCONDOTH * YCS(5)/TCSOTHER
          ENDIF
          CACIDAVE = CACIDAVE + CACID*DTNGC/DT

          CALL APM_MOVEBIN(NSO4,XN,XVA)
         ELSE
          CACIDAVE = CACIDAVE + CACID*DTNGC/DT
         ENDIF

         IF(IFNGC.EQ.1) THEN
! coag
          ZN2 =SUM(XN)
          DTCOAG1 = DTNGC
          ITYPE = 1 !sulfate
          CALL APM_COAG(ITYPE,NSO4,TK,PMB,DTCOAG1,RWETCM,DENWET,XN,XVA)
          NCOAG1 = 0  ! reset counter to 0 after coag is called

! Scavenging of sulfate particles by other types of aerosols
          IF((TCS-CSSULF).GT.1.E-5) THEN
           RLOSULF = MSULFLV/(SUMXVA1*DENSULF*1.d3)   !ratio of LO mass on SULF to total mass
           RLOSULF = MIN(RLOSULF, 1.0d0)

           CALL APM_COAGSCAV(
     &     TK,PMB,DTNGC,DENWET,RWETCM,YCS,RLOSULF,
     &     DENSALTWET,RSALTWETCM,XNSALT,
     &     DENBCWET,RBCWET,XNBC,
     &     DENOCWET,ROCWET,XNOC,
     &     DENAER,XNDST,
     &     XVA,MBCS,MOCS,MDSTS,MSALTS,MBCLV,MOCLV,MDSTLV,MSALTLV)
          ENDIF

!APM2+IFAG++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          IF(IFAG.EQ.1) THEN
! LVSOG condensation on secondary particles, when condensed, LVSOG becomes part of SULF (for now)
! Calculate condensational growth of LVSOG
           IF(CLVSOG.GT.CLVSOG0) THEN
            ICOND = 2
            XMCONDIN = XMLVSOG   !g/mol
            CACIDIN = CLVSOG
            PACIDIN = PLVSOG1  !#/cm3s
            TCSOTHERIN = TCSOTHER*SQRT(96./XMCONDIN)  !scale the CS
            SUMXVA0 = SUM(XVA)
            TOTCONDOTH = 0.

            CALL APM_GROW(ICOND,NSO4,TK,PRESS,CACIDIN,PACIDIN, 
     &      TCSOTHERIN,DTNGC,XN,XVA,RWETCM,TOTCONDOTH,XMCONDIN)   ! XN in #/cm3, XVA in cm3/cm3
                                                 ! TOTCONDOTH in # acid/cm3
            CALL APM_MOVEBIN(NSO4,XN,XVA)
! LVSOG MASS condensing on particles other than sulfate
            IF(TCSOTHERIN.GT.1.D-5) THEN
             MCONDOTH = TOTCONDOTH*V1LVSOG*DENSULF*1.d3  ! kg/m3
             MSALTLV = MSALTLV + MCONDOTH * YCS(2)/TCSOTHER
             MDSTLV = MDSTLV + MCONDOTH * YCS(3)/TCSOTHER
             MBCLV = MBCLV + MCONDOTH * YCS(4)/TCSOTHER
             MOCLV = MOCLV + MCONDOTH * YCS(5)/TCSOTHER
            ENDIF
            SUMXVA1 = MAX(SUM(XVA),1.D-30)
            DLVSOG = SUMXVA1 - SUMXVA0   ! amount of LVSOG condensed (cm3/cm3)
            MSULFLV = MSULFLV + DLVSOG*DENSULF*1.d3  ! kg/m3
            CLVSOG = CACIDIN
           ELSE
            IF(TCS.GT.1.d-6) THEN
! Allow the scaveging of residue LVSOG by particles
             TCSLV = TCS*SQRT(96./XMLVSOG)  !scale the CS
             TCSOTHERIN = TCSOTHER*SQRT(96./XMLVSOG)  !scale the CS
             CLVSOG1 = PLVSOG1/TCSLV+(CLVSOG - PLVSOG1/TCSLV) 
     &             *exp(-TCSLV*DTNGC)
            DLVSOG = CLVSOG + PLVSOG1*DTNGC - CLVSOG1
! Distribue all condensed LVSOG to primary particles for now
             MCONDOTH = DLVSOG*V1LVSOG*DENSULF*1.d3  ! kg/m3

             MSULFLV = MSULFLV + MCONDOTH * YCS(1)/TCS 
             XVA(NSO4-1) = XVA(NSO4-1)+MCONDOTH/
     &                   (DENSULF*1.d3) * YCS(1)/TCS  ! put in last bin of SP for now
             MSALTLV = MSALTLV + MCONDOTH * YCS(2)/TCS
             MDSTLV = MDSTLV + MCONDOTH * YCS(3)/TCS
             MBCLV = MBCLV + MCONDOTH * YCS(4)/TCS
             MOCLV = MOCLV + MCONDOTH * YCS(5)/TCS
            ELSE
             CLVSOG1 = CLVSOG + PLVSOG1*DTNGC
            ENDIF
            CLVSOG = CLVSOG1
           ENDIF
           RLOSULF = MSULFLV/(SUMXVA1*DENSULF*1.d3)   !ratio of LO mass on SULF to total mass
           RLOSULF = MIN(RLOSULF, 1.0d0)
          ENDIF
!APM2+ENDIF IFAG+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

         ENDIF

         IF(DTLEFT.GT.1.d-5) GOTO 80
      ELSE
!      ELSEIF(PACID.GT.1.d-3.and.TCS.GT.1.d-6) THEN
         !Luo20121201IF((PACID*DT).LT.5.d5) THEN
         IF((PACID*DT).LT.5.d5.or.TCS.LT.1.d-6) THEN
          CACID = CACID + PACID*DT
         ELSE
          CACID1 = PACID/TCS + (CACID - PACID/TCS)*exp(-TCS*DT)
          DACID = CACID + PACID*DT - CACID1
          if(DACID<0.d0)then
            write(*,*)'Luo warning DACID =',DACID,PACID,TCS,CACID,CACID1
            DACID=0.d0
            CACID1=CACID + PACID*DT
          endif

          MCONDOTH = DACID*V1ACID*DENSULF*1.d3  ! kg/m3

          XVA(NSO4-1) = XVA(NSO4-1)+MCONDOTH/
     &                (DENSULF*1.d3) * YCS(1)/TCS  ! put in NSO4-1 bin of SP for now
          MSALTS = MSALTS + MCONDOTH * YCS(2)/TCS
          MDSTS  = MDSTS + MCONDOTH * YCS(3)/TCS
          MBCS   = MBCS +   MCONDOTH * YCS(4)/TCS
          MOCS   = MOCS +   MCONDOTH * YCS(5)/TCS

          CACID = CACID1
         ENDIF
         CACIDAVE = CACID

         YJ = 1.E-20
      ENDIF

!APM2+IFAG++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      IF(IFAG.EQ.1) THEN
       IF(IFNGC.EQ.0) THEN  !Only when IFNGC=0 that LVSOG condensation has not yet been calculated above
! LVSOG condensation on secondary particles, when condensed, LVSOG becomes part of SULF
! Calculate condensational growth of LVSOG
        IF(CLVSOG.GT.CLVSOG0) THEN
         ICOND = 2
         XMCONDIN = XMLVSOG   !g/mol
!           WRITE(6,*)"V1LVSOG = ",V1LVSOG
!           V1LVSOG = XMLVSOG/(6.02d23*DENSULF)   ! assume SOA has same den as acid
!           WRITE(6,*)"V1LVSOG = ",V1LVSOG
         CACIDIN = CLVSOG
         PACIDIN = PLVSOG1  !#/cm3s
         TCSOTHERIN = TCSOTHER*SQRT(96./XMCONDIN)  !scale the CS
         SUMXVA0 = SUM(XVA)
         TOTCONDOTH = 0.
!
         CALL APM_GROW(ICOND,NSO4,TK,PRESS,CACIDIN,PACIDIN, 
     &      TCSOTHERIN,DT,XN,XVA,RWETCM,TOTCONDOTH,XMCONDIN)   ! XN in #/cm3, XVA in cm3/cm3
                                                 ! TOTCONDOTH in # acid/cm3

! LVSOG MASS condensing on particles other than sulfate
         IF(TCSOTHER.GT.1.D-5) THEN
          MCONDOTH = TOTCONDOTH*V1LVSOG*DENSULF*1.d3  ! kg/m3
          MSALTLV = MSALTLV + MCONDOTH * YCS(2)/TCSOTHER
          MDSTLV = MDSTLV + MCONDOTH * YCS(3)/TCSOTHER
          MBCLV = MBCLV + MCONDOTH * YCS(4)/TCSOTHER
          MOCLV = MOCLV + MCONDOTH * YCS(5)/TCSOTHER
         ENDIF
         SUMXVA1 = MAX(SUM(XVA),1.D-30)
         DLVSOG = SUMXVA1 - SUMXVA0   ! amount of LVSOG condensed (cm3/cm3)
         MSULFLV = MSULFLV + DLVSOG*DENSULF*1.d3  ! kg/m3
         CLVSOG = CACIDIN
        ELSE
! Allow the scaveging of residue LVSOG by particles
            IF(TCS.GT.1.d-6) THEN
! Allow the scaveging of residue LVSOG by particles
             TCSLV = TCS*SQRT(96./XMLVSOG)  !scale the CS
             TCSOTHERIN = TCSOTHER*SQRT(96./XMLVSOG)  !scale the CS
             CLVSOG1 = PLVSOG1/TCSLV+(CLVSOG - PLVSOG1/TCSLV) 
     &             *exp(-TCSLV*DTNGC)
             DLVSOG = CLVSOG + PLVSOG1*DTNGC - CLVSOG1
! Distribue all condensed LVSOG to primary particles for now
             MCONDOTH = DLVSOG*V1LVSOG*DENSULF*1.d3  ! kg/m3

             MSULFLV = MSULFLV + MCONDOTH * YCS(1)/TCS 
             XVA(NSO4-1) = XVA(NSO4-1)+MCONDOTH/
     &                   (DENSULF*1.d3) * YCS(1)/TCS  ! put in NSO4-1 bin of SP for now
             MSALTLV = MSALTLV + MCONDOTH * YCS(2)/TCS
             MDSTLV = MDSTLV + MCONDOTH * YCS(3)/TCS
             MBCLV = MBCLV + MCONDOTH * YCS(4)/TCS
             MOCLV = MOCLV + MCONDOTH * YCS(5)/TCS
            ELSE
             CLVSOG1 = CLVSOG + PLVSOG1*DTNGC
            ENDIF
            CLVSOG = CLVSOG1
        ENDIF

        SUMXVA1 = MAX(SUM(XVA),1.D-30)
        RLOSULF = MSULFLV/(SUMXVA1*DENSULF*1.d3)   !ratio of LO mass on SULF to total mass
        RLOSULF = MIN(RLOSULF, 1.0d0)

! Move particles across bins after growth
!         YN0 = sum(XN)
!         YV0 = sum(XVA)
        CALL APM_MOVEBIN(NSO4,XN,XVA)
!         YN1 = sum(XN)
!         YV1 = sum(XVA)
!         IF(abs(YN1-YN0).GT.(0.1*YN0).or.abs(YV1-YV0).GT.(0.1*YV0)) THEN
!            WRITE(7001,100)YN0,YN1,YV0,YV1
!            flush(7001)
!         ENDIF

       ENDIF
      ENDIF
!APM2+ENDIF IFAG+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


      IF(CACID.LT.0.0) THEN
         WRITE(6,*) "CACID < 0, set it to 1E2/cc"
        WRITE(6,99)II,JJ,LL,CACID,CACID0,PACID,YJAVE,CSSULF, 
     &   GFGAS,MMSA,MNIT,MNH4,MSOA,MSO4, 
     &   TCSOTHER,TOTCONDOTH,MCONDOTH,MDSTS
         CACID = 1.E2
      ENDIF

!******************************************************************************
! Coagulation  
      EPC = 0.1      ! ratio of coag timestep to coag half lifetime
      COAGN0 = 10.   ! only consider coag if substatial # of particle exits
!      NCOAGMAX = 240
      NCOAGMAX = 120

! Seasalt
      IF(ZTN(2).GT.COAGN0) THEN  ! only consider coag if substatial # of particle exits
        NCOAG(2) = MAX(1,INT(EPC/(0.5d-8*ZTN(2)*DT)))
        IF((NCOAG2+1).GE.NCOAG(2)
     &               .or.(NCOAG2+1).GE.NCOAGMAX) THEN   
         DTCOAG2 = float(NCOAG2+1)*DT
         ITYPE = 2 !seasalt
          CALL APM_COAG(ITYPE,NSEA,TK,PMB,DTCOAG2,
     &            RSALTWETCM,DENSALTWET,XNSALT,XVSALT)
          NCOAG2 = 0  ! reset counter to 0 after coag is called

          DO N = 1, NSEA
            XMSALT(N)=XVSALT(N)*DENAER(2)*1.d3  !XVA in cm3/cm3  XMA in kg/m3
            XM1D(NSO4+N)=XMSALT(N)
          ENDDO

        ELSE
         NCOAG2=NCOAG2 + 1   ! count step that coag is not call in the grid
        ENDIF
      ENDIF
!
! BCOC
      IF((ZTN(4)+ZTN(5)).GT.COAGN0) THEN  ! only consider coag if substatial # of particle exits
        NCOAG(4) = MAX(1,INT(EPC/(0.5D-8*(ZTN(4)+ZTN(5))*DT)))
        IF((NCOAG4+1).GE.NCOAG(4).or.(NCOAG4+1).GE.NCOAGMAX) THEN
         DTCOAG4 = float(NCOAG4+1)*DT
         ITYPE = 4 !BC
          CALL APM_COAG(ITYPE,NBCOC,TK,PMB,DTCOAG4,
     &                  RBCWET,DENBCWET,XNBC,XVBC)

         ITYPE = 5 !OC
          CALL APM_COAG(ITYPE,NBCOC,TK,PMB,DTCOAG4,
     &                  ROCWET,DENOCWET,XNOC,XVOC)

          NCOAG4 = 0  ! reset counter to 0 after coag is called
          DO N = 1, NBCOC
            XMBC(N)=XVBC(N)*DENAER(4)*1.d3  !XVA in cm3/cm3  XMA in kg/m3
            XMOC(N)=XVOC(N)*DENAER(5)*1.d3  !XVA in cm3/cm3  XMA in kg/m3
          ENDDO
        ELSE
         NCOAG4=NCOAG4 + 1   ! count step that coag is not call in the grid
        ENDIF
      ENDIF
! Sulfate
      TOTN = SUM(XN)
      IF(TOTN.GT.COAGN0.and.IFNGC.EQ.0) THEN  ! only consider coag if substatial # of particle exits
        NCOAG(1) = MAX(1,INT(EPC/(0.5d-8*TOTN*DT)))
        IF((NCOAG1+1).GE.NCOAG(1)
     &               .or.(NCOAG1+1).GE.NCOAGMAX) THEN   
         DTCOAG1 = float(NCOAG1+1)*DT

         ITYPE = 1 !sulfate
         CALL APM_COAG(ITYPE,NSO4,TK,PMB,DTCOAG1,RWETCM,DENWET,XN,XVA)
         NCOAG1 = 0  ! reset counter to 0 after coag is called
        ELSE
         NCOAG1=NCOAG1 + 1   ! count step that coag is not call in the grid
        ENDIF

! Scavenging of sulfate particles by other types of aerosols

        IF((TCS-CSSULF).GT.1.E-5) THEN 
         IF(IFAG.EQ.0) RLOSULF = 0. 
         CALL APM_COAGSCAV(
     &   TK,PMB,DT,DENWET,RWETCM,YCS,RLOSULF,
     &   DENSALTWET,RSALTWETCM,XNSALT,
     &   DENBCWET,RBCWET,XNBC,
     &   DENOCWET,ROCWET,XNOC,
     &   DENAER,XNDST,
     &   XVA,MBCS,MOCS,MDSTS,MSALTS,MBCLV,MOCLV,MDSTLV,MSALTLV)
        ENDIF

      ENDIF

! Here convert volume back to mass
      DO N = 1, NSO4
        XMA(N)=XVA(N)*DENSULF*1.d3  !XVA in cm3/cm3  XMA in kg/m3
        XM1D(N)=XMA(N)
! Update numb conc
        XN(N) =  XVA(N)/(1.E6*VDRY(N))   ! XVA in cm3/cm3, VDRY in m3,XN in #/cm3
      ENDDO

 99   FORMAT(I4,I4,I4,50(1PE9.2))
100   FORMAT(50(1PE9.2))

! tempout  
      IF(NTEMPOUT1.EQ.1) THEN
       TN10nm = SUM(XN(I10nm:NSO4))
       TNOTHER = SUM(ZTN(2:5))
       TEMPOUT1(1)=TN10nm + TNOTHER  ! total CN10 (#/cm3)
!APM2++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ELSEIF(NTEMPOUT1.GT.10) THEN
       ZCCN(1,1)=SUM(XN(IACT10:NSO4))
       ZCCN(1,2)=SUM(XN(IACT20:NSO4))
       ZCCN(1,3)=SUM(XN(IACT30:NSO4))
       ZTN(1) = SUM(XN(I3nm:NSO4))
       TN10nm = SUM(XN(I10nm:NSO4))

       TNOTHER = SUM(ZTN(2:5))
       TOTN = ZTN(1) + TNOTHER

       TCCN(1)=SUM(ZCCN(1:5,1))          !CCN0.8
       TCCN(2)=SUM(ZCCN(1:5,2))          !CCN0.4
       TCCN(3)=SUM(ZCCN(1:5,3))          !CCN0.2

       MWATERSP = SUM(MSP(1:5))*(1.0/1.8)*(GFWATERVOL-1.) ! total aerosol water (kg/m3)

       TEMPOUT1(1)=TCS
!       TEMPOUT1(2)=PACID
!       TEMPOUT1(3)=PLVSOG1  !LVSOG production rate
!       TEMPOUT1(3)=YYJAVE
!       TEMPOUT1(4)=XQ
        TEMPOUT1(2)=TK
        TEMPOUT1(3)=CACIDAVE
        TEMPOUT1(4)=CNH3
!       IF(IFNUCL.EQ.1.or.IFNUCL.EQ.2.or.IFNUCL.EQ.17) THEN
!        TEMPOUT1(3)=XNH3
!       ELSEIF(IFNUCL.GE.13) THEN
!        TEMPOUT1(3)=CHOM
!       ELSE
!        TEMPOUT1(3)=PLVSOG1  !LVSOG production rate
!       ENDIF
!       TEMPOUT1(4)=YYJAVE
       TEMPOUT1(5)=YJAVE

       TEMPOUT1(6)=TOTN  ! total >3nm particle # conc (#/cm3)
     &            * (1.013E5*TK)/(PRESS*273.15) ! normalize to standdard condition

       ATOM4N(:)=ATOM4N(:)
     &          * (1.013E5*TK)/(PRESS*273.15) ! normalize to standdard condition

       TEMPOUT1(7)=TN10nm  ! >~10 nm sulfate particle # conc (#/cm3)
       TEMPOUT1(8)=TNOTHER  ! other particle # conc (#/cm3)

       TEMPOUT1(9)=CHOM
!       TEMPOUT1(9)=ZCCN(1,1)  ! sulfate particle CCN0.8 # conc (#/cm3)
!       TEMPOUT1(9)=TCCN(2)  ! CCN0.4 # conc (#/cm3)
       TEMPOUT1(10)=TCCN(1)  ! total CCN0.8 (#/cm3)

       TEMPOUT1(11)=ZCCN(1,2) !  sulfate particle CCN0.4 # conc (#/cm3)
       TEMPOUT1(12)=ZCCN(2,2)  ! total seasalt CCN0.4 # (#/cm3)
       TEMPOUT1(13)=ZCCN(3,2) ! total dust CCN0.4# (#/cm3)
       TEMPOUT1(14)=ZCCN(4,2)  ! total BC CCN0.4 # (#/cm3)
       TEMPOUT1(15)=ZCCN(5,2)  ! total OC CCN0.4 # (#/cm3)

!       TEMPOUT1(11)=ZCDN(1) !  sulfate particle CCN0.4 # conc (#/cm3)
!       TEMPOUT1(12)=ZCDN(2)  ! total seasalt CCN0.4 # (#/cm3)
!       TEMPOUT1(13)=ZCDN(3) ! total dust CCN0.4# (#/cm3)
!       TEMPOUT1(14)=ZCDN(4)  ! total BC CCN0.4 # (#/cm3)
!       TEMPOUT1(15)=ZCDN(5)  ! total OC CCN0.4 # (#/cm3)

       TEMPOUT1(16)=ZCCN(1,3)  ! sulfate CCN0.2 (#/cm3)
       TEMPOUT1(17)=TCCN(3)    !total CCN0.2 

         TEMPOUT1(18)=ZTN(1)  ! >~3 nm sulfate particle # conc (#/cm3)
!      TEMPOUT1(19)=ZTN(2)  ! >total salt particle # conc (#/cm3)=TNOTHER - ZTN(1)-SUM(ZTN(3:5))
         TEMPOUT1(19)=ZTN(3)  ! >total dust particle # conc (#/cm3)
         TEMPOUT1(20)=ZTN(4)  ! >total BC particle # conc (#/cm3)
         TEMPOUT1(21)=ZTN(5)  ! >total OC particle # conc (#/cm3)

!         IF(IFOPT.EQ.1) THEN
!           TEMPOUT1(22)= ZBEXT(KWL1)  ! total ext coef    -- 550 nm
           TEMPOUT1(22)= sum(SURF(1:5))
!           TEMPOUT1(23)= ZW(KWL1)
!           TEMPOUT1(24)= ZG(KWL1)
!           TEMPOUT1(23)= YBEXT(1,KWL1) ! SP ext @
           TEMPOUT1(23)= RGM  !in m
!           TEMPOUT1(24)= SIG
           TEMPOUT1(24)= DUST500
!           TEMPOUT1(25)= YBEXT(1,KWL1) ! SP ext @
!           TEMPOUT1(26)= YBEXT(2,KWL1) ! salt ext @ 
!           TEMPOUT1(27)= YBEXT(3,KWL1) ! dust ext @
!!           TEMPOUT1(28)= YW(3,KWL1)
!           TEMPOUT1(28)= (1.-YW(3,KWL1))*YBEXT(3,KWL1)
!           TEMPOUT1(29)= YBEXT(4,KWL1) ! BC ext @ 
!!           TEMPOUT1(30)= YW(4,KWL1)
!           TEMPOUT1(30)= (1.-YW(4,KWL1))*YBEXT(4,KWL1)
!           TEMPOUT1(31)= YBEXT(5,KWL1) ! POC ext 
!         ENDIF
         TEMPOUT1(25)=YJATHNAVE 
         TEMPOUT1(26)=YJAVE_BH
         TEMPOUT1(27)=YJAVE_BIM 
         TEMPOUT1(28)=YJAVE_TH 
         TEMPOUT1(29)=YJAVE_OM 
         TEMPOUT1(30)=YJAVE_PON 
         TEMPOUT1(31)=YJAVE_POI 

         DO N =1, 5
          TEMPOUT1(31+N)=MSP(N)    !total SP species in each type
!          TEMPOUT1(31+N)=SURF(N)    !total SP species in each type
          TEMPOUT1(36+N)=MCORE(N)  !for SP, SO4 is the core
!          TEMPOUT1(41+N)=SURF(N)    !total SP species in each type
!          TEMPOUT1(41+N)=ZK(N)  !kapa
!          TEMPOUT1(41+N)=YBEXT(N,KWL1)   !
         ENDDO
         TEMPOUT1(42)=sum(YBEXT(1:5,KWL1))   !
         TEMPOUT1(43)=YJAVE_BN 
         TEMPOUT1(44)=YJAVE_BI 
         TEMPOUT1(45)=YJAVE_TN 
         TEMPOUT1(46)=YJAVE_TI 
         TEMPOUT1(47)=YJAVE_TIM 
         TEMPOUT1(48)=XQ 

! Parameterize BCOC e-folding time based on coat/core ratio
!         XBCLIFE= D0CONV/(DRDTBC*24.)     !days
!         XBCLIFE= MIN(1.0d1,XBCLIFE)  !max e-folding time 10 days
!         XOCLIFE= D0CONV/(DRDTOC*24.)     !days
!         XOCLIFE= MIN(1.0d1,XOCLIFE)  !max e-folding time 10 days

!         TEMPOUT1(49)= DRDTBC
!         TEMPOUT1(50)= DRDTOC

         IF(IFAEROCOMOUT==1)THEN

           AEROCOMOUT =0.d0

           TEMPOUT1(51:53)=0.d0
           I50nm=100
           I80nm=100
           I120nm=100

           DO N=1,NSO4
             IF(GFGAS*RDRY(N)*2.d0.GE.3.d-9) THEN
               TEMPOUT1(51)=TEMPOUT1(51)+XN(N)
             ENDIF
             IF(GFGAS*RDRY(N)*2.d0.GE.50.d-9) THEN
               TEMPOUT1(52)=TEMPOUT1(52)+XN(N)
               I50nm(1)=MIN(N,I50nm(1))
             ENDIF
             IF(GFGAS*RDRY(N)*2.d0.GE.120.d-9) THEN
               TEMPOUT1(53)=TEMPOUT1(53)+XN(N)
               I120nm(1)=MIN(N,I120nm(1))
             ENDIF
             IF(RRATIOBC*RDRY(N)*2.d0.GE.3.d-9) THEN
               TEMPOUT1(51)=TEMPOUT1(51)+XNBC(N)
             ENDIF
             IF(RRATIOBC*RDRY(N)*2.d0.GE.50.d-9) THEN
               TEMPOUT1(52)=TEMPOUT1(52)+XNBC(N)
               I50nm(4)=MIN(N,I50nm(4))
             ENDIF
             IF(RRATIOBC*RDRY(N)*2.d0.GE.120.d-9) THEN
               TEMPOUT1(53)=TEMPOUT1(53)+XNBC(N)
               I120nm(4)=MIN(N,I120nm(4))
             ENDIF
             IF(RRATIOOC*RDRY(N)*2.d0.GE.3.d-9) THEN
               TEMPOUT1(51)=TEMPOUT1(51)+XNOC(N)
             ENDIF
             IF(RRATIOOC*RDRY(N)*2.d0.GE.50.d-9) THEN
               TEMPOUT1(52)=TEMPOUT1(52)+XNOC(N)
               I50nm(5)=MIN(N,I50nm(5))
             ENDIF
             IF(RRATIOOC*RDRY(N)*2.d0.GE.120.d-9) THEN
               TEMPOUT1(53)=TEMPOUT1(53)+XNOC(N)
               I120nm(5)=MIN(N,I120nm(5))
             ENDIF
           ENDDO
           DO N=1,NSEA
             IF(RSALT(N)*2.d0.GE.3.d-9) THEN
               TEMPOUT1(51)=TEMPOUT1(51)+XNSALT(N)
             ENDIF
             IF(RSALT(N)*2.d0.GE.50.d-9) THEN
               TEMPOUT1(52)=TEMPOUT1(52)+XNSALT(N)
               I50nm(2)=MIN(N,I50nm(2))
             ENDIF
             IF(RSALT(N)*2.d0.GE.120.d-9) THEN
               TEMPOUT1(53)=TEMPOUT1(53)+XNSALT(N)
               I120nm(2)=MIN(N,I120nm(2))
             ENDIF
           ENDDO
           DO N=1,NDSTB
             IF(RDST(N)*2.d0.GE.3.d-9) THEN
               TEMPOUT1(51)=TEMPOUT1(51)+XNDST(N)
             ENDIF
             IF(RDST(N)*2.d0.GE.50.d-9) THEN
               TEMPOUT1(52)=TEMPOUT1(52)+XNDST(N)
               I50nm(3)=MIN(N,I50nm(3))
             ENDIF
             IF(RDST(N)*2.d0.GE.120.d-9) THEN
               TEMPOUT1(53)=TEMPOUT1(53)+XNDST(N)
               I120nm(3)=MIN(N,I120nm(3))
             ENDIF
           ENDDO

           TEMPOUT1(54)=MSULFT*32.d9/96.d0
           TEMPOUT1(55)=MCORE(4)*1.d9

!Luo&Yu Based on observation (Jimenez et al., 2009), O:C ratio in LV-SOA
!is 4:5 and in SV-SOA 2:5.
!Luo&Yu So use 5/9 as carbon fraction for LV-SOA and 5/7 for other SOA
!(with molecular weight of 150 g/mol).
!Luo&Yu Ignore H in the molecules for now.
           TEMPOUT1(56)=MCORE(5)*1.d9
     &   +(MSULFLV+MBCLV+MOCLV+MDSTLV+MSALTLV)*60.d9/124.d0

           TEMPOUT1(57)=MCORE(3)*1.d9
           TEMPOUT1(58)=MCORE(2)*1.d9

           IF(IFSITEOUT==1)THEN
           AEROCOMOUT(1)=TK
           AEROCOMOUT(2)=PRESS
           AEROCOMOUT(4)=TEMPOUT1(52)

           AEROCOMOUT(5)=0.d0
           DO N=1,NSO4
             IF(GFGAS*RDRY(N)*2.d0.GE.80.d-9) THEN
               AEROCOMOUT(5)=AEROCOMOUT(5)+XN(N)
               I80nm(1)=MIN(N,I80nm(1))
             ENDIF
             IF(RRATIOBC*RDRY(N)*2.d0.GE.80.d-9) THEN
               AEROCOMOUT(5)=AEROCOMOUT(5)+XNBC(N)
               I80nm(4)=MIN(N,I80nm(4))
             ENDIF
             IF(RRATIOOC*RDRY(N)*2.d0.GE.80.d-9) THEN
               AEROCOMOUT(5)=AEROCOMOUT(5)+XNOC(N)
               I80nm(5)=MIN(N,I80nm(5))
             ENDIF
           ENDDO
           DO N=1,NSEA
             IF(RSALT(N)*2.d0.GE.80.d-9) THEN
               AEROCOMOUT(5)=AEROCOMOUT(5)+XNSALT(N)
               I80nm(2)=MIN(N,I80nm(2))
             ENDIF
           ENDDO
           DO N=1,NDSTB
             IF(RDST(N)*2.d0.GE.80.d-9) THEN
               AEROCOMOUT(5)=AEROCOMOUT(5)+XNDST(N)
               I80nm(3)=MIN(N,I80nm(3))
             ENDIF
           ENDDO
           AEROCOMOUT(6)=TEMPOUT1(53)

!Need to be updated
           BINSIZE(1)=GFGAS*(RDRY(I50nm(1)+1)-RDRY(I50nm(1)))
           BINSIZE(2)=(RSALT(I50nm(2)+1)-RSALT(I50nm(2)))
           BINSIZE(3)=(RDST(I50nm(3)+1)-RDST(I50nm(3)))
           BINSIZE(4)=RRATIOBC*(RDRY(I50nm(4)+1)-RDRY(I50nm(4)))
           BINSIZE(5)=RRATIOOC*(RDRY(I50nm(5)+1)-RDRY(I50nm(5)))

           IF(MINVAL(BINSIZE)>1.D-30)THEN
           BINMASS(1)=MINVAL(BINSIZE)/BINSIZE(1)*XN(I50nm(1))*1.7d0
           BINMASS(2)=MINVAL(BINSIZE)/BINSIZE(2)*XNSALT(I50nm(2))*2.2d0
           BINMASS(3)=MINVAL(BINSIZE)/BINSIZE(3)*XNDST(I50nm(3))*2.65d0
           BINMASS(4)=MINVAL(BINSIZE)/BINSIZE(4)*XNBC(I50nm(4))*1.8d0
           BINMASS(5)=MINVAL(BINSIZE)/BINSIZE(5)*XNOC(I50nm(5))*1.8d0

!           AEROCOMOUT(7)=(0.9d0*BINMASS(1)+
!     &                    1.28d0*BINMASS(2)+
!     &                    0.1d0*BINMASS(5))/SUM(BINMASS)

           AEROCOMOUT(7)=(ZK(1)*BINMASS(1)+
     &      ZK(2)*BINMASS(2)+ZK(3)*BINMASS(3)+
     &      ZK(4)*BINMASS(4)+ZK(5)*BINMASS(5))/SUM(BINMASS)

           ELSE
           AEROCOMOUT(7)=0.9
           ENDIF

           BINSIZE(1)=GFGAS*(RDRY(I80nm(1)+1)-RDRY(I80nm(1)))
           BINSIZE(2)=(RSALT(I80nm(2)+1)-RSALT(I80nm(2)))
           BINSIZE(3)=(RDST(I80nm(3)+1)-RDST(I80nm(3)))
           BINSIZE(4)=RRATIOBC*(RDRY(I80nm(4)+1)-RDRY(I80nm(4)))
           BINSIZE(5)=RRATIOOC*(RDRY(I80nm(5)+1)-RDRY(I80nm(5)))

           IF(MINVAL(BINSIZE)>1.D-30)THEN
           BINMASS(1)=MINVAL(BINSIZE)/BINSIZE(1)*XN(I80nm(1))*1.7d0
           BINMASS(2)=MINVAL(BINSIZE)/BINSIZE(2)*XNSALT(I80nm(2))*2.2d0
           BINMASS(3)=MINVAL(BINSIZE)/BINSIZE(3)*XNDST(I80nm(3))*2.65d0
           BINMASS(4)=MINVAL(BINSIZE)/BINSIZE(4)*XNBC(I80nm(4))*1.8d0
           BINMASS(5)=MINVAL(BINSIZE)/BINSIZE(5)*XNOC(I80nm(5))*1.8d0

!           AEROCOMOUT(8)=(0.9d0*BINMASS(1)+
!     &                    1.28d0*BINMASS(2)+
!     &                    0.1d0*BINMASS(5))/SUM(BINMASS)
           AEROCOMOUT(8)=(ZK(1)*BINMASS(1)+
     &      ZK(2)*BINMASS(2)+ZK(3)*BINMASS(3)+
     &      ZK(4)*BINMASS(4)+ZK(5)*BINMASS(5))/SUM(BINMASS)
           ELSE
           AEROCOMOUT(8)=0.9
           ENDIF

           BINSIZE(1)=GFGAS*(RDRY(I120nm(1)+1)-RDRY(I120nm(1)))
           BINSIZE(2)=(RSALT(I120nm(2)+1)-RSALT(I120nm(2)))
           BINSIZE(3)=(RDST(I120nm(3)+1)-RDST(I120nm(3)))
           BINSIZE(4)=RRATIOBC*(RDRY(I120nm(4)+1)-RDRY(I120nm(4)))
           BINSIZE(5)=RRATIOOC*(RDRY(I120nm(5)+1)-RDRY(I120nm(5)))

           IF(MINVAL(BINSIZE)>1.D-30)THEN
           BINMASS(1)=MINVAL(BINSIZE)/BINSIZE(1)*XN(I120nm(1))*1.7d0
           BINMASS(2)=MINVAL(BINSIZE)/BINSIZE(2)*XNSALT(I120nm(2))*2.2d0
           BINMASS(3)=MINVAL(BINSIZE)/BINSIZE(3)*XNDST(I120nm(3))*2.65d0
           BINMASS(4)=MINVAL(BINSIZE)/BINSIZE(4)*XNBC(I120nm(4))*1.8d0
           BINMASS(5)=MINVAL(BINSIZE)/BINSIZE(5)*XNOC(I120nm(5))*1.8d0

!           AEROCOMOUT(9)=(0.9d0*BINMASS(1)+
!     &                    1.28d0*BINMASS(2)+
!     &                    0.1d0*BINMASS(5))/SUM(BINMASS)
           AEROCOMOUT(9)=(ZK(1)*BINMASS(1)+
     &      ZK(2)*BINMASS(2)+ZK(3)*BINMASS(3)+
     &      ZK(4)*BINMASS(4)+ZK(5)*BINMASS(5))/SUM(BINMASS)
           ELSE
           AEROCOMOUT(9)=0.9
           ENDIF

           DO M=10,22

           XSAT = XSATAEROCOM(M-9) ! assume 0.05%-1.0%
           CALL DACTS(TK,MSO4,MSULFLV,MMSAB,MNITB,MNH4B,MSOA,
     &                MCORE,MSP,XSAT,ZK,ZDACT)

           DO N=1,NSO4
             IF(GFGAS*RDRY(N)*2.d0.GE.ZDACT(1,4)) THEN
               AEROCOMOUT(M)=AEROCOMOUT(M)+XN(N)
             ENDIF
             IF(RRATIOBC*RDRY(N)*2.d0.GE.ZDACT(4,4)) THEN
               AEROCOMOUT(M)=AEROCOMOUT(M)+XNBC(N)
             ENDIF
             IF(RRATIOOC*RDRY(N)*2.d0.GE.ZDACT(5,4)) THEN
               AEROCOMOUT(M)=AEROCOMOUT(M)+XNOC(N)
             ENDIF
           ENDDO
           DO N=1,NSEA
             IF(RSALT(N)*2.d0.GE.ZDACT(2,4)) THEN
               AEROCOMOUT(M)=AEROCOMOUT(M)+XNSALT(N)
             ENDIF
           ENDDO
           DO N=1,NDSTB
             IF(RDST(N)*2.d0.GE.ZDACT(3,4)) THEN
               AEROCOMOUT(M)=AEROCOMOUT(M)+XNDST(N)
             ENDIF
           ENDDO

           ENDDO

           ENDIF

         ELSE
         DO N=1,8
          TEMPOUT1(50+N) = PM25(N)
         ENDDO
         ENDIF

!OPT+
!APM2++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ENDIF

      END SUBROUTINE APM_PHYS
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: dacts
!
! !DESCRIPTION: Subroutine DACTS determines SP activation sizes for 
!  the given composition.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE DACTS( TK,    MSO4,  MSULFLV, MMSAB, 
     &                  MNITB, MNH4B, MSOA,    MCORE,
     &                  MSP,   XSAT,  ZK,      ZDACT    )
!
! !USES:
!
      USE APM_INIT_MOD, ONLY : MTACT
      USE APM_INIT_MOD, ONLY : MKACT
      USE APM_INIT_MOD, ONLY : MSAT
      USE APM_INIT_MOD, ONLY : DACTTAB
      USE APM_INIT_MOD, ONLY : NTYP
!
! !INPUT PARAMETERS: 
!
      REAL*8  :: TK
!
! !INPUT/OUTPUT PARAMETERS: 
! 
      REAL*8  :: MSO4
      REAL*8  :: MSULFLV
      REAL*8  :: MMSAB
      REAL*8  :: MNITB
      REAL*8  :: MNH4B
      REAL*8  :: MSOA
      REAL*8  :: MCORE(NTYP)
      REAL*8  :: MSP(NTYP)
      REAL*8  :: XSAT
      REAL*8  :: ZK(NTYP)
      REAL*8  :: ZDACT(NTYP,4)
!
! !REVISION HISTORY: 
!  10 Jun 2010 - F. Yu       - Initial version  
!  08 Nov 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL*8  :: YK
      REAL*8  :: MSO4a,XNH4,XNIT,XSO4a,FNH4_SO4,FSO4_NH4,MTOTAL
      REAL*8  :: KAPA(NTYP)
      INTEGER :: IT, JK, N, IS
      DATA (KAPA(N),N=1,5)/0.90,1.28,0.0,0.0,0.1/

      MSO4a = MSO4 - MSULFLV  ! MSULFLV+MSO4a = MSO4 lumped
      MSO4a = MAX(1.0d-20,MSO4a)

      XNH4 = MNH4B/18.0
      XNIT = MNITB/62.0
      XSO4a = MSO4a/96.0

      IF(XNH4.LE.XNIT) THEN
        FNH4_SO4 = 0.   ! mole fraction of NH4 associated with SO4
        FSO4_NH4 = 0.   ! mole fraction of SO4 associated with NH4
      ELSE
        FNH4_SO4 = (XNH4-XNIT)/XNH4 
        FSO4_NH4 = 0.5*(XNH4-XNIT)/XSO4a ! assumed to be (NH4)2SO4
        FSO4_NH4 = MIN(1.0d0,FSO4_NH4)
      ENDIF

      MTOTAL = MSO4a + MSULFLV + MMSAB + MNITB + MNH4B + MSOA 
!
! Assumed to have same density, thus volume fraction = mass fraction
! need to modify when the densities are different
      YK = 1./MTOTAL *((MNITB+MNH4B*(1.-FNH4_SO4))*0.67 
     &       + MSOA*0.07 + MSULFLV*0.2
     &       + (MSO4a*FSO4_NH4 + MNH4B*FNH4_SO4)*0.61
     &       + ((1.-FSO4_NH4)*MSO4a + MMSAB)*0.9)

      ZK(1) = YK

! Assume SP coated on Primary has same Kapa as SP for now
      DO N=2,NTYP
         ZK(N)=(YK*MSP(N)+KAPA(N)*MCORE(N))/(MSP(N)+MCORE(N))
      ENDDO

      IT = INT(0.5*(TK-210.0)+0.5) + 1
      IT = MIN(IT,MTACT)
      IT = MAX(IT,1)

      IS = INT(XSAT/2.d-4)
      IS = MIN(IS,MSAT)
      IS = MAX(IS,1)

      DO N=1,NTYP
        JK = INT(DLOG10(ZK(N)/1.d-5)*30.) + 1
        JK = MIN(JK,MKACT)
        JK = MAX(JK,1)

        ZDACT(N,1) = DACTTAB(IT,JK,40)   !S=0.8%
        ZDACT(N,2) = DACTTAB(IT,JK,20)   !S=0.4%
        ZDACT(N,3) = DACTTAB(IT,JK,10)   !S=0.2%
        ZDACT(N,4) = DACTTAB(IT,JK,IS)   !S= input
      ENDDO

      END SUBROUTINE DACTS
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: dacts
!
! !DESCRIPTION: Subroutine DACTS determines SP activation sizes for 
!  the given composition.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE GETGAMMAN2O5(TK, RH_P, MSO4, MSULFLV, MMSAB, 
     &                  MNITB, MNH4B, MSOA, MCORE,
     &                  MSP, GAMMAN2O5 )
!
! !USES:
!
      USE APM_INIT_MOD, ONLY : MTACT
      USE APM_INIT_MOD, ONLY : MKACT
      USE APM_INIT_MOD, ONLY : MSAT
      USE APM_INIT_MOD, ONLY : DACTTAB
      USE APM_INIT_MOD, ONLY : NTYP
!
! !INPUT PARAMETERS: 
!
      REAL*8  :: TK
      REAL*8  :: RH_P
      REAL*8  :: RH_P1
!
! !INPUT/OUTPUT PARAMETERS: 
! 
      REAL*8  :: MSO4
      REAL*8  :: MSULFLV
      REAL*8  :: MMSAB
      REAL*8  :: MNITB
      REAL*8  :: MNH4B
      REAL*8  :: MSOA
      REAL*8  :: MCORE(NTYP)
      REAL*8  :: MSP(NTYP)
      REAL*8  :: GAMMAN2O5(NTYP)
!
! !REVISION HISTORY: 
!  10 Jun 2010 - F. Yu       - Initial version  
!  08 Nov 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL*8  :: YK,TTEMP,FACT
      REAL*8  :: MSO4a,XNH4,XNIT,XSO4a,FNH4_SO4,FSO4_NH4,MTOTAL
      REAL*8  :: GAMMA0(NTYP)
      INTEGER :: IT, JK, N

      GAMMA0=0.01d0

      MSO4a = MSO4 - MSULFLV  ! MSULFLV+MSO4a = MSO4 lumped
      MSO4a = MAX(1.0d-20,MSO4a)

      XNH4 = MNH4B/18.0
      XNIT = MNITB/62.0
      XSO4a = MSO4a/96.0

      IF(XNH4.LE.XNIT) THEN
        FNH4_SO4 = 0.   ! mole fraction of NH4 associated with SO4
        FSO4_NH4 = 0.   ! mole fraction of SO4 associated with NH4
      ELSE
        FNH4_SO4 = (XNH4-XNIT)/XNH4 
        FSO4_NH4 = 0.5*(XNH4-XNIT)/XSO4a ! assumed to be (NH4)2SO4
        FSO4_NH4 = MIN(1.0d0,FSO4_NH4)
      ENDIF

      MTOTAL = MSO4a + MSULFLV + MMSAB + MNITB + MNH4B + MSOA

      !===========================================================
      ! RH dependence from Kane et al., Heterogenous uptake of
      ! gaseous N2O5 by (NH4)2SO4, NH4HSO4 and H2SO4 aerosols
      ! J. Phys. Chem. A , 2001, 105, 6465-6470
      !===========================================================
      ! No RH dependence above 50.0% (lzh, 10/26/2011)
      ! According to Bertram and Thornton, ACP, 9, 8351-8363, 2009
      RH_P1  = MIN( RH_P, 50.d0 )

      GAMMA0(1) = 2.79e-4 + RH_P1*(  1.30e-4 +
     &                      RH_P1*( -3.43e-6 +
     &                      RH_P1*(  7.52e-8 ) ) )

      !===========================================================
      ! Temperature dependence factor (Cox et al, Cambridge UK)
      ! is of the form:
      !
      !          10^( LOG10( G294 ) - 0.04 * ( TTEMP - 294 ) )
      ! FACT = -------------------------------------------------
      !                     10^( LOG10( G294 ) )
      !
      ! Where G294 = 1e-2 and TTEMP is MAX( TEMP, 282 ).
      !
      ! For computational speed, replace LOG10( 1e-2 ) with -2
      ! and replace 10^( LOG10( G294 ) ) with G294
      !===========================================================
      TTEMP = MAX( TK, 282.d0 )
      FACT  = 10.e0**( -2e+0 - 4e-2
     &      *( TTEMP - 294.e+0 ) ) / 1e-2

      ! Apply temperature dependence
      GAMMA0(1) = GAMMA0(1) * FACT

      ! Based on IUPAC recomendation
      IF ( RH_P >= 62 ) THEN
        GAMMA0(2) = 0.03e+0
      ELSE
        GAMMA0(2) = 0.005e+0
      ENDIF
      GAMMA0(3) = 0.01e+0 ! Based on unpublished Crowley work
      GAMMA0(4) = 0.005e+0 ! Based on IUPAC recomendation
      !===========================================================
      ! Based on Thornton, Braban and Abbatt, 2003
      ! N2O5 hydrolysis on sub-micron organic aerosol: the effect
      ! of relative humidity, particle phase and particle size
      !===========================================================
      IF ( RH_P >= 57e+0 ) THEN
        GAMMA0(5) = 0.03e+0
      ELSE
        GAMMA0(5) = RH_P * 5.2e-4
      ENDIF

!
! Assumed to have same density, thus volume fraction = mass fraction
! need to modify when the densities are different
      YK = 1./MTOTAL *((MNITB+MNH4B*(1.-FNH4_SO4))*0.1*GAMMA0(1)
     &       + MSOA*GAMMA0(5) + MSULFLV*GAMMA0(5)
     &       + (MSO4a*FSO4_NH4 + MNH4B*FNH4_SO4)*GAMMA0(1)
     &       + ((1.-FSO4_NH4)*MSO4a + MMSAB)*GAMMA0(1))

      GAMMAN2O5(1) = YK

! Assume SP coated on Primary has same Kapa as SP for now
      DO N=2,NTYP
        GAMMAN2O5(N)=(YK*MSP(N)+GAMMA0(N)*MCORE(N))/(MSP(N)+MCORE(N))
      ENDDO

      END SUBROUTINE GETGAMMAN2O5
!EOC

! Binary & Ternary nucl parameterization based on Dunne et al., 2016
!
      subroutine JCLOUD(TK,CH2SO4,CNH3,XQ,XN0,XJBN,XJBI,XJTN,XJTI)
      REAL*8 :: TK,CH2SO4,CNH3,XQ,XN0,XJBN,XJBI,XJTN,XJTI
      REAL*8 :: XH2SO4, XNH3, XBI, XFN, XTI, ALPHA, X, CION

      REAL*8, PARAMETER   :: Pbn = 3.95
      REAL*8, PARAMETER   :: Ubn = 9.7
      REAL*8, PARAMETER   :: Vbn = 12.6
      REAL*8, PARAMETER   :: Wbn = -0.00707
      REAL*8, PARAMETER   :: Pbi = 3.37
      REAL*8, PARAMETER   :: Ubi = -11.5
      REAL*8, PARAMETER   :: Vbi = 25.5
      REAL*8, PARAMETER   :: Wbi = 0.181

      REAL*8, PARAMETER   :: Ptn = 2.89
      REAL*8, PARAMETER   :: Utn = 182.
      REAL*8, PARAMETER   :: Vtn = 1.2
      REAL*8, PARAMETER   :: Wtn = -4.19
      REAL*8, PARAMETER   :: Pan = 8.
      REAL*8, PARAMETER   :: An  = 1.6E-6
      REAL*8, PARAMETER   :: Pti = 3.14
      REAL*8, PARAMETER   :: Uti = -23.8
      REAL*8, PARAMETER   :: Vti = 37.
      REAL*8, PARAMETER   :: Wti = 0.227
      REAL*8, PARAMETER   :: Pai = 3.07
      REAL*8, PARAMETER   :: Ai  = 0.00485

! XH2SO4 and XNH3 in the unit of 1E6/cm3
      XH2SO4=CH2SO4/1.d6
      XNH3=CNH3/1.d6

      XJBN = EXP(Ubn-EXP(Vbn*(TK/1000.-Wbn)))*(XH2SO4**Pbn)

      XFN  = XNH3/(An*XH2SO4**(-Ptn)+XNH3**(-Pan))
      XJTN = XFN*EXP(Utn-EXP(Vtn*(TK/1000.-Wtn)))
      
      XBI  = EXP(Ubi-EXP(Vbi*(TK/1000.-Wbi)))*(XH2SO4**Pbi)
      XTI  = EXP(Uti-EXP(Vti*(TK/1000.-Wti)))*
     &         XNH3/(Ai*XH2SO4**(-Pti)+XNH3**(-Pai))
      
      ALPHA = 6.E-8*SQRT(300./TK)+6.E-26*XN0*(300./TK)**4.0
      X = 1./480. + XBI + XTI
      CION = (SQRT(X*X+4.*ALPHA*XQ)-X)/(2.*ALPHA)
      XJBI = XBI*CION
      XJTI = XTI*CION
!      WRITE(6,100) ALPHA, X, XBI,XTI,CION,XJTI
! 100  FORMAT(10(1PE9.2))
      
      RETURN
      END

      END MODULE APM_PHYS_MOD
!------------------------------------------------------------------------------
#endif
