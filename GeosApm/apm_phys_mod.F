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
     &                     NTEMPOUT1, PRESS,   TK,       RH,
     &                     XQ,        PLVSOG1, CACID,    PACID,
     &                     DT,        MMSA,    MNIT,     MNH4,
     &                     MBCS,      MOCS,    MDSTS,    MSALTS,
     &                     MBCOC8,    SOAT,    CLVSOG,   MSULFLV,
     &                     MBCLV,     MOCLV,   MDSTLV,   MSALTLV,
     &                     GFTOT1,    GFTOT2,  DENWET1,  DENWET2,
     &                     XM1D,      XN1D,    TEMPOUT1, XMDST,
     &                     FCLOUD1 )
!
! !USES:
!
      USE APM_INIT_MOD, ONLY : IFNUCL,IFAG
      USE APM_INIT_MOD, ONLY : NSO4,NSEA,NDSTB,NTYP
      USE APM_INIT_MOD, ONLY : RDRY, VDRY, RSALT, VSALT, YGF
      USE APM_INIT_MOD, ONLY : TOTNUMBC,TOTNUMOC, TOTAREABC,TOTAREAOC
      USE APM_INIT_MOD, ONLY : DACT1, DACT2, DACT3, DENSULF, XMLVSOG
      USE APM_INIT_MOD, ONLY : REFBCOC,DBCOC1, BCOCSIGMA1, FERF
      USE APM_INIT_MOD, ONLY : DBCOC2, BCOCSIGMA2
      USE APM_INIT_MOD, ONLY : V1LVSOG,V1ACID
      USE APM_INIT_MOD, ONLY : RDST,DENDST,VDST
      USE APM_INIT_MOD, ONLY : ONEPI
      USE APM_NUCL_MOD, ONLY : YUJIMN
      USE APM_COAG_MOD, ONLY : APM_COAG,APM_COAGSCAV
      USE APM_GROW_MOD, ONLY : APM_GROW,APM_MOVEBIN
      USE APM_INIT_MOD, ONLY : ISITES,JSITES
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
      INTEGER :: NCOAG1
      INTEGER :: NCOAG2
      INTEGER :: IACT10       ! bin index for cloud act 
      INTEGER :: IACT20       ! bin index for cloud act 
      INTEGER :: IACT30       ! bin index for cloud act 
      REAL*8  :: PRESS
      REAL*8  :: TK
      REAL*8  :: RH
      REAL*8  :: XQ
      REAL*8  :: PLVSOG1
      REAL*8  :: CACID
      REAL*8  :: PACID
      REAL*8  :: DT
      REAL*8  :: MMSA
      REAL*8  :: MNIT
      REAL*8  :: MNH4
      REAL*8  :: MBCS         ! mass of sulfate attached to primary particles
      REAL*8  :: MOCS         ! mass of sulfate attached to primary particles
      REAL*8  :: MDSTS        ! mass of sulfate attached to primary particles
      REAL*8  :: MSALTS       ! mass of sulfate attached to primary particles
      REAL*8  :: MBCOC8(8)
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
      REAL*8  :: XM1D(NSO4+NSEA)
      REAL*8  :: XN1D(NSO4)
      REAL*8  :: TEMPOUT1(NTEMPOUT1)
      REAL*8  :: XMDST(NDSTB)
      REAL*8  :: FCLOUD1(NSO4+4)
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
      INTEGER :: N,ITYPE, NI,NJ,NA
      INTEGER :: IRH
      REAL*8  :: PMB,XS,CACID1
      REAL*8  :: XN(NSO4),XMA(NSO4),XVA(NSO4)
      REAL*8  :: XNSALT(NSEA),XMSALT(NSEA),XVSALT(NSEA)
      REAL*8  :: XNDST(NDSTB),XVDST(NDSTB)
      REAL*8  :: TOTN
      REAL*8  :: YJ, DVJ,  DV
      REAL*8  :: VGAS,FREEP,TAREA,AREA,YKN,FCORR
      REAL*8  :: YCS(NTYP), YCSDUST, TCS
      REAL*8  :: CSSULF, AREASULF, XRCM, YJAVE, YRSTAR,CACID0
      REAL*8  :: MSO4,MSOA
      REAL*8  :: MSO4B,MNITB,MNH4B,MMSAB,MSULFT,TOTMP,FSO4B
      REAL*8  :: MSP(NTYP),MCORE(NTYP),ZK(NTYP)
      REAL*8  :: ZDACT(NTYP,3),ZCCN(NTYP,3),ZTN(NTYP),TCCN(3)
      REAL*8  :: DIAM1GAS,DIAM2GAS
      REAL*8  :: TNBC, TNOC,TNBC1,TNBC2,TNOC1,TNOC2
      REAL*8  :: GF2,GF3,GF4,FRACSOA,GFWATERVOL,GFWATER
      REAL*8  :: GFTOT,DENWET,DENWATER
      REAL*8  :: DENSALTWET,Y1
      REAL*8  :: TN10nm,TNOTHER
      REAL*8  :: RRATIO, RRATIOWET
      REAL*8  :: REFBCWET(2), REFOCWET(2),AREABC(2),AREAOC(2)
      REAL*8  :: RWET(NSO4),RSALTWET(NSEA),RWETCM(NSO4),RSALTWETCM(NSEA)
      REAL*8  :: GFGAS, CSSALT,AREASALT
      REAL*8  :: MCONDOTH,TCSOTHER, TOTCONDOTH
      REAL*8  :: YAREASALT(NSEA), MSALTTOT, YSALTS(NSEA), MSALTSOA(NSEA)
      REAL*8  :: RSALTGF(NSEA),TAREASALT
      REAL*8  :: YAREADST(NDSTB),CSDST,AREADST,RDSTWET(NDSTB)
      REAL*8  :: FBCOC1,FBCOC2
      REAL*8  :: EPC, DTCOAG1,DTCOAG2
      INTEGER :: NCOAG(2), NCOAGMAX, I3nm, IERF, I10nm
      REAL*8  :: YRHF, YAMOLF,YNTOT
      REAL*8  :: RLOSULF
      REAL*8  :: SUMXVA0,SUMXVA1,DLVSOG,DACID
      REAL*8  :: XMCONDIN,CACIDIN,PACIDIN,TCSOTHERIN
      REAL*8  :: DTLEFT, DTNGC, DNMAX
      INTEGER :: IFNGC
      INTEGER :: ISITE,JSITE
      REAL*8  :: ZN1,ZN2,ZN3
      REAL*8  :: TCSLV,CLVSOG1
      INTEGER :: ICOND

      REAL*8  :: DENAER(5)
      DATA (DENAER(NA),NA=1,5)/1.7,2.2,2.65,1.8,1.8/  ! density (g/cm3)
      DENWATER = 1.0  ! density of water (g/cm3)

      IF(PACID.LT.0.0) THEN
         WRITE(6,*) "PACID < 0, set it to 1.E-2"
         WRITE(6,99)II,JJ,LL,CACID,PACID
         PACID = 1.d-2
      ENDIF

      MCORE = 1.d-20
      MSP   = 1.d-21
      ZCCN  = 1.d-20

! The growth/movebin/coag subrountine is based on volume, convert mass to volume
! Sulfate
      MSO4 = 0.d0
      DO N=1,NSO4
        XN(N)=XN1D(N)
        XMA(N)=XM1D(N)
        MSO4 = MSO4 + XM1D(N)   ! total bin sulfate mass
        IF(XMA(N).LT.1.E-40)XMA(N)=1.E-40
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
        IF(XMSALT(N).LT.1.E-40)XMSALT(N)=1.E-40
        XVSALT(N)=XMSALT(N)*1.d-3/DENAER(2)   !XV in cm3/cm3  XM in kg/m3
        XNSALT(N) = XVSALT(N)/(1.E6*VSALT(N)) ! XN in #/cm3, VSALT in m3
      ENDDO
      ZTN(2) = SUM(XNSALT)
!
! Dust
      DO N=1,NDSTB
        MCORE(3)=MCORE(3) + XMDST(N)
        XVDST(N)=XMDST(N)*1.d-3/DENDST(N)   !XV in cm3/cm3,DENDST in g/cm3
        XNDST(N) = XVDST(N)/(1.E6*VDST(N)) ! XN in #/cm3, VDST in m3
        RDSTWET(N)=RDST(N)*100.   ! wet size in cm, use RDST for now
      ENDDO
      ZTN(3) = SUM(XNDST)

! BCOC
      MCORE(4) = MBCOC8(1)+MBCOC8(5)+MBCOC8(2)+MBCOC8(6)  !kg/m3
      MCORE(5) = MBCOC8(3)+MBCOC8(7)+MBCOC8(4)+MBCOC8(8)  !kg/m3

!******************************************************************************
! Sulfate particle dry size increase due to uptake of NIT, NH4, SOA via
! equilibrium/partition  (ratio same for all sizes for now)
! assume SO4, NIT, NH4, SOA have same density for now

! MSULFLV teated as a part of MSO4 for now but not involved in isoropia calculation
      MSO4B = MSO4-MSULFLV
      MSO4B = MAX(MSO4B,1.d-20)
      MCORE(1) = MSO4B

      DO N=1,NTYP
         IF(MCORE(N).LE.0.) THEN
           WRITE(6,*)"MCORE.LE.0", N, II, JJ,LL
           WRITE(6,100)(MBCOC8(NI),NI=1,8)
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
     &        +MCORE(5)*2.*SOAT/TOTMP
!
! determine cloud activation dry diameters at three S based on composition

      CALL DACTS(TK,MSO4,MSULFLV,MMSAB,MNITB,MNH4B,MSOA,
     &              MCORE,MSP,ZK,ZDACT)

      IF(MSO4.LE.0.) THEN
         write(6,99)II,JJ,LL,MSO4,MSULFLV,MMSA,MNIT,MNH4,MSOA
         GFGAS = 1.
      ELSE
! MSULFLV included in MSO4
         GFGAS=(MSP(1)/MSO4)**(1./3.) 
      ENDIF

! Growth due to uptake of water, need update for LV-SOA later
      IRH = INT(RH+0.5)
      IRH= MIN0(99,IRH)
      IRH= MAX0(1,IRH)
      GF2 = YGF(IRH,2)  ! growth factor of SO4+NIT+NH4 component (assuming ammonia bisulfate)
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
      DO N   = 2, NSO4
         DIAM1GAS = GFGAS*RDRY(N-1)*2.
         DIAM2GAS = GFGAS*RDRY(N)*2.
         IF(3.d-9.GE.DIAM1GAS.and.3.d-9.LT.DIAM2GAS) I3nm=N
         IF(1.d-8.GE.DIAM1GAS.and.1.d-8.LT.DIAM2GAS) I10nm=N
         IF(ZDACT(1,1).GE.DIAM1GAS.and.ZDACT(1,1).LT.DIAM2GAS) IACT10=N
         IF(ZDACT(1,2).GE.DIAM1GAS.and.ZDACT(1,2).LT.DIAM2GAS) IACT20=N
         IF(ZDACT(1,3).GE.DIAM1GAS.and.ZDACT(1,3).LT.DIAM2GAS) IACT30=N
      ENDDO
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
!
      ZCCN(1,1)=SUM(XN(IACT10:NSO4))
      ZCCN(1,2)=SUM(XN(IACT20:NSO4))
      ZCCN(1,3)=SUM(XN(IACT30:NSO4))
      ZTN(1) = SUM(XN(I3nm:NSO4))

! Get BCOC info based on assumed log-normal distributions and calculated mass.

! USe TOTNUMBC,TOTNUMOC (# per kg of BC or OC) to calculate
! TNBC, TNOC (#/cm3) from MBC, MOC (kg/m3)

      TNBC1 = (MBCOC8(1)+MBCOC8(5))*1.E-6*TOTNUMBC(1)     ! #/cm3--FF
      TNBC2 = (MBCOC8(2)+MBCOC8(6))*1.E-6*TOTNUMBC(2)   ! #/cm3--BB
      TNBC  = TNBC1 + TNBC2 
      ZTN(4) = TNBC

      TNOC1 = (MBCOC8(3)+MBCOC8(7))*1.E-6*TOTNUMOC(1)     ! #/cm3--FF
      TNOC2 = (MBCOC8(4)+MBCOC8(8))*1.E-6*TOTNUMOC(2)   ! #/cm3--BB
      TNOC  = TNOC1 + TNOC2
      ZTN(5) = TNOC

! Find BCOC CCNs
! only BCPI and OCPI can serve as CCN, the fraction of small mode depends on 
! the size of BCOC that takes into account other species associated.

      IF(MCORE(4).GT.0.D0)THEN
         RRATIO = (1. + MSP(4)/MCORE(4))**(1./3.)
         RRATIOWET=(1.+MSP(4)*GFWATER**3./MCORE(4))**(1./3.)
      ELSE
         RRATIO=1.D0
         RRATIOWET=1.D0
      ENDIF

      AREABC(1) = 1.d-2*((MBCOC8(1)+MBCOC8(5))*TOTAREABC(1))  !1.E-2 covert m2/m3 to cm2/cm3
     &  *RRATIOWET*RRATIOWET
      AREABC(2) = 1.d-2*((MBCOC8(2)+MBCOC8(6))*TOTAREABC(2))
     &   *RRATIOWET*RRATIOWET
      REFBCWET(1) = 100.*REFBCOC(1)*RRATIOWET   ! wet size for coag of sulf with BC
      REFBCWET(2) = 100.*REFBCOC(2)*RRATIOWET   ! in cm

      DO N =1, 3
       IERF = INT(LOG(ZDACT(4,N)/(RRATIO*DBCOC1))
     &    /(1.4142*LOG(BCOCSIGMA1))*100.) + 251
       IERF= MAX0(1,IERF)
       IERF= MIN0(501,IERF)
       FBCOC1 = 0.5 - 0.5*FERF(IERF)
       IERF = INT(LOG(ZDACT(4,N)/(RRATIO*DBCOC2))
     &    /(1.4142*LOG(BCOCSIGMA2))*100.) + 251
       IERF= MAX0(1,IERF)
       IERF= MIN0(501,IERF)
       FBCOC2 = 0.5 - 0.5*FERF(IERF)
       ZCCN(4,N) = MBCOC8(1)*1.E-6*TOTNUMBC(1)*FBCOC1     ! BCPI from FF
     &           + MBCOC8(2)*1.E-6*TOTNUMBC(2)*FBCOC2   ! BCPI from BB
      ENDDO

      IF(MCORE(5).GT.0.D0)THEN
         RRATIO = (1. + MSP(5)/MCORE(5))**(1./3.)
         RRATIOWET=(1.+MOCS*(GFGAS*GFWATER)**3./MCORE(5))**(1./3.)
      ELSE
         RRATIO=1.D0
         RRATIOWET=1.D0
      ENDIF

      AREAOC(1) = 1.d-2*((MBCOC8(3)+MBCOC8(7))*TOTAREAOC(1))
     &  *RRATIOWET*RRATIOWET
      AREAOC(2) = 1.d-2*((MBCOC8(4)+MBCOC8(8))*TOTAREAOC(2))
     &  *RRATIOWET*RRATIOWET

      REFOCWET(1) = 100.*REFBCOC(1)*RRATIOWET   ! wet size for coag of sulf with OC
      REFOCWET(2) = 100.*REFBCOC(2)*RRATIOWET   ! in cm

      DO N =1, 3
       IERF = INT(LOG(ZDACT(5,N)/(RRATIO*DBCOC1))
     &    /(1.4142*LOG(BCOCSIGMA1))*100.) + 251
       IERF= MAX0(1,IERF)
       IERF= MIN0(501,IERF)
       FBCOC1 = 0.5 - 0.5*FERF(IERF)
       IERF = INT(LOG(ZDACT(5,N)/(RRATIO*DBCOC2))
     &    /(1.4142*LOG(BCOCSIGMA2))*100.) + 251
       IERF= MAX0(1,IERF)
       IERF= MIN0(501,IERF)
       FBCOC2 = 0.5 - 0.5*FERF(IERF)

       ZCCN(5,N) = MBCOC8(3)*1.E-6*TOTNUMOC(1)*FBCOC1     ! OCPI from FF
     &           + MBCOC8(4)*1.E-6*TOTNUMOC(2)*FBCOC2   ! OCPI from BB
      ENDDO

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
              AREASULF = AREASULF + AREA
              YKN = FREEP/XRCM
              FCORR = YKN/(0.75+YKN) 
              CSSULF = CSSULF + 0.25*VGAS*AREA*FCORR  ! s-1
           ENDDO
           TAREA = TAREA + AREASULF
           YCS(N) = CSSULF
         ELSEIF(N.EQ.2) THEN   ! sea salt- use size distr to calculate CS
           CSSALT = 1.d-21
           AREASALT = 1.d-21
           DO NI = 1, NSEA
              XRCM = RSALTWET(NI)*100.   ! cm
              AREA = 4.*ONEPI*XRCM*XRCM*XNSALT(NI)  !cm2/cm3
              AREASALT = AREASALT + AREA
              YAREASALT(NI) = AREA
              YKN = FREEP/XRCM
              FCORR = YKN/(0.75+YKN) 
              CSSALT = CSSALT + 0.25*VGAS*AREA*FCORR  ! s-1
           ENDDO
           TAREA = TAREA + AREASALT
           YCS(N) = CSSALT
         ELSEIF(N.EQ.3) THEN   ! dust
           CSDST = 1.d-21
           AREADST = 1.d-21
           DO NI = 1, NDSTB
              XRCM = RDSTWET(NI)   ! cm
              AREA = 4.*ONEPI*XRCM*XRCM*XNDST(NI)  !cm2/cm3
              AREADST = AREADST + AREA
              YAREADST(NI) = AREA
              YKN = FREEP/XRCM
              FCORR = YKN/(0.75+YKN)
              CSDST = CSDST + 0.25*VGAS*AREA*FCORR  ! s-1
           ENDDO
           TAREA = TAREA + AREADST
           YCS(N) = CSDST
         ELSEIF(N.EQ.4) THEN   ! BC
           AREA = AREABC(1)+AREABC(2)
           TAREA = TAREA + AREA
           YKN = FREEP/REFBCWET(1)
           FCORR = YKN/(0.75+YKN)
           YCS(N) = 0.25*VGAS*AREABC(1)*FCORR  ! s-1 -- mode1

           YKN = FREEP/REFBCWET(2)
           FCORR = YKN/(0.75+YKN)
           YCS(N) = YCS(N)+0.25*VGAS*AREABC(2)*FCORR  ! s-1 add mode2
         ELSEIF(N.EQ.5) THEN   ! POC
           AREA = AREAOC(1)+AREAOC(2)
           TAREA = TAREA + AREA
           YKN = FREEP/REFOCWET(1)
           FCORR = YKN/(0.75+YKN)
           YCS(N) = 0.25*VGAS*AREAOC(1)*FCORR  ! s-1

           YKN = FREEP/REFOCWET(2)
           FCORR = YKN/(0.75+YKN)
           YCS(N) = YCS(N)+0.25*VGAS*AREAOC(2)*FCORR  ! s-1
         ENDIF
         TCS = TCS + YCS(N)   ! Total CS 
      ENDDO
      XS =  TAREA*1.E8   ! XS is total surface area in um2/cm3
      TCSOTHER = TCS - CSSULF  ! CS of particles other than sulfate

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

      DO N=1,NSEA  ! Sea-salt CCN at 3 S
         IF((2.*RSALTGF(N)).GE.ZDACT(2,1)) 
     &     ZCCN(2,1) = ZCCN(2,1) + XNSALT(N)
         IF((2.*RSALTGF(N)).GE.ZDACT(2,2)) 
     &     ZCCN(2,2) = ZCCN(2,2) + XNSALT(N)
         IF((2.*RSALTGF(N)).GE.ZDACT(2,3)) 
     &     ZCCN(2,3) = ZCCN(2,3) + XNSALT(N)
      ENDDO

! RDST below to be updated later
      DO N=1,NDSTB  
         IF((2.*RDST(N)).GE.ZDACT(3,1))
     &     ZCCN(3,1) = ZCCN(3,1) + XNDST(N)
         IF((2.*RDST(N)).GE.ZDACT(3,2))  ! dust CCN at S=~0.4%
     &     ZCCN(3,2) = ZCCN(3,2) + XNDST(N)
         IF((2.*RDST(N)).GE.ZDACT(3,3))
     &     ZCCN(3,3) = ZCCN(3,3) + XNDST(N)
      ENDDO

      ISITE = ISITES(4)
      JSITE = JSITES(4)
      IF(II.EQ.ISITE.and.JJ.EQ.JSITE.AND.LL.LE.1) THEN
        WRITE(31,100)MSO4,MSULFLV,MMSAB,MNITB,MNH4B,MSOA,FBCOC1,FBCOC2
        DO N =1, NTYP
         WRITE(31,100)MCORE(N),MSP(N),YCS(N),ZK(N),
     &       ZDACT(N,1),ZDACT(N,2),ZDACT(N,3),
     &       ZTN(N),ZCCN(N,1),ZCCN(N,2),ZCCN(N,3)
        ENDDO
        flush(31)
      ENDIF

! Fraction of BC and OC that serve as CCN0.4
      TCCN(2) = SUM(ZCCN(1:5,2))
      DO N=1,NSO4
         IF(N.LT.IACT20) THEN     ! average act diameter for sulfate aquous gain
            FCLOUD1(N) = 0.
         ELSE
            FCLOUD1(N) = XN(N)/TCCN(2)  ! parition based on number conc for now
         ENDIF
      ENDDO
      FCLOUD1(NSO4+1) = ZCCN(4,2)/TCCN(2)  ! partition based on number conc for now
      FCLOUD1(NSO4+2) = ZCCN(5,2)/TCCN(2)
      FCLOUD1(NSO4+3) = ZCCN(3,2)/TCCN(2)
      FCLOUD1(NSO4+4) = ZCCN(2,2)/TCCN(2)

!******************************************************************************
!******************************************************************************
! Reduced time step for nucleation/growth

      IFNGC = 0   ! timestep not reduced
      DNMAX = 500.  ! max change of XN1 due to nucl 

      CACID0 = CACID

      YJAVE = 0.
      DACID = 0.
      DLVSOG = 0.

      DTLEFT = DT
      IF(CACID.GT.5.d5)THEN   
 80      CONTINUE
! Nucl 
         IF(IFNUCL.EQ.1) THEN   ! IMN
            CALL YUJIMN(CACID,RH,TK,XQ,XS,YJ,YRSTAR)
         ELSEIF(IFNUCL.EQ.2) THEN   ! KBHN (XQ=1.E-20)
            CALL YUJIMN(CACID,RH,TK,XQ,XS,YJ,YRSTAR)
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
         XN(1) = XN(1) + YJ*DTNGC
         DVJ = YJ * DTNGC * VDRY(1)*1.E6   ! cm3/cm3
         XVA(1) = XVA(1) + DVJ   ! cm3/cm3
         CACID = CACID - DVJ/V1ACID

! Calculate condensational growth of H2SO4
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
         DACID = SUMXVA1 - SUMXVA0   ! amount condensed (cm3/cm3)
         CACID = CACIDIN
! SULFATE MASS condensing on particles other than sulfate 
         IF(TCSOTHER.GT.1.d-5) THEN
          MCONDOTH = TOTCONDOTH*V1ACID*DENSULF*1.d3  ! kg/m3
          MSALTS = MSALTS + MCONDOTH * YCS(2)/TCSOTHER
          MDSTS = MDSTS + MCONDOTH * YCS(3)/TCSOTHER
          MBCS   = MBCS +   MCONDOTH * YCS(4)/TCSOTHER
          MOCS   = MOCS +   MCONDOTH * YCS(5)/TCSOTHER
         ENDIF

         CALL APM_MOVEBIN(NSO4,XN,XVA)  

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
     &     REFBCWET,REFOCWET,TNBC1,TNBC2,TNOC1,TNOC2,DENAER,XNDST,
     &     XVA,MBCS,MOCS,MDSTS,MSALTS,MBCLV,MOCLV,MDSTLV,MSALTLV)
          ENDIF


         ENDIF

         IF(DTLEFT.GT.1.d-5) GOTO 80
       ELSE
         IF((PACID*DT).LT.5.E5) THEN
          CACID = CACID + PACID*DT
         ELSE
          CACID1 = PACID/TCS + (CACID - PACID/TCS)*exp(-TCS*DT)
          DACID = CACID + PACID*DT - CACID1
! Distribue all condensed acid to primary particles for now
          MCONDOTH = DACID*V1ACID*DENSULF*1.d3  ! kg/m3
          IF(TCSOTHER.GT.1.d-5) THEN
           MSALTS = MSALTS + MCONDOTH * YCS(2)/TCSOTHER
           MDSTS = MDSTS + MCONDOTH * YCS(3)/TCSOTHER
           MBCS   = MBCS +   MCONDOTH * YCS(4)/TCSOTHER
           MOCS   = MOCS +   MCONDOTH * YCS(5)/TCSOTHER
          ENDIF

          CACID = CACID1
         ENDIF
         
         YJ = 1.E-20
         YJAVE = 1.E-20
      ENDIF


      IF(CACID.LT.0.0) THEN
         WRITE(6,*) "CACID < 0, set it to 1E2/cc"
         WRITE(6,99)II,JJ,LL,CACID,CACID0,PACID,YJAVE,CSSULF,
     &   GFGAS,MMSA,MNIT,MNH4,MSOA,MSO4,
     &   TCSOTHER,TOTCONDOTH,MCONDOTH,MDSTS
         CACID = 1.E2
      ENDIF

!******************************************************************************
! Coagulation  
      EPC = 0.2      ! ratio of coag timestep to coag half lifetime
      NCOAGMAX = 240

! Seasalt
      IF(ZTN(2).GT.1.) THEN  ! only consider coag if substatial # of particle exits 
        NCOAG(2) = MAX(1,INT(EPC/(0.5d-8*ZTN(2)*DT)))
        IF((NCOAG2+1).GE.NCOAG(2)
     &               .or.(NCOAG2+1).GE.NCOAGMAX) THEN   
         DTCOAG2 = float(NCOAG2+1)*DT
         ITYPE = 4 !seasalt
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

! Sulfate
      TOTN = SUM(XN)
      IF(TOTN.GT.10.0.and.IFNGC.EQ.0) THEN  ! only consider coag if substatial # of particle exits 
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
     &   REFBCWET,REFOCWET, TNBC1,TNBC2,TNOC1,TNOC2,
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
103   FORMAT(I4,I3,I3,I3,I3,I3,I4,42(F6.3))



! tempout  
      IF(NTEMPOUT1.EQ.1) THEN
       TN10nm = SUM(XN(I10nm:NSO4))
       TNOTHER = SUM(ZTN(2:5))
       TEMPOUT1(1)=TN10nm + TNOTHER  ! total CN10 (#/cm3)
      ENDIF

      END SUBROUTINE APM_PHYS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
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
     &                  MSP,   ZK,    ZDACT           )
!
! !USES:
!
      USE APM_INIT_MOD, ONLY : MTACT
      USE APM_INIT_MOD, ONLY : MKACT
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
      REAL*8  :: ZK(NTYP)
      REAL*8  :: ZDACT(NTYP,3)
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
      INTEGER :: IT, JK, N
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

      DO N=1,NTYP
        JK = INT(DLOG10(ZK(N)/1.d-5)*30.) + 1
        JK = MIN(JK,MKACT)
        JK = MAX(JK,1)

        ZDACT(N,1) = DACTTAB(IT,JK,1)   !S=0.8%
        ZDACT(N,2) = DACTTAB(IT,JK,2)   !S=0.4%
        ZDACT(N,3) = DACTTAB(IT,JK,3)   !S=0.2%
      ENDDO

      END SUBROUTINE DACTS
!EOC
      END MODULE APM_PHYS_MOD
