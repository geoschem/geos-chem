#ifdef APM
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: apm_init_mod
!
! !DESCRIPTION: Module APM\_INIT\_MOD contains variables and routines for 
!  initializing the APM model.
!\\
!\\
! !INTERFACE:
!
      MODULE APM_INIT_MOD
!
! !USES:
!
      IMPLICIT NONE
      PRIVATE
!
! !PRIVATE MEMBER FUNCTIONS:
!
      PRIVATE :: HYDRGF1
      PRIVATE :: BINSULF
      PRIVATE :: BINSALT 
      PRIVATE :: BINBCOC
      PRIVATE :: BINDST7 
      PRIVATE :: BINDST15 
      PRIVATE :: SALTEMIT 
      PRIVATE :: SULF_EMIT
      PRIVATE :: BCOC_EMIT
!
! !PUBLIC MEMBER FUNCTIONS:
!  
      PUBLIC  :: APM_INIT
      PUBLIC  :: APM_NTRACERS
      PUBLIC  :: CLEANUP_APMARRAYS
!
! !PUBLIC DATA MEMBERS:
!
      PUBLIC  :: LAPM
      PUBLIC  :: IFNUCL,IFEMITBCOCS, FE0
      PUBLIC  :: WBAR     !for maxsat
      PUBLIC  :: IFEV    !0= no eruptive volcano
      PUBLIC  :: IFAG 
      PUBLIC  :: IFNEWCOD,IFNEWSF,IFNEWN2O5,IFONLINEMET
      PUBLIC  :: IFDOISRP
      PUBLIC  :: IFMODISOUT

      PUBLIC  :: NGCOND,NSO4,NSEA,NDSTB,NBCOC,NBCOCT
      PUBLIC  :: NCTSO4,NCTBC,NCTOC,NCTDST,NCTSEA,NAMINE
      PUBLIC  :: NTYP
      PUBLIC  :: DEDGE,RDST,DENDST,VDST

      PUBLIC  :: APMIDS

      PUBLIC  :: YGF                                ! HYDRGF1
      PUBLIC  :: RDRY, VDRY, COAGPAR                ! BINSULF
      PUBLIC  :: RSALT, VSALT, DFMSALT9,COAGPARSS   ! BINSALT
      PUBLIC  :: IACTSS1,IACTSS3                    ! BINSALT
      PUBLIC  :: CEMITSULF, FEMTBCOCS
      PUBLIC  :: CEMITSULF2
      PUBLIC  :: RBCOC, VBCOC, COAGPBCOC   ! BINBCOC

      PUBLIC :: CEMITBCOC1  !OPT+
      PUBLIC :: CEMITBCOC2  !OPT+

      PUBLIC  :: TOTNUMBC,TOTNUMOC, TOTAREABC,TOTAREAOC
      PUBLIC  :: DACT1, DACT2, DACT3, DENSULF
      PUBLIC  :: XMACID,XMLVSOG,V1ACID,V1LVSOG,M1ACID,M1LVSOG
      PUBLIC  :: APMTRACER_MW_G,APMTRACER_MW_KG
      PUBLIC  :: ONEPI,BK,AVG,RGAS
      PUBLIC  :: MAXSITE,MSITE,ISITES,JSITES,LOUT
      PUBLIC  :: IFSITE, IFSITEOUT,IFQANN,IFEMITH,IFPREIND
      PUBLIC  :: IFALBMODIS
      PUBLIC  :: IFCURRENT
      PUBLIC  :: IFAEROCOMOUT
      PUBLIC  :: IFRADF,IFOPT,IFCDNPDF
      PUBLIC  :: IFAMINEUP, AGAMA,ERATIO,OXRATE,IFATHN,NUCLAMINE,IFNH3
      PUBLIC  :: IFNUCLORG
      PUBLIC  :: IFATOM, IFSEABIRD,IFSSTSCALE
      PUBLIC  :: ATOMMONS,ATOMDAYS,ATOMMONE,ATOMDAYE
      PUBLIC  :: SITEOUT2D
!      PUBLIC  :: REFBCOC,DBCOC1, BCOCSIGMA1
      PUBLIC  :: REFBCOC!,DBCOC1, BCOCSIGMA1, FERF
      PUBLIC  :: DBCOC2, BCOCSIGMA2
      PUBLIC  :: MTACT,MKACT,MSAT
      PUBLIC  :: DACTTAB
      PUBLIC  :: FOOACRIT
      PUBLIC  :: DATA_DIR_1x1

      PUBLIC  :: BC_LIFEAPM, OC_LIFEAPM, FHBAPM, FHOAPM,D0CONV
      PUBLIC  :: IFNOBCOCFF
      PUBLIC  :: IFCOAT
      PUBLIC  :: IFCOATBC
      PUBLIC  :: DENBC,DENOC
      PUBLIC  :: IY
      PUBLIC  :: DTTEST
      PUBLIC  :: IFAPMSF, IFAPMAOD
!      PUBLIC  :: YRLIQ,YRICE

!
! !REVISION HISTORY: 
!  16 Mar 2010 - F. Yu       - Initial version
!  08 Nov 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
      !=================================================================
      ! MODULE VARIABLES
      !=================================================================
      LOGICAL             :: LAPM  ! ON/OFF switch for APM 
      INTEGER             :: IFEV
      INTEGER             :: IFAG
      INTEGER             :: IFNEWCOD,IFNEWSF,IFNEWN2O5,IFONLINEMET
      INTEGER             :: IFDOISRP
      INTEGER             :: IFMODISOUT
      INTEGER             :: NGCOND,NSO4,NSEA,NDSTB,NBCOC,NBCOCT 
      INTEGER             :: NCTSO4,NCTBC,NCTOC,NCTDST,NCTSEA,NAMINE

      TYPE APMidtype

      INTEGER :: id_SO4G,id_SO4BIN1,id_CTSEA,id_SEABIN1,id_DSTBIN1
      INTEGER :: id_CTSO4,id_CTBC,id_CTOC,id_CTDST
      INTEGER :: id_BCBIN1,id_OCBIN1
      INTEGER :: id_AMINE

      INTEGER :: id_NH4,  id_NIT, id_SO4, id_NH3
      INTEGER :: id_CO, id_NO,   id_NO2, id_NO3, id_HNO3
      INTEGER :: id_ISOP, id_MTPA
      INTEGER :: id_BCPO,id_BCPI,id_OCPO, id_OCPI
      INTEGER :: id_MSA
      INTEGER :: id_SALA, id_SALC
      INTEGER :: id_DST1,id_DST2,id_DST3,id_DST4

      ! mono+sesq products
      INTEGER :: id_TSOA1, id_TSOG1
      INTEGER :: id_TSOA2, id_TSOG2
      INTEGER :: id_TSOA3, id_TSOG3
      INTEGER :: id_TSOA0, id_TSOG0
      ! semivolpoa: add POA (hotp 2/26/09)
      INTEGER :: id_POA1,id_POA2
      INTEGER :: id_POG1, id_POG2
      ! semivolpoa4opoa: add OPOA, OPOG (hotp 3/18/09)
      INTEGER :: id_OPOA1, id_OPOG1
      INTEGER :: id_OPOA2, id_OPOG2
      ! lumped arom/IVOC (hotp 5/14/10)
      INTEGER :: id_ASOAN
      INTEGER :: id_ASOA1, id_ASOG1
      INTEGER :: id_ASOA2, id_ASOG2
      INTEGER :: id_ASOA3, id_ASOG3
      INTEGER :: id_SO2, id_H2O2, id_OH

      END TYPE APMidtype

      TYPE(APMidtype) :: APMIDS

      REAL*8, ALLOCATABLE :: YGF(:,:)
      REAL*8, ALLOCATABLE :: RDRY(:)              !m
      REAL*8, ALLOCATABLE :: VDRY(:)            !SULFATE dry volume (m3)
      REAL*8, ALLOCATABLE :: RSALT(:), RSALT80(:) !m
      REAL*8, ALLOCATABLE :: DFMSALT9(:) ! kg m-2 s-1
      REAL*8, ALLOCATABLE :: VSALT(:)            !SALT dry volume (m3)
      REAL*8, ALLOCATABLE :: COAGPAR(:,:,:)   ! for SULFATE
      REAL*8, ALLOCATABLE :: COAGPARSS(:,:,:)   ! for SEA SALT

      REAL*8, ALLOCATABLE :: CEMITSULF(:)
      REAL*8, ALLOCATABLE :: CEMITSULF2(:,:)
      REAL*8, ALLOCATABLE :: CEMITBCOC1(:,:)    !OPT+
      REAL*8, ALLOCATABLE :: CEMITBCOC2(:,:)    !OPT+
      REAL*8, ALLOCATABLE :: RBCOC(:)
      REAL*8, ALLOCATABLE :: VBCOC(:)
      REAL*8, ALLOCATABLE :: COAGPBCOC(:,:,:)   ! BINBCOC
!
! Dust for now
      REAL*8, ALLOCATABLE :: DEDGE(:),RDST(:)
      REAL*8, ALLOCATABLE :: DENDST(:),VDST(:)

! Cloud activation dry diameters
      REAL*8, PARAMETER   :: DACT1 = 3.3d-8   !convective precipitation
      REAL*8, PARAMETER   :: DACT2 = 5.7d-8   !average
      REAL*8, PARAMETER   :: DACT3 = 8.2d-8   !large scale precipitation

      REAL*8, PARAMETER   :: DENSULF = 1.7  ! density of sulfate (g/cm3)
      REAL*8, PARAMETER   :: XMACID = 98. ! Molecular mass of H2SO4 (g/mol)
      REAL*8, PARAMETER   :: XMLVSOG = 181. ! Molecular mass of LV_SOG (g/mol)

      REAL*8, PARAMETER   :: ONEPI = 3.1415926d0
      REAL*8, PARAMETER   :: AVG = 6.022d+23
      REAL*8, PARAMETER   :: BK = 1.3807d-16  !erg/K
      REAL*8, PARAMETER   :: RGAS = 8.3144d+7  !erg K / mol

      INTEGER :: IFNUCL, IFEMITBCOCS
      REAL*8  :: FEMTBCOCS, FE0
      REAL*8  :: WBAR
      REAL*8  :: Y_R1, Y_SIGMA1, FRAC1
      REAL*8  :: Y_R2, Y_SIGMA2, FRAC2
      INTEGER :: IACTSS1, IACTSS2, IACTSS3  ! bin index for cloud act 
                                            ! diameters of seasalt

      REAL*8  :: TOTNUMBC(2),TOTNUMOC(2) !total # of BCOC per kg of BCOC for
                                         ! mod1 (fossil fuel) and mod2 (biomass/biofuel)
      REAL*8  :: TOTAREABC(2),TOTAREAOC(2) !total surface area (m2) of BCOC per kg of BCOC for
                                         ! mod1 (fossil fuel) and mod2 (biomass/biofuel)
      REAL*8  :: REFBCOC(2)              ! Effective radius of BCOC (m)
!      REAL*8  :: DBCOC1, BCOCSIGMA1, FERF(501)
      REAL*8  :: DBCOC1, BCOCSIGMA1
      REAL*8  :: DBCOC2, BCOCSIGMA2

      CHARACTER(LEN=255)   :: DATA_DIR_1x1

!      INTEGER, PARAMETER  :: MHC      = 6
      INTEGER, PARAMETER  :: MHC      = 9
      INTEGER, PARAMETER  :: NPROD    = 3
      INTEGER, PARAMETER  :: MT       = 601  

      INTEGER, PARAMETER  :: MTACT    = 51  
      INTEGER, PARAMETER  :: MKACT    = 155  
      INTEGER, PARAMETER  :: MSAT   = 50  
      REAL*8  :: DACTTAB(MTACT,MKACT,MSAT)        

      INTEGER             :: IFSITE,IFQANN,IFEMITH,IFPREIND
      INTEGER             :: IFALBMODIS
      INTEGER             :: IFCURRENT
      INTEGER             :: IFAEROCOMOUT
      INTEGER             :: IFRADF, IFOPT, IFCDNPDF
      INTEGER        :: IFAMINEUP,IFATHN,NUCLAMINE,IFNH3,IFNUCLORG
      INTEGER             :: IFATOM,IFSEABIRD,IFSSTSCALE
      INTEGER             :: ATOMMONS,ATOMDAYS,ATOMMONE,ATOMDAYE
      REAL*8              :: AGAMA,ERATIO,OXRATE

!      INTEGER,PARAMETER   :: MAXSITE=47
      INTEGER,PARAMETER   :: MAXSITE=100
!      INTEGER,PARAMETER   :: LOUT=38      !# of layers for supersite output
      INTEGER,PARAMETER   :: LOUT=28      !# of layers for supersite output
!      INTEGER,PARAMETER   :: LOUT=70      !# of layers for supersite output
!      INTEGER,PARAMETER   :: LOUT=7      !# of layers for supersite output
      INTEGER             :: MSITE
      INTEGER             :: ISITES(MAXSITE),JSITES(MAXSITE)
      INTEGER             :: IFSITEOUT(MAXSITE)
!      REAL*8  :: SITEOUT2D(MAXSITE,22)
      REAL*8  :: SITEOUT2D(MAXSITE,56)

      INTEGER,PARAMETER   :: MAPMTRA=97  ! can make this dynamic later
      REAL*8  :: APMTRACER_MW_G(MAPMTRA),APMTRACER_MW_KG(MAPMTRA)

      REAL*8  :: FOOACRIT(NPROD,MHC, MT)
      REAL*8  :: V1ACID,V1LVSOG
      REAL*8  :: M1ACID,M1LVSOG

      REAL*8  :: BC_LIFEAPM, OC_LIFEAPM, FHBAPM, FHOAPM, D0CONV
      INTEGER :: IFNOBCOCFF
      INTEGER :: IFCOAT
      INTEGER :: IFCOATBC
      INTEGER :: IFAPMSF, IFAPMAOD
      REAL*8  :: Z_R1,Z_SIGMA1,Z_R2,Z_SIGMA2
      REAL*8  :: DENBC,DENOC
      INTEGER :: SEASCHEME
      REAL*8  :: FSALTSCAL
      INTEGER :: IY
      REAL*8  :: DTTEST
!      REAL*8  :: YRLIQ,YRICE

! Types of aerosols
      INTEGER, PARAMETER  :: NTYP = 5

      CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: apm_init 
!
! !DESCRIPTION: Subroutine APM\_INIT initializes the APM model.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE APM_INIT( DATA_DIR_1x1a )
!
! !USES:
!
      USE APM_ATHN_MOD, ONLY : READJATHN
      USE APM_TIMN_MOD, ONLY : READJTIMN
      USE APM_OPTI_MOD, ONLY : READOPTABLE  !OPT+
      USE APM_OPTI_MOD, ONLY : READOPTABLE_LW  !OPT+

!
! !INPUT PARAMETERS: 
!
      CHARACTER(LEN=255) :: DATA_DIR_1x1a   ! 1x1 data directory
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
      INTEGER :: IP,IC,IT
      INTEGER :: IACT10, IACT20, IACT30   ! bin index for cloud act 
                                          ! diameters corresponding to RDRY
      REAL*8  :: TK
      INTEGER ::  ITAB     !OPT+

      !=================================================================
      ! APM_INIT begins here!
      !=================================================================

      DATA_DIR_1x1 = DATA_DIR_1x1a

      ! First-time initialization
!      IF ( FIRST ) THEN
!
! Volume and mass of one acid and one LV-SOG molecule
      V1ACID = XMACID/(6.02252d+23*DENSULF)  ! volume (cm3) of one SO4 molecule
      V1LVSOG = XMLVSOG/(6.02d+23*DENSULF)   ! assume SOA has same den as acid
      WRITE(6,*)"V1ACID=",V1ACID,"V1LVSOG=",V1LVSOG

      M1ACID = XMACID*1.E-3/6.02252E+23  ! mass (kg) of one H2SO4 molecule
      M1LVSOG = XMLVSOG*1.E-3/6.02252E+23  ! mass (kg) of one LVSOG molecule
!
!
! READ IN APM CONTROLLING PARAMETERS
!
      CLOSE(14)
      OPEN(14,FILE='input.apm')
      READ(14,*)
      READ(14,*)BC_LIFEAPM,OC_LIFEAPM,FHBAPM,FHOAPM,D0CONV
      WRITE(6,*)'BC_LIFEAPM,OC_LIFEAPM,FHBAPM,FHOAPM,D0CONV in APM'
      WRITE(6,*)BC_LIFEAPM,OC_LIFEAPM,FHBAPM,FHOAPM,D0CONV
      READ(14,*)
      READ(14,*)IFNOBCOCFF
      WRITE(6,*)"IFNOBCOCFF=",IFNOBCOCFF
      READ(14,*)
      READ(14,*)IFCOAT,IFCOATBC
      WRITE(6,*)"*****IFCOAT= ", IFCOAT,"*****IFCOATBC= ",IFCOATBC
      READ(14,*)
      READ(14,*)IFALBMODIS
      WRITE(6,*)"*****IFALBMODIS= ", IFALBMODIS
      READ(14,*)
      READ(14,*)IFAEROCOMOUT
      WRITE(6,*)"*****IFAEROCOMOUT= ", IFAEROCOMOUT
      READ(14,*)
      READ(14,*)IFCURRENT
      WRITE(6,*)"*****IFCURRENT= ", IFCURRENT
      READ(14,*)
      READ(14,*)IFRADF,IFOPT
      IF(IFRADF.EQ.1.and.IFOPT.NE.1) IFOPT = 1   !force IFOPT = 1 if IFRADF=1
      WRITE(6,*)"*****IFRADF= ", IFRADF
      WRITE(6,*)"*****IFOPT= ", IFOPT

      READ(14,*)
      READ(14,*)Z_R1,Z_SIGMA1,Z_R2,Z_SIGMA2
      WRITE(6,*)"BCOC size parameters: Z_R1,Z_SIGMA1,Z_R2,Z_SIGMA2"
      WRITE(6,*)Z_R1,Z_SIGMA1,Z_R2,Z_SIGMA2
      READ(14,*)
      READ(14,*)DENBC
      WRITE(6,*)"DENBC (density of BC in kg/m3)"
      WRITE(6,*)DENBC

      READ(14,*)
      READ(14,*)SEASCHEME
      WRITE(6,*)"SEASCHEME (1=Clarke et al 2006; 2=Gong et al 2003)"
      WRITE(6,*)SEASCHEME

      READ(14,*)
      READ(14,*)FSALTSCAL,IFSSTSCALE
      WRITE(6,*)"FSALTSCALE (when SEASCHEME=2), scale with SST=1"
      WRITE(6,*)FSALTSCAL,IFSSTSCALE

      READ(14,*)
!      READ(14,*) IFAMINEUP, AGAMA, ERATIO,OXRATE
!      WRITE(6,*)"IFAMINEUP, AGAMA, ERATIO, OXRATE"
!      WRITE(6,*)IFAMINEUP, AGAMA, ERATIO,OXRATE
      READ(14,*) IFAMINEUP, AGAMA
      WRITE(6,*)"IFAMINEUP, AGAMA"
      WRITE(6,*)IFAMINEUP, AGAMA

      READ(14,*)
      READ(14,*)IFATHN,NUCLAMINE
      WRITE(6,*)"IFATHN,NUCLAMINE"
      WRITE(6,*)IFATHN,NUCLAMINE

      READ(14,*)
      READ(14,*)IFCDNPDF, IFNEWCOD
      WRITE(6,*)"*****IFCDNPDF= ", IFCDNPDF
      WRITE(6,*)"*****IFNEWCOD= ", IFNEWCOD
      IF(IFOPT<1)IFCDNPDF=0
      IF(IFCDNPDF==0)IFNEWCOD=0

      READ(14,*)
      READ(14,*)IFNEWSF
      WRITE(6,*)"*****IFNEWSF= ", IFNEWSF

      READ(14,*)
      READ(14,*)IFNEWN2O5
      WRITE(6,*)"*****IFNEWN2O5= ", IFNEWN2O5

      READ(14,*)
      READ(14,*)IFDOISRP
      WRITE(6,*)"*****IFDOISRP= ", IFDOISRP

      READ(14,*)
      READ(14,*)IFMODISOUT
      WRITE(6,*)"*****IFMODISOUT= ", IFMODISOUT

      READ(14,*)
      READ(14,*)IFONLINEMET
      WRITE(6,*)"*****IFONLINEMET= ", IFONLINEMET

      READ(14,*)
      READ(14,*)IFNH3
      WRITE(6,*)"*****IFNH3= ", IFNH3

      READ(14,*)
      READ(14,*)IFATOM, ATOMMONS, ATOMDAYS, ATOMMONE, ATOMDAYE
      WRITE(6,*)"*****IFATOM= ", IFATOM, ATOMMONS, ATOMDAYS,
     &          ATOMMONE, ATOMDAYE

      READ(14,*)
      READ(14,*)IFSEABIRD
      WRITE(6,*)"*****IFSEABIRD= ", IFSEABIRD


      READ(14,*)
      READ(14,*)IFNUCL, IFNUCLORG
      WRITE(6,*)"*****IFNUCLORG= ", IFNUCLORG


      CLOSE(14)

! ALLOCATE MODULE VARIABLES
         CALL INIT_APMARRAYS

! READ in IMN LOOKUP TABLES
         WRITE(6,*)"NUCLEATION SCHEME INDEX IFNUCL =", IFNUCL
         IF(IFNUCL.EQ.1.or.IFNUCL.EQ.2) CALL READJTIMN(DATA_DIR_1x1)
         IF(IFATHN.GE.1) CALL READJATHN(DATA_DIR_1x1)

!OPT+
      IF(IFOPT.EQ.1) THEN
        DO ITAB = 1, 3
         CALL  READOPTABLE(ITAB,DATA_DIR_1x1)
        ENDDO

        ITAB = 3   !dust LW
        CALL  READOPTABLE_LW(ITAB,DATA_DIR_1x1)
       ENDIF
!OPT+

! Calculate hygroscopic growth factor for ammonia bisulfate, seasalt, and SOA
         CALL HYDRGF1
         
!APM2++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        IF(IFAG.EQ.1) THEN
! READ in FOOACRIT TABLE (to be used in carbon_mod.f)
         CLOSE(49)

         OPEN(49,file=TRIM(DATA_DIR_1x1)//
     &       '/APM_data_201906/FSOGagMHC9DHDC.dat'
     &          ,status='old')
!         WRITE(6,*)"READ in FOOAaging.dat"
         WRITE(6,*)"READ in FSOGagMHC9DHDC.dat"


         DO IT = 1,MT     ! 260 K -320 K
           READ(49,99)TK,((FOOACRIT(IP,IC,IT),IC=1,MHC),IP=1,NPROD)
           IF(IT.EQ.1.and.TK.NE.260.0) THEN  ! make sure read in the right one
            WRITE(6,*)"STOP at apm_init_mod.f, need check FOOAaging.dat"
            STOP
           ENDIF
           IF(MOD(IT,10).EQ.1)
     &       WRITE(6,99)TK,((FOOACRIT(IP,IC,IT),IC=1,MHC),IP=1,NPROD)
         ENDDO
         CLOSE(49)
 99      FORMAT(F5.1,50(1PE10.3))
        ENDIF
!APM2++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


! Setup bin structures and emission parameterizations
!
! Primary sulfate emission parameters
!         FE0 = 0.    ! Fraction of sulfur emitted as sulfate (0 - 1)
         
! Parameters for Nucl Mode : 1.d6 5.0d-9 1.6d0 0.05
         Y_R1 = 5.0d-9   ! cm
         Y_SIGMA1 = 1.6d0
         FRAC1 = 0.05
! Parameters for Cond Mode   : 1.d6 3.5d-8 2.0d0 0.95
         Y_R2 = 3.5d-8   ! cm
         Y_SIGMA2 = 2.0d0
         FRAC2 = 1. - FRAC1
         
         IFEMITBCOCS = 1 ! Accumlation mode of primary sulfate added to BCOC
         CALL BINSULF
         CALL SULF_EMIT

         CALL BINSALT
         CALL SALTEMIT

         IF(NDSTB.EQ.7) CALL BINDST7
         IF(NDSTB.EQ.15) CALL BINDST15

         CALL BINBCOC
         CALL BCOC_EMIT

!         FIRST = .FALSE.
!      ENDIF

      END SUBROUTINE APM_INIT 
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: apm_ntracers
!
! !DESCRIPTION: Subroutine APM\_NTRACERS initializes the number of tracers 
!  associated with APM.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE APM_NTRACERS( State_Chm )

      USE State_Chm_Mod,  ONLY : ChmState
      USE State_Chm_Mod,  ONLY : Ind_
!
! !INPUT PARAMETERS: 
      TYPE(ChmState)           :: State_Chm   ! Chemistry State object
!
! !OUTPUT PARAMETERS:
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
      INTEGER            :: N, T

      IFEV = 0
      WRITE(6,*)"IFEV=",IFEV

      IFQANN = 1
      IFAG = 0
      IFEMITH = 1
      IFPREIND = 0
      IFSITE = 0  !  
        
      WBAR = 0.3  !m/s
!      WBAR = 0.15  !m/s
      WRITE(6,*)"WBAR (m/s) = ", WBAR



!APM2++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!         IFQANN = 2    ! =1: anuual mean; =2: monthly mean
!         IFAG = 1      ! =0: no SOG age;  =1: SOG age
!         IFEMITH = 2   ! =1: original;    =2: revised emit H
!         IFSITE = 1    !  supersites =0: none;=1: global sites; =2: nested EU

      CLOSE(199)
      OPEN(199,file='apm_tmp.dat',status='old')
      READ(199,*)
      READ(199,*)IFQANN,IFAG,IFEMITH,IFSITE,IFPREIND
      READ(199,*)
      READ(199,*)IY,DTTEST
      WRITE(6,*)"IY,DTTEST = ", IY,DTTEST
      CLOSE(199)

!      CLOSE(200)
!      OPEN(200,file='apm_rad.dat',status='old')
!      READ(200,*)
!      READ(200,*)YRLIQ,YRICE
!      WRITE(6,*)"YRLIQ= ", YRLIQ, "YRICE= ", YRICE
!      CLOSE(200)
 
!APM2++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      WRITE(6,*)"IFQANN etc =",IFQANN,IFAG,IFEMITH,IFSITE,IFPREIND


! # of various APM tracers, can be organized in a better way later
      
      NGCOND=1   ! Number of condensable gases
      NSO4=40    ! Number of bins for sulfate or secondary particles

      NCTSO4=0   ! Number of tracers associated with coating on SO4
      NCTBC=1    ! Number of tracers associated with coating on BC
      NCTOC=1    ! Number of tracers associated with coating on OC
      NCTDST=1   ! Number of tracers associated with coating on DUST
      NCTSEA=1   ! Number of tracers associated with coating on Sea salt
      
!APM2++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!20110728      IF(IFAG.EQ.1) THEN
        NGCOND=2   ! Number of condensable gases
        NCTSO4=1   ! Number of tracers associated with coating on SO4
        NCTBC=2    ! Number of tracers associated with coating on BC
        NCTOC=2    ! Number of tracers associated with coating on OC
        NCTDST=2   ! Number of tracers associated with coating on DUST
        NCTSEA=2   ! Number of tracers associated with coating on Sea salt
!20110728      ENDIF
!APM2++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


      NSEA=20    ! Number of bins for sea salt
      NDSTB=15    ! Number of bins for dust
      NBCOC=15    ! Number of bins for BC and POC

!      NBCPIF=1
!      NBCPIO=1
!      NOCPIF=1
!      NOCPIO=1
!      NBCPOF=1
!      NBCPOO=1
!      NOCPOF=1
!      NOCPOO=1
!      NBCOCT=NBCPIF+NBCPIO+NOCPIF+NOCPIO+NBCPOF+NBCPOO+NOCPOF+NOCPOO

      NBCOCT=NBCOC*2   !one for BC and the other for POC

      NAMINE=3

! set APM id_ #s
      APMIDS%id_SO4G = Ind_('APMH2SO4')
      APMIDS%id_SO4BIN1 = Ind_('APMSPBIN01')
      APMIDS%id_CTSEA = Ind_('APMCTSEA1')
      APMIDS%id_SEABIN1 = Ind_('APMSEABIN01')
      APMIDS%id_DSTBIN1 = Ind_('APMDSTBIN01')
      APMIDS%id_CTSO4 = Ind_('APMLVSOA')
      APMIDS%id_CTBC = Ind_('APMCTBC1')
      APMIDS%id_CTOC = Ind_('APMCTOC1')
      APMIDS%id_CTDST = Ind_('APMCTDST1')
      APMIDS%id_BCBIN1 = Ind_('APMBCBIN01')
      APMIDS%id_OCBIN1 = Ind_('APMOCBIN01')
      APMIDS%id_AMINE = Ind_('APMAMINE1')

      APMIDS%id_OH   = Ind_('OH'       )
      APMIDS%id_H2O2 = Ind_('H2O2')
      APMIDS%id_SO2 = Ind_('SO2')
      APMIDS%id_MSA = Ind_('MSA')
      APMIDS%id_SO4 = Ind_('SO4')
      APMIDS%id_HNO3 = Ind_('HNO3')
      APMIDS%id_NIT = Ind_('NIT')
      APMIDS%id_NH3 = Ind_('NH3')
      APMIDS%id_NH4 = Ind_('NH4')

      APMIDS%id_CO = Ind_('CO')
      APMIDS%id_NO = Ind_('NO')
      APMIDS%id_NO2 = Ind_('NO2')
      APMIDS%id_NO3 = Ind_('NO3')
      APMIDS%id_ISOP = Ind_('ISOP')
      APMIDS%id_MTPA = Ind_('MTPA')

      APMIDS%id_BCPO = Ind_('BCPO')
      APMIDS%id_BCPI = Ind_('BCPI')
      APMIDS%id_OCPO = Ind_('OCPO')
      APMIDS%id_OCPI = Ind_('OCPI')

      APMIDS%id_SALA = Ind_('SALA')
      APMIDS%id_SALC = Ind_('SALC')

      APMIDS%id_DST1 = Ind_('DST1')
      APMIDS%id_DST2 = Ind_('DST2')
      APMIDS%id_DST3 = Ind_('DST3')
      APMIDS%id_DST4 = Ind_('DST4')

      APMIDS%id_TSOA1 = Ind_('TSOA1')
      APMIDS%id_TSOG1 = Ind_('TSOG1')
      APMIDS%id_TSOA2 = Ind_('TSOA2')
      APMIDS%id_TSOG2 = Ind_('TSOG2')
      APMIDS%id_TSOA3 = Ind_('TSOA3')
      APMIDS%id_TSOG3 = Ind_('TSOG3')
      APMIDS%id_TSOA0 = Ind_('TSOA0')
      APMIDS%id_TSOG0 = Ind_('TSOG0')
      APMIDS%id_POA1 = Ind_('POA1')
      APMIDS%id_POA2 = Ind_('POA2')
      APMIDS%id_POG1 = Ind_('POG1')
      APMIDS%id_POG2 = Ind_('POG2')
      APMIDS%id_OPOA1 = Ind_('OPOA1')
      APMIDS%id_OPOG1 = Ind_('OPOG1')
      APMIDS%id_OPOA2 = Ind_('OPOA2')
      APMIDS%id_OPOG2 = Ind_('OPOG2')
      APMIDS%id_ASOAN = Ind_('ASOAN')
      APMIDS%id_ASOA1 = Ind_('ASOA1')
      APMIDS%id_ASOG1 = Ind_('ASOG1')
      APMIDS%id_ASOA2 = Ind_('ASOA2')
      APMIDS%id_ASOG2 = Ind_('ASOG2')
      APMIDS%id_ASOA3 = Ind_('ASOA3')
      APMIDS%id_ASOG3 = Ind_('ASOG3')

      END SUBROUTINE APM_NTRACERS 
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: 
!
! !DESCRIPTION: Subroutine HYDRGF1 calculates the hygroscopic growth factor 
!  table.
!\\
!\\
! !INTERFACE:
!     
      SUBROUTINE HYDRGF1
!
! !REMARKS:
!   YGF(IRH,ITYPE): Growth factor at IRH=1,99 for 
!     (NH4)2SO4 (ITYPE=1), 
!     (NH4)HSO4 (ITYPE=2), 
!     SOA (ITYPE=3), and 
!     Seasalt (ITYPE=4)
! 
! !REVISION HISTORY: 
!  26 Nov 2008 - F. Yu       - Initial version  
!  08 Nov 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
      ! Parameters  for seasalt
      REAL*8,  PARAMETER     :: C1 =  0.7674d0
      REAL*8,  PARAMETER     :: C2 =  3.079d0
      REAL*8,  PARAMETER     :: C3 =  2.573d-11
      REAL*8,  PARAMETER     :: C4 = -1.424d0
!     
! !LOCAL VARIABLES:
!
      REAL*8  :: RH
      REAL*8  :: XL(4),YL(4),ZL(4)
      REAL*8  :: RSALT(10),GFSALT(10)
      REAL*8  :: X,Y,Z,RCM,RWET,FAC1,FAC2,A,B,C,SOAGF

      INTEGER   :: J,K,IRH
      REAL*8  :: YGFSALT

! For (NH4)2SO4,  NH4HSO4, H2SO4: data from Li et al, JAS, 2001 
      data XL/-0.8082E-1, 0.5121, 0.4823E-2, 1.070/ ! (NH4)2SO4
      data YL/-0.1741E-1, 0.3702, 0.1367E-1, 1.116/ ! NH4HSO4
      data ZL/0.1842, 0.4416 , 0.3492E-2, 1.053/ ! H2SO4
      
      DO J=1,99
         RH = float(J)*0.01
         X =  exp(XL(1)+XL(2)*RH+XL(3)/(RH-XL(4))**2.) ! GF for (NH4)2SO4
         Y =  exp(YL(1)+YL(2)*RH+YL(3)/(RH-YL(4))**2.) ! GF for NH4HSO4
         Z =  exp(ZL(1)+ZL(2)*RH+ZL(3)/(RH-ZL(4))**2.) ! GF for H2SO4
!     write(6,100)RH,X,Y,Z
         X = MAX(X,1.01d0)
         Y = MAX(Y,1.01d0)
         YGF(J,1) = X
         YGF(J,2) = Y
      ENDDO
!
! SOA growth factor
! Varutbangkul et al., Atmos. Chem. Phys., 6, 2367-2388, 2006.
! GF values of the pure organic portion of the SOA at 85% RH are between 
! 1.09-1.16 for the C5-C8 cycloalkenes, 1.06-1.10 for the monoterpenes and 
! oxygenated terpenes, and 1.01-1.04 for the sesquiterpenes.
!
! GF = 1 + [(1-RH/100)**(-A) x B(RH/100)**C]
!
! We choose A=0.2817, B=0.0674, C=1.7026 
! which gives GF=1.087 at RH=85% to represent average GF
!
      A= 0.2817
      B= 0.0674
      C= 1.7026

      DO J=1,99
         RH = float(J)*0.01
         SOAGF = 1. + ((1.-RH)**(-A) * B * (RH)**C)
!         write(6,105)RH,SOAGF
         YGF(J,3) = SOAGF
      ENDDO

!  Sea salt growth with relative humidity in radius [m]
! (Gerber, 1985) (bec, 12/8/04)

      DO K=1,10                 ! check to see the effect of R on GF 
         RSALT(K) = 3.E-7*float(K)*float(K)*float(K)
      ENDDO
!     WRITE(6,*)"Sea Salt Growth Factor"
!     WRITE(6,110)(RSALT(K)*1.E4,K=1,10)
!      DO J=1,99
!         RH = float(J)*0.01
!         DO K=1,10      
!! Exponential factors
!            RCM = RSALT(K)      ! cm
!            FAC1 = C1 * ( RCM**C2 )
!            FAC2 = C3 * ( RCM**C4 )
!            RWET = (FAC1/(FAC2-DLOG(RH))+RCM**3.d0)**0.33d0 ! cm
!            GFSALT(K)=RWET/RCM
!         ENDDO
!!     write(6,100)RH,(GFSALT(K),K=1,10)
!!     write(6,100)RH,GFSALT(4)
!!     Particle dry size has some effect on GF when RH is very high
!!     here use GF @ RDRY=~200 nm as average GF
!!20110513         YGF(J,4) = GFSALT(4)
!         YGF(J,4) = 2.*0.54*(1. + 1./(1.-RH))**(1./3.)  !based on De Leeuw etal., ROG,2010 for ow.
!      ENDDO
!
!Yu+ 6/2/11 read from file
      CLOSE(99)
      OPEN(99,file=TRIM(DATA_DIR_1x1)//
     &            '/APM_data_201906/APMTABLES/HGFSEASALT.txt')

      DO J=1,99
         READ(99,102)IRH, YGFSALT
         YGF(J,4) = YGFSALT
      ENDDO
      CLOSE(99)
 102  FORMAT(I5,10(1PE10.3))

      WRITE(6,*)"Calculate hygroscopic growth factor tables"

      WRITE(6,*)"   RH   SULFATE    SOA   SEASALT"
      DO J=1,99,2
         write(6,100)FLOAT(J),YGF(J,1),YGF(J,3),YGF(J,4)
      ENDDO
      
 100  FORMAT(20(1PE9.2))
 105  FORMAT(20(1PE10.3))
 110  FORMAT("RSALT(um)",10(1PE9.2))

      END SUBROUTINE HYDRGF1
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: binsulf
!
! !DESCRIPTION: Subroutine BINSULF setup size-resolved sulfate sectional bin 
!  structure (allowing any variable bin resolution,i.e., variable bin volume 
!  ratio).
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE BINSULF
! 
! !REVISION HISTORY: 
!  16 Nov 2008 - F. Yu       - Initial version  
!  08 Nov 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER :: I,J,K,L, NMAX
      REAL*8  :: RMIN,VRAT,YVRATMAX,THIRD,VJK
      REAL*8  :: RDRY_TR, RDRY_TR1,RATIOVRAT
      REAL*8  :: YVRAT(NSO4)

      THIRD = 1.d0/3.d0
      RMIN = 0.6d-9    ! m
      VRAT = 1.60d0
      YVRATMAX  = 8.0d0
      RDRY_TR = 4.5d-8  ! m
      RDRY_TR1 = 3.d-7  ! m
      RATIOVRAT = 1.15d0

      NMAX      = NSO4
      DO I   = 1, NMAX
         IF(I.EQ.1) THEN
            RDRY(I) =RMIN
!         ELSEIF((DBLE(I)/DBLE(I-1)).GE.VRAT) THEN
!            YVRAT(I-1)=DBLE(I)/DBLE(I-1)
!            RDRY(I) = RDRY(I-1)*YVRAT(I-1)**THIRD
         ELSEIF(RDRY(I-1).LE.RDRY_TR) THEN
            YVRAT(I-1)= VRAT
            RDRY(I) = RDRY(I-1)*YVRAT(I-1)**THIRD
!         ELSEIF(RDRY(I-1).LE.RDRY_TR1) THEN
!            YVRAT(I-1)= 2.0d0
!            RDRY(I) = RDRY(I-1)*YVRAT(I-1)**THIRD
         ELSE
            YVRAT(I-1) = YVRAT(I-2) * RATIOVRAT
            IF(YVRAT(I-1).GT.YVRATMAX) YVRAT(I-1) = YVRATMAX
            RDRY(I) = RDRY(I-1)*YVRAT(I-1)**THIRD
         ENDIF
         VDRY(I)=4.d0/3.d0*3.1416d0*(RDRY(I)**3.d0)   ! m3
         IF(I.EQ.NMAX) YVRAT(I) = YVRAT(I-1)
      ENDDO
      WRITE(6,*)"SULF bin information"
      WRITE(6,*)"Bin# R(um)  Volume(m3)  VRAT"
      DO I   = 1, NMAX
         WRITE(6,55)I,RDRY(I)*1.E6,VDRY(I),YVRAT(I)
      ENDDO
 55   FORMAT(I4, 10(1PE11.3))

!------------------------------------------------------------------------------
      ! Coefficient for coagulation partition
      DO J=1,NMAX
         DO K=1, J
            VJK = VDRY(J) + VDRY(K)
            IF(VJK .GE. VDRY(NMAX)) THEN
              COAGPAR(j,k,NMAX) = 1.d0
              COAGPAR(k,j,NMAX) = 1.d0
            ELSE
              IF(J .LE. (NMAX-1))THEN
              DO I=J, NMAX-1
                IF((VJK .GE. VDRY(i)) .AND. (VJK .LT. VDRY(i+1))) THEN
                  COAGPAR(j,k,i)=((VDRY(i+1)-VJK)/(VDRY(i+1)-VDRY(i)))
     &                          *VDRY(i)/VJK
                  COAGPAR(j,k,i+1) = 1.- COAGPAR(j,k,i)
                  COAGPAR(k,j,i) = COAGPAR(j,k,i)
                  COAGPAR(k,j,i+1) = COAGPAR(j,k,i+1)
                ENDIF
              ENDDO
              ENDIF
            ENDIF
         ENDDO
      ENDDO
 100  FORMAT(I3,I3,I3,F5.2,F5.2)

      END SUBROUTINE BINSULF
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: 
!
! !DESCRIPTION: Subroutine BINSALT setups size-resolved sectional bins for 
!  sea salt (allowing any variable bin resolution,i.e., variable bin volume 
!  ratio)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE BINSALT
! 
! !REVISION HISTORY: 
!  16 Nov 2008 - F. Yu       - Initial version  
!  08 Nov 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!      
      INTEGER            :: I,J,K, NMAX
      REAL*8             :: RMIN,VRAT,YVRATMAX,THIRD,VJK
      REAL*8             :: RDRY_TR, RDRY_TR1,RATIOVRAT
      REAL*8             :: YVRAT(NSEA)
      REAL*8             :: RCM, FAC1, FAC2
!
! !DEFINED PARAMETERS:
!
      REAL*8,  PARAMETER :: C1 =  0.7674d0
      REAL*8,  PARAMETER :: C2 =  3.079d0
      REAL*8,  PARAMETER :: C3 =  2.573d-11
      REAL*8,  PARAMETER :: C4 = -1.424d0

      !=================================================================
      ! BINSALT begins here!
      !=================================================================

      THIRD = 1.d0/3.d0
      RMIN = 0.6d-8    ! m
      VRAT = 2.0d0
      YVRATMAX  = 8.0d0
      RDRY_TR = 0.5d-7  ! m
      RDRY_TR1 = 3.d-0  ! m
      RATIOVRAT = 1.138d0

      NMAX      = NSEA
      DO I   = 1, NMAX
         IF(I.EQ.1) THEN
            RSALT(I) =RMIN
         ELSEIF(I.LE.3) THEN
            YVRAT(I-1)= 2.5
            RSALT(I) = RSALT(I-1)*YVRAT(I-1)**THIRD
         ELSEIF((DBLE(I)/DBLE(I-1)).GE.VRAT) THEN
            YVRAT(I-1)=DBLE(I)/DBLE(I-1)
            RSALT(I) = RSALT(I-1)*YVRAT(I-1)**THIRD
         ELSEIF(RSALT(I-1).LE.RDRY_TR) THEN
            YVRAT(I-1)= VRAT
            RSALT(I) = RSALT(I-1)*YVRAT(I-1)**THIRD
         ELSE
            YVRAT(I-1) = YVRAT(I-2) * RATIOVRAT
            IF(YVRAT(I-1).GT.YVRATMAX) YVRAT(I-1) = YVRATMAX
            RSALT(I) = RSALT(I-1)*YVRAT(I-1)**THIRD
         ENDIF
         VSALT(I)=4.d0/3.d0*3.1416d0*(RSALT(I)**3.d0)   ! m3
         IF(I.EQ.NMAX) YVRAT(I) = YVRAT(I-1)
      ENDDO
      WRITE(6,*)"Sea Salt bin information"
      WRITE(6,*)"Bin# R(um)  Volume(m3)  VRAT"

      DO I   = 1, NMAX
      ! Sea salt radius [cm]
         RCM  = RSALT(I) * 100d0

!      ! Exponential factors
!         FAC1 = C1 * ( RCM**C2 )
!         FAC2 = C3 * ( RCM**C4 )
!
!! NOTE: DLOG is replaced by LOG (bmy, 2/16/11)
!!         RSALT80(I)=0.01d0*(FAC1/(FAC2-DLOG(0.8))+RCM**3.d0)**0.33d0
!         RSALT80(I)=0.01d0*(FAC1/(FAC2-LOG(0.8))+RCM**3.d0)**0.33d0
! de Leeuw, ROG, 2010, yfq,20110513
         RSALT80(I)= 2.*RSALT(I)
         WRITE(6,55)I,RSALT(I)*1.E6,RSALT80(I)*1.E6, VSALT(I),YVRAT(I)
      ENDDO
 55   FORMAT(I4, 10(1PE11.3))

      IACTSS1= 1
      IACTSS2= 1
      IACTSS3= 1
      DO I   = 2, NMAX
       IF(DACT1.GE.(2.*RSALT(I-1)).and.DACT1.LT.(2.*RSALT(I))) IACTSS1=I
       IF(DACT2.GE.(2.*RSALT(I-1)).and.DACT2.LT.(2.*RSALT(I))) IACTSS2=I
       IF(DACT3.GE.(2.*RSALT(I-1)).and.DACT3.LT.(2.*RSALT(I))) IACTSS3=I
      ENDDO

      WRITE(6,*)"IACTSS1=",IACTSS1,"IACTSS2=",IACTSS2,"IACTSS3=",IACTSS3
!------------------------------------------------------------------------------
      ! Coefficient for coagulation partition
      DO J=1,NMAX
         DO K=1, J
            VJK = VSALT(J) + VSALT(K)
            IF(VJK .GE. VSALT(NMAX)) THEN
              COAGPARSS(j,k,NMAX) = 1.d0
              COAGPARSS(k,j,NMAX) = 1.d0
            ELSE
              IF(J .LE. (NMAX-1))THEN
              DO I=J, NMAX-1
                IF((VJK .GE. VSALT(i)) .AND. (VJK .LT. VSALT(i+1))) THEN
                  COAGPARSS(j,k,i)=((VSALT(i+1)-VJK)/(VSALT(i+1)
     &                               -VSALT(i)))*VSALT(i)/VJK
                  COAGPARSS(j,k,i+1) = 1.- COAGPARSS(j,k,i)

                  COAGPARSS(k,j,i) = COAGPARSS(j,k,i)
                  COAGPARSS(k,j,i+1) = COAGPARSS(j,k,i+1)
                ENDIF
              ENDDO
              ENDIF
            ENDIF
         ENDDO
      ENDDO
 100  FORMAT(I3,I3,I3,F5.2,F5.2)

      END SUBROUTINE BINSALT 
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: bindst7
!
! !DESCRIPTION: Subroutine BINDST setups size-resolved sectional bins for 
!  dust particles (allowing any variable bin resolution,i.e., variable bin 
!  volume ratio)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE BINDST7
! 
! !REVISION HISTORY: 
!  13 Jun 2010 - F. Yu       - Initial version  
!  08 Nov 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!      
      INTEGER :: I,J,K, NMAX
!      REAL*8  :: RMIN,VRAT,YVRATMAX,THIRD,VJK
!      REAL*8  :: RDRY_TR, RDRY_TR1,RATIOVRAT
      REAL*8  :: YVRAT(NDSTB)
      REAL*8  :: RCM

      !=================================================================
      ! BINDST begins here!
      !=================================================================

      NMAX      = NDSTB
! based on GEOS-Chem size for now
      RDST(1) = 0.15d-6  !m
      RDST(2) = 0.25d-6  !m
      RDST(3) = 0.4d-6   !m
      RDST(4) = 0.8d-6   !m
      RDST(5) = 1.5d-6   !m
      RDST(6) = 2.5d-6   !m
      RDST(7) = 4.0d-6   !m
      DENDST(1) = 2.5    !g/cm3
      DENDST(2) = 2.5    !g/cm3
      DENDST(3) = 2.5    !g/cm3
      DENDST(4) = 2.5    !g/cm3
      DENDST(5) = 2.65   !g/cm3
      DENDST(6) = 2.65   !g/cm3
      DENDST(7) = 2.65   !g/cm3

      DO I   = 1, NMAX
        VDST(I)=4.d0/3.d0*3.1416d0*(RDST(I)**3.d0)   ! m3
      ENDDO
      DO I   = 1, NMAX-1
        YVRAT(I) = VDST(I+1)/VDST(I)
      ENDDO
      YVRAT(NMAX)=YVRAT(NMAX-1)
     
      WRITE(6,*)"Dust bin information"
      WRITE(6,*)"Bin# R(um)  Volume(m3)  VRAT"

      DO I   = 1, NMAX
      ! dust radius [cm]
         RCM  = RDST(I) * 100d0
         WRITE(6,55)I,RDST(I)*1.E6, VDST(I),YVRAT(I)
      ENDDO
 55   FORMAT(I4, 10(1PE11.3))

!------------------------------------------------------------------------------
!      ! Coefficient for coagulation partition
!                 --- Deal with this later
!      DO J=1,NMAX
!         DO K=1, J
!            VJK = VDST(J) + VDST(K)
!            IF(VJK .GE. VDST(NMAX)) THEN
!              COAGPARSS(j,k,NMAX) = 1.d0
!              COAGPARSS(k,j,NMAX) = 1.d0
!            ELSE
!              IF(J .LE. (NMAX-1))THEN
!              DO I=J, NMAX-1
!                IF((VJK .GE. VDST(i)) .AND. (VJK .LT. VDST(i+1))) THEN
!                  COAGPARSS(j,k,i)=((VDST(i+1)-VJK)/(VDST(i+1)
!     &                               -VDST(i)))*VDST(i)/VJK
!                  COAGPARSS(j,k,i+1) = 1.- COAGPARSS(j,k,i)
!
!                  COAGPARSS(k,j,i) = COAGPARSS(j,k,i)
!                  COAGPARSS(k,j,i+1) = COAGPARSS(j,k,i+1)
!                ENDIF
!              ENDDO
!              ENDIF
!            ENDIF
!         ENDDO
!      ENDDO
 100  FORMAT(I3,I3,I3,F5.2,F5.2)

      END SUBROUTINE BINDST7 
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: bindst15
!
! !DESCRIPTION: Subroutine BINDST setups size-resolved sectional bins for 
!  dust particles (allowing any variable bin resolution,i.e., variable bin 
!  volume ratio)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE BINDST15
! 
! !REVISION HISTORY: 
!  13 Jun 2010 - F. Yu       - Initial version  
!  08 Nov 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER :: I,J,K, NMAX
!      REAL*8  :: RMIN,VRAT,YVRATMAX,THIRD,VJK
!      REAL*8  :: RDRY_TR, RDRY_TR1,RATIOVRAT
      REAL*8  :: YVRAT(NDSTB),DEDGE1(16)
      REAL*8  :: RCM

!      DATA DEDGE1/0.05,0.09,0.18,0.35,0.6,1.0,1.55,2.5,   ! Edge diameter (um)
!     &            3.75,4.7,5.7,7.5,14.5,26.0,41.0,63.0/

      DATA DEDGE1/0.02,0.045,0.09,0.18,0.35,0.6,1.0,1.55,2.5,   ! Edge diameter (um)
     &            3.75,5.7,7.5,14.5,26.0,41.0,63.0/

      !=================================================================
      ! BINDST begins here!
      !=================================================================

      DEDGE = DEDGE1

      NMAX      = NDSTB

      DO I   = 1, NMAX
        RDST(I)=0.5*SQRT(DEDGE(I)*DEDGE(I+1))*1.d-6  ! m
        VDST(I)=4.d0/3.d0*3.1416d0*(RDST(I)**3.d0)   ! m3
        IF(RDST(I).LT.1.d-6) THEN
           DENDST(I)= 2.5 !g/cm3
        ELSE
           DENDST(I)= 2.65 !g/cm3
        ENDIF
      ENDDO
      DO I   = 1, NMAX-1
        YVRAT(I) = VDST(I+1)/VDST(I)
      ENDDO
      YVRAT(NMAX)=YVRAT(NMAX-1)
     
      WRITE(6,*)"Dust bin information"
      WRITE(6,*)"Bin# R(um)  Volume(m3)  VRAT"

      DO I   = 1, NMAX
      ! dust radius [cm]
         RCM  = RDST(I) * 100d0
         WRITE(6,55)I,RDST(I)*1.E6, VDST(I),YVRAT(I)
      ENDDO
 55   FORMAT(I4, 10(1PE11.3))

!------------------------------------------------------------------------------
!      ! Coefficient for coagulation partition
!                 --- Deal with this later
!      DO J=1,NMAX
!         DO K=1, J
!            VJK = VDST(J) + VDST(K)
!            IF(VJK .GE. VDST(NMAX)) THEN
!              COAGPARSS(j,k,NMAX) = 1.d0
!              COAGPARSS(k,j,NMAX) = 1.d0
!            ELSE
!              IF(J .LE. (NMAX-1))THEN
!              DO I=J, NMAX-1
!                IF((VJK .GE. VDST(i)) .AND. (VJK .LT. VDST(i+1))) THEN
!                  COAGPARSS(j,k,i)=((VDST(i+1)-VJK)/(VDST(i+1)
!     &                               -VDST(i)))*VDST(i)/VJK
!                  COAGPARSS(j,k,i+1) = 1.- COAGPARSS(j,k,i)
!
!                  COAGPARSS(k,j,i) = COAGPARSS(j,k,i)
!                  COAGPARSS(k,j,i+1) = COAGPARSS(j,k,i+1)
!                ENDIF
!              ENDDO
!              ENDIF
!            ENDIF
!         ENDDO
!      ENDDO
 100  FORMAT(I3,I3,I3,F5.2,F5.2)

      END SUBROUTINE BINDST15 
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE:
!
! !DESCRIPTION: Subroutine BINBCOC setups size-resolved sectional bins for
!  BCOC (allowing any variable bin resolution,i.e., variable bin volume
!  ratio)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE BINBCOC
!
! !REVISION HISTORY:
!  1 Feb 2016 - F. Yu       - Initial version
!  08 Nov 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER            :: I,J,K, NMAX
      REAL*8             :: RMIN,VRAT,YVRATMAX,THIRD,VJK
      REAL*8             :: RDRY_TR, RDRY_TR1,RATIOVRAT
      REAL*8             :: YVRAT(NBCOC)
      REAL*8             :: RCM, FAC1, FAC2
!
! !DEFINED PARAMETERS:
!
      REAL*8,  PARAMETER :: C1 =  0.7674d0
      REAL*8,  PARAMETER :: C2 =  3.079d0
      REAL*8,  PARAMETER :: C3 =  2.573d-11
      REAL*8,  PARAMETER :: C4 = -1.424d0

      !=================================================================
      ! BINBCOC begins here!
      !=================================================================

      THIRD = 1.d0/3.d0
      RMIN = 5.0d-9    ! m

      YVRAT = 2.0d0
      YVRAT(1)=8.0d0
      YVRAT(2)=4.0d0
      YVRAT(13)=4.0d0
      YVRAT(14)=8.0d0

      NMAX      = NBCOC
      DO I   = 1, NMAX
         IF(I.EQ.1) THEN
            RBCOC(I) =RMIN
         ELSE
            RBCOC(I) = RBCOC(I-1)*YVRAT(I-1)**THIRD
         ENDIF
         VBCOC(I)=4.d0/3.d0*3.1416d0*(RBCOC(I)**3.d0)   ! m3
         IF(I.EQ.NMAX) YVRAT(I) = YVRAT(I-1)
      ENDDO
      WRITE(6,*)"BCOC bin information"
      WRITE(6,*)"Bin# R(um)  Volume(m3)  VRAT"

      DO I   = 1, NMAX
      ! BCOC radius [um]
         WRITE(6,55)I,RBCOC(I)*1.E6,VBCOC(I),YVRAT(I)
      ENDDO
 55   FORMAT(I4, 10(1PE11.3))

!------------------------------------------------------------------------------
      ! Coefficient for coagulation partition
      DO J=1,NMAX
         DO K=1, J
            VJK = VBCOC(J) + VBCOC(K)
            IF(VJK .GE. VBCOC(NMAX)) THEN
              COAGPBCOC(j,k,NMAX) = 1.d0
              COAGPBCOC(k,j,NMAX) = 1.d0
            ELSE
              IF(J .LE. (NMAX-1))THEN
              DO I=J, NMAX-1
                IF((VJK .GE. VBCOC(i)) .AND. (VJK .LT. VBCOC(i+1))) THEN
                  COAGPBCOC(j,k,i)=((VBCOC(i+1)-VJK)/(VBCOC(i+1)
     &                                -VBCOC(i)))*VBCOC(i)/VJK
                  COAGPBCOC(j,k,i+1) = 1.- COAGPBCOC(j,k,i)

                  COAGPBCOC(k,j,i) = COAGPBCOC(j,k,i)
                  COAGPBCOC(k,j,i+1) = COAGPBCOC(j,k,i+1)
                ENDIF
              ENDDO
              ENDIF
            ENDIF
         ENDDO
      ENDDO
 100  FORMAT(I3,I3,I3,F5.2,F5.2)

      END SUBROUTINE BINBCOC
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: saltemit
!
! !DESCRIPTION: Sea salt emissions?
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE SALTEMIT
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
      REAL*8  :: YVRAT(NSEA)
      REAL*8  :: DR(NSEA),REDGE(NSEA+1),DLOGDp(NSEA)
      REAL*8  :: DFSALT9(NSEA),BETA3(6,3),BETA(6)
      REAL*8  :: VRHI,VRLOW,RUM,DFDR9,DFDLOGD9
              
      INTEGER :: I,J,NMAX

      DATA BETA3/-5.001D3,0.808D6,-1.98D7,2.188D8,-1.144D9,2.29D9
     &     ,3.854D3,1.168D4,-6.572D4,1.003D5,-6.407D4,1.493D4
     &     ,4.498D2,0.839D3,-5.394D2,1.218D2,-1.213D1,4.514D-1/

!
      NMAX = NSEA
!      SEASCHEME=1
      IF(SEASCHEME.EQ.1) THEN
         WRITE(6,*)"SEASCHEME: Clarke et al., 2006"
      ELSEIF(SEASCHEME.EQ.2) THEN
         WRITE(6,*)"SEASCHEME: Gong et al., 2003"
      ENDIF

      WRITE(6,*)"SEA SALT SIZE at RH=80%"

      DO I   = 1, NMAX-1
         YVRAT(I)=(RSALT80(I+1)/RSALT80(I))**3.0
      ENDDO
      YVRAT(NMAX) = YVRAT(NMAX-1)

      DO I = 1, NMAX
         IF(I.EQ.1) THEN
            VRHI = YVRAT(I)**(1./6.)
            VRLOW= 1./VRHI
            REDGE(I) = RSALT80(I)*VRLOW
            REDGE(I+1) = RSALT80(I)*VRHI
         ELSEIF(I.EQ.NMAX) THEN
            VRLOW= 1./(YVRAT(I-1)**(1./6.))
            VRHI = 1./VRLOW
            REDGE(I+1) = RSALT80(I)*VRHI
         ELSE
            VRHI = YVRAT(I)**(1./6.)
            VRLOW= 1./( YVRAT(I-1)**(1./6.))
            REDGE(I+1) = RSALT80(I)*VRHI
         ENDIF
         DR(I) = RSALT80(I)*(VRHI-VRLOW)
         WRITE(6,100)2.E6*RSALT80(I),2.E6*DR(I),
     &   2.E6*RSALT80(I)*VRLOW,2.E6*RSALT80(I)*VRHI,
     &             VRHI,VRLOW,YVRAT(I)
      ENDDO

      IF(SEASCHEME==2)THEN
!
! Size-resolved sea salt emission based on Gong et al., 2003
!
      DO I = 1, NMAX
         RUM = RSALT80(I)*1.D6   ! Radius @RH=80% in um
         DFDR9=1.373*9.d0**3.41*   ! dFdr (m-2 um-1 s-1)
     &    (RUM**(-4.7D0*((1.D0+30.D0*RUM)
     &   **(-0.017D0*(RUM**(-1.44D0))))))
     &   * (1.d0+0.057d0*(RUM**3.45d0))
     &   * (10.d0**(1.607d0*EXP(-1.0d0
     &   * ((0.433d0-DLOG10(RUM))/0.433d0)**2.d0)))

         DFSALT9(I) =DFDR9 * DR(I)*1.E6   ! Sea-salt flux dF (# m-2 s-1)

! Sea-salt mass flux dFM (kg m-2 s-1) at U10 = 9 m/s
         DFMSALT9(I)=DFSALT9(I)*2.2d3*4.0/3.0*3.1416*RSALT(I)**3.0 

! Yu+20111015 : Somehow Gong's scheme underpredict seasalt AOD by a factor 
! of ~2 in APM simulation when compared to MODIS and MISR data, double DFMSALT9 for now
! 20111016: x2 not enough --> x2.5
!         DFMSALT9(I)= DFMSALT9(I)*2.0
!         DFMSALT9(I)= DFMSALT9(I)*2.5
! use input: 2.5 for 4x5 and 1.7 for 2x2.5
         DFMSALT9(I)= DFMSALT9(I)*FSALTSCAL

      ENDDO
      ENDIF

! Dry size for plotting
      WRITE(6,*)"SEA SALT DRY SIZE"

      DO I   = 1, NMAX-1
         YVRAT(I)=(RSALT(I+1)/RSALT(I))**3.0
      ENDDO
      YVRAT(NMAX) = YVRAT(NMAX-1)

      DO I = 1, NMAX
         IF(I.EQ.1) THEN
            VRHI = YVRAT(I)**(1./6.)
            VRLOW= 1./VRHI
            REDGE(I) = RSALT(I)*VRLOW
            REDGE(I+1) = RSALT(I)*VRHI
         ELSEIF(I.EQ.NMAX) THEN
            VRLOW= 1./(YVRAT(I-1)**(1./6.))
            VRHI = 1./VRLOW
            REDGE(I+1) = RSALT(I)*VRHI
         ELSE
            VRHI = YVRAT(I)**(1./6.)
            VRLOW= 1./( YVRAT(I-1)**(1./6.))
            REDGE(I+1) = RSALT(I)*VRHI
         ENDIF
         DR(I) = RSALT(I)*(VRHI-VRLOW)
         DLOGDp(I) = DLOG10(VRHI/VRLOW)
         WRITE(6,100)2.E6*RSALT(I),2.E6*DR(I),
     &            2.E6*RSALT(I)*VRLOW,2.E6*RSALT(I)*VRHI,
     &             VRHI,VRLOW,YVRAT(I)
      ENDDO

      IF(SEASCHEME==1)THEN
!
! Size-resolved sea salt emission based on Clarke et al., 2006
!
      DO I = 1, NMAX

         RUM = RSALT(I)*2.D6   ! Diameter in um

         IF(RUM<0.132)THEN
           DO J=1,6
             BETA(J)=BETA3(J,1)
           ENDDO
         ELSE
           IF(RUM<1.2)THEN
             DO J=1,6
               BETA(J)=BETA3(J,2)
             ENDDO
           ELSE
             DO J=1,6
               BETA(J)=BETA3(J,3)
             ENDDO
           ENDIF
         ENDIF

         DFDR9=3.84d-2*(9.d0**3.41) ! dF (m-2 s-1)
     &        *(BETA(1)+BETA(2)*RUM+BETA(3)*RUM**2+BETA(4)*RUM**3
     &         +BETA(5)*RUM**4+BETA(6)*RUM**5)*DLOGDp(I)

         DFSALT9(I) =DFDR9 ! Sea-salt flux dF (# m-2 s-1)

! Clarke parameterization invalid for last bin which is 10 um
! Set to half of last bin
         IF(I.EQ.NMAX) DFSALT9(I) = 0.5*DFSALT9(I-1) 

! Sea-salt mass flux dFM (kg m-2 s-1) at U10 = 9 m/s
         DFMSALT9(I)=DFSALT9(I)*2.2d3*4.0/3.0*3.1416*RSALT(I)**3.0
      ENDDO
      ENDIF

      DO I = 1, NMAX
         DFDLOGD9= DFSALT9(I)/DLOGDp(I) ! Sea-salt flux dF/dlogDp (# m-2 s-1)

!         WRITE(28,100)RSALT(I)*2.E6,DFDLOGD9
!         WRITE(29,100)REDGE(I)*2.E6, 1.E-20
!         WRITE(29,100)REDGE(I)*2.E6,DFDLOGD9
!         WRITE(29,100)REDGE(I+1)*2.E6,DFDLOGD9
!         WRITE(29,100)REDGE(I+1)*2.E6,1.E-20
 100     FORMAT(10(1PE10.3))
      ENDDO
      !FLUSH(28)
      !FLUSH(29)

      END SUBROUTINE SALTEMIT 
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: sulf_emit
!
! !DESCRIPTION: Subroutine SULF\_EMIT computes size-resolved aerosol emissions
!  now for primary sulfate aerosols, can include other types later on.  
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE SULF_EMIT
!
! !REMARKS:
!  CEMITSULF2(NSO4,2) is the particle mass conc (kg/m3) in each bin of each 
!  mode per unit mass (kg/m3) of primary sulfate emitted.  Multiply 
!  CEMITSULF2(NSO4,2) with total mass of primary sulfate emission to obtain 
!  sulfate mass emitted into each bin.  CEMITSULF(NSO4) combines two modes 
!  based on the fractions and IFEMITBCOCS
! 
! !REVISION HISTORY: 
!  28 Aug 2008 - F. Yu       - Initial version  
!  08 Nov 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!      
      INTEGER :: I, IMOD, NMAX
      REAL*8  :: THIRD, CORPI 
      REAL*8  :: TOTNUM, TOTMASS
      REAL*8  :: VRHI,VRLOW,VMULT1,SIGG,VRAD,VR1,CSIG1,CSIG2
      REAL*8  :: RADRAT, CONEXP, SUMFRAC, XTEMP
      REAL*8  :: YVRAT(NSO4), CFRAC(NSO4)
      REAL*8  :: DR(NSO4), DLOGDp(NSO4), REDGE(NSO4+1)
      REAL*8  :: YVMULT1(2),YGSTANDEV(2),YVOLRAD(2)
      REAL*8  :: XMI, YTEMP,YTEMP2

      !=================================================================
      ! SULF_EMIT begins here!
      !=================================================================

! Initialize the particle number size distribution with two log-normal modes
!
      NMAX=NSO4

      WRITE(6,*) "Parameters for primary sulfate emission:"
      WRITE(6,*)"Mode 1:", Y_R1, Y_SIGMA1, FRAC1
      WRITE(6,*)"Mode 2:", Y_R2, Y_SIGMA2, FRAC2

      YVMULT1(1)  =  1000.0     ! can be any postive number, cancelled during norm
      YGSTANDEV(1) =  LOG(Y_SIGMA1)
      YVOLRAD(1)  =  Y_R1
!
      YVMULT1(2)  =  1000.0
      YGSTANDEV(2) =  LOG(Y_SIGMA2)
      YVOLRAD(2)  =  Y_R2
!
      THIRD = 1.d0/3.d0
      CORPI = 1.d0/SQRT(2.d0*3.1415926d0)


      DO 783 I   = 1, NMAX-1
         YVRAT(I)=(RDRY(I+1)/RDRY(I))**3.0
 783  CONTINUE
      YVRAT(NMAX) = YVRAT(NMAX-1)

      DO 786 I = 1, NMAX
!           VRHI = ( 2.0d0*YVRAT(I)/(1.0d0 + YVRAT(I)) )**THIRD
         IF(I.EQ.1) THEN
!           VRLOW= ( 2.0d0/(1.0d0 + YVRAT(I)) )**THIRD
            VRHI = YVRAT(I)**(1./6.)
            VRLOW= 1./VRHI
            REDGE(I) = RDRY(I)*VRLOW
            REDGE(I+1) = RDRY(I)*VRHI
         ELSEIF(I.EQ.NMAX) THEN
            VRLOW= 1./(YVRAT(I-1)**(1./6.))
            VRHI = 1./VRLOW
            REDGE(I+1) = RDRY(I)*VRHI
         ELSE
!           VRLOW= ( 2.0d0/(1.0d0 + YVRAT(I-1)) )**THIRD
            VRHI = YVRAT(I)**(1./6.)
            VRLOW= 1./( YVRAT(I-1)**(1./6.))
            REDGE(I+1) = RDRY(I)*VRHI
         ENDIF
         DR(I) = RDRY(I)*(VRHI-VRLOW)
         DLOGDp(I) = LOG10(VRHI/VRLOW)
         CEMITSULF2(I,1) = 1.d-50
         CEMITSULF2(I,2) = 1.d-50
!        WRITE(24,152)2.E6*RDRY(I),2.E6*DR(I),
!     &            2.E6*RDRY(I)*VRLOW,2.E6*RDRY(I)*VRHI,
!     &             VRHI,VRLOW,YVRAT(I)
 786  CONTINUE
!
! VMULT1 is total number.
!
      DO IMOD    = 1, 2
         VMULT1       = YVMULT1(IMOD)
         SIGG         = YGSTANDEV(IMOD)
         VRAD         = YVOLRAD(IMOD)
         
         VR1         = 1.d0  / VRAD
         CSIG1       = CORPI / SIGG
         CSIG2       = 0.5d0 / (SIGG * SIGG)
         SUMFRAC     = 0.d0
         DO I    = 1, NMAX
            RADRAT   = LOG( RDRY(I) * VR1 )*1.d0
            CONEXP   = RADRAT * RADRAT * CSIG2
            CFRAC(I) = DR(I) /RDRY(I) * CSIG1 * EXP(-CONEXP)
            SUMFRAC  = SUMFRAC + CFRAC(I)
         ENDDO
!
         TOTMASS      = 0.
         DO I        = 1, NMAX
            XTEMP        = CFRAC(I)*VMULT1/SUMFRAC
            CEMITSULF2(I,IMOD) = CEMITSULF2(I,IMOD) + XTEMP ! # per volume
            TOTMASS  = TOTMASS + XTEMP*VDRY(I)*DENSULF*1.d3 ! kg per volume
         ENDDO
!
! Scale to get number conc per unit mass of total emitted particles
         TOTNUM = 0.
         DO I        = 1, NMAX
            CEMITSULF2(I,IMOD) = CEMITSULF2(I,IMOD)/TOTMASS
            TOTNUM = TOTNUM + CEMITSULF2(I,IMOD) !total number conc per unit mass
            WRITE(6,152)2*RDRY(I)*1.E9,CEMITSULF2(I,IMOD)/DLOGDp(I)
         ENDDO
! Convert # conc to mass conc
         DO I        = 1, NMAX
            CEMITSULF2(I,IMOD)=CEMITSULF2(I,IMOD)*VDRY(I)*DENSULF*1.d3
            WRITE(6,152)2*RDRY(I)*1.E9,CEMITSULF2(I,IMOD)
         ENDDO

         WRITE(6,151)IMOD
         WRITE(6,*)" sigma   r_mode   Ntots(per kg)"
         WRITE(6,152)EXP(SIGG),VRAD,TOTNUM
 151     FORMAT("Parameters of background aerosols Mode ", I3)
 152     FORMAT(10(1PE10.3))
      ENDDO

!         WRITE(21,*)" Primary Sulfate Emission Parameterization Ni/kg"
!         WRITE(21,*)"  R     Mod1   Mod2  F1*Mod1+(1-F1)*Mod2"
!         WRITE(21,*)" (m)  (#/ug) (#/ug)   (#/ug)"
!         WRITE(22,*)" Primary Sulfate Emission Parameterization 
!     &                     dNi/dlogDp /ug"
!         WRITE(22,*)"  D F1*Mod1+(1-F1)*Mod2 Mod1  Mod2 "
!         WRITE(22,*)" (um)  (#/ug) (#/ug)   (#/ug)"
      DO I        = 1, NMAX
         IF(IFEMITBCOCS.EQ.1) THEN ! put accum mode primary sulfate in BCOC
            CEMITSULF(I)=CEMITSULF2(I,1)*FRAC1 
            FEMTBCOCS = FRAC2
         ELSE
            CEMITSULF(I)=CEMITSULF2(I,1)*FRAC1 + CEMITSULF2(I,2)*FRAC2
            FEMTBCOCS =  0.d0
         ENDIF
         
         XMI = VDRY(I)*DENSULF*1.E12 ! ug (VDRY in m3, DENSULF in g/cm3)
         YTEMP = XMI*DLOGDp(I)
!            WRITE(21,152)RDRY(I),CEMITSULF2(I,1)/XMI,
!     &                CEMITSULF2(I,2)/XMI, CEMITSULF(I)/XMI
!Output dN/dlogDp (#/kg-sulfate)
            YTEMP2 = CEMITSULF2(I,2)*FRAC2/YTEMP
!            WRITE(22,152)2*RDRY(I)*1.E6,CEMITSULF(I)/YTEMP,
!     &                YTEMP2
!     &                CEMITSULF2(I,1)/YTEMP,
!     &                CEMITSULF2(I,2)/YTEMP

!            WRITE(23,152)2.E6*REDGE(I),1.E-20,1.E-20
!            WRITE(23,152)2.E6*REDGE(I),CEMITSULF(I)/YTEMP,YTEMP2
!            WRITE(23,152)2.E6*REDGE(I+1),CEMITSULF(I)/YTEMP,YTEMP2
!            WRITE(23,152)2.E6*REDGE(I+1),1.E-20,1.E-20
!            WRITE(24,152)2.E6*RDRY(I),2.E6*DR(I),
!     &            2.E6*(REDGE(I+1)-REDGE(I)),
!     &             2.E6*REDGE(I),2.E6*REDGE(I+1),YVRAT(I)**(1./3.)
      ENDDO
      !FLUSH(21)
      !FLUSH(22)
      !FLUSH(23)
      !FLUSH(24)

      END SUBROUTINE SULF_EMIT
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: bcoc_emit 
!
! !DESCRIPTION: Subroutine BCOC\_EMIT computes size-resolved aerosol 
!  emissions for primary BCOC aerosols.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE BCOC_EMIT
!
! !REMARKS:
!  CEMITBCOC2(NMAX,2) is the # of BCOC particles distributed to each bin of 
!  per unit mass (kg) of primary BCOC. Multiple CEMITBCOC2(NNMAX,2) with total 
!  BC or OC mass conc in the grid to obtain # conc.  Before the size-resolved 
!  BCOC microphysics is implemented, here only use it for convert mass conc to 
!  number conc based on two assumed log-normal distributions:  Mod1 for fossil 
!  fuel; mod2 for biomass/biofuel.
! 
! !REVISION HISTORY: 
!  22 Jan 2009 - F. Yu       - Initial version  
!  08 Nov 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER :: I, IMOD, NMAX, N
      REAL*8  :: THIRD, CORPI 
      REAL*8  :: TOTMASS
      REAL*8  :: VRHI,VRLOW,VMULT1,SIGG,VRAD,VR1,CSIG1,CSIG2
      REAL*8  :: RADRAT, CONEXP, SUMFRAC, XTEMP
      REAL*8  :: YVRAT(NBCOC), CFRAC(NBCOC)
      REAL*8  :: DR(NBCOC), DLOGDp(NBCOC), REDGE(NBCOC+1)
      REAL*8  :: YVMULT1(2),YGSTANDEV(2),YVOLRAD(2)
      REAL*8  :: XMI, YTEMP
      REAL*8  :: x
      REAL    :: FV,SS,DACT,TK
      INTEGER :: IS,IT,ISOL,INTTK


! Read in pre-calculated table for the dry diameters of particles
! that can be activated at give supersaturation as a function 
! of T,k (hygroscopicity parameter)
!
      CLOSE(211)
      OPEN(211,file=TRIM(DATA_DIR_1x1)//
     &             '/APM_data_201906/APMTABLES/DACTS_50S.txt',
     &         status='old')

      DO I=1, 15
         READ(211,*)            ! skip lines
      ENDDO
      
      DO IS=1,MSAT
         DO IT=1,MTACT
            READ(211,101)SS,TK
            READ(211,101)(DACTTAB(IT,I,IS),I=1,MKACT)
         ENDDO
      ENDDO
!      flush(212)
!      close(212)
 101  FORMAT(500(1PE9.2))
      CLOSE(211)

      !=================================================================
      ! BCOC_EMIT begins here!
      !=================================================================

! Initialize the particle number size distribution with two log-normal modes
!
      NMAX=NBCOC

      WRITE(6,*) "Parameters for primary BCOC"

      DBCOC1 = Z_R1*2.0
      BCOCSIGMA1 = Z_SIGMA1
      DBCOC2 = Z_R2*2.0
      BCOCSIGMA2 = Z_SIGMA2

! Effective mode radius
      REFBCOC(1) = Z_R1*exp(2.*LOG(Z_SIGMA1)*LOG(Z_SIGMA1))
      REFBCOC(2) = Z_R2*exp(2.*LOG(Z_SIGMA2)*LOG(Z_SIGMA2))

      WRITE(6,*)"Mode 1:", Z_R1, Z_SIGMA1, REFBCOC(1)
      WRITE(6,*)"Mode 2:", Z_R2, Z_SIGMA2, REFBCOC(2)

!YuBCden        DENBC = 1000.0   !kg/m3
!20110712      DENBC = 1800.0            !kg/m3
!      DENBC = 1200.0            !kg/m3
      DENOC = 1800.0            !kg/m3

      YVMULT1(1)  =  1000.
      YGSTANDEV(1) =  LOG(Z_SIGMA1)
      YVOLRAD(1)  =  Z_R1
!     
      YVMULT1(2)  =  30.
      YGSTANDEV(2) =  LOG(Z_SIGMA2)
      YVOLRAD(2)  =  Z_R2
!
      THIRD = 1.d0/3.d0
      CORPI = 1.d0/SQRT(2.d0*3.1415926d0)

      DO 783 I   = 1, NMAX-1
         YVRAT(I)=(RBCOC(I+1)/RBCOC(I))**3.0
 783  CONTINUE
      YVRAT(NMAX) = YVRAT(NMAX-1)

      DO 786 I = 1, NMAX
!        VRHI = ( 2.0d0*YVRAT(I)/(1.0d0 + YVRAT(I)) )**THIRD
         IF(I.EQ.1) THEN
!           VRLOW= ( 2.0d0/(1.0d0 + YVRAT(I)) )**THIRD
            VRHI = YVRAT(I)**(1./6.)
            VRLOW= 1./VRHI
            REDGE(I) = RBCOC(I)*VRLOW
            REDGE(I+1) = RBCOC(I)*VRHI
         ELSEIF(I.EQ.NMAX) THEN
            VRLOW= 1./(YVRAT(I-1)**(1./6.))
            VRHI = 1./VRLOW
            REDGE(I+1) = RBCOC(I)*VRHI
         ELSE
!           VRLOW= ( 2.0d0/(1.0d0 + YVRAT(I-1)) )**THIRD
            VRHI = YVRAT(I)**(1./6.)
            VRLOW= 1./( YVRAT(I-1)**(1./6.))
            REDGE(I+1) = RBCOC(I)*VRHI
         ENDIF
         DR(I) = RBCOC(I)*(VRHI-VRLOW)
         DLOGDp(I) = LOG10(VRHI/VRLOW)
         CEMITBCOC1(I,1) = 1.d-50
         CEMITBCOC1(I,2) = 1.d-50
         CEMITBCOC2(I,1) = 1.d-50
         CEMITBCOC2(I,2) = 1.d-50
 786  CONTINUE
!
! VMULT1 is total number.
!
      DO IMOD    = 1, 2
         VMULT1       = YVMULT1(IMOD)
         SIGG         = YGSTANDEV(IMOD)
         VRAD         = YVOLRAD(IMOD)
         
         VR1         = 1.d0  / VRAD
         CSIG1       = CORPI / SIGG
         CSIG2       = 0.5d0 / (SIGG * SIGG)
         SUMFRAC     = 0.d0
         DO I    = 1, NMAX
            RADRAT   = LOG( RBCOC(I) * VR1 )*1.d0
            CONEXP   = RADRAT * RADRAT * CSIG2
            CFRAC(I) = DR(I) /RBCOC(I) * CSIG1 * EXP(-CONEXP)
            SUMFRAC  = SUMFRAC + CFRAC(I)
         ENDDO
!     
         TOTMASS      = 0.
         DO I        = 1, NMAX
            XTEMP        = CFRAC(I)*VMULT1/SUMFRAC
            CEMITBCOC1(I,IMOD) = CEMITBCOC1(I,IMOD)+XTEMP*VBCOC(I)*DENBC
            CEMITBCOC2(I,IMOD) = CEMITBCOC2(I,IMOD) + XTEMP !# per volume
            TOTMASS      = TOTMASS + XTEMP*VBCOC(I)*DENBC ! kg per volume
         ENDDO
!
! Scale to get number of BCOC per unit mass 
         TOTNUMBC(IMOD) = 0.
         TOTAREABC(IMOD) = 0.
         DO I        = 1, NMAX
            CEMITBCOC1(I,IMOD) = CEMITBCOC1(I,IMOD)/TOTMASS
            CEMITBCOC2(I,IMOD) = CEMITBCOC2(I,IMOD)/TOTMASS !# of BCOC particles per kg of BCOC
            TOTNUMBC(IMOD) = TOTNUMBC(IMOD) + CEMITBCOC2(I,IMOD) !total number of BC per kg of BC
            TOTAREABC(IMOD) = TOTAREABC(IMOD) + ! total area of BC per kg of BC
     &         CEMITBCOC2(I,IMOD)*4.*3.1416*RBCOC(I)*RBCOC(I) ! m2/kg BC
!               WRITE(6,152)2*RBCOC(I)*1.E9,CEMITBCOC2(I,IMOD)/DLOGDp(I)

         ENDDO
! Convert # conc to mass conc
!            DO I        = 1, NMAX
!               CEMITBCOC2(I,IMOD) = CEMITBCOC2(I,IMOD)*VBCOC(I)*DENBCOC
!               WRITE(6,152)2*RBCOC(I)*1.E9,CEMITBCOC2(I,IMOD)
!            ENDDO
         TOTNUMOC(IMOD) = TOTNUMBC(IMOD)* DENBC/DENOC !total number of OC per kg of OC
         TOTAREAOC(IMOD) = TOTAREABC(IMOD) * DENBC/DENOC ! total area of OC per kg of OC
            
         WRITE(6,151)IMOD
         WRITE(6,152)EXP(SIGG),VRAD, TOTNUMBC(IMOD), TOTNUMOC(IMOD),
     &         TOTAREABC(IMOD), TOTAREAOC(IMOD)
 151     FORMAT("Parameters of PRIMARY BCOC PARTICLES for mode ", I3)
 152     FORMAT(10(1PE10.3))
      ENDDO
         
      WRITE(6,*)" FosilFuel  Biomass/Biofuel"
      WRITE(6,*)"# & surface area (m2) of BC particles per kg of BC"
      WRITE(6,152)TOTNUMBC(1),TOTNUMBC(2),TOTAREABC(1),TOTAREABC(2)
      WRITE(6,*)"# & surface area (m2) of OC particles per kg of OC"
      WRITE(6,152)TOTNUMOC(1),TOTNUMOC(2),TOTAREAOC(1),TOTAREAOC(2)

      WRITE(6,*)"BC log-normal distributioc (dN/dlogDp per kg)"//
     &          " and Mass rate"
      WRITE(6,*)"Dia(nm)   mode1     mode2"
      DO I        = 1, NMAX
         WRITE(6,152)2*RBCOC(I)*1.E9,CEMITBCOC2(I,1)/DLOGDp(I),
     &         CEMITBCOC2(I,2)/DLOGDp(I),CEMITBCOC1(I,1),CEMITBCOC1(I,2)
      ENDDO

      END SUBROUTINE BCOC_EMIT
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_apmarrays
!
! !DESCRIPTION: Subroutine INIT\_APMARRAYS allocates and zeroes module arrays.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE INIT_APMARRAYS
! 
! !REVISION HISTORY: 
!  28 Aug 2008 - F. Yu       - Initial version  
!  08 Nov 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER :: AS

      !=================================================================
      ! INIT_AEROSOL begins here!
      !=================================================================

      ALLOCATE( RDRY( NSO4 ), STAT=AS )
      RDRY = 0d0

      ALLOCATE( VDRY( NSO4 ), STAT=AS )
      VDRY = 0d0

      ALLOCATE( COAGPAR(NSO4,NSO4,NSO4 ), STAT=AS )
      COAGPAR = 0d0

      ALLOCATE( CEMITSULF(NSO4), STAT=AS )
      CEMITSULF = 0d0

      ALLOCATE( CEMITSULF2(NSO4,2), STAT=AS )
      CEMITSULF2 = 0.

!OPT+
      ALLOCATE( RBCOC(NBCOC), STAT=AS )
      RBCOC = 0.

      ALLOCATE( VBCOC(NBCOC), STAT=AS )
      VBCOC = 0.

      ALLOCATE( COAGPBCOC(NBCOC,NBCOC,NBCOC), STAT=AS )
      COAGPBCOC = 0.

      ALLOCATE( CEMITBCOC1(NBCOC,2), STAT=AS )
      CEMITBCOC1 = 0.

      ALLOCATE( CEMITBCOC2(NBCOC,2), STAT=AS )
      CEMITBCOC2 = 0.

      ALLOCATE( RSALT( NSEA ), STAT=AS )
      RSALT = 0d0

      ALLOCATE( VSALT( NSEA ), STAT=AS )
      VSALT = 0d0

      ALLOCATE( RSALT80( NSEA ), STAT=AS )
      RSALT80 = 0d0

      ALLOCATE( DFMSALT9( NSEA ), STAT=AS )
      DFMSALT9 = 0d0

      ALLOCATE( COAGPARSS(NSEA,NSEA,NSEA ), STAT=AS )
      COAGPARSS = 0d0

      ALLOCATE( COAGPBCOC(NBCOC,NBCOC,NBCOC ), STAT=AS )
      COAGPBCOC = 0d0

      ALLOCATE( YGF(99,4), STAT=AS )
      YGF = 0d0

      ALLOCATE( DEDGE( NDSTB+1 ), STAT=AS )
      DEDGE = 0d0

      ALLOCATE( RDST( NDSTB ), STAT=AS )
      RDST = 0d0

      ALLOCATE( DENDST( NDSTB ), STAT=AS )
      DENDST = 0d0

      ALLOCATE( VDST( NDSTB ), STAT=AS )
      VDST = 0d0

      END SUBROUTINE INIT_APMARRAYS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cleanup_apmarrays 
!
! !DESCRIPTION: Subroutine CLEANUP_APMARRAYS deallocates all module arrays
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CLEANUP_APMARRAYS
! 
! !REVISION HISTORY: 
!  28 Aug 2008 - F. Yu       - Initial version  
!  08 Nov 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
      !=================================================================
      ! CLEANUP_AEROSOL begins here!
      !=================================================================
      IF ( ALLOCATED( RDRY        ) ) DEALLOCATE( RDRY        )
      IF ( ALLOCATED( VDRY        ) ) DEALLOCATE( VDRY        )
      IF ( ALLOCATED( COAGPAR     ) ) DEALLOCATE( COAGPAR     )
      IF ( ALLOCATED( RSALT       ) ) DEALLOCATE( RSALT       )
      IF ( ALLOCATED( VSALT       ) ) DEALLOCATE( VSALT       )
      IF ( ALLOCATED( RSALT80     ) ) DEALLOCATE( RSALT80     )
      IF ( ALLOCATED( DFMSALT9    ) ) DEALLOCATE( DFMSALT9    )
      IF ( ALLOCATED( COAGPARSS   ) ) DEALLOCATE( COAGPARSS   )
      IF ( ALLOCATED( COAGPBCOC   ) ) DEALLOCATE( COAGPBCOC   )
      IF ( ALLOCATED( CEMITSULF   ) ) DEALLOCATE( CEMITSULF   )
      IF ( ALLOCATED( CEMITSULF2  ) ) DEALLOCATE( CEMITSULF2  )
      IF ( ALLOCATED( CEMITBCOC1  ) ) DEALLOCATE( CEMITBCOC1  )   !OPT+
      IF ( ALLOCATED( CEMITBCOC2  ) ) DEALLOCATE( CEMITBCOC2  )   !OPT+
      IF ( ALLOCATED( YGF         ) ) DEALLOCATE( YGF         )

      END SUBROUTINE CLEANUP_APMARRAYS
!EOC
      END MODULE APM_INIT_MOD
#endif
