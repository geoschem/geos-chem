#ifdef APM
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: apm_driv_mod
!
! !DESCRIPTION: Module APM\_DRIV\_MOD contains variables and routines to drive
!  the Advanced Particle Microphysics (APM) model.  It serves as the
!  interface between APM module and the 3D model.
!\\
!\\
! !INTERFACE:
!
MODULE APM_DRIV_MOD

  USE PRECISION_MOD
  USE PhysConstants
!
! !USES:
!
  IMPLICIT NONE
  PRIVATE
!
! !DEFINED PARAMETERS:
!
  INTEGER, PARAMETER,   PUBLIC :: NTEMPOUT1 = 82
!
! !PUBLIC DATA MEMBERS:
!
  INTEGER,              PUBLIC :: NPOUTSTEPS
  INTEGER,              PUBLIC :: NTEMPOUT
  LOGICAL,              PUBLIC :: IFTEMPOUT

  REAL*4,  ALLOCATABLE, PUBLIC :: T3DAPM(:,:,:,:,:)
  REAL*4,  ALLOCATABLE, PUBLIC :: RH3DAPM(:,:,:,:,:)
  REAL*4,  ALLOCATABLE, PUBLIC :: PBLH2DAPM(:,:,:,:)

  REAL*8,  ALLOCATABLE, PUBLIC :: EMITNH3(:,:,:)
  REAL*8,  ALLOCATABLE, PUBLIC :: EMITSO2(:,:,:)
  REAL*8,  ALLOCATABLE, PUBLIC :: DRYDEP(:,:,:,:)
  REAL*8,  ALLOCATABLE, PUBLIC :: WETDEP(:,:,:,:)
  REAL*8,  ALLOCATABLE, PUBLIC :: CONDEP(:,:,:,:)
  REAL*8,  ALLOCATABLE, PUBLIC :: UPTAKE(:,:,:,:)
  REAL*8,  ALLOCATABLE, PUBLIC :: OXIDAT(:,:,:,:)
  REAL*8,  ALLOCATABLE, PUBLIC :: FCLOUD(:,:,:,:)
  REAL*8,  ALLOCATABLE, PUBLIC :: GFTOT3D(:,:,:,:)
  REAL*8,  ALLOCATABLE, PUBLIC :: DENWET3D(:,:,:,:)
  REAL*8,  ALLOCATABLE, PUBLIC :: MWSIZE3D(:,:,:,:)
  REAL*8,  ALLOCATABLE, PUBLIC :: SO2toSO4(:,:,:)
  REAL*8,  ALLOCATABLE, PUBLIC :: PLVSOG(:,:,:,:)
  REAL*8,  ALLOCATABLE, PUBLIC :: BCLIFE(:,:,:)
  REAL*8,  ALLOCATABLE, PUBLIC :: OCLIFE(:,:,:)
  REAL*8,  ALLOCATABLE, PUBLIC :: TEMPOUT(:,:,:,:)
  REAL*8,  ALLOCATABLE, PUBLIC :: SPGF(:,:,:)
  REAL*8,  ALLOCATABLE, PUBLIC :: TAONH3(:,:,:)
  REAL*8,  ALLOCATABLE, PUBLIC :: YCOD3D(:,:,:)
  REAL*8,  ALLOCATABLE, PUBLIC :: TCOD3D(:,:,:,:)
  REAL*8,  ALLOCATABLE, PUBLIC :: YHTRC3D(:,:,:)
  REAL*8,  ALLOCATABLE, PUBLIC :: YHTR3D(:,:,:)
  REAL*8,  ALLOCATABLE, PUBLIC :: YHTRC03D(:,:,:)
  REAL*8,  ALLOCATABLE, PUBLIC :: YHTR03D(:,:,:)
  REAL*8,  ALLOCATABLE, PUBLIC :: YUV(:,:,:,:)
  REAL*8,  ALLOCATABLE, PUBLIC :: AERAREA(:,:,:,:)
  REAL*8,  ALLOCATABLE, PUBLIC :: AERDRYR(:,:,:,:)
  REAL*8,  ALLOCATABLE, PUBLIC :: GAMMAPM(:,:,:,:)

  ! Bin index for cloud act diameters taking into account uptake of
  ! NIT, NH4, SOA
  INTEGER, ALLOCATABLE, PUBLIC :: IACT1(:,:,:)
  INTEGER, ALLOCATABLE, PUBLIC :: IACT2(:,:,:)
  INTEGER, ALLOCATABLE, PUBLIC :: IACT3(:,:,:)

  ! H2SO4 gas production rate
  REAL*8,  ALLOCATABLE, PUBLIC :: PSO4GAS(:,:,:)
  REAL*8,  ALLOCATABLE, PUBLIC :: XO3(:,:,:)

  !Yu+ seasalt SST correction
  REAL*8,  ALLOCATABLE, PUBLIC :: FSST(:,:)
  REAL*8,  ALLOCATABLE, PUBLIC :: XBEXT1k3D(:,:,:,:)

  !GanLuo+ add satellite local time
  INTEGER, ALLOCATABLE :: GOOD(:)
  REAL*8,  ALLOCATABLE, PUBLIC :: CODOUT(:,:,:,:)
  REAL*8,  ALLOCATABLE, PUBLIC :: CODOUTNUM(:,:,:)

  REAL*8,  ALLOCATABLE :: MODISOUT(:,:,:)
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC  :: APM_DRIV
  PUBLIC  :: AERONUM
  PUBLIC  :: INIT_APM3D
  PUBLIC  :: CLEANUP_APM3D
  PUBLIC  :: APM_RADFOUT
!
! !PRIVATE MEMBER FUNCTIONS:
!
  PRIVATE :: APM_RADFDRIV
  PRIVATE :: APM3DQGEOS
  PRIVATE :: GET_OH
!
! !REMARKS:
!  The APM model was designed and developed for implementation into GEOS-Chem
!  by Fangqun Yu and Gan Luo at State University of New York (SUNY) at Albany.
!  (Emails: yfq@asrc.cestm.albany.edu; ganluo@asrc.cestm.albany.edu)
!                                                                             .
!  Fangqun focused on overall strategy, aerosol structure (bins, compositions,
!  types, mixing states) design, computationally efficient schemes, particle
!  microphysics (nucleation, growth, coagulation), model evaluation and
!  improvement, and model application.
!  Gan focused on overall strategy, integration of APM with GEOS-Chem, model
!  input, emission, transport, size-resolved dry deposition and wet scavenging,
!  output visualization, and model application.
!                                                                             .
!  The major reference of the APM model implemented in GEOS-Chem is:
!                                                                             .
!  1. Yu, F., and G. Luo, Simulation of particle size distribution with a
!        global aerosol model: Contribution of nucleation to aerosol and CCN
!        number concentrations, Atmos. Chem. Phys., 9, 7691-7710, 2009.
!                                                                             .
!  The APM model is an advanced multi-type, multi-component, size-resolved
!  microphysics model developed for a wide range of applications. The current
!  APM model is the result of the past developments and validation efforts
!  aimed at explaining atmospheric particle observations, reported in a
!  number of previous publications including:
!                                                                             .
!  2. Turco, R. P., Hamill, P., Toon, O. B., Whitten, R. C., and Kiang, C. S.:
!        A one-dimensional model describing aerosol formation and evolution in
!        the stratosphere, Part I, Physical processes and mathematical analogs,
!        J. Atmos. Sci., 36, 699-717, 1979.
!  3. Toon, O. B., Turco, R. P., Westphal, D., Malone, R., and Liu, M. S.: A
!        multidimensional model for aerosols: Description of computational
!        analogs, J. Atmos. Sci., 45, 2123-2143, 1988.
!  4. Jacobson, M., Turco, R., Jensen, E. and Toon O.: Modeling coagulation
!        among particles of different composition and size, Atmos. Environ.,
!        28, 1327-1338, 1994.
!  5. Jacobson, M. Z., and Turco, R. P.: Simulating condensational growth,
!        evaporation and coagulation of aerosols using a combined moving and
!        stationary size grid, Aerosol Sci. Tech., 22, 73-92, 1995.
!  6. Yu, F., and R. P. Turco: The role of ions in the formation and evolution
!        of particles in aircraft plumes, Geophys. Res. Lett., 24, 1927-1930,
!        1997.
!  7. Yu, F.: A Study of the Formation and Evolution of Aerosols and Contrails
!        in Aircraft Wakes: Development, Validation and Application of an
!        Advanced Particle Microphysics (APM) Model, Doctoral Dissertation,
!        UCLA, 1998.
!  8. Yu, F., From molecular clusters to nanoparticles: Second-generation
!        ion-mediated nucleation model, Atmos. Chem. Phys.,6, 5193-5211, 2006.
!  9. Yu, F., and R. P. Turco, Case studies of particle formation events
!        observed in boreal forests: Implications for nucleation mechanisms,
!        Atmos. Chem. Phys., 8, 6085-6102, 2008.
!                                                                             .
!  While the core components of the present APM model implemented in GEOS-Chem
!  was written from scratch using fortran 95, some algorithms and ideas were
!  inherited from the above mentioned references.
!                                                                             .
!  For more information and updates, please check the APM aerosol microphysics
!  wiki page
!  http://wiki.seas.harvard.edu/geos-chem/index.php/APM_aerosol_microphysics.
!
! !REVISION HISTORY:
!  8/2008 - 10/2010 - F. Yu & G. Luo  - Initial versions
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE TYPES:
!
  REAL*8,        ALLOCATABLE :: XN4D(:,:,:,:)
  REAL*8,  SAVE, ALLOCATABLE :: XQ3D(:,:,:)
  INTEGER,       ALLOCATABLE :: NCOAG3D(:,:,:,:)
  INTEGER,       ALLOCATABLE :: IFOUTIJ(:,:)
  INTEGER,       ALLOCATABLE :: SITEID(:,:)

  REAL*8,        ALLOCATABLE :: RCLDL3D(:,:,:)
  REAL*8,        ALLOCATABLE :: CDN3D(:,:,:)
  REAL*8,        ALLOCATABLE :: RCLDI3D(:,:,:)
  REAL,          ALLOCATABLE :: RCLDI2D(:,:)
  REAL*8,        ALLOCATABLE :: IN3D(:,:,:)

  REAL*8,        ALLOCATABLE :: MASSISRP(:,:,:,:)
  REAL*8,        ALLOCATABLE :: MASSMESA(:,:,:,:)

  !INTEGER, PARAMETER :: MWL=5   ! 0.3, 0.55, 0.94, 1.785, 3.19 um
  INTEGER, PARAMETER :: MWL=16   !
  !DATA YWLS/0.23,0.30,0.39,0.50,0.53,0.55,0.70,1.01, &
  !          1.27,1.46,1.78,2.05,2.33,2.79,3.46,8.02/
  INTEGER, PARAMETER :: KOUT1 = 3   !390 nm
  !INTEGER, PARAMETER :: KOUT2 = 6   !550 nm
  !INTEGER, PARAMETER :: KOUT2 = 4   !500 nm
  INTEGER, PARAMETER :: KOUT2 = 5   !530 nm
  INTEGER, PARAMETER :: KOUT3 = 7   !700 nm

  REAL*8,  PARAMETER   :: MAIR         = 28.966d-3           ! kg/mol

  REAL*8  :: YWLS(MWL)
  INTEGER :: IWL,ITYP
  REAL*8  :: ZBEXT(MWL),ZW(MWL),ZG(MWL)

  REAL*8,        ALLOCATABLE :: ZBEXT3D(:,:,:,:)
  REAL*8,        ALLOCATABLE :: ZW3D(:,:,:,:)
  REAL*8,        ALLOCATABLE :: ZG3D(:,:,:,:)

  !OPT each aerosol species
  INTEGER, PARAMETER :: NTYP = 5
  REAL*8  :: YBEXT(NTYP,MWL),YW(NTYP,MWL),YG(NTYP,MWL)
  REAL*8  :: XBEXT1k(40,5)

  REAL*8,        ALLOCATABLE :: YBEXT3D(:,:,:,:,:)
  REAL*8,        ALLOCATABLE :: YW3D(:,:,:,:,:)
  REAL*8,        ALLOCATABLE :: YG3D(:,:,:,:,:)

  ! mxy +longwave
  INTEGER, PARAMETER :: NWL=9   ! 4.3,5.,6.,8.1,9.6,11.6,15.8,24.0,35. um
  REAL*8  :: ZBABS(NWL)

  REAL*8,        ALLOCATABLE :: ZBABS3D(:,:,:,:)

  ! for output
  !INTEGER, PARAMETER :: NBS = 4   !CCCMa
  INTEGER, PARAMETER    :: NBS = 14   !rrtmg
  INTEGER, PARAMETER    :: NBL = 9
  INTEGER               :: NWLS(NBS)
  REAL,    ALLOCATABLE  :: ALBDRR(:,:,:)
  REAL,    ALLOCATABLE  :: ZTCST(:,:,:)
  REAL,    ALLOCATABLE  :: ZTCSB(:,:,:)
  REAL,    ALLOCATABLE  :: ZTFST(:,:,:)
  REAL,    ALLOCATABLE  :: ZTFSB(:,:,:)
  REAL,    ALLOCATABLE  :: ZTFSA(:,:,:)
  REAL,    ALLOCATABLE  :: ZTCLT(:,:,:)
  REAL,    ALLOCATABLE  :: ZTCLB(:,:,:)
  REAL,    ALLOCATABLE  :: ZTFUL(:,:,:)
  REAL,    ALLOCATABLE  :: ZTFDL(:,:,:)
  REAL,    ALLOCATABLE  :: ZTFLA(:,:,:)
  REAL,    ALLOCATABLE  :: ZTAOD(:,:,:)
  REAL,    ALLOCATABLE  :: ZCLDF(:,:)
  REAL,    ALLOCATABLE  :: ZALB(:,:)
  REAL,    ALLOCATABLE  :: ZMALB(:,:,:)
  REAL,    ALLOCATABLE  :: ZCST(:,:)
  REAL,    ALLOCATABLE  :: ZCSB(:,:)
  REAL,    ALLOCATABLE  :: ZFST(:,:)
  REAL,    ALLOCATABLE  :: ZCLD(:,:)
  REAL,    ALLOCATABLE  :: ZCLD0(:,:)
  REAL,    ALLOCATABLE  :: ZFSB(:,:)
  REAL,    ALLOCATABLE  :: ZFSA(:,:)
  REAL,    ALLOCATABLE  :: ZCLT(:,:)
  REAL,    ALLOCATABLE  :: ZCLB(:,:)
  REAL,    ALLOCATABLE  :: ZFUL(:,:)
  REAL,    ALLOCATABLE  :: ZFDL(:,:)
  REAL,    ALLOCATABLE  :: ZFLA(:,:)
  REAL,    ALLOCATABLE  :: ZAOD(:,:)
  REAL,    ALLOCATABLE  :: ZCOD(:,:)
  REAL,    ALLOCATABLE  :: ZCODGC(:,:)
  REAL,    ALLOCATABLE  :: ZAODOUT1(:,:)
  REAL,    ALLOCATABLE  :: ZAODOUT3(:,:)
  REAL,    ALLOCATABLE  :: ZAOD25(:,:)
  REAL,    ALLOCATABLE  :: ZAOD50(:,:)
  REAL,    ALLOCATABLE  :: ZAOD75(:,:)
  INTEGER, ALLOCATABLE  :: NAOD25(:,:)
  INTEGER, ALLOCATABLE  :: NAOD50(:,:)
  INTEGER, ALLOCATABLE  :: NAOD75(:,:)
  REAL,    ALLOCATABLE  :: ZVIS(:,:)
  REAL,    ALLOCATABLE  :: THAZ(:,:)
  REAL,    ALLOCATABLE  :: TFOG(:,:)
  REAL,    ALLOCATABLE  :: ZAAOD(:,:)
  REAL,    ALLOCATABLE  :: ZAAODOUT1(:,:)
  REAL,    ALLOCATABLE  :: ZAAODOUT3(:,:)
  REAL,    ALLOCATABLE  :: ZABS(:,:)
  REAL,    ALLOCATABLE  :: ZWCL(:,:)
  REAL,    ALLOCATABLE  :: ZWCI(:,:)
  REAL,    ALLOCATABLE  :: AOD(:,:)
  REAL,    ALLOCATABLE  ::  AAOD(:,:)
  REAL,    ALLOCATABLE  :: AODOUT1(:,:)
  REAL,    ALLOCATABLE  :: AAODOUT1(:,:)
  REAL,    ALLOCATABLE  :: AODOUT3(:,:)
  REAL,    ALLOCATABLE  :: AAODOUT3(:,:)
  REAL,    ALLOCATABLE  :: CODGC(:,:)
  REAL,    ALLOCATABLE  :: TAOD(:,:,:)

  !DATA YWLS/ 0.34,0.38,0.443,0.469,0.5,0.554,0.645,0.675, &
  !           0.865,0.94,1.02,1.24,1.64,1.785,2.13,3.19/ ! um
  DATA YWLS/0.23,0.30,0.39,0.50,0.53,0.55,0.70,1.01, &
            1.27,1.46,1.78,2.05,2.33,2.79,3.46,8.02/

  DATA NWLS/15,14,13,12,11,10,9,8,7,5,3,2,1,16/  !position of RRTMG WL in LT

  REAL(fp),  ALLOCATABLE :: TCOSZ(:,:)
  REAL(fp),  ALLOCATABLE :: TTDAY(:,:)
  REAL(fp),  ALLOCATABLE :: COSZM(:,:)

  REAL(fp),  PARAMETER   :: XNUMOL_OH   = 6.022e+23_fp / 17e-3_fp
  REAL(fp),  PARAMETER   :: CM3PERM3    = 1.e6_fp

  !Yu for output
  REAL(f4),  POINTER     :: OH(:,:,:)  => NULL()

CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: apm_driv
!
! !DESCRIPTION: Subroutine APM\_DRIV is the interface between APM and
!  the GEOS-Chem model.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE APM_DRIV( Input_Opt,  State_Chm, State_Diag, &
                       State_Grid, State_Met, RC )
!
! !USES:
!
    USE TIME_MOD,     ONLY : GET_TS_CHEM, ITS_A_NEW_MONTH
    USE PRESSURE_MOD, ONLY : GET_PCENTER
    USE PRESSURE_MOD, ONLY : GET_PEDGE            !
    USE APM_INIT_MOD, ONLY : NGCOND,NSO4,NSEA,NDSTB
    USE APM_INIT_MOD, ONLY : NCTSO4,NCTBC,NCTOC,NCTDST,NCTSEA,NBCOC
    USE APM_INIT_MOD, ONLY : IFNOBCOCFF
    USE APM_INIT_MOD, ONLY : IY
    USE APM_INIT_MOD, ONLY : AGAMA

    USE TIME_MOD,     ONLY : GET_YEAR, GET_MONTH
    USE TIME_MOD,     ONLY : GET_DAY,GET_HOUR,GET_MINUTE

    USE APM_PHYS_MOD, ONLY : APM_PHYS
    USE APM_INIT_MOD, ONLY : IFNUCL,IFAG,IFCOATBC,NUCLAMINE,IFATHN
    USE APM_INIT_MOD, ONLY : IFOPT  !Yu + OPT switch
    USE APM_INIT_MOD, ONLY : IFRADF  !Yu + RADF switch
    USE APM_INIT_MOD, ONLY : IFDOISRP
    USE APM_INIT_MOD, ONLY : XMACID,XMLVSOG,M1ACID,M1LVSOG
    USE APM_INIT_MOD, ONLY : RDRY,RSALT,RDST
    USE APM_COAG_MOD, ONLY : READCK6DTABLE
    USE APM_INIT_MOD, ONLY : MAXSITE,MSITE,ISITES,JSITES,LOUT
    USE APM_INIT_MOD, ONLY : SITEOUT2D
    USE APM_INIT_MOD, ONLY : IFSITE,IFSITEOUT,IFQANN,IFATOM
    USE APM_INIT_MOD, ONLY : ATOMMONS,ATOMDAYS,ATOMMONE,ATOMDAYE
    USE APM_ICEN_MOD, ONLY : nucleati
    USE APM_INIT_MOD, ONLY : APMIDS

    use rrtmg_sw_init, only: rrtmg_sw_ini
    use rrtmg_sw_GCAPM, only : cldprop_swapm

    USE ErrCode_Mod
    USE Input_Opt_Mod,  ONLY : OptInput
    USE State_Met_Mod,  ONLY : MetState
    USE State_Chm_Mod,  ONLY : ChmState
    USE State_Chm_Mod,  ONLY : Ind_
    USE State_Grid_Mod, ONLY : GrdState
    USE State_Diag_Mod, ONLY : DgnState
    USE TIME_MOD,       ONLY : GET_LOCALTIME

    use parkind,        only : im => kind_im, rb => kind_rb

    USE IsorropiaII_Main_Mod, ONLY : Isorropia
    use module_mosaic_therm, only: mosaic
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
    TYPE(GrdState), INTENT(IN)    :: State_Grid  ! Grid State object
    TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
    TYPE(DgnState), INTENT(INOUT) :: State_Diag  ! Diagnostics State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REVISION HISTORY:
!  8/2008 - 10/2010 - F. Yu & G. Luo  - Initial versions
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
    real(kind=8), parameter :: cpdair = 1.004e3  ! Specific heat capacity
                                                 ! of dry air at constant
                                                 ! pressure at 273 K
                                                 ! (J kg-1 K-1)
!
! !LOCAL VARIABLES:
!
    INTEGER :: I,J,L,N,M,SIZENUM,MDAY

    REAL*8  :: PRESS,TK,RHIN,XQ,CACID,PACID,DT,DTAPM
    REAL*8  :: MSO4,MSO4BULK,MNIT,MNH4,MMSA,SOAT
    REAL*8  :: MBCS, MOCS   ! mass of sulfate attached to BC, OC
    REAL*8  :: MDSTS, MSALTS   ! mass of sulfate attached to dust,sea salt
    REAL*8  :: MSULFT   ! total sulfate
    REAL*8  :: MBC(NBCOC),MOC(NBCOC)
    REAL*8  :: XM1D(NSO4+NSEA),XN1D(NSO4),XNOLD,TEMPOUT1(NTEMPOUT1)
    REAL*8  :: ATOM4N(4)
    REAL*8  :: XMDST(NDSTB)

    REAL*8  :: MASS1, MASS2

    REAL*8  :: VOL, CLVSOG, PLVSOG01,PLVSOG1
    REAL*8  :: MSULFLV,MBCLV,MOCLV,MDSTLV,MSALTLV

    REAL*8  :: XOH,XSINK,XU,XV
    INTEGER :: KYEAR,KMON,KDAY,KHOUR,KMIN,ISITE,JSITE,NSITE
    REAL*8  :: TOP, TOPP
    INTEGER :: KKOUT
    REAL*8  :: XLAT, XLON

    REAL*8  :: CSO2,CNH3,XN0,CAMINE(3),CAMINEEMIT(3),YAMINEEMIT(3)
    REAL*8  :: NH3EMIT
    REAL*8  :: CSOG(14),CSOA(14)
    REAL*8  :: GFTOT1,GFTOT2,DENWET1,DENWET2
    INTEGER :: IACT10, IACT20, IACT30   ! bin index for cloud act
                                          ! diameters corresponding to RDRY
    INTEGER :: RACT1, RACT2, RACT3
    REAL*8  :: FCLOUD1(NSO4+4)
    REAL*8  :: AERAREA1(NSO4+NSEA+NDSTB+4)
    REAL*8  :: AERDRYR1(NSO4+NSEA+NDSTB+4)
    REAL*8  :: GAMMAPM1(NTYP)
    INTEGER :: NCOAG1,NCOAG2
    INTEGER, PARAMETER :: IFSITEADD = 1

    !REAL*8  :: SITEOUT(MAXSITE,LOUT,35),YTOP(MAXSITE,2)
    REAL*8  :: SITEOUT(MAXSITE,LOUT,59),YTOP(MAXSITE,2)
    REAL*8  :: SITEOUT1(MAXSITE,LOUT,NSO4+NSEA+NDSTB+NBCOC+NBCOC)
    REAL*8  :: SITEOUT2(MAXSITE,LOUT,NTEMPOUT1+6)
    REAL*8  :: SITEOUT3(MAXSITE,LOUT,11)
    REAL*8  :: AEROCOMOUT(MAXSITE,50,27),AEROCOMOUT1D(27)

    REAL*8  :: CCO, CNO,CNO2,CNO3,CHNO3,CISOP,CMTPA

    REAL*8  :: YSPGF,XBCLIFE,XOCLIFE  !Yu+   6/2/11
    REAL*8  :: YCCN,YCDN,YCDNSP,YCLDF,YCLDLIQ,YCLDICE,YRCLDL,VZ,YF,YC
    REAL*8  :: XCDN,XCDNSP
    REAL*8  :: wbar,relhum,yqc,yna(3),YB,YREI,YK
    !REAL*8  :: nuci, onihf, oniimm, onidep, onimey
    REAL*8  :: PRESS0, YSIGMA
    REAL*8  :: ACS,XAMINE,SCOS,LOCALTIME
    REAL*8  :: dumc, dumnc, pgam, lamc
    REAL(kind=8)  :: CCLD(State_Grid%NZ)      !!! CLOUD COVER
    REAL(kind=8)  :: CLDLIQ(State_Grid%NZ)    !!! CLOUD LIQUID WATER CONTENT
    REAL(kind=8)  :: CLDICE(State_Grid%NZ)    !!! CLOUD ICE WATER CONTENT
    REAL(kind=8)  :: REL(State_Grid%NZ)
    REAL(kind=8)  :: REI(State_Grid%NZ)
    REAL(kind=8)  :: taucloud(State_Grid%NZ,29)
    REAL(kind=8)  :: taucloudl(State_Grid%NZ,29)
    REAL(kind=8)  :: taucloudi(State_Grid%NZ,29)
    REAL(kind=8)  :: ssacloudl(State_Grid%NZ,29)
    REAL(kind=8)  :: ssacloudi(State_Grid%NZ,29)

    REAL*8  :: ERATIO(3)      ! ratio of NH3 emission to amine emission
    REAL*8  :: MWAMINE(3)
    REAL*8  :: OXRATE(3), XX0, XX1,DXX,XAREA,XCSNH3

    LOGICAL, SAVE    :: FIRST = .TRUE.
    LOGICAL, SAVE    :: FIRST1 = .TRUE. !mxy+
    LOGICAL, SAVE    :: FIRSTICE = .TRUE.
    LOGICAL, SAVE    :: FIRSTCOD = .TRUE.
    CHARACTER(LEN=4) :: YEAROUT='2000'
    CHARACTER(LEN=2) :: MONOUT='01'

    ! Make a pointer to the tracer array
    REAL*8, POINTER :: Spc(:,:,:,:)

    INTEGER :: vbs_nbin(1)=0
    REAL :: p_atm,cairclmdum,gas(4),aer(14),aerh2o,XM,VRATIO
    REAL :: gasold(4),aerold(14)
    REAL :: t_k4,rh4,dtchem,aersize

    ! Array dimensions
    INTEGER, PARAMETER       :: NOTHERA  =  9
    INTEGER, PARAMETER       :: NCTRLA   =  2
    INTEGER, PARAMETER       :: NCOMPA   =  8
    INTEGER, PARAMETER       :: NIONSA   = 10
    INTEGER, PARAMETER       :: NGASAQA  =  3
    INTEGER, PARAMETER       :: NSLDSA   = 19

    ! Concentration lower limit [mole/m3]
    REAL(fp),  PARAMETER       :: CONMIN = 1.0e-30_fp
!
! !LOCAL VARIABLES:
!
    REAL(fp)                   :: ANO3, GNO3
    REAL(f8)                   :: RHI, TEMPI
    REAL(fp)                   :: TCA,  TMG
    REAL(fp)                   :: TNA,  TCL,  TNH3, TNH4
    REAL(fp)                   :: TNIT, TNO3, TSO4
    REAL(f8)                   :: AERLIQ(NIONSA+NGASAQA+2)
    REAL(f8)                   :: AERSLD(NSLDSA)
    REAL(f8)                   :: GAS1(NGASAQA)
    REAL(f8)                   :: OTHER(NOTHERA)
    REAL(f8)                   :: WI(NCOMPA)
    REAL(f8)                   :: WT(NCOMPA)
    REAL(f8)                   :: CNTRL(NCTRLA)
    CHARACTER(LEN=255)       :: X
    CHARACTER(LEN=15)        :: SCASI
    REAL*8                   :: TSO4COAT,DNH3MAX
    REAL*8                   :: TH2O,TINORG   !Yu+ 6/1/11
    !--------------------------------------------------------------------------
    ! These do not appear to be used anymore (bmy, 6/18/19)
    !REAL,SAVE                :: SEABIRDEM(State_Grid%NX,State_Grid%NY,2)
    !REAL                     :: BIRDNH3(State_Grid%NX,State_Grid%NY)
    !--------------------------------------------------------------------------

    ! NOTE: We have removed GRID2x25, GRID4x5 switches (bmy, 6/18/19)
    INTEGER,SAVE :: ATOMGRID(2,400),GRIDNUM
    INTEGER :: RECID,ATOMGRIDIN(2),fileinfo
    CHARACTER(LEN=12) :: ATOMDATE

    WRITE(6,*)'    - APM calculation '

    !GanLuo+ local time satellite outputs
    CALL GET_LOCAL_TIME( State_Grid )

    ! Point to Spc
    Spc => State_Chm%Species

    ! Chemistry timestep [s]
    DT = GET_TS_CHEM()
    DTAPM = DT  !s

    ! Compute diurnal scaling for OH
    CALL OHNO3TIME( State_Grid )

    ! Calculate 3-D ionization rate
    IF(IFQANN.EQ.1) THEN
       IF(FIRST) THEN   ! only calculate Q once
          WRITE(6,*)"CALCULATE 3-D ANNUAL MEAN IONIZATION RATES"
          CALL APM3DQGEOS( Input_Opt, State_Grid, State_Met, IY)
          FIRST = .FALSE.
       ENDIF
    ELSE
       ! APM2++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
       IF(FIRST.or.ITS_A_NEW_MONTH()) THEN   ! only calculate Q once a month
          !IY = 1    ! Solar min, max Q
          !IY = -1   ! Solar max, min Q
          WRITE(6,*)"CALCULATE 3-D IONIZATION RATES", IY
          CALL APM3DQ( IY, Input_Opt, State_Grid, State_Met )
          FIRST = .FALSE.
       ENDIF
       ! APM2++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ENDIF

    ! TEMPOUT
    IF(FIRST1)THEN
       IF(IFATOM>0.5)THEN

          IF ( State_Grid%NX == 144 .and. State_Grid%NY == 91 ) THEN
             ATOMGRID=-999
             GRIDNUM=0
             CLOSE(100)
             OPEN(100,FILE='flightAtom-geos225.txt')
             DO
                READ(100,*,IOSTAT=fileinfo)(ATOMGRIDIN(I),I=1,2)
                if(fileinfo<0)exit
                GRIDNUM=GRIDNUM+1
                ATOMGRID(:,GRIDNUM)=ATOMGRIDIN
             ENDDO
             write(*,*)'ATOM2x2.5 outputs=',GRIDNUM

          ELSE IF ( State_Grid%NX == 72 .and. State_Grid%NY == 41 ) THEN
             ATOMGRID=-999
             GRIDNUM=0
             CLOSE(100)
             OPEN(100,FILE='flightAtom-geos45.txt')
             DO
                READ(100,*,IOSTAT=fileinfo)(ATOMGRIDIN(I),I=1,2)
                if(fileinfo<0)exit
                GRIDNUM=GRIDNUM+1
                ATOMGRID(:,GRIDNUM)=ATOMGRIDIN
             ENDDO
             write(*,*)'ATOM4x5 outputs=',GRIDNUM
          ELSE
             WRITE(*,*)'Can only run ATom with 2x2.5 or 4x5'
             STOP
          ENDIF
       ENDIF

       CALL  READCK6DTABLE

       ! Intitialize for first step useage
       GFTOT3D = 1.0
       DENWET3D = 2.0
       MWSIZE3D = 0.D0

       ! init FCLOUD
       FCLOUD=0.D0
       DO N=26,NSO4
       DO L=1,State_Grid%NZ
       DO J=1,State_Grid%NY
       DO I=1,State_Grid%NX
          FCLOUD(I,J,L,N)=1/15.D0
       ENDDO
       ENDDO
       ENDDO
       ENDDO

       ! APM2++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
       ! Get  (I, J) of selected sites
       IF(IFSITE.EQ.1) THEN
          CALL YSITESIJ_GL( State_Grid )             ! global sites
       ELSEIF(IFSITE.EQ.7.or.IFSITE.EQ.73) THEN
          CALL YSITESIJ_GL1( State_Grid )            ! global sites 1
       ELSEIF(IFSITE.EQ.9) THEN
          CALL YSITESIJ_BACCHUS( State_Grid )        ! for CCN intercomparison
       ELSEIF(IFSITE.EQ.2) THEN
          CALL YSITESIJ_EU( State_Grid )             ! nest EU sites
       ELSEIF(IFSITE.EQ.3) THEN
          CALL YSITESIJ_AC( State_Grid )             ! global AEROCOM sites
       ELSEIF(IFSITE.EQ.4) THEN
          CALL YSITESIJ_CH( State_Grid )             ! nest CH east asia sites
       ELSEIF(IFSITE.EQ.41) THEN
          CALL YSITESIJ_CH1( State_Grid )            ! nest CH east asia sites
       ELSEIF(IFSITE.EQ.5) THEN
          CALL YSITESIJ_NA( State_Grid )             ! nest NA sites
       ELSEIF(IFSITE.EQ.6) THEN
          CALL YSITESIJ_AOD_CN( State_Grid )         ! AOD_CN sites
       ENDIF
       ! APM2++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

       !define the variables you want to output
       IFTEMPOUT=.TRUE.
       IF(IFTEMPOUT)THEN
          NTEMPOUT=NTEMPOUT1 !number of species you want to output
          ALLOCATE( TEMPOUT( State_Grid%NX, State_Grid%NY, &
                             State_Grid%NZ, NTEMPOUT       ) )
          TEMPOUT = 0d0
          NPOUTSTEPS = 0
          FIRST1=.FALSE.
       ENDIF

    ENDIF

    MWAMINE(1) = 31  !  MA  -- CH3NH2 -- 31 g/mol
    MWAMINE(2) = 45. !  DMA -- (CH3)2NH -- 45 g/mol
    MWAMINE(3) = 59. !  TMA -- (CH3)3N -- 59 g/mol

    ERATIO(1) = 83./50000.
    ERATIO(2) = 33./50000.
    ERATIO(3) = 169./50000.

    SITEOUT   = 0.
    SITEOUT1  = 0.
    SITEOUT2  = 0.

    ZBEXT3D   = 0.
    ZW3D      = 1.
    ZG3D      = 0.
    YBEXT3D   = 0.
    XBEXT1k3D = 0.
    YW3D      = 1.
    YG3D      = 0.

    IF(FIRSTICE.or.ITS_A_NEW_MONTH()) THEN   ! only calculate Q once a month
       KYEAR = GET_YEAR()
       KMON = GET_MONTH()
       WRITE(YEAROUT,'(I4.4)')MIN(2015,MAX(2004,KYEAR))
       WRITE(MONOUT,'(I2.2)')KMON

       !Luo CLOSE(121)
       !Luo#if defined( GRID4x5 )
       !Luo OPEN(121,FILE='geosdata/GEOS_4x5/MODISIER/'//       &
       !Luo               'MODG45IER'//YEAROUT//MONOUT//'.bin', &
       !Luo               ACCESS='DIRECT',FORM='UNFORMATTED',   &
       !Luo               RECL=State_Grid%NX*State_Grid%NY)
       !Luo READ(121,REC=1)RCLDI2D
       !Luo#endif
       !Luo#if defined( GRID2x25 )
       !Luo OPEN(121,FILE='geosdata/GEOS_2x2.5/MODISIER/'//     &
       !Luo               'MODG45IER'//YEAROUT//MONOUT//'.bin', &
       !Luo               ACCESS='DIRECT',FORM='UNFORMATTED',   &
       !Luo               RECL=State_Grid%NX*State_Grid%NY)
       !Luo READ(121,REC=1)RCLDI2D
       !Luo#endif

       FIRSTICE=.FALSE.
    ENDIF

    IF(FIRSTCOD)THEN
       call rrtmg_sw_ini(cpdair)
       WRITE(6,*) "run rrtmg_sw_ini"
       FIRSTCOD = .FALSE.
    ENDIF

    !Doing MESA
    !$OMP PARALLEL DO        &
    !$OMP DEFAULT( SHARED )  &
    !$OMP PRIVATE( I, J, L ) &
    !$OMP SCHEDULE( DYNAMIC )
    DO L = 1, State_Grid%NZ
    DO J = 1, State_Grid%NY
    DO I = 1, State_Grid%NX
       UPTAKE(I,J,L,1)=Spc(I,J,L,APMIDS%id_NH3)
       SPGF(I,J,L) = 1.d0
       MASSISRP(I,J,L,1)=Spc(I,J,L,APMIDS%id_NH3)
       MASSISRP(I,J,L,2)=Spc(I,J,L,APMIDS%id_NH4)
       MASSISRP(I,J,L,3)=Spc(I,J,L,APMIDS%id_HNO3)
       MASSISRP(I,J,L,4)=Spc(I,J,L,APMIDS%id_NIT)

       MASSMESA(I,J,L,1)=Spc(I,J,L,APMIDS%id_NH3)
       MASSMESA(I,J,L,2)=Spc(I,J,L,APMIDS%id_NH4)
       MASSMESA(I,J,L,3)=Spc(I,J,L,APMIDS%id_HNO3)
       MASSMESA(I,J,L,4)=Spc(I,J,L,APMIDS%id_NIT)
    ENDDO
    ENDDO
    ENDDO
    !$OMP END PARALLEL DO

    if(IFDOISRP==1.or.IFDOISRP==2.or.IFDOISRP==3.or.IFDOISRP==4)then
       !******* DO ISORROPIA *******
       !$OMP PARALLEL DO                                                      &
       !$OMP DEFAULT( SHARED )                                                &
       !$OMP PRIVATE( I,    J,      L,       N,      WI,   WT,  GAS1,  TEMPI )&
       !$OMP PRIVATE( RHI,  VOL,    TSO4,    TNH3,   TNA,  TCL, ANO3, GNO3  ) &
       !$OMP PRIVATE( TCA,  TMG,    TK,      CNTRL,  SCASI                  ) &
       !$OMP PRIVATE( TNO3, AERLIQ, AERSLD,  OTHER,  TNH4, TNIT             ) &
       !$OMP PRIVATE( TSO4COAT ,DNH3MAX            )                          &
       !$OMP PRIVATE( TH2O, XM,VRATIO)                                        &
       !$OMP SCHEDULE( DYNAMIC )
       DO L = 1, State_Grid%NZ
       DO J = 1, State_Grid%NY
       DO I = 1, State_Grid%NX

          ! Temperature [K]
          TEMPI = MAX(235.D0,State_Met%T(I,J,L))

          if((Spc(I,J,L,APMIDS%id_NH3)  + &
              Spc(I,J,L,APMIDS%id_NH4)  + &
              Spc(I,J,L,APMIDS%id_HNO3) + &
              Spc(I,J,L,APMIDS%id_NIT)) > 1.d-12)then

             ! Initialize WI, WT
             DO N = 1, NCOMPA
                WI(N) = 0e+0_fp
                WT(N) = 0e+0_fp
             ENDDO

             ! Initialize GAS
             DO N = 1, NGASAQA
                GAS1(N) = 0e+0_fp
             ENDDO

             ! Relative humidity [unitless]
             RHI      = State_Met%RH(I,J,L) * 1.e-2_fp

             ! Force RH in the range 0.01 - 0.98
             RHI      = MAX( 0.01e+0_fp, RHI )
             RHI      = MIN( 0.98e+0_fp, RHI )

             ! Volume of grid box [m3]
             VOL      = State_Met%AIRVOL(I,J,L)

             !GanLuo+ Deduce SULFLO form TSO4,
             ! Add those sulfate coated on primary particles
             !GanLuo mv LVSOG         TSO4     = (SUM(Spc(I,J,L,APMIDS%id_SO4BIN1:(APMIDS%id_SO4BIN1+NSO4-1)))
             !GanLuo mv LVSOG     &            ) * 1.e+3_fp / ( 96.e+0_fp * VOL )
             TSO4 = (SUM( &
                    Spc(I,J,L,APMIDS%id_SO4BIN1:(APMIDS%id_SO4BIN1+NSO4-1))) &
                    - Spc(I,J,L,APMIDS%id_CTSO4))*1.e+3_fp / ( 96.e+0_fp * VOL )

             TSO4COAT = (Spc(I,J,L,APMIDS%id_CTBC) + &
                         Spc(I,J,L,APMIDS%id_CTOC) + &
                         Spc(I,J,L,APMIDS%id_CTDST)+ &
                         Spc(I,J,L,APMIDS%id_CTSEA)) * &
                         1.e+3_fp / ( 96.e+0_fp * VOL )
             TSO4 = TSO4 + TSO4COAT

             ! Total NH3 [mole/m3]
             IF(IFDOISRP==4)THEN
                DNH3MAX = Spc(I,J,L,APMIDS%id_NH3)
             ELSE
                DNH3MAX = Spc(I,J,L,APMIDS%id_NH3) * (1.-exp(-DT/TAONH3(I,J,L)))
             ENDIF

             TNH3 = Spc(I,J,L,APMIDS%id_NH4) * 1.e+3_fp / &
                    ( 18.e+0_fp * VOL ) + &
                    DNH3MAX * 1.e+3_fp / ( 17.e+0_fp * VOL )
             MASSISRP(I,J,L,1) = Spc(I,J,L,APMIDS%id_NH3)-DNH3MAX !remaining NH3

             ! Total Na+ (30.61% by weight of seasalt) [mole/m3]
             !TNA = SUM(Spc(I,J,L,APMIDS%id_SEABIN1:(APMIDS%id_SEABIN1+NSEA-1))) &
             !      * 0.3061e+0_fp * 1.e+3_fp /( 22.99e+0_fp  * VOL  )

             ! Total Cl- (55.04% by weight of seasalt) [mole/m3]
             !TCL = SUM(Spc(I,J,L,APMIDS%id_SEABIN1:(APMIDS%id_SEABIN1+NSEA-1))) &
             !      * 0.5504e+0_fp * 1.e+3_fp /( 22.99e+0_fp  * VOL  )

             !GLuo: Sea salt in ISORROPIA needs to be updated.
             TNA      = 0e+0_fp
             TCL      = 0e+0_fp

             TCA      = 0e+0_fp
             TK       = 0e+0_fp
             TMG      = 0e+0_fp

             !---------------------
             ! COUPLED SIMULATION
             !---------------------

             ! Compute gas-phase HNO3 [mole/m3] from HNO3 tracer
             GNO3  = Spc(I,J,L,APMIDS%id_HNO3)
             GNO3  = MAX( GNO3 * 1.e+3_fp / ( 63.e+0_fp * VOL ), CONMIN )

             ! Aerosol-phase NO3 [mole/m3]
             ANO3  = Spc(I,J,L,APMIDS%id_NIT) * 1.e+3_fp / ( 62.e+0_fp * VOL )

             ! Total NO3 [mole/m3]
             TNO3    = GNO3 + ANO3

             !---------------------------------
             ! Call ISORROPIA
             !---------------------------------

             ! set type of ISORROPIA call
             ! Forward problem, do not change this value
             ! 0e+0_fp represents forward problem
             CNTRL(1) = 0.0e+0_fp

             ! Metastable for now
             ! 1e+0_fp represents metastable problem
             CNTRL(2) = 1.0e+0_fp

             ! Insert concentrations [mole/m3] into WI & prevent underflow
             WI(1)    = MAX( TNA,  CONMIN )
             WI(2)    = MAX( TSO4, CONMIN )
             WI(3)    = MAX( TNH3, CONMIN )
             WI(4)    = MAX( TNO3, CONMIN )
             WI(5)    = MAX( TCL,  CONMIN )
             WI(6)    = MAX( TCA,  CONMIN )
             WI(7)    = MAX( TK,   CONMIN )
             WI(8)    = MAX( TMG,  CONMIN )

             ! Perform aerosol thermodynamic equilibrium
             ! ISORROPIA can be found in ISORROPIAIICODE.f
             ! inputs are WI, RHI, TEMPI, CNTRL
             CALL ISORROPIA( WI,    RHI,  TEMPI,  CNTRL,   &
                             WT,    GAS1,  AERLIQ, AERSLD, &
                             SCASI, OTHER)

             IF(IFDOISRP==1.or.IFDOISRP==3.or.IFDOISRP==4)THEN
                TH2O = AERLIQ(8)  ! aerosol water  (mole/m3)

                XM = WT(2)*96. + (WT(3)-GAS1(1))*18. +(WT(4)-GAS1(2))*62. &
                   + TNA*23. + TCL*35.45 !SO4, NH4,NIT,NACL (g/m3)
                XM = MAX(XM,CONMIN)
                VRATIO = 1. + (TH2O*(18.*1.0))/(XM/1.8) ! assume salt has density of 1.8 g/cm3

                SPGF(I,J,L) = VRATIO**(1./3.0)
             ENDIF

             TNH3 = MAX( 17.e-3_fp * VOL *   GAS1(1),           CONMIN )
             TNH4 = MAX( 18.e-3_fp * VOL * ( WT(3) - GAS1(1) ), CONMIN )
             TNIT = MAX( 62.e-3_fp * VOL * ( WT(4) - GAS1(2) ), CONMIN )

             ! Save tracers back into Spc array [kg]
             ! no longer save TSO4 back into Spc. SO4 is all aerosol phase
             ! (hotp 11/7/07)
             ! Spc(I,J,L,APMIDS%id_SO4) = TSO4
             MASSISRP(I,J,L,1) = MASSISRP(I,J,L,1) + TNH3
             MASSISRP(I,J,L,2) = TNH4
             MASSISRP(I,J,L,3) = MAX( 63.e-3_fp * VOL * GAS1(2), CONMIN )
             MASSISRP(I,J,L,4) = TNIT
             UPTAKE(I,J,L,1)=MASSISRP(I,J,L,1)-UPTAKE(I,J,L,1)  ! in kg

          endif

       ENDDO
       ENDDO
       ENDDO
       !$OMP END PARALLEL DO
    endif

    !******* MESA *******
    if(IFDOISRP==0.or.IFDOISRP==2.or.IFDOISRP==3)then

       DO L = 1, State_Grid%NZ
       DO J = 1, State_Grid%NY
       DO I = 1, State_Grid%NX

          if((Spc(I,J,L,APMIDS%id_NH3)  + &
              Spc(I,J,L,APMIDS%id_NH4)  + &
              Spc(I,J,L,APMIDS%id_HNO3) + &
              Spc(I,J,L,APMIDS%id_NIT)) > 1.d-16)then

             aersize=0.d0
             IF((SUM(Spc(I,J,L,APMIDS%id_SO4BIN1:(APMIDS%id_SO4BIN1+NSO4-1)))+ &
                 SUM(Spc(I,J,L,APMIDS%id_SEABIN1:(APMIDS%id_SEABIN1+NSEA-1)))) &
                 > 1.d-30)THEN
                DO N=1,NSO4
                   SIZENUM=APMIDS%id_SO4BIN1+N-1
                   aersize=aersize+Spc(I,J,L,SIZENUM)*RDRY(N)*GFTOT3D(I,J,L,1)
                ENDDO
                DO N=1,NSEA
                   SIZENUM=APMIDS%id_SEABIN1+N-1
                   aersize=aersize+Spc(I,J,L,SIZENUM)*RSALT(N)*GFTOT3D(I,J,L,2)
                ENDDO
                aersize=aersize/ &
                 (SUM(Spc(I,J,L,APMIDS%id_SO4BIN1:(APMIDS%id_SO4BIN1+NSO4-1)))+&
                  SUM(Spc(I,J,L,APMIDS%id_SEABIN1:(APMIDS%id_SEABIN1+NSEA-1))))
             ELSE
                aersize=50.D-9
             ENDIF

             t_k4  = max(235.D0,real(State_Met%T(I,J,L))) ! Temperature [K]
             PRESS = GET_PCENTER(I,J,L) * 1.d2 ! P at level center [Pa]
             p_atm = real(PRESS/101325.0)
             VOL   = State_Met%AIRVOL(I,J,L) ! Volume of grid box [m3]
             RHIN  = MIN( 98.d0, State_Met%RH(I,J,L) ) ! Cap RH at 99%
             RHIN  = MAX( 1.d0, RHIN ) ! Safety check
             cairclmdum=real(State_Met%AIRDEN(I,J,L)*1.e-6/MAIR)
             gas=0.
             aer=0.

             gas(2) = real(Spc(I,J,L,APMIDS%id_NH3)*MAIR/ &
                          (State_Met%AD(I,J,L)*17.0e-3))!NH3 mol/mol air
             gas(3) = real(Spc(I,J,L,APMIDS%id_HNO3)*MAIR/ &
                          (State_Met%AD(I,J,L)*63.0e-3))!HNO3

             aer(1) = real(((MAX(1.e-30, &
               (SUM(Spc(I,J,L,APMIDS%id_SO4BIN1:(APMIDS%id_SO4BIN1+NSO4-1))) &
                   -Spc(I,J,L,APMIDS%id_CTSO4)) )) &
                   +Spc(I,J,L,(APMIDS%id_CTBC))    &
                   +Spc(I,J,L,(APMIDS%id_CTOC))    &
                   +Spc(I,J,L,APMIDS%id_CTDST)+Spc(I,J,L,APMIDS%id_CTSEA))*MAIR&
                   /(State_Met%AD(I,J,L)*96.0e-3)) !SO4
             aer(2) = real(Spc(I,J,L,APMIDS%id_NH4)*MAIR/ &
                          (State_Met%AD(I,J,L)*18.0e-3))!NH4
             aer(3) = real(Spc(I,J,L,APMIDS%id_NIT)*MAIR/ &
                          (State_Met%AD(I,J,L)*62.0e-3))!NIT

             aer(4) = real(SUM( &
                      Spc(I,J,L,APMIDS%id_SEABIN1:(APMIDS%id_SEABIN1+NSEA-1))) &
                      *MAIR*0.3061/(State_Met%AD(I,J,L)*22.99e-3))!NA
             aer(5) = real(SUM( &
                      Spc(I,J,L,APMIDS%id_SEABIN1:(APMIDS%id_SEABIN1+NSEA-1))) &
                      *MAIR*0.5504/(State_Met%AD(I,J,L)*35.45e-3))!CL

             gasold=gas
             aerold=aer

             rh4=real(RHIN)

             aerh2o=0.D0
             TEMPOUT1(NTEMPOUT1-22)=aersize

             CALL mosaic( 1, 1, DTAPM, vbs_nbin, &
                          t_k4, p_atm, rh4, cairclmdum, &
                          gas, aer, aerh2o, aersize )

             MASSMESA(I,J,L,1)=MIN((Spc(I,J,L,APMIDS%id_NH3)+ &
                                    Spc(I,J,L,APMIDS%id_NH4)), &
                                    (gas(2)*(State_Met%AD(I,J,L)*17.0d-3)/MAIR))
             MASSMESA(I,J,L,2)=Spc(I,J,L,APMIDS%id_NH3)+ &
                               Spc(I,J,L,APMIDS%id_NH4)-MASSMESA(I,J,L,1)
             MASSMESA(I,J,L,3)=MIN((Spc(I,J,L,APMIDS%id_HNO3)+ &
                                    Spc(I,J,L,APMIDS%id_NIT)), &
                                   (gas(3)*(State_Met%AD(I,J,L)*63.0d-3)/MAIR))
             MASSMESA(I,J,L,4)=Spc(I,J,L,APMIDS%id_HNO3)+ &
                               Spc(I,J,L,APMIDS%id_NIT)-MASSMESA(I,J,L,3)

             !MASSMESA(I,J,L,1)=gas(2)*(State_Met%AD(I,J,L)*17.0e-3)/MAIR
             !MASSMESA(I,J,L,2)=aer(2)*(State_Met%AD(I,J,L)*18.0e-3)/MAIR
             !MASSMESA(I,J,L,3)=gas(3)*(State_Met%AD(I,J,L)*63.0e-3)/MAIR
             !MASSMESA(I,J,L,4)=aer(3)*(State_Met%AD(I,J,L)*62.0e-3)/MAIR

             IF(IFDOISRP==0.or.IFDOISRP==2)THEN
                XM = aer(1)*(State_Met%AD(I,J,L)*96.0)/MAIR/VOL &
                     + MASSMESA(I,J,L,2)*1.d3/VOL &
                     + MASSMESA(I,J,L,4)*1.d3/VOL !g/m3
                XM = MAX(XM,1.d-30)
                VRATIO = 1. + (aerh2o*1.d3)/(XM/1.8)   ! assume salt has density of 1.8 g/cm3

                SPGF(I,J,L) = VRATIO**(1./3.0)
             ENDIF

          endif
       ENDDO
       ENDDO
       ENDDO
    endif

    IF(IFDOISRP==1.or.IFDOISRP==4.or.IFDOISRP==5)THEN
       !$OMP PARALLEL DO       &
       !$OMP DEFAULT( SHARED ) &
       !$OMP PRIVATE( I, J, L) &
       !$OMP SCHEDULE( DYNAMIC )
       DO L = 1, State_Grid%NZ
       DO J = 1, State_Grid%NY
       DO I = 1, State_Grid%NX
          Spc(I,J,L,APMIDS%id_NH3) = MASSISRP(I,J,L,1)
          Spc(I,J,L,APMIDS%id_NH4) = MASSISRP(I,J,L,2)
          Spc(I,J,L,APMIDS%id_HNO3)= MASSISRP(I,J,L,3)
          Spc(I,J,L,APMIDS%id_NIT) = MASSISRP(I,J,L,4)
       ENDDO
       ENDDO
       ENDDO
       !$OMP END PARALLEL DO
    ENDIF

    if(IFDOISRP==0.or.IFDOISRP==2.or.IFDOISRP==3)then
       !$OMP PARALLEL DO       &
       !$OMP DEFAULT( SHARED ) &
       !$OMP PRIVATE( I, J, L) &
       !$OMP SCHEDULE( DYNAMIC )
       DO L = 1, State_Grid%NZ
       DO J = 1, State_Grid%NY
       DO I = 1, State_Grid%NX
          Spc(I,J,L,APMIDS%id_NH3) = MASSMESA(I,J,L,1)
          Spc(I,J,L,APMIDS%id_NH4) = MASSMESA(I,J,L,2)
          Spc(I,J,L,APMIDS%id_HNO3)= MASSMESA(I,J,L,3)
          Spc(I,J,L,APMIDS%id_NIT) = MASSMESA(I,J,L,4)
       ENDDO
       ENDDO
       ENDDO
       !$OMP END PARALLEL DO
    ENDIF

    !$OMP PARALLEL DO                                                       &
    !$OMP DEFAULT( SHARED )                                                 &
    !$OMP PRIVATE( I, J, L, N )                                             &
    !$OMP PRIVATE( SIZENUM, PRESS, TK, RHIN )                               &
    !$OMP PRIVATE( CACID,PACID )                                            &
    !$OMP PRIVATE( MSO4,MSO4BULK,MNIT,MNH4,SOAT)                            &
    !$OMP PRIVATE( MBCS, MOCS, MSULFT,MDSTS,MSALTS)                         &
    !$OMP PRIVATE( MBC, MOC, MMSA)                                          &
    !$OMP PRIVATE( XMDST)                                                   &
    !$OMP PRIVATE( MASS1, MASS2)                                            &
    !$OMP PRIVATE( CSO2,CNH3,XN0,CAMINE,CAMINEEMIT,YAMINEEMIT)              &
    !$OMP PRIVATE( CCO,CNO,CNO2,CNO3,CHNO3,CISOP,CMTPA)                     &
    !$OMP PRIVATE( NH3EMIT)                                                 &
    !$OMP PRIVATE( CSOG)                                                    &
    !$OMP PRIVATE( CSOA)                                                    &
    !$OMP PRIVATE( VOL)                                                     &
    !$OMP PRIVATE( CLVSOG,MSULFLV,MBCLV,MOCLV,MDSTLV,MSALTLV)               &
    !$OMP PRIVATE( XM1D,XN1D,XNOLD,TEMPOUT1,ATOM4N,AEROCOMOUT1D)            &
    !$OMP PRIVATE( XQ,PLVSOG01,PLVSOG1,GFTOT1,GFTOT2,DENWET1,DENWET2)       &
    !$OMP PRIVATE( IACT10,IACT20,IACT30,FCLOUD1,AERAREA1,AERDRYR1,GAMMAPM1) &
    !$OMP PRIVATE( RACT1,RACT2,RACT3)                                       &
    !$OMP PRIVATE( NCOAG1,NCOAG2)                                           &
    !$OMP PRIVATE( YSPGF,XBCLIFE,XOCLIFE,XCSNH3)                            &!Yu+
    !$OMP PRIVATE( XOH, XSINK,XAREA,XX0,XX1,DXX,ACS,XLAT, XLON,XAMINE)      &
    !$OMP PRIVATE( KYEAR,KMON,KDAY,KHOUR,KMIN,ISITE,JSITE,NSITE)            &
    !$OMP PRIVATE( XU,XV,TOP, TOPP)                                         &
    !$OMP PRIVATE( KKOUT)                                                   &
    !$OMP PRIVATE( ZBEXT,ZW,ZG)                                             &!OPT+
    !$OMP PRIVATE( ZBABS)                                                   &!OPT+
    !$OMP PRIVATE( YBEXT,XBEXT1k,YW,YG)                                     &!OPT+
    !$OMP PRIVATE( IWL)                                                     &!OPT+
    !$OMP PRIVATE( ITYP)                                                    &!mxy+
    !$OMP PRIVATE( YCCN,YCDN,YCDNSP,YCLDF,YCLDLIQ,YCLDICE,YRCLDL,VZ)        &!Yu+ 7/2012
    !$OMP PRIVATE( XCDN,XCDNSP )                                            &
    !$OMP PRIVATE( YF,YC,SCOS,LOCALTIME )                                   &
    !$OMP PRIVATE( PRESS0, YSIGMA )                                         &
    !$OMP PRIVATE( wbar,relhum,yqc,yna,YB,YREI,YK)                          &
    !$OMP PRIVATE( dumc, dumnc, pgam, lamc)                                 &
    !$OMP PRIVATE( CCLD,CLDLIQ,CLDICE )                                     &
    !$OMP PRIVATE( REL,REI)                                                 &
    !$OMP PRIVATE( taucloud, taucloudl, taucloudi, ssacloudl, ssacloudi )   &
    !!$OMP PRIVATE( nuci, onihf, oniimm, onidep, onimey)                    &
    !$OMP SCHEDULE( DYNAMIC )
    DO J = 1, State_Grid%NY
    DO I = 1, State_Grid%NX

       YF = 0.
       YC = 1.
       PRESS0 = GET_PEDGE(I,J,1) * 1.d2    ! surafce P in Pa

       ! Local time: 0-23.999
       LOCALTIME = GET_LOCALTIME( I, J, 1, State_Grid )

       SCOS   = COS((LOCALTIME/12.-1.)*3.1415926)

       DO L = 1, State_Grid%NZ
          PRESS = GET_PCENTER(I,J,L) * 1.d2           ! P at level center [Pa]
          YSIGMA= PRESS/PRESS0
          TK    = State_Met%T(I,J,L)                  ! Temperature [K]
          VOL   = State_Met%AIRVOL(I,J,L)             ! Volume of grid box [m3]
          XAREA = 1.d-6*VOL/State_Met%BXHEIGHT(I,J,L) ! Area of grid box [km2]
          RHIN  = MIN( 99.d0, State_Met%RH(I,J,L) )   ! Cap RH at 99%
          RHIN  = MAX( 0.d0, RHIN )                   ! Safety check

          YSPGF = SPGF(I,J,L)  !SP GF from ISOROPIA/MESA

          XOH = GET_OH( I, J, L, APMIDS%id_OH, Input_Opt, State_Chm, State_Met )

          !IF(IFNUCL.EQ.1.or.IFNUCL.EQ.6) THEN
          XQ = XQ3D(I,J,L)
          !ELSE
          !   XQ = 1.d-20   !ion-pairs/cm3s
          !ENDIF

          !kg/(box*timestep) to #/(cm3*s)
          PACID = XMACID/96.* PSO4GAS(I,J,L)/(VOL*M1ACID*DT)*1.d-6

          !kg/box to #/cm3
          CACID = Spc(I,J,L,APMIDS%id_SO4G)/(VOL*M1ACID)*1.d-6

          XN0 = 1.d-6*PRESS/(1.3807d-23*TK)   !#/cm3
          CNH3=Spc(I,J,L,APMIDS%id_NH3) &
               /(VOL*State_Chm%SpcData(APMIDS%id_NH3)%Info%emMW_g &
               *1.d-3)*6.02d17   ! #/cm3
          CNH3=CNH3*1.d12/XN0   ! CNH3 now in ppt
          CSO2=Spc(I,J,L,APMIDS%id_SO2) &
               /(VOL*State_Chm%SpcData(APMIDS%id_SO2)%Info%emMW_g &
               *1.d-3)*6.02d17   ! #/cm3
          CSO2=CSO2*1.d12/XN0   ! CSO2 now in ppt

          CCO=Spc(I,J,L,APMIDS%id_CO) &
              /(VOL*State_Chm%SpcData(APMIDS%id_CO)%Info%emMW_g &
              *1.d-3)*6.02d17   ! #/cm3
          CCO=CCO*1.d6/XN0   ! CCO now in ppm

          CNO=Spc(I,J,L,APMIDS%id_NO) &
              /(VOL*State_Chm%SpcData(APMIDS%id_NO)%Info%emMW_g &
              *1.d-3)*6.02d17   ! #/cm3
          CNO=CNO*1.d9/XN0   ! CNO now in ppb
          CNO2=Spc(I,J,L,APMIDS%id_NO2) &
               /(VOL*State_Chm%SpcData(APMIDS%id_NO2)%Info%emMW_g &
               *1.d-3)*6.02d17   ! #/cm3
          CNO2=CNO2*1.d9/XN0   ! CNO2 now in ppb
          CNO3=Spc(I,J,L,APMIDS%id_NO3) &
               /(VOL*State_Chm%SpcData(APMIDS%id_NO3)%Info%emMW_g &
               *1.d-3)*6.02d17   ! #/cm3
          CNO3=CNO3*1.d9/XN0   ! CNO3 now in ppb
          CHNO3=Spc(I,J,L,APMIDS%id_HNO3) &
                /(VOL*State_Chm%SpcData(APMIDS%id_HNO3)%Info%emMW_g &
                *1.d-3)*6.02d17   ! #/cm3
          CHNO3=CHNO3*1.d9/XN0   ! CHNO3 now in ppb

          CISOP=Spc(I,J,L,APMIDS%id_ISOP) &
                /(VOL*State_Chm%SpcData(APMIDS%id_ISOP)%Info%emMW_g &
                *1.d-3)*6.02d17   ! #/cm3
          CISOP=CISOP*1.d9/XN0   ! CISOP now in ppb
          CMTPA=Spc(I,J,L,APMIDS%id_MTPA) &
                /(VOL*State_Chm%SpcData(APMIDS%id_MTPA)%Info%emMW_g &
                *1.d-3)*6.02d17   ! #/cm3
          CMTPA=CMTPA*1.d9/XN0   ! CMTPA now in ppb

          NH3EMIT = MAX(0.d0,EMITNH3(I,J,L))*14./17.  ! kg N/box-sec  (NH3 emission)

          ! Consider diurnal variation in emissions
          NH3EMIT = NH3EMIT * (1.+0.7*SCOS)

          IF(IFATHN>0)THEN
             DO ITYP = 1, 3
                CAMINE(ITYP)=Spc(I,J,L,APMIDS%id_AMINE-1+ITYP)/(VOL* &
                   State_Chm%SpcData(APMIDS%id_NH3)%Info%emMW_g*1.d-3* &
                   MWAMINE(ITYP)/17.)*6.02d17   ! #/cm3
                CAMINE(ITYP)=CAMINE(ITYP)*1.d12/XN0   ! CAMINE now in ppt

                YAMINEEMIT(ITYP) = MWAMINE(ITYP)/14.*NH3EMIT*ERATIO(ITYP) ! convert to kg-amine/box-sec
                CAMINEEMIT(ITYP) = YAMINEEMIT(ITYP)/(VOL* &
                   State_Chm%SpcData(APMIDS%id_NH3)%Info%emMW_g*1.d-3* &
                   MWAMINE(ITYP)/17.)*6.02d17   ! # cm-3s-1
             ENDDO

             IF(NUCLAMINE.EQ.1) THEN
                XAMINE = CAMINE(1)
             ELSEIF(NUCLAMINE.EQ.2) THEN
                XAMINE = CAMINE(2)
             ELSEIF(NUCLAMINE.EQ.3) THEN
                XAMINE = CAMINE(3)
             ELSEIF(NUCLAMINE.EQ.4) THEN
                XAMINE = CAMINE(1)+CAMINE(2)+CAMINE(3)
             ELSE
                !XAMINE = 0.
                XAMINE = CAMINE(2)  !JATHN calculated anyway, switch in phy to include or not
             ENDIF
          ENDIF

          IF(IFAG.EQ.1) THEN
             PLVSOG1 = sum(PLVSOG(I,J,L,1:5))/(M1LVSOG*1.d+15)  !#/cm3s, PLVSOG in ug/m3s, M1LVSOG in kg
             PLVSOG01 = PLVSOG(I,J,L,1)/(M1LVSOG*1.d+15)  !#/cm3s, PLVSOG in ug/m3s, M1LVSOG in kg
             CLVSOG  = MAX(1.D-30,Spc(I,J,L,APMIDS%id_SO4G+1)/ &
                           (VOL*M1LVSOG)*1.d-6)
             ! coatd LV-SOA
             MSULFLV = Spc(I,J,L,APMIDS%id_CTSO4  )/VOL !kg/m3
             MBCLV   = Spc(I,J,L,APMIDS%id_CTBC +1)/VOL !kg/m3
             MOCLV   = Spc(I,J,L,APMIDS%id_CTOC +1)/VOL !kg/m3
             MDSTLV  = Spc(I,J,L,APMIDS%id_CTDST+1)/VOL !kg/m3
             MSALTLV = Spc(I,J,L,APMIDS%id_CTSEA+1)/VOL !kg/m3
          ELSE
             PLVSOG1 = 0.
             PLVSOG01= 0.
             CLVSOG  = 0.
             MSULFLV = 0.
             MBCLV   = 0.
             MOCLV   = 0.
             MDSTLV  = 0.
             MSALTLV = 0.
          ENDIF

          IF(IFAG.EQ.1) THEN
             CSOG(1)=Spc(I,J,L,APMIDS%id_TSOG0)/ &
                (VOL*State_Chm%SpcData(APMIDS%id_TSOG0)%Info%emMW_g* &
                1.d-3)*6.02d17
             CSOG(2)=Spc(I,J,L,APMIDS%id_TSOG1)/ &
                (VOL*State_Chm%SpcData(APMIDS%id_TSOG1)%Info%emMW_g* &
                1.d-3)*6.02d17
             CSOG(3)=Spc(I,J,L,APMIDS%id_TSOG2)/ &
                (VOL*State_Chm%SpcData(APMIDS%id_TSOG2)%Info%emMW_g* &
                1.d-3)*6.02d17
             CSOG(4)=Spc(I,J,L,APMIDS%id_TSOG3)/ &
                (VOL*State_Chm%SpcData(APMIDS%id_TSOG3)%Info%emMW_g* &
                1.d-3)*6.02d17

             ! Isoprene SOA from Pye et al. (2010) has been removed;
             ! the mechanistic isoprene SOA from Marais et al. is used instead.
             ! Set to low value to avoid breaking APM. (mps, 7/15/19)
             CSOG(5) = 1.0e-30
             CSOG(6) = 1.0e-30
             CSOG(7) = 1.0e-30

             CSOG(8) =Spc(I,J,L,APMIDS%id_ASOG1)/ &
                (VOL*State_Chm%SpcData(APMIDS%id_ASOG1)%Info%emMW_g* &
                1.d-3)*6.02d17
             CSOG(9) =Spc(I,J,L,APMIDS%id_ASOG2)/ &
                (VOL*State_Chm%SpcData(APMIDS%id_ASOG2)%Info%emMW_g* &
                1.d-3)*6.02d17
             CSOG(10)=Spc(I,J,L,APMIDS%id_ASOG3)/ &
                (VOL*State_Chm%SpcData(APMIDS%id_ASOG3)%Info%emMW_g* &
                1.d-3)*6.02d17

             CSOG(11)=Spc(I,J,L,APMIDS%id_POG1)/ &
                (VOL*State_Chm%SpcData(APMIDS%id_POG1)%Info%emMW_g* &
                1.d-3)*6.02d17
             CSOG(12)=Spc(I,J,L,APMIDS%id_POG2)/ &
                (VOL*State_Chm%SpcData(APMIDS%id_POG2)%Info%emMW_g* &
                1.d-3)*6.02d17

             CSOG(13)=Spc(I,J,L,APMIDS%id_OPOG1)/ &
                (VOL*State_Chm%SpcData(APMIDS%id_OPOG1)%Info%emMW_g* &
                1.d-3)*6.02d17
             CSOG(14)=Spc(I,J,L,APMIDS%id_OPOG2)/ &
                (VOL*State_Chm%SpcData(APMIDS%id_OPOG2)%Info%emMW_g* &
                1.d-3)*6.02d17

             CSOA(1)=Spc(I,J,L,APMIDS%id_TSOA0)/ &
                (VOL*State_Chm%SpcData(APMIDS%id_TSOA0)%Info%emMW_g* &
                1.d-3)*6.02d17
             CSOA(2)=Spc(I,J,L,APMIDS%id_TSOA1)/ &
                (VOL*State_Chm%SpcData(APMIDS%id_TSOA1)%Info%emMW_g* &
                1.d-3)*6.02d17
             CSOA(3)=Spc(I,J,L,APMIDS%id_TSOA2)/ &
                (VOL*State_Chm%SpcData(APMIDS%id_TSOA2)%Info%emMW_g* &
                1.d-3)*6.02d17
             CSOA(4)=Spc(I,J,L,APMIDS%id_TSOA3)/ &
                (VOL*State_Chm%SpcData(APMIDS%id_TSOA3)%Info%emMW_g* &
                1.d-3)*6.02d17

             ! Isoprene SOA from Pye et al. (2010) has been removed;
             ! the mechanistic isoprene SOA from Marais et al. is used instead.
             ! Set to low value to avoid breaking APM. (mps, 7/15/19)
             CSOA(5) = 1.0e-30
             CSOA(6) = 1.0e-30
             CSOA(7) = 1.0e-30

             CSOA(8) =Spc(I,J,L,APMIDS%id_ASOA1)/ &
                (VOL*State_Chm%SpcData(APMIDS%id_ASOA1)%Info%emMW_g* &
                1.d-3)*6.02d17
             CSOA(9) =Spc(I,J,L,APMIDS%id_ASOA2)/ &
                (VOL*State_Chm%SpcData(APMIDS%id_ASOA2)%Info%emMW_g* &
                1.d-3)*6.02d17
             CSOA(10)=Spc(I,J,L,APMIDS%id_ASOA3)/ &
                (VOL*State_Chm%SpcData(APMIDS%id_ASOA3)%Info%emMW_g* &
                1.d-3)*6.02d17

             CSOA(11)=Spc(I,J,L,APMIDS%id_POA1)/ &
                (VOL*State_Chm%SpcData(APMIDS%id_POA1)%Info%emMW_g* &
                1.d-3)*6.02d17
             CSOA(12)=Spc(I,J,L,APMIDS%id_POA2)/ &
                (VOL*State_Chm%SpcData(APMIDS%id_POA2)%Info%emMW_g* &
                1.d-3)*6.02d17

             CSOA(13)=Spc(I,J,L,APMIDS%id_OPOA1)/ &
                (VOL*State_Chm%SpcData(APMIDS%id_OPOA1)%Info%emMW_g* &
                1.d-3)*6.02d17
             CSOA(14)=Spc(I,J,L,APMIDS%id_OPOA2)/ &
                (VOL*State_Chm%SpcData(APMIDS%id_OPOA2)%Info%emMW_g* &
                1.d-3)*6.02d17


             SOAT= (Spc(I,J,L,APMIDS%id_TSOA1)+Spc(I,J,L,APMIDS%id_TSOA2)+  &
                    Spc(I,J,L,APMIDS%id_TSOA3)+Spc(I,J,L,APMIDS%id_TSOA0)+  &
                    Spc(I,J,L,APMIDS%id_ASOA1)+                             &
                    Spc(I,J,L,APMIDS%id_ASOA2)+Spc(I,J,L,APMIDS%id_ASOA3)+  &
                    Spc(I,J,L,APMIDS%id_POA1) +Spc(I,J,L,APMIDS%id_POA2) +  &
                    Spc(I,J,L,APMIDS%id_OPOA1)+Spc(I,J,L,APMIDS%id_OPOA2))/ &
                    VOL ! Total SV-&MV-SOA
          ELSE
             SOAT= 1.d-30
          ENDIF

          MSO4 = 0.d0
          DO N=1,NSO4
             SIZENUM=APMIDS%id_SO4BIN1+N-1
             XM1D(N)=Spc(I,J,L,SIZENUM)/VOL  !kg/m3
             MSO4 = MSO4 + XM1D(N)   ! total bin sulfate mass
             XN1D(N)=XN4D(I,J,L,N)
          ENDDO

          DO N=1,NSEA        ! sea salt
             SIZENUM=APMIDS%id_SEABIN1+N-1
             XM1D(NSO4+N)=Spc(I,J,L,SIZENUM)/VOL  !kg/m3
          ENDDO

          DO N=1,NDSTB      ! dust
             XMDST(N) = Spc(I,J,L,(APMIDS%id_DSTBIN1+N-1))/VOL   !kg/m3
             IF(XMDST(N).LT.1.d-20) THEN   !debug force
                XMDST(N) = 1.d-20
                Spc(I,J,L,(APMIDS%id_DSTBIN1+N-1)) = XMDST(N) * VOL
             ENDIF
             MWSIZE3D(I,J,L,3)=MWSIZE3D(I,J,L,3)+ &
                               Spc(I,J,L,(APMIDS%id_DSTBIN1+N-1))*RDST(N)
          ENDDO

          DO N= 1, NBCOC
             SIZENUM=APMIDS%id_BCBIN1+N-1
             MBC(N)=Spc(I,J,L,SIZENUM)/VOL  !kg/m3
             SIZENUM=APMIDS%id_OCBIN1+N-1
             MOC(N)=Spc(I,J,L,SIZENUM)/VOL  !kg/m3
          ENDDO

          MSO4BULK= Spc(I,J,L,APMIDS%id_SO4)/VOL   !kg/m3
          MNIT    = Spc(I,J,L,APMIDS%id_NIT)/VOL   !kg/m3
          MNH4    = Spc(I,J,L,APMIDS%id_NH4)/VOL   !kg/m3
          MMSA    = Spc(I,J,L,APMIDS%id_MSA)/VOL   !kg/m3

          ! coated sulfate
          MBCS   = Spc(I,J,L,(APMIDS%id_CTBC))/VOL ! SULF on BC kg/m3
          MOCS   = Spc(I,J,L,(APMIDS%id_CTOC))/VOL ! SULF on OC kg/m3
          MDSTS  = Spc(I,J,L,APMIDS%id_CTDST)/VOL ! SULF on DUST kg/m3
          MSALTS = Spc(I,J,L,APMIDS%id_CTSEA)/VOL ! SULF on SEA-SALT kg/m3

          IF(IFCOATBC.EQ.0)  MBCS = 1.d-30

          NCOAG1 = NCOAG3D(I,J,L,1)
          NCOAG2 = NCOAG3D(I,J,L,2)

          YCLDLIQ=  State_Met%QL(I,J,L)*State_Met%AIRDEN(I,J,L)*1E3  ! Cloud liquid water content [g/m3]

          VZ = -(State_Met%OMEGA(I,J,L))/(State_Met%AIRDEN(I,J,L)*9.8)  ! vertical velocity in m/s
          !VZ = 3.d0 ! vertical velocity in m/s

          CALL APM_PHYS(I,J,L, &
                 NCOAG1,NCOAG2,IACT10,IACT20,IACT30, &
                 RACT1,RACT2,RACT3,NTEMPOUT1, &
                 State_Met%AIRDEN(I,J,L), &
                 PRESS,YSIGMA,TK,RHIN,XQ,PLVSOG01,PLVSOG1,CACID,PACID,CNH3, &
                 DTAPM,MMSA,MNIT,MNH4,MBCS,MOCS,MDSTS,MSALTS,MBC,MOC, &
                 SOAT,CLVSOG,MSULFLV,MBCLV,MOCLV,MDSTLV,MSALTLV, &
                 GFTOT1,GFTOT2,DENWET1,DENWET2, &
                 YSPGF, XBCLIFE,XOCLIFE, &
!OPT+         &       XM1D,XN1D,TEMPOUT1,XMDST,FCLOUD1)
                 VZ,YCLDLIQ, XCDN,XCDNSP, &
                 XM1D,XN1D,TEMPOUT1,XMDST,FCLOUD1, &
                 ZBEXT,ZW,ZG,ZBABS,XBEXT1k, &       !LW absorption for dust only
                 YBEXT,YW,YG,ACS,XAMINE,AERAREA1,AERDRYR1, &
                 GAMMAPM1,AEROCOMOUT1D,IFOUTIJ(I,J), &
                 ATOM4N)

          TEMPOUT1(56)=TEMPOUT1(56)+ &
                       (Spc(I,J,L,APMIDS%id_TSOA1)+Spc(I,J,L,APMIDS%id_TSOA2)+ &
                        Spc(I,J,L,APMIDS%id_TSOA3)+Spc(I,J,L,APMIDS%id_TSOA0)+ &
                        +Spc(I,J,L,APMIDS%id_ASOA1)+ &
                        Spc(I,J,L,APMIDS%id_ASOA2)+Spc(I,J,L,APMIDS%id_ASOA3)) &
                       *60.d9/(92.d0*VOL)
          IF(IFAG.EQ.1) THEN
             TEMPOUT1(56)=TEMPOUT1(56)+ &
               (Spc(I,J,L,APMIDS%id_POA1)+Spc(I,J,L,APMIDS%id_POA2)+ &
                Spc(I,J,L,APMIDS%id_OPOA1)+Spc(I,J,L,APMIDS%id_OPOA2))*1.d9/VOL
          ENDIF

          TEMPOUT1(59)=State_Met%PHIS(I,J)+SUM(State_Met%BXHEIGHT(I,J,1:L))

          IF(IFOUTIJ(I,J).EQ.1)THEN
             AEROCOMOUT1D(23:27)=TEMPOUT1(54:58)
             AEROCOMOUT1D(3)=TEMPOUT1(59)
          ENDIF

          XBEXT1k3D(I,J,L,1)=SUM(XBEXT1k(1:2,3))*State_Met%BXHEIGHT(I,J,L)*100.
          XBEXT1k3D(I,J,L,2)=(XBEXT1k(3,3))*State_Met%BXHEIGHT(I,J,L)*100.
          XBEXT1k3D(I,J,L,3)=(XBEXT1k(4,3))*State_Met%BXHEIGHT(I,J,L)*100.
          XBEXT1k3D(I,J,L,4)=SUM(XBEXT1k(5:6,3))*State_Met%BXHEIGHT(I,J,L)*100.
          XBEXT1k3D(I,J,L,5)=(XBEXT1k(7,3))*State_Met%BXHEIGHT(I,J,L)*100.
          XBEXT1k3D(I,J,L,6)=SUM(XBEXT1k(8:9,3))*State_Met%BXHEIGHT(I,J,L)*100.
          XBEXT1k3D(I,J,L,7)=SUM(XBEXT1k(10:13,3))*State_Met%BXHEIGHT(I,J,L)*100.
          XBEXT1k3D(I,J,L,8)=SUM(XBEXT1k(:,1))*State_Met%BXHEIGHT(I,J,L)*100.
          XBEXT1k3D(I,J,L,9)=SUM(XBEXT1k(:,4))*State_Met%BXHEIGHT(I,J,L)*100.
          XBEXT1k3D(I,J,L,10)=SUM(XBEXT1k(:,5))*State_Met%BXHEIGHT(I,J,L)*100.

          XBEXT1k3D(I,J,L,11)=SUM(XBEXT1k(1:13,2))*State_Met%BXHEIGHT(I,J,L)*100.
          XBEXT1k3D(I,J,L,12)=SUM(XBEXT1k(14:18,2))*State_Met%BXHEIGHT(I,J,L)*100.

          !  CS is effective uptake CS from apm_phys,  oxidation rate from Carl and Crowley, 1998.
          !XSINK =  ACS + XOH*6.49d-11     ! oxid rate for DMA

          OXRATE(1) = 1.73d-11
          OXRATE(2) = 6.49d-11
          OXRATE(3) = 3.58d-11

          IF(IFATHN>0)THEN
             DO ITYP = 1, 3
                XSINK =  ACS + XOH*OXRATE(ITYP)     ! oxid rate for DMA

                IF(XSINK==0.D0)THEN
                   Spc(I,J,L,APMIDS%id_AMINE+ITYP-1)= &
                      Spc(I,J,L,APMIDS%id_AMINE+ITYP-1) &
                     +YAMINEEMIT(ITYP)*DTAPM
                   UPTAKE(I,J,L,1+ITYP) = 0.
                   OXIDAT(I,J,L,1+ITYP) = 0.
                ELSE
                   XX0 = Spc(I,J,L,APMIDS%id_AMINE+ITYP-1)+YAMINEEMIT(ITYP)*DTAPM
                   Spc(I,J,L,APMIDS%id_AMINE+ITYP-1)=YAMINEEMIT(ITYP)/XSINK- &
                    (YAMINEEMIT(ITYP)-XSINK*Spc(I,J,L,APMIDS%id_AMINE+ITYP-1))*&
                    EXP(-XSINK*DTAPM)/XSINK

                   DXX = (Spc(I,J,L,APMIDS%id_AMINE+ITYP-1)-XX0)*14./ &
                         MWAMINE(ITYP)   !kgN/box-step
                   DXX = DXX*8.64d4/DTAPM   !kgN/box-day
                   UPTAKE(I,J,L,1+ITYP) = DXX * ACS/XSINK
                   OXIDAT(I,J,L,1+ITYP) = DXX * (1.-ACS/XSINK)
                ENDIF

                ! DRYDEP from sulfate_mod.f in v/v-step  --> kgN/box-day
                DRYDEP(I,J,L,1+ITYP)=DRYDEP(I,J,L,1+ITYP)*14.d-3 &
                         *State_Met%AD(I,J,L) &   !AD -- kg dry air/box
                         *8.64d4/(MAIR*DT)

                ! CONDEP from sulfate_mod.f in v/v-step  --> kgN/box-day
                CONDEP(I,J,L,1+ITYP)=CONDEP(I,J,L,1+ITYP)*14.d-3 &
                         *State_Met%AD(I,J,L) &   !AD -- kg dry air/box
                         *8.64d4/(MAIR*DT)

                ! WETDEP from main.f in kg/box-step  --> kgN/box-day
                WETDEP(I,J,L,1+ITYP)=WETDEP(I,J,L,1+ITYP)*14./MWAMINE(ITYP) &
                         *8.64d4/DT

             ENDDO
          ENDIF

          ! DRYDEP and WETDEP for NH3
          DRYDEP(I,J,L,1)=DRYDEP(I,J,L,1)*14.d-3*State_Met%AD(I,J,L) &  !AD -- kg dry air/box
               *8.64d4/(MAIR*DT)
          CONDEP(I,J,L,1)=CONDEP(I,J,L,1)*14.d-3*State_Met%AD(I,J,L) &  !AD -- kg dry air/box
               *8.64d4/(MAIR*DT)
          WETDEP(I,J,L,1)=WETDEP(I,J,L,1)*14./17*8.64d4/DT

          ! NH3 oxidation
          XSINK =   XOH*1.6d-13     ! oxid rate for NH3 1.6d-13 (Atkinson, 1997)
          XX0 = Spc(I,J,L,APMIDS%id_NH3)
          Spc(I,J,L,APMIDS%id_NH3)=Spc(I,J,L,APMIDS%id_NH3)*EXP(-XSINK*DTAPM)
          DXX = (Spc(I,J,L,APMIDS%id_NH3)-XX0)*14./17.   !kgN/box-step
          DXX = DXX*8.64d4/DTAPM   !kgN/box-day
          OXIDAT(I,J,L,1) = DXX

          ! NH3 uptake lifetime
          XX1 = exp((9.02 - TK*0.0359)/(8.3144*1.d-3*TK/4.186))
          XCSNH3=XX1/(1.+XX1)*ACS/AGAMA*SQRT(98./17.)
          TAONH3(I,J,L)=1./XCSNH3    ! in s

          ! NH3 uptake from ISOROPIA
          UPTAKE(I,J,L,1) = UPTAKE(I,J,L,1)*14./17.*8.64d4/DT   !from kg NH3/box-step to kgN/box-day

          IF(IFOPT.EQ.1) THEN
             DO IWL = 1, MWL
                ZBEXT3D(I,J,L,IWL) = ZBEXT(IWL)
                ZW3D(I,J,L,IWL) = ZW(IWL)
                ZG3D(I,J,L,IWL) = ZG(IWL)
                DO ITYP = 1,NTYP
                   YBEXT3D(I,J,L,ITYP,IWL) = YBEXT(ITYP,IWL)
                   YW3D(I,J,L,ITYP,IWL) = YW(ITYP,IWL)
                   YG3D(I,J,L,ITYP,IWL) = YG(ITYP,IWL)
                ENDDO
             ENDDO

             DO IWL = 1, NWL
                ZBABS3D(I,J,L,IWL) = ZBABS(IWL)
             ENDDO
          ENDIF

          GFTOT3D(I,J,L,1) = GFTOT1
          GFTOT3D(I,J,L,2) = GFTOT2

          XNOLD=XN1D(1)
          DO N=2,NSO4
             IF(XN1D(N)>XNOLD )THEN
                GFTOT3D(I,J,L,3) = N*1.D0
                XNOLD=XN1D(N)
             ENDIF
          ENDDO

          DENWET3D(I,J,L,1) = DENWET1
          DENWET3D(I,J,L,2) = DENWET2
          IACT1(I,J,L) = RACT1-1
          IACT2(I,J,L) = RACT2-1
          IACT3(I,J,L) = RACT3-1
          DO N=1,NSO4+4
             FCLOUD(I,J,L,N) = FCLOUD1(N)
          ENDDO
          AERAREA(I,J,L,:) = AERAREA1(:)
          AERDRYR(I,J,L,:) = AERDRYR1(:)
          GAMMAPM(I,J,L,:) = GAMMAPM1(:)
          NCOAG3D(I,J,L,1)=NCOAG1 ! count step that coag is not call in the grid
          NCOAG3D(I,J,L,2)=NCOAG2 ! count step that coag is not call in the grid

          Spc(I,J,L,APMIDS%id_SO4G)=  CACID*VOL*M1ACID*1.d6

          DO N=1,NSO4
             SIZENUM=APMIDS%id_SO4BIN1+N-1
             Spc(I,J,L,SIZENUM)=XM1D(N)*VOL
             MWSIZE3D(I,J,L,1)=MWSIZE3D(I,J,L,1)+ &
                               Spc(I,J,L,SIZENUM)*RDRY(N)*GFTOT1
          ENDDO
          Spc(I,J,L,APMIDS%id_CTBC) = MBCS * VOL
          Spc(I,J,L,(APMIDS%id_CTOC)) = MOCS * VOL
          Spc(I,J,L,APMIDS%id_CTDST) = MDSTS * VOL
          Spc(I,J,L,APMIDS%id_CTSEA) = MSALTS * VOL

          IF(IFAG.EQ.1) THEN
             Spc(I,J,L,APMIDS%id_SO4G+1)=  CLVSOG*VOL*M1LVSOG*1.d6
             Spc(I,J,L,APMIDS%id_CTSO4)= MSULFLV * VOL   !kg
             Spc(I,J,L,APMIDS%id_CTBC+1) = MBCLV * VOL   !kg
             Spc(I,J,L,APMIDS%id_CTOC+1) = MOCLV * VOL   !kg
             Spc(I,J,L,APMIDS%id_CTDST+1) = MDSTLV * VOL   !kg
             Spc(I,J,L,APMIDS%id_CTSEA+1) = MSALTLV * VOL   !kg
          ENDIF

          DO N=1,NSEA        ! sea salt
             SIZENUM=APMIDS%id_SEABIN1+N-1
             Spc(I,J,L,SIZENUM)=XM1D(NSO4+N)*VOL  !kg
             MWSIZE3D(I,J,L,2)=MWSIZE3D(I,J,L,2)+ &
                               Spc(I,J,L,SIZENUM)*RSALT(N)*GFTOT2
          ENDDO

          ! no update for dust needed for now because of no coag yet (need update if coag)

          !BCLIFE(I,J,L) = TEMPOUT1(47)
          !OCLIFE(I,J,L) = TEMPOUT1(48)
          BCLIFE(I,J,L) = XBCLIFE
          OCLIFE(I,J,L) = XOCLIFE

          DO N= 1, NBCOC
             SIZENUM=APMIDS%id_BCBIN1+N-1
             Spc(I,J,L,SIZENUM)=MBC(N)*VOL  !kg/m3
             SIZENUM=APMIDS%id_OCBIN1+N-1
             Spc(I,J,L,SIZENUM)=MOC(N)*VOL  !kg/m3
          ENDDO

          ! Yu++ aerosol-cloud interaction
          ! CCN0.4

          !YCCN = TEMPOUT1(11)+TEMPOUT1(12) &
          ! +TEMPOUT1(13)+TEMPOUT1(14)+TEMPOUT1(15)
          !YCCN = TEMPOUT1(6)                  ! mxy+ YCCN here for CN3
          !YCDN=3.75D8 * (1.-exp(-2.5D-9*YCCN*1.d6))  ! YCCN in #/cm3, YCDN in #/m3

          YCDN = XCDN*1.D6     ! XCDN calculated based on S at a given updraft v
          YCDNSP = XCDNSP*1.D6 ! XCDN calculated based on S at a given updraft v
          ! test
          !YCDN = 2.*XCDN*1.D6  ! XCDN calculated based on S at a given updraft v

          !YCDN= MAX(YCDN, 1.d0)
          !YCDN= MAX(YCDN, 10.d0)
          !Luo YCDN= MAX(YCDN, 10.d6)
          CDN3D(I,J,L) = YCDN

          !Luo IF(L.GE.2) THEN
          !       YF = MIN(1.d0,State_Met%CMFMC(I,J,L)/1.d-2)
          !Luo    YF = MIN(1.d0,State_Met%CMFMC(I,J,L)/5.d-3)
          !Luo    YC = State_Met%AIRDEN(I,J,L)/State_Met%AIRDEN(I,J,L-1)
          !Luo    YCDN = CDN3D(I,J,L)*(1.-YF)+CDN3D(I,J,L-1)*YF*YC
          !Luo    CDN3D(I,J,L) = YCDN
          !Luo ENDIF

          YCLDF = State_Met%CLDF(I,J,L)

          IF(YCLDF.GT.1.D-3)THEN
             !! convert water content from kg/kg to g/m3
             YCLDLIQ= State_Met%QL(I,J,L)* &
                State_Met%AIRDEN(I,J,L)*1D3/YCLDF  ! Cloud liquid water content [g/m3]
             YCLDICE= State_Met%QI(I,J,L)* &
                State_Met%AIRDEN(I,J,L)*1D3/YCLDF  ! Cloud ice water content [g/m3]

             !IF(YCLDF.GT.1.d-5) THEN   ! RADF code use CUT=0.001
             ! IF(PRESS.GT.7.d4) THEN    !effect limit to cloud below 700 mb
             !Luo YRCLDL=(3.*YCLDLIQ/(4.*3.14*1.d6* YCDN))**(1./3.)*1.d6  ! in um
             !    YRCLDL = 1.12*YRCLDL   ! 1.12 to take into account relation between rvol and reff (Chen and Penner, ACP, 2005)
             !    YRCLDL = 1.2*YRCLDL   !test

             dumc = State_Met%QL(I,J,L)/YCLDF

             if(dumc>1.d-15)then
                dumnc = YCDN/State_Met%AIRDEN(I,J,L)

                call size_dist_param_liq(dumc, dumnc, &
                   0.d0, State_Met%AIRDEN(I,J,L), .true., &
                   pgam, lamc)

                YRCLDL = 0.7*gamma(pgam+4.)/gamma(pgam+3.)/ &
                   lamc*1.e6

             else
                !YRCLDL=(3.*YCLDLIQ/(4.*3.14*1.d6* YCDN))**(1./3.)*1.d6
                YRCLDL=55.d0
             endif

             !ELSE
             ! YRCLDL = 5.0
             !ENDIF

          ELSE
             YCLDLIQ=0.d0
             YCLDICE=0.d0
             YCLDF=1.d-3
             YRCLDL= 55.0
          ENDIF

          RCLDL3D(I,J,L)= MAX(2.5d0,YRCLDL)   !  ice particle radius
          RCLDL3D(I,J,L)= MIN(55.d0,YRCLDL)   !  ice particle radius

          YRCLDL = RCLDL3D(I,J,L)

          ! Ice Nuclei

          !YNA(1) = TEMPOUT1(16)  !SP in #/cm3
          !YNA(1) = TEMPOUT1(16)/3.  !SP in #/cm3
          YNA(1) = TEMPOUT1(11)  !SP in #/cm3
          YNA(2) = TEMPOUT1(13)  !dust
          !YNA(3) = 0.0  !soot  -- no soot IN for now
          !YNA(3) = TEMPOUT1(14)/100.  ! assume soot a factor 100 less effective compared to dust
          YNA(3) = TEMPOUT1(14)

          !IF(YCLDICE.GT.1.d-10.and.YCLDF.GT.1.d-5) THEN
          !
          !   wbar = 0.25
          !   relhum = RHIN/100.  ! (RH: 0-1)
          !   yqc   = State_Met%QL(I,J,L)   !kg/kg
          !
          !   call nucleati(wbar, TK, relhum, YCLDF, yqc, &
          !        YNA, nuci, onihf, oniimm, onidep, onimey)
          !
          !   IN3D(I,J,L) = nuci   !#/cm3
          !   IF(nuci.GT.1.d-4) THEN
          !    YREI=(3.*YCLDICE/(YCLDF*4.*3.14*0.92d6*nuci*1.d6)) &
          !                                   **(1./3.)*1.d6  ! in um
          !   ELSE
          !    YREI= 30.0
          !   ENDIF
          !
          !   IF(TK.LE.240.) THEN
          !     YB = -0.03387
          !   ELSE
          !     YB = 0.
          !   ENDIF
          !   IF(YCLDICE.LE.1.) THEN
          !     YC = -0.014738
          !   ELSE
          !     YC = 0.
          !   ENDIF
          !   YK = exp(-3.15393+YB*(TK-240.)+YC*log(YCLDICE))
          !
          !   YREI = YREI * (1./YK)**(1./3.)  ! relation between rvol and reff -- Liu et al., 2007
          !ELSE
          !  onihf = 0.
          !  oniimm = 0.
          !  onidep = 0.
          !  onimey = 0.
          !  IN3D(I,J,L) = 0.0
          !  YREI= 0.0
          !  YK = 1.
          !ENDIF
          !! RCLDI3D(I,J,L) =  YREI
          !! RCLDI3D(I,J,L) =  YREI/2.0            ! mxy+ increase cloud rf and palb
          !RCLDI3D(I,J,L) =  50.            ! test with fixed REI (in um)
          !RCLDI3D(I,J,L) =  100.*SQRT(YCLDICE/1.d-2)             ! test with fixed REI (in um), YCLDICE in g/m3
          !RCLDI3D(I,J,L) =  50.*YCLDICE/1.d-2             ! test with fixed REI (in um), YCLDICE in g/m3
          !RCLDI3D(I,J,L) =  30.*YCLDICE/1.d-2             ! test with fixed REI (in um), YCLDICE in g/m3
          !RCLDI3D(I,J,L) = 0.5*53.005*YCLDICE**0.06*EXP((TK-273.15)*0.013)! para based on Boudala et al., 2002
          !RCLDI3D(I,J,L) = 1.12*RCLDI3D(I,J,L)   ! 1.12 to take into account relation between rvol and reff (Chen and Penner, ACP, 2005)
          !RCLDI3D(I,J,L) = 1.2*RCLDI3D(I,J,L)   ! use 1.2 for now
          !RCLDI3D(I,J,L) = 1.25*RCLDI3D(I,J,L)   ! use 1.2 for now
          !Luo RCLDI3D(I,J,L) =  0.5*53.005*YCLDICE**0.06*EXP((TK-273.15)*0.013)   ! para based on Boudala et al., 2002
          !Luo RCLDI3D(I,J,L) = 2.0*RCLDI3D(I,J,L)   ! use 1.2 for now

          !Readin MODIS RCLDI
          IF(YCLDICE>0.D0)THEN
             !Luo RCLDI3D(I,J,L) = RCLDI2D(I,J)*1.d0
             YC = 0.124*exp(min(0.,max(-35.,(TK-273.15)))*0.038)
             YB = 0.114*exp(min(0.,max(-30.,(TK-273.15)))*0.054)

             YK = YC/(YB+YC)

             RCLDI3D(I,J,L) = 0.7* &
                (   YK *53.005*(max(0.01,(YCLDICE*YK))**(0.06))* &
                       EXP(min(0.,max(-35.,(TK-273.15)))*0.013) + &
                (1.-YK)*57.133*(max(0.01,(YCLDICE*(1.-YK)))**(-0.0313))* &
                       EXP(min(0.,max(-30.,(TK-273.15)))*0.011) )
          ELSE
             RCLDI3D(I,J,L) = 120.d0
          ENDIF

          ! set to RRTMG limit
          RCLDI3D(I,J,L)= MAX(5.0d0,RCLDI3D(I,J,L))   !  ice particle radius
          RCLDI3D(I,J,L)= MIN(140.d0,RCLDI3D(I,J,L))   !  ice particle radius

          YREI = RCLDI3D(I,J,L)

          !IF((MOD(I,10).EQ.1).and.(MOD(J,10).EQ.1).and.(MOD(L,5).EQ.1) &
          !   .and.L.LT.28)    THEN
          ! WRITE(111,108)I,J,L, PRESS, TK, &
          ! YCLDF,YCLDLIQ,YCCN,YCDN,YRCLDL, &
          ! (1./YK)**(1./3.), &
          ! YCLDICE,YNA(1),YNA(2),YNA(3),IN3D(I,J,L), &
          ! onihf,oniimm,onidep,onimey,YREI
          !ENDIF
          ! 108   FORMAT(I3,I3,I3,50(1PE9.2))
          !
          ! Temporally put OH in TEMPOUT1(4)
          !TEMPOUT1(4)=XOH
          !TEMPOUT1(4)=XSINK
          ! Temporally used the last four tempout1 to save CLDF,CLDLIQ,CDN,and RCLDL for output
          !TEMPOUT1(NTEMPOUT1-7)=YCLDLIQ
          !TEMPOUT1(NTEMPOUT1-6)=YCLDICE
          !TEMPOUT1(NTEMPOUT1-4)=YCLDF
          !TEMPOUT1(NTEMPOUT1-3)=YCDN*YCLDLIQ  !cloud liq weighted CDN
          !TEMPOUT1(NTEMPOUT1-2)=YRCLDL*YCLDLIQ   !cloud liq weighted RCLDL
          !TEMPOUT1(NTEMPOUT1-1)=IN3D(I,J,L)*YCLDICE  !cloud ice weighted CDN
          !TEMPOUT1(NTEMPOUT1)=YREI*YCLDICE   !cloud ice weighted RCLDL
          !
          !TEMPOUT1(19) = YUV(I,J,L,7)  !UV flux (289-320 nm)
          !TEMPOUT1(20) = YUV(I,J,L,8)  !UV flux (320-412 nm)
          !TEMPOUT1(21) = XO3(I,J,L)*1.d9  ! ppb
          !
          !TEMPOUT1(42) = YUV(I,J,L,1)  !NO2 phtolysis rate
          !TEMPOUT1(43) = YUV(I,J,L,2)  !HNO3 phtolysis rate
          !TEMPOUT1(44) = YUV(I,J,L,3)  !H2O2 phtolysis rate
          !TEMPOUT1(45) = YUV(I,J,L,4)  !CH2O phtolysis rate
          !TEMPOUT1(46) = YUV(I,J,L,5)  !O3 phtolysis rate
          !TEMPOUT1(47) = YUV(I,J,L,6)  !POH
          !! Heating profiles
          !!TEMPOUT1(57) = YHTRC03D(I,J,L)
          !!TEMPOUT1(58) = YHTRC3D(I,J,L)
          !!!TEMPOUT1(59) = YHTR03D(I,J,L)
          !!TEMPOUT1(60) = YHTR3D(I,J,L)
          !
          !TEMPOUT1(NTEMPOUT1-21)=YCDN*YCLDLIQ*YCLDF*YCLDF
          !TEMPOUT1(NTEMPOUT1-20)=YRCLDL*YCLDLIQ*YCLDF*YCLDF
          !IF(TK>233.D0)THEN
          !  TEMPOUT1(NTEMPOUT1-19)=State_Met%OPTDEP(I,J,L)
          !  TEMPOUT1(NTEMPOUT1-18)=0.D0
          !ELSE
          !  TEMPOUT1(NTEMPOUT1-19)=0.D0
          !  TEMPOUT1(NTEMPOUT1-18)=State_Met%OPTDEP(I,J,L)
          !ENDIF
          !
          !IF(TCOD3D(I,J,L,9)>1.D-3)THEN
          !  TEMPOUT1(NTEMPOUT1-17)=TCOD3D(I,J,L,4)/(TCOD3D(I,J,L,9)**1.5)
          !  TEMPOUT1(NTEMPOUT1-16)=TCOD3D(I,J,L,8)/(TCOD3D(I,J,L,9)**1.5)
          !ELSE
          !  TEMPOUT1(NTEMPOUT1-17)=0.D0
          !  TEMPOUT1(NTEMPOUT1-16)=0.D0
          !ENDIF
          !
          !IF(TK>233.D0)THEN
          !  TEMPOUT1(NTEMPOUT1-15)=State_Met%OPTDEP(I,J,L)*(YCLDF**1.5)
          !  TEMPOUT1(NTEMPOUT1-14)=0.D0
          !ELSE
          !  TEMPOUT1(NTEMPOUT1-15)=0.D0
          !  TEMPOUT1(NTEMPOUT1-14)=State_Met%OPTDEP(I,J,L)*(YCLDF**1.5)
          !ENDIF
          !TEMPOUT1(NTEMPOUT1-13)=TCOD3D(I,J,L,4)
          !TEMPOUT1(NTEMPOUT1-12)=TCOD3D(I,J,L,8)

          TEMPOUT1(NTEMPOUT1-21)=YCLDLIQ*YCLDF
          TEMPOUT1(NTEMPOUT1-20)=YCLDICE*YCLDF
          TEMPOUT1(NTEMPOUT1-19)=YCLDF
          TEMPOUT1(NTEMPOUT1-18)=YCDN*YCLDLIQ*YCLDF    !cloud liq weighted CDN
          TEMPOUT1(NTEMPOUT1-13)=YCDNSP*YCLDLIQ*YCLDF  !cloud liq weighted CDN
          TEMPOUT1(NTEMPOUT1-17)=YRCLDL*YCLDLIQ*YCLDF  !cloud liq weighted RCLDL
          TEMPOUT1(NTEMPOUT1-16)=IN3D(I,J,L)*YCLDICE*YCLDF !cloud ice weighted CDN
          TEMPOUT1(NTEMPOUT1-15)=YREI*YCLDICE*YCLDF    !cloud ice weighted RCLDL

          IF(YCLDF>0.8)THEN
             IF(State_Met%T(I,J,L)>268.)THEN
                IF(YCLDLIQ>1.D-20)THEN
                   MODISOUT(I,J,1)=MODISOUT(I,J,1)+YCDN
                   MODISOUT(I,J,2)=MODISOUT(I,J,2)+YRCLDL

                   MODISOUT(I,J,3)=MODISOUT(I,J,3)+YCDN*YCLDLIQ*YCLDF &
                                   *State_Met%BXHEIGHT(I,J,L)
                   MODISOUT(I,J,4)=MODISOUT(I,J,4)+YRCLDL*YCLDLIQ*YCLDF &
                                   *State_Met%BXHEIGHT(I,J,L)

                   MODISOUT(I,J,5)=MODISOUT(I,J,5)+1.
                   MODISOUT(I,J,6)=MODISOUT(I,J,6)+YCLDLIQ*YCLDF &
                        *State_Met%BXHEIGHT(I,J,L)
                ENDIF
             ENDIF
          ENDIF

          TEMPOUT1(NTEMPOUT1-14)=YUV(I,J,L,7)
          !TEMPOUT1(NTEMPOUT1-13)=YUV(I,J,L,8)
          !TEMPOUT1(NTEMPOUT1-12)=YUV(I,J,L,9)
          !TEMPOUT1(NTEMPOUT1-13)=XQ
          TEMPOUT1(NTEMPOUT1-12)=RHIN

          !TEMPOUT1(NTEMPOUT1-11)=YCLDLIQ*YCLDF*State_Met%BXHEIGHT(I,J,L)
          !TEMPOUT1(NTEMPOUT1-10)=YCLDICE*YCLDF*State_Met%BXHEIGHT(I,J,L)
          TEMPOUT1(NTEMPOUT1-11)=YSPGF
          TEMPOUT1(NTEMPOUT1-10)=0.D0   !reserved for later use?
          !TEMPOUT1(NTEMPOUT1-9)=State_Met%CLDF(I,J,L)
          !TEMPOUT1(NTEMPOUT1-9)=YSPGF
          TEMPOUT1(NTEMPOUT1-9)=EMITSO2(I,J,L)*32./64. &
              *8.64d4/XAREA ! SO2 emission -- kg S/km2-day

          !TEMPOUT1(NTEMPOUT1-8)=CACID
          !TEMPOUT1(NTEMPOUT1-7)=CLVSOG
          TEMPOUT1(NTEMPOUT1-8)=PACID
          TEMPOUT1(NTEMPOUT1-7)=PLVSOG1
          TEMPOUT1(NTEMPOUT1-6)=PLVSOG01
          TEMPOUT1(NTEMPOUT1-5)=XOH

          !TEMPOUT1(NTEMPOUT1-5)=YAMINEEMIT(2)*14./MWAMINE(2)!kg-N/box-sec (DMA)
          !TEMPOUT1(NTEMPOUT1-4)=CAMINEEMIT(2)   ! # cm-3s-1
          !TEMPOUT1(NTEMPOUT1-4)=EMITNH3(I,J,L)*14./17.
          TEMPOUT1(NTEMPOUT1-4)=NH3EMIT &
              *8.64d4/XAREA ! NH3 emission -- kg N/km2-day

!-------------------------------------------------------------------------------
! Prior to 6/28/19:
! Comment out problematic code.  Gan Luo says that we can set all of these
! to zero, because we do not simulate amines at this point.
! (bmy, 6/28/19)
!          TEMPOUT1(NTEMPOUT1-3)=YAMINEEMIT(2)*14./MWAMINE(2)   !kg-N/box-sec (DMA)
!          !TEMPOUT1(NTEMPOUT1-3)=CNH3
!          TEMPOUT1(NTEMPOUT1-2)=CAMINE(1)
!          TEMPOUT1(NTEMPOUT1-1)=CAMINE(2)
!          TEMPOUT1(NTEMPOUT1)  =CAMINE(3)
!------------------------------------------------------------------------------
          TEMPOUT1(NTEMPOUT1-3)=0.0d0
          !TEMPOUT1(NTEMPOUT1-3)=CNH3
          TEMPOUT1(NTEMPOUT1-2)=0.0d0
          TEMPOUT1(NTEMPOUT1-1)=0.0d0
          TEMPOUT1(NTEMPOUT1)  =0.0d0

          !TEMPOUT1(NTEMPOUT1-5) = YCDN*YCLDLIQ*YCLDF !cloud liq weighted CDN
          !TEMPOUT1(NTEMPOUT1-4) = YRCLDL*YCLDLIQ*YCLDF !cloud liq weighted RCLDL

          !TEMPOUT1(NTEMPOUT1-5) = State_Met%OPTDEP(I,J,L)
          !TEMPOUT1(NTEMPOUT1-4) = YCOD3D(I,J,L)/(TCOD3D(I,J,L,9)**1.5)
          !
          !TEMPOUT1(NTEMPOUT1-3) = State_Met%OPTDEP(I,J,L) * YCLDF
          !TEMPOUT1(NTEMPOUT1-2) = YCOD3D(I,J,L) &
          !                      * TCOD3D(I,J,L,9)/(TCOD3D(I,J,L,9)**1.5)
          !TEMPOUT1(NTEMPOUT1-1) = State_Met%OPTDEP(I,J,L) &
          !                      *(YCLDF**1.5)
          !TEMPOUT1(NTEMPOUT1)   = YCOD3D(I,J,L)

          !TEMPOUT1(NTEMPOUT1-24)=EMITNH3(I,J,L) ! kg NH3/box-sec
          !TEMPOUT1(NTEMPOUT1-24)=EMITNH3(I,J,L)*14./17.
          !TEMPOUT1(48)=EMITNH3(I,J,L)*14./17.
          !    *8.64d4/XAREA ! NH3 emission -- kg N/km2-day

          !Gan Luo+ Satellite output
          !CODOUT(I,J,L,1)=CODOUT(I,J,L,1)+GOOD(I)*YCLDF
          !
          !CODOUT(I,J,L,2)=CODOUT(I,J,L,2)+GOOD(I)*State_Met%OPTDEP(I,J,L)
          !
          !CODOUT(I,J,L,3)=CODOUT(I,J,L,3)+GOOD(I)* &
          !  YCOD3D(I,J,L)/(TCOD3D(I,J,L,9)**1.5)
          !CODOUT(I,J,L,4)=CODOUT(I,J,L,4)+GOOD(I)* &
          !  TCOD3D(I,J,L,4)/(TCOD3D(I,J,L,9)**1.5)
          !CODOUT(I,J,L,5)=CODOUT(I,J,L,5)+GOOD(I)* &
          !  TCOD3D(I,J,L,8)/(TCOD3D(I,J,L,9)**1.5)
          !
          !CODOUT(I,J,L,6)=CODOUT(I,J,L,6)+GOOD(I)*YCLDLIQ
          !CODOUT(I,J,L,7)=CODOUT(I,J,L,7)+GOOD(I)*YCLDICE
          !
          !CODOUT(I,J,L,8)=CODOUT(I,J,L,8)+GOOD(I)*YRCLDL*YCLDLIQ
          !CODOUT(I,J,L,9)=CODOUT(I,J,L,9)+GOOD(I)*YREI*YCLDICE
          !
          !CODOUT(I,J,L,10)=CODOUT(I,J,L,10)+GOOD(I)*YCDN*YCLDLIQ
          !
          !CODOUTNUM(I,J,L)=CODOUTNUM(I,J,L)+GOOD(I)*1.D0

          !TEMPOUT1(NTEMPOUT1-22)= &
          !  MASSISRP(I,J,L,1)*1.d9/(State_Met%AD(I,J,L)) !ppb
          !TEMPOUT1(NTEMPOUT1-21)= &
          !  MASSISRP(I,J,L,2)*1.d9/VOL !ug/m3
          !TEMPOUT1(NTEMPOUT1-20)= &
          !  MASSISRP(I,J,L,3)*1.d9/(State_Met%AD(I,J,L)) !ppb
          !TEMPOUT1(NTEMPOUT1-19)= &
          !  MASSISRP(I,J,L,4)*1.d9/VOL !ug/m3
          !
          !TEMPOUT1(NTEMPOUT1-18)= &
          !  MASSMESA(I,J,L,1)*1.d9/(State_Met%AD(I,J,L)) !ppb
          !TEMPOUT1(NTEMPOUT1-17)= &
          !  MASSMESA(I,J,L,2)*1.d9/VOL !ug/m3
          !TEMPOUT1(NTEMPOUT1-16)= &
          !  MASSMESA(I,J,L,3)*1.d9/(State_Met%AD(I,J,L)) !ppb
          !TEMPOUT1(NTEMPOUT1-15)= &
          !  MASSMESA(I,J,L,4)*1.d9/VOL !ug/m3

          !luogan temp output
          IF(IFTEMPOUT)THEN
             IF(IFATOM>0.5)THEN
                DO N=1,46
                   TEMPOUT(I,J,L,N)=TEMPOUT1(N)
                ENDDO
                TEMPOUT(I,J,L,47) =YUV(I,J,L,6)    !POH
                TEMPOUT(I,J,L,48) =XO3(I,J,L)*1.d9 ! ppb
                TEMPOUT(I,J,L,49) =CNO    ! ppb
                TEMPOUT(I,J,L,50) =CNO2   ! ppb
                TEMPOUT(I,J,L,51) =CNO3   ! ppb
                TEMPOUT(I,J,L,52) =CHNO3  ! ppb
                TEMPOUT(I,J,L,53) =CISOP  ! ppb
                TEMPOUT(I,J,L,54) =CMTPA  ! ppb
                TEMPOUT(I,J,L,55) =XX1/(1.+XX1)
                TEMPOUT(I,J,L,56) =TAONH3(I,J,L)
                TEMPOUT(I,J,L,57) =CCO    ! ppm
                TEMPOUT(I,J,L,58) =CSO2   ! ppb
                TEMPOUT(I,J,L,59) =XOH
                TEMPOUT(I,J,L,60) =PLVSOG1
                TEMPOUT(I,J,L,61) =CLVSOG
                TEMPOUT(I,J,L,62) =PACID
                TEMPOUT(I,J,L,63) =CACID
                TEMPOUT(I,J,L,64) =MBCS
                TEMPOUT(I,J,L,65) =MOCS
                TEMPOUT(I,J,L,66) =MSALTS
                TEMPOUT(I,J,L,67) =MDSTS
                TEMPOUT(I,J,L,68) =MSO4
                TEMPOUT(I,J,L,69) =MMSA
                TEMPOUT(I,J,L,70) =MNIT
                TEMPOUT(I,J,L,71) =MNH4
                TEMPOUT(I,J,L,72) =MSULFLV
                IF(IFATOM<0.5)THEN
                   TEMPOUT(I,J,L,73) =MBCLV
                   TEMPOUT(I,J,L,74) =MOCLV
                   TEMPOUT(I,J,L,75) =MDSTLV
                   TEMPOUT(I,J,L,76) =MSALTLV
                ELSE
                   TEMPOUT(I,J,L,73) =ATOM4N(1)
                   TEMPOUT(I,J,L,74) =ATOM4N(2)
                   TEMPOUT(I,J,L,75) =ATOM4N(3)
                   TEMPOUT(I,J,L,76) =ATOM4N(4)
                ENDIF
                TEMPOUT(I,J,L,77) =SOAT
             ELSE
                DO N=1,NTEMPOUT1
                   TEMPOUT(I,J,L,N)=TEMPOUT(I,J,L,N)+ TEMPOUT1(N)
                ENDDO
             ENDIF
          ENDIF

          !
          !APM2+IFSITE+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          ! output at selected sites
          IF(IFSITE.GT.0) THEN
             IF(IFOUTIJ(I,J).EQ.1.and.L.LE.LOUT) THEN
                NSITE = SITEID(I,J)

                AEROCOMOUT(NSITE,L,1:27)=AEROCOMOUT1D(1:27)

                DO N=1, NTEMPOUT1
                   SITEOUT2(NSITE,L,N)=TEMPOUT1(N)
                ENDDO
                !SITEOUT2(NSITE,L,NTEMPOUT1+1) = TK
                !SITEOUT2(NSITE,L,NTEMPOUT1+2) = RHIN
                !SITEOUT2(NSITE,L,NTEMPOUT1+3) = CACID
                !SITEOUT2(NSITE,L,NTEMPOUT1+4) = CLVSOG
                !SITEOUT2(NSITE,L,NTEMPOUT1+5) = State_Met%CLDF(I,J,L)
                !SITEOUT2(NSITE,L,NTEMPOUT1+6) = SPGF(I,J,L)

                SITEOUT3(NSITE,L,1)  = YUV(I,J,L,6)    !POH
                SITEOUT3(NSITE,L,2)  = XO3(I,J,L)*1.d9 ! ppb
                SITEOUT3(NSITE,L,3)  = CNO    ! ppb
                SITEOUT3(NSITE,L,4)  = CNO2   ! ppb
                SITEOUT3(NSITE,L,5)  = CNO3   ! ppb
                SITEOUT3(NSITE,L,6)  = CHNO3  ! ppb
                SITEOUT3(NSITE,L,7)  = CISOP  ! ppb
                SITEOUT3(NSITE,L,8)  = CMTPA  ! ppb
                SITEOUT3(NSITE,L,9)  = XX1/(1.+XX1)
                SITEOUT3(NSITE,L,10) = TAONH3(I,J,L)
                SITEOUT3(NSITE,L,11) = CCO    ! ppm

                IF(IFSITEADD.EQ.1) THEN   ! additional output for size dis and others
                   XU = State_Met%U(I,J,L)
                   XV = State_Met%V(I,J,L)

                   TOPP = State_Met%PBL_TOP_M(I,J) ! PBL top pressure (hpa)
                   TOP  = State_Met%PBL_TOP_L(I,J) ! PBL top model levels
                   YTOP(NSITE,1) = TOP
                   YTOP(NSITE,2) = TOPP
                   SITEOUT(NSITE,L,1)=PRESS/100.
                   SITEOUT(NSITE,L,2)=XU
                   SITEOUT(NSITE,L,3)=XV
                   SITEOUT(NSITE,L,4)=TK
                   SITEOUT(NSITE,L,5)=RHIN
                   SITEOUT(NSITE,L,6)=CSO2   ! ppb
                   SITEOUT(NSITE,L,7)=XOH
                   SITEOUT(NSITE,L,8)=PLVSOG1
                   SITEOUT(NSITE,L,9)=CLVSOG
                   SITEOUT(NSITE,L,10)=PACID
                   SITEOUT(NSITE,L,11)=CACID
                   SITEOUT(NSITE,L,12)=MBCS
                   SITEOUT(NSITE,L,13)=MOCS
                   SITEOUT(NSITE,L,14)=MSALTS
                   SITEOUT(NSITE,L,15)=MDSTS
                   SITEOUT(NSITE,L,16)=MSO4
                   SITEOUT(NSITE,L,17)=MMSA
                   SITEOUT(NSITE,L,18)=MNIT
                   SITEOUT(NSITE,L,19)=MNH4
                   SITEOUT(NSITE,L,20)=MSULFLV
                   SITEOUT(NSITE,L,21)=MBCLV
                   SITEOUT(NSITE,L,22)=MOCLV
                   SITEOUT(NSITE,L,23)=MDSTLV
                   SITEOUT(NSITE,L,24)=MSALTLV
                   SITEOUT(NSITE,L,25)=SOAT
                   DO N=1,14
                      SITEOUT(NSITE,L,25+N)=CSOA(N)
                      SITEOUT(NSITE,L,39+N)=CSOG(N)
                   ENDDO
                   SITEOUT(NSITE,L,54)=PLVSOG(I,J,L,1)/(M1LVSOG*1.d+15)  !#/cm3s
                   SITEOUT(NSITE,L,55)=PLVSOG(I,J,L,2)/(M1LVSOG*1.d+15)  !#/cm3s
                   SITEOUT(NSITE,L,56)=PLVSOG(I,J,L,3)/(M1LVSOG*1.d+15)  !#/cm3s
                   SITEOUT(NSITE,L,57)=PLVSOG(I,J,L,4)/(M1LVSOG*1.d+15)  !#/cm3s
                   SITEOUT(NSITE,L,58)=PLVSOG(I,J,L,5)/(M1LVSOG*1.d+15)  !#/cm3s
                   SITEOUT(NSITE,L,59)=MSO4BULK

                   !SITEOUT(NSITE,L,20)=CNH3/1000.   ! ppb
                   !SITEOUT(NSITE,L,35)=CNH3/1000.   ! ppb

                   DO N=1, NSO4+NSEA
                      SITEOUT1(NSITE,L,N)=XM1D(N)
                   ENDDO
                   DO N=1,NDSTB
                      SITEOUT1(NSITE,L,NSO4+NSEA+N)=XMDST(N)
                   ENDDO
                   DO N=1,NBCOC
                      SITEOUT1(NSITE,L,NSO4+NSEA+NDSTB+N)=MBC(N)
                      SITEOUT1(NSITE,L,NSO4+NSEA+NDSTB+NBCOC+N)=MOC(N)
                   ENDDO

                ENDIF

                !SITEOUT2(NSITE,L,1) = State_Met%UWND(I,J,L)
                !SITEOUT2(NSITE,L,1) = State_Met%UWND(I,J,L)
                !SITEOUT2(NSITE,L,2) = State_Met%VWND(I,J,L)
                !SITEOUT2(NSITE,L,1) = TK
                !SITEOUT2(NSITE,L,1) = State_Met%CMFMC(I,J,L+1)
                !SITEOUT2(NSITE,L,2) = VZ
                !SITEOUT2(NSITE,L,3) = TEMPOUT1(25) !SP ext
                !SITEOUT2(NSITE,L,4) = TEMPOUT1(26) !salt ext
                !SITEOUT2(NSITE,L,5) = TEMPOUT1(27) !dst ext
                !SITEOUT2(NSITE,L,6) = TEMPOUT1(29) !BC ext
                !SITEOUT2(NSITE,L,7) = TEMPOUT1(31) !POC ext
                !SITEOUT2(NSITE,L,8) = TEMPOUT1(6)  !CN3
                !SITEOUT2(NSITE,L,9) = TEMPOUT1(7)  !CN10_SP
                !SITEOUT2(NSITE,L,10) = TEMPOUT1(8)  !CN10_PP
                !SITEOUT2(NSITE,L,11) = TEMPOUT1(11)  !CC0.4_SP
                !SITEOUT2(NSITE,L,12) = SUM(TEMPOUT1(12:15))  !CC0.4_PP
                !SITEOUT2(NSITE,L,13) = RHIN
                !!SITEOUT2(NSITE,L,14)= SPGF(I,J,L)
                !SITEOUT2(NSITE,L,14) = State_Met%CLDF(I,J,L)
                !SITEOUT2(NSITE,L,15) = YCLDLIQ
                !SITEOUT2(NSITE,L,16)= YCLDICE
                !IF(State_Met%CLDF(I,J,L).GT.0.) THEN
                !   SITEOUT2(NSITE,L,17)=TEMPOUT1(NTEMPOUT1-1)/State_Met%CLDF(I,J,L)*1.d-6   !CDN in #/cm3
                !   SITEOUT2(NSITE,L,18) = TEMPOUT1(NTEMPOUT1)/State_Met%CLDF(I,J,L)         !RCLDL in um
                !ELSE
                !   SITEOUT2(NSITE,L,17)=0.
                !   SITEOUT2(NSITE,L,18)=0.
                !ENDIF
                !SITEOUT2(NSITE,L,19)= TEMPOUT1(NTEMPOUT1-2)
             ENDIF
          ENDIF
          !APM2+ENDIFSITE++++++++++++++++++++++++++++++++++++++++++++++++++++++

          CCLD(L)= MAX(1.d-3,State_Met%CLDF(I,J,L)) ! Cloud cover [Unitless]

          ! convert water content from kg/kg to g/m3 and then to g/m2
          CLDLIQ(L)=State_Met%QL(I,J,L)*State_Met%AIRDEN(I,J,L)*1.d3* &
                    State_Met%BXHEIGHT(I,J,L)/CCLD(L)
          !test GEOS5 CLDLIQ appears to be too small, increase 50% to enhance Cloud Forcing
          !CLDLIQ(L)=1.5*State_Met%QL(I,J,L)*State_Met%AIRDEN(I,J,L)*1.d3*State_Met%BXHEIGHT(I,J,L)
          CLDICE(L)=State_Met%QI(I,J,L)*State_Met%AIRDEN(I,J,L)*1.d3* &
                    State_Met%BXHEIGHT(I,J,L)/CCLD(L)

          REL(L) = RCLDL3D(I,J,L)   ! in um (RRTMG limit 2.5 - 60)
          REL(L)= MAX(2.5d0,REL(L)) !  liquid particle radius
          REL(L)= MIN(55.d0,REL(L)) !  liquid particle radius
          !REI(L) = 15.0            ! fixed to this size (um) for now
          !REI(L) = 20.0            ! fixed to this size (um) for now
          !REI(L) = 30.0            ! fixed to this size (um) for now
          REI(L) = RCLDI3D(I,J,L)
          ! set to RRTMG limit
          REI(L)= MAX(5.0d0,REI(L))    !  ice particle radius
          REI(L)= MIN(140.d0,REI(L))   !  ice particle radius

       ENDDO

       call cldprop_swapm(State_Grid%NZ, CCLD, CLDICE, CLDLIQ, REI, REL, &
                          taucloud, taucloudl, taucloudi, &
                          ssacloudl, ssacloudi)

       YCOD3D(I,J,:) = taucloud(:,23)
       TCOD3D(I,J,:,1) = ssacloudl(:,23)
       TCOD3D(I,J,:,2) = taucloudl(:,23)+taucloudi(:,23)
       TCOD3D(I,J,:,3) = taucloudl(:,24)
       TCOD3D(I,J,:,4) = taucloudl(:,23)
       TCOD3D(I,J,:,5) = ssacloudi(:,23)
       TCOD3D(I,J,:,6) = ssacloudl(:,23)+ssacloudi(:,23)
       TCOD3D(I,J,:,7) = taucloudi(:,24)
       TCOD3D(I,J,:,8) = taucloudi(:,23)
       TCOD3D(I,J,:,9) = CCLD
       ZCOD(I,J)= ZCOD(I,J)+ SUM(taucloud(:,23))

    ENDDO
    ENDDO
    !$OMP END PARALLEL DO

    write(*,*)'LuoSSA',sum(TCOD3D(:,:,:,1))/size(TCOD3D(:,:,:,1)), &
                       sum(TCOD3D(:,:,:,5))/size(TCOD3D(:,:,:,5))

    !
    ! radf calculation
    IF(IFRADF.EQ.1) THEN
       CALL APM_RADFDRIV( Input_Opt,  State_Chm, State_Diag, &
                          State_Grid, State_Met, RC )
    ENDIF

    !APM2+IFSITE+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    IF(IFSITE.GT.0) THEN
       KYEAR = GET_YEAR()
       KMON = GET_MONTH()
       KDAY = GET_DAY()
       KHOUR = GET_HOUR()
       KMIN = GET_MINUTE()
       DO NSITE=1,MSITE
          IF(IFSITEOUT(NSITE).EQ.1) THEN
             KKOUT = 100*NSITE
             WRITE(KKOUT+70,105)KYEAR,KMON,KDAY,KHOUR,KMIN, &
                  (SITEOUT2D(NSITE,N),N=1,33)
             DO L=1,LOUT
                !WRITE(KKOUT+70,104)L,(SITEOUT2(NSITE,L,N),N=1,NTEMPOUT1+6),
                WRITE(KKOUT+70,104)L,(SITEOUT2(NSITE,L,N),N=1,NTEMPOUT1), &
                     (SITEOUT3(NSITE,L,N),N=1,11)
             ENDDO
             !WRITE(KKOUT+30,104)1,(SITEOUT2(NSITE,1,N), &
             !      N=NTEMPOUT1-22,NTEMPOUT1-15)

             IF(IFSITEADD.EQ.1) THEN ! additional output for size dis and others
                WRITE(KKOUT+10,105)KYEAR,KMON,KDAY,KHOUR,KMIN, &
                      YTOP(NSITE,1),YTOP(NSITE,2)
                WRITE(KKOUT+40,105)KYEAR,KMON,KDAY,KHOUR,KMIN
                !WRITE(KKOUT+50,105)KYEAR,KMON,KDAY,KHOUR,KMIN
                DO L=1,LOUT
                   !WRITE(KKOUT+10,104)L,(SITEOUT(NSITE,L,N),N=1,35)
                   WRITE(KKOUT+10,104)L,(SITEOUT(NSITE,L,N),N=1,59)
                   WRITE(KKOUT+40,104)L, &
                        (SITEOUT1(NSITE,L,N),N=1,NSO4+NSEA+NDSTB+NBCOC+NBCOC)
                   !WRITE(KKOUT+50,104)L,(AEROCOMOUT(NSITE,L,N),N=1,27)
                ENDDO
             ENDIF
          ENDIF

       ENDDO
    ENDIF
104 FORMAT(I9,500(1PE9.2))
105 FORMAT(I4,I3,I3,I3,I3,500(1PE10.2))

    ! Global budget
    KYEAR = GET_YEAR()
    KMON  = GET_MONTH()
    KDAY  = GET_DAY()
    KHOUR = GET_HOUR()
    KMIN  = GET_MINUTE()

    XX0 = SUM(EMITNH3(:,:,:))*14./17.*8.64d4 ! NH3 emission -- kg N/day

    CONDEP = 0.
    OXIDAT = 0.
    UPTAKE = 0.
    DRYDEP = 0.
    WETDEP = 0.

    !APM2+ENDIFSITE++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !$OMP PARALLEL DO        &
    !$OMP DEFAULT( SHARED )  &
    !$OMP PRIVATE( I, J )    &
    !$OMP SCHEDULE( DYNAMIC )
    DO J = 1, State_Grid%NY
    DO I = 1, State_Grid%NX
       IF(MODISOUT(I,J,5)>0.D0)THEN
          TEMPOUT(I,J,1,NTEMPOUT1-10)=TEMPOUT(I,J,1,NTEMPOUT1-10) &
                                     +MODISOUT(I,J,1)/MODISOUT(I,J,5)

          TEMPOUT(I,J,2,NTEMPOUT1-10)=TEMPOUT(I,J,2,NTEMPOUT1-10) &
                                     +MODISOUT(I,J,2)/MODISOUT(I,J,5)

          TEMPOUT(I,J,3,NTEMPOUT1-10)=TEMPOUT(I,J,3,NTEMPOUT1-10) &
                                     +MODISOUT(I,J,5)
       ENDIF

       IF(MODISOUT(I,J,6)>0.D0)THEN
          TEMPOUT(I,J,4,NTEMPOUT1-10)=TEMPOUT(I,J,4,NTEMPOUT1-10) &
                                     +MODISOUT(I,J,3)/MODISOUT(I,J,6)

          TEMPOUT(I,J,5,NTEMPOUT1-10)=TEMPOUT(I,J,5,NTEMPOUT1-10) &
                                     +MODISOUT(I,J,4)/MODISOUT(I,J,6)

          TEMPOUT(I,J,6,NTEMPOUT1-10)=TEMPOUT(I,J,6,NTEMPOUT1-10) &
                                     +MODISOUT(I,J,6)
       ENDIF
       MODISOUT(I,J,:)=0.D0
    ENDDO
    ENDDO
    !$OMP END PARALLEL DO

    NPOUTSTEPS = NPOUTSTEPS + 1
    ! reset YUV to 0
    YUV = 0.

    !DO N=APMIDS%id_SO4G,(APMIDS%id_AMINE+2)
    !   WRITE(110,*)N,MAXVAL(Spc(:,:,:,N)), &
    !               SUM(Spc(:,:,:,N))/SIZE(Spc(:,:,:,N))
    !   FLUSH(110)
    !ENDDO

    IF(IFATOM>0.5)THEN

       KYEAR = GET_YEAR()
       KMON  = GET_MONTH()
       KDAY  = GET_DAY()
       KHOUR = GET_HOUR()
       KMIN  = GET_MINUTE()

       IF((KMON==ATOMMONS.AND.KDAY>ATOMDAYS) .OR. &
          (KMON>ATOMMONS.AND.KMON<ATOMMONE)  .OR. &
          (KMON==ATOMMONE.AND.KDAY<ATOMDAYE))THEN
          WRITE(ATOMDATE(1:4),'(I4.4)')KYEAR
          WRITE(ATOMDATE(5:6),'(I2.2)')KMON
          WRITE(ATOMDATE(7:8),'(I2.2)')KDAY
          WRITE(ATOMDATE(9:10),'(I2.2)')KHOUR
          WRITE(ATOMDATE(11:12),'(I2.2)')KMIN
          CLOSE(100)
          OPEN(100,FILE='ATOM/ATOM'//ATOMDATE//'.bin', &
               ACCESS='DIRECT',FORM='UNFORMATTED',     &
               RECL=(95+APMIDS%id_AMINE+4-APMIDS%id_SO4G))

          RECID=1
          DO J=1,GRIDNUM
             DO L=1,40
                PRESS=GET_PCENTER(ATOMGRID(1,J),ATOMGRID(2,J),L)
                WRITE(100,REC=RECID) &
                      REAL(ATOMGRID(1,J)),REAL(ATOMGRID(2,J)),REAL(L),        &
                      REAL(PRESS),                                            &
                      REAL(State_Met%T(ATOMGRID(1,J),ATOMGRID(2,J),L)),       &
                      REAL(State_Met%RH(ATOMGRID(1,J),ATOMGRID(2,J),L)),      &
                      REAL(State_Met%CLDF(ATOMGRID(1,J),ATOMGRID(2,J),L)),    &
                      REAL(State_Met%AIRVOL(ATOMGRID(1,J),ATOMGRID(2,J),L)),  &
                      REAL(State_Met%AIRDEN(ATOMGRID(1,J),ATOMGRID(2,J),L)),  &
                      REAL(State_Met%AD(ATOMGRID(1,J),ATOMGRID(2,J),L)),      &
                      REAL(State_Met%BXHEIGHT(ATOMGRID(1,J),ATOMGRID(2,J),L)),&
                      REAL(State_Met%U(ATOMGRID(1,J),ATOMGRID(2,J),L)),       &
                      REAL(State_Met%V(ATOMGRID(1,J),ATOMGRID(2,J),L)),       &
                     (REAL(TEMPOUT(ATOMGRID(1,J),ATOMGRID(2,J),L,N)),N=1,77), &
                     (REAL(Spc(ATOMGRID(1,J),ATOMGRID(2,J),L,M)),             &
                           M=APMIDS%id_SO4G,(APMIDS%id_AMINE+2)),             &
                      REAL(Spc(ATOMGRID(1,J),ATOMGRID(2,J),L,APMIDS%id_SO2)), &
                      REAL(Spc(ATOMGRID(1,J),ATOMGRID(2,J),L,APMIDS%id_NH3)), &
                      REAL(Spc(ATOMGRID(1,J),ATOMGRID(2,J),L,APMIDS%id_NH4)), &
                      REAL(Spc(ATOMGRID(1,J),ATOMGRID(2,J),L,APMIDS%id_HNO3)),&
                      REAL(Spc(ATOMGRID(1,J),ATOMGRID(2,J),L,APMIDS%id_NIT))
                      RECID=RECID+1
             ENDDO
          ENDDO
          CLOSE(100)
       ENDIF
    ENDIF

    !$OMP PARALLEL DO         &
    !$OMP DEFAULT( SHARED )   &
    !$OMP PRIVATE( I, J, L )  &
    !$OMP SCHEDULE( DYNAMIC )
    DO L = 1, State_Grid%NZ
    DO J = 1, State_Grid%NY
    DO I = 1, State_Grid%NX
       IF(SUM(Spc(I,J,L,APMIDS%id_SO4BIN1:(APMIDS%id_SO4BIN1+NSO4-1))) &
          >1.d-30)THEN
          MWSIZE3D(I,J,L,1)=MWSIZE3D(I,J,L,1)/ &
             SUM(Spc(I,J,L,APMIDS%id_SO4BIN1:(APMIDS%id_SO4BIN1+NSO4-1)))
       ELSE
          MWSIZE3D(I,J,L,1)=RDRY(39)*GFTOT3D(I,J,L,1)
       ENDIF
       IF(SUM(Spc(I,J,L,APMIDS%id_SEABIN1:(APMIDS%id_SEABIN1+NSEA-1))) &
          >1.d-30)THEN
          MWSIZE3D(I,J,L,2)=MWSIZE3D(I,J,L,2)/ &
             SUM(Spc(I,J,L,APMIDS%id_SEABIN1:(APMIDS%id_SEABIN1+NSEA-1)))
       ELSE
          MWSIZE3D(I,J,L,2)=RSALT(19)*GFTOT3D(I,J,L,2)
       ENDIF
       IF(SUM(Spc(I,J,L,APMIDS%id_DSTBIN1:(APMIDS%id_DSTBIN1+NDSTB-1))) &
          >1.d-30)THEN
          MWSIZE3D(I,J,L,3)=MWSIZE3D(I,J,L,3)/ &
             SUM(Spc(I,J,L,APMIDS%id_DSTBIN1:(APMIDS%id_DSTBIN1+NDSTB-1)))
       ELSE
          MWSIZE3D(I,J,L,3)=RDRY(14)
       ENDIF
    ENDDO
    ENDDO
    ENDDO
    !$OMP END PARALLEL DO

    ! Free pointer
    Spc => NULL()

  END SUBROUTINE APM_DRIV
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: apm_radfdriv
!
! !DESCRIPTION: Subroutine APM\_RADFDRIV is the interface between GEOS-Chem/APM
!  and 1-column radiative forcing module
!
!  Originally written for CCCMa RF module by X. Ma & F. Yu, SUNY-Albany. 2011
!
!  Modified to drive RRTMG by F. Yu, 08/2012
!
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE APM_RADFDRIV( Input_Opt,  State_Chm, State_Diag, &
                           State_Grid, State_Met, RC )
!
! !USES:
!
    USE APM_INIT_MOD, ONLY : IFSITE,SITEOUT2D
    USE APM_INIT_MOD, ONLY : IFALBMODIS
    USE TIME_MOD,     ONLY : GET_DAY_OF_YEAR
    USE TIME_MOD,     ONLY : GET_YEAR, GET_MONTH
    USE TIME_MOD,     ONLY : GET_DAY,GET_HOUR,GET_MINUTE
    USE TIME_MOD,     ONLY : GET_NHMS
    USE TIME_MOD,     ONLY : GET_LOCALTIME
    USE APM_ALB_MOD,  ONLY : APM_ALB         ! get surface albedo from MODIS
    USE APM_RADF_MOD, ONLY : APM_RADF        ! radiation
    USE PRESSURE_MOD, ONLY : GET_PEDGE       !
    USE PRESSURE_MOD, ONLY : GET_PCENTER

    use rrtmg_sw_GCAPM, only : rrtmg_sw
    use parkind, only : im => kind_im, rb => kind_rb

    USE ErrCode_Mod
    Use Input_Opt_Mod,  ONLY : OptInput
    USE State_Chm_Mod,  ONLY : ChmState
    Use State_Diag_Mod, ONLY : DgnState
    USE State_Grid_Mod, ONLY : GrdState
    USE State_Met_Mod,  ONLY : MetState
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
    TYPE(GrdState), INTENT(IN)    :: State_Grid  ! Grid State object
    TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
    TYPE(DgnState), INTENT(INOUT) :: State_Diag  ! Diagnostics State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REVISION HISTORY:
!  2-6/2011 - X. Ma & F. Yu, SUNY-Albany
!  8/2012 - F. Yu, SUNY-Albany
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: I,J,L,N

    INTEGER :: KYEAR,KMON,KDAY,KHOUR,KMIN,NSITE

    LOGICAL, SAVE    :: FIRST = .TRUE.
    LOGICAL, SAVE    :: FIRST1 = .TRUE.

    INTEGER :: DAY_YR
    INTEGER :: MINIMUM,DDAY,DAY_READ
    INTEGER, SAVE :: AFLAG(366)

    INTEGER :: K, LL, IJLOOP
    INTEGER :: NDAY
    INTEGER :: LEV

    INTEGER, PARAMETER :: KWL=7         !! mxy MODIS spectral bands
    INTEGER, PARAMETER :: IAE2 = 2
    INTEGER, PARAMETER :: IAE1 = 1
    INTEGER, PARAMETER :: IAE0 = 0

    INTEGER, PARAMETER :: icld = 2   ! for RRTMG
    INTEGER, PARAMETER :: mxmol = 7   ! for RRTMG

    REAL           :: MAS             !!! AIR MASS OF BOX
    REAL(kind=rb)  :: SCOS            !!! COSINE OF SOLAR ZENITH ANGLE
    !REAL          :: CLDT            !!! COLUMN CLOUD COVER
    REAL(kind=rb)  :: PSURF, TSURF

    REAL(kind=rb)  :: AIRD(State_Grid%NZ)
    REAL(kind=rb)  :: P1D (State_Grid%NZ)
    REAL(kind=rb)  :: PE1D (0:State_Grid%NZ)
    REAL(kind=rb)  :: PDP (State_Grid%NZ)
    REAL(kind=rb)  :: T1D (State_Grid%NZ)
    REAL(kind=rb)  :: TE1D (0:State_Grid%NZ)
    REAL(kind=rb)  :: WKL (mxmol,State_Grid%NZ)
    REAL(kind=rb)  :: COLDRY(State_Grid%NZ) !!! dry air column density (mol/cm2)
    REAL(kind=rb)  :: CCLD(State_Grid%NZ)   !!! CLOUD COVER
    REAL(kind=rb)  :: CLDLIQ(State_Grid%NZ) !!! CLOUD LIQUID WATER CONTENT
    REAL(kind=rb)  :: CLDICE(State_Grid%NZ) !!! CLOUD ICE WATER CONTENT
    REAL(kind=rb)  :: REL(State_Grid%NZ)
    REAL(kind=rb)  :: REI(State_Grid%NZ)

    REAL           :: YCOD(State_Grid%NZ)

    REAL           :: MALB(State_Grid%NX,State_Grid%NY,KWL)

    REAL(kind=rb)  :: TEXT  (State_Grid%NZ,NBS,NTYP)
    REAL(kind=rb)  :: TOMGA (State_Grid%NZ,NBS,NTYP)
    REAL(kind=rb)  :: TG    (State_Grid%NZ,NBS,NTYP)
    REAL(kind=rb)  :: TABS  (State_Grid%NZ,NBL,NTYP)

    REAL(kind=rb)  :: EXT  (State_Grid%NZ,NBS)
    REAL(kind=rb)  :: OMGA (State_Grid%NZ,NBS)
    REAL(kind=rb)  :: G    (State_Grid%NZ,NBS)
    REAL(kind=rb)  :: ABSA (State_Grid%NZ,NBL)

    REAL(kind=rb)  :: YHTRC(State_Grid%NZ),YHTR(State_Grid%NZ)    !heating profiles
    REAL(kind=rb)  :: YHTRC0(State_Grid%NZ),YHTR0(State_Grid%NZ)

    REAL(kind=rb)  :: SALB(NBS)

    REAL(kind=rb)  :: TCST(NTYP),TCSB(NTYP),TFST(NTYP),TFSB(NTYP)
    REAL(kind=rb)  :: TFSA(NTYP)
    REAL(kind=rb)  :: CST,CSB,FST,FSB,FSA
    REAL(kind=rb)  :: CST0,CSB0,FST0,FSB0,FSA0

    REAL  :: THRS(State_Grid%NZ,NTYP),THRL(State_Grid%NZ,NTYP)
    REAL  :: TCLT(NTYP),TCLB(NTYP),TFUL(NTYP),TFDL(NTYP),TFLA(NTYP)

    REAL  :: HRS  (State_Grid%NZ),HRL  (State_Grid%NZ)
    REAL  :: HRS0 (State_Grid%NZ),HRL0 (State_Grid%NZ)
    REAL  :: CLT,CLB,FUL,FDL,FLA
    REAL  :: CLT0,CLB0,FUL0,FDL0,FLA0
    REAL  :: WCL, WCI
    REAL  :: YYAOD,BEXTL13,YVIS,RHL1
    REAL  :: YYAOD1,YYAOD3, YYAODT

    REAL*8,  PARAMETER   :: MAIR         = 28.966d-3           ! kg/mol
    REAL*8,  PARAMETER   :: XNUMOL_AIR   = 6.022d23 / MAIR    ! molec/kg

    ! For RRTMG
    REAL*8,  PARAMETER   :: YCO2  = 3.55d-4        ! mol/mol
    REAL*8,  PARAMETER   :: YN2O  = 3.20d-7        ! mol/mol
    REAL*8,  PARAMETER   :: YCH4  = 1.65d-6        ! mol/mol
    REAL*8,  PARAMETER   :: YO2   = 2.09d-1        ! mol/mol

    ! Make a pointer to the tracer array
    REAL*8, POINTER :: Spc(:,:,:,:)

    INTEGER :: IDO3                 !
    INTEGER :: id_O3                ! O3 in g.mol ??
    INTEGER :: id_CO                !

    RC = GC_SUCCESS

    LEV = State_Grid%NZ+1

    ! Point to Spc
    Spc => State_Chm%Species

    WRITE(6,*)'    - APM RADF_DRIV  '

    IF(FIRST)THEN

       WRITE(6,*)"KOUT1=",KOUT1, YWLS(KOUT1)
       WRITE(6,*)"KOUT2=",KOUT2, YWLS(KOUT2)
       WRITE(6,*)"KOUT3=",KOUT3, YWLS(KOUT3)

       ! Need to pass state objects here (bmy, 6/18/19)
       CALL APM_RADFOUT( Input_Opt, State_Chm, State_Grid, State_Met, RC )

       AFLAG=0
       DO i=1,366
          IF(MOD(i,8).eq.1) AFLAG(i)=1
       ENDDO
180    format(366i1)

       FIRST = .FALSE.
    ENDIF

    ! mxy+ read in MODIS surface albedo for restart
    DAY_YR = GET_DAY_OF_YEAR()
    print*,'DAY_OF_YR',DAY_YR

    KYEAR = GET_YEAR()
    KMON = GET_MONTH()
    KDAY = GET_DAY()
    KHOUR = GET_HOUR()

    print*,'KYEAR',KYEAR
    print*,'KHOUR',KHOUR

    print*,'FIRST1',FIRST1
    print*,'AFLAG',AFLAG(DAY_YR)
    IF(FIRST1) THEN   ! read data when restart
       IF(IFALBMODIS.EQ.1) THEN
          WRITE(6,*)"NEED TO READ SURFACE ALBEDO DATA WHEN RESTART"
          IF(AFLAG(DAY_YR) == 1 .and. KHOUR == 0) THEN
             CALL APM_ALB(DAY_YR,KYEAR,State_Grid%NX,State_Grid%NY,KWL,MALB)
             CALL GETALBDRR(State_Grid%NX,State_Grid%NY,MALB,State_Met)
          ENDIF
          IF(AFLAG(DAY_YR) == 0) THEN
             MINIMUM = 366
             DO I = 1,366
                IF(AFLAG(I) == 1) THEN
                   DDAY=ABS(DAY_YR-I)
                   IF(DDAY.LT.MINIMUM) THEN
                      MINIMUM = DDAY
                      DAY_READ = I
                   ENDIF
                ENDIF
             ENDDO
             print*,'READ DATA FROM THIS DAY',DAY_READ
             IF( KHOUR == 0) THEN
                CALL APM_ALB(DAY_READ,KYEAR,State_Grid%NX,State_Grid%NY,KWL,MALB)
                CALL GETALBDRR(State_Grid%NX,State_Grid%NY,MALB,State_Met)
             ENDIF
          ENDIF
       ENDIF
       FIRST1 = .FALSE.
    ENDIF

    IF(IFALBMODIS.EQ.1) THEN
       ! read data every 8 days
       IF(AFLAG(DAY_YR) == 1 .and. KHOUR == 0) THEN
          CALL APM_ALB(DAY_YR,KYEAR,State_Grid%NX,State_Grid%NY,KWL,MALB)
          CALL GETALBDRR(State_Grid%NX,State_Grid%NY,MALB,State_Met)
       ENDIF
    ELSE
       DO IWL=1,NBS
          ALBDRR(:,:,IWL) = State_Met%ALBD(:,:)
       ENDDO
    ENDIF

    !$OMP PARALLEL DO                                   &
    !$OMP DEFAULT( SHARED )                             &
    !$OMP PRIVATE( I, J, L, LL, IJLOOP,IWL,ITYP,NSITE)  &
    !$OMP PRIVATE( TEXT, TOMGA,TG,TABS )                &
    !$OMP PRIVATE( EXT, OMGA,G,ABSA )                   &
    !$OMP PRIVATE( YHTRC,YHTR,YHTRC0,YHTR0)             &
    !$OMP PRIVATE( TCST, TCSB, TFST, TFSB, TFSA, THRS ) &
    !$OMP PRIVATE( TCLT, TCLB, TFUL, TFDL, TFLA, THRL ) &
    !$OMP PRIVATE( CST, CSB, FST, FSB, FSA, HRS )       &
    !$OMP PRIVATE( CST0, CSB0, FST0, FSB0, FSA0, HRS0 ) &
    !$OMP PRIVATE( CLT , CLB , FUL, FDL, FLA, HRL  )    &
    !$OMP PRIVATE( CLT0, CLB0, FUL0, FDL0, FLA0, HRL0 ) &
    !$OMP PRIVATE( WCL, WCI )                           &
    !$OMP PRIVATE( SALB )                               &
    !$OMP PRIVATE( YYAOD,BEXTL13,YVIS,RHL1 )            &
    !$OMP PRIVATE( YYAOD1 )                             &
    !$OMP PRIVATE( YYAODT )                             &
    !$OMP PRIVATE( PSURF,TSURF,SCOS )                   &
    !$OMP PRIVATE( AIRD,P1D,PE1D,PDP,T1D,TE1D)          &
    !$OMP PRIVATE( COLDRY,WKL)                          &
    !$OMP PRIVATE( CCLD,CLDLIQ,CLDICE )                 &
    !$OMP PRIVATE( REL,REI)                             &
    !$OMP PRIVATE( YCOD )                               &
    !$OMP SCHEDULE( DYNAMIC )
    DO J = 1, State_Grid%NY
    DO I = 1, State_Grid%NX

       ! Initialize parameters for 1-column rad calculation
       AOD(I,J)      =1.d-20
       TAOD(I,J,:)   =1.d-20

       AAOD(I,J)     =1.d-20
       AODOUT1(I,J)  =1.d-20
       AAODOUT1(I,J) =1.d-20
       AODOUT3(I,J)  =1.d-20
       AAODOUT3(I,J) =1.d-20

       ! Optical properties: RRTMG need AOD instead of ext, now EXT is layer AOD
       DO L = 1, State_Grid%NZ
          DO IWL = 1, NBS
             EXT(L,IWL) = ZBEXT3D(I,J,L,NWLS(IWL))* &
                          State_Met%BXHEIGHT(I,J,L)*100.
             OMGA(L,IWL)= ZW3D(I,J,L,NWLS(IWL))
             G(L,IWL)   = ZG3D(I,J,L,NWLS(IWL))
          ENDDO
          DO ITYP = 1,NTYP
             DO IWL = 1, NBS
                TEXT(L,IWL,ITYP) = YBEXT3D(I,J,L,ITYP,NWLS(IWL))* &
                                   State_Met%BXHEIGHT(I,J,L)*100.
                TOMGA(L,IWL,ITYP)= YW3D(I,J,L,ITYP,NWLS(IWL))
                TG(L,IWL,ITYP)   = YG3D(I,J,L,ITYP,NWLS(IWL))
             ENDDO
          ENDDO
          ! LW
          DO IWL=1,NWL
             ABSA(L,IWL)=ZBABS3D(I,J,L,IWL)
          ENDDO

          DO ITYP = 1,NTYP
             IF(ITYP.EQ.3)THEN !dust
                DO IWL=1,NBL
                   TABS(L,IWL,ITYP)=ZBABS3D(I,J,L,IWL)
                ENDDO
             ELSE
                DO IWL=1,NBL
                   TABS(L,IWL,ITYP)=0.
                ENDDO
             ENDIF
          ENDDO

          ! instantanous AOD
          YYAOD = State_Met%BXHEIGHT(I,J,L)*100.*ZBEXT3D(I,J,L,KOUT2)  !WL(KOUT2) nm

          AOD(I,J) =  AOD(I,J)+ YYAOD
          AAOD(I,J)= AAOD(I,J)+ YYAOD*(1.-ZW3D(I,J,L,KOUT2))  !WL(KOUT2) nm

          YYAOD1       = State_Met%BXHEIGHT(I,J,L)*100.*ZBEXT3D(I,J,L,KOUT1)  !
          AODOUT1(I,J) =  AODOUT1(I,J)+ YYAOD1
          AAODOUT1(I,J)= AAODOUT1(I,J)+ YYAOD1*(1.-ZW3D(I,J,L,KOUT1))  !

          YYAOD3       = State_Met%BXHEIGHT(I,J,L)*100.*ZBEXT3D(I,J,L,KOUT3)  !
          AODOUT3(I,J) =  AODOUT3(I,J)+ YYAOD3
          AAODOUT3(I,J)= AAODOUT3(I,J)+ YYAOD3*(1.-ZW3D(I,J,L,KOUT3))  !

          ! AOD
          DO ITYP = 1,NTYP
             YYAODT = State_Met%BXHEIGHT(I,J,L)*100.* &
                      YBEXT3D(I,J,L,ITYP,KOUT2)      !WL(KOUT2) nm
             TAOD(I,J,ITYP) = TAOD(I,J,ITYP)+ &
                              State_Met%BXHEIGHT(I,J,L)*100.* &
                              YBEXT3D(I,J,L,ITYP,KOUT2)      !WL(KOUT2) nm
             ZTAOD(I,J,ITYP) = ZTAOD(I,J,ITYP)+ &
                               State_Met%BXHEIGHT(I,J,L)*100.* &
                               YBEXT3D(I,J,L,ITYP,KOUT2)      !WL(KOUT2) nm
          ENDDO
          !cumulative AOD for average
          ZAOD(I,J)= ZAOD(I,J)+ YYAOD
          CODGC(I,J) = CODGC(I,J) + &
                       (State_Met%TAUCLW(I,J,L)+State_Met%TAUCLI(I,J,L)) &
                       *State_Met%CLDF(I,J,L)
          ZCODGC(I,J)= ZCODGC(I,J) + &
                       (State_Met%TAUCLW(I,J,L)+State_Met%TAUCLI(I,J,L)) &
                       *State_Met%CLDF(I,J,L)

          ZAODOUT1(I,J)= ZAODOUT1(I,J)+ YYAOD1  !380 nm
          ZAODOUT3(I,J)= ZAODOUT3(I,J)+ YYAOD3  !870 nm

          IF(State_Met%CLDFRC(I,J).LE.0.25) THEN
             ZAOD25(I,J)=ZAOD25(I,J) + YYAOD
             IF(L.EQ.1) NAOD25(I,J)=NAOD25(I,J) + 1
          ENDIF
          IF(State_Met%CLDFRC(I,J).LE.0.50) THEN
             ZAOD50(I,J)=ZAOD50(I,J) + YYAOD
             IF(L.EQ.1) NAOD50(I,J)=NAOD50(I,J) + 1
          ENDIF
          IF(State_Met%CLDFRC(I,J).LE.0.75) THEN
             ZAOD75(I,J)=ZAOD75(I,J) + YYAOD
             IF(L.EQ.1) NAOD75(I,J)=NAOD75(I,J) + 1
          ENDIF

          ZAAOD(I,J)= ZAAOD(I,J)+ State_Met%BXHEIGHT(I,J,L)*100.* &
                      ZBEXT3D(I,J,L,KOUT2)*(1.-ZW3D(I,J,L,KOUT2))  !WL(KOUT2) nm

          ZAAODOUT1(I,J)= ZAAODOUT1(I,J)+ State_Met%BXHEIGHT(I,J,L)*100.* &
                          ZBEXT3D(I,J,L,KOUT1)*(1.-ZW3D(I,J,L,KOUT1))  !380 nm

          ZAAODOUT3(I,J)= ZAAODOUT3(I,J)+ State_Met%BXHEIGHT(I,J,L)*100.* &
                          ZBEXT3D(I,J,L,KOUT3)*(1.-ZW3D(I,J,L,KOUT3))  !870 nm

          ZABS(I,J)= ZABS(I,J)+ State_Met%BXHEIGHT(I,J,L)*100.* &
                     ZBABS3D(I,J,L,6)  !LW, 1150 nm

       ENDDO   ! L=1,State_Grid%NZ

       ! visibility
       ! mean bext L1_L3 (H<387 m)
       BEXTL13 = (ZBEXT3D(I,J,1,6)*State_Met%BXHEIGHT(I,J,1) + &   ! 550 nm
                  ZBEXT3D(I,J,2,6)*State_Met%BXHEIGHT(I,J,2) + &
                  ZBEXT3D(I,J,3,6)*State_Met%BXHEIGHT(I,J,3))/ &
         (State_Met%BXHEIGHT(I,J,1)+State_Met%BXHEIGHT(I,J,2)+ &
          State_Met%BXHEIGHT(I,J,3))   ! in cm-1

       YVIS = 3.912/(3.9d-2 + BEXTL13*1.d5)  ! in km

       ZVIS(I,J) = ZVIS(I,J) + YVIS   ! accumulated for average

       RHL1 = State_Met%RH(I,J,1)/100. ! RH at surface layer
       IF(RHL1.LT.0.95.and.YVIS.LE.5.0) THEN
          THAZ(I,J)=THAZ(I,J) + 1.0
       ENDIF
       IF(RHL1.GT.0.95.and.YVIS.LE.2.0) THEN
          TFOG(I,J)=TFOG(I,J) + 1.0
       ENDIF

       ! Profiles for RF
       IJLOOP = ( (J-1) * State_Grid%NX ) + I
       SCOS   = State_Met%SUNCOS(I,J)

       PE1D(0) = GET_PEDGE(I,J,1)        ! P surface   [MB]
       TE1D(0) = State_Met%TSKIN(I,J)              ! T surface (K)
       TSURF   = TE1D(0)

       DO IWL = 1, NBS
          SALB(IWL) = ALBDRR(I,J,IWL)       ! Albedo surface
       ENDDO

       !IF(I.EQ.60.and.J.EQ.33) THEN
       !   WRITE(1002,99)DAY_YR,KHOUR,State_Grid%NZ,SCOS, &
       !      AOD(I,J),State_Met%CLDFRC(I,J),(SALB(IWL),IWL=1,NBS)
       !99 FORMAT(I4,I3,I3,F9.3,100(1PE9.2))
       !ENDIF

       DO L=1, State_Grid%NZ
          AIRD(L)= State_Met%AIRDEN(I,J,L)*XNUMOL_AIR* &
                   State_Met%BXHEIGHT(I,J,L)*1.d-4 !box column air mole #/cm2
                                                   ! AIRDEN unit: kg/m3
          P1D(L) = GET_PCENTER(I,J,L)      ! Pressure center    [MB]
          PE1D(L)= GET_PEDGE(I,J,L+1)     ! Pressure upper edge   [MB]

          PDP(L) = PE1D(L-1)-PE1D(L)

          T1D(L)= State_Met%T(I,J,L)     ! Temperature center [K]
          IF(L.EQ.State_Grid%NZ) THEN
             TE1D(L) = State_Met%T(I,J,L) + &
                       0.5*(State_Met%T(I,J,L)-State_Met%T(I,J,L-1))
          ELSE
             TE1D(L) = 0.5*(State_Met%T(I,J,L) + State_Met%T(I,J,L+1))
          ENDIF

          ! Gases:  in #/cm2 layer column density
          ! (CO2, N2O, CH4, and O2 use fixed mixing ratio for now)
          ! H2O  convert specific humidity from g/kg to g/g and then to volume
          ! mixing ratio and then to #/cm2
          WKL(1,L) = State_Met%SPHU(I,J,L)*1.d-3*MAIR/18.0d-3 * AIRD(L)

          WKL(2,L)=  YCO2 * AIRD(L)

          ! convert O3 and CO from kg/box to g/g and then to volume mixing ratio
          ! and then to #/cm2
          WKL(3,L)= Spc(I,J,L,id_O3)/State_Met%AD(I,J,L) &
                    *MAIR/48.0d-3 * AIRD(L)

          WKL(4,L)= YN2O * AIRD(L)
          ! CO
          WKL(5,L)= Spc(I,J,L,id_CO)/State_Met%AD(I,J,L) &
                    *MAIR/28.0d-3 * AIRD(L)

          WKL(6,L)= YCH4 * AIRD(L)

          WKL(7,L)= YO2  * AIRD(L)

          ! dry air column density (mol/cm2)
          COLDRY(L) = AIRD(L)- WKL(1,L)

          CCLD(L)= MAX(1.d-3,State_Met%CLDF(I,J,L))  ! Cloud cover [Unitless]

          ! convert water content from kg/kg to g/m3 and then to g/m2
          CLDLIQ(L)=State_Met%QL(I,J,L)*State_Met%AIRDEN(I,J,L)*1.d3* &
                    State_Met%BXHEIGHT(I,J,L)/CCLD(L)
          !test GEOS5 CLDLIQ appears to be too small, increase 50% to enhance Cloud Forcing
          !CLDLIQ(L)=1.5*State_Met%QL(I,J,L)*State_Met%AIRDEN(I,J,L)*1.d3*State_Met%BXHEIGHT(I,J,L)
          CLDICE(L)=State_Met%QI(I,J,L)*State_Met%AIRDEN(I,J,L)*1.d3* &
                    State_Met%BXHEIGHT(I,J,L)/CCLD(L)

          REL(L) = RCLDL3D(I,J,L)   ! in um (RRTMG limit 2.5 - 60)
          REL(L)= MAX(2.5d0,REL(L)) !  liquid particle radius
          REL(L)= MIN(55.d0,REL(L)) !  liquid particle radius
          !REI(L) = 15.0         ! fixed to this size (um) for now
          !REI(L) = 20.0         ! fixed to this size (um) for now
          !REI(L) = 30.0         ! fixed to this size (um) for now
          REI(L) = RCLDI3D(I,J,L)
          ! set to RRTMG limit
          REI(L)= MAX(5.0d0,REI(L))   !  ice particle radius
          REI(L)= MIN(140.d0,REI(L))  !  ice particle radius
       ENDDO

       !IF(I.EQ.60.and.J.EQ.33) THEN
       !   WRITE(1002,101) (P1D(L),L=1,State_Grid%NZ)
       !   WRITE(1002,101) (PE1D(L),L=0,State_Grid%NZ)
       !   WRITE(1002,101) (PDP(L),L=1,State_Grid%NZ)
       !   WRITE(1002,101) (T1D(L),L=1,State_Grid%NZ)
       !   WRITE(1002,101) (TE1D(L),L=0,State_Grid%NZ)
       !
       !   DO N = 1, mxmol
       !      WRITE(1002,101) (WKL(N,L),L=1,State_Grid%NZ)
       !   ENDDO
       !   WRITE(1002,101) (COLDRY(L),L=1,State_Grid%NZ)
       !
       !   WRITE(1002,101) (CCLD(L),L=1,State_Grid%NZ)
       !   WRITE(1002,101) (CLDLIQ(L),L=1,State_Grid%NZ)
       !   WRITE(1002,101) (CLDICE(L),L=1,State_Grid%NZ)
       !   WRITE(1002,101) (REL(L),L=1,State_Grid%NZ)
       !   WRITE(1002,101) (REI(L),L=1,State_Grid%NZ)
       !
       !   DO IWL=1,NBS
       !      WRITE(1002,101) (EXT(L,IWL),L=1,State_Grid%NZ)   !layer AOD
       !      WRITE(1002,101) (OMGA(L,IWL),L=1,State_Grid%NZ)
       !      WRITE(1002,101) (G(L,IWL),L=1,State_Grid%NZ)
       !   ENDDO
       !ENDIF
101    FORMAT(100(1PE9.2))
104    FORMAT(I3,100(1PE9.2))

       IF(SCOS.GT.0.) THEN
          CALL rrtmg_sw(I,J,2,icld,State_Grid%NZ,DAY_YR, SCOS,   &
                        PDP,P1D,T1D,PE1D,TE1D,TSURF,COLDRY, WKL, &
                        CCLD, CLDICE, CLDLIQ, REI, REL,          &
                        SALB,EXT, OMGA, G,                       &
                        YHTRC,YHTR,YHTRC0,YHTR0,                 &
                        CST,FST,CSB,FSB,CST0,FST0,CSB0,FSB0,     &
                        TEXT, TOMGA, TG, TCST,TCSB,TFST,TFSB)
       ELSE
          CST    = 0.
          FST    = 0.
          CSB    = 0.
          FSB    = 0.
          CST0   = 0.
          FST0   = 0.
          CSB0   = 0.
          FSB0   = 0.

          TCST   = 0.
          TFST   = 0.
          TCSB   = 0.
          TFSB   = 0.

          YHTRC  = 0.
          YHTR   = 0.
          YHTRC0 = 0.
          YHTR0  = 0.
       ENDIF

       FUL  = 0.
       FDL  = 0.
       FLA  = 0.
       HRL  = 0.
       FUL0 = 0.
       FDL0 = 0.
       FLA0 = 0.
       HRL0 = 0.

       FSA  = FST - FSB
       FSA0 = FST0 - FSB0

       TFSA(:) = TFST(:) - TFSB(:)

       ! Each aerosol type
       !DO ITYP = 1,NTYP
       !   CALL APM_RADF(J,I,State_Grid%NZ,LEV,NBS,NBL,                 &
       !              IAE2,TEXT(:,:,ITYP),TOMGA(:,:,ITYP),TG(:,:,ITYP), &
       !              TABS(:,:,ITYP),                                   &
       !              GT,SALB,SCOS,T1D,H2O,O3,P1D,AIRD,                 &
       !              CCLD,CLDLIQ,CLDICE,RCLDL,YCOD,                    &
       !              WCL,WCI,                                          &
       !              TCST(ITYP),TCSB(ITYP),TFST(ITYP),TFSB(ITYP),      &
       !              TFSA(ITYP),THRS(:,ITYP),                          &
       !              TCLT(ITYP),TCLB(ITYP),TFUL(ITYP),TFDL(ITYP),      &
       !              TFLA(ITYP),THRL(:,ITYP) )
       !ENDDO

       THRS = 0.
       TCLT = 0.
       TCLB = 0.
       TFUL = 0.
       TFDL = 0.
       TFLA = 0.
       THRL = 0.

       ! Total
       !CALL APM_RADF(J,I,State_Grid%NZ,LEV,NBS,NBL,                    &
       !              IAE1,EXT,OMGA,G,ABSA,                             &
       !              GT,SALB,SCOS,T1D,H2O,O3,P1D,AIRD,                 &
       !              CCLD,CLDLIQ,CLDICE,RCLDL,YCOD,                    &
       !              WCL, WCI,                                         &
       !              CST,CSB,FST,FSB,FSA,HRS,                          &
       !              CLT,CLB,FUL,FDL,FLA,HRL)
       !
       !print*,'ALBEDO,GT',IAE0,I,J,ALBEDO,GT,SCOS

       ! No aerosol
       !CALL APM_RADF(J,I,State_Grid%NZ,LEV,NBS,NBL,                    &
       !              IAE0,EXT0,OMGA0,G0,ABSA0,                         &
       !              GT,SALB,SCOS,T1D,H2O,O3,P1D,AIRD,                 &
       !              CCLD,CLDLIQ,CLDICE,RCLDL,YCOD,                    &
       !              WCL, WCI,                                         &
       !              CST0,CSB0,FST0,FSB0,FSA0,HRS0,                    &
       !              CLT0,CLB0,FUL0,FDL0,FLA0,HRL0)

       YHTRC3D(I,J,:) =YHTRC(:)    ! clear sky, with aerosol
       YHTR3D(I,J,:)  =YHTR(:)     ! all sky, with aerosol
       YHTRC03D(I,J,:)=YHTRC0(:)   ! clear sky, no aerosol
       YHTR03D(I,J,:) =YHTR0(:)    ! all sky, no aerosol

       DO ITYP = 1,NTYP
          ZTCST(I,J,ITYP)=ZTCST(I,J,ITYP) + (TCST(ITYP)-CST0)
          ZTCSB(I,J,ITYP)=ZTCSB(I,J,ITYP) + (TCSB(ITYP)-CSB0)
          ZTFST(I,J,ITYP)=ZTFST(I,J,ITYP) + (TFST(ITYP)-FST0)
          ZTFSB(I,J,ITYP)=ZTFSB(I,J,ITYP) + (TFSB(ITYP)-FSB0)
          ZTFSA(I,J,ITYP)=ZTFSA(I,J,ITYP) + (TFSA(ITYP)-FSA0)
       ENDDO

       ZCLDF(I,J)=ZCLDF(I,J) + State_Met%CLDFRC(I,J)
       ZALB(I,J)=ZALB(I,J) + ALBDRR(I,J,10)

       DO IWL = 1,NBS
          ZMALB(I,J,IWL)=ZMALB(I,J,IWL) + SALB(IWL)
       ENDDO

       ZCST(I,J) = ZCST(I,J) + (CST-CST0 )
       ZCSB(I,J) = ZCSB(I,J) + (CSB-CSB0 )
       ZFST(I,J) = ZFST(I,J) + (FST-FST0 )
       ZCLD(I,J) = ZCLD(I,J) + (FST-CST  )
       ZCLD0(I,J)= ZCLD0(I,J)+ (FST0-CST0)
       ZFSB(I,J) = ZFSB(I,J) + (FSB-FSB0 )
       ZFSA(I,J) = ZFSA(I,J) + (FSA-FSA0 )

       !if(i.eq.36.and.j.eq.25) then
       !   print*,'CST,CST0,CST-CST0',CST,CST0,CST-CST0
       !   print*,'FST,FST0,FST-FST0',FST,FST0,FST-FST0
       !   print*,'FSA,FSA0,FSA-FSA0',FSA,FSA0,FSA-FSA0
       !   print*,'FST-CST,FST0-CST0',FST-CST,FST0-CST0
       !endif
       if(i.eq.60.and.j.eq.33) then
          WRITE(1003,109) DAY_YR, AOD(I,J),State_Met%CLDFRC(I,J)
          WRITE(1003,107)'CST,CST0,CST-CST0',CST,CST0,CST-CST0
          WRITE(1003,107)'FST,FST0,FST-FST0',FST,FST0,FST-FST0
          WRITE(1003,107)'FSA,FSA0,FSA-FSA0',FSA,FSA0,FSA-FSA0
          WRITE(1003,107)'FST-CST,FST0-CST0',FST-CST,FST0-CST0
          DO ITYP = 1,NTYP
             WRITE(1003,107)'TCST-CST0,TFST-FST0', &
                             TCST(ITYP)-CST0,TFST(ITYP)-FST0
          ENDDO
109       FORMAT(I4,2x,10(F10.3))
107       FORMAT(A19,2x,10(F10.2))
       endif

       !ZCLT(I,J)=ZCLT(I,J) + (CLT-CLT0 )
       !ZCLB(I,J)=ZCLB(I,J) + (CLB-CLB0 )
       !ZFUL(I,J)=ZFUL(I,J) + (FUL0-FUL )
       !ZFDL(I,J)=ZFDL(I,J) + (FDL-FDL0 )
       !ZFLA(I,J)=ZFLA(I,J) + (FLA-FLA0 )

       ZWCL(I,J)=ZWCL(I,J) + WCL
       ZWCI(I,J)=ZWCI(I,J) + WCI

       IF((FST-FST0).gt.1000.0.or.(FST-FST0).lt.-2000.0) THEN
          write(1001,*) i,j
          !DO ITYP=1,NTYP
          ! write(1001,*) ITYP,TCST(ITYP),CST0
          ! DO L=1,State_Grid%NZ
          !  write(1001,*) L,TEXT(L,1,ITYP),TOMGA(L,1,ITYP),TG(L,1,ITYP)
          ! write(1001,*) L,TEXT(L,2,ITYP),TOMGA(L,2,ITYP),TG(L,2,ITYP)
          !  write(1001,*) L,TEXT(L,3,ITYP),TOMGA(L,3,ITYP),TG(L,3,ITYP)
          !  write(1001,*) L,TEXT(L,4,ITYP),TOMGA(L,4,ITYP),TG(L,4,ITYP)
          ! ENDDO
          !ENDDO

          write(1001,*) "EXT",EXT
          write(1001,*) "OMGA",OMGA
          write(1001,*) "G",G
          WRITE(1001,*) 'rf, claer sky toa'
          WRITE(1001,103)CST,CST0,CST-CST0
          WRITE(1001,*) '    clear sky surface'
          WRITE(1001,103)CSB,CSB0,CSB-CSB0
          WRITE(1001,*) '    all sky toa'
          WRITE(1001,103)FST,FST0,FST-FST0
          WRITE(1001,*) '    absorption'
          WRITE(1001,103)FSA,FSA0,FSA-FSA0
       ENDIF

102    FORMAT(100(1PE9.2))
103    FORMAT(3F9.4)

       !APM2+IFSITE++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
       ! output at selected sites
       IF(IFSITE.GT.0) THEN
          IF(IFOUTIJ(I,J).EQ.1) THEN
             NSITE = SITEID(I,J)
             SITEOUT2D(NSITE,1) = AODOUT1(I,J)
             SITEOUT2D(NSITE,2) = AOD(I,J)   !WL(KOUT2) AOD
             SITEOUT2D(NSITE,3) = AODOUT3(I,J)
             SITEOUT2D(NSITE,4) = -LOG10(AOD(I,J)/AODOUT1(I,J)) &
                                  /LOG10(YWLS(KOUT2)/YWLS(KOUT1))            !WL(KOUT1)_WL(KOUT2) Angstrom Exponent
             SITEOUT2D(NSITE,5) = -LOG10(AOD(I,J)/AODOUT3(I,J)) &
                                  /LOG10(YWLS(KOUT2)/YWLS(KOUT3))            !WL(KOUT2)_WL(KOUT3) Angstrom Exponent
             SITEOUT2D(NSITE,6) = AAODOUT1(I,J)
             SITEOUT2D(NSITE,7) = AAOD(I,J)  !WL(KOUT2) nm AAOD
             SITEOUT2D(NSITE,8) = AAODOUT3(I,J)
             SITEOUT2D(NSITE,9) = -LOG10(AAOD(I,J)/AAODOUT1(I,J)) &
                                  /LOG10(YWLS(KOUT2)/YWLS(KOUT1))            !WL(KOUT1)_WL(KOUT2) Angstrom Exponent
             SITEOUT2D(NSITE,10) = -LOG10(AAOD(I,J)/AAODOUT3(I,J)) &
                                   /LOG10(YWLS(KOUT2)/YWLS(KOUT3))           !WL(KOUT1)_WL(KOUT2) Angstrom Exponent

             SITEOUT2D(NSITE,11) = State_Met%CLDFRC(I,J)
             SITEOUT2D(NSITE,12) = CST-CST0  !clear sky aerosol forcing
             SITEOUT2D(NSITE,13) = FST-FST0  !full sky aerosol forcing
             DO ITYP=1,NTYP
                SITEOUT2D(NSITE,13+ITYP) = TFST(ITYP)-FST0   !full sky 5 type aerosol forcing
             ENDDO
             !SITEOUT2D(NSITE,13+NTYP+1)=YVIS
             !SITEOUT2D(NSITE,13+NTYP+2)=State_Met%PREACC(I,J)
             !SITEOUT2D(NSITE,19) = SUM(YCOD)   !cloud OD
             SITEOUT2D(NSITE,19) = CODGC(I,J)
             SITEOUT2D(NSITE,20) = FST-CST   !cloud forcing with aerosol
             SITEOUT2D(NSITE,22) = FST0-CST0   !cloud forcing without aerosol

             SITEOUT2D(NSITE,22)=State_Met%PRECTOT(I,J)
             SITEOUT2D(NSITE,23)=State_Met%PRECCON(I,J)
             SITEOUT2D(NSITE,24)=YVIS
             DO ITYP=1,NTYP
                SITEOUT2D(NSITE,24+ITYP)=TAOD(I,J,ITYP)
             ENDDO
             SITEOUT2D(NSITE,30)=FSB
             SITEOUT2D(NSITE,31)=FSB0
             SITEOUT2D(NSITE,32)=CSB
             SITEOUT2D(NSITE,33)=CSB0
          ENDIF
       ENDIF
       !APM2+ENDIFSITE+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ENDDO
    ENDDO
    !$OMP END PARALLEL DO

    FLUSH(890)
    FLUSH(891)
    FLUSH(1001)

    ! Clear the pointer
    Spc => NULL()

  END SUBROUTINE APM_RADFDRIV
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: apm_radfout
!
! !DESCRIPTION: Subroutine APM\_RADFOUT is the for radf output
!
!  Written by X. Ma & F. Yu, SUNY-Albany
!
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE APM_RADFOUT( Input_Opt, State_Chm, State_Grid, State_Met, RC )
!
! !USES:
!
    USE TIME_MOD,       ONLY : GET_YEAR, GET_MONTH
    USE TIME_MOD,       ONLY : GET_DAY,  GET_HOUR,  GET_MINUTE
    USE ErrCode_Mod
    Use Input_Opt_Mod,  ONLY : OptInput
    USE State_Chm_Mod,  ONLY : ChmState
    USE State_Grid_Mod, ONLY : GrdState
    USE State_Met_Mod,  ONLY : MetState
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
    TYPE(GrdState), INTENT(IN)    :: State_Grid  ! Grid State object
    TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REVISION HISTORY:
!  2-6/2011 - X. Ma & F. Yu, SUNY-Albany
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER          :: I,J,ITYP,IWL,RECID

    INTEGER          :: KYEAR,KMON,KDAY,KHOUR,KMIN

    LOGICAL, SAVE    :: FIRST = .TRUE.

    REAL*8           :: XCOUNT
    REAL*8           :: ZFSTAVE,ZCSTAVE,TAREA,YAREA
    REAL*8           :: ZCLDAVE,ZCLD0AVE
    REAL*8           :: ZFULAVE     ! longwave TOA
    REAL*8           :: ZFSAAVE, ZTFSTAVE(NTYP)
    ! fraction of time w/ CLDFR<25%, 50%, 75%
    REAL             :: FAOD25(State_Grid%NX,State_Grid%NY)
    REAL             :: FAOD50(State_Grid%NX,State_Grid%NY)
    REAL             :: FAOD75(State_Grid%NX,State_Grid%NY)
    REAL             :: ZAEEXT(State_Grid%NX,State_Grid%NY)
    REAL             :: ZAEABS(State_Grid%NX,State_Grid%NY)
    REAL             :: ZAEEXTOUT3(State_Grid%NX,State_Grid%NY)
    REAL             :: ZAEABSOUT3(State_Grid%NX,State_Grid%NY)
    CHARACTER(LEN=6) :: YYYYMM='200001'

    ! Initialize
    RC = GC_SUCCESS

    WRITE(6,*)'    - APM RADFOUT  '

    IF(FIRST)THEN
       ZTCST = 0.
       ZTCSB = 0.
       ZTFST = 0.
       ZTFSB = 0.
       ZTFSA = 0.
       ZTCLT = 0.
       ZTCLB = 0.
       ZTFUL = 0.
       ZTFDL = 0.
       ZTFLA = 0.

       ZTAOD = 0.

       ZMALB = 0.
       ZCLDF = 0.
       ZALB  = 0.

       ZCST  = 0.
       ZCSB  = 0.
       ZFST  = 0.
       ZCLD  = 0.
       ZCLD0 = 0.
       ZFSB  = 0.
       ZFSA  = 0.
       !ZCLT = 0.
       !ZCLB = 0.
       !ZFUL = 0.
       !ZFDL = 0.
       !ZFLA = 0.

       ZWCL  = 0.
       ZWCI  = 0.

       ZAOD  = 0.
       ZCOD  = 0.
       ZCODGC= 0.
       ZAAOD = 0.
       ZABS  = 0.

       ZAODOUT1  = 0.
       ZAAODOUT1 = 0.

       ZAOD25 = 0.
       ZAOD50 = 0.
       ZAOD75 = 0.

       ZVIS = 0.
       THAZ = 0.
       TFOG = 0.

       FIRST = .FALSE.
       RETURN
    ENDIF

    KYEAR = GET_YEAR()
    KMON  = GET_MONTH()
    KDAY  = GET_DAY()
    KHOUR = GET_HOUR()
    KMIN  = GET_MINUTE()

    WRITE(YYYYMM(1:4),'(I4.4)')KYEAR
    WRITE(YYYYMM(5:6),'(I2.2)')KMON

    CLOSE(686)
    CLOSE(687)
    CLOSE(688)
    CLOSE(689)
    CLOSE(690)

    CLOSE(785)

    CLOSE(786)
    CLOSE(787)
    CLOSE(788)
    CLOSE(789)
    CLOSE(790)

    CLOSE(811)
    CLOSE(812)
    CLOSE(813)
    CLOSE(814)
    CLOSE(815)
    CLOSE(816)

    CLOSE(821)
    CLOSE(822)
    CLOSE(823)

    CLOSE(883)
    CLOSE(884)
    CLOSE(885)

    CLOSE(886)
    CLOSE(887)
    CLOSE(888)
    CLOSE(588)
    CLOSE(589)
    CLOSE(889)
    CLOSE(890)
    CLOSE(590)
    CLOSE(591)
    CLOSE(891)
    CLOSE(892)
    CLOSE(893)
    CLOSE(894)
    CLOSE(895)

    OPEN(686,FILE='radf/radf'//YYYYMM//'CLT.out',ACCESS='DIRECT',    &
         FORM='UNFORMATTED', RECL=State_Grid%NX*State_Grid%NY)
    OPEN(687,FILE='radf/radf'//YYYYMM//'CLB.out',ACCESS='DIRECT',    &
         FORM='UNFORMATTED', RECL=State_Grid%NX*State_Grid%NY)
    OPEN(688,FILE='radf/radf'//YYYYMM//'FUL.out',ACCESS='DIRECT',    &
         FORM='UNFORMATTED', RECL=State_Grid%NX*State_Grid%NY)
    OPEN(689,FILE='radf/radf'//YYYYMM//'FDL.out',ACCESS='DIRECT',    &
         FORM='UNFORMATTED', RECL=State_Grid%NX*State_Grid%NY)
    OPEN(690,FILE='radf/radf'//YYYYMM//'FLA.out',ACCESS='DIRECT',    &
         FORM='UNFORMATTED', RECL=State_Grid%NX*State_Grid%NY)

    OPEN(785,FILE='radf/radf'//YYYYMM//'TAOD.out',ACCESS='DIRECT',   &
         FORM='UNFORMATTED', RECL=State_Grid%NX*State_Grid%NY)
    OPEN(786,FILE='radf/radf'//YYYYMM//'TCST.out',ACCESS='DIRECT',   &
         FORM='UNFORMATTED', RECL=State_Grid%NX*State_Grid%NY)
    OPEN(787,FILE='radf/radf'//YYYYMM//'TCSB.out',ACCESS='DIRECT',   &
         FORM='UNFORMATTED', RECL=State_Grid%NX*State_Grid%NY)
    OPEN(788,FILE='radf/radf'//YYYYMM//'TFST.out',ACCESS='DIRECT',   &
         FORM='UNFORMATTED', RECL=State_Grid%NX*State_Grid%NY)
    OPEN(789,FILE='radf/radf'//YYYYMM//'TFSA.out',ACCESS='DIRECT',   &
         FORM='UNFORMATTED', RECL=State_Grid%NX*State_Grid%NY)
    OPEN(790,FILE='radf/radf'//YYYYMM//'TFSB.out',ACCESS='DIRECT',   &
         FORM='UNFORMATTED', RECL=State_Grid%NX*State_Grid%NY)

    OPEN(883,FILE='radf/radf'//YYYYMM//'CLDF.out',ACCESS='DIRECT',   & ! cloud cover
         FORM='UNFORMATTED', RECL=State_Grid%NX*State_Grid%NY)
    OPEN(884,FILE='radf/radf'//YYYYMM//'MALB.out',ACCESS='DIRECT',   &
         FORM='UNFORMATTED', RECL=State_Grid%NX*State_Grid%NY)
    OPEN(885,FILE='radf/radf'//YYYYMM//'ALB.out',ACCESS='DIRECT',    &
         FORM='UNFORMATTED', RECL=State_Grid%NX*State_Grid%NY)
    OPEN(886,FILE='radf/radf'//YYYYMM//'CST.out',ACCESS='DIRECT',    &
         FORM='UNFORMATTED', RECL=State_Grid%NX*State_Grid%NY)
    OPEN(887,FILE='radf/radf'//YYYYMM//'CSB.out',ACCESS='DIRECT',    &
         FORM='UNFORMATTED', RECL=State_Grid%NX*State_Grid%NY)
    OPEN(888,FILE='radf/radf'//YYYYMM//'FST.out',ACCESS='DIRECT',    &
         FORM='UNFORMATTED', RECL=State_Grid%NX*State_Grid%NY)
    OPEN(588,FILE='radf/radf'//YYYYMM//'CLD.out',ACCESS='DIRECT',    &
         FORM='UNFORMATTED', RECL=State_Grid%NX*State_Grid%NY)
    OPEN(589,FILE='radf/radf'//YYYYMM//'CLD0.out',ACCESS='DIRECT',   &
         FORM='UNFORMATTED', RECL=State_Grid%NX*State_Grid%NY)
    OPEN(889,FILE='radf/radf'//YYYYMM//'FSA.out',ACCESS='DIRECT',    &
         FORM='UNFORMATTED', RECL=State_Grid%NX*State_Grid%NY)
    OPEN(890,FILE='radf/radf'//YYYYMM//'AOD.out',ACCESS='DIRECT',    &
         FORM='UNFORMATTED', RECL=State_Grid%NX*State_Grid%NY)
    OPEN(590,FILE='radf/radf'//YYYYMM//'COD.out',ACCESS='DIRECT',    &
         FORM='UNFORMATTED', RECL=State_Grid%NX*State_Grid%NY)
    OPEN(591,FILE='radf/radf'//YYYYMM//'CODGC.out',ACCESS='DIRECT',  &
         FORM='UNFORMATTED', RECL=State_Grid%NX*State_Grid%NY)
    OPEN(891,FILE='radf/radf'//YYYYMM//'AAOD.out',ACCESS='DIRECT',   &
         FORM='UNFORMATTED', RECL=State_Grid%NX*State_Grid%NY)
    OPEN(892,FILE='radf/radf'//YYYYMM//'ABS.out',ACCESS='DIRECT',    &
         FORM='UNFORMATTED', RECL=State_Grid%NX*State_Grid%NY)
    OPEN(893,FILE='radf/radf'//YYYYMM//'FSB.out',ACCESS='DIRECT',    &
         FORM='UNFORMATTED', RECL=State_Grid%NX*State_Grid%NY)
    OPEN(894,FILE='radf/radf'//YYYYMM//'WCL.out',ACCESS='DIRECT',    &
         FORM='UNFORMATTED', RECL=State_Grid%NX*State_Grid%NY)
    OPEN(895,FILE='radf/radf'//YYYYMM//'WCI.out',ACCESS='DIRECT',    &
         FORM='UNFORMATTED', RECL=State_Grid%NX*State_Grid%NY)

    OPEN(896,FILE='radf/radf'//YYYYMM//'AODOUT1.out',                &
         ACCESS='DIRECT',FORM='UNFORMATTED',                         &
         RECL=State_Grid%NX*State_Grid%NY)
    OPEN(897,FILE='radf/radf'//YYYYMM//'AAODOUT1.out',               &
         ACCESS='DIRECT',FORM='UNFORMATTED',                         &
         RECL=State_Grid%NX*State_Grid%NY)
    OPEN(898,FILE='radf/radf'//YYYYMM//'AEextOUT1.out',              &
         ACCESS='DIRECT',FORM='UNFORMATTED',                         &
         RECL=State_Grid%NX*State_Grid%NY)
    OPEN(899,FILE='radf/radf'//YYYYMM//'AEabsOUT1.out',              &
         ACCESS='DIRECT',FORM='UNFORMATTED',                         &
         RECL=State_Grid%NX*State_Grid%NY)

    OPEN(796,FILE='radf/radf'//YYYYMM//'AODOUT3.out',                &
         ACCESS='DIRECT',FORM='UNFORMATTED',                         &
         RECL=State_Grid%NX*State_Grid%NY)
    OPEN(797,FILE='radf/radf'//YYYYMM//'AAODOUT3.out',               &
         ACCESS='DIRECT',FORM='UNFORMATTED',                         &
         RECL=State_Grid%NX*State_Grid%NY)
    OPEN(798,FILE='radf/radf'//YYYYMM//'AEextOUT3.out',              &
         ACCESS='DIRECT',FORM='UNFORMATTED',                         &
         RECL=State_Grid%NX*State_Grid%NY)
    OPEN(799,FILE='radf/radf'//YYYYMM//'AEabsOUT3.out',              &
         ACCESS='DIRECT',FORM='UNFORMATTED',                         &
         RECL=State_Grid%NX*State_Grid%NY)

    OPEN(811,FILE='radf/radf'//YYYYMM//'AOD25.out',ACCESS='DIRECT',  &
         FORM='UNFORMATTED', RECL=State_Grid%NX*State_Grid%NY)
    OPEN(812,FILE='radf/radf'//YYYYMM//'AOD50.out',ACCESS='DIRECT',  &
         FORM='UNFORMATTED', RECL=State_Grid%NX*State_Grid%NY)
    OPEN(813,FILE='radf/radf'//YYYYMM//'AOD75.out',ACCESS='DIRECT',  &
         FORM='UNFORMATTED', RECL=State_Grid%NX*State_Grid%NY)
    OPEN(814,FILE='radf/radf'//YYYYMM//'FAOD25.out',ACCESS='DIRECT', &
         FORM='UNFORMATTED', RECL=State_Grid%NX*State_Grid%NY)
    OPEN(815,FILE='radf/radf'//YYYYMM//'FAOD50.out',ACCESS='DIRECT', &
         FORM='UNFORMATTED', RECL=State_Grid%NX*State_Grid%NY)
    OPEN(816,FILE='radf/radf'//YYYYMM//'FAOD75.out',ACCESS='DIRECT', &
         FORM='UNFORMATTED', RECL=State_Grid%NX*State_Grid%NY)

    OPEN(821,FILE='radf/radf'//YYYYMM//'VIS.out',ACCESS='DIRECT',    &
         FORM='UNFORMATTED', RECL=State_Grid%NX*State_Grid%NY)
    OPEN(822,FILE='radf/radf'//YYYYMM//'HAZE.out',ACCESS='DIRECT',   &
         FORM='UNFORMATTED', RECL=State_Grid%NX*State_Grid%NY)
    OPEN(823,FILE='radf/radf'//YYYYMM//'FOG.out',ACCESS='DIRECT',    &
         FORM='UNFORMATTED', RECL=State_Grid%NX*State_Grid%NY)

    !WRITE(1000,105) KYEAR,KMON,KDAY,KHOUR,KMIN
    WRITE(1001,105) KYEAR,KMON,KDAY,KHOUR,KMIN
    WRITE(1001,105) NPOUTSTEPS
    XCOUNT = float(NPOUTSTEPS)
105 FORMAT(I4,I3,I3,I3,I3,500(1PE10.2))

    !IF(NPOUTSTEPS.EQ.0.and.NPOUTSTEPS0.NE.0) THEN
    WRITE(985,*)KYEAR, KMON-1
    WRITE(987,*)KYEAR, KMON-1
    WRITE(988,*)KYEAR, KMON-1
    WRITE(989,*)KYEAR, KMON-1
    WRITE(992,*)KYEAR, KMON-1
    WRITE(991,*)KYEAR, KMON-1
    WRITE(993,*)KYEAR, KMON-1
    WRITE(994,*)KYEAR, KMON-1
    WRITE(995,*)KYEAR, KMON-1
    DO I = 1, State_Grid%NX
       WRITE(985,981)(ZCST(I,J)/XCOUNT,J=1,State_Grid%NY)
       !WRITE(987,981)(ZCSB(I,J)/XCOUNT,J=1,State_Grid%NY)
       WRITE(987,981)(ZCLD(I,J)/XCOUNT,J=1,State_Grid%NY)
       WRITE(988,981)(ZFST(I,J)/XCOUNT,J=1,State_Grid%NY)
       WRITE(989,981)(ZAOD(I,J)/XCOUNT,J=1,State_Grid%NY)
       !WRITE(992,981)(ZABS(I,J)/XCOUNT,J=1,State_Grid%NY)
       !WRITE(991,981)(ZFUL(I,J)/XCOUNT,J=1,State_Grid%NY)
       !WRITE(993,981)(ZFDL(I,J)/XCOUNT,J=1,State_Grid%NY)
       !WRITE(994,981)(ZCLT(I,J)/XCOUNT,J=1,State_Grid%NY)
       !WRITE(995,981)(ZCLB(I,J)/XCOUNT,J=1,State_Grid%NY)

       DO J=1,State_Grid%NY
          FAOD25(I,J)= float(NAOD25(I,J))/XCOUNT
          FAOD50(I,J)= float(NAOD50(I,J))/XCOUNT
          FAOD75(I,J)= float(NAOD75(I,J))/XCOUNT

          ZAOD25(I,J) = ZAOD25(I,J)/(1.d-20+float(NAOD25(I,J)))
          ZAOD50(I,J) = ZAOD50(I,J)/(1.d-20+float(NAOD50(I,J)))
          ZAOD75(I,J) = ZAOD75(I,J)/(1.d-20+float(NAOD75(I,J)))
       ENDDO
    ENDDO

    ZCLDF(:,:)=ZCLDF(:,:)/XCOUNT
    ZALB(:,:)=ZALB(:,:)/XCOUNT

    ZCST(:,:)=ZCST(:,:)/XCOUNT
    ZCSB(:,:)=ZCSB(:,:)/XCOUNT
    ZFST(:,:)=ZFST(:,:)/XCOUNT
    ZCLD(:,:)=ZCLD(:,:)/XCOUNT
    ZCLD0(:,:)=ZCLD0(:,:)/XCOUNT
    ZFSB(:,:)=ZFSB(:,:)/XCOUNT
    ZFSA(:,:)=ZFSA(:,:)/XCOUNT
    ZAOD(:,:)=ZAOD(:,:)/XCOUNT
    ZCOD(:,:)=ZCOD(:,:)/XCOUNT
    ZCODGC(:,:)=ZCODGC(:,:)/XCOUNT
    ZAODOUT1(:,:)=ZAODOUT1(:,:)/XCOUNT
    ZAODOUT3(:,:)=ZAODOUT3(:,:)/XCOUNT

    ZVIS(:,:)=ZVIS(:,:)/XCOUNT
    THAZ(:,:)=THAZ(:,:)/XCOUNT  ! fraction of time with haze
    TFOG(:,:)=TFOG(:,:)/XCOUNT  ! fraction of time with fog

    ZAAOD(:,:)=ZAAOD(:,:)/XCOUNT
    ZAAODOUT1(:,:)=ZAAODOUT1(:,:)/XCOUNT
    ZAAODOUT3(:,:)=ZAAODOUT3(:,:)/XCOUNT
    ZABS(:,:)=ZABS(:,:)/XCOUNT

    !ZCLT(:,:)=ZCLT(:,:)/XCOUNT
    !ZCLB(:,:)=ZCLB(:,:)/XCOUNT
    !ZFUL(:,:)=ZFUL(:,:)/XCOUNT
    !ZFDL(:,:)=ZFDL(:,:)/XCOUNT
    !ZFLA(:,:)=ZFLA(:,:)/XCOUNT

    ZWCL(:,:)=ZWCL(:,:)/XCOUNT
    ZWCI(:,:)=ZWCI(:,:)/XCOUNT

    ZMALB(:,:,:)=ZMALB(:,:,:)/XCOUNT

    DO I = 1, State_Grid%NX
       WRITE(991,981)(ZAOD25(I,J),J=1,State_Grid%NY)
       WRITE(992,981)(ZAOD75(I,J),J=1,State_Grid%NY)
       WRITE(993,981)(ZVIS(I,J),J=1,State_Grid%NY)
       WRITE(994,981)(THAZ(I,J),J=1,State_Grid%NY)
       WRITE(995,981)(TFOG(I,J),J=1,State_Grid%NY)

       DO J=1,State_Grid%NY
          ZAEEXT(I,J)=-LOG10(ZAOD(I,J)/ZAODOUT1(I,J)) &
                      /LOG10(YWLS(KOUT2)/YWLS(KOUT1))            !WL(KOUT1)_WL(KOUT2) Angstrom Exponent
          ZAEABS(I,J)=-LOG10(ZAAOD(I,J)/ZAAODOUT1(I,J)) &
                      /LOG10(YWLS(KOUT2)/YWLS(KOUT1))            !WL(KOUT1)_WL(KOUT2) Angstrom Exponent
          ZAEEXTOUT3(I,J)=-LOG10(ZAOD(I,J)/ZAODOUT3(I,J)) &
                          /LOG10(YWLS(KOUT2)/YWLS(KOUT3))            !WL(KOUT3)_WL(KOUT2) Angstrom Exponent
          ZAEABSOUT3(I,J)=-LOG10(ZAAOD(I,J)/ZAAODOUT3(I,J)) &
                          /LOG10(YWLS(KOUT2)/YWLS(KOUT3))            !WL(KOUT3)_WL(KOUT2) Angstrom Exponent
       ENDDO
       WRITE(996,981)(ZAODOUT1(I,J),J=1,State_Grid%NY)
       WRITE(997,981)(ZAAODOUT1(I,J),J=1,State_Grid%NY)
       WRITE(998,981)(ZAEEXT(I,J),J=1,State_Grid%NY)
       WRITE(999,981)(ZAEABS(I,J),J=1,State_Grid%NY)
    ENDDO

    !WRITE(686,REC=1)ZCLT
    !WRITE(687,REC=1)ZCLB
    !WRITE(688,REC=1)ZFUL
    !WRITE(689,REC=1)ZFDL
    !WRITE(690,REC=1)ZFLA

    WRITE(811,REC=1)ZAOD25
    WRITE(812,REC=1)ZAOD50
    WRITE(813,REC=1)ZAOD75

    WRITE(814,REC=1)FAOD25
    WRITE(815,REC=1)FAOD50
    WRITE(816,REC=1)FAOD75

    WRITE(821,REC=1)ZVIS
    WRITE(822,REC=1)THAZ
    WRITE(823,REC=1)TFOG

    WRITE(883,REC=1)ZCLDF
    WRITE(885,REC=1)ZALB

    WRITE(886,REC=1)ZCST
    WRITE(887,REC=1)ZCSB
    WRITE(888,REC=1)ZFST
    WRITE(588,REC=1)ZCLD
    WRITE(589,REC=1)ZCLD0
    WRITE(889,REC=1)ZFSA
    WRITE(890,REC=1)ZAOD
    WRITE(590,REC=1)ZCOD
    WRITE(591,REC=1)ZCODGC
    WRITE(891,REC=1)ZAAOD
    WRITE(892,REC=1)ZABS
    WRITE(893,REC=1)ZFSB

    WRITE(894,REC=1)ZWCL
    WRITE(895,REC=1)ZWCI

    WRITE(896,REC=1)ZAODOUT1
    WRITE(897,REC=1)ZAAODOUT1
    WRITE(898,REC=1)ZAEEXT
    WRITE(899,REC=1)ZAEABS

    WRITE(796,REC=1)ZAODOUT3
    WRITE(797,REC=1)ZAAODOUT3
    WRITE(798,REC=1)ZAEEXTOUT3
    WRITE(799,REC=1)ZAEABSOUT3

    RECID=0
    DO IWL=1,NBS
       RECID=RECID+1
       WRITE(884,REC=RECID) (ZMALB(:,:,IWL))
    ENDDO

    WRITE(1081,*)KYEAR, KMON-1
    WRITE(1082,*)KYEAR, KMON-1
    WRITE(1083,*)KYEAR, KMON-1
    WRITE(1084,*)KYEAR, KMON-1
    WRITE(1085,*)KYEAR, KMON-1

    DO I = 1, State_Grid%NX
       WRITE(1081,981)(ZTCST(I,J,1)/XCOUNT,J=1,State_Grid%NY)
       WRITE(1082,981)(ZTCST(I,J,2)/XCOUNT,J=1,State_Grid%NY)
       WRITE(1083,981)(ZTCST(I,J,3)/XCOUNT,J=1,State_Grid%NY)
       WRITE(1084,981)(ZTCST(I,J,4)/XCOUNT,J=1,State_Grid%NY)
       WRITE(1085,981)(ZTCST(I,J,5)/XCOUNT,J=1,State_Grid%NY)
    ENDDO

    ZTAOD(:,:,:)=ZTAOD(:,:,:)/XCOUNT

    ZTCST(:,:,:)=ZTCST(:,:,:)/XCOUNT
    ZTCSB(:,:,:)=ZTCSB(:,:,:)/XCOUNT
    ZTFST(:,:,:)=ZTFST(:,:,:)/XCOUNT
    ZTFSB(:,:,:)=ZTFSB(:,:,:)/XCOUNT
    ZTFSA(:,:,:)=ZTFSA(:,:,:)/XCOUNT

    RECID=0
    DO ITYP=1,NTYP
       RECID=RECID+1
       WRITE(785,REC=RECID) (ZTAOD(:,:,ITYP))

       WRITE(786,REC=RECID) (ZTCST(:,:,ITYP))
       WRITE(787,REC=RECID) (ZTCSB(:,:,ITYP))
       WRITE(788,REC=RECID) (ZTFST(:,:,ITYP))
       WRITE(789,REC=RECID) (ZTFSA(:,:,ITYP))
       WRITE(790,REC=RECID) (ZTFSB(:,:,ITYP))
    ENDDO

    ! AREA weighted average
    ZFSTAVE  = 0.
    ZCSTAVE  = 0.
    ZCLDAVE  = 0.
    ZCLD0AVE = 0.
    ZFSAAVE  = 0.
    !ZFULAVE = 0.     ! longwave

    ZTFSTAVE = 0.
    TAREA    = 0.
    DO J = 1, State_Grid%NY
       YAREA = State_Grid%AREA_M2(1,J)
       DO I = 1, State_Grid%NX
          ZFSTAVE = ZFSTAVE + ZFST(I,J)*YAREA
          ZCSTAVE = ZCSTAVE + ZCST(I,J)*YAREA
          ZCLDAVE = ZCLDAVE + ZCLD(I,J)*YAREA
          ZCLD0AVE = ZCLD0AVE + ZCLD0(I,J)*YAREA

          ZFSAAVE = ZFSAAVE + ZFSA(I,J)*YAREA
          !ZFULAVE = ZFULAVE + ZFUL(I,J)*YAREA
          TAREA = TAREA + YAREA
          DO ITYP=1,NTYP
             ZTFSTAVE(ITYP)=ZTFSTAVE(ITYP)+ZTFST(I,J,ITYP)*YAREA
          ENDDO
       ENDDO
    ENDDO
    ZFSTAVE = ZFSTAVE/TAREA
    ZCSTAVE = ZCSTAVE/TAREA
    ZCLDAVE = ZCLDAVE/TAREA
    ZCLD0AVE = ZCLD0AVE/TAREA

    ZFSAAVE = ZFSAAVE/TAREA
    !ZFULAVE = ZFULAVE/TAREA
    DO ITYP=1,NTYP
       ZTFSTAVE(ITYP)=ZTFSTAVE(ITYP)/TAREA
    ENDDO

    WRITE(990,106)KYEAR, KMON, KDAY, ZCLDAVE,ZCLD0AVE, &
          ZFSTAVE,ZCSTAVE, ZFSAAVE, &
          (ZTFSTAVE(ITYP),ITYP=1,NTYP)
106 FORMAT(I4,I3,I3,2(1PE11.3),20(1PE10.2))
    flush(990)
    !flush(1002)

    ZCLDF  = 0.
    ZALB   = 0.

    ZCST   = 0.
    ZCSB   = 0.
    ZFST   = 0.
    ZCLD   = 0.
    ZCLD0  = 0.
    ZFSB   = 0.
    ZFSA   = 0.
    ZAOD   = 0.
    ZCOD   = 0.
    ZCODGC = 0.
    ZAOD25 = 0.
    ZAOD50 = 0.
    ZAOD75 = 0.

    NAOD25 = 0
    NAOD50 = 0
    NAOD75 = 0

    ZVIS   = 0.
    THAZ   = 0.
    TFOG   = 0.

    ZAAOD  = 0.
    ZABS   = 0.
    !ZCLT  = 0.
    !ZCLB  = 0.
    !ZFUL  = 0.
    !ZFDL  = 0.
    !ZFLA  = 0.

    ZWCL   = 0.
    ZWCI   = 0.

    ZTAOD  = 0.

    ZTCST  = 0.
    ZTCSB  = 0.
    ZTFST  = 0.
    ZTFSB  = 0.
    ZTFSA  = 0.

981 FORMAT(200(1PE9.1))

    FLUSH(890)
    FLUSH(891)
    FLUSH(1001)

  END SUBROUTINE APM_RADFOUT
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: aeronum
!
! !DESCRIPTION: Subroutine AERONUM calculates aerosol number concentration
!  based on the mass.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE AERONUM( Input_Opt,  State_Chm, State_Diag, &
                      State_Grid, State_Met, RC )
!
! !USES:
!
    USE APM_INIT_MOD,   ONLY : APMIDS
    USE APM_INIT_MOD,   ONLY : NGCOND
    USE APM_INIT_MOD,   ONLY : NSO4
    USE APM_INIT_MOD,   ONLY : VDRY
    USE APM_INIT_MOD,   ONLY : DENSULF
    USE ErrCode_Mod
    USE Input_Opt_Mod,  ONLY : OptInput
    USE State_Diag_Mod, ONLY : DgnState
    USE State_Grid_Mod, ONLY : GrdState
    USE State_Chm_Mod,  ONLY : ChmState
    USE State_Met_Mod,  ONLY : MetState
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
    TYPE(GrdState), INTENT(IN)    :: State_Grid  ! Grid State object
    TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
    TYPE(DgnState), INTENT(INOUT) :: State_Diag  ! Diagnostics State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure
!
! !REVISION HISTORY:
!  28 Aug 2008 - F. Yu - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: I,J,L,N,SIZENUM
    REAL*8  :: XMATEMP

    ! Make a pointer to the tracer array
    REAL*8, POINTER :: Spc(:,:,:,:)

    ! Initialize
    RC = GC_SUCCESS

    ! Point to Spc
    Spc => State_Chm%Species

    XN4D = 0.D0

    !$OMP PARALLEL DO                             &
    !$OMP DEFAULT( SHARED                       ) &
    !$OMP SCHEDULE( DYNAMIC                     ) &
    !$OMP PRIVATE( I, J, L, N, SIZENUM, XMATEMP )
    DO L = 1, State_Grid%NZ
    DO J = 1, State_Grid%NY
    DO I = 1, State_Grid%NX

       ! Should we put the N loop on the outside?
       DO N = 1, NSO4
          SIZENUM =APMIDS%id_SO4BIN1+N-1 !Luodebug
          XMATEMP=Spc(I,J,L,SIZENUM) / State_Met%AIRVOL(I,J,L)
          XMATEMP= MAX( 1.d-40, XMATEMP)
          XN4D(I,J,L,N)=XMATEMP/(DENSULF*VDRY(N))*1.E-9 !XN4D in #/cm3, VDRY in m3
       ENDDO

    ENDDO
    ENDDO
    ENDDO
    !$OMP END PARALLEL DO

    ! Clear the pointer
    Spc => NULL()

  END SUBROUTINE AERONUM
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: apm3dqgeos
!
! !DESCRIPTION: Subroutine APM3DQGEOS ...
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE APM3DQGEOS( Input_Opt, State_Grid, State_Met, IY )
!
! !USES:
!
    USE PRESSURE_MOD,   ONLY : GET_PCENTER
    USE PRESSURE_MOD,   ONLY : GET_PEDGE
    USE APM_NUCL_MOD,   ONLY : IONRATE0
    USE Input_Opt_Mod,  ONLY : OptInput
    USE State_Grid_Mod, ONLY : GrdState
    USE State_Met_Mod,  ONLY : MetState
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN) :: Input_Opt   ! Input Options object
    TYPE(GrdState), INTENT(IN) :: State_Grid  ! Grid State object
    TYPE(MetState), INTENT(IN) :: State_Met   ! Meteorology State object
!
! !REVISION HISTORY:
!  17 Mar 2010 - F. Yu       - Initial version
!  08 Nov 2010 - R. Yantosca - Added ProTeX headers
!  01 Mar 2012 - R. Yantosca - Now use GET_AREA_CM2(I,J,L) from grid_mod.F90
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: IY
    INTEGER :: ISURF, I, J, L
    REAL*8  :: XLON,XLAT,YPR,YQ,YPSURF

    DO J = 1, State_Grid%NY
    DO I = 1, State_Grid%NX

       YPSURF = GET_PEDGE(I,J,1) !surface P in mb
       IF( State_Met%IsWater(I,J) .or. State_Met%IsIce(I,J) ) THEN
          ISURF = 0
       ELSE
          ISURF = 1
       ENDIF
       DO L = 1, State_Grid%NZ
          XLON = State_Grid%XMid( I, J ) ! Grid box longitude [degrees]
          XLAT = State_Grid%YMid( I, J ) ! Grid box latitude center [degree]
          YPR =  GET_PCENTER(I,J,L)      ! in mb
          CALL IONRATE0(Input_Opt%CHEM_INPUTS_DIR, &
                        IY,ISURF,YPSURF,XLON,XLAT,YPR,YQ)
          XQ3D(I,J,L)= YQ
       ENDDO
    ENDDO
    ENDDO

  END SUBROUTINE APM3DQGEOS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: apm3dq
!
! !DESCRIPTION: Function APM3DQ
!\\
!\\
! !INTERFACE:
!
  subroutine APM3DQ( IY, Input_Opt, State_Grid, State_Met )
!
! !USES:
!
    USE PRESSURE_MOD,   ONLY : GET_PCENTER
    USE TIME_MOD,       ONLY : GET_YEAR,   GET_MONTH
    USE APM_NUCL_MOD,   ONLY : IONRATE
    USE State_Met_Mod,  ONLY : MetState
    USE State_Grid_Mod, ONLY : GrdState
    USE Input_Opt_Mod,  ONLY : OptInput
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN) :: Input_Opt   ! Input Options object
    TYPE(GrdState), INTENT(IN) :: State_Grid  ! Grid State object
    TYPE(MetState), INTENT(IN) :: State_Met   ! Meteorology State object
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER   :: ISURF, I, J, L, IY, IYEAR, IMONTH
    REAL*8    :: XLON,XLAT,YPR,YQ

    IMONTH = GET_MONTH()
    IYEAR = GET_YEAR()   ! IY = 0, current modeling year
    IF(IY.EQ.1) IYEAR = 1997 ! solar min, max Q
    IF(IY.EQ.-1) IYEAR = 1990 ! solar max, min Q

    write(6,*)"GCR IONIZATION: YEAR =", IYEAR, "Month = ",IMONTH

    DO J = 1, State_Grid%NY
       XLAT = State_Grid%YMID(1,J) ! Grid box latitude center [degree]
       DO I = 1, State_Grid%NX
          XLON = State_Grid%XMID(I,J) !Grid box longitude [degrees]
          IF( State_Met%IsWater(I,J) .or. State_Met%IsIce(I,J) ) THEN
             ISURF = 0
          ELSE
             ISURF = 1
          ENDIF
          DO L = 1, State_Grid%NZ
             YPR =  GET_PCENTER(I,J,L) ! in mb
             CALL IONRATE(Input_Opt%CHEM_INPUTS_DIR, &
                          ISURF,XLON,XLAT,YPR,IYEAR,IMONTH,YQ)
             XQ3D(I,J,L)= YQ
          ENDDO
       ENDDO
    ENDDO

  end subroutine APM3DQ
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ysitesij_gl
!
! !DESCRIPTION: Function YSITESIJ_GL
!\\
!\\
! !INTERFACE:
!
  subroutine YSITESIJ_GL( State_Grid )
!
! !USES:
!
    USE APM_INIT_MOD,   ONLY : MAXSITE,MSITE,ISITES,JSITES,IFSITEOUT
    USE State_Grid_Mod, ONLY : GrdState
!
! !INPUT PARAMETERS:
!
    TYPE(GrdState), INTENT(IN) :: State_Grid
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER, PARAMETER :: MS   = 47
    REAL*8             :: ZLAT(MS),ZLON(MS), YLAT, XLON
    INTEGER            :: NSITE, ISITE, JSITE

    !  CN10LT.dat
    !  71.32   -156.6  231     !site1   Point Barrow
    !  67.97   24.12   802     !site2   Pallas
    !  67.77   29.58   823     !site3   Varri
    !  61.85   24.28   2016    !site4   Hyytiala
    !  59.78   21.38   2921    !site5   Uto
    !  58.77   17.4    2567    !site6   Aspvreten
    !  53.33   -9.9    1907    !site7   Mace Head
    !  51.5    12.9    4664    !site8   Melpitz
    !  43.9    -60     869     !site9   Stable Island
    !  43.11   -70.95  7039    !site10  Thomson Farm
    !  41.05   -124.15 918     !site11  Trinidad Head
    !  40.1    -88.3   5038    !site12  Bondville
    !  36.8    -97.5   5064    !site13  SGP
    !  -14.24  -170.57 270     !site14  American Samoa
    !  -25.54  25.75   2340    !site15  Botsalano, South Africa
    !  -35.66  148.15  1800    !site16  Bago State forest, Australia
    !  -40.682 144.688 1203.00 !site17  Cape Grim
    !  -70.65  -8.25   324     !site18  Neumayer
    !  -89.98  -24.8   156     !site19  South Pole
    !  29.44   79.62   2730    !site20  Mukteshwar, India
    !  36.28   100.9   2030    !site21  Mount Waliguan, China
    !
    !  CN10MT.dat
    !  47.6    8       803     !site22  Jungfraujoch
    !  27.96   86.82   900     !site23  Pyramid
    !  19.54   -155.58 363     !site24  Mauna Loa
    !
    !  Asian cities
    !  Urumqi        China  43.83   87.62  !site25
    !  Shenyang      China  41.81  123.43  !site26
    !  Beijing       China  39.91  116.41  !site27
    !  Taiyuan       China  37.87  112.55  !site28
    !  Seoul         Korea  37.57  126.98  !site29
    !  Anmyeon-do    Korea  36.37  126.32  !site30
    !  Lanzhou       China  36.06  103.83  !site31
    !  Tokoyo        Japan  35.69  139.69  !site32
    !  Shirahama     Japan  33.66  135.38  !site33
    !  Shanghai      China  31.23  121.47  !site34
    !  Chengdu       China  30.66  104.07  !site35
    !  Wuhan         China  30.59  114.31  !site36
    !  Lhasa         China  29.65   91.14  !site37
    !  New Delhi     India  28.64   77.22  !site38
    !  Taipei        China  25.09  121.56  !site39
    !  Kunming       China  25.04  102.72  !site40
    !  Guangzhou     China  23.13  113.26  !site41
    !
    !  AERONET -- Amazon
    !  Alta_Floresta        -9.87  -56.10  !site42
    !  Ji_Parana_SE        -10.93  -61.85  !site43
    !  CUIABA-MIRANDA      -15.73  -56.02  !site44
    !  Campo_Grande_SONDA  -20.43  -54.53  !site45
    !  Sao_Paulo           -23.56  -46.73  !site46
    !  !Cordoba-CETT       -31.52  -64.45  !site47  --> 3., 93.0 for IRF output

    DATA ZLAT/71.32, 67.97,67.77,61.85,  59.78, 58.77, 53.33,51.5,          &
              43.9, 43.11, 41.05, 40.1, 36.8, -14.24,-25.54,-35.66,         &
              -40.68,-70.65,-89.98,29.44,36.28,47.6,27.96,19.54,            &
              43.83,41.81,39.91,37.87,37.57,36.37,36.06,35.69,              &
              33.66,31.23,30.66,30.59,29.65,28.64,25.09,25.04,23.13,        &
              -9.87,-10.93,-15.73,-20.43,-23.56,3.0/

    DATA ZLON/-156.6,24.12,29.58,24.28,21.38,17.40, -9.9,12.9,              &
              -60.0,-70.95,-124.15,-88.3,-97.5,-170.57, 25.75,148.15,       &
              144.69, -8.25, -24.80,79.62,100.9,8.0,86.82,-155.58,          &
              87.62,123.43,116.41,112.55,126.98,126.32,103.83,139.69,       &
              135.38,121.47,104.07,114.31,91.14,77.22,121.56,102.72,113.26, &
              -56.10,-61.85,-56.02,-54.53,-46.73,93.0/

    MSITE = MS
    IF(MSITE.GT.MAXSITE) THEN
       WRITE(6,*)"MSITE>MAXSITE; NEED TO CHECK",MSITE,MAXSITE
       STOP
    ENDIF
    ISITES = 0
    JSITES = 0
    DO NSITE = 1, MSITE
       YLAT = ZLAT(NSITE)
       XLON = ZLON(NSITE)
       CALL SITEIJ(NSITE,XLON, YLAT,ISITE,JSITE,State_Grid)
       ISITES(NSITE) = ISITE
       JSITES(NSITE) = JSITE
    ENDDO

    ! output at selected sites
    IFSITEOUT = 0
    IFOUTIJ = 0
    SITEID  = 0
    DO NSITE=1,MSITE
       IF(NSITE.EQ.4.or.NSITE.EQ.13.or.NSITE.EQ.27) THEN !Hyytiala,SGP,Beijing
          !IF(NSITE.GE.45) THEN !No sites
          !IF(NSITE.EQ.1.or.NSITE.EQ.2.or.NSITE.EQ.4.or.    &
          !   NSITE.EQ.7.or.NSITE.EQ.8.or.NSITE.EQ.13.or.   &
          !   NSITE.EQ.15.or.NSITE.EQ.16.or.NSITE.EQ.18.or. &
          !   NSITE.EQ.20.or.NSITE.EQ.21.or.NSITE.EQ.27.or. &
          !   NSITE.EQ.34.or.NSITE.EQ.35.or.NSITE.EQ.36.or. &
          !   NSITE.GE.41) THEN
          !IF(NSITE.GE.2.and.NSITE.LE.8) THEN
          IFSITEOUT(NSITE) = 1
          ISITE = ISITES(NSITE)
          JSITE = JSITES(NSITE)
          IFOUTIJ(ISITE,JSITE)=1
          SITEID(ISITE,JSITE)=NSITE
       ENDIF
    ENDDO

    WRITE(6,*)"SITE LOCATION AND ID"
    DO NSITE = 1, MSITE
       WRITE(6,25)NSITE,ZLAT(NSITE),ZLON(NSITE), &
            ISITES(NSITE),JSITES(NSITE),IFSITEOUT(NSITE)
    ENDDO
25  FORMAT(I3,F8.2,F8.2,I4,I4,I2)

  end subroutine YSITESIJ_GL
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ysitesij_gl1
!
! !DESCRIPTION: Function YSITESIJ_GL1
!\\
!\\
! !INTERFACE:
!
  subroutine  YSITESIJ_GL1( State_Grid )
!
! !USES:
!
    USE APM_INIT_MOD,   ONLY : MAXSITE,MSITE,ISITES,JSITES,IFSITEOUT
    USE APM_INIT_MOD,   ONLY : IFSITE
    USE State_Grid_Mod, ONLY : GrdState
!
! !INPUT PARAMETERS:
!
    TYPE(GrdState), INTENT(IN) :: State_Grid
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER, PARAMETER :: MS   = 58
    REAL*8             :: ZLAT(MS),ZLON(MS), YLAT, XLON
    INTEGER            :: NSITE, ISITE, JSITE

    ! For size distributions
    ! site1 Point Barrow            (71.32, -156.6 )  231 m asl
    ! site2 Neumayer                (-70.65,  -8.25)
    ! site3 Troll                   (-72.02,   2.53)
    ! site4 Hyytiala                (61.85,   24.28)
    !
    ! EPA sopersites
    !site5 New York                 (40.7144,-74.01)
    !site6 Atlanta, Georgia         (33.75,  -84.39)
    !site7 Baltimore, Maryland      (39.29,  -76.61)
    !site8 Pittsburgh, Pennsylvania (40.44,  -80.0 )
    !!site9 Houston, Texas          (29.76,  -95.38)
    !
    ! Other sites
    !site9  APP, NC                 (36.213, -81.692)
    !site10 EGB, Canada             (44.23,  -79.783)
    !site11 DukeForest              (36.0,   -79.08 )
    !site12 Kent, Ohio              (41.15,  -81.36 )
    !site13 SGP                     (36.8    -97.5)
    !site14 Pinnacle State Park     (42.10,  -77.21)
    !site15 Whiteface Mtn           (44.4,   -73.9)
    !Site16 SPL                     (40.46 , -106.744 ) 3210 m asl
    !!site17  Boulder               (40.02,  -105.27)
    !site17  MMSF,Indiana           (39.32,  -86.42)
    !site18  Rochester              (43.16,  -77.61)
    !
    !site19 Manacapara, Brazil      (-3.21,  -60.60)
    !site20 CUIABA-MIRANDA          -15.73   -56.02
    !site21 Sao_Paulo               -23.56   -46.73
    !site22 Melpitz                 (51.5     12.9 )
    !site23 MexicoCity              (19.43,  -99.13)
    !site24 Shanghai China           31.23   121.47
    !site25 Wuhan China              30.59   114.31
    !site26 Chongqing China          29.57   106.52
    !site27 Zhengzhou China          34.77   113.62
    !site28 Beijing China            39.91   116.41
    !site28 Wuqing                  (39.39,  117.01) (same grid box as Beijing)
    !site29 UMBS                    (45.56N,  84.72W)
    !site30 NAAMES1                 (35N,     40W)
    !!site30 JINAN                  (36.67N, 116.95E)
    !site31 NAAMES2                 (46N,     40W)
    !site32 NAAMES3                 (57N,     40W)
    !!site32 Huangshan              (30.17N, 118.15E)
    !site33 LLN                     (23.47,  120.87 )
    !!site34 Chacaltaya            (-16.35,  -68.13 )
    !site34 Chacaltaya             (-16.20,  -68.10 )
    !site35 Puy de Dome             (45.767N,  2.966E,1465m)
    !
    !site36 Trinidad Head           (41.05, -124.15)
    !site37 Sacramento, CA          (38.645,-121.34)
    !site38 Los Angeles, California (34.05, -118.24)
    !site39 Bondville               (40.1,   -88.3 )
    !site40 Ozarks forest site      (38.74,  -92.2 )
    !site41 Brent,AL                (32.94,  -87.16)
    !
    !site42 CPR (Cape San Juan), PR (18.48,  -66.13)
    !site43 Aboa                   (-73.05,  -13.41)
    !Site44 Port-aux-Français      (-49.35,   70.219),French Southern and Antarctic Lands
    !Site45 Gual Pahari, India      (28.43 ,  77.15)
    !SIte46 Nepal Cli Obs           (27.96,   86.81)
    !site47 GuangzhouChina          (23.13,  113.26)
    !Site48 Taipei                  (25.09,  121.56)
    !Site49 Anmyeon-do, South Korea (36.54,  126.33)
    !Site50 DEM_Athens              (37.99,   23.82)
    !Site51 Izana                   (28.31,  -16.50)
    !Site52 Zeppelin mountain       (78.91,   11.89)
    !
    !Site53 Ascension Island        (-7.97,  -14.35)
    !Site54 GRW Graciosa Island, Azores, Portugal (39.09, -28.025)
    !Site55 SMO American Samoa     (-14.24, -170.57)
    !Site56 CGO Cape Grim          (-40.682, 144.688)
    !Site57 Chuuk, FSM               (7.44,  151.86)
    !site58  Albany, NY             (42.65,  -73.76)

    ! CPT?
    ! MLO
    ! ALT

    DATA ZLAT/71.32,-70.65,-72.02,61.85,40.71,33.75,39.29,40.44,            &
              !29.76,34.05,36.0,41.15,36.8,42.10,44.4,40.46,                &
              36.213,44.23,36.0,41.15,36.8,42.10,44.4,40.46,                &
              39.32,43.16,-3.21,-15.73,-23.56,51.5,19.43,31.23,             &
              !30.59,29.57,34.77,39.91,23.13/                               &
              30.59,29.57,34.77,39.91,45.56,35.,46.,57.,                    &
              !30.59,29.57,34.77,39.91,45.56,36.67,46.,30.17,               &
              23.47,-16.20,45.767,41.05,38.645,34.05,40.1,                  &
              38.74,32.94,18.48,-73.05,-49.35,28.43,27.96,                  &
              23.13,25.09,36.54,37.99,28.31,78.91,                          &
              -7.97,39.09,-14.24,-40.682,7.44,42.65/

    DATA ZLON/-156.6,-8.25,2.53,24.28,-74.01,-84.39,-76.61,-80.0,           &
              !-95.38,-118.24,-79.08,-81.36,-97.5,-77.21,-73.9,-106.744,    &
              -81.692,-79.783,-79.08,-81.36,-97.5,-77.21,-73.9,-106.744,    &
              -86.42,-77.61, -60.6,-56.02,-46.73,12.9,-99.13,121.47,        &
              !114.31,106.52,113.62,116.41,113.26/                          &
              114.31,106.52,113.62,116.41,-84.72,-40.,-40.,-40.,            &
              !114.31,106.52,113.62,116.41,-84.72,116.95,-40.,118.15,       &
              120.87,-68.10,2.966,-124.15,-121.34,-118.24,-88.3,-92.2,      &
              -87.16,-66.13,-13.41,70.219,77.15,86.81,113.26,121.56,        &
              126.33,23.82,-16.50,11.89,                                    &
              -14.35,-28.025,-170.57,144.688,151.86,-73.76/

    MSITE = MS
    IF(MSITE.GT.MAXSITE) THEN
       WRITE(6,*)"MSITE>MAXSITE; NEED TO CHECK",MSITE,MAXSITE
       STOP
    ENDIF
    ISITES = 0
    JSITES = 0
    DO NSITE = 1, MSITE
       YLAT = ZLAT(NSITE)
       XLON = ZLON(NSITE)
       CALL SITEIJ(NSITE,XLON, YLAT,ISITE,JSITE,State_Grid)
       ISITES(NSITE) = ISITE
       JSITES(NSITE) = JSITE
    ENDDO

    ! output at selected sites
    IFSITEOUT = 0
    IFOUTIJ   = 0
    SITEID    = 0
    DO NSITE=1,MSITE
       IF(IFSITE.EQ.73) THEN ! only output 3 sites
          IF(NSITE.EQ.4.or.NSITE.EQ.13.or.NSITE.EQ.27) THEN !Hyytiala,SGP, BJ
             IFSITEOUT(NSITE) = 1
             ISITE = ISITES(NSITE)
             JSITE = JSITES(NSITE)
             IFOUTIJ(ISITE,JSITE)=1
             SITEID(ISITE,JSITE)=NSITE
          ENDIF
       !ELSE  ! all sites
       ELSEIF((NSITE.GE.1.and.NSITE.LE.5).or.NSITE.EQ.7.or.NSITE.EQ.9 &
          .or.(NSITE.GE.12.and.NSITE.LE.22).or.NSITE.EQ.24            &
          .or. NSITE.EQ.28.or.(NSITE.GE.30.and.NSITE.LE.38)           &
          .or.(NSITE.GE.41.and.NSITE.LE.58)) THEN
          IFSITEOUT(NSITE) = 1
          ISITE = ISITES(NSITE)
          JSITE = JSITES(NSITE)
          IFOUTIJ(ISITE,JSITE)=1
          SITEID(ISITE,JSITE)=NSITE
       ENDIF
    ENDDO
    WRITE(6,*)"SITE LOCATION AND ID"
    DO NSITE = 1, MSITE
       WRITE(6,25)NSITE,ZLAT(NSITE),ZLON(NSITE), &
                  ISITES(NSITE),JSITES(NSITE),IFSITEOUT(NSITE)
    ENDDO
25  FORMAT(I3,F8.2,F8.2,I4,I4,I2)

  end subroutine YSITESIJ_GL1
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ysitesij_ac
!
! !DESCRIPTION: Function YSITESIJ_AC
!\\
!\\
! !INTERFACE:
!
  subroutine  YSITESIJ_AC( State_Grid )
!
! !USES:
!
    USE APM_INIT_MOD,   ONLY : MAXSITE,MSITE,ISITES,JSITES,IFSITEOUT
    USE APM_INIT_MOD,   ONLY : IFSITE
    USE State_Grid_Mod, ONLY : GrdState
!
! !INPUT PARAMETERS:
!
    TYPE(GrdState), INTENT(IN) :: State_Grid
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER, PARAMETER :: MS   = 47
    REAL*8             :: ZLAT(MS),ZLON(MS), YLAT, XLON
    INTEGER            :: NSITE, ISITE, JSITE

    ! From AEROCOM_diagnostics_091024.xls
    !STATIONID	Station Name	Lat	Lon	Height asl	Sampling Height	Mean Pressure	Country
    !1	Alert	82.450	-62.517	210.0	10.0	988.3	US
    !2	Aspvreten	58.800	17.380	20.0	-999.9	1010.9	SE
    !3	Auchencorth	55.800	-3.283	255.0	10.0	-999.9	UK
    !4	Barrow	71.323	-156.612	11.0	10.0	1012.3	US
    !5	Birkenes	58.383	8.250	190.0	5.0	990.7	NO
    !6	Bondville	40.050	-88.367	213.0	10.0	987.9	US
    !7	Cabauw	51.971	4.927	60.0	-999.9	1013.0	NL
    !8	Cape Grim	-40.682	144.688	94.0	10.0	1002.0	AS
    !9	Cape Point	-34.353	18.490	230.0	3.0	988.4	SA
    !10	Elandsfontein	-25.533	25.750	1424.0	10.0	-999.9	SA
    !11	Finokalia	35.333	25.670	250.0	2.0	990.7	GR
    !	Gual Pahari						IN
    !12	Harwell	51.608	-1.288	60.0	-999.9	1003.7	UK
    !13	Hohenpeissenberg	47.800	11.017	985.0	10.0	900.4	DE
    !14	Hyytiä	61.850	24.283	181.0	2.0	991.7	FI
    !15	Ispra	45.803	8.627	209.0	3.5	988.4	IT
    !16	Izana	28.309	-16.499	2373.0	8.0	759.4	ES
    !17	Jungfraujoch	46.548	7.987	3580.0	1.8	650.9	CH
    !18	Kosetice	49.583	15.083	534.0	-999.9	950.7	CR
    !19	K-Puszta	49.967	19.583	125.0	-999.9	998.3	HU
    !20	Mace Head	53.326	-9.899	5.0	10.0	1012.7	IR
    !21	Manaus	-2.583	60.200	100.0	40.0	-999.9	BR
    !22	Mauna Loa	19.539	-155.578	3397.0	12.2	666.4	US
    !23	Melpitz	51.533	12.900	87.0	-999.9	1002.9	DE
    !24	Monte Cimone	44.167	10.683	2165.0	5.0	778.8	IT
    !25	Mount Waliguan	36.283	100.900	3810.0	-999.9	631.8	PRC
    !26	Neumayer	-70.650	-8.250	42.0	7.5	1008.2	AQ
    !27	Pallas	67.974	24.116	560.0	7.0	948.0	FI
    !28	Payerne	46.817	6.950	490.0	-999.9	955.8	CH
    !29	Preila	55.210	21.040	5.0	-999.9	1012.9	LI
    !30	Puy de Dome	45.767	2.966	1465.0	12.0	849.2	FR
    !31	Sable Island	43.933	-60.017	4.0	10.0	1012.8	CA
    !32	Samoa	-14.232	-170.563	77.0	15.5	1008.2	US
    !33	Shang Dianzi	40.390	117.070	294.0	-999.9	978.4	PRC
    !34	Sonnblick	47.054	12.959	3106.0	-999.9	691.7	AU
    !35	South Pole	-89.997	-24.800	2841.0	13.3	715.4	AQ
    !36	Southern Great Plains	36.617	-97.500	320.0	10.0	975.4	US
    !37	Summit	72.580	-38.480	3238.0	-999.9	680.1	DK
    !38	Tahkuse	58.520	24.940	23.0	-999.9	1010.5	ET
    !39	Trinidad Head	41.050	-124.150	107.0	10.0	1000.5	US
    !40	Vavihill	56.017	13.150	172.0	5.0	992.8	SE
    !41	Zeppelin	78.908	11.881	474.0	7.0	947.8	NO
    !42	Zugspitze	47.417	10.983	2650.0	1.5	704.7	DE
    !43	Illmitz	47.767	16.767	117.0			Austra
    !44	Lille Valby	55.683	12.117	20.0			Denmark
    !45	Montseny	41.767	2.350	700.0			Spain
    !46	Backgarden	23.540	113.060				China
    !47	Cape San Juan	18.381	-65.618	65.0			USA

    DATA ZLAT/82.450, 58.800, 55.800, 71.323, 58.383, 40.050,           &
              51.971, -40.682, -34.353, -25.533, 35.333, 51.608,        &
              47.800, 61.850, 45.803, 28.309, 46.548, 49.583,           &
              49.967, 53.326, -2.583, 19.539, 51.533, 44.167,           &
              36.283, -70.650, 67.974, 46.817, 55.210, 45.767,          &
              43.933, -14.232, 40.390, 47.054, -89.997, 36.617,         &
              72.580, 58.520, 41.050, 56.017, 78.908, 47.417,           &
              47.767, 55.683, 41.767, 23.540, 18.381/


    DATA ZLON/-62.517, 17.380, -3.283, -156.612, 8.250, -88.367,        &
              4.927, 144.688, 18.490, 25.750, 25.670, -1.288,           &
              11.017, 24.283, 8.627, -16.499, 7.987, 15.083,            &
              19.583, -9.899, 60.200, -155.578, 12.900, 10.683,         &
              100.900, -8.250, 24.116, 6.950, 21.040, 2.966,            &
              -60.017, -170.563, 117.070, 12.959, -24.800, -97.500,     &
              -38.480, 24.940, -124.150, 13.150, 11.881, 10.983,        &
              16.767, 12.117, 2.350, 113.060, -65.618/

    MSITE = MS
    IF(MSITE.GT.MAXSITE) THEN
       WRITE(6,*)"MSITE>MAXSITE; NEED TO CHECK",MSITE,MAXSITE
       STOP
    ENDIF
    ISITES = 0
    JSITES = 0
    DO NSITE = 1, MSITE
       YLAT = ZLAT(NSITE)
       XLON = ZLON(NSITE)
       CALL SITEIJ(NSITE,XLON, YLAT,ISITE,JSITE,State_Grid)
       ISITES(NSITE) = ISITE
       JSITES(NSITE) = JSITE
    ENDDO

    ! output at selected sites
    IFSITEOUT = 0
    IFOUTIJ = 0
    SITEID  = 0
    DO NSITE=1,MSITE
       IFSITEOUT(NSITE) = 1
       ISITE = ISITES(NSITE)
       JSITE = JSITES(NSITE)
       IFOUTIJ(ISITE,JSITE)=1
       SITEID(ISITE,JSITE)=NSITE
    ENDDO

    WRITE(6,*)"SITE LOCATION AND ID"
    DO NSITE = 1, MSITE
       WRITE(6,25)NSITE,ZLAT(NSITE),ZLON(NSITE), &
                  ISITES(NSITE),JSITES(NSITE),IFSITEOUT(NSITE)
    ENDDO
25  FORMAT(I3,F8.2,F8.2,I4,I4,I2)

  end subroutine YSITESIJ_AC
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ysitesij_bacchus
!
! !DESCRIPTION: Function YSITESIJ_BACCHUS
!\\
!\\
! !INTERFACE:
!
  subroutine  YSITESIJ_BACCHUS( State_Grid )
!
! !USES:
!
    USE APM_INIT_MOD,   ONLY : MAXSITE,MSITE,ISITES,JSITES,IFSITEOUT
    USE APM_INIT_MOD,   ONLY : IFSITE
    USE State_Grid_Mod, ONLY : GrdState
!
! !INPUT PARAMETERS:
!
    TYPE(GrdState), INTENT(IN) :: State_Grid
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    !INTEGER, PARAMETER  :: MS   = 25
    INTEGER, PARAMETER :: MS   = 18
    REAL*8             :: ZLAT(MS),ZLON(MS), YLAT, XLON
    INTEGER            :: NSITE, ISITE, JSITE

    ! BACCHUS sites
    !51.95   4.93     !site1 Cabauw
    !35.33   25.67    !site2 Finokalia
    !46.53   7.98     !site3 Jungfraujoch
    !53.33   -9.9     !site4 Mace Head
    !51.54   12.93    !site5 Melpitz
    !37.75   137.60   !site6 NotoPeninsula
    !45.7    3.22     !site7 Puy de Dome
    !61.85   24.28    !site8 Hyytiala
    !56.02   13.15    !site9 Vavihill
    !
    ! Other sites
    !site10 Bondville (40.050, -88.367	213.0	10.0	987.9	US)
    !site11 SGP (36.8    -97.5)
    !site12 Shanghai China	31.23	121.47
    !site13 Beijing China	39.91	116.41
    !site14 Pinnacle State Park (42.10, -77.21)
    !site15 Chacaltaya  (-16.35, -68.13 )
    !Site16 SPL    (40.46 , -106.744 ) 3210 m asl
    !site17 LLN       (23.47,  120.87 )
    !site18 Rochester (43.16,  -77.61)

    DATA ZLAT/51.95,35.33,46.53,53.33,51.54,37.75,   &
              45.70,61.85,56.02,40.05,36.80,31.23,   &
              39.91,42.10,-16.35,40.46,23.47,43.16/

    DATA ZLON/4.93,25.67,7.98,-9.9,12.93,137.60,     &
              3.22,24.28,13.15,-88.367,-97.5,121.47, &
              116.41,-77.21,-68.13,-106.744,120.87,-77.61/

    MSITE = MS
    IF(MSITE.GT.MAXSITE) THEN
       WRITE(6,*)"MSITE>MAXSITE; NEED TO CHECK",MSITE,MAXSITE
       STOP
    ENDIF
    ISITES = 0
    JSITES = 0
    DO NSITE = 1, MSITE
       YLAT = ZLAT(NSITE)
       XLON = ZLON(NSITE)
       CALL SITEIJ(NSITE,XLON, YLAT,ISITE,JSITE,State_Grid)
       ISITES(NSITE) = ISITE
       JSITES(NSITE) = JSITE
    ENDDO

    ! output at selected sites
    IFSITEOUT = 0
    IFOUTIJ = 0
    SITEID  = 0
    DO NSITE=1,MSITE
       !IF(NSITE.GT.MS) THEN !No sites
       IF(NSITE.LE.MS) THEN !All sites
       !IF(NSITE.EQ.1.or.NSITE.EQ.4.or.NSITE.EQ.8. &
       !   or.NSITE.EQ.10.or.NSITE.EQ.12.
       !   or.NSITE.EQ.23.or.NSITE.EQ.24) THEN !Selected sites
          IFSITEOUT(NSITE) = 1
          ISITE = ISITES(NSITE)
          JSITE = JSITES(NSITE)
          IFOUTIJ(ISITE,JSITE)=1
          SITEID(ISITE,JSITE)=NSITE
       ENDIF
    ENDDO

    WRITE(6,*)"SITE LOCATION AND ID"
    DO NSITE = 1, MSITE
       WRITE(6,25)NSITE,ZLAT(NSITE),ZLON(NSITE), &
                  ISITES(NSITE),JSITES(NSITE),IFSITEOUT(NSITE)
    ENDDO
25  FORMAT(I3,F8.2,F8.2,I4,I4,I2)

  end subroutine YSITESIJ_BACCHUS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ysitesij_eu
!
! !DESCRIPTION: Function YSITESIJ_EU
!\\
!\\
! !INTERFACE:
!
  subroutine  YSITESIJ_EU( State_Grid )
!
! !USES:
!
    USE APM_INIT_MOD,   ONLY : MAXSITE,MSITE,ISITES,JSITES,IFSITEOUT
    USE APM_INIT_MOD,   ONLY : IFSITE
    USE State_Grid_Mod, ONLY : GrdState
!
! !INPUT PARAMETERS:
!
    TYPE(GrdState), INTENT(IN) :: State_Grid
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    !INTEGER, PARAMETER  :: MS   = 25
    INTEGER, PARAMETER  :: MS   = 24
    REAL*8         :: ZLAT(MS),ZLON(MS), YLAT, XLON
    INTEGER        :: NSITE, ISITE, JSITE

    ! EUCAARI sites
    !
    !Measurement site Site Geographic Altitude Site classification Surroundings
    !acronym coordinates (m a.s.l.)
    !Pallas, Finland PAL 67580 N, 24070 E 560 high-elevation, remote continental boreal forest, tundra
    !Hyyti¨al¨a, Finland HTL 61500 N, 24180 E 182 rural continental, background homogeneous boreal forest
    !Vavihill, Sweden VHL 56010 N, 139090 E 172 rural background, continental deciduous forest, field
    !Mace Head, Ireland MHD 53190 N, 09530W 5 marine background, coastal Atlantic ocean, bare land
    !Cabauw, Netherlands CBW 51570 N, 04530 E 0 clean marine, rural polluted ocean, field, urban
    !Melpitz, Germany MPZ 51320 N, 12540 E 87 rural polluted, continental pasture, suburban
    !Hohenpeissenberg, Germany HPB 47480 N, 11000 E 980 high-elevation, backgorund forest, meadows
    !K-Puszta, Hungary KPO 4664580 N, 1907350 E 125 rural continental, background field, deciduous forest
    !Jungfraujoch, Switzerland JFJ 4664320 N, 07570 E 3580 high altitude, background mountain
    !Puy de Dome, France PDD 45420 N, 03130 E 1465 high-elevation, background forest, mountain
    !San Pietro Capofiume, Italy SPC 44370 N, 11400 E 11 rural polluted, continental field, suburban
    !Finokalia, Greece FKL 35200 N, 2567400 E 250 marine background, coastal Mediterranean sea, dust

    ! Additional sites for GUAN
    !Bosel           rural lowland station                   53.00° N, 7.97° E
    !Waldhof         rural lowland station                   52.80° N 10.76° E
    !Melpitz         rural lowland station                   51.54° N 12.93° E (this is actually your site 8)
    !Schauinsland    rural mid-level mountain station        47.91° N  7.91° E
    !Hohenpeissenberg rural mid-level mountain station        47.80° N 11.00° E (this is actually your site 5)
    !Zugspitze        Alpine mountain station                 47.42° N 10.98° E
    !Kosetice         rural station in the Czech Republic     49° 35''N 15° 05''E (sorry, degrees here only!)

    !67.97   24.12    !site1 Pallas
    !56.02   13.15    !site2 Vavihill
    !51.95   4.88     !site3 Cabauw
    !61.85   24.28    !site4 Hyytiala
    !47.8    11.0     !site5 Hohenpeissenberg
    !46.97   19.60    !site6 K-Puszta
    !53.33   -9.9     !site7 Mace Head
    !51.54   12.93     !site8 Melpitz
    !46.53   7.95     !site9 Jungfraujoch
    !45.7    3.22     !site10 Puy de Dome
    !44.62   11.67    !site11 San Pietro Capofiume,
    !35.33   25.67    !site12 Finokalia
    !
    ! GUAN
    !50.0    7.97     !site13 Bosel
    !52.8    10.76    !site14 Waldhof
    !47.91   7.91     !site15 Schauinsland
    !47.42   10.98    !site16 Zugspitze
    !49.58   15.08    !site17 Kosetice

    ! EUSAAR
    !58.8, 17.38    !site 18 Aspvreten
    !58.38,8.25     !site 19 Birkenes
    !55.35,21.07    !site20 Preila
    !51.58,-1.32    !site21 Harwell
    !45.8, 8.63     !site 22 Ispra
    !44.18,10.7     !site23 Mt Cimone
    !42.17,23.58    !site24 BEO Moussala

    !52.80,10.76    !site21 Langenbrugge  !repeat 14

    DATA ZLAT/67.97,56.02,51.95,61.85,47.8,46.97,    &
              53.33,51.54,46.53,45.7,44.62,35.33,    &
              50.0,52.8,47.91,47.42,49.58,           &
              58.8,58.38,55.35,51.58,45.8,44.18,42.17/

    DATA ZLON/24.12,13.15,4.88,24.28,11.0,19.60,     &
              -9.9,12.93,7.95,3.22,11.67,25.67,      &
              7.97,10.76,7.91,10.98,15.08,           &
              17.38,8.25,21.07,-1.32,8.63,10.7,23.58/

    MSITE = MS
    IF(MSITE.GT.MAXSITE) THEN
       WRITE(6,*)"MSITE>MAXSITE; NEED TO CHECK",MSITE,MAXSITE
       STOP
    ENDIF
    ISITES = 0
    JSITES = 0
    DO NSITE = 1, MSITE
       YLAT = ZLAT(NSITE)
       XLON = ZLON(NSITE)
       CALL SITEIJ(NSITE,XLON, YLAT,ISITE,JSITE,State_Grid)
       ISITES(NSITE) = ISITE
       JSITES(NSITE) = JSITE
    ENDDO

    ! output at selected sites
    IFSITEOUT = 0
    IFOUTIJ = 0
    SITEID  = 0
    DO NSITE=1,MSITE
       !IF(NSITE.GT.MS) THEN !No sites
       !IF(NSITE.LE.MS) THEN !All sites
       IF(NSITE.EQ.1.or.NSITE.EQ.4.or.NSITE.EQ.8. &
          .or.NSITE.EQ.10.or.NSITE.EQ.12. &
          .or.NSITE.EQ.23.or.NSITE.EQ.24) THEN !Selected sites
          IFSITEOUT(NSITE) = 1
          ISITE = ISITES(NSITE)
          JSITE = JSITES(NSITE)
          IFOUTIJ(ISITE,JSITE)=1
          SITEID(ISITE,JSITE)=NSITE
       ENDIF
    ENDDO

    WRITE(6,*)"SITE LOCATION AND ID"
    DO NSITE = 1, MSITE
       WRITE(6,25)NSITE,ZLAT(NSITE),ZLON(NSITE), &
                  ISITES(NSITE),JSITES(NSITE),IFSITEOUT(NSITE)
    ENDDO
25  FORMAT(I3,F8.2,F8.2,I4,I4,I2)

  end subroutine YSITESIJ_EU
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ysitesij_ch
!
! !DESCRIPTION: Function YSITESIJ_CH
!\\
!\\
! !INTERFACE:
!
  subroutine  YSITESIJ_CH( State_Grid )
!
! !USES:
!
    USE APM_INIT_MOD,   ONLY : MAXSITE,MSITE,ISITES,JSITES,IFSITEOUT
    USE APM_INIT_MOD,   ONLY : IFSITE
    USE State_Grid_Mod, ONLY : GrdState
!
! !INPUT PARAMETERS:
!
    TYPE(GrdState), INTENT(IN) :: State_Grid
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER, PARAMETER :: MS   = 32
    REAL*8             :: ZLAT(MS),ZLON(MS), YLAT, XLON
    INTEGER            :: NSITE, ISITE, JSITE

    !Lat     Lon    Altitude
    !CN10LT.dat
    !36.28   100.9     !Mount Waliguan, China   !site1
    !29.44   79.62     !Mukteshwar, India       !site2
    !27.96   86.82   5079  !Pyramid, Himalayan Nepal !site3
    !Other sites
    !40.65  117.12        !Shang dian Zi of Beijing !site4
    !32.56  116.78   22.7 !DOE AMF Deployment, Shouxian, China !site5
    !28.43 77.156    243  !Gual Pahari, India (EUSSAR) !site6
    !28.6629,79.2563 2000 !Nainital, India (DOE AMF 2011) !site7

    ! Asia cities
    !Urumqi	 China	43.83	87.62  !site8
    !Shenyang China	41.81	123.43  !site9
    !Beijing China	39.91	116.41  !site10
    !TaiyuanChina	37.87	112.55  !site11
    !Seoul	 Korean	37.57	126.98  !site12
    !Anmyeon-doKorea36.37	126.32  !site13
    !Lanzhou China	36.06	103.83  !site14
    !Tokoyo	Japan	35.69	139.69  !site15
    !ShirahamaJapan	33.66	135.38  !site16
    !Shanghai China	31.23	121.47  !site17
    !Chengdu China	30.66	104.07  !site18
    !Wuhan	  China	30.59	114.31  !site19
    !Lhasa	 China	29.65	91.14  !site20
    !(New DelhiIndia	28.64	77.22  !site21 -- same as site 6)
    !Taipei	 China	25.09	121.56  !site21
    !Kunming China	25.04	102.72  !site22
    !GuangzhouChina	23.13	113.26  !site23
    !
    !Wuqing         39.39,  117.01   !site24
    !Nanjing        32.06,  118.80  !site25
    !Hongzhou       30.27,  120.16  !site26
    !Hefei          31.82,  117.23 !site27
    !Xian           34.34, 108.94  !site28
    !Jinan          36.65, 117.12  !site29
    !Shijiazhuang   38.04, 114.51   !site30
    !Changsha       28.23, 112.94   !site31
    !LLN            23.47,  120.87  !site32
    !

    DATA ZLAT/36.28,29.44,27.96,40.65,32.56,28.43,28.663,            &
              43.83,41.81,39.91,37.87,37.57,36.37,36.06,35.69,       &
              !33.66,31.23,30.66,30.59,29.65,28.64,25.09,25.04,23.13/
              33.66,31.23,30.66,30.59,29.65,25.09,25.04,23.13,       &
              39.39,32.06,30.27,31.82,34.34,36.65,38.04,28.23,23.47/

    DATA ZLON/100.9,79.62,86.82,117.12,116.78,77.156,79.256,                &
              87.62,123.43,116.41,112.55,126.98,126.32,103.83,139.69,       &
              !135.38,121.47,104.07,114.31,91.14,77.22,121.56,102.72,113.26/
              135.38,121.47,104.07,114.31,91.14,121.56,102.72,113.26,       &
              117.01,118.80,120.16,117.23,108.94,117.12,114.51,112.94,120.87/

    MSITE = MS
    IF(MSITE.GT.MAXSITE) THEN
       WRITE(6,*)"MSITE>MAXSITE; NEED TO CHECK",MSITE,MAXSITE
       STOP
    ENDIF
    ISITES = 0
    JSITES = 0
    DO NSITE = 1, MSITE
       YLAT = ZLAT(NSITE)
       XLON = ZLON(NSITE)
       CALL SITEIJ(NSITE,XLON, YLAT,ISITE,JSITE,State_Grid)
       ISITES(NSITE) = ISITE
       JSITES(NSITE) = JSITE
    ENDDO

    ! output at selected sites
    IFSITEOUT = 0
    IFOUTIJ = 0
    SITEID  = 0
    DO NSITE=1,MSITE
       !IF(NSITE.EQ.4.or.NSITE.EQ.13) THEN !Hyytiala,SGP
       !IF(NSITE.GE.45) THEN !No sites
       !IF(NSITE.EQ.1.or.NSITE.EQ.2.or.NSITE.EQ.4.or.  &
       ! NSITE.EQ.7.or.NSITE.EQ.8.or.NSITE.EQ.13.or.   &
       ! NSITE.EQ.15.or.NSITE.EQ.16.or.NSITE.EQ.18.or. &
       ! NSITE.EQ.20.or.NSITE.EQ.21.or.NSITE.EQ.27.or. &
       ! NSITE.EQ.34.or.NSITE.EQ.35.or.NSITE.EQ.36.or. &
       ! NSITE.EQ.41) THEN
       IF(NSITE.GE.1.and.NSITE.LE.MS) THEN
          IFSITEOUT(NSITE) = 1
          ISITE = ISITES(NSITE)
          JSITE = JSITES(NSITE)
          IFOUTIJ(ISITE,JSITE)=1
          SITEID(ISITE,JSITE)=NSITE
       ENDIF
    ENDDO

    WRITE(6,*)"SITE LOCATION AND ID"
    DO NSITE = 1, MSITE
       WRITE(6,25)NSITE,ZLAT(NSITE),ZLON(NSITE), &
                  ISITES(NSITE),JSITES(NSITE),IFSITEOUT(NSITE)
    ENDDO
25  FORMAT(I3,F8.2,F8.2,I4,I4,I2)

  end subroutine YSITESIJ_CH
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ysitesij_ch1
!
! !DESCRIPTION: Function YSITESIJ_CH1
!\\
!\\
! !INTERFACE:
!
  subroutine  YSITESIJ_CH1( State_Grid )
!
! !USES:
!
    USE APM_INIT_MOD,   ONLY : MAXSITE,MSITE,ISITES,JSITES,IFSITEOUT
    USE APM_INIT_MOD,   ONLY : IFSITE
    USE State_Grid_Mod, ONLY : GrdState
!
! !INPUT PARAMETERS:
!
    TYPE(GrdState), INTENT(IN) :: State_Grid
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER, PARAMETER :: MS   = 38
    REAL*8             :: ZLAT(MS),ZLON(MS), YLAT, XLON
    INTEGER            :: NSITE, ISITE, JSITE

    !Lat     Lon    Altitude
    !CN10LT.dat
    !36.28   100.9     !Mount Waliguan, China   !site1
    !29.44   79.62     !Mukteshwar, India       !site2
    !27.96   86.82   5079  !Pyramid, Himalayan Nepal !site3
    !Other sites
    !40.65  117.12        !Shang dian Zi of Beijing !site4
    !32.56  116.78   22.7 !DOE AMF Deployment, Shouxian, China !site5
    !28.43 77.156    243  !Gual Pahari, India (EUSSAR) !site6
    !28.6629,79.2563 2000 !Nainital, India (DOE AMF 2011) !site7

    ! Asia cities
    !Urumqi	 China	43.83	87.62  !site8
    !Shenyang China	41.81	123.43  !site9
    !Beijing China	39.91	116.41  !site10
    !TaiyuanChina	37.87	112.55  !site11
    !Seoul	 Korean	37.57	126.98  !site12
    !Anmyeon-doKorea36.37	126.32  !site13
    !Lanzhou China	36.06	103.83  !site14
    !Tokoyo	Japan	35.69	139.69  !site15
    !ShirahamaJapan	33.66	135.38  !site16
    !Shanghai China	31.23	121.47  !site17
    !Chengdu China	30.66	104.07  !site18
    !Wuhan	  China	30.59	114.31  !site19
    !Lhasa	 China	29.65	91.14  !site20
    !(New DelhiIndia	28.64	77.22  !site21 -- same as site 6)
    !!Taipei	25.09	121.56  !site21 -- same as site 35
    !Wufenshan	25.07	121.77  !site21
    !Kunming China	25.04	102.72  !site22
    !GuangzhouChina	23.13	113.26  !site23
    !
    !Wuqing         39.39,  117.01   !site24
    !Nanjing        32.06,  118.80  !site25
    !Hongzhou       30.27,  120.16  !site26
    !Hefei          31.82,  117.23 !site27
    !Xian           34.34, 108.94  !site28
    !Jinan          36.65, 117.12  !site29
    !Shijiazhuang   38.04, 114.51   !site30
    !Changsha       28.23, 112.94   !site31
    !LLN            23.47,  120.87  !site32
    !Pengchyu       25.63,  122.07  !site33
    !Fuguei Cape    25.30,  121.54   !site34
    !NCU            24.97,  121.19  !site35
    !HengChung      22.00,  120.71  !site36
    !Donshan        20.71,  116.73  !site37
    !Mt.Fuji(3776mASL) (35.36,138.73) !site38

    DATA ZLAT/36.28,29.44,27.96,40.65,32.56,28.43,28.663,            &
              43.83,41.81,39.91,37.87,37.57,36.37,36.06,35.69,       &
              !33.66,31.23,30.66,30.59,29.65,28.64,25.09,25.04,23.13/
              !33.66,31.23,30.66,30.59,29.65,25.09,25.04,23.13,      &
              33.66,31.23,30.66,30.59,29.65,25.07,25.04,23.13,       &
              39.39,32.06,30.27,31.82,34.34,36.65,38.04,28.23,23.47, &
              25.63, 25.30, 24.97, 22.00,20.71,35.36/

    DATA ZLON/100.9,79.62,86.82,117.12,116.78,77.156,79.256,                  &
              87.62,123.43,116.41,112.55,126.98,126.32,103.83,139.69,         &
              !135.38,121.47,104.07,114.31,91.14,77.22,121.56,102.72,113.26/
              !135.38,121.47,104.07,114.31,91.14,121.56,102.72,113.26,        &
              135.38,121.47,104.07,114.31,91.14,121.77,102.72,113.26,         &
              117.01,118.80,120.16,117.23,108.94,117.12,114.51,112.94,120.87, &
              122.07,121.54,121.19,120.71,116.73,138.73/

    MSITE = MS
    IF(MSITE.GT.MAXSITE) THEN
       WRITE(6,*)"MSITE>MAXSITE; NEED TO CHECK",MSITE,MAXSITE
       STOP
    ENDIF
    ISITES = 0
    JSITES = 0
    DO NSITE = 1, MSITE
       YLAT = ZLAT(NSITE)
       XLON = ZLON(NSITE)
       CALL SITEIJ(NSITE,XLON, YLAT,ISITE,JSITE,State_Grid)
       ISITES(NSITE) = ISITE
       JSITES(NSITE) = JSITE
    ENDDO

    ! output at selected sites
    IFSITEOUT = 0
    IFOUTIJ = 0
    SITEID  = 0
    DO NSITE=1,MSITE
       !IF(NSITE.EQ.4.or.NSITE.EQ.13) THEN !Hyytiala,SGP
       !IF(NSITE.GE.45) THEN !No sites
       !IF(NSITE.EQ.1.or.NSITE.EQ.2.or.NSITE.EQ.4.or.  &
       ! NSITE.EQ.7.or.NSITE.EQ.8.or.NSITE.EQ.13.or.   &
       ! NSITE.EQ.15.or.NSITE.EQ.16.or.NSITE.EQ.18.or. &
       ! NSITE.EQ.20.or.NSITE.EQ.21.or.NSITE.EQ.27.or. &
       ! NSITE.EQ.34.or.NSITE.EQ.35.or.NSITE.EQ.36.or. &
       ! NSITE.GE.32) THEN
       IF(NSITE.GE.1.and.NSITE.LE.MS) THEN
          IFSITEOUT(NSITE) = 1
          ISITE = ISITES(NSITE)
          JSITE = JSITES(NSITE)
          IFOUTIJ(ISITE,JSITE)=1
          SITEID(ISITE,JSITE)=NSITE
       ENDIF
    ENDDO

    WRITE(6,*)"SITE LOCATION AND ID"
    DO NSITE = 1, MSITE
       WRITE(6,25)NSITE,ZLAT(NSITE),ZLON(NSITE), &
                  ISITES(NSITE),JSITES(NSITE),IFSITEOUT(NSITE)
    ENDDO
25  FORMAT(I3,F8.2,F8.2,I4,I4,I2)

  end subroutine YSITESIJ_CH1
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ysitesij_na
!
! !DESCRIPTION: Function YSITESIJ_NA
!\\
!\\
! !INTERFACE:
!
  subroutine  YSITESIJ_NA( State_Grid )
!
! !USES:
!
    USE APM_INIT_MOD,   ONLY : MAXSITE,MSITE,ISITES,JSITES,IFSITEOUT
    USE APM_INIT_MOD,   ONLY : IFSITE
    USE State_Grid_Mod, ONLY : GrdState
!
! !INPUT PARAMETERS:
!
    TYPE(GrdState), INTENT(IN) :: State_Grid
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER, PARAMETER :: MS   = 31
    REAL*8             :: ZLAT(MS),ZLON(MS), YLAT, XLON
    INTEGER            :: NSITE, ISITE, JSITE

    ! For size distributions
    !site1 Bondville (40.1, -88.3  )
    !site2 Trinidad Head (41.05,-124.15  )
    !site3 WHI (Whistler), Canada (50.01, -122.95) 2182 masl
    !site4 ETL (East Trout Lake), Canada	(54.35, -104.99) 540 masl

    ! EPA sopersites
    !!site5 New York (40.7144,  -74.01 )
    !site5 Queens College, New York (40.74,  -73.82 )
    !site6 Atlanta, Georgia (33.75, -84.39)
    !!site7 Baltimore, Maryland (39.29,  -76.61)
    !site7 HUBC, Baltimore, Maryland (39.05,  -76.87)
    !site8 Pittsburgh, Pennsylvania (40.44,  -80.0 )
    !site9 APP, NC (36.213,  -81.692)
    !site10 EGB, Canada  (44.23,  -79.783)

    ! Other sites
    !site11 DukeForest (36.0, -79.08)
    !site12 Kent, Ohio (41.15, -81.36 )
    !site13 SGP (36.8    -97.5)
    !site14 Pinnacle State Park (42.10, -77.21)
    !site15 Whiteface Mtn  (44.4, -73.9)
    !Site16 SPL    (40.46 , -106.744 ) 3210 m asl
    !site17  MMSF,Indiana (39.32,  -86.42)
    !site18  Rochester (43.16,  -77.61)

    !site19 UMBS (45.56, -84.72)
    !site20  Boulder (40.02,  -105.27)
    !site21  Albany, NY (42.65, -73.76)

    !site22  Sweet Briar College (37.56,-79.08)
    !site23 Brent,AL (32.94,  -87.16)
    !site24 Houston, Texas (29.76,  -95.38)
    !site25 Los Angeles, California  (34.05,  -118.24)
    !site26 MexicoCity (19.43,  -99.13)
    !site27 Cool, CA (38.89, -120.97)
    !site28 Sacramento, CA (38.645, -121.34)
    !site29 Ozarks forest site (38.74, -92.2 )
    !site30 Brooklaven, NY (40.87, -72.89 )
    !site31 CPR (Cape San Juan), PR (18.48, -66.13)

    DATA ZLAT/40.1,41.05,50.01,54.35,40.74,33.75,39.05,   &
              40.44,36.213,44.23,36.0,41.15,36.8,42.10,   &
              44.4, 40.46,39.32,43.16,45.56,40.02,42.65,  &
              37.56,32.94,29.76,34.05,19.43,38.89,38.645, &
              38.74,40.87,18.48/

    DATA ZLON/-88.3,-124.15,-122.95,-104.99,-73.82,-84.39,-76.87,  &
              -80.0,-81.692,-79.783,-79.08,-81.36,-97.5,-77.21,    &
              -73.9,-106.744,-86.42,-77.61,-84.72,-105.27,-73.76,  &
              -79.08,-87.16,-95.38,-118.24,-99.13,-120.97,-121.34, &
              -92.2,-72.89,-66.13/

    MSITE = MS
    IF(MSITE.GT.MAXSITE) THEN
       WRITE(6,*)"MSITE>MAXSITE; NEED TO CHECK",MSITE,MAXSITE
       STOP
    ENDIF
    ISITES = 0
    JSITES = 0
    DO NSITE = 1, MSITE
       YLAT = ZLAT(NSITE)
       XLON = ZLON(NSITE)
       CALL SITEIJ(NSITE,XLON, YLAT,ISITE,JSITE,State_Grid)
       ISITES(NSITE) = ISITE
       JSITES(NSITE) = JSITE
    ENDDO

    ! output at selected sites
    IFSITEOUT = 0
    IFOUTIJ = 0
    SITEID  = 0
    DO NSITE=1,MSITE
       !IF(NSITE.EQ.4.or.NSITE.EQ.13.or.NSITE.EQ.27) THEN !
       IFSITEOUT(NSITE) = 1
       ISITE = ISITES(NSITE)
       JSITE = JSITES(NSITE)
       IFOUTIJ(ISITE,JSITE)=1
       SITEID(ISITE,JSITE)=NSITE
       !ENDIF
    ENDDO

    WRITE(6,*)"SITE LOCATION AND ID"
    DO NSITE = 1, MSITE
       WRITE(6,25)NSITE,ZLAT(NSITE),ZLON(NSITE), &
                  ISITES(NSITE),JSITES(NSITE),IFSITEOUT(NSITE)
    ENDDO
25  FORMAT(I3,F8.2,F8.2,I4,I4,I2)

  end subroutine YSITESIJ_NA
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ysitesij_aod_cn
!
! !DESCRIPTION: Function YSITESIJ_AOD_CN
!\\
!\\
! !INTERFACE:
!
  subroutine  YSITESIJ_AOD_CN( State_Grid )
!
! !USES:
!
    USE APM_INIT_MOD,   ONLY : MAXSITE,MSITE,ISITES,JSITES,IFSITEOUT
    USE APM_INIT_MOD,   ONLY : IFSITE
    USE State_Grid_Mod, ONLY : GrdState
!
! !INPUT PARAMETERS:
!
    TYPE(GrdState), INTENT(IN) :: State_Grid
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER, PARAMETER :: MS   = 100
    REAL*8             :: ZLAT(MS),ZLON(MS), YLAT, XLON
    INTEGER            :: NSITE, ISITE, JSITE

    ! Site info see C:\Users\yfq\Desktop\MyDocuments_Syn\Research\Data\AERONET\AOD_CN10_SITES.xls

    DATA ZLON/ -156.60, 24.28, 59.55, 37.51, 27.60, 20.79, 12.90,        &
               -104.70, 2.21, 28.82, -116.99, 12.51, -122.22, -79.75,    &
               4.88, 104.42, -70.95, 76.98, 12.65, 123.43, -88.37,       &
               -105.01, 116.38, -76.84, -7.91, -115.96, 54.27, -119.77,  &
               -97.49, 34.26, 3.22, 104.14, 25.28, 113.62, 69.22,        &
               -106.89, 58.68, 62.22, 135.36, 73.10, 44.23, -111.97,     &
               13.16, 37.13, 66.90, 120.22, -7.88, 46.08, 34.78,         &
               104.07, 114.31, 31.41, 47.52, 58.40, 79.52, -16.32,       &
               83.97, 86.81, 64.10, 75.81, 80.23, 84.13, 51.32,          &
               46.40, 67.03, 54.55, 58.28, 113.26, 72.63, 62.00,         &
               5.53, 39.77, -10.00, -155.58, -99.18, 73.81, -67.05,      &
               54.02, 47.10, -22.94, 15.00, -16.96, 2.66, -5.93,         &
               74.88, -1.40, 49.18, 4.34, -14.41, -56.10, 132.89,        &
               -170.57, 23.15, -70.31, -46.73, 31.59, 139.35, -58.50,    &
               144.69, -8.25/

    DATA ZLAT/ 71.32, 61.85, 57.04, 55.70, 53.92, 51.84, 51.50,          &
               50.28, 48.70, 47.00, 46.49, 45.31, 44.24, 44.23,          &
               43.93, 43.58, 43.11, 42.62, 41.84, 41.81, 40.05,          &
               40.05, 39.98, 38.99, 38.57, 38.50, 36.85, 36.79,          &
               36.61, 36.57, 36.32, 35.95, 35.33, 34.73, 34.55,          &
               34.35, 34.35, 34.22, 33.69, 33.62, 33.25, 33.07,          &
               32.67, 32.20, 32.10, 31.42, 31.21, 30.95, 30.86,          &
               30.66, 30.59, 30.12, 29.33, 29.10, 29.05, 28.48,          &
               28.15, 27.96, 26.97, 26.91, 26.51, 25.87, 25.20,          &
               24.91, 24.87, 24.25, 23.58, 23.13, 23.07, 23.00,          &
               22.79, 21.43, 20.00, 19.54, 19.33, 18.54, 17.97,          &
               17.67, 17.47, 16.73, 15.00, 14.39, 13.54, 13.28,          &
               12.92, 12.20, 11.28, 8.32, -7.98, -9.87, -12.66,          &
               -14.24, -15.25, -18.47, -23.56, -24.99, -25.90, -34.57,   &
               -40.68, -70.65/

    MSITE = MS
    IF(MSITE.GT.MAXSITE) THEN
       WRITE(6,*)"MSITE>MAXSITE; NEED TO CHECK",MSITE,MAXSITE
       STOP
    ENDIF
    ISITES = 0
    JSITES = 0
    DO NSITE = 1, MSITE
       YLAT = ZLAT(NSITE)
       XLON = ZLON(NSITE)
       CALL SITEIJ(NSITE,XLON, YLAT,ISITE,JSITE,State_Grid)
       ISITES(NSITE) = ISITE
       JSITES(NSITE) = JSITE
    ENDDO

    ! output at selected sites
    IFSITEOUT = 0
    IFOUTIJ = 0
    SITEID  = 0
    DO NSITE=1,MSITE
       !IF(NSITE.EQ.1) THEN  !SGP only
       IFSITEOUT(NSITE) = 1
       ISITE = ISITES(NSITE)
       JSITE = JSITES(NSITE)
       IFOUTIJ(ISITE,JSITE)=1
       SITEID(ISITE,JSITE)=NSITE
       !ENDIF
    ENDDO

    WRITE(6,*)"SITE LOCATION AND ID"
    DO NSITE = 1, MSITE
       WRITE(6,25)NSITE,ZLAT(NSITE),ZLON(NSITE), &
                  ISITES(NSITE),JSITES(NSITE),IFSITEOUT(NSITE)
    ENDDO
25  FORMAT(I3,F8.2,F8.2,I4,I4,I2)

  end subroutine YSITESIJ_AOD_CN
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: siteij
!
! !DESCRIPTION: Function SITEIJ finds I, J index for specific site (F. YU, 05/2009)
!\\
!\\
! !INTERFACE:
!
  subroutine SITEIJ(NSITE,XLON, YLAT,ISITE,JSITE,State_Grid)
!
! !USES:
!
    USE State_Grid_Mod, ONLY : GrdState
!
! !INPUT PARAMETERS:
!
    INTEGER,        INTENT(IN)  :: NSITE
    REAL*8,         INTENT(IN)  :: YLAT, XLON
    TYPE(GrdState), INTENT(IN)  :: State_Grid
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT) :: ISITE, JSITE
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES
!
    INTEGER :: I, J, IX, JY, L
    REAL*8  :: XEDGE1, XEDGE2, DLON1, YEDGE1, YEDGE2, DLAT1

    XEDGE1 = State_Grid%XEDGE(2,1)
    XEDGE2 = State_Grid%XEDGE(3,1)
    DLON1  = State_Grid%DX
    YEDGE1 = State_Grid%YEDGE(1,2)
    YEDGE2 = State_Grid%YEDGE(1,3)
    DLAT1  = State_Grid%DY

    ISITE = INT((XLON-XEDGE1)/DLON1)+2
    JSITE = INT((YLAT-YEDGE1)/DLAT1)+2
    IF(ISITE.LT.1.or.ISITE.GT.State_Grid%NX.or.JSITE.LT.1 &
       .or.JSITE.GT.State_Grid%NY) THEN
       WRITE(6,*)"SITE location out of nested boundary"
       IF(ISITE.LT.1) ISITE = 1
       IF(ISITE.GT.State_Grid%NX) ISITE = State_Grid%NX
       IF(JSITE.LT.1) JSITE = 1
       IF(JSITE.GT.State_Grid%NY) JSITE = State_Grid%NY
    ENDIF
    !WRITE(900+NSITE,520)ISITE, JSITE, XLON, YLAT,
    WRITE(6,520)ISITE, JSITE, XLON, YLAT,       &
                State_Grid%XEDGE(ISITE, JSITE), &
                State_Grid%XMID (ISITE, JSITE), &
                State_Grid%YEDGE(ISITE, JSITE), &
                State_Grid%YMID (ISITE, JSITE)
520 format(I4,I4,6F9.3, 10(1PE9.2))

  END SUBROUTINE SITEIJ
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get_local_time
!
! !DESCRIPTION: Function GET\_LOCAL\_TIME computes the local time and returns
!  an array of points where the local time is between two user-defined limits.
!  (bmy, 11/29/00, 12/10/08)
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GET_LOCAL_TIME( State_Grid )
!
! !USES:
!
    USE TIME_MOD,       ONLY : GET_LOCALTIME
    USE TIME_MOD,       ONLY : GET_TS_DYN
    USE State_Grid_Mod, ONLY : GrdState
!
! ! INPUT PARAMETERS:
!
    TYPE(GrdState), INTENT(IN) :: State_Grid
!
! !REVISION HISTORY:
!  29 Nov 2000 - R. Yantosca - Initial Version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: I
    REAL*8  :: LT, TS_DYN

    !=================================================================
    ! GET_LOCAL_TIME begins here!
    !=================================================================
    TS_DYN = GET_TS_DYN()
    TS_DYN = TS_DYN / 60e+0_fp

    GOOD = 0
    DO I = 1, State_Grid%NX

       ! Get local time
       LT = GET_LOCALTIME( I, 1, 1, State_Grid ) - TS_DYN
       IF ( LT < 0  ) LT = LT + 24e+0_fp

       ! GOOD indicates which boxes have local times between HR1 and HR2
       IF ( LT >= 13.d0 .and. LT <= 14.d0 ) THEN
          GOOD(I) = 1
       ENDIF
    ENDDO

  END SUBROUTINE GET_LOCAL_TIME
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_oh
!
! !DESCRIPTION: Function GET\_OH returns OH from State\_Chm%Species (for
!  coupled runs) or monthly mean OH (for offline runs).  Imposes a diurnal
!  variation on OH for offline simulations. (bmy, 12/16/02, 7/20/04)
!\\
!\\
! !INTERFACE:
!
  FUNCTION GET_OH( I, J, L, id_OH, Input_Opt, State_Chm, State_Met ) &
       RESULT( OH_MOLEC_CM3 )
!
! !USES:
!
    USE ERROR_MOD,          ONLY : ERROR_STOP
    USE Input_Opt_Mod,      ONLY : OptInput
    USE State_Chm_Mod,      ONLY : ChmState
    USE State_Met_Mod,      ONLY : MetState
    USE TIME_MOD,           ONLY : GET_TS_CHEM
!
! !INPUT PARAMETERS:
!
    INTEGER,        INTENT(IN) :: I, J, L     ! Lon, lat, level indices
    INTEGER,        INTENT(IN) :: id_OH
    TYPE(OptInput), INTENT(IN) :: Input_Opt   ! Input Options object
    TYPE(MetState), INTENT(IN) :: State_Met   ! Meteorology State object
    TYPE(ChmState), INTENT(IN) :: State_Chm   ! Chemistry State object
!
! !REVISION HISTORY:
!  16 Dec 2002 - R. Yantosca - Initial Version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    REAL(fp)            :: OH_MOLEC_CM3

    !=================================================================
    ! GET_OH begins here!
    !=================================================================

    IF ( Input_Opt%ITS_A_FULLCHEM_SIM ) THEN

       !---------------------
       ! Coupled simulation
       !---------------------

       ! Take OH from State_Chm%Species [molec/cm3]
       ! OH is defined only in the chemistry grid
       IF ( State_Met%InChemGrid(I,J,L) ) THEN
          OH_MOLEC_CM3 = State_Chm%Species(I,J,L,id_OH)
       ELSE
          OH_MOLEC_CM3 = 0e+0_fp
       ENDIF

    ELSE IF ( Input_Opt%ITS_AN_AEROSOL_SIM ) THEN

       !---------------------
       ! Offline simulation
       !---------------------

       ! Test for sunlight...
       IF ( State_Met%SUNCOS(I,J) > 0e+0_fp .and. TCOSZ(I,J) > 0e+0_fp ) THEN

          ! Impose a diurnal variation on OH during the day
          OH_MOLEC_CM3 = OH(I,J,L)                              * &
                         ( State_Met%SUNCOS(I,J) / TCOSZ(I,J) ) * &
                         ( 86400e+0_fp           / GET_TS_CHEM() )

          ! OH is in kg/m3 (from HEMCO), convert to molec/cm3 (mps, 9/18/14)
          OH_MOLEC_CM3 = OH_MOLEC_CM3 * XNUMOL_OH / CM3PERM3

          ! Make sure OH is not negative
          OH_MOLEC_CM3 = MAX( OH_MOLEC_CM3, 0e+0_fp )

       ELSE

          ! At night, OH goes to zero
          OH_MOLEC_CM3 = 0e+0_fp

       ENDIF

    ELSE

       !---------------------
       ! Invalid simulation
       !---------------------
       CALL ERROR_STOP( 'Invalid simulation!', 'GET_OH (apm_driv_mod.F90)')

    ENDIF

  END FUNCTION GET_OH
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ohno3time
!
! !DESCRIPTION: Subroutine OHNO3TIME computes the sum of cosine of the solar
!  zenith angle over a 24 hour day, as well as the total length of daylight.
!  This is needed to scale the offline OH and NO3 concentrations.
!  (rjp, bmy, 12/16/02, 3/30/04)
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE OHNO3TIME( State_Grid )
!
! !USES:
!
    USE TIME_MOD,       ONLY : GET_NHMSb,   GET_ELAPSED_SEC
    USE TIME_MOD,       ONLY : GET_TS_CHEM, GET_DAY_OF_YEAR, GET_GMT
    USE State_Grid_Mod, ONLY : GrdState
!
! !INPUT PARAMETERS:
!
    TYPE(GrdState), INTENT(IN) :: State_Grid
!
! !REVISION HISTORY:
!  16 Dec 2002 - R. Yantosca - Initial Version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    LOGICAL, SAVE       :: FIRST = .TRUE.
    INTEGER             :: I, J, L, N, NT, NDYSTEP
    REAL(fp)            :: A0, A1, A2, A3, B1, B2, B3
    REAL(fp)            :: LHR0, R, AHR, DEC, TIMLOC, YMID_R
    REAL(fp)            :: SUNTMP(State_Grid%NX,State_Grid%NY)

    !=================================================================
    ! OHNO3TIME begins here!
    !=================================================================

    !  Solar declination angle (low precision formula, good enough for us):
    A0 = 0.006918
    A1 = 0.399912
    A2 = 0.006758
    A3 = 0.002697
    B1 = 0.070257
    B2 = 0.000907
    B3 = 0.000148
    R  = 2.* PI * float( GET_DAY_OF_YEAR() - 1 ) / 365.

    DEC = A0 - A1*cos(  R) + B1*sin(  R) &
             - A2*cos(2*R) + B2*sin(2*R) &
             - A3*cos(3*R) + B3*sin(3*R)

    LHR0 = int(float( GET_NHMSb() )/10000.)

    ! Only do the following at the start of a new day
    IF ( FIRST .or. GET_GMT() < 1e-5 ) THEN

       ! Zero arrays
       TTDAY(:,:) = 0e+0_fp
       TCOSZ(:,:) = 0e+0_fp
       COSZM(:,:) = 0e+0_fp

       ! NDYSTEP is # of chemistry time steps in this day
       NDYSTEP = ( 24 - INT( GET_GMT() ) ) * 3600 / GET_TS_CHEM()

       ! NT is the elapsed time [s] since the beginning of the run
       NT = GET_ELAPSED_SEC()

       ! Loop forward through NDYSTEP "fake" timesteps for this day
       DO N = 1, NDYSTEP

          ! Zero SUNTMP array
          SUNTMP = 0e+0_fp

          ! Loop over surface grid boxes
          DO J = 1, State_Grid%NY
          DO I = 1, State_Grid%NX

             ! Grid box latitude center [radians]
             YMID_R = State_Grid%YMID_R(I,J)

             TIMLOC = real(LHR0) + real(NT)/3600.0 + State_Grid%XMID(I,J)/15.0

             DO WHILE (TIMLOC .lt. 0)
                TIMLOC = TIMLOC + 24.0
             ENDDO

             DO WHILE (TIMLOC .gt. 24.0)
                TIMLOC = TIMLOC - 24.0
             ENDDO

             AHR = abs(TIMLOC - 12.) * 15.0 * PI_180

             !===========================================================
             ! The cosine of the solar zenith angle (SZA) is given by:
             !
             !  cos(SZA) = sin(LAT)*sin(DEC) + cos(LAT)*cos(DEC)*cos(AHR)
             !
             ! where LAT = the latitude angle,
             !       DEC = the solar declination angle,
             !       AHR = the hour angle, all in radians.
             !
             ! If SUNCOS < 0, then the sun is below the horizon, and
             ! therefore does not contribute to any solar heating.
             !===========================================================

             ! Compute Cos(SZA)
             SUNTMP(I,J) = sin(YMID_R) * sin(DEC) + &
                           cos(YMID_R) * cos(DEC) * cos(AHR)

             ! TCOSZ is the sum of SUNTMP at location (I,J)
             ! Do not include negative values of SUNTMP
             TCOSZ(I,J) = TCOSZ(I,J) + MAX( SUNTMP(I,J), 0e+0_fp )

             ! COSZM is the peak value of SUMTMP during a day at (I,J)
             ! (rjp, bmy, 3/30/04)
             COSZM(I,J) = MAX( COSZM(I,J), SUNTMP(I,J) )

             ! TTDAY is the total daylight time at location (I,J)
             IF ( SUNTMP(I,J) > 0e+0_fp ) THEN
                TTDAY(I,J) = TTDAY(I,J) + DBLE( GET_TS_CHEM() )
             ENDIF
          ENDDO
          ENDDO

          !### Debug
          !PRINT*, '### IN OHNO3TIME'
          !PRINT*, '### N       : ', N
          !PRINT*, '### NDYSTEP : ', NDYSTEP
          !PRINT*, '### NT      : ', NT
          !PRINT*, '### JDAY    : ', JDAY
          !PRINT*, '### RLAT    : ', RLAT
          !PRINT*, '### XMID    : ', XMID
          !PRINT*, '### SUNTMP  : ', SUNTMP
          !PRINT*, '### TCOSZ   : ', MINVAL( TCOSZ ), MAXVAL( TCOSZ )
          !PRINT*, '### TTDAY   : ', MINVAL( TCOSZ ), MAXVAL( TCOSZ )

          ! Increment elapsed time [sec]
          NT = NT + GET_TS_CHEM()
       ENDDO

       ! Reset first-time flag
       FIRST = .FALSE.
    ENDIF

  END SUBROUTINE OHNO3TIME
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE:getalbdrr
!
! !DESCRIPTION: Function GETALBDDR interpolates or extrapolates to get albedo
!  at RRTMG wavelength from MODIS 7 bands albedo (F. YU, 08/2012)
!\\
!\\
! !INTERFACE:
!
  subroutine GETALBDRR(II,JJ,MALB,State_Met)
!
! !USES:
!
    USE State_Met_Mod, ONLY : MetState
!
! !INPUT PARAMETERS:
!
    TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
!
! !REMARKS:
! RRTMG WL(1-14): 3.46,2.79,2.33,2.05,1.78,1.46,1.27,
!                 1.01,0.70,0.53,0.39,0.30,0.23,8.02
!
! MODIS WL (1-7): 0.645, 0.859, 0.469, 0.550, 1.24, 1.64, 2.13
!
! !REVISION HISTORY:
!  Aug 2012 - F. Yu - Initial Version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER, PARAMETER :: KWL=7         !! mxy MODIS spectral bands
    INTEGER, INTENT(IN)    :: II,JJ
    REAL :: MALB(II,JJ,KWL)

    INTEGER :: I,J,IWL

    DO IWL=1,NBS
       ALBDRR(:,:,IWL) = State_Met%ALBD(:,:)
    ENDDO

    DO I=1,II
       DO J=1,JJ
          ! no extrapolation for now, should not have much effect
          IF(MALB(I,J,1).GT.0.0.and.MALB(I,J,2).GT.0.0) THEN
             ALBDRR(I,J,9)=MALB(I,J,1)+(MALB(I,J,2)-MALB(I,J,1))*0.055/0.214
          ENDIF

          IF(MALB(I,J,2).GT.0.0.and.MALB(I,J,5).GT.0.0) THEN
             ALBDRR(I,J,8)=MALB(I,J,2)+(MALB(I,J,5)-MALB(I,J,2))*0.142/0.381
          ENDIF

          IF(MALB(I,J,3).GT.0.) THEN
             ALBDRR(I,J,11) = MALB(I,J,3)
             ALBDRR(I,J,12) = MALB(I,J,3)
             ALBDRR(I,J,13) = MALB(I,J,3)
          ENDIF
          IF(MALB(I,J,4).GT.0.) ALBDRR(I,J,10)=MALB(I,J,4)

          IF(MALB(I,J,5).GT.0.) ALBDRR(I,J,7)=MALB(I,J,5)

          IF(MALB(I,J,5).GT.0.0.and.MALB(I,J,6).GT.0.0) THEN
             ALBDRR(I,J,6)=MALB(I,J,5)+(MALB(I,J,6)-MALB(I,J,5))*0.22/0.4
          ENDIF

          IF(MALB(I,J,6).GT.0.0.and.MALB(I,J,7).GT.0.0) THEN
             ALBDRR(I,J,5)=MALB(I,J,6)+(MALB(I,J,7)-MALB(I,J,6))*0.14/0.49
             ALBDRR(I,J,4)=MALB(I,J,6)+(MALB(I,J,7)-MALB(I,J,6))*0.41/0.49
          ENDIF
          IF(MALB(I,J,7).GT.0.) THEN
             ALBDRR(I,J,1) = MALB(I,J,7)
             ALBDRR(I,J,2) = MALB(I,J,7)
             ALBDRR(I,J,3) = MALB(I,J,7)
             ALBDRR(I,J,14) = MALB(I,J,7)
          ENDIF
          DO IWL=1,NBS
             IF(ALBDRR(I,J,IWL).LT.0.or.ALBDRR(I,J,IWL).GT.1.0) THEN
                WRITE(6,*)"STOP AT GETALBDRR: Need check"
                STOP
             ENDIF
          ENDDO
          !IF(MOD(I,5).EQ.1.and.MOD(J,5).EQ.1) THEN
          ! WRITE(1004,99) I,J,(ALBDRR(I,J,IWL),IWL=1,NBS)
          !ENDIF
          !99 FORMAT(I4,I4,30(F11.4))

       ENDDO
    ENDDO
    !flush(1004)

  end subroutine GETALBDRR
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_apm3d
!
! !DESCRIPTION: Subroutine INIT_APM3D allocates and zeroes module arrays.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE INIT_APM3D( Input_Opt, State_Grid, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod,       ONLY : OptInput
    USE State_Grid_Mod,      ONLY : GrdState
    USE APM_Init_Mod,        ONLY : NSO4, NSEA, NDSTB
    USE Module_Mosaic_Therm, ONLY : Load_Mosaic_Parameters
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)  :: Input_Opt    ! Input Options object
    TYPE(GrdState), INTENT(IN)  :: State_Grid   ! Grid State objectg
!
! ! OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT) :: RC           ! Success or failure?
!
! !REVISION HISTORY:
!  28 Aug 2008 - F. Yu       - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: NX, NY, NZ

    !=================================================================
    ! INIT_AEROSOL begins here!
    !=================================================================
    NX = State_Grid%NX
    NY = State_Grid%NY
    NZ = State_Grid%NZ

    ALLOCATE( T3DAPM(NX,NY,20,24,6), STAT=RC )
    CALL GC_CheckVar( 'apm_driv_mod.F%T3DAPM', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    T3DAPM = 0e0

    ALLOCATE( RH3DAPM(NX,NY,20,24,6), STAT=RC )
    CALL GC_CheckVar( 'apm_driv_mod.F%RH3DAPM', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    RH3DAPM = 0e0

    ALLOCATE( PBLH2DAPM(NX,NY,24,6), STAT=RC )
    CALL GC_CheckVar( 'apm_driv_mod.F%PBLH2DAPM', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    PBLH2DAPM = 0e0

    ALLOCATE( EMITNH3(NX,NY,NZ), STAT=RC )
    CALL GC_CheckVar( 'apm_driv_mod.F%EMITNH3', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    EMITNH3 = 0d0

    ALLOCATE( EMITSO2(NX,NY,NZ), STAT=RC )
    CALL GC_CheckVar( 'apm_driv_mod.F%EMITSO2', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    EMITSO2 = 0d0

    ALLOCATE( DRYDEP(NX,NY,NZ,4), STAT=RC )
    CALL GC_CheckVar( 'apm_driv_mod.F%DRYDEP', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    DRYDEP = 0d0

    ALLOCATE( WETDEP(NX,NY,NZ,4), STAT=RC )
    CALL GC_CheckVar( 'apm_driv_mod.F%WETDEP', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    WETDEP = 0d0

    ALLOCATE( CONDEP(NX,NY,NZ,4), STAT=RC )
    CALL GC_CheckVar( 'apm_driv_mod.F%CONDEP', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    CONDEP = 0d0

    ALLOCATE( UPTAKE(NX,NY,NZ,4), STAT=RC )
    CALL GC_CheckVar( 'apm_driv_mod.F%UPTAKE', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    UPTAKE = 0d0

    ALLOCATE( OXIDAT(NX,NY,NZ,4), STAT=RC )
    CALL GC_CheckVar( 'apm_driv_mod.F%OXIDAT', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    OXIDAT = 0d0

    ALLOCATE( FCLOUD(NX,NY,NZ,NSO4+4), STAT=RC )
    CALL GC_CheckVar( 'apm_driv_mod.F%FCLOUD', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    FCLOUD = 0d0

    ALLOCATE( AERAREA(NX,NY,NZ,NSO4+NSEA+NDSTB+4), STAT=RC )
    CALL GC_CheckVar( 'apm_driv_mod.F%AERAREA', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    AERAREA = 0d0

    ALLOCATE( AERDRYR(NX,NY,NZ,NSO4+NSEA+NDSTB+4), STAT=RC )
    CALL GC_CheckVar( 'apm_driv_mod.F%AERDRYR', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    AERDRYR = 0d0

    ALLOCATE( GAMMAPM(NX,NY,NZ,NTYP), STAT=RC )
    CALL GC_CheckVar( 'apm_driv_mod.F%GAMMAPM', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    GAMMAPM = 0d0

    ALLOCATE( XN4D(NX,NY,NZ,NSO4), STAT=RC )
    CALL GC_CheckVar( 'apm_driv_mod.F%XN4D', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    XN4D = 0d0

    ALLOCATE( XQ3D(NX,NY,NZ), STAT=RC )
    CALL GC_CheckVar( 'apm_driv_mod.F%XQ3D', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    XQ3D = 0d0

    ALLOCATE( RCLDL3D(NX,NY,NZ), STAT=RC )
    CALL GC_CheckVar( 'apm_driv_mod.F%RCLDL3D', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    RCLDL3D = 0d0

    ALLOCATE( CDN3D(NX,NY,NZ), STAT=RC )
    CALL GC_CheckVar( 'apm_driv_mod.F%CDN3D', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    CDN3D = 0d0

    ALLOCATE( RCLDI3D(NX,NY,NZ), STAT=RC )
    CALL GC_CheckVar( 'apm_driv_mod.F%RCLDI3D', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    RCLDI3D = 0d0

    ALLOCATE( RCLDI2D(NX,NY), STAT=RC )
    CALL GC_CheckVar( 'apm_driv_mod.F%RCLDL2D', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    RCLDI2D = 0d0

    ALLOCATE( IN3D(NX,NY,NZ), STAT=RC )
    CALL GC_CheckVar( 'apm_driv_mod.F%CDN3D', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    IN3D = 0d0

    ALLOCATE( MASSISRP(NX,NY,NZ,4), STAT=RC )
    CALL GC_CheckVar( 'apm_driv_mod.F%MASSISRP', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    MASSISRP = 0d0

    ALLOCATE( MASSMESA(NX,NY,NZ,4), STAT=RC )
    CALL GC_CheckVar( 'apm_driv_mod.F%MASSMESA', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    MASSMESA = 0d0

    ALLOCATE( GFTOT3D(NX,NY,NZ,3), STAT=RC )
    CALL GC_CheckVar( 'apm_driv_mod.F%GFTOT3D', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    GFTOT3D = 1d0

    ALLOCATE( DENWET3D(NX,NY,NZ,2), STAT=RC )
    CALL GC_CheckVar( 'apm_driv_mod.F%DENWET3D', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    DENWET3D = 2d0

    ALLOCATE( MWSIZE3D(NX,NY,NZ,3), STAT=RC )
    CALL GC_CheckVar( 'apm_driv_mod.F%MWSIZE3D', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    MWSIZE3D = 1.D-6

    ALLOCATE( SO2toSO4( NX, NY, NZ ), STAT=RC )
    CALL GC_CheckVar( 'apm_driv_mod.F%SO2toSO4', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    SO2toSO4 = 0d0

    ALLOCATE( NCOAG3D(NX,NY,NZ,2), STAT=RC )
    CALL GC_CheckVar( 'apm_driv_mod.F%NCOAG3D', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    NCOAG3D = 0d0

    ALLOCATE( IACT1(NX,NY,NZ), STAT=RC )
    CALL GC_CheckVar( 'apm_driv_mod.F%IACT1', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    IACT1 = 15

    ALLOCATE( IACT2(NX,NY,NZ), STAT=RC )
    CALL GC_CheckVar( 'apm_driv_mod.F%IACT2', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    IACT2 = 15

    ALLOCATE( IACT3(NX,NY,NZ), STAT=RC )
    CALL GC_CheckVar( 'apm_driv_mod.F%IACT3', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    IACT3 = 15

    ALLOCATE( PSO4GAS( NX, NY, NZ ), STAT=RC )
    CALL GC_CheckVar( 'apm_driv_mod.F%PSO4GAS', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    PSO4GAS = 0.d0

    ALLOCATE( XO3( NX, NY, NZ ), STAT=RC )
    CALL GC_CheckVar( 'apm_driv_mod.F%XO3', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    XO3 = 0.d0

    ALLOCATE( PLVSOG( NX, NY, NZ, 5 ), STAT=RC )
    CALL GC_CheckVar( 'apm_driv_mod.F%PLVSOG', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    PLVSOG = 0d0

    ALLOCATE( BCLIFE( NX, NY, NZ ), STAT=RC )
    CALL GC_CheckVar( 'apm_driv_mod.F%BCLIFE', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    PLVSOG = 0d0

    ALLOCATE( OCLIFE( NX, NY, NZ ), STAT=RC )
    CALL GC_CheckVar( 'apm_driv_mod.F%OCLIFE', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    PLVSOG = 0d0

    ALLOCATE( IFOUTIJ( NX, NY ) )
    CALL GC_CheckVar( 'apm_driv_mod.F%IFOUTIJ', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    IFOUTIJ = 0

    ALLOCATE( SITEID( NX, NY) )
    CALL GC_CheckVar( 'apm_driv_mod.F%SITEID', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    SITEID = 0

    !Yu+ seasalt SST correction 10/15/2011
    ALLOCATE( FSST( NX, NY), STAT=RC )
    CALL GC_CheckVar( 'apm_driv_mod.F%FSST', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    FSST = 1.d0

    !Yu+ 6/1/11
    ALLOCATE( SPGF( NX, NY, NZ ), STAT=RC )
    CALL GC_CheckVar( 'apm_driv_mod.F%SPGF', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    SPGF = 1.d0

    ALLOCATE( TAONH3( NX, NY, NZ ), STAT=RC )
    CALL GC_CheckVar( 'apm_driv_mod.F%TAONH3', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    TAONH3 = 1.d0

    ALLOCATE( YCOD3D( NX, NY, NZ ), STAT=RC )
    CALL GC_CheckVar( 'apm_driv_mod.F%YCOD3D', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    YCOD3D = 0.d0

    ALLOCATE( TCOD3D( NX, NY, NZ, 9 ), STAT=RC )
    CALL GC_CheckVar( 'apm_driv_mod.F%TCOD3D', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    TCOD3D = 0.d0
    TCOD3D(:,:,:,9) = 1.d-3

    ALLOCATE( YHTRC3D( NX, NY, NZ ), STAT=RC )
    CALL GC_CheckVar( 'apm_driv_mod.F%YHTRC3D', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    YHTRC3D = 0.d0

    ALLOCATE( YHTR3D( NX, NY, NZ ), STAT=RC )
    CALL GC_CheckVar( 'apm_driv_mod.F%YHTR3D', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    YHTR3D = 0.d0

    ALLOCATE( YHTRC03D( NX, NY, NZ ), STAT=RC )
    CALL GC_CheckVar( 'apm_driv_mod.F%YHTRC03D', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    YHTRC03D = 0.d0

    ALLOCATE( YHTR03D( NX, NY, NZ ), STAT=RC )
    CALL GC_CheckVar( 'apm_driv_mod.F%YHTR03D', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    YHTR03D = 0.d0

    ALLOCATE( ZBEXT3D(NX,NY,NZ,MWL), STAT=RC )
    CALL GC_CheckVar( 'apm_driv_mod.F%ZBEXT3D', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    ZBEXT3D = 0.d0

    ALLOCATE( ZW3D(NX,NY,NZ,MWL), STAT=RC )
    CALL GC_CheckVar( 'apm_driv_mod.F%ZW3D', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    ZW3D= 0.d0

    ALLOCATE( ZG3D(NX,NY,NZ,MWL), STAT=RC )
    CALL GC_CheckVar( 'apm_driv_mod.F%ZG3D', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    ZG3D= 1.d0

    ALLOCATE( YBEXT3D(NX,NY,NZ,NTYP,MWL), STAT=RC )
    CALL GC_CheckVar( 'apm_driv_mod.F%YBEXT3D', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    YBEXT3D= 0.d0

    ALLOCATE( XBEXT1k3D(NX,NY,NZ,12), STAT=RC )
    CALL GC_CheckVar( 'apm_driv_mod.F%XBEXT1K3D', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    XBEXT1k3D= 0.d0

    ALLOCATE( YW3D(NX,NY,NZ,NTYP,MWL), STAT=RC )
    CALL GC_CheckVar( 'apm_driv_mod.F%YW3D', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    YW3D= 0.d0

    ALLOCATE( YG3D(NX,NY,NZ,NTYP,MWL), STAT=RC )
    CALL GC_CheckVar( 'apm_driv_mod.F%YG3D', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    YG3D= 1.d0

    ALLOCATE( ZBABS3D(NX,NY,NZ,NWL), STAT=RC )
    CALL GC_CheckVar( 'apm_driv_mod.F%ZBABS3D', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    ZBABS3D= 0.d0

    ALLOCATE( YUV(NX,NY,NZ,9), STAT=RC )
    CALL GC_CheckVar( 'apm_driv_mod.F%YUV', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    YUV= 0.d0

    ! Array denoting where LT is between HR1 and HR2
    ALLOCATE( GOOD( NX ), STAT=RC )
    CALL GC_CheckVar( 'apm_driv_mod.F%GOOD', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    GOOD = 0

    ALLOCATE( CODOUT(NX,NY,NZ,10), STAT=RC )
    CALL GC_CheckVar( 'apm_driv_mod.F%CODOUT', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    CODOUT= 0.d0

    ALLOCATE( CODOUTNUM(NX,NY,NZ), STAT=RC )
    CALL GC_CheckVar( 'apm_driv_mod.F%CODOUTNUM', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    CODOUTNUM= 0.d0

    ALLOCATE( MODISOUT(NX,NY,6), STAT=RC )
    CALL GC_CheckVar( 'apm_driv_mod.F%MODISOUT', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    MODISOUT= 0.d0

    ALLOCATE( TCOSZ( NX, NY ), STAT=RC )
    CALL GC_CheckVar( 'apm_driv_mod.F%TCOSZ', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    TCOSZ = 0e+0_fp

    ALLOCATE( TTDAY( NX, NY ), STAT=RC )
    CALL GC_CheckVar( 'apm_driv_mod.F%TTDAY', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    TTDAY = 0e+0_fp

    ALLOCATE( COSZM( NX, NY ), STAT=RC )
    CALL GC_CheckVar( 'apm_driv_mod.F%COSZM', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    COSZM = 0e+0_fp

    ALLOCATE( ALBDRR(NX,NY,NBS), STAT=RC )
    CALL GC_CheckVar( 'apm_driv_mod.F%ALBDRR', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    ALBDRR = 0.0D0

    ALLOCATE( ZTCST(NX,NY,NTYP), STAT=RC )
    CALL GC_CheckVar( 'apm_driv_mod.F%ZTCST', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    ZTCST = 0.0D0

    ALLOCATE( ZTCSB(NX,NY,NTYP), STAT=RC )
    CALL GC_CheckVar( 'apm_driv_mod.F%ZTCSB', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    ZTCSB = 0.0D0

    ALLOCATE( ZTFST(NX,NY,NTYP), STAT=RC )
    CALL GC_CheckVar( 'apm_driv_mod.F%ZTFST', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    ZTFST = 0.0D0

    !ALLOCATE( YG3D, STAT=RC )
    !CALL GC_CheckVar( 'apm_driv_mod.F%YG3D', 0, RC )
    !IF ( RC /= GC_SUCCESS ) RETURN
    !YG3D = 0.0D0

    ALLOCATE( ZTFSA(NX,NY,NTYP), STAT=RC )
    CALL GC_CheckVar( 'apm_driv_mod.F%ZTFSA', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    ZTFSA = 0.0D0

    ALLOCATE( ZTCLT(NX,NY,NTYP), STAT=RC )
    CALL GC_CheckVar( 'apm_driv_mod.F%ZTCLT', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    ZTCLT = 0.0D0

    ALLOCATE( ZTCLB(NX,NY,NTYP), STAT=RC )
    CALL GC_CheckVar( 'apm_driv_mod.F%ZTCLB', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    ZTCLB = 0.0D0

    ALLOCATE( ZTFUL(NX,NY,NTYP), STAT=RC )
    CALL GC_CheckVar( 'apm_driv_mod.F%', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    ALBDRR = 0.0D0

    ALLOCATE( ZTFDL(NX,NY,NTYP), STAT=RC )
    CALL GC_CheckVar( 'apm_driv_mod.F%', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    ALBDRR = 0.0D0

    ALLOCATE( ZTFLA(NX,NY,NTYP), STAT=RC )
    CALL GC_CheckVar( 'apm_driv_mod.F%ZTFLA', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    ZTFLA = 0.0D0

    ALLOCATE( ZTAOD(NX,NY,NTYP), STAT=RC )
    CALL GC_CheckVar( 'apm_driv_mod.F%ZTAOD', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    ZTAOD = 0.0D0

    ALLOCATE( ZCLDF(NX,NY), STAT=RC )
    CALL GC_CheckVar( 'apm_driv_mod.F%ZCLDF', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    ZCLDF = 0.0D0

    ALLOCATE( ZALB(NX,NY), STAT=RC )
    CALL GC_CheckVar( 'apm_driv_mod.F%ZALB', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    ZALB = 0.0D0

    ALLOCATE( ZMALB(NX,NY,NBS), STAT=RC )
    CALL GC_CheckVar( 'apm_driv_mod.F%ZMALB', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    ZMALB = 0.0D0

    ALLOCATE( ZCST(NX,NY), STAT=RC )
    CALL GC_CheckVar( 'apm_driv_mod.F%ZCST', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    ZCST = 0.0D0

    ALLOCATE( ZCSB(NX,NY), STAT=RC )
    CALL GC_CheckVar( 'apm_driv_mod.F%ZCSB', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    ZCSB = 0.0D0

    ALLOCATE( ZFST(NX,NY), STAT=RC )
    CALL GC_CheckVar( 'apm_driv_mod.F%ZFST', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    ZFST = 0.0D0

    ALLOCATE( ZCLD(NX,NY), STAT=RC )
    CALL GC_CheckVar( 'apm_driv_mod.F%ZCLD', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    ZCLD = 0.0D0

    ALLOCATE( ZCLD0(NX,NY), STAT=RC )
    CALL GC_CheckVar( 'apm_driv_mod.F%ZCLD0', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    ZCLD0 = 0.0D0

    ALLOCATE( ZFSB(NX,NY), STAT=RC )
    CALL GC_CheckVar( 'apm_driv_mod.F%ZFSB', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    ZFSB = 0.0D0

    ALLOCATE( ZFSA(NX,NY), STAT=RC )
    CALL GC_CheckVar( 'apm_driv_mod.F%ZFSA', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    ZFSA = 0.0D0

    !ALLOCATE( ZCLT(NX,NY), STAT=RC )
    !CALL GC_CheckVar( 'apm_driv_mod.F%', 0, RC )
    !IF ( RC /= GC_SUCCESS ) RETURN
    !ALBDRR = 0.0D0

    !ALLOCATE( ZCLB(NX,NY), STAT=RC )
    !CALL GC_CheckVar( 'apm_driv_mod.F%', 0, RC )
    !IF ( RC /= GC_SUCCESS ) RETURN
    !ALBDRR = 0.0D0

    !ALLOCATE( ZFUL(NX,NY), STAT=RC )
    !CALL GC_CheckVar( 'apm_driv_mod.F%', 0, RC )
    !IF ( RC /= GC_SUCCESS ) RETURN
    !ALBDRR = 0.0D0

    !ALLOCATE( ZFDL(NX,NY), STAT=RC )
    !CALL GC_CheckVar( 'apm_driv_mod.F%', 0, RC )
    !IF ( RC /= GC_SUCCESS ) RETURN
    !ALBDRR = 0.0D0

    !ALLOCATE( ZFLA(NX,NY), STAT=RC )
    !CALL GC_CheckVar( 'apm_driv_mod.F%', 0, RC )
    !IF ( RC /= GC_SUCCESS ) RETURN
    !ALBDRR = 0.0D0

    ALLOCATE( ZAOD(NX,NY), STAT=RC )
    CALL GC_CheckVar( 'apm_driv_mod.F%ZAOD', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    ZAOD = 0.0D0

    ALLOCATE( ZCOD(NX,NY), STAT=RC )
    CALL GC_CheckVar( 'apm_driv_mod.F%ZCOD', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    ZCOD = 0.0D0

    ALLOCATE( ZCODGC(NX,NY), STAT=RC )
    CALL GC_CheckVar( 'apm_driv_mod.F%ZCODGC', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    ZCODGC = 0.0D0

    ALLOCATE( ZAODOUT1(NX,NY), STAT=RC )
    CALL GC_CheckVar( 'apm_driv_mod.F%ZADOUT1', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    ZAODOUT1 = 0.0D0

    ALLOCATE( ZAODOUT3(NX,NY), STAT=RC )
    CALL GC_CheckVar( 'apm_driv_mod.F%ZAODOUT3', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    ZAODOUT3 = 0.0D0

    ALLOCATE( ZAOD25(NX,NY), STAT=RC )
    CALL GC_CheckVar( 'apm_driv_mod.F%ZAOD25', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    ZAOD25 = 0.0D0

    ALLOCATE( ZAOD50(NX,NY), STAT=RC )
    CALL GC_CheckVar( 'apm_driv_mod.F%ZAOD50', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    ZAOD50 = 0.0D0

    ALLOCATE( ZAOD75(NX,NY), STAT=RC )
    CALL GC_CheckVar( 'apm_driv_mod.F%ZAOD75', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    ZAOD75 = 0.0D0

    ALLOCATE( NAOD25(NX,NY), STAT=RC )
    CALL GC_CheckVar( 'apm_driv_mod.F%NAOD25', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    NAOD25 = 0.0D0

    ALLOCATE( NAOD50(NX,NY), STAT=RC )
    CALL GC_CheckVar( 'apm_driv_mod.F%NAOD50', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    NAOD50 = 0.0D0

    ALLOCATE( NAOD75(NX,NY), STAT=RC )
    CALL GC_CheckVar( 'apm_driv_mod.F%NAOD75', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    NAOD75 = 0.0D0

    ALLOCATE( ZVIS(NX,NY), STAT=RC )
    CALL GC_CheckVar( 'apm_driv_mod.F%ZVIS', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    ZVIS = 0.0D0

    ALLOCATE( THAZ(NX,NY), STAT=RC )
    CALL GC_CheckVar( 'apm_driv_mod.F%THAZ', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    THAZ = 0.0D0

    ALLOCATE( TFOG(NX,NY), STAT=RC )
    CALL GC_CheckVar( 'apm_driv_mod.F%TFOG', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    TFOG = 0.0D0

    ALLOCATE( ZAAOD(NX,NY), STAT=RC )
    CALL GC_CheckVar( 'apm_driv_mod.F%ZAAOD', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    ZAAOD = 0.0D0

    ALLOCATE( ZAAODOUT1(NX,NY), STAT=RC )
    CALL GC_CheckVar( 'apm_driv_mod.F%ZAAODOUT1', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    ZAAODOUT1 = 0.0D0

    ALLOCATE( ZAAODOUT3(NX,NY), STAT=RC )
    CALL GC_CheckVar( 'apm_driv_mod.F%ZAAODOUT3', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    ZAAODOUT3 = 0.0D0

    ALLOCATE( ZABS(NX,NY), STAT=RC )
    CALL GC_CheckVar( 'apm_driv_mod.F%ZABS', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    ZABS = 0.0D0

    ALLOCATE( ZWCL(NX,NY), STAT=RC )
    CALL GC_CheckVar( 'apm_driv_mod.F%ZWCL', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    ZWCL = 0.0D0

    ALLOCATE( ZWCI(NX,NY), STAT=RC )
    CALL GC_CheckVar( 'apm_driv_mod.F%ZWCI', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    ZWCI = 0.0D0

    ALLOCATE( AOD(NX,NY), STAT=RC )
    CALL GC_CheckVar( 'apm_driv_mod.F%AOD', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    AOD = 0.0D0

    ALLOCATE( AAOD(NX,NY), STAT=RC )
    CALL GC_CheckVar( 'apm_driv_mod.F%AAOD', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    AAOD = 0.0D0

    ALLOCATE( AODOUT1(NX,NY), STAT=RC )
    CALL GC_CheckVar( 'apm_driv_mod.F%AODOUT1', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    AODOUT1 = 0.0D0

    ALLOCATE( AAODOUT1(NX,NY), STAT=RC )
    CALL GC_CheckVar( 'apm_driv_mod.F%AAODOUT1', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    AAODOUT1 = 0.0D0

    ALLOCATE( AODOUT3(NX,NY), STAT=RC )
    CALL GC_CheckVar( 'apm_driv_mod.F%', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    ALBDRR = 0.0D0

    ALLOCATE( AAODOUT3(NX,NY), STAT=RC )
    CALL GC_CheckVar( 'apm_driv_mod.F%AAODOUT3', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    AAODOUT3 = 0.0D0

    ALLOCATE( CODGC(NX,NY), STAT=RC )
    CALL GC_CheckVar( 'apm_driv_mod.F%CODGC', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    CODGC = 0.0D0

    ALLOCATE( TAOD(NX,NY,NTYP), STAT=RC )
    CALL GC_CheckVar( 'apm_driv_mod.F%TAOD', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    TAOD = 0.0D0

    ! sets up indices and other stuff once per simulation
    CALL load_mosaic_parameters()

  END SUBROUTINE INIT_APM3D
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cleanup_apm3d
!
! !DESCRIPTION: Subroutine CLEANUP\_APM3D deallocates all module arrays.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE CLEANUP_APM3D( Input_Opt, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod, ONLY : OptInput
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)  :: Input_Opt    ! Input Options object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT) :: RC           ! Success or failure?
!
! !REVISION HISTORY:
!  28 Aug 2008 - F. Yu       - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
    !=================================================================
    ! CLEANUP_APM3D begins here!
    !=================================================================
    IF ( ALLOCATED( T3DAPM ) ) THEN
       DEALLOCATE( T3DAPM, STAT=RC )
       CALL GC_CheckVar( 'apm_driv_mod.F%T3DAPM', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ALLOCATED( RH3DAPM ) ) THEN
       DEALLOCATE( RH3DAPM, STAT=RC )
       CALL GC_CheckVar( 'apm_driv_mod.F%RH3DAPM', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ALLOCATED( PBLH2DAPM ) ) THEN
       DEALLOCATE( PBLH2DAPM, STAT=RC )
       CALL GC_CheckVar( 'apm_driv_mod.F%PBLH2DAPM', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ALLOCATED( EMITNH3 ) ) THEN
       DEALLOCATE( EMITNH3, STAT=RC )
       CALL GC_CheckVar( 'apm_driv_mod.F%EMITNH3', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ALLOCATED( EMITSO2 ) ) THEN
       DEALLOCATE( EMITSO2, STAT=RC )
       CALL GC_CheckVar( 'apm_driv_mod.F%EMITSO2', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ALLOCATED( DRYDEP ) ) THEN
       DEALLOCATE( DRYDEP, STAT=RC )
       CALL GC_CheckVar( 'apm_driv_mod.F%DRYDEP', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ALLOCATED( WETDEP ) ) THEN
       DEALLOCATE( WETDEP, STAT=RC )
       CALL GC_CheckVar( 'apm_driv_mod.F%WETDEP', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ALLOCATED( CONDEP ) ) THEN
       DEALLOCATE( CONDEP, STAT=RC )
       CALL GC_CheckVar( 'apm_driv_mod.F%CONDEP', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ALLOCATED( UPTAKE ) ) THEN
       DEALLOCATE( UPTAKE, STAT=RC )
       CALL GC_CheckVar( 'apm_driv_mod.F%UPTAKE', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ALLOCATED( OXIDAT ) )  THEN
       DEALLOCATE( OXIDAT, STAT=RC )
       CALL GC_CheckVar( 'apm_driv_mod.F%OXIDAT', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ALLOCATED( FCLOUD ) ) THEN
       DEALLOCATE( FCLOUD, STAT=RC )
       CALL GC_CheckVar( 'apm_driv_mod.F%FCLOUD', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ALLOCATED( AERAREA ) ) THEN
       DEALLOCATE( AERAREA, STAT=RC )
       CALL GC_CheckVar( 'apm_driv_mod.F%AERAREA', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ALLOCATED( AERDRYR) ) THEN
       DEALLOCATE( AERDRYR, STAT=RC )
       CALL GC_CheckVar( 'apm_driv_mod.F%AERDRYR', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ALLOCATED( GAMMAPM ) ) THEN
       DEALLOCATE( GAMMAPM, STAT=RC )
       CALL GC_CheckVar( 'apm_driv_mod.F%GAMMAPM', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ALLOCATED( XN4D ) ) THEN
       DEALLOCATE( XN4D, STAT=RC )
       CALL GC_CheckVar( 'apm_driv_mod.F%XN4D', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ALLOCATED( XQ3D ) ) THEN
       DEALLOCATE( XQ3D, STAT=RC )
       CALL GC_CheckVar( 'apm_driv_mod.F%XQ3D', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ALLOCATED( RCLDL3D) ) THEN
       DEALLOCATE( RCLDL3D, STAT=RC )
       CALL GC_CheckVar( 'apm_driv_mod.F%RCLDL3D', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ALLOCATED( CDN3D) ) THEN
       DEALLOCATE( CDN3D, STAT=RC )
       CALL GC_CheckVar( 'apm_driv_mod.F%CDN3D', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ALLOCATED( RCLDI3D ) ) THEN
       DEALLOCATE( RCLDI3D, STAT=RC )
       CALL GC_CheckVar( 'apm_driv_mod.F%RCLDI3D', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ALLOCATED( RCLDI2D ) ) THEN
       DEALLOCATE( RCLDI2D, STAT=RC )
       CALL GC_CheckVar( 'apm_driv_mod.F%RCLDL2D', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ALLOCATED( IN3D ) ) THEN
       DEALLOCATE( IN3D, STAT=RC )
       CALL GC_CheckVar( 'apm_driv_mod.F%CDN3D', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ALLOCATED( MASSISRP ) ) THEN
       DEALLOCATE( MASSISRP, STAT=RC )
       CALL GC_CheckVar( 'apm_driv_mod.F%MASSISRP', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ALLOCATED( MASSMESA ) ) THEN
       DEALLOCATE( MASSMESA, STAT=RC )
       CALL GC_CheckVar( 'apm_driv_mod.F%MASSMESA', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ALLOCATED( GFTOT3D ) ) THEN
       DEALLOCATE( GFTOT3D, STAT=RC )
       CALL GC_CheckVar( 'apm_driv_mod.F%GFTOT3D', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ALLOCATED( DENWET3D ) ) THEN
       DEALLOCATE( DENWET3D, STAT=RC )
       CALL GC_CheckVar( 'apm_driv_mod.F%DENWET3D', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ALLOCATED( MWSIZE3D ) ) THEN
       DEALLOCATE( MWSIZE3D, STAT=RC )
       CALL GC_CheckVar( 'apm_driv_mod.F%MWSIZE3D', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ALLOCATED( SO2toSO4 ) ) THEN
       DEALLOCATE( SO2toSO4, STAT=RC )
       CALL GC_CheckVar( 'apm_driv_mod.F%SO2toSO4', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ALLOCATED( NCOAG3D ) ) THEN
       DEALLOCATE( NCOAG3D, STAT=RC )
       CALL GC_CheckVar( 'apm_driv_mod.F%NCOAG3D', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ALLOCATED( IACT1 ) ) THEN
       DEALLOCATE( IACT1, STAT=RC )
       CALL GC_CheckVar( 'apm_driv_mod.F%IACT1', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ALLOCATED( IACT2 ) ) THEN
       DEALLOCATE( IACT2, STAT=RC )
       CALL GC_CheckVar( 'apm_driv_mod.F%IACT2', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ALLOCATED( IACT3 ) ) THEN
       DEALLOCATE( IACT3, STAT=RC )
       CALL GC_CheckVar( 'apm_driv_mod.F%IACT3', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ALLOCATED( PSO4GAS) ) THEN
       DEALLOCATE( PSO4GAS, STAT=RC )
       CALL GC_CheckVar( 'apm_driv_mod.F%PSO4GAS', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ALLOCATED( XO3 ) ) THEN
       DEALLOCATE( XO3, STAT=RC )
       CALL GC_CheckVar( 'apm_driv_mod.F%XO3', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ALLOCATED( PLVSOG) ) THEN
       DEALLOCATE( PLVSOG, STAT=RC )
       CALL GC_CheckVar( 'apm_driv_mod.F%PLVSOG', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ALLOCATED( BCLIFE ) ) THEN
       DEALLOCATE( BCLIFE, STAT=RC )
       CALL GC_CheckVar( 'apm_driv_mod.F%BCLIFE', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ALLOCATED( OCLIFE ) ) THEN
       DEALLOCATE( OCLIFE, STAT=RC )
       CALL GC_CheckVar( 'apm_driv_mod.F%OCLIFE', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ALLOCATED( IFOUTIJ ) ) THEN
       DEALLOCATE( IFOUTIJ, STAT=Rc )
       CALL GC_CheckVar( 'apm_driv_mod.F%IFOUTIJ', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ALLOCATED( SITEID ) ) THEN
       DEALLOCATE( SITEID, STAT=RC )
       CALL GC_CheckVar( 'apm_driv_mod.F%SITEID', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ALLOCATED( FSST ) ) THEN
       DEALLOCATE( FSST, STAT=RC )
       CALL GC_CheckVar( 'apm_driv_mod.F%FSST', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ALLOCATED( SPGF ) ) THEN
       DEALLOCATE( SPGF, STAT=RC )
       CALL GC_CheckVar( 'apm_driv_mod.F%SPGF', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ALLOCATED( TAONH3 ) ) THEN
       DEALLOCATE( TAONH3, STAT=RC )
       CALL GC_CheckVar( 'apm_driv_mod.F%TAONH3', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ALLOCATED( YCOD3D ) ) THEN
       DEALLOCATE( YCOD3D, STAT=RC )
       CALL GC_CheckVar( 'apm_driv_mod.F%YCOD3D', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ALLOCATED( TCOD3D ) ) THEN
       DEALLOCATE( TCOD3D, STAT=RC )
       CALL GC_CheckVar( 'apm_driv_mod.F%TCOD3D', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ALLOCATED( YHTRC3D ) ) THEN
       DEALLOCATE( YHTRC3D, STAT=RC )
       CALL GC_CheckVar( 'apm_driv_mod.F%YHTRC3D', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ALLOCATED( YHTR3D ) ) THEN
       DEALLOCATE( YHTR3D, STAT=RC )
       CALL GC_CheckVar( 'apm_driv_mod.F%YHTR3D', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ALLOCATED( YHTRC03D ) ) THEN
       DEALLOCATE( YHTRC03D, STAT=RC )
       CALL GC_CheckVar( 'apm_driv_mod.F%YHTRC03D', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ALLOCATED( YHTR03D ) ) THEN
       DEALLOCATE( YHTR03D, STAT=RC )
       CALL GC_CheckVar( 'apm_driv_mod.F%YHTR03D', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ALLOCATED( ZBEXT3D ) ) THEN
       DEALLOCATE( ZBEXT3D, STAT=RC )
       CALL GC_CheckVar( 'apm_driv_mod.F%ZBEXT3D', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ALLOCATED( ZW3D ) ) THEN
       DEALLOCATE( ZW3D, STAT=RC )
       CALL GC_CheckVar( 'apm_driv_mod.F%ZW3D', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ALLOCATED( ZG3D ) ) THEN
       DEALLOCATE( ZG3D, STAT=RC )
       CALL GC_CheckVar( 'apm_driv_mod.F%ZG3D', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ALLOCATED( YBEXT3D ) ) THEN
       DEALLOCATE( YBEXT3D, STAT=RC )
       CALL GC_CheckVar( 'apm_driv_mod.F%YBEXT3D', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ALLOCATED( XBEXT1k3D ) ) THEN
       DEALLOCATE( XBEXT1k3D, STAT=RC )
       CALL GC_CheckVar( 'apm_driv_mod.F%XBEXT1K3D', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ALLOCATED( YW3D ) ) THEN
       DEALLOCATE( YW3D, STAT=RC )
       CALL GC_CheckVar( 'apm_driv_mod.F%YW3D', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ALLOCATED( YG3D ) ) THEN
       DEALLOCATE( YG3D, STAT=RC )
       CALL GC_CheckVar( 'apm_driv_mod.F%YG3D', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ALLOCATED( ZBABS3D ) ) THEN
       DEALLOCATE( ZBABS3D, STAT=RC )
       CALL GC_CheckVar( 'apm_driv_mod.F%ZBABS3D', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ALLOCATED( YUV ) ) THEN
       DEALLOCATE( YUV, STAT=RC )
       CALL GC_CheckVar( 'apm_driv_mod.F%YUV', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ALLOCATED( GOOD ) ) THEN
       DEALLOCATE( GOOD, STAT=RC )
       CALL GC_CheckVar( 'apm_driv_mod.F%GOOD', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ALLOCATED( CODOUT ) )  THEN
       DEALLOCATE( CODOUT, STAT=RC )
       CALL GC_CheckVar( 'apm_driv_mod.F%CODOUT', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ALLOCATED( CODOUTNUM ) ) THEN
       DEALLOCATE( CODOUTNUM, STAT=RC )
       CALL GC_CheckVar( 'apm_driv_mod.F%CODOUTNUM', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ALLOCATED( MODISOUT) ) THEN
       DEALLOCATE( MODISOUT, STAT=RC )
       CALL GC_CheckVar( 'apm_driv_mod.F%MODISOUT', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ALLOCATED( TCOSZ) ) THEN
       DEALLOCATE( TCOSZ, STAT=RC )
       CALL GC_CheckVar( 'apm_driv_mod.F%TCOSZ', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ALLOCATED( TTDAY) ) THEN
       DEALLOCATE( TTDAY, STAT=RC )
       CALL GC_CheckVar( 'apm_driv_mod.F%TTDAY', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ALLOCATED( COSZM) ) THEN
       DEALLOCATE( COSZM, STAT=RC )
       CALL GC_CheckVar( 'apm_driv_mod.F%COSZM', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ALLOCATED( ALBDRR ) ) THEN
       DEALLOCATE( ALBDRR, STAT=RC )
       CALL GC_CheckVar( 'apm_driv_mod.F%ALBDRR', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ALLOCATED( ZTCST ) ) THEN
       DEALLOCATE( ZTCST, STAT=RC )
       CALL GC_CheckVar( 'apm_driv_mod.F%ZTCST', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ALLOCATED( ZTCSB ) ) THEN
       DEALLOCATE( ZTCSB, STAT=RC )
       CALL GC_CheckVar( 'apm_driv_mod.F%ZTCSB', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ALLOCATED( ZTFST ) ) THEN
       DEALLOCATE( ZTFST, STAT=RC )
       CALL GC_CheckVar( 'apm_driv_mod.F%ZTFST', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ALLOCATED( ZTFSB ) ) THEN
       DEALLOCATE( ZTFSB, STAT=RC )
       CALL GC_CheckVar( 'apm_driv_mod.F%YG3D', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ALLOCATED( ZTFSA ) ) THEN
       DEALLOCATE( ZTFSA, STAT=RC )
       CALL GC_CheckVar( 'apm_driv_mod.F%ZTFSA', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ALLOCATED( ZTCLT ) ) THEN
       DEALLOCATE( ZTCLT, STAT=RC )
       CALL GC_CheckVar( 'apm_driv_mod.F%ZTCLT', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ALLOCATED( ZTCLB ) ) THEN
       DEALLOCATE( ZTCLB, STAT=RC )
       CALL GC_CheckVar( 'apm_driv_mod.F%ZTCLB', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ALLOCATED( ZTFUL ) ) THEN
       DEALLOCATE( ZTFUL, STAT=RC )
       CALL GC_CheckVar( 'apm_driv_mod.F%ZTFUL', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ALLOCATED( ZTFDL  ) ) THEN
       DEALLOCATE( ZTFDL, STAT=RC )
       CALL GC_CheckVar( 'apm_driv_mod.F%ZTFDL', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ALLOCATED( ZTFLA ) ) THEN
       DEALLOCATE( ZTFLA, STAT=RC )
       CALL GC_CheckVar( 'apm_driv_mod.F%ZTFLA', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ALLOCATED( ZTAOD ) ) THEN
       DEALLOCATE( ZTAOD, STAT=RC )
       CALL GC_CheckVar( 'apm_driv_mod.F%ZTAOD', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ALLOCATED( ZCLDF ) ) THEN
       DEALLOCATE( ZCLDF, STAT=RC )
       CALL GC_CheckVar( 'apm_driv_mod.F%ZCLDF', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ALLOCATED( ZALB ) ) THEN
       DEALLOCATE( ZALB, STAT=RC )
       CALL GC_CheckVar( 'apm_driv_mod.F%ZALB', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ALLOCATED( ZMALB ) ) THEN
       DEALLOCATE( ZMALB, STAT=RC )
       CALL GC_CheckVar( 'apm_driv_mod.F%ZMALB', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ALLOCATED( ZCST ) ) THEN
       DEALLOCATE( ZCST, STAT=RC )
       CALL GC_CheckVar( 'apm_driv_mod.F%ZCST', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ALLOCATED( ZCSB ) ) THEN
       DEALLOCATE( ZCSB, STAT=RC )
       CALL GC_CheckVar( 'apm_driv_mod.F%ZCSB', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ALLOCATED( ZFST ) ) THEN
       DEALLOCATE( ZFST, STAT=RC )
       CALL GC_CheckVar( 'apm_driv_mod.F%ZFST', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ALLOCATED( ZCLD ) ) THEN
       DEALLOCATE( ZCLD, STAT=RC )
       CALL GC_CheckVar( 'apm_driv_mod.F%ZCLD', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ALLOCATED( ZCLD0) ) THEN
       DEALLOCATE( ZCLD0, STAT=RC )
       CALL GC_CheckVar( 'apm_driv_mod.F%ZCLD0', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ALLOCATED( ZFSB ) ) THEN
       DEALLOCATE( ZFSB, STAT=RC )
       CALL GC_CheckVar( 'apm_driv_mod.F%ZFSB', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ALLOCATED( ZFSA ) ) THEN
       DEALLOCATE( ZFSA, STAT=RC )
       CALL GC_CheckVar( 'apm_driv_mod.F%ZFSA', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ALLOCATED( ZCLT ) ) THEN
       DEALLOCATE( ZCLT, STAT=RC )
       CALL GC_CheckVar( 'apm_driv_mod.F%', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ALLOCATED( ZCLB) ) THEN
       DEALLOCATE( ZCLB, STAT=RC )
       CALL GC_CheckVar( 'apm_driv_mod.F%', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ALLOCATED( ZFUL ) ) THEN
       DEALLOCATE( ZFUL, STAT=RC )
       CALL GC_CheckVar( 'apm_driv_mod.F%', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ALLOCATED( ZFDL ) ) THEN
       DEALLOCATE( ZFDL, STAT=RC )
       CALL GC_CheckVar( 'apm_driv_mod.F%', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ALLOCATED( ZFLA ) ) THEN
       DEALLOCATE( ZFLA, STAT=RC )
       CALL GC_CheckVar( 'apm_driv_mod.F%', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ALLOCATED( ZAOD) ) THEN
       DEALLOCATE( ZAOD, STAT=RC )
       CALL GC_CheckVar( 'apm_driv_mod.F%ZAOD', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ALLOCATED( ZCOD ) ) THEN
       DEALLOCATE( ZCOD, STAT=RC )
       CALL GC_CheckVar( 'apm_driv_mod.F%ZCOD', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ALLOCATED( ZCODGC ) ) THEN
       DEALLOCATE( ZCODGC, STAT=RC )
       CALL GC_CheckVar( 'apm_driv_mod.F%ZCODGC', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ALLOCATED( ZAODOUT1 ) ) THEN
       DEALLOCATE( ZAODOUT1, STAT=RC )
       CALL GC_CheckVar( 'apm_driv_mod.F%ZADOUT1', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ALLOCATED( ZAODOUT3 ) ) THEN
       DEALLOCATE( ZAODOUT3, STAT=RC )
       CALL GC_CheckVar( 'apm_driv_mod.F%ZAODOUT3', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ALLOCATED( ZAOD25 ) ) THEN
       DEALLOCATE( ZAOD25, STAT=RC )
       CALL GC_CheckVar( 'apm_driv_mod.F%ZAOD25', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ALLOCATED( ZAOD50 ) ) THEN
       DEALLOCATE( ZAOD50, STAT=RC )
       CALL GC_CheckVar( 'apm_driv_mod.F%ZAOD50', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ALLOCATED( ZAOD75 ) ) THEN
       DEALLOCATE( ZAOD75, STAT=RC )
       CALL GC_CheckVar( 'apm_driv_mod.F%ZAOD75', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ALLOCATED( NAOD25 ) ) THEN
       DEALLOCATE( NAOD25, STAT=RC )
       CALL GC_CheckVar( 'apm_driv_mod.F%NAOD25', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ALLOCATED( NAOD50 ) ) THEN
       DEALLOCATE( NAOD50, STAT=RC )
       CALL GC_CheckVar( 'apm_driv_mod.F%NAOD50', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ALLOCATED( NAOD75 ) ) THEN
       DEALLOCATE( NAOD75, STAT=RC )
       CALL GC_CheckVar( 'apm_driv_mod.F%NAOD75', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ALLOCATED( ZVIS) ) THEN
       DEALLOCATE( ZVIS, STAT=RC )
       CALL GC_CheckVar( 'apm_driv_mod.F%ZVIS', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ALLOCATED( THAZ ) ) THEN
       DEALLOCATE( THAZ, STAT=RC )
       CALL GC_CheckVar( 'apm_driv_mod.F%THAZ', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ALLOCATED( TFOG ) ) THEN
       DEALLOCATE( TFOG, STAT=RC )
       CALL GC_CheckVar( 'apm_driv_mod.F%TFOG', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ALLOCATED( ZAAOD ) ) THEN
       DEALLOCATE( ZAAOD, STAT=RC )
       CALL GC_CheckVar( 'apm_driv_mod.F%ZAAOD', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ALLOCATED( ZAAODOUT1 ) ) THEN
       DEALLOCATE( ZAAODOUT1, STAT=RC )
       CALL GC_CheckVar( 'apm_driv_mod.F%ZAAODOUT1', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ALLOCATED( ZAAODOUT3 ) ) THEN
       DEALLOCATE( ZAAODOUT3, STAT=RC )
       CALL GC_CheckVar( 'apm_driv_mod.F%ZAAODOUT3', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ALLOCATED( ZABS ) ) THEN
       DEALLOCATE( ZABS, STAT=RC )
       CALL GC_CheckVar( 'apm_driv_mod.F%ZABS', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ALLOCATED( ZWCL ) ) THEN
       DEALLOCATE( ZWCL, STAT=RC )
       CALL GC_CheckVar( 'apm_driv_mod.F%ZWCL', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ALLOCATED( ZWCI ) ) THEN
       DEALLOCATE( ZWCI, STAT=RC )
       CALL GC_CheckVar( 'apm_driv_mod.F%ZWCI', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ALLOCATED( AOD ) ) THEN
       DEALLOCATE( AOD, STAT=RC )
       CALL GC_CheckVar( 'apm_driv_mod.F%AOD', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ALLOCATED( AAOD ) ) THEN
       DEALLOCATE( AAOD, STAT=RC )
       CALL GC_CheckVar( 'apm_driv_mod.F%AAOD', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ALLOCATED( AODOUT1 ) ) THEN
       DEALLOCATE( AODOUT1, STAT=RC )
       CALL GC_CheckVar( 'apm_driv_mod.F%AODOUT1', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ALLOCATED( AAODOUT1 ) ) THEN
       DEALLOCATE( AAODOUT1, STAT=RC )
       CALL GC_CheckVar( 'apm_driv_mod.F%AAODOUT1', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ALLOCATED( AODOUT3 ) ) THEN
       DEALLOCATE( AODOUT3, STAT=RC )
       CALL GC_CheckVar( 'apm_driv_mod.F%', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ALLOCATED( AAODOUT3 ) ) THEN
       DEALLOCATE( AAODOUT3, STAT=RC )
       CALL GC_CheckVar( 'apm_driv_mod.F%AAODOUT3', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ALLOCATED( CODGC ) ) THEN
       DEALLOCATE( CODGC, STAT=RC )
       CALL GC_CheckVar( 'apm_driv_mod.F%CODGC', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ALLOCATED( TAOD) ) THEN
       DEALLOCATE( TAOD, STAT=RC )
       CALL GC_CheckVar( 'apm_driv_mod.F%TAOD', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

  END SUBROUTINE CLEANUP_APM3D
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: aeronum
!
! !DESCRIPTION: Subroutine SIZE\_DIST\_PARAM\_LIQ gets cloud droplet size
!  distribution parameters
!\\
!\\
! !INTERFACE:
!
  subroutine size_dist_param_liq(qcic, ncic, cdnl, rho, nadjflag, pgam, lamc)
!
! !INPUT/OUTPUT PARAMETERS:
!
    real*8, intent(in)    :: qcic
    real*8, intent(inout) :: ncic
    real*8, intent(in)    :: cdnl
    real*8, intent(in)    :: rho
    logical, intent(in)   :: nadjflag ! Whether to adjust number concentration to fall
                                      ! within certain bounds

    real*8, intent(out)   :: pgam
    real*8, intent(out)   :: lamc
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    real, parameter     :: pi=3.1415926, rhow=1000.
    real*8 :: dumgam1
    real*8 :: dumgam2
    real*8 :: lammin
    real*8 :: lammax

    if (qcic > 1.e-18) then

       if (nadjflag) then
          ! add upper limit to in-cloud number concentration to prevent numerical error
          ncic=min(ncic,qcic*1.e20)
          ! add lower limit to in-cloud number concentration
          ncic=max(ncic,cdnl/rho) ! sghan minimum in #/cm
       end if

       ! get pgam from fit to observations of martin et al. 1994
       pgam=0.0005714*(ncic/1.e6*rho)+0.2714
       pgam=1./(pgam**2)-1.
       pgam=max(pgam,2.d0)
       pgam=min(pgam,50.d0)

       ! calculate lamc
       dumgam1 = gamma(pgam+1.)
       dumgam2 = gamma(pgam+4.)

       lamc = (pi/6.*rhow*ncic*dumgam2/ &
              (qcic*dumgam1))**(1./3.)

       ! lammin, 50 micron diameter max mean size

       lammin = (pgam+1.)/50.e-6
       lammax = (pgam+1.)/2.e-6

       if (lamc < lammin) then
          lamc = lammin

          if (nadjflag) then
             ncic = 6. * lamc**3 * qcic * dumgam1/ &
                    (pi * rhow * dumgam2)
          end if
       else if (lamc > lammax) then
          lamc = lammax

          if (nadjflag) then
             ncic = 6. * lamc**3 * qcic * dumgam1/ &
                    (pi * rhow * dumgam2)
          end if
       end if

    else
       ! pgam not calculated in this case, so set it to a value likely to cause an error
       ! if it is accidentally used
       ! (gamma function undefined for negative integers)
       pgam = -100.
       lamc = 0.
    end if

  end subroutine size_dist_param_liq
!EOC
! ----------------------------------------------------------------------------
! function GAMMA
! ----------------------------------------------------------------------------
!BOC
  function gamma(x)
    implicit none
!
! Purpose:
!   Returns the gamma function
!
! Input:
!   [x]   value to compute gamma function of
!
! Returns:
!   gamma(x)
!
! Coded:
!   02/02/06  John Haynes (haynes@atmos.colostate.edu)
!   (original code of unknown origin)

! ----- INPUTS -----
    real*8, intent(in) :: x

! ----- OUTPUTS -----
    real*8 :: gamma

! ----- INTERNAL -----
    real*8 :: pi,ga,z,r,gr
    real*8 :: g(26)
    integer :: k,m1,m

    pi = acos(-1.)
    if (x ==int(x)) then
       if (x > 0.0) then
          ga=1.0
          m1=x-1
          do k=2,m1
             ga=ga*k
          enddo
       else
          ga=1.0+300
       endif
    else
       if (abs(x) > 1.0) then
          z=abs(x)
          m=int(z)
          r=1.0
          do k=1,m
             r=r*(z-k)
          enddo
          z=z-m
       else
          z=x
       endif
       data g/1.0,0.5772156649015329, &
              -0.6558780715202538, -0.420026350340952d-1, &
              0.1665386113822915,-.421977345555443d-1, &
              -.96219715278770d-2, .72189432466630d-2, &
              -.11651675918591d-2, -.2152416741149d-3, &
              .1280502823882d-3, -.201348547807d-4, &
              -.12504934821d-5, .11330272320d-5, &
              -.2056338417d-6, .61160950d-8, &
              .50020075d-8, -.11812746d-8, &
              .1043427d-9, .77823d-11, &
              -.36968d-11, .51d-12, &
              -.206d-13, -.54d-14, .14d-14, .1d-15/
       gr=g(26)
       do k=25,1,-1
          gr=gr*z+g(k)
       enddo
       ga=1.0/(gr*z)
       if (abs(x) > 1.0) then
          ga=ga*r
          if (x < 0.0) ga=-pi/(x*ga*sin(pi*x))
       endif
    endif
    gamma = ga
    return

  end function gamma
!EOC
END MODULE APM_DRIV_MOD
#endif
