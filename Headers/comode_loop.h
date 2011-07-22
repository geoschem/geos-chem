!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !INCLUDE: comode.h
!
! !DESCRIPTION: Header file COMODE contains common blocks and variables 
!  for the SMVGEAR II chemistry package.
!\\
!\\
! !REMARKS:
! *********************************************************************
! ************        WRITTEN BY MARK JACOBSON (1993)      ************
! ***             (C) COPYRIGHT, 1993 BY MARK Z. JACOBSON           *** 
! ***       U.S. COPYRIGHT OFFICE REGISTRATION NO. TXu 670-279      *** 
! ***                         (650) 723-6836                        *** 
! *********************************************************************
!
!         CCCCCCC  OOOOOOO  M     M  OOOOOOO  DDDDDD   EEEEEEE 
!         C        O     O  M M M M  O     O  D     D  E  
!         C        O     O  M  M  M  O     O  D     D  EEEEEEE  
!         C        O     O  M     M  O     O  D     D  E  
!         CCCCCCC  OOOOOOO  M     M  OOOOOOO  DDDDDD   EEEEEEE 
!
! *********************************************************************
! * THIS IS THE COMMON BLOCK FOR "SMVGEAR" AND "MIE," TWO ORDINARY    *
! * DIFFERENTIAL EQUATION SOLVERS. THE REFERENCE FOR THE CODES IS     *
! *                                                                   * 
! *   JACOBSON M. Z. AND TURCO R. P. (1993) SMVGEAR: A SPARSE-        * 
! *    MATRIX, VECTORIZED GEAR CODE FOR ATMOSPHERIC MODELS.           *
! *    SUBMITTED TO ATMOSPHERIC ENVIRONMENT, PART A. MAY 20, 1993     * 
! *                                                                   * 
! * COMODE.H SETS PARAMETER VALUES AND SERVES AS A COMMON BLOCK FOR   *
! * ALL DIMENSIONED AND NON-DIMENSIONED VARIABLES. COMODE.H ALSO      *
! * DEFINES EACH PARAMETER, BUT DATA FILE DEFINE.DAT EXPLAINS NON-    *
! * DIMENSIONED VARIABLES. INDIVIDUAL SUBROUTINES DEFINE DIMENSIONED  *
! * VARIABLES.                                                        *
! *********************************************************************
!
! !REVISION HISTORY:
!  22 Jul 2011 - R. Yantosca - Initial version, based on column code version
!
! !DEFINED PARAMETERS:
!
!      !---------------------------------------------------------------
!      ! Parameters for backwards compatibility with the older
!      ! "comode.h" file.  These can be removed later.
!      !---------------------------------------------------------------
!
!      INTEGER, PARAMETER :: ILAT    = JJPAR
!      INTEGER, PARAMETER :: ILONG   = IIPAR
!      INTEGER, PARAMETER :: IVERT   = LLTROP
!      INTEGER, PARAMETER :: IPLUME  = 0
!      INTEGER, PARAMETER :: IPVERT  = IVERT + IPLUME
!      INTEGER, PARAMETER :: ITLOOP  = ILAT * ILONG * IPVERT
!      INTEGER, PARAMETER :: IMLOOP  = ILAT * ILONG
!      INTEGER, PARAMETER :: ILAYER  = IVERT + 1
!      INTEGER, PARAMETER :: ILTLOOP = IMLOOP * ILAYER
!      INTEGER, PARAMETER :: MAXDAYS = 1000
!      INTEGER, PARAMETER :: MXBLOCK = 16 + ITLOOP/KBLOOP      
!
!      ! NEMPARA = max no. of anthropogenic emissions
!      ! NEMPARB = max no. of biogenic emissions
!      INTEGER IDEMS
!      COMMON /JTRCID/ IDEMS(NEMPARA+NEMPARB)

      !---------------------------------------------------------------
      ! Physical constants 
      !---------------------------------------------------------------

      ! Avogadro's #
      REAL*8,  PARAMETER :: AVG       = 6.02252d+23

      ! Boltzmann's constant [erg/K]
      REAL*8,  PARAMETER :: BK        = 1.38054d-16
      REAL*8,  PARAMETER :: BOLTG     = 1.38054d-16

      ! Condensation vapor pressure ?????? 
      REAL*8,  PARAMETER :: CONSVAP   = 6.1078d+03 / BOLTG
      
      ! PI (same value as in CMN_GCTM) and related quantities
      REAL*8,  PARAMETER :: ONEPI     = 3.14159265358979323d0 
      REAL*8,  PARAMETER :: EIGHTDPI  = 8.d0 / ONEPI

      ! Gas constant [erg/K/mole]
      REAL*8,  PARAMETER :: RGAS      = BOLTG * AVG
      
      ! Universal gas constant [g/cm2/s2/mole/K]
      REAL*8,  PARAMETER :: RSTARG    = 8.31450d+07

      ! Seconds per day
      REAL*8,  PARAMETER :: SCDAY     = 86400.0d0

      ! Molecular weight of air 
      REAL*8,  PARAMETER :: WTAIR     = 28.966d0

      !---------------------------------------------------------------
      ! Flags etc. mostly from "reader.f"
      !---------------------------------------------------------------

      ! Call chemistry solver
      INTEGER, PARAMETER :: IFSOLVE  = 1

      ! Chemistry mechanism goes in urban slot
      ! (NCSGAS is max # of slots that are used, set to 1)
      INTEGER, PARAMETER :: IFURBAN  = 1
      INTEGER, PARAMETER :: NCSURBAN = 1

      ! Nothing goes in trop slot
      INTEGER, PARAMETER :: IFTROP   = 0
      INTEGER, PARAMETER :: NCSTROP  = 0

      ! Nothing goes in strat slot
      INTEGER, PARAMETER :: IFSTRAT  = 0
      INTEGER, PARAMETER :: NCSSTRAT = 0

      ! Number of chemistry slots (set to NCSURBAN=1)
      INTEGER, PARAMETER :: NCSGAS   = NCSURBAN

      ! Max number of chemistry slots (set to NCSGAS=1)
      INTEGER, PARAMETER :: NCSALL   = NCSGAS
      INTEGER, PARAMETER :: NCSTRST  = NCSGAS

!------------------------------------------------------------------------------
! Prior to 7/22/11:
! For now, restore the NCS common block until we can rewrite the source code 
! to remove this from DO loops (bmy, 7/22/11)
!      ! Counter for the chemistry slot (NOTE: this was a DO loop variable
!      ! in a common block, very bad coding style.  We'll just make this a
!      ! parameter and rewrite the code to remove the do loops.
!      INTEGER, PARAMETER :: NCS      = NCSURBAN
!------------------------------------------------------------------------------
      INTEGER            :: NCS
      COMMON /TMPNCS/       NCS

      ! Print species and reactions in readchem.f
      INTEGER, PARAMETER :: IOSPEC   = 1
      INTEGER, PARAMETER :: IOREAC   = 1

      ! Print out info from "reader.f"
      INTEGER, PARAMETER :: IPREADER = 1

      ! We do not reorder boxes by stiffness
      INTEGER, PARAMETER :: ISREORD  = 0
     
      !---------------------------------------------------------------
      ! Small number tolerances
      !---------------------------------------------------------------
      
      REAL*8,  PARAMETER :: SMAL1    = 1d-06 
      REAL*8,  PARAMETER :: SMAL2    = 1.0d-99
      REAL*8,  PARAMETER :: SMAL3    = 1d-50

      !---------------------------------------------------------------
      ! Placeholder parameters to size some common blocks
      !---------------------------------------------------------------
!------------------------------------------------------------------------------
! Prior to 7/22/11:
! Use same value (KBLOOP=24) as in the current comode.h (bmy, 7/22/11)
!      ! NOTE: KBLOOP has to be dimensioned as large as the max # 
!      ! of vertical levels in a column.  The parameter MAX_COLUMN
!      ! is contained w/in include file "smv_dimension.h"
!      INTEGER, PARAMETER :: KBLOOP   = MAX_COLUMN
!------------------------------------------------------------------------------
      ! Number of individual grid boxes per grid block.  We pass this many 
      ! grid boxes to each CPU.  The convergence criterion is computed
      ! for the whole grid block (all KBLOOP grid boxes). (bmy, 7/22/11)
      INTEGER, PARAMETER :: KBLOOP   = 24

      !---------------------------------------------------------------
      ! Gas-phase parameters
      !---------------------------------------------------------------

      ! Maximum number of gases, active + inactive 
!------------------------------------------------------------------------------
! Prior to 7/22/11:
! Use same value (IGAS=300) as in the current comode.h (bmy, 7/22/11)
!      INTEGER, PARAMETER :: IGAS     = MAX_SPECIES
!------------------------------------------------------------------------------
      INTEGER, PARAMETER :: IGAS     = 300

      ! Maximum number of aqueous chemistry species 
      ! (set = 1 here since no aqueous chemistry is included)
      INTEGER, PARAMETER :: IAERTY   = 1  

      ! Maximum number of rates constants (max # reactions)
!------------------------------------------------------------------------------
! Prior to 7/22/11:
! Use same value (NMRATE=700) as in the current comode.h (bmy, 7/22/11)
!      INTEGER, PARAMETER :: NMRATE   = 700
!------------------------------------------------------------------------------
      INTEGER, PARAMETER :: NMRATE   = 700

      ! Maximum number of photo rates
!------------------------------------------------------------------------------
! Prior to 7/22/11:
! Use same value (IPHOT=85) as in the current comode.h (bmy, 7/22/11)
!      INTEGER, PARAMETER :: IPHOT    = 60
!------------------------------------------------------------------------------
      INTEGER, PARAMETER :: IPHOT    = 85

      ! Maximum number of kinetic + photo reactions
      INTEGER, PARAMETER :: NMTRATE  = NMRATE + IPHOT

      ! Maximum number of aqueous chemical reactions plus photo
      ! processes (set = 1 here since no aqueous chemistry included)
      INTEGER, PARAMETER :: NMQRATE  = 1           

      ! Maximum number of species in a reaction rate 
      INTEGER, PARAMETER :: NMRPROD  = 25

      ! Maximum number of dead species
      INTEGER, PARAMETER :: NMDEAD   = 100         

      ! Maximum number of gains / losses for each chemical species
!------------------------------------------------------------------------------
! Prior to 7/22/11:
! Use same value (MAXGL=1000) as in the current comode.h (bmy, 7/22/11)
!      INTEGER, PARAMETER :: MAXGL    = 430
!------------------------------------------------------------------------------
      INTEGER, PARAMETER :: MAXGL    = 1000

      ! An array dimension smaller than MAXGL
!------------------------------------------------------------------------------
! Prior to 7/22/11:
! Use same value (MAXGL2=90) as in the current comode.h (bmy, 7/22/11)
!      INTEGER, PARAMETER :: MAXGL2   = 50          
!------------------------------------------------------------------------------
      INTEGER, PARAMETER :: MAXGL2   = 90          

      ! An array dimension smaller than MAXGL2
!------------------------------------------------------------------------------
! Prior to 7/22/11:
! Use same value (MAXGL3=NNPAR) as in the current comode.h (bmy, 7/22/11)
!      ! NOTE: MAX_TRACERS is contained in include file "smv_dimension.h"
      INTEGER, PARAMETER :: MAXGL3   = MAX_TRACERS
!------------------------------------------------------------------------------
!      INTEGER, PARAMETER :: MAXGL3   = NNPAR

      ! An array dimension smaller than MAXGL2
      INTEGER, PARAMETER :: MAXGL4   = 10   

      ! Number OF TYPES OF CHEMISTRY (set to NCSGAS=1)
!------------------------------------------------------------------------------
! Prior to 7/22/11:
! Use same value (ICS=3) as in the current comode.h (bmy, 7/22/11)
!      INTEGER, PARAMETER :: ICS      = NCSGAS
!------------------------------------------------------------------------------
      INTEGER, PARAMETER :: ICS      = 3

      ! Types of chemistry * 2 (one for day, one for night)
      INTEGER, PARAMETER :: ICP      = 2*ICS       

      ! Maximum order for gear parameters for dimension purposes
      INTEGER, PARAMETER :: MORDER   = 7 

      ! Various other parameters
      INTEGER, PARAMETER :: IPHOT8   = IPHOT + 8      

      ! Dimension for XINP array
      INTEGER, PARAMETER :: IMISC    = 100 

      ! Dimension for MBCOMP, WTMB, JMBCOMP arrays
      INTEGER, PARAMETER :: IMASBAL  = 9          
      
      ! Dimension for ARRT, BRRT, KCRRT, etc. arrays
      INTEGER, PARAMETER :: MXCOF    = 5 

      !---------------------------------------------------------------
      ! Gas-phase parameters
      !---------------------------------------------------------------

      ! Larger of IGAS, IAERTY 
      INTEGER, PARAMETER :: MXGSAER  = IGAS
      
      ! Larger of NMTRATE, NMQRATE 
      INTEGER, PARAMETER :: MXRATE   = NMTRATE

      ! Larger of MXGSAER, MXARRAY
      INTEGER, PARAMETER :: MXARRAY  = 3000

      ! Maximum one-dimensional array-length of sparse matrix
      INTEGER, PARAMETER :: MXCC2    = MXARRAY

      ! Various other array dimensions
      INTEGER, PARAMETER :: MXCOUNT1 = MXGSAER * MAXGL3 * 3
      INTEGER, PARAMETER :: MXCOUNT2 = MXGSAER * MAXGL3 * 7  
      INTEGER, PARAMETER :: MXCOUNT3 = MXGSAER * 50
      INTEGER, PARAMETER :: MXCOUNT4 = MXGSAER * 20
      INTEGER, PARAMETER :: MXHOLD   = 250
!
! !PUBLIC DATA MEMBERS:
! 
      !---------------------------------------------------------------
      ! Character variables
      !---------------------------------------------------------------

      ! 2-character variables
      CHARACTER(LEN=2 )  :: SPECL(MXCOF)
      COMMON /CHAR2/        SPECL

      ! 14-character scalars
      CHARACTER(LEN=4 )  :: DINP
      CHARACTER(LEN=14)  :: JST 
      COMMON /CHAR14_0/     DINP, JST

      ! 14-character arrays #1
      CHARACTER(LEN=14)  :: XINP(IMISC)     
      CHARACTER(LEN=14)  :: RINP(IMISC)
      CHARACTER(LEN=14)  :: SINP(IMISC)     
      COMMON /CHAR14_1/     XINP, RINP, SINP

      ! 14-character arrays #2
      CHARACTER(LEN=14)  :: NAMD   (NMDEAD)
      CHARACTER(LEN=14)  :: NAMEMB (IMASBAL)
      CHARACTER(LEN=14)  :: CHEMTYP(ICS)
      COMMON /CHAR14_2/     NAMEMB, CHEMTYP, NAMD
  
      ! 14-character arrays #3
      CHARACTER(LEN=14)  :: NAMESPEC(0:IGAS,ICS)
      CHARACTER(LEN=14)  :: NAMEGAS (0:IGAS) 
      CHARACTER(LEN=14)  :: NAMENCS (0:MXGSAER,ICS)    
      COMMON /CHAR14_3/     NAMESPEC, NAMEGAS, NAMENCS

      ! 14-character arrays #4
      CHARACTER(LEN=14)  :: NAMEPHOT(NMRPROD,IPHOT)     
      COMMON /CHAR14_4/     NAMEPHOT

      ! 80-character variables
      CHARACTER(LEN=80)  :: HEADING, COMMENT
      COMMON /CHAR80/       HEADING, COMMENT

      !---------------------------------------------------------------
      ! INTEGER and REAL*8 scalars
      !---------------------------------------------------------------

      INTEGER ::         NTSPECGAS,   NMASBAL
      COMMON /CTLLOOP/   NTSPECGAS,   NMASBAL

      ! /CTLLOOP2/ needs to be declared THREADPRIVATE (bmy, 7/16/03)
      INTEGER ::         KTLOOP,      IFSUN
      COMMON /CTLLOOP2/  KTLOOP,      IFSUN
 
      REAL*8  ::         TINTERVAL,   CHEMINTV,  TIME,      OXYCONS
      REAL*8  ::         HMAXNIT,     FRACDEC
      COMMON /XYGRID/    TINTERVAL,   CHEMINTV,  TIME,      OXYCONS,    &
     &                   HMAXNIT,     FRACDEC

      INTEGER ::         IHOUR,       IRCHEM
      INTEGER ::         NPHOTALL,    IFDID,     IFNEVER,   IFNONE   
      COMMON /IXYGD/     IHOUR,       IRCHEM,                           &
     &                   NPHOTALL,    IFDID,     IFNEVER,   IFNONE   

      ! /IXYGD2/ needs to be held THREADPRIVATE.
      INTEGER ::         NCSP
      COMMON /IXYGD2/    NCSP

      REAL*8  ::         HMIN
      COMMON /DGEAR/     HMIN

      ! /DGEAR2/ needs to be held THREADPRIVATE (hamid, bmy, 7/16/03)
      REAL*8  ::         HMAX,        R1DELT,    DELT,      TIMREMAIN
      REAL*8  ::         XELAPS,      TOLD,      RDELT,     XELAPLAST
      REAL*8  ::         RMSERR
      COMMON /DGEAR2/    HMAX,        R1DELT,    DELT,      TIMREMAIN,  &
     &                   XELAPS,      TOLD,      RDELT,     XELAPLAST,  &
     &                   RMSERR

      INTEGER ::         NSFTOT,      NPDTOT,    NSTTOT
      INTEGER ::         IFAILTOT,    LFAILTOT,  NFAILTOT  
      COMMON /IGEAR/     NSFTOT,      NPDTOT,    NSTTOT,                &
     &                   IFAILTOT,    LFAILTOT,  NFAILTOT

      ! /IGEAR2/ has to be declared THREADPRIVATE (bmy, 7/16/03)
      INTEGER ::         NQQ,         NSUBFUN,   NPDERIV   
      INTEGER ::         NFAIL,       IFAIL,     LFAIL
      COMMON /IGEAR2/    NQQ,         NSUBFUN,   NPDERIV,               &
     &                   NFAIL,       IFAIL,     LFAIL

      INTEGER ::         NPHOT,       NPRODLO,   NPRODHI,   MSTEP
      INTEGER ::         MAXORD,      MBETWEEN,  IC3H8,     IC2H6
      COMMON /CHEM2/     NPHOT,       NPRODLO,   NPRODHI,   MSTEP,      &
     &                   MAXORD,      MBETWEEN,  IC3H8,     IC2H6

      ! /CHEM2A/ has to be held THREADPRIVATE (bmy, 7/16/03)
      INTEGER ::         ISCHAN,      NFDH3,     NFDL2,     NFDH2
      INTEGER ::         NFDL1,       NFDH1,     NFDREP,    NFDREP1   
      INTEGER ::         NFDL0,       NALLR 
      COMMON /CHEM2A/    ISCHAN,      NFDH3,     NFDL2,     NFDH2,      &      
     &                   NFDL1,       NFDH1,     NFDREP,    NFDREP1,    &
     &                   NFDL0,       NALLR

      INTEGER ::         NGAS,        NMREAC
      COMMON /CHEM3/     NGAS,        NMREAC

      ! Added NNADDG to /CHEM4/ for DMS+OH+O2 rxn (bdf, bmy, 4/18/03)
      ! Add NKN2O5 to /CHEM4/ to flag N2O5 hydrolysis rxn (mje, bmy, 8/7/03)
      INTEGER ::         NNADD1,      NNADDA,      NNADDB
      INTEGER ::         NNADDC,      NNADDD,      NNADDK
      INTEGER ::         NNADDV,      NNADDZ,      NKO3PHOT
      INTEGER ::         NNADDG,      NEMIS,       NDRYDEP
      INTEGER ::         NKHNO4,      NKN2O5      
      COMMON /CHEM4/     NNADD1,      NNADDA(ICS), NNADDB(  ICS),       &
     &                   NNADDC(ICS), NNADDD(ICS), NNADDK(  ICS),       &
     &                   NNADDV(ICS), NNADDZ,      NKO3PHOT(ICS),       &
     &                   NNADDG(ICS), NEMIS( ICS), NDRYDEP( ICS),       &
     &                   NKHNO4(ICS), NKN2O5

      INTEGER ::         IH2O,        IOXYGEN,   MB1,      MB2
      COMMON /SPECIES/   IH2O,        IOXYGEN,   MB1,      MB2

      ! Added for interannually-varying Methane (bnd, bmy, 7/1/03)
      INTEGER ::         ICH4
      COMMON /SPECIE2/   ICH4

      ! Added for interannually-varying Methane (bnd, bmy, 7/1/03)
      REAL*8  ::         C3090S,      C0030S,    C0030N,   C3090N
      COMMON /SPECIE3/   C3090S,      C0030S,    C0030N,   C3090N

      ! Added for tracking oxidation of ISOP by OH (dkh, bmy, 6/1/06)
      INTEGER ::         ILISOPOH
      COMMON /SPECIE4/   ILISOPOH

      INTEGER ::         IOUT,        KGLC,      KCPD,     IO93
      COMMON /FILES/     IOUT,        KGLC,      KCPD,     IO93

      !---------------------------------------------------------------
      ! INTEGER and REAL*8 arrays
      !---------------------------------------------------------------

      ! Add NKSPECG for DMS+OH+O2 rxn (bdf, bmy, 4/18/03)
      INTEGER NMOTH,NTSPEC,JPHOTRAT,ISGAINR,ISPORL,NOGAINE,NOUSE
      INTEGER NSPEC,NTRATES,ISGAINE,NTLOOPNCS,NSPCSOLV,ISCHANG,NRATES
      INTEGER NM3BOD,ITWOR,ITHRR,INOREP,NRATCUR,NSURFACE,NPRESM,NMAIR
      INTEGER NMO2,NMN2,NNEQ,NARR,NABR,NACR,NABC,NKSPECW,NKSPECX
      INTEGER NKSPECY,NKSPECZ,NKSPECV,ISLOSSR,NKSPECA,NKSPECB,NKSPECC
      INTEGER NKSPECD,NKSPECK,NKSPECG
      COMMON /IDICS/                                                     &
     &    NMOTH(   ICS), NTSPEC( ICS), JPHOTRAT(ICS),                    &
     &    ISGAINR( ICS), ISPORL( ICS), NOGAINE( ICS),   NOUSE(    ICS),  &
     &    NSPEC(   ICS), NTRATES(ICS), ISGAINE( ICS),   NTLOOPNCS(ICS),  & 
     &    NSPCSOLV(ICS), ISCHANG(ICS), NRATES(  ICS),   NM3BOD(   ICS),  &
     &    ITWOR(   ICS), ITHRR(  ICS), INOREP(  ICS),   NRATCUR(  ICS),  &
     &    NSURFACE(ICS), NPRESM( ICS), NMAIR(   ICS),   NMO2(     ICS),  &
     &    NMN2(    ICS), NNEQ(   ICS), NARR(    ICS),   NABR(     ICS),  &
     &    NACR(    ICS), NABC(   ICS), NKSPECW( ICS),   NKSPECX(  ICS),  &
     &    NKSPECY( ICS), NKSPECZ(ICS), NKSPECV(MAXGL2,ICS),ISLOSSR(ICS), &
     &    NKSPECA(  MAXGL3,ICS), NKSPECB(  MAXGL3,ICS),                  &
     &    NKSPECC(MAXGL3,ICS),NKSPECD(MAXGL3,ICS),NKSPECK(MAXGL3,ICS),   &
     &    NKSPECG(MAXGL2,ICS)

      ! re-define some nkspec* arrays for harvard chem mechanism (bdf)
      INTEGER ::        NOLOSP,         NGNFRAC,         NOLOSRAT
      INTEGER ::        IARRAY,         NALLRAT,         KZTLO
      INTEGER ::        KZTHI,          IONER,           NPLLO
      INTEGER ::        NPLHI,          NFRLO,           NFRHI
      INTEGER ::        NPDLO,          NPDHI,           IZLO
      INTEGER ::        JZLO,           JLLO,            JGLO
      INTEGER ::        IRMCT
      COMMON /IICP/     NOLOSP(ICP),    NGNFRAC(ICP),    NOLOSRAT(ICP),  &
     &                  IARRAY(ICP),    NALLRAT(ICP),    KZTLO(   ICP),  &
     &                  KZTHI( ICP),    IONER(  ICP),    NPLLO(   ICP),  &
     &                  NPLHI( ICP),    NFRLO(  ICP),    NFRHI(   ICP),  &
     &                  NPDLO( ICP),    NPDHI(  ICP),    IZLO (   ICP),  &
     &                  JZLO ( ICP),    JLLO(   ICP),    JGLO(    ICP),  &
     &                  IRMCT( ICP)

      REAL*8  ::        ABTOL,          ABST2
      REAL*8  ::        ERRMAX,         HMAXUSE,         TIMEINTV
      COMMON /DICS/     ABTOL(6,ICS),   ABST2(ICS),                      &
     &                  ERRMAX(ICS),    HMAXUSE(ICP),    TIMEINTV(ICS)    

      REAL*8  ::        WTGAS,          GQSUMINI
      REAL*8  ::        BSUMCHEM,       GQSUM,           QBKGAS
      COMMON /DIGAS/    WTGAS(IGAS),    GQSUMINI(IGAS),                  &
     &                  BSUMCHEM(IGAS), GQSUM(   IGAS),  QBKGAS(IGAS)       

      REAL*8  ::        CPREV,          CMODEL,          APORL
      COMMON /DMXGAER1/ CPREV(MXGSAER), CMODEL(MXGSAER), APORL(MXGSAER)

      INTEGER ::        IFPRGAS,        LGNUM,           NGMIX
      COMMON /IIGAS/    IFPRGAS(IGAS),  LGNUM(IGAS),     NGMIX(IGAS)

      REAL*8  ::        DEFPRAT
      COMMON /DIPHOT/   DEFPRAT(MXRATE,ICS)

      REAL*8  ::        ARRT,           BRRT
      REAL*8  ::        FCVT,           FCT1T,           FCT2T
      COMMON /DMXCOF/   ARRT(MXCOF),    BRRT( MXCOF),                    &
     &                  FCVT(MXCOF),    FCT1T(MXCOF),    FCT2T(MXCOF)

      INTEGER ::        KCRRT   
      COMMON /IMXCOF/   KCRRT(MXCOF)

      INTEGER NKARR,NKABR,NKACR,NKABC,IRORD
      COMMON /INMRAT2/                                                   &
        NKARR(NMTRATE,ICS),  NKABR(NMTRATE,ICS), NKACR(NMTRATE,ICS),     &
        NKABC(NMTRATE,ICS),  IRORD(NMTRATE,ICS)

      REAL*8 ARR,BRR,FCV,FCTEMP1,FCTEMP2
      COMMON /DNMTRATE/                                                  &
     &   ARR(    NMTRATE, ICS),  BRR(    NMTRATE, ICS),                  &
     &   FCV(    NMTRATE, ICS),  FCTEMP1(NMTRATE, ICS),                  &
     &   FCTEMP2(NMTRATE, ICS)  

      INTEGER IAPROD,NOLOSRN,NRUSE,NRREP,NPRODUC,IALLOSN,NCEQUAT
      INTEGER NEWFOLD,NKONER,NKTWOR,NKTHRR,IRMA,IRMB,KCRR,LSKIP,IRMC
      INTEGER JPHOTNK,IUSED,NOLDFNEW
      COMMON /INMTRATE/                                                  &
     &   IAPROD( NMTRATE,  ICS),  NOLOSRN( NMTRATE, ICS),                &
     &   NRUSE(  NMTRATE,  ICS),  NRREP(   NMTRATE, ICS),                &
     &   NPRODUC(NMTRATE,  ICS),  IALLOSN( MXRATE,  ICS),                &
     &   NCEQUAT(NMTRATE,  ICS),  NOLDFNEW(NMTRATE, ICS),                &
     &   NEWFOLD(NMTRATE*2,ICS),  NKONER(  NMTRATE, ICS),                &
     &   NKTWOR( NMTRATE,  ICS),  NKTHRR(  NMTRATE, ICS),                &
     &   KCRR(   NMTRATE,  ICS),  LSKIP(   MXRATE,  ICS),                &
     &   IRMC(   NMTRATE      ),  JPHOTNK( NMTRATE, ICS),                &
     &   IUSED(  MXRATE,   ICS)

      ! For COMPAQ, put IRMA, IRMB in /INMTRATE2/ (Q. Liang, bmy, 10/17/05)
      COMMON /INMTRATE2/                                                 &
     & IRMA(   NMTRATE      ),  IRMB(    NMTRATE     )

      INTEGER ::       NEWNK
      COMMON /IMAXGL3/ NEWNK(MAXGL)

      REAL*8  ::       FRACP 
      COMMON /DMAXGL2/ FRACP(MAXGL, ICS)

      INTEGER NREACOTH,LGASBINO,NKNLOSP,LOSINACP,IGNFRAC,NKGNFRAC
      INTEGER NREACAIR,NREAC3B,NREACEQ,NREQOTH,NREACN2,NREACO2,NREACPM
      INTEGER LGAS3BOD,NKSURF,NCOATG
      COMMON /IMAXGL2/                     &
     &   NREACOTH(MAXGL2,ICS), LGASBINO(MAXGL2,ICS),                     &
     &   NKNLOSP( MAXGL3,ICS), LOSINACP(MAXGL3,ICS),                     &
     &   IGNFRAC( MAXGL, ICS), NKGNFRAC(MAXGL, ICS),                     &
     &   NREACAIR(MAXGL3,ICS), NREAC3B( MAXGL3,ICS),                     &
     &   NREACEQ( MAXGL3,ICS), NREQOTH( MAXGL3,ICS),                     &
     &   NREACN2( MAXGL3,ICS), NREACO2( MAXGL3,ICS),                     &
     &   NREACPM( MAXGL3,ICS), LGAS3BOD(MAXGL3,ICS),                     &
     &   NKSURF(  MAXGL4    ), NCOATG(  MAXGL4    ) 

      INTEGER ::        MBCOMP,            MBTRACE
      COMMON /IIMASBAL/ MBCOMP(IMASBAL,2), MBTRACE(IMASBAL)

      ! /DKBLOOP2/ needs to be declared THREADPRIVATE
      REAL*8 CNEW,CEST,GLOSS,CHOLD,VDIAG,CBLK,DTLOS,EXPLIC,CONC
      REAL*8 RRATE,URATE,TRATE,CORIG
      COMMON /DKBLOOP2/              &
     &   CNEW(   KBLOOP,  MXGSAER),                                      &
     &   CEST(   KBLOOP,  MXGSAER),                                      &
     &   GLOSS(  KBLOOP,  MXGSAER),                                      &
     &   CHOLD(  KBLOOP,  MXGSAER),                                      &
     &   VDIAG(  KBLOOP,  MXGSAER),  CBLK(  KBLOOP,MXGSAER),             &
     &   DTLOS(  KBLOOP,  MXGSAER),  EXPLIC(KBLOOP,MXGSAER),             &
     &   CONC(   KBLOOP,MXGSAER*7),                                      &
     &   RRATE(  KBLOOP,  NMTRATE),                                      &
     &   URATE(  KBLOOP,NMTRATE,3),                                      &
     &   TRATE(  KBLOOP,NMTRATE*2),                                      &
     &   CORIG(  KBLOOP,  MXGSAER)

      ! /DKBLOOP0/ needs to be held THREADPRIVATE
      REAL*8 ::         CC2
      COMMON /DKBLOOP0/ CC2(KBLOOP,0:MXARRAY)  

      INTEGER NKPHOTRAT,NPPHOTRAT,NKNPHOTRT
      COMMON /DIPHOT2/                                                   &
     &   NKPHOTRAT(IPHOT,ICS),      NPPHOTRAT(IPHOT,ICS),                &
     &   NKNPHOTRT(IPHOT,ICS) 

      REAL*8 ::        FRACGAIN,              QBKCHEM
      COMMON /DIMXGS2/ FRACGAIN(MXGSAER,ICS), QBKCHEM( MXGSAER,ICS) 

      INTEGER NUMLOST,NUMGFRT,NUMLOSS,JPORL,NUMGAINT,NGAINE,NUMGAIN
      INTEGER IGAINR,IPORL,IGAINE,ISOLVSPC,INEWOLD,MAPPL,ISAPORL,NUMPORL
      INTEGER ISPARDER,JLOSST
      COMMON /IIMXGS2/                                                   &
     &   NUMLOST( MXGSAER,    ICS), NUMGFRT( MXGSAER,        ICS),       &
     &   NUMLOSS( MXGSAER,    ICP), JPORL(   MXGSAER,MAXGL,  ICS),       &
     &   NUMGAINT(MXGSAER,    ICS), NGAINE(  MXGSAER,        ICS),       &
     &   NUMGAIN( MXGSAER,    ICP), IGAINR(  MXGSAER,        ICS),       &
     &   IPORL(   MXGSAER,    ICS), IGAINE(  MXGSAER,        ICS),       &
     &   ISOLVSPC(MXGSAER,    ICS), INEWOLD( MXGSAER,        ICS),       &
     &   MAPPL(   MXGSAER,    ICS), ISAPORL( MXGSAER            ),       &
     &   NUMPORL( MXGSAER,    ICP), ISPARDER(MXGSAER,MXGSAER    ),       &
     &   JLOSST(  MXGSAER,MAXGL,ICS)

      INTEGER JZILCH,KZILCH,MZILCH
      COMMON /IGMXGLS/                                                   &
     &     JZILCH(MXGSAER),  KZILCH(MXGSAER),  MZILCH(MXGSAER) 

      INTEGER LZERO,JARRAYPT,IZILCH,JARRDIAG,JLOZ1,JHIZ1,IJTLO
      INTEGER IJTHI,IMZTOT,IFREPRO,IZLO1,IZLO2,IZHI0,IZHI1
      COMMON /IGMXGS2/                                                   &
     &   LZERO(   MXGSAER,MXGSAER), JARRAYPT(MXGSAER,MXGSAER),           &
     &   IZILCH(  MXGSAER,MXGSAER), JARRDIAG(MXGSAER,   ICP),            &
     &   JLOZ1(   MXGSAER,    ICP), JHIZ1(   MXGSAER,   ICP),            &
     &   IJTLO(   MXGSAER,    ICP), IJTHI(   MXGSAER,   ICP),            &
     &   IMZTOT(  MXGSAER,    ICP), IFREPRO( MXGSAER,MXRATE, ICS),       &
     &   IZLO1( MXGSAER,ICP),                                            &
     &   IZLO2( MXGSAER,ICP), IZHI0( MXGSAER,ICP), IZHI1( MXGSAER,ICP)    

      REAL*8 ::        FRACNFR,           FRACPL
      COMMON /DMXCOUN/ FRACNFR(MXCOUNT4), FRACPL(MXCOUNT2) 

      INTEGER JZERO,KZERO,MZERO,IZEROK,JZEROA,IKDECA,KJDECA,LOSSRA
      INTEGER IKDECB,KJDECB,LOSSRB,IKDECC,KJDECC,LOSSRC,IKDECD,KJDECD
      INTEGER LOSSRD,IKDECE,KJDECE,LOSSRE,KZEROA,MZEROA,KZEROB,MZEROB
      INTEGER KZEROC,MZEROC,KZEROD,MZEROD,KZEROE,MZEROE,IPOSPD,IIALPD
      INTEGER NKPDTERM,IJVAL,IKZTOT,JSPNPL,NKNFR,JSPCNFR
      COMMON /IMXCOUN/                                                   &
     &   JZERO( MXCOUNT3),  KZERO(  MXCOUNT3),  MZERO(   MXCOUNT3),      & 
     &   IZEROK(MXCOUNT2),  JZEROA( MXCOUNT3),                           &
     &   IKDECA(MXCOUNT3),  KJDECA( MXCOUNT3),  LOSSRA(  MXCOUNT4),      &
     &   IKDECB(MXCOUNT3),  KJDECB( MXCOUNT3),  LOSSRB(  MXCOUNT4),      &
     &   IKDECC(MXCOUNT3),  KJDECC( MXCOUNT3),  LOSSRC(  MXCOUNT4),      &
     &   IKDECD(MXCOUNT3),  KJDECD( MXCOUNT3),  LOSSRD(  MXCOUNT4),      &
     &   IKDECE(MXCOUNT3),  KJDECE( MXCOUNT3),  LOSSRE(  MXCOUNT4),      &
     &   KZEROA(MXCOUNT4),  MZEROA( MXCOUNT4),                           &
     &   KZEROB(MXCOUNT4),  MZEROB( MXCOUNT4),                           &
     &   KZEROC(MXCOUNT4),  MZEROC( MXCOUNT4),                           &
     &   KZEROD(MXCOUNT4),  MZEROD( MXCOUNT4),                           & 
     &   KZEROE(MXCOUNT4),  MZEROE( MXCOUNT4),                           & 
     &   IPOSPD(MXCOUNT2),  IIALPD( MXCOUNT2),  NKPDTERM(MXCOUNT2),      &
     &   IJVAL( MXCOUNT3),  IKZTOT( MXCOUNT4),  JSPNPL(  MXCOUNT4),      &
     &   NKNFR( MXCOUNT4),  JSPCNFR(MXCOUNT4) 

      INTEGER IDH5,IDH4,IDH3,IDH2,IDH1,IDL5,IDL4,IDL3,IDL2,IDL1,KBH5
      INTEGER KBH4,KBH3,KBH2,KBH1,KBL5,KBL4,KBL3,KBL2,KBL1,MBH5,MBH4
      INTEGER MBH3,MBH2,MBH1,MBL5,MBL4,MBL3,MBL2,MBL1,NPH5,NPH4,NPH3
      INTEGER NPH2,NPH1,NPL5,NPL4,NPL3,NPL2,NPL1
      COMMON /IMXCOU2/                                                   &
     &   IDH5(  MXCOUNT3),  IDH4(  MXCOUNT3),  IDH3(  MXCOUNT3),         &
     &   IDH2(  MXCOUNT3),  IDH1(  MXCOUNT3),  IDL5(  MXCOUNT3),         &
     &   IDL4(  MXCOUNT3),  IDL3(  MXCOUNT3),  IDL2(  MXCOUNT3),         &
     &   IDL1(  MXCOUNT3),                                               &
     &   KBH5(  MXCOUNT4),  KBH4(  MXCOUNT4),  KBH3(  MXCOUNT4),         &
     &   KBH2(  MXCOUNT4),  KBH1(  MXCOUNT4),  KBL5(  MXCOUNT4),         &
     &   KBL4(  MXCOUNT4),  KBL3(  MXCOUNT4),  KBL2(  MXCOUNT4),         &
     &   KBL1(  MXCOUNT4),                                               & 
     &   MBH5(  MXCOUNT4),  MBH4(  MXCOUNT4),  MBH3(  MXCOUNT4),         &
     &   MBH2(  MXCOUNT4),  MBH1(  MXCOUNT4),  MBL5(  MXCOUNT4),         &
     &   MBL4(  MXCOUNT4),  MBL3(  MXCOUNT4),  MBL2(  MXCOUNT4),         &
     &   MBL1(  MXCOUNT4),                                               &
     &   NPH5(  MXCOUNT4),  NPH4(  MXCOUNT4),  NPH3(  MXCOUNT4),         &
     &   NPH2(  MXCOUNT4),  NPH1(  MXCOUNT4),  NPL5(  MXCOUNT4),         &
     &   NPL4(  MXCOUNT4),  NPL3(  MXCOUNT4),  NPL2(  MXCOUNT4),         &
     &   NPL1(  MXCOUNT4)  

      REAL*8 ::          WTMB
      COMMON /DIMASBAL2/ WTMB(IMASBAL,MXGSAER,2)

      INTEGER ::         JMBCOMP
      COMMON /IIMASBAL2/ JMBCOMP(IMASBAL,MXGSAER,2) 

      REAL*8 ::          FKOEF
      REAL*8 ::          FK2 
      COMMON /DNMRPROD2/ FKOEF(NMRPROD,NMTRATE,ICS),                     &
     &                   FK2(  NMRPROD,NMTRATE,ICS)  

      INTEGER ::         IRM
      INTEGER ::         IRM2
      COMMON /INMRPROD2/ IRM( NMRPROD,NMTRATE,ICS),                      &
     &                   IRM2(NMRPROD,NMTRATE,ICS)                  

      REAL*8 ::          ASET,       PINP,     CVAR,     O3DOBS
      COMMON /DMISC/     ASET(10,8), PINP(20), CVAR(15), O3DOBS(12)  

      REAL*8 ::          ENQQ2,         ENQQ3,          CONPST
      REAL*8 ::          ENQQ1,         CONP15
      COMMON /IORDR/     ENQQ2(MORDER), ENQQ3( MORDER), CONPST(MORDER),  &
     &                   ENQQ1(MORDER), CONP15(MORDER)

      REAL*8 ::          PERTS2,           PERTST
      COMMON /DMORD/     PERTS2(MORDER,3), PERTST(MORDER,3)

      REAL*8 ::          XGDFCF, ASTKCF, RUARSL, RH100
      COMMON /XARSOL/    XGDFCF, ASTKCF, RUARSL, RH100

      INTEGER ::         IARSFA, MWARSL, MNTHARS
      COMMON /IARSOL/    IARSFA, MWARSL, MNTHARS

      INTEGER ::         NKEMIS,              NTEMIS
      INTEGER ::         NKDRY,               NTDEP
      COMMON /IMAXGL4/   NKEMIS(MAXGL3,ICS),  NTEMIS(MAXGL3,ICS),        &
     &                   NKDRY (MAXGL3,ICS),  NTDEP( MAXGL3)

      !=================================================================
      ! Declare some common blocks THREADPRIVATE for the OpenMP
      ! parallelization (bdf, bmy, 4/1/03)
      !=================================================================
!$OMP THREADPRIVATE( /CHEM2A/   )  
!$OMP THREADPRIVATE( /CTLLOOP2/ )
!$OMP THREADPRIVATE( /DGEAR2/   )
!$OMP THREADPRIVATE( /DKBLOOP0/ )
!$OMP THREADPRIVATE( /DKBLOOP2/ )
!$OMP THREADPRIVATE( /IGEAR2/   )
!$OMP THREADPRIVATE( /IXYGD2/   )
!EOP
!------------------------------------------------------------------------------
