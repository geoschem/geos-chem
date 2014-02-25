!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: comode_loop_mod
!
! !DESCRIPTION: Module file COMODE contains common blocks and variables for 
!  SMVGEAR II.
!\\
!\\
! !INTERFACE: 
!
      MODULE COMODE_LOOP_MOD
!
! !USES:
!
      USE CMN_SIZE_MOD

      IMPLICIT NONE
      PUBLIC
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC :: Init_COMODE_LOOP
      PUBLIC :: Cleanup_COMODE_LOOP
!
! !REMARKS:
!         CCCCCCC  OOOOOOO  M     M  OOOOOOO  DDDDDD   EEEEEEE 
!         C        O     O  M M M M  O     O  D     D  E  
!         C        O     O  M  M  M  O     O  D     D  EEEEEEE  
!         C        O     O  M     M  O     O  D     D  E  
!         CCCCCCC  OOOOOOO  M     M  OOOOOOO  DDDDDD   EEEEEEE 
!                                                                             .
! *********************************************************************
! * THIS IS THE COMMON BLOCK FOR "SMVGEAR" AND "MIE," TWO ORDINARY    *
! * DIFFERENTIAL EQUATION SOLVERS. THE REFERENCE FOR THE CODES IS     *
! *                                                                   * 
! *   JACOBSON M. Z. AND TURCO R. P. (1993) SMVGEAER: A SPARSE-       * 
! *    MATRIX, VECTORIZED GEAR CODE FOR ATMOSPHERIC MODELS.           *
! *    SUBMITTED TO ATMOSPHERIC ENVIRONMENT, PART A. MAY 20, 1993     * 
! *                                                                   * 
! * COMODE.H SETS PARAMETER VALUES AND SERVES AS A COMMON BLOCK FOR   *
! * ALL DIMENSIONED AND NON-DIMENSIONED VARIABLES. COMODE.H ALSO      *
! * DEFINES EACH PARAMETER, BUT DATA FILE DEFINE.DAT EXPLAINS NON-    *
! * DIMENSIONED VARIABLES. INDIVIDUAL SUBROUTINES DEFINE DIMENSIONED  *
! * VARIABLES.                                                        *
! *********************************************************************
!                                                                             .
! *********************************************************************
! *                       SET PARAMETERS                              *
! *********************************************************************
!                                                                             .
! ****************** COORDINATE-SYSTEM PARAMETERS *********************
!
! ILAT     = MAXIMUM NUMBER OF LATITUDE(ILAT)   GRID POINTS
! ILONG    = MAXIMUM NUMBER OF LONGITUDE(ILONG) GRID POINTS
! IMLOOP   = ILAT * ILONG - USED FOR MORE EFFICIENT ARRAYS
! IVERT    = MAXIMUM NUMBER OF LAYERS 
! ILAYER   = MAXIMUM OF LAYER BOUNDARIES
! KBLOOP   = MAXIMUM NUMBER OF GRID POINTS IN A VECTORIZED BLOCK 
!            SHOULD RANGE FROM 512 (BELOW WHICH VECTORIZATION DECREASES)
!            TO 1024 (ABOVE WHICH, ARRAY SPACE IS LIMITED)
! MXBLOCK  = MAXIMUM NUMBER OF GRID POINT BLOCKS
! MAXDAYS  = MAXIMUM NUMBER OF DAYS FOR THE MODEL TO RUN
!
! !REVISION HISTORY:
!  (M. Jacobson 1997; bdf, bmy, 4/23/03, 6/1/06)
!  (1 ) Removed many commented-out common blocks not needed for GEOS-CHEM.
!        Also updated comments.  Also make sure that MAXGL3 is dimensioned
!        for at least NNPAR tracers.  Add NNADDG and NKSPECG for DMS+OH+O2
!        rxn.  COEF12 and QRM2 are now obsolete for SMVGEAR II. (bmy, 4/23/03)
!  (2 ) Added ICH4 to the /SPECIE2/ common block for interannual-varying
!        CH4 concentration.  Added variables for latitude distribution of
!        CH4 to the /SPECIE3/ common block. (bmy, 7/1/03)
!  (3 ) Added ITS_NOT_A_ND65_FAMILY to the /LPL/ common block for the ND65
!        production/loss diagnostic.   Comment out counter variables, you can 
!        get the same info w/ a profiling run.  Updated comments, cosmetic 
!        changes. (bmy, 7/9/03)
!  (4 ) Removed the following variables from common blocks which are not needed
!        for GEOS-CHEM: COLENG, AERSURF, VHMET1, VHMET, VMET3, CINIT, RHO3,
!        GRIDVH, CSUMA1, XELRAT, T1BEG, T2BEG, T1FIN, T2FIN, DECLIN, RAGSUT,
!        SINDEC, COSDEC, SIGMAL, PRESSL, RHOA, DSIG_SMV, TEMPL, VMET, SIGDIF,
!        TMORN, PRESSC, XLAT, XLON, DMERIDUT, GRIDAREA, DSX, XLONUT, DSY, 
!        SINXLAT, COSXLAT, HMETT, HMET1, HMET2, RSET, RRIS, TZDIF, ZENRAT0,
!        ZENRAT1, MLOPJ, REORDER_SAVE, RHO3K, GRIDVH3K, FIELDXY, FIELDYZ,
!        FIELDXZ, RATMIX, GQSCHEM, C, QPRODA, QPRODB, QPRODC, QPRODD, QPROD,
!        CINP, NUMSDT, NKSDT, PRATE.  MONTHP, KYEAR, LDMONTH, ININT, ICLO, 
!        JCLO, FIELD1, MZLO, MZLO2, MZHI0, MZHI1, KZLO1, KZLO2, KZHI0, KZHI1, 
!        IHIZ1, IHIZ2, IHIZ3, PRESS5KM, KGRP, IABOVK, MROTAT1, MINROT1, 
!        NUMSUBS, LSPECEMIS, MROTAT2, MINROT2, MAXPOS, NOGAINR, NOLOSSR,
!        MAXSTEPS, YLOW, HMAXDAY, KPHT, KRDD, KMIX, KINS, KGCO, ABHSUMK, DX0, 
!        DY0, XU0, DTOUT, CONPSUR, DXLONG, DYLAT, SWLONDC, CONSTIM, SWLATDC, 
!        UTSECY, TOTSEC, FINHOUR, FINMIN, FINSEC, TFROMID, ZENFIXED, ZENITH, 
!        DENCONS, HALFDAY, GRAVC, FOURPI, TWOPI, REARTH, RPRIMB, AVOG1,
!        HALF, THIRD, THRPI2, PID180, PID2, SCTWOPI, AMRGAS, TWPISC.  
!        This should free up more memory for runs. (bmy, 7/16/03)
!  (5 ) Split off NOCC into the /CHEM3B/ common block, since it doesn't need 
!        to be held THREADPRIVATE.  Removed /DKBLOOP/ and /DKBLOOP5/, since
!        these contain variables which are used locally within either
!        "calcrate.f" or "smvgear.f".  Cosmetic changes. (bmy, 7/28/03)
!  (6 ) Add NKN2O5 to /CHEM4/ common block to flag N2O5 hydrolysis rxn.
!         (mje, bmy, 8/7/03)
!  (7 ) Eliminated SMALLCHEM cpp switch (bmy, 12/2/03)
!  (8 ) Now set MAXGL3 = NNPAR for new # of tracers (bmy, 4/6/04)
!  (9 ) Remove obsolete LGEOSCO and FULLCHEM Cpp switches (bmy, 6/24/05)
!  (10) For COMPAQ, put IRMA, IRMB in /INMTRATE2/ common block.  For COMPAQ, 
!        also declare /INMTRATE2/ THREADPRIVATE. (Q. Liang, bmy, 10/17/05)
!  (11) Now remove AVG, BOLTG, RGAS, SCDAY, BK, EIGHTDPI, RSTARG, WTAIR,
!        ONEPI, CONSVAP, SMAL1, SMAL2, SMAL3 from common blocks and declare 
!        these as parameters. (bec, bmy, 3/29/06)
!  (12) Added ILISOPOH, the index of ISOP lost to OH (dkh, bmy, 6/1/06)
!  (13) Added NKHO2 to /CHEM4/ common block to flag HO2 aerosol uptake
!        (jaegle 02/26/09)
!  (14) Add NNADDF and NNADDH to /CHEM4/ for HOC2H4O rxns
!       Add NKHOROI and NKHOROJ to /CHEM4/ for HOC2H4O rxns in EP photolysis
!       (tmf, 3/6/09)
!  (15) Added NKSPECF, NKSPECH to /IDICS/ for C2H4 chemistry (tmf, 3/6/09) 
!  (16) Increase IGAS, MAXGL, MAXGL2, NMRATE, IPHOT (tmf, 3/6/09)
!  (17) Add RRATE_FOR_KPP variable to DKBLOOP2 common block (phs,ks,dhk, 09/15/09)
!  (18) PINP(20) increased to PINP(IMISC) (FP 2/10)
!  (19) Added NKO3PHOTH2 (bhh, jmao, eam, 7/18/11)
!
!  !REVISION HISTORY:
!  23 Aug 2011 - M. Long     - Converted to Module from Header file
!  14 Aug 2012 - R. Yantosca - Add #include "define.h" to USES section
!  14 Aug 2012 - R. Yantosca - Now set KBLOOP=1 only if we are connecting to
!                              an external GCM (i.e. triggered by switches
!                              EXTERNAL_GRID or EXTERNAL_FORCING).
!  17 Oct 2012 - R. Yantosca - Do not define ILAYER if compiling w/ -DDEVEL
!  24 Oct 2012 - R. Yantosca - Rewrite GIGC #ifdef blocks for clarity
!  19 Nov 2012 - R. Yantosca - Now use funcitions INIT_COMODE_LOOP and
!                              CLEANUP_COMODE_LOOP to allocate/free arrays
!  13 Aug 2013 - M. Sulprizio- Add modifications for updated SOA (H. Pye)
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS: 
!
      INTEGER, PARAMETER :: IPLUME  = 0       ! Set to zero, we don't use this
      INTEGER, PARAMETER :: MAXDAYS = 1000    ! This is not used in G-C

#if defined( EXTERNAL_GRID ) || defined( EXTERNAL_FORCING )
      !-----------------------------------------------------------------
      !     %%%%% CONNECTING TO GEOS-5 GCM via ESMF INTERFACE %%%%%
      !
      ! KBLOOP is the # of boxes that SMVGEAR will process per CPU.
      ! Set KBLOOP=1 for connecting to an external GCM
      !
      ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      ! %%% NOTE: If you are using GEOS-Chem without ESMF, but with  %%%
      ! %%% the DEVEL=yes option (i.e. to test grid-independent      %%%
      ! %%% updates w/r/t a standard G-C simulation), then you must  %%%
      ! %%% make sure that KBLOOP is set to the same value in both   %%%
      ! %%% simulations.                                             %%%
      ! %%%                                                          %%%
      ! %%% The absolute and relative errors (which determine if the %%%
      ! %%% chemistry has converged to a solution) are computed over %%% 
      ! %%% all KBLOOP boxes at once.  Using different KBLOOP values %%%
      ! %%% in different simulations will cause slightly different   %%%
      ! %%% results in chemical concentrations (even after only one  %%%
      ! %%% timestep).                                               %%%
      ! %%%                                                          %%%
      ! %%% To this end, we now only set KBLOOP=1 if we are          %%%
      ! %%% connecting GEOS-Chem to an external GCM (i.e. if the Cpp %%%
      ! %%% switches EXTERNAL_GRID or EXTERNAL_FORCING are set in    %%%
      ! %%% define.h).                                               %%%
      ! %%%                                                          %%%
      ! %%%    -- Bob Yantosca (14 Aug 2012)                         %%%
      ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      !-----------------------------------------------------------------
      INTEGER, PARAMETER :: KBLOOP  = 1
#else
      !-----------------------------------------------------------------
      !              %%%%% TRADITIONAL GEOS-Chem %%%%%
      !
      ! KBLOOP is the # of boxes that SMVGEAR will process per CPU.
      ! For "traditional" G-C simulations, leave KBLOOP = 24
      !-----------------------------------------------------------------
      INTEGER, PARAMETER :: KBLOOP  = 24
#endif

      ! These values need to be defined as variables instead of
      ! parameters so that they can be defined dynamically.
      INTEGER            :: ILAT
      INTEGER            :: ILONG
      INTEGER            :: IVERT
      INTEGER            :: IPVERT
      INTEGER            :: ITLOOP
      INTEGER            :: IMLOOP
      INTEGER            :: ILAYER
      INTEGER            :: ILTLOOP
      INTEGER            :: MXBLOCK
!
! ************************* TRACER PARAMETERS ****************************
! IDEMS    = EMISSION ID NUMBER (WHICH SPECIES)
      INTEGER IDEMS(NEMPARA+NEMPARB)
      ! NEMPARA = max no. of anthropogenic emissions
      ! NEMPARB = max no. of biogenic emissions
!
! ************************* GAS-PHASE PARAMETERS **********************
!
! IGAS    = MAXIMUM NUMBER OF GASES, ACTIVE + INACTIVE 
! IAERTY  = MAXIMUM NUMBER OF AQUEOUS CHEMISTRY SPECIES (SET = 1
!           HERE SINCE NO AQUEOUS CHEMISTRY IS INCLUDED)
! NMRATE  = MAXIMUM NUMBER OF RATES CONSTANTS (MAX # REACTIONS)
! IPHOT   = MAXIMUM NUMBER OF PHOTO RATES
! NMTRATE = MAXIMUM NUMBER OF KINETIC + PHOTO REACTIONS 
! NMQRATE = MAXIMUM NUMBER OF AQUEOUS CHEMICAL REACTIONS PLUS PHOTO
!           PROCESSES (SET = 1 HERE SINCE NO AQUEOUS CHEMISTRY INCLUDED)
! NMRPROD = MAXIMUM NUMBER OF SPECIES IN A REACTION RATE 
! NMDEAD  = MAXIMUM NUMBER OF DEAD SPECIES
! MAXGL   = MAXIMUM NUMBER OF GAINS / LOSSES FOR EACH CHEMICAL SPECIES.
! MAXGL2  = AN ARRAY DIMENSION SMALLER THAN MAXGL
! MAXGL3  = AN ARRAY DIMENSION SMALLER THAN MAXGL2
! ICS     = NUMBER OF TYPES OF CHEMISTRY: up to 3 for gas phase
! ICP     = TYPES OF CHEMISTRY * 2 (ONE FOR DAY, ONE FOR NIGHT)
! MORDER  = MAXIMUM ORDER FOR GEAR PARAMETERS FOR DIMENSION PURPOSES
!
      INTEGER IGAS,IAERTY,NMRATE,IPHOT,NMTRATE,NMQRATE,NMRPROD,NMDEAD
      INTEGER MAXGL,MAXGL2,MAXGL3,MAXGL4,ICS,ICP,MORDER,IPHOT8,IMISC
      INTEGER IMASBAL,IALTS,MXCOF

      ! Updated for SMVGEAR II (bdf, bmy, 4/1/03)
      ! +4 for SEAC4RS FF passive tracers (jaf, 6/13/13)
      ! +4 for SEAC4RS BB passive tracers (jaf, 6/14/13)
      ! +1 for SEAC4RS ISOP passive tracer (jaf, 6/14/13)
      !PARAMETER ( IGAS    = 310,               IAERTY  = 1           )
      PARAMETER ( IGAS    = 320,               IAERTY  = 1           )
      ! +8 for SEAC4RS FF passive tracers (jaf, 6/13/13)
      ! +8 for SEAC4RS BB passive tracers (jaf, 6/14/13)
      ! +2 for SEAC4RS ISOP passive tracer (jaf, 6/14/13)
      !PARAMETER ( NMRATE  = 726,               IPHOT   = 92          )
      PARAMETER ( NMRATE  = 744,               IPHOT   = 92          )
      PARAMETER ( NMTRATE = NMRATE + IPHOT,    NMQRATE = 1           ) 
      PARAMETER ( NMRPROD = 25,                NMDEAD  = 100         )
      PARAMETER ( MAXGL   = 1000,              MAXGL2  = 90          )
      PARAMETER ( MAXGL3  = NNPAR,             MAXGL4  = 10          )
      PARAMETER ( ICS     = 3,                 ICP     = 2*ICS       ) 
      PARAMETER ( MORDER  = 7                                        ) 
      PARAMETER ( IPHOT8  = IPHOT + 8,         IMISC   = 100         )
      PARAMETER ( IMASBAL = 9,                 IALTS   = 22          )
      PARAMETER ( MXCOF   = 5                                        )
!
! ****************** PARAMETERS TO MINIMIZE ARRAY SPACE ***************
!
! MXGSAER    = LARGER OF IGAS, IAERTY 
! MXRATE     = LARGER OF NMTRATE, NMQRATE 
! MXCC2      = LARGER OF MXGSAER, MXARRAY  
! MXCOUNT1.. = ARRAYS SIZES USED TO MINIMIZE MATRIX SPACE
! MXARRAY    = MAXIMUM ONE-DIMENSIONAL ARRAY-LENGTH OF SPARSE MATRIX 
!
      INTEGER MXGSAER,MXRATE,MXARRAY,MXCC2,MXCOUNT1,MXCOUNT2,MXCOUNT3,
     &        MXCOUNT4,MXHOLD
      PARAMETER( MXGSAER  = IGAS)  
      PARAMETER( MXRATE   = NMTRATE)  
      PARAMETER( MXARRAY  = 3000)
      PARAMETER( MXCC2    = MXARRAY)  
      PARAMETER( MXCOUNT1 = MXGSAER * MAXGL3 * 3)  
      PARAMETER( MXCOUNT2 = MXGSAER * MAXGL3 * 7)  
      PARAMETER( MXCOUNT3 = MXGSAER * 50)
      PARAMETER( MXCOUNT4 = MXGSAER * 20)
      PARAMETER( MXHOLD   = 250)
!
!
! *********************************************************************
! *                    SET CHARACTER LENGTHS & DIMS                   *
! *********************************************************************
!
      CHARACTER*14 NAMESPEC(0:IGAS,ICS), NAMD(    NMDEAD),
     &             APGASA,   APGASB,   APGASC,   APGASD,   APGASE,
     &             APGASF,   APGASG,   APGASH,   IFSORM,     
     &             XINP(     IMISC),   RINP(      IMISC)
      CHARACTER*14 NAMEGAS(0:IGAS),    NAMEMB(IMASBAL),
     &             JST, NAMENCS(0:MXGSAER,ICS), ACORNER, SINP( IMISC), 
     &             NAMEPHOT(NMRPROD,IPHOT), CHEMTYP( ICS )
      CHARACTER*4  DINP,      DINPLAST
      CHARACTER*2  SPECIAL,   SPECL( MXCOF)
!
      CHARACTER*80 HEADING, COMMENT 
!
! *********************************************************************
! *           SET REAL AND INTEGER NON-ARRAY VARIABLES                *
! *********************************************************************
!
      !---------------------------------------------------------------
      ! Physical constants 
      ! (now make these PARAMETERS instead of COMMON block variables)
      !---------------------------------------------------------------

      ! Avogadro's #
      REAL*8, PARAMETER :: AVG      = 6.02252d+23

      ! Boltzmann's constant [erg/K]
      REAL*8, PARAMETER :: BK       = 1.38054d-16
      REAL*8, PARAMETER :: BOLTG    = 1.38054d-16

      ! Condensation vapor pressure ?????? 
      REAL*8, PARAMETER :: CONSVAP  = 6.1078d+03 / BOLTG
      
      ! PI (same value as in CMN_GCTM) and related quantities
      REAL*8, PARAMETER :: ONEPI    = 3.14159265358979323d0 
      REAL*8, PARAMETER :: EIGHTDPI = 8.d0 / ONEPI

      ! Gas constant [erg/K/mole]
      REAL*8, PARAMETER :: RGAS     = BOLTG * AVG
      
      ! Universal gas constant [g/cm2/s2/mole/K]
      REAL*8, PARAMETER :: RSTARG   = 8.31450d+07

      ! Seconds per day
      REAL*8, PARAMETER :: SCDAY    = 86400.0d0

      ! Molecular weight of air 
      REAL*8, PARAMETER :: WTAIR    = 28.966d0
     
      !---------------------------------------------------------------
      ! Small number tolerances
      ! (now make these PARAMETERS instead of COMMON block variables)
      !---------------------------------------------------------------
      
      REAL*8, PARAMETER :: SMAL1    = 1d-06 
      REAL*8, PARAMETER :: SMAL2    = 1.0d-99
      REAL*8, PARAMETER :: SMAL3    = 1d-50


      INTEGER ::         NLAT,      NLONG,     NLAYER,    NVERT
      INTEGER ::         NLOOP,     NTLOOP,    KULOOP,    NTSPECGAS
      ! FP added NREAD to indicate number of entries to read
      ! from globchem.dat
      ! NREAD may be 20 or 24 depending on globchem.dat
      ! NREAD is set in mglob.dat as part of CTLDIM (hotp 8/4/09)
      ! NREAD is used in readchem (hotp 7/31/09)
      INTEGER :: NREAD

      INTEGER ::         NMASBAL,   KSLOOP,    NTLOOPUSE, NPVERT
      INTEGER ::         NTTLOOP,   NIJLOOP

      ! /CTLLOOP2/ needs to be declared THREADPRIVATE (bmy, 7/16/03)
      INTEGER ::         KTLOOP,    JLOOPLO,   IFSUN
                         
      INTEGER ::         ICOORD,    IFPRAT,    INCVMIX,   IFSOLVE
      INTEGER ::         IFURBAN,   IFTROP,    IFSTRAT,   ISL
      INTEGER ::         IGLOBCHEM, ITESTGEAR, IFSIN,     IFBOX
                         
      INTEGER ::         IPRGASA,   IPRGASB,   IPRGASC,   IPRGASD
      INTEGER ::         IPRGASE,   IPRGASF,   IPRGASG,   IPRGASH
      INTEGER ::         IPRGASLO,  IPRGASHI,  NUMPRG,    IPGMTOT 
      INTEGER ::         IOXSEC,    IOSPEC,    IOREAC,    IPRTEMP
      INTEGER ::         IPRMANY,   IPREADER,  IPRMET,    IPMASBUD
      INTEGER ::         IFPR1,     IPONEND,   IPRATES,   IPRPRESS
      INTEGER ::         IUSRDUM,   IGRIDZ,    IPGASES,   INCXY
      INTEGER ::         INCXZ,     INCYZ,     IGRIDX,    IGRIDY
      INTEGER ::         LXOUT,     LYOUT,     LLOOP,     LLOOP2
      INTEGER ::         LZOUT

      REAL*8  ::         TINTERVAL,  CHEMINTV,  TIME,      OXYCONS
      REAL*8  ::         HMAXNIT,    FRACDEC
!
      INTEGER ::         IHOUR,       NCS,       NBLOCKS,   IRCHEM
      INTEGER ::         NCSGAS,      NCSURBAN,  NCSTROP,   NCSSTRAT
      INTEGER ::         NPHOTALL,    IFDID,     IFNEVER,   IFNONE   
      INTEGER ::         NCSALL,      NCSTRST

      ! /IXYGD2/ needs to be held THREADPRIVATE.  Also remove NSTEPS
      ! since this can be declared local w/in "smvgear.f" (bmy, 7/16/03)
      INTEGER ::         NCSP,        KBLK

      REAL*8  ::         HMIN,        PLOURB,    PLOTROP,   TSPMIDC     

      ! /DGEAR2/ needs to be held THREADPRIVATE (hamid, bmy, 7/16/03)
      REAL*8  ::         HMAX,        R1DELT,    DELT,      TIMREMAIN
      REAL*8  ::         XELAPS,      TOLD,      RDELT,     XELAPLAST
      REAL*8  ::         RMSERR

      ! /DGEAR3/ doesn't need to be held THREADPRIVATE (hamid, bmy, 7/16/03)
      REAL*8  ::         SUMAVGE,     SUMAVHI,   SUMRMSE,   SUMRMHI
      REAL*8  ::         TOTSTEP,     TOTIT,     TELAPS

      INTEGER ::         NSFTOT,      NPDTOT,    NSTTOT,    ISREORD
      INTEGER ::         IFREORD,     IFAILTOT,  LFAILTOT,  NFAILTOT  

      ! /IGEAR2/ has to be declared THREADPRIVATE (bmy, 7/16/03)
      INTEGER ::         NQQ,         NSUBFUN,   NPDERIV   
      INTEGER ::         NFAIL,       IFAIL,     LFAIL

      INTEGER ::         NPHOT,       NPRODLO,   NPRODHI,   MSTEP
      INTEGER ::         MAXORD,      MBETWEEN,  IC3H8,     IC2H6

      ! /CHEM2A/ has to be held THREADPRIVATE (bmy, 7/16/03)
      INTEGER ::         ISCHAN,      NFDH3,     NFDL2,     NFDH2
      INTEGER ::         NFDL1,       NFDH1,     NFDREP,    NFDREP1   
      INTEGER ::         NFDL0,       NALLR

      ! Split off from /CHEM2A/ (bmy, 7/28/03)
      INTEGER ::         NOCC
         
      INTEGER ::         NGAS,        NMREAC

      ! Added NNADDG to /CHEM4/ for DMS+OH+O2 rxn (bdf, bmy, 4/18/03)
      ! Add NKN2O5 to /CHEM4/ to flag N2O5 hydrolysis rxn (mje, bmy, 8/7/03)
      ! Add NKHO2 to /CHEM4/ to flag HO2 aerosol uptake (jaegle 02/26/09)
      ! Add NNADDF, NNADDH, NKHOROI and NKHOROJ to /CHEM4/ for HOC2H4O rxns
      ! (tmf, 3/6/09)
      ! Added NKGLYC and NKHAC NKMCO3 (hotp 7/31/09)
      ! Added NKBrNO3 to /CHEM4/ to flag BrNO3 hydrolysis rxn (jpp, 5/4/10)
      ! Added NKHOBr and NK1HBR to /CHEM4/ to flag the two HOBr + HBr
      !      pseudo reactions (in place of HBr + HOBr --het-->) (jpp, 3/22/11)
      ! Added NK2HOBr and NK2HBR to /CHEM4/ to flag the two HOBr + HBr
      !      pseudo reactions (in place of HBr + HOBr --ice-->) (jpp, 6/16/11)
      ! Add NKO3PHOTH2 to account for the H2 reaction of O1D (bhh, jmao,
      !  eam, 7/18/11)
      !Added NNRO2HO2(jmao, 11/17/2012)
      INTEGER            NNADD1,        NNADDA(ICS),  NNADDB(  ICS), 
     &                   NNADDC(ICS),   NNADDD(ICS),  NNADDK(  ICS), 
     &                   NNADDV(ICS),   NNADDZ,       
     &                   NKO3PHOT(ICS), NKO3PHOTH2(ICS),
     &                   NNADDG(ICS),   NEMIS( ICS),  NDRYDEP( ICS),
     &                   NNADDF(ICS),   NNADDH(ICS), 
     &                   NKHOROI(ICS),  NKHOROJ(ICS),
     &                   NKHNO4(ICS),   NKN2O5,       NKHO2, 
     &                   NKGLYC(ICS,2), NKHAC(ICS,2), NKMCO3(ICS,3),
     &                   NKBrNO3,       NKHOBr,       NK1HBr, 
     &                   NK2HOBr,       NK2HBr,       NNRO2HO2(ICS) 

      INTEGER ::         IH2O,        IOXYGEN,   MB1,      MB2

      ! Added for interannually-varying Methane (bnd, bmy, 7/1/03)
      INTEGER ::         ICH4

      ! Added for interannually-varying Methane (bnd, bmy, 7/1/03)
      REAL*8  ::         C3090S,      C0030S,    C0030N,   C3090N

      ! Added for tracking oxidation of ISOP by OH (dkh, bmy, 6/1/06)
      ! SOAupdate: Added for tracking oxidation of ISOP by NO3 (hotp 6/11/10)
      INTEGER ::         ILISOPOH, ILISOPNO3

      ! Added for tracking oxidation of aromatic RO2s (dkh, 10/06/06)  
      ! bug fix: change to INTEGER (dbm, dkh, 06/26/07)  
      ! SOAupdate: added ILNRO2H/N (hotp 7/21/10)
      INTEGER ::         ILBRO2H,    ILBRO2N,   ILTRO2H,  ILTRO2N
      INTEGER ::         ILXRO2H,    ILXRO2N,   ILNRO2H,  ILNRO2N

      ! Added option for zeroing CHBr3 and CH2Br2 inside calcrate,
      ! since zeroing their emissions doesn't prevent small
      ! conc's on the order of 0.01pptv (jpp, 4/15/10)
      INTEGER ::         ICHBr3, ICH2Br2

      INTEGER ::         IOUT,        KGLC,      KCPD,     IO93
!
! *********************************************************************
! *               SET REAL AND INTEGER ARRAY VARIABLES                *
! *********************************************************************
!
      INTEGER, ALLOCATABLE ::  JLOWVAR( :), KTLPVAR( :), 
     &                         JLOFIXED(:), JHIFIXED(:)

      INTEGER, ALLOCATABLE ::  JREORDER(:), LREORDER(:), 
     &                         ITWO(    :), NCSLOOP( :,:)

      ! Add NKSPECG for DMS+OH+O2 rxn (bdf, bmy, 4/18/03)
      ! Added NKSPECF, NKSPECH to /IDICS/ for C2H4 chemistry (tmf, 3/6/09) 
      ! Added NKSPECRO2HO2 for RO2+HO2 (jmao, 11/17/2012)
      INTEGER ::
     &  NMOTH(   ICS), NTSPEC( ICS), JPHOTRAT(ICS),
     &  ISGAINR( ICS), ISPORL( ICS), NOGAINE( ICS),   NOUSE(    ICS),
     &  NSPEC(   ICS), NTRATES(ICS), ISGAINE( ICS),   NTLOOPNCS(ICS),
     &  NSPCSOLV(ICS), ISCHANG(ICS), NRATES(  ICS),   NM3BOD(   ICS),
     &  ITWOR(   ICS), ITHRR(  ICS), INOREP(  ICS),   NRATCUR(  ICS),
     &  NSURFACE(ICS), NPRESM( ICS), NMAIR(   ICS),   NMO2(     ICS),
     &  NMN2(    ICS), NNEQ(   ICS), NARR(    ICS),   NABR(     ICS),
     &  NACR(    ICS), NABC(   ICS), NKSPECW( ICS),   NKSPECX(  ICS),
     &  NKSPECY( ICS), NKSPECZ(ICS), NKSPECV(MAXGL2,ICS),ISLOSSR(ICS),
     &  NKSPECA(  MAXGL3,ICS), NKSPECB(  MAXGL3,ICS),
     &  NKSPECC(MAXGL3,ICS),NKSPECD(MAXGL3,ICS),NKSPECK(MAXGL3,ICS),
     &  NKSPECG(MAXGL2,ICS),NKSPECRO2HO2(MAXGL3,ICS),
     &  NKSPECF(MAXGL3,ICS), NKSPECH(MAXGL3,ICS)

      ! re-define some nkspec* arrays for harvard chem mechanism (bdf)
      INTEGER ::        NOLOSP(ICP),    NGNFRAC(ICP),    NOLOSRAT(ICP), 
     &                  IARRAY(ICP),    NALLRAT(ICP),    KZTLO(   ICP), 
     &                  KZTHI( ICP),    IONER(  ICP),    NPLLO(   ICP),  
     &                  NPLHI( ICP),    NFRLO(  ICP),    NFRHI(   ICP),
     &                  NPDLO( ICP),    NPDHI(  ICP),    IZLO (   ICP), 
     &                  JZLO ( ICP),    JLLO(   ICP),    JGLO(    ICP), 
     &                  IRMCT( ICP)

      REAL*8  ::        ABTOL(6,ICS),   ABST2(ICS),  
     &                  ERRMAX(ICS),    HMAXUSE(ICP),    TIMEINTV(ICS)    

      REAL*8  ::        WTGAS(IGAS),    GQSUMINI(IGAS), 
     &                  BSUMCHEM(IGAS), GQSUM(   IGAS),  QBKGAS(IGAS)       

      REAL*8  ::        CPREV(MXGSAER), CMODEL(MXGSAER), APORL(MXGSAER)

      INTEGER ::        IFPRGAS(IGAS),  LGNUM(IGAS),     NGMIX(IGAS)

      REAL*8  ::        DEFPRAT(MXRATE,ICS)

      REAL*8  ::        ARRT(MXCOF),    BRRT( MXCOF),   
     &                  FCVT(MXCOF),    FCT1T(MXCOF),    FCT2T(MXCOF)

      INTEGER ::        KCRRT(MXCOF)

      INTEGER :: 
     &  NKARR(NMTRATE,ICS),  NKABR(NMTRATE,ICS), NKACR(NMTRATE,ICS), 
     &  NKABC(NMTRATE,ICS),  IRORD(NMTRATE,ICS)

      REAL*8  ::
     &  ARR(    NMTRATE, ICS),  BRR(    NMTRATE, ICS), 
     &  FCV(    NMTRATE, ICS),  FCTEMP1(NMTRATE, ICS),
     &  FCTEMP2(NMTRATE, ICS)  

      INTEGER ::
     &  IAPROD( NMTRATE,  ICS),  NOLOSRN( NMTRATE, ICS),
     &  NRUSE(  NMTRATE,  ICS),  NRREP(   NMTRATE, ICS),
     &  NPRODUC(NMTRATE,  ICS),  IALLOSN( MXRATE,  ICS),
     &  NCEQUAT(NMTRATE,  ICS),  NOLDFNEW(NMTRATE, ICS),
     &  NEWFOLD(NMTRATE*2,ICS),  NKONER(  NMTRATE, ICS),
     &  NKTWOR( NMTRATE,  ICS),  NKTHRR(  NMTRATE, ICS),
     &  KCRR(   NMTRATE,  ICS),  LSKIP(   MXRATE,  ICS),
     &  IRMC(   NMTRATE      ),  JPHOTNK( NMTRATE, ICS),
     &  IUSED(  MXRATE,   ICS)

      ! For COMPAQ, put IRMA, IRMB in /INMTRATE2/ (Q. Liang, bmy, 10/17/05)
      REAL*8  ::        
     & IRMA(   NMTRATE      ),  IRMB(    NMTRATE     )

      INTEGER :: NEWNK(MAXGL)

      REAL*8  :: FRACP(MAXGL, ICS)

      INTEGER ::
     &  NREACOTH(MAXGL2,ICS), LGASBINO(MAXGL2,ICS),
     &  NKNLOSP( MAXGL3,ICS), LOSINACP(MAXGL3,ICS),
     &  IGNFRAC( MAXGL, ICS), NKGNFRAC(MAXGL, ICS),
     &  NREACAIR(MAXGL3,ICS), NREAC3B( MAXGL3,ICS),
     &  NREACEQ( MAXGL3,ICS), NREQOTH( MAXGL3,ICS),
     &  NREACN2( MAXGL3,ICS), NREACO2( MAXGL3,ICS),
     &  NREACPM( MAXGL3,ICS), LGAS3BOD(MAXGL3,ICS),
     &  NKSURF(  MAXGL4    ), NCOATG(  MAXGL4    ) 

      INTEGER :: MBCOMP(IMASBAL,2), MBTRACE(IMASBAL)

      ! /DKBLOOP2/ needs to be declared THREADPRIVATE
      REAL*8 ::
     &  CNEW(   KBLOOP,  MXGSAER),  
     &  CEST(   KBLOOP,  MXGSAER),  
     &  GLOSS(  KBLOOP,  MXGSAER),  
     &  CHOLD(  KBLOOP,  MXGSAER), 
     &  VDIAG(  KBLOOP,  MXGSAER),  CBLK(  KBLOOP,MXGSAER),  
     &  DTLOS(  KBLOOP,  MXGSAER),  EXPLIC(KBLOOP,MXGSAER),
     &  CONC(   KBLOOP,MXGSAER*7),  
     &  RRATE(  KBLOOP,  NMTRATE), 
      !***************KPP_INTERFACE****************
     &  RRATE_FOR_KPP(  KBLOOP,  NMTRATE), 
      !********************************************  
     &  URATE(  KBLOOP,NMTRATE,3), 
     &  TRATE(  KBLOOP,NMTRATE*2), 
     &  CORIG(  KBLOOP,  MXGSAER)

      ! /DKBLOOP0/ needs to be held THREADPRIVATE
      REAL*8  :: CC2(KBLOOP,0:MXARRAY)  

      INTEGER, ALLOCATABLE :: MLOP(:,:), JLOP_SMV(:,:,:)

      INTEGER :: 
     &  NKPHOTRAT(IPHOT,ICS),      NPPHOTRAT(IPHOT,ICS),  
     &  NKNPHOTRT(IPHOT,ICS) 

      REAL*8  :: FRACGAIN(MXGSAER,ICS), QBKCHEM( MXGSAER,ICS) 

      INTEGER ::
     &  NUMLOST( MXGSAER,    ICS), NUMGFRT( MXGSAER,        ICS),  
     &  NUMLOSS( MXGSAER,    ICP), JPORL(   MXGSAER,MAXGL,  ICS),
     &  NUMGAINT(MXGSAER,    ICS), NGAINE(  MXGSAER,        ICS),  
     &  NUMGAIN( MXGSAER,    ICP), IGAINR(  MXGSAER,        ICS), 
     &  IPORL(   MXGSAER,    ICS), IGAINE(  MXGSAER,        ICS),  
     &  ISOLVSPC(MXGSAER,    ICS), INEWOLD( MXGSAER,        ICS),
     &  MAPPL(   MXGSAER,    ICS), ISAPORL( MXGSAER            ),
     &  NUMPORL( MXGSAER,    ICP), ISPARDER(MXGSAER,MXGSAER    ),
     &  JLOSST(  MXGSAER,MAXGL,ICS)

      INTEGER :: JZILCH(MXGSAER),  KZILCH(MXGSAER),  MZILCH(MXGSAER) 

      INTEGER ::
     &  LZERO(   MXGSAER,MXGSAER), JARRAYPT(MXGSAER,MXGSAER),
     &  IZILCH(  MXGSAER,MXGSAER), JARRDIAG(MXGSAER,   ICP),
     &  JLOZ1(   MXGSAER,    ICP), JHIZ1(   MXGSAER,   ICP), 
     &  IJTLO(   MXGSAER,    ICP), IJTHI(   MXGSAER,   ICP),
     &  IMZTOT(  MXGSAER,    ICP), IFREPRO( MXGSAER,MXRATE, ICS),
     &  IZLO1( MXGSAER,ICP),
     &  IZLO2( MXGSAER,ICP), IZHI0( MXGSAER,ICP), IZHI1( MXGSAER,ICP)

      REAL*8 :: FRACNFR(MXCOUNT4), FRACPL(MXCOUNT2) 

      INTEGER :: 
     &  JZERO( MXCOUNT3),  KZERO(  MXCOUNT3),  MZERO(   MXCOUNT3),
     &  IZEROK(MXCOUNT2),  JZEROA( MXCOUNT3),
     &  IKDECA(MXCOUNT3),  KJDECA( MXCOUNT3),  LOSSRA(  MXCOUNT4),  
     &  IKDECB(MXCOUNT3),  KJDECB( MXCOUNT3),  LOSSRB(  MXCOUNT4), 
     &  IKDECC(MXCOUNT3),  KJDECC( MXCOUNT3),  LOSSRC(  MXCOUNT4), 
     &  IKDECD(MXCOUNT3),  KJDECD( MXCOUNT3),  LOSSRD(  MXCOUNT4), 
     &  IKDECE(MXCOUNT3),  KJDECE( MXCOUNT3),  LOSSRE(  MXCOUNT4),
     &  KZEROA(MXCOUNT4),  MZEROA( MXCOUNT4), 
     &  KZEROB(MXCOUNT4),  MZEROB( MXCOUNT4), 
     &  KZEROC(MXCOUNT4),  MZEROC( MXCOUNT4), 
     &  KZEROD(MXCOUNT4),  MZEROD( MXCOUNT4), 
     &  KZEROE(MXCOUNT4),  MZEROE( MXCOUNT4), 
     &  IPOSPD(MXCOUNT2),  IIALPD( MXCOUNT2),  NKPDTERM(MXCOUNT2),  
     &  IJVAL( MXCOUNT3),  IKZTOT( MXCOUNT4),  JSPNPL(  MXCOUNT4),    
     &  NKNFR( MXCOUNT4),  JSPCNFR(MXCOUNT4) 

      INTEGER :: 
     &  IDH5(  MXCOUNT3),  IDH4(  MXCOUNT3),  IDH3(  MXCOUNT3),
     &  IDH2(  MXCOUNT3),  IDH1(  MXCOUNT3),  IDL5(  MXCOUNT3),
     &  IDL4(  MXCOUNT3),  IDL3(  MXCOUNT3),  IDL2(  MXCOUNT3),
     &  IDL1(  MXCOUNT3),  
     &  KBH5(  MXCOUNT4),  KBH4(  MXCOUNT4),  KBH3(  MXCOUNT4),
     &  KBH2(  MXCOUNT4),  KBH1(  MXCOUNT4),  KBL5(  MXCOUNT4),
     &  KBL4(  MXCOUNT4),  KBL3(  MXCOUNT4),  KBL2(  MXCOUNT4),
     &  KBL1(  MXCOUNT4),  
     &  MBH5(  MXCOUNT4),  MBH4(  MXCOUNT4),  MBH3(  MXCOUNT4),
     &  MBH2(  MXCOUNT4),  MBH1(  MXCOUNT4),  MBL5(  MXCOUNT4),
     &  MBL4(  MXCOUNT4),  MBL3(  MXCOUNT4),  MBL2(  MXCOUNT4),
     &  MBL1(  MXCOUNT4),  
     &  NPH5(  MXCOUNT4),  NPH4(  MXCOUNT4),  NPH3(  MXCOUNT4),
     &  NPH2(  MXCOUNT4),  NPH1(  MXCOUNT4),  NPL5(  MXCOUNT4),
     &  NPL4(  MXCOUNT4),  NPL3(  MXCOUNT4),  NPL2(  MXCOUNT4),
     &  NPL1(  MXCOUNT4)  

      REAL*8  :: GEARCONC(MXGSAER,0:MXHOLD,ICS)

      REAL*8  :: WTMB(IMASBAL,MXGSAER,2)

      INTEGER :: JMBCOMP(IMASBAL,MXGSAER,2) 

      REAL*8     FKOEF(NMRPROD,NMTRATE,ICS), 
     &           FK2(  NMRPROD,NMTRATE,ICS)  

      INTEGER :: IRM( NMRPROD,NMTRATE,ICS), 
     &           IRM2(NMRPROD,NMTRATE,ICS) 

      REAL*8 ::  ASET(10,8), PINP(IMISC), CVAR(15), O3DOBS(12)  

      REAL*8 ::  ENQQ2(MORDER), ENQQ3( MORDER), CONPST(MORDER),
     &           ENQQ1(MORDER), CONP15(MORDER)

      REAL*8 ::  PERTS2(MORDER,3), PERTST(MORDER,3)

      INTEGER, ALLOCATABLE :: JLLOW(:), KLHI(:)

      REAL*8 ::  XGDFCF, ASTKCF, RUARSL, RH100

      INTEGER :: IARSFA, MWARSL, MNTHARS

      INTEGER :: NKEMIS(MAXGL3,ICS),  NTEMIS(MAXGL3,ICS),
     &           NKDRY (MAXGL3,ICS),  NTDEP( MAXGL3)

      !=================================================================
      ! Common blocks for ND65 diagnostic (ljm, bmy, 7/8/03)
      !=================================================================
      INTEGER     ::     IFAM(MAXFAM), NFAMILIES

      CHARACTER*4 ::     PORL(MAXFAM)

      LOGICAL     ::     LFAMILY, ITS_NOT_A_ND65_FAMILY(IGAS)

      !=================================================================
      ! Declare some common blocks THREADPRIVATE for the OpenMP
      ! parallelization (bdf, bmy, 4/1/03)
      !=================================================================
!$OMP THREADPRIVATE( ISCHAN, NFDH3, NFDL2, NFDH2, NFDL1, NFDH1 )
!$OMP THREADPRIVATE( NFDREP, NFDREP1, NFDL0, NALLR )      
!$OMP THREADPRIVATE( KTLOOP,    JLOOPLO,   IFSUN )
!$OMP THREADPRIVATE( HMAX, R1DELT, DELT, TIMREMAIN, XELAPS )
!$OMP THREADPRIVATE( TOLD, RDELT, XELAPLAST, RMSERR )
!$OMP THREADPRIVATE( CC2 )
!$OMP THREADPRIVATE( CNEW, CEST, GLOSS, CHOLD, VDIAG, CBLK, DTLOS )
!$OMP THREADPRIVATE( EXPLIC, CONC, RRATE, RRATE_FOR_KPP )
!$OMP THREADPRIVATE( URATE, TRATE, CORIG )
!$OMP THREADPRIVATE( NQQ, NSUBFUN, NPDERIV, NFAIL, IFAIL, LFAIL )
!$OMP THREADPRIVATE( NCSP, KBLK )

#if   defined( COMPAQ )
! For COMPAQ, declare /INMTRATE2/ threadprivate (Q. Liang, bmy, 10/17/05)
!$OMP THREADPRIVATE( IRMA,IRMB )
#endif
!

! *********************************************************************
! ****************** END OF COMMON BLOCK COMODE.H *********************
! *********************************************************************
!
!BOC
      CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_comode_loop
!
! !DESCRIPTION: Subroutine INIT\_COMODE\_LOOP initializes size parameters
!  with the geospatial values obtained from the ESMF interface.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE Init_COMODE_LOOP( am_I_Root, Input_Opt, RC )
!
! !USES:
!
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
!
! !INPUT PARAMETERS: 
!
      LOGICAL,        INTENT(IN)  :: am_I_Root   ! Are we on the root CPU?
      TYPE(OptInput), INTENT(IN)  :: Input_Opt   ! Input Options object
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT) :: RC          ! Success or failure?
!
! !REMARKS:
!  Need to add error-checking on the allocation statements, so that we
!  exit the code upon error.
! 
! !REVISION HISTORY: 
!  17 Oct 2012 - R. Yantosca - Need to also set NLAT, NLONG, etc.
!  17 Oct 2012 - R. Yantosca - Added ProTeX headers
!  22 Oct 2012 - R. Yantosca - Now references gigc_errcode_mod.F90
!  27 Nov 2012 - R. Yantosca - Now pass Input_Opt via the argument list
!EOP
!------------------------------------------------------------------------------
!BOC
      ! Assume success
      RC     = GIGC_SUCCESS

      ! Set size variables (from comode_loop_mod.F)
      ILAT    = JJPAR
      ILONG   = IIPAR
      IVERT   = LLPAR
      IPVERT  = IVERT + IPLUME
      ILAYER  = LLPAR + 1
      ITLOOP  = ILAT * ILONG * IPVERT
      IMLOOP  = ILAT * ILONG
      ILTLOOP = IMLOOP * ILAYER
      MXBLOCK = 16 + ITLOOP/KBLOOP

      ! Set size variables (these declarations were in chemdr.F)
      NLAT    = JJPAR
      NLONG   = IIPAR
      NVERT   = IVERT 
      NPVERT  = NVERT
      NPVERT  = NVERT + IPLUME

      ! Allocate arrays
      ALLOCATE( MLOP( ILAT, ILONG ), STAT=RC )
      IF ( RC /= GIGC_SUCCESS ) THEN
         WRITE( 6, '(a)' ) 'Allocation error: MLOP'
         RETURN
      ENDIF

      ALLOCATE( JLOP_SMV( ILAT, ILONG, ILAYER ), STAT=RC )
      IF ( RC /= GIGC_SUCCESS ) THEN
         WRITE( 6, '(a)' ) 'Allocation error: JLOP_SMV'
         RETURN
      ENDIF

      ALLOCATE( JLLOW( MXBLOCK ), STAT=RC )
      IF ( RC /= GIGC_SUCCESS ) THEN
         WRITE( 6, '(a)' ) 'Allocation error: JLLOW'
         RETURN
      ENDIF

      ALLOCATE( KLHI( MXBLOCK ), STAT=RC )
      IF ( RC /= GIGC_SUCCESS ) THEN
         WRITE( 6, '(a)' ) 'Allocation error: KLHI'
         RETURN
      ENDIF

      ALLOCATE( JLOWVAR( MXBLOCK ), STAT=RC )
      IF ( RC /= GIGC_SUCCESS ) THEN
         WRITE( 6, '(a)' ) 'Allocation error: JLOWVAR'
         RETURN
      ENDIF

      ALLOCATE( KTLPVAR( MXBLOCK ), STAT=RC )
      IF ( RC /= GIGC_SUCCESS ) THEN
         WRITE( 6, '(a)' ) 'Allocation error: KLTPVAR'
         RETURN
      ENDIF

      ALLOCATE( JLOFIXED( MXBLOCK ), STAT=RC )
      IF ( RC /= GIGC_SUCCESS ) THEN
         WRITE( 6, '(a)' ) 'Allocation error: JLOFIXED'
         RETURN
      ENDIF

      ALLOCATE( JHIFIXED( MXBLOCK ), STAT=RC )
      IF ( RC /= GIGC_SUCCESS ) THEN
         WRITE( 6, '(a)' ) 'Allocation error: JHIFIXED'
         RETURN
      ENDIF

      ALLOCATE( JREORDER( ITLOOP ), STAT=RC )
      IF ( RC /= GIGC_SUCCESS ) THEN
         WRITE( 6, '(a)' ) 'Allocation error: JREORDER'
         RETURN
      ENDIF

      ALLOCATE( LREORDER( ITLOOP ), STAT=RC )
      IF ( RC /= GIGC_SUCCESS ) THEN
         WRITE( 6, '(a)' ) 'Allocation error: LREORDERP'
         RETURN
      ENDIF

      ALLOCATE( ITWO( ITLOOP ), STAT=RC )
      IF ( RC /= GIGC_SUCCESS ) THEN
         WRITE( 6, '(a)' ) 'Allocation error: ITWO'
         RETURN
      ENDIF

      ALLOCATE( NCSLOOP( ITLOOP, ICS ), STAT=RC )
      IF ( RC /= GIGC_SUCCESS ) THEN
         WRITE( 6, '(a)' ) 'Allocation error: ITLOOP'
         RETURN
      ENDIF

      IF ( Input_Opt%LPRT .and. am_I_Root ) THEN
         write( 6, * ) '###########################################'
         write( 6, * ) '### END of Init_COMODE_LOOP_MOD'
         write( 6, * ) '### IIPAR   : ', IIPAR
         write( 6, * ) '### JJPAR   : ', JJPAR
         write( 6, * ) '### LLPAR   : ', LLPAR
         write( 6, * ) '### LLTROP  : ', LLTROP
         write( 6, * ) '### ----------'
         write( 6, * ) '### ILAT    : ', ILAT   
         write( 6, * ) '### ILONG   : ', ILONG  
         write( 6, * ) '### IVERT   : ', IVERT  
         write( 6, * ) '### IPVERT  : ', IPVERT 
         write( 6, * ) '### ILAYER  : ', ILAYER 
         write( 6, * ) '### ITLOOP  : ', ITLOOP 
         write( 6, * ) '### IMLOOP  : ', IMLOOP 
         write( 6, * ) '### ILTLOOP : ', ILTLOOP
         write( 6, * ) '### MXBLOCK : ', MXBLOCK
         write( 6, * ) '### NLAT    : ', NLAT   
         write( 6, * ) '### NLONG   : ', NLONG  
         write( 6, * ) '### NVERT   : ', NVERT  
         write( 6, * ) '### NPVERT  : ', NPVERT 
         write( 6, * ) '###########################################'
      ENDIF

      END SUBROUTINE Init_COMODE_LOOP
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cleanup_comode_loop
!
! !DESCRIPTION: Subroutine CLEANUP\_COMODE\_LOOP deallocates all module arrays.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE Cleanup_COMODE_LOOP( am_I_Root, RC )
!
! !USES:
!
      USE GIGC_ErrCode_Mod

      IMPLICIT NONE
!
! !INPUT PARAMETERS: 
!
      LOGICAL, INTENT(IN)  :: am_I_Root   ! Are we on the root CPU?
!
! !OUTPUT PARAMETERS:
!
      INTEGER, INTENT(OUT) :: RC          ! Success or failure?
!
! !REMARKS:
!
! !REVISION HISTORY: 
!  19 Nov 2012 - R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
      ! Assume success
      RC = GIGC_SUCCESS

      ! Deallocate arrays
      IF ( ALLOCATED( MLOP     ) ) DEALLOCATE( MLOP     )  
      IF ( ALLOCATED( JLOP_SMV ) ) DEALLOCATE( JLOP_SMV )  
      IF ( ALLOCATED( JLLOW    ) ) DEALLOCATE( JLLOW    ) 
      IF ( ALLOCATED( KLHI     ) ) DEALLOCATE( KLHI     ) 
      IF ( ALLOCATED( JLOWVAR  ) ) DEALLOCATE( JLOWVAR  ) 
      IF ( ALLOCATED( KTLPVAR  ) ) DEALLOCATE( KTLPVAR  ) 
      IF ( ALLOCATED( JLOFIXED ) ) DEALLOCATE( JLOFIXED ) 
      IF ( ALLOCATED( JHIFIXED ) ) DEALLOCATE( JHIFIXED ) 
      IF ( ALLOCATED( JREORDER ) ) DEALLOCATE( JREORDER ) 
      IF ( ALLOCATED( LREORDER ) ) DEALLOCATE( LREORDER ) 
      IF ( ALLOCATED( ITWO     ) ) DEALLOCATE( ITWO     )
      IF ( ALLOCATED( NCSLOOP  ) ) DEALLOCATE( NCSLOOP  )

      END SUBROUTINE Cleanup_COMODE_LOOP
!EOC
      END MODULE COMODE_LOOP_MOD
