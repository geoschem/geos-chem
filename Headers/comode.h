! $Id: comode.h,v 1.2 2010/03/15 19:33:18 ccarouge Exp $
!
!******************************************************************************
!  Header file COMODE contains common blocks and variables for SMVGEAR II.
!  (M. Jacobson 1997; bdf, bmy, 4/23/03, 6/1/06)
!
!  NOTES:
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
!  28 Jun 2011 - M. Payer    - Add modifications for SOA + semivol POA (H. Pye)
!
!******************************************************************************
!
C         CCCCCCC  OOOOOOO  M     M  OOOOOOO  DDDDDD   EEEEEEE 
C         C        O     O  M M M M  O     O  D     D  E  
C         C        O     O  M  M  M  O     O  D     D  EEEEEEE  
C         C        O     O  M     M  O     O  D     D  E  
C         CCCCCCC  OOOOOOO  M     M  OOOOOOO  DDDDDD   EEEEEEE 
C
C *********************************************************************
C * THIS IS THE COMMON BLOCK FOR "SMVGEAR" AND "MIE," TWO ORDINARY    *
C * DIFFERENTIAL EQUATION SOLVERS. THE REFERENCE FOR THE CODES IS     *
C *                                                                   * 
C *   JACOBSON M. Z. AND TURCO R. P. (1993) SMVGEAER: A SPARSE-       * 
C *    MATRIX, VECTORIZED GEAR CODE FOR ATMOSPHERIC MODELS.           *
C *    SUBMITTED TO ATMOSPHERIC ENVIRONMENT, PART A. MAY 20, 1993     * 
C *                                                                   * 
C * COMODE.H SETS PARAMETER VALUES AND SERVES AS A COMMON BLOCK FOR   *
C * ALL DIMENSIONED AND NON-DIMENSIONED VARIABLES. COMODE.H ALSO      *
C * DEFINES EACH PARAMETER, BUT DATA FILE DEFINE.DAT EXPLAINS NON-    *
C * DIMENSIONED VARIABLES. INDIVIDUAL SUBROUTINES DEFINE DIMENSIONED  *
C * VARIABLES.                                                        *
C *********************************************************************
C
C *********************************************************************
C *                       SET PARAMETERS                              *
C *********************************************************************
C
C ****************** COORDINATE-SYSTEM PARAMETERS *********************
C
C ILAT     = MAXIMUM NUMBER OF LATITUDE(ILAT)   GRID POINTS
C ILONG    = MAXIMUM NUMBER OF LONGITUDE(ILONG) GRID POINTS
C IMLOOP   = ILAT * ILONG - USED FOR MORE EFFICIENT ARRAYS
C IVERT    = MAXIMUM NUMBER OF LAYERS 
C ILAYER   = MAXIMUM OF LAYER BOUNDARIES
C KBLOOP   = MAXIMUM NUMBER OF GRID POINTS IN A VECTORIZED BLOCK 
C            SHOULD RANGE FROM 512 (BELOW WHICH VECTORIZATION DECREASES)
C            TO 1024 (ABOVE WHICH, ARRAY SPACE IS LIMITED)
C MXBLOCK  = MAXIMUM NUMBER OF GRID POINT BLOCKS
C MAXDAYS  = MAXIMUM NUMBER OF DAYS FOR THE MODEL TO RUN
C

      INTEGER ILAT,ILONG,IVERT,IPLUME,IPVERT,ITLOOP,KBLOOP,MXBLOCK
      INTEGER IMLOOP,ILAYER,ILTLOOP,MAXDAYS
      PARAMETER (ILAT    = JJPAR                   )  
      PARAMETER (ILONG   = IIPAR                   )  

      ! LLTROP is the max number of tropospheric levels 
      PARAMETER (IVERT   = LLTROP                  )  

      ! GEOS-CHEM does not use plumes...set IPLUME=0 
      PARAMETER (IPLUME  = 0                       )
      PARAMETER (IPVERT  = IVERT + IPLUME          )
      PARAMETER (ITLOOP  = ILAT * ILONG * IPVERT   )

      ! Regular
      PARAMETER (KBLOOP  = 24                       )
      PARAMETER (IMLOOP  = ILAT * ILONG            )
      PARAMETER (ILAYER  = IVERT + 1               )
      PARAMETER (ILTLOOP = IMLOOP * ILAYER         )
      PARAMETER (MAXDAYS = 1000                    )
      PARAMETER (MXBLOCK = 16 + ITLOOP/KBLOOP      )
C
C ************************* TRACER PARAMETERS ****************************
C IDEMS    = EMISSION ID NUMBER (WHICH SPECIES)
      INTEGER IDEMS
      ! NEMPARA = max no. of anthropogenic emissions
      ! NEMPARB = max no. of biogenic emissions
      COMMON /JTRCID/ IDEMS(NEMPARA+NEMPARB)
C
C ************************* GAS-PHASE PARAMETERS **********************
C
C IGAS    = MAXIMUM NUMBER OF GASES, ACTIVE + INACTIVE 
C IAERTY  = MAXIMUM NUMBER OF AQUEOUS CHEMISTRY SPECIES (SET = 1
C           HERE SINCE NO AQUEOUS CHEMISTRY IS INCLUDED)
C NMRATE  = MAXIMUM NUMBER OF RATES CONSTANTS (MAX # REACTIONS)
C IPHOT   = MAXIMUM NUMBER OF PHOTO RATES
C NMTRATE = MAXIMUM NUMBER OF KINETIC + PHOTO REACTIONS 
C NMQRATE = MAXIMUM NUMBER OF AQUEOUS CHEMICAL REACTIONS PLUS PHOTO
C           PROCESSES (SET = 1 HERE SINCE NO AQUEOUS CHEMISTRY INCLUDED)
C NMRPROD = MAXIMUM NUMBER OF SPECIES IN A REACTION RATE 
C NMDEAD  = MAXIMUM NUMBER OF DEAD SPECIES
C MAXGL   = MAXIMUM NUMBER OF GAINS / LOSSES FOR EACH CHEMICAL SPECIES.
C MAXGL2  = AN ARRAY DIMENSION SMALLER THAN MAXGL
C MAXGL3  = AN ARRAY DIMENSION SMALLER THAN MAXGL2
C ICS     = NUMBER OF TYPES OF CHEMISTRY: up to 3 for gas phase
C ICP     = TYPES OF CHEMISTRY * 2 (ONE FOR DAY, ONE FOR NIGHT)
C MORDER  = MAXIMUM ORDER FOR GEAR PARAMETERS FOR DIMENSION PURPOSES
C
      INTEGER IGAS,IAERTY,NMRATE,IPHOT,NMTRATE,NMQRATE,NMRPROD,NMDEAD
      INTEGER MAXGL,MAXGL2,MAXGL3,MAXGL4,ICS,ICP,MORDER,IPHOT8,IMISC
      INTEGER IMASBAL,IALTS,MXCOF

      ! Updated for SMVGEAR II (bdf, bmy, 4/1/03)
!FP_ISOP 06/15/09 !IGAS, NMRATE and MAXGL increased
!      PARAMETER ( IGAS    = 200,               IAERTY  = 1           )
!      PARAMETER ( NMRATE  = 510,               IPHOT   = 85          )
      PARAMETER ( IGAS    = 300,               IAERTY  = 1           )
      PARAMETER ( NMRATE  = 700,               IPHOT   = 85          )
      PARAMETER ( NMTRATE = NMRATE + IPHOT,    NMQRATE = 1           ) 
      PARAMETER ( NMRPROD = 25,                NMDEAD  = 100         )
!      PARAMETER ( MAXGL   = 750,               MAXGL2  = 90          )
      PARAMETER ( MAXGL   = 1000,              MAXGL2  = 90          )
      PARAMETER ( MAXGL3  = NNPAR,             MAXGL4  = 10          )
      PARAMETER ( ICS     = 3,                 ICP     = 2*ICS       ) 
      PARAMETER ( MORDER  = 7                                        ) 
      PARAMETER ( IPHOT8  = IPHOT + 8,         IMISC   = 100         )
      PARAMETER ( IMASBAL = 9,                 IALTS   = 22          )
      PARAMETER ( MXCOF   = 5                                        )
C
C ****************** PARAMETERS TO MINIMIZE ARRAY SPACE ***************
C
C MXGSAER    = LARGER OF IGAS, IAERTY 
C MXRATE     = LARGER OF NMTRATE, NMQRATE 
C MXCC2      = LARGER OF MXGSAER, MXARRAY  
C MXCOUNT1.. = ARRAYS SIZES USED TO MINIMIZE MATRIX SPACE
C MXARRAY    = MAXIMUM ONE-DIMENSIONAL ARRAY-LENGTH OF SPARSE MATRIX 
C
      INTEGER MXGSAER,MXRATE,MXARRAY,MXCC2,MXCOUNT1,MXCOUNT2,MXCOUNT3,
     1        MXCOUNT4,MXHOLD
      PARAMETER( MXGSAER  = IGAS)  
      PARAMETER( MXRATE   = NMTRATE)  
      PARAMETER( MXARRAY  = 3000)
      PARAMETER( MXCC2    = MXARRAY)  
      PARAMETER( MXCOUNT1 = MXGSAER * MAXGL3 * 3)  
      PARAMETER( MXCOUNT2 = MXGSAER * MAXGL3 * 7)  
      PARAMETER( MXCOUNT3 = MXGSAER * 50)
      PARAMETER( MXCOUNT4 = MXGSAER * 20)
      PARAMETER( MXHOLD   = 250)
C
C
C *********************************************************************
C *                    SET CHARACTER LENGTHS                          *
C *********************************************************************
C
      CHARACTER*14 NAMESPEC, NAMD,
     1             APGASA,   APGASB,   APGASC,   APGASD,   APGASE,
     2             APGASF,   APGASG,   APGASH,   IFSORM,     
     3             XINP,     RINP
      CHARACTER*14 NAMEGAS,   NAMEMB,
     1             JST, NAMENCS, ACORNER, SINP, NAMEPHOT, CHEMTYP
      CHARACTER*4  DINP,      DINPLAST
      CHARACTER*2  SPECIAL,   SPECL
C
      CHARACTER*80 HEADING, COMMENT 
C
C *********************************************************************
C *                 SET CHARACTER DIMENSIONS                          *
C *********************************************************************
C
       COMMON / CHARAC /
     1  APGASA,   APGASB,    APGASC,   APGASD,   APGASE,
     2  APGASF,   APGASG,    APGASH,   IFSORM,
     3  DINP,     HEADING,   COMMENT,
     4  JST,      ACORNER,   SPECIAL,  DINPLAST

       COMMON / CHAR2 /
     1  XINP(      IMISC), RINP(       IMISC), SINP(       IMISC),
     2  NAMEMB(  IMASBAL), CHEMTYP(      ICS), NAMD(      NMDEAD),
     3  SPECL(     MXCOF) 

       COMMON / CHAR3 /
     1  NAMESPEC(0:IGAS,ICS), NAMEGAS(0:IGAS), NAMENCS(0:MXGSAER,ICS)

       COMMON / CHAR4 /
     2  NAMEPHOT(NMRPROD,IPHOT)
C
C *********************************************************************
C *           SET REAL AND INTEGER NON-ARRAY VARIABLES                *
C *********************************************************************
C
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
      COMMON /CTLLOOP/   NLAT,      NLONG,     NLAYER,    NVERT,     
     &                   NLOOP,     NTLOOP,    KULOOP,    NTSPECGAS, 
     &                   NMASBAL,   KSLOOP,    NTLOOPUSE, NPVERT,    
     &                   NTTLOOP,   NIJLOOP,   NREAD

      ! /CTLLOOP2/ needs to be declared THREADPRIVATE (bmy, 7/16/03)
      INTEGER ::         KTLOOP,    JLOOPLO,   IFSUN
      COMMON /CTLLOOP2/  KTLOOP,    JLOOPLO,   IFSUN
                         
      INTEGER ::         ICOORD,    IFPRAT,    INCVMIX,   IFSOLVE
      INTEGER ::         IFURBAN,   IFTROP,    IFSTRAT,   ISL
      INTEGER ::         IGLOBCHEM, ITESTGEAR, IFSIN,     IFBOX
      COMMON /CTLPROC/   ICOORD,    IFPRAT,    INCVMIX,   IFSOLVE,
     &                   IFURBAN,   IFTROP,    IFSTRAT,   ISL,
     &                   IGLOBCHEM, ITESTGEAR, IFSIN,     IFBOX
                         
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
      COMMON /CTLPRNT/   IPRGASA,   IPRGASB,   IPRGASC,   IPRGASD,   
     &                   IPRGASE,   IPRGASF,   IPRGASG,   IPRGASH,   
     &                   IPRGASLO,  IPRGASHI,  NUMPRG,    IPGMTOT,   
     &                   IOXSEC,    IOSPEC,    IOREAC,    IPRTEMP,   
     &                   IPRMANY,   IPREADER,  IPRMET,    IPMASBUD,  
     &                   IFPR1,     IPONEND,   IPRATES,   IPRPRESS,  
     &                   IUSRDUM,   IGRIDZ,    IPGASES,   INCXY,     
     &                   INCXZ,     INCYZ,     IGRIDX,    IGRIDY,
     &                   LXOUT,     LYOUT,     LLOOP,     LLOOP2,    
     &                   LZOUT

      REAL*8  ::         TINTERVAL,  CHEMINTV,  TIME,      OXYCONS
      REAL*8  ::         HMAXNIT,    FRACDEC
      COMMON /XYGRID/    TINTERVAL,  CHEMINTV,  TIME,      OXYCONS,    
     &                   HMAXNIT,    FRACDEC
C
      INTEGER ::         IHOUR,       NCS,       NBLOCKS,   IRCHEM
      INTEGER ::         NCSGAS,      NCSURBAN,  NCSTROP,   NCSSTRAT
      INTEGER ::         NPHOTALL,    IFDID,     IFNEVER,   IFNONE   
      INTEGER ::         NCSALL,      NCSTRST
      COMMON /IXYGD/     IHOUR,       NCS,       NBLOCKS,   IRCHEM,   
     &                   NCSGAS,      NCSURBAN,  NCSTROP,   NCSSTRAT, 
     &                   NPHOTALL,    IFDID,     IFNEVER,   IFNONE,    
     &                   NCSALL,      NCSTRST

      ! /IXYGD2/ needs to be held THREADPRIVATE.  Also remove NSTEPS
      ! since this can be declared local w/in "smvgear.f" (bmy, 7/16/03)
      INTEGER ::         NCSP,        KBLK
      COMMON /IXYGD2/    NCSP,        KBLK

      REAL*8  ::         HMIN,        PLOURB,    PLOTROP,   TSPMIDC     
      COMMON /DGEAR/     HMIN,        PLOURB,    PLOTROP,   TSPMIDC

      ! /DGEAR2/ needs to be held THREADPRIVATE (hamid, bmy, 7/16/03)
      REAL*8  ::         HMAX,        R1DELT,    DELT,      TIMREMAIN
      REAL*8  ::         XELAPS,      TOLD,      RDELT,     XELAPLAST
      REAL*8  ::         RMSERR
      COMMON /DGEAR2/    HMAX,        R1DELT,    DELT,      TIMREMAIN, 
     &                   XELAPS,      TOLD,      RDELT,     XELAPLAST, 
     &                   RMSERR

      ! /DGEAR3/ doesn't need to be held THREADPRIVATE (hamid, bmy, 7/16/03)
      REAL*8  ::         SUMAVGE,     SUMAVHI,   SUMRMSE,   SUMRMHI
      REAL*8  ::         TOTSTEP,     TOTIT,     TELAPS
      COMMON /DGEAR3/    SUMAVGE,     SUMAVHI,   SUMRMSE,   SUMRMHI, 
     &                   TOTSTEP,     TOTIT,     TELAPS

      INTEGER ::         NSFTOT,      NPDTOT,    NSTTOT,    ISREORD
      INTEGER ::         IFREORD,     IFAILTOT,  LFAILTOT,  NFAILTOT  
      COMMON /IGEAR/     NSFTOT,      NPDTOT,    NSTTOT,    ISREORD,   
     &                   IFREORD,     IFAILTOT,  LFAILTOT,  NFAILTOT

      ! /IGEAR2/ has to be declared THREADPRIVATE (bmy, 7/16/03)
      INTEGER ::         NQQ,         NSUBFUN,   NPDERIV   
      INTEGER ::         NFAIL,       IFAIL,     LFAIL
      COMMON /IGEAR2/    NQQ,         NSUBFUN,   NPDERIV, 
     &                   NFAIL,       IFAIL,     LFAIL

      INTEGER ::         NPHOT,       NPRODLO,   NPRODHI,   MSTEP
      INTEGER ::         MAXORD,      MBETWEEN,  IC3H8,     IC2H6
      COMMON /CHEM2/     NPHOT,       NPRODLO,   NPRODHI,   MSTEP,     
     &                   MAXORD,      MBETWEEN,  IC3H8,     IC2H6

      ! /CHEM2A/ has to be held THREADPRIVATE (bmy, 7/16/03)
      INTEGER ::         ISCHAN,      NFDH3,     NFDL2,     NFDH2
      INTEGER ::         NFDL1,       NFDH1,     NFDREP,    NFDREP1   
      INTEGER ::         NFDL0,       NALLR
      COMMON /CHEM2A/    ISCHAN,      NFDH3,     NFDL2,     NFDH2,       
     &                   NFDL1,       NFDH1,     NFDREP,    NFDREP1,     
     &                   NFDL0,       NALLR

      ! Split off from /CHEM2A/ (bmy, 7/28/03)
      INTEGER ::         NOCC
      COMMON /CHEM2B/    NOCC
         
      INTEGER ::         NGAS,        NMREAC
      COMMON /CHEM3/     NGAS,        NMREAC

      ! Added NNADDG to /CHEM4/ for DMS+OH+O2 rxn (bdf, bmy, 4/18/03)
      ! Add NKN2O5 to /CHEM4/ to flag N2O5 hydrolysis rxn (mje, bmy, 8/7/03)
      ! Add NKHO2 to /CHEM4/ to flag HO2 aerosol uptake (jaegle 02/26/09)
      ! Add NNADDF, NNADDH, NKHOROI and NKHOROJ to /CHEM4/ for HOC2H4O rxns
      ! (tmf, 3/6/09)
      !Added NKGLYC and NKHAC NKMCO3
      ! (hotp 7/31/09)
      INTEGER ::         NNADD1,      NNADDA,      NNADDB
      INTEGER ::         NNADDC,      NNADDD,      NNADDK
      INTEGER ::         NNADDV,      NNADDZ,      NKO3PHOT
      INTEGER ::         NNADDG,      NEMIS,       NDRYDEP
      INTEGER ::         NKHNO4,      NKN2O5,      NKHO2
      INTEGER ::         NNADDF,      NNADDH
      INTEGER ::         NKHOROI,     NKHOROJ
      INTEGER ::         NKGLYC,      NKHAC
      INTEGER ::         NKMCO3
      COMMON /CHEM4/     NNADD1,      NNADDA(ICS), NNADDB(  ICS), 
     &                   NNADDC(ICS), NNADDD(ICS), NNADDK(  ICS), 
     &                   NNADDV(ICS), NNADDZ,      NKO3PHOT(ICS),
     &                   NNADDG(ICS), NEMIS( ICS), NDRYDEP( ICS),
     &                   NNADDF(ICS), NNADDH(ICS), 
     &                   NKHOROI(ICS),NKHOROJ(ICS),
     &                   NKHNO4(ICS), NKN2O5,      NKHO2,
     &                   NKGLYC(ICS,2), NKHAC(ICS,2), NKMCO3(ICS,3)

      INTEGER ::         IH2O,        IOXYGEN,   MB1,      MB2
      COMMON /SPECIES/   IH2O,        IOXYGEN,   MB1,      MB2

      ! Added for interannually-varying Methane (bnd, bmy, 7/1/03)
      INTEGER ::         ICH4
      COMMON /SPECIE2/   ICH4

      ! Added for interannually-varying Methane (bnd, bmy, 7/1/03)
      REAL*8  ::         C3090S,      C0030S,    C0030N,   C3090N
      COMMON /SPECIE3/   C3090S,      C0030S,    C0030N,   C3090N

      ! Added for tracking oxidation of ISOP by OH (dkh, bmy, 6/1/06)
      ! SOAupdate: Added for tracking oxidation of ISOP by NO3 (hotp, mpayer,
      !  6/28/11)
      INTEGER ::         ILISOPOH, ILISOPNO3
      COMMON /SPECIE4/   ILISOPOH, ILISOPNO3

      ! Added for tracking oxidation of aromatic RO2s (dkh, 10/06/06)  
      ! bug fix: change to INTEGER (dbm, dkh, 06/26/07)
      ! SOAupdate: Added ILNRO2H and ILNRO2N (hotp, mpayer, 6/28/11) 
      INTEGER ::         ILBRO2H,    ILBRO2N,   ILTRO2H,  ILTRO2N
      INTEGER ::         ILXRO2H,    ILXRO2N,   ILNRO2H,  ILNRO2N
      COMMON /SPECIE5/   ILBRO2H,    ILBRO2N,   ILTRO2H,  ILTRO2N,
     &                   ILXRO2H,    ILXRO2N,   ILNRO2H,  ILNRO2N

      INTEGER ::         IOUT,        KGLC,      KCPD,     IO93
      COMMON /FILES/     IOUT,        KGLC,      KCPD,     IO93
C
C *********************************************************************
C *               SET REAL AND INTEGER ARRAY VARIABLES                *
C *********************************************************************
C
      INTEGER ::        JLOWVAR,           KTLPVAR
      INTEGER ::        JLOFIXED,          JHIFIXED
      COMMON /IMXBLOCK/ JLOWVAR( MXBLOCK), KTLPVAR( MXBLOCK), 
     &                  JLOFIXED(MXBLOCK), JHIFIXED(MXBLOCK)

      INTEGER ::        JREORDER,         LREORDER 
      INTEGER ::        ITWO,             NCSLOOP
      COMMON /IITLOOP/  JREORDER(ITLOOP), LREORDER(ITLOOP), 
     &                  ITWO(    ITLOOP), NCSLOOP( ITLOOP,ICS)

      ! Add NKSPECG for DMS+OH+O2 rxn (bdf, bmy, 4/18/03)
      ! Added NKSPECF, NKSPECH to /IDICS/ for C2H4 chemistry (tmf, 3/6/09) 
      INTEGER NMOTH,NTSPEC,JPHOTRAT,ISGAINR,ISPORL,NOGAINE,NOUSE
      INTEGER NSPEC,NTRATES,ISGAINE,NTLOOPNCS,NSPCSOLV,ISCHANG,NRATES
      INTEGER NM3BOD,ITWOR,ITHRR,INOREP,NRATCUR,NSURFACE,NPRESM,NMAIR
      INTEGER NMO2,NMN2,NNEQ,NARR,NABR,NACR,NABC,NKSPECW,NKSPECX
      INTEGER NKSPECY,NKSPECZ,NKSPECV,ISLOSSR,NKSPECA,NKSPECB,NKSPECC
      INTEGER NKSPECD,NKSPECK,NKSPECG
      INTEGER NKSPECF, NKSPECH

      COMMON /IDICS/
     1  NMOTH(   ICS), NTSPEC( ICS), JPHOTRAT(ICS),
     3  ISGAINR( ICS), ISPORL( ICS), NOGAINE( ICS),   NOUSE(    ICS),
     4  NSPEC(   ICS), NTRATES(ICS), ISGAINE( ICS),   NTLOOPNCS(ICS),
     5  NSPCSOLV(ICS), ISCHANG(ICS), NRATES(  ICS),   NM3BOD(   ICS),
     7  ITWOR(   ICS), ITHRR(  ICS), INOREP(  ICS),   NRATCUR(  ICS),
     8  NSURFACE(ICS), NPRESM( ICS), NMAIR(   ICS),   NMO2(     ICS),
     9  NMN2(    ICS), NNEQ(   ICS), NARR(    ICS),   NABR(     ICS),
     1  NACR(    ICS), NABC(   ICS), NKSPECW( ICS),   NKSPECX(  ICS),
     2  NKSPECY( ICS), NKSPECZ(ICS), NKSPECV(MAXGL2,ICS),ISLOSSR(ICS),
     3  NKSPECA(  MAXGL3,ICS), NKSPECB(  MAXGL3,ICS),
     4  NKSPECC(MAXGL3,ICS),NKSPECD(MAXGL3,ICS),NKSPECK(MAXGL3,ICS),
     5  NKSPECG(MAXGL2,ICS),
     6  NKSPECF(MAXGL3,ICS), NKSPECH(MAXGL3,ICS)

      ! re-define some nkspec* arrays for harvard chem mechanism (bdf)
      INTEGER ::        NOLOSP,         NGNFRAC,         NOLOSRAT
      INTEGER ::        IARRAY,         NALLRAT,         KZTLO
      INTEGER ::        KZTHI,          IONER,           NPLLO
      INTEGER ::        NPLHI,          NFRLO,           NFRHI
      INTEGER ::        NPDLO,          NPDHI,           IZLO
      INTEGER ::        JZLO,           JLLO,            JGLO
      INTEGER ::        IRMCT
      COMMON /IICP/     NOLOSP(ICP),    NGNFRAC(ICP),    NOLOSRAT(ICP), 
     &                  IARRAY(ICP),    NALLRAT(ICP),    KZTLO(   ICP), 
     &                  KZTHI( ICP),    IONER(  ICP),    NPLLO(   ICP),  
     &                  NPLHI( ICP),    NFRLO(  ICP),    NFRHI(   ICP),
     &                  NPDLO( ICP),    NPDHI(  ICP),    IZLO (   ICP), 
     &                  JZLO ( ICP),    JLLO(   ICP),    JGLO(    ICP), 
     &                  IRMCT( ICP)

      REAL*8  ::        ABTOL,          ABST2
      REAL*8  ::        ERRMAX,         HMAXUSE,         TIMEINTV
      COMMON /DICS/     ABTOL(6,ICS),   ABST2(ICS),  
     &                  ERRMAX(ICS),    HMAXUSE(ICP),    TIMEINTV(ICS)    

      REAL*8  ::        WTGAS,          GQSUMINI
      REAL*8  ::        BSUMCHEM,       GQSUM,           QBKGAS
      COMMON /DIGAS/    WTGAS(IGAS),    GQSUMINI(IGAS), 
     &                  BSUMCHEM(IGAS), GQSUM(   IGAS),  QBKGAS(IGAS)       

      REAL*8  ::        CPREV,          CMODEL,          APORL
      COMMON /DMXGAER1/ CPREV(MXGSAER), CMODEL(MXGSAER), APORL(MXGSAER)

      INTEGER ::        IFPRGAS,        LGNUM,           NGMIX
      COMMON /IIGAS/    IFPRGAS(IGAS),  LGNUM(IGAS),     NGMIX(IGAS)

      REAL*8  ::        DEFPRAT
      COMMON /DIPHOT/   DEFPRAT(MXRATE,ICS)

      REAL*8  ::        ARRT,           BRRT
      REAL*8  ::        FCVT,           FCT1T,           FCT2T
      COMMON /DMXCOF/   ARRT(MXCOF),    BRRT( MXCOF),   
     &                  FCVT(MXCOF),    FCT1T(MXCOF),    FCT2T(MXCOF)

      INTEGER ::        KCRRT   
      COMMON /IMXCOF/   KCRRT(MXCOF)

      INTEGER NKARR,NKABR,NKACR,NKABC,IRORD
      COMMON /INMRAT2/
     1  NKARR(NMTRATE,ICS),  NKABR(NMTRATE,ICS), NKACR(NMTRATE,ICS), 
     2  NKABC(NMTRATE,ICS),  IRORD(NMTRATE,ICS)

      REAL*8 ARR,BRR,FCV,FCTEMP1,FCTEMP2
      COMMON /DNMTRATE/
     1  ARR(    NMTRATE, ICS),  BRR(    NMTRATE, ICS), 
     2  FCV(    NMTRATE, ICS),  FCTEMP1(NMTRATE, ICS),
     3  FCTEMP2(NMTRATE, ICS)  

      INTEGER IAPROD,NOLOSRN,NRUSE,NRREP,NPRODUC,IALLOSN,NCEQUAT
      INTEGER NEWFOLD,NKONER,NKTWOR,NKTHRR,IRMA,IRMB,KCRR,LSKIP,IRMC
      INTEGER JPHOTNK,IUSED,NOLDFNEW
      COMMON /INMTRATE/
     2  IAPROD( NMTRATE,  ICS),  NOLOSRN( NMTRATE, ICS),
     3  NRUSE(  NMTRATE,  ICS),  NRREP(   NMTRATE, ICS),
     4  NPRODUC(NMTRATE,  ICS),  IALLOSN( MXRATE,  ICS),
     5  NCEQUAT(NMTRATE,  ICS),  NOLDFNEW(NMTRATE, ICS),
     6  NEWFOLD(NMTRATE*2,ICS),  NKONER(  NMTRATE, ICS),
     7  NKTWOR( NMTRATE,  ICS),  NKTHRR(  NMTRATE, ICS),
     9  KCRR(   NMTRATE,  ICS),  LSKIP(   MXRATE,  ICS),
     1  IRMC(   NMTRATE      ),  JPHOTNK( NMTRATE, ICS),
     2  IUSED(  MXRATE,   ICS)

      ! For COMPAQ, put IRMA, IRMB in /INMTRATE2/ (Q. Liang, bmy, 10/17/05)
      COMMON /INMTRATE2/
     & IRMA(   NMTRATE      ),  IRMB(    NMTRATE     )

      INTEGER ::       NEWNK
      COMMON /IMAXGL3/ NEWNK(MAXGL)

      REAL*8  ::       FRACP 
      COMMON /DMAXGL2/ FRACP(MAXGL, ICS)

      INTEGER NREACOTH,LGASBINO,NKNLOSP,LOSINACP,IGNFRAC,NKGNFRAC
      INTEGER NREACAIR,NREAC3B,NREACEQ,NREQOTH,NREACN2,NREACO2,NREACPM
      INTEGER LGAS3BOD,NKSURF,NCOATG
      COMMON /IMAXGL2/
     1  NREACOTH(MAXGL2,ICS), LGASBINO(MAXGL2,ICS),
     2  NKNLOSP( MAXGL3,ICS), LOSINACP(MAXGL3,ICS),
     3  IGNFRAC( MAXGL, ICS), NKGNFRAC(MAXGL, ICS),
     4  NREACAIR(MAXGL3,ICS), NREAC3B( MAXGL3,ICS),
     5  NREACEQ( MAXGL3,ICS), NREQOTH( MAXGL3,ICS),
     6  NREACN2( MAXGL3,ICS), NREACO2( MAXGL3,ICS),
     7  NREACPM( MAXGL3,ICS), LGAS3BOD(MAXGL3,ICS),
     8  NKSURF(  MAXGL4    ), NCOATG(  MAXGL4    ) 

      INTEGER ::        MBCOMP,            MBTRACE
      COMMON /IIMASBAL/ MBCOMP(IMASBAL,2), MBTRACE(IMASBAL)

      ! /DKBLOOP2/ needs to be declared THREADPRIVATE
      REAL*8 CNEW,CEST,GLOSS,CHOLD,VDIAG,CBLK,DTLOS,EXPLIC,CONC
      REAL*8 RRATE,URATE,TRATE,CORIG
      !***************KPP_INTERFACE****************
      REAL*8 RRATE_FOR_KPP
      !********************************************
      COMMON /DKBLOOP2/
     2  CNEW(   KBLOOP,  MXGSAER),  
     3  CEST(   KBLOOP,  MXGSAER),  
     4  GLOSS(  KBLOOP,  MXGSAER),  
     5  CHOLD(  KBLOOP,  MXGSAER), 
     6  VDIAG(  KBLOOP,  MXGSAER),  CBLK(  KBLOOP,MXGSAER),  
     7  DTLOS(  KBLOOP,  MXGSAER),  EXPLIC(KBLOOP,MXGSAER),
     1  CONC(   KBLOOP,MXGSAER*7),  
     2  RRATE(  KBLOOP,  NMTRATE), 
      !***************KPP_INTERFACE****************
     2  RRATE_FOR_KPP(  KBLOOP,  NMTRATE), 
      !********************************************  
     3  URATE(  KBLOOP,NMTRATE,3), 
     4  TRATE(  KBLOOP,NMTRATE*2), 
     7  CORIG(  KBLOOP,  MXGSAER)

      ! /DKBLOOP0/ needs to be held THREADPRIVATE
      REAL*8 ::         CC2
      COMMON /DKBLOOP0/ CC2(KBLOOP,0:MXARRAY)  

      INTEGER ::      MLOP,             JLOP_SMV  
      COMMON /IILAT2/ MLOP(ILAT,ILONG), JLOP_SMV(ILAT,ILONG,ILAYER)

      INTEGER NKPHOTRAT,NPPHOTRAT,NKNPHOTRT
      COMMON /DIPHOT2/ 
     1  NKPHOTRAT(IPHOT,ICS),      NPPHOTRAT(IPHOT,ICS),  
     2  NKNPHOTRT(IPHOT,ICS) 

      REAL*8 ::        FRACGAIN,              QBKCHEM
      COMMON /DIMXGS2/ FRACGAIN(MXGSAER,ICS), QBKCHEM( MXGSAER,ICS) 

      INTEGER NUMLOST,NUMGFRT,NUMLOSS,JPORL,NUMGAINT,NGAINE,NUMGAIN
      INTEGER IGAINR,IPORL,IGAINE,ISOLVSPC,INEWOLD,MAPPL,ISAPORL,NUMPORL
      INTEGER ISPARDER,JLOSST
      COMMON /IIMXGS2/
     1  NUMLOST( MXGSAER,    ICS), NUMGFRT( MXGSAER,        ICS),  
     2  NUMLOSS( MXGSAER,    ICP), JPORL(   MXGSAER,MAXGL,  ICS),
     3  NUMGAINT(MXGSAER,    ICS), NGAINE(  MXGSAER,        ICS),  
     4  NUMGAIN( MXGSAER,    ICP), IGAINR(  MXGSAER,        ICS), 
     9  IPORL(   MXGSAER,    ICS), IGAINE(  MXGSAER,        ICS),  
     2  ISOLVSPC(MXGSAER,    ICS), INEWOLD( MXGSAER,        ICS),
     3  MAPPL(   MXGSAER,    ICS), ISAPORL( MXGSAER            ),
     7  NUMPORL( MXGSAER,    ICP), ISPARDER(MXGSAER,MXGSAER    ),
     8  JLOSST(  MXGSAER,MAXGL,ICS)

      INTEGER JZILCH,KZILCH,MZILCH
      COMMON /IGMXGLS/
     &     JZILCH(MXGSAER),  KZILCH(MXGSAER),  MZILCH(MXGSAER) 

      INTEGER LZERO,JARRAYPT,IZILCH,JARRDIAG,JLOZ1,JHIZ1,IJTLO
      INTEGER IJTHI,IMZTOT,IFREPRO,IZLO1,IZLO2,IZHI0,IZHI1
      COMMON /IGMXGS2/
     1  LZERO(   MXGSAER,MXGSAER), JARRAYPT(MXGSAER,MXGSAER),
     2  IZILCH(  MXGSAER,MXGSAER), JARRDIAG(MXGSAER,   ICP),
     3  JLOZ1(   MXGSAER,    ICP), JHIZ1(   MXGSAER,   ICP), 
     4  IJTLO(   MXGSAER,    ICP), IJTHI(   MXGSAER,   ICP),
     5  IMZTOT(  MXGSAER,    ICP), IFREPRO( MXGSAER,MXRATE, ICS),
     6  IZLO1( MXGSAER,ICP),
     7  IZLO2( MXGSAER,ICP), IZHI0( MXGSAER,ICP), IZHI1( MXGSAER,ICP)

      REAL*8 ::        FRACNFR,           FRACPL
      COMMON /DMXCOUN/ FRACNFR(MXCOUNT4), FRACPL(MXCOUNT2) 

      INTEGER JZERO,KZERO,MZERO,IZEROK,JZEROA,IKDECA,KJDECA,LOSSRA
      INTEGER IKDECB,KJDECB,LOSSRB,IKDECC,KJDECC,LOSSRC,IKDECD,KJDECD
      INTEGER LOSSRD,IKDECE,KJDECE,LOSSRE,KZEROA,MZEROA,KZEROB,MZEROB
      INTEGER KZEROC,MZEROC,KZEROD,MZEROD,KZEROE,MZEROE,IPOSPD,IIALPD
      INTEGER NKPDTERM,IJVAL,IKZTOT,JSPNPL,NKNFR,JSPCNFR
      COMMON /IMXCOUN/
     1  JZERO( MXCOUNT3),  KZERO(  MXCOUNT3),  MZERO(   MXCOUNT3),
     2  IZEROK(MXCOUNT2),  JZEROA( MXCOUNT3),
     3  IKDECA(MXCOUNT3),  KJDECA( MXCOUNT3),  LOSSRA(  MXCOUNT4),  
     4  IKDECB(MXCOUNT3),  KJDECB( MXCOUNT3),  LOSSRB(  MXCOUNT4), 
     5  IKDECC(MXCOUNT3),  KJDECC( MXCOUNT3),  LOSSRC(  MXCOUNT4), 
     6  IKDECD(MXCOUNT3),  KJDECD( MXCOUNT3),  LOSSRD(  MXCOUNT4), 
     7  IKDECE(MXCOUNT3),  KJDECE( MXCOUNT3),  LOSSRE(  MXCOUNT4),
     8  KZEROA(MXCOUNT4),  MZEROA( MXCOUNT4), 
     9  KZEROB(MXCOUNT4),  MZEROB( MXCOUNT4), 
     1  KZEROC(MXCOUNT4),  MZEROC( MXCOUNT4), 
     2  KZEROD(MXCOUNT4),  MZEROD( MXCOUNT4), 
     3  KZEROE(MXCOUNT4),  MZEROE( MXCOUNT4), 
     4  IPOSPD(MXCOUNT2),  IIALPD( MXCOUNT2),  NKPDTERM(MXCOUNT2),  
     5  IJVAL( MXCOUNT3),  IKZTOT( MXCOUNT4),  JSPNPL(  MXCOUNT4),    
     7  NKNFR( MXCOUNT4),  JSPCNFR(MXCOUNT4) 

      INTEGER IDH5,IDH4,IDH3,IDH2,IDH1,IDL5,IDL4,IDL3,IDL2,IDL1,KBH5
      INTEGER KBH4,KBH3,KBH2,KBH1,KBL5,KBL4,KBL3,KBL2,KBL1,MBH5,MBH4
      INTEGER MBH3,MBH2,MBH1,MBL5,MBL4,MBL3,MBL2,MBL1,NPH5,NPH4,NPH3
      INTEGER NPH2,NPH1,NPL5,NPL4,NPL3,NPL2,NPL1
      COMMON /IMXCOU2/
     1  IDH5(  MXCOUNT3),  IDH4(  MXCOUNT3),  IDH3(  MXCOUNT3),
     2  IDH2(  MXCOUNT3),  IDH1(  MXCOUNT3),  IDL5(  MXCOUNT3),
     3  IDL4(  MXCOUNT3),  IDL3(  MXCOUNT3),  IDL2(  MXCOUNT3),
     4  IDL1(  MXCOUNT3),  
     5  KBH5(  MXCOUNT4),  KBH4(  MXCOUNT4),  KBH3(  MXCOUNT4),
     6  KBH2(  MXCOUNT4),  KBH1(  MXCOUNT4),  KBL5(  MXCOUNT4),
     7  KBL4(  MXCOUNT4),  KBL3(  MXCOUNT4),  KBL2(  MXCOUNT4),
     8  KBL1(  MXCOUNT4),  
     9  MBH5(  MXCOUNT4),  MBH4(  MXCOUNT4),  MBH3(  MXCOUNT4),
     1  MBH2(  MXCOUNT4),  MBH1(  MXCOUNT4),  MBL5(  MXCOUNT4),
     2  MBL4(  MXCOUNT4),  MBL3(  MXCOUNT4),  MBL2(  MXCOUNT4),
     3  MBL1(  MXCOUNT4),  
     4  NPH5(  MXCOUNT4),  NPH4(  MXCOUNT4),  NPH3(  MXCOUNT4),
     5  NPH2(  MXCOUNT4),  NPH1(  MXCOUNT4),  NPL5(  MXCOUNT4),
     6  NPL4(  MXCOUNT4),  NPL3(  MXCOUNT4),  NPL2(  MXCOUNT4),
     7  NPL1(  MXCOUNT4)  

      REAL*8 ::          GEARCONC 
      COMMON /DIMXG2/    GEARCONC(MXGSAER,0:MXHOLD,ICS)

      REAL*8 ::          WTMB
      COMMON /DIMASBAL2/ WTMB(IMASBAL,MXGSAER,2)

      INTEGER ::         JMBCOMP
      COMMON /IIMASBAL2/ JMBCOMP(IMASBAL,MXGSAER,2) 

      REAL*8 ::          FKOEF
      REAL*8 ::          FK2 
      COMMON /DNMRPROD2/ FKOEF(NMRPROD,NMTRATE,ICS), 
     &                   FK2(  NMRPROD,NMTRATE,ICS)  

      INTEGER ::         IRM
      INTEGER ::         IRM2
      COMMON /INMRPROD2/ IRM( NMRPROD,NMTRATE,ICS), 
     &                   IRM2(NMRPROD,NMTRATE,ICS) 

      REAL*8 ::          ASET,       PINP,     CVAR,     O3DOBS
!      COMMON /DMISC/     ASET(10,8), PINP(20), CVAR(15), O3DOBS(12)  
      COMMON /DMISC/     ASET(10,8), PINP(IMISC), CVAR(15), O3DOBS(12)  

      REAL*8 ::          ENQQ2,         ENQQ3,          CONPST
      REAL*8 ::          ENQQ1,         CONP15
      COMMON /IORDR/     ENQQ2(MORDER), ENQQ3( MORDER), CONPST(MORDER),
     &                   ENQQ1(MORDER), CONP15(MORDER)

      REAL*8 ::          PERTS2,           PERTST
      COMMON /DMORD/     PERTS2(MORDER,3), PERTST(MORDER,3)

      INTEGER ::         JLLOW,          KLHI
      COMMON /IMXBLOCK/  JLLOW(MXBLOCK), KLHI(MXBLOCK)

      REAL*8 ::          XGDFCF, ASTKCF, RUARSL, RH100
      COMMON /XARSOL/    XGDFCF, ASTKCF, RUARSL, RH100

      INTEGER ::         IARSFA, MWARSL, MNTHARS
      COMMON /IARSOL/    IARSFA, MWARSL, MNTHARS

      INTEGER ::         NKEMIS,              NTEMIS
      INTEGER ::         NKDRY,               NTDEP
      COMMON /IMAXGL4/   NKEMIS(MAXGL3,ICS),  NTEMIS(MAXGL3,ICS),
     &                   NKDRY (MAXGL3,ICS),  NTDEP( MAXGL3)

C BDF COULD BE SOME MISSING STUFF FROM HERE.
C      INTEGER KZERO,IZEROA,IZEROB,IZEROC,IZEROD,IZER2A,IZER2B,
C     1        IZER2C,IZER2D,JZEROA,JZEROB,IRMSEC,IRMTHD,NKARRY,
C     2        LOSSRA,JPRODA,LOSSRB,LOSSRC,LOSSRD,JPRODB,JPRODC,
C     3        JPRODD
C      COMMON /IMXCOUN/
C     1  KZERO( MXARRAY,ICP), IZEROA( MXARRAY   ), 
C     7  IZEROB( MXARRAY),    IZEROC( MXARRAY   ),  IZEROD( MXARRAY),
C     8  IZER2A( MXARRAY),    IZER2B( MXARRAY   ),
C     9  IZER2C( MXARRAY),    IZER2D( MXARRAY   ),
C     1  JZEROA( MXARRAY),    JZEROB( MXARRAY   ),
C     3  IRMSEC(MXCOUNT3),    IRMTHD(MXCOUNT3   ),  NKARRY(MXCOUNT3),
C     2  LOSSRA(MXCOUNT4),    JPRODA(MXCOUNT4   ),
C     4  LOSSRB(MXCOUNT4),    LOSSRC(MXCOUNT4   ),  LOSSRD(MXCOUNT4),
C     5  JPRODB(MXCOUNT4),    JPRODC(MXCOUNT4   ),
C     6  JPRODD(MXCOUNT4),    JPRODT(  MXGSAER,MAXGL,ICS),
C     5  NKSINGL( MXGSAER,ICS,    2), NKNUMSL( MXGSAER,MAXGL,ICS),
C     6  NKDOUBL( MXGSAER,ICS      ), NKNUMDL( MXGSAER,MAXGL,ICS),
C     7  NKTRIPL( MXGSAER,ICS      ), NKNUMTL( MXGSAER,MAXGL,ICS),
C     8  LOSSLEFT(MXGSAER,ICS      ), LOSSREM( MXGSAER,MAXGL,ICS),
C     9  ILOSSR(  MXGSAER,ICS      ), NLOSSR(  MXGSAER,      ICS),
C     1  NGAINR(  MXGSAER,ICS      ), ICHANSPC(MXGSAER,      ICS),
C     4  NUML1(   MXGSAER,ICP      ), NUML2(   MXGSAER,      ICP),
C     5  NUMP1(   MXGSAER,ICP      ), NUMP2(   MXGSAER,      ICP),
C     6  JHIZ2(   MXGSAER,      ICP)

      !=================================================================
      ! Common blocks for ND65 diagnostic (ljm, bmy, 7/8/03)
      !=================================================================
      INTEGER     ::     IFAM,         NFAMILIES
      COMMON /IPL/       IFAM(MAXFAM), NFAMILIES

      CHARACTER*4 ::     PORL
      COMMON /CPL/       PORL(MAXFAM)

      LOGICAL     ::     LFAMILY, ITS_NOT_A_ND65_FAMILY
      COMMON /LPL/       LFAMILY, ITS_NOT_A_ND65_FAMILY(IGAS)

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

#if   defined( COMPAQ )
! For COMPAQ, declare /INMTRATE2/ threadprivate (Q. Liang, bmy, 10/17/05)
!$OMP THREADPRIVATE( /INMTRATE2/ )
#endif

C
C *********************************************************************
C ****************** END OF COMMON BLOCK COMODE.H *********************
C *********************************************************************
C
