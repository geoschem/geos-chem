! $Id: comode.h,v 1.4 2003/07/21 15:09:25 bmy Exp $
!
!******************************************************************************
!  Header file COMODE contains common blocks and variables for SMVGEAR II.
!  (M. Jacobson 1997; bdf, bmy, 4/23/03, 7/16/03)
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
      !PARAMETER (KBLOOP  = 64                      )
      PARAMETER (KBLOOP  = 24                      )
      PARAMETER (IMLOOP  = ILAT * ILONG            )
      PARAMETER (ILAYER  = IVERT + 1               )
      PARAMETER (ILTLOOP = IMLOOP * ILAYER         )
      PARAMETER (MAXDAYS = 1000                    )
C Debug
C     PARAMETER (KBLOOP  = 1                       )
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

#if   defined ( FULLCHEM )
      ! Updated for SMVGEAR II (bdf, bmy, 4/1/03)
      PARAMETER ( IGAS    = 125,               IAERTY  = 1           )
      PARAMETER ( NMRATE  = 360,               IPHOT   = 60          )
      PARAMETER ( NMTRATE = NMRATE + IPHOT,    NMQRATE = 1           ) 
      PARAMETER ( NMRPROD = 25,                NMDEAD  = 100         )
      PARAMETER ( MAXGL   = 430,               MAXGL2  = 50          )
      PARAMETER ( MAXGL3  = 32,                MAXGL4  = 10          )
      PARAMETER ( ICS     = 3,                 ICP     = 2*ICS       ) 
      PARAMETER ( MORDER  = 7                                        ) 
      PARAMETER ( IPHOT8  = IPHOT + 8,         IMISC   = 100         )
      PARAMETER ( IMASBAL = 9,                 IALTS   = 22          )
      PARAMETER ( MXCOF   = 5                                        )

#elif defined ( SMALLCHEM ) || defined ( LGEOSCO )
      ! New settings for small chemistry to save array space (bmy, 1/5/98)
      ! Need these to also be defined for LGEOSCO run (bmy, 10/3/00)
      PARAMETER ( IGAS    = 35,                IAERTY  = 1           )
      PARAMETER ( NMRATE  = 70,                IPHOT   = 15          )
      PARAMETER ( NMTRATE = NMRATE + IPHOT,    NMQRATE = 1           )
      PARAMETER ( NMRPROD = 25,                NMDEAD  = 15          )
      PARAMETER ( MAXGL   = 130,               MAXGL2  = 50          )
      PARAMETER ( MAXGL3  = 25                                       )
      PARAMETER ( ICS     = 1,                 ICP     = 2*ICS       )
      PARAMETER ( MORDER  = 7                                        )
      PARAMETER ( IPHOT8  = IPHOT + 8,         IMISC   = 100         )
      PARAMETER ( IMASBAL = 9,                 IALTS   = 22          )
      PARAMETER ( MXCOF   = 5                                        )

#endif
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
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Prior to 7/16/03:
! Remove HALFDAY, GRAVC, FOURPI, TWOPI, REARTH, RPRIMB, AVOG1, HALF, THIRD,
! THRPI2, PID180, PID2, SCTWOPI, AMRGAS, TWPISC -- they are obsolete.
! (bmy, 7/16/03)
!      REAL*8 AVG,BOLTG,HALFDAY,GRAVC
!      REAL*8 FOURPI,    TWOPI,     RGAS,      REARTH,    RPRIMB
!      REAL*8 SCDAY,     AVOG1,     HALF,      THIRD,     BK
!      REAL*8 THRPI2,    PID180,    PID2,      EIGHTDPI,  RSTARG
!      REAL*8 SCTWOPI,   WTAIR,     ONEPI,     AMRGAS,    TWPISC
!      COMMON /BSCPARM /
!     1  AVG,       BOLTG,     HALFDAY,   GRAVC,
!     2  FOURPI,    TWOPI,     RGAS,      REARTH,    RPRIMB,
!     3  SCDAY,     AVOG1,     HALF,      THIRD,     BK,
!     4  THRPI2,    PID180,    PID2,      EIGHTDPI,  RSTARG,
!     5  SCTWOPI,   WTAIR,     ONEPI,     AMRGAS,    TWPISC
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      REAL*8  ::         AVG,       BOLTG,     RGAS
      REAL*8  ::         SCDAY,     BK,        EIGHTDPI  
      REAL*8  ::         RSTARG,    WTAIR,     ONEPI
      COMMON /BSCPARM/   AVG,       BOLTG,     RGAS, 
     &                   SCDAY,     BK,        EIGHTDPI,  
     &                   RSTARG,    WTAIR,     ONEPI
                         
      INTEGER ::         NLAT,      NLONG,     NLAYER,    NVERT
      INTEGER ::         NLOOP,     NTLOOP,    KULOOP,    NTSPECGAS
      INTEGER ::         NMASBAL,   KSLOOP,    NTLOOPUSE, NPVERT
      INTEGER ::         NTTLOOP,   NIJLOOP
      COMMON /CTLLOOP/   NLAT,      NLONG,     NLAYER,    NVERT,     
     &                   NLOOP,     NTLOOP,    KULOOP,    NTSPECGAS, 
     &                   NMASBAL,   KSLOOP,    NTLOOPUSE, NPVERT,    
     &                   NTTLOOP,   NIJLOOP

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

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Prior to 7/16/03:
! Remove DX0, DY0, XU0, DTOUT, CONPSUR, DXLONG, DYLAT, SWLONDC, CONSTIM, 
! SWLATDC, UTSECY, TOTSEC, FINHOUR, FINMIN, FINSEC, TFROMID, ZENFIXED,
! ZENITH, DENCONS from common block -- they're obsolete. (bmy, 7/16/03)
!      REAL*8 DX0,       DY0,       TINTERVAL,
!     2     XU0,       DTOUT,     CONPSUR,   CHEMINTV,
!     3     DXLONG,    DYLAT,     SWLONDC,   CONSTIM,   SWLATDC,
!     4     FRACDEC,   UTSECY,    TIME,      CONSVAP,
!     5     TOTSEC,    FINHOUR,   FINMIN,    FINSEC,
!     6     TFROMID,   ZENFIXED,  ZENITH,
!     7     OXYCONS,   DENCONS,   HMAXNIT,
!     8     SMAL1,     SMAL2,     SMAL3
!      COMMON /XYGRID /
!     1  DX0,       DY0,       TINTERVAL,
!     2  XU0,       DTOUT,     CONPSUR,   CHEMINTV,
!     3  DXLONG,    DYLAT,     SWLONDC,   CONSTIM,   SWLATDC,
!     4  FRACDEC,   UTSECY,    TIME,      CONSVAP,
!     5  TOTSEC,    FINHOUR,   FINMIN,    FINSEC,
!     6  TFROMID,   ZENFIXED,  ZENITH,
!     7  OXYCONS,   DENCONS,   HMAXNIT,
!     8  SMAL1,     SMAL2,     SMAL3
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      REAL*8  ::        TINTERVAL,  CHEMINTV,  TIME,      CONSVAP
      REAL*8  ::        OXYCONS,    HMAXNIT,   SMAL1,     SMAL2
      REAL*8  ::        SMAL3,      FRACDEC
      COMMON /XYGRID/   TINTERVAL,  CHEMINTV,  TIME,      CONSVAP, 
     &                  OXYCONS,    HMAXNIT,   SMAL1,     SMAL2,  
     &                  SMAL3,      FRACDEC
C
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Prior to 7/16/03:
! Remove IMIN, ISEC, IDAY_SMV, IMONTH, IYEAR_SMV, IRUN, NRUN, NDAYS, LEAP,
! NSOUT, IDAYR, LASTCHEM, LOTEMP1 from the common block...these are not used.
! (bmy, 7/16/03)
!      INTEGER IHOUR,     IMIN,      ISEC,      IDAY_SMV,  IMONTH,
!     2        IYEAR_SMV, IRUN,      NRUN,      NDAYS,     LEAP,
!     3        NSOUT,     IDAYR,     NCS,       NCSP,      LASTCHEM,
!     4        NSTEPS,    KBLK,      NBLOCKS,   IRCHEM,    NCSGAS,
!     5        NCSURBAN,  NCSTROP,   NCSSTRAT,  NPHOTALL,  IFDID,
!     6        IFNEVER,   IFNONE,    LOTEMP1,   NCSALL,    NCSTRST
!      COMMON /IXYGD /
!     1  IHOUR,     IMIN,      ISEC,      IDAY_SMV,  IMONTH,
!     2  IYEAR_SMV, IRUN,      NRUN,      NDAYS,     LEAP,
!     3  NSOUT,     IDAYR,     NCS,       LASTCHEM,
!     4  NBLOCKS,   IRCHEM,    NCSGAS,
!     5  NCSURBAN,  NCSTROP,   NCSSTRAT,  NPHOTALL,  IFDID,
!     6  IFNEVER,   IFNONE,    LOTEMP1,   NCSALL,    NCSTRST
!C
! Prior to 7/16/03:
! Remove NSTEPS from /IXYGD2/ since that can be declared as a local variable
! w/in "smvgear.f".  Also tidy up the declarations below. (bmy, 7/16/03)
!      COMMON /IXYGD2 /
!     1  NCSP,  NSTEPS,  KBLK
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      INTEGER ::       IHOUR,       NCS,       NBLOCKS,   IRCHEM
      INTEGER ::       NCSGAS,      NCSURBAN,  NCSTROP,   NCSSTRAT
      INTEGER ::       NPHOTALL,    IFDID,     IFNEVER,   IFNONE   
      INTEGER ::       NCSALL,      NCSTRST
      COMMON /IXYGD/   IHOUR,       NCS,       NBLOCKS,   IRCHEM,   
     &                 NCSGAS,      NCSURBAN,  NCSTROP,   NCSSTRAT, 
     &                 NPHOTALL,    IFDID,     IFNEVER,   IFNONE,    
     &                 NCSALL,      NCSTRST

      ! /IXYGD2/ needs to be held THREADPRIVATE.  Also remove NSTEPS
      ! since this can be declared local w/in "smvgear.f" (bmy, 7/16/03)
      INTEGER ::       NCSP,        KBLK
      COMMON /IXYGD2/  NCSP,        KBLK

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Prior to 7/16/03:
! Now declare variables above the common block in which they belong.
! (bmy, 7/16/03)
!      REAL*8 HMAX ,     HMIN,
!     2     PLOURB,    PLOTROP,   TSPMIDC,
!     3     R1DELT,    DELT,      TIMREMAIN, XELAPS,    TOLD,
!     4     RDELT,     XELAPLAST, SUMAVGE,   SUMAVHI,   SUMRMSE,
!     5     SUMRMHI,   TOTSTEP,   TOTIT,     TELAPS,    RMSERR
!     COMMON /DGEAR/
!    1  HMIN,
!    2  PLOURB,    PLOTROP,   TSPMIDC
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      REAL*8  ::       HMIN,        PLOURB,    PLOTROP,   TSPMIDC     
      COMMON /DGEAR/   HMIN,        PLOURB,    PLOTROP,   TSPMIDC

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Prior to 7/16/03:
! Split /DGEAR2/ into two common blocks, since some of these variables need
! to be held THREADPRIVATE while others do not. (hamid, bmy, 7/16/03)
!      COMMON /DGEAR2/
!     3  HMAX,      R1DELT,    DELT,      TIMREMAIN, XELAPS,    TOLD,
!     4  RDELT,     XELAPLAST, SUMAVGE,   SUMAVHI,   SUMRMSE,
!     5  SUMRMHI,   TOTSTEP,   TOTIT,     TELAPS,    RMSERR
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      ! /DGEAR2/ needs to be held THREADPRIVATE (hamid, bmy, 7/16/03)
      REAL*8  ::       HMAX,        R1DELT,    DELT,      TIMREMAIN
      REAL*8  ::       XELAPS,      TOLD,      RDELT,     XELAPLAST
      REAL*8  ::       RMSERR
      COMMON /DGEAR2/  HMAX,        R1DELT,    DELT,      TIMREMAIN, 
     &                 XELAPS,      TOLD,      RDELT,     XELAPLAST, 
     &                 RMSERR

      ! /DGEAR3/ doesn't need to be held THREADPRIVATE (hamid, bmy, 7/16/03)
      REAL*8  ::       SUMAVGE,     SUMAVHI,   SUMRMSE,   SUMRMHI
      REAL*8  ::       TOTSTEP,     TOTIT,     TELAPS
      COMMON /DGEAR3/  SUMAVGE,     SUMAVHI,   SUMRMSE,   SUMRMHI, 
     &                 TOTSTEP,     TOTIT,     TELAPS
C
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Prior to 7/16/03:
! MONTHP and KYEAR are not needed for SMVGEAR II, so remove these from the
! common block. (bmy, 7/16/03)
!      INTEGER NSFTOT,    NPDTOT,    NSTTOT,
!     3     ISREORD,   IFREORD,
!     4     IFAILTOT,  LFAILTOT,  NFAILTOT,  MONTHP,    KYEAR
!
!      COMMON /IGEAR/
!     2  NSFTOT,    NPDTOT,    NSTTOT,
!     3  ISREORD,   IFREORD,
!     4  IFAILTOT,  LFAILTOT,  NFAILTOT,  MONTHP,    KYEAR
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      INTEGER ::       NSFTOT,      NPDTOT,    NSTTOT,    ISREORD
      INTEGER ::       IFREORD,     IFAILTOT,  LFAILTOT,  NFAILTOT  
      COMMON /IGEAR/   NSFTOT,      NPDTOT,    NSTTOT,    ISREORD,   
     &                 IFREORD,     IFAILTOT,  LFAILTOT,  NFAILTOT

      ! /IGEAR2/ has to be declared THREADPRIVATE (bmy, 7/16/03)
      INTEGER ::       NQQ,         NSUBFUN,   NPDERIV   
      INTEGER ::       NFAIL,       IFAIL,     LFAIL
      COMMON /IGEAR2/  NQQ,         NSUBFUN,   NPDERIV, 
     &                 NFAIL,       IFAIL,     LFAIL

      INTEGER ::       NPHOT,       NPRODLO,   NPRODHI,   MSTEP
      INTEGER ::       MAXORD,      MBETWEEN,  IC3H8,     IC2H6
      COMMON /CHEM2/   NPHOT,       NPRODLO,   NPRODHI,   MSTEP,     
     &                 MAXORD,      MBETWEEN,  IC3H8,     IC2H6

      ! /CHEM2A/ has to be held THREADPRIVATE (bmy, 7/16/03)
      INTEGER ::       ISCHAN,      NOCC,      NFDH3,     NFDL2   
      INTEGER ::       NFDH2,       NFDL1,     NFDH1,     NFDREP
      INTEGER ::       NFDREP1,     NFDL0,     NALLR
      COMMON /CHEM2A/  ISCHAN,      NOCC,      NFDH3,     NFDL2,
     &                 NFDH2,       NFDL1,     NFDH1,     NFDREP, 
     &                 NFDREP1,     NFDL0,     NALLR
                       
      INTEGER ::       NGAS,        NMREAC
      COMMON /CHEM3/   NGAS,        NMREAC

      ! Added NNADDG to /CHEM4/ for DMS+OH+O2 rxn (bdf, bmy, 4/18/03)
      INTEGER ::       NNADD1,      NNADDA,      NNADDB
      INTEGER ::       NNADDC,      NNADDD,      NNADDK
      INTEGER ::       NNADDV,      NNADDZ,      NKO3PHOT
      INTEGER ::       NNADDG,      NEMIS,       NDRYDEP
      INTEGER ::       NKHNO4      
      COMMON /CHEM4/   NNADD1,      NNADDA(ICS), NNADDB(  ICS), 
     &                 NNADDC(ICS), NNADDD(ICS), NNADDK(  ICS), 
     &                 NNADDV(ICS), NNADDZ,      NKO3PHOT(ICS),
     &                 NNADDG(ICS), NEMIS( ICS), NDRYDEP( ICS),
     &                 NKHNO4(ICS)

      INTEGER ::       IH2O,        IOXYGEN,   MB1,      MB2
      COMMON /SPECIES/ IH2O,        IOXYGEN,   MB1,      MB2

      ! Added for interannually-varying Methane (bnd, bmy, 7/1/03)
      INTEGER ::       ICH4
      COMMON /SPECIE2/ ICH4

      ! Added for interannually-varying Methane (bnd, bmy, 7/1/03)
      REAL*8  ::       C3090S,      C0030S,    C0030N,   C3090N
      COMMON /SPECIE3/ C3090S,      C0030S,    C0030N,   C3090N

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Prior to 7/16/03:
! Remove KPHT, KRDD, KMIX, KINS, KGCO -- these are obsolete. (bmy, 7/16/03)
!      INTEGER IOUT,KGLC,  KPHT,  KRDD,  KMIX,  KCPD,
!     2     KINS,KGCO
!      COMMON /FILES/ 
!     1  IOUT,      KGLC,      KPHT,      KRDD,      KMIX,      KCPD,
!     2  KINS,      KGCO
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      INTEGER ::       IOUT,        KGLC,      KCPD,     IO93
      COMMON /FILES/   IOUT,        KGLC,      KCPD,     IO93
C
C *********************************************************************
C *               SET REAL AND INTEGER ARRAY VARIABLES                *
C *********************************************************************
C
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Prior to 7/16/03:
! GEOS-CHEM does not need these variables which compute photorates.  
! GEOS-CHEM photorates are computed either in FAST-J or SLOW-J. (bmy, 7/16/03)
!      REAL*8 DECLIN,RAGSUT,SINDEC,COSDEC
!      COMMON /IMXDAY/
!     1  DECLIN(MAXDAYS),   RAGSUT(MAXDAYS),  SINDEC(MAXDAYS),
!     2  COSDEC(MAXDAYS)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C
      INTEGER ::        JLOWVAR,           KTLPVAR
      INTEGER ::        JLOFIXED,          JHIFIXED
      COMMON /IMXBLOCK/ JLOWVAR( MXBLOCK), KTLPVAR( MXBLOCK), 
     &                  JLOFIXED(MXBLOCK), JHIFIXED(MXBLOCK)
C
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Prior to 7/16/03:
! GEOS-CHEM implementation of SMVGEAR II does not need these variables,
! so comment out this common block. (bmy, 7/16/03)
!      REAL*8 SIGMAL,PRESSL,RHOA,DSIG_SMV
!      COMMON /DILAYER/
!     1  DSIG_SMV( ILAYER), SIGMAL(   ILAYER),
!     2  PRESSL(   ILAYER), RHOA(     ILAYER)
!C
! Prior to 7/16/03:
! GEOS-CHEM implementation of SMVGEAR II does not need these variables,
! so comment out this common block. (bmy, 7/16/03)
!      REAL*8 TEMPL,VMET,SIGDIF,TMORN,PRESSC
!      COMMON /DIVERT/
!     1  TEMPL(     IVERT), VMET(      IVERT), SIGDIF(   IVERT), 
!     2  TMORN(     IVERT), PRESSC(    IVERT)    
!C
! Prior to 7/16/03:
! GEOS-CHEM implementation of SMVGEAR II does not need these variables,
! so comment out this common block. (bmy, 7/16/03)
!      REAL*8 XLAT,XLON,DMERIDUT,GRIDAREA,DSX,XLONUT
!      REAL*8 DSY,SINXLAT,COSXLAT,HMETT,HMET1,HMET2,RSET,RRIS,TZDIF
!      REAL*8 ZENRAT0,ZENRAT1
!      COMMON /DIMLOOP/
!     2  XLAT(     IMLOOP), XLON(     IMLOOP), DMERIDUT(IMLOOP), 
!     3  GRIDAREA( IMLOOP), DSX(      IMLOOP), XLONUT(  IMLOOP), 
!     4  DSY(      IMLOOP), SINXLAT(  IMLOOP), COSXLAT( IMLOOP), 
!     5  HMETT(    IMLOOP), HMET1(    IMLOOP), HMET2(   IMLOOP),
!     6  RSET(     IMLOOP), RRIS(     IMLOOP), TZDIF(   IMLOOP), 
!     7  ZENRAT0(  IMLOOP), ZENRAT1(  IMLOOP)  
!C
! Prior to 7/16/03:
! The GEOS-CHEM implementation of SMVGEAR II does not use COLENG, AERSURF,
! VHMET1, VHMET, VMET3, RHO3, GRIDVH, and CSUMA1, so remove them from the 
! common block.  (bmy, 7/16/03)
!      REAL*8 RHO3,VHMET1,VHMET,VMET3,CSUMA1,CSUMA,CSUMC,GRIDVH
!      REAL*8 COLENG,ERRMX2,AERSURF
!      COMMON /DITLOOP/
!     1  RHO3(     ITLOOP), VHMET1(   ITLOOP), VHMET(   ITLOOP), 
!     2  VMET3(    ITLOOP),
!     3  CSUMA1(   ITLOOP), CSUMA(    ITLOOP), CSUMC(   ITLOOP),  
!     4  GRIDVH(   ITLOOP), COLENG(   ITLOOP), ERRMX2(  ITLOOP),
!     5  AERSURF(  ITLOOP)   
!
! Prior to 7/16/03:
! Now made these arrays allocatable in "comode_mod.f".  This will save
! memory for non-fullchem runs. (bmy, 7/18/03)
!      REAL*8 ::        CSUMA,         CSUMC,         ERRMX2
!      COMMON /DITLOOP/ CSUMA(ITLOOP), CSUMC(ITLOOP), ERRMX2(ITLOOP)   
!
! Prior to 7/16/03:
! The GEOS-CHEM implementation of SMVGEAR II does not use MLOPJ and
! REORDER_SAVE, so remove them from the common block. (bmy, 7/16/03) 
!      INTEGER MLOPJ,JREORDER,LREORDER,ITWO,REORDER_SAVE
!     3       ,NCSLOOP
!      COMMON /IITLOOP/
!     1  MLOPJ(    ITLOOP), JREORDER( ITLOOP), LREORDER( ITLOOP), 
!     2  ITWO(     ITLOOP), REORDER_SAVE(ITLOOP,24)
!     3 ,NCSLOOP(ITLOOP,ICS)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      INTEGER ::       JREORDER,         LREORDER 
      INTEGER ::       ITWO,             NCSLOOP
      COMMON /IITLOOP/ JREORDER(ITLOOP), LREORDER(ITLOOP), 
     &                 ITWO(    ITLOOP), NCSLOOP( ITLOOP,ICS)
C
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Prior to 7/16/03:
! The GEOS-CHEM implementation of SMVGEAR II does not use RHO3K and GRIDVHK,
! so remove them from the common block. (bmy, 7/16/03) 
!      REAL*8 DENAIR,CONCO2,CONCN2,T3I,TEMP1,RHO3K,T3K,GRIDVHK
!      REAL*8 PRESSK,DELY,ERRHOLD
!      COMMON /DKBLOOP/
!     1  DENAIR(   KBLOOP), CONCO2(   KBLOOP), CONCN2(   KBLOOP),
!     2  T3I(      KBLOOP), TEMP1(    KBLOOP), RHO3K(    KBLOOP),    
!     3  T3K(      KBLOOP), GRIDVHK(  KBLOOP),  
!     4  PRESSK(   KBLOOP), DELY(     KBLOOP), ERRHOLD(  KBLOOP)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      REAL*8 ::        DENAIR,         CONCO2,         CONCN2
      REAL*8 ::        T3I,            TEMP1,          T3K
      REAL*8 ::        PRESSK,         DELY,           ERRHOLD
      COMMON /DKBLOOP/ DENAIR(KBLOOP), CONCO2(KBLOOP), CONCN2( KBLOOP),
     &                 T3I(   KBLOOP), TEMP1( KBLOOP), T3K(    KBLOOP),  
     &                 PRESSK(KBLOOP), DELY(  KBLOOP), ERRHOLD(KBLOOP)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Prior to 7/16/03:
! The GEOS-CHEM implementation of SMVGEAR II does not need XELRAT, T1BEG,
! T2BEG, T1FIN, T2FIN, so remove them from the common block. (bmy, 7/16/03)
!      REAL*8 XELRAT,T1BEG,T2BEG,YABST,T1FIN,T2FIN,AREAXT
!      COMMON /DKBLOOP5/
!     1  XELRAT(   KBLOOP), T1BEG(    KBLOOP), T2BEG(    KBLOOP),    
!     2  YABST(    KBLOOP), T1FIN(    KBLOOP), T2FIN(    KBLOOP),
!     3  AREAXT(   KBLOOP) 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      REAL*8 ::         YABST
      COMMON /DKBLOOP5/ YABST(KBLOOP)
C
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Prior to 7/16/03:
! IABOVK is only ever used in "smvgear.f"; therefore declare it as a local 
! variable in "smvgear.f" and remove this common block. (bmy, 7/16/03)
!      INTEGER IABOVK
!      COMMON /IKBLOOP2/
!     1  IABOVK(   KBLOOP) 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C
      ! Add NKSPECG for DMS+OH+O2 rxn (bdf, bmy, 4/18/03)
      INTEGER NMOTH,NTSPEC,JPHOTRAT,ISGAINR,ISPORL,NOGAINE,NOUSE
      INTEGER NSPEC,NTRATES,ISGAINE,NTLOOPNCS,NSPCSOLV,ISCHANG,NRATES
      INTEGER NM3BOD,ITWOR,ITHRR,INOREP,NRATCUR,NSURFACE,NPRESM,NMAIR
      INTEGER NMO2,NMN2,NNEQ,NARR,NABR,NACR,NABC,NKSPECW,NKSPECX
      INTEGER NKSPECY,NKSPECZ,NKSPECV,ISLOSSR,NKSPECA,NKSPECB,NKSPECC
      INTEGER NKSPECD,NKSPECK,NKSPECG
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
     5  NKSPECG(MAXGL2,ICS)

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

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Prior to 7/16/03:
! IPORD is obsolete, so remove this common block. (bmy, 7/16/03)
!      INTEGER IPORD
!      COMMON /IIPHOT/
!     1  IPORD(   IPHOT8)    
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      REAL*8  ::        ARRT,           BRRT
      REAL*8  ::        FCVT,           FCT1T,           FCT2T
      COMMON /DMXCOF/   ARRT(MXCOF),    BRRT( MXCOF),   
     &                  FCVT(MXCOF),    FCT1T(MXCOF),    FCT2T(MXCOF)

      INTEGER ::        KCRRT   
      COMMON /IMXCOF/   KCRRT(MXCOF)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Prior to 7/16/03:
! The GEOS-CHEM implementation of SMVGEAR II does not need PRATE, 
! so comment out this common block.  (bmy, 7/16/03)
!      REAL*8 PRATE
!      COMMON /DPHOT2/
!     1  PRATE(IVERT,IPHOT) 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C
      INTEGER NKARR,NKABR,NKACR,NKABC,IRORD
      COMMON /INMRAT2/
     1  NKARR(NMTRATE,ICS),  NKABR(NMTRATE,ICS), NKACR(NMTRATE,ICS), 
     2  NKABC(NMTRATE,ICS),  IRORD(NMTRATE,ICS)
C
      REAL*8 ARR,BRR,FCV,FCTEMP1,FCTEMP2
      COMMON /DNMTRATE/
     1  ARR(    NMTRATE, ICS),  BRR(    NMTRATE, ICS), 
     2  FCV(    NMTRATE, ICS),  FCTEMP1(NMTRATE, ICS),
     3  FCTEMP2(NMTRATE, ICS)  
C
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
     8  IRMA(   NMTRATE      ),  IRMB(    NMTRATE     ),
     9  KCRR(   NMTRATE,  ICS),  LSKIP(   MXRATE,  ICS),
     1  IRMC(   NMTRATE      ),  JPHOTNK( NMTRATE, ICS),
     2  IUSED(  MXRATE,   ICS)
C
      INTEGER ::       NEWNK
      COMMON /IMAXGL3/ NEWNK(MAXGL)
C
      REAL*8  ::       FRACP 
      COMMON /DMAXGL2/ FRACP(MAXGL, ICS)
C
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
C
      INTEGER ::        MBCOMP,            MBTRACE
      COMMON /IIMASBAL/ MBCOMP(IMASBAL,2), MBTRACE(IMASBAL)
C
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Prior to 7/16/03:
! PRESS5KM is not used anywhere, so delete the common block (bmy, 7/16/03)
!      REAL*8 PRES5KM
!      COMMON /DIALTS/
!     1  PRES5KM(  IALTS)
!C
! Prior to 7/16/03:
! CINIT, PRATK1, PRATKD are not needed in the GEOS_CHEM implementation of 
! SMVGEAR II, so remove them from the common block (bmy, 7/16/03)
!      REAL*8 CINIT,CNEW,CEST,GLOSS,CHOLD,VDIAG,CBLK,DTLOS,EXPLIC,CONC
!      REAL*8 RRATE,URATE,TRATE,PRATK1,PRATKD,CORIG
!      COMMON /DKBLOOP2/
!     1  CINIT(  KBLOOP,  MXGSAER),  
!     2  CNEW(   KBLOOP,  MXGSAER),  
!     3  CEST(   KBLOOP,  MXGSAER),  
!     4  GLOSS(  KBLOOP,  MXGSAER),  
!     5  CHOLD(  KBLOOP,  MXGSAER), 
!     6  VDIAG(  KBLOOP,  MXGSAER),  CBLK(  KBLOOP,MXGSAER),  
!     7  DTLOS(  KBLOOP,  MXGSAER),  EXPLIC(KBLOOP,MXGSAER),
!     1  CONC(   KBLOOP,MXGSAER*7),  
!     2  RRATE(  KBLOOP,  NMTRATE),  
!     3  URATE(  KBLOOP,NMTRATE,3), 
!     4  TRATE(  KBLOOP,NMTRATE*2), 
!     5  PRATK1( KBLOOP,    IPHOT),  
!     6  PRATKD( KBLOOP,    IPHOT), 
!     7  CORIG(  KBLOOP,  MXGSAER)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      REAL*8 CNEW,CEST,GLOSS,CHOLD,VDIAG,CBLK,DTLOS,EXPLIC,CONC
      REAL*8 RRATE,URATE,TRATE,PRATK1,PRATKD,CORIG
      COMMON /DKBLOOP2/
     2  CNEW(   KBLOOP,  MXGSAER),  
     3  CEST(   KBLOOP,  MXGSAER),  
     4  GLOSS(  KBLOOP,  MXGSAER),  
     5  CHOLD(  KBLOOP,  MXGSAER), 
     6  VDIAG(  KBLOOP,  MXGSAER),  CBLK(  KBLOOP,MXGSAER),  
     7  DTLOS(  KBLOOP,  MXGSAER),  EXPLIC(KBLOOP,MXGSAER),
     1  CONC(   KBLOOP,MXGSAER*7),  
     2  RRATE(  KBLOOP,  NMTRATE),  
     3  URATE(  KBLOOP,NMTRATE,3), 
     4  TRATE(  KBLOOP,NMTRATE*2), 
     7  CORIG(  KBLOOP,  MXGSAER)
C
      REAL*8 ::         CC2
      COMMON /DKBLOOP0/ CC2(KBLOOP,0:MXARRAY)  
C
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Prior to 7/16/03:
! The GEOS-CHEM implementation of SMVGEAR II does not need C.  
! Comment out the common block. (bmy, 7/16/03)
!      REAL*8 C
!      COMMON /DITLOOP3/
!     1  C(      ITLOOP,   IGAS)
!C
! Prior to 7/16/03:
! KGRP is only ever used w/in "smvgear.f", so declare it there as a local
! variable and remove it from "comode.h". (bmy, 7/16/03)
!      INTEGER KGRP
!      COMMON /IKBLOOP/
!     1  KGRP(   KBLOOP,        5)
!C
! Prior to 7/16/03:
! The GEOS-CHEM implementation of SMVGEAR II does not need FIELDXY, FIELDYZ,
! and FIELDXZ.  Comment out the common block. (bmy, 7/16/03)
!      REAL*8 FIELDXY,FIELDYZ,FIELDXZ
!      COMMON /FIELD/
!     1  FIELDXY(ILONG, ILAT), FIELDYZ(ILAT,IVERT),
!     2  FIELDXZ(ILONG,IVERT)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C
      INTEGER ::      MLOP,             JLOP_SMV  
      COMMON /IILAT2/ MLOP(ILAT,ILONG), JLOP_SMV(ILAT,ILONG,ILAYER)
C
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Prior to 7/16/03:
! The GEOS-CHEM implementation of SMVGEAR II does not need RATMIX and GQSCHEM.
! Comment out the common block. (bmy, 7/16/03)
!      REAL*8 RATMIX,GQSCHEM
!      COMMON /DIGAS2/
!     1  RATMIX(       IGAS,  IALTS), GQSCHEM(IGAS,ICS)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C
      INTEGER NKPHOTRAT,NPPHOTRAT,NKNPHOTRT
      COMMON /DIPHOT2/ 
     1  NKPHOTRAT(IPHOT,ICS),      NPPHOTRAT(IPHOT,ICS),  
     2  NKNPHOTRT(IPHOT,ICS) 
C
      REAL*8 ::        FRACGAIN,              QBKCHEM
      COMMON /DIMXGS2/ FRACGAIN(MXGSAER,ICS), QBKCHEM( MXGSAER,ICS) 
C
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
C
      INTEGER JZILCH,KZILCH,MZILCH
      COMMON /IGMXGLS/
     &     JZILCH(MXGSAER),  KZILCH(MXGSAER),  MZILCH(MXGSAER) 
C
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
C
      REAL*8 ::        FRACNFR,           FRACPL
      COMMON /DMXCOUN/ FRACNFR(MXCOUNT4), FRACPL(MXCOUNT2) 
C
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
C
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
      COMMON /DMISC/     ASET(10,8), PINP(20), CVAR(15), O3DOBS(12)  

      REAL*8 ::          ENQQ2,         ENQQ3,          CONPST
      REAL*8 ::          ENQQ1,         CONP15
      COMMON /IORDR/     ENQQ2(MORDER), ENQQ3( MORDER), CONPST(MORDER),
     &                   ENQQ1(MORDER), CONP15(MORDER)

      REAL*8 ::          PERTS2,           PERTST
      COMMON /DMORD/     PERTS2(MORDER,3), PERTST(MORDER,3)
C
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Prior to 7/16/03:
! We don't need LDMONTH, ININT, so remove the common block.  
! ININT may be declared local w/in "readchem.f" (bmy, 7/16/03)
!      INTEGER LDMONTH,ININT
!      COMMON /IMISC2/
!     1  LDMONTH(  12),   ININT(    10)   
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      INTEGER ::         JLLOW,          KLHI
      COMMON /IMXBLOCK/  JLLOW(MXBLOCK), KLHI(MXBLOCK)

      REAL*8 ::          XGDFCF, ASTKCF, RUARSL, RH100
      COMMON /XARSOL/    XGDFCF, ASTKCF, RUARSL, RH100

      INTEGER ::         IARSFA, MWARSL, MNTHARS
      COMMON /IARSOL/    IARSFA, MWARSL, MNTHARS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Prior to 7/16/03:
! ABSHUMK is only ever used w/in "calcrate.f".  Remove the /DKBLOOP4/ common
! block and THREADPRIVATE declaration, and declare ABSHUMK as a local variable
! w/in "calcrate.f" (bmy, 7/16/03)
!      REAL*8 ::         ABSHUMK
!      COMMON /DKBLOOP4/ ABSHUMK(KBLOOP)
!
! Prior to 7/16/03:
! These are not used anywhere, so remove the common block (bmy, 7/16/03)
!      INTEGER MROTAT1,MINROT1,NUMSUBS,LSPECEMIS,MROTAT2,MINROT2,
!     1        MAXPOS,NOGAINR,NOLOSSR,MAXSTEPS
!      COMMON /IDICS0/
!     1  MROTAT1( ICS),  MINROT1( ICS), NUMSUBS(ICS), LSPECEMIS(ICS), 
!     2  MROTAT2( ICS),  MINROT2( ICS), MAXPOS( ICS), 
!     3  NOGAINR(ICS), NOLOSSR(  ICS),
!     4  MAXSTEPS(ICS)
C
! Prior to 7/16/03:
! These are not used anywhere, so remove the common block (bmy, 7/16/03)
!      ! Split these off from common block /IDICS0/ to avoid mixing
!      ! INTEGER and REAL types in the same common block (bmy, 4/20/99)
!      REAL*8 YLOW,HMAXDAY
!      COMMON /IDICSR/
!     1  YLOW(   ICS),  
!     2  HMAXDAY(  ICS)
C
! Prior to 7/16/03:
! ICLO, JCLO are not used; so we can delete the common block (bmy, 7/16/03)
!      INTEGER ICLO,JCLO
!      COMMON /IICP0J/   ICLO(     ICP),  JCLO(    ICP) 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      INTEGER ::         NKEMIS,              NTEMIS
      INTEGER ::         NKDRY,               NTDEP
      COMMON /IMAXGL4/   NKEMIS(MAXGL3,ICS),  NTEMIS(MAXGL3,ICS),
     &                   NKDRY (MAXGL3,ICS),  NTDEP( MAXGL3)
C
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Prior to 7/16/03:
! FIELD1 is not used, so we can delete the common block. (bmy, 7/16/03)
!      ! Split FIELD1 off into a separate common block, to avoid mixing
!      ! INTEGER and REAL types (bmy, 4/20/99)
!      REAL*8 FIELD1
!      COMMON /IILAT3/ FIELD1(ILONG, ILAT)
!
! Prior to 7/16/03:
! None of these are needed for SMVGEAR II, remove the common block.
! (bmy, 7/16/03)
!      INTEGER MZLO1,MZLO2,MZHI0,MZHI1,KZLO1,
!     1        KZLO2,KZHI0,KZHI1,IHIZ1,IHIZ2,IHIZ3
!      COMMON /JGMXGS2/ 
!     2  MZLO1( MXGSAER,ICP), MZLO2( MXGSAER,ICP), MZHI0( MXGSAER,ICP),
!     3  MZHI1( MXGSAER,ICP), KZLO1( MXGSAER,ICP), KZLO2( MXGSAER,ICP),
!     4  KZHI0( MXGSAER,ICP), KZHI1( MXGSAER,ICP),
!     5  IHIZ1(  MXGSAER,MXGSAER,ICP), IHIZ2(   MXGSAER,MXGSAER,ICP),
!     6  IHIZ3(  MXGSAER,MXGSAER,ICP)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Prior to 7/16/03:
! The GEOS-CHEM implementation of SMVGEAR II does not require QPRODA, QPRODB,
! QPRODC, QPRODD, QPROD.  Comment out this common block. (bmy, 7/16/03) 
!      ! Split these off from common block /IMXCOUN/ to avoid mixing
!      ! INTEGER and REAL*8 types in the same common block (bmy, 4/20/99)
!      REAL*8 QPRODA,QPRODB,QPRODC,QPRODD,QPROD
!      COMMON /IMXCOUN2/
!     &  QPRODA(MXCOUNT4),    QPRODB(MXCOUNT4),    QPRODC(MXCOUNT4   ),  
!     &  QPRODD(MXCOUNT4),    QPROD(MXGSAER,MAXGL,ICS)
!C
! Prior to 7/16/03:
! The GEOS-CHEM implementation of SMVGEAR II does not require CINP.
! Comment out this common block. (bmy, 7/16/03)
!      ! Split this off from common block /IMISC/ to avoid mixing
!      ! INTEGER and REAL*8 types in the same common block (bmy, 4/20/99)
!      REAL*8 CINP      
!      COMMON /IMISC2/ CINP(NMRPROD)
!C
! Prior to 7/16/03:
! The GEOS-CHEM implementation of SMVGEAR II does not require NUMSDT or
! NKSDT.  Comment out this common block. (bmy, 7/16/03)
!      INTEGER NUMSDT,NKSDT
!      COMMON /IDIC20/ NUMSDT(3, ICS,2),  NKSDT( 3,MXRATE,ICS)
!C
! Prior to 7/16/03:
! Now add this to the /FILES/ common block (bmy, 7/16/03)
!      INTEGER IO93
!      COMMON /IOUNIT/ IO93
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C
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
!$OMP THREADPRIVATE( /DKBLOOP/  )
!$OMP THREADPRIVATE( /DKBLOOP0/ )
!$OMP THREADPRIVATE( /DKBLOOP2/ )
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Prior to 7/16/03:
! Remove THREADPRIVATE declarations for /DKBLOOP4/, since these contained 
! variables which are now declared locally w/in "calcrate.f"
!!$OMP THREADPRIVATE( /DKBLOOP4/ )
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!$OMP THREADPRIVATE( /DKBLOOP5/ )
!$OMP THREADPRIVATE( /IGEAR2/   )
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Prior to 7/16/03:
! Remove THREADPRIVATE declarations for IKBLOOP and IKBLOOP2, since these 
! contained variables which are now declared locally w/in "smvgear.f"
!!$OMP THREADPRIVATE( /IKBLOOP/  )
!!$OMP THREADPRIVATE( /IKBLOOP2/ )
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!$OMP THREADPRIVATE( /IXYGD2/   )
C
C *********************************************************************
C ****************** END OF COMMON BLOCK COMODE.H *********************
C *********************************************************************
C
