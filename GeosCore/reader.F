! $Id: reader.f,v 1.2 2010/03/15 19:33:21 ccarouge Exp $
      SUBROUTINE READER( FIRSTCHEM, am_I_Root )
!
!******************************************************************************
!  Subroutine READER reads on/off switches and other settings for SMVGEAR II.
!  (M. Jacobson 1997; bdf, bmy, 4/18/03, 10/16/06)
!
!  NOTES:
!  (1 ) Now force double-precision values with the "D" exponent.  Also use
!        consistent physical constant values w/ GEOS-CHEM.  Now use GEOS-CHEM
!        unit IU_FILE number to read the "mglob.dat" file.  Now references
!        GEOS_CHEM_STOP from "error_mod.f".  Now force double-precision with
!        the "D" exponent.  Set KGLC = IU_CHEMDAT = 7 from "file_mod.f" 
!        (bmy, 4/18/03)
!  (2 ) Remove obsolete variables AERSURF, MLOPJ, REARTHC, DENCONS, HALFDAY, 
!        GRAVC, FOURPI, TWOPI, REARTH, RPRIMB, AVOG1, HALF, THIRD, THRPI2, 
!        PID180, PID2, SCTWOPI, AMRGAS, TWPISC, REARTH. these aren't used w/in 
!        "reader.f" anymore.  Use F90-style variable declarations.  Also 
!        remove obsolete variables from documentation. (bmy, 7/16/03)
!  (3 ) Redefine CHEMINTV [s] to the value in "input.geos" so that we don't 
!        have a discrepancy with the value in "mglob.dat".  SLOW-J is now 
!        obsolete; remove LSLOWJ #ifdef blocks (bmy, 6/23/05)
!  (4 ) Physical constants and some error tolerances are now defined as 
!        parameters in "comode.h".  In this way, their values will be defined 
!        before the first call to READER for the offline aerosol simulation. 
!        (bec, bmy, 3/29/06)
!  (5 ) Increase max # of products that a reaction can have from 12 to 14.  
!        This coincides w/ the new globchem.dat. (bmy, 8/9/06)
!  (6 ) At the end of this subrouitne, now set NCS=NCSURBAN (=1) instead of 
!        hardwiring it. (dbm, bmy, 10/16/06)
!  01 Aug 2012 - R. Yantosca - Add reference to findFreeLUN from inqure_mod.F90
!  03 Aug 2012 - R. Yantosca - Move calls to findFreeLUN out of DEVEL block
!  07 Aug 2012 - R. Yantosca - Now print LUN used to open file
!  19 Mar 2013 - R. Yantosca - When using ESMF interface to GEOS-5, append
!                              ".rc" to input.geos (instead of ___.rc)
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!******************************************************************************
!
      USE ERROR_MOD,  ONLY : GEOS_CHEM_STOP
      USE TIME_MOD,   ONLY : GET_TS_CHEM
      USE inquireMod, ONLY : findFreeLUN
      USE CMN_SIZE_MOD                       ! Size parameters
      USE COMODE_LOOP_MOD                    ! SMVGEAR II arrays
      USE CMN_GCTM_MOD                       ! Re, PI

      IMPLICIT NONE

!
! *********************************************************************
! *  THIS SUBROUTINE OPENS ALL DATA FILES AND READS DATA FROM m.dat ***
! *  FOR DEFINITIONS OF THE PARAMETERS READ IN HERE, SEE define.dat ***
! *********************************************************************
!
!        RRRRRR   EEEEEEE     A     DDDDDD   EEEEEEE  RRRRRR 
!        R     R  E          A A     D    D  E        R     R 
!        RRRRRR   EEEEE     A   A    D    D  EEEEE    RRRRRR
!        R  R     E        AAAAAAA   D    D  E        R  R  
!        R   R    EEEEEEE A       A DDDDDD   EEEEEEE  R   R
!
!
! *********************************************************************
! *              NAMELIST DATA FOR DATA FILE m.dat                    *
! *********************************************************************
!
! *********************************************************************
!                           MAIN SWITCHES
! *********************************************************************
! IFSOLVE   = 1: SOLVE CHEMICAL EQUATIONS WITH SMVGEAR 
!             0: DO NOT SOLVE ANY CHEMICAL EQUATIONS (mglob.dat)
! IFPRAT    = 1: USE DEFAULT PHOTORATES FROM photrate.dat;
!           = 0: USE DEFAULT PHOTORATES FROM globchem.dat
! INCVMIX   = 1: INTERPOLATE MIXING RATIO PROFILES FROM DATA IN MIXRATIO.DAT
! ITESTGEAR = 1: CREATE EXACT SOLUTION TO COMPARE OTHER GEAR SOLUTIONS AGAINST
!           = 2: COMPARE CURRENT SOLUTION TO EXACT SOLUTION
!
! IFURBAN  IFTROP  IFSTRAT            TYPE OF CHEMISTRY SOLVED  
!                             (U=URBAN, T=TROPOSPHERIC, S=STRATOSPHERIC)
! ----------------------------------------------------------------------
!    0        0       0      DO NOT SOLVE CHEMISTRY 
!    1        0       0      SOLVE U EVERYWHERE 
!    0        1       0      SOLVE T EVERYWHERE 
!    0        0       1      SOLVE S EVERYWHERE 
!    1        1       1      SOLVE U BELOW PLOURB, T BETWEEN PLOURB,
!                            PLOTROP, AND S ABOVE PLOTROP  
!    0        2       2      SOLVE T/S   CHEMISTRY EVERYWHERE   
!    2        2       2      SOLVE U/T/S CHEMISTRY EVERYWHERE   
!
      LOGICAL, INTENT(IN) :: FIRSTCHEM
      LOGICAL, INTENT(IN) :: am_I_Root   ! Is this the root CPU?
      INTEGER             :: K,       M2,      M1,       MLOOP,   KLOOP 
      INTEGER             :: JLOOP,   IAVBLOK, IAVGSIZE, IREMAIN, JADD
      INTEGER             :: IFCHEM,  I,       NALLREAC, NMPROD,  I1
      INTEGER             :: J,      NK
!      INTEGER             :: IU_CHEMDAT
      INTEGER             :: IU_FILE
      INTEGER             :: IU_SMV2LOG

      REAL*8              :: ERRMAXU, YLOWU, YHIU, HMAXDAYU
      REAL*8              :: ERRMAXR, YLOWR, YHIR, HMAXDAYR
      REAL*8              :: ERRMAXS, YLOWS, YHIS, HMAXDAYS
      REAL*8              :: ABHI,    ABLO
      CHARACTER(LEN=255)  :: FILENAME


      NAMELIST /CTLFLG/ IFSOLVE, ITESTGEAR,   
     1                  IFURBAN, IFTROP,  IFSTRAT
!
! *********************************************************************
!                            DIMENSIONS
! *********************************************************************
! NLAT      = # SOUTH-NORTH GRID CELLS
! NLONG     = # WEST-EAST GRID CELLS
! NVERT     = # VERTICAL LAYERS
! KULOOP    = MAXIMUM ACTUAL # OF GRID CELLS IN A GRID BLOCK
! LYOUT     = SPECIFIC SOUTH-NORTH CELL FOR PRINTING
! LXOUT     = SPECIFIC WEST-EAST CELL FOR PRINTING
! LZOUT     = SPECIFIC VERTICAL LAYER FOR PRINTING
!
! Now include NREAD (the number of entries in globchem.dat)
! which had been increased to allow for more species by FP
! (hotp 8/1/09)
      NAMELIST /CTLDIM/ KULOOP,   
     1                  LYOUT, LXOUT,  LZOUT, 
     1                  NREAD
!
! *********************************************************************
!            SWITCHES FOR TIME, TIME-STEPS, AND OUTPUT
! *********************************************************************
! CHEMINTV   = TIME STEP FOR GAS AND RADIATIVE PROCESS CALCULATIONS  
!
      NAMELIST /CTLTIM/ CHEMINTV
!
! *********************************************************************
!                         SWITCHES FOR OUTPUT
! *********************************************************************
! IPRATES    = 1: PRINT CHEMICAL RATE COEFFICIENT DATA IN UPDATE.F 
! IPREADER   = 1: PRINT INPUT DATA READ IN READER.F
! IOREAC     = 1: PRINT LIST OF REACTIONS IN READCHEM.F
! APGASA..H  = GASES FOR WHICH OUTPUT ARE PRINTED. OVERRIDES IPRMANY
!
      NAMELIST /CTLPRT/ IPRATES,  IPREADER, 
     1                  IOSPEC,   IOREAC,         
     3                  APGASA,   APGASB,   APGASC,      
     4                  APGASD,   APGASE,   APGASF,   
     5                  APGASG,   APGASH
!
! *********************************************************************
!                      SWITCHES FOR CHEMISTRY
! *********************************************************************
! IFREORD    = 1: REORDER GRID CELLS BY STIFFNESS DURING CHEMISTRY
! FRACDEC    = FRACTION THE TIME STEP IS DECREASED IF CONVERGENCE FAILS
! PLOTROP    = PRESSURE (MB) ABOVE WHICH STRAT CHEM IS SOLVED   
! PLOURB     = PRESSURE (MB) BELOW WHICH URBAN CHEM IS SOLVED  
! ERRMAXU    = RELATIVE ERROR TOLERANCE (FRACTION) FOR URBAN CHEMISTRY 
! ERRMAXR    = RELATIVE ERROR TOLERANCE (FRACTION) FOR TROPOSPHERIC CHEMISTRY 
! ERRMAXS    = RELATIVE ERROR TOLERANCE (FRACTION) FOR STRATOSPHERIC CHEMISTRY 
! YLOWU,YHIU = LOW /HIGH ABS. ERROR TOLERANCES (MOLEC. CM-3) FOR URBAN CHEM 
! YLOWR,YHIR = LOW /HIGH ABS. ERROR TOLERANCES (MOLEC. CM-3) FOR TROP. CHEM 
! YLOWS,YHIS = LOW /HIGH ABS. ERROR TOLERANCES (MOLEC. CM-3) FOR STRAT. CHEM 
! HMAXDAYU   = MAXIMUM TIME STEP FOR DAYTIME URBAN CHEMISTRY (S)
! HMAXDAYR   = MAXIMUM TIME STEP FOR DAYTIME TROP. CHEMISTRY (S)
! HMAXDAYS   = MAXIMUM TIME STEP FOR DAYTIME STRAT. CHEMISTRY (S)
! HMAXNIT    = MAXIMUM TIME STEP FOR NIGHTTIME CHEMISTRY EVERYWHERE
!
      NAMELIST /CLGEAR/ IFREORD,  FRACDEC,  
     2                  PLOURB,   PLOTROP, 
     3                  ERRMAXU,  YLOWU, YHIU,    HMAXDAYU,  
     4                  ERRMAXR,  YLOWR, YHIR,    HMAXDAYR,
     5                  ERRMAXS,  YLOWS, YHIS,    HMAXDAYS,  
     8                  HMAXNIT   
!
! *********************************************************************
! *********************** OPEN CONTROL INPUT FILE ********************* 
! *********************************************************************
!

      ! Find a free file LUN
      IU_FILE  = findFreeLUN()

      ! Name of file to open
      FILENAME = "mglob.dat"

#if defined( ESMF_ )
      ! When connecting to GEOS-5 via ESMF, we need to append '.rc'
      ! to "mglob.dat" so that GEOS-5 will copy the file to the proper
      ! run directory (bmy, 3/19/13)
      FILENAME = TRIM( FILENAME ) // ".rc"
#endif

      OPEN( UNIT=IU_FILE, FILE=TRIM( FILENAME ), ACTION="read",
     &      FORM="formatted", STATUS="old", ACCESS="sequential" )
  
      READ( IU_FILE, 100    ) HEADING 
      READ( IU_FILE, 100    ) COMMENT
      READ( IU_FILE, CTLFLG )
      READ( IU_FILE, CTLDIM )
      READ( IU_FILE, CTLTIM )
      READ( IU_FILE, CTLPRT )
      READ( IU_FILE, CLGEAR )
      CLOSE( UNIT=IU_FILE )

      IF ( am_I_Root ) THEN 
         WRITE( 6, 300 ) TRIM( FILENAME )
 300     FORMAT( 'GEOSCHEM::READER: Successfully read ' , a )
      ENDIF

      ! NOTE: Redefine CHEMINTV [s] to the value in "input.geos" so 
      ! that we don't have a discrepancy with the value in "mglob.dat" 
      ! (bmy, 5/10/05)
      CHEMINTV = GET_TS_CHEM() * 60d0
!
! *********************************************************************
! *                     DEFINE SOME GRID PARAMETERS                   *
! *********************************************************************
! NLOOP     = NUMBER OF GRID-CELLS IN A VERTICAL LAYER
! NTLOOP    = NUMBER OF GRID-CELLS IN THE ENTIRE GRID-DOMAIN
! NLAYER    = NVERT + 1
! LX,Y,ZOUT = IDENTIFY GRID POINT WHERE OUTPUT IS PRINTED
! INCVMIX   = 1: INITIALIZE MIXING RATIOS FROM mixratio.dat
! IFPRAT    = 1: USE DEFAULT PHOTORATES FROM photrate.dat 
! ICOORD    = 1: RECTANGULAR; 2: SPHERICAL; 3: GLOBAL SPHERICAL 
! IFBOX     = 1: SETS UP BOX MODEL TO SOLVE URBAN/TROP/STRAT CHEM TOGETHER  
!                USING DEFAULT PHOTORATES
! ITESTGEAR = 1: SETS UP BOX MODEL TO COMPARE URBAN/TROP/STRAT  
!                CHEMISTRY TO EXACT SOLUTION
!           = 2: SETS UP BOX MODEL TO CREATE URBAN/TROP/STRAT 
!                CHEMISTRY EXACT SOLUTION
!
      IF (ITESTGEAR.GT.0) THEN
         NLAT    = 1 
         NLONG   = 1 
         NVERT   = 1 
         ICOORD  = 1 
         LXOUT   = 1
         LYOUT   = 1
         LZOUT   = 1
      ENDIF 
!
      ! nlat and nlong are defined in chemdr.f (bdf, 4/1/03)
      !NLOOP     = NLAT  * NLONG
      !NTLOOP    = NLOOP * NVERT

      ! needed in reader.f for kuloop (bdf, 4/1/03)
      NTLOOP    = IIPAR*JJPAR*NVERT  
!
      NLAYER    = LLTROP
      LXOUT     = MIN0(LXOUT,NLONG)
      LYOUT     = MIN0(LYOUT,NLAT)
      LZOUT     = MIN0(LZOUT,NVERT)
!
! *********************************************************************
!                            OPEN MORE FILES 
! *********************************************************************
!
!----> Transferred to READCHEM.F
!      ! Find a free file LUN
!      IU_CHEMDAT = findFreeLUN()
!
!      ! Name of file to open
!      FILENAME   ='globchem.dat' 
!
!#if defined( ESMF_ )
!      ! We need to add the ".rc" extension to the "globchem.dat" file
!      ! so that GEOS-5 will copy it to the proper run directory.
!      FILENAME   = TRIM( FILENAME ) // 'rc'
!#endif
!      KGLC       = IU_CHEMDAT
!
!      IF ( am_I_Root ) THEN 
!         WRITE( 6, 310 ) TRIM( FILENAME ), KGLC
! 310     FORMAT( 'GEOSCHEM::READER: Opening ', a, ' on unit ', i4  )
!      ENDIF
!
!      ! Open chemistry mechanism file
!      OPEN( KGLC, FILE=TRIM( FILENAME ) )
!------------------------------------------------------------------

      ! Open "smv2.log" for echoback output as unit #93
      IF ( am_I_Root ) THEN

         ! Find a free file LUN
         IU_SMV2LOG = findFreeLUN()

         WRITE( 6, 320 ) IU_SMV2LOG
 320     FORMAT( 'GEOSCHEM::READER: Opening smv2.log on unit', i4 )

         IO93 = IU_SMV2LOG
         OPEN( IO93, FILE='smv2.log', STATUS='UNKNOWN' )
!
! *********************************************************************
! *                    PRINT INFORMATION FROM m.dat                   *
! *********************************************************************
!
         IF ( IPREADER.EQ. 1 .AND. FIRSTCHEM .and. am_I_Root ) THEN
            WRITE( IO93, '(a)'   ) REPEAT( '=', 79 )
            WRITE( IO93, '(a,/)' ) 'SMV2.LOG -- SMVGEAR II information' 
            WRITE( IO93, '(a)'   ) 'Switches in mglob.dat'
            WRITE( IO93, '(a)'   ) REPEAT( '=', 79 )
            WRITE( IO93, *       ) 'IFSOLVE   = ', IFSOLVE
            WRITE( IO93, *       ) 'ITESTGEAR = ', ITESTGEAR
            WRITE( IO93, *       ) 'IFURBAN   = ', IFURBAN
            WRITE( IO93, *       ) 'IFTROP    = ', IFTROP
            WRITE( IO93, *       ) 'IFSTRAT   = ', IFSTRAT
            WRITE( IO93, *       ) 'KULOOP    = ', KULOOP
            WRITE( IO93, *       ) 'LYOUT     = ', LYOUT
            WRITE( IO93, *       ) 'LXOUT     = ', LXOUT
            WRITE( IO93, *       ) 'LZOUT     = ', LZOUT
            WRITE( IO93, *       ) 'CHEMINTV  = ', CHEMINTV
            WRITE( IO93, *       ) 'IPRATES   = ', IPRATES
            WRITE( IO93, *       ) 'IPREADER  = ', IPREADER
            WRITE( IO93, *       ) 'IOSPEC    = ', IOSPEC
            WRITE( IO93, *       ) 'IOREAC    = ', IOREAC
            WRITE( IO93, *       ) 'APGASA    = ', APGASA
            WRITE( IO93, *       ) 'APGASB    = ', APGASB
            WRITE( IO93, *       ) 'APGASC    = ', APGASC
            WRITE( IO93, *       ) 'APGASD    = ', APGASD
            WRITE( IO93, *       ) 'APGASE    = ', APGASE
            WRITE( IO93, *       ) 'APGASF    = ', APGASF
            WRITE( IO93, *       ) 'APGASG    = ', APGASG
            WRITE( IO93, *       ) 'IFREORD   = ', IFREORD 
            WRITE( IO93, *       ) 'FRACDEC   = ', FRACDEC
            WRITE( IO93, *       ) 'PLOURB    = ', PLOURB 
            WRITE( IO93, *       ) 'PLOTROP   = ', PLOTROP
            WRITE( IO93, *       ) 'ERRMAXU   = ', ERRMAXU
            WRITE( IO93, *       ) 'YLOWU     = ', YLOWU
            WRITE( IO93, *       ) 'YHIU      = ', YHIU 
            WRITE( IO93, *       ) 'HMAXDAYU  = ', HMAXDAYU 
            WRITE( IO93, *       ) 'ERRMAXR   = ', ERRMAXR
            WRITE( IO93, *       ) 'YLOWR     = ', YLOWR          
            WRITE( IO93, *       ) 'YHIR      = ', YHIR 
            WRITE( IO93, *       ) 'HMAXDAYR  = ', HMAXDAYR 
            WRITE( IO93, *       ) 'ERRMAXS   = ', ERRMAXS
            WRITE( IO93, *       ) 'YLOWS     = ', YLOWS 
            WRITE( IO93, *       ) 'YHIS      = ', YHIS
            WRITE( IO93, *       ) 'HMAXDAYS  = ', HMAXDAYS 
            WRITE( IO93, *       ) 'HMAXNIT   = ', HMAXNIT
            WRITE( IO93, '(/,a)' ) 'Using U.C.I. Fast-J photolysis'  
            WRITE(    6, '(  a)' ) 'Using U.C.I. Fast-J photolysis'     

            ! Write spacer line to "smv2.log
            WRITE( IO93, '(a)' )
         ENDIF
      ENDIF
!
! *********************************************************************
! *******             THE VALUES OF BASIC PARAMETERS            *******
! *********************************************************************
! BOLTG   = BOLTZMANN"S CONSTANT, 1.381E-16 ERG DEG K**-1 = RGAS / AVG 
!         = (1 J = 10**7 ERG = 1 N-M = 1 KG M2 S-2)
! RSTARG  = UNIVERSAL GAS CONSTANT = 8.3145E+07 G CM2 S-2 MOLE-1 K-1
! AVG     = AVOGADRO"S NUMBER,MOL**-1
! WTAIR   = MOLECULAR WEIGHT OF AIR; 
! RGAS    = GAS CONSTANT             (ERG DEG K-1 MOL-1) 
!           1 ERG = 1 DYNE-CM = 10**-7 J
!           1 ATM = 1.013 BAR = 10**5 PA. 1PA = 1 N M-2 = 10 DYNES CM-2.
! SCDAY   = SECONDS PER DAY
!
!
      NMASBAL     = 9 
      NAMEMB(  1) = 'SULFUR ATOMS'
      NAMEMB(  2) = 'NITROGEN NO3'
      NAMEMB(  3) = 'NITROGEN NH4'
      NAMEMB(  4) = 'CARBON ATOMS'
      NAMEMB(  5) = 'CHLORINE ATOMS'
      NAMEMB(  6) = 'BROMINE ATOMS'
      NAMEMB(  7) = 'FLOURINE ATOMS'
      NAMEMB(  8) = 'HYDROGEN ATOMS'
      NAMEMB(  9) = 'OXYGEN ATOMS'
!
! *********************************************************************
!
      IF (NLAT.GT.ILAT.OR.NLONG.GT.ILONG.OR.NVERT.GT.IVERT) THEN
       IF ( am_I_Root ) THEN
          WRITE(6,*)'READER: NLAT, NLONG, OR NVERT TOO BIG'
       ENDIF
       CALL GEOS_CHEM_STOP
      END IF
!
! *********************************************************************
! *                  SETUP LOOP-LOCATING ARRAYS                       * 
! *********************************************************************
!
! VALUE OF JLOOP CORRESPONDING TO EACH GRID-CELL FOR GRID
! OF      NLAT = 3, NLONG = 5, NVERT = 2.   
!
!         LAYER 1 (TOP)                  LAYER NVERT = 2 (BOTTOM)         
! M1                                 M1       
! 3 |  11  12  13  14  15             3 |  26  27  28  29  30  
! 2 |   6   7   8   9  10             2 |  21  22  23  24  25  
! 1 |   1   2   3   4   5             1 |  16  17  18  19  20    
!      -------------------                 -------------------
!       1   2   3   4   5   M2             1   2   3   4   5   M2    
!
      DO 210 M2        = 1, NLONG
       DO 210 M1       = 1, NLAT
        MLOOP          = (M1 - 1) * NLONG + M2  
 210    MLOP(M1,M2)    = MLOOP
!
      DO 220 K         = 1, NLAYER
       KLOOP           = (K - 1) * NLOOP
       DO 220 M2       = 1, NLONG
        DO 220 M1      = 1, NLAT
         MLOOP         = MLOP(M1,M2)
         JLOOP         = MLOOP + KLOOP

         ! JLOP set differently in ruralbox (bdf, 4/1/03)
         JLOP_SMV(M1,M2,K) = JLOOP
 220  CONTINUE
!
      LLOOP            = JLOP_SMV(LYOUT,LXOUT,LZOUT)
!
! *********************************************************************
!           DETERMINE HOW MANY PROCESSES SOLVED FOR IN SMVGEAR
! *********************************************************************
!
! IFURBAN  IFTROP  IFSTRAT            TYPE OF CHEMISTRY SOLVED  
!                             (U=URBAN, T=TROPOSPHERIC, S=STRATOSPHERIC)
! ----------------------------------------------------------------------
!    0        0       0      DO NOT SOLVE CHEMISTRY 
!    1        0       0      SOLVE U EVERYWHERE 
!    0        1       0      SOLVE T EVERYWHERE 
!    0        0       1      SOLVE S EVERYWHERE 
!    1        1       1      SOLVE U BELOW PLOURB, T BETWEEN PLOURB,
!                            PLOTROP, AND S ABOVE PLOTROP  
!    0        2       2      SOLVE T/S   CHEMISTRY EVERYWHERE   
!    2        2       2      SOLVE U/T/S CHEMISTRY EVERYWHERE   
! 
! IGLOBCHEM = -2 --> SOLVE ALL GAS CHEMISTRY WITH COMBINATION OF U/R/S SETS
!           = -1 --> SOLVE ALL GAS CHEMISTRY WITH COMBINATION OF R/S SETS
!           = 0  --> SOLVE ALL GAS CHEMISTRY WITH EITHER U, R, OR S SETS
!           = 1  --> SOLVE EACH REGION SEPARATELY WITH U, R, OR S SET 
!
      IF (IFURBAN.EQ.2.AND.IFTROP.EQ.2.AND.IFSTRAT.EQ.2) THEN 
       IGLOBCHEM             = -2 
       NCSALL                = 1 
       NCSTRST               = 0 
       NCSURBAN              = 0 
       NCSTROP               = 0  
       NCSSTRAT              = 0 
       NCSGAS                = 1 
       NTLOOPNCS(NCSGAS)     = NTLOOP 
      ELSEIF (IFURBAN.EQ.0.AND.IFTROP.EQ.2.AND.IFSTRAT.EQ.2) THEN 
       IGLOBCHEM             = -1 
       NCSALL                = 0 
       NCSTRST               = 1 
       NCSURBAN              = 0 
       NCSTROP               = 0  
       NCSSTRAT              = 0 
       NCSGAS                = 1 
       NTLOOPNCS(NCSGAS)     = NTLOOP 
      ELSEIF (IFURBAN.EQ.1.AND.IFTROP.EQ.1.AND.IFSTRAT.EQ.1) THEN 
       IGLOBCHEM             = 1 
       NCSALL                = 0 
       NCSTRST               = 0 
       NCSURBAN              = 1 
       NCSTROP               = 2  
       NCSSTRAT              = 3 
       NCSGAS                = 3 
      ELSE 
       IGLOBCHEM             = 0 
       NCSALL                = 0 
       NCSTRST               = 0 
       NCSURBAN              = 0 
       NCSTROP               = 0  
       NCSSTRAT              = 0 
       NCSGAS                = 1  
       IF (IFURBAN.EQ.1.AND.IFTROP.EQ.0.AND.IFSTRAT.EQ.0) THEN 
        NTLOOPNCS(NCSGAS)    = NTLOOP 
        NCSURBAN             = 1  
       ELSEIF (IFURBAN.EQ.0.AND.IFTROP.EQ.1.AND.IFSTRAT.EQ.0) THEN 
        NTLOOPNCS(NCSGAS)    = NTLOOP 
        NCSTROP              = 1  
       ELSEIF (IFURBAN.EQ.0.AND.IFTROP.EQ.0.AND.IFSTRAT.EQ.1) THEN 
        NTLOOPNCS(NCSGAS)    = NTLOOP
        NCSSTRAT             = 1  
       ELSEIF (IFURBAN.EQ.0.AND.IFTROP.EQ.0.AND.IFSTRAT.EQ.0) THEN 
        IFCHEM               = 0
        IFSOLVE              = 0
       ELSE 
        IF ( am_I_Root ) WRITE(6,265)  
        CALL GEOS_CHEM_STOP
       ENDIF 
      ENDIF 
!
 265  FORMAT('READER: NEED IFURBAN, IFSTRAT, IFTROP ALL = 1 OR JUST ',
     1       'ONE = 1')
!
! ITESTGEAR = 1: TEST SMVGEAR TO ACCURATE SOLUTION FOUND IN compare.dat
! ITESTGEAR = 2: GENERATE SMVGEAR ACCURATE SOLUTION AND WRITE TO compare.dat
!
      IF (ITESTGEAR.EQ.2) THEN
       ERRMAXU = 1.00d-09
       ERRMAXR = 1.00d-09
       ERRMAXS = 1.00d-09
!
       YLOWU   = 1.00d-10  
       YLOWR   = 1.00d-10  
       YLOWS   = 1.00d-10   
!
       YHIU    = 1.00d-10  
       YHIR    = 1.00d-10  
       YHIS    = 1.00d-10   
      ENDIF 
!
      DO 269 NCS     = 1, ICS
       ABTOL(1,NCS)  = 0.d0
       ABTOL(6,NCS)  = 0.d0
 269  CONTINUE 
!
! URBAN / REGIONAL / STRATOSPHERIC CHEMISTRY TOGETHER
!
      IF (NCSALL.GT.0) THEN 
       NCS            = NCSALL 
       NCSP           = NCS + ICS 
       CHEMTYP(  NCS) = 'URB/REG/STR' 
       ERRMAX(   NCS) = ERRMAXU  
       ABTOL(1,  NCS) = YHIU
       ABTOL(6,  NCS) = YLOWU
       TIMEINTV( NCS) = CHEMINTV 
       ABST2(    NCS) = 1. / (CHEMINTV * CHEMINTV) 
       HMAXUSE(  NCS) = HMAXDAYU    
       HMAXUSE( NCSP) = HMAXNIT
      ENDIF 
!
! REGIONAL / STRATOSPHERIC CHEMISTRY TOGETHER
!
      IF (NCSTRST.GT.0) THEN 
       NCS            = NCSTRST
       NCSP           = NCS + ICS 
       CHEMTYP(  NCS) = 'REG/STR' 
       ERRMAX(   NCS) = ERRMAXR   
       ABTOL(1,  NCS) = YHIR 
       ABTOL(6,  NCS) = YLOWR 
       TIMEINTV( NCS) = CHEMINTV 
       ABST2(    NCS) = 1. / (CHEMINTV * CHEMINTV) 
       HMAXUSE(  NCS) = HMAXDAYR   
       HMAXUSE( NCSP) = HMAXNIT
      ENDIF 
!
! URBAN CHEMISTRY 
!
      IF (NCSURBAN.GT.0) THEN 
       NCS            = NCSURBAN 
       NCSP           = NCS + ICS 
       CHEMTYP( NCS)  = 'URBAN' 
       ERRMAX(  NCS)  = ERRMAXU  
       ABTOL(1, NCS)  = YHIU  
       ABTOL(6, NCS)  = YLOWU  
       TIMEINTV(NCS)  = CHEMINTV 
       ABST2(   NCS)  = 1. / (CHEMINTV * CHEMINTV) 
       HMAXUSE( NCS)  = HMAXDAYU    
       HMAXUSE(NCSP)  = HMAXNIT
      ENDIF 
!
! TROPOSPHERIC CHEMISTRY 
!
      IF (NCSTROP.GT.0) THEN 
       NCS            = NCSTROP
       NCSP           = NCS + ICS 
       CHEMTYP( NCS)  = 'TROPOSPHERIC'  
       ERRMAX(  NCS)  = ERRMAXR
       ABTOL(1, NCS)  = YHIR   
       ABTOL(6, NCS)  = YLOWR   
       TIMEINTV(NCS)  = CHEMINTV 
       ABST2(   NCS)  = 1.d0 / (CHEMINTV * CHEMINTV) 
       HMAXUSE( NCS)  = HMAXDAYR     
       HMAXUSE(NCSP)  = HMAXNIT
      ENDIF 
!
! STRATOSPHERIC CHEMISTRY 
!
      IF (NCSSTRAT.GT.0) THEN 
       NCS            = NCSSTRAT  
       NCSP           = NCS + ICS 
       CHEMTYP( NCS)  = 'STRATOSPHERIC'  
       ERRMAX(  NCS)  = ERRMAXS 
       ABTOL(1, NCS)  = YHIS    
       ABTOL(6, NCS)  = YLOWS    
       TIMEINTV(NCS)  = CHEMINTV 
       ABST2(   NCS)  = 1.d0 / (CHEMINTV * CHEMINTV) 
       HMAXUSE( NCS)  = HMAXDAYS      
       HMAXUSE(NCSP)  = HMAXNIT
      ENDIF 
!
! CALCULATE ALL POSSIBLE REMAINING ABSOLUTE ERROR TOLERANCES
!
      DO 272 NCS      = 1, NCSGAS
       ABHI           = LOG10(ABTOL(1,NCS))
       ABLO           = LOG10(ABTOL(6,NCS))
!
       IF (ABHI.LT.ABLO) THEN
        IF ( am_I_Root ) THEN
           WRITE(6,*)'READER: ABHI < ABLO - INCREASE UPPER BOUND OF', 
     &               'ABSOLUTE ERROR TOLERANCE FOR NCS = ',NCS, 
     &                ABTOL(1,NCS),ABTOL(6,NCS) 
        ENDIF
        CALL GEOS_CHEM_STOP
       ENDIF
!
       DO 270 I       = 2, 5 
        ABTOL(I,NCS)  = 10.d0**(ABLO + (ABHI - ABLO) *FLOAT(6-I) / 5.d0)
 270   CONTINUE 
 272  CONTINUE
!
! *********************************************************************
!
      !FP_ISOP
      ! update NMPROD FP (hotp 8/1/09)
      NMREAC         = 3
      NALLREAC       = 4 
      !NMPROD         = 14
      NMPROD         = NREAD-NALLREAC
      NPRODLO        = NALLREAC + 1 
      NPRODHI        = NALLREAC + NMPROD
      IFDID          = 0
      IFNEVER        = 0
      IFNONE         = 0 
      NSFTOT         = 0
      NPDTOT         = 0
      NSTTOT         = 0
      IFAILTOT       = 0
      LFAILTOT       = 0
      NFAILTOT       = 0
      NOCC           = 0
      SUMAVGE        = 0.d0
      SUMAVHI        = 0.d0
      SUMRMSE        = 0.d0
      SUMRMHI        = 0.d0
      TOTSTEP        = 0.d0
      TOTIT          = 0.d0
      TELAPS         = 0.d0
      RMSERR         = 0.d0
!
      MB1            = 1 
      MB2            = 2 
      DO 660 I       = 1, IMASBAL
       MBCOMP(I,MB1) = 0.d0
       MBCOMP(I,MB2) = 0.d0
 660  CONTINUE
!
      DO 705 NCS       = 1, ICS
       NAMENCS(0,NCS) = ' '
       NMOTH(    NCS) = 0
       NTSPEC(   NCS) = 0
       JPHOTRAT( NCS) = 0 
       ISGAINR(  NCS) = 0 
       ISPORL(   NCS) = 0 
       NOGAINE(  NCS) = 0 
       NOUSE(    NCS) = 0
       NSPEC(    NCS) = 0
       NTRATES(  NCS) = 0  
       ISGAINE(  NCS) = 0 
       NSPCSOLV( NCS) = 0 
       ISCHANG(  NCS) = 0 
       NRATES(   NCS) = 0
       NM3BOD(   NCS) = 0
       ITWOR(    NCS) = 0
       ITHRR(    NCS) = 0 
       INOREP(   NCS) = 0
       NRATCUR(  NCS) = 0
       NSURFACE( NCS) = 0 
       NPRESM(   NCS) = 0 
       NMAIR(    NCS) = 0 
       NMO2(     NCS) = 0 
       NMN2(     NCS) = 0 
       NNEQ(     NCS) = 0 
       NARR(     NCS) = 0 
       NABR(     NCS) = 0 
       NACR(     NCS) = 0 
       NABC(     NCS) = 0 
       NKSPECW(  NCS) = 0
       NKSPECX(  NCS) = 0
       NKSPECY(  NCS) = 0
       NKSPECZ(  NCS) = 0
 705  CONTINUE

      ! Zero out entire nkspecv array (bdf, 4/1/03)
      NKSPECV = 0d0
!
      DO 710 NCS       = 1, ICP 
       NOLOSP(   NCS) = 0  
       NGNFRAC(  NCS) = 0
       NOLOSRAT( NCS) = 0
       IARRAY(   NCS) = 0
       NALLRAT(  NCS) = 0
       KZTLO(    NCS) = 0
       KZTHI(    NCS) = 0
       IONER(    NCS) = 0
       NPLLO(    NCS) = 0
       NPLHI(    NCS) = 0
       NFRLO(    NCS) = 0
       NFRHI(    NCS) = 0
       NPDLO(    NCS) = 0
       NPDHI(    NCS) = 0
 710  CONTINUE  
!
      DO 715 NCS        = 1, ICS
       DO 714 I         = 1, MAXGL
        FRACP(   I,NCS) = 0  
        IGNFRAC( I,NCS) = 0
        NKGNFRAC(I,NCS) = 0  
 714   CONTINUE
 715  CONTINUE
!
      DO 720 NCS        = 1, ICS
       DO 719 I         = 1, MAXGL2  
        NREACOTH(I,NCS) = 0
        LGASBINO(I,NCS) = 0
 719   CONTINUE
 720  CONTINUE
!
      DO 725 NCS        = 1, ICS
       DO 724 I         = 1, MAXGL3 
        NKNLOSP( I,NCS) = 0
        LOSINACP(I,NCS) = 0
        NREACAIR(I,NCS) = 0
        NREAC3B( I,NCS) = 0
        NREACEQ( I,NCS) = 0
        NREQOTH( I,NCS) = 0
        NREACN2( I,NCS) = 0
        NREACO2( I,NCS) = 0
        NREACPM( I,NCS) = 0
        LGAS3BOD(I,NCS) = 0
 724   CONTINUE
 725  CONTINUE
!
      DO 735 NCS        = 1, ICS
       DO 734 I         = 1, MXGSAER 
        NAMENCS( I,NCS) = ' '
        FRACGAIN(I,NCS) = 0.d0
        QBKCHEM( I,NCS) = 0.d0
        NUMLOST( I,NCS) = 0        
        NUMGFRT( I,NCS) = 0 
        NUMGAINT(I,NCS) = 0
        NGAINE(  I,NCS) = 0
        IGAINR(  I,NCS) = 0
        IPORL(   I,NCS) = 0
        IGAINE(  I,NCS) = 0
        ISOLVSPC(I,NCS) = 0
        INEWOLD( I,NCS) = 0
        MAPPL(   I,NCS) = 0
 734   CONTINUE 
 735  CONTINUE 
!
      DO 740 NCS        = 1, ICP 
       DO 739 I         = 1, MXGSAER 
        NUMLOSS( I,NCS) = 0        
        NUMGAIN( I,NCS) = 0
        NUMPORL( I,NCS) = 0        
 739   CONTINUE 
 740  CONTINUE 
!
      DO 745 NCS        = 1, ICS
       DO 744 I         = 1, NMTRATE
        I1              = NMTRATE + I
        ARR(     I,NCS) = 0.d0  
        BRR(     I,NCS) = 0.d0
        FCV(     I,NCS) = 0.d0
        FCTEMP1( I,NCS) = 0.d0
        FCTEMP2( I,NCS) = 0.d0
        NKARR(   I,NCS) = 0
        NKABR(   I,NCS) = 0
        NKACR(   I,NCS) = 0
        NKABC(   I,NCS) = 0
        IRORD(   I,NCS) = 0 
        IAPROD(  I,NCS) = 0
        NOLOSRN( I,NCS) = 0
        NRUSE(   I,NCS) = 0
        NRREP(   I,NCS) = 0
        NPRODUC( I,NCS) = 0
        NCEQUAT( I,NCS) = 0 
        NOLDFNEW(I,NCS) = 0
        NEWFOLD( I,NCS) = 0
        NEWFOLD(I1,NCS) = 0
        NKONER(  I,NCS) = 0
        NKTWOR(  I,NCS) = 0
        NKTHRR(  I,NCS) = 0
        KCRR(    I,NCS) = 0 
        JPHOTNK( I,NCS) = 0
 744   CONTINUE
 745  CONTINUE
!
      DO 755 NCS         = 1, ICS
       DO 754 J          = 1, IPHOT
        NKPHOTRAT(J,NCS) = 0 
        NPPHOTRAT(J,NCS) = 0 
        NKNPHOTRT(J,NCS) = 0 
 754   CONTINUE
 755  CONTINUE
!
      DO 765 NCS         = 1, ICP 
       DO 764 I          = 1, MXGSAER
        JARRDIAG(I,NCS)  = 0
        JLOZ1(   I,NCS)  = 0
        JHIZ1(   I,NCS)  = 0
        IJTLO(   I,NCS)  = 0
        IJTHI(   I,NCS)  = 0
        IMZTOT(  I,NCS)  = 0
 764   CONTINUE
 765  CONTINUE

      DO 770 NCS         = 1, ICS
       DO 769 NK         = 1, NMTRATE
        DO 768 I         = 1, NMRPROD
         IRM(  I,NK,NCS) = 0
         IRM2( I,NK,NCS) = 0
         FKOEF(I,NK,NCS) = 0.d0
         FK2(  I,NK,NCS) = 0.d0
 768    CONTINUE
 769   CONTINUE
 770  CONTINUE
!
      DO 775 NCS         = 1, ICS
       DO 774 J          = 1, MAXGL
        DO 773 I         = 1, MXGSAER
         JPORL(I,J,NCS)  = 0
 773    CONTINUE
 774   CONTINUE
 775  CONTINUE

      ! Set NCS=NCSURBAN here since we have defined our tropospheric
      ! chemistry mechanism in the urban slot of SMVGEAR II
      NCS = NCSURBAN
!     
! *********************************************************************
! ********************** END OF SUBROUTINE READER *********************
! *********************************************************************
!
 100  FORMAT(A72)
 110  FORMAT(32X,'SMVGEAR II')
 115  FORMAT(//,35X,'*****                  MAIN SWITCHES',
     1              '                  *****',/)
!
      RETURN
      END SUBROUTINE READER
