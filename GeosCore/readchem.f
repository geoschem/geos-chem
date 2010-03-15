! $Id: readchem.f,v 1.2 2010/03/15 19:33:21 ccarouge Exp $
      SUBROUTINE READCHEM 
!
!******************************************************************************
!  Subroutine READCHEM reads species 2names, chemical rxns, and photolysis 
!  reactions from the "globchem.dat" chemistry mechanism file for SMVGEAR II.  
!  (M. Jacobson 1997; bdf, bmy, 5/9/03, 8/9/06)
!
!  NOTES:
!  (1 ) Added space in FORMAT strings for more products.  Also now references
!        MAXDEP from "drydep_mod.f".  Now also writes species and reactions
!        to the "smv2.log" output file as unit #93.  Now call GEOS_CHEM_STOP
!        to deallocate all arrays and stop the run safely.  Add NNADDG and 
!        NKSPECG for DMS+OH+O2 rxn.  Now also force double precision with
!        the "D" exponent.  Now call SETPL before JSPARSE so that the prod/loss
!        families will be computed correctly. (bmy, 5/9/03)
!  (2 ) Now initialize ICH4 -- the location of CH4 in the CSPEC array.  Now 
!        define lookup table ITS_NOT_A_ND65_FAMILY, which is used to exclude
!        ND65 prod/loss families from modifying the SMVGEAR II convergence
!        criteria.  (bnd, bmy, 7/9/03)
!  (3 ) Now declare ININT as a local variable instead of being declared w/in 
!        "comode.h".  Remove reference to IPORD. (bmy, 7/16/03)
!  (4 ) Now flag the N2O5 hydrolysis rxn for later use. (mje, bmy, 8/7/03)
!  (5 ) Now references SETJFAM & SETPL from "diag_pl_mod.f" (bmy, 7/20/04)
!  (6 ) Now look up ILISOPOH, the index of ISOP lost to OH (dkh, bmy, 6/1/06)
!  (7 ) Increase FORMAT 510 so that it has space for 14 products (bmy, 8/9/06)
!  (8 ) Now flag the HO2 heterogeneous uptake rxn for later use  
!        (jaegle, 02/26/09)
!  (9 ) Added identifier to mark branching rxns for HOC2H4O (tmf, 1/7/09)
!          HOC2H4O ------> HO2 + 2CH2O         : marked as 'F' in P column
!          HOC2H4O --O2--> HO2 + GLYC          : marked as 'H' in P column
!
!       The same branching rxns are also implemented in EP photolysis
!          HOC2H4O ------> HO2 + 2CH2O         : marked as 'I' in P column
!          HOC2H4O --O2--> HO2 + GLYC          : marked as 'J' in P column
!******************************************************************************
!
      ! References to F90 modules
      USE DRYDEP_MOD,  ONLY : MAXDEP
      USE ERROR_MOD,   ONLY : GEOS_CHEM_STOP
      USE DIAG_PL_MOD, ONLY : SETJFAM, SETPL

      !FP_ISOP (6/2009)
      USE ERROR_MOD,   ONLY : DEBUG_MSG
      USE LOGICAL_MOD, ONLY : LPRT

      IMPLICIT NONE

#     include "CMN_SIZE"  ! Size parameters
#     include "comode.h"  ! SMVGEAR II arrays
C
C *********************************************************************
C ************        WRITTEN BY MARK JACOBSON (1990-4)    ************
C ***                         (650) 723-6836                        *** 
C *********************************************************************
C
C  RRRRRR   EEEEEEE     A     DDDDDD  CCCCCCC  H     H  EEEEEEE M     M  
C  R     R  E          A A    D     D C        H     H  E       M M M M 
C  RRRRRR   EEEE      A   A   D     D C        HHHHHHH  EEEEEEE M  M  M  
C  R  R     E        AAAAAAA  D     D C        H     H  E       M     M
C  R   R    EEEEEEE A       A DDDDDD  CCCCCCC  H     H  EEEEEEE M     M 
C   
C *********************************************************************
C *  THIS IS THE SETUP ROUTINE FOR GAS-PHASE CHEMISTRY. IT READS      *
C *  SPECIES NAMES, CHEMICAL REACTIONS, AND PHOTOPROCESSES FROM AN    *
C *  INPUT DATA SET. IT THEN PLACES ALL NECESSARY INFORMATION INTO    *
C *  ARRAYS AND PRINTS OUT THE INPUT INFORMATION.                     *
C *                                                                   *
C * HOW TO CALL SUBROUTINE:                                           *
C * ----------------------                                            *
C *  CALL READCHEM.F FROM MAIN.F                                      *
C *                                                                   *
C *********************************************************************
C
C *********************************************************************
C *                  SOME PARAMETER DEFINITIONS                       * 
C *********************************************************************
C
C DEFPRAT   = DEFAULT PHOTORATE (SEC-1)
C IGAS      = DIMENSION OF MAXIMUM NUMBER OF GAS SPECIES, ACTIVE + INACTIVE.
C IPHOT     = MAXIMUM NUMBER OF RADIATIVELY ACTIVE SPECIES 
C IPORD     = ORDINAL # OF PHOTOPROCESS (USED TO IDENTIFY REACTION)
C IRORD     = ORDINAL # OF KINET. REACT. (USED TO IDENTIFY REACTION) 
C JMBCOMP   = SPECIES NUMBER FOR EACH SPECIES IN A MASS BAL. GROUP
C MBCOMP    = COUNTS THE MASS BALANCE SPECIES IN EACH M.B. GROUP (I.E.
C             SULFUR IS A M.B. GROUP.
C NACTIVE   = NUMBER OF ACTIVE SPECIES READ IN -- (A) IN COLUMN ONE
C             OF INPUT DATA SET, CONVERTED TO NSPEC LATER 
C NALLREAC  = TOTAL NUMBER OF REACTANT POSITIONS IN A REACTION (BUT
C             NUMBER OF ACTIVE POSITIONS IN NMREAC) 
C NAMD      = NAMES OF SPECIES WHICH MAY APPEAR IN REACTIONS BUT WHICH ARE
C            "DEAD" WITH RESPECT TO THE PHOTOCHEMISTRY AND THUS ARE NOT
C             PRINTED OUT.
C NAMEGAS   = CHARACTER ARRAY OF SPECIES NAMES. 
C NAMENCS   = CHARACTER ARRAY OF SPECIES NAMES. 
C NCS       = 1..NCSGAS FOR GAS CHEMISTRY                            
C NGAS      = NSPEC, THE NUMBER OF ACTIVE SPECIES
C NINAC     = NUMBER OF INACTIVE SPECIES READ -- (I) IN COLUMN ONE 
C NMAIR     = # REACTIONS WHERE THE SPECIES IN THE THIRD POSITION
C             IS 'M' = 'O2 + N2' 
C NMASBAL   = NUMBER OF MASS BALANCE GROUPS (E.G. S, N, C ARE GROUPS)
C NMN2      = # REACTIONS WHERE SPECIES IN THE THIRD POSITION IS N2
C NMO2      = # REACTIONS WHERE SPECIES IN THE THIRD POSITION IS O2
C NM3BOD    = # REACTIONS WHERE SPECIES IN THE THIRD POSITION
C             IS ANY OTHER SPECIES:  IRM(3,NK) = -SPECIES NUMBER
C NMPROD    = MAXIMUM NUMBER OF ACTIVE PRODUCTS IN A REACTION (READER.F)
C NMREAC    = MAXIMUM NUMBER OF ACTIVE REACTANTS IN A REACTION
C             OF INPUT DATA SET. (DEFINED IN READER.F)   
C NNEQ      = # THERMALLY DISSOCIATING EQUILIBRIUM REACTIONS. PREVIOUS
C             EQUATION MUST BE PRESSURE-DEPENDENT.
C NPHOTALL  = NUMBER OF ACTIVE GAS PHOTOPROCESSES
C NPRESM    = # PRESSURE DEPENDENT 3-BODY REACTIONS  
C NPRODHI   = HIGHEST PRODUCT POSITON IN A REACTION
C NPRODLO   = LOWEST PRODUCT POSITON IN A REACTION
C NRATES    = NUMBER OF ACTIVE REACTIONS (EXCLUDING PHOTPROCESSES)
C NSDEAD    = NUMBER OF DEAD SPECIES READ IN -- (D) IN COLUMN ONE OF
C             INPUT DATA SET.
C NSPEC     = TOTAL NUMBER OF ACTIVE SPECIES.
C NTRATES   = NUMBER OF ACTIVE KINETIC REACTIONS PLUS PHOTOPROCESSES
C NTSPEC    = ACTUAL NUMBER OF ACTIVE + INACTIVE (BUT NOT DEAD) SPECIES. 
C QBKGAS    = DEFAULT BACKGROUND CONCENTRATION (VOL MIXING RATIO) 
C RINP      = CHARACTER ARRAY FOR READING IN INFORMATION FROM DATA SETS.
C WTMB      = MASS BALANCE WEIGHT FOR EACH M. B. SPECIES
C XINP      = CHARACTER ARRAY FOR READING IN INFORMATION FROM DATA SETS.
C
C *********************************************************************
C * SET INITIAL VALUES AND READ INITIAL COMMENTS FROM INPUT DATA SET  *
C *********************************************************************
C

      INTEGER NINAC,NACTIVE,NSDEAD,NOTHGS,IDOPHOT,I,NMBGAS,NM,JGAS
      INTEGER MB,MBP,INACT1,JGAS1,J,IORD,NCOF,JORD,NDUM,NK,NAR,NK1
      INTEGER JPR,JNUM,ITHIRDB,NM2,NR,NN,JGAS2,JGAS3,NA,N,NS2

      REAL*8 C1,CSTRAT,CTROPL,CTROPS,CURBAN,QTHERMG

      ! ININT used to be defined w/in "comode.h", but it is only ever used
      ! w/in "readchem.f".  Declare here as a local variable. (bmy, 7/16/03)
      INTEGER :: ININT(10)  

      !=================================================================
      ! READCHEM begins here!
      !=================================================================
      NINAC          = 0
      NACTIVE        = 0 
      NSDEAD         = 0 
      NOTHGS         = 0 
      NPHOTALL       = 0 
      IDOPHOT        = 0
C
      NAMEGAS(0)     = ' '
C
      ! Initialize flag for N2O5 hydrolysis rxn (bmy, 8/7/03)
      NKN2O5         = 0

      ! Initialize flag for HO2 hydrolysis rxn (jaegle, 02/26/09)
      NKHO2          = 0

      DO 44 I            = 1, NMDEAD
       NAMD(I)           = ' '
 44   CONTINUE
C
      DO 46 I            = 1, IGAS
       NAMEGAS(I)        = ' '
       WTGAS(  I)        = 0.d0
       QBKGAS( I)        = 0.d0
 46   CONTINUE
C
      DO 47 I            = 1, MXGSAER
       CPREV( I)         = 0.d0
       CMODEL(I)         = 0.d0
 47   CONTINUE
C
      DO 48 I            = 1, MAXGL4  
       NKSURF(I)         = 0
       NCOATG(I)         = 0
 48   CONTINUE
C
      DO 49 I            = 1, IPHOT
       DEFPRAT(I,:)        = 0.d0
 49   CONTINUE
C
      READ(KGLC,21)  HEADING 
      !WRITE(6,21) HEADING 
 13   READ(KGLC,21)  HEADING 
      IF (HEADING.NE.'BEGIN') GOTO 13 
 21   FORMAT(A76)
C
C *********************************************************************
C                   READ IN MASS BALANCE GROUPS USED 
C *********************************************************************
C
      NMBGAS = 9 
      READ(KGLC,59) (RINP(I), I = 1, NMBGAS)  
      READ(KGLC,61) (SINP(I), I = 1, NMBGAS)   
C
      DO 36 I         = 1, NMBGAS
       MBTRACE(I)     = 0 
       IF (SINP(I).EQ.'A') THEN 
        DO 34 NM       = 1, NMASBAL 
         IF((NAMEMB(NM).EQ.'SULFUR ATOMS'  .AND.RINP(I).EQ.'SUL').OR. 
     1      (NAMEMB(NM).EQ.'NITROGEN NO3'  .AND.RINP(I).EQ.'NO3').OR.
     2      (NAMEMB(NM).EQ.'NITROGEN NH4'  .AND.RINP(I).EQ.'NH4').OR.
     3      (NAMEMB(NM).EQ.'CARBON ATOMS'  .AND.RINP(I).EQ.'CAR').OR.
     4      (NAMEMB(NM).EQ.'CHLORINE ATOMS'.AND.RINP(I).EQ.'CHL').OR.
     5      (NAMEMB(NM).EQ.'BROMINE ATOMS' .AND.RINP(I).EQ.'BRO').OR.
     6      (NAMEMB(NM).EQ.'FLOURINE ATOMS'.AND.RINP(I).EQ.'FLO').OR.
     7      (NAMEMB(NM).EQ.'HYDROGEN ATOMS'.AND.RINP(I).EQ.'HYD').OR.
     8      (NAMEMB(NM).EQ.'OXYGEN ATOMS'  .AND.RINP(I).EQ.'OXY'))THEN
          MBTRACE(I) = NM 
          GOTO 36 
         ENDIF 
 34     CONTINUE 
        WRITE(6,33) RINP(I) 
        CALL GEOS_CHEM_STOP
       ENDIF 
C
 36   CONTINUE
C
 59   FORMAT(20X,A3,8(1X,A3)) 
 61   FORMAT(21X,A1,8(3X,A1)) 
 33   FORMAT('READCHEM: MASS BALANCE GROUP ',A14,' NOT SET') 
C
C *********************************************************************
C * READ IN THE SPECIES AND OTHER DATA FOR THIS RUN FROM INPUT DATA   *
C *********************************************************************
C
C                    ITEMS IN THE FIRST READ STATEMENT
C                    ---------------------------------
C
C A/I/D 
C         D = SPECIES IS DEAD AND NOT USED
C         I = INACTIVE BUT USED (THESE SPECIES MUST ALSO BE INITIALIZED) 
C         A = SPECIES USED WHEN IFURBAN, IFTROP, OR IFSTRAT > 0
C             (URBAN, TROPSOSPHERIC AND STRATOSPHERIC SETS) 
C         U = SPECIES USED WHEN IFURBAN > 0
C         S = SPECIES USED WHEN IFSTRAT > 0
C         T = SPECIES USED WHEN IFTROP  > 0
C         R = SPECIES USED WHEN IFURBAN OR IFTROP  > 0
C         H = SPECIES USED WHEN IFTROP  OR IFSTRAT > 0 
C SPEC      = NAME OF THE SPECIES, 
C AB          TELLS WHETHER SPECIES ABSORBS RADIATION (THE SPECIES
C             DOES NOT NECESSARILY PHOTOLYZE)
C MW        = ATOMIC MASS IN AMU;  
C IFSTRAT   = 1: SOLVE STRATOSPHERIC CHEMISTRY
C IFTROP    = 1: SOLVE FREE TROPOSPHERIC CHEMISTRY
C IFURBAN   = 1: SOLVE URBAN CHEMISTRY
C INITCONC  = DEF'T BACKGROUND CONC. AT LOWEST LEVEL (VOL MIXING RATIO); 
C CSTRAT    = DEFAULT VOL MIX RATIO (FRACTION) IN STRATOSPHERE (25 KM) 
C CTROPL    = DEFAULT VOL MIX RATIO (FRACTION) IN FREE TROP OVER LAND (0 KM) 
C CTROPS    = DEFAULT VOL MIX RATIO (FRACTION) IN FREE TROP OVER SEA  (0 KM) 
C CURBAN    = DEFAULT VOL MIX RATIO (FRACTION) IN URBAN REGIONS (0 KM) 
C
C *********************************************************************
C ******************** READ IN SPECIES INFORMATION ********************
C *********************************************************************
C
C            FORMAT OF ITEMS IN THE SPECIES-LIST READ STATEMENT
C
C A/I/D SPEC   AB   MW  CSTRAT CTROPL CTROPS  CURBAN 
C A1,1X, A14,A2,1X,F6.0,E9.2,  E9.2,  E9.2    E9.2   
C
C        S     NO3  NH4     C    CL    BR    F     H     O    
C 20X,   I2,2X,I2,2X,I2,2X,I2,2X,I2,2X,I2,2X,I2,2X,I2,2X,I2 
C

      ! Read 1st line of species list
 10   READ(KGLC,11) (XINP(I), I=1,3),C1,CSTRAT,CTROPL,CTROPS,CURBAN
 11   FORMAT(A1,1X,A14,A2,1X,0PF6.2,4(1PE10.3))

      ! Test for "END" here (bmy, 4/7/03)
      IF (XINP(2).EQ.'END') GOTO 15

      ! Read 2nd line of species list 
      READ(KGLC,14) (ININT(I),I = 1, NMBGAS)
 14   FORMAT(20X,I2,8(2X,I2))
C
C *********************************************************************
C *  COUNT ACTIVE, INACTIVE, AND DEAD SPECIES. ALSO, SET UP ARRAYS    * 
C *  FOR OTHER INFORMATION.                                           *
C *********************************************************************

      IF (XINP(2).EQ.'END') GOTO 15
C
      IF (XINP(1).EQ.'U'.OR.XINP(1).EQ.'T'.OR.XINP(1).EQ.'S'.OR.
     1    XINP(1).EQ.'R'.OR.XINP(1).EQ.'H') THEN
       IF ((XINP(1).EQ.'U'.AND.IFURBAN.EQ.0).OR.
     1     (XINP(1).EQ.'T'.AND.IFTROP .EQ.0).OR.
     2     (XINP(1).EQ.'S'.AND.IFSTRAT.EQ.0).OR.
     3     (XINP(1).EQ.'R'.AND.IFURBAN.EQ.0.AND.IFTROP.EQ.0).OR.
     4     (XINP(1).EQ.'H'.AND.IFSTRAT.EQ.0.AND.IFTROP.EQ.0)) THEN
        XINP(1)              = 'D'
       ELSE 
        XINP(1)              = 'A'
       ENDIF 
      ENDIF 
C
      IF (XINP(1).EQ.'D') THEN 
       NSDEAD                = NSDEAD + 1
       NAMD(NSDEAD)          = XINP(2) 
       GOTO 10
      ELSEIF (XINP(1).EQ.'I') THEN 
       NINAC                 = NINAC   + 1
       JGAS                  = IGAS    - NINAC + 1
      ELSEIF (XINP(1).EQ.'A') THEN 
       NACTIVE               = NACTIVE + 1
       JGAS                  = NACTIVE 
C
       DO 41 I                 = 1, NMBGAS 
        MB                     = MBTRACE(I) 
        IF (MB.GT.0.AND.ININT(I).GT.0) THEN 
         MBCOMP(MB,MB1)       = MBCOMP(MB,MB1) + 1 
         MBP                  = MBCOMP(MB,MB1)
         JMBCOMP(MB,MBP,MB1)  = NACTIVE  
         WTMB(MB,NACTIVE,MB1) = ININT(I)
        ENDIF
 41    CONTINUE
      ELSE
       WRITE(6,19) XINP(2), XINP(1)
       CALL GEOS_CHEM_STOP
      ENDIF
C
      NAMEGAS(JGAS)       = XINP(2) 
      WTGAS(  JGAS)       = C1
C
      IF (IFSTRAT.EQ.1.AND.IFTROP.EQ.0.AND.IFURBAN.EQ.0) THEN
       QBKGAS( JGAS)      = CSTRAT
      ELSEIF (IFSTRAT.EQ.0.AND.IFTROP.EQ.1.AND.IFURBAN.EQ.0) THEN
       QBKGAS( JGAS)      = CTROPL 
      ELSEIF (IFSTRAT.EQ.0.AND.IFTROP.EQ.0.AND.IFURBAN.EQ.1) THEN
       QBKGAS( JGAS)      = CURBAN
      ELSE
       QBKGAS( JGAS)      = CTROPL 
      ENDIF
C
      GOTO 10
C
C *********************************************************************
C *  SET NSPEC AS NUMBER OF ACTIVE SPECIES - 1 SINCE JUST INCREASED   *
C *  NACTIVE BEFORE THE 'END' STATEMENT. ALSO, CHECK SOME DIMENSIONS. *
C *********************************************************************
C
 15   CONTINUE
      NGAS         = NACTIVE
      NTSPECGAS    = NGAS + NINAC

      !=================================================================
      ! Chemical prod-loss diagnostic (bdf, bmy, 4/18/03)
      !=================================================================
      IF ( LFAMILY ) THEN

         ! Find species and rxns for ND65 diagnostic families
         CALL SETJFAM( NACTIVE, NINAC )

         ! Reset quantities after SETJFAM 
         NSPEC(NCS)  = NACTIVE - 1
         NGAS        = NSPEC(NCS)
         NTSPECGAS   = NGAS + NINAC
         NTSPEC(NCS) = NGAS + NINAC
      ENDIF
C
      IF (NTSPECGAS.GT.IGAS   .OR.  NSDEAD.GT. NMDEAD) THEN
        WRITE(6,18) IGAS,   NTSPECGAS,  NMDEAD,   NSDEAD   
        CALL GEOS_CHEM_STOP
      ENDIF
C
 18   FORMAT('READCHEM: ONE OF THE DIMENSIONS BELOW IS TOO SMALL:',/,
     1       'DIMENSION: IGAS     = ',I3,' VARIABLE: NTSPECGS = ',I3,/,  
     2       'DIMENSION: NMDEAD   = ',I3,' VARIABLE: NSDEAD   = ',I3)
C
C *********************************************************************
C *  RE-ARRANGE INACTIVE GAS ARRAYS SO THAT INFORMATION OF INACTIVE   *
C *  SPECIES APPEARS IMMEDIATELY AFTER INFORMATION OF ACTIVE SPECIES  * 
C *********************************************************************
C

      IF (NINAC.GT.0) THEN 
       INACT1             = IGAS - NINAC
       DO 26 N            = 1, NINAC
        JGAS              = NGAS   + N
        JGAS1             = INACT1 + N
        NAMEGAS(JGAS)     = NAMEGAS(JGAS1)
        QBKGAS( JGAS)     = QBKGAS( JGAS1)
        WTGAS(  JGAS)     = WTGAS(  JGAS1)
 26    CONTINUE
      END IF
C
C *********************************************************************
C * PRINT SPECIES INFORMATION IF IOSPEC = 1 + PRINT MASS BALANCE INFO * 
C *********************************************************************
C
      IF (IOSPEC.EQ.1) THEN 

       ! Write species header
       WRITE( IO93, '(/,a)' ) REPEAT( '=', 79 )
       WRITE( IO93, 23      )
       WRITE( IO93, '(a,/)' ) REPEAT( '=', 79 )
       WRITE( IO93, 22      )

       ! Write species to "smv2.log"
       DO 25 N = 1, NGAS
          WRITE(IO93,24) N, NAMEGAS(N), WTGAS(  N), QBKGAS(N)
 25    CONTINUE
       IF (NINAC.GT.0)  WRITE(IO93,28)(NAMEGAS(NGAS+N),N=1,NINAC)
       IF (NSDEAD.GT.0) WRITE(IO93,31)(NAMD(N),N=1,NSDEAD)
      END IF
C
      WRITE(6,*)
      DO 77 I = 1, NMASBAL
       IF (MBCOMP(I,MB1).GT.0) THEN 
        WRITE(6,98) NAMEMB(I) 
        WRITE(6,99)(WTMB( I,JMBCOMP(I,J,MB1),MB1),
     1              NAMEGAS(JMBCOMP(I,J,MB1)), J = 1, MBCOMP(I,MB1))
       ENDIF 
 77   CONTINUE 
C
 19   FORMAT('SPECIES ACTIVITY IS UNDEFINED ',A14,' TYPE = ',A2   )
 23   FORMAT('SPECIES FOR THIS RUN.  PHYSICAL CONSTS AND BOUNDARY', 
     1 ' CONDITIONS ALSO GIVEN.')
 22   FORMAT( 'NBR', 1X, 'NAME', 12X, 'MW', 1X, 'BKGAS(VMRAT)' )
 24   FORMAT(I3,1X,   A14,F7.2,1PE9.2)
 28   FORMAT(/'INACTIVE SPECIES FOR THIS RUN ARE:'//4(1X,A14))
 31   FORMAT(/'THE DEAD SPECIES FOR THIS RUN ARE:'//4(1X,A14)) 
 98   FORMAT('WEIGHTS AND SPECIES FOR MASS BALANCE EQUATION # ',A14)
 99   FORMAT(4(0PF5.1,1X,A14))  
C
C *********************************************************************
C *  SEARCH FOR SPECIFIC SPECIES NUMBERS USED IN OTHER SUBROUTINES    * 
C *********************************************************************
C
      ! Initialize for safety's sake (bmy, 7/7/03)
      IOXYGEN  = 0
      IH2O     = 0
      ICH4     = 0
      ILISOPOH =0
      ! add these as well (dkh, 10/06/06)  
      ILBRO2H = 0
      ILBRO2N = 0
      ILTRO2H = 0
      ILTRO2N = 0
      ILXRO2H = 0
      ILXRO2N = 0

      ! Locate positions of O2, H2O, CH4, LISOPOH in CSPEC array
      DO I = 1, NTSPECGAS
         SELECT CASE ( TRIM( NAMEGAS(I) ) )
            CASE( 'O2'      )
               IOXYGEN  = I
            CASE( 'H2O'     )
               IH2O     = I
            CASE( 'CH4'     )
               ICH4     = I
            CASE( 'LISOPOH' )
               ILISOPOH = I
            ! Add definitions (dkh, 10/06/06)  
            CASE( 'LBRO2H' )
               ILBRO2H = I
            CASE( 'LBRO2N' )
               ILBRO2N = I
            CASE( 'LTRO2H' )
               ILTRO2H = I
            CASE( 'LTRO2N' )
               ILTRO2N = I
            CASE( 'LXRO2H' )
               ILXRO2H = I
            CASE( 'LXRO2N' )
               ILXRO2N = I
            CASE DEFAULT
               ! Nothing
         END SELECT
      ENDDO
C
C *********************************************************************
C *****************    READ IN REACTION RATES   *********************** 
C *********************************************************************
C
C HERE, WE CAN HAVE 3 REACTANTS (EACH WITH COEFFICIENT OF 1) AND 12 
C             PRODUCTS [EACH WITH ANY REAL COEFFICIENT]. 
C FOR A 3-BODY M REACTION, PLACE AN M IN THE THIRD REACTANT POSITION 
C             [NO [+] BEFORE IT] 
C FOR A 3-BODY OTHER SPECIES REACTION, PLACE THE SPECIES NAME IN THE THIRD
C             REACTANT POSITION [NO [+] BEFORE IT] 
C FOR A 3RD REACTANT, PLACE THE SPECIES NAME IN THE THIRD REACTANT POSITION
C             WITH A PLUS BEFORE IT. 
C FOR A REACTANT NOT INCLUDED IN THE REACTION RATE [I.E.02] PLACE THE
C             SPECIES NAME IN THE FOURTH REACTANT POSITION [NO [+]  
C             BEFORE IT]. THE SPECIES MAY HAVE A COEFFICIENT PRECEDING IT.
C A PRODUCT MAY EITHER BE LISTED TWICE [OR MORE TIMES] OR 
C             HAVE A COEFFICIENT [I.E. 2 OR 3, 0.34] IMMEDIATELY BEFORE IT.
C N COLUMN:   NUMBER OF RATE COEFFICIENT LINES FOLLOWING TOP LINE
C P COLUMN:
C   P         = REACTION IS PRESSURE DEPENDENT 3-BODY REACTION. 
C               THE FIRST COEFFICIENT IS A 3-BODY COEF. THE SECOND IS 2-BODY.
C   S         = IDENTIFIES A ONE-BODY SURFACE REACTION 
C   E         = IDENTIFIES REVERSE EQUILIBRIUM REACTION
C   V         = IDENTIFIES CH3SCH3 + OH --> CH3S(OH)CH3
C   W         = IDENTIFIES O(1D) + N2 OR O2 
C   X         = IDENTIFIES OH  + HNO3 
C   Y         = IDENTIFIES OH  + CO   
C   Z         = IDENTIFIES HO2 + HO2  
C   G         = IDENTIFIES DMS + OH + O2
C   K         = IDENTIFIES WETDEP or HYDROLYSIS REACTIONS
C
C Fc COLUMN   = VALUE OF Fc FOR THREE-BODY RATE REACTIONS (SEE REF 9, P.1145)
C Fc(T)       = Fc CALCULATED AS EXP(-T(K)/Fc(T)) 
C
C REACTION RATES HAVE THE FORM K = A * (300 / T)**B * EXP(C / T)
C
C *********************************************************************
C *                     READ PRELIMINARY COMMENTS                     * 
C *********************************************************************
C
C                   ----- REACTION RATE FORMAT -----
C
C A/D ORD AR BR CR N P Fc FcT COM    X  +Y  +Z  IV  =aA  +bB  +cC  +dD  +...  
C
C A/D
C   D       = REACTION IS DEAD. SKIP THIS REACTION. 
C   A       = REACTION ACTIVE AND INCLUDED IN ALL CHEMISTRY SETS
C             (URBAN, TROPSOSPHERIC AND STRATOSPHERIC SETS) 
C   U       = REACTION IN URBAN CHEMISTRY SET  
C   S       = REACTION IN STRATOSPHERIC CHEMISTRY SET  
C   T       = REACTION IN TROPOSPHERIC CHEMISTRY SET  
C   R       = REACTION IN TROPOSPHERIC AND URBAN CHEMISTRY SETS  
C   H       = REACTION IN TROPOSPHERIC AND STRATOSPHERIC CHEMISTRY SETS  
C ORD       = ORDINAL NUMBER OF REACTION
C AR,BR,CR  = RATE COEFFICIENTS: K = AR * (300/T) * BR * EXP( CR / T)
C AR        = DEFAULT PHOTORATE (S-1) FOR PHOTOPROCESSES
C NCOF      = DESCRIBED IN 'N COLUMN' ABOVE
C P         = DESCRIBED IN 'P COLUMN' ABOVE
C FCVT       = CHARACTERIZES FALLOFF CURVE IN PRESSURE-DEPENDENT REACTION
C FCT1T,2T  = EXPONENTS GIVING TEMPERATURE DEPENDENCE OF FCVT
C             FCVT = EXP(-T/FCT1)     OR  
C             FCVT = EXP(-T/FCT1) + EXP(-FCT2/T)   
C COM       = A9 AT THE END IS CURRENTLY FOR COMMENTS.
C X,Y,Z     = REACTANTS
C Z         = REACTANT OR 3RD BODY 'M' OR OTHER THIRD BODY
C I         = COEFFICIENT (INTEGER) FOR V 
C V         = REACT NOT INCLUDED IN REACT. RATE, BUT WHICH IS LOST IN REACTION.
C a,b,c,d.  = COEFFICIENTS FOR PRODUCTS (1,2,3,0.45,1.32, ETC (>=0.))
C A,B,C,D.. = PRODUCTS
C
C *********************************************************************
C *                           READ REACTIONS                          * 
C *********************************************************************
C
 102  READ(KGLC,21)  HEADING 
      IF (HEADING.NE.'BEGIN') GOTO 102  
C
 310  READ(KGLC,330) DINP,IORD,ARRT(1),BRRT(1),KCRRT(1),NCOF,SPECL(1),
     1               FCVT(1),FCT1T(1),FCT2T(1),COMMENT

!FP_ISOP - for debug (6/2009)

      IF (LPRT) THEN

      WRITE(*,*) DINP,IORD,ARRT(1),BRRT(1),KCRRT(1),NCOF,SPECL(1),
     &               FCVT(1),FCT1T(1),FCT2T(1),COMMENT

      ENDIF

      IF (NCOF+1.GT.MXCOF) THEN
       WRITE(6,155) NCOF+1, MXCOF, IORD
       CALL GEOS_CHEM_STOP
      ENDIF
C
      DO 350 I = 2, NCOF + 1
       READ(KGLC,330) JST,JORD,ARRT(I),BRRT(I),KCRRT(I),NDUM,SPECL(I),
     1                FCVT(I),FCT1T(I),FCT2T(I),COMMENT
 350  CONTINUE
C
      ! Now read 20 entries instead of 16 (bdf, bmy, 4/1/03)
      !FP_ISOP 24: now read 20 or 24 depending on input file (6/2009)
      ! NREAD is in comode.h (8/1/09)
      !READ(KGLC,332) (RINP(I),PINP(I),XINP(I),I=1,20)
      READ(KGLC,332) (RINP(I),PINP(I),XINP(I),I=1,NREAD)
C
 155  FORMAT('READCHEM: NCOF + 1 > MXCOF IN GLOBCHEM.DAT',3I4)
 330  FORMAT(A1,1X,I4,1X,ES8.2,1X,ES8.1,1X,I6,1X,I1,1X,A2,F6.2,1X,
     1       2(F6.0,1X),A20) 

      ! Increase format string to 14 products (bdf, 4/1/03)
 332  FORMAT(4(A1,0PF5.3,A14)/4(A1,0PF5.3,A14)/
     1     4(A1,0PF5.3,A14)/4(A1,0PF5.3,A14)/4(A1,0PF5.3,A14))
C
      IF (DINP.NE.'A'.AND.DINP.NE.'U'.AND.DINP.NE.'T'.AND.
     1    DINP.NE.'S'.AND.DINP.NE.'R'.AND.DINP.NE.'H') DINP = 'D'
C
      IF ((DINP.EQ.'U'.AND.IFURBAN.EQ.0).OR.
     1    (DINP.EQ.'T'.AND.IFTROP .EQ.0).OR.
     2    (DINP.EQ.'S'.AND.IFSTRAT.EQ.0).OR.
     3    (DINP.EQ.'R'.AND.IFURBAN.EQ.0.AND.IFTROP.EQ.0).OR.
     4    (DINP.EQ.'H'.AND.IFSTRAT.EQ.0.AND.IFTROP.EQ.0)) DINP = 'D'
C
      IF (XINP(1).EQ.'END KINETIC') THEN
C
       DO 323 NCS   = 1, NCSGAS 
        NRATES(NCS) = NTRATES(NCS)
 323   CONTINUE
C
       IDOPHOT     = 1
       GOTO 102 
      ELSEIF (XINP(1).EQ.'END PHOTOLYSIS') THEN
       GOTO 660 
      ELSEIF (DINP.EQ.'D') THEN
       GOTO 310 
      ENDIF
C
C *********************************************************************
C *      UPDATE REACTION NUMBER FOR REACTIONS THAT ARE USED           *
C *********************************************************************
C NRATCUR = CURRENT REACTION RATE NUMBER FOR A SET OF RATE COEFFICIENTS 
C NTRATES = CURRENT RATE COEFFICIENT NUMBER
C NALLRAT = COUNTS THE NUMBER OF ACTUAL REACTIONS
C SKIP URBAN         ('A', 'U', OR 'R')      REACTIONS IF NOT USED
C SKIP TROPOSPHERIC  ('A', 'T', 'R', OR 'H') REACTIONS IF NOT USED
C SKIP STRATOSPHERIC ('A', 'S', OR 'H')      REACTIONS IF NOT USED
C
      DO 325 NCS      = 1, NCSGAS
C
       NOUSE(NCS)     = 1
       IF (NCS.EQ.NCSALL  .AND.(DINP.EQ.'A'.OR.DINP.EQ.'U'.OR.
     1     DINP.EQ.'R'.OR.DINP.EQ.'S'.OR.DINP.EQ.'T'.OR.
     2     DINP.EQ.'H'))                NOUSE(NCS) = 0
       IF (NCS.EQ.NCSTRST .AND.(DINP.EQ.'A'.OR.DINP.EQ.'R'.OR.
     1     DINP.EQ.'T'.OR.DINP.EQ.'S'.OR.DINP.EQ.'H'))
     2                                  NOUSE(NCS) = 0
       IF (NCS.EQ.NCSURBAN.AND.(DINP.EQ.'A'.OR.DINP.EQ.'U'.OR.
     1     DINP.EQ.'R'))                NOUSE(NCS) = 0
       IF (NCS.EQ.NCSTROP .AND.(DINP.EQ.'A'.OR.DINP.EQ.'T'.OR.
     1     DINP.EQ.'R'.OR.DINP.EQ.'H')) NOUSE(NCS) = 0
       IF (NCS.EQ.NCSSTRAT.AND.(DINP.EQ.'A'.OR.DINP.EQ.'S'.OR.
     1     DINP.EQ.'H'))                NOUSE(NCS) = 0
C
       IF (NOUSE(NCS).EQ.0) THEN
        NK               = NTRATES(NCS) +  1
        NRATCUR(NCS)     = NK
        NALLRAT(NCS)     = NALLRAT(NCS) + 1
        NAR              = NALLRAT(NCS)
        NCEQUAT(NAR,NCS) = NK
C
        DO 320 I          = 1, NCOF + 1
         NTRATES(NCS)     = NTRATES( NCS) + 1
         NK1              = NTRATES( NCS)
         IRORD(  NK1,NCS) = IORD
         ARR(    NK1,NCS) = ARRT(I)
         BRR(    NK1,NCS) = BRRT(I)
         KCRR(   NK1,NCS) = KCRRT(I)
         FCV(    NK1,NCS) = FCVT( I)
         FCTEMP1(NK1,NCS) = FCT1T(I)
         FCTEMP2(NK1,NCS) = FCT2T(I)
 320    CONTINUE
       ENDIF
 325  CONTINUE

C *********************************************************************
C       SET UP A DEFAULT PHOTORATE (SEC-1), STORE ORDINAL NUMBER
C *********************************************************************
C
      IF (IDOPHOT.EQ.1) THEN
       NPHOTALL             = NPHOTALL      + 1
       ! record photalysis numbers for harvard-geos code (bdf, 4/18/03)
       NPHOT = NPHOTALL
       !DEFPRAT(NPHOTALL)    = ARRT(1)

       DO 640 NCS           = 1, NCSGAS
        IF (NOUSE(NCS).EQ.0) THEN
         NK                 = NRATCUR(NCS)
         DEFPRAT(NK,NCS)    = ARRT(1)
         JPHOTRAT(NCS)      = JPHOTRAT(NCS) + 1
         JPR                = JPHOTRAT(NCS)
         NKPHOTRAT(JPR,NCS) = NK
         NPPHOTRAT(JPR,NCS) = NPHOTALL
         JPHOTNK(  NK, NCS) = JPR
        ENDIF
 640   CONTINUE
      ENDIF
C
C *********************************************************************
C *  CHECK WHETHER EACH SPECIES SPOT IN INPUT REACTION SET IS ACTIVE, *
C *  INACTIVE, BLANK, DEAD, OR 'M'. STOP IF THE SPECIES IS NONE       *
C *  JNUM = -J = NON 'M' THIRD BODY IN PRESSURE-DEPENDENT REACTION.   *
C *********************************************************************
C
      DO 360 I = 1, NPRODHI
       IF (XINP(I).NE.' ') THEN
C
        IF (I.LE.NMREAC.AND.PINP(I).NE.0.) THEN
         WRITE(6,450) IORD
         CALL GEOS_CHEM_STOP
        ENDIF
C
        JNUM                      = 0
        IF (I.EQ.NMREAC.AND.RINP(3).EQ.' ') THEN
         IF (XINP(I).EQ.'M' ) JNUM = -9999
         IF (XINP(I).EQ.'O2') JNUM = -9998
         IF (XINP(I).EQ.'N2') JNUM = -9997
         IF (JNUM.LT.0)      GOTO 380
        ENDIF

        DO 370 J = 1, NTSPECGAS
         IF(XINP(I).EQ.NAMEGAS(J)) THEN
          IF (I.EQ.NMREAC.AND.RINP(3).EQ.' ') THEN
           JNUM   = -J
          ELSE
           JNUM   = J
          ENDIF
          GOTO 380
         ENDIF
 370    CONTINUE
C
        IF (I.GT.NMREAC) THEN
         DO 390 J = 1, NSDEAD 
          IF (XINP(I).EQ.NAMD(J)) GOTO 360 
 390     CONTINUE
        ENDIF
C
        WRITE(6,400) IORD, XINP(I)
        CALL GEOS_CHEM_STOP
C
 380    DO 410 NCS        = 1, NCSGAS
         IF (NOUSE(NCS).EQ.0) THEN 
          NK               = NRATCUR(NCS) 
          IRM(I,NK,NCS)    = JNUM 
          NPRODUC(NK,NCS)  = I 
          IF (PINP(I).EQ.0.) THEN
           FKOEF(I,NK,NCS) = 1.0d0
          ELSE 
           FKOEF(I,NK,NCS) = PINP(I)
          ENDIF 
         ENDIF 
 410    CONTINUE 
C
        IF (IDOPHOT.EQ.1) NAMEPHOT(I,NPHOTALL) = XINP(I)
       ENDIF
 360  CONTINUE
C
      DO 415 NCS           = 1, NCSGAS 
       IF (NOUSE(NCS).EQ.0) THEN 
        NK                 = NRATCUR(NCS) 
        IF (IRM(1,NK,NCS).EQ.0.OR.(IRM(3,NK,NCS).GT.0.AND.
     1      IRM(2,NK,NCS).EQ.0)) THEN
         WRITE(6,430) IORD
         CALL GEOS_CHEM_STOP
        ENDIF
       ENDIF
 415  CONTINUE
C
C *********************************************************************
C *             PLACE SPECIAL-RATE INFORMATION INTO ARRAYS            *
C *********************************************************************
C
C NMAIR   = # REACTIONS WHERE THE SPECIES IN THE THIRD POSITION
C           IS 'M' = 'O2 + N2': IRM(3,NK) = -9999 
C NMO2    = # REACTIONS WHERE THE SPECIES IN THE THIRD POSITION
C           IS O2:  IRM(3,NK) = -9998  
C NMN2    = # REACTIONS WHERE THE SPECIES IN THE THIRD POSITION
C           IS N2:  IRM(3,NK) = -9997  
C NM3BOD  = # REACTIONS WHERE THE SPECIES IN THE THIRD POSITION
C           IS ANY OTHER SPECIES:  IRM(3,NK) = -SPECIES NUMBER
C NNEQ    = # THERMALLY DISSOCIATING EQUILIBRIUM REACTIONS. PREVIOUS
C             EQUATION MUST BE PRESSURE-DEPENDENT. (SPECIAL = 'E')
C             THUS: REQUIRES 3 REACTIONS TOTAL (2 FOR PRESS. DEP, 1 EQ.)
C NPRESM  = # PRESSURE DEPENDENT 3-BODY REACTIONS  (SPECIAL = 'P')
C

! bdf smv2 count number of emission and drydep reactions. Used in calcrate
!  to put emissions into the chemistry, and drydep out of chemistry.


      DO NCS=1,NCSGAS
         IF (NOUSE(NCS) .EQ. 0 ) THEN
            NK                 = NRATCUR(NCS)
            IF (XINP(1).EQ.'EMISSION') THEN
               NEMIS(NCS)           = NEMIS(NCS) + 1
               IF (NEMIS(NCS)   .GT. MAXGL3) THEN
                  WRITE(*,*) 'ERROR NEMIS   Greater then MAXGL3',
     x                 ' NEMIS(NCS)   = ',NEMIS(NCS), 'MAXGL3 = ',MAXGL3
                  WRITE( 6, '(a)' ) 'STOP 124'
                  CALL GEOS_CHEM_STOP
               ENDIF
               NKEMIS(NEMIS(NCS),NCS)   = NK
            ENDIF
         ENDIF
C
C  NDRYDEP  = number of dry deposition reactions read in
C  NKDRY    = reaction numbers of dry deposition reactions
         IF (NOUSE(NCS) .EQ. 0) THEN
            NK                 = NRATCUR(NCS)
            IF (XINP(NPRODLO).EQ.'DRYDEP') THEN
               NDRYDEP(NCS)           = NDRYDEP(NCS) + 1
               IF (NDRYDEP(NCS) .GT. MAXDEP) THEN
                  WRITE(*,*) 'ERROR NDRYDEP Greater then MAXDEP',
     x                 ' NDRYDEP(NCS)=',NDRYDEP(NCS),'MAXDEP=',MAXDEP
                  WRITE( 6, '(a)' ) 'STOP 125'
                  CALL GEOS_CHEM_STOP
               ENDIF
               NKDRY(NDRYDEP(NCS),NCS)   = NK
            ENDIF
         ENDIF

! bdf smv2 use Q to flag O3 photolysis, code is not confused by 'A''s
         IF (NOUSE(NCS) .EQ. 0) THEN
            NK                 = NRATCUR(NCS)
            IF (SPECL(1).EQ.'Q') NKO3PHOT(NCS)=NK !Flag O3 photolysis
            IF (SPECL(1).EQ.'T') NKHNO4(NCS)  =NK !Flag HNO4 photolysis (gcc)
            IF (SPECL(1).EQ.'I') NKHOROI(NCS) = NK !Flag CH2O-producing branch in EP photolysis
            IF (SPECL(1).EQ.'J') NKHOROJ(NCS) = NK !Flag GLYC-producing branch in EP photolysis

!FLAG HACET AND GLYC reaction to modify yield according to temperature
! (fp, 8/09)
            IF (SPECL(1).EQ.'N') NKGLYC(NCS,1)=NK !Flag GLYC 
            IF (SPECL(1).EQ.'O') NKGLYC(NCS,2)=NK !Flag GLYC 
            IF (SPECL(1).EQ.'F') NKHAC(NCS,1)  =NK !Flag HAC photolysis
            IF (SPECL(1).EQ.'L') NKHAC(NCS,2)  =NK !Flag HAC photolysis

!FP_ISOP FLAG FOR MCO3
            IF (SPECL(1).EQ.'DA') NKMCO3(NCS,1)=NK 
            IF (SPECL(1).EQ.'DB') NKMCO3(NCS,2)=NK 
            IF (SPECL(1).EQ.'DC') NKMCO3(NCS,3)=NK 

         ENDIF
      ENDDO

      IF (IDOPHOT.EQ.0) THEN
         IF ((SPECL(1).EQ.'V'.AND.NCOF.NE.1).OR.
     1        (SPECL(1).EQ.'W'.AND.NCOF.NE.1).OR.
     2        (SPECL(1).EQ.'X'.AND.NCOF.NE.2).OR.  
     3        (SPECL(1).EQ.'Y'.AND.NCOF.NE.0).OR.  
     4        (SPECL(1).EQ.'Z'.AND.NCOF.NE.1).OR.  
     5        (SPECL(1).EQ.'P'.AND.NCOF.NE.1).OR.  
     6        (SPECL(1).EQ.'E'.AND.NCOF.NE.0).OR.
     7        (SPECL(1).EQ.'S'.AND.NCOF.NE.0)) THEN
            WRITE(6,440) IORD, SPECL(1), NCOF
            CALL GEOS_CHEM_STOP
         ENDIF 
C
         DO 420 NCS           = 1, NCSGAS 
            IF (NOUSE(NCS).EQ.0) THEN 
               NK                 = NRATCUR(NCS) 
C
               ITHIRDB            = IRM(3,NK,NCS)
C
               IF (ITHIRDB.EQ.-9999) THEN
                  NMAIR(NCS)        = NMAIR(NCS) + 1
                  NM2               = NMAIR(NCS) 
                  NREACAIR(NM2,NCS) = NK 
               ELSEIF (ITHIRDB.EQ.-9998) THEN  
                  NMO2(NCS)         = NMO2(NCS)  + 1
                  NM2               = NMO2(NCS) 
                  NREACO2(NM2,NCS)  = NK 
               ELSEIF (ITHIRDB.EQ.-9997) THEN
                  NMN2(NCS)         = NMN2(NCS)  + 1
                  NM2               = NMN2(NCS) 
                  NREACN2(NM2,NCS)  = NK 
               ELSEIF (ITHIRDB.LT.0) THEN
                  NM3BOD(NCS)       = NM3BOD(NCS) + 1
                  NM2               = NM3BOD(NCS)
                  NREAC3B(NM2,NCS)  = NK
                  LGAS3BOD(NM2,NCS) = -ITHIRDB 
               ENDIF
C
               IF (SPECL(1).EQ.'P') THEN
                  NPRESM(NCS)       = NPRESM(NCS) + 1  
                  NR                = NPRESM(NCS) 
                  NREACPM(NR,NCS)   = NK
               ELSEIF (SPECL(1).EQ.'E') THEN
                  NNEQ(NCS)         = NNEQ(NCS)   + 1
                  NN                = NNEQ(NCS)
                  NREACEQ(NN,NCS)   = NK
C
C EQUILIBRIUM REACTIONS USE THE PREVIOUS REACTION AS PART OF THE 
C RATE CALCULATION (SEE CALCRATE.F). THE PREVIOUS REACTION MAY BE
C PRESSURE DEPENDENT.
C
                  NREQOTH(NN,NCS)   = NCEQUAT(NALLRAT(NCS)-1,NCS)
               ENDIF
C
C NKSPECV = SPECIAL REACTION CH3SCH3 + OH = CH3S(OH)CH3 (SPECL = 'V') 
C NKSPECW = SPECIAL REACTION O(1D)   + O2,N2            (SPECL = 'W') 
C NKSPECX = SPECIAL REACTION OH      + HNO3             (SPECL = 'X')
C NKSPECY = SPECIAL REACTION OH      + CO               (SPECL = 'Y')
C NKSPECZ = SPECIAL REACTION HO2     + HO2              (SPECL = 'Z') 
C
               ! bdf smv2 'V' reaction has a special rate.  
               ! More than one reaction of this type
               IF (SPECL(1).EQ.'V') THEN
                  NNADDV(NCS)                = NNADDV(NCS)+1
                  NKSPECV( NNADDV(NCS),NCS ) = NK
               ENDIF

               ! Added for DMS+OH+O2 rxn (bdf, bmy, 4/18/03)
               IF (SPECL(1).EQ.'G') THEN
                  NNADDG(NCS)                = NNADDG(NCS)+1
                  NKSPECG( NNADDG(NCS),NCS ) = NK
               ENDIF

               ! add flag for wet dep reaction (bdf, bmy, 4/18/03)
               IF (SPECL(1).EQ.'K')  THEN
                  NNADDK(NCS)               = NNADDK(NCS) + 1
                  NKSPECK( NNADDK(NCS),NCS) = NK

                  ! Also denote N2O5 hydrolysis rxn (mje, bmy, 8/7/03)
                  IF ( XINP(1) == 'N2O5' ) THEN
                     NKN2O5 = NK
                  ENDIF

                  ! Same for HO2 hydrolysis rxn (jaegle, 02/26/09)
                  IF ( XINP(1) == 'HO2' ) THEN
                     NKHO2 = NK
                  ENDIF

               ENDIF
C

               ! DA, DB, and DC are now in use by FP
               ! This part applies only to 'D' (hotp 8/1/09)
               IF (SPECL(1).EQ.'D' .AND. SPECL(2).NE.'A' .AND.
     &             SPECL(2).NE.'B' .AND. SPECL(2).NE.'C')  THEN
                  NNADDD(NCS)               = NNADDD(NCS) + 1
                  NKSPECD( NNADDD(NCS),NCS) = NK
               ENDIF
C
               IF (SPECL(1).EQ.'A')  THEN
                  NNADDA(NCS)               = NNADDA(NCS) + 1
                  NKSPECA( NNADDA(NCS),NCS) = NK
               ENDIF
C
               IF (SPECL(1).EQ.'B')  THEN
                   NNADDB(NCS)              = NNADDB(NCS) + 1
                  NKSPECB( NNADDB(NCS),NCS) = NK
               ENDIF
C
               IF (SPECL(1).EQ.'C')  THEN
                   NNADDC(NCS)               = NNADDC(NCS) + 1
                   NKSPECC( NNADDC(NCS),NCS) = NK
               ENDIF

               ! F: HOC2H4O ------> HO2 + 2CH2O
               IF (SPECL(1).EQ.'F') THEN
                   NNADDF(NCS)               = NNADDF(NCS) + 1
                   NKSPECF( NNADDF(NCS),NCS) = NK
               ENDIF

               ! H: HOC2H4O --O2--> HO2 + GLYC
               IF (SPECL(1).EQ.'H') THEN
                   NNADDH(NCS)               = NNADDH(NCS) + 1
                   NKSPECH( NNADDH(NCS),NCS) = NK
               ENDIF

               IF (SPECL(1).EQ.'W') THEN
                  NKSPECW(NCS) = NK
               ENDIF
               IF (SPECL(1).EQ.'X') THEN
                  NKSPECX(NCS) = NK
               ENDIF
               IF (SPECL(1).EQ.'Y') THEN
                  NKSPECY(NCS) = NK
               ENDIF
               IF (SPECL(1).EQ.'Z') THEN
                  NKSPECZ(NCS) = NK
               ENDIF
C
C *********************************************************************
C *                       SURFACE REACTIONS                           * 
C *********************************************************************
C ARR(INIT)         = REACTION PROBABILITY 
C ARR(FINAL)        = REACTION PROBABILITY * QTHERMG  
C QTHERMG * SQRT(T) = (1/4) * THERMAL VELOCITY OF GAS (CM S-1)
C
               IF (SPECL(1).EQ.'S') THEN
                  NSURFACE(NCS)    = NSURFACE(NCS) + 1
                  NS2              = NSURFACE(NCS) 
                  JGAS1            = IRM(1,NK,NCS) 
                  JGAS2            = IRM(2,NK,NCS) 
                  JGAS3            = IRM(3,NK,NCS) 
                  QTHERMG  = 0.25d0*SQRT(EIGHTDPI*RSTARG/WTGAS(JGAS1)) 
                  ARR(NK,NCS)      = ARR(NK,NCS) * QTHERMG 
                  NKSURF(NS2)      = NK 
                  NCOATG(NS2)      = JGAS2  
C
                  IF (JGAS3.NE.0) THEN
                     WRITE(6,470) NK 
                     CALL GEOS_CHEM_STOP
                  ENDIF  
               ENDIF 
C
C *********************************************************************
C *       SET ARRAYS FOR CALCULATING REACTION RATES EFFICIENTLY       * 
C *********************************************************************
C NARR = NUMBER OF REACTIONS OF THE FORM K = A
C NABR = NUMBER OF REACTIONS OF THE FORM K = A * (300 / T)**B
C NACR = NUMBER OF REACTIONS OF THE FORM K = A * EXP(C / T)
C NABC = NUMBER OF REACTIONS OF THE FORM K = A * (300 / T)**B * EXP(C / T)
C NKARR, NKBRR, NKACR, NKABC = REACTION RATE NUMBERS OF EACH
C NARR,  NABR,  NACR,  NABC  REACTION 
C


               NK1              = NK - 1 

               DO 425 I         = 1, NCOF + 1
                  NK1             = NK1 + 1
                  IF (KCRR(NK1,NCS).EQ.0) THEN
                     IF (BRR(NK1,NCS).EQ.0.) THEN
                        NARR(NCS)     = NARR(NCS) + 1
                        NA            = NARR(NCS) 
                        NKARR(NA,NCS) = NK1 
                     ELSE
                        NABR(NCS)     = NABR(NCS) + 1
                        NA            = NABR(NCS) 
                        NKABR(NA,NCS) = NK1 
                     ENDIF
                  ELSE
                     IF (BRR(NK1,NCS).EQ.0.) THEN
                        NACR(NCS)     = NACR(NCS) + 1
                        NA            = NACR(NCS) 
                        NKACR(NA,NCS) = NK1 
                     ELSE
                        NABC(NCS)     = NABC(NCS) + 1
                        NA            = NABC(NCS) 
                        NKABC(NA,NCS) = NK1 
                     ENDIF
                  ENDIF
 425           CONTINUE
C
            ENDIF
C       ENDIF NOUSE(NCS).EQ.0
C
 420     CONTINUE
      ENDIF
C     ENDIF IDOPHOT.EQ.0
C
      GOTO 310 
C
 400  FORMAT('INVALID REACT',I4,' W UNRECOGNIZABLE OR DEAD SPEC ',A14,
     1       'ALL REACTANTS MUST BE ACTIVE/INACTIVE. PRODS CAN BE DEAD') 
 430  FORMAT('READCHEM:REACT ',I3,' 1ST SPOT EMPTY OR 3RD SPOT FILLED ', 
     1       ' BUT 2ND EMPTY')
 440  FORMAT('READCHEM: REACT ',I3,'OR BEFORE: SPECIAL REACTION WITH ',
     1       'DELIMETER ',A2,' HAD INCORRECT # OF REACTIONS ',I5)
 450  FORMAT('READCHEM: ORD# REACT ',I3,' CANT HAVE COEFF > 1')
 470  FORMAT('READCHEM: SURFACE REACTION ',I5,'HAS THREE REACTANTS ')
 510  FORMAT(I3,1X,ES7.1,1X,ES7.1,I6,1X,0PF3.2,1X,
     1       A6,2(A1,A6),14(A1,0PF3.1,A6))  
 520  FORMAT( 'KINETIC REACTIONS FOR ', A,' CHEMISTRY',/,   
     1        'RATE CONSTANTS HAVE FORM K = A * (300/T)**B * EXP(C/T).')
 521  FORMAT( 'NMBR   A       B     C     Fv       REACTION' )
 525  FORMAT( 'PHOTOPROCESS REACTIONS FOR ', A,' CHEMISTRY' )   
 526  FORMAT( 'NMBR   DEFP (S-1)                        REACTION' )
C
C *********************************************************************
C *                  PRINT OUT REACTION INFORMATION                   *
C *********************************************************************
C
 660  IF (IOREAC.EQ.1) THEN 
       DO 502 NCS    = 1, NCSGAS 

        ! Write reaction header
        WRITE( IO93, '(/,a)' ) REPEAT( '=', 79 )
        WRITE( IO93, 520     ) TRIM( CHEMTYP(NCS) )
        WRITE( IO93, '(a,/)' ) REPEAT( '=', 79 )
        WRITE( IO93, 521     )

        DO 500 NK    = 1, NTRATES(NCS) 

         ! Write photo rxn header
         IF ( NK .EQ. NRATES(NCS)+1 ) THEN
            WRITE( IO93, '(/,a)' ) REPEAT( '=', 79 ) 
            WRITE( IO93, 525     ) TRIM( CHEMTYP(NCS) )
            WRITE( IO93, '(a,/)' ) REPEAT( '=', 79 ) 
            WRITE( IO93, 526     )
         ENDIF

         DO 490 I    = 1, NPRODHI
          RINP(I)    = '+'
          PINP(I)    = FKOEF(I,NK,NCS)  
          JGAS       = IRM(I,NK,NCS)
          IF (JGAS.GE.0)                     XINP(I) = NAMEGAS(JGAS) 
          IF (JGAS.EQ.-9999)                 XINP(I) = 'M' 
          IF (JGAS.EQ.-9998)                 XINP(I) = 'O2' 
          IF (JGAS.EQ.-9997)                 XINP(I) = 'N2' 
          IF (JGAS.LT.0.AND.JGAS.GT.-NTSPECGAS) XINP(I) = NAMEGAS(-JGAS) 
 490     CONTINUE
C
         RINP(5)    = '='
         WRITE(IO93,510) NK,ARR(NK,NCS),BRR(NK,NCS),KCRR(NK,NCS),
     1                FCV(NK,NCS),XINP(1),'+',XINP(2),
     2                '+',XINP(3),(RINP(I),PINP(I),XINP(I),
     3                I = 5,NPRODUC(NK,NCS))
 500    CONTINUE 
 502   CONTINUE 
      ENDIF

C
C *********************************************************************
C ***********************  CHECK SOME DIMENSIONS  *********************
C *********************************************************************
C
      DO 670 NCS = 1, NCSGAS  
       IF (NTRATES(NCS) .GT. NMTRATE  .OR. NPHOTALL   .GT. IPHOT   .OR. 
     1     NTSPECGAS    .GT. IGAS     .OR. NSDEAD     .GT. NMDEAD  .OR.
     2     NPRODHI      .GT. NMRPROD .OR.
     3     NMAIR(NCS)   .GT. MAXGL3   .OR. NMO2(NCS)  .GT. MAXGL3  .OR.
     4     NMN2(NCS)    .GT. MAXGL2   .OR. NPRESM(NCS).GT. MAXGL2  .OR.
     5     NSURFACE(NCS).GT. MAXGL4   .OR. NM3BOD(NCS).GT. MAXGL3) THEN
C
        WRITE(6,680)
     1  NMTRATE,NTRATES(NCS), IPHOT  , NPHOTALL,     IGAS  ,NTSPECGAS,
     2  NMDEAD ,NSDEAD,       NMRPROD, NPRODHI,
     3  MAXGL3 ,NMAIR(NCS),   MAXGL3 , NMO2(NCS),    MAXGL2,NMN2(NCS),
     4  MAXGL2 ,NPRESM(NCS),  MAXGL4 , NSURFACE(NCS),MAXGL3,NM3BOD(NCS)  
        CALL GEOS_CHEM_STOP
       ENDIF
 670  CONTINUE 
C
 680  FORMAT('ONE OF THE DIMENSIONS BELOW IS TOO SMALL:',/,
     1       'DIMENSION: NMTRATE  =  ',I4,' VARIABLE: NTRATES  = ',I4/  
     2       'DIMENSION: IPHOT    =  ',I4,' VARIABLE: NPHOTALL = ',I4/  
     3       'DIMENSION: IGAS     =  ',I4,' VARIABLE: NTSPECGS = ',I4/  
     4       'DIMENSION: NMDEAD   =  ',I4,' VARIABLE: NSDEAD   = ',I4/  
     6       'DIMENSION: NMRPROD  =  ',I4,' VARIABLE: NPRODHI  = ',I4/
     7       'DIMENSION: MAXGL3   =  ',I4,' VARIABLE: NMAIR    = ',I4/ 
     8       'DIMENSION: MAXGL3   =  ',I4,' VARIABLE: NMO2     = ',I4/
     9       'DIMENSION: MAXGL2   =  ',I4,' VARIABLE: NMN2     = ',I4/ 
     1       'DIMENSION: MAXGL2   =  ',I4,' VARIABLE: NPRESM   = ',I4/ 
     2       'DIMENSION: MAXGL4   =  ',I4,' VARIABLE: NSURFACE = ',I4/ 
     3       'DIMENSION: MAXGL3   =  ',I4,' VARIABLE: NM3BOD   = ',I4)
C
C *********************************************************************
C ************************** SET KEY PARAMETERS ***********************
C *********************************************************************
C
      DO 702 NCS          = 1, NCSGAS 
       NSPEC(NCS)         = NGAS
       NTSPEC(NCS)        = NTSPECGAS
C
       DO 700 JGAS        = 1, NTSPECGAS
        NAMENCS(JGAS,NCS) = NAMEGAS(JGAS)  
        QBKCHEM(JGAS,NCS) = QBKGAS( JGAS)
 700   CONTINUE 
 702  CONTINUE 

!---smv2-s
! Update (gcc)
C bdf smv2, put this in here for now.
C *********************************************************************
C ****************   READ INFO FOR AEROSOL REACTIONS   ****************
C *********************************************************************

C astkcf  -- sticking coefficient (no unit), order of 0.1
C xgdfcf  -- gas phase diffusion coefficient (cm2/s), order of 0.1
C iarsfa  -- fortran unit number for reading sulfate abundance file
C mwarsl  -- aerosol molecular wright (g/mol)    [H2SO4=98]
C ruarsl  -- density of aerosol (g/cc)
C RH100   -- deliquescence point, relative humidity below which we
C            have no wet aerosols
C
	OPEN(7,FILE='chemga.dat',FORM='FORMATTED',STATUS='OLD')
	READ(7,*)
	READ(7,610) ASTKCF
	READ(7,*)
        READ(7,610)
	READ(7,*)
	READ(7,620) MWARSL
	READ(7,610) RUARSL
	READ(7,630) RH100
	READ(7,620) IARSFA
	CLOSE(7)
 610	FORMAT(E10.3)
 620  FORMAT(I10)
 630  FORMAT(F10.2)

C
C *********************************************************************
C ***** CALL JSPARSE TO SET ARRAYS FOR SOLVING CHEMICAL EQUATIONS ***** 
C *********************************************************************
C
      ! Call SETPL to setup ND65 prod/loss diagnostic
      ! SETPL must be called before JSPARSE (ljm, bmy, 5/9/03)
      IF ( LFAMILY ) CALL SETPL

      ! IFSOLVE = 1 means we are calling the chemistry solver
      IF ( IFSOLVE .EQ. 1 ) THEN

         ! Loop over chemistry regimes (for now NCSGAS=NCSURBAN=1)
         DO NCS = 1, NCSGAS

            ! Set up sparse matrix stuff
            CALL JSPARSE

            !===========================================================
            ! Determine which species are ND65 families and which are 
            ! not.  Do this once (after JSPARSE) & store in the lookup
            ! table ITS_NOT_A_ND65_FAMILY. (bmy, 7/9/03)
            !===========================================================

            ! Loop over all species
            DO J = 1, ISCHANG(NCS)
               
               ! Initialize lookup table
               ITS_NOT_A_ND65_FAMILY(J) = .TRUE.

               ! Test if species J is a ND65 prodloss family
               ! MAPPL is the reordered species index after JSPARSE
               DO N = 1, NFAMILIES
                  IF ( J == MAPPL(IFAM(N),NCS) ) THEN
                     ITS_NOT_A_ND65_FAMILY(J) = .FALSE.
                     EXIT
                  ! dkh 
                  ELSEIF ( J == MAPPL(ILBRO2H,NCS) .or.
     &                     J == MAPPL(ILBRO2N,NCS) .or.
     &                     J == MAPPL(ILTRO2H,NCS) .or.
     &                     J == MAPPL(ILTRO2N,NCS) .or.
     &                     J == MAPPL(ILXRO2H,NCS) .or.
     &                     J == MAPPL(ILXRO2N,NCS) ) THEN
                     ITS_NOT_A_ND65_FAMILY(J) = .FALSE.
                     EXIT
                  ENDIF
               ENDDO  
            ENDDO
         ENDDO

      ELSE

         ! If we are not calling the chemistry solver, then
         ! set number of active gas photoprocesses to zero
         NPHOTALL = 0

      ENDIF
C
C *********************************************************************
C *******************   END OF SUBROUTINE READCHEM   ****************** 
C *********************************************************************
C     

      RETURN
      END SUBROUTINE READCHEM
