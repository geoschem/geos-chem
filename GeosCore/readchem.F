      SUBROUTINE READCHEM( am_I_Root, Input_Opt, RC )
!
!******************************************************************************
!  Subroutine READCHEM reads species 2names, chemical rxns, and photolysis 
!  reactions from the "globchem.dat" chemistry mechanism file for SMVGEAR II.  
!  (M. Jacobson 1997; bdf, bmy, 5/9/03, 6/1/10)
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
!  (10) In "510 FORMAT", the format of 'B' in kinetic reactions section of 
!        smv2.log does not match variable type and causes code to crash when 
!        run with ifort -check all flag. (hotp, bmy, 6/1/10)
!  13 Jun 2011 - R. Yantosca - Bug fix: prevent out-of-bounds errors in MAPPL
!                              array.  Test species flags to make sure they
!                              are nonzero before using them.
!  30 Jul 2012 - R. Yantosca - Now accept am_I_Root as an argument even when
!                              called from the traditional driver "main.F"
!  03 Aug 2012 - R. Yantosca - Move calls to findFreeLUN out of DEVEL block
!  07 Aug 2012 - R. Yantosca - Now print LUN used to open file
!  19 Mar 2013 - R. Yantosca - When using ESMF interface to GEOS-5, append
!                              ".rc" to input.geos (instead of ___.rc)
!  13 Aug 2013 - M. Sulprizio- Add modifications for updated SOA and SOA +
!                              semivolatile POA simulations (H. Pye)
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!******************************************************************************
!
! !USES:
!
      USE DRYDEP_MOD,         ONLY : MAXDEP
      USE ERROR_MOD,          ONLY : GEOS_CHEM_STOP
      USE DIAG_PL_MOD,        ONLY : SETJFAM, SETPL
      USE ERROR_MOD,          ONLY : DEBUG_MSG
      USE CMN_SIZE_MOD
      USE COMODE_LOOP_MOD
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE inquireMod,         ONLY : findFreeLUN

!#if defined( DEVEL ) || defined ( ESMF_ )
!      USE GIGC_MPI_WRAP
!endif

      IMPLICIT NONE
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)  :: am_I_Root   ! Is this the root CPU?
      TYPE(OptInput), INTENT(IN)  :: Input_Opt   ! Input Options object
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT) :: RC          ! Success or failure?
!
! *********************************************************************
! ************        WRITTEN BY MARK JACOBSON (1990-4)    ************
! ***                         (650) 723-6836                        *** 
! *********************************************************************
!
!  RRRRRR   EEEEEEE     A     DDDDDD  CCCCCCC  H     H  EEEEEEE M     M  
!  R     R  E          A A    D     D C        H     H  E       M M M M 
!  RRRRRR   EEEE      A   A   D     D C        HHHHHHH  EEEEEEE M  M  M  
!  R  R     E        AAAAAAA  D     D C        H     H  E       M     M
!  R   R    EEEEEEE A       A DDDDDD  CCCCCCC  H     H  EEEEEEE M     M 
!   
! *********************************************************************
! *  THIS IS THE SETUP ROUTINE FOR GAS-PHASE CHEMISTRY. IT READS      *
! *  SPECIES NAMES, CHEMICAL REACTIONS, AND PHOTOPROCESSES FROM AN    *
! *  INPUT DATA SET. IT THEN PLACES ALL NECESSARY INFORMATION INTO    *
! *  ARRAYS AND PRINTS OUT THE INPUT INFORMATION.                     *
! *                                                                   *
! * HOW TO CALL SUBROUTINE:                                           *
! * ----------------------                                            *
! *  CALL READCHEM.F FROM MAIN.F                                      *
! *                                                                   *
! *********************************************************************
!
! *********************************************************************
! *                  SOME PARAMETER DEFINITIONS                       * 
! *********************************************************************
!
! DEFPRAT   = DEFAULT PHOTORATE (SEC-1)
! IGAS      = DIMENSION OF MAXIMUM NUMBER OF GAS SPECIES, ACTIVE + INACTIVE.
! IPHOT     = MAXIMUM NUMBER OF RADIATIVELY ACTIVE SPECIES 
! IPORD     = ORDINAL # OF PHOTOPROCESS (USED TO IDENTIFY REACTION)
! IRORD     = ORDINAL # OF KINET. REACT. (USED TO IDENTIFY REACTION) 
! JMBCOMP   = SPECIES NUMBER FOR EACH SPECIES IN A MASS BAL. GROUP
! MBCOMP    = COUNTS THE MASS BALANCE SPECIES IN EACH M.B. GROUP (I.E.
!             SULFUR IS A M.B. GROUP.
! NACTIVE   = NUMBER OF ACTIVE SPECIES READ IN -- (A) IN COLUMN ONE
!             OF INPUT DATA SET, CONVERTED TO NSPEC LATER 
! NALLREAC  = TOTAL NUMBER OF REACTANT POSITIONS IN A REACTION (BUT
!             NUMBER OF ACTIVE POSITIONS IN NMREAC) 
! NAMD      = NAMES OF SPECIES WHICH MAY APPEAR IN REACTIONS BUT WHICH ARE
!            "DEAD" WITH RESPECT TO THE PHOTOCHEMISTRY AND THUS ARE NOT
!             PRINTED OUT.
! NAMEGAS   = CHARACTER ARRAY OF SPECIES NAMES. 
! NAMENCS   = CHARACTER ARRAY OF SPECIES NAMES. 
! NCS       = 1..NCSGAS FOR GAS CHEMISTRY                            
! NGAS      = NSPEC, THE NUMBER OF ACTIVE SPECIES
! NINAC     = NUMBER OF INACTIVE SPECIES READ -- (I) IN COLUMN ONE 
! NMAIR     = # REACTIONS WHERE THE SPECIES IN THE THIRD POSITION
!             IS 'M' = 'O2 + N2' 
! NMASBAL   = NUMBER OF MASS BALANCE GROUPS (E.G. S, N, C ARE GROUPS)
! NMN2      = # REACTIONS WHERE SPECIES IN THE THIRD POSITION IS N2
! NMO2      = # REACTIONS WHERE SPECIES IN THE THIRD POSITION IS O2
! NM3BOD    = # REACTIONS WHERE SPECIES IN THE THIRD POSITION
!             IS ANY OTHER SPECIES:  IRM(3,NK) = -SPECIES NUMBER
! NMPROD    = MAXIMUM NUMBER OF ACTIVE PRODUCTS IN A REACTION (READER.F)
! NMREAC    = MAXIMUM NUMBER OF ACTIVE REACTANTS IN A REACTION
!             OF INPUT DATA SET. (DEFINED IN READER.F)   
! NNEQ      = # THERMALLY DISSOCIATING EQUILIBRIUM REACTIONS. PREVIOUS
!             EQUATION MUST BE PRESSURE-DEPENDENT.
! NPHOTALL  = NUMBER OF ACTIVE GAS PHOTOPROCESSES
! NPRESM    = # PRESSURE DEPENDENT 3-BODY REACTIONS  
! NPRODHI   = HIGHEST PRODUCT POSITON IN A REACTION
! NPRODLO   = LOWEST PRODUCT POSITON IN A REACTION
! NRATES    = NUMBER OF ACTIVE REACTIONS (EXCLUDING PHOTPROCESSES)
! NSDEAD    = NUMBER OF DEAD SPECIES READ IN -- (D) IN COLUMN ONE OF
!             INPUT DATA SET.
! NSPEC     = TOTAL NUMBER OF ACTIVE SPECIES.
! NTRATES   = NUMBER OF ACTIVE KINETIC REACTIONS PLUS PHOTOPROCESSES
! NTSPEC    = ACTUAL NUMBER OF ACTIVE + INACTIVE (BUT NOT DEAD) SPECIES. 
! QBKGAS    = DEFAULT BACKGROUND CONCENTRATION (VOL MIXING RATIO) 
! RINP      = CHARACTER ARRAY FOR READING IN INFORMATION FROM DATA SETS.
! WTMB      = MASS BALANCE WEIGHT FOR EACH M. B. SPECIES
! XINP      = CHARACTER ARRAY FOR READING IN INFORMATION FROM DATA SETS.
!
! *********************************************************************
! * SET INITIAL VALUES AND READ INITIAL COMMENTS FROM INPUT DATA SET  *
! *********************************************************************
!
      INTEGER NINAC,NACTIVE,NSDEAD,NOTHGS,IDOPHOT,I,NMBGAS,NM,JGAS
      INTEGER MB,MBP,INACT1,JGAS1,J,IORD,NCOF,JORD,NDUM,NK,NAR,NK1
      INTEGER JPR,JNUM,ITHIRDB,NM2,NR,NN,JGAS2,JGAS3,NA,N,NS2

      REAL*8 C1,CSTRAT,CTROPL,CTROPS,CURBAN,QTHERMG

      ! ININT used to be defined w/in "comode.h", but it is only ever used
      ! w/in "readchem.f".  Declare here as a local variable. (bmy, 7/16/03)
      INTEGER            :: ININT(10)  

      CHARACTER(LEN=255) :: FILENAME
      LOGICAL            :: LPRT

      INTEGER             :: IU_CHEMDAT ! GLOBCHEM.DAT LUN

!#if defined( DEVEL )
!      INTEGER, ALLOCATABLE :: INTVAL(:)
!#endif

      !=================================================================
      ! READCHEM begins here!
      !=================================================================

      ! Assume success
      RC             = GIGC_SUCCESS

      ! Copy logical fields from INPUT_OPT to local variables for use below
      LPRT           = Input_Opt%LPRT

      NINAC          = 0
      NACTIVE        = 0 
      NSDEAD         = 0 
      NOTHGS         = 0 
      NPHOTALL       = 0 
      IDOPHOT        = 0
!
      NAMEGAS(0)     = ' '
!
      ! Initialize flag for N2O5 hydrolysis rxn (bmy, 8/7/03)
      NKN2O5         = 0

      ! Initialize flag for HO2 hydrolysis rxn (jaegle, 02/26/09)
      NKHO2          = 0

      ! Initialize flag for BrNO3 hydrolysis rxn (jpp, 5/4/10)
      NKBrNO3        = 0
      ! Initialize flags for HOBr + HBr het. pseudo-rxns (jpp, 3/22/11)
      NKHOBr         = 0
      NK1HBr         = 0

      DO 44 I            = 1, NMDEAD
       NAMD(I)           = ' '
 44   CONTINUE
!
      DO 46 I            = 1, IGAS
       NAMEGAS(I)        = ' '
       WTGAS(  I)        = 0.d0
       QBKGAS( I)        = 0.d0
 46   CONTINUE
!
      DO 47 I            = 1, MXGSAER
       CPREV( I)         = 0.d0
       CMODEL(I)         = 0.d0
 47   CONTINUE
!
      DO 48 I            = 1, MAXGL4  
       NKSURF(I)         = 0
       NCOATG(I)         = 0
 48   CONTINUE
!
      DO 49 I            = 1, IPHOT
       DEFPRAT(I,:)        = 0.d0
 49   CONTINUE
!
!
! *********************************************************************
!                            OPEN GLOBCHEM
!      MSL - 01-Feb-13 : Transferred from READER.F
! *********************************************************************
!
      !IF ( am_I_Root ) PRINT *,"READCHEM: KGLC=",KGLC
      !IF ( am_I_Root ) PRINT *,"READCHEM: IO93=",IO93
      ! Find a free file LUN
      IU_CHEMDAT = findFreeLUN()

      ! Name of file to open
      FILENAME   ='globchem.dat' 

#if defined( ESMF_ ) 
      ! We need to add the ".rc" extension to the "globchem.dat" file
      ! so that GEOS-5 will copy it to the proper run directory.
      FILENAME   = TRIM( FILENAME ) // '.rc'
#endif
      KGLC       = IU_CHEMDAT

      IF ( am_I_Root ) THEN 
         WRITE( 6, 301 ) TRIM( FILENAME ), KGLC
 301     FORMAT( 'GEOSCHEM::READER: Opening ', a, ' on unit ', i4  )
      ENDIF

      ! Open chemistry mechanism file
      OPEN( KGLC, FILE=TRIM( FILENAME ) )

      READ(KGLC,21)  HEADING 
      !WRITE(6,21) HEADING 
 13   READ(KGLC,21)  HEADING 
      IF (HEADING.NE.'BEGIN') GOTO 13 
 21   FORMAT(A76)
!
! *********************************************************************
!                   READ IN MASS BALANCE GROUPS USED 
! *********************************************************************
!
      NMBGAS = 9 
      READ(KGLC,59) (RINP(I), I = 1, NMBGAS)  
      READ(KGLC,61) (SINP(I), I = 1, NMBGAS)   
!
!#if defined( ESMF_ )
!
!      ALLOCATE(INTVAL(1))
!      INTVAL(1)  = NMASBAL
!      CALL GIGC_BCAST_CHAR( RINP,   14*NMBGAS,  RC )
!      CALL GIGC_BCAST_CHAR( NAMEMB, 14*IMASBAL, RC )
!      CALL GIGC_BCAST_INT( INTVAL, 1, RC )
!      NMASBAL = INTVAL(1)
!      INTVAL(1)  = NCSURBAN
!      CALL GIGC_BCAST_INT( INTVAL, 1, RC )
!      NCSURBAN = INTVAL(1)
!      DEALLOCATE(INTVAL)
!#endif

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
        IF ( am_I_Root ) WRITE(6,33) RINP(I) 
        CALL GEOS_CHEM_STOP
       ENDIF 
!
 36   CONTINUE
!
 59   FORMAT(20X,A3,8(1X,A3)) 
 61   FORMAT(21X,A1,8(3X,A1)) 
 33   FORMAT('READCHEM: MASS BALANCE GROUP ',A14,' NOT SET') 
!
! *********************************************************************
! * READ IN THE SPECIES AND OTHER DATA FOR THIS RUN FROM INPUT DATA   *
! *********************************************************************
!
!                    ITEMS IN THE FIRST READ STATEMENT
!                    ---------------------------------
!
! A/I/D 
!         D = SPECIES IS DEAD AND NOT USED
!         I = INACTIVE BUT USED (THESE SPECIES MUST ALSO BE INITIALIZED) 
!         A = SPECIES USED WHEN IFURBAN, IFTROP, OR IFSTRAT > 0
!             (URBAN, TROPSOSPHERIC AND STRATOSPHERIC SETS) 
!         U = SPECIES USED WHEN IFURBAN > 0
!         S = SPECIES USED WHEN IFSTRAT > 0
!         T = SPECIES USED WHEN IFTROP  > 0
!         R = SPECIES USED WHEN IFURBAN OR IFTROP  > 0
!         H = SPECIES USED WHEN IFTROP  OR IFSTRAT > 0 
! SPEC      = NAME OF THE SPECIES, 
! AB          TELLS WHETHER SPECIES ABSORBS RADIATION (THE SPECIES
!             DOES NOT NECESSARILY PHOTOLYZE)
! MW        = ATOMIC MASS IN AMU;  
! IFSTRAT   = 1: SOLVE STRATOSPHERIC CHEMISTRY
! IFTROP    = 1: SOLVE FREE TROPOSPHERIC CHEMISTRY
! IFURBAN   = 1: SOLVE URBAN CHEMISTRY
! INITCONC  = DEF'T BACKGROUND CONC. AT LOWEST LEVEL (VOL MIXING RATIO); 
! CSTRAT    = DEFAULT VOL MIX RATIO (FRACTION) IN STRATOSPHERE (25 KM) 
! CTROPL    = DEFAULT VOL MIX RATIO (FRACTION) IN FREE TROP OVER LAND (0 KM) 
! CTROPS    = DEFAULT VOL MIX RATIO (FRACTION) IN FREE TROP OVER SEA  (0 KM) 
! CURBAN    = DEFAULT VOL MIX RATIO (FRACTION) IN URBAN REGIONS (0 KM) 
!
! *********************************************************************
! ******************** READ IN SPECIES INFORMATION ********************
! *********************************************************************
!
!            FORMAT OF ITEMS IN THE SPECIES-LIST READ STATEMENT
!
! A/I/D SPEC   AB   MW  CSTRAT CTROPL CTROPS  CURBAN 
! A1,1X, A14,A2,1X,F6.0,E9.2,  E9.2,  E9.2    E9.2   
!
!        S     NO3  NH4     C    CL    BR    F     H     O    
! 20X,   I2,2X,I2,2X,I2,2X,I2,2X,I2,2X,I2,2X,I2,2X,I2,2X,I2 
!

      ! Read 1st line of species list
 10   READ(KGLC,11) (XINP(I), I=1,3),C1,CSTRAT,CTROPL,CTROPS,CURBAN
 11   FORMAT(A1,1X,A14,A2,1X,0PF6.2,4(1PE10.3))

      ! Test for "END" here (bmy, 4/7/03)
      IF (XINP(2).EQ.'END') GOTO 15

      ! Read 2nd line of species list 
      READ(KGLC,14) (ININT(I),I = 1, NMBGAS)
 14   FORMAT(20X,I2,8(2X,I2))
!
! *********************************************************************
! *  COUNT ACTIVE, INACTIVE, AND DEAD SPECIES. ALSO, SET UP ARRAYS    * 
! *  FOR OTHER INFORMATION.                                           *
! *********************************************************************

      IF (XINP(2).EQ.'END') GOTO 15
!
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
!
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
!
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
       IF ( am_I_Root ) WRITE(6,19) XINP(2), XINP(1)
       CALL GEOS_CHEM_STOP
      ENDIF
!
      NAMEGAS(JGAS)       = XINP(2) 
      WTGAS(  JGAS)       = C1
!
      IF (IFSTRAT.EQ.1.AND.IFTROP.EQ.0.AND.IFURBAN.EQ.0) THEN
       QBKGAS( JGAS)      = CSTRAT
      ELSEIF (IFSTRAT.EQ.0.AND.IFTROP.EQ.1.AND.IFURBAN.EQ.0) THEN
       QBKGAS( JGAS)      = CTROPL 
      ELSEIF (IFSTRAT.EQ.0.AND.IFTROP.EQ.0.AND.IFURBAN.EQ.1) THEN
       QBKGAS( JGAS)      = CURBAN
      ELSE
       QBKGAS( JGAS)      = CTROPL 
      ENDIF
!
      GOTO 10
!
! *********************************************************************
! *  SET NSPEC AS NUMBER OF ACTIVE SPECIES - 1 SINCE JUST INCREASED   *
! *  NACTIVE BEFORE THE 'END' STATEMENT. ALSO, CHECK SOME DIMENSIONS. *
! *********************************************************************
!
 15   CONTINUE
      NGAS         = NACTIVE
      NTSPECGAS    = NGAS + NINAC

      !=================================================================
      ! Chemical prod-loss diagnostic (bdf, bmy, 4/18/03)
      !=================================================================
      IF ( LFAMILY ) THEN

         ! Find species and rxns for ND65 diagnostic families
         CALL SETJFAM( NACTIVE, NINAC, am_I_Root )

         ! Reset quantities after SETJFAM 
         NSPEC(NCS)  = NACTIVE - 1
         NGAS        = NSPEC(NCS)
         NTSPECGAS   = NGAS + NINAC
         NTSPEC(NCS) = NGAS + NINAC
      ENDIF
!
      IF (NTSPECGAS.GT.IGAS   .OR.  NSDEAD.GT. NMDEAD) THEN
        IF ( am_I_Root ) WRITE(6,18) IGAS, NTSPECGAS, NMDEAD, NSDEAD   
        CALL GEOS_CHEM_STOP
      ENDIF
!
 18   FORMAT('READCHEM: ONE OF THE DIMENSIONS BELOW IS TOO SMALL:',/,
     1       'DIMENSION: IGAS     = ',I3,' VARIABLE: NTSPECGS = ',I3,/,  
     2       'DIMENSION: NMDEAD   = ',I3,' VARIABLE: NSDEAD   = ',I3)
!
! *********************************************************************
! *  RE-ARRANGE INACTIVE GAS ARRAYS SO THAT INFORMATION OF INACTIVE   *
! *  SPECIES APPEARS IMMEDIATELY AFTER INFORMATION OF ACTIVE SPECIES  * 
! *********************************************************************
!

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
!
! *********************************************************************
! * PRINT SPECIES INFORMATION IF IOSPEC = 1 + PRINT MASS BALANCE INFO * 
! *********************************************************************
!
      IF ( am_I_Root ) THEN
         IF (IOSPEC.EQ.1) THEN 

            ! Write species header
            WRITE( IO93, '(/,a)' ) REPEAT( '=', 79 )
            WRITE( IO93, 23      )
            WRITE( IO93, '(a,/)' ) REPEAT( '=', 79 )
            WRITE( IO93, 22      )
            
            ! Write species to "smv2.log"
            DO N = 1, NGAS
               WRITE(IO93,24) N, NAMEGAS(N), WTGAS(  N), QBKGAS(N)
            END DO
            IF (NINAC.GT.0)  WRITE(IO93,28)(NAMEGAS(NGAS+N),N=1,NINAC)
            IF (NSDEAD.GT.0) WRITE(IO93,31)(NAMD(N),N=1,NSDEAD)
         ENDIF
!
         PRINT *,' '
         DO 77 I = 1, NMASBAL
            IF (MBCOMP(I,MB1).GT.0) THEN 
               WRITE(6,98) NAMEMB(I) 
               WRITE(6,99)(WTMB( I,JMBCOMP(I,J,MB1),MB1),
     1              NAMEGAS(JMBCOMP(I,J,MB1)), J = 1, MBCOMP(I,MB1))
            ENDIF 
 77      CONTINUE
      ENDIF

 19   FORMAT('SPECIES ACTIVITY IS UNDEFINED ',A14,' TYPE = ',A2   )
 23   FORMAT('SPECIES FOR THIS RUN.  PHYSICAL CONSTS AND BOUNDARY', 
     1 ' CONDITIONS ALSO GIVEN.')
 22   FORMAT( 'NBR', 1X, 'NAME', 12X, 'MW', 1X, 'BKGAS(VMRAT)' )
 24   FORMAT(I3,1X,   A14,F7.2,1PE9.2)
 28   FORMAT(/'INACTIVE SPECIES FOR THIS RUN ARE:'//4(1X,A14))
 31   FORMAT(/'THE DEAD SPECIES FOR THIS RUN ARE:'//4(1X,A14)) 
 98   FORMAT('WEIGHTS AND SPECIES FOR MASS BALANCE EQUATION # ',A14)
 99   FORMAT(4(0PF5.1,1X,A14))  
!
! *********************************************************************
! *  SEARCH FOR SPECIFIC SPECIES NUMBERS USED IN OTHER SUBROUTINES    * 
! *********************************************************************
!
      ! Initialize for safety's sake (bmy, 7/7/03)
      IOXYGEN  = 0
      IH2O     = 0
      ICH4     = 0
      ILISOPOH = 0
      ILISOPNO3= 0 ! SOAupdate hotp 6/11/10
      ! add these as well (dkh, 10/06/06)  
      ILBRO2H  = 0
      ILBRO2N  = 0
      ILTRO2H  = 0
      ILTRO2N  = 0
      ILXRO2H  = 0
      ILXRO2N  = 0
      ! and for the new VSL index holders (jpp, 4/15/10)
      ICH2Br2  = 0
      ICHBr3   = 0
      ! SOAupdate: add NAP (hotp 6/25/09)
      ILNRO2H  = 0 
      ILNRO2N  = 0 

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
            CASE( 'LISOPNO3' )  ! SOAupdate: hotp 6/11/10
               ILISOPNO3= I
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
            CASE( 'CHBr3'   ) ! adding selection of CHBr3
               ICHBr3   = I
            CASE( 'CH2Br2'  ) ! adding selection of CH2Br2
               ICH2Br2  = I
            ! SOAupdate: add NAP (hotp 6/25/09)
            CASE( 'LNRO2H' )
               ILNRO2H   = I
            CASE( 'LNRO2N' )
               ILNRO2N   = I
            CASE DEFAULT
               ! Nothing
         END SELECT

      ENDDO
!
! *********************************************************************
! *****************    READ IN REACTION RATES   *********************** 
! *********************************************************************
!
! HERE, WE CAN HAVE 3 REACTANTS (EACH WITH COEFFICIENT OF 1) AND 12 
!             PRODUCTS [EACH WITH ANY REAL COEFFICIENT]. 
! FOR A 3-BODY M REACTION, PLACE AN M IN THE THIRD REACTANT POSITION 
!             [NO [+] BEFORE IT] 
! FOR A 3-BODY OTHER SPECIES REACTION, PLACE THE SPECIES NAME IN THE THIRD
!             REACTANT POSITION [NO [+] BEFORE IT] 
! FOR A 3RD REACTANT, PLACE THE SPECIES NAME IN THE THIRD REACTANT POSITION
!             WITH A PLUS BEFORE IT. 
! FOR A REACTANT NOT INCLUDED IN THE REACTION RATE [I.E.02] PLACE THE
!             SPECIES NAME IN THE FOURTH REACTANT POSITION [NO [+]  
!             BEFORE IT]. THE SPECIES MAY HAVE A COEFFICIENT PRECEDING IT.
! A PRODUCT MAY EITHER BE LISTED TWICE [OR MORE TIMES] OR 
!             HAVE A COEFFICIENT [I.E. 2 OR 3, 0.34] IMMEDIATELY BEFORE IT.
! N COLUMN:   NUMBER OF RATE COEFFICIENT LINES FOLLOWING TOP LINE
! P COLUMN:
!   P         = REACTION IS PRESSURE DEPENDENT 3-BODY REACTION. 
!               THE FIRST COEFFICIENT IS A 3-BODY COEF. THE SECOND IS 2-BODY.
!   S         = IDENTIFIES A ONE-BODY SURFACE REACTION 
!   E         = IDENTIFIES REVERSE EQUILIBRIUM REACTION
!   V         = IDENTIFIES CH3SCH3 + OH --> CH3S(OH)CH3
!   W         = IDENTIFIES O(1D) + N2 OR O2 
!   X         = IDENTIFIES OH  + HNO3 
!   Y         = IDENTIFIES OH  + CO   
!   Z         = IDENTIFIES HO2 + HO2  
!   G         = IDENTIFIES DMS + OH + O2
!   K         = IDENTIFIES WETDEP or HYDROLYSIS REACTIONS
!
! Fc COLUMN   = VALUE OF Fc FOR THREE-BODY RATE REACTIONS (SEE REF 9, P.1145)
! Fc(T)       = Fc CALCULATED AS EXP(-T(K)/Fc(T)) 
!
! REACTION RATES HAVE THE FORM K = A * (300 / T)**B * EXP(C / T)
!
! *********************************************************************
! *                     READ PRELIMINARY COMMENTS                     * 
! *********************************************************************
!
!                   ----- REACTION RATE FORMAT -----
!
! A/D ORD AR BR CR N P Fc FcT COM    X  +Y  +Z  IV  =aA  +bB  +cC  +dD  +...  
!
! A/D
!   D       = REACTION IS DEAD. SKIP THIS REACTION. 
!   A       = REACTION ACTIVE AND INCLUDED IN ALL CHEMISTRY SETS
!             (URBAN, TROPSOSPHERIC AND STRATOSPHERIC SETS) 
!   U       = REACTION IN URBAN CHEMISTRY SET  
!   S       = REACTION IN STRATOSPHERIC CHEMISTRY SET  
!   T       = REACTION IN TROPOSPHERIC CHEMISTRY SET  
!   R       = REACTION IN TROPOSPHERIC AND URBAN CHEMISTRY SETS  
!   H       = REACTION IN TROPOSPHERIC AND STRATOSPHERIC CHEMISTRY SETS  
! ORD       = ORDINAL NUMBER OF REACTION
! AR,BR,CR  = RATE COEFFICIENTS: K = AR * (300/T) * BR * EXP( CR / T)
! AR        = DEFAULT PHOTORATE (S-1) FOR PHOTOPROCESSES
! NCOF      = DESCRIBED IN 'N COLUMN' ABOVE
! P         = DESCRIBED IN 'P COLUMN' ABOVE
! FCVT       = CHARACTERIZES FALLOFF CURVE IN PRESSURE-DEPENDENT REACTION
! FCT1T,2T  = EXPONENTS GIVING TEMPERATURE DEPENDENCE OF FCVT
!             FCVT = EXP(-T/FCT1)     OR  
!             FCVT = EXP(-T/FCT1) + EXP(-FCT2/T)   
! COM       = A9 AT THE END IS CURRENTLY FOR COMMENTS.
! X,Y,Z     = REACTANTS
! Z         = REACTANT OR 3RD BODY 'M' OR OTHER THIRD BODY
! I         = COEFFICIENT (INTEGER) FOR V 
! V         = REACT NOT INCLUDED IN REACT. RATE, BUT WHICH IS LOST IN REACTION.
! a,b,c,d.  = COEFFICIENTS FOR PRODUCTS (1,2,3,0.45,1.32, ETC (>=0.))
! A,B,C,D.. = PRODUCTS
!
! *********************************************************************
! *                           READ REACTIONS                          * 
! *********************************************************************
!
 102  READ(KGLC,21)  HEADING 
      IF (HEADING.NE.'BEGIN') GOTO 102  
!
 310  READ(KGLC,330) DINP,IORD,ARRT(1),BRRT(1),KCRRT(1),NCOF,SPECL(1),
     1               FCVT(1),FCT1T(1),FCT2T(1),COMMENT

!FP_ISOP - for debug (6/2009)

!#if !defined( ESMF_ )
!     IF ( LPRT .and. am_I_Root ) THEN
!         WRITE(*,*) DINP,IORD,ARRT(1),BRRT(1),KCRRT(1),NCOF,SPECL(1),
!     &              FCVT(1),FCT1T(1),FCT2T(1),COMMENT
!      ENDIF
!#endif

      IF (NCOF+1.GT.MXCOF) THEN
       IF ( am_I_Root ) WRITE(6,155) NCOF+1, MXCOF, IORD
       CALL GEOS_CHEM_STOP
      ENDIF
!
      DO 350 I = 2, NCOF + 1
       READ(KGLC,330) JST,JORD,ARRT(I),BRRT(I),KCRRT(I),NDUM,SPECL(I),
     1                FCVT(I),FCT1T(I),FCT2T(I),COMMENT
 350  CONTINUE
!
      ! Now read 20 entries instead of 16 (bdf, bmy, 4/1/03)
      !FP_ISOP 24: now read 20 or 24 depending on input file (6/2009)
      ! NREAD is in comode.h (8/1/09)
      !READ(KGLC,332) (RINP(I),PINP(I),XINP(I),I=1,20)
      READ(KGLC,332) (RINP(I),PINP(I),XINP(I),I=1,NREAD)
!
 155  FORMAT('READCHEM: NCOF + 1 > MXCOF IN GLOBCHEM.DAT',3I4)
 330  FORMAT(A1,1X,I4,1X,ES8.2,1X,ES8.1,1X,I6,1X,I1,1X,A2,F6.2,1X,
     1       2(F6.0,1X),A20) 

      ! Increase format string to 14 products (bdf, 4/1/03)
 332  FORMAT(4(A1,0PF5.3,A14)/4(A1,0PF5.3,A14)/
     1     4(A1,0PF5.3,A14)/4(A1,0PF5.3,A14)/4(A1,0PF5.3,A14))
!
      IF (DINP.NE.'A'.AND.DINP.NE.'U'.AND.DINP.NE.'T'.AND.
     1    DINP.NE.'S'.AND.DINP.NE.'R'.AND.DINP.NE.'H') DINP = 'D'
!
      IF ((DINP.EQ.'U'.AND.IFURBAN.EQ.0).OR.
     1    (DINP.EQ.'T'.AND.IFTROP .EQ.0).OR.
     2    (DINP.EQ.'S'.AND.IFSTRAT.EQ.0).OR.
     3    (DINP.EQ.'R'.AND.IFURBAN.EQ.0.AND.IFTROP.EQ.0).OR.
     4    (DINP.EQ.'H'.AND.IFSTRAT.EQ.0.AND.IFTROP.EQ.0)) DINP = 'D'
!
      IF (XINP(1).EQ.'END KINETIC') THEN
!
       DO 323 NCS   = 1, NCSGAS 
        NRATES(NCS) = NTRATES(NCS)
 323   CONTINUE
!
       IDOPHOT     = 1
       GOTO 102 
      ELSEIF (XINP(1).EQ.'END PHOTOLYSIS') THEN
       GOTO 660 
      ELSEIF (DINP.EQ.'D') THEN
       GOTO 310 
      ENDIF

!#if defined( DEVEL ) || defined ( ESMF_ ) 
!      CALL GIGC_READCHEM_BCAST( C1,    CSTRAT, CTROPL, CTROPS, CURBAN, 
!    &                       ININT, IORD,   NCOF,   RC              )
!#endif
!
! *********************************************************************
! *      UPDATE REACTION NUMBER FOR REACTIONS THAT ARE USED           *
! *********************************************************************
! NRATCUR = CURRENT REACTION RATE NUMBER FOR A SET OF RATE COEFFICIENTS 
! NTRATES = CURRENT RATE COEFFICIENT NUMBER
! NALLRAT = COUNTS THE NUMBER OF ACTUAL REACTIONS
! SKIP URBAN         ('A', 'U', OR 'R')      REACTIONS IF NOT USED
! SKIP TROPOSPHERIC  ('A', 'T', 'R', OR 'H') REACTIONS IF NOT USED
! SKIP STRATOSPHERIC ('A', 'S', OR 'H')      REACTIONS IF NOT USED
!
      DO 325 NCS      = 1, NCSGAS
!
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
!
       IF (NOUSE(NCS).EQ.0) THEN
        NK               = NTRATES(NCS) +  1
        NRATCUR(NCS)     = NK
        NALLRAT(NCS)     = NALLRAT(NCS) + 1
        NAR              = NALLRAT(NCS)
        NCEQUAT(NAR,NCS) = NK
!
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

! *********************************************************************
!       SET UP A DEFAULT PHOTORATE (SEC-1), STORE ORDINAL NUMBER
! *********************************************************************
!
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
!
! *********************************************************************
! *  CHECK WHETHER EACH SPECIES SPOT IN INPUT REACTION SET IS ACTIVE, *
! *  INACTIVE, BLANK, DEAD, OR 'M'. STOP IF THE SPECIES IS NONE       *
! *  JNUM = -J = NON 'M' THIRD BODY IN PRESSURE-DEPENDENT REACTION.   *
! *********************************************************************
!
      DO 360 I = 1, NPRODHI
       IF (XINP(I).NE.' ') THEN
!
        IF (I.LE.NMREAC.AND.PINP(I).NE.0.) THEN
         IF ( am_I_Root ) WRITE(6,450) IORD
         CALL GEOS_CHEM_STOP
        ENDIF
!
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
!
        IF (I.GT.NMREAC) THEN
         DO 390 J = 1, NSDEAD 
          IF (XINP(I).EQ.NAMD(J)) GOTO 360 
 390     CONTINUE
        ENDIF
!
        IF ( am_I_Root ) WRITE(6,400) IORD, XINP(I)
        CALL GEOS_CHEM_STOP
!
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
!
        IF (IDOPHOT.EQ.1) NAMEPHOT(I,NPHOTALL) = XINP(I)
       ENDIF
 360  CONTINUE
!
      DO 415 NCS           = 1, NCSGAS 
       IF (NOUSE(NCS).EQ.0) THEN 
        NK                 = NRATCUR(NCS) 
        IF (IRM(1,NK,NCS).EQ.0.OR.(IRM(3,NK,NCS).GT.0.AND.
     1      IRM(2,NK,NCS).EQ.0)) THEN
         IF ( am_I_Root ) WRITE(6,430) IORD
         CALL GEOS_CHEM_STOP
        ENDIF
       ENDIF
 415  CONTINUE
!
! *********************************************************************
! *             PLACE SPECIAL-RATE INFORMATION INTO ARRAYS            *
! *********************************************************************
!
! NMAIR   = # REACTIONS WHERE THE SPECIES IN THE THIRD POSITION
!           IS 'M' = 'O2 + N2': IRM(3,NK) = -9999 
! NMO2    = # REACTIONS WHERE THE SPECIES IN THE THIRD POSITION
!           IS O2:  IRM(3,NK) = -9998  
! NMN2    = # REACTIONS WHERE THE SPECIES IN THE THIRD POSITION
!           IS N2:  IRM(3,NK) = -9997  
! NM3BOD  = # REACTIONS WHERE THE SPECIES IN THE THIRD POSITION
!           IS ANY OTHER SPECIES:  IRM(3,NK) = -SPECIES NUMBER
! NNEQ    = # THERMALLY DISSOCIATING EQUILIBRIUM REACTIONS. PREVIOUS
!             EQUATION MUST BE PRESSURE-DEPENDENT. (SPECIAL = 'E')
!             THUS: REQUIRES 3 REACTIONS TOTAL (2 FOR PRESS. DEP, 1 EQ.)
! NPRESM  = # PRESSURE DEPENDENT 3-BODY REACTIONS  (SPECIAL = 'P')
!

! bdf smv2 count number of emission and drydep reactions. Used in calcrate
!  to put emissions into the chemistry, and drydep out of chemistry.


      DO NCS=1,NCSGAS
         IF (NOUSE(NCS) .EQ. 0 ) THEN
            NK                 = NRATCUR(NCS)
            IF (XINP(1).EQ.'EMISSION') THEN
               NEMIS(NCS)           = NEMIS(NCS) + 1
               IF (NEMIS(NCS)   .GT. MAXGL3) THEN
                  IF ( am_I_Root ) THEN
                     WRITE(*,*) 'ERROR NEMIS   Greater then MAXGL3',
     x                ' NEMIS(NCS)   = ',NEMIS(NCS), 'MAXGL3 = ',MAXGL3
                     WRITE( 6, '(a)' ) 'STOP 124'
                  ENDIF
                  CALL GEOS_CHEM_STOP
               ENDIF
               NKEMIS(NEMIS(NCS),NCS)   = NK
            ENDIF
         ENDIF
!
!  NDRYDEP  = number of dry deposition reactions read in
!  NKDRY    = reaction numbers of dry deposition reactions
         IF (NOUSE(NCS) .EQ. 0) THEN
            NK                 = NRATCUR(NCS)
            IF (XINP(NPRODLO).EQ.'DRYDEP') THEN
               NDRYDEP(NCS)           = NDRYDEP(NCS) + 1
               IF (NDRYDEP(NCS) .GT. MAXDEP) THEN
                  WRITE(*,*) 'ERROR NDRYDEP Greater then MAXDEP',
     x                 ' NDRYDEP(NCS)=',NDRYDEP(NCS),'MAXDEP=',MAXDEP
                  IF ( am_I_Root ) WRITE( 6, '(a)' ) 'STOP 125'
                  CALL GEOS_CHEM_STOP
               ENDIF
               NKDRY(NDRYDEP(NCS),NCS)   = NK
            ENDIF
         ENDIF

! bdf smv2 use Q to flag O3 photolysis, code is not confused by 'A''s
         IF (NOUSE(NCS) .EQ. 0) THEN
            NK                 = NRATCUR(NCS)
            IF (SPECL(1).EQ.'Q') NKO3PHOT(NCS)=NK !Flag O3 photolysis
            IF (SPECL(1).EQ.'QB') NKO3PHOTH2(NCS)=NK !Flag O3 photolysis (bhh, jmao, eam, 7/18/11)
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
            IF ( am_I_Root ) WRITE(6,440) IORD, SPECL(1), NCOF
            CALL GEOS_CHEM_STOP
         ENDIF 
!
         DO 420 NCS           = 1, NCSGAS 
            IF (NOUSE(NCS).EQ.0) THEN 
               NK                 = NRATCUR(NCS) 
!
               ITHIRDB            = IRM(3,NK,NCS)
!
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
!
               IF (SPECL(1).EQ.'P') THEN
                  NPRESM(NCS)       = NPRESM(NCS) + 1  
                  NR                = NPRESM(NCS) 
                  NREACPM(NR,NCS)   = NK
               ELSEIF (SPECL(1).EQ.'E') THEN
                  NNEQ(NCS)         = NNEQ(NCS)   + 1
                  NN                = NNEQ(NCS)
                  NREACEQ(NN,NCS)   = NK
!
! EQUILIBRIUM REACTIONS USE THE PREVIOUS REACTION AS PART OF THE 
! RATE CALCULATION (SEE CALCRATE.F). THE PREVIOUS REACTION MAY BE
! PRESSURE DEPENDENT.
!
                  NREQOTH(NN,NCS)   = NCEQUAT(NALLRAT(NCS)-1,NCS)
               ENDIF
!
! NKSPECV = SPECIAL REACTION CH3SCH3 + OH = CH3S(OH)CH3 (SPECL = 'V') 
! NKSPECW = SPECIAL REACTION O(1D)   + O2,N2            (SPECL = 'W') 
! NKSPECX = SPECIAL REACTION OH      + HNO3             (SPECL = 'X')
! NKSPECY = SPECIAL REACTION OH      + CO               (SPECL = 'Y')
! NKSPECZ = SPECIAL REACTION HO2     + HO2              (SPECL = 'Z') 
!
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

                  ! Same for BrNO3 hydrolysis rxn (jpp, 5/4/10)
                  IF ( XINP(1) == 'BrNO3' ) THEN
                     NKBrNO3 = NK
                  ENDIF

                  ! Same for HOBr psuedo-rxn (jpp, 03/22/11)
                  !   ( HOBr + HBr ---het---> Br2 + H2O )
                  IF ( (XINP(1) == 'HOBr') .and. 
     &                 (trim(adjustl(comment)) == 'h1') ) THEN
                     NKHOBr = NK
                  ENDIF

                  ! Same for 1st HBr psuedo-rxn (jpp, 03/22/11)
                  !   ( HOBr + HBr ---het---> Br2 + H2O )
                  IF ( XINP(1) == 'HBr' .and. 
     &                 (trim(adjustl(comment)) == 'h1') ) THEN
                     NK1HBr = NK
                  ENDIF

                  ! Now the HOBr ice psuedo-rxn (jpp, 06/16/11)
                  !   ( HOBr + HBr ---ice---> Br2 + H2O )
                  IF ( (XINP(1) == 'HOBr') .and. 
     &                 (trim(adjustl(comment)) == 'h2') ) THEN
                     NK2HOBr = NK
                  ENDIF

                  ! Now the HOBr ice psuedo-rxn (jpp, 06/16/11)
                  !   ( HOBr + HBr ---ice---> Br2 + H2O )
                  IF ( (XINP(1) == 'HBr') .and. 
     &                 (trim(adjustl(comment)) == 'h2') ) THEN
                     NK2HBr = NK
                  ENDIF

               ENDIF
!

               IF (SPECL(1).EQ.'HR')  THEN    !modification of the rate for HO2+RO2
                  NNRO2HO2(NCS)               = NNRO2HO2(NCS) + 1
                  NKSPECRO2HO2( NNRO2HO2(NCS),NCS) = NK
               ENDIF

               ! DA, DB, and DC are now in use by FP
               ! This part applies only to 'D' (hotp 8/1/09)
               IF (SPECL(1).EQ.'D' .AND. SPECL(2).NE.'A' .AND.
     &             SPECL(2).NE.'B' .AND. SPECL(2).NE.'C')  THEN
                  NNADDD(NCS)               = NNADDD(NCS) + 1
                  NKSPECD( NNADDD(NCS),NCS) = NK
               ENDIF
!
               IF (SPECL(1).EQ.'A')  THEN
                  NNADDA(NCS)               = NNADDA(NCS) + 1
                  NKSPECA( NNADDA(NCS),NCS) = NK
               ENDIF
!
               IF (SPECL(1).EQ.'B')  THEN
                   NNADDB(NCS)              = NNADDB(NCS) + 1
                  NKSPECB( NNADDB(NCS),NCS) = NK
               ENDIF
!
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
!
! *********************************************************************
! *                       SURFACE REACTIONS                           * 
! *********************************************************************
! ARR(INIT)         = REACTION PROBABILITY 
! ARR(FINAL)        = REACTION PROBABILITY * QTHERMG  
! QTHERMG * SQRT(T) = (1/4) * THERMAL VELOCITY OF GAS (CM S-1)
!
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
!
                  IF (JGAS3.NE.0) THEN
                     IF ( am_I_Root ) WRITE(6,470) NK 
                     CALL GEOS_CHEM_STOP
                  ENDIF  
               ENDIF 
!
! *********************************************************************
! *       SET ARRAYS FOR CALCULATING REACTION RATES EFFICIENTLY       * 
! *********************************************************************
! NARR = NUMBER OF REACTIONS OF THE FORM K = A
! NABR = NUMBER OF REACTIONS OF THE FORM K = A * (300 / T)**B
! NACR = NUMBER OF REACTIONS OF THE FORM K = A * EXP(C / T)
! NABC = NUMBER OF REACTIONS OF THE FORM K = A * (300 / T)**B * EXP(C / T)
! NKARR, NKBRR, NKACR, NKABC = REACTION RATE NUMBERS OF EACH
! NARR,  NABR,  NACR,  NABC  REACTION 
!


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
!
            ENDIF
!       ENDIF NOUSE(NCS).EQ.0
!
 420     CONTINUE
      ENDIF
!     ENDIF IDOPHOT.EQ.0
!
      GOTO 310 
!
 400  FORMAT('INVALID REACT',I4,' W UNRECOGNIZABLE OR DEAD SPEC ',A14,
     1       'ALL REACTANTS MUST BE ACTIVE/INACTIVE. PRODS CAN BE DEAD') 
 430  FORMAT('READCHEM:REACT ',I3,' 1ST SPOT EMPTY OR 3RD SPOT FILLED ', 
     1       ' BUT 2ND EMPTY')
 440  FORMAT('READCHEM: REACT ',I3,'OR BEFORE: SPECIAL REACTION WITH ',
     1       'DELIMETER ',A2,' HAD INCORRECT # OF REACTIONS ',I5)
 450  FORMAT('READCHEM: ORD# REACT ',I3,' CANT HAVE COEFF > 1')
 470  FORMAT('READCHEM: SURFACE REACTION ',I5,'HAS THREE REACTANTS ')
 510  FORMAT(I3,1X,ES7.1,1X,F7.2,1X,I6,1X,0PF3.2,1X,
     1       A6,2(A1,A6),14(A1,0PF3.1,A6)) 
 520  FORMAT( 'KINETIC REACTIONS FOR ', A,' CHEMISTRY',/,   
     1        'RATE CONSTANTS HAVE FORM K = A * (300/T)**B * EXP(C/T).')
 521  FORMAT( 'NMBR   A       B     C     Fv       REACTION' )
 525  FORMAT( 'PHOTOPROCESS REACTIONS FOR ', A,' CHEMISTRY' )   
 526  FORMAT( 'NMBR   DEFP (S-1)                        REACTION' )
!
! *********************************************************************
! *                  PRINT OUT REACTION INFORMATION                   *
! *********************************************************************
!
 660  IF (IOREAC.EQ.1) THEN 
       DO 502 NCS    = 1, NCSGAS 

        ! Write reaction header
	IF ( am_I_Root ) THEN
           WRITE( IO93, '(/,a)' ) REPEAT( '=', 79 )
           WRITE( IO93, 520     ) TRIM( CHEMTYP(NCS) )
           WRITE( IO93, '(a,/)' ) REPEAT( '=', 79 )
           WRITE( IO93, 521     )
        ENDIF

        DO 500 NK    = 1, NTRATES(NCS) 

         ! Write photo rxn header
         IF ( NK .EQ. NRATES(NCS)+1 .AND. am_I_Root) THEN
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
!
         RINP(5)    = '='
	 IF ( am_I_Root ) THEN
            WRITE(IO93,510) NK,ARR(NK,NCS),BRR(NK,NCS),KCRR(NK,NCS),
     1           FCV(NK,NCS),XINP(1),'+',XINP(2),
     2           '+',XINP(3),(RINP(I),PINP(I),XINP(I),
     3           I = 5,NPRODUC(NK,NCS))
         ENDIF
 500    CONTINUE 
 502   CONTINUE 
      ENDIF

!
! *********************************************************************
! ***********************  CHECK SOME DIMENSIONS  *********************
! *********************************************************************
!
      DO 670 NCS = 1, NCSGAS  
       IF (NTRATES(NCS) .GT. NMTRATE  .OR. NPHOTALL   .GT. IPHOT   .OR. 
     1     NTSPECGAS    .GT. IGAS     .OR. NSDEAD     .GT. NMDEAD  .OR.
     2     NPRODHI      .GT. NMRPROD .OR.
     3     NMAIR(NCS)   .GT. MAXGL3   .OR. NMO2(NCS)  .GT. MAXGL3  .OR.
     4     NMN2(NCS)    .GT. MAXGL2   .OR. NPRESM(NCS).GT. MAXGL2  .OR.
     5     NSURFACE(NCS).GT. MAXGL4   .OR. NM3BOD(NCS).GT. MAXGL3) THEN
!
        IF ( am_I_Root ) THEN
           WRITE(6,680)
     &          NMTRATE,NTRATES(NCS), IPHOT  , NPHOTALL,     IGAS,
     &          NTSPECGAS,NMDEAD ,NSDEAD,       NMRPROD, NPRODHI,
     &          MAXGL3 ,NMAIR(NCS),   MAXGL3 , NMO2(NCS),    
     &          MAXGL2,NMN2(NCS), MAXGL2 ,NPRESM(NCS),  MAXGL4 , 
     &          NSURFACE(NCS),MAXGL3,NM3BOD(NCS)  
        ENDIF
        CALL GEOS_CHEM_STOP
       ENDIF
 670  CONTINUE 
!
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
!
! *********************************************************************
! ************************** SET KEY PARAMETERS ***********************
! *********************************************************************
!
      DO 702 NCS          = 1, NCSGAS 
       NSPEC(NCS)         = NGAS
       NTSPEC(NCS)        = NTSPECGAS
!
       DO 700 JGAS        = 1, NTSPECGAS
        NAMENCS(JGAS,NCS) = NAMEGAS(JGAS)  
        QBKCHEM(JGAS,NCS) = QBKGAS( JGAS)
 700   CONTINUE 
 702  CONTINUE 

!#if defined( ESMF_ )
!      CALL GIGC_BCAST_INT(NSPEC, ICS, RC)
!      CALL GIGC_BCAST_INT(NTSPEC, ICS, RC)
!!      CALL GIGC_BCAST_REAL8(NTSPEC, ICS, RC)
!#endif
!---smv2-s
! Update (gcc)
! bdf smv2, put this in here for now.
! *********************************************************************
! ****************   READ INFO FOR AEROSOL REACTIONS   ****************
! *********************************************************************

! astkcf  -- sticking coefficient (no unit), order of 0.1
! xgdfcf  -- gas phase diffusion coefficient (cm2/s), order of 0.1
! iarsfa  -- fortran unit number for reading sulfate abundance file
! mwarsl  -- aerosol molecular wright (g/mol)    [H2SO4=98]
! ruarsl  -- density of aerosol (g/cc)
! RH100   -- deliquescence point, relative humidity below which we
!            have no wet aerosols

      ! Find a free file LUN
      I        = findFreeLUN()

      ! Name of file to open
      FILENAME = 'chemga.dat'

#if defined( ESMF_)
      ! We need to add the ".rc" extension to the "chemga.dat" file
      ! so that GEOS-5 will copy it to the proper run directory.
      FILENAME = TRIM( FILENAME ) // '.rc'
#endif

      IF ( am_I_Root ) THEN
         WRITE( 6, 300 ) TRIM( FILENAME ), I
 300     FORMAT( 'GEOSCHEM::READCHEM: Opening ', a, ' on unit ', i4 )
      ENDIF

      OPEN(I, FILE=TRIM( FILENAME ), FORM='FORMATTED', STATUS='OLD' )
      READ(I,*)
      READ(I,610) ASTKCF
      READ(I,*)
      READ(I,610)
      READ(I,*)
      READ(I,620) MWARSL
      READ(I,610) RUARSL
      READ(I,630) RH100
      READ(I,620) IARSFA
      CLOSE(I)
 610  FORMAT(E10.3)
 620  FORMAT(I10)
 630  FORMAT(F10.2)
!
! *********************************************************************
! ***** CALL JSPARSE TO SET ARRAYS FOR SOLVING CHEMICAL EQUATIONS ***** 
! *********************************************************************
!
      ! Call SETPL to setup ND65 prod/loss diagnostic
      ! SETPL must be called before JSPARSE (ljm, bmy, 5/9/03)
      IF ( LFAMILY ) CALL SETPL( am_I_Root, Input_Opt, RC )

      ! IFSOLVE = 1 means we are calling the chemistry solver
      IF ( IFSOLVE .EQ. 1 ) THEN

         ! Loop over chemistry regimes (for now NCSGAS=NCSURBAN=1)
         DO NCS = 1, NCSGAS

            ! Set up sparse matrix stuff
            CALL JSPARSE( am_I_Root )

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
                  ENDIF

                  ! Avoid out-of-bounds-errors by making sure that each
                  ! IL* index is nonzero before doing the IF test for J. 
                  ! (bmy, 6/13/11)
                  IF ( ILBRO2H > 0 ) THEN 
                     IF ( J == MAPPL(ILBRO2H,NCS) ) THEN
                        ITS_NOT_A_ND65_FAMILY(J) = .FALSE.
                        EXIT
                     ENDIF
                  ENDIF

                  IF ( ILBRO2N > 0 ) THEN 
                     IF ( J == MAPPL(ILBRO2N,NCS) ) THEN
                        ITS_NOT_A_ND65_FAMILY(J) = .FALSE.
                        EXIT
                     ENDIF 
                  ENDIF

                  IF ( ILTRO2H > 0 ) THEN 
                     IF ( J == MAPPL(ILTRO2H,NCS) ) THEN
                        ITS_NOT_A_ND65_FAMILY(J) = .FALSE.
                        EXIT
                     ENDIF 
                  ENDIF

                  IF ( ILTRO2N > 0 ) THEN 
                     IF ( J == MAPPL(ILTRO2N,NCS) ) THEN
                        ITS_NOT_A_ND65_FAMILY(J) = .FALSE.
                        EXIT
                     ENDIF 
                  ENDIF

                  IF ( ILXRO2H > 0 ) THEN 
                     IF ( J == MAPPL(ILXRO2H,NCS)  ) THEN
                        ITS_NOT_A_ND65_FAMILY(J) = .FALSE.
                        EXIT
                     ENDIF 
                  ENDIF

                  IF ( ILXRO2N > 0 ) THEN 
                     IF ( J == MAPPL(ILXRO2N,NCS) ) THEN
                        ITS_NOT_A_ND65_FAMILY(J) = .FALSE.
                        EXIT
                     ENDIF
                  ENDIF

                  ! SOAupdate: add NAP (hotp 6/25/09)
                  IF ( ILNRO2H > 0 ) THEN 
                     IF ( J == MAPPL(ILNRO2H,NCS) ) THEN
                        ITS_NOT_A_ND65_FAMILY(J) = .FALSE.
                        EXIT
                     ENDIF
                  ENDIF

                  IF ( ILNRO2N > 0 ) THEN 
                     IF ( J == MAPPL(ILNRO2N,NCS) ) THEN
                        ITS_NOT_A_ND65_FAMILY(J) = .FALSE.
                        EXIT
                     ENDIF
                  ENDIF

               ENDDO  
            ENDDO
         ENDDO

      ELSE

         ! If we are not calling the chemistry solver, then
         ! set number of active gas photoprocesses to zero
         NPHOTALL = 0

      ENDIF
!
! *********************************************************************
! *******************   END OF SUBROUTINE READCHEM   ****************** 
! *********************************************************************
!
      RETURN
      END SUBROUTINE READCHEM
