! $Id: calcrate.f,v 1.10 2005/05/09 14:33:56 bmy Exp $
      SUBROUTINE CALCRATE( SUNCOS )
!
!******************************************************************************
!  Subroutine CALCRATE computes reaction rates before passing them to the
!  SMVGEAR solver.  (M. Jacobson 1997; gcc, bdf, bmy, 4/1/03, 2/17/05)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) SUNCOS (REAL*8) : Array of COSINE( solar zenith angle )
!  
!  NOTES:
!  (1 ) For GEOS-CHEM we had to remove several arrays from "comode.h" and 
!        declare these allocatable in "comode_mod.f".  This allows us to only 
!        allocate these if we are doing a fullchem run.  Now also references
!        routines from "diag_mod.f", "drydep_mod.f", "error_mod", and 
!        "planeflight_mod.f".  Also, CMN_SAV has now been eliminated. 
!        Also modified ND22 FAST-J diagnostics accordingly for SMVGEAR II. 
!        Now added special rxn for DMS+OH+O2.  Force double precision with
!        "D" exponents. (gcc, bdf, bmy, 4/1/03)
!  (2 ) Now implement interannually-varying CH4 field.  Now reference GET_YMID
!        from "grid_mod.f".  Now reference AIRDENS array from "comode_mod.f". 
!        Added YLAT variable for grid-box latitude.  Cosmetic changes.
!        (bnd, bmy, 7/1/03)
!  (3 ) Comment out AREAXT, this is not needed.  Also comment out sections 
!        which compute surface rxns and 3-body rxns, since these are not
!        applicable to GEOS-CHEM.  Declare ABSHUMK as a local variables since 
!        it is only ever used w/in "smvgear.f".  Remove obsolete variables 
!        from documentation.  Now call ARCHIVE_RXNS_FOR_PF to save rxn rates
!        for the ND40 planeflight diagnostic before exiting. (bmy, 7/16/03)
!  (4 ) Now apply dry deposition throughout the entire PBL, in order to prevent
!        short-lived species such as HNO3 from being depleted too much in
!        the shallow GEOS-3 surface layer.  Now reference PBLFRAC from
!        "drydep_mod.f".  Now declare DENAIR, CONCO2, CONCN2, T3I, TEMP1, T3K
!        and PRESSK as local variables, since these are only used w/in 
!        this routine and nowhere else -- also remove these from /DKBLOOP/ in
!        "comode.h".  (rjp, bmy, 7/30/03)
!  (5 ) Now references GEOS_CHEM_STOP from "error_mod.f".  Added internal
!        function N2O5 to compute the GAMMA "stickiness" parameter for N2O5
!        hydrolysis, which is a function of aerosol type.  Now also pass N2O5
!        reaction rate to ARCHIVE_RXNS_FOR_PF. (bmy, 8/8/03)
!  (6 ) Updated loss rate for O(1D) with H2O according to new rate measurement
!        from JPL (mje, bmy, 5/26/04)
!  (7 ) Now use GET_FRAC_UNDER_PBLTOP from "pbl_mix_mod.f" instead of
!        PBLFRAC from "drydep_mod.f" (bmy, 2/17/05)
!******************************************************************************
!
      ! References to F90 modules 
      USE COMODE_MOD,      ONLY : ABSHUM, AIRDENS, ERADIUS, IXSAVE, 
     &                            IYSAVE, IZSAVE,  JLOP,    PRESS3,  
     &                            REMIS,  T3,      TAREA
      USE DIAG_MOD,        ONLY : AD22,   LTJV
      USE DRYDEP_MOD,      ONLY : DEPSAV
      USE ERROR_MOD,       ONLY : ERROR_STOP, GEOS_CHEM_STOP
      USE GRID_MOD,        ONLY : GET_YMID
      USE PBL_MIX_MOD,     ONLY : GET_FRAC_UNDER_PBLTOP
      USE PLANEFLIGHT_MOD, ONLY : ARCHIVE_RXNS_FOR_PF

      IMPLICIT NONE

#     include "CMN_SIZE"  ! Size parameters
#     include "CMN"       ! STT, etc
#     include "comode.h"  ! SMVGEAR II arrays
#     include "CMN_DIAG"  ! ND22, LD22
      
      ! Local variables
      INTEGER :: KLOOP,JLOOP,I,NK,JOLD2,JOLD,NK1,NKN,NH,MLOOP,J
      INTEGER :: NP,K,IFNC,IBRCH,IX,IY,IZ,IJWINDOW,INDEX,NN,ii,jj
      INTEGER :: PHOTVAL,KSUN

      REAL*8  :: ARRNK,FCVNK,FCT1,FCT2,FCT,XYRAT,BLOG,FEXP,RATE3M
      REAL*8  :: CONSTC,RATE2AIR,RATE3H2O,RIS,RST,TBEGIN,TFINISH
      REAL*8  :: PBEG,PFIN,PBEGNEW,PFINNEW,TOFDAYB,TOFDAYE,HOURANGB
      REAL*8  :: HOURANGE,SINFUNCB,SINFUNCE,XAREA,XRADIUS,XSQM,RRATE2
      REAL*8  :: XSTKCF,GMU,SUNCOS(MAXIJ),DUMMY(KBLOOP),XDENA,XSTK
      REAL*8  :: TK,     CONSEXP, VPRESH2O

      ! External functions
      REAL*8, EXTERNAL :: RTFUNC, FYRNO3,  ARSL1K,  FJFUNC

      CHARACTER*4      :: SPECNAME

      ! Added for heterogeneous chemistry (bmy, 11/15/01, 8/7/03)
      LOGICAL          :: HETCHEM
      INTEGER          :: N
      REAL*8           :: SUMAREA, TOTAREA, DUMMY2(KBLOOP)

      ! For grid-box latitude (bnd, bmy, 7/1/03)
      REAL*8           :: YLAT

      ! Variables from "comode.h" which are only ever used in "calcrate.f"
      ! Remove them from "comode.h" and the THREADPRIVATE declarations
      ! (bmy, 7/28/03) 
      REAL*8           :: ABSHUMK(KBLOOP), DENAIR(KBLOOP)
      REAL*8           :: CONCO2(KBLOOP),  CONCN2(KBLOOP)
      REAL*8           :: T3I(KBLOOP),     TEMP1(KBLOOP)
      REAL*8           :: T3K(KBLOOP),     PRESSK(KBLOOP) 

#if   defined( LSLOWJ )
      ! Include SLOW-J header file if FAST-J is turned off (bmy, 9/30/99)
#     include "comsol.h"
#endif

      ! FAST-J: Zero out the dummy array (bmy, 9/30/99)
      DUMMY = 0d0
C
C *********************************************************************
C ************        WRITTEN BY MARK JACOBSON (1993)      ************
C ***             (C) COPYRIGHT, 1993 BY MARK Z. JACOBSON           *** 
C ***       U.S. COPYRIGHT OFFICE REGISTRATION NO. TXu 670-279      *** 
C ***                         (650) 723-6836                        *** 
C *********************************************************************
C
C CCCCCCC     A     L       CCCCCCC  RRRRRRR     A     TTTTTTT  EEEEEEE
C C          A A    L       C        R     R    A A       T     E 
C C         A   A   L       C        RRRRRRR   A   A      T     EEEEEEE
C C        AAAAAAA  L       C        R  R     AAAAAAA     T     E 
C CCCCCCC A       A LLLLLLL CCCCCCC  R    R  A       A    T     EEEEEEE
C
C *********************************************************************
C * THIS SUBROUTINE CALCULATES KINETIC REACTION AND PHOTORATES        *
C * (S-1, CM3 S-1, OR CM6 S-2) AND PRESSURE AND TEMPERATURE-          *
C * DEPENDENCE FOR GAS-PHASE CHEMICAL REACTIONS.                      *
C *                                                                   *
C * HOW TO CALL SUBROUTINE:                                           *
C * ----------------------                                            *
C *  CALL CALCRATE.F FROM PHYSPROC.F WITH                             * 
C *     NCS  = 1..NCSGAS FOR GAS CHEMISTRY                            *
C *********************************************************************
C
C *********************************************************************
C *********************  GAS-PHASE CHEMISTRY  *************************
C *********************************************************************
C
      IF (NCS.LE.NCSGAS) THEN 
C
C *********************************************************************
C *            CALCULATE CONCENTRATIONS OF FIXED SPECIES              *
C *********************************************************************
C AERSURF  = PARTICLE SURFACE AREA (CM2 CM-3)
C KTLOOP   = NUMBER OF GRID-CELLS IN A GRID-BLOCK
C AM       = MOLECULAR WEIGHT OF AIR (28.966 G MOLE-1)
C AVG      = AVOGADRO'S NUMBER (6.02252E+23  # MOLE-1)
C BOLTG    = BOLTZMANN'S CONSTANT (1.38054E-16 ERG K-1) 
C RHO3     = DENSITY OF AIR         (G CM-3) = WTAIR*P(DYN CM-2)/(RSTARG*T) 
C RSTARG   = BOLTG * AVG
C DENAIR   = DENSITY OF AIR         (# CM-3) = P(DYNCM-2) / (T * BOLTG) 
C CONCO2   = OXYGEN CONCENTRATION   (# CM-3)
C CONCN2   = NITROGEN CONCENTRATION (# CM-3)
C RRATE    = RATE CONST (EITHER S-1, CM**3-AIR #-1 S-1, CM**6 #-2 S-1,
C            OR CM**9 #-3 S-1.
C PRESS3   = AIR PRESSURE AT VERTICAL CENTER OF LAYER (MB)
C
C --------------------------- AIR, O2, N2 -----------------------------
C

!   be sure to check out these values!!!
         KSUN=0
         DO 20 KLOOP        = 1, KTLOOP
            JLOOP             = LREORDER(JLOOPLO+KLOOP)

            ! Add DENAIR here instead of in physproc.f, so that we
            ! can eliminate the /DKBLOOP/ common block (bmy, 7/28/03)
            DENAIR(KLOOP)     = AIRDENS(JLOOP)
     
            PRESSK(KLOOP)     = PRESS3(JLOOP)
            T3K(KLOOP)        = T3(JLOOP)
            T3I(KLOOP)        = 1.d0/T3(JLOOP)
            ABSHUMK(KLOOP)    = ABSHUM(JLOOP)
            TEMP1(KLOOP)      = 300.d0    / T3K(KLOOP)
            CONCO2(KLOOP)     = 0.2095d0  * DENAIR(KLOOP)
            CONCN2(KLOOP)     = 0.7808d0  * DENAIR(KLOOP)
C
C   Check if sun is up anywhere in this block of grid-boxes.
C   IFSUN gets used in CALCRATE
C   Get the right index for SUNCOS, which is calculated
C   outside of chemistry module.
C   (This works for LEMBED= .TRUE. or .FALSE.)
C
            IX                = IXSAVE(JLOOP)
            IY                = IYSAVE(JLOOP)
            IJWINDOW          = (IY-1)*IIPAR + IX
            GMU               = SUNCOS(IJWINDOW)
            IF(GMU .GT. 0.D0) KSUN = 1
 20      CONTINUE
         IFSUN = 2-KSUN
C
C ---------------------------     H2O     -----------------------------
C
      IF (IH2O.GT.NGAS) THEN
         DO KLOOP      = 1, KTLOOP 
            TK               = T3K(KLOOP) 
            CONSEXP          = 17.2693882D0 * (TK - 273.16D0) / 
     x           (TK - 35.86D0) 
            VPRESH2O         = CONSVAP * EXP(CONSEXP) * T3I(KLOOP) 
            CBLK(KLOOP,IH2O) = ABSHUMK(KLOOP) * VPRESH2O  
         ENDDO
      END IF
C
C -----------------  SET O2 TO CONCO2 IF O2 INACTIVE ------------------
C
      IF (IOXYGEN.GT.NGAS) THEN
         DO KLOOP          = 1, KTLOOP
            CBLK(KLOOP,IOXYGEN) = CONCO2(KLOOP)
         ENDDO
      ENDIF
C
C *********************************************************************
C *     INTERANNUALLY-VARYING CH4 CONCENTRATION (bnd, bmy, 7/1/03)    *
C *********************************************************************
C
      ! Test if CH4 is defined as an inert SMVGEAR II species
      IF ( ICH4 > NGAS ) THEN

         ! Loop over boxes per grid block
         DO KLOOP = 1, KTLOOP 

            ! 1-D grid box index
            JLOOP = KLOOP + JLOOPLO

            ! Grid-box latitude index
            YLAT  = GET_YMID( IYSAVE(JLOOP) )
            
            ! Pick the CH4 concentration [ppbv] for the proper lat bin
            ! CH4 values are read in "chemdr.f" (outside the parallel loop)
            IF ( YLAT < -30d0 ) THEN
               CBLK(KLOOP,ICH4) = C3090S
            ELSE IF ( YLAT >= -30d0 .and. YLAT < 0d0  ) THEN
               CBLK(KLOOP,ICH4) = C0030S
            ELSE IF ( YLAT >=   0d0 .and. YLAT < 30d0 ) THEN
               CBLK(KLOOP,ICH4) = C0030N
            ELSE
               CBLK(KLOOP,ICH4) = C3090N
            ENDIF

            ! Convert from [ppbv CH4] to [molec CH4/cm3]
            CBLK(KLOOP,ICH4) = CBLK(KLOOP,ICH4) *1d-9 * AIRDENS(JLOOP)
         ENDDO
      ENDIF
C
C *********************************************************************
C *   CALCULATE KINETIC REACTION RATES USING ARRHENIUS PARAMETERS     * 
C *********************************************************************
C REACTION RATES HAVE THE FORM K = A * (300 / T)**B * EXP(C / T)
C
C NARR = NUMBER OF REACTIONS OF THE FORM K = A
C NABR = NUMBER OF REACTIONS OF THE FORM K = A * (300 / T)**B
C NACR = NUMBER OF REACTIONS OF THE FORM K = A                * EXP(C / T)
C NABC = NUMBER OF REACTIONS OF THE FORM K = A * (300 / T)**B * EXP(C / T)
C NKARR, NKBRR, NKACR, NKABC = REACTION RATE NUMBERS OF EACH
C NARR,  NABR,  NACR,  NABC  REACTION 
C
         DO 37 I           = 1, NARR(NCS)
            NK               = NKARR(I,NCS)
            ARRNK            = ARR(NK,NCS)
            DO 35 KLOOP      = 1, KTLOOP
               RRATE(KLOOP,NK) = ARRNK
 35         CONTINUE
 37      CONTINUE
         DO 42 I           = 1, NABR(NCS)
            NK               = NKABR(I,NCS)
            ARRNK            = ARR(NK,NCS)
            DO 40 KLOOP      = 1, KTLOOP
               RRATE(KLOOP,NK) = ARRNK * TEMP1(KLOOP)**BRR(NK,NCS)
 40         CONTINUE
 42      CONTINUE
C
         DO 47 I           = 1, NACR(NCS) 
            NK               = NKACR(I,NCS)
            ARRNK            = ARR(NK,NCS)  
            DO 45 KLOOP      = 1, KTLOOP
               RRATE(KLOOP,NK) = ARRNK * EXP(KCRR(NK,NCS) / T3K(KLOOP)) 
 45         CONTINUE
 47      CONTINUE
C
         DO 52 I           = 1, NABC(NCS) 
            NK               = NKABC(I,NCS)
            ARRNK            = ARR(NK,NCS)  
            DO 50 KLOOP      = 1, KTLOOP
               RRATE(KLOOP,NK) = ARRNK * TEMP1(KLOOP)**BRR(NK,NCS)
     1              * EXP(KCRR(NK,NCS) / T3K(KLOOP))
 50         CONTINUE
 52      CONTINUE

C
C *********************************************************************
C ******                   SET EMISSION RATES                    ******
C *********************************************************************
C
      DO I = 1,NEMIS(NCS)
C get tracer number corresponding to emission species I
         NN = IDEMS(I)
         IF (NN.NE.0) THEN
C find reaction number for emission of tracer NN
            NK = NTEMIS(NN,NCS)
            IF (NK.NE.0) THEN
               DO KLOOP = 1,KTLOOP
                  JLOOP = LREORDER(KLOOP+JLOOPLO)
                  RRATE(KLOOP,NK) = REMIS(JLOOP,I)
               ENDDO
            ENDIF
         ENDIF
      ENDDO
C
C *********************************************************************
C ******                SET DRY DEPOSITION RATES                 ******
C ******                                                         ******
C ******   NOTE: Now compute drydep throughout the mixed layer   ******
C ******   (a.k.a. PBL) in order to prevent short-lived species  ******
C ******   such as HNO3 from being depleted in the shallow       ******
C ******   surface layer. (rjp, bmy, 7/30/03)                    ******   
C *********************************************************************
C
      DO I = 1,NDRYDEP(NCS)
         NK = NTDEP(I)
         IF (NK.NE.0) THEN
            DO KLOOP = 1,KTLOOP

               ! 1-D grid box index (accounts for reordering)
               JLOOP = LREORDER(KLOOP+JLOOPLO)

               ! 3-D grid box index
               IX    = IXSAVE(JLOOP)
               IY    = IYSAVE(JLOOP)
               IZ    = IZSAVE(JLOOP)
               
               ! Now compute drydep throughout the entire PBL
               ! GET_FRAC_UNDER_PBLTOP returns the fraction of layer
               ! (IX, IY, IZ) that is beneath the PBL top
               RRATE(KLOOP,NK) = DEPSAV(IX,IY,I) * 
     &                           GET_FRAC_UNDER_PBLTOP( IX, IY, IZ )
            ENDDO
         ENDIF
      ENDDO
C
C *********************************************************************
C ********  MULTIPLY RATES BY CONSTANT SPECIES CONCENTRATIONS  ********
C *       (EITHER M, O2, N2, OR ANY ACTIVE OR INACTIVE SPECIES)       *   
C *********************************************************************
C NMAIR    = # REACTIONS WHERE THE SPECIES IN THE THIRD POSITION IS 
C              IS 'M' = 'O2 + N2'
C NMO2     = # REACTIONS WHERE THE SPECIES IN THE THIRD POSITION IS O2 
C NMN2     = # REACTIONS WHERE THE SPECIES IN THE THIRD POSITION IS N2 
C NMOTH    = # OCCURENCES OF SPECIES IN THIRD POSITION THAN ARE NOT
C              O2, N2, OR M, OR OF SPECIES IN ANY POSITION THAT ARE
C              INACTIVE.
C LGASBINO = JGAS (SET IN READCHEM AND GASCONC)
C
         DO 72 I           = 1, NMAIR(NCS)   
            NK               = NREACAIR(I,NCS)
            DO 70 KLOOP      = 1, KTLOOP 
               RRATE(KLOOP,NK) = RRATE(KLOOP,NK) * DENAIR(KLOOP) 
 70         CONTINUE
 72      CONTINUE
C
         DO 82 I           = 1, NMO2(NCS) 
            NK               = NREACO2(I,NCS)
            DO 80 KLOOP      = 1, KTLOOP 
               RRATE(KLOOP,NK) = RRATE(KLOOP,NK) * CONCO2(KLOOP) 
 80         CONTINUE
 82      CONTINUE
C        
         DO 92 I           = 1, NMN2(NCS)  
            NK               = NREACN2(I,NCS)
            DO 90 KLOOP      = 1, KTLOOP 
               RRATE(KLOOP,NK) = RRATE(KLOOP,NK) * CONCN2(KLOOP) 
 90         CONTINUE
 92      CONTINUE
C
C *********************************************************************
C *                   PRESSURE-DEPENDENT EFFECTS                      * 
C * ADD THE THIRD BODY EFFECT FOR PRESSURE DEPENDENCE OF RATE         *
C * COEFFICIENTS. THE REACTIONS WERE READ IN IN PAIRS (LOW AND HIGH   *
C * PRESSURE LIMITS) WITH THE SPECIFIC INDICATOR, 'P'  IN THE INPUT   *
C * DATA SET. SEE DEMORE ET AL. (1990) JPL 90-1 FOR MORE DETAILS      *
C *********************************************************************
C NPRESM  = # PRESSURE DEPENDENT 3-BODY REACTIONS  
C FCV     = CHARACTERIZES FALLOFF CURVE (SEE ATKINSON ET. AL (1992)
C           J. PHYS. CHEM. REF. DATA 21, P. 1145). USUALLY = 0.6 
C           HOWEVER, TWO TEMPERATURE-DEPENDENT EXPRESSIONS ARE ALSO USED: 
C             FCV = EXP(-T/FCT1)  OR 
C             FCV = EXP(-T/FCT1)+EXP(-FCT2/T)
C RATE(NK)   = K(0,T)[M], WHERE K(0,T) = 3-BODY, LOW PRESSURE LIMIT COEF. 
C RATE(NK+1) = K(INF,T) = 2-BODY, HIGH PRESSURE LIMIT COEF. 
C
         DO 165 I           = 1, NPRESM(NCS) 
            NK                = NREACPM(I,NCS)
            FCVNK             = FCV(    NK,NCS) 
            FCT1              = FCTEMP1(NK,NCS)
            FCT2              = FCTEMP2(NK,NCS)
            IF (FCT2.NE.0) THEN 
               DO 150 KLOOP     = 1, KTLOOP 
                  FCT            = EXP(-T3K(KLOOP) / FCT1)
     1                 + EXP(-FCT2       / T3K(KLOOP)) 
                  XYRAT          = RRATE(KLOOP,NK) / RRATE(KLOOP,NK+1) 
                  BLOG           = LOG10(XYRAT)
                  FEXP           = 1.d0 / (1.d0 + BLOG * BLOG)
                  RRATE(KLOOP,NK)= RRATE(KLOOP,NK)*FCT**FEXP/(1d0+XYRAT) 
 150           CONTINUE
            ELSEIF (FCT1.NE.0.) THEN 
               DO 155 KLOOP     = 1, KTLOOP 
                  FCT            = EXP(-T3K(KLOOP) / FCT1)
                  XYRAT          = RRATE(KLOOP,NK) / RRATE(KLOOP,NK+1) 
                  BLOG           = LOG10(XYRAT)
                  FEXP           = 1.d0 / (1.d0 + BLOG * BLOG)
                  RRATE(KLOOP,NK)= RRATE(KLOOP,NK)*FCT**FEXP/(1d0+XYRAT)  
 155           CONTINUE 
            ELSE
               DO 160 KLOOP     = 1, KTLOOP 
                  XYRAT          = RRATE(KLOOP,NK) / RRATE(KLOOP,NK+1) 
                  BLOG           = LOG10(XYRAT)
                  FEXP           = 1.d0 / (1.d0 + BLOG * BLOG)
                 RRATE(KLOOP,NK)=RRATE(KLOOP,NK)*FCVNK**FEXP/(1d0+XYRAT)  
 160           CONTINUE
            ENDIF
 165     CONTINUE
C
C *********************************************************************
C * SET THE RATES OF ALL THERMALLY DISSOCIATING SPECIES. SEE DEMORE   *
C * ET AL. (1990). CHEMICAL KINETICS AND PHOTOCHEMICAL DATA FOR USE   *
C * IN STRATOSPHERIC MODELING. JPL. 90-1, P. 93. THE RATE HAS THE     *
C * FORM Kf / [A EXP (C / T)], WHERE Kf IS THE REACTION IN THE        *
C * REVERSE DIRECTION.                                                * 
C *********************************************************************
C NNEQ      = # THERMALLY DISSOCIATING EQUILIBRIUM REACTIONS. PREVIOUS
C             EQUATION MUST BE PRESSURE-DEPENDENT. 
C RATE(NK1) = CM3 MOLEC.-1 S-1 (BIMOLECULAR RATE FROM PRESSURE-DEPEND)
C RATE(NK)  = CM3 MOLEC.-1 (EQUILIBRIUM CONSTANT) (BEFORE CALCULATION)
C RATE(NK)  = S-1 (UNIMOLECULAR RATE AFTER CALCULATION) 
C
         DO 182 I           = 1, NNEQ(NCS) 
            NK                = NREACEQ(I,NCS)
            NK1               = NREQOTH(I,NCS)
            DO 180 KLOOP      = 1, KTLOOP
               RRATE(KLOOP,NK)  = RRATE(KLOOP,NK1) / RRATE(KLOOP,NK)
 180        CONTINUE
 182     CONTINUE 
C
C *********************************************************************
C     MULTIPLY RATE COEFFICIENT BY OTHER INACTIVE CONCENTRATIONS
C *********************************************************************
C THIS LOOP MUST OCCUR AFTER EQUILIBRIUM REACTIONS 
C
         DO 183 I          = 1, NMOTH(NCS) 
            NK               = NREACOTH(I,NCS)
            JOLD             = LGASBINO(I,NCS)
            DO 181 KLOOP     = 1, KTLOOP 
               RRATE(KLOOP,NK) = RRATE(KLOOP,NK) * CBLK(KLOOP,JOLD) 
 181        CONTINUE
 183     CONTINUE
C
C *********************************************************************
C *                       SET SPECIAL RATES                           *
C *********************************************************************
C
C ---  K = K1 + K2  ---- 
         IF (NKSPECW(NCS).GT.0) THEN
            NK         = NKSPECW( I )
            DO KLOOP   = 1, KTLOOP
               RRATE(KLOOP,NK) = RRATE(KLOOP,NK) + RRATE(KLOOP,NK+1)
            ENDDO
         ENDIF
C
C ---  K = K1*FYRNO3(K2,M,T)  ---   addition branch of RO2+NO
         DO I          = 1, NNADDA(NCS)
            NK         = NKSPECA( I,NCS )
            DO KLOOP   = 1, KTLOOP
               RRATE(KLOOP,NK) = RRATE(KLOOP,NK) * 
     +              FYRNO3(RRATE(KLOOP,NK+1),DENAIR(KLOOP),T3K(KLOOP))
            ENDDO
         ENDDO
C
C ---  K = K1*(1-FYRNO3(K2,M,T))  ---  abstraction branch of RO2+NO
         DO I          = 1, NNADDB(NCS)
            NK         = NKSPECB( I,NCS )
            DO KLOOP   = 1, KTLOOP
               RRATE(KLOOP,NK) = RRATE(KLOOP,NK) *
     $              (1.D0 - FYRNO3(RRATE(KLOOP,NK+1), DENAIR(KLOOP), 
     $              T3K(KLOOP)))
            ENDDO
         ENDDO
C
C ---  K = K1*([O2]+3.5D18)/(2*[O2]+3.5D18) --- HO2+2*CO branch of GLYX+OH/NO3
         DO I          = 1, NNADDC(NCS)
            NK         = NKSPECC( I,NCS )
            DO KLOOP   = 1, KTLOOP
               RRATE(KLOOP,NK) = RRATE(KLOOP,NK) * 
     +              (CONCO2(KLOOP)+3.5D18)/(2.D0*CONCO2(KLOOP)+3.5D18)
            ENDDO
         ENDDO
C
C ---  K = K1*[O2]/(2*[O2]+3.5D18)  --- GLCO3 branch of GLYX+OH/NO3
         DO I          = 1, NNADDD(NCS)
            NK         = NKSPECD( I,NCS )
            DO KLOOP   = 1, KTLOOP
               RRATE(KLOOP,NK) = RRATE(KLOOP,NK) * 
     +              (CONCO2(KLOOP))/(2.D0*CONCO2(KLOOP)+3.5D18)
            ENDDO
         ENDDO
C
C ---  OH + HNO3:   K = K0 + K3[M] / (1 + K3[M]/K2)  ------
         IF (NKSPECX(NCS).GT.0) THEN
            NK               = NKSPECX(NCS)
            DO KLOOP     = 1, KTLOOP
            RRATE2=RRATE(KLOOP,NK+2)*DENAIR(KLOOP)
            RRATE(KLOOP,NK) = RRATE(KLOOP,NK) + RRATE2 /
     1           (1.D0 + RRATE2 / RRATE(KLOOP,NK+1))    
            ENDDO
         ENDIF
C
C ---    OH + CO: K = K0(1+0.6 Patm)  ------------ 
C    CONSTC includes a factor to convert PRESS3 from (dyn cm-2) to (atm)
         IF (NKSPECY(NCS).GT.0) THEN
            NK           = NKSPECY(NCS)
            CONSTC       = 0.6D0 * 9.871D-07
            DO KLOOP     = 1, KTLOOP
               JLOOP           = LREORDER(JLOOPLO + KLOOP)
               RRATE(KLOOP,NK) = RRATE(KLOOP,NK) *
     1              (1.D0 + CONSTC*PRESS3(JLOOP))
            ENDDO
         ENDIF
C
C ---    MCO3+MO2:  K = K1 / (1+K2)   ------------
C  temperature-dependent branching ratio
         DO I          = 1,NNADDV(NCS)
            NK         = NKSPECV( I,NCS )
            DO KLOOP   = 1, KTLOOP
               RRATE(KLOOP,NK)=RRATE(KLOOP,NK)/(1.d0+RRATE(KLOOP,NK+1))
            ENDDO
         ENDDO
C
         ! Add special reaction for DMS + OH + O2 (bdf, bmy, 4/18/03)
         DO I          = 1,NNADDG(NCS)
            NK         = NKSPECG( I,NCS )
            DO KLOOP   = 1, KTLOOP
               RRATE(KLOOP,NK)=RRATE(KLOOP,NK)/
     &              (1.d0+RRATE(KLOOP,NK+1)*CONCO2(KLOOP))
            ! SMVGEARII doesn't have structure to multiply rate(nk+1) by [O2]
            ENDDO
         ENDDO
C
C ---  HO2/NO3 + HO2:  K = (K1 + K2)*(1+1.4E-21*[H2O]*EXP(2200/T))  --- 
C  dependence of HO2/NO3 + HO2 on water vapor
         IF (NKSPECZ(NCS).GT.0) THEN
            NK         = NKSPECZ(NCS)
            DO KLOOP   = 1, KTLOOP
               RRATE(KLOOP,NK) =
     +          (RRATE(KLOOP,NK)+RRATE(KLOOP,NK+1)*DENAIR(KLOOP)) * 
     +          (1.D0+1.4D-21*CBLK(KLOOP,IH2O)*EXP(2200.D0/T3K(KLOOP)))
            ENDDO
         ENDIF

      !=================================================================
      ! Perform loss on wet aerosol 
      !=================================================================

      ! Set HETCHEM = T to perform het chem on aerosols
      HETCHEM = .TRUE.

      DO KLOOP = 1, KTLOOP

         ! 1-D grid box index
         JLOOP = LREORDER(JLOOPLO+KLOOP)

         IF ( HETCHEM ) THEN

            !===========================================================
            ! Perform heterogeneous chemistry on sulfate aerosol
            ! plus each of the NDUST dust size bins from FAST-J
            !===========================================================
            XDENA   = DENAIR(KLOOP)
            XSTK    = SQRT(T3K(KLOOP))

            DO I       = 1, NNADDK(NCS)
               NK      = NKSPECK(I,NCS)
               XSQM    = SQRT(ARR(NK,NCS))

               ! Initialize
               RRATE(KLOOP,NK) = 0d0
               DUMMY2(KLOOP)   = 0d0

               ! Sum up total surface area for all aerosol types
               ! so that we can use it for the planeflight diagnostic
               ! (mje, bmy, 8/7/03)
               IF ( NK == NKN2O5 ) THEN
                  TOTAREA = 0.d0
                  DO N = 1, NDUST + NAER
                     TOTAREA = TOTAREA + TAREA(JLOOP,N)
                  ENDDO
               ENDIF

               ! Loop over sulfate and other aerosols
               DO N = 1, NDUST + NAER

                  ! Surface area of aerosol [cm2 aerosol/cm3 air]
                  XAREA = TAREA(JLOOP,N) 
                  
                  ! Test if N2O5 hydrolysis rxn
                  IF ( NK == NKN2O5 ) THEN
                     
                     ! Get GAMMA for N2O5 hydrolysis, which is
                     ! a function of aerosol type, temp, and RH
                     XSTKCF = N2O5( N, T3K(KLOOP), ABSHUMK(KLOOP) )

                     ! Archive N2O5 hydrolysis for ND40 diagnostic
                     DUMMY2(KLOOP) = DUMMY2(KLOOP) +
     &                               ( XAREA / TOTAREA * XSTKCF )
                    
                  ELSE
                     
                     ! Get GAMMA for species other than N2O5
                     XSTKCF = BRR(NK,NCS)

                  ENDIF

                  ! Radius for dust size bin N
                  XRADIUS = ERADIUS(JLOOP,N) 

                  ! Reaction rate for dust size bin N
                  RRATE(KLOOP,NK) = RRATE(KLOOP,NK) + 
     $                 ARSL1K(XAREA,XRADIUS,XDENA,XSTKCF,XSTK,XSQM)
               ENDDO
            ENDDO
         ELSE

            !===========================================================
            ! Don't perform heterogeneous chemistry at all
            !===========================================================
            XAREA   = TAREA(JLOOP,1)
            XRADIUS = ERADIUS(JLOOP,1)
            XDENA   = DENAIR(KLOOP)
            XSTK    = SQRT(T3K(KLOOP))
            DO I       = 1, NNADDK(NCS)
               NK      = NKSPECK(I,NCS)
               XSQM    = SQRT(ARR(NK,NCS))
               XSTKCF  = BRR(NK,NCS)
               RRATE(KLOOP,NK) =
     &              ARSL1K(XAREA,XRADIUS,XDENA,XSTKCF,XSTK,XSQM)
            ENDDO
         ENDIF
      ENDDO

      ENDIF
C     ENDIF NCS.EQ.1 OR 2
C
C *********************************************************************
C *********************************************************************
C *             REORDER RRATE ARRAY THEN CALL SOLVER                  * 
C *                                                                   *
C * NOTE: If after this point you want to reference a SMVGEAR rxn #,  *
C * then you must use NOLDFNEW(NK,1) instead of NK. (bmy, 8/8/03)     *
C *********************************************************************
C
      NFDH3              = ITHRR(NCS) 
      NFDL2              = NFDH3  + 1 
      NFDREP             = INOREP(NCS)
      NFDREP1            = NFDREP + 1
      NFDH2              = NFDH3  + ITWOR(NCS) 
      NFDL1              = NFDH2  + 1
      NFDH1              = NFDH2  + IONER(NCS) 
      NFDL0              = NFDH1  + 1 
      NALLR              = NALLRAT(NCS) 

C
      DO 730 NKN         = 1, NALLR
         NK                = NOLDFNEW(NKN,NCS)
         IRMA(NKN)         = IRM2(1,NK,NCS) 
         IRMB(NKN)         = IRM2(2,NK,NCS) 
         IRMC(NKN)         = IRM2(3,NK,NCS) 
 730  CONTINUE 
C
C *********************************************************************
C                        REORDER RRATE ARRAY 
C *********************************************************************
C                 TRATE USED HERE AS A DUMMY ARRAY 
C *********************************************************************
C
C
      DO 745 NK          = 1, NTRATES(NCS)
         DO 740 KLOOP      = 1, KTLOOP
            TRATE(KLOOP,NK)  = RRATE(KLOOP,NK)
 740     CONTINUE
 745  CONTINUE
C
      DO 755 NKN         = 1, NALLR
         NK                = NOLDFNEW(NKN,NCS)
         DO 750 KLOOP      = 1, KTLOOP
            RRATE(KLOOP,NKN) = TRATE(KLOOP,NK)
 750     CONTINUE
 755  CONTINUE
C
C *********************************************************************
C REPLACE INACTIVE REACTION RATE COEFFICIENT ARRAY WITH NEW ARRAY 
C THESE REACTIONS HAVE NO ACTIVE LOSS TERMS. PHOTORATE TERMS HERE
C ARE REPLACED IN UPDATE.F .
C *********************************************************************
C                 TRATE USED HERE AS A REAL ARRAY 
C *********************************************************************
C
      DO 765 NKN          = NFDL0, NALLR
         NH                 = NKN + NALLR
         DO 760 KLOOP       = 1, KTLOOP
            TRATE(KLOOP,NKN)  =  RRATE(KLOOP,NKN)
            TRATE(KLOOP,NH)   = -RRATE(KLOOP,NKN)
 760     CONTINUE
 765  CONTINUE
C
C *********************************************************************
C              Photorates for Harvard Geos Code
C *********************************************************************
C PRATE           = PHOTORATE (SEC-1) IF SUN IS DOWN, PRATE = 0.
C RRATE           = RATE COEFFICIENT (SEC-1)
C NRATES          = NUMBER OF KINETIC REACTION RATES.
C
C
      IF(IFSUN.EQ.1) THEN
         DO I                = 1, JPHOTRAT(NCS)
            NK               = NRATES(NCS) + I
            NKN              = NKNPHOTRT(I,NCS)
            SPECNAME         = NAMEGAS(IRM(1,NK,NCS))
            IFNC             = DEFPRAT(NK,NCS) + 0.01D0
            IBRCH            = 10.D0*(DEFPRAT(NK,NCS)-IFNC) + 0.5D0

#if   defined( LSLOWJ )
            ! ISPEC is only needed for SLOW-J photolysis
            ISPEC            = INAME(I)
#endif

            DO KLOOP            = 1, KTLOOP 
               JLOOP            = LREORDER(KLOOP+JLOOPLO)

               ! Translate 1-D to 3-D grid box indices
               IX               = IXSAVE(JLOOP)
               IY               = IYSAVE(JLOOP)
               IZ               = IZSAVE(JLOOP)                  

               ! Get cosine( SZA ) using 1-D array index
               IJWINDOW         = (IY-1)*IIPAR + IX
               GMU              = SUNCOS(IJWINDOW)

               ! For daylight boxes...
               IF(GMU.GT. 0.D0) THEN

#if   defined( LFASTJ )

                  ! For FAST-J, get photorates from fjfunc.f
                  RRATE(KLOOP,NKN)  = FJFUNC(IX,IY,IZ,I,IBRCH,SPECNAME)

#elif defined( LSLOWJ )

                  ! For SLOW-J, get photorates from rtfunc.f
                  TT               = T3K(KLOOP)
                  PRESMB           = PRESS3(JLOOP) * 1.D-3
                  RRATE(KLOOP,NKN) = RTFUNC(GMU,ISPEC,IFNC,IBRCH,TT
     &                                     ,PRESMB,JLOOP)

#endif

!### Debug: warn if there are negative J-values, for either 
!### FAST-J or SLOW-J photolysis (bmy, 10/1/98)
!###                  IF ( RRATE(KLOOP,NK) < 0 ) THEN
!###                     PRINT*, 'CALCRATE.F: J-Value < 0: ', 
!###     &                  IX, IY, IZ, IBRCH, SPECNAME, KLOOP, NK,
!###     &                  RRATE(KLOOP,NK)
!###                  ENDIF

               ELSE

                  ! Nighttime: photorates are zero
                  RRATE(KLOOP,NKN)  = 0.D0

               ENDIF
            ENDDO
         ENDDO

         !==============================================================
         ! HARDWIRE addition of 1e-5 s-1 photolysis rate to 
         ! HNO4 -> HO2+NO2 to account for HNO4 photolysis in near-IR -- 
         ! see Roehl et al. 'Photodissociation of peroxynitric acid in 
         ! the near-IR', 2002. (amf, bmy, 1/7/02)
         !
         ! Add NCS index to NKHNO4 for SMVGEAR II (gcc, bmy, 4/1/03)
         !==============================================================
         IF ( NKHNO4(NCS) > 0 ) THEN

            ! Put J(HNO4) in correct spot for SMVGEAR II
            PHOTVAL = NKHNO4(NCS) - NRATES(NCS)
            NKN     = NKNPHOTRT(PHOTVAL,NCS)

            DO KLOOP=1,KTLOOP
               RRATE(KLOOP,NKN)=RRATE(KLOOP,NKN) + 1d-5
            ENDDO
         ENDIF

         !==============================================================
         ! SPECIAL TREATMENT FOR O3+hv -> OH+OH
         ! [O1D]ss=J[O3]/(k[H2O]+k[N2]+k[O2])
         ! SO, THE EFFECTIVE J-VALUE IS J*k[H2O]/(k[H2O]+k[N2]+k[O2])
         !
         ! Add NCS index to NKHNO4 for SMVGEAR II (gcc, bmy, 4/1/03)
         !==============================================================
         IF ( NKO3PHOT(NCS) > 0 ) THEN

            ! Put J(O3) in correct spot for SMVGEAR II
            PHOTVAL = NKO3PHOT(NCS) - NRATES(NCS)
            NKN     = NKNPHOTRT(PHOTVAL,NCS)

            DO KLOOP = 1, KTLOOP

               ! Save old value of J-O3 in a diagnostic array 
               ! (gcc, bmy, 4/1/03)
               DUMMY(KLOOP) = RRATE(KLOOP,NKN)

               !========================================================
               ! Change rate of O(1D)+ N2 to be 3.1e-11 at 298K rather
               ! than 2.6e-11.  The temperature dependence remains the
               ! same, so the constant changes from 1.8e-11 to 2.14e-11
               ! according to Heard, pers. comm.,2002. (amf, bmy, 1/7/02)
               !========================================================
               ! Change the rate of O(1D)+H2O from 2.2e-10 to 1.45e-10*
               ! exp(89/temp) on the basis of Dunlea and Ravishankara
               ! 'Measurement of the Rate coefficient for the reaction 
               ! of O(1D) with H2O and re-evaluation of the atmospheric
               ! OH Production Rate'.  One of the RSC Journals
               ! (mje 4/5/04)
               !========================================================
               RRATE(KLOOP,NKN) = RRATE(KLOOP,NKN) *
     $            1.45d-10 * EXP( 89.d0*T3I(KLOOP)) * CBLK(KLOOP,IH2O) /
     $          ( 1.45d-10 * EXP( 89.d0*T3I(KLOOP)) * CBLK(KLOOP,IH2O) +
     $            2.14d-11 * EXP(110.d0*T3I(KLOOP)) * CONCN2(KLOOP)    +
     $            3.20d-11 * EXP( 70.d0*T3I(KLOOP)) * CONCO2(KLOOP)    )


            ENDDO
         ENDIF
      ELSEIF(IFSUN.EQ.2) THEN
         DO I          = 1, JPHOTRAT(NCS)
            NKN        = NKNPHOTRT(I,NCS)
            DO KLOOP   = 1, KTLOOP
               RRATE(KLOOP,NKN)  = 0.D0
            ENDDO
         ENDDO
      ELSE
         ! ERROR IN IFSUN
         CALL ERROR_STOP( 'ERROR in IFSUN -- STOP 0345', 'calcrate.f' )
      ENDIF

      !=================================================================
      ! Store J-values for 5 rxns + POH in ND22 diagnostic 
      !=================================================================
      IF ( ND22 > 0 ) THEN
         DO I  = 1, JPHOTRAT(NCS)
            NK  = NRATES(NCS) + I
            NKN = NKNPHOTRT(I,NCS)

            ! Name of species being photolyzed
            SPECNAME = NAMEGAS(IRM(1,NK,NCS))

            SELECT CASE ( TRIM( SPECNAME ) )
               CASE ( 'NO2' )
                  INDEX = 1
               CASE ( 'HNO3' )
                  INDEX = 2
               CASE ( 'H2O2' )
                  INDEX = 3
               CASE ( 'CH2O' )
               !CASE ( 'ACET' )  ! for testing (bey, 1/7/99)
                  INDEX = 4
               CASE ( 'O3'   )
                  INDEX = 5
               CASE DEFAULT
                  INDEX = 0
            END SELECT

            ! If this is not one of the 5 J-values, go to next reaction
            IF ( INDEX == 0 ) CYCLE

            ! Loop over I-J-L boxes
            DO KLOOP = 1, KTLOOP
               JLOOP = LREORDER( KLOOP + JLOOPLO )

               ! I-J-L indices
               IX = IXSAVE(JLOOP)
               IY = IYSAVE(JLOOP)
               IZ = IZSAVE(JLOOP)

               ! Save J-values for 2PM diagnostic boxes
               ! Use AD22 array for J-value diagnostic (bmy, 9/30/99)
               IF ( LTJV(IX,IY) > 0 .and. IZ <= LD22 ) THEN
                  IF ( INDEX == 5 ) THEN

                     ! Store unadjusted J-O3 as AD22(:,:,:,5)
                     AD22(IX,IY,IZ,5) =
     &                    AD22(IX,IY,IZ,5) + DUMMY(KLOOP)

                     ! Store POH as AD22(:,:,:,6)
                     AD22(IX,IY,IZ,6) =
     &                    AD22(IX,IY,IZ,6) + RRATE(KLOOP,NKN)
                  ELSE
                     ! Store other J-Values in their appropriate slots
                     AD22(IX,IY,IZ,INDEX) =
     &                    AD22(IX,IY,IZ,INDEX) + RRATE(KLOOP,NKN)
                  ENDIF
               ENDIF
            ENDDO
         ENDDO
      ENDIF
C
C *********************************************************************
C                         RESET NCSP  
C *********************************************************************
C NCS       = 1..NCSGAS FOR GAS CHEMISTRY                           
C NCSP      = NCS       FOR DAYTIME   GAS CHEM            
C           = NCS + ICS FOR NIGHTTIME GAS CHEM           
C
      NCSP            = (IFSUN - 1) * ICS + NCS
C
C *********************************************************************
C                ARCHIVE FOR PLANE-FOLLOWING DIAGNOSTIC   
C *********************************************************************
C
      ! Pass JO1D and N2O5 to "planeflight_mod.f" (mje, bmy, 8/7/03)
      CALL ARCHIVE_RXNS_FOR_PF( DUMMY, DUMMY2 )
C
C *********************************************************************
C                     RETURN TO CALLING PROGRAM
C *********************************************************************
C
      RETURN
C
C *********************************************************************
C                       INTERNAL SUBROUTINES 
C *********************************************************************
C
      CONTAINS

      FUNCTION N2O5( AEROTYPE, TEMP, RH ) RESULT( GAMMA )

      !=================================================================
      ! Internal function N2O5 computes the GAMMA sticking factor
      ! for N2O5 hydrolysis. (mje, bmy, 8/7/030
      ! 
      ! Arguments as Input:
      ! ----------------------------------------------------------------
      ! (1 ) AEROTYPE (INTEGER) : # denoting aerosol type (cf FAST-J)
      ! (2 ) TEMP     (REAL*8 ) : Temperature [K]
      ! (3 ) RH       (REAL*8 ) : Relative Humidity [fraction]
      !
      ! NOTES:
      !=================================================================
      
      ! Arguments
      INTEGER, INTENT(IN) :: AEROTYPE
      REAL*8,  INTENT(IN) :: TEMP, RH

      ! Local variables
      REAL*8              :: RH_P, FACT, TTEMP

      ! Function return value
      REAL*8              :: GAMMA
      
      !=================================================================
      ! N2O5 begins here!
      !=================================================================

      ! Convert RH to % (max = 100%)
      RH_P  = MIN( RH * 100d0, 100d0 )

      ! Default value
      GAMMA = 0.01d0

      ! Special handling for various aerosols
      SELECT CASE ( AEROTYPE )

         !----------------
         ! Dust 
         !----------------
         CASE ( 1, 2, 3, 4, 5, 6, 7 )      
                                
            ! Based on unpublished Crowley work
            GAMMA = 0.01d0

         !----------------
         ! Sulfate
         !----------------
         CASE ( 8 )            
    
            !===========================================================
            ! RH dependence from Kane et al., Heterogenous uptake of 
            ! gaseous N2O5 by (NH4)2SO4, NH4HSO4 and H2SO4 aerosols
            ! J. Phys. Chem. A , 2001, 105, 6465-6470 
            !===========================================================
            GAMMA = 2.79d-4 + RH_P*(  1.30d-4 + 
     &                        RH_P*( -3.43d-6 + 
     &                        RH_P*(  7.52d-8 ) ) )

            !===========================================================
            ! Temperature dependence factor (Cox et al, Cambridge UK) 
            ! is of the form:
            !
            !          10^( LOG10( G294 ) - 0.04 * ( TTEMP - 294 ) )
            ! FACT = -------------------------------------------------
            !                     10^( LOG10( G294 ) )
            !
            ! Where G294 = 1e-2 and TTEMP is MAX( TEMP, 282 ).
            ! 
            ! For computational speed, replace LOG10( 1e-2 ) with -2
            ! and replace 10^( LOG10( G294 ) ) with G294 
            !===========================================================
            TTEMP = MAX( TEMP, 282d0 )
            FACT  = 10.d0**( -2d0 - 4d-2*( TTEMP - 294.d0 ) ) / 1d-2

            ! Apply temperature dependence
            GAMMA = GAMMA * FACT

         !----------------
         ! Black Carbon
         !----------------
         CASE ( 9 )  

             ! From IUPAC
             GAMMA = 0.005d0

         !----------------
         ! Organic Carbon
         !----------------           
         CASE ( 10 )          

            !===========================================================
            ! Based on Thornton, Braban and Abbatt, 2003
            ! N2O5 hydrolysis on sub-micron organic aerosol: the effect
            ! of relative humidity, particle phase and particle size
            !===========================================================
            IF ( RH_P >= 57d0 ) THEN
               GAMMA = 0.03d0
            ELSE
               GAMMA = RH_P * 5.2d-4
            ENDIF

         !----------------
         ! Sea salt
         ! accum & coarse
         !----------------
         CASE ( 11, 12 )        
          
            ! Based on IUPAC recomendation
            IF ( RH_P >= 62 ) THEN 
               GAMMA = 0.03d0
            ELSE
               GAMMA = 0.005d0
            ENDIF

         !----------------         
         ! Default
         !----------------
         CASE DEFAULT
            WRITE (6,*) 'Not a suitable aerosol surface '
            WRITE (6,*) 'for N2O5 hydrolysis'
            WRITE (6,*) 'AEROSOL TYPE =',AEROTYPE
            CALL GEOS_CHEM_STOP

      END SELECT   
         
      ! Return to CALCRATE
      END FUNCTION N2O5
C
C *********************************************************************
C ******************** END OF SUBROUTINE CALCRATE *********************
C *********************************************************************
C
      END SUBROUTINE CALCRATE
