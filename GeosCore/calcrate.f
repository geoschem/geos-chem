! $Id: calcrate.f,v 1.6 2010/03/15 19:33:25 ccarouge Exp $
      SUBROUTINE CALCRATE( SUNCOS )
!
!******************************************************************************
!  Subroutine CALCRATE computes reaction rates before passing them to the
!  SMVGEAR solver.  (M. Jacobson 1997; gcc, bdf, bmy, 4/1/03, 11/19/08)
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
!  (8 ) SLOW-J is now obsolete; remove LSLOWJ #ifdef blocks (bmy, 6/23/05)
!  (9 ) Now use NUMDEP instead of NDRYDEP(NCS) for the loop limit over drydep
!        species.  NDRYDEP is the # of rxns in "globchem.dat", and NUMDEP is
!        the # of drydep species in GEOS-Chem.  The two values may not be the 
!        same. (dbm, phs, 11/19/08)
!  (10) Now use new gamma(HO2) based on Thornton, Jaegle, and McNeill
!       (JGR, 2008) (jaegle, 02/26/09)
!  (11) Added branching ratio for C2H4 oxidation and photolysis (tmf, 12/14/06)
!  (12) Outputs GLYX and MGLY J-values. (tmf, 1/31/06)
!  (13) Modified the dry deposition rate reference, such that gas tracers 
!         which appear after the depositing aerosols will be referenced 
!         correctly. (tmf, 11/08/06)
!  (14) Updated OH+CO and O(1D)+H2O rates. (jmao, 4/20/09)
!  (15) Add options for emissions and depositions for non-local PBL scheme.
!       (ccc, 5/21/09)
!  (16) Added interface to get rx rates for solver from kpp (phs,ks,dhk, 09/15/09)  
!  (17) FastJX has taken near-IR photolysis into account with cross section at 
!       574nm, so we don't need to add the 1e-5 to the JHNO4 rate anymore.
!       (jmao, bmy, 10/27/09)
!  (18) Add change yield for HAC. Make CONSTC a parameter. (fp, 2/2/10)
!  (10) Updated O3 photolysis reaction rate and DMS oxidation reaction rates. 
!       See http://wiki.seas.harvard.edu/geos-chem/index.php/Updating_standard_chemistry_with_JPL_10-6. (bhh, jmao, eam, 7/18/11)
!******************************************************************************
!
      ! References to F90 modules 
      ! Add CSPEC to extract HO2 concentration (jaegle 2/26/09)
      USE COMODE_MOD,      ONLY : ABSHUM, AIRDENS, ERADIUS, IXSAVE, 
     &                            IYSAVE, IZSAVE,  JLOP,    PRESS3,  
     &                            REMIS,  T3,      TAREA,    CSPEC
      ! Add AD52 for gamma_ho2 diagnostic (jaegle 02/26/09)
      USE DIAG_MOD,        ONLY : AD22,   LTJV, AD52
      USE DRYDEP_MOD,      ONLY : DEPSAV, NUMDEP, DEPNAME
      ! Add CHECK_VALUE (jaegle 2/26/09)
      USE ERROR_MOD,       ONLY : ERROR_STOP, GEOS_CHEM_STOP,
     &                            CHECK_VALUE,
     &                            DEBUG_MSG
      USE GRID_MOD,        ONLY : GET_YMID
      ! Add GET_PBL_TOP_L (jaegle 02/26/09)
      USE PBL_MIX_MOD,     ONLY : GET_FRAC_UNDER_PBLTOP,GET_PBL_TOP_L
      USE PLANEFLIGHT_MOD, ONLY : ARCHIVE_RXNS_FOR_PF
      ! Add IDHO2 for HO2 concentration and IS_ICE (jaegle 02/26/09)
      USE TRACERID_MOD,    ONLY : IDHO2
      USE DAO_MOD,         ONLY : IS_ICE
      ! KPP interface (phs,ks,dhk, 09/15/09)
      USE GCKPP_GLOBAL,    ONLY : IND
      USE LOGICAL_MOD,     ONLY : LNLPBL, LKPP, LPRT ! (Lin, 03/31/09)


      IMPLICIT NONE

#     include "CMN_SIZE"  ! Size parameters
#     include "CMN"       ! STT, etc
#     include "comode.h"  ! SMVGEAR II arrays
      ! added LD52 and ND52 (jaegle 2/26/09)
#     include "CMN_DIAG"  ! ND22, LD22, ND52, LD52
      ! added FRCLND (jaegle 02/26/09)
#     include "CMN_DEP"  ! FRCLND
      
      ! Local variables
      INTEGER :: KLOOP,JLOOP,I,NK,JOLD2,JOLD,NK1,NKN,NH,MLOOP,J
      INTEGER :: NKNH2O, NKNH2    ! (bhh, jmao, eam, 7/18/11)
      LOGICAL :: LJO1DH2, LJO1DH2O, LPCELL     ! (bhh, jmao, eam, 7/18/11)
      INTEGER :: NP,K,IFNC,IBRCH,IX,IY,IZ,IJWINDOW,INDEX,NN,ii,jj
      INTEGER :: PHOTVAL,KSUN

      REAL*8  :: ARRNK,FCVNK,FCT1,FCT2,FCT,XYRAT,BLOG,FEXP,RATE3M
      REAL*8  :: RO1D, RO1DplH2O, RO1DplH2   ! (bhh, jmao, eam, 7/18/11)
C      REAL*8  :: CONSTC,RATE2AIR,RATE3H2O,RIS,RST,TBEGIN,TFINISH
      ! CONSTC declared as a parameter by FP (hotp 7/31/09)
      REAL*8  :: RATE2AIR,RATE3H2O,RIS,RST,TBEGIN,TFINISH
      REAL*8  :: PBEG,PFIN,PBEGNEW,PFINNEW,TOFDAYB,TOFDAYE,HOURANGB
      REAL*8  :: HOURANGE,SINFUNCB,SINFUNCE,XAREA,XRADIUS,XSQM,RRATE2
      REAL*8  :: XSTKCF,GMU,SUNCOS(MAXIJ),DUMMY(KBLOOP),XDENA,XSTK
      REAL*8  :: TK,     CONSEXP, VPRESH2O
      ! New variables for new reations (jmao, 4/20/09)
      REAL*8  :: KHI1,KLO1,XYRAT1,BLOG1,FEXP1,KHI2,KLO2,XYRAT2,BLOG2
      REAL*8  :: FEXP2,KCO1,KCO2,KCO
      !FP_ISOP
      ! FP variables for isoprene chemistry (hotp 7/31/09)
      REAL*8   :: HAC_FRAC      
      REAL*8   :: GLYC_FRAC   

      ! External functions
      REAL*8, EXTERNAL :: RTFUNC, FYRNO3,  ARSL1K,  FJFUNC, FYHORO

      CHARACTER*4      :: SPECNAME

      ! Added for heterogeneous chemistry (bmy, 11/15/01, 8/7/03)
      LOGICAL          :: HETCHEM
      INTEGER          :: N
      REAL*8           :: SUMAREA, TOTAREA, DUMMY2(KBLOOP)

      ! Added for HO2 het uptake (jaegle, 2/26/09)
      REAL*8           :: HO2_MOLEC_CM3, DUMMY3(KBLOOP)
      INTEGER          :: CONTINENTAL_PBL

      ! For grid-box latitude (bnd, bmy, 7/1/03)
      REAL*8           :: YLAT

      ! Variables from "comode.h" which are only ever used in "calcrate.f"
      ! Remove them from "comode.h" and the THREADPRIVATE declarations
      ! (bmy, 7/28/03) 
      REAL*8           :: ABSHUMK(KBLOOP), DENAIR(KBLOOP)
      REAL*8           :: CONCO2(KBLOOP),  CONCN2(KBLOOP)
      REAL*8           :: CONCH2(KBLOOP)     ! (bhh, jmao, eam, 7/18/11)
      REAL*8           :: T3I(KBLOOP),     TEMP1(KBLOOP)
      REAL*8           :: T3K(KBLOOP),     PRESSK(KBLOOP) 
C FP_ISOP
      ! Variables from Fabien (hotp 7/31/09)
      REAL*8, PARAMETER :: ONE_SIXTY = 1d0 / 60d0
      REAL*8, PARAMETER :: ONE_SEVENTYTHREE = 1d0 / 73d0
      REAL*8, PARAMETER :: CONSTC= 9.871D-07
      REAL*8 :: YA,YB,YC
      REAL*8 :: PTEMP
  

      ! FAST-J: Zero out the dummy array (bmy, 9/30/99)
      DUMMY = 0d0

      ! Zero rate arrays for this column, so that we don't have any
      ! "leftovers" from the last iteration.  This facilitates 
      ! comparing rates directly with the column code. (bmy, 5/22/09)
      RRATE = 0d0
      TRATE = 0d0
      URATE = 0d0
      GLOSS = 0d0
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
            ! Added H2 concentration (bhh, jmao, eam, 7/18/11)
            ! Seasonal variability of H2 may be important, 
            ! but not included in this update (bhh, jmao, eam, 7/18/11)
            CONCH2(KLOOP)     = 0.5000d-6 * DENAIR(KLOOP) 
C
C   Check if sun is up anywhere in this block of grid-boxes.
C   IFSUN gets used in CALCRATE
C   Get the right index for SUNCOS, which is calculated
C   outside of chemistry module.
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
! Add option for non-local PBL mixing (Lin, 03/31/09)

      IF (LNLPBL) THEN

         DO I = 1,NEMIS(NCS)
C get tracer number corresponding to emission species I
            NN = IDEMS(I)
            IF (NN.NE.0) THEN
C find reaction number for emission of tracer NN
               NK = NTEMIS(NN,NCS)
               IF (NK.NE.0) THEN
                  DO KLOOP = 1,KTLOOP
                     JLOOP = LREORDER(KLOOP+JLOOPLO)
                     RRATE(KLOOP,NK) = 0.d0
                     ! Surface emissions of gases are constrained to the 
                     ! lowest model layer, and are considered in the 
                     ! PBL mixing module vdiff_mod, not anymore in SMVGEAR.
                     ! Emissions at higher levels (e.g. aircraft) are still
                     ! managed by SMVGEAR. 
                     !*** As of 05/02/08, REMIS only contains emissions of 
                     ! gases for the SMVGEAR mechanism. (Lin, 05/02/08)
                     !RRATE(KLOOP,NK) = REMIS(JLOOP,I)
                     IZ    = IZSAVE(JLOOP)
                     if (IZ .EQ. 1) then
                        RRATE(KLOOP,NK) = 0.D0
                     else
                        RRATE(KLOOP,NK) = REMIS(JLOOP,I)
                     endif
                  ENDDO
               ENDIF
            ENDIF
         ENDDO

      ELSE

         DO I = 1,NEMIS(NCS)
C get tracer number corresponding to emission species I
            NN = IDEMS(I)
            IF (NN.NE.0) THEN
C find reaction number for emission of tracer NN
               NK = NTEMIS(NN,NCS)
               IF (NK.NE.0) THEN
                  DO KLOOP = 1,KTLOOP
                     RRATE(KLOOP,NK) = 0.d0
                     JLOOP = LREORDER(KLOOP+JLOOPLO)
                     RRATE(KLOOP,NK) = REMIS(JLOOP,I)
                  ENDDO
               ENDIF
            ENDIF
         ENDDO

      ENDIF

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
! Add option for non-local PBL mixing (Lin, 03/31/09)
      IF (LNLPBL) THEN      

         DO I = 1,NUMDEP
            NK = NTDEP(I)
            IF (NK.NE.0) THEN
               DO KLOOP = 1,KTLOOP

                  ! 1-D grid box index (accounts for reordering)
                  JLOOP = LREORDER(KLOOP+JLOOPLO)

                  ! 3-D grid box index
                  IX    = IXSAVE(JLOOP)
                  IY    = IYSAVE(JLOOP)
                  IZ    = IZSAVE(JLOOP)
               
                  ! constrain the dry deposition to the lowest model layer
                  ! (Lin, 04/28/08)
                  if (IZ .NE. 1) then
                     RRATE(KLOOP,NK) = 0.D0 
                  else
                     SELECT CASE ( DEPNAME(I) )
                        CASE ( 'DST1', 'DST2', 'DST3', 'DST4', 'SALA', 
     &                         'SALC' )
                           ! dusts and sea salts
                           RRATE(KLOOP,NK) = DEPSAV(IX,IY,I)
                        CASE DEFAULT
                           ! gases + aerosols for full chemistry
                           RRATE(KLOOP,NK) = 0.D0
                     END SELECT
                  endif
               ENDDO
            ENDIF
         ENDDO

      ELSE
  
         DO I = 1,NUMDEP
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
     &                              GET_FRAC_UNDER_PBLTOP( IX, IY, IZ )

               ENDDO
            ENDIF
         ENDDO

      ENDIF

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
C Add branching ratio for HOC2H4O for C2H4 oxidation (tmf, 12/14/06) 
C
C ---  KF = K*(1-FYHORO(M,T))  ---  HOC2H4O ------> HO2 + 2CH2O 
         DO I          = 1, NNADDF(NCS)
            NK         = NKSPECF( I,NCS )
            DO KLOOP   = 1, KTLOOP
               RRATE(KLOOP,NK) = RRATE(KLOOP,NK) *
     $              (1.D0 - FYHORO(DENAIR(KLOOP), T3K(KLOOP)))
            ENDDO
         ENDDO
C
C ---  KH = K*FYHORO(M,T)  ---  HOC2H4O --O2--> HO2 + GLYC
         DO I          = 1, NNADDH(NCS)
            NK         = NKSPECH( I,NCS )
            DO KLOOP   = 1, KTLOOP
               RRATE(KLOOP,NK) = RRATE(KLOOP,NK) * 
     +              FYHORO(DENAIR(KLOOP), T3K(KLOOP))

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
            ! FP defined CONSTC = 9.871d-7 as a parameter (hotp 7/31/09)
C            CONSTC       = 0.6D0 * 9.871D-07
            DO KLOOP     = 1, KTLOOP
               JLOOP           = LREORDER(JLOOPLO + KLOOP)
               RRATE(KLOOP,NK) = RRATE(KLOOP,NK) *
C     1              (1.D0 + CONSTC*PRESS3(JLOOP))
     1              (1.D0 + 0.6d0*CONSTC*PRESS3(JLOOP))
c new OH+CO rate from JPL2006.
C Watch out! KCO1 and KCO2 have different form!!!!!!!!!!!!!!!(jmao,02/26/09)
               KLO1=5.9D-33*(300*T3I(KLOOP))**(1.4D0) 
               KHI1=1.1D-12*(300*T3I(KLOOP))**(-1.3D0)
               XYRAT1=KLO1*DENAIR(KLOOP)/KHI1
               BLOG1=LOG10(XYRAT1)
               FEXP1=1.D0/(1.D0+BLOG1*BLOG1)
               KCO1=KLO1*DENAIR(KLOOP)*0.6**FEXP1/(1.d0+XYRAT1)
               KLO2=1.5D-13*(300*T3I(KLOOP))**(-0.6D0)
               KHI2=2.1D09 *(300*T3I(KLOOP))**(-6.1D0)
               XYRAT2=KLO2*DENAIR(KLOOP)/KHI2
               BLOG2=LOG10(XYRAT2)
               FEXP2=1.D0/(1.D0+BLOG2*BLOG2)
               KCO2=KLO2*0.6**FEXP2/(1.d0+XYRAT2)
               KCO=KCO1+KCO2
               RRATE(KLOOP,NK)=KCO
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
               ! Changed reaction rate calculation from concO2 to 0.2095d0
               ! (bhh, jmao, eam, 7/18/11)
               RRATE(KLOOP,NK)=RRATE(KLOOP,NK)/
     &              (1.d0+RRATE(KLOOP,NK+1)*0.2095d0)
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
C FP_ISOP 
C Change yield for HAC  

         ! Questions for Fabien: (hotp 7/31/09)
         ! What are HAC, MCO3, and GLYC? Can you provide some comments for
         ! this new code?

         IF (NKHAC(NCS,1) .GT. 0 .AND. NKHAC(NCS,2) .GT. 0) THEN
            DO KLOOP   = 1, KTLOOP
               HAC_FRAC=1d0-23.7d0*EXP(-ONE_SIXTY*T3K(KLOOP)) !FRAC MGLY 

               IF (HAC_FRAC<0d0) HAC_FRAC=0d0

              RRATE(KLOOP,NKHAC(NCS,1))=RRATE(KLOOP,NKHAC(NCS,1))
     &                                  *HAC_FRAC
              RRATE(KLOOP,NKHAC(NCS,2))=RRATE(KLOOP,NKHAC(NCS,2))
     &                                  *(1d0-HAC_FRAC)

            ENDDO
         ENDIF

         IF (NKMCO3(NCS,1) .GT. 0 .AND. NKMCO3(NCS,2) .GT. 0 .AND. 
     &       NKMCO3(NCS,3) .GT. 0) THEN

            DO KLOOP   = 1, KTLOOP
              
              JLOOP           = LREORDER(JLOOPLO + KLOOP)
              PTEMP=CONSTC*PRESS3(JLOOP)

              YA=-4.262d-1*PTEMP**2.925d0+7.42d-1
              YB=1.834d-1*PTEMP**1.396d0+5.069d-3
              YC=1d0-YA-YB

              ! YA,YB, and YC capitalized (hotp 7/31/09)
              IF (YA<0d0) YA=0d0
              IF (YB<0d0) YB=0d0
              IF (YC<0d0) YC=0d0
              IF (YA>1d0) YA=1d0
              IF (YB<0d0) YB=1d0
              IF (YC<0d0) YC=1d0

              !WRITE(*,*) 'P=',PTEMP
              !WRITE(*,*) 'Ya=',Ya
              !WRITE(*,*) 'Yb=',Yb
              !WRITE(*,*) 'Yc=',Yc

              RRATE(KLOOP,NKMCO3(NCS,1))=RRATE(KLOOP,NKMCO3(NCS,1))
     &                                  *YA
              RRATE(KLOOP,NKMCO3(NCS,2))=RRATE(KLOOP,NKMCO3(NCS,2))
     &                                  *YB
              RRATE(KLOOP,NKMCO3(NCS,3))=RRATE(KLOOP,NKMCO3(NCS,2))
     &                                  *YC


            ENDDO
         ENDIF

         IF (NKGLYC(NCS,1) .GT. 0 .AND. NKGLYC(NCS,2) .GT. 0) THEN
            DO KLOOP   = 1, KTLOOP
               GLYC_FRAC=1d0-11.0729d0*EXP(-ONE_SEVENTYTHREE*T3K(KLOOP)) !FRAC HCOOH

              IF (GLYC_FRAC<0d0) GLYC_FRAC=0d0

              RRATE(KLOOP,NKGLYC(NCS,1))=RRATE(KLOOP,NKGLYC(NCS,1))
     &                                  *GLYC_FRAC
              RRATE(KLOOP,NKGLYC(NCS,2))=RRATE(KLOOP,NKGLYC(NCS,2))
     &                                  *(1d0-GLYC_FRAC)

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

         ! Added I-J-L indices to archive diagnostic AD52
         ! jaegle (2/26/09)
         ! I-J-L indices
         IX = IXSAVE(JLOOP)
         IY = IYSAVE(JLOOP)
         IZ = IZSAVE(JLOOP)

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
               ! Initialize DUMMY3 (jaegle, 2/26/09)
               DUMMY3(KLOOP)   = 0d0

               ! Sum up total surface area for all aerosol types
               ! so that we can use it for the planeflight diagnostic
               ! (mje, bmy, 8/7/03)
               IF ( NK == NKN2O5 .or. NK == NKHO2 ) THEN
                  TOTAREA = 0.d0
                  DO N = 1, NDUST + NAER
                     TOTAREA = TOTAREA + TAREA(JLOOP,N)
                  ENDDO
               ENDIF

               ! Loop over sulfate and other aerosols
               DO N = 1, NDUST + NAER

                  ! Surface area of aerosol [cm2 aerosol/cm3 air]
                  XAREA = TAREA(JLOOP,N) 

                  ! Radius for aerosol size bin N (jaegle 2/26/09)
                  XRADIUS = ERADIUS(JLOOP,N) 
                  
                  ! Test if N2O5 hydrolysis rxn
                  IF ( NK == NKN2O5 ) THEN
                     
                     ! Get GAMMA for N2O5 hydrolysis, which is
                     ! a function of aerosol type, temp, and RH
                     XSTKCF = N2O5( N, T3K(KLOOP), ABSHUMK(KLOOP) )

                     ! Archive N2O5 hydrolysis for ND40 diagnostic
                     DUMMY2(KLOOP) = DUMMY2(KLOOP) +
     &                               ( XAREA / TOTAREA * XSTKCF )
                    
                  ! Test if HO2 het uptake reaction
                  ELSE IF ( NK == NKHO2 ) THEN
                     ! Calculate GAMMA for HO2 self-reaction on aerosols, 
                     ! which is a function of aerosol type, radius, 
                     ! temperature, air density, and HO2 concentration 
                     ! (jaegle - 02/26/09)

                     HO2_MOLEC_CM3 = CSPEC(JLOOP,IDHO2)

                     ! Find out whether we are in the continental
                     ! boundary layer and set the CONTINENTAL_PBL
                     ! flag to 1 (also assume that there is no ice/snow)
                     IF (  IZ <= GET_PBL_TOP_L( IX , IY ) .and.
     &                     FRCLND(IX,IY) >= 0.5 .and.
     &                     (.not. IS_ICE(IX,IY) ) ) THEN 
                        CONTINENTAL_PBL=1 
                     ELSE
                        CONTINENTAL_PBL=0
                     ENDIF

                     XSTKCF = HO2( XRADIUS, T3K(KLOOP), XDENA, XSQM,
     &                             HO2_MOLEC_CM3, N , CONTINENTAL_PBL)

                     ! Now call CHECK_VALUE to make sure that XSTKCF is 
                     ! not a NaN or an infinity
                     CALL CHECK_VALUE( XSTKCF, (/KLOOP,0,0,0/),
     &                                 'GAMMA_HO2', 'at calcrate')

                     ! Archive gamma HO2 for ND52 diagnostic
                     DUMMY3(KLOOP) = DUMMY3(KLOOP) +
     &                               ( XAREA / TOTAREA * XSTKCF )


                  ELSE 

                     ! Get GAMMA for species other than N2O5
                     XSTKCF = BRR(NK,NCS)

                  ENDIF

                  ! Reaction rate for dust size bin N
                  RRATE(KLOOP,NK) = RRATE(KLOOP,NK) + 
     $                 ARSL1K(XAREA,XRADIUS,XDENA,XSTKCF,XSTK,XSQM)
               ENDDO
               IF ( ND52 > 0 ) THEN
                  ! Archive gamma HO2 in AD52
                  IF ( IZ <= LD52 ) THEN
                     AD52(IX,IY,IZ) =
     &                    AD52(IX,IY,IZ) + DUMMY3(KLOOP)
                  ENDIF
               ENDIF
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

                  ! For FAST-J, get photorates from fjfunc.f
                  RRATE(KLOOP,NKN)  = FJFUNC(IX,IY,IZ,I,IBRCH,SPECNAME)

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
         ! HARDWIRE the effect of branching ratio of HOC2H4O in EP photolysis
         !   HOC2H4O ------> HO2 + 2CH2O    : marked as I in P column of 
         !                                    'globchem.dat'
         !   HOC2H4O --O2--> HO2 + GLYC     : marked as J in P column of 
         !                                    'globchem.dat'
         !
         ! Add NCS index to NKHOROI and HKHOROJ for SMVGEAR II (tmf, 12/16/06)
         !==============================================================
         IF ( NKHOROI(NCS) > 0 ) THEN

            ! Put J(EP) in correct spot for SMVGEAR II
            PHOTVAL = NKHOROI(NCS) - NRATES(NCS)
            NKN     = NKNPHOTRT(PHOTVAL,NCS)

            DO KLOOP=1,KTLOOP
               RRATE(KLOOP,NKN)=RRATE(KLOOP,NKN) *
     +            ( 1.D0-FYHORO(DENAIR(KLOOP), T3K(KLOOP)) )
            ENDDO
         ENDIF

         IF ( NKHOROJ(NCS) > 0 ) THEN

            ! Put J(EP) in correct spot for SMVGEAR II
            PHOTVAL = NKHOROJ(NCS) - NRATES(NCS)
            NKN     = NKNPHOTRT(PHOTVAL,NCS)

            DO KLOOP=1,KTLOOP
               RRATE(KLOOP,NKN)=RRATE(KLOOP,NKN) *
     +            FYHORO(DENAIR(KLOOP), T3K(KLOOP)) 
            ENDDO
         ENDIF

         !==============================================================
         ! SPECIAL TREATMENT FOR O3+hv -> OH+OH
         ! [O1D]ss=J[O3]/(k[H2O]+k[N2]+k[O2])
         ! SO, THE EFFECTIVE J-VALUE IS J*k[H2O]/(k[H2O]+k[N2]+k[O2])
         !
         ! Add NCS index to NKHNO4 for SMVGEAR II (gcc, bmy, 4/1/03)
         !==============================================================
         ! Now also consider O3 + hv --> HO2 + OH (bhh, jmao, eam, 7/18/11)
         IF ( (NKO3PHOT(NCS) > 0)  .OR. (NKO3PHOTH2(NCS) > 0) ) THEN
            ! Put J(O3) in correct spot for SMVGEAR II
            LJO1DH2 = .false.
            LJO1DH2O = .false.
            IF ( NKO3PHOT(NCS) > 0 ) THEN
              PHOTVAL = NKO3PHOT(NCS) - NRATES(NCS)
              NKNH2O     =  NKNPHOTRT(PHOTVAL,NCS)
              LJO1DH2O = .true.
            ENDIF            
            IF ( NKO3PHOTH2(NCS) > 0 ) THEN
              PHOTVAL = NKO3PHOTH2(NCS) - NRATES(NCS)
              NKNH2     =  NKNPHOTRT(PHOTVAL,NCS)
              LJO1DH2 = .true.
            ENDIF

            DO KLOOP = 1, KTLOOP

               ! Save old value of J-O3 in a diagnostic array 
               ! (gcc, bmy, 4/1/03)
               DUMMY(KLOOP) = RRATE(KLOOP,NKNH2O)

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
c Updated from JPL2006, the difference is pretty small.(jmao,02/26/2009)
C Additional update from JPL 10-6 for reaction of O3 + hv --> HO2 + OH and 
C O1D + H2:
C Includes calculation of k[O3], where
C k[O3]  = J[O3]*1.2e-10/k[O1D], where
C k[O1D] = ([O1D]*[H2]+[O1D]*[H2O])/([O1D]*[H2]+[O1D]*[H2O]+[O1D][N2]+[O1D][O2]) (bhh, jmao, eam, 7/18/11)
               RO1DplH2O = 1.63d-10 * EXP( 60.d0*T3I(KLOOP)) *
     $                     CBLK(KLOOP,IH2O)
               RO1DplH2 = 1.2e-10 * CONCH2(KLOOP)
               RO1D = (RO1DplH2O + RO1DplH2 +
     $             2.15d-11 * EXP(110.d0*T3I(KLOOP)) * CONCN2(KLOOP)   +
     $             3.30d-11 * EXP( 55.d0*T3I(KLOOP)) * CONCO2(KLOOP))
               IF (LJO1DH2O) THEN
                  RRATE(KLOOP,NKNH2O) = RRATE(KLOOP,NKNH2O) *
     $                                 RO1DplH2O / RO1D
               ENDIF
               IF (LJO1DH2) THEN
                  RRATE(KLOOP,NKNH2) = RRATE(KLOOP,NKNH2) *
     $                                  RO1DplH2 / RO1D
               ENDIF
c               RRATE(KLOOP,NKN) = RRATE(KLOOP,NKN) *
c     $            1.45d-10 * EXP( 89.d0*T3I(KLOOP)) * CBLK(KLOOP,IH2O) /
c     $          ( 1.45d-10 * EXP( 89.d0*T3I(KLOOP)) * CBLK(KLOOP,IH2O) +
c     $            2.14d-11 * EXP(110.d0*T3I(KLOOP)) * CONCN2(KLOOP)    +
c     $            3.20d-11 * EXP( 70.d0*T3I(KLOOP)) * CONCO2(KLOOP)    )

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
               CASE ( 'GLYX' )
                  INDEX = 7
               CASE ( 'MGLY' )
                  INDEX = 8
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


      !*********KPP_INTERFACE************ (phs,ks,dhk, 09/15/09)
      IF ( LKPP ) THEN 
         I = 1
         DO NK          = 1, NTRATES(NCS)
            DO KLOOP    = 1, KTLOOP
               IF ( NEWFOLD(NK,NCS) > 0 ) THEN
                  IF ( KLOOP.eq.1 ) THEN
                     IND(I) = NK
                     I      = I +1
                  ENDIF
                  RRATE_FOR_KPP(KLOOP,NK) = RRATE(KLOOP,NEWFOLD(NK,NCS)) ! Saving rates for KPP
               ENDIF
            ENDDO
         ENDDO
      ENDIF
      !********************************************
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

C *********************************************************************

      FUNCTION HO2( RADIUS, TEMP, DENAIR, SQM, HO2DENS,
     &              AEROTYPE, CONTINENTAL_PBL ) RESULT( GAMMA )

      !=================================================================
      ! Internal function HO2 computes the GAMMA reaction probability
      ! for HO2 loss in aerosols based on the recommendation of 
      ! Thornton, Jaegle, and McNeill, 
      ! "Assessing Known Pathways For HO2 Loss in Aqueous Atmospheric
      !  Aerosols: Regional and Global Impacts on Tropospheric Oxidants"
      !  J. Geophys. Res.,  doi:10.1029/2007JD009236, 2008  
      !
      ! gamma(HO2) is a function of aerosol type, radius, temperature
      !
      ! jaegle 01/22/2008
      ! 
      ! Arguments as Input:
      ! ----------------------------------------------------------------
      ! (1 ) RADIUS   (REAL*8 ) : Aerosol radius [cm]
      ! (2 ) TEMP     (REAL*8 ) : Temperature [K]
      ! (3 ) DENAIR   (REAL*8 ) : Air Density [molec/cm3]
      ! (4 ) HO2DENS  (REAL*8 ) : HO2 Number Density [molec/cm3]
      ! (5 ) SQM      (REAL*8 ) : Square root of molecular weight [g/mole]
      ! (6 ) AEROTYPE (INTEGER) : # denoting aerosol type (cf FAST-J)
      ! (7 ) CONTINENTAL_PBL (INTEGER)  : Flag set to 1 if the
      !         box is located in the continenal boundary layer,
      !         otherwise it is zero. Also check for ICE/SNOW (to
      !         disable this at high latitudes)
      !
      ! NOTES:
      !=================================================================
      

      IMPLICIT NONE

      ! Arguments
      REAL*8,  INTENT(IN) :: RADIUS, TEMP, DENAIR, HO2DENS, SQM
      INTEGER, INTENT(IN) :: AEROTYPE, CONTINENTAL_PBL

      ! Local variables
      REAL*8              :: ALPHA
      REAL*8              :: delG, Keq, w, H_eff
      REAL*8              :: A1, B1, k1, k2, A, B, C
      REAL*8              :: kaq, kmt, o2_ss, fluxrxn, DFKG
      REAL*8              :: TEST


      ! Avogadro's number
      REAL*8,  PARAMETER   :: Na = 6.022d23

      ! Ideal gas constant [atm cm3/mol/K], Raq
      REAL*8,  PARAMETER   :: Raq=82.d0

      ! Function return value
      REAL*8              :: GAMMA
      
      !=================================================================
      ! HO2 begins here!
      !=================================================================

      ! Default value
      GAMMA = 0.0d0

      ! Special handling for various aerosols
      SELECT CASE ( AEROTYPE )

         !----------------
         ! Dust 
         !----------------
         CASE ( 1, 2, 3, 4, 5, 6, 7 )      
                                
            ! Assume default gamma=0.1 on dust aerosols
            ! This is tentative as no lab measurements presently exist
            ! for gamma(HO2) on dust aerosols. We assume the rate to
            ! be fast on dust aerosols as transition metal ion induced
            ! chemistry is likely to occur in a thin aqueous surface layer.
            GAMMA = 0.1d0

         !----------------
         ! For Sulfate(8), Black Carbon (9), Organic Carbon (10),
         ! Sea-salt accum & coarse (11,12) calculate the 
         ! reaction probability due to self reaction 
         ! by using the algebraic expression in Thornton et al.  (2008)
         ! (equation 7) which is a function of temperature, aerosol radius,
         ! air density and HO2 concentration. 
         !
         ! Transition metal ions (such as copper and iron) in sea-salt and 
         ! carbonaceous aerosols are complexed to ligands and/or exist at 
         ! a concentration too low to catalyze HO2 loss efficiently, so we 
         ! apply the HO2 self reaction expression directly for these aerosols.
         ! 
         ! In the case of sulfate aerosol, the aerosols likely
         ! contain copper in the continental boundary layer and
         ! HO2 uptake proceeds rapidly. To account for the metal catalyzed
         ! uptake, we assume gamma(HO2)=0.07 (in the mid-range of the recommended
         ! 0.04-0.1 by Thornton et al, based on observed copper concentrations
         ! in the US boundary layer). Outside the continental boundary layer, we
         ! use the HO2-only algebraic expression.
         !
         !----------------
         CASE ( 8, 9, 10, 11, 12)  

            ! Mean molecular speed [cm/s]
            w = 14550.5d0 * sqrt(TEMP/(SQM*SQM))

            ! DFKG = Gas phase diffusion coeff [cm2/s]
            DFKG  = 9.45D17/DENAIR * SQRT(TEMP) * 
     &              SQRT(3.472D-2 + 1.D0/(SQM*SQM))

            !calculate T-dependent solubility and aq. reaction rate constants
            ! hydronium ion concentration
            ! A1 = 1.+(Keq/hplus) 
            ! with Keq = 2.1d-5 [M], Equilibrium constant for 
            ! HO2aq = H+ + O2- (Jacob, 2000)
            !      hplus=10.d0^(-pH), with pH = 5
            ! B1 = Req * TEMP
            ! with Req = 1.987d-3 [kcal/K/mol], Ideal gas constant
            ! Note that we assume a constant pH of 5.
            A1 = 1.+ (2.1d-5 / (10.d0**(-5) ) )
            B1 = 1.987d-3 * TEMP

            ! Free energy change for solvation of HO2 (Hanson 1992, Golden 1991)
            ! in [kcal/mol]:
            ! delG = -4.9-(TEMP-298d0)*delS
            ! with delS=-0.023  [kcal/mol/K],  Entropy change for solvation of HO2
            delG  = -4.9d0 - (TEMP-298.d0) * (-0.023)
            H_eff = exp( -delG / B1 ) * A1

            ! Estimated temp dependent value for HO2 + O2- (k1) and 
            ! HO2+HO2 (see Jacob 1989)
            k1  =   1.58d10 * exp( -3. / B1 )
            k2  =   2.4d9   * exp( -4.7 / B1 )
            kaq = ( k1 * (A1 - 1.d0) + k2) / (A1**2)

            ! Calculate the mass transfer rate constant and s.s. conc. of 
            ! total HO2 in the aqueous phase:
            ! kmt = (RADIUS/DFKG + 4d0/w/alpha)^(-1)
            ! with alpha = mass accomodation coefficient, assumed 
            ! to be 1 (Thornton et al.)
            kmt = 1.d0/( RADIUS/DFKG + 4d0/w/1. )

            !use quadratic formula to obtain [O2-] in particle of radius RADIUS
            A = -2d0 * kaq
            B = -3d0 * kmt / RADIUS / (H_eff * 0.082 * TEMP)
            C =  3d0 * kmt * HO2DENS * 1000d0 / RADIUS / Na

            ! Error check that B^2-(4d0*A*C) is not negative
            TEST= B**2-(4d0*A*C)
            IF ( TEST < 0d0 ) THEN
                GAMMA = 0d0
            ELSE
                ! Calculate the concentration of O2- in the aerosol
                o2_ss= ( -B  -sqrt(B**2-(4d0*A*C)) )/(2d0*A)

                ! Calculate the reactive flux
                fluxrxn = kmt*HO2DENS - o2_ss*Na*kmt/H_eff/Raq/TEMP

                IF ( fluxrxn <= 0d0 ) THEN
                   GAMMA = 0d0
                ELSE
                   ! Gamma for HO2 at TEMP, ho2, and RADIUS given
                   GAMMA = 1./( ( ( HO2DENS/fluxrxn ) - 
     &                            ( RADIUS/DFKG ) ) * w / 4.d0 )
                ENDIF
            ENDIF
            ! For sulfate aerosols, check whether we are in
            ! the continental boundary layer, in which case
            ! copper catalyzed HO2 uptake likely dominates and
            ! speeds up the reaction: we assume gamma=0.07,
            ! which is in the middle of the 0.04-0.1 range recommended
            ! by Thornton et al. (2008)
            !
            IF ( AEROTYPE == 8 .and. CONTINENTAL_PBL == 1) THEN
                GAMMA = 0.07
            ENDIF 

         !----------------
         ! Default
         !----------------
         CASE DEFAULT
            WRITE (6,*) 'Not a suitable aerosol surface '
            WRITE (6,*) 'for HO2 uptake'
            WRITE (6,*) 'AEROSOL TYPE =',AEROTYPE
            CALL GEOS_CHEM_STOP

      END SELECT
     
      ! If negative value is calculated, set it to zero
      IF ( GAMMA  <= 0d0 ) GAMMA = 0d0

      ! Return to CALCRATE
      END FUNCTION HO2
C
C *********************************************************************
C ******************** END OF SUBROUTINE CALCRATE *********************
C *********************************************************************
C
      END SUBROUTINE CALCRATE
