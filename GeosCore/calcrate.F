!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: 
!
! !DESCRIPTION: Subroutine CALCRATE computes reaction rates before passing 
!  them to the SMVGEAR solver.
!\\
!\\
! !INTERFACE:
!
       SUBROUTINE CALCRATE( am_I_Root, Input_Opt, 
     &                      State_Met, State_Chm, RC )
!
! !USES:
!
      USE CMN_SIZE_MOD
      USE COMODE_LOOP_MOD
      USE CMN_DIAG_MOD
      USE CMN_O3_MOD
      USE COMODE_MOD,         ONLY : ABSHUM, AIRDENS, ERADIUS, IXSAVE, 
     &                               IYSAVE, IZSAVE,  JLOP,    PRESS3,  
     &                               REMIS,  T3,      TAREA,    CSPEC
      USE DAO_MOD,            ONLY : IS_ICE
      USE DIAG_MOD,           ONLY : AD22,   LTJV, AD52
      USE DRYDEP_MOD,         ONLY : DEPSAV, NUMDEP, DEPNAME  !SHIPO3DEP
      USE DRYDEP_MOD,         ONLY : NTRAIND
      USE ERROR_MOD,          ONLY : is_safe_div
      USE ERROR_MOD,          ONLY : ERROR_STOP, GEOS_CHEM_STOP
      USE ERROR_MOD,          ONLY : CHECK_VALUE
      USE GCKPP_GLOBAL,       ONLY : IND
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE GIGC_State_Chm_Mod, ONLY : ChmState
      USE GRID_MOD,           ONLY : GET_YMID
      USE PBL_MIX_MOD,        ONLY : GET_FRAC_UNDER_PBLTOP,GET_PBL_TOP_L
      USE PLANEFLIGHT_MOD,    ONLY : ARCHIVE_RXNS_FOR_PF
      USE TRACERID_MOD,       ONLY : IDTHBr,  IDTHOBr
      USE TRACERID_MOD,       ONLY : IDHO2
      USE TRACERID_MOD,       ONLY : IDTSO4, IDTNIT

      ! Also IDHCl and IDHBr for PSCs (SDE 04/17/13)
      USE TRACERID_MOD,       ONLY : IDHCl,   IDHBr,   IDHOCl, IDHOBr
      USE TRACERID_MOD,       ONLY : IDBrNO3, IDClNO3, IDN2O5
      ! Online CH4 (SDE 03/25/13, derived from mpayer 01/24/12)
      USE TRACERID_MOD,       ONLY : IDTCH4
      USE FAST_JX_MOD,        ONLY : FJXFUNC
      USE UCX_MOD,            ONLY : SO4_PHOTFRAC, CHECK_NAT
      USE UCX_MOD,            ONLY : KHETI_SLA
      USE ERROR_MOD,          ONLY : IT_IS_FINITE

      ! To get HEMCO deposition values
      USE HCOI_GC_MAIN_MOD,   ONLY : GetHcoVal

      ! Geos-Chem Flexible Precision
      USE PRECISION_MOD

      IMPLICIT NONE
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU?
      TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
      TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS:
!
      TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT)   :: RC          ! Success or failure
!
! !REVISION HISTORY: 
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
!  (19) Updated O3 photolysis reaction rate and DMS oxidation reaction rates. 
!       See http://wiki.seas.harvard.edu/geos-chem/index.php/
!       Updating_standard_chemistry_with_JPL_10-6. (bhh, jmao, eam, 7/18/11)
!  01 Mar 2012 - R. Yantosca - Now use new functions from grid_mod.F90:
!                              GET_YMID(I,J,L) and GET_AREA_CM2(I,J,L)
!  20 Jun 2012 - R. Yantosca - Under some circumstances it seems that INT_OPE
!                              computed by PARANOX can be negative.  This
!                              would make EMIS_BL and REMIS negative.  If any
!                              negative emissions occur, reset to zero.
!  24 Sep 2012 - R. Yantosca - Bug fix: Dimension QICE and QLIQ with KBLOOP
!                              instead of (IIPAR,JJPAR,LLPAR).  This should
!                              lead to better scalability.
!  09 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!  27 Nov 2012 - R. Yantosca - Replace SUNCOS with State_Met%SUNCOSmid
!  05 Mar 2013 - R. Yantosca - Now get logical, XNMOL fields from Input_Opt
!  05 Mar 2013 - R. Yantosca - Pass am_I_Root, Input_Opt, RC arguments
!  25 Mar 2013 - M. Payer    - Now pass State_Chm via the argument list
!  14 Mar 2013 - M. Payer    - Replace NOx and Ox emissions with NO and O3
!                              emissions as part of removal of NOx-Ox
!                              partitioning
!  16 May 2013 - M. Payer    - Add nitrate effect (L. Zhang)
!  17 May 2013 - M. Payer    - Add improved HO2 uptake (J. Mao)
!  21 Mar 2013 - C. Holmes   - Treat O3 destruction in the PARANOX ship plume 
!                              parameterization as dry deposition. These 
!                              "negative emissions" were previously set to zero.
!                              Now we calculate an effective dry deposition 
!                              frequency that will remove O3 at the same rate 
!                              that the plume parameterization calculates.
!  06 Jun 2013 - M. Payer    - Move ship NOx plume emissions to routine
!                              SETEMIS (C. Holmes)
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!  26 Sep 2013 - R. Yantosca - Renamed GEOS_57 Cpp switch to GEOS_FP
!  25 Apr 2014 - M. Sulprizio- Now store J-values for ACET
!  27 May 2014 - M. Sulprizio- Now store J-values for ALD2, MVK, MACR, HAC, GLYC
!  23 Jul 2014 - R. Yantosca - Remove reference to obsolete CMN_MOD
!  25 Jul 2014 - R. Yantosca - Now inlined functions ARSLIK, FYRNO3, FYHORO
!  25 Jul 2014 - R. Yantosca - Added ProTeX headers
!  13 Oct 2014 - C. Keller   - Added HEMCO deposition values for full-mixing
!                              PBL. This includes O3 deposition from PARANOX.
!                              Further, exclude now SALA and SALC deposition
!                              if non-local PBL scheme is used since those are
!                              already applied in vdiff_mod.F90
!  06 Nov 2014 - M. Yannetti - Added PRECISION_MOD
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! jpp, rate constant for BrNO3 cloud recycling
      REAL(fp) :: CLD_BrNO3_RC

      ! Local variables
      INTEGER :: KLOOP,JLOOP,I,NK,JOLD2,JOLD,NK1,NKN,NH,MLOOP,J
      INTEGER :: NKNH2O, NKNH2               ! (bhh, jmao, eam, 7/18/11)
      LOGICAL :: LJO1DH2, LJO1DH2O, LPCELL   ! (bhh, jmao, eam, 7/18/11)
      INTEGER :: NP,K,IFNC,IBRCH,IX,IY,IZ,IJWINDOW,INDEX,NN,ii,jj
      INTEGER :: PHOTVAL,KSUN
      INTEGER :: NKN_2, PHOTVAL_2

      REAL(fp)  :: ARRNK,FCVNK,FCT1,FCT2,FCT,XYRAT,BLOG,FEXP,RATE3M
      REAL(fp)  :: RO1D, RO1DplH2O, RO1DplH2   ! (bhh, jmao, eam, 7/18/11)
!      REAL(fp)  :: CONSTC,RATE2AIR,RATE3H2O,RIS,RST,TBEGIN,TFINISH
      ! CONSTC declared as a parameter by FP (hotp 7/31/09)
      REAL(fp)  :: RATE2AIR,RATE3H2O,RIS,RST,TBEGIN,TFINISH
      REAL(fp)  :: PBEG,PFIN,PBEGNEW,PFINNEW,TOFDAYB,TOFDAYE,HOURANGB
      REAL(fp)  :: HOURANGE,SINFUNCB,SINFUNCE,XAREA,XRADIUS,XSQM,RRATE2
      REAL(fp)  :: XSTKCF,GMU,DUMMY(KBLOOP),XDENA,XSTK
      REAL(fp)  :: DUMMYO3P(KBLOOP)
      REAL(fp)  :: TK,     CONSEXP, VPRESH2O
      ! New variables for new reations (jmao, 4/20/09)
      REAL(fp)  :: KHI1,KLO1,XYRAT1,BLOG1,FEXP1,KHI2,KLO2,XYRAT2,BLOG2
      REAL(fp)  :: FEXP2,KCO1,KCO2,KCO
      !FP_ISOP
      ! FP variables for isoprene chemistry (hotp 7/31/09)
      REAL(fp)   :: HAC_FRAC      
      REAL(fp)   :: GLYC_FRAC   

      CHARACTER*8      :: SPECNAME

      ! Added for heterogeneous chemistry (bmy, 11/15/01, 8/7/03)
      LOGICAL          :: HETCHEM
      INTEGER          :: N
      REAL(fp)           :: SUMAREA, TOTAREA, DUMMY2(KBLOOP)

      ! Added for HO2 het uptake (jaegle, 2/26/09)
      REAL(fp)           :: HO2_MOLEC_CM3, DUMMY3(KBLOOP)
      INTEGER          :: CONTINENTAL_PBL

      ! For grid-box latitude (bnd, bmy, 7/1/03)
      REAL(fp)           :: YLAT

      ! Variables from "comode.h" which are only ever used in "calcrate.f"
      ! Remove them from "comode.h" and the THREADPRIVATE declarations
      ! (bmy, 7/28/03) 
      REAL(fp)           :: ABSHUMK(KBLOOP), DENAIR(KBLOOP)
      REAL(fp)           :: CONCO2(KBLOOP),  CONCN2(KBLOOP)
      REAL(fp)           :: CONCH2(KBLOOP)     ! (bhh, jmao, eam, 7/18/11)
      REAL(fp)           :: T3I(KBLOOP),     TEMP1(KBLOOP)
      REAL(fp)           :: T3K(KBLOOP),     PRESSK(KBLOOP) 

! FP_ISOP
      ! Variables from Fabien (hotp 7/31/09)
      REAL(fp), PARAMETER :: ONE_SIXTY = 1e+0_fp / 60e+0_fp
      REAL(fp), PARAMETER :: ONE_SEVENTYTHREE = 1e+0_fp / 73e+0_fp
      REAL(fp), PARAMETER :: CONSTC= 9.871e-07_fp
      REAL(fp) :: YA,YB,YC
      REAL(fp) :: PTEMP
  
      !%%% FOR BROMINE CHEMISTRY %%%

      ! True/False variables for treating rates for the
      ! psuedo-reactions that account for bromine heterogeneous
      ! chemistry involving 2 gas phase species. (jpp, 4/12/10)
      LOGICAL :: yn_nkhobr, yn_nk1hbr
      REAL(fp)  :: br_aer1_rc ! temp storage for HBr+HOBr+aerosol rate constant

      REAL(fp)  :: hobr_rtemp, hbr_rtemp
      REAL(fp)  :: ki_hobr, ki_hbr
      REAL(fp)  :: darea ! dummy variable for ice surface area (jpp, 7/5/2011)
      LOGICAL :: yn_div_safe

      ! Local variables for ice and water mixing ratios
      REAL(fp)  :: QLIQ(KBLOOP)
      REAL(fp)  :: QICE(KBLOOP)

      ! Add nitrate effect (lzh, 10/25/2011)
      REAL(fp)           :: TMP1, TMP2

      ! SDE: Strat-chem
      REAL(fp)  :: DUMMYO2(KBLOOP)
      REAL(fp)  :: EDUCTCONC, LIMITCONC
      INTEGER :: PSCIDX
      REAL(fp)  :: PSCEDUCTCONC(11,2)
      LOGICAL :: NATSURFACE
      REAL(fp)  :: ADJUSTEDRATE
      LOGICAL :: KII_KI, PSCBOX, STRATBOX
      REAL(fp), PARAMETER  :: PSCMAXRATE=1.e-13_fp

      ! Minimum lifetime of a species for PSC heterogeneous chemistry
      REAL(fp), PARAMETER  :: PSCMINLIFE=1.e-3_fp
      CHARACTER(LEN=255) :: PSCMESSAGE

      ! Avoid array temporaries in CHECK_VALUE
      INTEGER            :: ERR_LOC(4)
      CHARACTER(LEN=255) :: ERR_VAR
      CHARACTER(LEN=255) :: ERR_MSG

      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      !%%% NEED TO GET THESE FIELDS FROM Input_Opt SO THAT     %%%
      !%%% THE MPI PARALLELIZATION WILL WORK PROPERLY FOR ESMF %%%
      !%%% (bmy, 3/5/13)                                       %%%
      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      LOGICAL :: LKPP
      LOGICAL :: LNLPBL
      LOGICAL :: LVISTAS
      !LOGICAL :: LWARWICK_VSLS
      LOGICAL :: LDICARB
      REAL(fp)  :: XNUMOL(Input_Opt%N_TRACERS)

      ! For HEMCO drydep values (ckeller, 10/13/14)
      INTEGER :: TRCID
      REAL*8  :: Dep 
      LOGICAL :: Fnd

      ! Pointers
      REAL(fp), POINTER :: STT(:,:,:,:)
      REAL(fp), POINTER :: FRCLND(:,:)
      REAL(fp), POINTER :: QI(:,:,:)
      REAL(fp), POINTER :: QL(:,:,:)
      REAL(fp), POINTER :: SUNCOSmid(:,:)
      REAL(fp), POINTER :: T(:,:,:)

      !=================================================================
      ! CALCRATE begins here!
      !=================================================================

      ! Assume success
      RC            = GIGC_SUCCESS

      ! Save fields from Input_Opt into local variables
      LKPP          = Input_Opt%LKPP
      LNLPBL        = Input_Opt%LNLPBL
      LDICARB       = Input_Opt%LDICARB
      XNUMOL        = Input_Opt%XNUMOL

      ! Initialize pointers
      FRCLND       => State_Met%FRCLND
      QI           => State_Met%QI
      QL           => State_Met%QL
      SUNCOSmid    => State_Met%SUNCOSmid
      T            => State_Met%T

      ! Initialize
      darea         = 0e+0_fp

      ! initialize the true/false selectors (jpp, 4/12/10)
      yn_nkhobr     = .false.; yn_nk1hbr = .false.
      br_aer1_rc    = 0.e+0_fp

      ! FAST-JX: Zero out the dummy arrays (bmy, 9/30/99)
      DUMMY   = 0e+0_fp
      DUMMYO3P= 0e+0_fp
      DUMMYO2 = 0e+0_fp

      ! Zero rate arrays for this column, so that we don't have any
      ! "leftovers" from the last iteration.  This facilitates 
      ! comparing rates directly with the column code. (bmy, 5/22/09)
      RRATE         = 0e+0_fp
      TRATE         = 0e+0_fp
      URATE         = 0e+0_fp
      GLOSS         = 0e+0_fp

      ! Zero QICE, QLIQ arrays for safety's sake (bmy, 9/24/12)
      QICE          = 0e+0_fp
      QLIQ          = 0e+0_fp
!
! *********************************************************************
! ************        WRITTEN BY MARK JACOBSON (1993)      ************
! ***             (C) COPYRIGHT, 1993 BY MARK Z. JACOBSON           *** 
! ***       U.S. COPYRIGHT OFFICE REGISTRATION NO. TXu 670-279      *** 
! ***                         (650) 723-6836                        *** 
! *********************************************************************
!
! CCCCCCC     A     L       CCCCCCC  RRRRRRR     A     TTTTTTT  EEEEEEE
! C          A A    L       C        R     R    A A       T     E 
! C         A   A   L       C        RRRRRRR   A   A      T     EEEEEEE
! C        AAAAAAA  L       C        R  R     AAAAAAA     T     E 
! CCCCCCC A       A LLLLLLL CCCCCCC  R    R  A       A    T     EEEEEEE
!
! *********************************************************************
! * THIS SUBROUTINE CALCULATES KINETIC REACTION AND PHOTORATES        *
! * (S-1, CM3 S-1, OR CM6 S-2) AND PRESSURE AND TEMPERATURE-          *
! * DEPENDENCE FOR GAS-PHASE CHEMICAL REACTIONS.                      *
! *                                                                   *
! * HOW TO CALL SUBROUTINE:                                           *
! * ----------------------                                            *
! *  CALL CALCRATE.F FROM PHYSPROC.F WITH                             * 
! *     NCS  = 1..NCSGAS FOR GAS CHEMISTRY                            *
! *********************************************************************
!
! *********************************************************************
! *********************  GAS-PHASE CHEMISTRY  *************************
! *********************************************************************
!
      IF (NCS.LE.NCSGAS) THEN 
!
! *********************************************************************
! *            CALCULATE CONCENTRATIONS OF FIXED SPECIES              *
! *********************************************************************
! AERSURF  = PARTICLE SURFACE AREA (CM2 CM-3)
! KTLOOP   = NUMBER OF GRID-CELLS IN A GRID-BLOCK
! AM       = MOLECULAR WEIGHT OF AIR (28.966 G MOLE-1)
! AVG      = AVOGADRO'S NUMBER (6.02252E+23  # MOLE-1)
! BOLTG    = BOLTZMANN'S CONSTANT (1.38054E-16 ERG K-1) 
! RHO3     = DENSITY OF AIR         (G CM-3) = WTAIR*P(DYN CM-2)/(RSTARG*T) 
! RSTARG   = BOLTG * AVG
! DENAIR   = DENSITY OF AIR         (# CM-3) = P(DYNCM-2) / (T * BOLTG) 
! CONCO2   = OXYGEN CONCENTRATION   (# CM-3)
! CONCN2   = NITROGEN CONCENTRATION (# CM-3)
! RRATE    = RATE CONST (EITHER S-1, CM**3-AIR #-1 S-1, CM**6 #-2 S-1,
!            OR CM**9 #-3 S-1.
! PRESS3   = AIR PRESSURE AT VERTICAL CENTER OF LAYER (MB)
!
! --------------------------- AIR, O2, N2 -----------------------------
!

!   be sure to check out these values!!!
         KSUN=0
         DO 20 KLOOP        = 1, KTLOOP
      
            JLOOP             = LREORDER(JLOOPLO+KLOOP)

            ! Add DENAIR here instead of in physproc.f, so that we
            ! can eliminate the /DKBLOOP/ common block (bmy, 7/28/03)
            DENAIR(KLOOP)     = AIRDENS(JLOOP)
     
            PRESSK(KLOOP)     = PRESS3(JLOOP)
            T3K(KLOOP)        = T3(JLOOP)
            T3I(KLOOP)        = 1.e+0_fp/T3(JLOOP)
            ABSHUMK(KLOOP)    = ABSHUM(JLOOP)
            TEMP1(KLOOP)      = 300.e+0_fp    / T3K(KLOOP)
            CONCO2(KLOOP)     = 0.2095e+0_fp  * DENAIR(KLOOP)
            CONCN2(KLOOP)     = 0.7808e+0_fp  * DENAIR(KLOOP)
            ! Added H2 concentration (bhh, jmao, eam, 7/18/11)
            ! Seasonal variability of H2 may be important, 
            ! but not included in this update (bhh, jmao, eam, 7/18/11)
            CONCH2(KLOOP)     = 0.5000e-6_fp * DENAIR(KLOOP)

            ! Check if sun is up anywhere in this block of grid-boxes.
            ! IFSUN gets used in CALCRATE
            ! Get the right index for SUNCOS, which is calculated
            ! outside of chemistry module
            IX                = IXSAVE(JLOOP)
            IY                = IYSAVE(JLOOP)
            IJWINDOW          = (IY-1)*IIPAR + IX
            !GMU               = SUNCOS(IJWINDOW)
            GMU               = SUNCOSmid(IX,IY)
            IF(GMU .GT. 0.E+0_fp) KSUN = 1

            !==========================================================
            ! Updates for bromine chemistry.  Moved this section here 
            ! from the  DO loop that had been added at the top of 
            ! CALCRATE.  This should improve scalability, and prevent
            ! doing IIPAR*JJPAR*LLPAR operations for each grid block
            ! of KBLOOP grid boxes. (bmy, 9/24/12)
            !==========================================================

            ! Vertical index
            IZ                = IZSAVE(JLOOP)

#if   defined( GEOS_5 ) || defined( MERRA ) || defined( GEOS_FP )
            
            ! GEOS-5 / MERRA / GEOS-5.7 have QI and QL defined as met fields
            ! so use these to define the QICE, QLIQ arrays.  Also, QICE and
            ! QLIQ only need to be dimensioned with KBLOOP instead of
            ! IIPAR,JJPAR,LLPAR. (bmy, 9/24/12)
            QICE(KLOOP)       = QI(IX,IY,IZ)
            QLIQ(KLOOP)       = QL(IX,IY,IZ)

#else

            ! Otherwise, compute QLIQ as a function of temperature ...
            IF ( T(IX,IY,IZ) .LE. 248e+0_fp ) THEN
               QLIQ(KLOOP)    = 0e+0_fp
            ELSE IF ( T(IX,IY,IZ) .GE. 268e+0_fp ) THEN
               QLIQ(KLOOP)    = 1e-6_fp
            ELSE
               QLIQ(KLOOP)    = 1e-6_fp *
     &                      ( ( T(IX,IY,IZ) - 248e+0_fp ) / 20e+0_fp)
            ENDIF

            ! ... and compute QICE from QLIQ (bmy, 9/24/12)
            QICE(KLOOP)       = 1e-6_fp - QLIQ(KLOOP)

#endif
 20      CONTINUE
         IFSUN = 2-KSUN
!
! ---------------------------     H2O     -----------------------------
!
      IF (IH2O.GT.NGAS) THEN
         DO KLOOP      = 1, KTLOOP 
            TK               = T3K(KLOOP) 
            CONSEXP          = 17.2693882E+0_fp * (TK - 273.16E+0_fp) / 
     x           (TK - 35.86E+0_fp) 
            VPRESH2O         = CONSVAP * EXP(CONSEXP) * T3I(KLOOP) 
            CBLK(KLOOP,IH2O) = ABSHUMK(KLOOP) * VPRESH2O  
         ENDDO
      END IF
!
! -----------------  SET O2 TO CONCO2 IF O2 INACTIVE ------------------
!
      IF (IOXYGEN.GT.NGAS) THEN
         DO KLOOP          = 1, KTLOOP
            CBLK(KLOOP,IOXYGEN) = CONCO2(KLOOP)
         ENDDO
      ENDIF
!
! *********************************************************************
! *     INTERANNUALLY-VARYING CH4 CONCENTRATION (bnd, bmy, 7/1/03)    *
! *********************************************************************
!
      ! Test if CH4 is defined as an inert SMVGEAR II species
      IF ( ( IDTCH4 == 0) .and. ( ICH4 > NGAS ) ) THEN

         ! Loop over boxes per grid block
         DO KLOOP = 1, KTLOOP 

            ! 1-D grid box index
            JLOOP = KLOOP + JLOOPLO

            ! Grid-box latitude index
            YLAT  = GET_YMID( IXSAVE(JLOOP), 
     &                        IYSAVE(JLOOP), 
     &                        IZSAVE(JLOOP)  )
            
            ! Pick the CH4 concentration [ppbv] for the proper lat bin
            ! CH4 values are read in "chemdr.f" (outside the parallel loop)
            IF ( YLAT < -30e+0_fp ) THEN
               CBLK(KLOOP,ICH4) = C3090S
            ELSE IF ( YLAT >= -30e+0_fp .and. YLAT < 0e+0_fp  ) THEN
               CBLK(KLOOP,ICH4) = C0030S
            ELSE IF ( YLAT >=   0e+0_fp .and. YLAT < 30e+0_fp ) THEN
               CBLK(KLOOP,ICH4) = C0030N
            ELSE
               CBLK(KLOOP,ICH4) = C3090N
            ENDIF

            ! Convert from [ppbv CH4] to [molec CH4/cm3]
            CBLK(KLOOP,ICH4) = CBLK(KLOOP,ICH4) *1e-9_fp 
     &                       * AIRDENS(JLOOP)
         ENDDO
      ENDIF

      ! ---------------------------------------------------------
      ! Test whether we are setting the VSL bromine source gases
      ! equal to zero or not: (jpp, 4/15/10)
      ! ---------------------------------------------------------
      IF ( ( ICHBr3 > 0 ) .and. ( ICH2Br2 > 0 ) ) THEN

! TODO: testing only
! --- comment for now, not sure if this is important (ckeller) ---
!         IF ( .not. LWARWICK_VSLS ) THEN
!            ! Loop over boxes per grid block
!            DO KLOOP = 1, KTLOOP 
!               CBLK(KLOOP, ICHBr3)  = 0.e+0_fp
!               CBLK(KLOOP, ICH2Br2) = 0.e+0_fp
!            ENDDO
!         ENDIF
      ENDIF

!
! *********************************************************************
! *   CALCULATE KINETIC REACTION RATES USING ARRHENIUS PARAMETERS     * 
! *********************************************************************
! REACTION RATES HAVE THE FORM K = A * (300 / T)**B * EXP(C / T)
!
! NARR = NUMBER OF REACTIONS OF THE FORM K = A
! NABR = NUMBER OF REACTIONS OF THE FORM K = A * (300 / T)**B
! NACR = NUMBER OF REACTIONS OF THE FORM K = A                * EXP(C / T)
! NABC = NUMBER OF REACTIONS OF THE FORM K = A * (300 / T)**B * EXP(C / T)
! NKARR, NKBRR, NKACR, NKABC = REACTION RATE NUMBERS OF EACH
! NARR,  NABR,  NACR,  NABC  REACTION 
!
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
!
         DO 47 I           = 1, NACR(NCS) 
            NK               = NKACR(I,NCS)
            ARRNK            = ARR(NK,NCS)  
            DO 45 KLOOP      = 1, KTLOOP
               RRATE(KLOOP,NK) = ARRNK * EXP(KCRR(NK,NCS) / T3K(KLOOP)) 
 45         CONTINUE
 47      CONTINUE
!
         DO 52 I           = 1, NABC(NCS) 
            NK               = NKABC(I,NCS)
            ARRNK            = ARR(NK,NCS)  
            DO 50 KLOOP      = 1, KTLOOP
               RRATE(KLOOP,NK) = ARRNK * TEMP1(KLOOP)**BRR(NK,NCS)
     1              * EXP(KCRR(NK,NCS) / T3K(KLOOP))
 50         CONTINUE
 52      CONTINUE

!
! *********************************************************************
!              Photorates for Harvard Geos Code
! *********************************************************************
! PRATE           = PHOTORATE (SEC-1) IF SUN IS DOWN, PRATE = 0.
! RRATE           = RATE COEFFICIENT (SEC-1)
! NRATES          = NUMBER OF KINETIC REACTION RATES.
!
!
      IF(IFSUN.EQ.1) THEN
         DO I                = 1, JPHOTRAT(NCS)
            NK               = NRATES(NCS) + I
            NKN              = NKNPHOTRT(I,NCS)
            SPECNAME         = NAMEGAS(IRM(1,NK,NCS))
            IFNC             = DEFPRAT(NK,NCS) + 0.01e+0_fp
            IBRCH            = 10.E+0_fp*(DEFPRAT(NK,NCS)-IFNC) + 0.5E+0_fp

            DO KLOOP            = 1, KTLOOP 
               JLOOP            = LREORDER(KLOOP+JLOOPLO)

               ! Translate 1-D to 3-D grid box indices
               IX               = IXSAVE(JLOOP)
               IY               = IYSAVE(JLOOP)
               IZ               = IZSAVE(JLOOP)                  

               ! Get cosine( SZA ) using 1-D array index
               IJWINDOW         = (IY-1)*IIPAR + IX
               GMU              = SUNCOSmid(IX,IY)

               ! For daylight boxes...
               IF(GMU.GT. 0.E+0_fp) THEN

                  ! For FAST-JX, get photorates from fjfunc.f
                  RRATE(KLOOP,NKN)  = FJXFUNC(IX,IY,IZ,I,IBRCH,SPECNAME)

!### Debug: warn if there are negative J-values, for either 
!### FAST-J or SLOW-J photolysis (bmy, 10/1/98)
!###                  IF ( RRATE(KLOOP,NK) < 0 ) THEN
!###                     PRINT*, 'CALCRATE.F: J-Value < 0: ', 
!###     &                  IX, IY, IZ, IBRCH, SPECNAME, KLOOP, NK,
!###     &                  RRATE(KLOOP,NK)
!###                  ENDIF

               ELSE

                  ! Nighttime: photorates are zero
                  RRATE(KLOOP,NKN)  = 0.E+0_fp

               ENDIF
            ENDDO
         ENDDO
      ENDIF

!
! *********************************************************************
! ********  MULTIPLY RATES BY CONSTANT SPECIES CONCENTRATIONS  ********
! *       (EITHER M, O2, N2, OR ANY ACTIVE OR INACTIVE SPECIES)       *   
! *********************************************************************
! NMAIR    = # REACTIONS WHERE THE SPECIES IN THE THIRD POSITION IS 
!              IS 'M' = 'O2 + N2'
! NMO2     = # REACTIONS WHERE THE SPECIES IN THE THIRD POSITION IS O2 
! NMN2     = # REACTIONS WHERE THE SPECIES IN THE THIRD POSITION IS N2 
! NMOTH    = # OCCURENCES OF SPECIES IN THIRD POSITION THAN ARE NOT
!              O2, N2, OR M, OR OF SPECIES IN ANY POSITION THAT ARE
!              INACTIVE.
! LGASBINO = JGAS (SET IN READCHEM AND GASCONC)
!
         DO 72 I           = 1, NMAIR(NCS)   
            NK               = NREACAIR(I,NCS)
            DO 70 KLOOP      = 1, KTLOOP 
               RRATE(KLOOP,NK) = RRATE(KLOOP,NK) * DENAIR(KLOOP) 
 70         CONTINUE
 72      CONTINUE
!
         DO 82 I           = 1, NMO2(NCS) 
            NK               = NREACO2(I,NCS)
            DO 80 KLOOP      = 1, KTLOOP 
               RRATE(KLOOP,NK) = RRATE(KLOOP,NK) * CONCO2(KLOOP) 
 80         CONTINUE
 82      CONTINUE
!        
         DO 92 I           = 1, NMN2(NCS)  
            NK               = NREACN2(I,NCS)
            DO 90 KLOOP      = 1, KTLOOP 
               RRATE(KLOOP,NK) = RRATE(KLOOP,NK) * CONCN2(KLOOP) 
 90         CONTINUE
 92      CONTINUE
!
! *********************************************************************
! *                   PRESSURE-DEPENDENT EFFECTS                      * 
! * ADD THE THIRD BODY EFFECT FOR PRESSURE DEPENDENCE OF RATE         *
! * COEFFICIENTS. THE REACTIONS WERE READ IN IN PAIRS (LOW AND HIGH   *
! * PRESSURE LIMITS) WITH THE SPECIFIC INDICATOR, 'P'  IN THE INPUT   *
! * DATA SET. SEE DEMORE ET AL. (1990) JPL 90-1 FOR MORE DETAILS      *
! *********************************************************************
! NPRESM  = # PRESSURE DEPENDENT 3-BODY REACTIONS  
! FCV     = CHARACTERIZES FALLOFF CURVE (SEE ATKINSON ET. AL (1992)
!           J. PHYS. CHEM. REF. DATA 21, P. 1145). USUALLY = 0.6 
!           HOWEVER, TWO TEMPERATURE-DEPENDENT EXPRESSIONS ARE ALSO USED: 
!             FCV = EXP(-T/FCT1)  OR 
!             FCV = EXP(-T/FCT1)+EXP(-FCT2/T)
! RATE(NK)   = K(0,T)[M], WHERE K(0,T) = 3-BODY, LOW PRESSURE LIMIT COEF. 
! RATE(NK+1) = K(INF,T) = 2-BODY, HIGH PRESSURE LIMIT COEF. 
!
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
                  FEXP           = 1.e+0_fp / (1.e+0_fp + BLOG * BLOG)
                  RRATE(KLOOP,NK)= RRATE(KLOOP,NK)*FCT**FEXP
     &                           /(1e+0_fp+XYRAT) 
 150           CONTINUE
            ELSEIF (FCT1.NE.0.) THEN 
               DO 155 KLOOP     = 1, KTLOOP 
                  FCT            = EXP(-T3K(KLOOP) / FCT1)
                  XYRAT          = RRATE(KLOOP,NK) / RRATE(KLOOP,NK+1) 
                  BLOG           = LOG10(XYRAT)
                  FEXP           = 1.e+0_fp / (1.e+0_fp + BLOG * BLOG)
                  RRATE(KLOOP,NK)= RRATE(KLOOP,NK)*FCT**FEXP
     &                           /(1e+0_fp+XYRAT)  
 155           CONTINUE 
            ELSE
               DO 160 KLOOP     = 1, KTLOOP 
                  XYRAT          = RRATE(KLOOP,NK) / RRATE(KLOOP,NK+1) 
                  BLOG           = LOG10(XYRAT)
                  FEXP           = 1.e+0_fp / (1.e+0_fp + BLOG * BLOG)
                 RRATE(KLOOP,NK)=RRATE(KLOOP,NK)*FCVNK**FEXP
     &                          /(1e+0_fp+XYRAT)  
 160           CONTINUE
            ENDIF
 165     CONTINUE
!
! *********************************************************************
! * SET THE RATES OF ALL THERMALLY DISSOCIATING SPECIES. SEE DEMORE   *
! * ET AL. (1990). CHEMICAL KINETICS AND PHOTOCHEMICAL DATA FOR USE   *
! * IN STRATOSPHERIC MODELING. JPL. 90-1, P. 93. THE RATE HAS THE     *
! * FORM Kf / [A EXP (C / T)], WHERE Kf IS THE REACTION IN THE        *
! * REVERSE DIRECTION.                                                * 
! *********************************************************************
! NNEQ      = # THERMALLY DISSOCIATING EQUILIBRIUM REACTIONS. PREVIOUS
!             EQUATION MUST BE PRESSURE-DEPENDENT. 
! RATE(NK1) = CM3 MOLEC.-1 S-1 (BIMOLECULAR RATE FROM PRESSURE-DEPEND)
! RATE(NK)  = CM3 MOLEC.-1 (EQUILIBRIUM CONSTANT) (BEFORE CALCULATION)
! RATE(NK)  = S-1 (UNIMOLECULAR RATE AFTER CALCULATION) 
!
         DO 182 I           = 1, NNEQ(NCS) 
            NK                = NREACEQ(I,NCS)
            NK1               = NREQOTH(I,NCS)
            DO 180 KLOOP      = 1, KTLOOP
               RRATE(KLOOP,NK)  = RRATE(KLOOP,NK1) / RRATE(KLOOP,NK)
 180        CONTINUE
 182     CONTINUE 
!
! *********************************************************************
!     MULTIPLY RATE COEFFICIENT BY OTHER INACTIVE CONCENTRATIONS
! *********************************************************************
! THIS LOOP MUST OCCUR AFTER EQUILIBRIUM REACTIONS 
!
         DO 183 I          = 1, NMOTH(NCS) 
            NK               = NREACOTH(I,NCS)
            JOLD             = LGASBINO(I,NCS)
            DO 181 KLOOP     = 1, KTLOOP 
               RRATE(KLOOP,NK) = RRATE(KLOOP,NK) * CBLK(KLOOP,JOLD) 
 181        CONTINUE
 183     CONTINUE
!
! *********************************************************************
! *                       SET SPECIAL RATES                           *
! *********************************************************************
!

! FGAP  25/02/2011 
! add carbon dependence of RO2+HO2
         DO I  = 1, NNRO2HO2(NCS)
            NK         = NKSPECRO2HO2( I,NCS )
            DO KLOOP   = 1, KTLOOP
               RRATE(KLOOP,NK) = RRATE(KLOOP,NK) * 
     &              FCRO2HO2(RRATE(KLOOP,NK+1))
           ENDDO
         ENDDO


! ---  K = K1 + K2  ---- 
         IF (NKSPECW(NCS).GT.0) THEN
            NK         = NKSPECW( I )
            DO KLOOP   = 1, KTLOOP
               RRATE(KLOOP,NK) = RRATE(KLOOP,NK) + RRATE(KLOOP,NK+1)
            ENDDO
         ENDIF
!
! ---  K = K1*FYRNO3(K2,M,T)  ---   addition branch of RO2+NO
         DO I          = 1, NNADDA(NCS)
            NK         = NKSPECA( I,NCS )
            DO KLOOP   = 1, KTLOOP
               RRATE(KLOOP,NK) = RRATE(KLOOP,NK) * 
     +              FYRNO3(RRATE(KLOOP,NK+1),DENAIR(KLOOP),T3K(KLOOP))
            ENDDO
         ENDDO
!
! ---  K = K1*(1-FYRNO3(K2,M,T))  ---  abstraction branch of RO2+NO
         DO I          = 1, NNADDB(NCS)
            NK         = NKSPECB( I,NCS )
            DO KLOOP   = 1, KTLOOP
               RRATE(KLOOP,NK) = RRATE(KLOOP,NK) *
     $              (1.E+0_fp - FYRNO3(RRATE(KLOOP,NK+1), DENAIR(KLOOP),
     $              T3K(KLOOP)))
            ENDDO
         ENDDO
!
! ---  K = K1*([O2]+3.5D18)/(2*[O2]+3.5D18) --- HO2+2*CO branch of GLYX+OH/NO3
         DO I          = 1, NNADDC(NCS)
            NK         = NKSPECC( I,NCS )
            DO KLOOP   = 1, KTLOOP
               RRATE(KLOOP,NK) = RRATE(KLOOP,NK) * 
     +              (CONCO2(KLOOP)+3.5E+18_fp)/(2.E+0_fp*CONCO2(KLOOP)
     +              +3.5E+18_fp)
            ENDDO
         ENDDO
!
! ---  K = K1*[O2]/(2*[O2]+3.5D18)  --- GLCO3 branch of GLYX+OH/NO3
         DO I          = 1, NNADDD(NCS)
            NK         = NKSPECD( I,NCS )
            DO KLOOP   = 1, KTLOOP
               RRATE(KLOOP,NK) = RRATE(KLOOP,NK) * 
     +              (CONCO2(KLOOP))/(2.E+0*CONCO2(KLOOP)+3.5E+18)
            ENDDO
         ENDDO
! Add branching ratio for HOC2H4O for C2H4 oxidation (tmf, 12/14/06) 
!
! ---  KF = K*(1-FYHORO(M,T))  ---  HOC2H4O ------> HO2 + 2CH2O 
         DO I          = 1, NNADDF(NCS)
            NK         = NKSPECF( I,NCS )
            DO KLOOP   = 1, KTLOOP
               RRATE(KLOOP,NK) = RRATE(KLOOP,NK) *
     $              (1.E+0_fp - FYHORO(DENAIR(KLOOP), T3K(KLOOP)))
            ENDDO
         ENDDO
!
! ---  KH = K*FYHORO(M,T)  ---  HOC2H4O --O2--> HO2 + GLYC
         DO I          = 1, NNADDH(NCS)
            NK         = NKSPECH( I,NCS )
            DO KLOOP   = 1, KTLOOP
               RRATE(KLOOP,NK) = RRATE(KLOOP,NK) * 
     +              FYHORO(DENAIR(KLOOP), T3K(KLOOP))

            ENDDO
         ENDDO
!
! ---  OH + HNO3:   K = K0 + K3[M] / (1 + K3[M]/K2)  ------
         IF (NKSPECX(NCS).GT.0) THEN
            NK               = NKSPECX(NCS)
            DO KLOOP     = 1, KTLOOP
            RRATE2=RRATE(KLOOP,NK+2)*DENAIR(KLOOP)
            RRATE(KLOOP,NK) = RRATE(KLOOP,NK) + RRATE2 /
     1           (1.E+0_fp + RRATE2 / RRATE(KLOOP,NK+1))    
            ENDDO
         ENDIF
!
! ---    OH + CO: K = K0(1+0.6 Patm)  ------------ 
!    CONSTC includes a factor to convert PRESS3 from (dyn cm-2) to (atm)
         IF (NKSPECY(NCS).GT.0) THEN
            NK           = NKSPECY(NCS)
            ! FP defined CONSTC = 9.871d-7 as a parameter (hotp 7/31/09)
!            CONSTC       = 0.6D0 * 9.871D-07
            DO KLOOP     = 1, KTLOOP
               JLOOP           = LREORDER(JLOOPLO + KLOOP)
               RRATE(KLOOP,NK) = RRATE(KLOOP,NK) *
!     1              (1.D0 + CONSTC*PRESS3(JLOOP))
     1              (1.E+0_fp + 0.6e+0_fp*CONSTC*PRESS3(JLOOP))
c new OH+CO rate from JPL2006.
! Watch out! KCO1 and KCO2 have different form!!!!!!!!!!!!!!!(jmao,02/26/09)
               KLO1=5.9E-33_fp*(300*T3I(KLOOP))**(1.4E+0_fp) 
               KHI1=1.1E-12_fp*(300*T3I(KLOOP))**(-1.3E0_fp)
               XYRAT1=KLO1*DENAIR(KLOOP)/KHI1
               BLOG1=LOG10(XYRAT1)
               FEXP1=1.E+0_fp/(1.E+0_fp+BLOG1*BLOG1)
               KCO1=KLO1*DENAIR(KLOOP)*0.6**FEXP1/(1.e+0_fp+XYRAT1)
               KLO2=1.5E-13_fp*(300*T3I(KLOOP))**(-0.6E+0_fp)
               KHI2=2.1e+09_fp *(300*T3I(KLOOP))**(-6.1E+0_fp)
               XYRAT2=KLO2*DENAIR(KLOOP)/KHI2
               BLOG2=LOG10(XYRAT2)
               FEXP2=1.E+0_fp/(1.E+0_fp+BLOG2*BLOG2)
               KCO2=KLO2*0.6**FEXP2/(1.e+0_fp+XYRAT2)
               KCO=KCO1+KCO2
               RRATE(KLOOP,NK)=KCO
            ENDDO
         ENDIF
!
! ---    MCO3+MO2:  K = K1 / (1+K2)   ------------
!  temperature-dependent branching ratio
         DO I          = 1,NNADDV(NCS)
            NK         = NKSPECV( I,NCS )
            DO KLOOP   = 1, KTLOOP
               RRATE(KLOOP,NK)=RRATE(KLOOP,NK)/(1.e+0_fp
     &                        +RRATE(KLOOP,NK+1))
            ENDDO
         ENDDO
!
         ! Add special reaction for DMS + OH + O2 (bdf, bmy, 4/18/03)
         DO I          = 1,NNADDG(NCS)
            NK         = NKSPECG( I,NCS )
            DO KLOOP   = 1, KTLOOP
               ! Changed reaction rate calculation from concO2 to 0.2095e+0_fp
               ! (bhh, jmao, eam, 7/18/11)
               RRATE(KLOOP,NK)=RRATE(KLOOP,NK)/
     &              (1.e+0_fp+RRATE(KLOOP,NK+1)*0.2095e+0_fp)
            ! SMVGEARII doesn't have structure to multiply rate(nk+1) by [O2]
            ENDDO
         ENDDO
!
! ---  HO2/NO3 + HO2:  K = (K1 + K2)*(1+1.4E-21*[H2O]*EXP(2200/T))  --- 
!  dependence of HO2/NO3 + HO2 on water vapor
         IF (NKSPECZ(NCS).GT.0) THEN
            NK         = NKSPECZ(NCS)
            DO KLOOP   = 1, KTLOOP
#if defined( UCX )
               JLOOP = LREORDER(JLOOPLO+KLOOP)
               RRATE(KLOOP,NK) =
     +          (RRATE(KLOOP,NK)+RRATE(KLOOP,NK+1)*DENAIR(KLOOP)) * 
     +          (1.D0+1.4E-21_fp*CSPEC(JLOOP,IH2O)*EXP(2200.E+0_fp
     +          /T3K(KLOOP)))
#else
               RRATE(KLOOP,NK) =
     +          (RRATE(KLOOP,NK)+RRATE(KLOOP,NK+1)*DENAIR(KLOOP)) * 
     +          (1.D0+1.4E-21_fp*CBLK(KLOOP,IH2O)*EXP(2200.e+0_fp
     +          /T3K(KLOOP)))
#endif
            ENDDO
         ENDIF
! FP_ISOP 
! Change yield for HAC  

         ! Questions for Fabien: (hotp 7/31/09)
         ! What are HAC, MCO3, and GLYC? Can you provide some comments for
         ! this new code?

         IF (NKHAC(NCS,1) .GT. 0 .AND. NKHAC(NCS,2) .GT. 0) THEN
            DO KLOOP   = 1, KTLOOP
               HAC_FRAC=1e+0_fp-23.7e+0_fp*EXP(-ONE_SIXTY*T3K(KLOOP)) !FRAC MGLY 

               IF (HAC_FRAC<0e+0_fp) HAC_FRAC=0e+0_fp

              RRATE(KLOOP,NKHAC(NCS,1))=RRATE(KLOOP,NKHAC(NCS,1))
     &                                  *HAC_FRAC
              RRATE(KLOOP,NKHAC(NCS,2))=RRATE(KLOOP,NKHAC(NCS,2))
     &                                  *(1e+0_fp-HAC_FRAC)

            ENDDO
         ENDIF

         IF (NKMCO3(NCS,1) .GT. 0 .AND. NKMCO3(NCS,2) .GT. 0 .AND. 
     &       NKMCO3(NCS,3) .GT. 0) THEN

            DO KLOOP   = 1, KTLOOP
              
              JLOOP           = LREORDER(JLOOPLO + KLOOP)
              PTEMP=CONSTC*PRESS3(JLOOP)

              YA=-4.262e-1_fp*PTEMP**2.925e+0_fp+7.42e-1_fp
              YB=1.834e-1_fp*PTEMP**1.396e+0_fp+5.069e-3_fp
              YC=1e+0_fp-YA-YB

              ! YA,YB, and YC capitalized (hotp 7/31/09)
              IF (YA<0e+0_fp) YA=0e+0_fp
              IF (YB<0e+0_fp) YB=0e+0_fp
              IF (YC<0e+0_fp) YC=0e+0_fp
              IF (YA>1e+0_fp) YA=1e+0_fp
              IF (YB<0e+0_fp) YB=1e+0_fp
              IF (YC<0e+0_fp) YC=1e+0_fp

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
               GLYC_FRAC=1e+0_fp-11.0729e+0_fp*EXP(-ONE_SEVENTYTHREE
     &                  *T3K(KLOOP)) !FRAC HCOOH

              IF (GLYC_FRAC<0e+0_fp) GLYC_FRAC=0e+0_fp

              RRATE(KLOOP,NKGLYC(NCS,1))=RRATE(KLOOP,NKGLYC(NCS,1))
     &                                  *GLYC_FRAC
              RRATE(KLOOP,NKGLYC(NCS,2))=RRATE(KLOOP,NKGLYC(NCS,2))
     &                                  *(1e+0_fp-GLYC_FRAC)

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

         ! Initialize, for safety's sake
         PSCBOX     = .FALSE.
         STRATBOX   = .FALSE.
         NATSURFACE = .FALSE.

         IF ( HETCHEM ) THEN

            !===========================================================
            ! Perform heterogeneous chemistry on sulfate aerosol
            ! plus each of the NDUST dust size bins from FAST-J
            !===========================================================
            XDENA   = DENAIR(KLOOP)
            XSTK    = SQRT(T3K(KLOOP))
                  
#if defined( UCX )
            ! Check surface type of PSCs (SDE 04/17/13)
            CALL CHECK_NAT( IX, IY, IZ, NATSURFACE, PSCBOX, STRATBOX,
     &                      Input_Opt,  State_Met,  State_Chm )
#endif

            DO I       = 1, NNADDK(NCS)
               NK      = NKSPECK(I,NCS)
               XSQM    = SQRT(ARR(NK,NCS))

               ! Initialize
               RRATE(KLOOP,NK) = 0e+0_fp
               DUMMY2(KLOOP)   = 0e+0_fp
               ! Initialize DUMMY3 (jaegle, 2/26/09)
               DUMMY3(KLOOP)   = 0e+0_fp

               ! Sum up total surface area for all aerosol types
               ! so that we can use it for the planeflight diagnostic
               ! (mje, bmy, 8/7/03)
               IF ( NK == NKN2O5 .or. NK == NKHO2 ) THEN
                  TOTAREA = 0.e+0_fp
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

                  ! Divide by educt concentration
                  KII_KI = .FALSE.

                  ! Test if N2O5 hydrolysis rxn
                  IF ( NK == NKN2O5 ) THEN 

                     KII_KI = .TRUE.

                     ! Get GAMMA for N2O5 hydrolysis, which is
                     ! a function of aerosol type, temp, and RH
                     IF (N.eq.14) THEN
                        IF (NATSURFACE) THEN
                           XSTKCF = 4.0e-4_fp ! NAT
                        ELSE
                           XSTKCF = 0.02e+0_fp ! Ice
                        ENDIF
                     ELSEIF (N.eq.13) THEN
                        ! Stratospheric aerosol
                        XSTKCF = KHETI_SLA(IX,IY,IZ,1)
                     ELSE
                        ! In UCX, ABSHUMK will have been set by
                        ! STT(I,J,L,IDTH2O)
                        XSTKCF = N2O5( N, T3K(KLOOP), ABSHUMK(KLOOP) )
                     ENDIF

                     ! Nitrate effect; reduce the gamma on nitrate by a
                     ! factor of 10 (lzh, 10/25/2011)
                     IF ( N == 8 ) THEN
                        TMP1 = State_Chm%Tracers(IX,IY,IZ,IDTSO4) +
     &                         State_Chm%Tracers(IX,IY,IZ,IDTNIT)
                        TMP2 = State_Chm%tracers(IX,IY,IZ,IDTNIT)
                        IF ( TMP1 .GT. 0.0 ) THEN
                           XSTKCF = XSTKCF * ( 1.0e+0_fp - 0.9e+0_fp
     &                            *TMP2/TMP1 )
                        ENDIF
                     ENDIF

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
     &                     (.not. IS_ICE( IX,IY,State_Met) ) ) THEN 

                        CONTINENTAL_PBL=1 
                     ELSE
                        CONTINENTAL_PBL=0
                     ENDIF

                     IF (N.gt.12) THEN
                        XSTKCF = TINY(1e+0_fp)
                     ELSE
                        XSTKCF = HO2( XRADIUS, T3K(KLOOP),    XDENA,
     &                                XSQM,    HO2_MOLEC_CM3, N,
     &                                CONTINENTAL_PBL, Input_Opt )
                     ENDIF

                     ! Now call CHECK_VALUE to make sure that XSTKCF is 
                     ! not a NaN or an infinity
                     ERR_VAR = 'GAMMA_HO2'
                     ERR_MSG = 'at calcrate'
                     ERR_LOC = (/ KLOOP,0,0,0 /)
                     CALL CHECK_VALUE( XSTKCF,  ERR_LOC, 
     &                                 ERR_VAR, ERR_MSG  )

                     ! Archive gamma HO2 for ND52 diagnostic
                     DUMMY3(KLOOP) = DUMMY3(KLOOP) +
     &                               ( XAREA / TOTAREA * XSTKCF )


                  ELSE IF ( NK == NKBrNO3 ) THEN

                     ! Only apply adjustment if at high altitude
                     KII_KI = STRATBOX

                     ! Get the aerosol type
                     ! If it's sulfate then use 0.8 for alpha, following
                     !  JPL 2006 kinetics evaluation... holds for many
                     !  temperatures and percent weights of sulfate.
                     ! If not, then use the IUPAC recommendation of
                     !  0.3, which is an input in globchem.dat
                     ! (jpp, 5/4/10)
                     IF ( N == 8 ) THEN
                        ! sulfate aerosol
                        XSTKCF = 0.8e+0_fp
                     ELSE IF ( (N == 11) .OR. ( N == 12) ) THEN
                        ! 2 modes of sea-salt
                        XSTKCF = BRR(NK,NCS)
                     ELSE IF ( N == 13 ) THEN
                        ! SSA/STS
                        XSTKCF = KHETI_SLA(IX,IY,IZ,6)
                     ELSE IF ( N == 14 ) THEN 
                        ! Ice/NAT PSC
                        IF (NATSURFACE) THEN 
                           XSTKCF = 0.001e+0_fp
                        ELSE
                           XSTKCF = 0.3e+0_fp
                        ENDIF
                     ELSE
                        XSTKCF = 0e+0_fp
                     ENDIF

                  ELSE IF ( NK == NKHOBr .OR. NK == NK1HBr) THEN
                     ! jpp, 3/22/11: set the sticking coefficient to 
                     !  ~0 for aerosol types we don't want reactions on
                     !  for the HBr and HOBr surface reaction

                     KII_KI = STRATBOX

                     ! Select proper aerosol type
                     IF ( (N == 8) .OR. (N == 11) .OR. (N == 12)) THEN
                        ! sulfate, 2 modes of sea-salt
                        XSTKCF = BRR(NK,NCS)
                     ELSEIF ( N == 13 ) THEN
                        XSTKCF = KHETI_SLA(IX,IY,IZ,11)
                     ELSEIF ( N == 14 ) THEN
                        XSTKCF = 0.1e+0_fp
                     ELSE
                        XSTKCF = 0e+0_fp
                     ENDIF

                  ! Only consider PSC reactions in strat
                  ELSE IF (STRATBOX) THEN
                     IF ( NK == NKPSC(NCS,2 ) ) THEN
                        ! N2O5(g)    + HCl(l,s)
                        KII_KI = .TRUE.
                        IF (N.eq.8) THEN
                           XSTKCF = 0.1e-4_fp ! Sulfate
                        ELSEIF (N.eq.13) THEN
                           XSTKCF = KHETI_SLA(IX,IY,IZ,2)
                        ELSEIF (N.eq.14) THEN
                           IF (NATSURFACE) THEN
                              XSTKCF = 0.003e+0_fp ! NAT
                           ELSE
                              XSTKCF = 0.03e+0_fp ! Ice
                           ENDIF
                        ELSE
                           XSTKCF = 0e+0_fp
                        ENDIF
                     ELSE IF ( NK == NKPSC(NCS,3 ) ) THEN
                        ! ClNO3(g)   + H2O(l,s)
                        KII_KI = .TRUE.

                        IF (N.eq.8) THEN
                           XSTKCF = 0.1e-3_fp ! Sulfate
                        ELSEIF (N.eq.13) THEN
                           XSTKCF = KHETI_SLA(IX,IY,IZ,3)
                        ELSEIF (N.eq.14) THEN
                           IF (NATSURFACE) THEN
                              XSTKCF = 0.004e+0_fp ! NAT
                           ELSE
                              XSTKCF = 0.3e+0_fp ! Ice
                           ENDIF
                        ELSE
                           XSTKCF = 0e+0_fp
                        ENDIF
                     ELSE IF ( NK == NKPSC(NCS,4 ) ) THEN
                        ! ClNO3(g)   + HCl(l,s)
                        KII_KI = .TRUE.
                        IF (N.eq.8) THEN
                           XSTKCF = 0.1e-4_fp ! Sulfate
                        ELSEIF (N.eq.13) THEN
                           XSTKCF = KHETI_SLA(IX,IY,IZ,4)
                        ELSEIF (N.eq.14) THEN
                           IF (NATSURFACE) THEN
                              XSTKCF = 0.2e+0_fp ! NAT
                           ELSE
                              XSTKCF = 0.3e+0_fp ! Ice
                           ENDIF
                        ELSE
                           XSTKCF = 0e+0_fp
                        ENDIF
                     ELSE IF ( NK == NKPSC(NCS,5 ) ) THEN
                        ! ClNO3(g)   + HBr(l,s)
                        KII_KI = .TRUE.
                        IF (N.eq.8) THEN
                           XSTKCF = 0.e+0_fp ! Sulfate
                        ELSEIF (N.eq.13) THEN
                           XSTKCF = KHETI_SLA(IX,IY,IZ,5)
                        ELSEIF (N.eq.14) THEN
                           IF (NATSURFACE) THEN
                              XSTKCF = 0.3e+0_fp ! NAT
                           ELSE
                              XSTKCF = 0.3e+0_fp ! Ice
                           ENDIF
                        ELSE
                           XSTKCF = 0e+0_fp
                        ENDIF
                     ELSE IF ( NK == NKPSC(NCS,7 ) ) THEN
                        ! BrNO3(g)   + HCl(l,s)
                        KII_KI = .TRUE.
                        IF (N.eq.8) THEN
                           XSTKCF = 0.9e+0_fp ! Sulfate
                        ELSEIF (N.eq.13) THEN
                           XSTKCF = KHETI_SLA(IX,IY,IZ,7)
                        ELSEIF (N.eq.14) THEN
                           IF (NATSURFACE) THEN
                              XSTKCF = 0.3e+0_fp ! NAT
                           ELSE
                              XSTKCF = 0.3e+0_fp ! Ice
                           ENDIF
                        ELSE
                           XSTKCF = 0e+0_fp
                        ENDIF
                     ELSE IF ( NK == NKPSC(NCS,8 ) ) THEN
                        ! HOCl(g)    + HCl(l,s)
                        KII_KI = .TRUE.
                        IF (N.eq.8) THEN
                           XSTKCF = 0.8e+0_fp ! Sulfate
                        ELSEIF (N.eq.13) THEN
                           XSTKCF = KHETI_SLA(IX,IY,IZ,8)
                        ELSEIF (N.eq.14) THEN
                           IF (NATSURFACE) THEN
                              XSTKCF = 0.1e+0_fp ! NAT
                           ELSE
                              XSTKCF = 0.2e+0_fp ! Ice
                           ENDIF
                        ELSE
                           XSTKCF = 0e+0_fp
                        ENDIF
                     ELSE IF ( NK == NKPSC(NCS,9 ) ) THEN
                        ! HOCl(g)    + HBr(l,s)
                        KII_KI = .TRUE.
                        IF (N.eq.8) THEN
                           XSTKCF = 0.8e+0_fp ! Sulfate
                        ELSEIF (N.eq.13) THEN
                           XSTKCF = KHETI_SLA(IX,IY,IZ,9)
                        ELSEIF (N.eq.14) THEN
                           IF (NATSURFACE) THEN
                              XSTKCF = 0.3e+0_fp ! NAT
                           ELSE
                              XSTKCF = 0.3e+0_fp ! Ice
                           ENDIF
                        ELSE
                           XSTKCF = 0e+0_fp
                        ENDIF
                     ELSE IF ( NK == NKPSC(NCS,10) ) THEN
                        ! HOBr(g)    + HCl(l,s)
                        KII_KI = .TRUE.
                        IF (N.eq.8) THEN
                           XSTKCF = 0.8e+0_fp ! Sulfate
                        ELSEIF (N.eq.13) THEN
                           XSTKCF = KHETI_SLA(IX,IY,IZ,10)
                        ELSEIF (N.eq.14) THEN
                           IF (NATSURFACE) THEN
                              XSTKCF = 0.1e+0_fp ! NAT
                           ELSE
                              XSTKCF = 0.3e+0_fp ! Ice
                           ENDIF
                        ELSE
                           XSTKCF = 0e+0_fp
                        ENDIF
                     ENDIF
                  ELSE

                     ! Get GAMMA for species other than N2O5
                     XSTKCF = BRR(NK,NCS)

                  ENDIF

                  IF (N.eq.13) THEN
                     ! Calculate for stratospheric liquid aerosol
                     ! Note that XSTKCF is actually a premultiplying
                     ! factor in this case, including c-bar
                     ADJUSTEDRATE = XAREA * XSTKCF
                  ELSE
                     ! Reaction rate for surface of aerosol
                     ADJUSTEDRATE=
     &                 ARSL1K(XAREA,XRADIUS,XDENA,XSTKCF,XSTK,XSQM)
                  ENDIF

                  IF (KII_KI .and. N.gt.12) THEN
                     ! PSC reaction - prevent excessive reaction rate
                     IF (ADJUSTEDRATE.gt.(1.e+0_fp/PSCMINLIFE)) THEN
                        ADJUSTEDRATE = 1.e+0_fp/PSCMINLIFE
                     ENDIF
                  ENDIF

                  ! Add to overall reaction rate
                  RRATE(KLOOP,NK) = RRATE(KLOOP,NK) + ADJUSTEDRATE

               ENDDO

               IF ( ND52 > 0 ) THEN
                  ! Archive gamma HO2 in AD52
                  IF ( IZ <= LD52 ) THEN
                     AD52(IX,IY,IZ) =
     &                    AD52(IX,IY,IZ) + DUMMY3(KLOOP)
                  ENDIF
               ENDIF
            ENDDO

            ! ---------------------------------------------------
            ! Setting the HOBr and HBr pseudo-rxns equal to
            ! whichever rate is smaller/limiting. (jpp, 4/12/10)
            !
            ! Make sure that both rates have been calculated,
            ! then, take the lowest rate and apply it to
            ! both pseudo-rxns (rate limiting). (jpp, 4/12/10)
            ! -------------------------------------------------

            ! Initialize GEOS-Chem tracer array [kg] from Chemistry State object
            ! (mpayer, 12/6/12)
            STT => State_Chm%Tracers

            IF ( ( NK1HBr > 0) .and. ( NKHOBr  > 0 ) ) THEN

               ! select the min of the two rates
               hbr_rtemp  = RRATE(KLOOP, NK1HBr) *
     &                      STT(IX,IY,IZ, IDTHBr)
               hobr_rtemp = RRATE(KLOOP, NKHOBr) * 
     &                      STT(IX,IY,IZ, IDTHOBr)

               ! ---------------------------------------------
               ! Kludging the rates to be equal to one another
               ! to avoid having to keep setting equality in
               ! SMVGEAR solver. (jpp, 5/10/2011)
               ! ---------------------------------------------
               IF ( hbr_rtemp > hobr_rtemp ) THEN

                  ! 1. is it safe to divide?
                  yn_div_safe = is_safe_div( 
     &                          RRATE(KLOOP, NKHOBr) *
     &                          STT(IX,IY,IZ, IDTHOBr), 
     &                          STT(IX,IY,IZ, IDTHBr) )

                  IF (yn_div_safe) THEN
                     ! 2. if it is safe, then go ahead
                     RRATE(KLOOP,NK1HBr) = RRATE(KLOOP, NKHOBr)   *
     &                                     STT(IX,IY,IZ, IDTHOBr) /
     &                                     STT(IX,IY,IZ, IDTHBr)
                  ELSE
                     ! if not, then set rates really small...
                     ! b/c the largest contributor is very small.
                     RRATE(KLOOP, NK1HBr) = TINY(1.e+0_fp)
                     RRATE(KLOOP, NKHOBr) = TINY(1.e+0_fp)
                  ENDIF

               ! if HOBr rate is larger than HBr rate
               ELSE

                  ! 1. is it safe to divide?
                  yn_div_safe = is_safe_div( 
     &                          RRATE(KLOOP, NK1HBr) *
     &                          STT(IX,IY,IZ, IDTHBr), 
     &                          STT(IX,IY,IZ, IDTHOBr) )

                  IF (yn_div_safe) THEN
                     ! 2. if it is safe, then go ahead
                     RRATE(KLOOP,NKHOBr) = RRATE(KLOOP, NK1HBr)  *
     &                                     STT(IX,IY,IZ, IDTHBr) /
     &                                     STT(IX,IY,IZ, IDTHOBr)
                  ELSE
                     ! if not, then set rates really small...
                     ! b/c the largest contributor is very small.
                     RRATE(KLOOP, NK1HBr) = TINY(1.e+0_fp)
                     RRATE(KLOOP, NKHOBr) = TINY(1.e+0_fp)
                  ENDIF
               ENDIF
            ENDIF

            ! ----------------------------------------------
            !  Add rate for cloud heterogeneous chemistry
            !  to RRATE(KLOOP, NKBrNO3). (jpp, 2/28/2011)
            ! ----------------------------------------------
            IF ( NKBrNO3 > 0 ) THEN
               IF (.not.PSCBOX) THEN
                  cld_brno3_rc = CLD1K_BrNO3(IX,          IY,
     &                                       IZ,          XDENA,
     &                                       QLIQ(KLOOP), State_Met )

                  RRATE(KLOOP, NKBrNO3) = RRATE(KLOOP, NKBrNO3) + 
     &                                    cld_brno3_rc
               ENDIF
            ENDIF

            ! ----------------------------------------------
            !  Calculate rates for HOBr + HBr + ice --> Br2
            !  for cold and mixed clouds. (jpp, 6/16/2011)
            ! ----------------------------------------------
            IF ( ( NK2HBr > 0 ) .and. ( NK2HOBr > 0 ) ) THEN
               IF (.not.PSCBOX) THEN
                  CALL cldice_hbrhobr_rxn( IX, IY, IZ, 
     &                                     XDENA,
     &                                     QICE(KLOOP),
     &                                     STT(IX,IY,IZ,IDTHBr),
     &                                     STT(IX,IY,IZ,IDTHOBr),
     &                                     ki_hbr,
     &                                     ki_hobr,
     &                                     DAREA,
     &                                     State_Met )
               ELSE
                  ! For PSCs, het chem already accounted for in
                  ! aerosol code
                  ki_hbr = 0e+0_fp
                  ki_hobr = 0e+0_fp
               ENDIF

               RRATE(KLOOP, NK2HBr)  = ki_hbr
               RRATE(KLOOP, NK2HOBr) = ki_hobr
            ENDIF

            ! Free pointer
            NULLIFY( STT )

#if defined( UCX )
            ! SDE 05/30/13: Loop through PSC het. chem. and limit rates
            ! to prevent solver failure

            ! Set PSC educt concentrations (SDE 04/24/13)
            PSCEDUCTCONC( 1,1) = CSPEC(JLOOP,IDN2O5)
            PSCEDUCTCONC( 1,2) = CSPEC(JLOOP,IH2O)

            PSCEDUCTCONC( 2,1) = CSPEC(JLOOP,IDN2O5)
            PSCEDUCTCONC( 2,2) = CSPEC(JLOOP,IDHCl)
             
            PSCEDUCTCONC( 3,1) = CSPEC(JLOOP,IDClNO3)
            PSCEDUCTCONC( 3,2) = CSPEC(JLOOP,IH2O)
            
            PSCEDUCTCONC( 4,1) = CSPEC(JLOOP,IDClNO3)
            PSCEDUCTCONC( 4,2) = CSPEC(JLOOP,IDHCl)
            
            PSCEDUCTCONC( 5,1) = CSPEC(JLOOP,IDClNO3)
            PSCEDUCTCONC( 5,2) = CSPEC(JLOOP,IDHBr)
            
            PSCEDUCTCONC( 6,1) = CSPEC(JLOOP,IDBrNO3)
            PSCEDUCTCONC( 6,2) = CSPEC(JLOOP,IH2O)
            
            PSCEDUCTCONC( 7,1) = CSPEC(JLOOP,IDBrNO3)
            PSCEDUCTCONC( 7,2) = CSPEC(JLOOP,IDHCl)
            
            PSCEDUCTCONC( 8,1) = CSPEC(JLOOP,IDHOCl)
            PSCEDUCTCONC( 8,2) = CSPEC(JLOOP,IDHCl)
            
            PSCEDUCTCONC( 9,1) = CSPEC(JLOOP,IDHOCl)
            PSCEDUCTCONC( 9,2) = CSPEC(JLOOP,IDHBr)
            
            PSCEDUCTCONC(10,1) = CSPEC(JLOOP,IDHOBr)
            PSCEDUCTCONC(10,2) = CSPEC(JLOOP,IDHCl)

            ! This is still pseudo-first-order - ignore            
            PSCEDUCTCONC(11,1) = CSPEC(JLOOP,IDHOBr)
            PSCEDUCTCONC(11,2) = CSPEC(JLOOP,IDHBr)

            DO PSCIDX=1,10
               NK = NKPSC(NCS,PSCIDX)
               ! Pseudo-first-order reactions - divide by number-conc
               ! of aerosol-phase educt to yield 2nd-order constant
               IF (NK.gt.0) THEN
                  EDUCTCONC = PSCEDUCTCONC(PSCIDX,2)
                  LIMITCONC = PSCEDUCTCONC(PSCIDX,1)
                  ADJUSTEDRATE = RRATE(KLOOP,NK)

                  ! ---SAFETY-CHECK REACTION---
                  ! Definition of 2nd order reaction rate:
                  ! k[A][B] = -d[A]/dt = -d[B]/dt
                  ! 
                  ! However, here we are using a pseudo-first order
                  ! reaction rate, ki, and assuming that [B] is
                  ! abundant. To get k, we will therefore perform:
                  ! k = ki/[B]
                  !
                  ! This will yield the following when solved:
                  ! -d[A]/dt = ki[A] = -d[B]/dt
                  !
                  ! This has some problems, especially for small [B]!
                  ! To get around this, we run the following checks:
                  !
                  ! 1. The lifetime of [A] is 1/ki. If this is below
                  !    PSCMINLIFE, limit reaction rate to yield the
                  !    specified lifetimedepletion (ki = 1/60)
                  ! 2. The depletion time of [B] is [B]/(ki[A]). If
                  !    this is below PSCMINLIFE, limit reaction rate
                  !    (ki = [B]/(T*[A])
                  ! 3. If [B] is < 100 molec/cm3, or ki/[B] yields
                  !    a Nan, force k = 0.
                  !
                  ! If all these checks are passed, we set k = ki/[B].
                  ! Rxn 11 is first-order - ignore
                  IF (PSCIDX.eq.1) THEN
                     ! Convert from 1st-order to 2nd-order
                     IF (IS_SAFE_DIV(EDUCTCONC,LIMITCONC)) THEN
                        ! Temporarily store [B]/(T*[A])
                        LIMITCONC = EDUCTCONC/(PSCMINLIFE*LIMITCONC)
                        IF (ADJUSTEDRATE.gt.LIMITCONC) THEN
                           ADJUSTEDRATE = LIMITCONC
                        ENDIF
                     ELSE
                        ADJUSTEDRATE = 0e+0_fp
                     ENDIF
                     IF ((EDUCTCONC.gt.1.d2).and.
     &                      (IS_SAFE_DIV(ADJUSTEDRATE,EDUCTCONC))) THEN
                        ADJUSTEDRATE = ADJUSTEDRATE/EDUCTCONC
                     ELSE
                        ADJUSTEDRATE = 0e+0_fp
                     ENDIF
                  ELSEIF (PSCIDX.ne.11) THEN
                     IF (ADJUSTEDRATE.gt.(1.e+0_fp/PSCMINLIFE)) THEN
                        ADJUSTEDRATE = 1.e+0_fp/PSCMINLIFE
                     ENDIF
                     ! Convert from 1st-order to 2nd-order
                     IF (IS_SAFE_DIV(EDUCTCONC,LIMITCONC)) THEN
                        ! Temporarily store [B]/(T*[A])
                        LIMITCONC = EDUCTCONC/(PSCMINLIFE*LIMITCONC)
                        IF (ADJUSTEDRATE.gt.LIMITCONC) THEN
                           ADJUSTEDRATE = LIMITCONC
                        ENDIF
                     ELSE
                        ADJUSTEDRATE = 0e+0_fp
                     ENDIF
                     IF ((EDUCTCONC.gt.1.e+2_fp).and.
     &                      (IS_SAFE_DIV(ADJUSTEDRATE,EDUCTCONC))) THEN
                        ADJUSTEDRATE = ADJUSTEDRATE/EDUCTCONC
                     ELSE
                        ADJUSTEDRATE = 0e+0_fp
                     ENDIF
                  ENDIF
                  RRATE(KLOOP,NK) = ADJUSTEDRATE
               ENDIF
            ENDDO
#endif

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
!     ENDIF NCS.EQ.1 OR 2

      ! SDE 05/14/13: Moved initial photorate calculation to be earlier
      ! so that O1D/O3P/H/N steady-states can be found
      ! Special photorates
      IF (IFSUN.EQ.1) THEN

         !==============================================================
         ! HARDWIRE the effect of branching ratio of HOC2H4O in EP photolysis
         !   HOC2H4O ------> HO2 + 2CH2O    : marked as I in P column of 
         !                                    'globchem.dat'
         !   HOC2H4O --O2--> HO2 + GLYC     : marked as J in P column of 
         !                                    'globchem.dat'
         !
         ! Add NCS index to NKHOROI and HKHOROJ for SMVGEARII (tmf, 12/16/06)
         !==============================================================

         ! Not yet modified this for compatibility with Fast-JX v7.0.
         ! (SDE 04/01/13)

         ! Now only do the following if using the dicarbonyls mechanism
         ! (sde, mps, 5/28/14)
         IF ( LDICARB ) THEN
            IF ( NKHOROI(NCS) > 0 ) THEN

               ! Put J(EP) in correct spot for SMVGEAR II
               PHOTVAL = NKHOROI(NCS) - NRATES(NCS)
               NKN     = NKNPHOTRT(PHOTVAL,NCS)

               DO KLOOP=1,KTLOOP
                  RRATE(KLOOP,NKN)=RRATE(KLOOP,NKN) *
     +               ( 1.E+0_fp-FYHORO(DENAIR(KLOOP), T3K(KLOOP)) )
               ENDDO
            ENDIF

            IF ( NKHOROJ(NCS) > 0 ) THEN

               ! Put J(EP) in correct spot for SMVGEAR II
               PHOTVAL = NKHOROJ(NCS) - NRATES(NCS)
               NKN     = NKNPHOTRT(PHOTVAL,NCS)

               DO KLOOP=1,KTLOOP
                  RRATE(KLOOP,NKN)=RRATE(KLOOP,NKN) *
     +               FYHORO(DENAIR(KLOOP), T3K(KLOOP)) 
               ENDDO
            ENDIF
         ENDIF

         !==============================================================
         ! SPECIAL TREATMENT FOR H2SO4+hv -> SO2 + 2OH
         !
         ! Only allow photolysis of H2SO4 when gaseous (SDE 04/11/13)
         !
         !==============================================================

#if defined( UCX )       
         IF ( NKSO4PHOT(NCS) > 0 ) THEN
            ! Calculate if H2SO4 expected to be gaseous or aqueous
            ! Only allow photolysis above 6 hPa
            PHOTVAL = NKSO4PHOT(NCS)-NRATES(NCS)
            NKN     = NKNPHOTRT(PHOTVAL,NCS)
            
            DO KLOOP = 1, KTLOOP 
               JLOOP            = LREORDER(KLOOP+JLOOPLO)

               ! Find location
               IX               = IXSAVE(JLOOP)
               IY               = IYSAVE(JLOOP)
               IZ               = IZSAVE(JLOOP)                  
               RRATE(KLOOP,NKN)=RRATE(KLOOP,NKN)
     &                            *SO4_PHOTFRAC(IX,IY,IZ)
            ENDDO
         ENDIF
#endif

         !==============================================================
         ! SPECIAL TREATMENT FOR O2+hv -> O+O
         !
         ! O2 is an inactive species, so we need to premultiply RRATE
         ! by the O2 concentration prior to sending to the solver.
         ! (SDE 04/05/13)
         !==============================================================
         IF ( NKO2PHOT(NCS) > 0) THEN
         
            !GET CORRECT NKN LOCATION OF O2 PHOT
            PHOTVAL = NKO2PHOT(NCS)-NRATES(NCS)
            NKN     = NKNPHOTRT(PHOTVAL,NCS)

            !Loop through boxes and multiply by O2 concentration
            DO KLOOP = 1, KTLOOP
               ! Store actual rate in dummy array         
               DUMMYO2(KLOOP) = RRATE(KLOOP,NKN)

               RRATE(KLOOP,NKN) = RRATE(KLOOP,NKN)*CONCO2(KLOOP)
            ENDDO
         ENDIF

         !==============================================================
         ! SPECIAL TREATMENT FOR O3+hv -> OH+OH (trop-only simulation)
         !                    or O3+hv -> O+O2  (UCX simulation)
         !
         ! [O1D]ss=J[O3]/(k[H2O]+k[N2]+k[O2])
         ! SO, THE EFFECTIVE J-VALUE IS J*k[H2O]/(k[H2O]+k[N2]+k[O2])
         !
         ! Add NCS index to NKHNO4 for SMVGEAR II (gcc, bmy, 4/1/03)
         !
         ! We don't want to do this if strat-chem is in use, as all
         ! the intermediate reactions are included - this would be
         ! double-counting (SDE 04/01/13)
         !==============================================================
         ! Now also consider O3 + hv --> HO2 + OH (bhh, jmao, eam, 7/18/11)
         IF ( (NKO3PHOT(NCS) > 0)  .OR. (NKO3PHOTH2(NCS) > 0) ) THEN
            ! Put J(O3) in correct spot for SMVGEAR II
            LJO1DH2     = .false.
            LJO1DH2O    = .false.
            IF ( NKO3PHOT(NCS) > 0 ) THEN
               PHOTVAL  = NKO3PHOT(NCS) - NRATES(NCS)
               NKNH2O   = NKNPHOTRT(PHOTVAL,NCS)
               LJO1DH2O = .true.
            ENDIF
            IF ( NKO3PHOTH2(NCS) > 0 ) THEN
               PHOTVAL  = NKO3PHOTH2(NCS) - NRATES(NCS)
               NKNH2    = NKNPHOTRT(PHOTVAL,NCS)
               LJO1DH2  = .true.
            ENDIF

#if defined( UCX )
            IF ( NKO3PPHOT(NCS) > 0 ) THEN
               PHOTVAL_2 = NKO3PPHOT(NCS) - NRATES(NCS)
               NKN_2     = NKNPHOTRT(PHOTVAL_2,NCS)
            ENDIF
#endif

            DO KLOOP = 1, KTLOOP

#if defined( UCX )
               ! Need to subtract O3->O1D from rate
               RRATE(KLOOP,NKN_2)=RRATE(KLOOP,NKN_2)-RRATE(KLOOP,NKNH2O)
               DUMMYO3P(KLOOP) = RRATE(KLOOP,NKN_2)
#endif

               ! Save old value of J-O3 in a diagnostic array 
               ! (gcc, bmy, 4/1/03)
               DUMMY(KLOOP) = RRATE(KLOOP,NKNH2O)

#if !defined( UCX )
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
               ! Updated from JPL2006, the difference is pretty small.
               ! (jmao,02/26/2009)
               !========================================================
               ! Additional update from JPL 10-6 for reaction of
               ! O3 + hv --> HO2 + OH and O1D + H2:
               ! Includes calculation of k[O3], where
               ! k[O3]  = J[O3]*1.2e-10/k[O1D], where
               ! k[O1D] = ([O1D]*[H2]+[O1D]*[H2O])/
               !          ([O1D]*[H2]+[O1D]*[H2O]+[O1D][N2]+[O1D][O2])
               ! (bhh, jmao, eam, 7/18/11)
               !========================================================
               RO1DplH2O = 1.63e-10_fp * EXP( 60.e+0_fp*T3I(KLOOP)) *
     $                     CBLK(KLOOP,IH2O)
               RO1DplH2  = 1.2e-10 * CONCH2(KLOOP)
               RO1D      = (RO1DplH2O + RO1DplH2 +
     $             2.15e-11_fp * EXP(110.e+0_fp*T3I(KLOOP)) 
     $             * CONCN2(KLOOP)   +
     $             3.30e-11_fp * EXP( 55.e+0_fp*T3I(KLOOP)) 
     $             * CONCO2(KLOOP))
               IF (LJO1DH2O) THEN
                  RRATE(KLOOP,NKNH2O) = RRATE(KLOOP,NKNH2O) *
     $                                  RO1DplH2O / RO1D
               ENDIF
               IF (LJO1DH2) THEN
                  RRATE(KLOOP,NKNH2)  = RRATE(KLOOP,NKNH2) *
     $                                  RO1DplH2 / RO1D
               ENDIF
!               RRATE(KLOOP,NKN) = RRATE(KLOOP,NKN) *
!     $            1.45d-10 * EXP( 89.e+0_fp*T3I(KLOOP)) * CBLK(KLOOP,IH2O) /
!     $          ( 1.45d-10 * EXP( 89.e+0_fp*T3I(KLOOP)) * CBLK(KLOOP,IH2O) +
!     $            2.14d-11 * EXP(110.e+0_fp*T3I(KLOOP)) * CONCN2(KLOOP)    +
!     $            3.20d-11 * EXP( 70.e+0_fp*T3I(KLOOP)) * CONCO2(KLOOP)    )
#endif

            ENDDO
         ENDIF
      ELSEIF(IFSUN.EQ.2) THEN
         DO I          = 1, JPHOTRAT(NCS)
            NKN        = NKNPHOTRT(I,NCS)
            DO KLOOP   = 1, KTLOOP
               RRATE(KLOOP,NKN)  = 0.E+0_fp
            ENDDO
         ENDDO
      ELSE
         ! ERROR IN IFSUN
         CALL ERROR_STOP( 'ERROR in IFSUN -- STOP 0345', 'calcrate.f' )
      ENDIF

      !=================================================================
      ! Store J-values for J(NO2) and J(O1D) for ship emissions 
      !=================================================================
      DO I  = 1, JPHOTRAT(NCS)
         NK  = NRATES(NCS) + I
         NKN = NKNPHOTRT(I,NCS)

         ! Name of species being photolyzed
         SPECNAME = NAMEGAS(IRM(1,NK,NCS))

         SELECT CASE ( TRIM( SPECNAME ) )
            CASE ( 'NO2'   )
               INDEX = 1
            CASE ( 'HNO3'  )
               INDEX = 2
            CASE ( 'H2O2'  )
               INDEX = 3
            CASE ( 'CH2O'  )
               INDEX = 4
            CASE ( 'O3'    )
               INDEX = 5
            !  INDEX = 6 for POH (trop-only) or JO3_O3P (UCX)
            CASE ( 'GLYX'  )
               INDEX = 7
            CASE ( 'MGLY'  )
               INDEX = 8
            CASE ( 'BrO'   )
               INDEX = 9
            CASE ( 'HOBr'  )
               INDEX = 10
            CASE ( 'BrNO2' )
               INDEX = 11
            CASE ( 'BrNO3' )
               INDEX = 12
            CASE ( 'CHBr3' )
               INDEX = 13
            CASE ( 'Br2'   )
               INDEX = 14
            CASE ( 'O2'    )
               INDEX = 15
            CASE ( 'N2O'   )
               INDEX = 16
            CASE ( 'NO'    )
               INDEX = 17
            CASE ( 'NO3'   )
               INDEX = 18
            CASE ( 'CFC11' )
               INDEX = 19
            CASE ( 'CFC12' )
               INDEX = 20
            CASE ( 'CCl4'  )
               INDEX = 21
            CASE ( 'CH3Cl' )
               INDEX = 22
            CASE ( 'ACET'  )
               INDEX = 23
            CASE ( 'ALD2'  )
               INDEX = 24
            CASE ( 'MVK'  )
               INDEX = 25
            CASE ( 'MACR'  )
               INDEX = 26
            CASE ( 'HAC'  )
               INDEX = 27
            CASE ( 'GLYC'  )
               INDEX = 28
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
             
            !===========================================================
            ! Store J-values in ND22 diagnostic 
            !===========================================================
            IF ( ND22 > 0 ) THEN
               ! Save J-values for 2PM diagnostic boxes
               ! Use AD22 array for J-value diagnostic (bmy, 9/30/99)
               IF ( LTJV(IX,IY) > 0 .and. IZ <= LD22 ) THEN
                  IF ( INDEX == 5 ) THEN
#if !defined( UCX )
                     ! Store unadjusted total J-O3 as AD22(:,:,:,5)
                     AD22(IX,IY,IZ,5) = AD22(IX,IY,IZ,5) +  
     &                                     DUMMY(KLOOP)
                     ! Store POH as AD22(:,:,:,6)
                     AD22(IX,IY,IZ,6) = AD22(IX,IY,IZ,6) +
     &                                     RRATE(KLOOP,NKN)
#else
                     ! Store O1D/O3P rates separately
                     AD22(IX,IY,IZ,5) = AD22(IX,IY,IZ,5) + 
     &                                     DUMMY(KLOOP)
                     AD22(IX,IY,IZ,6) = AD22(IX,IY,IZ,6) +
     &                                     DUMMYO3P(KLOOP)
#endif
                  ELSEIF ( INDEX == 15 ) THEN
                     ! Store unadjusted O2
                     AD22(IX,IY,IZ,15) = AD22(IX,IY,IZ,15) +
     &                                     DUMMYO2(KLOOP)
                  ELSE
                     ! Store other J-Values in their appropriate slots
                     AD22(IX,IY,IZ,INDEX) = AD22(IX,IY,IZ,INDEX) +
     &                                        RRATE(KLOOP,NKN)
                  ENDIF
               ENDIF
            ENDIF

         ENDDO
      ENDDO

!
! *********************************************************************
! ******                   SET EMISSION RATES                    ******
! *********************************************************************
!
! Add option for non-local PBL mixing (Lin, 03/31/09)

      IF (LNLPBL) THEN

         DO I = 1,NEMIS(NCS)
! get tracer number corresponding to emission species I
            NN = IDEMS(I)
            IF (NN.NE.0) THEN
! find reaction number for emission of tracer NN
               NK = NTEMIS(NN,NCS)
               IF (NK.NE.0) THEN
                  DO KLOOP = 1,KTLOOP
                     JLOOP = LREORDER(KLOOP+JLOOPLO)
                     RRATE(KLOOP,NK) = 0.e+0_fp
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
                        RRATE(KLOOP,NK) = 0.E+0_fp
                     else
                        RRATE(KLOOP,NK) = REMIS(JLOOP,I)
                     endif
                  ENDDO
               ENDIF
            ENDIF
         ENDDO

      ELSE

         DO I = 1,NEMIS(NCS)
! get tracer number corresponding to emission species I
            NN = IDEMS(I)
            IF (NN.NE.0) THEN
! find reaction number for emission of tracer NN
               NK = NTEMIS(NN,NCS)
               IF (NK.NE.0) THEN
                  DO KLOOP = 1,KTLOOP
                     RRATE(KLOOP,NK) = 0.e+0_fp
                     JLOOP = LREORDER(KLOOP+JLOOPLO)
                     RRATE(KLOOP,NK) = REMIS(JLOOP,I)
                  ENDDO
               ENDIF
            ENDIF
         ENDDO

      ENDIF
!
! *********************************************************************
! ******                SET DRY DEPOSITION RATES                 ******
! ******                                                         ******
! ******   NOTE: Now compute drydep throughout the mixed layer   ******
! ******   (a.k.a. PBL) in order to prevent short-lived species  ******
! ******   such as HNO3 from being depleted in the shallow       ******
! ******   surface layer. (rjp, bmy, 7/30/03)                    ******   
! *********************************************************************
!
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
                     RRATE(KLOOP,NK) = 0.E+0_fp 
                  else

                     ! This used to include dust and sea salt aerosols,
                     ! but sea salt deposition is already done in PBL 
                     ! code. So exclude now from here 
                     ! No need to check here for HEMCO drydep rates of
                     ! dust (that shouldn't exist anyways!) since those
                     ! would have been added to the tracers array in 
                     ! DUSTMIX already (ckeller, 10/13/14).
                     SELECT CASE ( DEPNAME(I) )
!                        CASE ( 'DST1', 'DST2', 'DST3', 'DST4', 'SALA', 
!     &                         'SALC' )
                        CASE ( 'DST1', 'DST2', 'DST3', 'DST4' )
                           ! dusts 
                           RRATE(KLOOP,NK) = DEPSAV(IX,IY,I)
                        CASE DEFAULT
                           ! gases + aerosols for full chemistry
                           RRATE(KLOOP,NK) = 0.E+0_fp
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

                  ! Also add deposition values from HEMCO. (ckeller, 10/13/14)
                  TRCID = NTRAIND(N)
                  IF ( TRCID > 0 ) THEN
                     CALL GetHcoVal( TRCID, IX, IY, 1, Fnd, Dep=Dep )
                     IF ( Fnd ) THEN
                       RRATE(KLOOP,NK) = RRATE(KLOOP,NK)
     &                                 + Dep
     &                                 * GET_FRAC_UNDER_PBLTOP(IX,IY,IZ)
                     ENDIF          
                  ENDIF          

                  ! Special case for O3. Increase the deposition frequency (SHIPO3DEP)
                  ! when there is O3 destruction in subgrid ship plume 
                  ! parameterization. This is roughly equivalent to negative
                  ! emissions, which were used previously by PARANOX,
                  ! but caused instability in the chemical solver
                  ! (cdh, 3/21/2013)
                  ! 
                  ! This is now handled through HEMCO.
                  ! (ckeller, 10/13/14)
!                  SELECT CASE ( TRIM(DEPNAME(I)) )
!                     CASE ( 'O3' )
!                        RRATE(KLOOP,NK) = RRATE(KLOOP,NK) + 
!     &                              SHIPO3DEP(IX,IY) * 
!     &                              GET_FRAC_UNDER_PBLTOP( IX, IY, IZ )
!                     CASE DEFAULT  ! Do nothing
!                  END SELECT
               ENDDO
            ENDIF
         ENDDO

      ENDIF
!
! *********************************************************************
! *********************************************************************
! *             REORDER RRATE ARRAY THEN CALL SOLVER                  * 
! *                                                                   *
! * NOTE: If after this point you want to reference a SMVGEAR rxn #,  *
! * then you must use NOLDFNEW(NK,1) instead of NK. (bmy, 8/8/03)     *
! *********************************************************************
!
      NFDH3              = ITHRR(NCS) 
      NFDL2              = NFDH3  + 1 
      NFDREP             = INOREP(NCS)
      NFDREP1            = NFDREP + 1
      NFDH2              = NFDH3  + ITWOR(NCS) 
      NFDL1              = NFDH2  + 1
      NFDH1              = NFDH2  + IONER(NCS) 
      NFDL0              = NFDH1  + 1 
      NALLR              = NALLRAT(NCS) 

!
      DO 730 NKN         = 1, NALLR
         NK                = NOLDFNEW(NKN,NCS)
         IRMA(NKN)         = IRM2(1,NK,NCS) 
         IRMB(NKN)         = IRM2(2,NK,NCS) 
         IRMC(NKN)         = IRM2(3,NK,NCS) 
 730  CONTINUE 
!
! *********************************************************************
!                        REORDER RRATE ARRAY 
! *********************************************************************
!                 TRATE USED HERE AS A DUMMY ARRAY 
! *********************************************************************
!
!
      DO 745 NK          = 1, NTRATES(NCS)
         DO 740 KLOOP      = 1, KTLOOP
            TRATE(KLOOP,NK)  = RRATE(KLOOP,NK)
 740     CONTINUE
 745  CONTINUE
!
      DO 755 NKN         = 1, NALLR
         NK                = NOLDFNEW(NKN,NCS)
         DO 750 KLOOP      = 1, KTLOOP
            RRATE(KLOOP,NKN) = TRATE(KLOOP,NK)
 750     CONTINUE
 755  CONTINUE
!
! *********************************************************************
! REPLACE INACTIVE REACTION RATE COEFFICIENT ARRAY WITH NEW ARRAY 
! THESE REACTIONS HAVE NO ACTIVE LOSS TERMS. PHOTORATE TERMS HERE
! ARE REPLACED IN UPDATE.F .
! *********************************************************************
!                 TRATE USED HERE AS A REAL ARRAY 
! *********************************************************************
!
      DO 765 NKN          = NFDL0, NALLR
         NH                 = NKN + NALLR
         DO 760 KLOOP       = 1, KTLOOP
            TRATE(KLOOP,NKN)  =  RRATE(KLOOP,NKN)
            TRATE(KLOOP,NH)   = -RRATE(KLOOP,NKN)
 760     CONTINUE
 765  CONTINUE
!
! *********************************************************************
!                         RESET NCSP  
! *********************************************************************
! NCS       = 1..NCSGAS FOR GAS CHEMISTRY                           
! NCSP      = NCS       FOR DAYTIME   GAS CHEM            
!           = NCS + ICS FOR NIGHTTIME GAS CHEM           
!
      NCSP            = (IFSUN - 1) * ICS + NCS
!
! *********************************************************************
!                ARCHIVE FOR PLANE-FOLLOWING DIAGNOSTIC   
! *********************************************************************
!
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

      ! Free pointers
      NULLIFY( FRCLND    )
      NULLIFY( QI        )
      NULLIFY( QL        )
      NULLIFY( SUNCOSmid )
      NULLIFY( T         )

!
! *********************************************************************
!                     RETURN TO CALLING PROGRAM
! *********************************************************************
!
      RETURN
!
! *********************************************************************
!                       INTERNAL SUBROUTINES 
! *********************************************************************
!
      CONTAINS

      FUNCTION N2O5( AEROTYPE, TEMP, RH ) RESULT( GAMMA )

      !=================================================================
      ! Internal function N2O5 computes the GAMMA sticking factor
      ! for N2O5 hydrolysis. (mje, bmy, 8/7/030
      ! 
      ! Arguments as Input:
      ! ----------------------------------------------------------------
      ! (1 ) AEROTYPE (INTEGER) : # denoting aerosol type (cf FAST-J)
      ! (2 ) TEMP     (REAL(fp) ) : Temperature [K]
      ! (3 ) RH       (REAL(fp) ) : Relative Humidity [fraction]
      !
      ! NOTES:
      !=================================================================
      
      ! Arguments
      INTEGER, INTENT(IN) :: AEROTYPE
      REAL(fp),  INTENT(IN) :: TEMP, RH

      ! Local variables
      REAL(fp)              :: RH_P, FACT, TTEMP

      ! Function return value
      REAL(fp)              :: GAMMA
      
      !=================================================================
      ! N2O5 begins here!
      !=================================================================

      ! Convert RH to % (max = 100%)
      RH_P  = MIN( RH * 100e+0_fp, 100e0_fp )

      ! Default value
      GAMMA = 0.01e+0_fp

      ! Special handling for various aerosols
      SELECT CASE ( AEROTYPE )

         !----------------
         ! Dust 
         !----------------
         CASE ( 1, 2, 3, 4, 5, 6, 7 )      
                                
            ! Based on unpublished Crowley work
            GAMMA = 0.01e+0_fp

         !----------------
         ! Sulfate
         !----------------
         CASE ( 8 )            
    
            !===========================================================
            ! RH dependence from Kane et al., Heterogenous uptake of 
            ! gaseous N2O5 by (NH4)2SO4, NH4HSO4 and H2SO4 aerosols
            ! J. Phys. Chem. A , 2001, 105, 6465-6470 
            !===========================================================
            ! No RH dependence above 50.0% (lzh, 10/26/2011)
            ! According to Bertram and Thornton, ACP, 9, 8351-8363, 2009
            RH_P  = MIN( RH_P, 50e+0_fp )

            GAMMA = 2.79e-4_fp + RH_P*(  1.30e-4_fp + 
     &                        RH_P*( -3.43e-6_fp + 
     &                        RH_P*(  7.52e-8_fp ) ) )

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
            TTEMP = MAX( TEMP, 282e0_fp )
            FACT  = 10.e0_fp**( -2e+0_fp - 4e-2_fp
     &            *( TTEMP - 294.e+0_fp ) ) / 1e-2_fp

            ! Apply temperature dependence
            GAMMA = GAMMA * FACT

         !----------------
         ! Black Carbon
         !----------------
         CASE ( 9 )  

             ! From IUPAC
             GAMMA = 0.005e+0_fp

         !----------------
         ! Organic Carbon
         !----------------           
         CASE ( 10 )          

            !===========================================================
            ! Based on Thornton, Braban and Abbatt, 2003
            ! N2O5 hydrolysis on sub-micron organic aerosol: the effect
            ! of relative humidity, particle phase and particle size
            !===========================================================
            IF ( RH_P >= 57e+0_fp ) THEN
               GAMMA = 0.03e+0_fp
            ELSE
               GAMMA = RH_P * 5.2e-4_fp
            ENDIF

         !----------------
         ! Sea salt
         ! accum & coarse
         !----------------
         CASE ( 11, 12 )        
          
            ! Based on IUPAC recomendation
            IF ( RH_P >= 62 ) THEN 
               GAMMA = 0.03e+0_fp
            ELSE
               GAMMA = 0.005e+0_fp
            ENDIF

         !----------------
         ! Strat. aerosols
         !----------------
         CASE ( 13, 14 )
       
            ! These are handled outside this routine - something
            ! is wrong if AEROTYPE=13 or 14 reaches this point
            WRITE (6,*) 'Stratospheric aerosols should not '
            WRITE (6,*) 'be passed to general N2O5 het. '
            WRITE (6,*) 'chem. subroutine'
            WRITE (6,*) 'AEROSOL TYPE =',AEROTYPE
            CALL GEOS_CHEM_STOP

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

! *********************************************************************

      FUNCTION HO2( RADIUS, TEMP, DENAIR, SQM, HO2DENS,
     &              AEROTYPE, CONTINENTAL_PBL, Input_Opt )
     &              RESULT( GAMMA )

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
      ! (1 ) RADIUS   (REAL(fp) ) : Aerosol radius [cm]
      ! (2 ) TEMP     (REAL(fp) ) : Temperature [K]
      ! (3 ) DENAIR   (REAL(fp) ) : Air Density [molec/cm3]
      ! (4 ) HO2DENS  (REAL(fp) ) : HO2 Number Density [molec/cm3]
      ! (5 ) SQM      (REAL(fp) ) : Square root of molecular weight [g/mole]
      ! (6 ) AEROTYPE (INTEGER) : # denoting aerosol type (cf FAST-J)
      ! (7 ) CONTINENTAL_PBL (INTEGER)  : Flag set to 1 if the
      !         box is located in the continenal boundary layer,
      !         otherwise it is zero. Also check for ICE/SNOW (to
      !         disable this at high latitudes)
      !
      ! References:
      ! ----------------------------------------------------------------
      ! (1) Jacob, D.J., Heterogeneous chemistry and tropospheric ozone,
      !     Atmos. Environ., 34, 2131-2159, 2000. [full text (pdf)]
      ! (2) J. Mao, Fan, S., Jacob, D. J., and Travis, K. R.: Radical
      !     loss in the atmosphere from Cu-Fe redox coupling in aerosols,
      !     Atmos. Chem. Phys., 13, 509-519, doi:10.5194/acp-13-509-2013,
      !     2013.
      !
      ! REVISION HISTORY:
      !  17 May 2013 - M. Payer    - Add improved HO2 uptake (J. Mao)
      !  22 May 2013 - M. Payer    - Added option to read GAMMA_HO2 from
      !                              input.geos. Recommended value is 0.2
      !                              based on Jacob et al (2000) and Mao
      !                              et al. (2013).
      !  20 Aug 2013 - R. Yantosca - Removed "define.h", this is 
      !                              now obsolete
      !=================================================================
      
      ! References to F90 modules 
      USE GIGC_Input_Opt_Mod, ONLY : OptInput

      IMPLICIT NONE

      ! Arguments
      REAL(fp),  INTENT(IN) :: RADIUS, TEMP, DENAIR, HO2DENS, SQM
      INTEGER, INTENT(IN) :: AEROTYPE, CONTINENTAL_PBL
      TYPE(OptInput), INTENT(IN) :: Input_Opt   ! Input options

      ! Local variables
      REAL(fp)              :: ALPHA
      REAL(fp)              :: delG, Keq, w, H_eff
      REAL(fp)              :: A1, B1, k1, k2, A, B, C
      REAL(fp)              :: kaq, kmt, o2_ss, fluxrxn, DFKG
      REAL(fp)              :: TEST


      ! Avogadro's number
      REAL(fp),  PARAMETER   :: Na = 6.022e+23_fp

      ! Ideal gas constant [atm cm3/mol/K], Raq
      REAL(fp),  PARAMETER   :: Raq=82.e+0_fp

      ! Function return value
      REAL(fp)              :: GAMMA
      
      !=================================================================
      ! HO2 begins here!
      !=================================================================

      ! Default value
      GAMMA = 0.0e+0_fp

      ! Error check
      IF (RADIUS.le.1e-30_fp) THEN
         RETURN
      ENDIF

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
            GAMMA = 0.1e+0_fp

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
         ! SDE 04/18/13: Added stratospheric sulfur aerosols
         !
         !----------------
         CASE ( 8, 9, 10, 11, 12, 13 )  

            ! Mean molecular speed [cm/s]
            w = 14550.5e+0_fp * sqrt(TEMP/(SQM*SQM))

            ! DFKG = Gas phase diffusion coeff [cm2/s]
            DFKG  = 9.45E+17_fp/DENAIR * SQRT(TEMP) * 
     &              SQRT(3.472E-2_fp + 1.E+0_fp/(SQM*SQM))

            !calculate T-dependent solubility and aq. reaction rate constants
            ! hydronium ion concentration
            ! A1 = 1.+(Keq/hplus) 
            ! with Keq = 2.1d-5 [M], Equilibrium constant for 
            ! HO2aq = H+ + O2- (Jacob, 2000)
            !      hplus=10.e+0_fp^(-pH), with pH = 5
            ! B1 = Req * TEMP
            ! with Req = 1.987d-3 [kcal/K/mol], Ideal gas constant
            ! Note that we assume a constant pH of 5.
            A1 = 1.+ (2.1e-5_fp / (10.e+0_fp**(-5) ) )
            B1 = 1.987e-3_fp * TEMP

            ! Free energy change for solvation of HO2 (Hanson 1992, Golden 1991)
            ! in [kcal/mol]:
            ! delG = -4.9-(TEMP-298e+0_fp)*delS
            ! with delS=-0.023  [kcal/mol/K],  Entropy change for solvation of HO2
            delG  = -4.9e+0_fp - (TEMP-298.e+0_fp) * (-0.023)
            H_eff = exp( -delG / B1 ) * A1

            ! Estimated temp dependent value for HO2 + O2- (k1) and 
            ! HO2+HO2 (see Jacob 1989)
            k1  =   1.58e+10_fp * exp( -3. / B1 )
            k2  =   2.4e+9_fp   * exp( -4.7 / B1 )
            kaq = ( k1 * (A1 - 1.e+0_fp) + k2) / (A1**2)

            ! Calculate the mass transfer rate constant and s.s. conc. of 
            ! total HO2 in the aqueous phase:
            ! kmt = (RADIUS/DFKG + 4e+0_fp/w/alpha)^(-1)
            ! with alpha = mass accomodation coefficient, assumed 
            ! to be 1 (Thornton et al.)
            kmt = 1.e+0_fp/( RADIUS/DFKG + 4e+0_fp/w/1. )

            !use quadratic formula to obtain [O2-] in particle of radius RADIUS
            A = -2e+0_fp * kaq
            B = -3e+0_fp * kmt / RADIUS / (H_eff * 0.082 * TEMP)
            C =  3e+0_fp * kmt * HO2DENS * 1000e+0_fp / RADIUS / Na

            ! Error check that B^2-(4e+0_fp*A*C) is not negative
            TEST= B**2-(4e+0_fp*A*C)
            IF ( TEST < 0e+0_fp ) THEN
                GAMMA = 0e+0_fp
            ELSE
                ! Calculate the concentration of O2- in the aerosol
                o2_ss= ( -B  -sqrt(B**2-(4e+0_fp*A*C)) )/(2e+0_fp*A)

                ! Calculate the reactive flux
                fluxrxn = kmt*HO2DENS - o2_ss*Na*kmt/H_eff/Raq/TEMP

                IF ( fluxrxn <= 0e0_fp ) THEN
                   GAMMA = 0e+0_fp
                ELSE
                   ! Gamma for HO2 at TEMP, ho2, and RADIUS given
                   GAMMA = 1./( ( ( HO2DENS/fluxrxn ) - 
     &                            ( RADIUS/DFKG ) ) * w / 4.e+0_fp )
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
         ! NAT/ice (SDE 04/18/13)
         !----------------
         CASE ( 14 )
       
            GAMMA = 0.e+0_fp

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
      IF ( GAMMA  <= 0e+0_fp ) GAMMA = 0e+0_fp

      ! This is for the improved HO2 uptake (J. Mao)
      GAMMA = Input_Opt%GAMMA_HO2

      ! Return to CALCRATE
      END FUNCTION HO2
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cld1k_brno3
!
! !DESCRIPTION: Function CLD1K_BrNO3 calculates the rate constant for
!  heterogeneous cycling of BrNO3 off of cloud particles.
!\\
!\\
! !INTERFACE:
!
      FUNCTION CLD1K_BrNO3( I,      J,        L,   DENAIR, 
     &                      QL,     State_Met )    RESULT( cld1k )

!
! !USES:
!
      USE DAO_MOD,            ONLY : IS_LAND, IS_ICE
      USE GIGC_State_Met_Mod, ONLY : MetState
      IMPLICIT NONE
!
! !INPUT PARAMETERS:
!
      INTEGER,        INTENT(IN) :: I         ! Longitude index
      INTEGER,        INTENT(IN) :: J         ! Latitude  index
      INTEGER,        INTENT(IN) :: L         ! Altitude  index
      REAL(fp),         INTENT(IN) :: DENAIR    ! Density of air [#/cm3]
      REAL(fp),         INTENT(IN) :: QL        ! Cloud water mixing ratio [kg/kg]
      TYPE(MetState), INTENT(IN) :: State_Met ! Meteorology State object
!
! !RETURN VALUE:
!
      REAL(fp)              :: cld1k    ! Rate constant for heterogeneous cycling
                                      ! of BrNO3 off of cloud particles
!
! !REMARKS:
!  The rate constant for heterogeneous cycling of BrNO3 off of cloud particles
!  is calculated assuming:
!                                                                             .
!    1. A sticking coefficient of 0.3 [Yang et al. 2005]
!    2. uniform cloud droplet size for 2 types of clouds
!       - continental warm clouds: r =  6d-4 [cm]
!       - marine warm clouds:      r = 10d-4 [cm]
!       * no distributions are assumed
!
!  ** Calculation of a 1st order rate constent barrowed from the
!     subroutine arsl1k.f. Below are comments from that code:
!                                                                             .
!       The 1st-order loss rate on wet aerosol (Dentener's Thesis, p. 14)
!       is computed as:
!                                                                             .
!         ARSL1K [1/s] = area / [ radius/dfkg + 4./(stkcf * nu) ]        
!                                                                             .
!       where nu   = Mean molecular speed [cm/s] = sqrt(8R*TK/pi/M) for Maxwell
!             DFKG = Gas phase diffusion coeff [cm2/s] (order of 0.1)
!
! !REVISION HISTORY:
!  27 Feb 2011 - J. Parrella - Initial version
!  22 May 2012 - M. Payer    - Added ProTeX headers
!  09 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
      ! Cloud droplet radius in continental warm clouds [cm]
      REAL(fp), PARAMETER :: XCLDR_CONT =  6.e-4_fp

      ! Cloud droplet radius in marine warm clouds [cm]
      REAL(fp), PARAMETER :: XCLDR_MARI = 10.e-4_fp

      REAL(fp), PARAMETER :: R = 8.314472                  ! [J/mol/K]
      REAL(fp), PARAMETER :: mw_brno3 = 0.142              ! [kg/mol]
      rEAL(fp), PARAMETER :: pi = 3.14159265358979323846e+0_fp ! [unitless]
      REAL(fp), PARAMETER :: alpha = 0.3                   ! sticking coefficient
      REAL(fp), PARAMETER :: dens_h2o = 0.001e+0_fp            ! [kg/cm3]
!
! !LOCAL VARIABLES:
!
      REAL(fp)            :: nu         ! Mean molecular speed
      REAL(fp)            :: RADIUS     ! Radius of cloud droplet      [cm]
      REAL(fp)            :: SQM        ! Square root of molec. weight [g/mol]
      REAL(fp)            :: STK        ! Square root of temperature   [K]
      REAL(fp)            :: AREA       ! Surface area                 [m2]
      REAL(fp)            :: DFKG       ! Gas diffusion coefficient    [cm2/s]
      REAL(fp)            :: Vc         ! Volume of the cloud          [cm3]
      REAL(fp)            :: XAIRM3     ! Volume of air                [m3]
      LOGICAL           :: yn_continue
   
      ! Pointers
      REAL(fp), POINTER   :: AD(:,:,:)
      REAL(fp), POINTER   :: AIRVOL(:,:,:)
      REAL(fp), POINTER   :: CLDF(:,:,:)
      REAL(fp), POINTER   :: FRLAND(:,:)
      REAL(fp), POINTER   :: FROCEAN(:,:)
      REAL(fp), POINTER   :: T(:,:,:)

      !=================================================================
      ! CLD1K_BrNO3 begins here!
      !=================================================================

      ! Initialize pointers
      AD      => State_Met%AD
      AIRVOL  => State_Met%AIRVOL
      CLDF    => State_Met%CLDF
      FRLAND  => State_Met%FRLAND
      FROCEAN => State_Met%FROCEAN
      T       => State_Met%T

      ! ----------------------------------------------
      ! 1.
      !   calculate the mean molecular speed of the
      !   molecules given the temperature.
      ! ----------------------------------------------
      nu   = sqrt( 8.e+0_fp * R * T(I,J,L) / (mw_brno3 * pi) )

      ! ----------------------------------------------
      ! Test conditions to see if we want to continue
      ! or set the cloud rate equal to zero.
      ! ----------------------------------------------

      ! continental or marine clouds only...
#if defined( GEOS_5 ) || defined( MERRA ) || defined( GEOS_FP )
      IF ( (FRLAND (I,J) > 0) .or. (FROCEAN(I,J) > 0) ) THEN
#else
      ! Above line is to skip over land ice (Greenland and Antartica). This
      ! should do the same (and also work for GEOS-5, but leave above for now).
      IF ( IS_LAND( I, J, State_Met )    .and. .not.
     &   ( IS_ICE ( I, J, State_Met )) ) THEN
#endif
         ! do we have clouds? and do we have warm temperatures?
         IF ( ( CLDF(L,I,J) > 0    )   .and.
     &        ( T(I,J,L)    > 258.0) ) THEN
            yn_continue = .TRUE.
         ELSE
            yn_continue = .FALSE.
         ENDIF
      ELSE
         yn_continue = .FALSE.
      ENDIF

      ! test
      IF ( .not. yn_continue ) THEN
         ! nothing to calculate...
         cld1k = 0.e+0_fp
         RETURN
      ENDIF


      ! ----------------------------------------------
      ! 2.
      !   calculate the surface area of cloud droplets
      !   in the given grid box, assuming 1 of 2
      !   conditions:
      !     a. marine warm cloud
      !       or
      !     b. continental warm cloud
      !
      !
      !   * Calculation for area is derived follows,
      !     assuming that RADIUS is constant:
      !
      !                         4/3 (pi) (RADIUS)**3
      !  1) FC = Vc / Vb = N  -------------------------
      !                                  Vb
      !
      !
      !       where N      = number of cloud droplets
      !             RADIUS = radius of cloud droplet
      !             Vc     = volumn of the cloud
      !             Vb     = volumn of the box = AIRVOL (in GEOS-Chem)
      !
      !
      !                     Vb
      !  2) N = FC --------------------
      !            4/3 (pi) (RADIUS)**3
      !
      !
      !  So the surface area [m2] is calculated as
      !
      !  3) total surface A = N * 4 * (pi) * (RADIUS)**2
      !
      !                  3*Vb
      !          = FC ----------
      !                 RADIUS
      !
      !  4) for this routine though we want
      !     AREA in [cm2/cm3], surface area to volume air:
      !
      !                   3
      !     AREA = FC ---------
      !                RADIUS (in cm)
      !
      !
      !    or    
      !                   3 x Vc
      !     AREA =  -----------------
      !              AIRVOL x RADIUS      (in cm)
      ! ----------------------------------------------
#if defined( GEOS_5 ) || defined( MERRA ) || defined( GEOS_FP )
      IF ( FRLAND(I,J) > FROCEAN(I,J) ) THEN
#else
      IF ( IS_LAND(I,J, State_Met) ) THEN
#endif
         ! Continental cloud droplet radius [cm]
         RADIUS = XCLDR_CONT
      ELSE
         ! Marine cloud droplet radius [cm]
         RADIUS = XCLDR_MARI
      ENDIF

      ! store the volume of air [m3]
      XAIRM3 = AIRVOL(I,J,L)
      ! convert to [cm3]
      XAIRM3 = XAIRM3 * (100.e+0_fp)**3

      ! get the volume of cloud [cm3]
#if defined( GEOS_5 ) || defined( MERRA ) || defined( GEOS_FP )
      ! QL is [g/g]
      Vc = CLDF(L,I,J) * QL * AD(I,J,L) / dens_h2o
#else
      ! QL is [cm3/cm3]
      Vc = CLDF(L,I,J) * QL * XAIRM3
#endif

      ! now calculate the cloud droplet surface area
      AREA    = 3.e+0_fp * (Vc/XAIRM3) / (RADIUS) ! keep Radius in [cm]

      ! ----------------------------------------------------
      ! 3.
      !   Now finish calculating the 1st order rate
      !   constant for BrNO3 hydrolysis.
      !
      !   (a) calculate the gas phase diffusion coefficient;
      !
      !   (b) calculate the hydrolysis rxn rate.
      ! ----------------------------------------------------
      SQM = sqrt(mw_brno3 * 1.e+3_fp)    ! square root of molar mass [g/mole]
      STK = sqrt(T(I,J,L)) ! square root of temperature [K]

      ! DFKG = Gas phase diffusion coeff [cm2/s] (order of 0.1)
      DFKG  = 9.45E+17_fp/DENAIR * STK * SQRT(3.472E-2_fp 
     &     + 1.E+0_fp/(SQM*SQM))

      ! Compute ARSL1K according to the formula listed above
      cld1k = AREA / ( RADIUS/DFKG + 2.749064E-4 
     &     * SQM/(alpha*STK) )

      ! Free Pointers
      NULLIFY( AD      )
      NULLIFY( AIRVOL  )
      NULLIFY( CLDF    )
      NULLIFY( FRLAND  )
      NULLIFY( FROCEAN )
      NULLIFY( T       )

      END FUNCTION CLD1K_BrNO3
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: fcro2ho2
!
! !DESCRIPTION: !fgap, based on saunder 2003 k14
!\\
!\\
! !INTERFACE:
!
      REAL(fp) FUNCTION FCRO2HO2( XCARBN )
!
! !REVISION HISTORY:
!  24 Jul 2014 - R. Yantosca - Now inlined to calcrate.F
!EOP
!------------------------------------------------------------------------------
!BOC
!
      IMPLICIT NONE

      ! Arguments
      REAL(fp), INTENT(IN) :: XCARBN

      FCRO2HO2 = 1E+0_fp-EXP(-0.245E+0_fp*XCARBN)
     
      ! Return to calling program
      END FUNCTION FCRO2HO2
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !FUNCTION: FYHORO
!
! !DESCRIPTION: \subsection*{Overview}
!  Function FYHORO returns returns the branching ratio between 
!  HOC2H4O oxidation and dissociation:
!  (1) HOC2H4 --O2--> HO2 + GLYC
!  (2) HOC2H4 ------> HO2 + 2CH2O
!
!\subsection*{References}
!  \begin{enumerate}
!  \item Orlando et al., 1998: \emph{Laboratory and theoretical study of the 
!         oxyradicals in the OH- and Cl-initiated oxidation of ethene}, 
!        \underline{J. Phys. Chem. A}, \textbf{102}, 8116-8123.
!  \item Orlando et al., 2003: \emph{The atmospheric chemistry of alkoxy 
!         radicals}, \underline{Chem. Rev.}, \textbf{103}, 4657-4689.
!  \end{enumerate}
!
!\\
!\\
! !INTERFACE: 
!
      REAL(fp) FUNCTION FYHORO( ZDNUM, TT )
! 
! !USES:
!
      IMPLICIT NONE
!
! !INPUT PARAMETERS: 
!
      REAL(fp), INTENT(IN) :: ZDNUM   ! Air density   [molec/cm3 ]
      REAL(fp), INTENT(IN) :: TT      ! Temperature   [K         ]
!
! !REVISION HISTORY:
!  (1 ) Branching ratio calculation (tmf, 2/6/05).
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!  25 Jul 2014 - R. Yantosca - Now inlined into calcrate.F
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL(fp) :: K1, K2, O2DNUM

      !=================================================================
      ! FYHORO begins here!
      !=================================================================
      O2DNUM = ZDNUM * 0.21E+0_fp
      K1     = 6.0E-14_fp * EXP(-550.E+0_fp/TT) * O2DNUM
      K2     = 9.5E+13_fp * EXP(-5988.E+0_fp/TT) 

      FYHORO = K1 / (K1 + K2)
     
      END FUNCTION FYHORO
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: fyrno3
!
! !DESCRIPTION: Function FYRNO3 returns organic nitrate yields 
!  YN = RKA/(RKA+RKB) from RO2+NO reactions as a function of the number 
!  N of carbon atoms.
!\\
!\\
! !INTERFACE:
!
      REAL(fp) FUNCTION FYRNO3( XCARBN, ZDNUM, TT )
!
! !USES:
!
      IMPLICIT NONE
!
! !INPUT PARAMETERS: 
!
      REAL(fp), INTENT(IN) :: XCARBN   ! Number of C atoms in RO2
      REAL(fp), INTENT(IN) :: ZDNUM    ! Air density   [molec/cm3 ]
      REAL(fp), INTENT(IN) :: TT       ! Temperature   [K         ]
! 
! !REVISION HISTORY: 
!  (1 ) Original code from Larry Horowitz, Jinyou Liang, Gerry Gardner,
!        and Daniel Jacob circa 1989/1990.
!  (2 ) Updated following Atkinson 1990.
!  (3 ) Change yield from Isoprene Nitrate (ISN2) from 0.44% to 12%,
!        according to Sprengnether et al., 2002. (amf, bmy, 1/7/02)
!  (4 ) Eliminate obsolete code from 1/02 (bmy, 2/27/02)
!  (5 ) Updated comment description of XCARBN (bmy, 6/26/03)
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!  25 Jul 2014 - R. Yantosca - Now inlined into calcrate.F
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL(fp) :: YYYN, XXYN,  AAA,  RARB, ZZYN
      REAL(fp) :: XF,   ALPHA, Y300, BETA, XMINF, XM0

      ! Initialize static variables
      DATA   Y300,ALPHA,BETA,XM0,XMINF,XF/.826,1.94E-22,.97,0.,8.1,.411/

      !=================================================================
      ! FYRNO3 begins here!
      !=================================================================
      XXYN   = ALPHA*EXP(BETA*XCARBN)*ZDNUM*((300./TT)**XM0)
      YYYN   = Y300*((300./TT)**XMINF)
      AAA    = LOG10(XXYN/YYYN)
      ZZYN   = 1./(1.+ AAA*AAA )
      RARB   = (XXYN/(1.+ (XXYN/YYYN)))*(XF**ZZYN)
      FYRNO3 = RARB/(1. + RARB)
     
      END FUNCTION FYRNO3
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: arsl1k
!
! !DESCRIPTION: Function ARSL1K calculates the 1st-order loss rate of species 
!  on wet aerosol surface.
!\\
!\\
! !INTERFACE:
!
      REAL(fp) FUNCTION ARSL1K( AREA, RADIUS, DENAIR, STKCF, STK, SQM )
!
! !USES:
!
      IMPLICIT NONE
!
! !INPUT PARAMETERS: 
!
      ! Surface  area of wet aerosols/volume of air [cm2/cm3]
      REAL(fp), INTENT(IN) :: AREA     

      ! Radius of wet aerosol [cm], order of 0.01-10 um;
      ! Note that radius here is Rd, not Ro
      REAL(fp), INTENT(IN) :: RADIUS 
  
      ! Density of air [#/cm3]
      REAL(fp), INTENT(IN) :: DENAIR  
 
      ! Sticking coefficient [unitless], order of 0.1
      REAL(fp), INTENT(IN) :: STKCF  
  
      ! Square root of temperature [K]
      REAL(fp), INTENT(IN) :: STK  
    
      ! Square root of molecular weight [g/mole]
      REAL(fp), INTENT(IN) :: SQM      
!
! !REMARKS:
!  The 1st-order loss rate on wet aerosol (Dentener's Thesis, p. 14)
!  is computed as:
!                                                                             .
!      ARSL1K [1/s] = area / [ radius/dfkg + 4./(stkcf * xmms) ]        
!                                                                             .
!  where XMMS = Mean molecular speed [cm/s] = sqrt(8R*TK/pi/M) for Maxwell 
!        DFKG = Gas phase diffusion coeff [cm2/s] (order of 0.1)

! !REVISION HISTORY:
!  01 Jul 1994 - lwh, jyl, gmg, djj - Initial version 
!  04 Apr 2003 - R. Yantosca - Updated comments, cosmetic changes
!  07 Apr 2004 - R. Yantosca - Now return w/ default value if RADIUS is zero 
!                              (i.e. is smaller than a very small number)
!  03 Dec 2009 - R. Yantosca - Prevent div-by-zero errors by returning the
!                              default value if any of the args are zero 
!  03 Dec 2009 - R. Yantosca - Added ProTeX Header
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL(fp) :: DFKG

      !=================================================================
      ! ARSL1K begins here!
      !=================================================================
      IF ( AREA < 0e0_fp   .or. DENAIR < 1e-30_fp .or. 
     &     RADIUS < 1e-30_fp  .or.
     &     SQM  < 1e-30_fp .or. STK    < 1e-30_fp .or. 
     &     STKCF  < 1e-30_fp ) THEN

         ! Use default value if any of the above values are zero
         ! This will prevent div-by-zero errors in the eqns below
         ! Value changed from 1d-3 to 1d-30 (bhh, jmao, eam, 7/18/2011)
         ARSL1K = 1.E-30_fp

      ELSE

         ! DFKG = Gas phase diffusion coeff [cm2/s] (order of 0.1)
         DFKG  = 9.45E+17_fp/DENAIR * STK * SQRT(3.472E-2_fp + 1.E0_fp/
     &    (SQM*SQM))

         ! Compute ARSL1K according to the formula listed above
         ARSL1K = AREA / ( RADIUS/DFKG + 2.749064E-4*SQM/(STKCF*STK) )

      ENDIF

      ! Return to calling program
      END FUNCTION ARSL1K
!EOC
! *********************************************************************
! ******************** END OF SUBROUTINE CALCRATE *********************
! *********************************************************************
!
      END SUBROUTINE CALCRATE
