!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: tagged_co_mod
!
! !DESCRIPTION: Module TAGGED\_CO\_MOD contains variables and routines 
!  used for the  geographically tagged CO simulation.
!\\
!\\
! !INTERFACE: 
!
      MODULE TAGGED_CO_MOD
!
! !USES:
!
      USE PRECISION_MOD    ! For GEOS-Chem Precision (fp)

      IMPLICIT NONE
      PRIVATE
!
! !PRIVATE MEMBER FUNCTIONS:
!
      PRIVATE :: CALC_DIURNAL
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC  :: CHEM_TAGGED_CO
      PUBLIC  :: GET_ALPHA_ISOP
      PUBLIC  :: READ_PCO_LCO_STRAT
      PUBLIC  :: GET_PCO_LCO_STRAT
      PUBLIC  :: INIT_TAGGED_CO
      PUBLIC  :: CLEANUP_TAGGED_CO
!
! !PUBLIC DATA MEMBERS:
!
      INTEGER, PUBLIC, ALLOCATABLE :: BB_REGION(:,:)   ! BB CO regions
      INTEGER, PUBLIC, ALLOCATABLE :: FF_REGION(:,:)   ! FF CO regions
      REAL(fp),  PUBLIC, ALLOCATABLE :: EMACET(:,:)      ! Acetone emissions
!
! !REMARKS:
!  Tagged CO Tracers (you can modify these as needs be!)
!  ============================================================================
!  (1 ) Total CO
!  (2 ) CO from North American fossil fuel 
!  (3 ) CO from European fossil fuel 
!  (4 ) CO from Asian fossil fuel 
!  (5 ) CO from fossil fuel from everywhere else
!  (6 ) CO from South American biomass burning 
!  (7 ) CO from African biomass burning 
!  (8 ) CO from Southeast Asian biomass burning
!  (9 ) CO from Oceania biomass burning
!  (10) CO from European biomass burning
!  (11) CO from North American biomass burning
!  (12) CO chemically produced from Methane
!  (13) CO from Biofuel burning (whole world)
!  (14) CO chemically produced from Isoprene 
!  (15) CO chemically produced from Monoterpenes
!  (16) CO chemically produced from Methanol (CH3OH)
!  (17) CO chemically produced from Acetone
!
! !REVISION HISTORY:
!  28 Jul 2000- R. Yantosca - Initial version
!  (1 ) Removed obsolete code from CHEM_TAGGED_CO (bmy, 12/21/00)
!  (2 ) Added CO sources from oxidation of biofuel VOC's, biomass burning
!        VOC's, fossil fuel VOC's, and natural VOC's (bnd, bmy, 1/2/01)
!  (3 ) Added chemical P(CO) from CH3OH and MONOTERPENES (bnd, bmy, 1/2/01)
!  (4 ) Now cap SCALEYEAR at 1997 in "emiss_tagged_co" (bnd, bmy, 4/6/01) 
!  (5 ) Removed obsolete commented-out code (bmy, 4/23/01)
!  (6 ) Added new module variables SUMACETCO, EMACET, CO_PRODS, CO_LOSSS,
!        ISOP96, MONO96, and MEOH96.  Also added new module routines 
!        GET_ALPHA_ISOP, READ_PCO_LCO_STRAT, GET_PCO_LCO_STRAT, 
!        READ_ACETONE, and READ_BIOG_FOR_GEOS3. (bnd, bmy, 6/14/01)
!  (7 ) Now read files from DATA_DIR/tagged_CO_200106/ (bmy, 6/19/01) 
!  (8 ) Removed ISOP96, MONO96, and CH3OH96 since we now use the new GEOS-3
!        fields and no longer have to correct for the surface temperature.
!        (bmy, 8/21/01)
!  (9 ) Bug fix: don't call GLOBAL_NOX_MOD in routine CHEM_TAGGED_CO 
!        unless a logical switch is set (bmy, 8/28/01)
!  (10) Updated comments (bmy, 9/6/01)
!  (11) Deleted obsolete code for 1998 GEOS-3 fix.  Also now archive ND46
!        diagnostic as [atoms C/cm2/s] (bmy, 9/13/01)
!  (12) Bug fix in CHEM_TAGGED_CO: now save CO sources/sinks into the
!        Total CO tracer (N=1). (qli, bmy, 9/21/01)
!  (13) Resize arrays of (IIPAR,JJPAR) to (IIPAR,JJPAR) (bmy, 9/28/01)
!  (14) Removed obsolete code from 9/28/01 (bmy, 10/23/01)
!  (15) Updated comments (bmy, 2/15/02)
!  (16) Removed double definition of SUMCH3OHCO (bmy, 3/20/02)
!  (17) Now use P(I,J) + PTOP instead of PS(I,J) (bmy, 4/11/02)
!  (18) Now divide module header into MODULE PRIVATE, MODULE VARIABLES, and
!        MODULE ROUTINES sections.  Updated comments (bmy, 5/28/02)  
!  (19) Now references "pressure_mod.f" (dsa, bdf, bmy, 8/21/02)
!  (20) Now reference AD, BXHEIGHT, T and SUNCOS from "dao_mod.f".  Also 
!        removed obsolete code from various routines.  Now references
!        ERROR_STOP from "error_mod.f". (bmy, 10/15/02)
!  (21) Now references "grid_mod.f" and the new "time_mod.f". (bmy, 2/3/03)
!  (22) Bug fix for NTAU in EMISS_TAGGED_CO.  Bug fix for FILENAME in routine
!        READ_PCO_LCO_STRAT. (ave, bnd, bmy, 6/3/03)
!  (23) Updated arg list in call to EMISOP in EMISS_TAGGED_CO (bmy, 12/9/03)
!  (24) Now references "directory_mod.f", "logical_mod.f", "tracer_mod.f".
!        Now remove IJLOOP_CO. (bmy, 7/20/04)
!  (25) Fixed bug in CHEM_TAGGED_CO (bmy, 3/7/05)
!  (26) Now reads data from both GEOS and GCAP grids.  Now also references
!        "tropopause_mod.f". (bmy, 8/16/05)
!  (27) Now modified for new "biomass_mod.f" (bmy, 4/5/06)
!  (28) BIOMASS(:,:,IDBCO) from "biomass_mod.f" is now in units of 
!        [molec CO/cm2/s].  Adjust unit conversion accordingly. (bmy, 9/27/06)
!  (29) Routines GET_ALPHA_ISOP, GET_PCO_LCO_STRAT, READ_PCO_LCO_STRAT,
!        READ_ACETONE and INIT_TAGGED_CO are public now. Variable
!        EMACET is public now.  (phs, 9/18/07)
!  13 Aug 2010 - R. Yantosca - Add modifications for MERRA (treat like GEOS-5)
!  08 Feb 2012 - R. Yantosca - Add modifications for GEOS-5.7.2 met
!  01 Mar 2012 - R. Yantosca - Now reference new grid_mod.F90
!  23 Oct 2012 - R. Yantosca - Update prod/loss for new GMI strat chem
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!  10 Apr 2014 - R. Yantosca - Now change CO_PRODS and CO_LOSSS to REAL*4
!  11 Jun 2014 - J.A. Fisher - Replace SUNCOS with SUNCOSmid for MEGAN
!  11 Jun 2014 - J.A. Fisher - Treat monoterpenes as in full chemistry
!  12 Jun 2014 - J.A. Fisher - Add option for non-local PBL mixing
!  10 Jul 2014 - J.A. Fisher - Add diurnal cycle for OH
!  10 Jul 2014 - J.A. Fisher - Change FIRSTEMISS to FIRSTTIME to allow
!                              initialisation when emissions are off
!  20 Nov 2014 - M. Yannetti - Added PRECISION_MOD
!  11 Mar 2015 - R. Yantosca - Add OH pointer for HEMCO data input
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE TYPES:
!
      REAL(f4), ALLOCATABLE :: CO_PRODS  (:,:,:)   ! Strat P(CO) from CH4
      REAL(f4), ALLOCATABLE :: CO_LOSSS  (:,:,:)   ! Strat L(CO) from CO + OH
      REAL(fp), ALLOCATABLE :: SUMACETCO (:,:  )   ! P(CO) from Acetone
      REAL(fp), ALLOCATABLE :: SUMCH3OHCO(:,:  )   ! P(CO) from CH3OH
      REAL(fp), ALLOCATABLE :: SUMISOPCO (:,:  )   ! P(CO) from Isoprene
      REAL(fp), ALLOCATABLE :: SUMMONOCO (:,:  )   ! P(CO) from Monoterpenes
      REAL(fp), ALLOCATABLE :: TCOSZ     (:,:  )   ! Daily sum COS(SZA)

      ! Pointers to fields in the HEMCO data structure.
      ! These need to be declared REAL(f4), aka REAL*4.
      REAL(f4), POINTER     :: OH(:,:,:) => NULL()
!
! !DEFINED PARAMETERS:
!
      ! FMOL_CO     - kg CO / mole CO
      ! XNUMOL_CO   - molecules CO / kg CO
      REAL(fp),  PARAMETER   :: FMOL_CO     = 28e-3_fp
      REAL(fp),  PARAMETER   :: XNUMOL_CO   = 6.022e+23_fp/FMOL_CO

      ! Molecules C / kg C
      REAL(fp),  PARAMETER   :: XNUMOL_C    = 6.022e+23_fp/12e-3_fp 

      CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: chem_tagged_co
!
! !DESCRIPTION: Subroutine CHEM\_TAGGED\_CO performs CO chemistry on 
!  geographically "tagged" CO tracers.  Loss is via reaction with OH. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CHEM_TAGGED_CO( am_I_Root, Input_Opt,
     &                           State_Met, State_Chm, RC )
!
! !USES:
!
      USE CHEMGRID_MOD,       ONLY : ITS_IN_THE_STRAT
      USE CMN_SIZE_MOD
      USE CMN_DIAG_MOD
      USE DIAG_PL_MOD,        ONLY : AD65
      USE ERROR_MOD,          ONLY : CHECK_VALUE
      USE ERROR_MOD,          ONLY : ERROR_STOP
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GIGC_State_Chm_Mod, ONLY : ChmState
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE HCO_EmisList_Mod,   ONLY : HCO_GetPtr
!      USE HCOI_GC_MAIN_MOD,   ONLY : GetHcoVal
      USE GRID_MOD,           ONLY : GET_YMID
      USE PRESSURE_MOD,       ONLY : GET_PCENTER
      USE TIME_MOD,           ONLY : GET_TS_CHEM,     GET_TS_EMIS
      USE TIME_MOD,           ONLY : GET_MONTH,       GET_YEAR
      USE TIME_MOD,           ONLY : ITS_A_NEW_MONTH, ITS_A_NEW_YEAR
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
      INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
! 
! !REVISION HISTORY: 
!  19 Oct 1999 - Q. Li, B. Duncan, B. Field - Initial version
!  (1 ) Now do chemistry all the way to the model top. 
!  (2 ) Use monthly mean OH fields for oxidation.
!  (3 ) Now reference the monthly mean OH array and the routine which reads
!        it from disk in "global_oh_mod.f" (bmy, 7/28/00)
!  (4 ) Removed obsolete code from 10/6/00 (bmy, 12/21/00)
!  (5 ) Added P(CO) from CH3OH and MONOTERPENES.  Also account for the 
!        variation of CH4 conc. w/ latitude and year (bnd, bmy, 1/2/01)
!  (6 ) Removed obsolete commented-out code (bmy, 4/23/01)
!  (7 ) Updated with new stuff from bnd.  Also references module routines
!        & arrays from "global_nox_mod.f" and "error_mod.f".  Removed
!        THISMONTH as an argument since we can use the MONTH value
!        from the "CMN" include file. (bmy, 6/14/01)
!  (8 ) Remove GEOS-3 fix for biogenic fields (ISOP96, MONO96, CH3OH96),
!        since we now use updated met fields w/o the surface temperature
!        problem (bmy, 8/21/01)
!  (9 ) Now only call GLOBAL_NOX_MOD to define the BNOX array, if switch
!        ALPHA_ISOP_FROM_NOX is set.  The NOx concentrations used to compute
!        the CO yield from isoprene are currently only at 4x5 (bmy, 8/28/01)
!  (10) Bug fix: now make sure to add CO production and loss into the 
!        STT(:,:,:,1), which is the Total CO tracer. (qli, bmy, 9/21/01)
!  (11) Updated comments.  Bug fix: multiply CO_OH (after using it to update 
!        tagged tracers) by GCO (the initial value of STT in molec/cm3) to 
!        convert it to an amount of CO lost by OH [molec/cm3]. (bmy, 2/19/02)
!  (12) Removed PS as an argument; use P(I,J) + PTOP instead of PS, in order
!        to ensure that we use P and AD computed from the same preesure
!        by AIRQNT. (bmy, 4/11/02)
!  (13) Now use GET_PCENTER from "pressure_mod.f" to compute the pressure
!        at the midpoint of box (I,J,L).  Also deleted obsolete, commented-out
!        code. (dsa, bdf, bmy, 8/21/02)
!  (14) Now reference AD and T from "dao_mod.f".  Now make FIRSTCHEM a local
!        SAVEd variable. (bmy, 11/15/02)
!  (15) Now replace YLMID(J) with routine GET_YMID of "grid_mod.f".  Now
!        uses functions GET_TS_CHEM, GET_MONTH, GET_YEAR from the new
!        "time_mod.f".  (bmy, 2/10/03) 
!  (16) Now reference STT & N_TRACERS from "tracer_mod.f".  Now references
!        LSPLIT from "logical_mod.f".  Now references AD65 from 
!        "diag_pl_mod.f".  Updated comments. (bmy, 7/20/04)
!  (17) Bug fix: re-insert ELSE between (1a-1) and (1a-2); it appears to have
!        been mistakenly deleted. (bmy, 3/7/05)
!  (18) Now references ITS_IN_THE_STRAT from "tropopause_mod.f".  Now remove
!        reference to "CMN", it's obsolete. (bmy, 8/22/05)
!  (19) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (20) Remove reference to "global_ch4_mod.f" (bmy, 5/31/06)
!  (21) Use newest JPL 2006 rate constant for CO+OH (jaf, jmao, 3/4/09)
!  01 Mar 2012 - R. Yantosca - Now use GET_AREA_CM2(I,J,L) from grid_mod.F90
!  01 Mar 2012 - R. Yantosca - Now use GET_YMID(I,J,L) from grid_mod.F90
!  22 Oct 2012 - R. Yantosca - Now pass am_I_Root=.TRUE. to GET_GLOBAL_CH4
!  23 Oct 2012 - R. Yantosca - Added ProTeX headers
!  09 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!  25 Mar 2013 - R. Yantosca - Now accept am_I_Root, Input_Opt, State_Chm, RC
!  11 Apr 2013 - R. Yantosca - Now pass Input_Opt to GET_GLOBAL_NOX
!  10 Jun 2013 - R. Yantosca - Avoid array temporaries in CHECK_VALUE
!  23 Oct 2013 - R. Yantosca - Now pass objects to GET_GLOBAL_OH routine
!  10 Apr 2014 - R. Yantosca - Now use pointers to fields in State_Met
!  10 Apr 2014 - R. Yantosca - Bug fix: DTSRCE should be DTCHEM
!  10 Jul 2014 - J.A. Fisher - Add diurnal cycle for OH
!  10 Jul 2014 - J.A. Fisher - Restore assignment of DTCHEM
!  10 Jul 2014 - J.A. Fisher - Remove unused DTSRCE, PCO, and AREA_CM2
!  10 Jul 2014 - J.A. Fisher - Moved ALPHA definitions out of code
!  27 Aug 2014 - R. Yantosca - Bug fix: pass Input_Opt to GET_GLOBAL_CH4
!  21 Oct 2014 - C. Keller   - Now use GetHcoVal instead ot Trac_Tend array
!  11 Mar 2015 - R. Yantosca - Remove references to global_nox_mod.F
!  11 Mar 2015 - R. Yantosca - Now get OH via the HEMCO data structure
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL                :: LSPLIT
      INTEGER                :: N_TRACERS
      INTEGER                :: I,         J,          L
      INTEGER                :: N,         MONTH,      YEAR
      REAL(fp)               :: ALPHA_ISOP
      REAL(fp)               :: DTCHEM,    GCO
      REAL(fp)               :: STTCO,     KRATE,      CH4
      REAL(fp)               :: CO_CH4,    CO_ISOP,    CO_MONO
      REAL(fp)               :: CO_CH3OH,  CO_OH,      CO_ACET
      REAL(fp)               :: CH4RATE,   DENS
      REAL(fp)               :: CORATE,    YMID
      REAL(fp)               :: BOXVL

      ! Variables for new JPL 2006 rate constants (jaf, jmao, 3/4/09)
      REAL(fp)               :: KHI1,      KLO1,       XYRAT1
      REAL(fp)               :: BLOG1,     FEXP1,      KHI2
      REAL(fp)               :: KLO2,      XYRAT2,     BLOG2
      REAL(fp)               :: FEXP2,     KCO1,       KCO2

      ! For saving CH4 latitudinal gradient
      REAL(fp),  SAVE        :: A3090S, A0030S, A0030N, A3090N

      ! For OH diurnal cycle (jaf, 3/12/14)
      REAL(fp)               :: OH_MOLEC_CM3, FAC_DIURNAL, SUNCOS

      ! Avoid array temporaries in CHECK_VALUE
      INTEGER                :: ERR_LOC(4)
      CHARACTER(LEN=255)     :: ERR_VAR
      CHARACTER(LEN=255)     :: ERR_MSG

      ! Pointers to fields in State_Met and State_Chm
      REAL(fp),  POINTER     :: AD    (:,:,:  )
      REAL(fp),  POINTER     :: AIRVOL(:,:,:  )
      REAL(fp),  POINTER     :: STT   (:,:,:,:)
      REAL(fp),  POINTER     :: T     (:,:,:  )

      ! To get HEMCO values
      REAL(fp)               :: Emis
      LOGICAL                :: FOUND

      ! Strings
      CHARACTER(LEN=255)     :: LOC='CHEM_TAGGED_CO (tagged_co_mod.F)'
!
! !DEFINED PARAMETERS:
!
      ! Switch to scale yield of isoprene from NOx concentration or not
      LOGICAL,   PARAMETER   :: ALPHA_ISOP_FROM_NOX = .FALSE.

      ! MW of air (g/mol)
      REAL(fp),  PARAMETER   :: WTAIR = 28.966e+0_fp

      ! Yield of CO from CH4
      REAL(fp),  PARAMETER   :: ALPHA_CH4  = 1e+0_fp

      ! Yield of CO from monoterpenes
      REAL(fp),  PARAMETER   :: ALPHA_MONO = 2e-1_fp

      ! Yield of CO from acetone
      REAL(fp),  PARAMETER   :: ALPHA_ACET = 2e+0_fp/3e+0_fp

      !=================================================================
      ! CHEM_TAGGED_CO begins here!
      !
      ! Do the following on the first call to CHEM_TAGGED_CO...
      !=================================================================

      ! Assume success
      RC        = GIGC_SUCCESS

      ! Get fields from Input_Opt
      LSPLIT    = Input_Opt%LSPLIT
      N_TRACERS = Input_Opt%N_TRACERS

      ! DTCHEM is the number of seconds per chemistry timestep
      DTCHEM    = GET_TS_CHEM() * 60e+0_fp

      ! Point to fields in the State_Met and State_Chm objects (bmy, 4/10/14)
      AD        => State_Met%AD
      AIRVOL    => State_Met%AIRVOL
      STT       => State_Chm%Tracers
      T         => State_Met%T

      ! Compute diurnal cycle for OH every day (check for new day inside
      ! subroutine) - jaf 7/10/14
      CALL CALC_DIURNAL

      !=================================================================
      ! Read in OH, NOx, P(CO), and L(CO) fields for the current month
      !=================================================================
      IF ( ITS_A_NEW_MONTH() ) THEN
         
         ! Get current month
         MONTH = GET_MONTH()

         ! Get a pointer to the OH field from the HEMCO list (bmy, 3/11/15)
         CALL HCO_GetPtr( am_I_Root, 'GLOBAL_OH', OH,   RC )
         IF ( RC /= GIGC_SUCCESS ) THEN
            CALL ERROR_STOP ( 'Cannot get pointer to GLOBAL_OH', LOC )
         ENDIF

         ! Read in the loss/production of CO in the stratosphere.
         CALL READ_PCO_LCO_STRAT( am_I_Root, Input_Opt, MONTH, RC )
            
      ENDIF

      !=================================================================
      ! Get the yearly and latitudinal gradients for CH4
      ! This only needs to be called once per year
      !
      ! NOTE: If you are going to run w/ future emissions you must 
      !       pass the future emissions year to GET_GLOBAL_CH4.  
      !       See the modification that was made in "chemdr.f" 
      !       (bmy, 1/24/08)
      !=================================================================
      IF ( ITS_A_NEW_YEAR() ) THEN 

         ! Get current year
         YEAR = GET_YEAR()

         ! Get the semi-hemispheric CH4 gradients
         CALL GET_GLOBAL_CH4( YEAR,  .TRUE.,  A3090S,    A0030S, 
     &                        A0030N, A3090N, am_I_Root, Input_Opt )
      ENDIF

      ! The following is now done in mixing_mod.F90 (ckeller, 3/4/15)
      !=================================================================
      ! Add emissions above the PBL to the tracer array. 
      ! For more information about emissions and PBL mixing, see:
      ! http://wiki.geos-chem.org/Distributing_emissions_in_the_PBL
      !=================================================================
!      DO N = 1, N_TRACERS
!      DO L = 1, LLPAR
!      DO J = 1, JJPAR
!      DO I = 1, IIPAR
!
!         ! If we are using the non-local PBL mixing scheme, then
!         ! only add emissions above the PBL top.  All emissions within
!         ! the PBL are handled in routine VDIFFDR in vdiff_mod.F90.
!         ! Skip if we have PBL mixing completely turned off, or if
!         ! we are using the TURBDAY full-mixing scheme.
!         IF ( Input_Opt%LNLPBL ) THEN
!            IF ( L < State_Met%PBL_TOP_L(I,J) + 1 ) CYCLE
!         ENDIF
!
!         ! Add emissions into the tracer array
!         ! Convert from [kg/m2/s] to [kg]
!         ! Now use GetHcoVal (ckeller, 10/21/2014)
!         CALL GetHcoVal( N, I, J, L, FOUND, Emis=Emis )
!         IF ( .NOT. FOUND ) CYCLE
!         State_Chm%TRACERS(I,J,L,N) = State_Chm%TRACERS(I,J,L,N) 
!     &                              + ( Emis
!     &                              *   State_Met%AREA_M2(I,J,1)
!     &                              *   DTCHEM                  )
!
!      ENDDO
!      ENDDO
!      ENDDO
!      ENDDO

      !=================================================================
      ! Do tagged CO chemistry -- Put everything within a large 
      ! DO-loop over all grid boxes to facilitate parallelization
      !=================================================================
!$OMP PARALLEL DO 
!$OMP+DEFAULT( SHARED                                             ) 
!$OMP+PRIVATE( I, J, L, N, STTCO, GCO, DENS, CH4RATE, CO_CH4, CH4 )
!$OMP+PRIVATE( KRATE, CO_ISOP, CO_CH3OH, CO_MONO, CO_ACET, CORATE )
!$OMP+PRIVATE( CO_OH, YMID, ALPHA_ISOP                            )
!$OMP+PRIVATE( KHI1, KLO1, XYRAT1, BLOG1, FEXP1, KHI2, KLO2       )
!$OMP+PRIVATE( XYRAT2, BLOG2, FEXP2, KCO1, KCO2                   )
!$OMP+PRIVATE( OH_MOLEC_CM3, FAC_DIURNAL, SUNCOS                  )
!$OMP+PRIVATE( ERR_LOC, ERR_VAR, ERR_MSG, BOXVL )
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         
         ! Latitude of grid box
         YMID     = GET_YMID( I, J, L )

         ! Grid box volume [cm3]
         BOXVL    = State_Met%AIRVOL(I,J,L) * 1e+6_fp

         !==============================================================
         ! (0) Define useful quantities
         !==============================================================

         ! STTCO [molec CO/cm3/kg CO] converts [kg CO] --> [molec CO/cm3]
         ! kg CO/box * box/cm3 * mole/0.028 kg CO * Avog.#/mole
         STTCO = 1e+0_fp     / AIRVOL(I,J,L) / 1e+6_fp /
     &           FMOL_CO * 6.023e+23_fp

         ! GCO is CO concentration in [molec CO/cm3]
         GCO   = STT(I,J,L,1) * STTCO
         
         ! DENS is the number density of air [molec air/cm3]
         DENS  = AD(I,J,L) * 1000.e+0_fp   /
     &           BOXVL     * 6.023e+23_fp  / WTAIR

         ! Cosine of the solar zenith angle [unitless]
         SUNCOS = State_Met%SUNCOSmid(I,J)

         ! Scaling factor for diurnal cycles - zero at night
         IF ( SUNCOS > 0e+0_fp .and. TCOSZ(I,J) > 0e+0_fp ) THEN
            FAC_DIURNAL = ( SUNCOS     / TCOSZ(I,J)    ) *
     &                    ( 1440e+0_fp / GET_TS_CHEM() )
         ELSE
            FAC_DIURNAL = 0e+0_fp
         ENDIF

         ! Now impose a diurnal cycle on OH.
         ! This is done in other offline simulations but was
         ! missing from tagged CO (jaf, 3/12/14)
         OH_MOLEC_CM3 = OH(I,J,L) * FAC_DIURNAL

         ! Make sure OH is not negative
         OH_MOLEC_CM3 = MAX( OH_MOLEC_CM3, 0e+0_fp )

         !==============================================================
         ! (1a) Production of CO by reaction with CH4
         !==============================================================

         ! Initialize
         CO_CH4 = 0e+0_fp

         ! Test level for stratosphere or troposphere
         IF ( ITS_IN_THE_STRAT( I, J, L, State_Met ) ) THEN

            !===========================================================
            ! (1a-1) Production of CO from CH4 in the stratosphere 
            !===========================================================

            ! Call GET_PCO_LCO_STRAT to get the P(CO) rate from CH4
            CH4RATE = GET_PCO_LCO_STRAT( .TRUE., I, J, L )
            
            ! Convert units of CH4RATE from [v/v/s] to [molec CO/cm3]
            CO_CH4 = CH4RATE * DTCHEM * DENS

         ELSE
            
            !===========================================================
            ! (1a-2) Production of CO from CH4 in the troposphere 
            !===========================================================

            ! CH4 concentration [ppbv] for the given latitude band 
            ! (bmy, 1/2/01)
            CH4 = A3090S
            IF ( YMID >= -30.0 .and. YMID < 0.0  ) CH4 = A0030S
            IF ( YMID >=   0.0 .and. YMID < 30.0 ) CH4 = A0030N
            IF ( YMID >=  30.0                   ) CH4 = A3090N
            
            ! Convert CH4 from [ppbv] to [molec CH4/cm3]
            CH4 = CH4 * 1e-9_fp * DENS

            ! Calculate updated rate constant [s-1] (bnd, bmy, 1/2/01)
            KRATE = 2.45e-12_fp * EXP( -1775.e+0_fp / T(I,J,L) ) 

            ! Production of CO from CH4 = alpha * k * [CH4] * [OH] * dt 
            ! Units are [molec CO/cm3]
            CO_CH4 = ALPHA_CH4 * KRATE * CH4 * OH_MOLEC_CM3 * DTCHEM 

         ENDIF

         ! Check CO_CH4 for NaN or Infinity
         ERR_LOC = (/ I, J, L, 0 /)
         ERR_VAR = 'CO_CH4'
         ERR_MSG = 'STOP at tagged_co_mod:1'
         CALL CHECK_VALUE( CO_CH4, ERR_LOC, ERR_VAR, ERR_MSG )

         !==============================================================
         ! (1b) Production of CO from ISOPRENE and METHANOL (CH3OH)
         !==============================================================

         ! Initialize
         CO_ISOP  = 0e+0_fp
         CO_CH3OH = 0e+0_fp

         ! Isoprene is emitted only into the surface layer 
         IF ( L == 1 ) THEN 

            !===========================================================
            ! Yield of CO from ISOP: 30%, from Miyoshi et al., 1994.
            ! They estimate globally 105 Tg C/yr of CO is produced 
            ! from isoprene oxidation. 
            !-----------------------------------------------------------
            ! We need to scale the Isoprene flux to get the CH3OH 
            ! (methanol) flux.  Currently, the annual isoprene flux in 
            ! GEOS-CHEM is ~ 397 Tg C.
            !
            ! Daniel Jacob recommends a flux of 100 Tg/yr CO from CH3OH 
            ! oxidation based on Singh et al. 2000 [JGR 105, 3795-3805] 
            ! who estimate a global methanol source of 122 Tg yr-1, of 
            ! which most (75 Tg yr-1) is "primary biogenic".  He also 
            ! recommends for now that the CO flux from CH3OH oxidation 
            ! be scaled to monthly mean isoprene flux.
            !
            ! To get CO from METHANOL oxidation, we must therefore 
            ! multiply the ISOPRENE flux by the following scale factor:
            !  ( 100 Tg CO / 397 Tg C ) * ( 12 g C/mole / 28 g CO/mole )
            !-----------------------------------------------------------
            ! We now call GET_ALPHA_ISOP to get the yield factor of
            ! CO produced from isoprene, as a function of NOx, or
            ! as a constant. (bnd, bmy, 6/14/01)
            !===========================================================

            ! Get CO yield from isoprene
            ALPHA_ISOP = GET_ALPHA_ISOP( .FALSE. )

            ! P(CO) from Isoprene Flux = ALPHA_ISOP * Flux(ISOP)
            ! Convert from [molec ISOP/box] to [molec CO/cm3]
            !
            ! Units of SUMISOPCO are [atoms C/box/time step]. 
            ! Division by 5 is necessary to convert to 
            ! [molec ISOP/box/timestep].
            !
            ! Units of ALPHA_ISOP are [molec CO/molec ISOP]
            ! Units of CO_ISOP are [molec CO/cm3]
            CO_ISOP = SUMISOPCO(I,J) / BOXVL      /
     &                5.e+0_fp           * ALPHA_ISOP

            ! P(CO) from CH3OH is scaled to Isoprene Flux (see above)
            ! Units are [molec CO/cm3]
            CO_CH3OH = ( SUMISOPCO(I,J) / BOXVL ) * 
     &                 ( 100e+0_fp          / 397e+0_fp ) * 
     &                 ( 12e+0_fp           / 28e+0_fp  )
 
            ! Zero SUMISOPCO and SUMCH3OHCO for the next emission step
            SUMISOPCO(I,J)  = 0e+0_fp
            SUMCH3OHCO(I,J) = 0e+0_fp

            ! Check CO_ISOP for NaN or Infinity
            ERR_LOC = (/ I, J, L, 0 /)
            ERR_VAR = 'CO_ISOP'
            ERR_MSG = 'STOP at tagged_co_mod:2'
            CALL CHECK_VALUE( CO_ISOP,  ERR_LOC, ERR_VAR, ERR_MSG )

            ! Check CO_CH3OH for NaN or Infinity
            ERR_VAR = 'CO_CH3OH'
            ERR_MSG = 'STOP at tagged_co_mod:3'
            CALL CHECK_VALUE( CO_CH3OH, ERR_LOC, ERR_VAR, ERR_MSG )

         ENDIF

         !==============================================================
         ! (1c) Production of CO from MONOTERPENE oxidation
         !==============================================================

         ! Initialize
         CO_MONO = 0.e+0_fp

         ! Monoterpenes are emitted only into the surface layer
         IF ( L == 1 ) THEN

            !===========================================================
            ! Assume the production of CO from monoterpenes is 
            ! instantaneous even though the lifetime of intermediate 
            ! species may be on the order of hours or days.  This 
            ! assumption will likely cause CO from monoterpene 
            ! oxidation to be too high in the box in which the 
            ! monoterpene is emitted.
            !-----------------------------------------------------------
            ! The CO yield here is taken from:
            !   Hatakeyama et al. JGR, Vol. 96, p. 947-958 (1991)
            !   Vinckier et al. Fresenius Env. Bull., Vol. 7, p.361-368 
            !     (1998)
            !
            ! Hatakeyama:  "The ultimate yield of CO from the 
            !   tropospheric oxidation of terpenes (including both O3 
            !   and OH reactions) was estimated to be 20% on the carbon 
            !   number basis."  They studied ALPHA- & BETA-pinene.
            !
            ! Vinckier  :  "R(CO)=1.8+/-0.3" : 1.8/10 is about 20%.
            !-----------------------------------------------------------
            ! Calculate source of CO per time step from monoterpene 
            ! flux (assume lifetime very short) using the C number basis:
            !
            !   CO [molec CO/cm3] = Flux [atoms C from MONO/box] /
            !                       Grid Box Volume [cm^-3]       *
            !                       ALPHA_MONO 
            !
            ! where ALPHA_MONO = 0.2 as explained above.
            !===========================================================

            ! P(CO) from Monoterpene Flux =  alpha * Flux(Mono)
            ! Units are [molec CO/cm3]
            CO_MONO = ( SUMMONOCO(I,J) / BOXVL ) *
     &                ALPHA_MONO

            ! Zero SUMMONOCO for the next emission step
            SUMMONOCO(I,J) = 0e+0_fp

            ! Check CO_MONO for NaN or Infinity
            ERR_LOC = (/ I, J, L, 0 /)
            ERR_VAR = 'CO_MONO'
            ERR_MSG = 'STOP at tagged_co_mod:4'
            CALL CHECK_VALUE( CO_MONO, ERR_LOC, ERR_VAR, ERR_MSG )

         ENDIF

         !==============================================================
         ! (1d) Production of CO from oxidation of ACETONE 
         !
         ! ALPHA_ACET = 2/3 to get a yield for CO.  This accounts 
         ! for acetone loss from reaction with OH And photolysis.  
         ! The acetone sources taken into account are:
         ! 
         ! (a) Primary emissions of acetone from biogenic sources
         ! (b) Secondary production of acetone from monoterpene 
         !      oxidation
         ! (c) Secondary production of acetone from ALK4 and 
         !      propane oxidation
         ! (d) Direct emissions of acetone from biomass burning and 
         !      fossil fuels
         ! (e) direct emissions from ocean
         !
         ! Calculate source of CO per time step from biogenic acetone 
         ! # molec CO/cc = ALPHA * ACET Emission Rate * dt
         !==============================================================

         ! Initialize
         CO_ACET = 0.e+0_fp

         ! Biogenic acetone sources are emitted only into the surface layer
         IF ( L == 1 ) THEN

            ! Units are [molec CO/cc]
            CO_ACET = SUMACETCO(I,J) / BOXVL * ALPHA_ACET

            ! Zero SUMACETCO for the next emission step
            SUMACETCO(I,J) = 0e+0_fp

            ! Check CO_ACET for NaN or Infinity
            ERR_LOC = (/ I, J, L, 0 /)
            ERR_VAR = 'CO_ACET'
            ERR_MSG = 'STOP at tagged_co_mod:5'
            CALL CHECK_VALUE( CO_ACET, ERR_LOC, ERR_VAR, ERR_MSG )

         ENDIF

         !==============================================================
         ! (1e) Add production of CO into the following tagged tracers:
         !
         ! (a) Tracer #12: CO produced from CH4
         ! (b) Tracer #14: CO produced from ISOPRENE
         ! (c) Tracer #15: CO produced from MONOTERPENES
         ! (d) Tracer #16: CO produced from METHANOL (CH3OH)
         ! (e) Tracer #17: CO produced from ACETONE
         !
         ! %%% NOTE: If you are modifying the tagged CO simulation, 
         ! %%% and your simulation has less than 12 tracers, then 
         ! %%% then comment out this section.  If you don't you can
         ! %%% get an array-out-of-bounds error (bmy, 6/11/08)
         !==============================================================
         IF ( LSPLIT ) THEN
            STT(I,J,L,12) = STT(I,J,L,12) + CO_CH4   / STTCO
            STT(I,J,L,14) = STT(I,J,L,14) + CO_ISOP  / STTCO 
            STT(I,J,L,15) = STT(I,J,L,15) + CO_MONO  / STTCO
            STT(I,J,L,16) = STT(I,J,L,16) + CO_CH3OH / STTCO
            STT(I,J,L,17) = STT(I,J,L,17) + CO_ACET  / STTCO
         ENDIF

         !==============================================================
         ! (2a) Loss of CO due to chemical reaction w/ OH
         !==============================================================

         ! Select out tropospheric or stratospheric boxes
         IF ( ITS_IN_THE_STRAT( I, J, L, State_Met ) ) THEN

            !===========================================================
            ! (2a-1) Stratospheric loss of CO due to chemical rxn w/ OH
            !===========================================================

            ! Get the L(CO) rate in the stratosphere in [s-1]
            CORATE = GET_PCO_LCO_STRAT( .FALSE., I, J, L )

            ! CO_OH is the fraction of CO lost to OH [unitless]
            CO_OH = CORATE * DTCHEM

            ! Check CO_OH for NaN or Infinity
            ERR_LOC = (/ I, J, L, 0 /)
            ERR_VAR = 'CO_OH'
            ERR_MSG = 'STOP at tagged_co_mod:6'
            CALL CHECK_VALUE( CO_OH, ERR_LOC, ERR_VAR, ERR_MSG )

            ! Handle strat loss by OH for regional CO tracers
            IF ( LSPLIT ) THEN

               ! Loop over regional CO tracers
               DO N = 2, N_TRACERS

                  ! Loss
                  STT(I,J,L,N) = STT(I,J,L,N) * ( 1e+0_fp - CO_OH )

                  ! STT shouldn't be less than zero
                  IF ( STT(I,J,L,N) < 0e+0_fp ) STT(I,J,L,N) = 0e+0_fp

                  ! Error check
                  ERR_LOC = (/ I, J, L, N /)
                  ERR_VAR = 'STT'
                  ERR_MSG = 'STOP at tagged_co_mod:7'
                  CALL CHECK_VALUE( STT(I,J,L,N), ERR_LOC,
     &                              ERR_VAR,      ERR_MSG )

                  ! ND65 diagnostic -- loss of CO by OH for "tagged" tracers
                  IF ( ND65 > 0 .and. L <= LD65 ) THEN
                     AD65(I,J,L,N) = AD65(I,J,L,N) +
     &                               ( CORATE * STT(I,J,L,N) * STTCO )
                  ENDIF
               ENDDO
            ENDIF

            ! CO_OH above is just the fraction of CO lost by OH.  Here
            ! we multiply it by GCO (the initial value of STT in molec/cm3)
            ! to convert it to an amount of CO lost by OH [molec/cm3]
            ! (bmy, 2/19/02)
            CO_OH = GCO * CO_OH

         ELSE
   
            !===========================================================
            ! (2a-2) Tropospheric loss of CO due to chemical rxn w/ OH
            !
            !  DECAY RATE
            !  The decay rate (KRATE) is calculated by:
            !
            !     No change from JPL '97 to JPL '03 (jaf, 2/27/09)
            !     OH + CO -> products (JPL '03)
            !     k = (1 + 0.6Patm) * 1.5E-13
            !
            !  KRATE has units of [ molec^2 CO / cm6 / s ]^-1, 
            !  since this is a 2-body reaction.
            !
            ! Updated rate constant from JPL 2006; now more complicated.
            ! From JPL 2006: "The  reaction between HO and CO to yield
            ! H + CO2 akes place on a potential energy surface that 
            ! contains the radical HOCO.  The yield of H and CO2 is 
            ! diminished as the pressure rises.  The loss of reactants 
            ! is thus the sum of two processes, an association to yield 
            ! HOCO and the chemical activation process yielding H and
            ! CO2." So we now need two complicated reactions. The code
            ! is more or less copied from calcrate.f as implemented by
            ! jmao (jaf, 3/4/09)
            !===========================================================

            ! Decay rate - new JPL 2006 version (jaf, 3/4/09)
            ! KLO1 = k_0(T) from JPL Data Eval (page 2-1)
            KLO1   = 5.9e-33_fp * ( T(I,J,L) / 300 )**(-1.4e+0_fp)
            ! KHI1 = k_inf(T) from JPL Data Eval (page 2-1)
            KHI1   = 1.1e-12_fp * ( T(I,J,L) / 300 )**(1.3e+0_fp)
            XYRAT1 = KLO1 * DENS / KHI1
            BLOG1  = LOG10(XYRAT1)
            FEXP1  = 1.e+0_fp / ( 1.e+0_fp + BLOG1 * BLOG1 )
            ! KCO1 = k_f([M],T) from JPL Data Eval (page 2-1)
            KCO1   = KLO1 * DENS * 0.6**FEXP1 / ( 1.e+0_fp + XYRAT1 )
            ! KLO2 = k_0(T) from JPL Data Eval (page 2-1)
            KLO2   = 1.5e-13_fp * ( T(I,J,L) / 300 )**(0.6e+0_fp)
            ! KHI2 = k_inf(T) from JPL Data Eval (page 2-1)
            KHI2   = 2.1e+09_fp  * ( T(I,J,L) / 300 )**(6.1e+0_fp)
            XYRAT2 = KLO2 * DENS / KHI2
            BLOG2  = LOG10(XYRAT2)
            FEXP2  = 1.e+0_fp / ( 1.e+0_fp + BLOG2 * BLOG2 )
            ! KCO2 = k_f^ca([M],T) from JPL Data Eval (page 2-2)
            KCO2   = KLO2 * 0.6**FEXP2 / ( 1.e+0_fp + XYRAT2 )

            ! KRATE is the sum of the two.
            KRATE  = KCO1 + KCO2

            ! CO_OH = Tropospheric loss of CO by OH [molec/cm3]
            ! Now use OH_MOLEC_CM3, which includes a diurnal cycle.
            CO_OH = KRATE * GCO * OH_MOLEC_CM3 * DTCHEM

            ! Handle trop loss by OH for regional CO tracers
            IF ( LSPLIT ) THEN

               ! Loop over regional CO tracers
               DO N = 2, N_TRACERS

                  ! Use tropospheric rate constant 
                  STT(I,J,L,N) = STT(I,J,L,N) *
     &                 ( 1e+0_fp - KRATE * OH_MOLEC_CM3 * DTCHEM )

                  ! Error check
                  ERR_LOC = (/ I, J, L, N /)
                  ERR_VAR = 'STT'
                  ERR_MSG = 'STOP at tagged_co_mod:8'
                  CALL CHECK_VALUE( STT(I,J,L,N), ERR_LOC,
     &                              ERR_VAR,      ERR_MSG )

                  ! ND65 diagnostic -- loss of CO by OH for "tagged" tracers
                  IF ( ND65 > 0 .and. L <= LD65 ) THEN
                     AD65(I,J,L,N) = AD65(I,J,L,N) + 
     &                   ( KRATE * OH_MOLEC_CM3 * STT(I,J,L,N) * STTCO )
                  ENDIF
               ENDDO
            ENDIF

         ENDIF
         
         !==============================================================
         ! Save the total chemical production from various sources
         ! into the total CO tracer STT(I,J,L,1)
         !==============================================================

         ! GCO is the total CO before chemistry was applied [molec CO/cm3]
         ! Add to GCO the sources and sinks listed above
         GCO = GCO + CO_CH4   + CO_MONO + CO_ACET + 
     &               CO_CH3OH + CO_ISOP - CO_OH

         ! Convert net CO from [molec CO/cm3] to [kg] and store in STT
         STT(I,J,L,1) = GCO / STTCO

         !==============================================================
         ! Archive ND65 diagnostics -- Production & Loss of CO
         ! Also save P(CO) from CH3OH and MONOTERPENES (bmy, 1/2/01)
         !==============================================================
         IF ( ND65 > 0 .and. L <= LD65 ) THEN

            ! Loss of CO by OH (global) [s-1]
            N             = 1
            AD65(I,J,L,N) = AD65(I,J,L,N) + ( CO_OH / DTCHEM )

            ! Production of CO from Isoprene [molec CO/cm3/s]
            N             = N_TRACERS + 1
            AD65(I,J,L,N) = AD65(I,J,L,N) + ( CO_ISOP / DTCHEM )

            ! Production of CO from CH4 [molec CO/cm3/s]
            N             = N_TRACERS + 2
            AD65(I,J,L,N) = AD65(I,J,L,N) + ( CO_CH4 / DTCHEM )

            ! Production of CO from CH3OH [molec CO/cm3/s]            
            N             = N_TRACERS + 3
            AD65(I,J,L,N) = AD65(I,J,L,N) + ( CO_CH3OH / DTCHEM )

            ! Production of CO from MONO [molec CO/cm3/s]
            N             = N_TRACERS + 4
            AD65(I,J,L,N) = AD65(I,J,L,N) + ( CO_MONO / DTCHEM )    

            ! Production of CO from ACET [molec CO/cm3/s]
            N             = N_TRACERS + 5 
            AD65(I,J,L,N) = AD65(I,J,L,N) + ( CO_ACET / DTCHEM )
         ENDIF
      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Free pointers
      NULLIFY( AD, AIRVOL, STT, T )

      END SUBROUTINE CHEM_TAGGED_CO
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: calc_diurnal
!
! !DESCRIPTION: Subroutine CALC\_DIRUNAL computes the sume of the cosine
!  of the solar zenith angle over a 24 hour day as well as the total
!  length of daylight to scale the offline OH concentrations.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CALC_DIURNAL
!
! !USES:
!
      USE CMN_GCTM_MOD
      USE CMN_SIZE_MOD
      USE GRID_MOD, ONLY : GET_YMID_R
      USE TIME_MOD, ONLY : ITS_A_NEW_DAY, GET_MINUTE,      GET_HOUR
      USE TIME_MOD, ONLY : GET_TS_CHEM,   GET_DAY_OF_YEAR, GET_LOCALTIME
! 
! !REVISION HISTORY: 
!  12 Mar 2014 - J. Fisher   - Copied from OHNO3TIME in carbon_mod and
!                              COSSZA in dao_mod
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL, SAVE       :: FIRST = .TRUE.
      INTEGER             :: I, J, N, NDYSTEP
      INTEGER             :: MINUTE,  TS_SUN
      REAL*8              :: GMT_MID, TIMLOC, FACTOR
      REAL*8              :: R,       AHR,    DEC
      REAL*8              :: YMID_R,  SUNTMP_MID
!
! !DEFINED PARAMETERS:
!
      ! Coefficients for solar declination angle
      REAL*8,  PARAMETER :: A0 = 0.006918d0
      REAL*8,  PARAMETER :: A1 = 0.399912d0
      REAL*8,  PARAMETER :: A2 = 0.006758d0
      REAL*8,  PARAMETER :: A3 = 0.002697d0
      REAL*8,  PARAMETER :: B1 = 0.070257d0
      REAL*8,  PARAMETER :: B2 = 0.000907d0
      REAL*8,  PARAMETER :: B3 = 0.000148d0
 
      !=================================================================
      ! CALC_DIURNAL begins here!
      !=================================================================

      ! Only do at the start of a new day
      IF ( FIRST .or. ITS_A_NEW_DAY() ) THEN 
 
         ! Zero array
         TCOSZ = 0d0

         ! Get time for central chemistry timestep
         TS_SUN = GET_TS_CHEM()                     ! Chemistry interval
         MINUTE = GET_MINUTE()                      ! Current minutes
         FACTOR = MINUTE / TS_SUN                   ! Multiplying factor

         ! GMT at the midpoint of the chemistry time interval for first
         ! timestep of the day
         GMT_MID  = ( DBLE( GET_HOUR()        )      )
     &            + ( DBLE( TS_SUN * FACTOR ) / 60d0 )
     &            + ( DBLE( TS_SUN / 2      ) / 60d0 )

         ! Solar declination angle (low precision formula):
         ! Path length of earth's orbit traversed since Jan 1 [radians]
         R = ( 2d0 * PI / 365d0 ) * FLOAT( GET_DAY_OF_YEAR() - 1 )
         DEC = A0 - A1*COS(    R) + B1*SIN(    R)
     &            - A2*COS(2d0*R) + B2*SIN(2d0*R)
     &            - A3*COS(3d0*R) + B3*SIN(3d0*R)

         ! NDYSTEP is # of chemistry time steps
         NDYSTEP = INT( 24d0 * 60d0 / GET_TS_CHEM() )

         ! Loop forward through NDYSTEP "fake" timesteps for this day
         DO N = 1, NDYSTEP
            
            ! Increment GMT (hours) to midpoint of next timestep
            IF ( N > 1 ) GMT_MID = GMT_MID + TS_SUN / 60d0

            ! Loop over surface grid boxes
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, YMID_R, TIMLOC, AHR, SUNTMP_MID )
            DO J = 1, JJPAR
            DO I = 1, IIPAR

               ! Zero SUNTMP_MID
               SUNTMP_MID = 0d0

               ! Grid box latitude center [radians]
               YMID_R = GET_YMID_R( I, J, 1 )

               ! Local time at box (I,J) [hours]
               TIMLOC = GET_LOCALTIME( I, J, 1, GMT=GMT_MID )

               ! Hour angle at box (I,J) [radians]
               AHR = ABS( TIMLOC - 12d0 ) * 15d00 * PI_180

               !===========================================================
               ! The cosine of the solar zenith angle (SZA) is given by:
               !
               !  cos(SZA) = sin(LAT)*sin(DEC) + cos(LAT)*cos(DEC)*cos(AHR)
               !
               ! where LAT = the latitude angle,
               !       DEC = the solar declination angle,
               !       AHR = the hour angle, all in radians.
               !
               ! If SUNCOS < 0, then the sun is below the horizon, and
               ! therefore does not contribute to any solar heating.
               !===========================================================

               ! Compute Cos(SZA)
               SUNTMP_MID = sin(YMID_R) * sin(DEC) +
     &                      cos(YMID_R) * cos(DEC) * cos(AHR)

               ! TCOSZ is the sum of SUNTMP_MID at location (I,J)
               ! Do not include negative values of SUNTMP_MID
               TCOSZ(I,J) = TCOSZ(I,J) + MAX( SUNTMP_MID, 0d0 )

            ENDDO
            ENDDO
!$OMP END PARALLEL DO
         ENDDO

         ! Reset first-time flag
         FIRST = .FALSE.
      ENDIF

      ! Return to calling program
      END SUBROUTINE CALC_DIURNAL

!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_alpha_isop
!
! !DESCRIPTION: Function GET\_ALPHA\_ISOP returns the CO yield from 
!  Isoprene (ALPHA\_ISOP) either as a function of NOx or as a constant. 
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_ALPHA_ISOP( FROM_NOX, NOX ) RESULT( ALPHA_ISOP )
!
! !USES:
!
      USE ERROR_MOD, ONLY : ERROR_STOP
!
! !INPUT PARAMETERS: 
!
      LOGICAL, INTENT(IN)           :: FROM_NOX    ! If =T, will take ALPHA_ISOP
                                                   ! as f(NOx).  If =F, will set
                                                   ! ALPHA_ISOP to a constant
      REAL(fp),  INTENT(IN), OPTIONAL :: NOX         ! NOx concentration [ppbv]
!
! !RETURN VALUE:
!
      REAL(fp)                        :: ALPHA_ISOP  ! CO yield from ISOP
! 
! !REVISION HISTORY: 
!  13 Jun 2001 - B. Duncan - Initial version
!  01 Oct 1995 - R. Yantosca - Initial version
!  (1 ) Now make NOx an optional argument (bmy, 8/28/01)
!  (2 ) Now reference ERROR_STOP from "error_mod.f" (bmy, 10/15/02)
!  (3 ) Updated comments (bmy, 7/20/04)
!  23 Oct 2012 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
      !=================================================================
      ! If FROM_NOX = T, then we are taking the CO yield from Isoprene 
      ! as a function of NOx.  Otherwise, we set it to a constant.
      !=================================================================
      IF ( FROM_NOX ) THEN

         ! Make sure NOX is passed if FROM_NOX = TRUE
         IF ( PRESENT( NOX ) ) THEN

            ! The CO yield here is taken from acf's calculation 
            ! (from group meeting (1-20-99).  Assumming linearity, 
            ! ALPHA=0.8*[NOx]+0.6, with an upper and lower limit of 2.1 
            ! and 0.8, respectively. [NOx] is concentration in ppbv !!
            ALPHA_ISOP = ( 0.8e+0_fp * NOX ) + 0.6e+0_fp

            ! Force lower & upper limits
            IF ( NOX < 0.5e+0_fp ) ALPHA_ISOP = 0.8e+0_fp
            IF ( NOX > 1.8e+0_fp ) ALPHA_ISOP = 2.1e+0_fp
         
         ELSE

            ! Stop w/ error message
            CALL ERROR_STOP( 'NOx argument not passed!', 
     &                       'GET_ALPHA_ISOP (tagged_co_mod.f)' )
           
         ENDIF
   
      ELSE

         ! Otherwise, use a 30% yield from Miyoshi et al., 1994.
         ! They estimate globally 105 Tg C/yr of CO is produced
         ! from isoprene oxidation.
         ! ALPHA_ISOP = (0.3 molec CO/atoms C) x (5 atoms C/molec ISOP)
         ALPHA_ISOP = 1.5e+0_fp

      ENDIF
      
      END FUNCTION GET_ALPHA_ISOP
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_pco_lco_strat
!
! !DESCRIPTION: Subroutine READ\_PCO\_LCO\_STRAT reads production and 
!  destruction rates for CO in the stratosphere. (bnd, bmy, 9/13/00, 8/13/10) 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE READ_PCO_LCO_STRAT( am_I_Root, Input_Opt, 
     &                               THISMONTH, RC         ) 
!
! !USES:
!
      ! GEOS-Chem modules
      USE CMN_SIZE_MOD   
      USE BPCH2_MOD,          ONLY : GET_NAME_EXT, GET_RES_EXT
      USE TRANSFER_MOD,       ONLY : TRANSFER_3D
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput

      ! netCDF modules
      USE m_netcdf_io_open
      USE m_netcdf_io_read
      USE m_netcdf_io_close
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)  :: am_I_Root   ! Are we on the root CPU?
      TYPE(OptInput), INTENT(IN)  :: Input_Opt   ! Input Options object
      INTEGER,        INTENT(IN)  :: THISMONTH   ! Month to read
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT) :: RC          ! Success or failure?!
! 
! !REMARKS:
!  ##########################################################################
!  #####    NOTE: BINARY PUNCH INPUT IS BEING PHASED OUT.  THIS DATA    #####
!  #####    WILL EVENTUALLY BE READ IN FROM netCDF FILES VIA HEMCO!     #####
!  #####       -- Bob Yantosca (05 Mar 2015)                            #####
!  ##########################################################################

! !REVISION HISTORY: 
!  13 Sep 2000 - B. Duncan   - Initial version
!  (1 ) Now use IOS /= 0 to trap both I/O errors and EOF. (bmy, 9/13/00)
!  (2 ) Added to module "tagged_co_mod.f" (bmy, 6/13/01)
!  (3 ) ARRAY needs to be of size (IIPAR,JJPAR).  Use TRANSFER_ZONAL from
!        "transfer_mod.f" to cast data from REAL*4 to REAL(fp), and also to
!        copy data to an array of size (JJPAR,LLPAR).  Use 3 arguments (M/D/Y)
!        in call to GET_TAU0.  Use JJPAR,LGLOB in call to READ_BPCH2. 
!        (bmy, 9/28/01)
!  (4 ) Removed obsolete code from 9/28/01 (bmy, 10/22/01)
!  (5 ) Updated comments (bmy, 2/15/02)
!  (6 ) Update FILENAME so that it looks in the "pco_lco_200203" subdirectory
!        of DATA_DIR. (bnd, bmy, 6/30/03)
!  (7 ) Now references DATA_DIR from "directory_mod.f" (bmy, 7/20/04)
!  (8 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  13 Aug 2010 - R. Yantosca - Treat MERRA in the same way as GEOS-5
!  08 Feb 2012 - R. Yantosca - Treat GEOS-5.7.2 in the same way as MERRA
!  08 Dec 2009 - R. Yantosca - Added ProTeX headers
!  20 Jun 2014 - R. Yantosca - Now accept am_I_Root, Input_Opt, RC
!  05 Mar 2015 - R. Yantosca - Add Input_Opt%RES_DIR to data path
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER             :: fileID
      CHARACTER(LEN=255)  :: FILENAME
      REAL*4              :: ARRAY(IIPAR,JJPAR,LGLOB) 

      !=================================================================
      ! READ_PCO_LCO_STRAT begins here!
      !=================================================================

      ! Assume success
      RC        =  GIGC_SUCCESS

      ! Initialize some variables
      ARRAY    = 0e0
      CO_PRODS = 0e+0_fp
      CO_LOSSS = 0e+0_fp

      ! Construct filename and open file
      FILENAME = 'strat_chem_201206/gmi.clim.CO.' // 
     &                GET_NAME_EXT() // '.' // GET_RES_EXT() // '.nc'
      FILENAME = TRIM( Input_Opt%DATA_DIR ) // 
     &           TRIM( Input_Opt%RES_DIR  ) // TRIM( FILENAME )

      call NcOp_Rd( fileID, TRIM( FILENAME ) )

      !=================================================================
      ! Read in CO production rates [v/v/s]
      !=================================================================

      ! Read P(CO)
      call NcRd( array, fileID, 'prod',
     &           (/     1,     1,     1, THISMONTH /),   ! Start
     &           (/ IIPAR, JJPAR, LGLOB,         1 /)  ) ! Count

      ! Cast from REAL*4 to REAL(fp) and resize to 1:LLPAR
      CALL TRANSFER_3D( ARRAY, CO_PRODS )

      !=================================================================
      ! Read in CO loss frequency [s^-1]
      !=================================================================

      ! Read P(CO)
      call NcRd( array, fileID, 'loss',
     &           (/     1,     1,     1, THISMONTH /),   ! Start
     &           (/ IIPAR, JJPAR, LGLOB,         1 /)  ) ! Count

      ! Cast from REAL*4 to REAL(fp) and resize to 1:LLPAR
      CALL TRANSFER_3D( ARRAY, CO_LOSSS )

      ! Close file
      call NcCl( fileID )

      END SUBROUTINE READ_PCO_LCO_STRAT
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_pco_lco_strat
!
! !DESCRIPTION: Function GET\_CO\_STRAT\_RATE uses production and loss rates
!  for CO to calculate net production of CO in the stratosphere.  The 
!  purpose of this SR is to prevent high CO concentrations from building 
!  up in the stratosphere; in these layers only transport is simulated 
!  (i.e., no chemistry).  For a long simulation, a buildup of high 
!  concentrations could occur causing the stratosphere to become a 
!  significant source of CO. 
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_PCO_LCO_STRAT( IS_PROD, I, J, L ) RESULT( RATE )
!
! !USES:
!
      USE CMN_SIZE_MOD                  ! Size parameters

!
! !INPUT PARAMETERS: 
!
      LOGICAL, INTENT(IN)  :: IS_PROD   ! If = T, return CO production
                                        ! If = F, return CO loss
      INTEGER, INTENT(IN)  :: I, J, L   ! GEOS-Chem grid box indices
!
! !RETURN VALUE:
!
      REAL(fp)               :: RATE      ! CO prod [v/v/s] or loss [1/s]
! 
! !REVISION HISTORY: 
!  13 Jun 2001 - B. Duncan - Initial version
!  (1 ) Production (mixing ratio/sec) and loss (1/sec) rates provided by 
!        Dylan Jones.  Only production by CH4+OH and destruction by CO+OH 
!        are considered.
!  (2 ) The annual mean tropopause is stored in the LPAUSE array 
!       (from header file "CMN").  LPAUSE is defined such that: 
!       Levels            1 <= L <= LPAUSE(I,J) - 1 are tropospheric
!               LPAUSE(I,J) <= L <= LLPAR           are stratospheric
!  (3 ) LPAUSE_MIN = minimun tropopause height.  Start L-loop from the 
!        lowest stratospheric level!
!  (4 ) Added to module "tagged_co_mod.f" (bmy, 6/18/01)
!  (5 ) Updated comments (bmy, 2/19/02)
!  (6 ) Removed reference to CMN, it's not needed (bmy, 7/20/04)
!  23 Oct 2012 - R. Yantosca - Added ProTeX headers
!  23 Oct 2012 - J. Fisher   - Now use (I,J,L) to reference CO_PRODS 
!                              and CO_LOSSS arrays for GMI strat chem
!EOP
!------------------------------------------------------------------------------
!BOC
!
      !=================================================================
      ! Pick P(CO) or L(CO) depending on IS_PROD
      !=================================================================
      IF ( IS_PROD ) THEN
         RATE = CO_PRODS(I,J,L)   ! P(CO) from CH4 + OH in [v/v/s]
      ELSE
         RATE = CO_LOSSS(I,J,L)   ! L(CO) from CO + OH  in [s^-1]
      ENDIF
    
      END FUNCTION GET_PCO_LCO_STRAT
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_tagged_co
!
! !DESCRIPTION: Subroutine INIT\_TAGGED\_CO allocates memory to module arrays.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE INIT_TAGGED_CO( am_I_Root, Input_Opt, RC )
!
! !USES:
!
      USE CMN_SIZE_MOD
      USE ERROR_MOD,          ONLY : ALLOC_ERR
      USE GIGC_ErrCode_MOd
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU?
      TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?      
! 
! !REVISION HISTORY: 
!  19 Jul 2000 - R. Yantosca - Initial version
!  (1 ) Added ISOP96, MONO96, CH3OH96 for GEOS-3 (bnd, bmy, 6/14/01)
!  (2 ) Removed ISOP96, MONO96, CH3OH96 for GEOS-3, since the new GEOS-3
!        fields make these no longer necessary (bmy, 8/21/09)
!  (3 ) Now allocate BB_REGION, FF_REGION as (IIPAR,JJPAR) (bmy, 9/28/01)
!  (4 ) Removed obsolete code from 9/28/01 (bmy, 10/22/01)
!  (5 ) Now references ALLOC_ERR from "error_mod.f" (bmy, 10/15/02)
!  (6 ) Now remove IJLOOP_CO (bmy, 7/20/04)
!  (7 ) Now public. Now references ITS_A_H2HD_SIM from "tracer_mod.f". 
!        Allocate needed variables if H2/HD simulation (phs, 9/18/07)
!  23 Oct 2012 - R. Yantosca - Added ProTeX headers
!  23 Oct 2012 - J. Fisher   - Dimension CO_PRODS and CO_LOSSS with
!                              (IIPAR,JJPAR,LLPAR) for new GMI strat chem
!  10 Apr 2014 - R. Yantosca - Now call DEFINE* routines to define the
!                              fossil fuel and biofuel regions
!  10 Apr 2014 - R. Yantosca - Now accept am_I_Root, Input_Opt, RC arguments
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER :: AS, I, J, IJLOOP
      
      !=================================================================
      ! INIT_TAGGED_CO begins here!
      !=================================================================

      ! Assume success until otherwise known
      RC = GIGC_SUCCESS

      ! Allocate SUMACETCO -- array for CO from isoprene
      ALLOCATE( EMACET( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EMACET' )         
      EMACET = 0e+0_fp 

      ! Allocate and initialize IJLOOP_CO -- 1-D array index
      ALLOCATE( CO_PRODS( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'CO_PRODS' )    
!      CO_PRODS = 0e+0_fp
      CO_PRODS = 0e0

      ! Allocate and initialize IJLOOP_CO -- 1-D array index
      ALLOCATE( CO_LOSSS( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'CO_LOSSS' )    
!      CO_LOSSS = 0e+0_fp
      CO_PRODS = 0e0

      ! Allocate only what is needed
      IF ( Input_Opt%ITS_A_H2HD_SIM ) RETURN

      ! Allocate BB_REGION -- array for biomass burning regions
      ALLOCATE( BB_REGION( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'BB_REGION' )         
      BB_REGION = 0

      ! Allocate FF_REGION -- array for fossil fuel regions
      ALLOCATE( FF_REGION( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'FF_REGION' )
      FF_REGION = 0
      
      ! Allocate SUMISOPCO -- array for CO from isoprene
      ALLOCATE( SUMISOPCO( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'SUMISOPCO' )         
      SUMISOPCO = 0e+0_fp 

      ! Allocate SUMISOPCO -- array for CO from isoprene
      ALLOCATE( SUMMONOCO( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'SUMMONOCO' )         
      SUMMONOCO = 0e+0_fp 

      ! Allocate SUMISOPCO -- array for CO from isoprene
      ALLOCATE( SUMCH3OHCO( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'SUMCH3OHCO' )         
      SUMCH3OHCO = 0e+0_fp 

      ! Allocate SUMACETCO -- array for CO from isoprene
      ALLOCATE( SUMACETCO( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'SUMACETCO' )         
      SUMACETCO = 0e+0_fp 

      ! Allocate TCOSZ -- array for sum of COS(SZA)
      ALLOCATE( TCOSZ( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'TCOSZ' )         
      TCOSZ = 0d0 

!      ! Define geographic regions for both fossil fuel & biomass burning
!      CALL DEFINE_FF_CO_REGIONS( FF_REGION )         
!      CALL DEFINE_BB_CO_REGIONS( BB_REGION )

      END SUBROUTINE INIT_TAGGED_CO
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cleanup_tagged_co
!
! !DESCRIPTION: Subroutine CLEANUP\_TAGGED\_CO deallocates memory from 
!  previously allocated module arrays.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CLEANUP_TAGGED_CO
! 
! !REVISION HISTORY: 
!  19 Jul 2000 - R. Yantosca - Initial version
!  (1 ) Added ISOP96, MONO96, CH3OH96 for GEOS-3 (bnd, bmy, 6/14/01)
!  (2 ) Removed ISOP96, MONO96, CH3OH96 for GEOS-3, since the new GEOS-3
!        fields make these no longer necessary (bmy, 8/21/09)
!  (3 ) Now remove IJLOOP_CO (bmy, 7/20/04)
!  23 Oct 2012 - R. Yantosca - Added ProTeX headers
!  11 Mar 2015 - R. Yantosca - Now free the OH pointer
!EOP
!------------------------------------------------------------------------------
!BOC
!
      IF ( ALLOCATED( BB_REGION  ) ) DEALLOCATE( BB_REGION  )
      IF ( ALLOCATED( FF_REGION  ) ) DEALLOCATE( FF_REGION  )
      IF ( ALLOCATED( SUMISOPCO  ) ) DEALLOCATE( SUMISOPCO  )
      IF ( ALLOCATED( SUMMONOCO  ) ) DEALLOCATE( SUMMONOCO  )
      IF ( ALLOCATED( SUMCH3OHCO ) ) DEALLOCATE( SUMCH3OHCO )
      IF ( ALLOCATED( SUMACETCO  ) ) DEALLOCATE( SUMACETCO  )
      IF ( ALLOCATED( EMACET     ) ) DEALLOCATE( EMACET     )
      IF ( ALLOCATED( CO_PRODS   ) ) DEALLOCATE( CO_PRODS   )
      IF ( ALLOCATED( TCOSZ      ) ) DEALLOCATE( TCOSZ      )
      IF ( ALLOCATED( CO_LOSSS   ) ) DEALLOCATE( CO_LOSSS   )

      ! Free OH pointer
      OH => NULL()

      END SUBROUTINE CLEANUP_TAGGED_CO
!EOC
      END MODULE TAGGED_CO_MOD
