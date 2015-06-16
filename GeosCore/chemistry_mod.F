!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: chemistry_mod
!
! !DESCRIPTION: Module CHEMISTRY\_MOD is used to call the proper chemistry 
!  subroutine for the various GEOS-Chem simulations. 
!\\
!\\
! !INTERFACE:
!
      MODULE CHEMISTRY_MOD
!
! !USES:
!
      USE PRECISION_MOD    ! For GEOS-Chem Precision (fp)

      IMPLICIT NONE
      PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC :: INIT_CHEMISTRY
      PUBLIC :: DO_CHEMISTRY
#if defined( USE_REAL8 )
      PUBLIC :: GCKPP_DRIVER
#endif
      PUBLIC :: RECOMPUTE_OD
!
! !REVISION HISTORY: 
!  (1 ) Bug fix in DO_CHEMISTRY (bnd, bmy, 4/14/03)
!  (2 ) Now references DEBUG_MSG from "error_mod.f" (bmy, 8/7/03)
!  (3 ) Now references "tagged_ox_mod.f"(bmy, 8/18/03)
!  (4 ) Now references "Kr85_mod.f" (jsw, bmy, 8/20/03)
!  (5 ) Bug fix: Now also call OPTDEPTH for GEOS-4 (bmy, 1/27/04)
!  (6 ) Now references "carbon_mod.f" and "dust_mod.f" (rjp, tdf, bmy, 4/5/04)
!  (7 ) Now references "seasalt_mod.f" (rjp, bec, bmy, 4/20/04)
!  (8 ) Now references "logical_mod.f", "tracer_mod.f", "diag20_mod.f", and
!        "diag65_mod.f", and "aerosol_mod." (bmy, 7/20/04)
!  (9 ) Now references "mercury_mod.f" (bmy, 12/7/04)
!  (10) Updated for SO4s, NITs chemistry (bec, bmy, 4/13/05)
!  (11) Now call CHEM_HCN_CH3CN from "hcn_ch3cn_mod.f".  Also remove all
!        references to the obsolete CO-OH param simulation. (xyp, bmy, 6/24/05)
!  (12) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (13) Now call MAKE_RH from "main.f" (bmy, 3/16/06)
!  (14) Updated for SOA from isoprene (dkh, bmy, 6/1/06)
!  (15) Remove support for GEOS-1 and GEOS-STRAT met fields (bmy, 8/4/06)
!  (16) For now, replace use RPMARES instead of ISORROPIA. (bmy, 4/2/08)
!  (17) Added KPP chemistry driver subroutine (phs,ks,dhk, 09/15/09)
!  (18) Added public member function recompute_OD (skim, 02/03/11)
!  17 Dec 2009 - R. Yantosca - Added ProTeX headers
!  28 Jan 2010 - C. Carouge, R. Yantosca - Modified for ISORROPIA II
!  08 Aug 2012 - R. Yantosca - Now align IF statements better
!  10 Aug 2012 - R. Yantosca - Cosmetic changes
!  25 Mar 2013 - M. Payer    - Now pass State_Chm to several routines
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!  19 May 2014 - C. Keller   - Added INIT_CHEMISTRY
!  15 Dec 2014 - M. Yannetti - KPP code is commented out unless compiling KPP
!  08 Jan 2015 - M. Sulprizio- Now restrict KPP to REAL*8 to allow for KPP code
!                              to compile properly
!EOP
!------------------------------------------------------------------------------
!BOC
      CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: do_chemistry
!
! !DESCRIPTION: Subroutine DO\_CHEMISTRY is the driver routine which calls 
!  the appropriate chemistry subroutine for the various GEOS-Chem simulations.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE DO_CHEMISTRY( am_I_Root, Input_Opt, 
     &                         State_Chm, State_Met, RC )
!
! !USES:
!
      USE AEROSOL_MOD,        ONLY : AEROSOL_CONC
      USE AEROSOL_MOD,        ONLY : AEROSOL_RURALBOX
      USE AEROSOL_MOD,        ONLY : RDAER
      USE AEROSOL_MOD,        ONLY : SOILDUST
      USE C2H6_MOD,           ONLY : CHEMC2H6
      USE CARBON_MOD,         ONLY : CHEMCARBON
      USE CMN_DIAG_MOD        
      USE CMN_SIZE_MOD        
      USE COMODE_LOOP_MOD     
      USE DUST_MOD,           ONLY : CHEMDUST
      USE DUST_MOD,           ONLY : RDUST_ONLINE
      USE ERROR_MOD,          ONLY : DEBUG_MSG
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GIGC_State_Chm_Mod, ONLY : ChmState
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE GLOBAL_CH4_MOD,     ONLY : CHEMCH4
      USE ISOROPIAII_MOD,     ONLY : DO_ISOROPIAII
      USE MERCURY_MOD,        ONLY : CHEMMERCURY
      USE OPTDEPTH_MOD,       ONLY : OPTDEPTH
      USE POPS_MOD,           ONLY : CHEMPOPS
      USE RnPbBe_MOD,         ONLY : CHEMRnPbBe
      USE RPMARES_MOD,        ONLY : DO_RPMARES
      USE SEASALT_MOD,        ONLY : CHEMSEASALT
      USE SULFATE_MOD,        ONLY : CHEMSULFATE
      USE STRAT_CHEM_MOD,     ONLY : DO_STRAT_CHEM
      USE TAGGED_CO_MOD,      ONLY : CHEM_TAGGED_CO
      USE TAGGED_OX_MOD,      ONLY : CHEM_TAGGED_OX
      USE TIME_MOD,           ONLY : GET_ELAPSED_MIN
      USE TIME_MOD,           ONLY : GET_TS_CHEM
      USE TRACERID_MOD,       ONLY : IDTISOP
      USE TRACERID_MOD,       ONLY : IDTDST1
#if defined (UCX )
      USE UCX_MOD,            ONLY : CALC_STRAT_AER ! (SDE 04/20/13)
      USE UCX_MOD,            ONLY : READ_PSC_FILE
      USE UCX_MOD,            ONLY : WRITE_STATE_PSC
#endif
#if   defined( TOMAS )
      USE TOMAS_MOD,          ONLY : DO_TOMAS  !(win, 7/14/09)
      USE TRACERID_MOD,       ONLY : IDTNK1    !(win, 7/14/09)
#endif
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)    :: am_I_Root   ! Is this the root CPU?
      TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
!
! !INPUT/OUTPUT PARAMETERS:
!
      TYPE(MetState), INTENT(INOUT) :: State_Met   ! Meteorology State object
      TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT)   :: RC          ! Success or failure
!
! !REMARKS:
!  (1) State_Chm%TRACERS needs to enter this routine w/ units of [kg].
!  (2) State_Chm%SPECIES needs to enter this routine w/ units of [molec/cm3].
!  (3) As of 25 Oct 2012, we have not replaced the met-field arrays from
!       dao_mod.F with the fields in State_Met.  This will be done in a
!       piecemeal fashion so as not to disrupt the existing G-C workflow.
!
! !REVISION HISTORY: 
!  (1 ) Now reference DELP, T from "dao_mod.f" since we need to pass this
!        to OPTDEPTH for GEOS-1 or GEOS-STRAT met fields (bnd, bmy, 4/14/03)
!  (2 ) Now references DEBUG_MSG from "error_mod.f" (bmy, 8/7/03)
!  (3 ) Removed call to CHEMO3, it's obsolete.  Now calls CHEM_TAGGED_OX !
!        from "tagged_ox_mod.f" when NSRCX==6.  Now calls Kr85 chemistry if 
!        NSRCX == 12 (jsw, bmy, 8/20/03)
!  (4 ) Bug fix: added GEOS-4 to the #if block in the call to OPTDEPTH.
!        (bmy, 1/27/04)
!  (5 ) Now calls CHEMCARBON and CHEMDUST to do carbon aerosol & dust 
!        aerosol chemistry (rjp, tdf, bmy, 4/2/04)
!  (6 ) Now calls CHEMSEASALT to do seasalt aerosol chemistry 
!        (rjp, bec, bmy, 4/20/04)
!  (7 ) Now references "logical_mod.f" & "tracer_mod.f".  Now references
!        AEROSOL_CONC, AEROSOL_RURALBOX, and RDAER from "aerosol_mod.f".  
!        Now includes "CMN_DIAG" and "comode.h".  Also call READER, READCHEM, 
!        and INPHOT to initialize the FAST-J arrays so that we can save out !
!        AOD's to the ND21 diagnostic for offline runs. (bmy, 7/20/04)
!  (8 ) Now call routine CHEMMERCURY from "mercury_mod.f" for an offline
!        Hg0/Hg2/HgP simulation. (eck, bmy, 12/7/04)
!  (9 ) Now do not call DO_RPMARES if we are doing an offline aerosol run
!        with crystalline sulfur & aqueous tracers (cas, bmy, 1/7/05)
!  (10) Now use ISOROPIA for aer thermodyn equilibrium if we have seasalt 
!        tracers defined, or RPMARES if not.  Now call CHEMSEASALT before
!        CHEMSULFATE.  Now do aerosol thermodynamic equilibrium before
!        aerosol chemistry for offline aerosol runs.  Now also reference 
!        CLDF from "dao_mod.f" (bec, bmy, 4/20/05)
!  (11) Now modified for GCAP met fields.  Now call CHEM_HCN_CH3CN from 
!        "hcn_ch3cn_mod.f".  Also remove allreferences to the obsolete 
!         CO-OH param simulation. (xyp, bmy, 6/23/05)
!  (12) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (13) Now call MAKE_RH from "main.f" (bmy, 3/16/06)
!  (14) Removed ISOP_PRIOR as a local variable (dkh, bmy, 6/1/06)
!  (15) Remove support for GEOS-1 and GEOS-STRAT met fields (bmy, 8/4/06)
!  (16) Now use DRYFLXH2HD and CHEM_H2_HD for H2/HD sim (lyj, phs, 9/18/07)
!  (17) Bug fix: now hardwired to use RPMARES since ISORROPIA can return very
!        unphysical values at low RH.  Wait for ISORROPIA II. (bmy, 4/2/08)
!  (18) The dry deposition diagnostic (ND44) is done in vdiff_mod if using non-
!        local PBL (lin, ccc, 5/29/09)
!  (19) Now calls CHEMPOPS from "pops_mod.f" for an offline POPs simulation
!       (eck, 9/20/10)
!  17 Dec 2009 - R. Yantosca - Added ProTeX headers
!  25 Jan 2010 - R. Yantosca - Now call DO_TOMAS for TOMAS microphysics
!  28 Jan 2010 - C. Carouge, R. Yantosca - Modified for ISORROPIA II
!  19 Mar 2012 - R. Yantosca - Add C-preprocessor switch to shut off 
!                              ISORROPIA to facilitate debugging
!  30 Jul 2012 - R. Yantosca - Now accept am_I_Root as an argument, and pass
!                              this down to lower-level chem routines for GIGC
!  08 Aug 2012 - R. Yantosca - Now align IF statements better
!  10 Aug 2012 - R. Yantosca - Cosmetic changes
!  18 Oct 2012 - R. Yantosca - Rename GC_MET argument to State_Met
!  18 Oct 2012 - R. Yantosca - Rename CHEM_STATE argument to State_Chem
!  19 Oct 2012 - R. Yantosca - Now reference gigc_state_chm_mod.F90
!  19 Oct 2012 - R. Yantosca - Now reference gigc_state_met_mod.F90
!  25 Oct 2012 - R. Yantosca - Add comments for GIGC #ifdefs
!  25 Oct 2012 - R. Yantosca - Add the RC output argument for the GIGC
!  08 Nov 2012 - R. Yantosca - Now pass Input_Opt argument for the GIGC and
!                              use fields of Input_Opt to replace logicals
!  15 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!  26 Nov 2012 - R. Yantosca - Now pass Input_Opt, State_Chm, RC to routine
!                              DO_STRAT_CHEM (in GeosCore/strat_chem_mod.F90)
!  11 Dec 2012 - R. Yantosca - Remove NI, NJ, NL, NCNST arguments; these are
!                              now obtained either from CMN_SIZE_mod.F or
!                              from the Input_Opt object
!  05 Mar 2013 - R. Yantosca - Now pass am_I_Root, Input_Opt, RC to DRYFLX
!  28 Mar 2013 - S.D. Eastham- Updated to use FAST_JX_MOD
!  31 May 2013 - R. Yantosca - Now pass Input_Opt, State_Chm to DO_TOMAS
!  19 May 2014 - C. Keller   - Removed call for acetone ocean sink - now done
!                              in HEMCO.
!  06 Nov 2014 - M. Yannetti - Added PRECISION_MOD
!  08 May 2015 - C. Keller   - Added WRITE_STATE_PSC.
!  18 May 2015 - R. Yantosca - Remove DIAG_STATE_PSC, that is not used anymore
!  15 Jun 2015 - R. Yantosca - Removed calls to DRYFLXRnPbBe, that's obsolete
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
! 
      INTEGER       :: N_TRACERS
      INTEGER       :: N_TROP
      INTEGER       :: MONTH
      INTEGER       :: YEAR
      INTEGER       :: WAVELENGTH
      LOGICAL       :: IT_IS_A_C2H6_SIM
      LOGICAL       :: IT_IS_A_CH3I_SIM
      LOGICAL       :: IT_IS_A_CH4_SIM
      LOGICAL       :: IT_IS_A_FULLCHEM_SIM
      LOGICAL       :: IT_IS_A_H2HD_SIM
      LOGICAL       :: IT_IS_A_HCN_SIM
      LOGICAL       :: IT_IS_A_MERCURY_SIM
      LOGICAL       :: IT_IS_A_RnPbBe_SIM
      LOGICAL       :: IT_IS_A_TAGCO_SIM
      LOGICAL       :: IT_IS_A_TAGOX_SIM
      LOGICAL       :: IT_IS_AN_AEROSOL_SIM
      LOGICAL       :: IT_IS_NOT_COPARAM_OR_CH4
      LOGICAL       :: IT_IS_A_POPS_SIM
      LOGICAL       :: LCARB
      LOGICAL       :: LCHEM
      LOGICAL       :: LCRYST
      LOGICAL       :: LDUST
      LOGICAL       :: LSCHEM
      LOGICAL       :: LPRT
      LOGICAL       :: LSSALT
      LOGICAL       :: LSULF
      LOGICAL       :: LSOA
      LOGICAL       :: LNLPBL

      ! SAVEd variables
      LOGICAL, SAVE :: FIRST = .TRUE.

      ! Pointers
      ! We need to define local arrays to hold corresponding values 
      ! from the Chemistry State (State_Chm) object. (mpayer, 12/6/12)
      REAL(fp), POINTER :: STT(:,:,:,:)

      !=================================================================
      ! DO_CHEMISTRY begins here!
      !=================================================================

      ! Assume success
      RC                       = GIGC_SUCCESS

      ! Copy fields from INPUT_OPT to local variables for use below
      N_TRACERS                = Input_Opt%N_TRACERS
      LCARB                    = Input_Opt%LCARB                        
      LCHEM                    = Input_Opt%LCHEM
      LCRYST                   = Input_Opt%LCRYST
      LDUST                    = Input_Opt%LDUST
      LSCHEM                   = Input_Opt%LSCHEM
      LPRT                     = Input_Opt%LPRT
      LSSALT                   = Input_Opt%LSSALT
      LSULF                    = Input_Opt%LSULF
      LSOA                     = Input_Opt%LSOA
      LNLPBL                   = Input_Opt%LNLPBL
      IT_IS_A_C2H6_SIM         = Input_Opt%ITS_A_C2H6_SIM
      IT_IS_A_CH3I_SIM         = Input_Opt%ITS_A_CH3I_SIM
      IT_IS_A_CH4_SIM          = Input_Opt%ITS_A_CH4_SIM 
      IT_IS_A_FULLCHEM_SIM     = Input_Opt%ITS_A_FULLCHEM_SIM
      IT_IS_A_H2HD_SIM         = Input_Opt%ITS_A_H2HD_SIM
      IT_IS_A_HCN_SIM          = Input_Opt%ITS_A_HCN_SIM
      IT_IS_A_MERCURY_SIM      = Input_Opt%ITS_A_MERCURY_SIM
      IT_IS_A_RnPbBe_SIM       = Input_Opt%ITS_A_RnPbBe_SIM
      IT_IS_A_TAGCO_SIM        = Input_Opt%ITS_A_TAGCO_SIM
      IT_IS_A_TAGOX_SIM        = Input_Opt%ITS_A_TAGOX_SIM
      IT_IS_A_POPS_SIM         = Input_Opt%ITS_A_POPS_SIM
      IT_IS_AN_AEROSOL_SIM     = Input_Opt%ITS_AN_AEROSOL_SIM
      IT_IS_NOT_COPARAM_OR_CH4 = Input_Opt%ITS_NOT_COPARAM_OR_CH4

      !-----------------------------------------------------------------
      !         %%%%%%% GEOS-Chem HP (with ESMF & MPI) %%%%%%%
      !
      ! NOTE: State_Chm%TRACERS needs to enter this routine with
      ! units of [kg].  If we are calling this routine from the 
      ! traditional "main.F" driver, then the tracers are already in
      ! [kg] and no further unit conversion is necessary.  If on the
      ! other hand we are calling this routine from the ESMF interface,
      ! then we need to convert State_Chm%TRACERS from [v/v] to [kg]
      ! in the calling routine GIGC_Do_Chem (in ESMF/gigc_chemdr.F90).
      ! (bmy, 4/4/13) 
      !-----------------------------------------------------------------

      ! Initialize GEOS-Chem tracer array [kg] from Chemistry State object
      ! (mpayer, 12/6/12)
      STT                      => State_Chm%Tracers

      ! Compute optical depths (except for CH4 simulation)
      IF ( .not. IT_IS_A_CH4_SIM ) THEN
         CALL OPTDEPTH( am_I_Root, Input_Opt, State_Met, RC )
      ENDIF

      !=================================================================
      ! If LCHEM=T then call the chemistry subroutines
      !=================================================================

      IF ( LCHEM ) THEN

         !---------------------------------
         ! NOx-Ox-HC (w/ or w/o aerosols) 
         !---------------------------------
         IF ( IT_IS_A_FULLCHEM_SIM ) THEN 

#if defined( UCX )
            ! We need to get the STATE_PSC from the HEMCO restart file.
            ! NOTE: Do this before calling CHEMDR to avoid reading this
            ! data from w/in an !$OMP PARALLEL region. (bmy, 4/8/15)
            ! Now call READ_PSC_FILE on every time step to accomodate 
            ! replay simulations in GEOS-5 (ckeller, 5/8/15)
            CALL READ_PSC_FILE( am_I_Root, RC )
            IF ( RC /= GIGC_SUCCESS ) RETURN
#endif
            
            ! Call gas-phase chemistry
            CALL CHEMDR
     &         ( am_I_Root, Input_Opt, State_Met, State_Chm, RC )

            ! Call linearized stratospheric scheme
            IF ( LSCHEM ) THEN 
               CALL DO_STRAT_CHEM
     &            ( am_I_Root, Input_Opt, State_Met, State_Chm, RC )
            ENDIF

            ! Do seasalt aerosol chemistry
            IF ( LSSALT ) THEN
               CALL CHEMSEASALT
     &            ( am_I_Root, Input_Opt, State_Met, State_Chm, RC )
            ENDIF

#if defined( UCX )
            ! Recalculate PSC properties
            CALL CALC_STRAT_AER
     &         ( am_I_Root, Input_Opt, State_Met, State_Chm, RC )
#endif

            ! Also do sulfate chemistry
            IF ( LSULF ) THEN

               ! Do sulfate chemistry
               CALL CHEMSULFATE
     &            ( am_I_Root, Input_Opt, State_Met, State_Chm, RC )

               ! Do aerosol thermodynamic equilibrium
               IF ( LSSALT ) THEN

#if   !defined( NO_ISORROPIA )
                  ! ISOROPIA takes Na+, Cl- into account
                  CALL DO_ISOROPIAII
     &               ( am_I_Root, Input_Opt, State_Met, State_Chm, RC )
#endif

               ELSE

                  ! RPMARES does not take Na+, Cl- into account
                  CALL DO_RPMARES
     &               ( am_I_Root, Input_Opt, State_Met, State_Chm, RC )

               ENDIF

            ENDIF

            ! Do carbonaceous aerosol chemistry
            IF ( LCARB ) THEN
               CALL CHEMCARBON
     &            ( am_I_Root, Input_Opt, State_Met, State_Chm, RC )
            ENDIF


            ! Do dust aerosol chemistry
            IF ( LDUST .AND. IDTDST1 > 0 ) THEN
               CALL CHEMDUST
     &            ( am_I_Root, Input_Opt, State_Met, State_Chm, RC )
            ENDIF
 
#if   defined( TOMAS )
            ! Do TOMAS aerosol microphysics and dry dep
            IF ( IDTNK1 > 0 ) THEN 
               CALL DO_TOMAS
     &             ( am_I_Root, Input_Opt, State_Met, State_Chm, RC )
            ENDIF
#endif

            ! ND43 chemical production
            CALL DIAGOH()

            ! Write STATE_PSC to diagnostics. This is only of relevance in
            ! an ESMF environment, where the restart variables (STATE_PSC)
            ! need to be written to the internal state object on every time
            ! step (ckeller, 5/8/15). 
#if defined( UCX )
            CALL WRITE_STATE_PSC( am_I_Root, RC )
            IF ( RC /= GIGC_SUCCESS ) RETURN
#endif
 
         !---------------------------------
         ! Offline aerosol simulation
         !---------------------------------
         ELSE IF ( IT_IS_AN_AEROSOL_SIM ) THEN

            ! Define loop index and other SMVGEAR arrays
            ! N_TROP, the # of trop boxes, is returned
            CALL AEROSOL_RURALBOX( N_TROP,    am_I_Root, 
     &                             Input_Opt, State_Met  )
 
            ! Initialize FAST-JX quantities for computing AOD's
            IF ( FIRST ) THEN

               ! Reset NTLOOP and NTTLOOP after call to READER
               ! with the actual # of boxes w/in the ann mean trop
               NTLOOP  = N_TROP
               NTTLOOP = N_TROP
 
               ! Reset first-time flag
               FIRST = .FALSE.
            ENDIF
 
            ! Compute aerosol & dust concentrations [kg/m3]
            ! (NOTE: SOILDUST in "aerosol_mod.f" is computed here)
            CALL AEROSOL_CONC
     &         ( am_I_Root, Input_Opt, State_Met, State_Chm, RC )

            ! Compute AOD's and surface areas at 999 nm
            MONTH      = 0
            YEAR       = 0
            WAVELENGTH = 0
            CALL RDAER( am_I_Root, Input_Opt, State_Met,  RC,
     &                  MONTH,     YEAR,      WAVELENGTH      )
            !*** AEROSOL THERMODYNAMIC EQUILIBRIUM ***
            IF ( LSULF ) THEN
               IF ( LSSALT ) THEN

#if   !defined( NO_ISORROPIA )
                  ! ISOROPIA takes Na+, Cl- into account
                  CALL DO_ISOROPIAII
     &               ( am_I_Root, Input_Opt, State_Met, State_Chm, RC )
#endif

               ELSE

                  ! RPMARES does not take Na+, Cl- into account
                  ! (skip for crystalline & aqueous offline run)
                  IF ( .not. LCRYST ) THEN
                     CALL DO_RPMARES( am_I_Root, Input_Opt, 
     &                                State_Met, State_Chm, RC )
                  ENDIF
               ENDIF
            ENDIF

            !*** SEASALT AEROSOLS ***
            IF ( LSSALT ) THEN
               CALL CHEMSEASALT
     &            ( am_I_Root, Input_Opt, State_Met, State_Chm, RC )
            ENDIF

            !*** SULFATE AEROSOLS ***
            IF ( LSULF .or. LCRYST ) THEN
 
               ! Do sulfate chemistry
               CALL CHEMSULFATE
     &            ( am_I_Root, Input_Opt, State_Met, State_Chm, RC )

            ENDIF

            !*** CARBON AND 2NDARY ORGANIC AEROSOLS ***
            IF ( LCARB ) THEN
               CALL CHEMCARBON
     &            ( am_I_Root, Input_Opt, State_Met, State_Chm, RC )
            ENDIF

            !*** MINERAL DUST AEROSOLS ***
            IF ( LDUST ) THEN 
               ! Do dust aerosol chemsitry
               CALL CHEMDUST
     &            ( am_I_Root, Input_Opt, State_Met, State_Chm, RC )

               ! Compute dust OD's & surface areas
               WAVELENGTH = 0
               CALL RDUST_ONLINE
     &             ( am_I_Root, Input_Opt,  State_Met, 
     &               SOILDUST,  WAVELENGTH, RC         )
            ENDIF
 
         !---------------------------------
         ! Rn-Pb-Be
         !---------------------------------                 
         ELSE IF ( IT_IS_A_RnPbBe_SIM ) THEN
            CALL CHEMRnPbBe 
     &         (  am_I_Root, Input_Opt, State_Met, State_Chm, RC )

         !---------------------------------
         ! Tagged O3
         !---------------------------------
         ELSE IF ( IT_IS_A_TAGOX_SIM ) THEN 
            CALL CHEM_TAGGED_OX
     &         ( am_I_Root, Input_Opt, State_Met, State_Chm, RC )

            ! Call linearized stratospheric scheme
            IF ( LSCHEM ) THEN 
               CALL DO_STRAT_CHEM
     &            ( am_I_Root, Input_Opt, State_met, State_Chm, RC )
            ENDIF

         !---------------------------------
         ! Tagged CO
         !---------------------------------
         ELSE IF ( IT_IS_A_TAGCO_SIM ) THEN
            CALL CHEM_TAGGED_CO
     &         ( am_I_Root, Input_Opt, State_Met, State_Chm, RC )
 
         !---------------------------------
         ! C2H6
         !---------------------------------
         ELSE IF ( IT_IS_A_C2H6_SIM ) THEN
            CALL CHEMC2H6
     &         ( am_I_Root, Input_Opt, State_Met, State_Chm, RC )
 
         !---------------------------------
         ! CH4
         !---------------------------------
         ELSE IF ( IT_IS_A_CH4_SIM ) THEN
 
            ! Only call after the first 24 hours
            IF ( GET_ELAPSED_MIN() >= GET_TS_CHEM() ) THEN
               CALL CHEMCH4
     &            ( am_I_Root, Input_Opt, State_Met, State_Chm, RC )
            ENDIF

         !---------------------------------
         ! Mercury
         !---------------------------------
         ELSE IF ( IT_IS_A_MERCURY_SIM ) THEN
 
            ! Do Hg chemistry
            CALL CHEMMERCURY
     &         ( am_I_Root, Input_Opt, State_Met, State_Chm, RC )

         !---------------------------------
         ! POPs
         !---------------------------------
         ELSE IF ( IT_IS_A_POPS_SIM ) THEN
 
            ! Do POPS chemistry
            CALL CHEMPOPS
     &         ( am_I_Root, Input_Opt, State_Met, State_Chm, RC )

         ENDIF

         !### Debug
         IF ( LPRT .and. am_I_Root ) THEN
            CALL DEBUG_MSG( '### MAIN: a CHEMISTRY' )
         ENDIF
      ENDIF
      
      !-----------------------------------------------------------------
      !         %%%%%%% GEOS-Chem HP (with ESMF & MPI) %%%%%%%
      !
      ! ALSO NOTE: State_Chm%TRACERS needs to leave this routine with
      ! units of [kg].  If we are calling this routine from the 
      ! traditional "main.F" driver, then the tracers are already in
      ! [kg] and no further unit conversion is necessary.  If on the
      ! other hand we are calling this routine from the ESMF interface,
      ! then we need to convert State_Chm%TRACERS from [kg] to [v/v]
      ! in the calling routine GIGC_Do_Chem (in ESMF/gigc_chemdr.F90).
      ! (bmy, 4/4/13) 
      !-----------------------------------------------------------------

      ! Free pointer
      NULLIFY( STT )
    
      END SUBROUTINE DO_CHEMISTRY
!EOC
#if defined( USE_REAL8 )
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: gckpp_driver
!
! !DESCRIPTION: Subroutine GCKPP\_DRIVER is the driver routine to perform 
!  integration with the full KPP chemistry mechanism.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE GCKPP_DRIVER( KTLOOP, JLOOPLO, R_KPP, NSPEC_GC,
     &                         Input_Opt ) 
!
! !USES:
!
      USE COMODE_MOD,           ONLY : JLOP,   CSPEC
      USE COMODE_MOD,           ONLY : IXSAVE, IYSAVE,    IZSAVE
      USE GCKPP_COMODE_MOD,     ONLY : HSAVE_KPP 
      USE TIME_MOD,             ONLY : GET_TS_CHEM
      USE GIGC_Input_Opt_Mod,   ONLY : OptInput
      USE GCKPP_UTIL,           ONLY : SHUFFLE_KPP2USER
      USE GCKPP_INITIALIZE,     ONLY : INITIALIZE
      USE GCKPP_MODEL
      USE GCKPP_GLOBAL    
      USE GCKPP_RATES,          ONLY : UPDATE_RCONST
      USE GCKPP_MONITOR,        ONLY : SPC_NAMES
      USE GCKPP_FUNCTION
      USE ERROR_MOD,            ONLY : ERROR_STOP
      USE GCKPP_INTEGRATOR,     ONLY : NHNEW, NHEXIT, INTEGRATE
      USE CMN_SIZE_MOD,         ONLY : IIPAR, JJPAR
!
! !INPUT PARAMETERS: 
!
      INTEGER,        INTENT(IN) :: KTLOOP       ! Local loop index
      INTEGER,        INTENT(IN) :: JLOOPLO      ! JLOOPLO + KLOOP = JLOOP
      REAL(fp),       INTENT(IN) :: R_KPP(:,:)   ! Array of reaction rates
      INTEGER,        INTENT(IN) :: NSPEC_GC     ! # of active chemical species
      TYPE(OptInput), INTENT(IN) :: Input_Opt  ! Input Options object
!
! !REMARKS:
!  ############################################################################
!  ###   NOTE: THIS ROUTINE WILL BE REPLACED WHEN FLEXCHEM IS ADDED TO    #####
!  ###   GEOS-CHEM IN THE VERSION FOLLOWING GEOS-CHEM v10-01!             #####
!  ############################################################################
!
!  Variables used to pass the last/first step size b/w call 
!                                                                             .
!  For Rosenbrock:
!  ----------------
!      Nhexit=2, Nhnew = 3
!  OUT    
!      RSTATUS(2)  -> Hexit, last accepted step before exit
!      RSTATUS(3)  -> Hnew, last predicted step (not yet taken)
!      For multiple restarts, use Hnew as Hstart in the subsequent run
!  IN
!      RCNTRL(3)   -> Hstart, starting value for the integration step size
!                                                                             .
!                                                                             .
!  For LSODE:
!  ------------
!  OUT
!      RSTATUS(1)  -> Texit, the time corresponding to the
!                     computed Y upon return
!      RSTATUS(2)  -> Hexit, last predicted step before exit
!      For multiple restarts, use Hexit as Hstart in the following run
!
!  IN      
!      RCNTRL(3)  -> Hstart, starting value for the integration step size
!                                                                             .
!                                                                             .
!  For RADAU5:
!  ------------
!  OUT
!      RSTATUS(1)  -> final time     
!  IN      
!      RCNTRL(3)   -> not used
!                                                                             .
!                                                                             .
! For RUNGE-KUTTA
! ---------------
! OUT
!     same as Rosenbrock
! 
! !REVISION HISTORY: 
!  24 Jan 2008 - Kumaresh    - Based on Daven Henze's GCKPP_DRIVER.
!  16 Sep 2009 - R. Yantosca - Commented, and updated to call various 
!  03 Dec 2009 - C. Carouge  - Use CSPEC instead of CSPEC_FOR_KPP 
!                              to save memory space
!  17 Dec 2009 - R. Yantosca - Added ProTeX headers
!  20 Jan 2010 - C. Carouge  - Now call GCKPP_DRIVER from physproc.f to save 
!                              memory. 
!  20 Jan 2010 - C. Carouge  - Now use the # of active species from GC to 
!                              update CSPEC and not the of variable species 
!                              from KPP.
!  12 Apr 2013 - R. Yantosca - If -DDEVEL is used, when the chemistry can't
!                              converge, we shall increase the tolerances
!                              by a factor of 2 and try again.  This is often
!                              needed to run GEOS-Chem in the GEOS-5 GCM.
!  01 Apr 2015 - R. Yantosca - Skip chemistry in the nested-grid buffer zone
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL(fp)        :: T, TIN, TOUT
      INTEGER       :: ICNTRL(20)
      REAL(kind=dp) :: RCNTRL(20)
      INTEGER       :: ISTATUS(20)
      INTEGER       :: I, J, L, N, KLOOP, IERR
      REAL(kind=dp) :: RSTATE(20)
      LOGICAL, SAVE :: FIRST = .TRUE.
      CHARACTER(LEN=255) :: ERR_MSG

      !=================================================================

      STEPMIN = 0.0e+0_fp
      STEPMAX = 0.0e+0_fp

      ! Suggested tolerance for Rosenbrock solvers
      DO i=1,NVAR
#if defined (UCX )
         RTOL(i) = 2.e-2_fp
#else
         RTOL(i) = 2.e-1_fp
#endif
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%% NOTE: This update was made in the SEAC4RS code.  Leave here in case we
!%%% need to restore this. (bmy, 9/25/13)
!         RTOL(i) = 5.d-2
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         ATOL(i) = 1.e-2_fp
      END DO
      
      ! for LSODE
!      DO i=1,NVAR
!         RTOL(i) = 1.0d-3
!         ATOL(i) = 1.0d6
!      END DO

      
      ! Set parameters to default. See comments in Integrator module for
      ! a list of the defaults.
      ICNTRL(:) = 0
      RCNTRL(:) = 0.e+0_fp

      ! Change some parameters from the default to new values
      ICNTRL(1) = 1             ! Autonomous
      ICNTRL(2) = 0             ! Nonautonomous

      ! Select Integrator
      !    ICNTRL(3)  -> selection of a particular method.
      !                  For Rosenbrock, options are:
      !        = 0 :  default method is Rodas3
      !        = 1 :  method is  Ros2
      !        = 2 :  method is  Ros3 
      !        = 3 :  method is  Ros4 
      !        = 4 :  method is  Rodas3
      !        = 5:   method is  Rodas4
      ICNTRL(3) = 4    
      ICNTRL(7) = 1 ! No adjoint

      IF(FIRST)THEN

         RSTATE(Nhexit) = 0e+0_fp

         FIRST = .FALSE. 

      ENDIF

      ! GET TS_CHEM and convert it to seconds. 
      DT = GET_TS_CHEM() * 60e+0_fp

      ! Set time parameters. 
      T = 0e+0_fp
      TIN = T
      TOUT = T + DT

      !=================================================================
      ! Solve Chemistry
      !=================================================================
 100  format('No. of function calls:', i6, /,
     &       'No. of jacobian calls:', i6, /,
     &       'No. of steps:         ', i6, /,
     &       'No. of accepted steps:', i6, /,
     &       'No. of rejected steps ', i6, /,
     $       '       (except at very beginning)',          /,
     &       'No. of LU decompositions:             ', i6, /,
     &       'No. of forward/backward substitutions:', i6, /,
     &       'No. of singular matrix decompositions:', i6, /,
     &       /,
     &       'Texit, the time corresponding to the      ',        /,
     &       '       computed Y upon return:            ', f11.4, /,
     &       'Hexit, last accepted step before exit:    ', f11.4, /,
     &       'Hnew, last predicted step (not yet taken):', f11.4 )


      DO KLOOP = 1, KTLOOP

         ! 1D grid box index on the global grid
         JLOOP         = JLOOPLO+KLOOP
           
         ! Get 3D coords from SMVGEAR's 1D coords
         I = IXSAVE(JLOOP)
         J = IYSAVE(JLOOP)
         L = IZSAVE(JLOOP)

         ! Skipping buffer zone (lzh, 08/10/2014)
         IF ( Input_Opt%ITS_A_NESTED_GRID ) THEN
           IF ( J <=         Input_Opt%NESTED_J0W ) CYCLE
           IF ( J >  JJPAR - Input_Opt%NESTED_J0E ) CYCLE
           IF ( I <=         Input_Opt%NESTED_I0W ) CYCLE
           IF ( I >  IIPAR - Input_Opt%NESTED_I0E ) CYCLE
         ENDIF

         ! Pass tracer concentrations from CSPEC to KPP working vector V_CSPEC
         DO N =1, NVAR
            V_CSPEC(N) = CSPEC(JLOOP,N)
         END DO

         ! Pass tracer concentrations from V_CSPEC to KPP working vector VAR.
         ! This also initializes the constant rate constants FIX.
         CALL Initialize()


         ! starting value for integration time step
         RCNTRL(3) = HSAVE_KPP(I,J,L)

         ! 1D grid box index in the KTLOOP subset.
         ! R_KPP is only defined on KTLOOP boxes (ccc, 12/3/09)
         JLOOP = KLOOP
         CALL Update_RCONST(R_KPP)
         
         ! 1D grid box index on the global grid
         JLOOP = JLOOPLO+KLOOP


         ! Integrate FWD -- phs --
         call integrate(TIN, TOUT, ICNTRL, RCNTRL, ISTATUS,
     $        RSTATE, IERR)
         
         ! Try another time if it failed
         IF ( IERR < 0 ) THEN
#if defined( KPP_SOLVE_ALWAYS )
            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            !%%% For KPP tolerance fix (mlong, 4/12/13) %%%
            !%%% Re-entry point for non-convergence     %%%
            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 123        CONTINUE
#endif
            write(6,*) ''
            write(6,100) istatus(1:8),rstate(1:3)
            write(6,*) ''
            write(ERR_MSG,'(a, i3)') 'Integrator error code :',IERR
            write(6,*)  'JLOOP, I, J, L ', JLOOP, I, J, L

            ! Reset first time step and start concentrations
            RCNTRL(3)  = 0e+0_fp
            CALL Initialize( )  ! v2.1 

            ! 1D grid box index in the KTLOOP subset.
            ! R_KPP is only defined on KTLOOP boxes (ccc, 12/3/09)
            JLOOP = KLOOP
            CALL Update_RCONST(R_KPP)

            ! 1D grid box index on the global grid
            JLOOP = JLOOPLO+KLOOP

#if defined( KPP_SOLVE_ALWAYS )
            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            !%%% FOR KPP TOLERANCE FIX (mlong, 4/13/12) %%%
            !%%% Increase tolerances by factor of 2     %%%
            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           !RTOL = 2d0 * RTOL    ! Relative
            ATOL = 2e+0_fp * ATOL    ! Absolute
#endif

            call Integrate( TIN,    TOUT,    ICNTRL, 
     &                      RCNTRL, ISTATUS, RSTATE, IERR )
            IF ( IERR < 0 ) THEN 
#if defined( KPP_SOLVE_ALWAYS )
               !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
               !%%% FOR KPP TOLERANCE FIX (mlong, 4/13/12) %%%
               !%%% Try to solve again; don't quit the run %%%
               !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
               GOTO 123
#endif
               print*, 'failed twice !!! '
               CALL ERROR_STOP(ERR_MSG, 'INTEGRATE_KPP')
            ENDIF
                            
         ENDIF

!--------- using the adjoint integrator (limited to Rosenbrock)           
!adj!           ! using adjoint integrator           
!adj!           CALL INTEGRATE_ADJ( 1, VAR, VAR_ADJ, TIN, TOUT,ATOL, 
!adj!     &        RTOL, ICNTRL, RCNTRL, ISTATUS, RSTATE)
!adj!
!adj!
!adj!           ! Update the printout and functionality (dkh, 07/28/09)
!adj!           IF ( ISTATUS(NIERR) < 0 ) THEN 
!adj!              print*, 'KPP Integrator failed.  Trying again' 
!adj!              print*, 'IERR = ', RSTATE(20)
!adj!              print*, 'RSTAT = ', RSTATE
!adj!              print*, 'ISTAT = ', ISTATUS
!adj!              print*, 'RCNTRL = ', RCNTRL
!adj!              print*, 'ICNTRL = ', ICNTRL
!adj!              print*, 'JLOOP, I, J, L ', JLOOP, I, J, L
!adj!              rcntrl(3)  = 0d0
!adj!              CALL Initialize( ) ! v2.1 
!adj!              CALL Update_RCONST()
!adj!              CALL INTEGRATE_ADJ( 1, VAR, VAR_ADJ, TIN, TOUT,
!adj!     &             ATOL, RTOL, ICNTRL, RCNTRL, ISTATUS, RSTATE)
!adj!              IF ( ISTATUS(NIERR) < 0 ) THEN
!adj!                 print*, 'failed twice !!! '
!adj!                 CALL ERROR_STOP('IERR < 0 ', 'INTEGRATE_ADJ')
!adj!              ENDIF
!adj!           ENDIF
!--------- using the adjoint integrator (limited to Rosenbrock)           


         ! Set negative values to SMAL2
         DO N = 1, NVAR
            VAR(N) = MAX(VAR(N),SMAL2)
         ENDDO

         ! save next integration time step
         HSAVE_KPP(I,J,L) = RSTATE(Nhnew)

         ! Pass KPP concentrations from VAR to V_CSPEC
         CALL Shuffle_kpp2user(VAR,V_CSPEC)  

         ! Pass KPP concentrations V_CSPEC to geos-chem CSPEC
         DO N =1, NSPEC_GC
            CSPEC(JLOOP,N) = V_CSPEC(N)
         ENDDO

      ENDDO

      ! Return to calling program
      END SUBROUTINE GCKPP_DRIVER
#endif
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: recompute_od
!
! !DESCRIPTION: Subroutine RECOMPUTE\_OD will update the optical depth values 
!  before accumulating or writing the diagnostics.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE RECOMPUTE_OD( am_I_Root, Input_Opt, 
     &                         State_Met, State_Chm, RC )
!
! !USES:
!
      ! References to F90 modules
      USE AEROSOL_MOD,        ONLY : AEROSOL_CONC
      USE AEROSOL_MOD,        ONLY : RDAER
      USE AEROSOL_MOD,        ONLY : SOILDUST
      USE DUST_MOD,           ONLY : RDUST_ONLINE
      USE DUST_MOD,           ONLY : RDUST_OFFLINE
      USE ERROR_MOD,          ONLY : DEBUG_MSG
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GIGC_State_Chm_Mod, ONLY : ChmState
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE TIME_MOD,           ONLY : GET_MONTH
      USE TIME_MOD,           ONLY : GET_YEAR
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)    :: am_I_Root   ! Is this the root CPU?
      TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
      TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
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
!  03 Fev 2011 - Adapted from chemdr.f by skim
!  30 Jul 2012 - R. Yantosca - Now accept am_I_Root as an argument when
!                              running with the traditional driver main.F
!  13 Nov 2012 - R. Yantosca - Now pass Input_Opt and RC arguments for GIGC
!  15 Nov 2012 - M. Payer    - Now pass all met fields via State_Met
!  25 Mar 2013 - R. Yantosca - Now accept am_I_Root, Input_Opt, State_Chm, RC
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL :: IT_IS_A_FULLCHEM_SIM
      LOGICAL :: IT_IS_AN_AEROSOL_SIM
      LOGICAL :: LCARB, LCHEM,  LCRYST,     LDUST
      LOGICAL :: LPRT,  LSSALT, LSULF,      LSOA
      INTEGER :: MONTH, YEAR,   WAVELENGTH

      !=================================================================
      ! RECOMPUTE_OD begins here!
      !=================================================================

      ! Get month and year
      MONTH                = GET_MONTH()
      YEAR                 = GET_YEAR()

      ! Assume success
      RC                   = GIGC_SUCCESS

      ! Copy fields from INPUT_OPT to local variables for use below
      LCARB                = Input_Opt%LCARB 
      LCHEM                = Input_Opt%LCHEM
      LCRYST               = Input_Opt%LCRYST
      LDUST                = Input_Opt%LDUST
      LPRT                 = Input_Opt%LPRT
      LSSALT               = Input_Opt%LSSALT
      LSULF                = Input_Opt%LSULF
      LSOA                 = Input_Opt%LSOA
      IT_IS_A_FULLCHEM_SIM = Input_Opt%ITS_A_FULLCHEM_SIM
      IT_IS_AN_AEROSOL_SIM = Input_Opt%ITS_AN_AEROSOL_SIM    

      ! First make sure chemistry is turned on
      IF ( LCHEM ) THEN

         ! Then make sure that the simulations use aerosol tracers
         IF ( IT_IS_A_FULLCHEM_SIM .or. IT_IS_AN_AEROSOL_SIM ) THEN

            ! And then make sure that the aersol tracers are defined
            IF ( LSULF .or. LCARB .or. LDUST .or. LSSALT ) THEN

               ! Skip this section if all of these are turned off
               CALL AEROSOL_CONC
     &            ( am_I_Root, Input_Opt, State_Met, State_Chm, RC )

               !===============================================================
               ! Call RDAER -- computes aerosol optical depths
               !===============================================================

               ! Calculate the AOD at the wavelength specified in jv_spec_aod
               WAVELENGTH = 1
              write(6,*) 'before RDAER, AOD'
               CALL RDAER( am_I_Root, Input_Opt, State_Met,  RC,
     &                     MONTH,     YEAR,      WAVELENGTH      )
               write(6,*) 'after RDAER, AOD'
               !### Debug
               IF ( LPRT .and. am_I_Root ) THEN 
                  CALL DEBUG_MSG( '### RECOMPUTE_OD: after RDAER' )
               ENDIF

               !===============================================================
               ! If LDUST is turned on, then we have online dust aerosol in
               ! GEOS-CHEM...so just pass SOILDUST to RDUST_ONLINE in order to
               ! compute aerosol optical depth for FAST-JX, etc.
               !
               ! If LDUST is turned off, then we don't have online dust aerosol
               ! in GEOS-CHEM...so read monthly-mean dust files from disk.
               ! (rjp, tdf, bmy, 4/1/04)
               !===============================================================
               IF ( LDUST ) THEN
                  CALL RDUST_ONLINE
     &               ( am_I_Root, Input_Opt,  State_Met,
     &                 SOILDUST,  WAVELENGTH, RC         )
#if  !defined( TOMAS )
               ELSE
                  CALL RDUST_OFFLINE
     &               ( am_I_Root, Input_Opt, State_Met,
     &                 MONTH,     YEAR,      WAVELENGTH, RC )
#endif
               ENDIF

               !### Debug
               IF ( LPRT .and. am_I_Root ) THEN
                  CALL DEBUG_MSG( '### RECOMPUTE_OD: after RDUST' )
               ENDIF

            ENDIF
         ENDIF

      ENDIF

      END SUBROUTINE RECOMPUTE_OD
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: INIT_CHEMISTRY
!
! !DESCRIPTION: Subroutine INIT\_CHEMISTRY initializes chemistry
! variables.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE INIT_CHEMISTRY( am_I_Root, Input_Opt, State_Chm, RC )
!
! !USES:
!
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE CMN_SIZE_MOD,       ONLY : LLTROP
      USE TRACERID_MOD,       ONLY : SETTRACE
      USE GIGC_State_Chm_Mod, ONLY : ChmState
      USE FAST_JX_MOD,        ONLY : INIT_FJX
      USE COMODE_LOOP_MOD
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)     :: am_I_Root   ! Is this the root CPU?
      TYPE(OptInput), INTENT(IN)     :: Input_Opt   ! Input Options object
!
! !INPUT/OUTPUT PARAMETERS: 
!
      TYPE(ChmState), INTENT(INOUT)  :: State_Chm   ! Chemistry State object
      INTEGER,        INTENT(INOUT)  :: RC          ! Success or failure?
!
! !REVISION HISTORY: 
!  19 May 2014 - C. Keller - Initial version (stripped from do_chemistry
!                            and chemdr.F)
!  20 Jun 2014 - R. Yantosca - Now pass Input_Opt to INIT_FJX
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL, SAVE   :: FIRST = .TRUE.

      !=================================================================
      ! INIT_CHEMISTRY begins here!
      !=================================================================

      ! Assume success
      RC                   = GIGC_SUCCESS

      IF ( FIRST ) THEN

         ! Init chem mechanism and photolysis
         CALL READER  ( FIRST, am_I_Root                           )
         CALL READCHEM(        am_I_Root, Input_Opt,            RC )
         CALL INIT_FJX(        am_I_Root, Input_Opt,            RC )
         CALL SETTRACE(        am_I_Root, Input_Opt, State_Chm, RC )

         ! Reset NCS with NCSURBAN
         NCS     = NCSURBAN

         ! Adjust first flag
         FIRST = .FALSE.
      ENDIF

      END SUBROUTINE INIT_CHEMISTRY
!EOC
      END MODULE CHEMISTRY_MOD
