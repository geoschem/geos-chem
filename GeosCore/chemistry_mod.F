!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: chemistry_mod.F
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
      USE GEOS_TIMERS_MOD  ! For GEOS-Chem timers (optional)

      IMPLICIT NONE
      PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC  :: INIT_CHEMISTRY
      PUBLIC  :: DO_CHEMISTRY
      PUBLIC  :: RECOMPUTE_OD
      PRIVATE :: CHEM_PASSIVE_TRACER
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
!  13 Aug 2015 - E. Lundgren - Tracer units are now kg/kg and converted to
!                              kg within DO_CHEMISTRY
!  03 Nov 2016 - C. Keller   - Added wrapper routine for passive tracers.
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
      USE AEROSOL_MOD,        ONLY : RDAER
      USE AEROSOL_MOD,        ONLY : SOILDUST
      USE C2H6_MOD,           ONLY : CHEMC2H6
      USE CARBON_MOD,         ONLY : CHEMCARBON
      USE CMN_DIAG_MOD        
      USE CMN_SIZE_MOD        
      USE DUST_MOD,           ONLY : CHEMDUST
      USE DUST_MOD,           ONLY : RDUST_ONLINE
      USE ErrCode_Mod
      USE ERROR_MOD
      USE FlexChem_Mod,       ONLY : Do_FlexChem
      USE GLOBAL_CH4_MOD,     ONLY : CHEMCH4
      USE Input_Opt_Mod,      ONLY : OptInput
      USE ISOROPIAII_MOD,     ONLY : DO_ISOROPIAII
      USE MERCURY_MOD,        ONLY : CHEMMERCURY
      USE POPS_MOD,           ONLY : CHEMPOPS
      USE RnPbBe_MOD,         ONLY : CHEMRnPbBe
      USE RPMARES_MOD,        ONLY : DO_RPMARES
      USE SEASALT_MOD,        ONLY : CHEMSEASALT
      USE SULFATE_MOD,        ONLY : CHEMSULFATE
      USE State_Chm_Mod,      ONLY : ChmState
      USE State_Chm_Mod,      ONLY : Ind_
      USE State_Met_Mod,      ONLY : MetState
      USE STRAT_CHEM_MOD,     ONLY : DO_STRAT_CHEM
      USE TAGGED_CO_MOD,      ONLY : CHEM_TAGGED_CO
      USE TAGGED_O3_MOD,      ONLY : CHEM_TAGGED_O3
      USE TIME_MOD,           ONLY : GET_ELAPSED_MIN
      USE TIME_MOD,           ONLY : GET_TS_CHEM
#if defined( USE_TEND )
      USE TENDENCIES_MOD
#endif
      USE UnitConv_Mod

#if defined (UCX )
      USE UCX_MOD,            ONLY : CALC_STRAT_AER ! (SDE 04/20/13)
      USE UCX_MOD,            ONLY : READ_PSC_FILE
      USE UCX_MOD,            ONLY : WRITE_STATE_PSC
#endif
#if defined( TOMAS )
      USE TOMAS_MOD,          ONLY : DO_TOMAS  !(win, 7/14/09)
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
!  25 Mar 2013 - H. Amos     - merged C. Friedman's PAH code into v9-01-03
!  28 Mar 2013 - S.D. Eastham- Updated to use FAST_JX_MOD
!  31 May 2013 - R. Yantosca - Now pass Input_Opt, State_Chm to DO_TOMAS
!  19 May 2014 - C. Keller   - Removed call for acetone ocean sink - now done
!                              in HEMCO.
!  06 Nov 2014 - M. Yannetti - Added PRECISION_MOD
!  08 May 2015 - C. Keller   - Added WRITE_STATE_PSC.
!  18 May 2015 - R. Yantosca - Remove DIAG_STATE_PSC, that is not used anymore
!  15 Jun 2015 - R. Yantosca - Removed calls to DRYFLXRnPbBe, that's obsolete
!  04 Sep 2015 - C. Keller   - Added passive tracer call.
!  17 Mar 2016 - M. Sulprizio- Remove call to OPTDEPTH. The optical depth fields
!                              are now saved into State_Met%OPTD in the routines
!                              that read the met fields from disk.
!  16 May 2016 - M. Sulprizio- Remove call to AEROSOL_RURALBOX. The FlexChem
!                              implementation has rendered the routine obsolete.
!  16 Jun 2016 - C. Miller   - Now use Ind_ function to define species ID's
!  17 Jun 2016 - R. Yantosca - Now define species ID's only on first call
!  17 Jun 2016 - R. Yantosca - Now reset first-time flag at end of routine
!  30 Jun 2016 - R. Yantosca - Remove instances of STT.
!  19 Jul 2016 - R. Yantosca - Now bracket DO_TEND calls with #ifdef USE_TEND
!  10 Aug 2016 - R. Yantosca - Remove temporary tracer-removal code
!  11 Aug 2016 - R. Yantosca - Clean up calls to error subroutines
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
! 
      ! Scalars
      INTEGER            :: N_TROP, N
      INTEGER            :: MONTH
      INTEGER            :: YEAR
      INTEGER            :: WAVELENGTH
      LOGICAL            :: IT_IS_A_C2H6_SIM
      LOGICAL            :: IT_IS_A_CH3I_SIM
      LOGICAL            :: IT_IS_A_CH4_SIM
      LOGICAL            :: IT_IS_A_FULLCHEM_SIM
      LOGICAL            :: IT_IS_A_H2HD_SIM
      LOGICAL            :: IT_IS_A_HCN_SIM
      LOGICAL            :: IT_IS_A_MERCURY_SIM
      LOGICAL            :: IT_IS_A_RnPbBe_SIM
      LOGICAL            :: IT_IS_A_TAGCO_SIM
      LOGICAL            :: IT_IS_A_TAGO3_SIM
      LOGICAL            :: IT_IS_AN_AEROSOL_SIM
      LOGICAL            :: IT_IS_NOT_COPARAM_OR_CH4
      LOGICAL            :: IT_IS_A_POPS_SIM
      LOGICAL            :: LCARB
      LOGICAL            :: LCHEM
      LOGICAL            :: LCRYST
      LOGICAL            :: LDUST
      LOGICAL            :: LSCHEM
      LOGICAL            :: LPRT
      LOGICAL            :: LSSALT
      LOGICAL            :: LSULF
      LOGICAL            :: LSOA
      LOGICAL            :: LNLPBL

      ! SAVEd scalars
      LOGICAL, SAVE      :: FIRST = .TRUE.
      INTEGER, SAVE      :: id_DST1
      INTEGER, SAVE      :: id_NK1
#if defined( USE_TEND )
      INTEGER, SAVE      :: id_O3 
      INTEGER, SAVE      :: id_CO
      REAL(fp)           :: DT_TEND
#endif

      ! Strings
      CHARACTER(LEN=255) :: MSG, LOC

      !=================================================================
      ! DO_CHEMISTRY begins here!
      !=================================================================

      ! Set the location for error messages
      LOC = 'DO_CHEMISTRY (in module GeosCore/chemistry_mod.F)'

      ! Assume success
      RC                       = GC_SUCCESS

      ! Copy fields from INPUT_OPT to local variables for use below
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
      IT_IS_A_TAGO3_SIM        = Input_Opt%ITS_A_TAGO3_SIM
      IT_IS_A_POPS_SIM         = Input_Opt%ITS_A_POPS_SIM
      IT_IS_AN_AEROSOL_SIM     = Input_Opt%ITS_AN_AEROSOL_SIM
      IT_IS_NOT_COPARAM_OR_CH4 = Input_Opt%ITS_NOT_COPARAM_OR_CH4

      ! Save species ID"s on first call
      IF ( FIRST ) THEN
         id_DST1 = Ind_('DST1')
         id_NK1  = Ind_('NK1' )
#if defined( USE_TEND )
         id_O3   = Ind_('O3'  )
         id_CO   = Ind_('CO'  )
#endif
      ENDIF

      !-----------------------------------------------------------------
      !         %%%%%%% GEOS-Chem HP (with ESMF & MPI) %%%%%%%
      !
      ! DO_CHEMISTRY is passed species concentration in 
      ! units of kg/kg and must be converted to kg at the start of the
      ! routine. Units are then converted back to kg/kg at the end. To work
      ! in an ESMF environment, GIGC_Do_Chem (in ESMF/gigc_chemdr.F90)
      ! will need to be updated to convert State_Chm%Species from [v/v]
      ! to [kg/kg dry air]. (ewl, 8/13/15)
      !-----------------------------------------------------------------

      ! Convert species units from [kg/kg] to [kg] for chemistry (ewl, 8/12/15)
      CALL ConvertSpc_KgKgDry_to_Kg( am_I_Root, State_Met, 
     &                               State_Chm, RC         ) 
      IF ( RC /= GC_SUCCESS ) THEN
         MSG = 'Unit conversion error!'
         CALL GC_Error( MSG, RC, LOC )
         RETURN
      ENDIF  

#if defined( USE_TEND )
      !=================================================================
      ! Archive species concentrations for tendencies (ckeller,7/15/2015)
      !=================================================================
      CALL TEND_STAGE1( am_I_Root, Input_Opt, State_Met, 
     &                  State_Chm, 'CHEM', RC )
#endif

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
            CALL READ_PSC_FILE( am_I_Root, State_Chm, RC )
            IF ( RC /= GC_SUCCESS ) RETURN
#endif
#if defined( USE_TIMERS )
            CALL GEOS_Timer_Start( "=> Gas-phase chem", RC )
#endif

            ! Call gas-phase chemistry
            CALL Do_FlexChem
     &         ( am_I_Root, Input_Opt, State_Met, State_Chm, RC )

            ! Check units (ewl, 10/5/15)
            IF ( TRIM( State_Chm%Spc_Units ) /= 'kg' ) THEN
               MSG = 'Incorrect species units after FLEX_CHEMDR!'
               CALL ERROR_STOP( MSG, LOC )
            ENDIF  

#if defined( USE_TIMERS )
            CALL GEOS_Timer_End( "=> Gas-phase chem", RC )
#endif

            ! Call linearized stratospheric scheme
            IF ( LSCHEM ) THEN 

#if defined( USE_TIMERS )
               CALL GEOS_Timer_Start( "=> Strat chem", RC )
#endif

               CALL DO_STRAT_CHEM
     &            ( am_I_Root, Input_Opt, State_Met, State_Chm, RC )

               ! Check units (ewl, 10/5/15)
               IF ( TRIM( State_Chm%Spc_Units ) /= 'kg' ) THEN
                  MSG = 'Incorrect species units after DO_STRAT_CHEM!'
                  CALL ERROR_STOP( MSG, LOC )
               ENDIF  

#if defined( USE_TIMERS )
               CALL GEOS_Timer_End( "=> Strat chem", RC )
#endif

            ENDIF

#if defined( USE_TIMERS )
            CALL GEOS_Timer_Start( "=> All aerosol chem", RC )
#endif

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

               ! Check units (ewl, 10/5/15)
               IF ( TRIM( State_Chm%Spc_Units ) /= 'kg' ) THEN
                  MSG =  'Incorrect species units after CHEMSULFATE!'
                  CALL ERROR_STOP( MSG, LOC )
               ENDIF  

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
            IF ( LDUST .AND. id_DST1 > 0 ) THEN
               CALL CHEMDUST
     &            ( am_I_Root, Input_Opt, State_Met, State_Chm, RC )
            ENDIF
 
#if   defined( TOMAS )
            ! Do TOMAS aerosol microphysics and dry dep
            IF ( id_NK1 > 0 ) THEN 
               CALL DO_TOMAS
     &             ( am_I_Root, Input_Opt, State_Met, State_Chm, RC )
               IF ( RC /= GC_SUCCESS ) RETURN
            ENDIF

            ! Check units (ewl, 10/5/15)
            IF ( TRIM( State_Chm%Spc_Units ) /= 'kg' ) THEN
               MSG = 'Incorrect species units after DO_TOMAS!' 
               CALL ERROR_STOP( MSG, LOC )
            ENDIF  

#endif

#if defined( USE_TIMERS )
            CALL GEOS_Timer_End( "=> All aerosol chem", RC )
#endif

            ! ND43 chemical production
            CALL DIAGOH()

            ! Write STATE_PSC to diagnostics. This is only of relevance in
            ! an ESMF environment, where the restart variables (STATE_PSC)
            ! need to be written to the internal state object on every time
            ! step (ckeller, 5/8/15). 
#if defined( UCX )
            CALL WRITE_STATE_PSC( am_I_Root, State_Chm, RC )
            IF ( RC /= GC_SUCCESS ) RETURN
#endif
 
         !---------------------------------
         ! Offline aerosol simulation
         !---------------------------------
         ELSE IF ( IT_IS_AN_AEROSOL_SIM ) THEN

#if defined( USE_TIMERS )
            CALL GEOS_Timer_Start( "=> All aerosol chem", RC )
#endif

            ! Compute aerosol & dust concentrations [kg/m3]
            ! (NOTE: SOILDUST in "aerosol_mod.f" is computed here)
            CALL AEROSOL_CONC
     &         ( am_I_Root, Input_Opt, State_Met, State_Chm, RC )

            ! Check units (ewl, 10/5/15)
            IF ( TRIM( State_Chm%Spc_Units ) /= 'kg' ) THEN
               MSG = 'Incorrect species units after AEROSOL_CONC!'             
               CALL ERROR_STOP( MSG, LOC )
            ENDIF  

            ! Compute AOD's and surface areas at 999 nm
            MONTH      = 0
            YEAR       = 0
            WAVELENGTH = 0
            CALL RDAER( am_I_Root, Input_Opt, State_Met, State_Chm,
     &                  RC,        MONTH,     YEAR,      WAVELENGTH )

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
               ! Do dust aerosol chemistry
               CALL CHEMDUST
     &            ( am_I_Root, Input_Opt, State_Met, State_Chm, RC )

               ! Compute dust OD's & surface areas
               WAVELENGTH = 0
               CALL RDUST_ONLINE
     &             ( am_I_Root, Input_Opt,  State_Met, State_Chm,
     &               SOILDUST,  WAVELENGTH, RC )
            ENDIF

#if defined( USE_TIMERS )
            CALL GEOS_Timer_End( "=> All aerosol chem", RC )
#endif

         !---------------------------------
         ! Rn-Pb-Be
         !---------------------------------                 
         ELSE IF ( IT_IS_A_RnPbBe_SIM ) THEN
            CALL CHEMRnPbBe 
     &         (  am_I_Root, Input_Opt, State_Met, State_Chm, RC )

         !---------------------------------
         ! Tagged O3
         !---------------------------------
         ELSE IF ( IT_IS_A_TAGO3_SIM ) THEN 
            CALL CHEM_TAGGED_O3
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

         !-----------------------------------------------------
         ! Passive tracers
         ! This performs a simple loss chemistry on all passive 
         ! tracers. Call this routine for all simulation types
         ! since passive tracers can be defined for various 
         ! simulations (as additional tracers to the default
         ! ones). ckeller, 09/04/15
         ! CHEM_PASSIVE_TRACER is defined below (ckeller, 11/3/16)
         !-----------------------------------------------------
         CALL CHEM_PASSIVE_TRACER
     &       ( am_I_Root, Input_Opt, State_Met, State_Chm, RC )
         IF ( RC /= GC_SUCCESS ) RETURN

         !### Debug
         IF ( LPRT .and. am_I_Root ) THEN
            CALL DEBUG_MSG( '### MAIN: a CHEMISTRY' )
         ENDIF

      ENDIF
     
      !-----------------------------------------------------------------
      !         %%%%%%% GEOS-Chem HP (with ESMF & MPI) %%%%%%%
      !
      ! DO_CHEMISTRY is passed species concentration in units of 
      ! kg/kg and is converted to kg at the start of the routine.
      ! Units are then converted back to kg/kg at the end. To work in
      ! an ESMF environment, GIGC_Do_Chem (in ESMF/gigc_chemdr.F90)
      ! will need to be updated to convert State_Chm%Species from [kg]
      ! to [kg/kg dry air] at the end (not [v/v]). (ewl, 8/13/15)
      !-----------------------------------------------------------------
      ! Convert species units back to [kg/kg] after chemistry (ewl, 8/12/15)
      CALL ConvertSpc_Kg_to_KgKgDry( am_I_Root, State_Met,  
     &                               State_Chm, RC         ) 
      IF ( RC /= GC_SUCCESS ) THEN
         MSG = 'Unit conversion error!'
         CALL GC_Error( MSG, RC, LOC )
         RETURN
      ENDIF  

#if defined( USE_TEND )
      !=================================================================
      ! Calculate tendencies and write to diagnostics 
      ! (ckeller,7/15/2015)
      !=================================================================

      ! Chemistry timestep [s]
      DT_TEND = GET_TS_CHEM() * 60.0_fp

      ! Compute tendencies
      CALL TEND_STAGE2( am_I_Root, Input_Opt, State_Met, 
     &                  State_Chm, 'CHEM', DT_TEND, RC )
#endif
    
      END SUBROUTINE DO_CHEMISTRY
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
      USE ErrCode_Mod
      USE ERROR_MOD,          ONLY : DEBUG_MSG
      USE Input_Opt_Mod,      ONLY : OptInput
      USE State_Chm_Mod,      ONLY : ChmState
      USE State_Met_Mod,      ONLY : MetState
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
!  12 Aug 2015 - E. Lundgren  - Input tracer units are now [kg/kg] and 
!                               are converted to [kg] for recomputing OD
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
      RC                   = GC_SUCCESS

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

         ! Then make sure that the simulations use aerosol species
         IF ( IT_IS_A_FULLCHEM_SIM .or. IT_IS_AN_AEROSOL_SIM ) THEN

            ! And then make sure that the aersol species are defined
            IF ( LSULF .or. LCARB .or. LDUST .or. LSSALT ) THEN

               ! Skip this section if all of these are turned off
               CALL AEROSOL_CONC
     &            ( am_I_Root, Input_Opt, State_Met, State_Chm, RC )

               !===============================================================
               ! Call RDAER -- computes aerosol optical depths
               !===============================================================

               ! Calculate the AOD at the wavelength specified in jv_spec_aod
               WAVELENGTH = 1
               CALL RDAER( am_I_Root, Input_Opt, State_Met, State_Chm,
     &                     RC,        MONTH,     YEAR,      WAVELENGTH )

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
     &               ( am_I_Root, Input_Opt,  State_Met, State_Chm,
     &                 SOILDUST,  WAVELENGTH, RC )
#if  !defined( TOMAS )
               ELSE
                  CALL RDUST_OFFLINE
     &               ( am_I_Root, Input_Opt, State_Met,  State_Chm,
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
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: chem_passive_tracer
!
! !DESCRIPTION: Subroutine RUN\_PASSIVE\_TRACER performs loss chemistry 
!  on all passive tracers.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CHEM_PASSIVE_TRACER ( am_I_Root, Input_Opt,  
     &                                 State_Met, State_Chm, RC ) 
!
! !USES:
!
      USE ErrCode_Mod
      USE ERROR_MOD,          ONLY : GC_Error
      USE Input_Opt_Mod,      ONLY : OptInput
      USE State_Chm_Mod,      ONLY : ChmState
      USE State_Met_Mod,      ONLY : MetState
      USE State_Chm_Mod,      ONLY : ind_ 
      USE PASSIVE_TRACER_MOD, ONLY : PASSIVE_TRACER_GETRATE
      USE PASSIVE_TRACER_MOD, ONLY : NPASSIVE, PASSIVE_NAME 
      USE TIME_MOD,           ONLY : GET_TS_CHEM
!
! !INPUT PARAMETERS:
!
      LOGICAL,          INTENT(IN   )  :: am_I_Root  ! root CPU?
      TYPE(OptInput),   INTENT(IN   )  :: Input_Opt  ! Input options object
      TYPE(MetState),   INTENT(IN   )  :: State_Met  ! Meteorology state object
!
! !INPUT/OUTPUT PARAMETERS:
!
      TYPE(ChmState),   INTENT(IN   )  :: State_Chm  ! Chemistry state object 
      INTEGER,          INTENT(INOUT)  :: RC         ! Failure or success
!
! !REMARKS
!
! !REVISION HISTORY: 
!  04 Sep 2015 - C. Keller    - Initial version 
!  03 Nov 2016 - C. Keller    - Moved to chemistry_mod
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER            :: N,  GCID
      REAL(fp)           :: DT, Rate
      CHARACTER(LEN=255) :: MSG
      CHARACTER(LEN=255) :: LOC='CHEM_PASSIVE_TRACER (chemistry_mod.F)'

      DT = GET_TS_CHEM() * 60.0_fp

      ! Assume success
      RC = GC_SUCCESS

      ! Return if there are no passive tracers
      IF ( NPASSIVE <= 0 ) RETURN

      ! Prosa
      IF ( am_I_Root ) THEN
         WRITE(*,*) '### Passive tracer chemistry: '
      ENDIF

      ! Do for every passive tracer
      DO N = 1, NPASSIVE

         ! Sanity check (GCID should never be negative) 
         GCID = ind_( TRIM(PASSIVE_NAME(N)) )
         IF ( GCID <= 0 ) CYCLE

         ! Calculate decay rate (unitless)
         CALL PASSIVE_TRACER_GETRATE( am_I_Root, PASSIVE_NAME(N), 
     &                                DT, RATE, RC )
         IF ( RC /= GC_SUCCESS ) THEN
            MSG = 'Passive tracer rate error: '//TRIM(PASSIVE_NAME(N))
            CALL GC_Error( MSG, RC, LOC )
            RETURN
         ENDIF

         ! Apply loss
         State_Chm%Species(:,:,:,GCID) = 
     &      State_Chm%Species(:,:,:,GCID) * Rate

         ! Prosa
         IF ( am_I_Root ) THEN
            write(*,*) '- Pass. tracer, Tracer ID, loss rate: ', 
     &         TRIM(PASSIVE_NAME(N)), GCID, Rate
         ENDIF

      ENDDO

      ! Return w/ success 
      RC = GC_SUCCESS
 
      END SUBROUTINE CHEM_PASSIVE_TRACER
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_chemistry
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
      USE CMN_SIZE_MOD,       ONLY : LLTROP
      USE ErrCode_Mod
      USE FAST_JX_MOD,        ONLY : INIT_FJX
      USE FlexChem_Mod,       ONLY : Init_FlexChem
      USE Input_Opt_Mod,      ONLY : OptInput
      USE State_Chm_Mod,      ONLY : ChmState
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)     :: am_I_Root   ! Is this the root CPU?
!
! !INPUT/OUTPUT PARAMETERS: 
!
      TYPE(OptInput), INTENT(INOUT)  :: Input_Opt   ! Input Options object
      TYPE(ChmState), INTENT(INOUT)  :: State_Chm   ! Chemistry State object
      INTEGER,        INTENT(INOUT)  :: RC          ! Success or failure?
!
! !REVISION HISTORY: 
!  19 May 2014 - C. Keller   - Initial version (stripped from do_chemistry
!                              and chemdr.F)
!  20 Jun 2014 - R. Yantosca - Now pass Input_Opt to INIT_FJX
!  23 Jun 2016 - R. Yantosca - Remove call to SETTRACE, it's obsolete
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL, SAVE :: FIRST = .TRUE.

      !=================================================================
      ! INIT_CHEMISTRY begins here!
      !=================================================================

      ! Assume success
      RC = GC_SUCCESS

      IF ( FIRST ) THEN

         ! Initialize FlexChem
         CALL Init_FlexChem( am_I_Root, Input_Opt, RC )

         ! Initialize Fast-JX photolysis
         CALL INIT_FJX(      am_I_Root, Input_Opt, RC )
         
         ! Adjust first flag
         FIRST = .FALSE.
      ENDIF

      END SUBROUTINE INIT_CHEMISTRY
!EOC
      END MODULE CHEMISTRY_MOD
