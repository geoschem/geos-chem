!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: chemistry_mod.F90
!
! !DESCRIPTION: Module CHEMISTRY\_MOD is used to call the proper chemistry 
!  subroutine for the various GEOS-Chem simulations. 
!\\
!\\
! !INTERFACE:
!
MODULE Chemistry_Mod
!
! !USES:
!
  USE Precision_Mod    ! For GEOS-Chem Precision (fp)
  USE Geos_Timers_Mod  ! For GEOS-Chem timers (optional)

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC  :: INIT_CHEMISTRY
  PUBLIC  :: DO_CHEMISTRY
  PUBLIC  :: RECOMPUTE_OD
!
! !PRIVATE MEMBER FUNCTIONS:
!
  PRIVATE :: CHEM_PASSIVE_SPECIES
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
!  17 Nov 2017 - R. Yantosca - Now in F90 format; added Diag_OH_HO2_O1D_O3P
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE TYPES:
!

  INTEGER :: id_DST1, id_NK1   ! Species ID flags

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
  SUBROUTINE DO_CHEMISTRY( am_I_Root, Input_Opt,  State_Met,                 &
                           State_Chm, State_Diag, RC                        )
!
! !USES:
!
    USE AEROSOL_MOD,     ONLY : AEROSOL_CONC
    USE AEROSOL_MOD,     ONLY : RDAER
    USE AEROSOL_MOD,     ONLY : SOILDUST
    USE C2H6_MOD,        ONLY : CHEMC2H6
    USE CARBON_MOD,      ONLY : CHEMCARBON
#if defined( BPCH_DIAG )
    USE CMN_DIAG_MOD  
#endif
    USE CMN_SIZE_MOD
#if defined( NC_DIAG )
    USE Diagnostics_Mod, ONLY : Compute_Column_Mass
    USE Diagnostics_Mod, ONLY : Compute_Budget_Diagnostics
#endif
    USE DUST_MOD,        ONLY : CHEMDUST
    USE DUST_MOD,        ONLY : RDUST_ONLINE
    USE ErrCode_Mod      
    USE ERROR_MOD        
    USE FlexChem_Mod,    ONLY : Do_FlexChem
    USE GLOBAL_CH4_MOD,  ONLY : CHEMCH4
    USE Input_Opt_Mod,   ONLY : OptInput
    USE ISOROPIAII_MOD,  ONLY : DO_ISOROPIAII
    USE MERCURY_MOD,     ONLY : CHEMMERCURY
    USE POPS_MOD,        ONLY : CHEMPOPS
    USE RnPbBe_MOD,      ONLY : CHEMRnPbBe
    USE RPMARES_MOD,     ONLY : DO_RPMARES
    USE SEASALT_MOD,     ONLY : CHEMSEASALT
    USE SULFATE_MOD,     ONLY : CHEMSULFATE
    USE State_Chm_Mod,   ONLY : ChmState
    USE State_Chm_Mod,   ONLY : Ind_
    USE State_Diag_Mod,  ONLY : DgnState
    USE State_Met_Mod,   ONLY : MetState
    USE STRAT_CHEM_MOD,  ONLY : DO_STRAT_CHEM
    USE TAGGED_CO_MOD,   ONLY : CHEM_TAGGED_CO
    USE TAGGED_O3_MOD,   ONLY : CHEM_TAGGED_O3
    USE TIME_MOD,        ONLY : GET_ELAPSED_SEC
    USE TIME_MOD,        ONLY : GET_TS_CHEM
#if defined( USE_TEND )  
    USE TENDENCIES_MOD   
#endif                   
#if defined( TOMAS )     
    USE TOMAS_MOD,       ONLY : DO_TOMAS  !(win, 7/14/09)
#endif                   
    USE UCX_MOD,         ONLY : CALC_STRAT_AER
    USE UnitConv_Mod,    ONLY : Convert_Spc_Units
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
    TYPE(DgnState), INTENT(INOUT) :: State_Diag  ! Diagnostics State object
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
!  09 Mar 2017 - C. Keller   - Bug fix: call TEND_STAGE1 before unit conversion
!  28 Sep 2017 - E. Lundgren - Simplify unit conversions with wrapper routine
!  03 Oct 2017 - E. Lundgren - Pass State_Diag as argument
!  09 Nov 2017 - R. Yantosca - Reorder arguments for consistency: Input_Opt,
!                              State_Met, State_Chm, State_Diag
!  20 Nov 2017 - R. Yantosca - Move ND43 diagnostics to flexchem_mod.F90
!  03 Jan 2018 - M. Sulprizio- Replace UCX CPP switch with Input_Opt%LUCX
!  06 Feb 2018 - E. Lundgren - Change GET_ELAPSED_MIN to GET_ELAPSED_SEC to
!                              match new timestep unit of seconds
!  28 Aug 2018 - E. Lundgren - Implement budget diagnostics
!  26 Oct 2018 - M. Sulprizio- Remove calls to read and write STATE_PSC
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
    LOGICAL            :: LDUST
    LOGICAL            :: LSCHEM
    LOGICAL            :: LPRT
    LOGICAL            :: LSSALT
    LOGICAL            :: LSULF
    LOGICAL            :: LSOA
    LOGICAL            :: LNLPBL
    LOGICAL            :: LUCX
#if defined( USE_TEND ) || defined( NC_DIAG )
    REAL(fp)           :: DT_Chem
#endif

    ! SAVEd scalars
    LOGICAL, SAVE      :: FIRST = .TRUE.

    ! Strings
    CHARACTER(LEN=63)  :: OrigUnit
    CHARACTER(LEN=255) :: ErrMsg,  ThisLoc

    !=======================================================================
    ! DO_CHEMISTRY begins here!
    !=======================================================================

    ! Initialize
    RC      = GC_SUCCESS
    ErrMsg  = ''
    ThisLoc = ' -> at Do_Chemistry  (in module GeosCore/chemistry_mod.F90)'

    ! Copy fields from INPUT_OPT to local variables for use below
    LCARB                    = Input_Opt%LCARB                        
    LCHEM                    = Input_Opt%LCHEM
    LDUST                    = Input_Opt%LDUST
    LSCHEM                   = Input_Opt%LSCHEM
    LPRT                     = Input_Opt%LPRT
    LSSALT                   = Input_Opt%LSSALT
    LSULF                    = Input_Opt%LSULF
    LSOA                     = Input_Opt%LSOA
    LNLPBL                   = Input_Opt%LNLPBL
    LUCX                     = Input_Opt%LUCX
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
    
    ! Save species ID"s on first call
    IF ( FIRST ) THEN
       id_DST1 = Ind_('DST1')
       id_NK1  = Ind_('NK1' )
    ENDIF

#if defined( NC_DIAG )
    !----------------------------------------------------------
    ! Chemistry budget diagnostics - Part 1 of 2
    !----------------------------------------------------------
    IF ( State_Diag%Archive_BudgetChemistry ) THEN
       ! Get initial column masses
       CALL Compute_Column_Mass( am_I_Root,                              & 
                                 Input_Opt, State_Met, State_Chm,        &
                                 State_Chm%Map_Advect,                   &
                                 State_Diag%Archive_BudgetChemistryFull, &
                                 State_Diag%Archive_BudgetChemistryTrop, &
                                 State_Diag%Archive_BudgetChemistryPBL,  &
                                 State_Diag%BudgetMass1,                 &
                                 RC ) 
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Chemistry budget diagnostics error 1'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDIF
#endif

#if defined( USE_TEND )
    !=======================================================================
    ! Archive species concentrations for tendencies (ckeller,7/15/2015)
    !=======================================================================
    CALL Tend_Stage1( am_I_Root, Input_Opt, State_Met,                       &
                      State_Chm, 'CHEM', RC                                 )
#endif

    !=======================================================================
    ! Convert species units to [kg] for chemistry (ewl, 8/12/15)
    !=======================================================================
    CALL Convert_Spc_Units( am_I_Root,        Input_Opt, State_Met,          &
                            State_Chm,        'kg',      RC,                 &
                            OrigUnit=OrigUnit                               )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Unit conversion error (kg/kg dry -> kg)'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !=======================================================================
    ! If LCHEM=T then call the chemistry subroutines
    !=======================================================================
    IF ( LCHEM ) THEN

       !====================================================================
       ! Full-chemistry simulations:
       !
       ! (1) Benchmark; (2) Standard; (3) SimpleSOA; (4) complexSOA, 
       ! (5) complexSOA-SVPOA; (6) aciduptake; (7) marinePOA
       !====================================================================
       IF ( IT_IS_A_FULLCHEM_SIM ) THEN 
             
#if defined( USE_TIMERS )
          CALL GEOS_Timer_Start( "=> Gas-phase chem", RC )
#endif

          !----------------------------------------
          ! Dry-run sulfate chem to get cloud pH
          !----------------------------------------
          IF ( LSULF ) THEN

             ! Dry run only
             CALL ChemSulfate( am_I_Root, Input_Opt,  State_Met,            &
                               State_Chm, State_Diag, .FALSE.,  RC         )

             ! Trap potential errors
             IF ( RC /= GC_SUCCESS ) THEN
                ErrMsg = 'Error encountered in "ChemSulfate"!'
                CALL GC_Error( ErrMsg, RC, ThisLoc )
                RETURN
             ENDIF
          ENDIF

          !---------------------------
          ! Call gas-phase chemistry
          !---------------------------
          CALL Do_FlexChem( am_I_Root, Input_Opt,  State_Met,               &
                            State_Chm, State_Diag, RC                      )

          ! Check units (ewl, 10/5/15)
          IF ( TRIM( State_Chm%Spc_Units ) /= 'kg' ) THEN
             ErrMsg = 'Incorrect species units after FLEX_CHEMDR!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF

          ! Trap potential errors
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = 'Error encountered in "Do_FlexChem"!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF

#if defined( USE_TIMERS )
          CALL GEOS_Timer_End( "=> Gas-phase chem", RC )
#endif

          !----------------------------------------
          ! Call linearized stratospheric scheme
          !----------------------------------------
          IF ( LSCHEM ) THEN 

#if defined( USE_TIMERS )
             CALL GEOS_Timer_Start( "=> Strat chem", RC )
#endif

             ! Do linearized chemistry for the stratosphere (tropchem)
             ! or the mesosphere (UCX)
             CALL Do_Strat_Chem( am_I_Root, Input_Opt, State_Met,            &
                                 State_Chm, RC                              )

             ! Check units (ewl, 10/5/15)
             IF ( TRIM( State_Chm%Spc_Units ) /= 'kg' ) THEN
                ErrMsg = 'Incorrect species units after DO_STRAT_CHEM!'
                CALL GC_Error( ErrMsg, RC, ThisLoc )
             ENDIF

             ! Trap potential errors
             IF ( RC /= GC_SUCCESS ) THEN
                ErrMsg = 'Error encountered in ""!'
                CALL GC_Error( ErrMsg, RC, ThisLoc )
                RETURN
             ENDIF

#if defined( USE_TIMERS )
             CALL GEOS_Timer_End( "=> Strat chem", RC )
#endif

          ENDIF

#if defined( USE_TIMERS )
          CALL GEOS_Timer_Start( "=> All aerosol chem", RC )
#endif

          !--------------------------------
          ! Do seasalt aerosol chemistry
          !--------------------------------
          IF ( LSSALT ) THEN
             CALL ChemSeaSalt( am_I_Root, Input_Opt,  State_Met,             &
                               State_Chm, State_Diag, RC                    )

             ! Trap potential errors
             IF ( RC /= GC_SUCCESS ) THEN
                ErrMsg = 'Error encountered in "ChemSeaSalt"!'
                CALL GC_Error( ErrMsg, RC, ThisLoc )
                RETURN
             ENDIF
          ENDIF

          !-------------------------------
          ! Recalculate PSC properties
          !-------------------------------
          IF ( LUCX ) THEN

#if defined( USE_TIMERS )
             CALL GEOS_Timer_End  ( "=> All aerosol chem", RC )
             CALL GEOS_Timer_Start( "=> Strat chem",       RC )
#endif
             
             ! Recalculate PSC
             CALL Calc_Strat_Aer( am_I_Root, Input_Opt, State_Met,           &
                                  State_Chm, RC )

             ! Trap potential errors
             IF ( RC /= GC_SUCCESS ) THEN
                ErrMsg = 'Error encountered in "Calc_Strat_Aer"!'
                CALL GC_Error( ErrMsg, RC, ThisLoc )
                RETURN
             ENDIF

#if defined( USE_TIMERS )
             CALL GEOS_Timer_End  ( "=> Strat chem",       RC )
             CALL GEOS_Timer_Start( "=> All aerosol chem", RC )
#endif

          ENDIF

          !--------------------------------
          ! Also do sulfate chemistry
          !--------------------------------
          IF ( LSULF ) THEN

             ! Do sulfate chemistry
             CALL ChemSulfate( am_I_Root, Input_Opt,  State_Met,             &
                               State_Chm, State_Diag, .TRUE.,    RC         )

             ! Check units (ewl, 10/5/15)
             IF ( TRIM( State_Chm%Spc_Units ) /= 'kg' ) THEN
                ErrMsg =  'Incorrect species units after CHEMSULFATE!'
                CALL GC_Error( ErrMsg, RC, ThisLoc )
             ENDIF

             ! Trap potential errors
             IF ( RC /= GC_SUCCESS ) THEN
                ErrMsg = 'Error encountered in "ChemSulfate"!'
                CALL GC_Error( ErrMsg, RC, ThisLoc )
                RETURN
             ENDIF
       
             !-----------------------------------------
             ! Do aerosol thermodynamic equilibrium
             !-----------------------------------------
             IF ( LSSALT ) THEN

#if   !defined( NO_ISORROPIA )
                ! ISOROPIA takes Na+, Cl- into account
                CALL Do_IsoropiaII( am_I_Root, Input_Opt,  State_Met,        &
                                    State_Chm, State_Diag, RC               )

                ! Trap potential errors
                IF ( RC /= GC_SUCCESS ) THEN
                   ErrMsg = 'Error encountered in "Do_ISOROPIAII"!'
                   CALL GC_Error( ErrMsg, RC, ThisLoc )
                   RETURN
                ENDIF
#endif

             ELSE

                ! RPMARES does not take Na+, Cl- into account
                CALL Do_RPMARES( am_I_Root, Input_Opt, State_Met,            &
                                 State_Chm, RC                              )

             ENDIF

          ENDIF

          !-----------------------------------
          ! Do carbonaceous aerosol chemistry
          !-----------------------------------
          IF ( LCARB ) THEN
             CALL ChemCarbon( am_I_Root, Input_Opt,  State_Met,              &
                              State_Chm, State_Diag, RC                     )

             ! Trap potential errors
             IF ( RC /= GC_SUCCESS ) THEN
                ErrMsg = 'Error encountered in "ChemCarbon"!'
                CALL GC_Error( ErrMsg, RC, ThisLoc )
                RETURN
             ENDIF
          ENDIF

          !------------------------------------
          ! Do dust aerosol chemistry/removal
          !------------------------------------
          IF ( LDUST .AND. id_DST1 > 0 ) THEN
             CALL ChemDust( am_I_Root, Input_Opt,  State_Met,                &
                            State_Chm, State_Diag, RC                       )

             ! Trap potential errors
             IF ( RC /= GC_SUCCESS ) THEN
                ErrMsg = 'Error encountered in "ChemDust"!'
                CALL GC_Error( ErrMsg, RC, ThisLoc )
                RETURN
             ENDIF
          ENDIF
 
#if   defined( TOMAS )
          !--------------------------------------------
          ! Do TOMAS aerosol microphysics and dry dep
          !--------------------------------------------
          IF ( id_NK1 > 0 ) THEN 
             CALL Do_TOMAS( am_I_Root, Input_Opt,  State_Met,               &
                            State_Chm, State_Diag, RC                       )

             ! Check units (ewl, 10/5/15)
             IF ( TRIM( State_Chm%Spc_Units ) /= 'kg' ) THEN
                ErrMsg = 'Incorrect species units after DO_TOMAS!' 
                CALL GC_Error( ErrMsg, RC, ThisLoc )
             ENDIF

             ! Trap potential errors
             IF ( RC /= GC_SUCCESS ) THEN
                ErrMsg = 'Error encountered in "Do_TOMAS"!'
                CALL GC_Error( ErrMsg, RC, ThisLoc )
                RETURN
             ENDIF
          ENDIF
#endif

#if defined( USE_TIMERS )
          CALL GEOS_Timer_End( "=> All aerosol chem", RC )
#endif
          
       !====================================================================
       ! Aerosol-only simulation
       !====================================================================
       ELSE IF ( IT_IS_AN_AEROSOL_SIM ) THEN

#if defined( USE_TIMERS )
          CALL GEOS_Timer_Start( "=> All aerosol chem", RC )
#endif

          !-------------------------------------------------------
          ! Compute aerosol & dust concentrations [kg/m3]
          ! (NOTE: SOILDUST in "aerosol_mod.f" is computed here)
          !-------------------------------------------------------
          CALL Aerosol_Conc( am_I_Root, Input_Opt,  State_Met,               &
                             State_Chm, State_Diag, RC                      )

          ! Check units (ewl, 10/5/15)
          IF ( TRIM( State_Chm%Spc_Units ) /= 'kg' ) THEN
             ErrMsg = 'Incorrect species units after AEROSOL_CONC!'             
             CALL GC_Error( ErrMsg, RC, ThisLoc )
          ENDIF

          ! Trap potential errors
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = 'Error encountered in "Aerosol_Conc"!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF

          !-------------------------------------------
          ! Compute AOD's and surface areas at 999 nm
          !-------------------------------------------
          MONTH      = 0
          YEAR       = 0
          WAVELENGTH = 0
          CALL RdAer( am_I_Root,  Input_Opt, State_Met, State_Chm,           &
                      State_Diag, RC,        MONTH,     YEAR,                &
                      WAVELENGTH                                            )

          ! Trap potential errors
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = 'Error encountered in "RdAer"!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF

          !--------------------------------------------
          ! Aerosol Thermodynamic Equilibrium
          !--------------------------------------------
          IF ( LSULF ) THEN
             IF ( LSSALT ) THEN

#if   !defined( NO_ISORROPIA )
                ! ISOROPIA takes Na+, Cl- into account
                CALL Do_IsoropiaII( am_I_Root, Input_Opt,  State_Met,        &
                                    State_Chm, State_Diag, RC               )
#endif

                ! Trap potential errors
                IF ( RC /= GC_SUCCESS ) THEN
                   ErrMsg = 'Error encountered in "Do_IsoropiaII"!'
                   CALL GC_Error( ErrMsg, RC, ThisLoc )
                   RETURN
                ENDIF

             ELSE

                ! RPMARES does not take Na+, Cl- into account
                ! (skip for crystalline & aqueous offline run)
                CALL Do_RPMARES( am_I_Root, Input_Opt,                    &
                                 State_Met, State_Chm, RC )

                ! Trap potential errors
                IF ( RC /= GC_SUCCESS ) THEN
                   ErrMsg = 'Error encountered in "Do_RPMARES"!'
                   CALL GC_Error( ErrMsg, RC, ThisLoc )
                   RETURN
                ENDIF
             ENDIF
          ENDIF

          !-----------------------------
          ! Seasalt Aerosols
          !-----------------------------
          IF ( LSSALT ) THEN
             CALL ChemSeaSalt( am_I_Root, Input_Opt,  State_Met,             &
                               State_Chm, State_Diag, RC                    )

             ! Trap potential errors
             IF ( RC /= GC_SUCCESS ) THEN
                ErrMsg = 'Error encountered in "ChemSeaSalt"!'
                CALL GC_Error( ErrMsg, RC, ThisLoc )
                RETURN
             ENDIF
          ENDIF

          !-------------------
          ! Sulfate aerosols
          !-------------------
          IF ( LSULF ) THEN
 
             ! Do sulfate chemistry
             CALL ChemSulfate( am_I_Root, Input_Opt,  State_Met,             &
                               State_Chm, State_Diag, .TRUE.,    RC         )

             ! Trap potential errors
             IF ( RC /= GC_SUCCESS ) THEN
                ErrMsg = 'Error encountered in "ChemSulfate"!'
                CALL GC_Error( ErrMsg, RC, ThisLoc )
                RETURN
             ENDIF
          ENDIF
            
          !-----------------------------------------
          ! Carbon and Secondary Organic Aerosols
          !-----------------------------------------
          IF ( LCARB ) THEN
             CALL ChemCarbon( am_I_Root, Input_Opt,  State_Met,              &
                              State_Chm, State_Diag, RC                     )

             ! Trap potential errors
             IF ( RC /= GC_SUCCESS ) THEN
                ErrMsg = 'Error encountered in ""!'
                CALL GC_Error( ErrMsg, RC, ThisLoc )
                RETURN
             ENDIF
          ENDIF

          !------------------------
          ! Mineral Dust Aerosols
          !------------------------
          IF ( LDUST ) THEN 

             ! Do dust aerosol chemistry
             CALL ChemDust( am_I_Root, Input_Opt,  State_Met,                &
                            State_Chm, State_Diag, RC                       )

             ! Trap potential errors
             IF ( RC /= GC_SUCCESS ) THEN
                ErrMsg = 'Error encountered in "ChemDust"!'
                CALL GC_Error( ErrMsg, RC, ThisLoc )
                RETURN
             ENDIF

             ! Compute dust OD's & surface areas
             WAVELENGTH = 0
             CALL Rdust_Online( am_I_Root,  Input_Opt,  State_Met,           &
                                State_Chm,  State_Diag, SOILDUST,            &
                                WAVELENGTH, RC                              )

             ! Trap potential errors
             IF ( RC /= GC_SUCCESS ) THEN
                ErrMsg = 'Error encountered in "Rdust_Online"!'
                CALL GC_Error( ErrMsg, RC, ThisLoc )
                RETURN
             ENDIF
          ENDIF

#if defined( USE_TIMERS )
          CALL GEOS_Timer_End( "=> All aerosol chem", RC )
#endif

       !====================================================================
       ! Rn-Pb-Be
       !====================================================================
       ELSE IF ( IT_IS_A_RnPbBe_SIM ) THEN

#if defined( USE_TIMERS )
          CALL GEOS_Timer_Start( "=> Gas-phase chem", RC )
#endif

          ! Do Rn-Pb-Be chemistry
          CALL ChemRnPbBe( am_I_Root, Input_Opt,  State_Met,                 &
                           State_Chm, State_Diag, RC                        )

          ! Trap potential errors
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = 'Error encountered in "ChemRnPbBe"!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF

#if defined( USE_TIMERS )
          CALL GEOS_Timer_End( "=> Gas-phase chem", RC )
#endif

       !====================================================================
       ! Tagged O3
       !====================================================================
       ELSE IF ( IT_IS_A_TAGO3_SIM ) THEN 

#if defined( USE_TIMERS )
          CALL GEOS_Timer_Start( "=> Gas-phase chem", RC )
#endif

          !-----------------------------------------------
          ! Do Tagged O3 chemistry
          !-----------------------------------------------
          CALL Chem_Tagged_O3( am_I_Root, Input_Opt,  State_Met,             &
                               State_Chm, State_Diag, RC                    )

          ! Trap potential errors
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = 'Error encountered in "Chem_Tagged_O3"!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF
          
#if defined( USE_TIMERS )
          CALL GEOS_Timer_End( "=> Gas-phase chem", RC )
#endif

          !-----------------------------------------------
          ! Call linearized stratospheric scheme (LINOZ)
          !-----------------------------------------------
          IF ( LSCHEM ) THEN 

#if defined( USE_TIMERS )
             CALL GEOS_Timer_Start( "=> Strat chem", RC )
#endif

             ! Do LINOZ for Ozone
             CALL Do_Strat_Chem( am_I_Root, Input_Opt, State_Met,            &
                                 State_Chm, RC                              )

             ! Trap potential errors
             IF ( RC /= GC_SUCCESS ) THEN
                ErrMsg = 'Error encountered in "Do_Strat_Chem"!'
                CALL GC_Error( ErrMsg, RC, ThisLoc )
                RETURN
             ENDIF

#if defined( USE_TIMERS )
             CALL GEOS_Timer_End( "=> Strat chem", RC )
#endif

          ENDIF

       !====================================================================
       ! Tagged CO
       !====================================================================
       ELSE IF ( IT_IS_A_TAGCO_SIM ) THEN

#if defined( USE_TIMERS )
          CALL GEOS_Timer_Start( "=> Gas-phase chem", RC )
#endif

          ! Do tagged CO chemistry
          CALL Chem_Tagged_CO( am_I_Root, Input_Opt,  State_Met,             &
                               State_Chm, State_Diag, RC                    )

          ! Trap potential errors
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = 'Error encountered in "Chem_Tagged_CO"!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF

#if defined( USE_TIMERS )
          CALL GEOS_Timer_End( "=> Gas-phase chem", RC )
#endif

       !====================================================================
       ! C2H6
       !====================================================================
       ELSE IF ( IT_IS_A_C2H6_SIM ) THEN
          CALL ChemC2H6( am_I_Root, Input_Opt, State_Met, State_Chm, RC )
 
          ! Trap potential errors
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = 'Error encountered in "ChemC2H6"!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF

       !====================================================================
       ! CH4
       !====================================================================
       ELSE IF ( IT_IS_A_CH4_SIM ) THEN

#if defined( USE_TIMERS )
          CALL GEOS_Timer_Start( "=> Gas-phase chem", RC )
#endif 

          CALL ChemCh4( am_I_Root, Input_Opt,  State_Met,                 &
                        State_Chm, State_Diag, RC                        )

          ! Trap potential errors
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = 'Error encountered in "ChemCh4"!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF

#if defined( USE_TIMERS )
          CALL GEOS_Timer_End( "=> Gas-phase chem", RC )
#endif

       !====================================================================
       ! Mercury
       !====================================================================
       ELSE IF ( IT_IS_A_MERCURY_SIM ) THEN
 
#if defined( USE_TIMERS )
          CALL GEOS_Timer_Start( "=> Gas-phase chem", RC )
#endif

          ! Do Hg chemistry
          CALL ChemMercury( am_I_Root, Input_Opt,  State_Met,                &
                            State_Chm, State_Diag, RC                       )

          ! Trap potential errors
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = 'Error encountered in "ChemMercury"!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF

#if defined( USE_TIMERS )
          CALL GEOS_Timer_End( "=> Gas-phase chem", RC )
#endif

       !====================================================================
       ! POPs
       !====================================================================
       ELSE IF ( IT_IS_A_POPS_SIM ) THEN
 
#if defined( USE_TIMERS )
          CALL GEOS_Timer_Start( "=> Gas-phase chem", RC )
#endif

          ! Do POPS chemistry
          CALL ChemPOPs( am_I_Root, Input_Opt,  State_Met,                   &
                         State_Chm, State_Diag, RC                          )

          ! Trap potential errors
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = 'Error encountered in "ChemPOPs"!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF

#if defined( USE_TIMERS )
          CALL GEOS_Timer_End( "=> Gas-phase chem", RC )
#endif
       ENDIF

       !====================================================================
       ! PASSIVE SPECIES
       !
       ! This performs a simple loss chemistry on passive species.  Call 
       ! this routine for all simulation types since passive species can 
       ! be defined for various simulations (as additional species to the 
       ! default! ones). ckeller, 09/04/15
       !
       ! NOTE: To speed up execution, only call Chem_Passive_Species
       ! if there is at least one passive species with a finite 
       ! atmospheric lifetime.  There is no reason to apply a loss rate
       ! of unity to those passive species whose lifetime is infinity.  
       ! This will speed up GEOS-Chem simulations. (bmy, 12/13/17)
       !====================================================================
       IF ( Input_Opt%NPassive_Decay > 0 ) THEN

          ! Apply loss rate to passive species with finite lifetimes
          CALL Chem_Passive_Species( am_I_Root, Input_Opt,                   & 
                                     State_Met, State_Chm, RC               )

          ! Trap potential errors
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = 'Error encountered in "Chem_Passive_Species"!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF

          !### Debug
          IF ( LPRT .and. am_I_Root ) THEN
             CALL Debug_Msg( '### MAIN: a CHEMISTRY' )
          ENDIF

       ENDIF

    ENDIF
     
    !=======================================================================
    ! Convert species units back to original unit (ewl, 8/12/15)
    !=======================================================================
    CALL Convert_Spc_Units( am_I_Root, Input_Opt, State_Met,                 &
                            State_Chm, OrigUnit,  RC                        )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Unit conversion error'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

#if defined( USE_TEND ) || defined( NC_DIAG )
    ! Chemistry timestep [s]
    DT_Chem = Get_Ts_Chem()
#endif

#if defined( USE_TEND )
    !=======================================================================
    ! Calculate tendencies and write to diagnostics (ckeller,7/15/2015)
    !=======================================================================

    ! Compute tendencies
    CALL Tend_Stage2( am_I_Root, Input_Opt, State_Met,                       &
                      State_Chm, 'CHEM',    DT_Chem,   RC                   ) 

    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered in ""!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
#endif

#if defined( NC_DIAG )
    !----------------------------------------------------------
    ! Chemistry budget diagnostics - Part 2 of 2
    !----------------------------------------------------------
    IF ( State_Diag%Archive_BudgetChemistry ) THEN
       ! Get final column masses and compute diagnostics
       CALL Compute_Column_Mass( am_I_Root,                              &
                                 Input_Opt, State_Met, State_Chm,        &
                                 State_Chm%Map_Advect,                   &
                                 State_Diag%Archive_BudgetChemistryFull, &
                                 State_Diag%Archive_BudgetChemistryTrop, &
                                 State_Diag%Archive_BudgetChemistryPBL,  &
                                 State_Diag%BudgetMass2,                 &
                                 RC )       
       CALL Compute_Budget_Diagnostics( am_I_Root,                           &
                                     State_Chm%Map_Advect,                   &
                                     DT_Chem,                                &
                                     State_Diag%Archive_BudgetChemistryFull, &
                                     State_Diag%Archive_BudgetChemistryTrop, &
                                     State_Diag%Archive_BudgetChemistryPBL,  &
                                     State_Diag%BudgetChemistryFull,         &
                                     State_Diag%BudgetChemistryTrop,         &
                                     State_Diag%BudgetChemistryPBL,          &
                                     State_Diag%BudgetMass1,                 &
                                     State_Diag%BudgetMass2,                 &
                                     RC )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Chemistry budget diagnostics error 2'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDIF
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
  SUBROUTINE RECOMPUTE_OD( am_I_Root, Input_Opt,  State_Met,                &
                           State_Chm, State_Diag, RC                       )
!
! !USES:
!
    ! References to F90 modules
    USE AEROSOL_MOD,    ONLY : AEROSOL_CONC
    USE AEROSOL_MOD,    ONLY : RDAER
    USE AEROSOL_MOD,    ONLY : SOILDUST
    USE DUST_MOD,       ONLY : RDUST_ONLINE
    USE DUST_MOD,       ONLY : RDUST_OFFLINE
    USE ErrCode_Mod
    USE ERROR_MOD,      ONLY : Debug_Msg
    USE Input_Opt_Mod,  ONLY : OptInput
    USE State_Chm_Mod,  ONLY : ChmState
    USE State_Diag_Mod, ONLY : DgnState
    USE State_Met_Mod,  ONLY : MetState
    USE TIME_MOD,       ONLY : GET_MONTH
    USE TIME_MOD,       ONLY : GET_YEAR
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
    TYPE(DgnState), INTENT(INOUT) :: State_Diag  ! Diagnostics State object
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
!  03 Nov 2017 - R. Yantosca - Now accept State_Diag as an argument
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    LOGICAL            :: IT_IS_A_FULLCHEM_SIM
    LOGICAL            :: IT_IS_AN_AEROSOL_SIM
    LOGICAL            :: LCARB, LCHEM,  LDUST
    LOGICAL            :: LPRT,  LSSALT, LSULF,      LSOA
    INTEGER            :: MONTH, YEAR,   WAVELENGTH

    ! Strings
    CHARACTER(LEN=255) :: ErrMsg, ThisLoc

    !=======================================================================
    ! RECOMPUTE_OD begins here!
    !=======================================================================

    ! Initialize
    RC      = GC_SUCCESS
    ErrMsg  = ''
    ThisLoc = ' -> at Recompute_OD  (in module GeosCore/chemistry_mod.F)'

    ! Get month and year
    MONTH                = GET_MONTH()
    YEAR                 = GET_YEAR()

    ! Copy fields from INPUT_OPT to local variables for use below
    LCARB                = Input_Opt%LCARB 
    LCHEM                = Input_Opt%LCHEM
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
             CALL AEROSOL_CONC( am_I_Root, Input_Opt,  State_Met,            &
                                State_Chm, State_Diag, RC                   )

             !==============================================================
             ! Call RDAER -- computes aerosol optical depths
             !==============================================================

             ! Calculate the AOD at the wavelength specified in jv_spec_aod
             WAVELENGTH = 1
             CALL RDAER( am_I_Root, Input_Opt,  State_Met,                   &
                         State_Chm, State_Diag, RC,                          &
                         MONTH,     YEAR,       WAVELENGTH                  )

             ! Trap potential errors
             IF ( RC /= GC_SUCCESS ) THEN
                ErrMsg = 'Error encountered in "RdAer"!'
                CALL GC_Error( ErrMsg, RC, ThisLoc )
                RETURN
             ENDIF

             !### Debug
             IF ( LPRT .and. am_I_Root ) THEN 
                CALL Debug_Msg( '### RECOMPUTE_OD: after RDAER' )
             ENDIF

             !==============================================================
             ! If LDUST is turned on, then we have online dust aerosol in
             ! GEOS-CHEM...so just pass SOILDUST to RDUST_ONLINE in order 
             ! to compute aerosol optical depth for FAST-JX, etc.
             !
             ! If LDUST is turned off, then we don't have online dust 
             ! aerosol in GEOS-CHEM...so read monthly-mean dust files
             ! from disk. (rjp, tdf, bmy, 4/1/04)
             !==============================================================
             IF ( LDUST ) THEN
                CALL RDUST_ONLINE( am_I_Root,  Input_Opt,  State_Met,        &
                                   State_Chm,  State_Diag, SOILDUST,         &
                                   WAVELENGTH, RC                           )

                ! Trap potential errors
                IF ( RC /= GC_SUCCESS ) THEN
                   ErrMsg = 'Error encountered in "Rdust_Online"!'
                   CALL GC_Error( ErrMsg, RC, ThisLoc )
                   RETURN
                ENDIF

#if  !defined( TOMAS )
             ELSE
                CALL RDUST_OFFLINE( am_I_Root, Input_Opt,  State_Met,        &
                                    State_Chm, State_Diag, MONTH,            &
                                    YEAR,      WAVELENGTH, RC               )

                ! Trap potential errors
                IF ( RC /= GC_SUCCESS ) THEN
                   ErrMsg = 'Error encountered in "Rdust_Offline"!'
                   CALL GC_Error( ErrMsg, RC, ThisLoc )
                   RETURN
                ENDIF
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
! !IROUTINE: chem_passive_species
!
! !DESCRIPTION: Subroutine RUN\_PASSIVE\_SPECIES performs loss chemistry 
!  on passive species with finite atmospheric lifetimes.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Chem_Passive_Species( am_I_Root, Input_Opt,                     &
                                   State_Met, State_Chm, RC                 ) 
!
! !USES:
!
    USE CMN_SIZE_Mod,   ONLY : IIPAR, JJPAR, LLPAR
    USE ErrCode_Mod
    USE Input_Opt_Mod,  ONLY : OptInput
    USE State_Chm_Mod,  ONLY : ChmState
    USE State_Met_Mod,  ONLY : MetState
    USE State_Chm_Mod,  ONLY : ind_ 
    USE Time_Mod,       ONLY : Get_Ts_Chem
!
! !INPUT PARAMETERS:
!
    LOGICAL,         INTENT(IN   )  :: am_I_Root   ! root CPU?
    TYPE(OptInput),  INTENT(IN   )  :: Input_Opt   ! Input options object
    TYPE(MetState),  INTENT(IN   )  :: State_Met   ! Meteorology state object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState),  INTENT(IN   )  :: State_Chm   ! Chemistry state object 
!
! !OUTPUT PARAMETERS:
!
    INTEGER,         INTENT(INOUT)  :: RC          ! Failure or success
!
! !REMARKS:
!
! !REVISION HISTORY: 
!  04 Sep 2015 - C. Keller   - Initial version 
!  03 Nov 2016 - C. Keller   - Moved to chemistry_mod
!  26 Jun 2017 - R. Yantosca - GC_ERROR is now contained in errcode_mod.F90
!  14 Jul 2017 - E. Lundgren - Remove dependency on passive_species_mod.F90
!  02 Aug 2017 - R. Yantosca - Turn off debug print unless ND70 is activated
!  13 Dec 2017 - R. Yantosca - Now apply decay only to those passive species
!                              with finite atmospheric lifetimes
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    LOGICAL             :: prtDebug
    INTEGER             :: I,       J,      L
    INTEGER             :: N,       GCId,   Id
    REAL(fp)            :: DT,      Decay,  Rate

    ! SAVEd scalars
    LOGICAL,  SAVE      :: First = .TRUE.

    ! Strings
    CHARACTER(LEN=255)  :: ErrMsg,  ThisLoc
!
! !DEFINED PARAMETERS:
!   
    REAL(fp), PARAMETER :: ln2 = 0.693147181E+00_fp

    !=======================================================================
    ! Chem_Passive_Species begins here!
    !=======================================================================

    ! Initialize
    RC       = GC_SUCCESS
    prtDebug = ( am_I_Root .and. Input_Opt%LPRT )
    ErrMsg   = ''
    ThisLoc  = &
       ' -> at Chem_Passive_Species (in module GeosCore/chemistry_mod.F)'

    DT       = GET_TS_CHEM() ! timestep in seconds

    !=======================================================================
    ! Apply decay loss rate only to those passive species that have a
    ! finite atmospheric lifetime (this speeds up execution)
    !=======================================================================

    ! Loop over all decaying passive species
    DO N = 1, Input_Opt%NPassive_Decay

       !----------------------------------
       ! Find the GEOS-Chem species Id
       !----------------------------------

       ! Get the Id of the species in the passive decay menu
       Id   = Input_Opt%Passive_DecayID(N)

       ! Convert this to a GEOS-Chem species Id number
       GcId = Ind_( TRIM( Input_Opt%PASSIVE_NAME(Id) ) )

       ! Make sure the model ID is valid
       IF ( GcId < 0 ) THEN
          ErrMsg = 'Could not find the GEOS-Chem species ID # '        // &
                   'for passive species : '                            // &
                   TRIM( Input_Opt%PASSIVE_NAME(Id) )
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       !----------------------------------
       ! Compute the decay rate
       !----------------------------------

       ! Compute the decay rate for each passive species
       Decay = ln2 / Input_Opt%PASSIVE_TAU(Id)
       Rate  = EXP( - DT * Decay )

       !### Debug output
       IF ( First ) THEN
          IF ( prtDebug ) THEN
             WRITE( 6,100 ) ADJUSTL( Input_Opt%PASSIVE_NAME(Id) ),           &
                            GcId, Rate
 100         FORMAT( '     -  Pass. species name, Id, loss rate:',           &
                      a15, i5, 1x, es13.6 )
          ENDIF
          First = .FALSE.
       ENDIF

       !----------------------------------
       ! Apply loss
       !----------------------------------

       !$OMP PARALLEL DO                  &
       !$OMP DEFAULT( SHARED            ) &
       !$OMP PRIVATE( I, J, L           )
       DO L = 1, LLPAR
       DO J = 1, JJPAR
       DO I = 1, IIPAR
          State_Chm%Species(I,J,L,GcId) = State_Chm%Species(I,J,L,GcId)      &
                                        * Rate
       ENDDO
       ENDDO
       ENDDO
       !$OMP END PARALLEL DO

    ENDDO
 
  END SUBROUTINE Chem_Passive_Species
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
  SUBROUTINE Init_Chemistry( am_I_Root, Input_Opt,                           &
                             State_Chm, State_Diag, RC                      ) 
!
! !USES:
!
    USE ErrCode_Mod
    USE FAST_JX_MOD,    ONLY : Init_FJX
    USE FlexChem_Mod,   ONLY : Init_FlexChem
    USE Input_Opt_Mod,  ONLY : OptInput
    USE State_Chm_Mod,  ONLY : ChmState
    USE State_Chm_Mod,  ONLY : Ind_
    USE State_Diag_Mod, ONLY : DgnState
!
! !INPUT PARAMETERS:
!
    LOGICAL,        INTENT(IN)     :: am_I_Root   ! Is this the root CPU?
!
! !INPUT/OUTPUT PARAMETERS: 
!
    TYPE(OptInput), INTENT(INOUT)  :: Input_Opt   ! Input Options object
    TYPE(ChmState), INTENT(INOUT)  :: State_Chm   ! Chemistry State object
    TYPE(DgnState), INTENT(INOUT)  :: State_Diag  ! Diagnostics State object
    INTEGER,        INTENT(INOUT)  :: RC          ! Success or failure?
!
! !REVISION HISTORY: 
!  19 May 2014 - C. Keller   - Initial version (stripped from do_chemistry
!                              and chemdr.F)
!  20 Jun 2014 - R. Yantosca - Now pass Input_Opt to INIT_FJX
!  23 Jun 2016 - R. Yantosca - Remove call to SETTRACE, it's obsolete
!  03 Nov 2017 - R. Yantosca - Now accept State_Diag as an argument
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:

    ! SAVEd scalars
    LOGICAL, SAVE      :: FIRST = .TRUE.

    ! Strings
    CHARACTER(LEN=255) :: ErrMsg, ThisLoc

    !=======================================================================
    ! INIT_CHEMISTRY begins here!
    !=======================================================================
    
    ! Initialize
    RC       = GC_SUCCESS
    ErrMsg   = ''
    ThisLoc  = ' -> at Init_Chemistry  (in module GeosCore/chemistry_mod.F)'

    ! Skip if we have already done this
    IF ( FIRST ) THEN

       ! Adjust first flag
       FIRST  = .FALSE.

       ! Define species ID's
       id_DST1 = Ind_( 'DST1' )
       id_NK1  = Ind_( 'NK1'  )

       !--------------------------------------------------------------------
       ! Initialize FlexChem
       !--------------------------------------------------------------------
       CALL Init_FlexChem( am_I_Root, Input_Opt, State_Chm, State_Diag, RC  )

       ! Trap potential errors
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Error encountered in "Init_FlexChem"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! Initialize Fast-JX photolysis
       !--------------------------------------------------------------------
       CALL Init_FJX( am_I_Root, Input_Opt, State_Chm, State_Diag, RC       )

       ! Trap potential errors
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Error encountered in "Init_FJX"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

    ENDIF

  END SUBROUTINE Init_Chemistry
!EOC
END MODULE Chemistry_Mod
