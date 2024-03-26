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
  USE Timers_Mod       ! For GEOS-Chem timers (optional)

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC  :: Do_Chemistry
  PUBLIC  :: Recompute_OD
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE TYPES:
!
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
  SUBROUTINE Do_Chemistry( Input_Opt,  State_Chm, State_Diag,                &
                           State_Grid, State_Met, RC                        )
!
! !USES:
!
    USE AEROSOL_MOD,      ONLY : AEROSOL_CONC
    USE AEROSOL_MOD,      ONLY : RDAER
    USE CARBON_MOD,       ONLY : CHEMCARBON
    USE Carbon_Gases_Mod, ONLY : Chem_Carbon_Gases
    USE Diagnostics_Mod,  ONLY : Compute_Budget_Diagnostics
    USE DUST_MOD,         ONLY : CHEMDUST
    USE DUST_MOD,         ONLY : RDUST_ONLINE
    USE ErrCode_Mod
    USE ERROR_MOD
    USE FullChem_Mod,     ONLY : Do_FullChem
    USE GLOBAL_CH4_MOD,   ONLY : CHEMCH4
    USE Input_Opt_Mod,    ONLY : OptInput
    USE AEROSOL_THERMODYNAMICS_MOD,  ONLY : DO_ATE
    USE LINEAR_CHEM_MOD,  ONLY : DO_LINEAR_CHEM
    USE MERCURY_MOD,      ONLY : CHEMMERCURY
    USE POPS_MOD,         ONLY : CHEMPOPS
    USE RnPbBe_MOD,       ONLY : CHEMRnPbBe
    USE RPMARES_MOD,      ONLY : DO_RPMARES
    USE SEASALT_MOD,      ONLY : CHEMSEASALT
    USE SULFATE_MOD,      ONLY : CHEMSULFATE
    USE State_Chm_Mod,    ONLY : ChmState
    USE State_Chm_Mod,    ONLY : Ind_
    USE State_Diag_Mod,   ONLY : DgnState
    USE State_Grid_Mod,   ONLY : GrdState
    USE State_Met_Mod,    ONLY : MetState
    USE TAGGED_CO_MOD,    ONLY : CHEM_TAGGED_CO
    USE TAGGED_O3_MOD,    ONLY : CHEM_TAGGED_O3
    USE TIME_MOD,         ONLY : GET_TS_CHEM
    USE Tracer_Mod,       ONLY : Tracer_Sink_Phase
    USE UCX_MOD,          ONLY : CALC_STRAT_AER
    USE UnitConv_Mod
#ifdef APM
    USE APM_INIT_MOD,     ONLY : APMIDS
    USE APM_DRIV_MOD,     ONLY : PSO4GAS
    USE APM_DRIV_MOD,     ONLY : AERONUM
    USE APM_DRIV_MOD,     ONLY : APM_DRIV
#endif
#ifdef TOMAS
    USE TOMAS_MOD,        ONLY : DO_TOMAS  !(win, 7/14/09)
#endif

!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
    TYPE(GrdState), INTENT(IN)    :: State_Grid  ! Grid State object
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
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! SAVEd scalars
    INTEGER, SAVE      :: id_DST1, id_NK01, id_CO2   ! Species ID flags

    ! Scalars
    INTEGER            :: N_TROP, N
    INTEGER            :: MONTH
    INTEGER            :: YEAR
    INTEGER            :: WAVELENGTH
    INTEGER            :: TS_Chem
    REAL(f8)           :: DT_Chem, sDTFC, fDTFC
#ifdef APM
    INTEGER            :: I,J,L
    REAL*8             :: CONCTMPSO4(State_Grid%NX,State_Grid%NY,State_Grid%NZ)
#endif

    ! SAVEd scalars
    LOGICAL, SAVE      :: FIRST = .TRUE.

    ! Strings
    INTEGER            :: previous_units
    CHARACTER(LEN=255) :: ErrMsg,  ThisLoc

    !=======================================================================
    ! DO_CHEMISTRY begins here!
    !=======================================================================

    ! Initialize
    RC        = GC_SUCCESS
    ErrMsg    = ''
    ThisLoc   = ' -> at Do_Chemistry  (in module GeosCore/chemistry_mod.F90)'

    ! Save species ID"s on first call
    IF ( FIRST ) THEN
       id_DST1 = Ind_('DST1')
       id_NK01 = Ind_('NK01')
       id_CO2  = Ind_('CO2' )
    ENDIF

    !========================================================================
    ! Chemistry budget diagnostics - Part 1 of 2
    !========================================================================
    IF ( State_Diag%Archive_BudgetChemistry ) THEN

       ! Get initial column masses (full, trop, PBL)
       CALL Compute_Budget_Diagnostics(                                      &
            Input_Opt   = Input_Opt,                                         &
            State_Chm   = State_Chm,                                         &
            State_Diag  = State_Diag,                                        &
            State_Grid  = State_Grid,                                        &
            State_Met   = State_Met,                                         &
            isFull      = State_Diag%Archive_BudgetChemistryFull,            &
            diagFull    = NULL(),                                            &
            mapDataFull = State_Diag%Map_BudgetChemistryFull,                &
            isTrop      = State_Diag%Archive_BudgetChemistryTrop,            &
            diagTrop    = NULL(),                                            &
            mapDataTrop = State_Diag%Map_BudgetChemistryTrop,                &
            isPBL       = State_Diag%Archive_BudgetChemistryPBL,             &
            diagPBL     = NULL(),                                            &
            mapDataPBL  = State_Diag%Map_BudgetChemistryPBL,                 &
            isLevs      = State_Diag%Archive_BudgetChemistryLevs,            &
            diagLevs    = NULL(),                                            &
            mapDataLevs = State_Diag%Map_BudgetChemistryLevs,                &
            colMass     = State_Diag%BudgetColumnMass,                       &
            before_op   = .TRUE.,                                            &
            msg         = 'Compute budget diag (1)'//TRIM(ThisLoc),          &
            RC          = RC                                                )

       ! Trap potential errors
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Chemistry budget diagnostics error 1'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDIF

    !========================================================================
    ! Convert species units to [kg] for chemistry (ewl, 8/12/15)
    !========================================================================

    ! Here, units are still in mol/mol dry.   For fullchem-simulation only,
    ! set CO2 to 421 ppm (or 421e-6 mol/mol dry) since this is the global
    ! average value. This is necessary to reduce the error norm in KPP.
    ! See https://github.com/geoschem/geos-chem/issues/1529.
    IF ( Input_Opt%ITS_A_FULLCHEM_SIM ) THEN
       State_Chm%Species(id_CO2)%Conc = 421.0e-6_fp
    ENDIF

    ! Halt "All chemistry" timer (so that diags can be timed separately)
    IF ( Input_Opt%useTimers ) THEN
       CALL Timer_End( "All chemistry", RC )
    ENDIF

    ! Convert units from mol/mol dry to kg
    CALL Convert_Spc_Units(                                                  &
         Input_Opt      = Input_Opt,                                         &
         State_Chm      = State_Chm,                                         &
         State_Grid     = State_Grid,                                        &
         State_Met      = State_Met,                                         &
         new_units      = KG_SPECIES,                                        &
         previous_units = previous_units,                                    &
         RC             = RC                                                )

    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Unit conversion error (kg/kg dry -> kg)'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    ! Start "All chemistry" timer again
    IF ( Input_Opt%useTimers ) THEN
       CALL Timer_Start( "All chemistry", RC )
    ENDIF

    !========================================================================
    ! If Input_Opt%LCHEM=T then call the chemistry subroutines
    !========================================================================
    IF ( Input_Opt%LCHEM ) THEN

       !=====================================================================
       ! Full-chemistry simulations:
       !
       ! (1) Benchmark; (2) Standard; (3) SimpleSOA; (4) complexSOA,
       ! (5) complexSOA-SVPOA; (6) aciduptake; (7) marinePOA
       !=====================================================================
       IF ( Input_Opt%ITS_A_FULLCHEM_SIM ) THEN

          IF ( Input_Opt%useTimers ) THEN
             CALL Timer_Start  ( "=> Aerosol chem", RC )
          ENDIF

          !------------------------------------------------------------------
          ! Dry-run sulfate chem to get cloud pH
          !------------------------------------------------------------------
          IF ( Input_Opt%LSULF ) THEN

             ! Calculate stratospheric aerosol properties
             CALL Calc_Strat_Aer( Input_Opt  = Input_Opt,                    &
                                  State_Chm  = State_Chm,                    &
                                  State_Grid = State_Grid,                   &
                                  State_Met  = State_Met,                    &
                                  RC         = RC                           )

             ! Trap potential errors
             IF ( RC /= GC_SUCCESS ) THEN
                ErrMsg = 'Error encountered in "Calc_Strat_Aer"!'
                CALL GC_Error( ErrMsg, RC, ThisLoc )
                RETURN
             ENDIF

             ! Compute aerosol concentrations (needed for AOD computations)
             CALL Aerosol_Conc( Input_Opt  = Input_Opt,                      &
                                State_Chm  = State_Chm,                      &
                                State_Diag = State_Diag,                     &
                                State_Grid = State_Grid,                     &
                                State_Met  = State_Met,                      &
                                RC         = RC                             )

             ! Trap potential errors
             IF ( RC /= GC_SUCCESS ) THEN
                ErrMsg = 'Error encountered in "AEROSOL_CONC"!'
                CALL GC_Error( ErrMsg, RC, ThisLoc )
                RETURN
             ENDIF

          ENDIF

          !------------------------------------------------------------------
          ! Call RDAER
          !------------------------------------------------------------------
          waveLength = 0
          CALL RdAer( Input_Opt  = Input_Opt,                                &
                      State_Chm  = State_Chm,                                &
                      State_Diag = State_Diag,                               &
                      State_Grid = State_Grid,                               &
                      State_Met  = State_Met,                                &
                      month      = month,                                    &
                      year       = year,                                     &
                      odSwitch   = waveLength,                               &
                      RC         = RC                                       )

          ! Trap potential errors
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = 'Error encountered in "RDAER"!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF

          !==================================================================
          ! If LDUST is turned on, then we have online dust aerosol in
          ! GEOS-CHEM.
          ! If LDUST is turned off, then we do not have online dust aerosol
          ! in GEOS-CHEM...so read monthly-mean dust files from disk.
          ! (rjp, tdf, bmy, 4/1/04)
          !==================================================================
          IF ( Input_Opt%LDUST ) THEN
             CALL RDust_Online( Input_Opt  = Input_Opt,                      &
                                State_Chm  = State_Chm,                      &
                                State_Diag = State_Diag,                     &
                                State_Grid = State_Grid,                     &
                                State_Met  = State_Met,                      &
                                odSwitch   = waveLength,                     &
                                RC         = RC                             )

             ! Trap potential errors
             IF ( RC /= GC_SUCCESS ) THEN
                ErrMsg = 'Error encountered in "RDUST_ONLINE"!'
                CALL GC_Error( ErrMsg, RC, ThisLoc )
                RETURN
             ENDIF
          ENDIF

          !------------------------------------------------------------------
          ! Dry-run sulfate chem to get cloud pH
          !------------------------------------------------------------------
          IF ( Input_Opt%LSULF ) THEN

             ! Dry run only
             CALL ChemSulfate( Input_Opt  = Input_Opt,                       &
                               State_Chm  = State_Chm,                       &
                               State_Diag = State_Diag,                      &
                               State_Grid = State_Grid,                      &
                               State_Met  = State_Met,                       &
                               FullRun    = .FALSE.,                         &
                               RC         = RC                              )

             ! Trap potential errors
             IF ( RC /= GC_SUCCESS ) THEN
                ErrMsg = 'Error encountered in "ChemSulfate"!'
                CALL GC_Error( ErrMsg, RC, ThisLoc )
                RETURN
             ENDIF

             !---------------------------------------------------------------
             ! Do aerosol thermodynamic equilibrium
             !---------------------------------------------------------------
             IF ( Input_Opt%LSSALT ) THEN

#ifndef APM
                ! ISORROPIA/HETP take Na+, Cl- into account
                CALL Do_ATE( Input_Opt  = Input_Opt,                 &
                             State_Chm  = State_Chm,                 &
                             State_Diag = State_Diag,                &
                             State_Grid = State_Grid,                &
                             State_Met  = State_Met,                 &
                             RC         = RC                        )

                ! Trap potential errors
                IF ( RC /= GC_SUCCESS ) THEN
                   ErrMsg = 'Error encountered in "Do_ATE"!'
                   CALL GC_Error( ErrMsg, RC, ThisLoc )
                   RETURN
                ENDIF
#endif

             ELSE

#ifdef APM
                ! Exit with error if RPMARES + APM is selected
                ErrMsg = 'Warning: APM does not want to use DO_RPMARES'
                CALL GC_Error( ErrMsg, RC, ThisLoc )
                RETURN
#endif

                ! RPMARES does not take Na+, Cl- into account
                CALL Do_RPMARES( Input_Opt  = Input_Opt,                     &
                                 State_Chm  = State_Chm,                     &
                                 State_Grid = State_Grid,                    &
                                 State_Met  = State_Met,                     &
                                 RC         = RC                            )
             ENDIF

          ENDIF

#ifdef APM
          ! Save SO4 concentration before chemistry
          N          = APMIDS%id_SO4
          CONCTMPSO4 = State_Chm%Species(N)%Conc

          CALL AERONUM( Input_Opt  = Input_Opt,                              &
                        State_Chm  = State_Chm,                              &
                        State_Diag = State_Diag,                             &
                        State_Grid = State_Grid,                             &
                        State_Met  = State_Met,                              &
                        RC         = RC                                     )

          ! Trap potential errors
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = 'Error encountered in "AERONUM"!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF
#endif

          IF ( Input_Opt%useTimers ) THEN
             CALL Timer_End( "=> Aerosol chem", RC )
          ENDIF

          !------------------------------------------------------------------
          ! Call gas-phase chemistry
          !------------------------------------------------------------------
          IF ( Input_Opt%useTimers ) THEN
             CALL Timer_Start( "=> Gas-phase chem", RC )
          ENDIF

          ! Solve the KPP-generated mechanism
          CALL Do_FullChem( Input_Opt  = Input_Opt,                          &
                            State_Chm  = State_Chm,                          &
                            State_Diag = State_Diag,                         &
                            State_Grid = State_Grid,                         &
                            State_Met  = State_Met,                          &
                            RC         = RC                                 )

          ! Trap potential errors
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = 'Error encountered in "Do_FullChem"!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF

          IF ( Input_Opt%useTimers ) THEN
             CALL Timer_End( "=> Gas-phase chem", RC )
          ENDIF

          !------------------------------------------------------------------
          ! Linearized chemistry above chemistry grid
          !------------------------------------------------------------------
          IF ( Input_Opt%LINEAR_CHEM ) THEN

             IF ( Input_Opt%useTimers ) THEN
                CALL Timer_Start( "=> Linearized chem", RC )
             ENDIF

             ! Do linearized chemistry for the mesosphere
             CALL Do_Linear_Chem( Input_Opt  = Input_Opt,                    &
                                  State_Chm  = State_Chm,                    &
                                  State_Grid = State_Grid,                   &
                                  State_Met  = State_Met,                    &
                                  errCode    = RC                           )

             ! Trap potential errors
             IF ( RC /= GC_SUCCESS ) THEN
                ErrMsg = 'Error encountred in "Do_LinearChem"!'
                CALL GC_Error( ErrMsg, RC, ThisLoc )
                RETURN
             ENDIF

             ! Make sure all units are still in kg
             IF ( .not. Check_Units( State_Chm, KG_SPECIES ) ) THEN
                ErrMsg = 'Incorrect species after calling "Do_Linear_Chem"!'
                CALL GC_Error( ErrMsg, RC, ThisLoc )
                RETURN
             ENDIF

             IF ( Input_Opt%useTimers ) THEN
                CALL Timer_End( "=> Linearized chem", RC )
             ENDIF

          ENDIF

#ifdef APM
          ! Obtain SO4 production after chemistry
          N = APMIDS%id_SO4
          !$OMP PARALLEL DO         &
          !$OMP DEFAULT( SHARED   ) &
          !$OMP PRIVATE( I, J, L  ) &
          !$OMP SCHEDULE( DYNAMIC )
          DO L = 1, State_Grid%NZ
          DO J = 1, State_Grid%NY
          DO I = 1, State_Grid%NX
             IF ( State_Chm%Species(N)%Conc(I,J,L) > CONCTMPSO4(I,J,L) ) THEN
                PSO4GAS(I,J,L) = State_Chm%Species(N)%Conc(I,J,L)                  &
                               - CONCTMPSO4(I,J,L)
             ELSE
                PSO4GAS(I,J,L) = 0.D0
             ENDIF
          ENDDO
          ENDDO
          ENDDO
          !$OMP END PARALLEL DO
#endif

          IF ( Input_Opt%useTimers ) THEN
             CALL Timer_Start( "=> Aerosol chem", RC )
          ENDIF

          !------------------------------------------------------------------
          ! Do seasalt aerosol chemistry
          !------------------------------------------------------------------
          IF ( Input_Opt%LSSALT ) THEN
             CALL ChemSeaSalt( Input_Opt  = Input_Opt,                       &
                               State_Chm  = State_Chm,                       &
                               State_Diag = State_Diag,                      &
                               State_Grid = State_Grid,                      &
                               State_Met  = State_Met,                       &
                               RC         = RC                              )

             ! Trap potential errors
             IF ( RC /= GC_SUCCESS ) THEN
                ErrMsg = 'Error encountered in "ChemSeaSalt"!'
                CALL GC_Error( ErrMsg, RC, ThisLoc )
                RETURN
             ENDIF
          ENDIF

          !------------------------------------------------------------------
          ! Recalculate PSC properties
          !------------------------------------------------------------------
          CALL Calc_Strat_Aer( Input_Opt  = Input_Opt,                       &
                               State_Chm  = State_Chm,                       &
                               State_Grid = State_Grid,                      &
                               State_Met  = State_Met,                       &
                               RC         = RC                              )

          ! Trap potential errors
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = 'Error encountered in "Calc_Strat_Aer"!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF

          !------------------------------------------------------------------
          ! Also do sulfate chemistry
          !------------------------------------------------------------------
          IF ( Input_Opt%LSULF ) THEN

             ! Do sulfate chemistry
             CALL ChemSulfate( Input_Opt  = Input_Opt,                       &
                               State_Chm  = State_Chm,                       &
                               State_Diag = State_Diag,                      &
                               State_Grid = State_Grid,                      &
                               State_Met  = State_Met,                       &
                               FullRun    = .TRUE.,                          &
                               RC         = RC                              )

             ! Trap potential errors
             IF ( RC /= GC_SUCCESS ) THEN
                ErrMsg = 'Error encountered after calling "ChemSulfate"!'
                CALL GC_Error( ErrMsg, RC, ThisLoc )
                RETURN
             ENDIF

             ! Make sure all units are still in kg
             IF ( .not. Check_Units( State_Chm, KG_SPECIES ) ) THEN
                ErrMsg = 'Incorrect species after calling "ChemSulfate"!'
                CALL GC_Error( ErrMsg, RC, ThisLoc )
                RETURN
             ENDIF

          ENDIF

          !------------------------------------------------------------------
          ! Do carbonaceous aerosol chemistry
          !------------------------------------------------------------------
          IF ( Input_Opt%LCARB ) THEN
             CALL ChemCarbon( Input_Opt  = Input_Opt,                        &
                              State_Chm  = State_Chm,                        &
                              State_Diag = State_Diag,                       &
                              State_Grid = State_Grid,                       &
                              State_Met  = State_Met,                        &
                              RC         = RC                               )

             ! Trap potential errors
             IF ( RC /= GC_SUCCESS ) THEN
                ErrMsg = 'Error encountered in "ChemCarbon"!'
                CALL GC_Error( ErrMsg, RC, ThisLoc )
                RETURN
             ENDIF
          ENDIF

          !------------------------------------------------------------------
          ! Do dust aerosol chemistry/removal
          !------------------------------------------------------------------
          IF ( Input_Opt%LDUST .AND. id_DST1 > 0 ) THEN
             CALL ChemDust( Input_Opt  = Input_Opt,                          &
                            State_Chm  = State_Chm,                          &
                            State_Diag = State_Diag,                         &
                            State_Grid = State_Grid,                         &
                            State_Met  = State_Met,                          &
                            RC         = RC                                 )

             ! Trap potential errors
             IF ( RC /= GC_SUCCESS ) THEN
                ErrMsg = 'Error encountered in "ChemDust"!'
                CALL GC_Error( ErrMsg, RC, ThisLoc )
                RETURN
             ENDIF
          ENDIF

#ifdef APM
          !------------------------------------------------------------------
          ! Do APM aerosol microphysics
          !------------------------------------------------------------------
          CALL APM_Driv( Input_Opt  = Input_Opt,                             &
                         State_Chm  = State_Chm,                             &
                         State_Diag = State_Diag,                            &
                         State_Grid = State_Grid,                            &
                         State_Met  = State_Met,                             &
                         RC         = RC                                    )

          ! Trap potential errors
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = 'Error encountered in routine "APM_DRIV"!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF
#endif

#ifdef TOMAS
          !------------------------------------------------------------------
          ! Do TOMAS aerosol microphysics and dry dep
          !------------------------------------------------------------------
          IF ( id_NK01 > 0 ) THEN
             CALL Do_TOMAS( Input_Opt  = Input_Opt,                          &
                            State_Chm  = State_Chm,                          &
                            State_Diag = State_Diag,                         &
                            State_Grid = State_Grid,                         &
                            State_Met  = State_Met,                          &
                            RC         = RC                                 )

             ! Trap potential errors
             IF ( RC /= GC_SUCCESS ) THEN
                ErrMsg = 'Error encountered in "Do_TOMAS"!'
                CALL GC_Error( ErrMsg, RC, ThisLoc )
                RETURN
             ENDIF

             ! Check units (ewl, 10/5/15)
             IF ( .not. Check_Units( State_Chm, KG_SPECIES ) ) THEN
                ErrMsg = 'Not all species have units "kg"!'
                CALL GC_Error( ErrMsg, RC, ThisLoc )
             ENDIF
          ENDIF
#endif

          IF ( Input_Opt%useTimers ) THEN
             CALL Timer_End( "=> Aerosol chem", RC )
          ENDIF

       !=====================================================================
       ! Aerosol-only simulation
       !=====================================================================
       ELSE IF ( Input_Opt%ITS_AN_AEROSOL_SIM ) THEN

          IF ( Input_Opt%useTimers ) THEN
             CALL Timer_Start( "=> Aerosol chem", RC )
          ENDIF

          !------------------------------------------------------------------
          ! Compute aerosol & dust concentrations [kg/m3]
          ! (NOTE: SOILDUST in "aerosol_mod.F90" is computed here)
          !------------------------------------------------------------------
          CALL Aerosol_Conc( Input_Opt  = Input_Opt,                         &
                             State_Chm  = State_Chm,                         &
                             State_Diag = State_Diag,                        &
                             State_Grid = State_Grid,                        &
                             State_Met  = State_Met,                         &
                             RC         = RC                                )

          ! Trap potential errors
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = 'Error encountered in "Aerosol_Conc"!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF

          !------------------------------------------------------------------
          ! Compute AOD's and surface areas at 999 nm
          !------------------------------------------------------------------
          month      = 0
          year       = 0
          waveLength = 0
          CALL RdAer( Input_Opt  = Input_Opt,                                &
                      State_Chm  = State_Chm,                                &
                      State_Diag = State_Diag,                               &
                      State_Grid = State_Grid,                               &
                      State_Met  = State_Met,                                &
                      month      = month,                                    &
                      year       = year,                                     &
                      ODswitch   = Wavelength,                               &
                      RC         = RC                                       )

          ! Trap potential errors
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = 'Error encountered in "RdAer"!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF

          !------------------------------------------------------------------
          ! Aerosol Thermodynamic Equilibrium
          !------------------------------------------------------------------
          IF ( Input_Opt%LSULF ) THEN
             IF ( Input_Opt%LSSALT ) THEN

#ifndef APM
                ! ISORROPIA/HETP take Na+, Cl- into account
                CALL Do_ATE( Input_Opt  = Input_Opt,                 &
                             State_Chm  = State_Chm,                 &
                             State_Diag = State_Diag,                &
                             State_Grid = State_Grid,                &
                             State_Met  = State_Met,                 &
                             RC         = RC                        )
#endif

                ! Trap potential errors
                IF ( RC /= GC_SUCCESS ) THEN
                   ErrMsg = 'Error encountered in "Do_ATE"!'
                   CALL GC_Error( ErrMsg, RC, ThisLoc )
                   RETURN
                ENDIF

             ELSE

#ifdef APM
                ! Exit with error if RPMARES + APM is selected
                ErrMsg = 'Warning: APM does not want to use DO_RPMARES'
                CALL GC_Error( ErrMsg, RC, ThisLoc )
                RETURN
#endif

                ! RPMARES does not take Na+, Cl- into account
                ! (skip for crystalline & aqueous offline run)
                CALL Do_RPMARES( Input_Opt  = Input_Opt,                     &
                                 State_Chm  = State_Chm,                     &
                                 State_Grid = State_Grid,                    &
                                 State_Met  = State_Met,                     &
                                 RC         = RC                            )

                ! Trap potential errors
                IF ( RC /= GC_SUCCESS ) THEN
                   ErrMsg = 'Error encountered in "Do_RPMARES"!'
                   CALL GC_Error( ErrMsg, RC, ThisLoc )
                   RETURN
                ENDIF
             ENDIF
          ENDIF

          !------------------------------------------------------------------
          ! Seasalt Aerosols
          !------------------------------------------------------------------
          IF ( Input_Opt%LSSALT ) THEN
             CALL ChemSeaSalt( Input_Opt  = Input_Opt,                       &
                               State_Chm  = State_Chm,                       &
                               State_Diag = State_Diag,                      &
                               State_Grid = State_Grid,                      &
                               State_Met  = State_Met,                       &
                               RC         = RC                              )

             ! Trap potential errors
             IF ( RC /= GC_SUCCESS ) THEN
                ErrMsg = 'Error encountered in "ChemSeaSalt"!'
                CALL GC_Error( ErrMsg, RC, ThisLoc )
                RETURN
             ENDIF
          ENDIF

          !------------------------------------------------------------------
          ! Sulfate aerosols
          !------------------------------------------------------------------
          IF ( Input_Opt%LSULF ) THEN

             ! Do sulfate chemistry
             CALL ChemSulfate( Input_Opt  = Input_Opt,                       &
                               State_Chm  = State_Chm,                       &
                               State_Diag = State_Diag,                      &
                               State_Grid = State_Grid,                      &
                               State_Met  = State_Met,                       &
                               FullRun    = .TRUE.,                          &
                               RC         = RC                              )

             ! Trap potential errors
             IF ( RC /= GC_SUCCESS ) THEN
                ErrMsg = 'Error encountered in "ChemSulfate"!'
                CALL GC_Error( ErrMsg, RC, ThisLoc )
                RETURN
             ENDIF
          ENDIF

          !------------------------------------------------------------------
          ! Carbon and Secondary Organic Aerosols
          !------------------------------------------------------------------
          IF ( Input_Opt%LCARB ) THEN
             CALL ChemCarbon( Input_Opt  = Input_Opt,                        &
                              State_Chm  = State_Chm,                        &
                              State_Diag = State_Diag,                       &
                              State_Grid = State_Grid,                       &
                              State_Met  = State_Met,                        &
                              RC         = RC                               )

             ! Trap potential errors
             IF ( RC /= GC_SUCCESS ) THEN
                ErrMsg = 'Error encountered in ""!'
                CALL GC_Error( ErrMsg, RC, ThisLoc )
                RETURN
             ENDIF
          ENDIF

          !------------------------------------------------------------------
          ! Mineral Dust Aerosols
          !------------------------------------------------------------------
          IF ( Input_Opt%LDUST ) THEN

             ! Do dust aerosol chemistry
             CALL ChemDust( Input_Opt  = Input_Opt,                          &
                            State_Chm  = State_Chm,                          &
                            State_Diag = State_Diag,                         &
                            State_Grid = State_Grid,                         &
                            State_Met  = State_Met,                          &
                            RC         = RC                                 )

             ! Trap potential errors
             IF ( RC /= GC_SUCCESS ) THEN
                ErrMsg = 'Error encountered in "ChemDust"!'
                CALL GC_Error( ErrMsg, RC, ThisLoc )
                RETURN
             ENDIF

             ! Compute dust OD's & surface areas
             WAVELENGTH = 0
             CALL Rdust_Online( Input_Opt  = Input_Opt,                   &
                                State_Chm  = State_Chm,                   &
                                State_Diag = State_Diag,                  &
                                State_Grid = State_Grid,                  &
                                State_Met  = State_Met,                   &
                                ODswitch   = waveLength,                  &
                                RC         = RC                          )

             ! Trap potential errors
             IF ( RC /= GC_SUCCESS ) THEN
                ErrMsg = 'Error encountered in "Rdust_Online"!'
                CALL GC_Error( ErrMsg, RC, ThisLoc )
                RETURN
             ENDIF
          ENDIF

          IF ( Input_Opt%useTimers ) THEN
             CALL Timer_End( "=> Aerosol chem", RC )
          ENDIF

       !=====================================================================
       ! Transport Tracers
       !=====================================================================
       ELSE IF ( Input_Opt%ITS_A_TRACER_SIM ) THEN

          ! Do Rn-Pb-Be chemistry
          CALL ChemRnPbBe( Input_Opt  = Input_Opt,                           &
                           State_Chm  = State_Chm,                           &
                           State_Diag = State_Diag,                          &
                           State_Grid = State_Grid,                          &
                           State_Met  = State_Met,                           &
                           RC         = RC                                  )

          ! Trap potential errors
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = 'Error encountered in "ChemRnPbBe"!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF

          ! Apply sinks for transport tracers
          CALL Tracer_Sink_Phase( Input_Opt  = Input_Opt,                 &
                                  State_Chm  = State_Chm,                 &
                                  State_Grid = State_Grid,                &
                                  State_Met  = State_Met,                 &
                                  RC         = RC                        )

          ! Trap potential errors
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = 'Error encountered in "Tracer_Sink_Phase"!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF

       !=====================================================================
       ! Tagged O3
       !=====================================================================
       ELSE IF ( Input_Opt%ITS_A_TAGO3_SIM ) THEN

          !------------------------------------------------------------------
          ! Do Tagged O3 chemistry
          !------------------------------------------------------------------
          CALL Chem_Tagged_O3( Input_Opt  = Input_Opt,                       &
                               State_Chm  = State_Chm,                       &
                               State_Diag = State_Diag,                      &
                               State_Grid = State_Grid,                      &
                               State_Met  = State_Met,                       &
                               RC         = RC                              )

          ! Trap potential errors
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = 'Error encountered in "Chem_Tagged_O3"!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF

          !------------------------------------------------------------------
          ! Linearized chemistry
          !------------------------------------------------------------------
          IF ( Input_Opt%LINEAR_CHEM ) THEN

             IF ( Input_Opt%useTimers ) THEN
                CALL Timer_Start( "=> Linearized chem", RC )
             ENDIF

             ! Do LINOZ for Ozone
             CALL Do_Linear_Chem( Input_Opt  = Input_Opt,                    &
                                  State_Chm  = State_Chm,                    &
                                  State_Grid = State_Grid,                   &
                                  State_Met  = State_Met,                    &
                                  errCode    = RC                           )

             ! Trap potential errors
             IF ( RC /= GC_SUCCESS ) THEN
                ErrMsg = 'Error encountered in "Do_Linear_Chem"!'
                CALL GC_Error( ErrMsg, RC, ThisLoc )
                RETURN
             ENDIF

             IF ( Input_Opt%useTimers ) THEN
                CALL Timer_End( "=> Linearized chem", RC )
             ENDIF

          ENDIF

       !=====================================================================
       ! Tagged CO
       !=====================================================================
       ELSE IF ( Input_Opt%ITS_A_TAGCO_SIM ) THEN

          ! Do tagged CO chemistry
          CALL Chem_Tagged_CO( Input_Opt  = Input_Opt,                       &
                               State_Chm  = State_Chm,                       &
                               State_Diag = State_Diag,                      &
                               State_Grid = State_Grid,                      &
                               State_Met  = State_Met,                       &
                               RC         = RC                              )

          ! Trap potential errors
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = 'Error encountered in "Chem_Tagged_CO"!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF

       !=====================================================================
       ! CH4
       !=====================================================================
       ELSE IF ( Input_Opt%ITS_A_CH4_SIM ) THEN

          ! Do CH4 chemistry
          CALL ChemCh4( Input_Opt  = Input_Opt,                              &
                        State_Chm  = State_Chm,                              &
                        State_Diag = State_Diag,                             &
                        State_Grid = State_Grid,                             &
                        State_Met  = State_Met,                              &
                        RC         = RC                                     )

          ! Trap potential errors
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = 'Error encountered in "ChemCh4"!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF

       !=====================================================================
       ! Carbon gases (configure with -DMECH=carbon)
       !=====================================================================
       ELSE IF ( Input_Opt%ITS_A_CARBON_SIM ) THEN

          IF ( Input_Opt%useTimers ) THEN
             CALL Timer_Start( "=> Gas-phase chem", RC )
          ENDIF

          ! Do carbon chemistry
          CALL Chem_Carbon_Gases( Input_Opt  = Input_Opt,                    &
                                  State_Met  = State_Met,                    &
                                  State_Chm  = State_Chm,                    &
                                  State_Grid = State_Grid,                   &
                                  State_Diag = State_Diag,                   &
                                  RC         = RC                           )

          ! Trap potential errors
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = 'Error encountered in "Chem_Carbon_Gases"!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF

          IF ( Input_Opt%useTimers ) THEN
             CALL Timer_End( "=> Gas-phase chem", RC )
          ENDIF

       !====================================================================
       ! Mercury (configure with -DMECH=Hg)
       !=====================================================================
       ELSE IF ( Input_Opt%ITS_A_MERCURY_SIM ) THEN

          IF ( Input_Opt%useTimers ) THEN
             CALL Timer_Start( "=> Gas-phase chem", RC )
          ENDIF

          ! Do Hg chemistry
          CALL ChemMercury( Input_Opt  = Input_Opt,                          &
                            State_Chm  = State_Chm,                          &
                            State_Diag = State_Diag,                         &
                            State_Grid = State_Grid,                         &
                            State_Met  = State_Met,                          &
                            RC         = RC                                 )

          ! Trap potential errors
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = 'Error encountered in "ChemMercury"!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF

          IF ( Input_Opt%useTimers ) THEN
             CALL Timer_End( "=> Gas-phase chem", RC )
          ENDIF

       !=====================================================================
       ! POPs
       !=====================================================================
       ELSE IF ( Input_Opt%ITS_A_POPS_SIM ) THEN

          ! Do POPS chemistry
          CALL ChemPOPs( Input_Opt  = Input_Opt,                             &
                         State_Chm  = State_Chm,                             &
                         State_Diag = State_Diag,                            &
                         State_Grid = State_Grid,                            &
                         State_Met  = State_Met,                             &
                         RC         = RC                                    )

          ! Trap potential errors
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = 'Error encountered in "ChemPOPs"!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF

       ENDIF

       !### Debug
       IF ( Input_Opt%Verbose ) THEN
          CALL Debug_Msg( '### MAIN: a CHEMISTRY' )
       ENDIF

    ENDIF

    !========================================================================
    ! Convert species units back to original unit (ewl, 8/12/15)
    !========================================================================

    ! Halt "All chemistry" timer (so unitconv+diags can be timed separately)
    IF ( Input_Opt%useTimers ) THEN
       CALL Timer_End( "All chemistry", RC )
    ENDIF

    ! Convert units
    CALL Convert_Spc_Units(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Chm  = State_Chm,                                             &
         State_Grid = State_Grid,                                            &
         State_Met  = State_Met,                                             &
         new_units  = previous_units,                                        &
         RC         = RC                                                    )

    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Unit conversion error'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    ! Start "All chemistry" timer again
    IF ( Input_Opt%useTimers ) THEN
       CALL Timer_Start( "All chemistry", RC )
    ENDIF

    !========================================================================
    ! Chemistry budget diagnostics - Part 2 of 2
    !========================================================================
    IF ( State_Diag%Archive_BudgetChemistry ) THEN

       ! Chemistry timestep [s]
       TS_Chem = Get_Ts_Chem()
       DT_Chem = DBLE( Ts_Chem )

       ! Compute change in column masses (after chemistry - before chemistry)
       ! and store in diagnostic arrays.  Units are [kg/s].
       CALL Compute_Budget_Diagnostics(                                      &
            Input_Opt   = Input_Opt,                                         &
            State_Chm   = State_Chm,                                         &
            State_Diag  = State_Diag,                                        &
            State_Grid  = State_Grid,                                        &
            State_Met   = State_Met,                                         &
            isFull      = State_Diag%Archive_BudgetChemistryFull,            &
            diagFull    = State_Diag%BudgetChemistryFull,                    &
            mapDataFull = State_Diag%Map_BudgetChemistryFull,                &
            isTrop      = State_Diag%Archive_BudgetChemistryTrop,            &
            diagTrop    = State_Diag%BudgetChemistryTrop,                    &
            mapDataTrop = State_Diag%Map_BudgetChemistryTrop,                &
            isPBL       = State_Diag%Archive_BudgetChemistryPBL,             &
            diagPBL     = State_Diag%BudgetChemistryPBL,                     &
            mapDataPBL  = State_Diag%Map_BudgetChemistryPBL,                 &
            isLevs      = State_Diag%Archive_BudgetChemistryLevs,            &
            diagLevs    = State_Diag%BudgetChemistryLevs,                    &
            mapDataLevs = State_Diag%Map_BudgetChemistryLevs,                &
            colMass     = State_Diag%BudgetColumnMass,                       &
            timeStep    = DT_Chem,                                           &
            msg         = 'Compute budget diag (2)'//TRIM(ThisLoc),          &
            RC          = RC                                                )

       ! Trap potential errors
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Chemistry budget diagnostics error 2'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDIF

  END SUBROUTINE Do_Chemistry
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
  SUBROUTINE RECOMPUTE_OD( Input_Opt,  State_Chm, State_Diag,                &
                           State_Grid, State_Met, RC                        )
!
! !USES:
!
    USE AEROSOL_MOD,    ONLY : AEROSOL_CONC
    USE AEROSOL_MOD,    ONLY : RDAER
    USE DUST_MOD,       ONLY : RDUST_ONLINE
    USE ErrCode_Mod
    USE ERROR_MOD,      ONLY : Debug_Msg
    USE Input_Opt_Mod,  ONLY : OptInput
    USE State_Chm_Mod,  ONLY : ChmState
    USE State_Diag_Mod, ONLY : DgnState
    USE State_Grid_Mod, ONLY : GrdState
    USE State_Met_Mod,  ONLY : MetState
    USE TIME_MOD,       ONLY : GET_MONTH
    USE TIME_MOD,       ONLY : GET_YEAR
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
    TYPE(GrdState), INTENT(IN)    :: State_Grid  ! Grid State object
    TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
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
!  See https://github.com/geoschem/geos-chem for complete history
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
    LOGICAL            :: LSSALT, LSULF, LSOA
    INTEGER            :: MONTH, YEAR,   WAVELENGTH

    ! Strings
    CHARACTER(LEN=255) :: ErrMsg, ThisLoc

    !========================================================================
    ! RECOMPUTE_OD begins here!
    !========================================================================

    ! Initialize
    RC       = GC_SUCCESS
    MONTH    = GET_MONTH()
    YEAR     = GET_YEAR()
    ErrMsg   = ''
    ThisLoc  = ' -> at Recompute_OD  (in module GeosCore/chemistry_mod.F90)'

    ! Get month and year

    ! First make sure chemistry is turned on
    IF ( Input_Opt%LCHEM ) THEN

       ! Then make sure that the simulations use aerosol species
       IF (  Input_Opt%ITS_A_FULLCHEM_SIM   .or.                             &
             Input_Opt%ITS_AN_AEROSOL_SIM ) THEN

          ! And then make sure that the aersol species are defined
          IF ( Input_Opt%LSULF .or. Input_Opt%LCARB    .or.                  &
               Input_Opt%LDUST .or. Input_Opt%LSSALT ) THEN

             ! Skip this section if all of these are turned off
             CALL AEROSOL_CONC( Input_Opt,  State_Chm, State_Diag, &
                                State_Grid, State_Met, RC )

             !===============================================================
             ! Call RDAER -- computes aerosol optical depths
             !===============================================================

             ! Calculate the AOD at the wavelength specified in jv_spec_aod
             WAVELENGTH = 1
             CALL RDAER( Input_Opt  = Input_Opt,                             &
                         State_Chm  = State_Chm,                             &
                         State_Diag = State_Diag,                            &
                         State_Grid = State_Grid,                            &
                         State_Met  = State_Met,                             &
                         month      = month,                                 &
                         year       = year,                                  &
                         ODswitch   = waveLength,                            &
                         RC         = RC                                    )

             ! Trap potential errors
             IF ( RC /= GC_SUCCESS ) THEN
                ErrMsg = 'Error encountered in "RdAer"!'
                CALL GC_Error( ErrMsg, RC, ThisLoc )
                RETURN
             ENDIF

             !### Debug
             IF ( Input_Opt%Verbose ) THEN
                CALL Debug_Msg( '### RECOMPUTE_OD: after RDAER' )
             ENDIF

             !===============================================================
             ! If LDUST is turned on, then we have online dust aerosol in
             ! GEOS-CHEM.
             ! If LDUST is turned off, then we don't have online dust
             ! aerosol in GEOS-CHEM...so read monthly-mean dust files
             ! from disk. (rjp, tdf, bmy, 4/1/04)
             !===============================================================
             IF ( Input_Opt%LDUST ) THEN
                CALL Rdust_Online( Input_Opt  = Input_Opt,                   &
                                   State_Chm  = State_Chm,                   &
                                   State_Diag = State_Diag,                  &
                                   State_Grid = State_Grid,                  &
                                   State_Met  = State_Met,                   &
                                   ODswitch   = waveLength,                  &
                                   RC         = RC                          )

                ! Trap potential errors
                IF ( RC /= GC_SUCCESS ) THEN
                   ErrMsg = 'Error encountered in "Rdust_Online"!'
                   CALL GC_Error( ErrMsg, RC, ThisLoc )
                   RETURN
                ENDIF
             ENDIF

             !### Debug
             IF ( Input_Opt%Verbose ) THEN
                CALL DEBUG_MSG( '### RECOMPUTE_OD: after RDUST' )
             ENDIF
          ENDIF
       ENDIF
    ENDIF

  END SUBROUTINE RECOMPUTE_OD
!EOC
END MODULE Chemistry_Mod
