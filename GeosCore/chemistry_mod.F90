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
  PUBLIC  :: Init_Chemistry
  PUBLIC  :: Do_Chemistry
  PUBLIC  :: Recompute_OD
!
! !PRIVATE MEMBER FUNCTIONS:
!
  PRIVATE :: Chem_Passive_Species
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
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
    USE ISORROPIAII_MOD,  ONLY : DO_ISORROPIAII
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
    USE TIME_MOD,         ONLY : GET_ELAPSED_SEC
    USE TIME_MOD,         ONLY : GET_TS_CHEM
    USE UCX_MOD,          ONLY : CALC_STRAT_AER
    USE UnitConv_Mod,     ONLY : Convert_Spc_Units
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
    ! Scalars
    INTEGER            :: N_TROP, N
    INTEGER            :: MONTH
    INTEGER            :: YEAR
    INTEGER            :: WAVELENGTH
    LOGICAL            :: prtDebug
    INTEGER            :: TS_Chem
    REAL(f8)           :: DT_Chem, sDTFC, fDTFC
#ifdef APM
    INTEGER            :: I,J,L
    REAL*8             :: CONCTMPSO4(State_Grid%NX,State_Grid%NY,State_Grid%NZ)
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
    RC        = GC_SUCCESS
    ErrMsg    = ''
    prtDebug  = ( Input_Opt%LPRT .and. Input_Opt%amIRoot )
    ThisLoc   = ' -> at Do_Chemistry  (in module GeosCore/chemistry_mod.F90)'

    ! Save species ID"s on first call
    IF ( FIRST ) THEN
       id_DST1 = Ind_('DST1')
       id_NK1  = Ind_('NK1' )
    ENDIF

    !========================================================================
    ! Chemistry budget diagnostics - Part 1 of 2
    !========================================================================
    IF ( State_Diag%Archive_BudgetChemistry ) THEN

       ! Get initial column masses (full, trop, PBL)
       CALL Compute_Budget_Diagnostics(                                      &
            Input_Opt   = Input_Opt,                                         &
            State_Chm   = State_Chm,                                         &
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
            colMass     = State_Diag%BudgetColumnMass,                       &
            before_op   = .TRUE.,                                            &
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
    CALL Convert_Spc_Units( Input_Opt  = Input_Opt,                          &
                            State_Chm  = State_Chm,                          &
                            State_Grid = State_Grid,                         &
                            State_Met  = State_Met,                          &
                            OutUnit    = 'kg',                               &
                            OrigUnit   = OrigUnit,                           &
                            RC         = RC                                 )

    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Unit conversion error (kg/kg dry -> kg)'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
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
                ! ISORROPIA takes Na+, Cl- into account
                CALL Do_IsorropiaII( Input_Opt  = Input_Opt,                 &
                                     State_Chm  = State_Chm,                 &
                                     State_Diag = State_Diag,                &
                                     State_Grid = State_Grid,                &
                                     State_Met  = State_Met,                 &
                                     RC         = RC                        )

                ! Trap potential errors
                IF ( RC /= GC_SUCCESS ) THEN
                   ErrMsg = 'Error encountered in "Do_ISORROPIAII"!'
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
             CALL Timer_Start( "=> FlexChem", RC )
          ENDIF

          ! Solve the KPP-generated mechanism
          CALL Do_FullChem( Input_Opt  = Input_Opt,                          &
                            State_Chm  = State_Chm,                          &
                            State_Diag = State_Diag,                         &
                            State_Grid = State_Grid,                         &
                            State_Met  = State_Met,                          &
                            RC         = RC                                 )

          ! Check units (ewl, 10/5/15)
          IF ( TRIM( State_Chm%Spc_Units ) /= 'kg' ) THEN
             ErrMsg = 'Incorrect species units after Do_FullChem!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF

          ! Trap potential errors
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = 'Error encountered in "Do_FlexChem"!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF

          IF ( Input_Opt%useTimers ) THEN
             CALL Timer_End( "=> FlexChem", RC )
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

             ! Check units (ewl, 10/5/15)
             IF ( TRIM( State_Chm%Spc_Units ) /= 'kg' ) THEN
                ErrMsg = 'Incorrect species units after DO_LINEARCHEM!'
                CALL GC_Error( ErrMsg, RC, ThisLoc )
             ENDIF

             ! Trap potential errors
             IF ( RC /= GC_SUCCESS ) THEN
                ErrMsg = 'Error encountered in ""!'
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
          IF ( id_NK1 > 0 ) THEN
             CALL Do_TOMAS( Input_Opt  = Input_Opt,                          &
                            State_Chm  = State_Chm,                          &
                            State_Diag = State_Diag,                         &
                            State_Grid = State_Grid,                         &
                            State_Met  = State_Met,                          &
                            RC         = RC                                 )


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

          ! Check units (ewl, 10/5/15)
          IF ( TRIM( State_Chm%Spc_Units ) /= 'kg' ) THEN
             ErrMsg = 'Incorrect species units after AEROSOL_CONC'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
          ENDIF

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
                ! ISORROPIA takes Na+, Cl- into account
                CALL Do_IsorropiaII( Input_Opt  = Input_Opt,                 &
                                     State_Chm  = State_Chm,                 &
                                     State_Diag = State_Diag,                &
                                     State_Grid = State_Grid,                &
                                     State_Met  = State_Met,                 &
                                     RC         = RC                        )
#endif

                ! Trap potential errors
                IF ( RC /= GC_SUCCESS ) THEN
                   ErrMsg = 'Error encountered in "Do_IsorropiaII"!'
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
       ! Rn-Pb-Be
       !=====================================================================
       ELSE IF ( Input_Opt%ITS_A_RnPbBe_SIM ) THEN

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

       !====================================================================
       ! Mercury (configure with -DMECH=Hg)
       !=====================================================================
       ELSE IF ( Input_Opt%ITS_A_MERCURY_SIM ) THEN

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

       !=====================================================================
       ! POPs (only used when compiled with BPCH_DIAG=y)
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

       !=====================================================================
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
       !=====================================================================
       IF ( Input_Opt%NPassive_Decay > 0 ) THEN

          ! Apply loss rate to passive species with finite lifetimes
          CALL Chem_Passive_Species( Input_Opt  = Input_Opt,                 &
                                     State_Chm  = State_Chm,                 &
                                     State_Grid = State_Grid,                &
                                     State_Met  = State_Met,                 &
                                     RC         = RC                        )

          ! Trap potential errors
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = 'Error encountered in "Chem_Passive_Species"!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF

          !### Debug
          IF ( prtDebug ) THEN
             CALL Debug_Msg( '### MAIN: a CHEMISTRY' )
          ENDIF

       ENDIF

    ENDIF

    !========================================================================
    ! Convert species units back to original unit (ewl, 8/12/15)
    !========================================================================
    CALL Convert_Spc_Units( Input_Opt  = Input_Opt,                          &
                            State_Chm  = State_Chm,                          &
                            State_Grid = State_Grid,                         &
                            State_Met  = State_Met,                          &
                            OutUnit    = OrigUnit,                           &
                            RC         = RC                                 )

    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Unit conversion error'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
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
            colMass     = State_Diag%BudgetColumnMass,                       &
            timeStep    = DT_Chem,                                           &
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
    LOGICAL            :: prtDebug
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
    prtDebug = ( Input_Opt%LPRT .and. Input_Opt%amIRoot )
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
             IF ( prtDebug ) THEN
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
             IF ( prtDebug ) THEN
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
! !IROUTINE: chem_passive_species
!
! !DESCRIPTION: Subroutine CHEM\_PASSIVE\_SPECIES performs loss chemistry
!  on passive species with finite atmospheric lifetimes.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Chem_Passive_Species( Input_Opt, State_Chm, State_Grid,         &
                                   State_Met, RC                            )
!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod,  ONLY : OptInput
    USE State_Chm_Mod,  ONLY : ChmState
    USE State_Chm_Mod,  ONLY : ind_
    USE State_Grid_Mod, ONLY : GrdState
    USE State_Met_Mod,  ONLY : MetState
    USE Time_Mod,       ONLY : Get_Ts_Chem
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN   ) :: Input_Opt   ! Input options object
    TYPE(GrdState), INTENT(IN   ) :: State_Grid  ! Grid state object
    TYPE(MetState), INTENT(IN   ) :: State_Met   ! Meteorology state object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry state object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Failure or success
!
! !REMARKS:
!
! !REVISION HISTORY:
!  04 Sep 2015 - C. Keller   - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    LOGICAL             :: prtDebug
    LOGICAL             :: Is_HalfLife
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

    !========================================================================
    ! Chem_Passive_Species begins here!
    !========================================================================

    ! Initialize
    RC       = GC_SUCCESS
    prtDebug = ( Input_Opt%LPRT .and. Input_Opt%amIRoot )
    ErrMsg   = ''
    ThisLoc  = &
       ' -> at Chem_Passive_Species (in module GeosCore/chemistry_mod.F90)'

    DT       = GET_TS_CHEM() ! timestep in seconds

    ! For now, always compute decay using e-folding time
    Is_HalfLife = .FALSE.

    !========================================================================
    ! Apply decay loss rate only to those passive species that have a
    ! finite atmospheric lifetime (this speeds up execution)
    !========================================================================

    ! Loop over all decaying passive species
    DO N = 1, Input_Opt%NPassive_Decay

       !---------------------------------------------------------------------
       ! Find the GEOS-Chem species Id
       !---------------------------------------------------------------------

       ! Get the Id of the species in the passive decay menu
       Id   = Input_Opt%Passive_DecayID(N)

       ! Convert this to a GEOS-Chem species Id number
       GcId = Ind_( TRIM( Input_Opt%PASSIVE_NAME(Id) ) )

       ! Make sure the model ID is valid
       IF ( GcId < 0 ) THEN
          ErrMsg = 'Could not find the GEOS-Chem species ID # '           // &
                   'for passive species : '                               // &
                   TRIM( Input_Opt%PASSIVE_NAME(Id) )
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       !---------------------------------------------------------------------
       ! Compute the decay rate
       !---------------------------------------------------------------------

       ! Compute the decay rate for each passive species
       IF ( Is_HalfLife ) THEN
          Decay = ln2 / Input_Opt%PASSIVE_TAU(Id)
       ELSE
          Decay = 1.0 / Input_Opt%PASSIVE_TAU(Id)
       ENDIF
       Rate  = EXP( - DT * Decay )

       !### Debug output
       IF ( First ) THEN
          IF ( prtDebug ) THEN
             WRITE( 6,100 ) ADJUSTL( Input_Opt%PASSIVE_NAME(Id) ),           &
                            GcId, Rate
 100         FORMAT( '     -  Pass. species name, Id, loss rate: ',          &
                      a15, i5, 1x, es13.6 )
          ENDIF
       ENDIF

       !---------------------------------------------------------------------
       ! Apply loss
       !---------------------------------------------------------------------

       !$OMP PARALLEL DO                                                     &
       !$OMP DEFAULT( SHARED                                                )&
       !$OMP PRIVATE( I, J, L                                               )&
       !$OMP COLLAPSE( 3                                                    )
       DO L = 1, State_Grid%NZ
       DO J = 1, State_Grid%NY
       DO I = 1, State_Grid%NX
          State_Chm%Species(GcId)%Conc(I,J,L) =                &
                    State_Chm%Species(GcId)%Conc(I,J,L) * Rate
       ENDDO
       ENDDO
       ENDDO
       !$OMP END PARALLEL DO

    ENDDO

    ! Reset after the first time
    IF ( First ) First = .FALSE.

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
  SUBROUTINE Init_Chemistry( Input_Opt, State_Chm, State_Diag, State_Grid, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE FAST_JX_MOD,       ONLY : Init_FJX
    USE FullChem_Mod,      ONLY : Init_FullChem
    USE Input_Opt_Mod,     ONLY : OptInput
    USE State_Chm_Mod,     ONLY : ChmState
    USE State_Chm_Mod,     ONLY : Ind_
    USE State_Diag_Mod,    ONLY : DgnState
    USE State_Grid_Mod,    ONLY : GrdState
!
! !INPUT PARAMETERS:
!
    TYPE(GrdState), INTENT(IN)    :: State_Grid  ! Grid State object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(INOUT) :: Input_Opt   ! Input Options object
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
    TYPE(DgnState), INTENT(INOUT) :: State_Diag  ! Diagnostics State object
    INTEGER,        INTENT(INOUT) :: RC          ! Success or failure?
!
! !REMARKS:
!  We initialize relevant fullchem and carbon KPP mechanism variables
!  here in order to use values from the Species Database.  When the other
!  modules are initialized (most of which are done in GC_Init_Extra), at
!  that point the Species Database has not been read from the YAML file,
!  so we must call Init_FullChem and Init_Carbon_Gases here.
!
! !REVISION HISTORY:
!  19 May 2014 - C. Keller   - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
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
    ThisLoc  = ' -> at Init_Chemistry  (in module GeosCore/chemistry_mod.F90)'

    ! Skip if we have already done this
    IF ( FIRST ) THEN

       ! Adjust first flag
       FIRST  = .FALSE.

       !--------------------------------------------------------------------
       ! Initialize Fast-JX photolysis
       ! (except for carbon, which has no photolysis)
       !
       ! NOTE: we need to call this for a dry-run so that we can get
       ! a list of all of the lookup tables etc. that FAST-JX reads
       !--------------------------------------------------------------------
       IF ( .not. Input_Opt%ITS_A_CARBON_SIM ) THEN
          CALL Init_FJX( Input_Opt, State_Chm, State_Diag, State_Grid, RC )
       ENDIF

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
