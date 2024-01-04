!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: emissions_mod.F90
!
! !DESCRIPTION: Module emissions\_mod.F90 is a wrapper module to interface
! GEOS-Chem and HEMCO. It basically just calls the GEOS-Chem - HEMCO interface
! routines. For some specialty sims, a few additional steps are required that
! are also executed here.
!\\
!\\
! !INTERFACE:
!
MODULE Emissions_Mod
!
! !USES:
!
  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: Emissions_Init
  PUBLIC :: Emissions_Run
  PUBLIC :: Emissions_Final
!
! !REVISION HISTORY:
!  27 Aug 2014 - C. Keller   - Initial version.
!  See https://github.com/geoschem/geos-chem for complete history
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
! !IROUTINE: emissions_init
!
! !DESCRIPTION: Subroutine EMISSIONS\_INIT calls the HEMCO - GEOS-Chem
! interface initialization routines.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Emissions_Init( Input_Opt, State_Chm, State_Grid, State_Met,  &
                             RC,        HcoConfig )
!
! !USES:
!
    USE ErrCode_Mod
    USE HCO_Interface_GC_Mod, ONLY : HCoi_GC_Init
    USE HCO_Types_Mod,        ONLY : ConfigObj
    USE Input_Opt_Mod,        ONLY : OptInput
    USE State_Chm_Mod,        ONLY : ChmState
    USE State_Chm_Mod,        ONLY : Ind_
    USE State_Grid_Mod,       ONLY : GrdState
    USE State_Met_Mod,        ONLY : MetState
!
! !INPUT PARAMETERS:
!
    TYPE(ChmState),  INTENT(IN   )          :: State_Chm  ! Chemistry state
    TYPE(GrdState),  INTENT(IN   )          :: State_Grid ! Grid state
    TYPE(MetState),  INTENT(IN   )          :: State_Met  ! Met state
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(OptInput),  INTENT(INOUT)          :: Input_Opt  ! Input opts
    INTEGER,         INTENT(INOUT)          :: RC         ! Failure or success
    TYPE(ConfigObj), POINTER,      OPTIONAL :: HcoConfig  ! HEMCO config object
!
! !REVISION HISTORY:
!  27 Aug 2014 - C. Keller   - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC

    ! Strings
    CHARACTER(LEN=255) :: ErrMsg, ThisLoc

    !=================================================================
    ! EMISSIONS_INIT begins here!
    !=================================================================

    ! Initialize
    RC       = GC_SUCCESS
    ErrMsg   = ''
    ThisLoc  = ' -> at Emissions_Init (in module GeosCore/emissions_mod.F90)'

    ! Initialize the HEMCO environment for this GEOS-Chem run.
    CALL HCOI_GC_Init( Input_Opt, State_Chm, State_Grid, State_Met, &
                       RC,        HcoConfig=HcoConfig )

    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered in "HCOI_GC_Init"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

  END SUBROUTINE Emissions_Init
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: emissions_run
!
! !DESCRIPTION: Subroutine EMISSIONS\_RUN calls the HEMCO - GEOS-Chem
! interface run routines.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Emissions_Run( Input_Opt, State_Chm, State_Diag, State_Grid,    &
                            State_Met, EmisTime,  Phase,      RC            )
!
! !USES:
!
    USE CARBON_MOD,            ONLY : EmissCarbon
    USE Carbon_Gases_Mod,      ONLY : CO2_Production
    USE CO2_MOD,               ONLY : EmissCO2
    USE ErrCode_Mod
    USE HCO_Interface_GC_Mod,  ONLY : HCOI_GC_Run
    USE Input_Opt_Mod,         ONLY : OptInput
    USE Mercury_Mod,           ONLY : EmissMercury
    USE Precision_Mod
    USE State_Chm_Mod,         ONLY : ChmState
    USE State_Diag_Mod,        ONLY : DgnState
    USE State_Grid_Mod,        ONLY : GrdState
    USE State_Met_Mod,         ONLY : MetState
    USE Time_Mod,              ONLY : Get_Ts_Emis
    USE Tracer_Mod,            ONLY : Tracer_Source_Phase
    Use SfcVmr_Mod,            ONLY : FixSfcVmr_Run
#ifdef TOMAS
    USE CARBON_MOD,            ONLY : EmissCarbonTomas
    USE SULFATE_MOD,           ONLY : EmissSulfateTomas
#endif
!
! !INPUT PARAMETERS:
!
    LOGICAL,        INTENT(IN   )  :: EmisTime   ! Emissions in this time step
    INTEGER,        INTENT(IN   )  :: Phase      ! Run phase
    TYPE(GrdState), INTENT(IN   )  :: State_Grid ! Grid State object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT)  :: State_Chm  ! Chemistry State object
    TYPE(DgnState), INTENT(INOUT)  :: State_Diag ! Diagnostics State object
    TYPE(MetState), INTENT(INOUT)  :: State_Met  ! Meteorology State object
    TYPE(OptInput), INTENT(INOUT)  :: Input_Opt  ! Input Options object
    INTEGER,        INTENT(INOUT)  :: RC         ! Failure or success
!
! !REVISION HISTORY:
!  27 Aug 2014 - C. Keller   - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! LOCAL VARIABLES:
!
    ! Strings
    CHARACTER(LEN=255) :: ErrMsg, ThisLoc

    !=================================================================
    ! EMISSIONS_RUN begins here!
    !=================================================================

    ! Initialize
    RC       = GC_SUCCESS
    ErrMsg   = ''
    ThisLoc  = ' -> at Emissions_Run (in module GeosCore/emissions_mod.F90)'

    ! Run HEMCO. Phase 1 will only update the HEMCO clock and the
    ! HEMCO data list, phase 2 will perform the emission calculations.
    CALL HCOI_GC_Run( Input_Opt, State_Chm,  State_Grid, &
                      State_Met, State_Diag, EmisTime, Phase, RC )
    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered in "HCOI_GC_Run"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    ! Exit if Phase 0 or 1, or if it is a GEOS-Chem dry-run
    IF ( ( Phase == 1 ) .or. ( Phase == 0 ) .or. Input_Opt%DryRun ) RETURN

    ! Call carbon emissions module to make sure that sesquiterpene
    ! emissions calculated in HEMCO (SESQ) are passed to the internal
    ! species array in carbon, as well as to ensure that POA emissions
    ! are correctly treated.
    IF ( Input_Opt%ITS_A_FULLCHEM_SIM   .or. &
         Input_Opt%ITS_AN_AEROSOL_SIM ) THEN
       CALL EmissCarbon( Input_Opt, State_Chm, State_Grid, State_Met, RC )

       ! Trap potential errors
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Error encountered in "EmissCarbon"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDIF

#ifdef TOMAS
    ! Call TOMAS emission routines (JKodros 6/2/15)
    CALL EmissCarbonTomas( Input_Opt, State_Chm, State_Grid, State_Met, &
               State_Diag,  RC )

    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered in "EmissCarbonTomas"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    CALL EmissSulfateTomas( Input_Opt, State_Chm, State_Grid, State_Met, RC )

    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered in "EmissSulfateTomas"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
#endif

    ! For the CO2 simulation, we manually add the chemical production of CO2
    ! from CO oxidation (which is listed as a non-chemical source in HEMCO)
    ! to State_Chm%Species.  This is done in EmissCO2.  All other CO2 emissions
    ! (as of GEOS-Chem 12.0.2) are now added via HEMCO, and diagnostics for
    ! these quantities are saved out via HEMCO diagnostics. (bmy, 10/18/18)
    IF ( Input_Opt%ITS_A_CO2_SIM ) THEN
       CALL EmissCO2( Input_Opt, State_Chm, State_Diag, State_Grid, &
                      State_Met, RC )

       ! Trap potential errors
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Error encountered in "EmissCO2"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDIF

    ! For transport tracer simulation
    IF ( Input_Opt%ITS_A_TRACER_SIM ) THEN
       CALL Tracer_Source_Phase( Input_Opt, State_Chm, State_Grid, &
                                 State_Met, RC )

       ! Trap potential errors
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Error encountered in "Tracer_Source_Phase"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDIF

    ! Carbon simulation (e.g. CO2-CO-CH4-OCS via KPP)
    !
    ! Computes CO2 production from CO oxidation
    IF ( Input_Opt%ITS_A_CARBON_SIM ) THEN
       CALL CO2_Production( Input_Opt,  State_Chm, State_Diag,            &
                            State_Grid, State_Met, RC                    )

       ! Trap potential errors
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Error encountered in "Emiss_Carbon_Gases"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDIF

    ! For mercury, use old emissions code for now
    IF ( Input_Opt%ITS_A_MERCURY_SIM ) THEN
       CALL EmissMercury( Input_Opt,  State_Chm, State_Diag,                 &
                          State_Grid, State_Met, RC                         )

       ! Trap potential errors
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Error encountered in "EmissMercury"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDIF

    ! Prescribe some concentrations if needed
    IF ( Input_Opt%ITS_A_FULLCHEM_SIM ) THEN
       ! Set other (non-UCX) fixed VMRs
       If ( Input_Opt%DoEmissions ) Then
          CALL FixSfcVMR_Run( Input_Opt, State_Chm, State_Grid, State_Met, RC )

          ! Trap potential errors
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = 'Error encountered in "FixSfcVmr_Run"!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF

       ENDIF

    ENDIF

   END SUBROUTINE EMISSIONS_RUN
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: emissions_final
!
! !DESCRIPTION: Subroutine EMISSIONS\_FINAL calls the HEMCO - GEOS-Chem
! interface finalization routines.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Emissions_Final( Error, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE HCO_Interface_GC_Mod, ONLY : HCOI_GC_Final
!
! !INPUT PARAMETERS:
!
    LOGICAL, INTENT(IN)  :: Error      ! Cleanup arrays after crash?
!
! !OUTPUT PARAMETERS:
!
    INTEGER, INTENT(OUT) :: RC         ! Success or failure?
!
! !REVISION HISTORY:
!  27 Aug 2014 - C. Keller    - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC

    ! Strings
    CHARACTER(LEN=255) :: ErrMsg, ThisLoc

    !=================================================================
    ! EMISSIONS_FINAL begins here!
    !=================================================================

    ! Initialize
    RC      = GC_SUCCESS
    ErrMsg  = ''
    ThisLoc = ' -> at HCOI_GC_Final (in module GeosCore/hcoi_gc_final_mod.F90)'

    ! Finalize HEMCO
    CALL HCOI_GC_Final( Error, RC )

    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered in "HCOI_GC_Final"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
  END SUBROUTINE Emissions_Final
!EOC
END MODULE Emissions_Mod
