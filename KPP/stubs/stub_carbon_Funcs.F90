MODULE carbon_Funcs
  ! Stub module for KPP/carbon/carbon_Funcs.F90 needed to satisfy
  ! compile-time dependencies for non-carbon chemistry mechanisms.

  USE gckpp_Precision
  USE gckpp_Parameters
  USE gckpp_Global
  USE Precision_Mod,   ONLY : fp
  USE rateLawUtilFuncs

CONTAINS

  !------------------------------------------------------------
  SUBROUTINE carbon_ConvertKgToMolecCm3( I, J, L, State_Chm, State_Met )
    USE State_Chm_Mod, ONLY : ChmState
    USE State_Met_Mod, ONLY : MetState
    INTEGER,        INTENT(IN) :: I, J, L
    TYPE(MetState), INTENT(IN) :: State_Met
    TYPE(ChmState), INTENT(IN) :: State_Chm
  END SUBROUTINE carbon_ConvertKgToMolecCm3

  !------------------------------------------------------------
  SUBROUTINE carbon_ComputeRateConstants(                       &
             I,             J,                 L,               &
             ConcClMnd,     ConcOHMnd,         LCH4_in_Strat,   &
             LCO_in_Strat,  OHdiurnalFac,      PCO_fr_CH4_use,  &
             PCO_fr_CH4,    PCO_fr_NMVOC_use,  PCO_fr_NMVOC,    &
             PCO_in_Strat,  dtChem,            State_Chm,       &
             State_Met                                         )
    USE State_Chm_Mod, ONLY : ChmState
    USE State_Met_Mod, ONLY : MetState
    INTEGER,        INTENT(IN) :: I, J, L
    REAL(fp),       INTENT(IN) :: ConcClMnd
    REAL(fp),       INTENT(IN) :: ConcOHmnd
    REAL(fp),       INTENT(IN) :: LCH4_in_Strat
    REAL(fp),       INTENT(IN) :: LCO_in_Strat
    REAL(fp),       INTENT(IN) :: OHdiurnalFac
    LOGICAL,        INTENT(IN) :: PCO_fr_CH4_use
    REAL(fp),       INTENT(IN) :: PCO_fr_CH4
    LOGICAL,        INTENT(IN) :: PCO_fr_NMVOC_use
    REAL(fp),       INTENT(IN) :: PCO_fr_NMVOC
    REAL(fp),       INTENT(IN) :: PCO_in_Strat
    REAL(fp),       INTENT(IN) :: dtChem
    TYPE(ChmState), INTENT(IN) :: State_Chm
    TYPE(MetState), INTENT(IN) :: State_Met
  END SUBROUTINE carbon_ComputeRateConstants
  
  !------------------------------------------------------------
  SUBROUTINE carbon_ConvertMolecCm3ToKg( I, J, L, State_Chm, State_Met )
    USE State_Chm_Mod, ONLY : ChmState
    USE State_Met_Mod, ONLY : MetState
    TYPE(MetState), INTENT(IN)    :: State_Met
    TYPE(ChmState), INTENT(INOUT) :: State_Chm
  END SUBROUTINE carbon_ConvertMolecCm3ToKg
  
  !------------------------------------------------------------
  SUBROUTINE carbon_InitCarbonKPPFuncs( kgmolec_CH4, kgmolec_CO, kgmolec_CO2, RC )
    USE ErrCode_Mod
    USE State_Chm_Mod,  ONLY : Ind_
    REAL(fp), INTENT(IN)  :: kgmolec_CH4
    REAL(fp), INTENT(IN)  :: kgmolec_CO
    REAL(fp), INTENT(IN)  :: kgmolec_CO2
    INTEGER,  INTENT(OUT) :: RC
  END SUBROUTINE carbon_InitCarbonKPPFuncs

  !------------------------------------------------------------
  SUBROUTINE carbon_CleanupCarbonKPPFuncs( RC )
    USE ErrCode_Mod
    INTEGER,  INTENT(OUT) :: RC
  END SUBROUTINE carbon_CleanupCarbonKPPFuncs

  !------------------------------------------------------------
  FUNCTION carbon_Get_COfromCH4_Flux( dtChem ) RESULT ( flux )
    REAL(dp), INTENT(IN) :: dtChem
    REAL(dp)             :: flux
  END FUNCTION carbon_Get_COfromCH4_Flux

  !------------------------------------------------------------
  FUNCTION carbon_Get_COfromNMVOC_Flux( dtChem ) RESULT ( flux )
    REAL(dp), INTENT(IN) :: dtChem
    REAL(dp)             :: flux
  END FUNCTION carbon_Get_COfromNMVOC_Flux

END MODULE carbon_Funcs
