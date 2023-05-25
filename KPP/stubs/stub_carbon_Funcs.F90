MODULE carbon_Funcs
  !
  ! Stub module for KPP/carbon/carbon_Funcs.F90
  !
  USE gckpp_Precision
  USE gckpp_Parameters
  USE gckpp_Global
  USE Precision_Mod,   ONLY : fp
  USE rateLawUtilFuncs
  !
CONTAINS
  !
  SUBROUTINE carbon_ConvertKgToMolecCm3(                                &
             I,        J,          L,          id_CH4,    id_CO,             &
             id_CO2,   xnumol_CH4, xnumol_CO2, xnumol_CO, State_Chm,         &
             State_Met                                                      )
    !
    ! Stub routine for carbon_ConvertKgToMolecCm3,
    ! needed to satisfy compile-time dependencies
    !
    USE State_Chm_Mod, ONLY : ChmState
    USE State_Met_Mod, ONLY : MetState
    !
    INTEGER,        INTENT(IN) :: I, J, L
    INTEGER,        INTENT(IN) :: id_CH4
    INTEGER,        INTENT(IN) :: id_CO
    INTEGER,        INTENT(IN) :: id_CO2
    REAL(fp),       INTENT(IN) :: xnumol_CH4
    REAL(fp),       INTENT(IN) :: xnumol_CO
    REAL(fp),       INTENT(IN) :: xnumol_CO2
    TYPE(MetState), INTENT(IN) :: State_Met
    TYPE(ChmState), INTENT(IN) :: State_Chm
  END SUBROUTINE carbon_ConvertKgToMolecCm3
  !
  SUBROUTINE carbon_ComputeRateConstants(                               &
             I,             J,                 L,                            &
             ConcClMnd,     ConcOHMnd,         LCH4_by_OH,                   &
             LCO_in_Strat,  OHdiurnalFac,      PCO_fr_CH4_use,               &
             PCO_fr_CH4,    PCO_fr_NMVOC_use,  PCO_fr_NMVOC,                 &
             PCO_in_Strat,  dtChem,            State_Chm,                    &
             State_Met                                                      )
    !
    ! Stub for routine carbon_ComputeRateConstants,
    ! needed to satisfy compile-time dependencies
    !
    USE State_Chm_Mod, ONLY : ChmState
    USE State_Met_Mod, ONLY : MetState
    !
    INTEGER,        INTENT(IN) :: I, J, L
    REAL(fp),       INTENT(IN) :: ConcClMnd
    REAL(fp),       INTENT(IN) :: ConcOHmnd
    REAL(fp),       INTENT(IN) :: LCH4_by_OH
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
  !
  SUBROUTINE carbon_ConvertMolecCm3ToKg(                                &
             I,         J,          L,         id_CH4,     id_CO,            &
             id_COch4,  id_COnmvoc, id_CO2,    xnumol_CH4, xnumol_CO2,       &
             xnumol_CO, State_Chm,  State_Met                               )
    !
    ! Stub for carbon_ConvertMolecCm3ToKg,
    ! needed to satisfy compile-time dependencies
    !
    USE State_Chm_Mod, ONLY : ChmState
    USE State_Met_Mod, ONLY : MetState
    !
    INTEGER,        INTENT(IN)    :: id_CH4
    INTEGER,        INTENT(IN)    :: id_CO
    INTEGER,        INTENT(IN)    :: id_COch4
    INTEGER,        INTENT(IN)    :: id_COnmvoc
    INTEGER,        INTENT(IN)    :: id_CO2
    REAL(fp),       INTENT(IN)    :: xnumol_CH4
    REAL(fp),       INTENT(IN)    :: xnumol_CO
    REAL(fp),       INTENT(IN)    :: xnumol_CO2
    TYPE(MetState), INTENT(IN)    :: State_Met
    TYPE(ChmState), INTENT(INOUT) :: State_Chm
  END SUBROUTINE carbon_ConvertMolecCm3ToKg

  FUNCTION carbon_Get_COfromCH4_Flux( dtChem ) RESULT ( flux )
    !
    ! Stub for carbon_Get_CO_CH4_Flux
    ! needed to satisfy compile-time dependencies
    !
    REAL(dp), INTENT(IN) :: dtChem
    REAL(dp)             :: flux
  END FUNCTION carbon_Get_COfromCH4_Flux

  FUNCTION carbon_Get_COfromNMVOC_Flux( dtChem ) RESULT ( flux )
    !
    ! Stub for carbon_Get_CO_NMVOC_Flux
    ! needed to satisfy compile-time dependencies
    !
    REAL(dp), INTENT(IN) :: dtChem
    REAL(dp)             :: flux
  END FUNCTION carbon_Get_COfromNMVOC_Flux

  FUNCTION  carbon_Get_CO2fromOH_Flux( dtChem ) RESULT ( flux )
    !
    ! Stub for carbon_Get_CO_NMVOC_Flux
    ! needed to satisfy compile-time dependencies
    !
    REAL(dp), INTENT(IN) :: dtChem
    REAL(dp)             :: flux
  END FUNCTION carbon_Get_CO2fromOH_Flux

  FUNCTION carbon_Get_FixedOH_Flux( dtChem ) RESULT ( flux )
    !
    ! Stub for ccarbon_Get_OH_E_Flux
    ! needed to satisfy compile-time dependencies
    !
    REAL(dp), INTENT(IN) :: dtChem
    REAL(dp)             :: flux
  END FUNCTION carbon_Get_FixedOH_Flux

END MODULE carbon_Funcs
