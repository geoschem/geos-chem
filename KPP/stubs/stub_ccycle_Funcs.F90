MODULE ccycle_Funcs
  !
  ! Stub module for KPP/ccycle/ccycle_Funcs.F90
  !
  USE gckpp_Precision
  USE gckpp_Parameters
  USE gckpp_Global
  USE Precision_Mod,   ONLY : fp
  USE rateLawUtilFuncs
  !
CONTAINS
  !
  SUBROUTINE ccycle_ConvertKgToMolecCm3( I,          J,          L,          &
                                         id_CH4,     id_CO,      id_CO2,     &
                                         xnumol_CH4, xnumol_CO2, xnumol_CO,  &
                                         State_Chm,  State_Met              )
    !
    ! Stub routine for ccycle_ConvertKgToMolecCm3,
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
  END SUBROUTINE ccycle_ConvertKgToMolecCm3
  !
  SUBROUTINE ccycle_ComputeRateConstants(                                    &
             I,        J,        L,           bAirDens,    bCl,              &
             bOH,      CH4loss,  GMI_Prod_CO, GMI_Loss_CO, PCO_nmVOC,        &
             PCO_CH4,  LPCO_CH4, dtChem,      tCosZ,       State_Chm,        &
             State_Met                                                      )
    !
    ! Stub for routine ccycle_ComputeRateConstants,
    ! needed to satisfy compile-time dependencies
    !
    USE State_Chm_Mod, ONLY : ChmState
    USE State_Met_Mod, ONLY : MetState
    !
    INTEGER,        INTENT(IN)    :: I, J, L
    REAL(fp),       INTENT(IN)    :: bAirDens
    REAL(fp),       INTENT(IN)    :: bCl
    REAL(fp),       INTENT(IN)    :: bOH
    REAL(fp),       INTENT(IN)    :: CH4loss
    REAL(fp),       INTENT(IN)    :: GMI_Prod_CO
    REAL(fp),       INTENT(IN)    :: GMI_Loss_CO
    REAL(fp),       INTENT(IN)    :: PCO_NMVOC
    LOGICAL,        INTENT(IN)    :: LPCO_CH4
    REAL(fp),       INTENT(IN)    :: PCO_CH4
    REAL(fp),       INTENT(IN)    :: dtChem
    REAL(fp),       INTENT(IN)    :: tCosZ
    TYPE(ChmState), INTENT(IN)    :: State_Chm
    TYPE(MetState), INTENT(IN)    :: State_Met
  END SUBROUTINE ccycle_ComputeRateConstants
  !
  SUBROUTINE ccycle_ConvertMolecCm3ToKg( I,          J,         L,           &
                                         id_CH4,     id_CO,     id_COch4,    &
                                         id_COnmvoc, id_CO2,    xnumol_CH4,  &
                                         xnumol_CO2, xnumol_CO, State_Chm,   &
                                         State_Met                          )
    !
    ! Stub for ccycle_ConvertMolecCm3ToKg,
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
  END SUBROUTINE ccycle_ConvertMolecCm3ToKg

END MODULE ccycle_Funcs
