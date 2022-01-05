MODULE fullchem_SulfurChemFuncs
  !
  ! Stub module for GeosChem/fullchem_mod.F90
  !
  IMPLICIT NONE
  PRIVATE
  !
  PUBLIC :: fullchem_SulfurAqChem
  PUBLIC :: fullchem_SulfurCldChem
  PUBLIC :: fullchem_InitSulfurChem
  !
CONTAINS
  !
  SUBROUTINE fullchem_SulfurAqChem( I,         J,          L,                &
                                    Input_Opt, State_Chm,  State_Grid,       &
                                    State_Met, RC                           )
    !
    ! Stub routine for GeosChem/fullchem_mod.F90
    !
    USE Input_Opt_Mod,    ONLY : OptInput
    USE State_Chm_Mod,    ONLY : ChmState
    USE State_Met_Mod,    ONLY : MetState
    USE State_Grid_Mod,   ONLY : GrdState
    !
    INTEGER,        INTENT(IN)    :: I, J, L    ! Lon, lat, level indices
    TYPE(MetState), INTENT(IN)    :: State_Met  ! Meteorology State object
    TYPE(GrdState), INTENT(IN)    :: State_Grid  ! Grid State object
    TYPE(OptInput), INTENT(IN)    :: Input_Opt  ! Input Options object
    TYPE(ChmState), INTENT(INOUT) :: State_Chm  ! Chemistry State object
    INTEGER,        INTENT(OUT)   :: RC         ! Success or failure
    !
    RC = 0
  END SUBROUTINE fullchem_SulfurAqChem
  !
  SUBROUTINE fullchem_SulfurCldChem( I,         J,         L,                &
                                    Input_Opt,  State_Chm, State_Diag,       &
                                    State_Grid, State_Met, RC               )
    !
    ! Stub routine for GeosChem/fullchem_mod.F90
    !
    USE Input_Opt_Mod,      ONLY : OptInput
    USE State_Chm_Mod,      ONLY : ChmState
    USE State_Diag_Mod,     ONLY : DgnState
    USE State_Grid_Mod,     ONLY : GrdState
    USE State_Met_Mod,      ONLY : MetState
    !
    INTEGER,        INTENT(IN)    :: I, J, L
    TYPE(OptInput), INTENT(IN)    :: Input_Opt
    TYPE(GrdState), INTENT(IN)    :: State_Grid
    TYPE(MetState), INTENT(IN)    :: State_Met
    TYPE(ChmState), INTENT(INOUT) :: State_Chm
    TYPE(DgnState), INTENT(INOUT) :: State_Diag
    INTEGER,        INTENT(OUT)   :: RC
    !
    RC = 0
  END SUBROUTINE fullchem_SulfurCldChem
  !
  SUBROUTINE fullchem_InitSulfurChem( RC )
    !
    ! Stub routine for GeosChem/fullchem_mod.F90
    !    
    INTEGER, INTENT(OUT) :: RC
    !
    RC = 0
  END SUBROUTINE fullchem_InitSulfurChem
  !
END MODULE fullchem_SulfurChemFuncs
