MODULE fullchem_SulfurChemFuncs
  !
  ! Stub module for GeosCore/fullchem_mod.F90
  ! This allows us to satisfy compilation dependencies for the "fullchem"
  ! mechanism when compiling other mechanisms (e.g. Hg)
  !
  IMPLICIT NONE
  PUBLIC
  !
CONTAINS
  !
  SUBROUTINE fullchem_ConvertAlkToEquiv()
    !
    ! Stub for GeosCore/fullchem_mod.F90:fullchem_ConvertAlkToEquiv
    !
  END SUBROUTINE fullchem_ConvertAlkToEquiv
  !
  SUBROUTINE fullchem_ConvertEquivToAlk()
    !
    ! Stub for GeosCore/fullchem_mod.F90:fullchem_ConvertEquivToAlk
    !
  END SUBROUTINE fullchem_ConvertEquivToAlk
  !
  SUBROUTINE fullchem_HetDropChem( I,         J,         L,         SR,      &
                                   Input_Opt, State_Met, State_Chm          )
    !
    ! Stub for GeosCore/fullchem_mod.F90:fullchem_HetDropChem
    !
    USE gckpp_Precision, ONLY : dp
    USE Input_Opt_Mod,   ONLY : OptInput
    USE State_Chm_Mod,   ONLY : ChmState
    USE State_Met_Mod,   ONLY : MetState
    !
    INTEGER,        INTENT(IN)    :: I, J, L
    TYPE(OptInput), INTENT(IN)    :: Input_Opt
    TYPE(MetState), INTENT(IN)    :: State_Met
    TYPE(ChmState), INTENT(INOUT) :: State_Chm
    REAL(dp),       INTENT(OUT)   :: SR
  END SUBROUTINE fullchem_HetDropChem
  !
  SUBROUTINE fullchem_InitSulfurChem( RC )
    !
    ! Stub for GeosCore/fullchem_mod.F90:fullchem_InitSulfurChem
    !
    INTEGER, INTENT(OUT) :: RC
    !
    RC = 0
  END SUBROUTINE fullchem_InitSulfurChem
  !
  SUBROUTINE fullchem_SulfurAqChem( I,         J,          L,                &
                                    Input_Opt, State_Chm,  State_Grid,       &
                                    State_Met, RC                           )
    !
    ! Stub for GeosCore/fullchem_mod.F90:fullchem_SulfurAqChem
    !
    USE Input_Opt_Mod,  ONLY : OptInput
    USE State_Chm_Mod,  ONLY : ChmState
    USE State_Met_Mod,  ONLY : MetState
    USE State_Grid_Mod, ONLY : GrdState
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
                                     Input_Opt,  State_Chm, State_Diag,      &
                                     State_Grid, State_Met, size_res,        &
                                     RC                                     )
    !
    ! Stub for GeosCore/fullchem_mod.F90:fullchem_SulfurCldChem
    !
    USE Input_Opt_Mod,      ONLY : OptInput
    USE State_Chm_Mod,      ONLY : ChmState
    USE State_Diag_Mod,     ONLY : DgnState
    USE State_Grid_Mod,     ONLY : GrdState
    USE State_Met_Mod,      ONLY : MetState
    !
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
    TYPE(GrdState), INTENT(IN)    :: State_Grid  ! Grid State object
    TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
    INTEGER,        INTENT(IN)    :: I, J, L
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
    TYPE(DgnState), INTENT(INOUT) :: State_Diag  ! Diagnostics State object
    LOGICAL,        INTENT(OUT)   :: size_res    ! Should we call HetDropChem?
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
    !
    !
    RC       = 0
    size_res = .FALSE.
  END SUBROUTINE fullchem_SulfurCldChem
  !
  SUBROUTINE fullchem_UpdateHSO3mAndSO3mm( I,         J,         L,          &
                                           State_Chm, State_Het, State_Met  )
    !
    ! Stub for GeosCore/fullchem_mod.F90:fullchem_UpdateHSO3mAndSO3mm
    !
    USE gckpp_Global,  ONLY : HetState
    USE State_Chm_Mod, ONLY : ChmState
    USE State_Met_Mod, ONLY : MetState
    !
    INTEGER,        INTENT(IN) :: I, J, L
    TYPE(ChmState), INTENT(IN) :: State_Chm
    TYPE(HetState), INTENT(IN) :: State_Het
    TYPE(MetState), INTENT(IN) :: State_Met
  END SUBROUTINE fullchem_UpdateHSO3mAndSO3mm
  !
END MODULE fullchem_SulfurChemFuncs
