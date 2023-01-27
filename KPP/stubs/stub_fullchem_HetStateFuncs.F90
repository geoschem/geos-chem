MODULE fullchem_HetStateFuncs
  !
  ! Contains stub routines to satisfy compilation requirements in
  ! fullchem_mod.F90 when building other KPP-based mechanisms.
  !
  IMPLICIT NONE
  PRIVATE
  !
  PUBLIC :: fullchem_SetStateHet
  !
CONTAINS
  !
  SUBROUTINE fullchem_SetStateHet( I,         J,         L,                  &
                                   Input_Opt, State_Chm, State_Met,          &
                                   H,         RC                            )
    !
    ! Stub routine to avoid compilation errors
    !
    USE GcKpp_Global,     ONLY : HetState
    USE GcKpp_Precision
    USE Input_Opt_Mod,    ONLY : OptInput
    USE State_Chm_Mod,    ONLY : ChmState
    USE State_Met_Mod,    ONLY : MetState
    !
    INTEGER,        INTENT(IN)    :: I
    INTEGER,        INTENT(IN)    :: J
    INTEGER,        INTENT(IN)    :: L
    TYPE(OptInput), INTENT(IN)    :: Input_Opt
    TYPE(ChmState), INTENT(IN)    :: State_Chm
    TYPE(MetState), INTENT(IN)    :: State_Met
    TYPE(HetState), INTENT(INOUT) :: H
    INTEGER,        INTENT(OUT)   :: RC
    !
    RC = 0
  END SUBROUTINE fullchem_SetStateHet
  !
END MODULE fullchem_HetStateFuncs
