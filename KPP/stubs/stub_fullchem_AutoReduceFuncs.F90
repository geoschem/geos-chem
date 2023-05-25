!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: stub_fullchem_AutoReduceFuncs
!
! !DESCRIPTION: Stub routines corresponding to fullchem_AutoReduceFuncs.F90.
!  This allows us to satisfy compilation requirements for fullchem_mod.F90
!  when building other KPP mechanisms (e.g. Hg, carbon, etc.)
!\\
!\\
! !INTERFACE:
!
MODULE fullchem_AutoReduceFuncs
!
! !USES:
!
  IMPLICIT NONE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: fullchem_AR_KeepHalogensActive
  PUBLIC :: fullchem_AR_SetKeepActive
  PUBLIC :: fullchem_AR_UpdateKppDiags
!
!EOP
!------------------------------------------------------------------------------
!BOC
CONTAINS
!EOC
  SUBROUTINE fullchem_AR_KeepHalogensActive( doPrint )
    !
    ! Sets halogen species to "fast" for the rosenbrock_autoreduce integrator.
    !
    LOGICAL, INTENT(IN) :: doPrint
    !
  END SUBROUTINE fullchem_AR_KeepHalogensActive
  !
  SUBROUTINE fullchem_AR_SetKeepActive( option ) 
    !
    ! Abstracts setting the rosenbrock_autoreduce keepActive flag
    ! out of fullchem_mod.F90
    !
    LOGICAL, INTENT(IN) :: option
    !
  END SUBROUTINE fullchem_AR_SetKeepActive
  !
  SUBROUTINE fullchem_AR_UpdateKppDiags( I, J, L, RSTATE, State_Diag )
    !
    ! Updates KPP diagnostics for the rosenbrock_autoreduce solver
    ! 
    USE gckpp_Precision
    USE State_Diag_Mod, ONLY : DgnState
    !
    INTEGER,        INTENT(IN)    :: I, J, L
    REAL(dp),       INTENT(IN)    :: RSTATE(20)
    TYPE(DgnState), INTENT(INOUT) :: State_Diag
    !
  END SUBROUTINE fullchem_AR_UpdateKppDiags

  SUBROUTINE fullchem_AR_SetIntegratorOptions( Input_Opt, State_Chm,         &
                                               State_Met, FirstChem,         &
                                               I,         J,         L,      &
                                               ICNTRL,    RCNTRL )
    !
    USE gckpp_Parameters
    USE gckpp_Precision
    USE Input_Opt_Mod, ONLY : OptInput
    USE State_Chm_Mod, ONLY : ChmState
    USE State_Met_Mod, ONLY : MetState
    !
    TYPE(OptInput), INTENT(IN)    :: Input_Opt
    TYPE(ChmState), INTENT(IN)    :: State_Chm
    TYPE(MetState), INTENT(IN)    :: State_Met
    LOGICAL,        INTENT(IN)    :: FirstChem
    INTEGER,        INTENT(IN)    :: I, J, L
    INTEGER,        INTENT(INOUT) :: ICNTRL(20)
    REAL(dp),       INTENT(INOUT) :: RCNTRL(20)
  END SUBROUTINE fullchem_AR_SetIntegratorOptions
!EOC
END MODULE fullchem_AutoReduceFuncs
