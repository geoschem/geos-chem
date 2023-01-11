!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: stub_fullchem_AutoReduceFuncs
!
! !DESCRIPTION: Stub routines corresponding to fullchem_AutoReduceFuncs.F90.
!  This allows us to satisfy compilation requirements for fullchem_mod.F90
!  when building other KPP mechanisms (e.g. Hg, carboncycle, etc.)
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
!EOC
END MODULE fullchem_AutoReduceFuncs