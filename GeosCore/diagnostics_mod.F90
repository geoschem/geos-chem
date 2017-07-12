!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: diagnostics_mod.F90
!
! !DESCRIPTION: Module diagnostics\_mod.F90 is the first crack of the new
! diagnostics package for GEOS-Chem.
!
! !INTERFACE:
!
MODULE Diagnostics_Mod
!
! !USES:
!
  USE CMN_SIZE_Mod
  USE ErrCode_Mod
  USE Error_Mod,          ONLY : Error_Stop
  USE HCO_Error_Mod
  USE HCO_INTERFACE_MOD,  ONLY : HcoState
  USE Precision_Mod

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC  :: Diagnostics_Write
!
! !REVISION HISTORY:
!  09 Jan 2015 - C. Keller   - Initial version. 
!  14 Jan 2016 - E. Lundgren - Add several GEOS-Chem diagnostics
!  29 Jan 2016 - E. Lundgren - Update diagnostics for recent HEMCO updates
!  20 Jul 2016 - R. Yantosca - Replace #ifdef DEVEL with #ifdef DIAG_DEVEL
!  29 Nov 2016 - R. Yantosca - grid_mod.F90 is now gc_grid_mod.F90
!  12 Jul 2017 - R. Yantosca - Remove everything in #if defined( NC_DIaG )
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE TYPES:
!
  ! Species ID flags
  INTEGER :: id_Rn, id_Pb, id_Be7
!BOC
CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Diagnostics_Write
!
! !DESCRIPTION: Subroutine Diagnostics\_Write writes the GEOS-Chem diagnostics
! to disk. If the variable RESTART is set to true, all GEOS-Chem diagnostics
! are passed to the restart file.
!\\
!\\
! The Diagnostics\_Write routine is called from main.F, at the end of the time
! loop and during cleanup (to write the restart file).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Diagnostics_Write ( am_I_Root, Input_Opt, State_Chm, RESTART, RC ) 
!
! !USES:
!
    USE HCO_INTERFACE_MOD,  ONLY : HcoState, GetHcoID
    USE HCOI_GC_MAIN_MOD,   ONLY : HCOI_GC_WriteDiagn 
    USE HCOIO_Diagn_Mod,    ONLY : HCOIO_Diagn_WriteOut
    USE Input_Opt_Mod,      ONLY : OptInput
    USE State_Chm_Mod,      ONLY : ChmState
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   )  :: am_I_Root  ! Are we on the root CPU?
    LOGICAL,          INTENT(IN   )  :: RESTART    ! Write restart file? 
    TYPE(OptInput),   INTENT(IN   )  :: Input_Opt  ! Input Options object
    TYPE(ChmState),   INTENT(IN   )  :: State_Chm  ! Chemistry State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(OUT  )  :: RC         ! Failure or success
!
! !REVISION HISTORY: 
!  09 Jan 2015 - C. Keller   - Initial version
!  15 Jan 2015 - R. Yantosca - Now accept Input_Opt via the arg list
!  29 Apr 2016 - R. Yantosca - Don't initialize pointers in declaration stmts
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    CHARACTER(LEN=255)       :: LOC = 'Diagnostics_Write (diagnostics_mod.F90)'

    !=======================================================================
    ! Diagnostics_Write begins here 
    !=======================================================================

    ! Write HEMCO diagnostics
    CALL HCOI_GC_WriteDiagn( am_I_Root, Input_Opt, RESTART, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    ! Leave w/ success
    RC = GC_SUCCESS

  END SUBROUTINE Diagnostics_Write
!EOC
END MODULE Diagnostics_Mod
