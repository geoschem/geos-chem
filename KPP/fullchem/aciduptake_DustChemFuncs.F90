!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: aciduptake_DustChemFuncs
!
! !DESCRIPTION: Stub routine for GeosCore/fullchem_mod.F90,
!\\
!\\
! !INTERFACE:

MODULE aciduptake_DustChemFuncs
!
! !USES:
!
  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: aciduptake_InitDustChem

CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: aciduptake_InitDustChem
!
! !DESCRIPTION: Stub routine for GeosCore/fullchem_mod.F90,
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE aciduptake_InitDustChem( RC )
!
! !USES:
!
    USE ErrCode_Mod
!
! !INPUT/OUTPUT PARAMETERS: 
!
    INTEGER, INTENT(INOUT) :: RC   ! Success or failure
!EOP
!------------------------------------------------------------------------------
!BOC
    RC = GC_SUCCESS
  END SUBROUTINE aciduptake_InitDustChem
!EOC
END MODULE aciduptake_DustChemFuncs
