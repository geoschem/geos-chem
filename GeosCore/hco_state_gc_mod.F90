!------------------------------------------------------------------------------
!                    Harmonized Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: hco_state_gc_mod.F90
!
! !DESCRIPTION: Module hco\_state\_gc\_mod.F90 holds the HEMCO state objects
!  in GEOS-Chem to avoid circular dependencies.
!
! !INTERFACE:
!
MODULE HCO_State_GC_Mod
!
! !USES:
!
  USE HCOX_State_Mod, ONLY : Ext_State
  USE HCO_State_Mod,  ONLY : HCO_State
  USE State_Grid_Mod, ONLY : GrdState
  IMPLICIT NONE
  PRIVATE
!
! !REMARKS:
!  Exists to avoid circular dependencies. Also holds the HEMCO intermediate grid
!  description, State\_Grid\_HCO, for GC-Classic.
!
!  The intermediate grid descriptor uses GrdState to avoid code duplication. It will
!  be initialized by HCO\_Interface\_GC\_Mod.
!
! !REVISION HISTORY:
!  13 Mar 2020 - H.P. Lin    - Initial version.
!  04 Jun 2020 - H.P. Lin    - Now holding HEMCO intermediate grid object.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PUBLIC MODULE VARIABLES:
!
  !--------------------------
  ! %%% Pointers %%%
  !--------------------------

  ! HEMCO state
  TYPE(HCO_State), POINTER, PUBLIC :: HcoState => NULL()

  ! HEMCO extensions state
  TYPE(Ext_State), POINTER, PUBLIC :: ExtState => NULL()

#if defined ( MODEL_CLASSIC )
  !--------------------------
  ! %%% HEMCO Intermediate Grid %%%
  !--------------------------
  TYPE(GrdState), PUBLIC           :: State_Grid_HCO
#endif
!EOC
END MODULE HCO_State_GC_Mod