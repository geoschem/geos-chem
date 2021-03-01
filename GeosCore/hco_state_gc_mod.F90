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
  IMPLICIT NONE
  PRIVATE
!
! !REMARKS:
!  Exists to avoid circular dependencies.
!
! !REVISION HISTORY:
!  13 Mar 2020 - H.P. Lin    - Initial version.
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
!EOC
END MODULE HCO_State_GC_Mod