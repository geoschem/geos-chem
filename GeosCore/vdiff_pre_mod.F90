!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: vdiff_pre_mod.F90
!
! !DESCRIPTION: Module VDIFF\_PRE\_MOD contains variables used in VDIFF\_MOD.
!\\
!\\
! !INTERFACE:
!
MODULE VDIFF_PRE_MOD
!
! !USES:
!
  IMPLICIT NONE

  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: Set_VDIFF_VALUES
!
! !PUBLIC DATA MEMBERS:
!
  LOGICAL, PUBLIC      :: prtDebug           ! Passes prtDebug to vdiff_mod
  LOGICAL, PUBLIC      :: LTURB              ! Passes LTURB to vdiff_mod
  INTEGER, PUBLIC      :: PCNST              ! Passes N_TRACERS to vdiff_mod
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: set_vdiff_values
!
! !DESCRIPTION: Subroutine SET\_VDIFF\_VALUES initializes the PCNST value, which
!  is the number of advected species.  This is needed in vdiff\_mod.F90.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Set_VDIFF_VALUES( Input_Opt, State_Chm, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod,      ONLY : OptInput
    USE State_Chm_Mod,      ONLY : ChmState
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)  :: Input_Opt   ! Input Options object
    TYPE(ChmState), INTENT(IN)  :: State_Chm   ! Chemistry State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT) :: RC          ! Success or failure?
!
! !REMARKS:
!  This routine has to be called after routine READ_INPUT_FILE.
!
! !REVISION HISTORY:
!  24 Jun 2014 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
    ! Assume success
    RC    = GC_SUCCESS

    !=====================================================================
    ! The following fields of Input_Opt have to be set in module
    ! variables.  This will allow us to pass these to vdiff_mod.F90,
    ! now that logical_mod.F and tracer_mod.F have been retired.
    !=====================================================================

    ! Number of advected species
    PCNST = State_Chm%nAdvect

    ! Debug print?
    prtDebug = ( Input_Opt%LPRT .and. Input_Opt%amIRoot )

    ! Use PBL mixing?
    LTURB = Input_Opt%LTURB

  END SUBROUTINE SET_VDIFF_VALUES
!EOC
END MODULE VDIFF_PRE_MOD

