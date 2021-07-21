!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: fjx_gc_interface_mod.F90
!
! !DESCRIPTION: Module fjx\_gc\_interface\_mod.F90 contains routines and
! variables to interface GEOS-Chem and FAST-JX. It contains the FAST-JX state
! state object (FjxState) as well as init-run-finalize driver routines
! to run FAST-JX within GEOS-Chem.
!\\
!\\
! !INTERFACE:
!
MODULE FJX_GC_Interface_Mod
!
! !USES:
!
  USE Precision_Mod
  USE Error_Mod

  ! Import the FJX states and their types from the state container
!  USE HCOX_State_Mod, ONLY : Ext_State

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC  :: Set_State_FJX
!
! !PRIVATE MEMBER FUNCTIONS:
!
!
! !REMARKS:
!
! !REVISION HISTORY:
!  20 Jul 2021 - E. Lundgren   - Initial version.
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
! !IROUTINE: set_state_fjx
!
! !DESCRIPTION: Subroutine SET\_STATE\_FJX allocates and initializes
!  fields of the FAST-JX object.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE SET_STATE_FJX( Input_Opt,  State_Chm, State_Grid, State_Met, &
                            State_Diag, State_FJX, RC )
!
! !USES:
!
    USE ERRCODE_MOD
    USE Fast_jx_mod,        ONLY : FjxState, Init_State_FJX
    USE Input_Opt_Mod,      ONLY : OptInput
    USE State_Chm_Mod,      ONLY : ChmState
    USE State_Chm_Mod,      ONLY : Ind_
    USE State_Diag_Mod,     ONLY : DgnState
    USE State_Grid_Mod,     ONLY : GrdState
    USE State_Met_Mod,      ONLY : MetState
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input options
    TYPE(ChmState), INTENT(IN)    :: State_Chm   ! Chemistry State object
    TYPE(GrdState), INTENT(IN)    :: State_Grid  ! Grid State object
    TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(DgnState), INTENT(INOUT) :: State_Diag  ! Diagnostics State object

!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(FjxState), INTENT(INOUT) :: State_Fjx ! FJX state object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT) :: RC          ! Success or failure?
!
! !REVISION HISTORY:
!  20 Jul 2021 - E. Lundgren
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    CHARACTER(LEN=255) :: ErrMsg, ThisLoc

    !=================================================================
    ! SET_STATE_FJX begins here!
    !=================================================================

    ! Initialize
    RC          = GC_SUCCESS
    ErrMsg      = ''
    ThisLoc     = ' -> at Set_FJX_State (in module GeosCore/fast_jx_mod.F90)'

    ! Nullify or zero all State_FJX variables
    CALL Init_State_FJX( State_FJX, RC )

    ! Set parameters - do later

    ! Set values from state_chm
    ! others from state_chm:

    ! input configuration
    state_fjx%config%DryRun            = Input_Opt%DryRun
    state_fjx%config%amIRoot           = Input_Opt%amIRoot
    state_fjx%config%LPRT              = Input_Opt%LPRT
    state_fjx%config%USE_ONLINE_O3     = Input_Opt%USE_ONLINE_O3
    state_fjx%config%LBRC              = Input_Opt%LBRC
    state_fjx%config%LUCX              = Input_Opt%LUCX
    state_fjx%config%FAST_JX_DIR       = Input_Opt%FAST_JX_DIR
    state_fjx%config%CHEM_INPUTS_DIR   = Input_Opt%CHEM_INPUTS_DIR
    state_fjx%config%hvAerNIT          = Input_Opt%hvAerNIT
    state_fjx%config%hvAerNIT_JNITs    = Input_Opt%hvAerNIT_JNITs
    state_fjx%config%hvAerNIT_JNIT     = Input_Opt%hvAerNIT_JNIT
    state_fjx%config%JNITChanA         = Input_Opt%JNITChanA
    state_fjx%config%JNITChanB         = Input_Opt%JNITChanB
#if defined( MODEL_GEOS )
    state_fjx%config%FJX_EXTRAL_ITERMAX= Input_Opt%FJX_EXTRAL_ITERMAX
    state_fjx%config%FJX_EXTRAL_ERR    = Input_Opt%FJX_EXTRAL_ERR
#endif
    state_fjx%config%NWVSELECT         = Input_Opt%NWVSELECT
   
    ALLOCATE( state_fjx%config%WVSELECT( Input_Opt%NWVSELECT ), STAT=RC )
    CALL GC_CheckVar( 'state_fjx%config%WVSELECT', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    state_fjx%config%WVSELECT(1:state_fjx%config%NWVSELECT) = Input_Opt%WVSELECT(1:Input_Opt%NWVSELECT)

    ! grid quantities
    state_fjx%grid%NX         = State_Grid%NX
    state_fjx%grid%NY         = State_Grid%MaxChemLev

    ! met quantities
    state_fjx%met%PEDGE     => State_Met%PEDGE
    state_fjx%met%T         => State_Met%T
    state_fjx%met%SuncosMid => State_Met%SuncosMid
    state_fjx%met%UVALBEDO  => State_Met%UVALBEDO
    state_fjx%met%OPTD      => State_Met%OPTD
    state_fjx%met%CLDF      => State_Met%CLDF


    ! State_Diag vals needed (move elsewhere, last)

  END SUBROUTINE Set_STATE_FJX
!EOC
END MODULE FJX_GC_Interface_Mod
