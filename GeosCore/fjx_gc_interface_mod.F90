!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: FJX_GC_Interface_Mod.F90
!
! !DESCRIPTION: Module FJX\_GC\_Interface\_Mod.F90 contains routines and
! variables to interface GEOS-Chem and FAST-JX. It sets the FAST-JX
! state object (FjxState) from other GEOS-Chem state objects.
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

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC  :: FjxState_GC_Set
  PUBLIC  :: FjxState_Print
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
! !IROUTINE: FjxState_GC_Set
!
! !DESCRIPTION: Subroutine FjxState\_GC\_Set sets the FAST-JX state object
!  values and pointers from other GEOS-Chem state objects.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE FjxState_GC_Set( Input_Opt,  State_Chm, State_Grid, State_Met, &
                              State_Diag, FjxState, RC )
!
! !USES:
!
    USE ERRCODE_MOD
    USE Fast_JX_Mod,        ONLY : Fjx_State, FjxState_Init
    USE Input_Opt_Mod,      ONLY : OptInput
    USE State_Chm_Mod,      ONLY : ChmState
    USE State_Chm_Mod,      ONLY : Ind_
    USE State_Diag_Mod,     ONLY : DgnState
    USE State_Grid_Mod,     ONLY : GrdState
    USE State_Met_Mod,      ONLY : MetState
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput),  INTENT(IN)    :: Input_Opt   ! Input options
    TYPE(ChmState),  INTENT(IN)    :: State_Chm   ! Chemistry State object
    TYPE(GrdState),  INTENT(IN)    :: State_Grid  ! Grid State object
    TYPE(MetState),  INTENT(IN)    :: State_Met   ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(DgnState),  INTENT(INOUT) :: State_Diag  ! Diagnostics State object
    TYPE(Fjx_State), INTENT(INOUT) :: FjxState    ! FJX state object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,         INTENT(OUT)   :: RC          ! Success or failure?
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
    ! FjxState_GC_Set begins here!
    !=================================================================

    ! Initialize
    RC          = GC_SUCCESS
    ErrMsg      = ''
    ThisLoc     = ' -> at FjxState_GC_Set (in module GeosCore/fjx_gc_interface_mod.F90)'

    ! Initialize FjxState object
    CALL FjxState_Init( FjxState, RC )

    ! Set values from state_chm: todo

    ! Set Config from GEOS-Chem object Input_Opt
    FjxState%Config%DryRun            = Input_Opt%DryRun
    FjxState%Config%amIRoot           = Input_Opt%amIRoot
    FjxState%Config%LPRT              = Input_Opt%LPRT
    FjxState%Config%USE_ONLINE_O3     = Input_Opt%USE_ONLINE_O3
    FjxState%Config%LBRC              = Input_Opt%LBRC
    FjxState%Config%LUCX              = Input_Opt%LUCX
    FjxState%Config%FAST_JX_DIR       = Input_Opt%FAST_JX_DIR
    FjxState%Config%CHEM_INPUTS_DIR   = Input_Opt%CHEM_INPUTS_DIR
    FjxState%Config%hvAerNIT          = Input_Opt%hvAerNIT
    FjxState%Config%hvAerNIT_JNITs    = Input_Opt%hvAerNIT_JNITs
    FjxState%Config%hvAerNIT_JNIT     = Input_Opt%hvAerNIT_JNIT
    FjxState%Config%JNITChanA         = Input_Opt%JNITChanA
    FjxState%Config%JNITChanB         = Input_Opt%JNITChanB
#if defined( MODEL_GEOS )
    FjxState%Config%FJX_EXTRAL_ITERMAX= Input_Opt%FJX_EXTRAL_ITERMAX
    FjxState%Config%FJX_EXTRAL_ERR    = Input_Opt%FJX_EXTRAL_ERR
#endif
    FjxState%Config%NWVSELECT         = Input_Opt%NWVSELECT
   
    ALLOCATE( FjxState%Config%WVSELECT( Input_Opt%NWVSELECT ), STAT=RC )
    CALL GC_CheckVar( 'FjxState%Config%WVSELECT', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    FjxState%Config%WVSELECT(1:FjxState%Config%NWVSELECT) = Input_Opt%WVSELECT(1:Input_Opt%NWVSELECT)

    ! Set Grid from GEOS-Chem object State_Grid
    FjxState%Grid%NX         =  State_Grid%NX
    FjxState%Grid%NY         =  State_Grid%NY
    FjxState%Grid%NZ         =  State_Grid%NZ
    FjxState%Grid%MaxChemLev =  State_Grid%MaxChemLev
    FjxState%Grid%YMID       => State_Grid%YMID

    ! Set Met from GEOS-Chem object State_Met
    FjxState%Met%ChemGridLev => State_Met%ChemGridLev
    FjxState%Met%SuncosMid   => State_Met%SuncosMid
    FjxState%Met%UVALBEDO    => State_Met%UVALBEDO
    FjxState%Met%PEDGE       => State_Met%PEDGE
    FjxState%Met%T           => State_Met%T
    FjxState%Met%OPTD        => State_Met%OPTD
    FjxState%Met%CLDF        => State_Met%CLDF

    ! State_Diag vals needed: todo

  END SUBROUTINE FjxState_GC_Set
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: FjxState_Print
!
! !DESCRIPTION: Subroutine FjxState\_Print prints information about the
!  FAST-JX object.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE FjxState_Print( Input_Opt, FjxState, RC )
!
! !USES:
!
    USE ERRCODE_MOD
    USE Fast_jx_mod,        ONLY : Fjx_State
    USE Input_Opt_Mod,      ONLY : OptInput
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput),  INTENT(IN)  :: Input_Opt   ! Input options
    TYPE(Fjx_State), INTENT(IN)  :: FjxState    ! FJX state object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,         INTENT(OUT) :: RC          ! Success or failure?
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
    ! FjxState_Print begins here!
    !=================================================================

    ! Initialize
    RC      = GC_SUCCESS
    ErrMsg  = ''
    ThisLoc = ' -> at FjxState_Print (in module GeosCore/fjx_gc_interface_mod.F90)'

    IF ( input_opt%AmIRoot ) THEN

       ! input configuration
       print *, ' '
       print *, '=== Printing FAST-JX State Object ==='
       print *, 'FjxState%Config%amIRoot: ',         FjxState%Config%amIRoot
       print *, 'FjxState%Config%DryRun: ',          FjxState%Config%DryRun
       print *, 'FjxState%Config%LPRT: ',            FjxState%Config%LPRT
       print *, 'FjxState%Config%USE_ONLINE_O3: ',   FjxState%Config%USE_ONLINE_O3
       print *, 'FjxState%Config%LUCX: ',            FjxState%Config%LUCX
       print *, 'FjxState%Config%LBRC: ',            FjxState%Config%LBRC
       print *, 'FjxState%Config%FAST_JX_DIR: ',     TRIM(FjxState%Config%FAST_JX_DIR)
       print *, 'FjxState%Config%CHEM_INPUTS_DIR: ', TRIM(FjxState%Config%CHEM_INPUTS_DIR)
       print *, 'FjxState%Config%hvAerNIT: ',        FjxState%Config%hvAerNIT
       print *, 'FjxState%Config%hvAerNIT_JNITs: ',  FjxState%Config%hvAerNIT_JNITs
       print *, 'FjxState%Config%hvAerNIT_JNIT: ',   FjxState%Config%hvAerNIT_JNIT
       print *, 'FjxState%Config%JNITChanA: ',       FjxState%Config%JNITChanA
       print *, 'FjxState%Config%JNITChanB: ',       FjxState%Config%JNITChanB
       print *, 'FjxState%Config%NWVSELECT: ',       FjxState%Config%NWVSELECT       
       print *, 'FjxState%Config%WVSELECT: ',        FjxState%Config%WVSELECT
#if defined( MODEL_GEOS )
       print *, 'FjxState%Config%FJX_EXTRAL_ITERMAX: ', &
                                                     FjxState%Config%FJX_EXTRAL_ITERMAX
       print *, 'FjxState%Config%FJX_EXTRAL_ERR: ', FjxState%Config%FJX_EXTRAL_ERR
#endif
       
       ! grid quantities
       print *, 'FjxState%Grid%NX: ',         FjxState%Grid%NX
       print *, 'FjxState%Grid%NY: ',         FjxState%Grid%NY
       print *, 'FjxState%Grid%NZ: ',         FjxState%Grid%NZ
       print *, 'FjxState%Grid%MaxChemLev: ', FjxState%Grid%MaxChemLev
       print *, 'FjxState%Grid%YMID max: ',   MAXVAL(FjxState%Grid%YMID)
       
       ! met quantities
       print *, 'FjxState%Met%ChemGridLev max: ',MAXVAL(FjxState%Met%ChemGridLev)
       print *, 'FjxState%Met%SuncosMid max: ',  MAXVAL(FjxState%Met%SuncosMid)
       print *, 'FjxState%Met%UVALBEDO max: ',   MAXVAL(FjxState%Met%UVALBEDO)
       print *, 'FjxState%Met%PEDGE max: ',      MAXVAL(FjxState%Met%PEDGE)
       print *, 'FjxState%Met%T max: ',          MAXVAL(FjxState%Met%T)
       print *, 'FjxState%Met%OPTD max: ',       MAXVAL(FjxState%Met%OPTD)
       print *, 'FjxState%Met%CLDF max: ',       MAXVAL(FjxState%Met%CLDF)
       
       
       ! State_Diag vals needed (move elsewhere, last)
      
       print *, '======'

    ENDIF

  END SUBROUTINE FjxState_Print
!EOC
END MODULE FJX_GC_Interface_Mod
