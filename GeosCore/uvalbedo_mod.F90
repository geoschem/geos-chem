!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: uvalbedo_mod.F90
!
! !DESCRIPTION: Module UVALBEDO\_MOD contains variables and routines for
!  reading the UV Albedo data.  This data is required by the FAST-JX photolysis
!  module.  UV albedo data will now be obtained from the HEMCO data structure.
!\\
!\\
! !INTERFACE:
!
MODULE UValbedo_Mod
!
! !USES:
!
  USE Precision_Mod    ! For GEOS-Chem Precision (fp)

  IMPLICIT NONE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: Get_UValbedo
!
! !REMARKS:
!  References:
!  ============================================================================
!  Herman, J.R and Celarier, E.A., "Earth surface reflectivity climatology
!    at 340-380 nm from TOMS data", __J. Geophys. Res__, Vol. 102, D23,
!    pp. 28003-28011, Dec 20, 1997.
!
! !REVISION HISTORY:
!  06 Jan 2015 - R. Yantosca - Initial version
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
! !IROUTINE: get_uvalbedo
!
! !DESCRIPTION: Copies the UV Albedo data from the HEMCO data structure
!  into the State\_Met derived type object.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Get_UValbedo( Input_Opt, State_Grid, State_Met, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE HCO_Utilities_GC_Mod, ONLY : HCO_GC_EvalFld
    USE Input_Opt_Mod,        ONLY : OptInput
    USE State_Met_Mod,        ONLY : MetState
    USE State_Grid_Mod,       ONLY : GrdState
!
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
    TYPE(GrdState), INTENT(IN)    :: State_Grid
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(MetState), INTENT(INOUT) :: State_Met   ! Meteorology State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REVISION HISTORY:
!  06 Jan 2015 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Strings
    CHARACTER(LEN=255) :: ErrMsg
    CHARACTER(LEN=255) :: ThisLoc

    !=======================================================================
    ! READ_UVALBEDO begins here!
    !=======================================================================

    ! Assume success
    RC      = GC_SUCCESS
    ErrMsg  = ''
    ThisLoc = ' -> at Get_UValbedo (in GeosCore/uvalbedo_mod.F90)'

    ! Skip unless we are doing a fullchem or aerosol-only simulation
    IF ( ( .not. Input_Opt%ITS_A_FULLCHEM_SIM ) .and. &
         ( .not. Input_Opt%ITS_AN_AEROSOL_SIM ) ) THEN
       RETURN
    ENDIF

    ! Evalulate the UV albedo from HEMCO
    CALL HCO_GC_EvalFld( Input_Opt, State_Grid, 'UV_ALBEDO', State_Met%UVALBEDO, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Could not find UV_ALBEDO in HEMCO data list!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

  END SUBROUTINE Get_UValbedo
!EOC
END MODULE UValbedo_Mod
