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
!   19 Apr 2002 - R. Yantosca - Initial version
!  (1 ) Now read uvalbedo file directly from DATA_DIR/uvalbedo_200111
!        subdirectory.  (bmy, 4/2/02)
!  (2 ) Now divide module header into MODULE PRIVATE, MODULE VARIABLES, and
!        MODULE ROUTINES sections. (bmy, 5/28/02)
!  (3 ) Now references "error_mod.f" (bmy, 10/15/02)
!  (4 ) Minor modification in READ_UVALBEDO (bmy, 3/14/03)
!  (5 ) Now references "directory_mod.f" (bmy, 7/20/04)
!  (6 ) Bug fix for GCAP grid in READ_UVALBEDO (bmy, 8/16/05)
!  (7 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!  24 Nov 2014 - M. Yannetti - Added PRECISION_MOD
!  17 Dec 2014 - R. Yantosca - Leave time/date variables as 8-byte
!  12 Jan 2015 - R. Yantosca - Remove CLEANUP_UVALBEDO routine
!  04 Mar 2015 - R. Yantosca - UV albedo data now comes via HEMCO
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
  SUBROUTINE Get_UValbedo( am_I_root, Input_Opt, State_Met, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE HCO_INTERFACE_MOD,  ONLY : HcoState
    USE HCO_EmisList_Mod,   ONLY : HCO_GetPtr 
    USE Input_Opt_Mod,      ONLY : OptInput
    USE State_Met_Mod,      ONLY : MetState
!
!
! !INPUT PARAMETERS:
!
    LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU?
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
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
!  29 Apr 2016 - R. Yantosca - Don't initialize pointers in declaration stmts
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!   
    ! Scalars
    LOGICAL :: FND

    ! Pointers
    REAL(f4), POINTER :: Ptr2D(:,:)

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

    ! Nullify pointer
    Ptr2D => NULL()

    ! Get the pointer to the UV albedo data in the HEMCO data structure
    CALL HCO_GetPtr( am_I_Root, HcoState, 'UV_ALBEDO', Ptr2D, RC, FOUND=FND )

    ! Trap potential errors
    IF ( RC /= GC_SUCCESS .or. ( .not. FND ) ) THEN
       ErrMsg = 'Could not find UV_ALBEDO in HEMCO data list!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN       
    ENDIF

    ! Add to State_Met
    State_Met%UVALBEDO = Ptr2D(:,:) 

    ! Free the pointer
    Ptr2d => NULL()

  END SUBROUTINE Get_UValbedo
!EOC
END MODULE UValbedo_Mod
