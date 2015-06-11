!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: hcox_tools_mod.F90
!
! !DESCRIPTION: Module HCOX\_Tools\_Mod contains a collection of helper
! routines for the HEMCO extensions. 
!\\
!\\
! !INTERFACE: 
!
MODULE HCOX_TOOLS_MOD
!
! !USES:
!
  USE HCO_ERROR_MOD

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: HCOX_Mask 
!
! !MODULE VARIABLES:
!
  CHARACTER(LEN=31), PARAMETER, PUBLIC  :: HCOX_NOMASK = 'none' 
!
! !PRIVATE MEMBER FUNCTIONS:
!
! !REVISION HISTORY:
!  11 Jun 2015 - C. Keller   - Initial version
!EOP
!-----------------------------------------------------------------------------
!BOC
!
! !MODULE INTERFACES: 
!
  INTERFACE HCOX_Mask 
     MODULE PROCEDURE HCOX_Mask_sp2D 
     MODULE PROCEDURE HCOX_Mask_sp3D 
     MODULE PROCEDURE HCOX_Mask_dp2D 
     MODULE PROCEDURE HCOX_Mask_dp3D 
  END INTERFACE HCOX_Mask 

CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: HCOX_Mask_sp2D 
!
! !DESCRIPTION: Applies mask `MaskName` to the passed 2D sp field. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCOX_Mask_sp2D( am_I_Root, HcoState, Arr, MaskName, RC )
!
! !USES:
!
    USE HCO_CALC_MOD,   ONLY : HCO_MaskFld
    USE HCO_STATE_MOD,  ONLY : HCO_State
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   )  :: am_I_Root  ! Root CPU?
    TYPE(HCO_STATE),  POINTER        :: HcoState   ! HcoState obj 
    CHARACTER(LEN=*), INTENT(IN   )  :: MaskName   ! Mask to be used
!
! !INPUT/OUTPUT PARAMETERS:
!
    REAL(sp),         INTENT(INOUT)  :: Arr(:,:)   ! Array to be scaled 
    INTEGER,          INTENT(INOUT)  :: RC         ! Success or failure?
!
! !REVISION HISTORY:
!  11 Jun 2013 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    REAL(sp)            :: MASK(HcoState%NX,HcoState%NY)

    !======================================================================
    ! HCOX_Mask_sp2D begins here
    !======================================================================

    IF ( TRIM(MaskName) /= HCOX_NOMASK ) THEN
       
       ! Get mask field
       CALL HCO_MaskFld ( am_I_Root, HcoState, TRIM(MaskName), MASK, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN

       ! Set array to zero outside of mask region
       Arr = Arr * MASK
    ENDIF

    ! Return w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE HCOX_Mask_sp2D 
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: HCOX_Mask_sp3D 
!
! !DESCRIPTION: Applies mask `MaskName` to the passed 3D sp field. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCOX_Mask_sp3D( am_I_Root, HcoState, Arr, MaskName, RC )
!
! !USES:
!
    USE HCO_CALC_MOD,   ONLY : HCO_MaskFld
    USE HCO_STATE_MOD,  ONLY : HCO_State
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   )  :: am_I_Root  ! Root CPU?
    TYPE(HCO_STATE),  POINTER        :: HcoState   ! HcoState obj 
    CHARACTER(LEN=*), INTENT(IN   )  :: MaskName   ! Mask to be used
!
! !INPUT/OUTPUT PARAMETERS:
!
    REAL(sp),         INTENT(INOUT)  :: Arr(:,:,:) ! Array to be scaled 
    INTEGER,          INTENT(INOUT)  :: RC         ! Success or failure?
!
! !REVISION HISTORY:
!  11 Jun 2013 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    REAL(sp)            :: MASK(HcoState%NX,HcoState%NY)
    INTEGER             :: I, NZ

    !======================================================================
    ! HCOX_Mask_sp3D begins here
    !======================================================================

    IF ( TRIM(MaskName) /= HCOX_NOMASK ) THEN
       
       ! Get mask field
       CALL HCO_MaskFld ( am_I_Root, HcoState, TRIM(MaskName), MASK, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN 

       ! Number of levels
       NZ = SIZE(Arr,3)

       DO I = 1, NZ
          Arr(:,:,I) = Arr(:,:,1) * MASK
       ENDDO  

    ENDIF

    ! Return w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE HCOX_Mask_sp3D 
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: HCOX_Mask_dp2D 
!
! !DESCRIPTION: Applies mask `MaskName` to the passed 2D dp field. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCOX_Mask_dp2D( am_I_Root, HcoState, Arr, MaskName, RC )
!
! !USES:
!
    USE HCO_CALC_MOD,   ONLY : HCO_MaskFld
    USE HCO_STATE_MOD,  ONLY : HCO_State
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   )  :: am_I_Root  ! Root CPU?
    TYPE(HCO_STATE),  POINTER        :: HcoState   ! HcoState obj 
    CHARACTER(LEN=*), INTENT(IN   )  :: MaskName   ! Mask to be used
!
! !INPUT/OUTPUT PARAMETERS:
!
    REAL(dp),         INTENT(INOUT)  :: Arr(:,:)   ! Array to be scaled 
    INTEGER,          INTENT(INOUT)  :: RC         ! Success or failure?
!
! !REVISION HISTORY:
!  11 Jun 2013 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    REAL(sp)            :: MASK(HcoState%NX,HcoState%NY)

    !======================================================================
    ! HCOX_Mask_dp2D begins here
    !======================================================================

    IF ( TRIM(MaskName) /= HCOX_NOMASK ) THEN
       
       ! Get mask field
       CALL HCO_MaskFld ( am_I_Root, HcoState, TRIM(MaskName), MASK, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN 

       ! Set array to zero outside of mask region
       Arr = Arr * MASK

    ENDIF

    ! Return w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE HCOX_Mask_dp2D 
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: HCOX_Mask_dp3D 
!
! !DESCRIPTION: Applies mask `MaskName` to the passed 3D dp field. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCOX_Mask_dp3D( am_I_Root, HcoState, Arr, MaskName, RC )
!
! !USES:
!
    USE HCO_CALC_MOD,   ONLY : HCO_MaskFld
    USE HCO_STATE_MOD,  ONLY : HCO_State
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   )  :: am_I_Root  ! Root CPU?
    TYPE(HCO_STATE),  POINTER        :: HcoState   ! HcoState obj 
    CHARACTER(LEN=*), INTENT(IN   )  :: MaskName   ! Mask to be used
!
! !INPUT/OUTPUT PARAMETERS:
!
    REAL(dp),         INTENT(INOUT)  :: Arr(:,:,:) ! Array to be scaled 
    INTEGER,          INTENT(INOUT)  :: RC         ! Success or failure?
!
! !REVISION HISTORY:
!  11 Jun 2013 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    REAL(sp)            :: MASK(HcoState%NX,HcoState%NY)
    INTEGER             :: I, NZ

    !======================================================================
    ! HCOX_Mask_dp3D begins here
    !======================================================================

    IF ( TRIM(MaskName) /= HCOX_NOMASK ) THEN
      
       ! Get mask field
       CALL HCO_MaskFld ( am_I_Root, HcoState, TRIM(MaskName), MASK, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN 

       ! Number of levels
       NZ = SIZE(Arr,3)
       DO I = 1, NZ
          Arr(:,:,I) = Arr(:,:,1) * MASK
       ENDDO  
    ENDIF

    ! Return w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE HCOX_Mask_dp3D 
!EOC
END MODULE HCOX_TOOLS_MOD
