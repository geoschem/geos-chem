!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: hco_geotools_mod.F90
!
! !DESCRIPTION: Module HCO\_GeoTools\_Mod contains a collection of 
! helper routines for extracting geographical information. These 
! routines are based upon GEOS-5 data and may need to be revised
! for other met. fields! 
! \\
! !INTERFACE: 
!
MODULE HCO_GeoTools_Mod
!
! !USES:
!
  USE HCO_Error_Mod

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: HCO_LandType
  INTERFACE HCO_LandType
     MODULE PROCEDURE HCO_LandType_Dp
     MODULE PROCEDURE HCO_LandType_Sp
  END INTERFACE HCO_LandType
!
! !PRIVATE MEMBER FUNCTIONS:
!
  PRIVATE:: HCO_LandType_Dp
  PRIVATE:: HCO_LandType_Sp
!
! !REVISION HISTORY:
!  18 Dec 2013 - C. Keller   - Initialization
!  01 Jul 2014 - R. Yantosca - Cosmetic changes in ProTeX headers
!  01 Jul 2014 - R. Yantosca - Now use F90 free-format indentation
!EOP
!------------------------------------------------------------------------------
!BOC
CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !FUNCTION: HCO_LandType_Sp
!
! !DESCRIPTION: Function HCO\_LANDTYPE returns the land type based upon 
!  the land water index (0=water,1=land,2=ice) and the surface albedo.
!  Inputs are in single precision.
!\\
!\\
! !INTERFACE:
!
  FUNCTION HCO_LandType_Sp( WLI, Albedo ) Result ( LandType )
!
! !INPUT PARAMETERS:
!
    REAL(sp), INTENT(IN) :: WLI       ! Land type: 0=water,1=land,2=ice
    REAL(sp), INTENT(IN) :: Albedo    ! Surface albedo
!
! !RETURN VALUE
!
    INTEGER              :: LandType  ! Land type: 0=water,1=land,2=ice
!
! !REMARKS:
!  This function is largely based on the GEOS-Chem functions in dao_mod.F. 
!
! !REVISION HISTORY:
!  18 Dec 2013 - C. Keller - Initialization!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
    ! Set threshold for albedo over ice. All surfaces w/ an albedo above 
    ! this value will be classified as ice.
    REAL(sp), PARAMETER :: albd_ice = 0.695_sp

    !--------------------------
    ! HCO_LANDTYPE begins here
    !--------------------------

    ! Water:
    IF ( NINT(WLI) == 0 .AND. Albedo < albd_ice ) THEN
       LandType = 0

    ! Land: 
    ELSEIF ( NINT(WLI) == 1 .AND. Albedo < albd_ice ) THEN
       LandType = 1

    ! Ice:
    ELSEIF ( NINT(WLI) == 2 .OR. Albedo >= albd_ice ) THEN
       LandType = 2
    ENDIF

  END FUNCTION HCO_LandType_Sp
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !FUNCTION: HCO_LandType_Dp 
!
! !DESCRIPTION: Function HCO\_LandType\_Dp returns the land type based upon 
! the land water index (0=water,1=land,2=ice) and the surface albedo.
! Inputs are in double precision.
!\\
!\\
! !INTERFACE:
!
  FUNCTION HCO_LandType_Dp( WLI, Albedo ) Result ( LandType )
!
! !INPUT PARAMETERS:
!
    REAL(dp), INTENT(IN) :: WLI       ! Land type: 0=water,1=land,2=ice
    REAL(dp), INTENT(IN) :: Albedo    ! Surface albedo
!
! !RETURN VALUE:
!
    INTEGER              :: LandType  ! Land type: 0=water,1=land,2=ice
!
! !REMARKS:
!  This function is largely based on the GEOS-Chem functions in dao_mod.F. 
!
! !REVISION HISTORY:
!  18 Dec 2013 - C. Keller - Initialization
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS::
!
    ! Set threshold for albedo over ice. All surfaces w/ an albedo above 
    ! this value will be classified as ice.
    REAL(dp), PARAMETER :: albd_ice = 0.695_dp

    !--------------------------
    ! HCO_LANDTYPE begins here
    !--------------------------

    ! Water:
    IF ( NINT(WLI) == 0 .AND. Albedo < albd_ice ) THEN
       LandType = 0
         
    ! Land: 
    ELSEIF ( NINT(WLI) == 1 .AND. Albedo < albd_ice ) THEN
       LandType = 1

    ! Ice:
    ELSEIF ( NINT(WLI) == 2 .OR. Albedo >= albd_ice ) THEN
       LandType = 2
    ENDIF

  END FUNCTION HCO_LandType_Dp
!EOC
END MODULE HCO_GeoTools_MOd
!EOM
