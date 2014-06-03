!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!
! !MODULE: hco_geotools_mod 
!
! !DESCRIPTION: Module HCO\_GEOTOOLS\_MOD contains a collection of 
! helper routines for extracting geographical information. These 
! routines are based upon GEOS-5 data and may need to be revised
! for other met. fields! 
! \\
! !INTERFACE: 
!
      MODULE HCO_GEOTOOLS_MOD
!
! !USES:
!
      USE HCO_ERROR_MOD

      IMPLICIT NONE
      PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC :: HCO_LANDTYPE
!
! !REVISION HISTORY:
!  18 Dec 2013 - C. Keller - Initialization
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE MEMBER FUNCTIONS:
!
!
! !MODULE INTERFACES: 
!
      INTERFACE HCO_LANDTYPE
         MODULE PROCEDURE HCO_LANDTYPE_DP
         MODULE PROCEDURE HCO_LANDTYPE_SP
      END INTERFACE

      CONTAINS
!
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !FUNCTION: HCO_LANDTYPE 
!
! !DESCRIPTION: Function HCO\_LANDTYPE returns the land type based upon 
! the land water index (0=water,1=land,2=ice) and the surface albedo.
!\\
!\\
! !INTERFACE:
!
      FUNCTION HCO_LANDTYPE_SP( WLI, Albedo ) Result ( LandType )
!
! !ARGUMENTS:
!
      REAL(sp), INTENT(IN   ) :: WLI       ! land-water index: 0=water,1=land,2=ice
      REAL(sp), INTENT(IN   ) :: Albedo    ! Surface albedo
      INTEGER                 :: LandType  ! Land type: 0=water,1=land,2=ice
!
! ! NOTES: This function is largely based on the GEOS-Chem functions 
! in dao_mod.F. 
!
! !REVISION HISTORY:
!  18 Dec 2013 - C. Keller - Initialization
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
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

      END FUNCTION HCO_LANDTYPE_SP 
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !FUNCTION: HCO_LANDTYPE 
!
! !DESCRIPTION: Function HCO\_LANDTYPE returns the land type based upon 
! the land water index (0=water,1=land,2=ice) and the surface albedo.
!\\
!\\
! !INTERFACE:
!
      FUNCTION HCO_LANDTYPE_DP( WLI, Albedo ) Result ( LandType )
!
! !ARGUMENTS:
!
      REAL(dp), INTENT(IN   ) :: WLI       ! land-water index: 0=water,1=land,2=ice
      REAL(dp), INTENT(IN   ) :: Albedo    ! Surface albedo
      INTEGER                 :: LandType  ! Land type: 0=water,1=land,2=ice
!
! ! NOTES: This function is largely based on the GEOS-Chem functions 
! in dao_mod.F. 
!
! !REVISION HISTORY:
!  18 Dec 2013 - C. Keller - Initialization
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
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

      END FUNCTION HCO_LANDTYPE_DP 
!EOC
      END MODULE HCO_GEOTOOLS_MOD
!EOM
