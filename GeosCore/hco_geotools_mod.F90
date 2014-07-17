!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: hco_geotools_mod.F90
!
! !DESCRIPTION: Module HCO\_GeoTools\_Mod contains a collection of 
! helper routines for extracting geographical information. 
! \\
! Note that some of these routines are based upon GEOS-5 data and may 
! need to be revised for other met. fields! 
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
  PUBLIC :: HCO_ValidateLon

  INTERFACE HCO_LandType
     MODULE PROCEDURE HCO_LandType_Dp
     MODULE PROCEDURE HCO_LandType_Sp
  END INTERFACE HCO_LandType

  INTERFACE HCO_ValidateLon
     MODULE PROCEDURE HCO_ValidateLon_Dp
     MODULE PROCEDURE HCO_ValidateLon_Sp
  END INTERFACE HCO_ValidateLon
!
! !PRIVATE MEMBER FUNCTIONS:
!
  PRIVATE:: HCO_LandType_Dp
  PRIVATE:: HCO_LandType_Sp
  PRIVATE:: HCO_ValidateLon_Dp
  PRIVATE:: HCO_ValidateLon_Sp
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
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !FUNCTION: HCO_ValidateLon_Sp
!
! !DESCRIPTION: Subroutine HCO\_ValidateLon\_Sp ensures that the passed 
! single precision longitude axis LON is steadily increasing.
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_ValidateLon_Sp ( NLON, LON, RC )
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,  INTENT(IN   ) :: NLON        ! # of lons
!
! !INPUT/OUTPUT PARAMETERS:
!
    REAL(sp), INTENT(INOUT) :: LON(NLON)   ! longitude axis
    INTEGER,  INTENT(INOUT) :: RC          ! Return code
!
! !REVISION HISTORY:
!  16 Jul 2014 - C. Keller - Initialization
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER             :: I, CNT
    LOGICAL             :: REDO
    INTEGER, PARAMETER  :: MAXIT = 10

    !--------------------------
    ! HCO_ValidateLon_Sp begins here
    !--------------------------

    REDO = .TRUE.
    CNT  = 0

    DO WHILE ( REDO )

       ! Exit w/ error after 10 iterations
       CNT = CNT + 1
       IF ( CNT > MAXIT ) THEN
          CALL HCO_ERROR ( '>10 iterations', RC, &
                           THISLOC='HCO_ValidateLon (HCO_GEOTOOLS_MOD.F90)' )
          RETURN
       ENDIF

       DO I = 1, NLON

          ! If we reach the last grid box, all values are steadily
          ! increasing (otherwise, the loop would have been exited).
          IF ( I == NLON ) THEN
             REDO = .FALSE.
             EXIT
          ENDIF

          ! Check if next lon value is lower, in which case we subtract
          ! a value of 360 (degrees) from all longitude values up to
          ! this point. Then repeat the lookup (from the beginning).
          IF ( LON(I+1) < LON(I) ) THEN
             LON(1:I) = LON(1:I) - 360.0_sp
             EXIT
          ENDIF

       ENDDO !I
    ENDDO ! REDO

    ! Leave w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE HCO_ValidateLon_Sp
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !FUNCTION: HCO_ValidateLon_Dp
!
! !DESCRIPTION: Subroutine HCO\_ValidateLon\_Sp ensures that the passed 
! double precision longitude axis LON is steadily increasing.
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_ValidateLon_Dp ( NLON, LON, RC )
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,  INTENT(IN   ) :: NLON        ! # of lons
!
!
! !INPUT/OUTPUT PARAMETERS:
!
    REAL(dp), INTENT(INOUT) :: LON(NLON)   ! longitude axis
    INTEGER,  INTENT(INOUT) :: RC          ! Return code
!
! !REVISION HISTORY:
!  16 Jul 2014 - C. Keller - Initialization
!EOP
!------------------------------------------------------------------------------
!BOC
    INTEGER             :: I, CNT
    LOGICAL             :: REDO
    INTEGER, PARAMETER  :: MAXIT = 10

    !--------------------------
    ! HCO_ValidateLon_Dp begins here
    !--------------------------

    REDO = .TRUE.
    CNT  = 0

    DO WHILE ( REDO )

       ! Exit w/ error after 10 iterations
       CNT = CNT + 1
       IF ( CNT > MAXIT ) THEN
          CALL HCO_ERROR ( '>10 iterations', RC, &
                           THISLOC='HCO_ValidateLon (HCO_GEOTOOLS_MOD.F90)' )
          RETURN
       ENDIF

       DO I = 1, NLON

          ! If we reach the last grid box, all values are steadily
          ! increasing (otherwise, the loop would have been exited).
          IF ( I == NLON ) THEN
             REDO = .FALSE.
             EXIT
          ENDIF

          ! Check if next lon value is lower, in which case we subtract
          ! a value of 360 (degrees) from all longitude values up to
          ! this point. Then repeat the lookup (from the beginning).
          IF ( LON(I+1) < LON(I) ) THEN
             LON(1:I) = LON(1:I) - 360.0_dp
             EXIT
          ENDIF

       ENDDO !I
    ENDDO ! REDO

    ! Leave w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE HCO_ValidateLon_Dp
!EOC
END MODULE HCO_GeoTools_Mod
!EOM
