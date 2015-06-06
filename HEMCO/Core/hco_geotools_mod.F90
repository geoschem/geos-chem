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
  PUBLIC :: HCO_ValidateLon
  PUBLIC :: HCO_GetSUNCOS
  PUBLIC :: HCO_GetHorzIJIndex

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
!  16 Jul 2014 - C. Keller   - Added HCO_ValidateLon
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
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !SUBROUTINE: HCO_GetSUNCOS
!
! !DESCRIPTION: Subroutine HCO\_GetSUNCOS calculates the solar zenith angle
! for the given date.
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_GetSUNCOS( am_I_Root, HcoState, SUNCOS, DT, RC )
!
! !USES
!
    USE HCO_STATE_MOD,   ONLY : HCO_STATE
    USE HCO_CLOCK_MOD,   ONLY : HcoClock_Get
    USE HCO_CLOCK_MOD,   ONLY : HcoClock_GetLocal
!
! !INPUT/OUTPUT PARAMETERS:
!
    LOGICAL,         INTENT(IN   )  :: am_I_Root      ! Root CPU? 
    TYPE(HCO_State), POINTER        :: HcoState       ! HEMCO state object
    INTEGER,         INTENT(IN   )  :: DT             ! Time shift relative to current date [hrs]
!
! !OUTPUT PARAMETERS:
!
    REAL(hp),        INTENT(  OUT)  :: SUNCOS(HcoState%NX,HcoState%NY)
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,         INTENT(INOUT)  :: RC             ! Return code
!
! !REVISION HISTORY:
!  22 May 2015 - C. Keller - Initial version, based on GEOS-Chem's dao_mod.F.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER              :: I, J,   DOY, HOUR
    LOGICAL              :: ERR
    REAL(hp)             :: YMID_R, S_YMID_R,  C_YMID_R
    REAL(hp)             :: R,      DEC
    REAL(hp)             :: S_DEC,  C_DEC 
    REAL(hp)             :: SC,     LHR
    REAL(hp)             :: AHR

    ! Coefficients for solar declination angle
    REAL(hp),  PARAMETER :: A0 = 0.006918e+0_hp
    REAL(hp),  PARAMETER :: A1 = 0.399912e+0_hp
    REAL(hp),  PARAMETER :: A2 = 0.006758e+0_hp
    REAL(hp),  PARAMETER :: A3 = 0.002697e+0_hp
    REAL(hp),  PARAMETER :: B1 = 0.070257e+0_hp
    REAL(hp),  PARAMETER :: B2 = 0.000907e+0_hp
    REAL(hp),  PARAMETER :: B3 = 0.000148e+0_hp

    !-------------------------------
    ! HCO_GetSUNCOS starts here! 
    !-------------------------------

    ! Get current time information
    CALL HcoClock_Get( cDOY=DOY, cH=HOUR, RC=RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Add time adjustment 
    HOUR = HOUR + DT

    ! Make sure HOUR is within valid range (0-24)
    IF ( HOUR < 0 ) THEN
       HOUR = HOUR + 24
       DOY  = DOY  - 1
    ELSEIF ( HOUR > 23 ) THEN
       HOUR = HOUR - 24
       DOY  = DOY  + 1
    ENDIF

    ! Make sure DOY is within valid range of 1 to 365
    DOY = MAX(MIN(DOY,365),1)

    ! Path length of earth's orbit traversed since Jan 1 [radians]
    R        = ( 2e+0_hp * HcoState%Phys%PI / 365e+0_hp ) * DBLE( DOY - 1 )

    ! Solar declination angle (low precision formula) [radians]
    DEC      = A0 - A1*COS(         R ) + B1*SIN(         R ) &
                  - A2*COS( 2e+0_hp*R ) + B2*SIN( 2e+0_hp*R ) &
                  - A3*COS( 3e+0_hp*R ) + B3*SIN( 3e+0_hp*R )

    ! Pre-compute sin & cos of DEC outside of DO loops (for efficiency)
    S_DEC    = SIN( DEC )
    C_DEC    = COS( DEC )

    ! Init
    ERR = .FALSE.

    !=================================================================
    ! Compute cosine of solar zenith angle
    !=================================================================
!$OMP PARALLEL DO                                         &
!$OMP DEFAULT( SHARED )                                   &
!$OMP PRIVATE( I,      J,   YMID_R, S_YMID_R,  C_YMID_R ) &
!$OMP PRIVATE( LHR,    AHR, SC,     RC                  )
    DO J = 1, HcoState%NY 
    DO I = 1, HcoState%NX

         ! Latitude of grid box [radians]
         YMID_R     = HcoState%Grid%YMID%Val(I,J) * HcoState%Phys%PI_180 

         ! Pre-compute sin & cos of DEC outside of I loop (for efficiency)
         S_YMID_R   = SIN( YMID_R )
         C_YMID_R   = COS( YMID_R )

         !==============================================================
         ! Compute cosine of SZA at the midpoint of the chem timestep
         ! Required for photolysis, chemistry, emissions, drydep
         !==============================================================

         ! Local time [hours] at box (I,J) at the midpt of the chem timestep
         CALL HcoClock_GetLocal ( HcoState, I, J, cH=LHR, RC=RC )
         IF ( RC /= HCO_SUCCESS ) THEN
            ERR = .TRUE.
            EXIT
         ENDIF

         ! Adjust for time shift
         LHR = LHR + DT
         IF ( LHR <   0.0_hp ) LHR = LHR + 24.0_hp
         IF ( LHR >= 24.0_hp ) LHR = LHR - 24.0_hp

         ! Hour angle at box (I,J) [radians]
         AHR = ABS( LHR - 12.0_hp ) * 15.0_hp * HcoState%Phys%PI_180
         
         ! Corresponding cosine( SZA ) at box (I,J) [unitless]
         SC = ( S_YMID_R * S_DEC              ) &
            + ( C_YMID_R * C_DEC * COS( AHR ) )

         ! COS(SZA) at the current time
         SUNCOS(I,J) = SC
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

    ! Check error status
    IF ( ERR ) THEN 
       CALL HCO_ERROR ( 'Cannot calculate SZA', RC, &
          THISLOC='HCO_GetSUNCOS (hco_geotools_mod.F90)' )
       RETURN
    ENDIF

    ! Leave w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE HCO_GetSUNCOS
!EOC
#if defined(ESMF_)
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !SUBROUTINE: HCO_GetHorzIJIndex
!
! !DESCRIPTION: Function HCO\_GetHorzIJIndex returns the grid box index for 
!  the given longitude (deg E, -180...180), and latitude (deg N, -90...90).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_GetHorzIJIndex( am_I_Root, HcoState, N, Lon, Lat, idx, jdx, RC ) 
!
! !USES
!
#include "MAPL_Generic.h"
    USE ESMF
    USE MAPL_Mod
    USE HCO_STATE_MOD,   ONLY : HCO_STATE
!
! !INPUT PARAMETERS:
!
    LOGICAL,         INTENT(IN   )  :: am_I_Root      ! Root CPU? 
    TYPE(HCO_State), POINTER        :: HcoState       ! HEMCO state object
    INTEGER,         INTENT(IN   )  :: N           
    REAL(hp),        INTENT(IN   )  :: Lon(N)
    REAL(hp),        INTENT(IN   )  :: Lat(N)
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,         INTENT(INOUT)  :: RC
!
! !OUTPUT PARAMETERS:
!
    INTEGER,         INTENT(  OUT)  :: IDX(N), JDX(N)
! 
! !REVISION HISTORY:
!  04 Jun 2015 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    REAL             :: LonR(N), LatR(N)
    REAL, PARAMETER  :: radToDeg = 57.2957795
    TYPE(ESMF_Grid)  :: Grid

    ! Defined Iam and STATUS
    __Iam__("HCO_GetHorzIJIndex (hco_geotools_mod.F90)")

    !-------------------------------
    ! HCO_GetHorzIJIndex begins here
    !-------------------------------

    ! Get grid
    ASSERT_(ASSOCIATED(HcoState%GridComp))
    CALL ESMF_GridCompGet( HcoState%GridComp, Grid=Grid, __RC__ ) 

    ! Shadow variables
    LonR(:) = Lon / radToDeg 
    LatR(:) = Lat / radToDeg 

    ! Get indeces
    CALL MAPL_GetHorzIJIndex(LonR,LatR,N,Grid,idx,jdx,__RC__)

    ! Return w/ success
    RC =  HCO_SUCCESS

  END SUBROUTINE HCO_GetHorzIJIndex 
!EOC
#else
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !SUBROUTINE: HCO_GetHorzIJIndex
!
! !DESCRIPTION: Function HCO\_GetHorzIJIndex returns the grid box index for 
!  the given longitude (deg E, -180...180), and latitude (deg N, -90...90).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_GetHorzIJIndex( am_I_Root, HcoState, N, Lon, Lat, idx, jdx, RC ) 
!
! !USES
!
    USE HCO_STATE_MOD,   ONLY : HCO_STATE
!
! !INPUT PARAMETERS:
!
    LOGICAL,         INTENT(IN   )  :: am_I_Root      ! Root CPU? 
    TYPE(HCO_State), POINTER        :: HcoState       ! HEMCO state object
    INTEGER,         INTENT(IN   )  :: N           
    REAL(hp),        INTENT(IN   )  :: Lon(N)
    REAL(hp),        INTENT(IN   )  :: Lat(N)
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,         INTENT(INOUT)  :: RC
!
! !OUTPUT PARAMETERS:
!
    INTEGER,         INTENT(  OUT)  :: IDX(N), JDX(N)
! 
! !REVISION HISTORY:
!  04 Jun 2015 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER  :: I, J, L, FOUND 
    REAL(hp) :: iLon1, iLon2
    REAL(hp) :: iLat1, iLat2
    REAL(hp) :: delta

    !-------------------------------
    ! HCO_GetHorzIJIndex begins here
    !-------------------------------

    ! Initialize
    IDX(:) = -1
    JDX(:) = -1
    FOUND  =  0

    ! do for every grid box 
    DO J = 1, HcoState%NY 
    DO I = 1, HcoState%NX 

       ! Get grid edges for this box
 
       ! Longitude edges 
       IF ( ASSOCIATED(HcoState%Grid%XEDGE%Val) ) THEN 
          iLon1 = HcoState%Grid%XEDGE%Val(I,  J)
          iLon2 = HcoState%Grid%XEDGE%Val(I+1,J)

       ! Approximate from mid points
       ELSE
          iLon1 = HcoState%Grid%XMID%Val(I,J)
          IF ( I < HcoState%NX ) THEN
             delta = HcoState%Grid%XMID%Val(I+1,J)-iLon1
          ELSE
             delta = iLon1 - HcoState%Grid%XMID%Val(I-1,J)
          ENDIF
          iLon2 = iLon1 + (delta / 2.0_hp)
          iLon1 = iLon1 - (delta / 2.0_hp)
       ENDIF

       ! Latitude edges
       IF ( ASSOCIATED(HcoState%Grid%YEDGE%Val) ) THEN 
          iLat1 = HcoState%Grid%YEDGE%Val(I,J)
          iLat2 = HcoState%Grid%YEDGE%Val(I,J+1)

       ! Approximate from mid points
       ELSE
          iLat1 = HcoState%Grid%YMID%Val(I,J)
          IF ( J < HcoState%NY ) THEN
             delta = HcoState%Grid%YMID%Val(I,J+1)-iLat1
          ELSE
             delta = iLat1 - HcoState%Grid%YMID%Val(I,J-1)
          ENDIF
          iLat2 = iLat1 + (delta / 2.0_hp)
          iLat1 = iLat1 - (delta / 2.0_hp)
       ENDIF

       ! Check if it's within this box
       DO L = 1, N
          IF ( IDX(L) > 0 ) CYCLE 

          IF ( Lon(L) >= HcoState%Grid%XEDGE%Val(I,  J  ) .AND. &
               Lon(L) <= HcoState%Grid%XEDGE%Val(I+1,J  ) .AND. &
               Lat(L) >= HcoState%Grid%YEDGE%Val(I  ,J  ) .AND. &
               Lat(L) <= HcoState%Grid%YEDGE%Val(I  ,J+1)        ) THEN
             IDX(L) = I
             JDX(L) = J
             FOUND  = FOUND + 1
             IF ( FOUND == N ) EXIT
          ENDIF
       ENDDO 
    ENDDO 
    ENDDO 

    ! Return w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE HCO_GetHorzIJIndex 
!EOC
#endif
END MODULE HCO_GeoTools_Mod
!EOM
