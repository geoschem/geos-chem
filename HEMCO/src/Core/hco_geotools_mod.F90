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
!\\
!\\
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
  PUBLIC :: HCO_CalcVertGrid
  PUBLIC :: HCO_SetPBLm
!  PUBLIC :: HCO_CalcPBLlev

  INTERFACE HCO_LandType
     MODULE PROCEDURE HCO_LandType_Dp
     MODULE PROCEDURE HCO_LandType_Sp
  END INTERFACE HCO_LandType

  INTERFACE HCO_ValidateLon
     MODULE PROCEDURE HCO_ValidateLon_Dp
     MODULE PROCEDURE HCO_ValidateLon_Sp
  END INTERFACE HCO_ValidateLon

!  INTERFACE HCO_CalcPBLlev
!     MODULE PROCEDURE HCO_CalcPBLlev2D
!     MODULE PROCEDURE HCO_CalcPBLlev3D
!  END INTERFACE HCO_CalcPBLlev
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
! !IROUTINE: HCO_LandType_Sp
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
! !IROUTINE: HCO_LandType_Dp
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
! !IROUTINE: HCO_ValidateLon_Sp
!
! !DESCRIPTION: Subroutine HCO\_ValidateLon\_Sp ensures that the passed
! single precision longitude axis LON is steadily increasing.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_ValidateLon_Sp ( HcoState, NLON, LON, RC )
!
! !USES:
!
    USE HCO_STATE_MOD,   ONLY : HCO_STATE
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(HCO_State), POINTER       :: HcoState       ! HEMCO state object
    INTEGER,         INTENT(IN   ) :: NLON        ! # of lons
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
          CALL HCO_ERROR ( HcoState%Config%Err, '>10 iterations', RC, &
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
! !IROUTINE: HCO_ValidateLon_Dp
!
! !DESCRIPTION: Subroutine HCO\_ValidateLon\_Sp ensures that the passed
! double precision longitude axis LON is steadily increasing.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_ValidateLon_Dp ( HcoState, NLON, LON, RC )
!
! !USES:
!
    USE HCO_STATE_MOD,   ONLY : HCO_STATE
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(HCO_State), POINTER       :: HcoState       ! HEMCO state object
    INTEGER,         INTENT(IN   ) :: NLON        ! # of lons
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
          CALL HCO_ERROR ( HcoState%Config%Err, '>10 iterations', RC, &
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
! !IROUTINE: HCO_GetSUNCOS
!
! !DESCRIPTION: Subroutine HCO\_GetSUNCOS calculates the solar zenith angle
! for the given date.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_GetSUNCOS( HcoState, SUNCOS, DT, RC )
!
! !USES
!
    USE HCO_STATE_MOD,   ONLY : HCO_STATE
    USE HCO_CLOCK_MOD,   ONLY : HcoClock_Get
    USE HCO_CLOCK_MOD,   ONLY : HcoClock_GetLocal
!
! !INPUT/OUTPUT PARAMETERS:
!
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
!  22 May 2015 - C. Keller   - Initial version, based on GEOS-Chem's dao_mod.F.
!  10 Jul 2015 - R. Yantosca - Corrected issues in ProTeX header
!  02 Mar 2017 - R. Yantosca - Now compute local time as UTC + Longitude/15,
!                              so as to avoid using Voronoi TZ's for SUNCOS
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
    CALL HcoClock_Get( HcoState%Clock, cDOY=DOY, cH=HOUR, RC=RC )
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

!-----------------------------------------------------------------------------
! Prior to 3/2/17:
! Seb Eastham suggested to comment out the call to HcoClock_GetLocal.  If
! the Voronoi timezones are used, this will compute the timezones on political
! boundaries and not strictly on longitude.  This will cause funny results.
! Replace this with a strict longitudinal local time. (bmy, 3/27/17)
!         ! Local time [hours] at box (I,J) at the midpt of the chem timestep
!         CALL HcoClock_GetLocal ( HcoState, I, J, cH=LHR, RC=RC )
!
!         IF ( RC /= HCO_SUCCESS ) THEN
!            ERR = .TRUE.
!            EXIT
!         ENDIF
!-----------------------------------------------------------------------------
! Prior to 3/2/17:
! Seb Eastham says that HOUR (in the new formula below) already contains DT.
! so we need to comment this out and just use HOUR + LONGITUDE/15.
! (bmy, 3/2/17)
!         ! Adjust for time shift
!         LHR = LHR + DT
!----------------------------------------------------------------------------

         ! Compute local time as UTC + longitude/15 (bmy, 3/2/17)
         LHR = HOUR + ( HcoState%Grid%XMid%Val(I,J) / 15.0_hp )

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
       CALL HCO_ERROR ( HcoState%Config%Err, &
         'Cannot calculate SZA', RC, &
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
! !IROUTINE: HCO_GetHorzIJIndex
!
! !DESCRIPTION: Function HCO\_GetHorzIJIndex returns the grid box index for
!  the given longitude (deg E, -180...180), and latitude (deg N, -90...90).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_GetHorzIJIndex( HcoState, N, Lon, Lat, idx, jdx, RC )
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
!  10 Jul 2015 - R. Yantosca - Corrected issues in ProTeX header
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
    CALL MAPL_GetHorzIJIndex( npts=N,   II=idx,   JJ=jdx,    &
                              lon=LonR, lat=LatR, Grid=Grid, &
                              __RC__)

!!! old version of MAPL:
!    CALL MAPL_GetHorzIJIndex(N,idx,jdx,LonR,LatR,Grid=Grid,__RC__)

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
! !IROUTINE: HCO_GetHorzIJIndex
!
! !DESCRIPTION: Function HCO\_GetHorzIJIndex returns the grid box index for
!  the given longitude (deg E, -180...180), and latitude (deg N, -90...90).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_GetHorzIJIndex( HcoState, N, Lon, Lat, idx, jdx, RC )
!
! !USES:
!
    USE HCO_STATE_MOD,   ONLY : HCO_STATE
!
! !INPUT PARAMETERS:
!
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
!  10 Jul 2015 - R. Yantosca - Corrected issues in ProTeX header
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
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCO_CalcVertGrid
!
! !DESCRIPTION: Function HCO\_CalcVertGrid calculates the vertical grid
!  quantities surface pressure PSFC [Pa], surface geopotential height ZSFC
!  [m], grid box height BXHEIGHT [m], and pressure edges PEDGE [Pa]. Any of
!  these fields can be passed explicitly to the routine, in which case these
!  fields are being used. If not passed through the routine (i.e. if the
!  corresponding input argument pointer is nullified), the field is searched
!  in the HEMCO configuration file. If not found in the configuration file,
!  the field is approximated from other quantities (if possible). For example,
!  if surface pressures are provided (either passed as argument or in the
!  HEMCO configuration file as field PSFC), pressure edges are calculated
!  from PSFC and the vertical grid coordinates (Ap and Bp for a hybrid sigma
!  coordinate system). The temperature field TK [K] is needed to approximate
!  box heights and/or geopotential height (via the hydrostatic equation).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_CalcVertGrid ( HcoState, PSFC, ZSFC, TK, BXHEIGHT, PEDGE, RC )
!
! !USES
!
    USE HCO_Arr_Mod,      ONLY : HCO_ArrAssert
    USE HCO_STATE_MOD,    ONLY : HCO_STATE
    USE HCO_CALC_MOD,     ONLY : HCO_EvalFld
!
! !INPUT PARAMETERS:
!
    TYPE(HCO_State), POINTER        :: HcoState        ! HEMCO state object
    REAL(hp),        POINTER        :: PSFC(:,:)       ! surface pressure (Pa)
    REAL(hp),        POINTER        :: ZSFC(:,:)       ! surface geopotential height (m)
    REAL(hp),        POINTER        :: TK  (:,:,:)     ! air temperature (K)
    REAL(hp),        POINTER        :: BXHEIGHT(:,:,:) ! grid box height (m)
    REAL(hp),        POINTER        :: PEDGE(:,:,:)    ! pressure edges (Pa)
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,         INTENT(INOUT)  :: RC
!
! !REVISION HISTORY:
!  28 Sep 2015 - C. Keller - Initial version
!  26 Oct 2016 - R. Yantosca - Don't nullify local ptrs in declaration stmts
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER                       :: NX, NY, NZ
    INTEGER                       :: I,  J,  L
    LOGICAL                       :: Verb
    LOGICAL                       :: FoundPSFC
    LOGICAL                       :: FoundZSFC
    LOGICAL                       :: FoundTK
    LOGICAL                       :: FoundPEDGE
    LOGICAL                       :: FoundBXHEIGHT
    LOGICAL                       :: ERRBX, ERRZSFC
    REAL(hp)                      :: P1, P2
    REAL(hp), ALLOCATABLE, TARGET :: TmpTK(:,:,:)
    REAL(hp), POINTER             :: ThisTK(:,:,:)
    CHARACTER(LEN=255)            :: MSG
    CHARACTER(LEN=255)            :: LOC = 'HCO_CalcVertGrid (hco_geotools_mod.F90)'

    LOGICAL, SAVE                 :: FIRST         = .TRUE.
    LOGICAL, SAVE                 :: EVAL_PSFC     = .TRUE.
    LOGICAL, SAVE                 :: EVAL_ZSFC     = .TRUE.
    LOGICAL, SAVE                 :: EVAL_TK       = .TRUE.
    LOGICAL, SAVE                 :: EVAL_PEDGE    = .TRUE.
    LOGICAL, SAVE                 :: EVAL_BXHEIGHT = .TRUE.

    !-------------------------------
    ! HCO_CalcVertGrid begins here
    !-------------------------------

    ! Init
    Verb          = .FALSE.
    FoundPSFC     = .FALSE.
    FoundZSFC     = .FALSE.
    FoundTK       = .FALSE.
    FoundPEDGE    = .FALSE.
    FoundBXHEIGHT = .FALSE.
    ThisTK        => NULL()

    ! Verbose statements
    IF ( HcoState%amIRoot .AND. FIRST .AND. &
         HCO_IsVerb(HcoState%Config%Err,2) ) THEN
       Verb = .TRUE.
    ENDIF
    IF ( Verb ) THEN
       MSG = 'Details about vertical grid calculations (only shown on first time step):'
       CALL HCO_MSG(HcoState%Config%Err,MSG,SEP1='-')
       MSG = '1. Input data availability: '
       CALL HCO_MSG(HcoState%Config%Err,MSG,SEP1=' ')
    ENDIF

    ! ------------------------------------------------------------------
    ! TK
    ! ------------------------------------------------------------------

    ! If associated, make sure that array size is correct
    ! and pass to HEMCO surface pressure field
    IF ( ASSOCIATED(TK) ) THEN
       NX = SIZE(TK,1)
       NY = SIZE(TK,2)
       NZ = SIZE(TK,3)
       IF ( NX /= HcoState%NX .OR. NY /= HcoState%NY .OR. &
            NZ /= HcoState%NZ ) THEN
          WRITE(MSG,*) 'Wrong TK array size: ', NX, NY, NZ, &
                       '; should be: ', HcoState%NX, HcoState%NY, HcoState%NZ
          CALL HCO_ERROR( HcoState%Config%Err, MSG, RC, THISLOC=LOC )
          RETURN
       ENDIF

       ! TK is not a field in grid, so don't pass
       ThisTK  => TK
       FoundTK = .TRUE.

       ! Verbose
       IF ( Verb ) THEN
          WRITE(MSG,*) ' - Temperature field TK obtained from model interface (min,max): ', MINVAL(ThisTK), MAXVAL(ThisTK)
          CALL HCO_MSG(HcoState%Config%Err,MSG)
       ENDIF

    ELSEIF ( EVAL_TK ) THEN
       ALLOCATE(TmpTK(HcoState%NX,HcoState%NY,HcoState%NZ))
       CALL HCO_EvalFld ( HcoState, 'TK', TmpTK, RC, FOUND=FoundTK )
       IF ( RC /= HCO_SUCCESS ) RETURN
       EVAL_TK = FoundTK
       IF ( FoundTK ) ThisTK => TmpTk

       ! Verbose
       IF ( Verb ) THEN
          IF ( FoundTK ) THEN
             WRITE(MSG,*) ' - Temperature field TK [K] obtained from configuration file (min,max): ', MINVAL(ThisTK), MAXVAL(ThisTK)
             CALL HCO_MSG(HcoState%Config%Err,MSG)
          ELSE
             WRITE(MSG,*) ' - No temperature field TK found - some vertical grid calculations may not be performed...'
             CALL HCO_MSG(HcoState%Config%Err,MSG)
          ENDIF
       ENDIF
    ENDIF

    ! ------------------------------------------------------------------
    ! PSFC
    ! ------------------------------------------------------------------
    CALL HCO_ArrAssert( HcoState%Grid%PSFC, HcoState%NX, HcoState%NY, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! If associated, make sure that array size is correct
    ! and pass to HEMCO surface pressure field
    IF ( ASSOCIATED(PSFC) ) THEN
       NX = SIZE(PSFC,1)
       NY = SIZE(PSFC,2)
       IF ( NX /= HcoState%NX .OR. NY /= HcoState%NY ) THEN
          WRITE(MSG,*) 'Wrong PSFC array size: ', NX, NY, &
                       '; should be: ', HcoState%NX, HcoState%NY
          CALL HCO_ERROR( HcoState%Config%Err, MSG, RC, THISLOC=LOC )
          RETURN
       ENDIF

       ! Pass to HcoState array
       HcoState%Grid%PSFC%Val = PSFC
       FoundPSFC              = .TRUE.

       ! Verbose
       IF ( Verb ) THEN
          WRITE(MSG,*) ' - Surface pressure PSFC [Pa] obtained from model interface (min, max): ', MINVAL(HcoState%Grid%PSFC%Val), MAXVAL(HcoState%Grid%PSFC%VAL)
          CALL HCO_MSG(HcoState%Config%Err,MSG)
       ENDIF

    ! Otherwise, try to read from HEMCO configuration file
    ELSEIF ( EVAL_PSFC ) THEN
       CALL HCO_EvalFld ( HcoState, 'PSFC', HcoState%Grid%PSFC%Val, RC, FOUND=FoundPSFC )
       IF ( RC /= HCO_SUCCESS ) RETURN
       EVAL_PSFC = FoundPSFC

       ! Verbose
       IF ( Verb ) THEN
          IF ( FoundPSFC ) THEN
             WRITE(MSG,*) ' - Surface pressure PSFC [Pa] obtained from configuration file (min, max): ', MINVAL(HcoState%Grid%PSFC%Val), MAXVAL(HcoState%Grid%PSFC%VAL)
             CALL HCO_MSG(HcoState%Config%Err,MSG)
          ELSE
             MSG = ' - Surface pressure PSFC not found. Will attempt to calculate it.'
             CALL HCO_MSG(HcoState%Config%Err,MSG)
          ENDIF
       ENDIF
    ENDIF

    ! ------------------------------------------------------------------
    ! ZSFC
    ! ------------------------------------------------------------------
    CALL HCO_ArrAssert( HcoState%Grid%ZSFC, HcoState%NX, HcoState%NY, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! If associated, make sure that array size is correct
    ! and pass to HEMCO surface pressure field
    IF ( ASSOCIATED(ZSFC) ) THEN
       NX = SIZE(ZSFC,1)
       NY = SIZE(ZSFC,2)
       IF ( NX /= HcoState%NX .OR. NY /= HcoState%NY ) THEN
          WRITE(MSG,*) 'Wrong ZSFC array size: ', NX, NY, &
                       '; should be: ', HcoState%NX, HcoState%NY
          CALL HCO_ERROR( HcoState%Config%Err, MSG, RC, THISLOC=LOC )
          RETURN
       ENDIF

       ! Pass to HcoState array
       HcoState%Grid%ZSFC%Val = ZSFC
       FoundZSFC              = .TRUE.

       ! Verbose
       IF ( Verb ) THEN
          WRITE(MSG,*) ' - Surface geopotential height ZSFC [m] obtained from model interface (min, max): ', MINVAL(HcoState%Grid%ZSFC%Val), MAXVAL(HcoState%Grid%ZSFC%VAL)
          CALL HCO_MSG(HcoState%Config%Err,MSG)
       ENDIF

    ! Otherwise, try to read from HEMCO configuration file
    ELSEIF ( EVAL_ZSFC ) THEN
       CALL HCO_EvalFld ( HcoState, 'ZSFC', HcoState%Grid%ZSFC%Val, RC, FOUND=FoundZSFC )
       IF ( RC /= HCO_SUCCESS ) RETURN
       EVAL_ZSFC = FoundZSFC

       ! Verbose
       IF ( Verb ) THEN
          IF ( FoundZSFC ) THEN
             WRITE(MSG,*) ' - Surface geopotential height ZSFC [m] obtained from configuration file (min, max): ', MINVAL(HcoState%Grid%ZSFC%Val), MAXVAL(HcoState%Grid%ZSFC%VAL)
             CALL HCO_MSG(HcoState%Config%Err,MSG)
          ELSE
             MSG = ' - Surface geopotential height ZSFC not found. Will attempt to calculate it.'
             CALL HCO_MSG(HcoState%Config%Err,MSG)
          ENDIF
       ENDIF
    ENDIF

    ! ------------------------------------------------------------------
    ! PEDGE
    ! ------------------------------------------------------------------
    CALL HCO_ArrAssert( HcoState%Grid%PEDGE, HcoState%NX, &
                        HcoState%NY,         HcoState%NZ+1, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    IF ( ASSOCIATED( PEDGE ) ) THEN

       NX = SIZE(PEDGE,1)
       NY = SIZE(PEDGE,2)
       NZ = SIZE(PEDGE,3)
       IF ( NX /= HcoState%NX .OR. NY /= HcoState%NY .OR. &
            NZ /= (HcoState%NZ + 1) ) THEN
          WRITE(MSG,*) 'Wrong PEDGE array size: ', NX, NY, NZ, &
                       '; should be: ', HcoState%NX, HcoState%NY, HcoState%NZ+1
          CALL HCO_ERROR( HcoState%Config%Err, MSG, RC, THISLOC=LOC )
          RETURN
       ENDIF

       HcoState%Grid%PEDGE%Val = PEDGE
       FoundPEDGE = .TRUE.

       ! Verbose
       IF ( Verb ) THEN
          WRITE(MSG,*) ' - Pressure edges PEDGE obtained from model interface (min, max): ', MINVAL(HcoState%Grid%PEDGE%Val), MAXVAL(HcoState%Grid%PEDGE%VAL)
          CALL HCO_MSG(HcoState%Config%Err,MSG)
       ENDIF

    ELSEIF ( EVAL_PEDGE ) THEN
       CALL HCO_EvalFld ( HcoState, 'PEDGE', &
                          HcoState%Grid%PEDGE%Val, RC, FOUND=FoundPEDGE )
       IF ( RC /= HCO_SUCCESS ) RETURN
       EVAL_PEDGE = FoundPEDGE

       ! Verbose
       IF ( Verb ) THEN
          IF ( FoundPEDGE ) THEN
             WRITE(MSG,*) ' - Pressure edges PEDGE obtained from configuration file (min, max): ', MINVAL(HcoState%Grid%PEDGE%Val), MAXVAL(HcoState%Grid%PEDGE%VAL)
             CALL HCO_MSG(HcoState%Config%Err,MSG)
          ELSE
             MSG = ' - Pressure edges PEDGE not found. Will attempt to calculate it.'
             CALL HCO_MSG(HcoState%Config%Err,MSG)
          ENDIF
       ENDIF
    ENDIF

    ! ------------------------------------------------------------------
    ! BXHEIGHT
    ! ------------------------------------------------------------------
    CALL HCO_ArrAssert( HcoState%Grid%BXHEIGHT_M, HcoState%NX, &
                        HcoState%NY,              HcoState%NZ, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    IF ( ASSOCIATED( BXHEIGHT ) ) THEN

       NX = SIZE(BXHEIGHT,1)
       NY = SIZE(BXHEIGHT,2)
       NZ = SIZE(BXHEIGHT,3)
       IF ( NX /= HcoState%NX .OR. NY /= HcoState%NY .OR. &
            NZ /= HcoState%NZ ) THEN
          WRITE(MSG,*) 'Wrong BXHEIGHT array size: ', NX, NY, NZ, &
                       '; should be: ', HcoState%NX, HcoState%NY, HcoState%NZ
          CALL HCO_ERROR( HcoState%Config%Err, MSG, RC, THISLOC=LOC )
          RETURN
       ENDIF

       HcoState%Grid%BXHEIGHT_M%Val = BXHEIGHT
       FoundBXHEIGHT = .TRUE.

       ! Verbose
       IF ( Verb ) THEN
          WRITE(MSG,*) ' - Boxheights BXHEIGHT_M obtained from model interface (min, max): ', MINVAL(HcoState%Grid%BXHEIGHT_M%Val), MAXVAL(HcoState%Grid%BXHEIGHT_M%VAL)
          CALL HCO_MSG(HcoState%Config%Err,MSG)
       ENDIF

    ! Otherwise, try to read from HEMCO configuration file
    ELSEIF ( EVAL_BXHEIGHT ) THEN
       CALL HCO_EvalFld ( HcoState, 'BXHEIGHT_M', &
                          HcoState%Grid%BXHEIGHT_M%Val, RC, FOUND=FoundBXHEIGHT )
       IF ( RC /= HCO_SUCCESS ) RETURN
       EVAL_BXHEIGHT = FoundBXHEIGHT

       ! Verbose
       IF ( Verb ) THEN
          IF ( FoundBXHEIGHT ) THEN
             WRITE(MSG,*) ' - Boxheights BXHEIGHT_M obtained from configuration file (min, max): ', MINVAL(HcoState%Grid%BXHEIGHT_M%Val), MAXVAL(HcoState%Grid%BXHEIGHT_M%VAL)
             CALL HCO_MSG(HcoState%Config%Err,MSG)
          ELSE
             MSG = ' - Boxheights BXHEIGHT_M not found. Will attempt to calculate it.'
             CALL HCO_MSG(HcoState%Config%Err,MSG)
          ENDIF
       ENDIF
    ENDIF

    ! ------------------------------------------------------------------
    ! Calculate various quantities: the goal is to have the following
    ! quantities defined: ZSFC, PSFC, PEDGE, BXHEIGHT_M.
    ! - If PEDGE is not yet defined, it is calculated from surface
    !   pressure (PSFC).
    ! - If BXHEIGHT is not yet defined, it is calculated from PEDGE
    !   and TK.
    ! - If PSFC is not yet defined, it is set to the first level of
    !   PEDGE (if defined) or uniformly set to 101325 Pa.
    ! - If ZSFC is not yet defined, it is calculated from PSFC and TK.
    ! - If TK is not defined, no attempt is made to initialize it to
    !   a useful quantity. Calculations that require TK are omitted.
    ! ------------------------------------------------------------------

    ! Verbose
    IF ( Verb ) THEN
       MSG = '2. Grid calculations: '
       CALL HCO_MSG(HcoState%Config%Err,MSG,SEP1=' ')
    ENDIF

    ! Set PSFC
    IF ( .NOT. FoundPSFC ) THEN
       IF ( FoundPEDGE ) THEN
          HcoState%Grid%PSFC%Val(:,:) = HcoState%Grid%PEDGE%Val(:,:,1)

          ! Verbose
          IF ( Verb ) THEN
             MSG = ' - Surface pressure set to surface pressure edge.'
             CALL HCO_MSG(HcoState%Config%Err,MSG)
          ENDIF
       ELSE
          HcoState%Grid%PSFC%Val(:,:) = 101325.0_hp
          IF ( FIRST .AND. HcoState%amIRoot ) THEN
             MSG = 'Surface pressure PSFC uniformly set to 101325 Pa! ' // &
                   'This may affect the accuracy of vertical grid '     // &
                   'quantities. It is recommended you provide PSFC via '// &
                   'the model-HEMCO interface or the HEMCO configuration file!'
             CALL HCO_WARNING( HcoState%Config%Err,MSG, RC, THISLOC=LOC, WARNLEV=1 )
          ENDIF

          ! Verbose
          IF ( Verb ) THEN
             MSG = ' - Surface pressure uniformly set to 101325.0 Pa.'
             CALL HCO_MSG(HcoState%Config%Err,MSG)
          ENDIF
       ENDIF
       FoundPSFC = .TRUE.
    ENDIF

    ! Set PEDGE
    IF ( .NOT. FoundPEDGE ) THEN
!$OMP PARALLEL DO                                                      &
!$OMP DEFAULT( SHARED )                                                &
!$OMP PRIVATE( I, J, L )                                               &
!$OMP SCHEDULE( DYNAMIC )
       DO L = 1, HcoState%NZ+1
       DO J = 1, HcoState%NY
       DO I = 1, HcoState%NX
          HcoState%Grid%PEDGE%Val(I,J,L) &
           = HcoState%Grid%zGrid%AP(L)   &
           + ( HcoState%Grid%zGrid%BP(L) &
             * HcoState%Grid%PSFC%Val(I,J) )
       ENDDO
       ENDDO
       ENDDO
!$OMP END PARALLEL DO
       FoundPEDGE = .TRUE.

       ! Verbose
       IF ( Verb ) THEN
          WRITE(MSG,*) ' - PEDGE calculated from PSFC, Ap, and Bp (min, max): ', &
             MINVAL(HcoState%Grid%PEDGE%Val), MAXVAL(HcoState%Grid%PEDGE%Val)
          CALL HCO_MSG(HcoState%Config%Err,MSG)
       ENDIF
    ENDIF

    ! Set surface height and/or grid box height
    IF ( .NOT. FoundZSFC .OR. .NOT. FoundBXHEIGHT ) THEN
       IF ( FoundTK .AND. FoundPEDGE ) THEN

          ! Initialize error flags
          ERRZSFC = .FALSE.
          ERRBX   = .FALSE.

!$OMP PARALLEL DO                                                      &
!$OMP DEFAULT( SHARED )                                                &
!$OMP PRIVATE( I, J, L, P1, P2 )                                       &
!$OMP SCHEDULE( DYNAMIC )
          DO L = 1, HcoState%NZ
          DO J = 1, HcoState%NY
          DO I = 1, HcoState%NX

             ! BOXHEIGHT (hydrostatic equation)
             IF ( .NOT. FoundBXHEIGHT ) THEN
                P1 = HcoState%Grid%PEDGE%Val(I,J,L)
                P2 = HcoState%Grid%PEDGE%Val(I,J,L+1)
                IF ( P2 == 0.0_hp ) THEN
                   ERRBX = .TRUE.
                ELSE
                   HcoState%Grid%BXHEIGHT_M%Val(I,J,L) = HcoState%Phys%Rdg0 &
                                                       * ThisTK(I,J,1)      &
                                                       * LOG( P1 / P2 )
                ENDIF
             ENDIF

             ! ZSFC
             IF ( L == 1 .AND. .NOT. FoundZSFC ) THEN
                P1 = 101325.0_hp
                P2 = HcoState%Grid%PEDGE%Val(I,J,1)
                IF ( P2 == 0.0_hp ) THEN
                   ERRZSFC = .TRUE.
                ELSE
                   HcoState%Grid%ZSFC%Val(I,J) = HcoState%Phys%Rdg0 &
                                               * ThisTK(I,J,1)      &
                                               * LOG( P1 / P2 )
                ENDIF
             ENDIF
          ENDDO
          ENDDO
          ENDDO
!$OMP END PARALLEL DO

          IF ( ERRZSFC ) THEN
             MSG = 'Cannot calculate surface geopotential heights - at least one ' // &
                   'surface pressure value is zero! You can either provide an '    // &
                   'updated pressure edge field (PEDGE) or add a field with the '  // &
                   'surface geopotential height to your configuration file (ZSFC)'
             CALL HCO_ERROR( HcoState%Config%Err, MSG, RC, THISLOC=LOC )
             RETURN
          ELSE
             FoundZSFC = .TRUE.

             ! Verbose
             IF ( Verb ) THEN
                WRITE(MSG,*) ' - ZSFC calculated from PSFC and T (min, max): ', &
                   MINVAL(HcoState%Grid%ZSFC%Val), MAXVAL(HcoState%Grid%ZSFC%Val)
                CALL HCO_MSG(HcoState%Config%Err,MSG)
             ENDIF
          ENDIF

          IF ( ERRBX ) THEN
             MSG = 'Cannot calculate grid box heights - at least one ' // &
                   'pressure value is zero! You can either provide an '    // &
                   'updated pressure edge field (PEDGE) or add a field with the '  // &
                   'grid box heights to your configuration file (BOXHEIGHT_M)'
             CALL HCO_ERROR( HcoState%Config%Err, MSG, RC, THISLOC=LOC )
             RETURN
          ELSE
             FoundZSFC = .TRUE.

             ! Verbose
             IF ( Verb ) THEN
                WRITE(MSG,*) ' - Boxheights calculated from PEDGE and T (min, max): ', &
                   MINVAL(HcoState%Grid%BXHEIGHT_M%Val), MAXVAL(HcoState%Grid%BXHEIGHT_M%Val)
                CALL HCO_MSG(HcoState%Config%Err,MSG)
             ENDIF
          ENDIF

       ! PEDGE and/or TK not defined
       ELSE
          IF ( .NOT. FoundZSFC .AND. FIRST .AND. HcoState%amIRoot ) THEN
             MSG = 'Cannot set surface height ZSFC. This may cause '      // &
                   'some extensions to fail. HEMCO tries to calculate '   // &
                   'ZSFC from surface pressure and air temperature, but ' // &
                   'at least one of these variables seem to be missing.'
             CALL HCO_WARNING( HcoState%Config%Err,MSG, RC, THISLOC=LOC, WARNLEV=1 )
          ENDIF
          IF ( .NOT. FoundBXHEIGHT .AND. FIRST .AND. HcoState%amIRoot ) THEN
             MSG = 'Cannot set boxheights BXHEIGHT_M. This may cause '      // &
                   'some extensions to fail. HEMCO tries to calculate '     // &
                   'BXHEIGHT from pressure edges and air temperature, but ' // &
                   'at least one of these variables seem to be missing.'
             CALL HCO_WARNING( HcoState%Config%Err,MSG, RC, THISLOC=LOC, WARNLEV=1 )
          ENDIF
       ENDIF
    ENDIF

    ! ------------------------------------------------------------------
    ! Wrap up and leave
    ! ------------------------------------------------------------------

    ! Verbose
    IF ( Verb ) THEN
       WRITE(MSG,*) 'Vertical grid calculations done.'
       CALL HCO_MSG(HcoState%Config%Err,MSG,SEP2='-')
    ENDIF

    ! Cleanup
    ThisTK => NULL()
    IF ( ALLOCATED(TmpTK) ) DEALLOCATE(TmpTK)

    ! Update first flag
    FIRST = .FALSE.

    ! Return w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE HCO_CalcVertGrid
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCO_SetPBLm
!
! !DESCRIPTION: Subroutine HCO\_SetPBLm sets the HEMCO PBL mixing height in
! meters. It first tries to read it from field 'FldName' (from the HEMCO data
! list), then to fill it from field 'PBLM', and then assigns the default value
! 'DefVal' to it.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_SetPBLm( HcoState, FldName, PBLM, DefVal, RC )
!
! !USES
!
    USE HCO_Arr_Mod,      ONLY : HCO_ArrAssert
    USE HCO_STATE_MOD,    ONLY : HCO_STATE
    USE HCO_CALC_MOD,     ONLY : HCO_EvalFld
!
! !INPUT PARAMETERS:
!
    TYPE(HCO_State),  POINTER                  :: HcoState    ! HEMCO state object
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN   )  :: FldName     ! field name
    REAL(hp),         OPTIONAL, POINTER        :: PBLM(:,:)   ! pbl mixing height
    REAL(hp),         OPTIONAL, INTENT(IN   )  :: DefVal      ! default value
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(INOUT)  :: RC
!
! !REVISION HISTORY:
!  28 Sep 2015 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER                       :: NX, NY
    LOGICAL                       :: FOUND
    CHARACTER(LEN=255)            :: MSG
    CHARACTER(LEN=255)            :: LOC = 'HCO_SetPBLm (hco_geotools_mod.F90)'

    !-------------------------------
    ! HCO_SetPBLm begins here
    !-------------------------------

    ! Init
    FOUND = .FALSE.

    ! Try to read from file first
    IF ( PRESENT( FldName ) ) THEN
       CALL HCO_EvalFld ( HcoState, FldName, &
          HcoState%Grid%PBLHEIGHT%Val, RC, FOUND=FOUND )
       IF ( RC /= HCO_SUCCESS ) RETURN

       ! Verbose
       IF ( HcoState%amIRoot .AND. HCO_IsVerb(HcoState%Config%Err,2) ) THEN
          IF ( FOUND ) THEN
             WRITE(MSG,*) 'HEMCO PBL heights obtained from field ',TRIM(FldName)
             CALL HCO_MSG(HcoState%Config%Err,MSG,SEP2='-')
          ENDIF
       ENDIF
    ENDIF

    ! Pass 2D field if available
    IF ( .not. FOUND ) THEN
       IF ( PRESENT( PBLM ) ) THEN
          IF ( ASSOCIATED(PBLM) ) THEN
             NX = SIZE(PBLM,1)
             NY = SIZE(PBLM,2)
             IF ( NX /= HcoState%NX .OR. NY /= HcoState%NY ) THEN
                WRITE(MSG,*) 'Wrong PBLM array size: ', NX, NY, &
                             '; should be: ', HcoState%NX, HcoState%NY
                CALL HCO_ERROR( HcoState%Config%Err, MSG, RC, THISLOC=LOC )
                RETURN
             ENDIF

             ! Make sure size is ok
             CALL HCO_ArrAssert( HcoState%Grid%PBLHEIGHT, &
                                 HcoState%NX, HcoState%NY, RC )
             IF ( RC /= HCO_SUCCESS ) RETURN

             ! Pass data
             HcoState%Grid%PBLHEIGHT%Val = PBLM
             FOUND                       = .TRUE.

             ! Verbose
             IF ( HcoState%amIRoot .AND. HCO_IsVerb(HcoState%Config%Err,2) ) THEN
                WRITE(MSG,*) 'HEMCO PBL heights obtained from provided 2D field.'
                CALL HCO_MSG(HcoState%Config%Err,MSG,SEP2='-')
             ENDIF
          ENDIF
       ENDIF
    ENDIF

    ! Finally, assign default value if field not yet set
!    IF ( .NOT. FOUND .AND. PRESENT(DefVal) ) THEN
    IF ( .NOT. FOUND ) THEN
       IF ( PRESENT(DefVal) ) THEN
          ! Make sure size is ok
          CALL HCO_ArrAssert( HcoState%Grid%PBLHEIGHT, &
                              HcoState%NX, HcoState%NY, RC )
          IF ( RC /= HCO_SUCCESS ) RETURN

          ! Pass data
          HcoState%Grid%PBLHEIGHT%Val = DefVal
          FOUND                       = .TRUE.

          ! Verbose
          IF ( HcoState%amIRoot .AND. HCO_IsVerb(HcoState%Config%Err,2) ) THEN
             WRITE(MSG,*) 'HEMCO PBL heights uniformly set to ', DefVal
             CALL HCO_MSG(HcoState%Config%Err,MSG,SEP2='-')
          ENDIF
       ENDIF
    ENDIF

    ! Error check
    IF ( .NOT. FOUND ) THEN
       WRITE(MSG,*) 'Cannot set PBL height: a valid HEMCO data field, ', &
          'an explicit 2D field or a default value must be provided!'
       CALL HCO_ERROR( HcoState%Config%Err, MSG, RC, THISLOC=LOC )
       RETURN
    ENDIF

    ! Return w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE HCO_SetPBLm
!EOC
!!------------------------------------------------------------------------------
!!                  Harvard-NASA Emissions Component (HEMCO)                   !
!!------------------------------------------------------------------------------
!!BOP
!!
!! !SUBROUTINE: HCO_CalcPBLlev3D
!!
!! !DESCRIPTION:
!!\\
!!\\
!! !INTERFACE:
!!
!  SUBROUTINE HCO_CalcPBLlev3D ( HcoState, PBLFRAC, RC )
!!
!! !USES
!!
!    USE HCO_Arr_Mod,      ONLY : HCO_ArrAssert
!    USE HCO_STATE_MOD,    ONLY : HCO_STATE
!!
!! !INPUT PARAMETERS:
!!
!    TYPE(HCO_State), POINTER        :: HcoState        ! HEMCO state object
!    REAL(hp),        POINTER        :: PBLFRAC(:,:,:)  ! planetary PBL fraction
!!
!! !INPUT/OUTPUT PARAMETERS:
!!
!    INTEGER,         INTENT(INOUT)  :: RC
!!
!! !REVISION HISTORY:
!!  05 May 2016 - C. Keller - Initial version
!!EOP
!!------------------------------------------------------------------------------
!!BOC
!!
!! !LOCAL VARIABLES:
!!
!    INTEGER             :: I,  J,  L
!    CHARACTER(LEN=255)  :: MSG
!    CHARACTER(LEN=255)  :: LOC = 'HCO_CalcPBLlev3D (hco_geotools_mod.F90)'
!
!    !-------------------------------
!    ! HCO_CalcPBLlev3D begins here
!    !-------------------------------
!
!    ! Check input array size
!    IF ( SIZE(PBLFRAC,1) /= HcoState%NX .OR. &
!         SIZE(PBLFRAC,2) /= HcoState%NY .OR. &
!         SIZE(PBLFRAC,3) /= HcoState%NZ        ) THEN
!       WRITE(MSG,*) 'Input array PBLFRAC has wrong horiz. dimensions: ', &
!                     SIZE(PBLFRAC,1),SIZE(PBLFRAC,2),SIZE(PBLFRAC,3)
!       CALL HCO_ERROR( HcoState%Config%Err, MSG, RC, THISLOC=LOC )
!       RETURN
!    ENDIF
!
!    ! Make sure array is associated
!    CALL HCO_ArrAssert( HcoState%Grid%PBL, HcoState%NX, HcoState%NY, RC )
!    IF ( RC /= HCO_SUCCESS ) RETURN
!
!    ! Initialize values
!    HcoState%Grid%PBL%Val = 1.0
!
!!$OMP PARALLEL DO                                                      &
!!$OMP DEFAULT( SHARED )                                                &
!!$OMP PRIVATE( I, J, L )                                               &
!!$OMP SCHEDULE( DYNAMIC )
!    DO J = 1, HcoState%NY
!    DO I = 1, HcoState%NX
!       ! Search for first level where PBL fraction is zero
!       DO L = 1, HcoState%NZ
!          IF ( PBLFRAC(I,J,L) > 0.0_hp ) CYCLE
!          HcoState%Grid%PBL%Val(I,J) = MAX(L-1,1)
!          EXIT
!       ENDDO
!    ENDDO
!    ENDDO
!!$OMP END PARALLEL DO
!
!    ! Return w/ success
!    RC = HCO_SUCCESS
!
!  END SUBROUTINE HCO_CalcPBLlev3D
!!EOC
!!------------------------------------------------------------------------------
!!                  Harvard-NASA Emissions Component (HEMCO)                   !
!!------------------------------------------------------------------------------
!!BOP
!!
!! !SUBROUTINE: HCO_CalcPBLlev2D
!!
!! !DESCRIPTION:
!!\\
!!\\
!! !INTERFACE:
!!
!  SUBROUTINE HCO_CalcPBLlev2D ( HcoState, PBLlev, RC )
!!
!! !USES
!!
!    USE HCO_Arr_Mod,      ONLY : HCO_ArrAssert
!    USE HCO_STATE_MOD,    ONLY : HCO_STATE
!!
!! !INPUT PARAMETERS:
!!
!    TYPE(HCO_State), POINTER        :: HcoState        ! HEMCO state object
!    INTEGER,         POINTER        :: PBLlev(:,:)     ! planetary PBL lev
!!
!! !INPUT/OUTPUT PARAMETERS:
!!
!    INTEGER,         INTENT(INOUT)  :: RC
!!
!! !REVISION HISTORY:
!!  05 May 2016 - C. Keller - Initial version
!!EOP
!!------------------------------------------------------------------------------
!!BOC
!!
!! !LOCAL VARIABLES:
!!
!    INTEGER             :: I,  J,  L
!    CHARACTER(LEN=255)  :: MSG
!    CHARACTER(LEN=255)  :: LOC = 'HCO_CalcPBLlev2D (hco_geotools_mod.F90)'
!
!    !-------------------------------
!    ! HCO_CalcPBLlev2D begins here
!    !-------------------------------
!
!    ! Make sure array is associated
!    CALL HCO_ArrAssert( HcoState%Grid%PBL, HcoState%NX, HcoState%NY, RC )
!    IF ( RC /= HCO_SUCCESS ) RETURN
!
!    ! Check input array size
!    IF ( SIZE(PBLlev,1) /= HcoState%NX .OR. SIZE(PBLlev,2) /= HcoState%NY ) THEN
!       WRITE(MSG,*) 'Input array PBLlev has wrong horiz. dimensions: ', &
!              SIZE(PBLlev,1),SIZE(PBLlev,2)
!       CALL HCO_ERROR( HcoState%Config%Err, MSG, RC, THISLOC=LOC )
!       RETURN
!    ENDIF
!
!    ! Set values from PBLlev
!!$OMP PARALLEL DO                                                      &
!!$OMP DEFAULT( SHARED )                                                &
!!$OMP PRIVATE( I, J )                                                  &
!!$OMP SCHEDULE( DYNAMIC )
!    DO J = 1, HcoState%NY
!    DO I = 1, HcoState%NX
!       HcoState%Grid%PBL%Val(I,J) = MAX(PBLlev(I,J),1)
!    ENDDO
!    ENDDO
!!$OMP END PARALLEL DO
!
!    ! Return w/ success
!    RC = HCO_SUCCESS
!
!  END SUBROUTINE HCO_CalcPBLlev2D
!!EOC
END MODULE HCO_GeoTools_Mod
!EOM
