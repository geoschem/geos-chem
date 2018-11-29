!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: gc_grid_mod.F90
!
! !DESCRIPTION: Module GC\_GRID\_MOD contains variables and routines which are 
!  used to specify the parameters of a GEOS-Chem horizontal grid. Grid 
!  parameters are computed as 3D arrays, which are required for interfacing
!  with a GCM.
!\\  
!\\
! !INTERFACE: 
!
MODULE GC_Grid_Mod
! 
! !USES:
!
  USE Error_Mod                        ! Error-handling routines
  USE Precision_Mod                    ! For GEOS-Chem Precision (fp)
  USE PhysConstants                    ! Physical constants
  USE Registry_Mod, ONLY : MetaRegItem

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC  :: Cleanup_Grid
  PUBLIC  :: Compute_Grid
  PUBLIC  :: DoGridComputation
  PUBLIC  :: Get_Area_m2
  PUBLIC  :: Get_Area_cm2
  PUBLIC  :: Get_Bounding_Box
  PUBLIC  :: Get_xEdge
  PUBLIC  :: Get_xMid
  PUBLIC  :: Get_yEdge
  PUBLIC  :: Get_yEdge_r
  PUBLIC  :: Get_yMid
  PUBLIC  :: Get_yMid_r
  PUBLIC  :: Get_yMid_r_w
  PUBLIC  :: Get_ySin
  PUBLIC  :: Get_xOffSet
  PUBLIC  :: Get_yOffSet
  PUBLIC  :: Init_Grid
  PUBLIC  :: Its_A_Nested_Grid
  PUBLIC  :: Set_xOffSet
  PUBLIC  :: Set_yOffSet
  PUBLIC  :: SetGridFromCtr
#if defined ( MODEL_WRF )
  PUBLIC  :: SetGridFromCtrEdges
#endif
  PUBLIC  :: GET_IJ

! Make some arrays public
  PUBLIC  :: XMID, YMID, XEDGE, YEDGE, YSIN, AREA_M2
!
! !REVISION HISTORY:
!  23 Feb 2012 - R. Yantosca - Initial version, based on grid_mod.F
!  01 Mar 2012 - R. Yantosca - Validated for nested grids
!  03 Apr 2012 - M. Payer    - Added ySin for map_a2a regrid (M. Cooper)
!  04 Dec 2012 - R. Yantosca - Modified for GIGC running in ESMF environment
!  26 Feb 2013 - R. Yantosca - Fixed bug in computation of lons & lats when
!                              connecting GEOS-Chem to the GEOS-5 GCM
!  19 May 2013 - C. Keller   - Added wrapper routine DoGridComputation so that
!                              module can also be used by HEMCO.
!  02 Dec 2014 - M. Yannetti - Added PRECISION_MOD
!  26 Mar 2015 - R. Yantosca - Removed obsolete, commented-out code
!  29 Nov 2016 - R. Yantosca - Renamed to gc_grid_mod.F90 to avoid namespace
!                              conflicts when interfacing GCHP to the BCC model
!  09 Aug 2017 - R. Yantosca - Now register lon (slice of XMID), lat (slice of
!                              YMID) and AREA_M2.  Added registry routines etc.
!  18 Aug 2017 - R. Yantosca - Move roundoff routines to roundoff_mod.F90
!  23 Aug 2017 - R. Yantosca - Registry is now moved to grid_registry_mod.F90
!  11 Nov 2018 - H.P. Lin    - Added SetGridFromCtrEdges for running in WRF environment
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE TYPES:
!
  ! Scalars
  LOGICAL                       :: IS_NESTED        ! Nested grid (T/F)?
  INTEGER                       :: I0               ! Nested grid X-offset
  INTEGER                       :: J0               ! Nested-grid Y-ofset

  ! Arrays
  REAL(fp), ALLOCATABLE, TARGET :: XMID     (:,:,:) ! Lon centers [degrees]
  REAL(fp), ALLOCATABLE, TARGET :: XEDGE    (:,:,:) ! Lon edges   [degrees]
  REAL(fp), ALLOCATABLE, TARGET :: YMID     (:,:,:) ! Lat centers [degrees]
  REAL(fp), ALLOCATABLE, TARGET :: YEDGE    (:,:,:) ! Lat edges   [degrees]
  REAL(fp), ALLOCATABLE, TARGET :: YMID_R   (:,:,:) ! Lat centers [radians]
  REAL(fp), ALLOCATABLE, TARGET :: YEDGE_R  (:,:,:) ! Lat edges   [radians]
  REAL(fp), ALLOCATABLE, TARGET :: YSIN     (:,:,:) ! SIN( lat edges )
  REAL(fp), ALLOCATABLE, TARGET :: YMID_R_W (:,:,:) ! Lat ctrs  nest grid [rad]
  REAL(fp), ALLOCATABLE, TARGET :: YEDGE_R_W(:,:,:) ! Lat edges nest grid [rad]
  REAL(fp), ALLOCATABLE, TARGET :: AREA_M2  (:,:,:) ! Grid box areas [m2]
  
CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Compute_Grid
!
! !DESCRIPTION: Subroutine COMPUTE\_GRID is the wrapper routine to 
! initializes the longitude, latitude and surface area arrays. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Compute_Grid( am_I_Root,                          &
                           I1, I2, J1,   J2,   JSP,  JNP,      &
                           L1, L2, DLON, DLAT, I_LO, J_LO, RC )
!
! !USES:
!
    USE ErrCode_Mod
!
! !INPUT PARAMETERS:
!
    LOGICAL,  INTENT(IN)  :: am_I_Root                     ! Root CPU?

    ! Variables with local CPU indices
    INTEGER,  INTENT(IN)  :: I1,  I2                       ! Min lon index
    INTEGER,  INTENT(IN)  :: J1,  J2                       ! Local lat indices
    INTEGER,  INTENT(IN)  :: JSP, JNP                      ! Polar lat indices
    INTEGER,  INTENT(IN)  :: L1,  L2                       ! Local lev indices
    REAL(fp), INTENT(IN)  :: DLON(I2-I1+1,J2-J1+1,L2-L1+1) ! Delta lon [deg]
    REAL(fp), INTENT(IN)  :: DLAT(I2-I1+1,J2-J1+1,L2-L1+1) ! Delta lat [deg]

    ! Variables with global CPU indices
    INTEGER,  INTENT(IN)  :: I_LO                          ! Min global lon
    INTEGER,  INTENT(IN)  :: J_LO                          ! Min global lat
!
! !OUTPUT PARAMETERS:
!  
    INTEGER,  INTENT(OUT) :: RC                            ! Success/failure?
!
! !REMARKS:
!  (1) Lon/lat loop indices IG, JG are global indices.
!  (2) Lon/lat loop indices I,  J  are local to each CPU.
!  (3) We do not need to have global loop indices for vertical levels,
!       because we will always decompose the grid for MPI parallelization
!       in longitude and/or latitude.  All vertical levels must be present
!       on each CPU for the grid-independent GEOS-Chem to function properly.
!
! !REVISION HISTORY:
!  19 May 2014 - C. Keller   - Initial version: now wrapper routine that
!                              calls DoGridComputation.
!EOP
!------------------------------------------------------------------------------
!BOC

  !======================================================================
  ! Compute_Grid begins here!
  !======================================================================

  ! Assume success
  RC = GC_SUCCESS

  ! Compute the GEOS-Chem grid specifications
  Call DoGridComputation( am_I_Root,                                     &
                          I1,      I2,       J1,        J2,      JSP,    &
                          JNP,     L1,       L2,        DLON,    DLAT,   & 
                          I_LO,    J_LO,     I0,        J0,      XMID,   &
                          XEDGE,   YMID,     YEDGE,     YSIN,    YMID_R, & 
                          YEDGE_R, YMID_R_W, YEDGE_R_W, AREA_M2, RC       )

  END SUBROUTINE COMPUTE_GRID
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: DoGridComputation 
!
! !DESCRIPTION: Subroutine DoGridComputation initializes the longitude, 
!  latitude and surface area arrays. This used to be subroutine COMPUTE\_GRID.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE DoGridComputation( am_I_Root,                               &
                                I1,   I2,    J1,    J2,    JSP,   JNP,   &
                                L1,   L2,    DLON,  DLAT,  I_LO,  J_LO,  &
                                IOFF, JOFF,  XMD,   XDG,   YMD,   YDG,   & 
                                YSN,  YMDR,  YDGR,  YMDRW, YDGRW, AM2, RC )
!
! !USES:
!
    USE ErrCode_Mod
!
! !INPUT PARAMETERS:
!
    LOGICAL, INTENT(IN)  :: am_I_Root                      ! Root CPU?

    ! Variables with local CPU indices
    INTEGER,  INTENT(IN)  :: I1,  I2                       ! Min lon index
    INTEGER,  INTENT(IN)  :: J1,  J2                       ! Local lat indices
    INTEGER,  INTENT(IN)  :: JSP, JNP                      ! Polar lat indices
    INTEGER,  INTENT(IN)  :: L1,  L2                       ! Local lev indices
    REAL(fp), INTENT(IN)  :: DLON(I2-I1+1,J2-J1+1,L2-L1+1) ! Delta lon [deg]
    REAL(fp), INTENT(IN)  :: DLAT(I2-I1+1,J2-J1+1,L2-L1+1) ! Delta lat [deg]

    ! Variables with global CPU indices
    INTEGER,  INTENT(IN)  :: I_LO                          ! Min global lon
    INTEGER,  INTENT(IN)  :: J_LO                          ! Min global lat

    ! Offsets (for nested grids) 
    INTEGER,  INTENT(IN)  :: IOFF
    INTEGER,  INTENT(IN)  :: JOFF
!
! !OUTPUT PARAMETERS:
! 
    REAL(fp), INTENT(OUT) :: XMD  (:,:,:)                  ! Lon centers [deg]
    REAL(fp), INTENT(OUT) :: XDG  (:,:,:)                  ! Lon edges [deg]
    REAL(fp), INTENT(OUT) :: YMD  (:,:,:)                  ! Lat centers [deg]
    REAL(fp), INTENT(OUT) :: YDG  (:,:,:)                  ! Lat edges [deg]
    REAL(fp), INTENT(OUT) :: YSN  (:,:,:)                  ! SIN( lat edges )
    REAL(fp), INTENT(OUT) :: YMDR (:,:,:)                  ! Lat centers [rad]
    REAL(fp), INTENT(OUT) :: YDGR (:,:,:)                  ! Lat edges [rad]
    REAL(fp), INTENT(OUT) :: YMDRW(:,:,:)                  ! window lat centers
    REAL(fp), INTENT(OUT) :: YDGRW(:,:,:)                  ! Window lat edes
    REAL(fp), INTENT(OUT) :: AM2  (:,:,:)                  ! Area [m2]
!
! !OUTPUT PARAMETERS:
!  
    INTEGER, INTENT(OUT) :: RC                             ! Success or failure?
!
! !REMARKS:
!  (1) Lon/lat loop indices IG, JG are global indices.
!  (2) Lon/lat loop indices I,  J  are local to each CPU.
!  (3) We do not need to have global loop indices for vertical levels,
!       because we will always decompose the grid for MPI parallelization
!       in longitude and/or latitude.  All vertical levels must be present
!       on each CPU for the grid-independent GEOS-Chem to function properly.
!
! !REVISION HISTORY:
!  23 Feb 2012 - R. Yantosca - Initial version, based on grid_mod.F
!  30 Jul 2012 - R. Yantosca - Now accept am_I_Root as an argument when
!                              running with the traditional driver main.F
!  03 Dec 2012 - R. Yantosca - Add RC to argument list
!  04 Dec 2012 - R. Yantosca - Now define arrays with local CPU lon/lat bounds
!  07 Dec 2012 - R. Yantosca - Bug fix: make sure the last longitude edge is
!                              computed properly.  Test for IG==I2, not I==I2.
!  07 Dec 2012 - R. Yantosca - Also do not apply half-polar boxes when running
!                              in ESMF environment
!  26 Feb 2013 - R. Yantosca - Bug fix: now compute IND_X and IND_Y properly
!                              when connecting GEOS-Chem to the GEOS-5 GCM
!  21 Mar 2013 - R. Yantosca - Add fix to prevent zero surface area at poles
!  21 Mar 2013 - R. Yantosca - Rename loop indices to prevent confusion
!  06 Jun 2013 - M. Payer    - Add fix to compute sine of last latitude edge
!                              for MAP_A2A regridding (C. Keller)
!  19 May 2014 - C. Keller   - Renamed from Compute_grid to DoGridComputation.
!  06 Nov 2014 - C. Keller   - Now use LBOUND to get leftmost index of YMDRW.
!  26 Mar 2015 - R. Yantosca - Fix apparent optimization error by using 
!                              scalars in call to the SIN function
!  26 Mar 2015 - R. Yantosca - Cosmetic changes; improve indentation
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
! 
    ! Scalars
    INTEGER  :: I,     J,       L,        IG,        JG,      LBND
    REAL(fp) :: SIN_N, SIN_S,   SIN_DIFF, YEDGE_VAL, YSIN_VAL

    ! Arrays
    REAL(fp) :: IND_X( I1:I2+1 )
    REAL(fp) :: IND_Y( J1:J2+1 )

    !======================================================================
    ! Initialization
    !======================================================================

    ! Index array for longitudes
    DO I = I1, I2+1
       IND_X(I) = ( ( I + IOFF - 1 ) * 1e+0_fp ) + ( I_LO - 1 )
    ENDDO

    ! Index array for latitudes
    DO J = J1, J2+1
       IND_Y(J) = ( ( J + JOFF - 1 ) * 1e+0_fp ) + ( J_LO - 1 )
    ENDDO

    ! Left bound of YMDRW. Since we now pass YMID_R_W as an argument to 
    ! this routine, we cannot know for sure that the 2nd subscript starts
    ! at index 0 (ckeller, 11/06/14).
    LBND = LBOUND(YMDRW,2)

    !======================================================================
    ! Compute longitude and latitude arrays
    !======================================================================
    
    ! We can set L=1 instead of looping over vertical levels
    L = 1
       
    !----------------------------------------------------------------------
    ! Longitude center and edge arrays
    !----------------------------------------------------------------------

    ! Loop over local latitudes
    DO J = J1, J2

       ! Loop over local longitudes
       DO I = I1, I2

          ! Longitude centers
          XMD(I,J,L)  = ( DLON(I,J,L) * IND_X(I) ) - 180e+0_fp
          
          ! Longitude edges
          XDG(I,J,L) = XMD(I,J,L) - ( DLON(I,J,L) * 0.5e+0_fp )

          ! Compute the last longitude edge
          IF ( I == I2 ) THEN
             XDG(I+1,J,L) = XDG(I,J,L) + DLON(I,J,L)
          ENDIF
             
       ENDDO
    ENDDO

    !----------------------------------------------------------------------
    ! Latitude center and edge arrays
    !-------------------------------------------------------------------===

    ! Loop over local latitudes
    DO J = J1, J2

       ! Global latitude index
       JG = J + ( J_LO - 1 )

       ! Loop over local longitudes
       DO I = I1, I2

          !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          !%%%%%%   LATITUDE CENTERS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if defined( ESMF_ ) || defined( EXTERNAL_GRID ) || defined( EXTERNAL_FORCING )
          !----------------------------------------------------------------
          !         %%%%%%% GEOS-Chem HP (with ESMF & MPI) %%%%%%%
          !
          ! Do not define half-sized polar boxes (bmy, 3/21/13)
          !----------------------------------------------------------------
#else
          !----------------------------------------------------------------
          !         %%%%%%% GEOS-Chem CLASSIC (with OpenMP) %%%%%%%
          !
          ! Current practice in the standard GEOS-Chem is to make
          ! the polar grid boxes (for non-GCAP global grids only) be
          ! half the size of other grid boxes.  This lets us make +90
          ! degrees and -90 degrees be the edges of the grid.
          ! (bmy, 7/2/13)
          !----------------------------------------------------------------

          ! Lat centers (degrees)
          YMD(I,J,L)     = ( DLAT(I,J,L) * IND_Y(J) ) - 90e+0_fp

          IF ( JG == JSP ) THEN
             YMD(I,J,L)  = -90e+0_fp + ( 0.5e+0_fp * DLAT(I,J,L) )   ! S pole
          ELSE IF ( JG == JNP ) THEN
             YMD(I,J,L)  = +90e+0_fp - ( 0.5e+0_fp * DLAT(I,J,L) )   ! N pole
          ENDIF
#endif

          ! Lat centers (radians)
          YMDR(I,J,L)   = ( PI_180 * YMD(I,J,L)  )

#if defined( NESTED_CH ) || defined( NESTED_EU ) || defined( NESTED_NA ) || defined( NESTED_AS )

          !----------------------------------------------------------------
          !              %%%%% FOR NESTED GRIDS ONLY %%%%%
          !----------------------------------------------------------------

          ! Lat centers (radians), for nested grid window array
          YMDRW(I,LBND+J,L) = YMDR(I,J,L)

          ! Compute YMID_R_W at edges of nested region
          IF ( J == J1 ) THEN
             YMDRW(I,LBND+J-1,1) = YMDR(I,J,L) - ( DLAT(I,J,L) * PI_180 )
          ELSE IF ( J == J2 ) THEN
             YMDRW(I,LBND+J2+1,1) = YMDR(I,J,L) + ( DLAT(I,J,L) * PI_180 )
          ENDIF
                    
#endif

          !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          !%%%%%%   LATITUDE EDGES   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

          ! Lat edges (degrees and radians)
          YDG(I,J,L)    = YMD(I,J,L) - ( DLAT(I,J,L) * 0.5e+0_fp )
            
#if defined( ESMF_ ) || defined( EXTERNAL_GRID ) || defined( EXTERNAL_FORCING )
          !----------------------------------------------------------------
          !         %%%%%%% GEOS-Chem HP (with ESMF & MPI) %%%%%%%
    
          ! Do not define half-sized polar boxes (bmy, 3/21/13)
          !----------------------------------------------------------------
#else
          !----------------------------------------------------------------
          !         %%%%%%% GEOS-Chem CLASSIC (with OpenMP) %%%%%%%
          !
          ! Current practice in the standard GEOS-Chem is to force
          ! the northern edge of grid boxes along the SOUTH POLE to
          ! be -90 degrees latitude. (bmy, 3/21/13)
          !----------------------------------------------------------------
          IF ( JG == JSP ) THEN
             YDG(I,J,L) = -90e+0_fp                             ! S pole
          ENDIF
#endif          

          ! Lat edges (radians)
          YDGR(I,J,L)   = ( PI_180  * YDG(I,J,L) )
          
          ! mjc - Compute sine of latitude edges (needed for map_a2a regrid)
          YEDGE_VAL     = YDGR(I,J,L)           ! Lat edge in radians
          YSIN_VAL      = SIN( YEDGE_VAL )      ! SIN( lat edge )
          YSN(I,J,L)    = YSIN_VAL              ! Store in YSN array

       ENDDO
    ENDDO

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%   LATITUDE EDGES (the last edge)   %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
    ! Test for North Pole
    IF ( J2 == JNP ) THEN
          
       ! North pole case (global grids only)
       DO I = I1, I2
#if defined( ESMF_ ) || defined( EXTERNAL_GRID ) || defined( EXTERNAL_FORCING )
          !----------------------------------------------------------------
          !         %%%%%%% GEOS-Chem HP (with ESMF & MPI) %%%%%%%
          !
          ! Do not define half-sized polar boxes (bmy, 3/21/13)
          !----------------------------------------------------------------
          YDG  (I,J2+1,L)  = YDG(I,J2,L)   + DLAT(I,J2,L)
#else
          !----------------------------------------------------------------
          !         %%%%%%% GEOS-Chem CLASSIC (with OpenMP) %%%%%%%
          !
          ! Current practice in the standard GEOS-Chem is to force
          ! the northern edge of grid boxes along the NORTH POLE to
          ! be +90 degrees latitude. (bmy, 3/21/13)
          !----------------------------------------------------------------
          YDG (I,J2+1,L)   = +90e+0_fp
#endif
          YDGR(I,J2+1,L)   = YDG(I,J2+1,L) * PI_180

          ! Also compute sine of last latitude edge! (ckeller, 02/13/12)
          YEDGE_VAL        = YDGR(I,J2+1,L)     ! Lat edge in radians
          YSIN_VAL         = SIN( YEDGE_VAL )   ! SIN( lat edge )
          YSN(I,J2+1,L)    = YSIN_VAL           ! Store in YSN array

       ENDDO
          
    ELSE
          
       !-------------------------------------------------------------------
       !                %%%%% FOR NESTED GRIDS ONLY %%%%%
       !-------------------------------------------------------------------

       ! No north pole (nested grids only)
       DO I = I1, I2
          YDG (I,J2+1,L)  = YDG(I,J2,L  ) + DLAT(I,J2,L)
          YDGR(I,J2+1,L)  = YDG(I,J2+1,L) * PI_180

          ! Also compute sine of last latitude edge! (ckeller, 02/13/12)
          YEDGE_VAL       = YDGR(I,J2+1,L)
          YSIN_VAL        = SIN( YEDGE_VAL )
          YSN(I,J2+1,L)   = YSIN_VAL

       ENDDO
    ENDIF

    !======================================================================
    ! Compute grid box surface areas (algorithm from old "input.f")
    !
    ! The surface area of a grid box is derived as follows:
    ! 
    !    Area = dx * dy
    !
    ! Where:
    !
    !    dx is the arc length of the box in longitude
    !    dy is the arc length of the box in latitude
    !  
    ! Which are computed as:
    !  
    !    dx = r * delta-longitude
    !       = ( Re * cos[ YMID[J] ] ) * ( 2 * PI / IIIPAR )
    !
    !    dy = r * delta-latitude
    !       = Re * ( YEDGE[J+1] - YEDGE[J] )
    !  
    ! Where:
    !    
    !    Re         is the radius of the earth
    !    YMID[J]    is the latitude at the center of box J
    !    YEDGE[J+1] is the latitude at the N. Edge of box J
    !    YEDGE[J]   is the latitude at the S. Edge of box J
    !
    ! So, the surface area is thus:
    ! 
    !    Area = ( Re * cos( YMID[J] ) * ( 2 * PI / IIIPAR ) *
    !             Re * ( YEDGE[J+1] - YEDGE[J] )
    !
    !    2*PI*Re^2    {                                            }      
    ! = ----------- * { cos( YMID[J] ) * ( YEDGE[J+1] - YEDGE[J] ) }
    !     IIIPAR      {                                            }
    !
    ! And, by using the trigonometric identity:
    !
    !    d sin(x) = cos x * dx
    !
    ! The following term:
    !
    !    cos( YMID[J] ) * ( YEDGE[J+1] - YEDGE[J] ) 
    !
    ! May also be written as a difference of sines:
    !
    !    sin( YEDGE[J+1] ) - sin( YEDGE[J] ) 
    ! 
    ! So the final formula for surface area of a grid box is:
    ! 
    !            2*PI*Re^2    {                                     }
    !    Area = ----------- * { sin( YEDGE[J+1] ) - sin( YEDGE[J] ) }
    !              IIIPAR     {                                     }
    !
    !
    ! NOTES:
    ! (1) The formula with sines is more numerically stable, and will 
    !      yield identical global total surface areas for all grids.
    ! (2) The units are determined by the radius of the earth Re.
    !      if you use Re [m], then surface area will be in [m2], or
    !      if you use Re [cm], then surface area will be in [cm2], etc.
    ! (3) The grid box surface areas only depend on latitude, as they
    !      are symmetric in longitude.  To compute the global surface
    !      area, multiply the surface area arrays below by the number
    !      of longitudes (e.g. IIIPAR).
    ! (4) At present, assumes that GEOS-Chem will work on a
    !      Cartesian grid.
    !
    ! (bmy, 4/20/06, 2/24/12)
    !====================================================================== 

    ! Loop over local latitudes
    DO J = J1, J2

       ! Global latitude index
       JG  = J + ( J_LO - 1 )
          
       ! Loop over local longitudes
       DO I = I1, I2

          ! Sine of latitudes at N and S edges of grid box (I,J,L)
          SIN_N       = SIN( YDGR(I,J+1,L) )
          SIN_S       = SIN( YDGR(I,J,  L) )

          ! Difference of sin(latitude) at N and S edges of grid box
          SIN_DIFF    = SIN_N - SIN_S

#if defined( ESMF_ ) || defined( EXTERNAL_GRID ) || defined( EXTERNAL_FORCING )
          !-------------------------------------------------------------
          !      %%%%%%% GEOS-Chem HP (with ESMF & MPI) %%%%%%%
          !
          ! The GEOS-5 GCM is a grid-point model, with the polar
          ! boxes (+90 & -90 degrees latitude) being grid box
          ! centers.  But since GEOS-Chem also needs to know the
          ! latitudes at the north & south edges of each grid box,
          ! we have to make a kludge.
          !
          ! For all grid boxes along the NORTH POLE (+90 lat), we
          ! let the northern edge be greater than 90 degrees.  For
          ! example, if the grid spacing is 1 degree in latitude,
          ! then for all longitudes at the North Pole:
          !
          !   * The N. EDGE of each grid box is +90.5 degrees
          !   * The CENTER  of each grid box is +90   degrees
          !   * The S. EDGE of each grid box is +89.5 degrees.
          !
          ! Similarly, at for all grid boxes along the SOUTH POLE,
          ! we have this condition (also assuming a grid spacing
          ! of 1 degree in latitude):
          !
          !   * The N. EDGE of each grid box is -89.5 degrees
          !   * The CENTER  of each grid box is -90   degrees
          !   * The S. EDGE of each grid box is -90.5 degrees.
          !
          ! Therefore, at the poles, the latitudes at the northern
          ! and southern edges of each grid box are symmetric around
          ! either +90 degrees or -90 degrees.  When you take the
          ! difference of the sine of the latitudes at the north and
          ! south edges of a polar grid box, the terms will cancel
          ! each other out, resulting in a grid box surface area
          ! that is zero.
          !
          ! We can take advantage of this symmetry around +90 and
          ! -90 degrees to make a simple fix:
          !
          ! (1) AT THE SOUTH POLE: Subtract the sine of the latitude
          !     at the north edge of the grid box from the sine
          !     of -90 degrees (which is -1) and then multiply by 2.
          !
          ! (2) AT THE NORTH POLE: Subtract the sine of +90 degrees
          !     (which is 1) from the sine of the latitude at the
          !     south edge of the grid box, and then multiply by 2.
          !
          ! This fix avoids having polar grid boxes with zero area.
          !    -- Bob Yantosca (21 Mar 2013)
          !--------------------------------------------------------------

          ! South pole kludge
          IF ( JG == JSP ) THEN
             SIN_DIFF = 2e+0_fp * ( SIN_N - ( -1e+0_fp ) )
          ENDIF

          ! North pole kludge
          IF ( JG == JNP ) THEN
             SIN_DIFF = 2e+0_fp * ( 1e+0_fp - SIN_S )
          ENDIF
#endif

          ! Grid box surface areas [m2]
          AM2(I,J,L) = ( DLON(I,J,L) * PI_180 ) * ( Re**2 ) * SIN_DIFF

       ENDDO
    ENDDO

#if defined( ESMF_ ) || defined( EXTERNAL_GRID ) || defined( EXTERNAL_FORCING )
    WRITE( 6, '(''Call incorrectly made for lat-lon grid comp'')' )
    WRITE( 6, '(''Should be using external grids only'')' )
    RC = GC_FAILURE
#else
    ! Return successfully
    RC = GC_SUCCESS
#endif

    !======================================================================
    ! Echo info to stdout
    !======================================================================
    IF ( am_I_Root ) THEN
       WRITE( 6, '(''Nested-Grid X-offset [boxes]:'', i4 )' ) IOFF
       WRITE( 6, '(''Nested-Grid Y-offset [boxes]:'', i4 )' ) JOFF
       WRITE( 6, '(a)' )
       WRITE( 6, '(''Grid box longitude centers [degrees]: '')' )
       WRITE( 6, '(8(f8.3,1x))' ) ( XMD(I,1,1),  I=1,I2-I1+1 )
       WRITE( 6, '(a)' )
       WRITE( 6, '(''Grid box longitude edges [degrees]: '')' )
       WRITE( 6, '(8(f8.3,1x))' ) ( XDG(I,1,1), I=1,I2+1-I1+1 )
       WRITE( 6, '(a)' )
       WRITE( 6, '(''Grid box latitude centers [degrees]: '')' )
       WRITE( 6, '(8(f8.3,1x))' ) ( YMD(1,J,1),  J=1,J2-J1+1 )
       WRITE( 6, '(a)' )
       WRITE( 6, '(''Grid box latitude edges [degrees]: '')' )
       WRITE( 6, '(8(f8.3,1x))' ) ( YDG(1,J,1), J=1,J2+1-J1+1 )
       WRITE( 6, '(a)' )
       WRITE( 6, '(''SIN( grid box latitude edges )'')' )
       WRITE( 6, '(8(f8.3,1x))' ) ( YSN(1,J,1), J=1,J2+1-J1+1 )
    ENDIF

  END SUBROUTINE DoGridComputation 
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: SetGridFromCtr 
!
! !DESCRIPTION: Subroutine SetGridFromCtr sets the grid based upon the passed
! mid-points. This routine is primarily intented to provide an interface to 
! GEOS-5 in an ESMF-environment.
!\\
!\\
! This routine does not update the grid box areas (AREA\_M2) of grid\_mod.F90.
! These need to be updated manually. We cannot do this within this routine
! since in GEOS-5, the grid box areas are not yet available during the 
! initialization phase (they are imported from superdynamics).
! !INTERFACE:
!
  SUBROUTINE SetGridFromCtr( am_I_Root, NX, NY, lonCtr, latCtr, RC ) 
!
! USES
!
    USE ErrCode_Mod
    USE Roundoff_Mod
!
! !INPUT PARAMETERS: 
!
    LOGICAL,  INTENT(IN)   :: am_I_Root      ! Root CPU?
    INTEGER,  INTENT(IN)   :: NX             ! # of lons
    INTEGER,  INTENT(IN)   :: NY             ! # of lats
    REAL(f4), INTENT(IN)   :: lonCtr(NX,NY)  ! Lon ctrs [rad]
    REAL(f4), INTENT(IN)   :: latCtr(NX,NY)  ! Lat ctrs [rad]
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER, INTENT(INOUT) :: RC 
!
! !REVISION HISTORY:
!  02 Jan 2014 - C. Keller   - Initial version
!  26 Mar 2015 - R. Yantosca - Fix apparent optimization error by using 
!                              scalars in call to the SIN function
!  03 Sep 2015 - C. Keller   - Bug fix: need to explicitly calculate adjacent
!                              mid-point to calculate first xedge and yedge. 
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER            :: I,         J,        L
    INTEGER            :: NI,        NJ,       NL
    REAL(fp)           :: YEDGE_VAL, YSIN_VAL, TMP

    ! Strings
    CHARACTER(LEN=255) :: ErrMsg
    CHARACTER(LEN=255) :: ThisLoc

    !======================================================================
    ! SetGridFromCtr begins here! 
    !======================================================================

    ! Initialize
    RC      = GC_SUCCESS
    ErrMsg  = ''
    ThisLoc = ' -> at SetGridFromCtr (in module GeosUtil/gc_grid_mod.F90)'

    ! Get array size
    NI = SIZE(XMID,1)
    NJ = SIZE(XMID,2)
    NL = SIZE(XMID,3)

    ! Horizontal dimensions must agree
    IF ( NX /= NI .OR. NY /= NJ ) THEN
       WRITE(ErrMsg,*) 'Grid dimension mismatch: ',NX, '/=',NI,              &
                       ' and/or ',NY,'/=',NJ
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    ! Loop over all grid boxes
    DO L = 1, NL 
    DO J = 1, NJ
    DO I = 1, NI

       ! Mid points: get directly from passed value
       XMID(I,J,L)      = RoundOff( lonCtr(I,J) / PI_180, 4 )
       YMID(I,J,L)      = RoundOff( latCtr(I,J) / PI_180, 4 )
       YMID_R(I,J,L)    = YMID(I,J,L) * PI_180 

       IF ( ALLOCATED(YMID_R_W) ) THEN
          YMID_R_W(I,J,L)  = YMID_R(I,J,L)
       ENDIF

       ! Edges: approximate from neighboring mid points.
       IF ( I == 1 ) THEN
          TMP           = RoundOff( lonCtr(I+1,J) / PI_180, 4 )
          XEDGE(I,J,L)  = XMID(I,J,L) - ( ( TMP - XMID(I,J,L) ) / 2.0_f4 )
       ELSE
          XEDGE(I,J,L) = ( XMID(I,J,L) + XMID(I-1,J,L) ) / 2.0_f4
       ENDIF

       IF ( J == 1 ) THEN
          TMP          = RoundOff( latCtr(I,J+1) / PI_180, 4 )
          YEDGE(I,J,L) = YMID(I,J,L) - ( ( TMP - YMID(I,J,L) ) / 2.0_f4 )
       ELSE
          YEDGE(I,J,L) = ( YMID(I,J,L) + YMID(I,J-1,L) ) / 2.0_f4
       ENDIF

       ! Special treatment at uppermost edge
       IF ( I == NI ) THEN
          XEDGE(I+1,J,L) = XMID(I,J,L) + ( ( XMID(I,J,L) - XMID(I-1,J,L) ) / 2.0_f4 )
       ENDIF
       IF ( J == NJ ) THEN
          YEDGE(I,J+1,L) = YMID(I,J,L) + ( ( YMID(I,J,L) - YMID(I,J-1,L) ) / 2.0_f4 )
       ENDIF

       ! Special quantities directly derived from YEDGE
       YEDGE_R(I,J,L)   = YEDGE(I,J,L) * PI_180
       YEDGE_VAL        = YEDGE_R(I,J,L)           ! Lat edge [radians]
       YSIN_VAL         = SIN( YEDGE_VAL)          ! SIN( lat edge )
       YSIN(I,J,L)      = YSIN_VAL                 ! Store in YSIN array

    ENDDO
    ENDDO
    ENDDO 

    ! Return w/ success
    RC = GC_SUCCESS

  END SUBROUTINE SetGridFromCtr
!EOC
#if defined ( MODEL_WRF )
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: SetGridFromCtrEdges
!
! !DESCRIPTION: Subroutine SetGridFromCtrEdges sets the grid based upon the passed
! mid-points and edge-points given an external grid. This interface is primarily
! used for GEOS-Chem to interface with the WRF model.
!\\
!\\
! This routine does not update the grid box areas (AREA\_M2) of grid\_mod.F90.
! These need to be updated manually from State_Met%AREA_M2 to maintain consistency
! with the GEOS-Chem interface to GEOS-5.
! !INTERFACE:
!
  SUBROUTINE SetGridFromCtrEdges( am_I_Root, NX, NY, lonCtr, latCtr, lonEdge, latEdge, RC )
!
! USES
!
    USE ErrCode_Mod
    USE Roundoff_Mod
!
! !INPUT PARAMETERS:
!
    LOGICAL,  INTENT(IN)   :: am_I_Root           ! Root CPU?
    INTEGER,  INTENT(IN)   :: NX                  ! # of lons
    INTEGER,  INTENT(IN)   :: NY                  ! # of lats
    REAL(f4), INTENT(IN)   :: lonCtr (NX,NY)      ! Lon ctrs [rad]
    REAL(f4), INTENT(IN)   :: latCtr (NX,NY)      ! Lat ctrs [rad]
    REAL(f4), INTENT(IN)   :: lonEdge(NX+1,NY)    ! Lon edges [rad]
    REAL(f4), INTENT(IN)   :: latEdge(NX,NY+1)    ! Lat edges [rad]
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER, INTENT(INOUT) :: RC
!
! !REVISION HISTORY:
!  11 Nov 2018 - H.P. Lin    - Initial version based on SetGridFromCtr
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER            :: I,         J,        L
    INTEGER            :: NI,        NJ,       NL
    REAL(fp)           :: YEDGE_VAL, YSIN_VAL

    ! Strings
    CHARACTER(LEN=255) :: ErrMsg
    CHARACTER(LEN=255) :: ThisLoc

    !======================================================================
    ! SetGridFromCtrEdges begins here!
    !======================================================================

    ! Initialize
    RC      = GC_SUCCESS
    ErrMsg  = ''
    ThisLoc = ' -> at SetGridFromCtrEdges (in module GeosUtil/gc_grid_mod.F90)'

    ! Get array size
    NI = SIZE(XMID,1)
    NJ = SIZE(XMID,2)
    NL = SIZE(XMID,3)

    ! Horizontal dimensions must agree
    IF ( NX /= NI .OR. NY /= NJ ) THEN
       WRITE(ErrMsg,*) 'Grid dimension mismatch: ',NX, '/=',NI,              &
                       ' and/or ',NY,'/=',NJ
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    ! Loop over all grid boxes
    DO L = 1, NL
    DO J = 1, NJ
    DO I = 1, NI

       ! Mid points: get directly from passed value
       XMID(I,J,L)      = RoundOff( lonCtr(I,J) / PI_180, 4 )
       YMID(I,J,L)      = RoundOff( latCtr(I,J) / PI_180, 4 )
       YMID_R(I,J,L)    = YMID(I,J,L) * PI_180

       IF ( ALLOCATED(YMID_R_W) ) THEN
          YMID_R_W(I,J,L)  = YMID_R(I,J,L)
       ENDIF

       ! Edges: get directly from passed value
       XEDGE(I,J,L)     = RoundOff( lonEdge(I,J) / PI_180, 4 )
       YEDGE(I,J,L)     = RoundOff( latEdge(I,J) / PI_180, 4 )

       ! Special treatment at uppermost edge
       IF ( I == NI ) THEN
          XEDGE(I+1,J,L) = RoundOff( lonEdge(I+1,J) / PI_180, 4 )
       ENDIF
       IF ( J == NJ ) THEN
          YEDGE(I,J+1,L)   = RoundOff( latEdge(I,J+1) / PI_180, 4 )
          YEDGE_R(I,J+1,L) = YEDGE(I,J+1,L) * PI_180
          YSIN(I,J+1,L)    = SIN( YEDGE_R(I,J+1,L) )
       ENDIF

       ! Special quantities directly derived from YEDGE
       YEDGE_R(I,J,L)   = YEDGE(I,J,L) * PI_180
       YEDGE_VAL        = YEDGE_R(I,J,L)           ! Lat edge [radians]
       YSIN_VAL         = SIN( YEDGE_VAL)          ! SIN( lat edge )
       YSIN(I,J,L)      = YSIN_VAL                 ! Store in YSIN array

    ENDDO
    ENDDO
    ENDDO

    ! Return w/ success
    RC = GC_SUCCESS

  END SUBROUTINE SetGridFromCtrEdges
!EOC
#endif
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Set_Xoffset
!
! !DESCRIPTION: Function SET\_XOFFSET initializes the nested-grid longitude 
!  offset variable I0.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Set_xOffSet( X_OFFSET )
!
! !INPUT PARAMETERS: 
!
    INTEGER, INTENT(IN) :: X_OFFSET  ! Value to assign to I0
!
! !REVISION HISTORY:
!  24 Feb 2012 - R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
    I0 = X_OFFSET

  END SUBROUTINE Set_xOffSet
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Set_Yoffset
!
! !DESCRIPTION: Function SET\_YOFFSET initializes the nested-grid latitude 
!  offset variable J0.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Set_yOffSet( Y_OFFSET )
!
! !INPUT PARAMETERS: 
!
    INTEGER, INTENT(IN) :: Y_OFFSET  ! Value to assign to J0
!
! !REVISION HISTORY:
!  24 Feb 2012 - R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC

    J0 = Y_OFFSET

  END SUBROUTINE Set_yOffSet
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get_Xoffset
!
! !DESCRIPTION: Function GET\_XOFFSET returns the nested-grid longitude 
!  offset to the calling program.
!\\
!\\
! !INTERFACE:
!
  FUNCTION Get_xOffSet( GLOBAL ) RESULT( X_OFFSET )
!
! !INPUT PARAMETERS: 
!
    ! If GLOBAL is passed, then return the actual window offset.
    ! This is necessary for certain instances (e.g. diagnostics)
    LOGICAL, INTENT(IN), OPTIONAL :: GLOBAL
!
! !RETURN VALUE:
!
    INTEGER                       :: X_OFFSET
!
! !REVISION HISTORY:
!  24 Feb 2012 - R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
    IF ( PRESENT( GLOBAL ) ) THEN

       ! If GLOBAL is passed, then return the actual window offset.
       ! This is necessary for certain instances (e.g. diagnostics)
       X_OFFSET = I0

    ELSE

       ! Otherwise, if we have a nested grid, then all of the met
       ! fields have been cut down to size already.  Return 0.
       X_OFFSET = 0
       
    ENDIF

  END FUNCTION Get_xOffSet
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get_Yoffset
!
! !DESCRIPTION: Function GET\_XOFFSET returns the nested-grid longitude 
!  offset to the calling program.
!\\
!\\
! !INTERFACE:
!
  FUNCTION Get_yOffSet( GLOBAL ) RESULT( Y_OFFSET )
!
! !INPUT PARAMETERS: 
!
    ! If GLOBAL is passed, then return the actual window offset.
    ! This is necessary for certain instances (e.g. diagnostics)
    LOGICAL, INTENT(IN), OPTIONAL :: GLOBAL
!
! !RETURN VALUE:
!
    INTEGER                       :: Y_OFFSET
!
! !REVISION HISTORY:
!  24 Feb 2012 - R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
    IF ( PRESENT( GLOBAL ) ) THEN 

       ! If GLOBAL is passed, then return the actual window offset.
       ! This is necessary for certain instances (e.g. diagnostics)
       Y_OFFSET = J0

    ELSE

       ! Otherwise, if we have a nested grid, then all of the met
       ! fields have been cut down to size already.  Return 0.
       Y_OFFSET = 0

    ENDIF

  END FUNCTION Get_yOffSet
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get_Xmid
!
! !DESCRIPTION: Function GET\_XMID returns the longitude in degrees at the 
!  center of a GEOS-Chem grid box.
!\\
!\\
! !INTERFACE:
!
  FUNCTION Get_xMid( I, J, L ) RESULT( X )
!
! !INPUT PARAMETERS: 
!
    INTEGER, INTENT(IN) :: I   ! Longitude index
    INTEGER, INTENT(IN) :: J   ! Latitude index
    INTEGER, INTENT(IN) :: L   ! Level index
!
! !RETURN VALUE:
!
    REAL(fp)              :: X   ! Corresponding lon value @ grid box ctr
!
! !REVISION HISTORY:
!  24 Feb 2012 - R. Yantosca - Initial version
!  09 Jun 2015 - C. Keller   - Now ensure that -180.0 <= x < +180.0.
!EOP
!------------------------------------------------------------------------------
!BOC

    X = XMID(I,J,1)

    ! Make sure value is between -180 and +180
    IF ( X >= 180.0_fp ) X = X - 360.0_fp

  END FUNCTION Get_xMid
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get_Xedge
!
! !DESCRIPTION: Function GET\_XEDGE returns the longitude in degrees at the 
!  western edge of a GEOS-Chem grid box.
!\\
!\\
! !INTERFACE:
!
  FUNCTION Get_xEdge( I, J, L ) RESULT( X )
!
! !INPUT PARAMETERS: 
!
    INTEGER, INTENT(IN) :: I   ! Longitude index
    INTEGER, INTENT(IN) :: J   ! Latitude index
    INTEGER, INTENT(IN) :: L   ! Level index
!
! !RETURN VALUE:
!
    REAL(fp)              :: X   ! Corresponding lon value @ W edge of grid box
!
! !REVISION HISTORY:
!  24 Feb 2012 - R. Yantosca - Initial version
!  09 Jun 2015 - C. Keller   - Now ensure that -180.0 <= x < +180.0.
!EOP
!------------------------------------------------------------------------------
!BOC

    X = XEDGE(I,J,1)

    ! Make sure value is between -180 and +180
    IF ( X >= 180.0_fp ) X = X - 360.0_fp

  END FUNCTION Get_xEdge
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get_Ymid
!
! !DESCRIPTION: Function GET\_YMID returns the latitude in degrees at the 
!  center of a GEOS-Chem grid box.
!\\
!\\
! !INTERFACE:
!
  FUNCTION Get_yMid( I, J, L ) RESULT( Y )
!
! !INPUT PARAMETERS: 
!
    INTEGER, INTENT(IN) :: I   ! Longitude index
    INTEGER, INTENT(IN) :: J   ! Latitude index
    INTEGER, INTENT(IN) :: L   ! Level index
!
! !RETURN VALUE:
!
    REAL(fp)              :: Y   ! Latitude value at @ grid box ctr [degrees]
!
! !REVISION HISTORY:
!  24 Feb 2012 - R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC

    Y = YMID(I,J,1)
      
  END FUNCTION Get_yMid
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get_Yedge
!
! !DESCRIPTION: Function GET\_YEDGE returns the latitude in degrees at the 
!  southern edge of a GEOS-Chem grid box.
!\\
!\\
! !INTERFACE:
!
  FUNCTION Get_yEdge( I, J, L ) RESULT( Y )
!
! !INPUT PARAMETERS: 
!
!
    INTEGER, INTENT(IN) :: I   ! Longitude index
    INTEGER, INTENT(IN) :: J   ! Latitude index
    INTEGER, INTENT(IN) :: L   ! Level index
!
! !RETURN VALUE:
!
    REAL(fp)              :: Y   ! Latitude value @ S edge of grid box [degrees]
!
! !REVISION HISTORY:
!  24 Feb 2012 - R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC

    Y = YEDGE(I,J,1)

  END FUNCTION Get_yEdge
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get_Ymid_R
!
! !DESCRIPTION: Function GET\_YMID\_R returns the latitude in radians at 
!  the center of a GEOS-Chem grid box.
!\\
!\\
! !INTERFACE:
!
  FUNCTION Get_yMid_R( I, J, L ) RESULT( Y )
!
! !INPUT PARAMETERS: 
!
    INTEGER, INTENT(IN) :: I   ! Longitude index
    INTEGER, INTENT(IN) :: J   ! Latitude index
    INTEGER, INTENT(IN) :: L   ! Level index
!
! !RETURN VALUE:
!
    REAL(fp)              :: Y   ! Latitude value at @ grid box ctr [radians]
!
! !REVISION HISTORY:
!  24 Feb 2012 - R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC

    Y = YMID_R(I,J,1)

  END FUNCTION Get_yMid_R
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get_Ymid_R_W
!
! !DESCRIPTION: Function GET\_YMID3\_R\_W returns the latitude in radians at 
!  the center of a GEOS-Chem grid box for the GEOS-5 nested grid.
!\\
!\\
! !INTERFACE:
!
  FUNCTION Get_yMid_R_W( I, J, L ) RESULT( Y )
!
! !INPUT PARAMETERS: 
!
    INTEGER, INTENT(IN) :: I   ! Longitude index
    INTEGER, INTENT(IN) :: J   ! Latitude index
    INTEGER, INTENT(IN) :: L   ! Level index
!
! !RETURN VALUE:
!
    REAL(fp)              :: Y   ! Latitude value at @ grid box ctr [radians]
!
! !REVISION HISTORY:
!  24 Feb 2012 - R. Yantosca - Initial version
!  01 Mar 2012 - R. Yantosca - Bracket with #ifdef for nested grids
!EOP
!------------------------------------------------------------------------------
!BOC

#if defined( NESTED_CH ) || defined( NESTED_EU ) || defined( NESTED_NA ) || defined( NESTED_AS )

    ! For nested grids, return the latitude center of the window
    ! region (in radians)
    Y = YMID_R_W(I,J,1)

#else

    ! Otherwise return a fake value
    Y = -1e+0_fp

#endif

  END FUNCTION Get_yMid_R_W
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get_Yedge_R
!
! !DESCRIPTION: Function GET\_YEDGE\_R returns the latitude in radians at 
!  the southern edge of a GEOS-Chem grid box.
!\\
!\\
! !INTERFACE:
!
  FUNCTION Get_yEdge_R( I, J, L ) RESULT( Y )
!
! !INPUT PARAMETERS: 
!
    INTEGER, INTENT(IN) :: I   ! Longitude index
    INTEGER, INTENT(IN) :: J   ! Latitude index
    INTEGER, INTENT(IN) :: L   ! Level index
!
! !RETURN VALUE:
!
    REAL(fp)              :: Y   ! Latitude value @ S edge of grid box [radians]
!
! !REVISION HISTORY:
!  24 Feb 2012 - R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC

    Y = YEDGE_R(I,J,1)

  END FUNCTION Get_yEdge_R
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get_Ysin
!
! !DESCRIPTION: Function GET\_YSIN returns the sine of the southern edge 
! of a GEOS-Chem grid box. 
!\\
!\\
! !INTERFACE:
!
      FUNCTION Get_ySin( I, J, L ) RESULT( Y ) 
!   
! !INPUT PARAMETERS:
!                  
    INTEGER, INTENT(IN) :: I   ! Longitude index
    INTEGER, INTENT(IN) :: J   ! Latitude index
    INTEGER, INTENT(IN) :: L   ! Level index
!
! !RETURN VALUE:
!
    REAL(fp)              :: Y   ! Sine of Latitude value @ S edge of grid box
!
! !REVISION HISTORY:
!  03 Apr 2012 - M. Payer    -  Initial version (M. Cooper)
!EOP
!------------------------------------------------------------------------------
!BOC

    Y = YSIN(I,J,1)

  END FUNCTION Get_ySin
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get_Area_m2
!
! !DESCRIPTION: Function GET\_AREA\_M2 returns the surface area [m2] of a 
!  GEOS-Chem grid box.  
!\\
!\\
! !INTERFACE:
!
  FUNCTION Get_Area_m2( I, J, L ) RESULT( A )
!
! !INPUT PARAMETERS: 
!
    INTEGER, INTENT(IN) :: I   ! Longitude index
    INTEGER, INTENT(IN) :: J   ! Latitude index
    INTEGER, INTENT(IN) :: L   ! Level index
!
! !RETURN VALUE:
!
    REAL(fp)              :: A   ! Grid box surface area [m2]
!
! !REVISION HISTORY:
!  24 Feb 2012 - R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC

    A = AREA_M2(I,J,1)

  END FUNCTION Get_Area_m2
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get_Area_cm2
!
! !DESCRIPTION: Function GET\_AREA\_CM2 returns the surface area [cm2] of a 
!  GEOS-Chem grid box.  Works for nested grids too.
!\\
!\\
! !INTERFACE:
!
  FUNCTION Get_area_cm2( I, J, L ) RESULT( A )
!
! !INPUT PARAMETERS: 
!
    INTEGER, INTENT(IN) :: I   ! Longitude index
    INTEGER, INTENT(IN) :: J   ! Latitude index
    INTEGER, INTENT(IN) :: L   ! Level index
!
! !RETURN VALUE:
!
    REAL(fp)              :: A  ! Grid box surface area [cm2]
!
! !REVISION HISTORY:
!  24 Feb 2012 - R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC

    A = AREA_M2(I,J,1) * 1e+4_fp
    
  END FUNCTION Get_Area_cm2
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get_Bounding_Box
!
! !DESCRIPTION: Subroutine GET\_BOUNDING\_BOX returns the indices which 
!  specify the lower left (LL) and upper right (UR) corners of a rectangular 
!  region, given the corresponding longitude and latitude values. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Get_Bounding_Box( I1, I2, J1, J2, L, COORDS, INDICES )
!
! !INPUT PARAMETERS: 
!
    INTEGER, INTENT(IN)  :: I1, I2      ! Lon indices
    INTEGER, INTENT(IN)  :: J1, J2      ! Lat indices
    INTEGER, INTENT(IN)  :: L
    REAL(fp),  INTENT(IN)  :: COORDS(4)   ! (/LON_LL, LAT_LL, LON_UR, LAT_UR/)
!
! !INPUT/OUTPUT PARAMETERS: 
!
    INTEGER, INTENT(OUT) :: INDICES(4)  ! (/I_LL, J_LL, I_UR, J_UR/)
!
! !REMARKS:
!  For now, this only works with the surface layer (which is OK since this
!  routine is mostly just called to find a window for surface emissions) 
! 
! !REVISION HISTORY:
!  24 Feb 2012 - R. Yantosca - Initial version
!  01 Mar 2012 - R. Yantosca - Modified for  grids, added input parameters
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
! 
    INTEGER            :: I, J
    CHARACTER(LEN=255) :: LOCATION

    !=================================================================
    ! Initialization
    !=================================================================
      
    ! Location
    LOCATION = 'GET_BOUNDING_BOX (grid_mod.f)'

    ! Initialize
    INDICES(:) = 0

    !=================================================================
    ! Longitude search
    !=================================================================
    DO J = 1,  1
    DO I = I1, I2

       ! Locate index corresponding to the lower-left longitude
       IF ( COORDS(1) >  XEDGE(I,  J,1)   .and.          &
            COORDS(1) <= XEDGE(I+1,J,1) ) INDICES(1) = I

         ! Locate index corresponding to upper-right longitude
       IF ( COORDS(3) >  XEDGE(I,  J,1)   .and.          &
            COORDS(3) <= XEDGE(I+1,J,1) ) INDICES(3) = I

    ENDDO
    ENDDO

    ! Error check lower-left longitude
    IF ( INDICES(1) == 0 ) THEN
       CALL ERROR_STOP( 'Invalid lower-left lon index!',  LOCATION )
    ENDIF
    
    ! Error check upper-right longitude
    IF ( INDICES(3) == 0 ) THEN
       CALL ERROR_STOP( 'Invalid upper-right lon index!', LOCATION )
    ENDIF
      
    !=================================================================
    ! Latitude search
    !=================================================================
    DO J = J1, J2
    DO I = 1,  1

       ! Locate index corresponding to the lower-left latitude
       IF ( COORDS(2) >  YEDGE(I,J,  1)   .and.          &
            COORDS(2) <= YEDGE(I,J+1,1) ) INDICES(2) = J

       ! Locate index corresponding to the upper-right latitude
       IF ( COORDS(4) >  YEDGE(I,J,  1)   .and.          &
            COORDS(4) <= YEDGE(I,J+1,1) ) INDICES(4) = J

    ENDDO
    ENDDO

    ! Error check lower-left longitude
    IF ( INDICES(2) == 0 ) THEN
       CALL ERROR_STOP( 'Invalid lower-left lat index!',  LOCATION )
    ENDIF

    ! Error check upper-right longitude
    IF ( INDICES(4) == 0 ) THEN
       CALL ERROR_STOP( 'Invalid upper-right lat index!', LOCATION )
    ENDIF

  END SUBROUTINE Get_Bounding_Box
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Its_A_Nested_Grid
!
! !DESCRIPTION: Function GET\_AREA\_CM2 returns the surface area [cm2] of a 
!  GEOS-Chem grid box.  Works for nested grids too.
!\\
!\\
! !INTERFACE:
!
  FUNCTION ITS_A_NESTED_GRID() RESULT( IT_IS_NESTED )
!
! !RETURN VALUE:
!
    LOGICAL :: IT_IS_NESTED   ! =T if it's a nested grid; =F otherwise
!
! !REVISION HISTORY:
!  24 Feb 2012 - R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
    IT_IS_NESTED = IS_NESTED
    
  END FUNCTION ITS_A_NESTED_GRID
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get_IJ
!
! !DESCRIPTION: Function GET\_IJ returns I and J index for a LON, LAT
!  coordinate (dkh, 11/16/06). Updated to support nested domains and made much
!  simpler (zhe, 1/19/11).
!\\
!\\
! !INTERFACE:
!
  FUNCTION GET_IJ( LON, LAT ) RESULT ( IIJJ )
!
! !USES:
!
    USE CMN_SIZE_MOD
    USE ERROR_MOD,    ONLY : ERROR_STOP
!
! !INPUT PARAMETERS:
!
    REAL*4, INTENT(IN)  :: LAT, LON
!
! !RETURN VALUE:
!
    INTEGER               :: IIJJ(2)
!
! !REVISION HISTORY:
!  16 Jun 2017 - M. Sulprizio- Initial version based on routine from adjoint
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL(fp) :: TLON, TLAT
      REAL(fp) :: I0,   J0

      !=================================================================
      ! GET_IJ begins here!
      !=================================================================
      I0 = GET_XOFFSET( GLOBAL=.TRUE. )
      J0 = GET_YOFFSET( GLOBAL=.TRUE. )

      TLON = INT( ( LON + 180e+0_fp ) / DISIZE + 1.5e+0_fp )
      TLAT = INT( ( LAT +  90e+0_fp ) / DJSIZE + 1.5e+0_fp )

#if defined( NESTED_CH ) || defined( NESTED_NA ) || defined( NESTED_CA )
      TLON = TLON - I0
      TLAT = TLAT - J0
      IF ( TLAT < 1 .or. TLAT > JJPAR ) THEN
         CALL ERROR_STOP('Beyond the nested window', 'GET_IJ')
      ENDIF
#else
      IF ( TLON > IIPAR ) TLON = TLON - IIPAR

      ! Check for impossible values
      IF ( TLON > IIPAR .or. TLAT > JJPAR .or. &
           TLON < 1     .or. TLAT < 1          ) THEN
         CALL ERROR_STOP('Error finding grid box', 'GET_IJ')
      ENDIF

#endif

      IIJJ(1) = TLON
      IIJJ(2) = TLAT

    END FUNCTION GET_IJ
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Init_Grid
!
! !DESCRIPTION: Subroutine INIT\_GRID initializes variables and allocates
!  module arrays.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Init_Grid( am_I_Root, Input_Opt, IM, JM, LM, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod, ONLY : OptInput
    USE Registry_Mod,  ONLY : Registry_AddField   
!
! !INPUT PARAMETERS:
!
    LOGICAL,        INTENT(IN)  :: am_I_Root   ! Are we on the root CPU 
    TYPE(OptInput), INTENT(IN)  :: Input_Opt   ! Input Options object
    INTEGER,        INTENT(IN)  :: IM          ! # of longitudes on this CPU
    INTEGER,        INTENT(IN)  :: JM          ! # of latitudes  on this CPU
    INTEGER,        INTENT(IN)  :: LM          ! # of levels     on this CPU
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT) :: RC          ! Success or failure?
!
! !REVISION HISTORY:
!  24 Feb 2012 - R. Yantosca - Initial version, based on grid_mod.F
!  01 Mar 2012 - R. Yantosca - Now define IS_NESTED based on Cpp flags
!  03 Dec 2012 - R. Yantosca - Add am_I_Root, RC to argument list
!  04 Dec 2012 - R. Yantosca - Now dimension arrays with IM, JM, LM instead 
!                              of I1, J1, L1, I2, J2, L2. 
!  01 Apr 2015 - R. Yantosca - Now accept Input_Opt as an argument 
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
! 
    ! Scalars
    INTEGER            :: L,        AS

    REAL(fp), SAVE     :: TIME_IND(1)

    ! Strings
    CHARACTER(LEN=255) :: Variable, Desc, Units

    !======================================================================
    ! Initialize module variables
    !======================================================================

    ! Assume success
    RC        = GC_SUCCESS

    ! IS_NESTED is now a local shadow variable (bmy, 4/1/15)
    IS_NESTED = Input_Opt%ITS_A_NESTED_GRID

    ! We only need one level (ckeller, 01/02/15).
    L         = 1

    ALLOCATE( XMID   ( IM,   JM,   L ), STAT=RC )
    IF ( RC /= 0 ) CALL ALLOC_ERR( 'XMID' )
    XMID    = 0e+0_fp
    
    ALLOCATE( XEDGE  ( IM+1, JM,   L ), STAT=RC )
    IF ( RC /= 0 ) CALL ALLOC_ERR( 'XEDGE' )
    XEDGE   = 0e+0_fp
    
    ALLOCATE( YMID   ( IM,   JM,   L ), STAT=RC )
    IF ( RC /= 0 ) CALL ALLOC_ERR( 'YMID' )
    YMID    = 0e+0_fp
    
    ALLOCATE( YEDGE  ( IM,   JM+1, L ), STAT=RC )
    IF ( RC /= 0 ) CALL ALLOC_ERR( 'YEDGE' )
    YEDGE   = 0e+0_fp

    ALLOCATE( YSIN   ( IM,   JM+1, L ), STAT=RC )
    IF ( RC /= 0 ) CALL ALLOC_ERR( 'YSIN' )
    YSIN    = 0e+0_fp

    ALLOCATE( YMID_R ( IM,   JM,   L ), STAT=RC )               
    IF ( RC /= 0 ) CALL ALLOC_ERR( 'YMID_R' )
    YMID_R  = 0e+0_fp
   
    ALLOCATE( YEDGE_R( IM,   JM+1, L ), STAT=RC )
    IF ( RC /= 0 ) CALL ALLOC_ERR( 'YEDGE_R' )
    YEDGE_R = 0e+0_fp

    ALLOCATE( AREA_M2( IM,   JM,   L ), STAT=RC )
    IF ( RC /= 0 ) CALL ALLOC_ERR( 'AREA_M2' )
    AREA_M2 = 0e+0_fp

    !======================================================================
    ! Special settings for nested-grid simulations only
    !======================================================================

    ! Allocate nested-grid window array of lat centers (radians)
    IF ( Input_Opt%ITS_A_NESTED_GRID ) THEN
       ALLOCATE( YMID_R_W( IM, 0:JM+1, L ), STAT=AS ) 
       IF ( RC /= 0 ) CALL ALLOC_ERR( 'YMID_R_W' )
       YMID_R_W = 0e+0_fp
    ENDIF

  END SUBROUTINE Init_Grid
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Cleanup_Grid
!
! !DESCRIPTION: Subroutine CLEANUP\_GRID deallocates all module arrays.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Cleanup_Grid( am_I_Root, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE Registry_Mod, ONLY : Registry_Destroy
!
! !INPUT PARAMETERS:
!
    LOGICAL, INTENT(IN)  :: am_I_Root
!
! !OUTPUT PARAMETERS:
!
    INTEGER, INTENT(OUT) :: RC
!
! !REVISION HISTORY:
!  24 Feb 2012 - R. Yantosca - Initial version, based on grid_mod.F
!EOP
!------------------------------------------------------------------------------
!BOC
!
! ! LOCAL VARIABLES
!
    ! Strings
    CHARACTER(LEN=255) :: ErrMsg, ThisLoc

    !=======================================================================
    ! Initialize
    !=======================================================================
    RC      = GC_SUCCESS
    ErrMsg  = ''
    ThisLoc = ' -> at Cleanup_Grid (in GeosUtil/gc_grid_mod.F90)'

    !=======================================================================
    ! Deallocate module variables
    !=======================================================================
    IF ( ALLOCATED( XMID      ) ) DEALLOCATE( XMID      )
    IF ( ALLOCATED( XEDGE     ) ) DEALLOCATE( XEDGE     )
    IF ( ALLOCATED( YMID      ) ) DEALLOCATE( YMID      )
    IF ( ALLOCATED( YEDGE     ) ) DEALLOCATE( YEDGE     )
    IF ( ALLOCATED( YSIN      ) ) DEALLOCATE( YSIN      ) 
    IF ( ALLOCATED( YMID_R    ) ) DEALLOCATE( YMID_R    )
    IF ( ALLOCATED( YMID_R_W  ) ) DEALLOCATE( YMID_R_W  )  
    IF ( ALLOCATED( YEDGE_R   ) ) DEALLOCATE( YEDGE_R   )
    IF ( ALLOCATED( AREA_M2   ) ) DEALLOCATE( AREA_M2   )

  END SUBROUTINE Cleanup_Grid
!EOC
END MODULE GC_Grid_Mod
