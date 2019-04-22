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
  USE ErrCode_Mod
  USE Error_Mod
  USE Precision_Mod
  USE PhysConstants
  USE Registry_Mod,   ONLY : MetaRegItem
  USE State_Grid_Mod, ONLY : GrdState

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC  :: Compute_Grid
  PUBLIC  :: DoGridComputation
  PUBLIC  :: Get_Bounding_Box
  PUBLIC  :: SetGridFromCtr
#if defined ( MODEL_WRF )
  PUBLIC  :: SetGridFromCtrEdges
#endif
  PUBLIC  :: GET_IJ
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
!  10 Mar 2019 - M. Sulprizio- Move module arrays to fields in State_Grid
!EOP
!------------------------------------------------------------------------------
!BOC
!
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
  SUBROUTINE Compute_Grid( am_I_Root, Input_Opt, State_Grid, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod, ONLY : OptInput
!
! !INPUT PARAMETERS:
!
    LOGICAL,        INTENT(IN)    :: am_I_Root         ! Root CPU?
    TYPE(OptInput), INTENT(IN)    :: Input_Opt         ! Input Options
!
! !INPUT/OUTPUT PARAMETERS:
!  
    TYPE(GrdState), INTENT(INOUT) :: State_Grid        ! Grid State object
!
! !OUTPUT PARAMETERS:
!  
    INTEGER,        INTENT(OUT)   :: RC                ! Success/failure?
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
  Call DoGridComputation( am_I_Root, Input_Opt, State_Grid, RC )

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
  SUBROUTINE DoGridComputation( am_I_Root, Input_Opt, State_Grid, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod,  ONLY : OptInput
    USE State_Grid_Mod, ONLY : GrdState
!
! !INPUT PARAMETERS:
!
    LOGICAL,        INTENT(IN)    :: am_I_Root         ! Root CPU?
    TYPE(OptInput), INTENT(IN)    :: Input_Opt         ! Input Options
!
! !INPUT/OUTPUT PARAMETERS:
!  
    TYPE(GrdState), INTENT(INOUT) :: State_Grid        ! Grid State object
!
! !OUTPUT PARAMETERS:
!  
    INTEGER,        INTENT(OUT)   :: RC                ! Success or failure?
!
! !REMARKS:
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
    INTEGER  :: I, J, L, IG, JG
    REAL(fp) :: YEDGE_VAL, YSIN_VAL
    REAL(fp) :: SIN_N, SIN_S, SIN_DIFF

    ! Arrays
    REAL(fp) :: Ind_X(State_Grid%NX+1)
    REAL(fp) :: Ind_Y(State_Grid%NY+1)

    ! Strings
    CHARACTER(LEN=255) :: ErrMsg
    CHARACTER(LEN=255) :: ThisLoc

    !======================================================================
    ! Initialization
    !======================================================================

    ! Assume success
    RC = GC_SUCCESS
    ThisLoc = 'DoGridComputation (gc_grid_mod.F90)'

    !======================================================================
    ! Vertical Grid
    !======================================================================

    ! Both GEOS-FP and MERRA-2 have native vertical resolution of 72 levels
    State_Grid%NativeNZ = 72

    ! Hardcode maximum number of levels below tropopause and stratopause
    ! (formerly set in CMN_SIZE_mod.F)
    IF ( State_Grid%NZ == 47 ) THEN
       State_Grid%MaxTropLev  = 38
       State_Grid%MaxStratLev = 44
    ELSE IF ( State_Grid%NZ == 72 ) THEN
       State_Grid%MaxTropLev  = 40
       State_Grid%MaxStratLev = 59
    ENDIF

    ! Set maximum number of levels in the chemistry grid
    IF ( Input_Opt%LUCX ) THEN
       State_Grid%MaxChemLev  = State_Grid%MaxStratLev
    ELSE
       State_Grid%MaxChemLev  = State_Grid%MaxTropLev
    ENDIF

    !======================================================================
    ! Global Horizontal Grid
    !
    ! First, we need to compute the XMid and YMid values on the global
    ! grid at the specified resolution to that we can compute X and Y
    ! offsets.
    !======================================================================

    ! Compute number of grid boxes on global grid
    State_Grid%GlobalNX =   360.0_fp / State_Grid%DX
    State_Grid%GlobalNY = ( 180.0_fp / State_Grid%DY ) + 1

    !----------------------------------------
    ! Arrays of delta-X and delta-Y values
    !----------------------------------------

    ! Define global array of delta-X and delta-Y values
    ! Loop over horizontal grid
    DO J = 1, State_Grid%GlobalNY
    DO I = 1, State_Grid%GlobalNX
       State_Grid%DeltaX(I,J) = State_Grid%DX
       State_Grid%DeltaY(I,J) = State_Grid%DY

       ! If using half-sized polar boxes on a global grid,
       ! multiply delta-Y by 1/2 at the South and North Poles
       IF ( (J == 1 .or. J == State_Grid%NY) .and. State_Grid%HalfPolar ) THEN
          State_Grid%DeltaY(I,J) = State_Grid%DY * 0.5_fp
       ENDIF

    ENDDO
    ENDDO

    !----------------------------------------------------------------------
    ! Calculate grid box centers on global grid
    !----------------------------------------------------------------------

    ! Allocate arrays
    ALLOCATE( State_Grid%GlobalXMid(State_Grid%GlobalNX,State_Grid%GlobalNY),&
              STAT=RC )
    CALL GC_CheckVar( 'State_Grid%GlobalXMid', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Grid%GlobalXMid = 0e+0_fp
    
    ALLOCATE( State_Grid%GlobalYMid(State_Grid%GlobalNX,State_Grid%GlobalNY),&
              STAT=RC )
    CALL GC_CheckVar( 'State_Grid%GlobalYMid', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Grid%GlobalYMid = 0e+0_fp

    ! Loop over horizontal grid
    DO J = 1, State_Grid%GlobalNY
    DO I = 1, State_Grid%GlobalNX

       !--------------------------------
       ! Longitude centers [degrees]
       !--------------------------------
       State_Grid%GlobalXMid(I,J) = ( State_Grid%DeltaX(I,J) * (I-1) ) - &
                                      180e+0_fp

       !--------------------------------
       ! Latitude centers [degrees]
       !--------------------------------
       State_Grid%GlobalYMid(I,J) = ( State_Grid%DeltaY(I,J) * (J-1) ) - &
                                      90e+0_fp

       ! If using half-sized polar boxes, multiply DY by 1/2 at poles
       IF ( State_Grid%HalfPolar .and. J == 1)  THEN
          ! South Pole
          State_Grid%GlobalYMid(I,J) =  &
             -90e+0_fp + ( 0.5e+0_fp * State_Grid%DeltaY(I,J) )
       ENDIF
       IF ( State_Grid%HalfPolar .and. J == State_Grid%GlobalNY ) THEN
          ! North Pole
          State_Grid%GlobalYMid(I,J) = &
             +90e+0_fp - ( 0.5e+0_fp * State_Grid%DeltaY(I,J) )
       ENDIF

    ENDDO
    ENDDO

    !======================================================================
    ! User-defined Horizontal Grid
    !======================================================================

    ! Determine X offsets based on global grid
    DO I = 1, State_Grid%GlobalNX
       IF ( State_Grid%GlobalXMid(I,1) >= State_Grid%XMin ) THEN
          State_Grid%XMinOffset = I-1
          EXIT
       ENDIF
    ENDDO
    DO I = 1, State_Grid%GlobalNX
       IF ( State_Grid%GlobalXMid(I,1)+State_Grid%DX >= State_Grid%XMax ) THEN
          State_Grid%XMaxOffset = I
          EXIT
       ENDIF
    ENDDO

    ! Determine Y offsets based on global grid
    DO J = 1, State_Grid%GlobalNY
       IF ( State_Grid%GlobalYMid(1,J) >= State_Grid%YMin ) THEN
          State_Grid%YMinOffset = J-1
          EXIT
       ENDIF
    ENDDO
    DO J = 1, State_Grid%GlobalNY
       IF ( State_Grid%GlobalYMid(1,J)+State_Grid%DY >= State_Grid%YMax ) THEN
          State_Grid%YMaxOffset = J
          EXIT
       ENDIF
    ENDDO

    !----------------------------------------------------------------------
    ! Calculate grid box centers and edges on global grid
    !----------------------------------------------------------------------   

    ! Loop over horizontal grid
    DO J = 1, State_Grid%NY
    DO I = 1, State_Grid%NX

       ! Index value for user-defined grid on the global grid
       IG = I + ( State_Grid%XMinOffset - 1 )
       JG = J + ( State_Grid%YMinOffset - 1 )
       
       !--------------------------------
       ! Longitude centers [degrees]
       !--------------------------------
       State_Grid%XMid(I,J) = ( State_Grid%DeltaX(I,J) * IG ) - 180e+0_fp
          
       !--------------------------------
       ! Longitude edges [degrees]
       !--------------------------------
       State_Grid%XEdge(I,J) = State_Grid%XMid(I,J) - &
                             ( State_Grid%DeltaX(I,J) * 0.5e+0_fp )

       ! Compute the last longitude edge
       IF ( I == State_Grid%NX ) THEN
          State_Grid%XEdge(I+1,J) = State_Grid%XEdge(I,J) + &
                                    State_Grid%DeltaX(I,J)
       ENDIF

       !--------------------------------
       ! Latitude centers [degrees]
       !--------------------------------
       State_Grid%YMid(I,J) = ( State_Grid%DeltaY(I,J) * JG ) - 90e+0_fp

       ! If using half-sized polar boxes on a global grid,
       ! multiply DY by 1/2 at poles
       IF ( State_Grid%HalfPolar .and. .not. State_Grid%NestedGrid .and. &
            J == 1 ) THEN
          ! South Pole
          State_Grid%YMid(I,J) = -90e+0_fp + &
                                 ( 0.5e+0_fp * State_Grid%DeltaY(I,J) )
       ENDIF
       IF ( State_Grid%HalfPolar .and. .not. State_Grid%NestedGrid .and. &
            J == State_Grid%NY ) THEN
          ! North Pole
          State_Grid%YMid(I,J) = +90e+0_fp - &
                                 ( 0.5e+0_fp * State_Grid%DeltaY(I,J) )
       ENDIF

       !--------------------------------
       ! Latitude centers [radians]
       !--------------------------------
       State_Grid%YMid_R(I,J) = ( PI_180 * State_Grid%YMid(I,J)  )

       !--------------------------------
       ! Latitude edges [degrees]
       !--------------------------------
       State_Grid%YEdge(I,J) = State_Grid%YMid(I,J) - &
                               ( State_Grid%DeltaY(I,J) * 0.5e+0_fp )

       ! If using half-sized polar boxeson a global grid,
       ! force the northern edge of grid boxes along the SOUTH POLE
       ! to be -90 degrees latitude
       IF ( State_Grid%HalfPolar .and. .not. State_Grid%NestedGrid .and. &
            J == 1 ) THEN
          State_Grid%YEdge(I,J) = -90e+0_fp
       ENDIF

       !--------------------------------
       ! Lat edges [radians]
       !--------------------------------
       State_Grid%YEdge_R(I,J) = ( PI_180  * State_Grid%YEdge(I,J) )
          
       ! mjc - Compute sine of latitude edges (needed for map_a2a regrid)
       YEDGE_VAL = State_Grid%YEdge_R(I,J) ! Lat edge in radians
       YSIN_VAL  = SIN( YEDGE_VAL )          ! SIN( lat edge )
       State_Grid%YSIN(I,J) = YSIN_VAL     ! Store in YSIN array

       ! Test for North Pole if using global grid
       IF ( .not. State_Grid%NestedGrid .and. J == State_Grid%NY ) THEN

#if defined( ESMF_ ) || defined( EXTERNAL_GRID ) || defined( EXTERNAL_FORCING )
          !----------------------------------------------------------------
          !         %%%%%%% GEOS-Chem HP (with ESMF & MPI) %%%%%%%
          !
          ! Do not define half-sized polar boxes (bmy, 3/21/13)
          !----------------------------------------------------------------
          State_Grid%YEdge(I,State_Grid%NY+1) = &
          State_Grid%YEdge(I,State_Grid%NY  ) + &
          State_Grid%DeltaY(I,State_Grid%NY)
#else
          !----------------------------------------------------------------
          !         %%%%%%% GEOS-Chem CLASSIC (with OpenMP) %%%%%%%
          !
          ! Current practice in the standard GEOS-Chem is to force
          ! the northern edge of grid boxes along the NORTH POLE to
          ! be +90 degrees latitude. (bmy, 3/21/13)
          !----------------------------------------------------------------
          State_Grid%YEdge(I,State_Grid%NY+1) = +90e+0_fp
#endif
          State_Grid%YEdge_R(I,State_Grid%NY+1) = &
          State_Grid%YEdge  (I,State_Grid%NY+1) * PI_180

          ! Also compute sine of last latitude edge! (ckeller, 02/13/12)
          YEDGE_VAL = State_Grid%YEdge_R(I,State_Grid%NY+1)
          YSIN_VAL  = SIN( YEDGE_VAL )
          State_Grid%YSIN(I,State_Grid%NY+1) = YSIN_VAL

       ENDIF

    ENDDO
    ENDDO

    !======================================================================
    ! Compute grid box surface areas
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

    ! Loop over horizontal grid
    DO J = 1, State_Grid%NY
    DO I = 1, State_Grid%NX

       ! Sine of latitudes at N and S edges of grid box (I,J)
       SIN_N = SIN( State_Grid%YEdge_R(I,J+1) )
       SIN_S = SIN( State_Grid%YEdge_R(I,J  ) )

       ! Difference of sin(latitude) at N and S edges of grid box
       SIN_DIFF = SIN_N - SIN_S

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
       IF ( J == 1 ) THEN
          SIN_DIFF = 2e+0_fp * ( SIN_N - ( -1e+0_fp ) )
       ENDIF

       ! North pole kludge
       IF ( .not. State_Grid%NestedGrid .and. J == State_Grid%NY ) THEN
          SIN_DIFF = 2e+0_fp * ( 1e+0_fp - SIN_S )
       ENDIF
#endif

       ! Grid box surface areas [m2]
       State_Grid%Area_M2(I,J) = ( State_Grid%DeltaX(I,J) * PI_180 ) * &
                                   ( Re**2 ) * SIN_DIFF

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
       WRITE( 6, '(a)' )
       WRITE( 6, '(''%%%%%%%%%%%%%%% GLOBAL GRID %%%%%%%%%%%%%%%'')' )
       WRITE( 6, '(a)' )
       WRITE( 6, '(''Grid box longitude centers [degrees]: '')' )
       WRITE( 6, '(8(f8.3,1x))' ) ( State_Grid%GlobalXMid(I,1), &
                                    I=1,State_Grid%GlobalNX )
       WRITE( 6, '(a)' )
       WRITE( 6, '(''Grid box latitude centers [degrees]: '')' )
       WRITE( 6, '(8(f8.3,1x))' ) ( State_Grid%GlobalYMid(1,J), &
                                    J=1,State_Grid%GlobalNY )
       WRITE( 6, '(a)' )
       WRITE( 6, '(''%%%%%%%%%%%% USER-DEFINED GRID %%%%%%%%%%%%'')' )
       WRITE( 6, '(a)' )
       WRITE( 6, * ) ' XMinOffset : ', State_Grid%XMinOffset
       WRITE( 6, * ) ' XMaxOffset : ', State_Grid%XMaxOffset
       WRITE( 6, * ) ' YMinOffset : ', State_Grid%YMinOffset
       WRITE( 6, * ) ' YMaxOffset : ', State_Grid%YMaxOffset
       WRITE( 6, '(a)' )
       WRITE( 6, '(a)' )
       WRITE( 6, '(''Grid box longitude centers [degrees]: '')' )
       WRITE( 6, '(8(f8.3,1x))' ) ( State_Grid%XMid(I,1), I=1,State_Grid%NX )
       WRITE( 6, '(a)' )
       WRITE( 6, '(''Grid box longitude edges [degrees]: '')' )
       WRITE( 6, '(8(f8.3,1x))' ) ( State_Grid%XEdge(I,1), I=1,State_Grid%NX+1 )
       WRITE( 6, '(a)' )
       WRITE( 6, '(''Grid box latitude centers [degrees]: '')' )
       WRITE( 6, '(8(f8.3,1x))' ) ( State_Grid%YMid(1,J), J=1,State_Grid%NY )
       WRITE( 6, '(a)' )
       WRITE( 6, '(''Grid box latitude edges [degrees]: '')' )
       WRITE( 6, '(8(f8.3,1x))' ) ( State_Grid%YEdge(1,J), J=1,State_Grid%NY+1 )
       WRITE( 6, '(a)' )
       WRITE( 6, '(''SIN( grid box latitude edges )'')' )
       WRITE( 6, '(8(f8.3,1x))' ) ( State_Grid%YSIN(1,J), J=1,State_Grid%NY+1 )
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
  SUBROUTINE SetGridFromCtr( am_I_Root,  NX, NY, lonCtr, latCtr, &
                             State_Grid, RC ) 
!
! USES
!
    USE ErrCode_Mod
    USE Roundoff_Mod
!
! !INPUT PARAMETERS: 
!
    LOGICAL,        INTENT(IN)    :: am_I_Root      ! Root CPU?
    INTEGER,        INTENT(IN)    :: NX             ! # of lons
    INTEGER,        INTENT(IN)    :: NY             ! # of lats
    REAL(f4),       INTENT(IN)    :: lonCtr(NX,NY)  ! Lon ctrs [rad]
    REAL(f4),       INTENT(IN)    :: latCtr(NX,NY)  ! Lat ctrs [rad]
    TYPE(GrdState), INTENT(IN)    :: State_Grid     ! Grid State object
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(INOUT) :: RC 
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
    NI = SIZE(State_Grid%XMid,1)
    NJ = SIZE(State_Grid%XMid,2)

    ! Horizontal dimensions must agree
    IF ( State_Grid%NX /= NI .OR. State_Grid%NY /= NJ ) THEN
       WRITE(ErrMsg,*) 'Grid dimension mismatch: ',State_Grid%NX, '/=',NI, &
                       ' and/or ',State_Grid%NY,'/=',NJ
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    ! Loop over all grid boxes
    DO J = 1, NJ
    DO I = 1, NI

       ! Mid points: get directly from passed value
       State_Grid%XMid(I,J)   = RoundOff( lonCtr(I,J) / PI_180, 4 )
       State_Grid%YMid(I,J)   = RoundOff( latCtr(I,J) / PI_180, 4 )
       State_Grid%YMid_R(I,J) = State_Grid%YMid(I,J) * PI_180 

       ! Edges: approximate from neighboring mid points.
       IF ( I == 1 ) THEN
          TMP = RoundOff( lonCtr(I+1,J) / PI_180, 4 )
          State_Grid%XEdge(I,J) = State_Grid%XMid(I,J) - &
                                  ( ( TMP - State_Grid%XMid(I,J) ) / 2.0_f4 )
       ELSE
          State_Grid%XEdge(I,J) = ( State_Grid%XMid(I,J) + &
                                      State_Grid%XMid(I-1,J) ) / 2.0_f4
       ENDIF

       IF ( J == 1 ) THEN
          TMP = RoundOff( latCtr(I,J+1) / PI_180, 4 )
          State_Grid%YEdge(I,J) = State_Grid%YMid(I,J) - &
                                  ( ( TMP - State_Grid%YMid(I,J) ) / 2.0_f4 )
       ELSE
          State_Grid%YEdge(I,J) = ( State_Grid%YMid(I,J) + &
                                      State_Grid%YMid(I,J-1) ) / 2.0_f4
       ENDIF

       ! Special treatment at uppermost edge
       IF ( I == NI ) THEN
          State_Grid%XEdge(I+1,J) = State_Grid%XMid(I,J) + &
             ( ( State_Grid%XMid(I,J) - State_Grid%XMid(I-1,J) ) / 2.0_f4 )
       ENDIF
       IF ( J == NJ ) THEN
          State_Grid%YEdge(I,J+1) = State_Grid%YMid(I,J) + &
             ( ( State_Grid%YMid(I,J) - State_Grid%YMid(I,J-1) ) / 2.0_f4 )
       ENDIF

       ! Special quantities directly derived from State_Grid%YEdge
       State_Grid%YEdge_R(I,J) = State_Grid%YEdge(I,J) * PI_180
       YEDGE_VAL                 = State_Grid%YEdge_R(I,J) ! Lat edge[radians]
       YSIN_VAL                  = SIN( YEDGE_VAL)           ! SIN( lat edge )
       State_Grid%YSIN(I,J)    = YSIN_VAL                  ! Store in YSIN

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
  SUBROUTINE SetGridFromCtrEdges( am_I_Root, NX, NY, lonCtr, latCtr, &
                                  lonEdge, latEdge, State_Grid, RC )
!
! USES
!
    USE ErrCode_Mod
    USE Roundoff_Mod
!
! !INPUT PARAMETERS:
!
    LOGICAL,        INTENT(IN)    :: am_I_Root         ! Root CPU?
    INTEGER,        INTENT(IN)    :: NX                ! # of lons
    INTEGER,        INTENT(IN)    :: NY                ! # of lats
    REAL(f4),       INTENT(IN)    :: lonCtr (NX,NY)    ! Lon ctrs [rad]
    REAL(f4),       INTENT(IN)    :: latCtr (NX,NY)    ! Lat ctrs [rad]
    REAL(f4),       INTENT(IN)    :: lonEdge(NX+1,NY)  ! Lon edges [rad]
    REAL(f4),       INTENT(IN)    :: latEdge(NX,NY+1)  ! Lat edges [rad]
    TYPE(GrdState), INTENT(IN)    :: State_Grid        ! Grid State object
!!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(INOUT) :: RC
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
    NI = SIZE(State_Grid%XMid,1)
    NJ = SIZE(State_Grid%XMid,2)

    ! Horizontal dimensions must agree
    IF ( State_Grid%NX /= NI .OR. State_Grid%NY /= NJ ) THEN
       WRITE(ErrMsg,*) 'Grid dimension mismatch: ',State_Grid%NX, '/=',NI, &
                       ' and/or ',State_Grid%NY,'/=',NJ
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    ! Loop over all grid boxes
    DO J = 1, NJ
    DO I = 1, NI

       ! Mid points: get directly from passed value
       State_Grid%XMid(I,J)      = RoundOff( lonCtr(I,J) / PI_180, 4 )
       State_Grid%YMid(I,J)      = RoundOff( latCtr(I,J) / PI_180, 4 )
       State_Grid%YMid_R(I,J)    = State_Grid%YMid(I,J) * PI_180

       ! Edges: get directly from passed value
       State_Grid%XEdge(I,J)     = RoundOff( lonEdge(I,J) / PI_180, 4 )
       State_Grid%YEdge(I,J)     = RoundOff( latEdge(I,J) / PI_180, 4 )

       ! Special treatment at uppermost edge
       IF ( I == NI ) THEN
          XEDGE(I+1,J) = RoundOff( lonEdge(I+1,J) / PI_180, 4 )
       ENDIF
       IF ( J == NJ ) THEN
          State_Grid%YEDdge(I,J+1)   = RoundOff( latEdge(I,J+1) / PI_180, 4 )
          State_Grid%YEdge_R(I,J+1) = State_Grid%YEdge(I,J+1) * PI_180
          State_Grid%YSIN(I,J+1)    = SIN( State_Grid%YEdge_R(I,J+1) )
       ENDIF

       ! Special quantities directly derived from YEDGE
       State_Grid%YEdge_R(I,J)   = State_Grid%YEdge(I,J) * PI_180
       State_Grid%YEdge_VAL      = State_Grid%YEdge_R(I,J) ! Lat edge [radians]
       State_Grid%YSIN_VAL       = SIN( YEDGE_VAL)         ! SIN( lat edge )
       State_Grid%YSIN(I,J)      = YSIN_VAL                ! Store in YSIN array

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
! !IROUTINE: Get_Bounding_Box
!
! !DESCRIPTION: Subroutine GET\_BOUNDING\_BOX returns the indices which 
!  specify the lower left (LL) and upper right (UR) corners of a rectangular 
!  region, given the corresponding longitude and latitude values. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Get_Bounding_Box( I1, I2, J1, J2, L, COORDS, State_Grid, &
                               INDICES )
!
! !INPUT PARAMETERS: 
!
    INTEGER,        INTENT(IN)  :: I1, I2     ! Lon indices
    INTEGER,        INTENT(IN)  :: J1, J2     ! Lat indices
    INTEGER,        INTENT(IN)  :: L
    REAL(fp),       INTENT(IN)  :: COORDS(4)  ! (LON_LL,LAT_LL,LON_UR,LAT_UR)
    TYPE(GrdState), INTENT(IN)  :: State_Grid ! Grid State object
!
! !INPUT/OUTPUT PARAMETERS: 
!
    INTEGER,        INTENT(OUT) :: INDICES(4) ! (I_LL,J_LL,I_UR,J_UR)
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
       IF ( COORDS(1) >  State_Grid%XEdge(I,  J)   .and.          &
            COORDS(1) <= State_Grid%XEdge(I+1,J) ) INDICES(1) = I

         ! Locate index corresponding to upper-right longitude
       IF ( COORDS(3) >  State_Grid%XEdge(I,  J)   .and.          &
            COORDS(3) <= State_Grid%XEdge(I+1,J) ) INDICES(3) = I

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
       IF ( COORDS(2) >  State_Grid%YEdge(I,J  )   .and.          &
            COORDS(2) <= State_Grid%YEdge(I,J+1) ) INDICES(2) = J

       ! Locate index corresponding to the upper-right latitude
       IF ( COORDS(4) >  State_Grid%YEdge(I,J  )   .and.          &
            COORDS(4) <= State_Grid%YEdge(I,J+1) ) INDICES(4) = J

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
! !IROUTINE: Get_IJ
!
! !DESCRIPTION: Function GET\_IJ returns I and J index for a LON, LAT
!  coordinate (dkh, 11/16/06). Updated to support nested domains and made much
!  simpler (zhe, 1/19/11).
!\\
!\\
! !INTERFACE:
!
  FUNCTION GET_IJ( LON, LAT, State_Grid ) RESULT ( IIJJ )
!
! !USES:
!
    USE CMN_SIZE_MOD
    USE ERROR_MOD,      ONLY : ERROR_STOP
    USE State_Grid_Mod, ONLY : GrdState
!
! !INPUT PARAMETERS:
!
    REAL*4,         INTENT(IN)  :: LAT, LON
    TYPE(GrdState), INTENT(IN)  :: State_Grid  ! Grid State object
!
! !RETURN VALUE:
!
    INTEGER                     :: IIJJ(2)
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
      I0 = State_Grid%XMinOffset
      J0 = State_Grid%YMinOffset

      TLON = INT( ( LON + 180e+0_fp ) / State_Grid%DX + 1.5e+0_fp )
      TLAT = INT( ( LAT +  90e+0_fp ) / State_Grid%DY + 1.5e+0_fp )

      IF ( State_Grid%NestedGrid ) THEN
         TLON = TLON - I0
         TLAT = TLAT - J0
         IF ( TLAT < 1 .or. TLAT > State_Grid%NX ) THEN
            CALL ERROR_STOP('Beyond the nested window', 'GET_IJ')
         ENDIF
      ELSE
         IF ( TLON > State_Grid%NX ) TLON = TLON - State_Grid%NX

         ! Check for impossible values
         IF ( TLON > State_Grid%NX .or. TLAT > State_Grid%NY .or. &
              TLON < 1     .or. TLAT < 1          ) THEN
            CALL ERROR_STOP('Error finding grid box', 'GET_IJ')
         ENDIF
      ENDIF

      IIJJ(1) = TLON
      IIJJ(2) = TLAT

    END FUNCTION GET_IJ
!EOC
END MODULE GC_Grid_Mod
