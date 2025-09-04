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
#if !(defined( EXTERNAL_GRID ) || defined( EXTERNAL_FORCING ))
  PUBLIC  :: Compute_Grid
  PUBLIC  :: Compute_Scaled_Grid
#endif
  PUBLIC  :: GET_IJ
  PUBLIC  :: SetGridFromCtr
#if defined ( MODEL_WRF ) || defined( MODEL_CESM )
  PUBLIC  :: SetGridFromCtrEdges
#endif
!
! !REVISION HISTORY:
!  23 Feb 2012 - R. Yantosca - Initial version, based on grid_mod.F90
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
CONTAINS
!EOC
#if !(defined( EXTERNAL_GRID ) || defined( EXTERNAL_FORCING ))
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Compute_Grid
!
! !DESCRIPTION: Subroutine COMPUTE\_GRID initializes the longitude, latitude,
!  and surface area arrays.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Compute_Grid( Input_Opt, State_Grid, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod, ONLY : OptInput
!
! !INPUT PARAMETERS:
!
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
!
! !REVISION HISTORY:
!  22 May 2019 - M. Sulprizio- Initial version: Consolidated Compute_Grid and
!                              DoGridComputation into single routine that
!                              computes fields in State_Grid.
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    LOGICAL  :: is_GMAO
    INTEGER  :: I, J, L, IG, JG
    REAL(fp) :: YEDGE_VAL, YSIN_VAL
    REAL(fp) :: SIN_N, SIN_S, SIN_DIFF

    ! Strings
    CHARACTER(LEN=255) :: ErrMsg
    CHARACTER(LEN=255) :: ThisLoc

    !======================================================================
    ! Initialization
    !======================================================================

    ! Assume success
    RC = GC_SUCCESS
    ThisLoc = 'Compute_Grid (gc_grid_mod.F90)'

    !======================================================================
    ! Vertical Grid
    !======================================================================

    ! Both GEOS-FP and MERRA-2 have native vertical resolution of 72 levels
    State_Grid%NativeNZ = 72

    ! Hardcode maximum number of levels below tropopause and stratopause
    IF ( State_Grid%NZ == 47 ) THEN
       State_Grid%MaxTropLev  = 38
       State_Grid%MaxStratLev = 44
    ELSE IF ( State_Grid%NZ == 72 ) THEN
       State_Grid%MaxTropLev  = 40
       State_Grid%MaxStratLev = 59
    ELSE IF ( State_Grid%NZ == 40 ) THEN
       State_Grid%NativeNZ = 40
       State_Grid%MaxTropLev  = 28
       State_Grid%MaxStratLev = 40
    ELSE IF ( State_Grid%NZ == 74 ) THEN
       State_Grid%NativeNZ = 102
       State_Grid%MaxTropLev  = 60
       State_Grid%MaxStratLev = 72
    ELSE IF ( State_Grid%NZ == 102 ) THEN
       State_Grid%NativeNZ = 102
       State_Grid%MaxTropLev  = 60
       State_Grid%MaxStratLev = 91
    ELSE
       ErrMsg = 'State_Grid%GridRes = ' // Trim( State_Grid%GridRes)// &
                ' does not have MaxTropLev and MaxStratLev defined.'// &
                ' Please add these definitions in gc_grid_mod.F90.'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    ! Set maximum number of levels in the chemistry grid
    State_Grid%MaxChemLev  = State_Grid%MaxStratLev

    ! Set a flag to denote if this is a GMAO met field
    ! (which has half-sized polar grid boxes)
    SELECT CASE( TRIM( Input_Opt%MetField ) )
       CASE( 'MERRA2', 'GEOSFP' )
          is_GMAO = .TRUE.
       CASE DEFAULT
          is_GMAO = .FALSE.
    END SELECT

    !======================================================================
    ! Global Horizontal Grid
    !
    ! First, we need to compute the XMid and YMid values on the global
    ! grid at the specified resolution to that we can compute X and Y
    ! offsets.
    !======================================================================

    ! Compute number of grid boxes on global grid
    State_Grid%GlobalNX =   360.0_fp / State_Grid%DX
    IF ( State_Grid%HalfPolar .or. is_GMAO ) THEN
       State_Grid%GlobalNY = ( 180.0_fp / State_Grid%DY ) + 1
    ELSE
       State_Grid%GlobalNY = ( 180.0_fp / State_Grid%DY )
    ENDIF

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

    ALLOCATE( State_Grid%GlobalXEdge(State_Grid%GlobalNX+1,State_Grid%GlobalNY),&
              STAT=RC )
    CALL GC_CheckVar( 'State_Grid%GlobalXEdge', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Grid%GlobalXEdge = 0e+0_fp

    ALLOCATE( State_Grid%GlobalYEdge(State_Grid%GlobalNX,State_Grid%GlobalNY+1),&
              STAT=RC )
    CALL GC_CheckVar( 'State_Grid%GlobalYEdge', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Grid%GlobalYEdge = 0e+0_fp

    !======================================================================
    ! Global Horizontal Grid: Edges
    !======================================================================
    DO J = 1, State_Grid%GlobalNY
    DO I = 1, State_Grid%GlobalNX

       ! Latitude edges
       IF ( J == 1 ) THEN

          ! South Pole
          State_Grid%GlobalYEdge(I,J) = -90.0_fp

       ELSE

          IF ( State_Grid%HalfPolar .or. is_GMAO ) THEN
             State_Grid%GlobalYEdge(I,J) = -90.0_fp                          &
                                         - ( State_Grid%DY * 0.5_fp    )     &
                                         + ( State_Grid%DY * ( J - 1 ) )
          ELSE
             State_Grid%GlobalYEdge(I,J) = -90.0_fp                          &
                                         + ( State_Grid%DY * ( J - 1 ) )
          ENDIF

          IF ( J == State_Grid%GlobalNY ) THEN

             ! North pole
             State_Grid%GlobalYEdge(I,J+1) = +90.0_fp

          ENDIF

       ENDIF

       ! Longitude edges
       IF ( I == 1 ) THEN

          IF ( State_Grid%Center180 ) THEN
             State_Grid%GlobalXEdge(I,J) = -180.0_fp                         &
                                         - ( State_Grid%DX * 0.5_fp )
          ELSE
             State_Grid%GlobalXEdge(I,J) = -180.0_fp
          ENDIF

       ELSE

          State_Grid%GlobalXEdge(I,J) = State_Grid%GlobalXEdge(1,J)          &
                                      + ( State_Grid%DX * ( I - 1 ) )

          IF ( I == State_Grid%GlobalNX ) THEN
             State_Grid%GlobalXEdge(I+1,J) = State_Grid%GlobalXEdge(1,J)     &
                                           + 360.0_fp
          ENDIF

       ENDIF

    ENDDO
    ENDDO

    !======================================================================
    ! Global Horizontal Grid: Centers
    !======================================================================
    DO J = 1, State_Grid%GlobalNY
    DO I = 1, State_Grid%GlobalNX

       ! Latitude centers
       State_Grid%GlobalYMid(I,J) = 0.5_fp *                                 &
            ( State_Grid%GlobalYEdge(I,J) + State_Grid%GlobalYEdge(I,J+1) )

       ! Longitude centers
       State_Grid%GlobalXMid(I,J) = 0.5_fp *                                 &
            ( State_Grid%GlobalXEdge(I,J) + State_Grid%GlobalXEdge(I+1,J) )

    ENDDO
    ENDDO

    !======================================================================
    ! User-defined Horizontal Grid
    !======================================================================

    ! Determine X offsets based on global grid
    DO I = 1, State_Grid%GlobalNX
       IF ( State_Grid%GlobalXMid(I,1) >= State_Grid%XMin ) THEN
          State_Grid%XMinOffset = I - 1
          EXIT
       ENDIF
    ENDDO
    DO I = State_Grid%GlobalNX, 1, -1
       IF ( State_Grid%GlobalXMid(I,1) <= State_Grid%XMax ) THEN
          State_Grid%XMaxOffset = I - 1
          EXIT
       ENDIF
    ENDDO

    ! Determine Y offsets based on global grid
    DO J = 1, State_Grid%GlobalNY
       IF ( State_Grid%GlobalYMid(1,J) >= State_Grid%YMin ) THEN
          State_Grid%YMinOffset = J - 1
          EXIT
       ENDIF
    ENDDO
    DO J = State_Grid%GlobalNY, 1, -1
       IF ( State_Grid%GlobalYMid(1,J) <= State_Grid%YMax ) THEN
          State_Grid%YMaxOffset = J - 1
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
       IG = I + State_Grid%XMinOffset
       JG = J + State_Grid%YMinOffset

       !--------------------------------
       ! Longitude centers [degrees]
       !--------------------------------
       State_Grid%XMid(I,J) = State_Grid%GlobalXMid( IG, JG )

       !--------------------------------
       ! Longitude edges [degrees]
       !--------------------------------
       State_Grid%XEdge(I,J) = State_Grid%GlobalXEdge( IG, JG )

       ! Compute the last longitude edge
       IF ( I == State_Grid%NX ) THEN
          State_Grid%XEdge(I+1,J) = State_Grid%GlobalXEdge(IG+1,J)
       ENDIF

       !--------------------------------
       ! Latitude centers [degrees]
       !--------------------------------
       State_Grid%YMid(I,J) = State_Grid%GlobalYMid( IG, JG )

       !--------------------------------
       ! Latitude centers [radians]
       !--------------------------------
       State_Grid%YMid_R(I,J) = ( PI_180 * State_Grid%YMid(I,J)  )

       !-----------------------------------------
       ! Latitude edges [degrees; rad; sin(lat)]
       !-----------------------------------------
       State_Grid%YEdge(I,J)   = State_Grid%GlobalYEdge( IG, JG )
       State_Grid%YEdge_R(I,J) = ( PI_180  * State_Grid%YEdge(I,J) )
       State_Grid%YSIN(I,J)    = SIN( State_Grid%YEdge_R(I,J) )

       ! Compute the last latitude edge
       IF ( J == State_Grid%NY ) THEN
          State_Grid%YEdge(I,J+1)   = State_Grid%GlobalYEdge(I,JG+1)
          State_Grid%YEdge_R(I,J+1) = ( PI_180  * State_Grid%YEdge(I,J+1) )
          State_Grid%YSIN(I,J+1)    = SIN( State_Grid%YEdge_R(I,J+1) )
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

       ! Grid box surface areas [m2]
       State_Grid%Area_M2(I,J) = ( State_Grid%DX * PI_180 ) * &
                                   ( Re**2 ) * SIN_DIFF

    ENDDO
    ENDDO

    ! Return successfully
    RC = GC_SUCCESS

    !======================================================================
    ! Echo info to stdout
    !======================================================================
    IF ( Input_Opt%amIRoot ) THEN
       WRITE( 6, '(a)' )
       WRITE( 6, '(''%%%%%%%%%%%%%%% GLOBAL GRID %%%%%%%%%%%%%%%'')' )
       WRITE( 6, '(a)' )
       WRITE( 6, '(''Grid box longitude centers [degrees]: '')' )
       WRITE( 6, '(7(f10.5,1x))' ) ( State_Grid%GlobalXMid(I,1), &
                                    I=1,State_Grid%GlobalNX )
       WRITE( 6, '(a)' )
       WRITE( 6, '(''Grid box latitude centers [degrees]: '')' )
       WRITE( 6, '(7(f10.5,1x))' ) ( State_Grid%GlobalYMid(1,J), &
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
       WRITE( 6, '(7(f10.5,1x))') ( State_Grid%XMid(I,1), I=1,State_Grid%NX )
       WRITE( 6, '(a)' )
       WRITE( 6, '(''Grid box longitude edges [degrees]: '')' )
       WRITE( 6, '(7(f10.5,1x))') ( State_Grid%XEdge(I,1), I=1,State_Grid%NX+1 )
       WRITE( 6, '(a)' )
       WRITE( 6, '(''Grid box latitude centers [degrees]: '')' )
       WRITE( 6, '(7(f10.5,1x))') ( State_Grid%YMid(1,J), J=1,State_Grid%NY )
       WRITE( 6, '(a)' )
       WRITE( 6, '(''Grid box latitude edges [degrees]: '')' )
       WRITE( 6, '(7(f10.5,1x))') ( State_Grid%YEdge(1,J), J=1,State_Grid%NY+1 )
       WRITE( 6, '(a)' )
       WRITE( 6, '(''SIN( grid box latitude edges )'')' )
       WRITE( 6, '(7(f10.5,1x))') ( State_Grid%YSIN(1,J), J=1,State_Grid%NY+1 )
    ENDIF

  END SUBROUTINE Compute_Grid
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Compute_Scaled_Grid
!
! !DESCRIPTION: Subroutine COMPUTE\_SCALED\_GRID populates a secondary Grid
!  State object ("Destination") by performing a linear scaling refinement
!  of the primary ("Source") Grid. e.g. a scale of 2 will yield
!  2 x 2.5 -> 1 x 1.25.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Compute_Scaled_Grid( Input_Opt, State_Grid, State_Grid_Dst,      &
                                  XScale,    YScale,      RC                 )
!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod,  ONLY : OptInput
    USE State_Grid_Mod, ONLY : Init_State_Grid, Allocate_State_Grid
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt       ! Input Options
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(GrdState), INTENT(IN   ) :: State_Grid      ! Grid State object, orig
    TYPE(GrdState), INTENT(INOUT) :: State_Grid_Dst  ! Grid State object, scaled
    INTEGER,        INTENT(IN   ) :: XScale          ! Long dim scale (>=1)
    INTEGER,        INTENT(IN   ) :: YScale          ! Lat  dim scale (>=1)
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC                ! Success/failure?
!
! !REMARKS:
!  Unlike Compute_Grid, Compute_Scaled_Grid takes care of allocating
!  the State_Grid_Dst derived type object.
!  Only works with rectilinear lat-lon (assumptions given)
!
! !REVISION HISTORY:
!  15 Jun 2020 - H.P. Lin    - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER                       :: I, J, L, II, JJ
    REAL(fp)                      :: delta

    ! Strings
    CHARACTER(LEN=255) :: ErrMsg
    CHARACTER(LEN=255) :: ThisLoc

    ! Code begins here!
    RC = GC_SUCCESS
    ErrMsg  = ''
    ThisLoc = ' -> at Compute_Scaled_Grid (in module GeosUtil/gc_grid_mod.F90)'

    ! Verify the scales are correct
    IF ( XScale .lt. 1 .or. YScale .lt. 1 ) THEN
      RC = GC_FAILURE
      ErrMsg = 'Cannot scale Compute_Scaled_Grid to a negative factor!'
      CALL GC_Error( ErrMsg, RC, ThisLoc )
    ENDIF

    ! Initialize the destination object for safety sake
    CALL Init_State_Grid ( Input_Opt, State_Grid_Dst, RC )

    ! First, copy generic properties of the grid state information
    State_Grid_Dst%HalfPolar  = State_Grid%HalfPolar
    State_Grid_Dst%NestedGrid = State_Grid%NestedGrid

    ! The vertical direction is NOT scaled!
    State_Grid_Dst%NZ   = State_Grid%NZ

    ! Compute the scaling. First NX, NY will be x number the scaled grid boxes
    State_Grid_Dst%NX   = State_Grid%NX * XScale
    State_Grid_Dst%NY   = State_Grid%NY * YScale

    State_Grid_Dst%XMin = State_Grid%XMin
    State_Grid_Dst%XMax = State_Grid%XMax
    State_Grid_Dst%YMin = State_Grid%YMin
    State_Grid_Dst%YMax = State_Grid%YMax

    ! The delta is also divided. Note however that the polar boxes
    ! will also be evenly divided so you may have more than 1 1/2-sized polar box.
    !
    ! LATITUDE EDGES example:
    !   [deg_Src]                   [deg_Dst]
    !    90.0                        90.0
    !                                89.5*   * additional polar edge
    !    89.0                        89.0
    !                                88.0
    !    87.0                        87.0
    !     ...  in 2x2.5 grid          ...  in 1xLon grid (2x scaling)
    State_Grid_Dst%DX   = State_Grid%DX / XScale
    State_Grid_Dst%DY   = State_Grid%DY / YScale

    ! Now, allocate the derived type object as we have NX, NY and NZ now!
    CALL Allocate_State_Grid ( Input_Opt, State_Grid_Dst, RC )

    IF ( Input_Opt%amIRoot ) THEN
      WRITE(6,*) "GC_GRID_MOD: Scaled grid NX, NY, DX, DY", State_Grid_Dst%NX, State_Grid_Dst%NY, &
                 State_Grid_Dst%DX, State_Grid_Dst%DY
    ENDIF

    ! Now, fill in the appropriate fields. Note that we do not fill
    ! all available fields in State_Grid, only the ones that are used
    ! by HEMCO for now. This is just to save development time (hplin, 6/15/20) <-- person to blame
    ! XMid, YMid, XEdge, YEdge, YSin, Area_M2

    ! First, construct the edges. The end goals is that the edges
    ! of the destination grid will be equal as the original grid, plus the refined edges.
    ! So, fill in the original edges, then interpolate
    DO I = 1, State_Grid%NX+1         ! For each original longitude EDGE...
      II = 1 + XScale * (I - 1)
      State_Grid_Dst%XEdge(II,:) = State_Grid%XEdge(I,1)

      ! Fill in the gaps above, except if there isnt one
      IF ( II .lt. State_Grid_Dst%NX ) THEN
        Delta = (State_Grid%XEdge(I+1,1) - State_Grid%XEdge(I,1)) / XScale
        DO JJ = 1, XScale
          State_Grid_Dst%XEdge(II+JJ,:) = State_Grid_Dst%XEdge(II,:) + Delta * JJ
        ENDDO
      ENDIF
    ENDDO

    ! For latitude, do the same
    DO I = 1, State_Grid%NY+1
      II = 1 + YScale * (I - 1)
      State_Grid_Dst%YEdge(:,II) = State_Grid%YEdge(1,I)

      ! Fill in the gaps above, except if there isnt one
      IF ( II .lt. State_Grid_Dst%NY ) THEN
        Delta = (State_Grid%YEdge(1,I+1) - State_Grid%YEdge(1,I)) / YScale
        DO JJ = 1, YScale
          State_Grid_Dst%YEdge(:,II+JJ) = State_Grid_Dst%YEdge(:,II) + Delta * JJ
        ENDDO
      ENDIF
    ENDDO

    DO J = 1, State_Grid_Dst%NY+1   ! Including edges!
    DO I = 1, State_Grid_Dst%NX
      State_Grid_Dst%YEdge_R(I,J) = ( PI_180 * State_Grid_Dst%YEdge(I,J) )
      State_Grid_Dst%YSIN(I,J) = SIN( State_Grid_Dst%YEdge_R(I,J) )
    ENDDO
    ENDDO

    ! Compute surface areas and midpoints
    DO J = 1, State_Grid_Dst%NY
    DO I = 1, State_Grid_Dst%NX
      State_Grid_Dst%XMid(I,J) = (State_Grid_Dst%XEdge(I,J) + State_Grid_Dst%XEdge(I+1,J)) / 2.0
      State_Grid_Dst%YMid(I,J) = (State_Grid_Dst%YEdge(I,J) + State_Grid_Dst%YEdge(I,J+1)) / 2.0

      ! Grid box surface areas [m2]
      State_Grid_Dst%Area_M2(I,J) = ( State_Grid_Dst%DX * PI_180 ) * ( Re**2 ) * &
                                    ( SIN( State_Grid_Dst%YEdge_R(I,J+1) ) - SIN( State_Grid_Dst%YEdge_R(I,J  ) ))
    ENDDO
    ENDDO

    IF ( Input_Opt%amIRoot ) THEN
      WRITE( 6, '(''%%%%%%%%%%%%%%% SCALED (HEMCO) GRID %%%%%%%%%%%%%%%'')' )
      WRITE( 6, '(''Grid box longitude centers [degrees]: '')' )
      WRITE( 6, '(7(f10.5,1x))' ) ( State_Grid_Dst%XMid(I,1), I=1,State_Grid_Dst%NX )
      WRITE( 6, '(a)' )
      WRITE( 6, '(''Grid box longitude edges [degrees]: '')' )
      WRITE( 6, '(7(f10.5,1x))' ) ( State_Grid_Dst%XEdge(I,1), I=1,State_Grid_Dst%NX+1 )
      WRITE( 6, '(a)' )
      WRITE( 6, '(''Grid box latitude centers [degrees]: '')' )
      WRITE( 6, '(7(f10.5,1x))' ) ( State_Grid_Dst%YMid(1,J), J=1,State_Grid_Dst%NY )
      WRITE( 6, '(a)' )
      WRITE( 6, '(''Grid box latitude edges [degrees]: '')' )
      WRITE( 6, '(7(f10.5,1x))' ) ( State_Grid_Dst%YEdge(1,J), J=1,State_Grid_Dst%NY+1 )
      WRITE( 6, '(a)' )
      WRITE( 6, '(''SIN( grid box latitude edges )'')' )
      WRITE( 6, '(7(f10.5,1x))' ) ( State_Grid_Dst%YSIN(1,J), J=1,State_Grid_Dst%NY+1 )
    ENDIF
  END SUBROUTINE Compute_Scaled_Grid
!EOC
#endif
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
  SUBROUTINE SetGridFromCtr( Input_Opt, State_Grid, lonCtr, latCtr, RC )
!
! USES
!
    USE ErrCode_Mod
    USE Input_Opt_Mod, ONLY : OptInput
    USE Roundoff_Mod
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt      ! Input Options object
    REAL(f4),       INTENT(IN)    :: lonCtr(:,:)    ! Lon ctrs [rad]
    REAL(f4),       INTENT(IN)    :: latCtr(:,:)    ! Lat ctrs [rad]
    TYPE(GrdState), INTENT(IN)    :: State_Grid     ! Grid State object
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(INOUT) :: RC
!
! !REVISION HISTORY:
!  02 Jan 2014 - C. Keller   - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER            :: I, J
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

    ! Loop over all grid boxes
    DO J = 1, State_Grid%NY
    DO I = 1, State_Grid%NX

       ! Mid points: get directly from passed value
       State_Grid%XMid(I,J)   = RoundOff( lonCtr(I,J) / PI_180, 4 )
       State_Grid%YMid(I,J)   = RoundOff( latCtr(I,J) / PI_180, 4 )
       State_Grid%YMid_R(I,J) = State_Grid%YMid(I,J) * PI_180

       ! Edges: approximate from neighboring mid points.
#if defined( MODEL_CESM )
       ! If using CESM, prevent out-of-bound errors when running with
       ! NX or NY equal to 1
       IF ( State_Grid%NX > 1 ) THEN
#endif
           IF ( I == 1 ) THEN
              TMP = RoundOff( lonCtr(I+1,J) / PI_180, 4 )
              State_Grid%XEdge(I,J) = State_Grid%XMid(I,J) - &
                                      ( ( TMP - State_Grid%XMid(I,J) ) / 2.0_f4 )
           ELSE
              State_Grid%XEdge(I,J) = ( State_Grid%XMid(I,J) + &
                                          State_Grid%XMid(I-1,J) ) / 2.0_f4
           ENDIF
#if defined( MODEL_CESM )
       ELSE
           State_Grid%XEdge(I,J) = State_Grid%XMid(I,J)
       ENDIF

       IF ( State_Grid%NY > 1 ) THEN
#endif
           IF ( J == 1 ) THEN
              TMP = RoundOff( latCtr(I,J+1) / PI_180, 4 )
              State_Grid%YEdge(I,J) = State_Grid%YMid(I,J) - &
                                      ( ( TMP - State_Grid%YMid(I,J) ) / 2.0_f4 )
           ELSE
              State_Grid%YEdge(I,J) = ( State_Grid%YMid(I,J) + &
                                          State_Grid%YMid(I,J-1) ) / 2.0_f4
           ENDIF
#if defined( MODEL_CESM )
       ELSE
           State_Grid%YEdge(I,J) = State_Grid%YMid(I,J)
       ENDIF
#endif

       ! Special treatment at uppermost edge
#if defined( MODEL_CESM )
       IF ( State_Grid%NX > 1 ) THEN
#endif
           IF ( I == State_Grid%NX ) THEN
              State_Grid%XEdge(I+1,J) = State_Grid%XMid(I,J) + &
                 ( ( State_Grid%XMid(I,J) - State_Grid%XMid(I-1,J) ) / 2.0_f4 )
           ENDIF
#if defined( MODEL_CESM )
       ENDIF
       IF ( State_Grid%NY > 1 ) THEN
#endif
           IF ( J == State_Grid%NY ) THEN
              State_Grid%YEdge(I,J+1) = State_Grid%YMid(I,J) + &
                 ( ( State_Grid%YMid(I,J) - State_Grid%YMid(I,J-1) ) / 2.0_f4 )
           ENDIF
#if defined( MODEL_CESM )
       ENDIF
#endif

       ! Special quantities directly derived from State_Grid%YEdge
       State_Grid%YEdge_R(I,J) = State_Grid%YEdge(I,J) * PI_180
       YEDGE_VAL               = State_Grid%YEdge_R(I,J) ! Lat edge[radians]
       YSIN_VAL                = SIN( YEDGE_VAL)         ! SIN( lat edge )
       State_Grid%YSIN(I,J)    = YSIN_VAL                ! Store in YSIN

    ENDDO
    ENDDO

    ! Return w/ success
    RC = GC_SUCCESS

  END SUBROUTINE SetGridFromCtr
!EOC
#if defined ( MODEL_WRF ) || defined( MODEL_CESM )
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: SetGridFromCtrEdges
!
! !DESCRIPTION: Subroutine SetGridFromCtrEdges sets the grid based upon the
!  passed mid-points and edge-points given an external grid. This interface
!  is primarily used for GEOS-Chem to interface with the WRF and CESM models.
!\\
!\\
! This routine does not update the grid box areas (AREA\_M2) of grid\_mod.F90.
! These need to be updated manually from State\_Grid%AREA\_M2 to maintain
! consistency with the GEOS-Chem interface to GEOS-5.
! !INTERFACE:
!
  SUBROUTINE SetGridFromCtrEdges( Input_Opt, State_Grid, lonCtr, latCtr, &
                                  lonEdge, latEdge, RC )
!
! USES
!
    USE ErrCode_Mod
    USE Input_Opt_Mod, ONLY : OptInput
    USE Roundoff_Mod
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt      ! Input Options object
    REAL(f4),       INTENT(IN)    :: lonCtr (:,:)   ! Lon ctrs [rad]
    REAL(f4),       INTENT(IN)    :: latCtr (:,:)   ! Lat ctrs [rad]
    REAL(f4),       INTENT(IN)    :: lonEdge(:,:)   ! Lon edges [rad]
    REAL(f4),       INTENT(IN)    :: latEdge(:,:)   ! Lat edges [rad]
    TYPE(GrdState), INTENT(IN)    :: State_Grid     ! Grid State object
!!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(INOUT) :: RC
!
! !REVISION HISTORY:
!  11 Nov 2018 - H.P. Lin    - Initial version based on SetGridFromCtr
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER            :: I, J
    REAL(fp)           :: YEDGE_VAL, YSIN_VAL, TMP

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

    ! Loop over all grid boxes
    DO J = 1, State_Grid%NY
    DO I = 1, State_Grid%NX

       ! Mid points: get directly from passed value
       State_Grid%XMid(I,J)      = RoundOff( lonCtr(I,J) / PI_180, 4 )
       State_Grid%YMid(I,J)      = RoundOff( latCtr(I,J) / PI_180, 4 )
       State_Grid%YMid_R(I,J)    = State_Grid%YMid(I,J) * PI_180

       ! Edges: get directly from passed value
       State_Grid%XEdge(I,J)     = RoundOff( lonEdge(I,J) / PI_180, 4 )
       State_Grid%YEdge(I,J)     = RoundOff( latEdge(I,J) / PI_180, 4 )

       ! Special treatment at uppermost edge
       IF ( I == State_Grid%NX ) THEN
          State_Grid%XEdge(I+1,J) = RoundOff( lonEdge(I+1,J) / PI_180, 4 )
       ENDIF
       IF ( J == State_Grid%NY ) THEN
          State_Grid%YEdge(I,J+1)   = RoundOff( latEdge(I,J+1) / PI_180, 4 )
          State_Grid%YEdge_R(I,J+1) = State_Grid%YEdge(I,J+1) * PI_180
          State_Grid%YSIN(I,J+1)    = SIN( State_Grid%YEdge_R(I,J+1) )
       ENDIF

       ! Special quantities directly derived from State_Grid%YEdge
       State_Grid%YEdge_R(I,J) = State_Grid%YEdge(I,J) * PI_180
       YEDGE_VAL               = State_Grid%YEdge_R(I,J) ! Lat edge[radians]
       YSIN_VAL                = SIN( YEDGE_VAL )        ! SIN( lat edge )
       State_Grid%YSIN(I,J)    = YSIN_VAL                ! Store in YSIN

    ENDDO
    ENDDO

    ! Return w/ success
    RC = GC_SUCCESS

  END SUBROUTINE SetGridFromCtrEdges
#endif
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
!  See https://github.com/geoschem/geos-chem for complete history
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
         IF ( TLAT < 1 .or. TLAT > State_Grid%NY ) THEN
            CALL ERROR_STOP('Beyond the nested window', 'GET_IJ')
         ENDIF
      ELSE
         IF ( TLON > State_Grid%NX ) TLON = TLON - State_Grid%NX

         ! Check for impossible values
         IF ( TLON > State_Grid%NX .or. TLAT > State_Grid%NY .or. &
              TLON < 1             .or. TLAT < 1          ) THEN
            CALL ERROR_STOP('Error finding grid box', 'GET_IJ')
         ENDIF
      ENDIF

      IIJJ(1) = TLON
      IIJJ(2) = TLAT

    END FUNCTION GET_IJ
!EOC
END MODULE GC_Grid_Mod
