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
#endif
  PUBLIC  :: GET_IJ
  PUBLIC  :: SetGridFromCtr
#ifdef MODEL_WRF
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
       State_Grid%GlobalXMid(I,J) = ( State_Grid%DX * (I-1) ) - 180e+0_fp

       !--------------------------------
       ! Latitude centers [degrees]
       !--------------------------------
       State_Grid%GlobalYMid(I,J) = ( State_Grid%DY * (J-1) ) - 90e+0_fp

       ! If using half-sized polar boxes, multiply DY by 1/4 at poles
       IF ( State_Grid%HalfPolar .and. J == 1)  THEN
          ! South Pole
          State_Grid%GlobalYMid(I,J) = -90e+0_fp + (0.25e+0_fp * State_Grid%DY)
       ENDIF
       IF ( State_Grid%HalfPolar .and. J == State_Grid%GlobalNY ) THEN
          ! North Pole
          State_Grid%GlobalYMid(I,J) = +90e+0_fp - (0.25e+0_fp * State_Grid%DY)
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
       State_Grid%XMid(I,J) = ( State_Grid%DX * IG ) - 180e+0_fp

       !--------------------------------
       ! Longitude edges [degrees]
       !--------------------------------
       State_Grid%XEdge(I,J) = State_Grid%XMid(I,J) - &
                             ( State_Grid%DX * 0.5e+0_fp )

       ! Compute the last longitude edge
       IF ( I == State_Grid%NX ) THEN
          State_Grid%XEdge(I+1,J) = State_Grid%XEdge(I,J) + State_Grid%DX
       ENDIF

       !--------------------------------
       ! Latitude centers [degrees]
       !--------------------------------
       State_Grid%YMid(I,J) = ( State_Grid%DY * JG ) - 90e+0_fp

       ! If using half-sized polar boxes on a global grid,
       ! multiply DY by 1/4 at poles
       IF ( State_Grid%HalfPolar .and. .not. State_Grid%NestedGrid .and. &
            J == 1 ) THEN
          ! South Pole
          State_Grid%YMid(I,J) = -90e+0_fp + ( 0.25e+0_fp * State_Grid%DY )
       ENDIF
       IF ( State_Grid%HalfPolar .and. .not. State_Grid%NestedGrid .and. &
            J == State_Grid%NY ) THEN
          ! North Pole
          State_Grid%YMid(I,J) = +90e+0_fp - ( 0.25e+0_fp * State_Grid%DY )
       ENDIF

       !--------------------------------
       ! Latitude centers [radians]
       !--------------------------------
       State_Grid%YMid_R(I,J) = ( PI_180 * State_Grid%YMid(I,J)  )

       !--------------------------------
       ! Latitude edges [degrees]
       !--------------------------------
       State_Grid%YEdge(I,J) = State_Grid%YMid(I,J) - &
                               ( State_Grid%DY * 0.5e+0_fp )

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

       ! Compute last latitude edge
       IF ( J == State_Grid%NY ) THEN

          ! Test for North Pole if using global grid
          IF ( .not. State_Grid%NestedGrid ) THEN

             ! Force the northern edge of grid boxes along the NORTH POLE to
             ! be +90 degrees latitude
             State_Grid%YEdge(I,J+1) = +90e+0_fp

             ! Make adjustment for second-to-last latitude edge
             State_Grid%YEdge(I,J  ) = State_Grid%YEdge(I,J+1) - &
                                     ( State_Grid%DY * 0.5e+0_fp )
             State_Grid%YEdge_R(I,J  ) = State_Grid%YEdge(I,J  ) * PI_180
             YEDGE_VAL = State_Grid%YEdge_R(I,J)
             YSIN_VAL  = SIN( YEDGE_VAL )
             State_Grid%YSIN(I,J) = YSIN_VAL
          ELSE

             !----------------------------------------------------------------
             !                %%%%% FOR NESTED GRIDS ONLY %%%%%
             !----------------------------------------------------------------
             State_Grid%YEdge(I,J+1) = State_Grid%YEdge(I,J) + State_Grid%DY

          ENDIF

          ! Last latitude edge [radians]
          State_Grid%YEdge_R(I,J+1) = State_Grid%YEdge(I,J+1) * PI_180

          ! Also compute sine of last two latitude edges! (ckeller, 02/13/12)
          YEDGE_VAL = State_Grid%YEdge_R(I,J+1)
          YSIN_VAL  = SIN( YEDGE_VAL )
          State_Grid%YSIN(I,J+1) = YSIN_VAL

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

  END SUBROUTINE Compute_Grid
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
       IF ( I == State_Grid%NX ) THEN
          State_Grid%XEdge(I+1,J) = State_Grid%XMid(I,J) + &
             ( ( State_Grid%XMid(I,J) - State_Grid%XMid(I-1,J) ) / 2.0_f4 )
       ENDIF
       IF ( J == State_Grid%NY ) THEN
          State_Grid%YEdge(I,J+1) = State_Grid%YMid(I,J) + &
             ( ( State_Grid%YMid(I,J) - State_Grid%YMid(I,J-1) ) / 2.0_f4 )
       ENDIF

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
#ifdef MODEL_WRF
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: SetGridFromCtrEdges
!
! !DESCRIPTION: Subroutine SetGridFromCtrEdges sets the grid based upon the
!  passed mid-points and edge-points given an external grid. This interface
!  is primarily used for GEOS-Chem to interface with the WRF model.
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
