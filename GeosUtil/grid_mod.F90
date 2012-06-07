!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: grid_mod.F90
!
! !DESCRIPTION: Module GRID\_MOD contains variables and routines which are 
!  used to specify the parameters of a GEOS-Chem horizontal grid. Grid 
!  parameters are computed as 3D arrays, which are required for interfacing
!  with a GCM.
!\\  
!\\
! !INTERFACE: 
!
MODULE Grid_Mod
! 
! !USES:
!
  USE CMN_GCTM_Mod             ! Physical constants
  USE Error_Mod                ! Error-handling routines

  IMPLICIT NONE
# include "define.h"
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC  :: Cleanup_Grid
  PUBLIC  :: Compute_Grid
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

#if defined( DEVEL )
!      PUBLIC  :: AREA_M2 ! Permit setting this externally
      PUBLIC  :: YMID, XMID, YEDGE, XEDGE
#endif
!
! !REVISION HISTORY:
!  23 Feb 2012 - R. Yantosca - Initial version, based on grid_mod.F
!  01 Mar 2012 - R. Yantosca - Validated for nested grids
!  03 Apr 2012 - M. Payer    - Added ySin for map_a2a regrid (M. Cooper)
!EOP
!------------------------------------------------------------------------------
!BOC

  !=======================================================================
  ! MODULE VARIABLES:
  !
  ! IS_NESTED  : =T if we are using a nested-grid 
  ! I0         : Nested-grid offset in longitude (X) dimension
  ! J0         : Nested-grid offset in latitude  (Y) dimension
  ! XMID       : Array of grid-box lon centers   [degrees]
  ! XEDGE      : Array of grid-box lon edges     [degrees]
  ! YMID       : Array of grid-box lat centers   [degrees]
  ! YEDGE      : Array of grid-box lat edges     [degrees]
  ! YMID_R     : Array of grid-box lat centers   [radians]
  ! YEDGE_R    : Array of grid-box lat edges     [radians]
  ! AREA_M2    : Array of grid-box surface areas [m2     ]
  ! AREA_CM2   : Array of grid-box surface areas [cm2    ]
  !=======================================================================

  ! Scalars
  LOGICAL              :: IS_NESTED
  INTEGER              :: I0
  INTEGER              :: J0

  ! Arrays
  REAL*8,  ALLOCATABLE :: XMID     (:,:,:)
  REAL*8,  ALLOCATABLE :: XEDGE    (:,:,:)
  REAL*8,  ALLOCATABLE :: YMID     (:,:,:)
  REAL*8,  ALLOCATABLE :: YEDGE    (:,:,:)
  REAL*8,  ALLOCATABLE :: YSIN     (:,:,:)
  REAL*8,  ALLOCATABLE :: YMID_R   (:,:,:)
  REAL*8,  ALLOCATABLE :: YEDGE_R  (:,:,:)
  REAL*8,  ALLOCATABLE :: YMID_R_W (:,:,:)
  REAL*8,  ALLOCATABLE :: YEDGE_R_W(:,:,:)
  REAL*8,  ALLOCATABLE :: AREA_M2  (:,:,:)
  
CONTAINS
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: compute_grid
!
! !DESCRIPTION: Subroutine COMPUTE\_GRID initializes the longitude, 
!  latitude and surface area arrays. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Compute_Grid( I1, I2, J1, J2, JSP, JNP, L1, L2, DLON, DLAT )
!
! !INPUT PARAMETERS:
!
    INTEGER, INTENT(IN) :: I1,  I2                  ! Local CPU lon idx bounds
    INTEGER, INTENT(IN) :: J1,  J2                  ! Local CPU lat idx bounds
    INTEGER, INTENT(IN) :: JSP, JNP                 ! Polar lat indices
    INTEGER, INTENT(IN) :: L1,  L2                  ! Local CPU lev idx bounds
    REAL*8,  INTENT(IN) :: DLON(I1:I2,J1:J2,L1:L2)  ! Delta lon [degrees]
    REAL*8,  INTENT(IN) :: DLAT(I1:I2,J1:J2,L1:L2)  ! Delta lat [degrees]
!
! !REVISION HISTORY:
!  23 Feb 2012 - R. Yantosca - Initial version, based on grid_mod.F
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
! 
    ! Scalars
    INTEGER :: I, J, L

    ! Arrays
    REAL*8  :: IND_X( I1:I2+1 )
    REAL*8  :: IND_Y( J1:J2+1 ),  TMP

    !=================================================================
    ! Initialization
    !=================================================================

    ! Index array for longitudes
    DO I = I1, I2+1
       IND_X(I) = ( I + I0 - 1 ) * 1d0 
    ENDDO

    ! Index array for latitudes
    DO J = J1, J2+1
       IND_Y(J) = ( J + J0 - 1 ) * 1d0
    ENDDO

    !=================================================================
    ! Compute longitude and latitude arrays
    !=================================================================
    
    ! Loop over levels
    DO L = L1, L2
       
       !--------------------------------------------------------------
       ! Longitude center and edge arrays
       !--------------------------------------------------------------
       DO J = J1, J2
       DO I = I1, I2
             
          ! Longitude centers
          XMID(I,J,L)  = ( DLON(I,J,L) * IND_X(I) ) - 180d0
          
          ! Longitude edges
          XEDGE(I,J,L) = XMID(I,J,L) - ( DLON(I,J,L) * 0.5d0 )

          ! Compute the last longitude edge
          IF ( I == I2 ) THEN
             XEDGE(I2+1,J,L) = XEDGE(I2,J,L) + DLON(I2,J,L)
          ENDIF

       ENDDO
       ENDDO

       !--------------------------------------------------------------
       ! Latitude center and edge arrays
       !---------------------------------------==---------------------
       DO J = J1, J2
       DO I = I1, I2

          !%%%%%% LATITUDE CENTERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

          ! Lat centers (degrees)
          YMID(I,J,L)     = ( DLAT(I,J,L) * IND_Y(J) ) - 90d0
          
          ! Adjust for halfpolar boxes (global grids only)
          IF ( J == JSP ) THEN
             YMID(I,J,L)  = -90d0 + ( 0.5d0 * DLAT(I,J,L) )   ! S pole
          ELSE IF ( J == JNP ) THEN
             YMID(I,J,L)  = +90d0 - ( 0.5d0 * DLAT(I,J,L) )   ! N pole
          ENDIF

          ! Lat centers (radians)
          YMID_R(I,J,L)   = ( PI_180 * YMID(I,J,L)  )

#if defined( NESTED_CH ) || defined( NESTED_EU ) || defined( NESTED_NA ) || defined( SEAC4RS )

          ! Lat centers (radians), for nested grid window array
          YMID_R_W(I,J,L) = YMID_R(I,J,L)

          ! Compute YMID_R_W at edges of nested region
          IF ( J == J1 ) THEN
             YMID_R_W(I,J1-1,1) = YMID_R(I,J1,L) - ( DLAT(I,J1,L) * PI_180 )
          ELSE IF ( J == J2 ) THEN
             YMID_R_W(I,J2+1,1) = YMID_R(I,J2,L) + ( DLAT(I,J2,L) * PI_180 )
          ENDIF
                    
#endif

          !%%%%%% LATITUDE EDGES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

          ! Lat edges (degrees and radians)
          YEDGE(I,J,L)    = YMID(I,J,L) - ( DLAT(I,J,L) * 0.5d0 )
            
          ! Adjust for halfpolar boxes
          IF ( J == JSP ) THEN
             YEDGE(I,J,L) = -90d0                             ! S pole
          ENDIF
          
          ! Lat edges (radians)
          YEDGE_R(I,J,L)  = ( PI_180  * YEDGE(I,J,L) )

          ! mjc - Compute sine of latitude edges (needed for map_a2a regrid)
          YSIN(I,J,L) = SIN ( YEDGE_R(I,J,L) )

       ENDDO
       ENDDO

       !%%%%%% LATITUDE EDGES (the last edge) %%%%%%%%%%%%%%%%%%%%%%%%
       
       ! Test for North Pole
       IF ( J2 == JNP ) THEN

          ! North pole case (global grids only)
          DO I = I1, I2
             YEDGE  (I,JNP+1,L) = +90d0
             YEDGE_R(I,JNP+1,L) = YEDGE(I,JNP+1,L) * PI_180
          ENDDO
          
       ELSE
          
          ! No north pole (nested grids only)
          DO I = I1, I2 
             YEDGE  (I,J2+1,L)  = YEDGE(I,J2,L  ) + DLAT(I,J2,L)
             YEDGE_R(I,J2+1,L)  = YEDGE(I,J2+1,L) * PI_180
          ENDDO
       ENDIF

       !=================================================================
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
       !=================================================================  
       DO J = J1, J2
       DO I = I1, I2

          ! Grid box surface areas [m2]
          AREA_M2(I,J,L) = ( DLON(I,J,L) * PI_180 ) * Re**2                &
                         * ( SIN( YEDGE_R(I,J+1,L) ) - SIN( YEDGE_R(I,J,L) ) )

       ENDDO
       ENDDO

    ENDDO

    !=================================================================
    ! Echo info to stdout
    !=================================================================
    WRITE( 6, '(''Nested-Grid X-offset [boxes]:'', i4 )' ) I0
    WRITE( 6, '(''Nested-Grid Y-offset [boxes]:'', i4 )' ) J0
    WRITE( 6, '(a)' )
    WRITE( 6, '(''Grid box longitude centers [degrees]: '')' )
    WRITE( 6, '(8(f8.3,1x))' ) ( XMID(I,1,1),  I=I1,I2 )
    WRITE( 6, '(a)' )
    WRITE( 6, '(''Grid box longitude edges [degrees]: '')' )
    WRITE( 6, '(8(f8.3,1x))' ) ( XEDGE(I,1,1), I=I1,I2+1 )
    WRITE( 6, '(a)' )
    WRITE( 6, '(''Grid box latitude centers [degrees]: '')' )
    WRITE( 6, '(8(f8.3,1x))' ) ( YMID(1,J,1),  J=J1,J2 )
    WRITE( 6, '(a)' )
    WRITE( 6, '(''Grid box latitude edges [degrees]: '')' )
    WRITE( 6, '(8(f8.3,1x))' ) ( YEDGE(1,J,1), J=J1,J2+1 )

  END SUBROUTINE Compute_Grid
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!-----------------------------<------------------------------------------------
!BOP
!
! !IROUTINE: set_xoffset
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
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: set_yoffset
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
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_xoffset
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
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_xoffset
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
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_xmid
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
    REAL*8              :: X   ! Corresponding lon value @ grid box ctr
!
! !REVISION HISTORY:
!  24 Feb 2012 - R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
    X = XMID(I,J,L)

  END FUNCTION Get_xMid
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_xedge
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
    REAL*8              :: X   ! Corresponding lon value @ W edge of grid box
!
! !REVISION HISTORY:
!  24 Feb 2012 - R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
    X = XEDGE(I,J,L)

  END FUNCTION Get_xEdge
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_ymid
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
    REAL*8              :: Y   ! Latitude value at @ grid box ctr [degrees]
!
! !REVISION HISTORY:
!  24 Feb 2012 - R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
    Y = YMID(I,J,L)
      
  END FUNCTION Get_yMid
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_yedge
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
    REAL*8              :: Y   ! Latitude value @ S edge of grid box [degrees]
!
! !REVISION HISTORY:
!  24 Feb 2012 - R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
    Y = YEDGE(I,J,L)

  END FUNCTION Get_yEdge
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_ymid_r
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
    REAL*8              :: Y   ! Latitude value at @ grid box ctr [radians]
!
! !REVISION HISTORY:
!  24 Feb 2012 - R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
    Y = YMID_R(I,J,L)

  END FUNCTION Get_yMid_R
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_ymid_r_w
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
    REAL*8              :: Y   ! Latitude value at @ grid box ctr [radians]
!
! !REVISION HISTORY:
!  24 Feb 2012 - R. Yantosca - Initial version
!  01 Mar 2012 - R. Yantosca - Bracket with #ifdef for nested grids
!EOP
!------------------------------------------------------------------------------
!BOC

#if defined( NESTED_CH ) || defined( NESTED_EU ) || defined( NESTED_NA ) || defined( SEAC4RS )

    ! For nested grids, return the latitude center of the window
    ! region (in radians)
    Y = YMID_R_W(I,J,L)

#else

    ! Otherwise return a fake value
    Y = -1d0

#endif

  END FUNCTION Get_yMid_R_W
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_yedge_r
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
    REAL*8              :: Y   ! Latitude value @ S edge of grid box [radians]
!
! !REVISION HISTORY:
!  24 Feb 2012 - R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC

    Y = YEDGE_R(I,J,L)

  END FUNCTION Get_yEdge_R
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_ysin
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
    REAL*8              :: Y   ! Sine of Latitude value @ S edge of grid box
!
! !REVISION HISTORY:
!  03 Apr 2012 - M. Payer    -  Initial version (M. Cooper)
!EOP
!------------------------------------------------------------------------------
!BOC

      Y = YSIN(I,J,L)

      END FUNCTION Get_ySin
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_area_m2
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
    REAL*8              :: A   ! Grid box surface area [m2]
!
! !REVISION HISTORY:
!  24 Feb 2012 - R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
    A = AREA_M2(I,J,L)

  END FUNCTION Get_Area_m2
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_area_cm2
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
    REAL*8              :: A  ! Grid box surface area [cm2]
!
! !REVISION HISTORY:
!  24 Feb 2012 - R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
    A = AREA_M2(I,J,L) * 1d4
    
  END FUNCTION Get_Area_cm2
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_bounding_box
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
    REAL*8,  INTENT(IN)  :: COORDS(4)   ! (/LON_LL, LAT_LL, LON_UR, LAT_UR/)
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
       IF ( COORDS(1) >  XEDGE(I,  J,L)   .and.          &
            COORDS(1) <= XEDGE(I+1,J,L) ) INDICES(1) = I

         ! Locate index corresponding to upper-right longitude
       IF ( COORDS(3) >  XEDGE(I,  J,L)   .and.          &
            COORDS(3) <= XEDGE(I+1,J,L) ) INDICES(3) = I

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
       IF ( COORDS(2) >  YEDGE(I,J,  L)   .and.          &
            COORDS(2) <= YEDGE(I,J+1,L) ) INDICES(2) = J

       ! Locate index corresponding to the upper-right latitude
       IF ( COORDS(4) >  YEDGE(I,J,  L)   .and.          &
            COORDS(4) <= YEDGE(I,J+1,L) ) INDICES(4) = J

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
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: its_a_nested_grid
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
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_grid
!
! !DESCRIPTION: Subroutine INIT\_GRID initializes variables and allocates
!  module arrays.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Init_Grid( I1, I2, J1, J2, L1, L2 )
!
! !INPUT PARAMETERS:
!
    INTEGER, INTENT(IN) :: I1, I2   ! Local CPU lon index bounds
    INTEGER, INTENT(IN) :: J1, J2   ! Local CPU lat index bounds
    INTEGER, INTENT(IN) :: L1, L2   ! Local CPU lev index bounds
!
! !REVISION HISTORY:
!  24 Feb 2012 - R. Yantosca - Initial version, based on grid_mod.F
!  01 Mar 2012 - R. Yantosca - Now define IS_NESTED based on Cpp flags
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
! 
    INTEGER :: AS

    !======================================================================
    ! Initialize module variables
    !======================================================================

    ! First assume that we are doing a global simulation
    IS_NESTED = .FALSE.

    ALLOCATE( XMID( I1:I2, J1:J2, L1:L2 ), STAT=AS )
    IF ( AS /= 0 ) CALL ALLOC_ERR( 'XMID' )
    XMID = 0
    
    ALLOCATE( XEDGE( I1:I2+1, J1:J2, L1:L2 ), STAT=AS )
    IF ( AS /= 0 ) CALL ALLOC_ERR( 'XEDGE' )
    XEDGE = 0d0
    
    ALLOCATE( YMID( I1:I2, J1:J2, L1:L2 ), STAT=AS )
    IF ( AS /= 0 ) CALL ALLOC_ERR( 'YMID' )
    YMID = 0d0
    
    ALLOCATE( YEDGE( I1:I2, J1:J2+1, L1:L2 ), STAT=AS )
    IF ( AS /= 0 ) CALL ALLOC_ERR( 'YEDGE' )
    YEDGE = 0d0

    ALLOCATE( YSIN( I1:I2, J1:J2+1, L1:L2 ), STAT=AS )
    IF ( AS /= 0 ) CALL ALLOC_ERR( 'YSIN' )
    YSIN = 0d0

    ALLOCATE( YMID_R( I1:I2, J1:J2, L1:L2 ), STAT=AS )               
    IF ( AS /= 0 ) CALL ALLOC_ERR( 'YMID_R' )
    YMID_R = 0d0
   
    ALLOCATE( YEDGE_R( I1:I2, J1:J2+1, L1:L2 ), STAT=AS )
    IF ( AS /= 0 ) CALL ALLOC_ERR( 'YEDGE_R' )
    YEDGE_R = 0d0

    ALLOCATE( AREA_M2( I1:I2, J1:J2, L1:L2 ), STAT=AS )
    IF ( AS /= 0 ) CALL ALLOC_ERR( 'AREA_M2' )
    AREA_M2 = 0d0

#if defined( NESTED_CH ) || defined( NESTED_EU ) || defined( NESTED_NA ) || defined( SEAC4RS )

    !======================================================================
    ! Special settings for nested-grid simulations only
    !======================================================================

    ! Denote that this is a nested-grid simulation
    IS_NESTED = .TRUE.
    
    ! Allocate nested-grid window array of lat centers (radians)
    ALLOCATE( YMID_R_W( I1:I2, 0:J2+1, L1:L2 ), STAT=AS ) 
    IF ( AS /= 0 ) CALL ALLOC_ERR( 'YMID_R_W' )
    YMID_R_W = 0d0

#endif

  END SUBROUTINE Init_Grid
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cleanup_grid
!
! !DESCRIPTION: Subroutine CLEANUP\_GRID deallocates all module arrays.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Cleanup_Grid
!
! !REVISION HISTORY:
!  24 Feb 2012 - R. Yantosca - Initial version, based on grid_mod.F
!EOP
!------------------------------------------------------------------------------
!BOC
    IF ( ALLOCATED( XMID       ) ) DEALLOCATE( XMID       )
    IF ( ALLOCATED( XEDGE      ) ) DEALLOCATE( XEDGE      )
    IF ( ALLOCATED( YMID       ) ) DEALLOCATE( YMID       )
    IF ( ALLOCATED( YEDGE      ) ) DEALLOCATE( YEDGE      )
    IF ( ALLOCATED( YSIN       ) ) DEALLOCATE( YSIN       )
    IF ( ALLOCATED( YMID_R     ) ) DEALLOCATE( YMID_R     )
    IF ( ALLOCATED( YMID_R_W   ) ) DEALLOCATE( YMID_R_W   )  
    IF ( ALLOCATED( YEDGE_R    ) ) DEALLOCATE( YEDGE_R    )
    IF ( ALLOCATED( AREA_M2    ) ) DEALLOCATE( AREA_M2    )
    
  END SUBROUTINE Cleanup_Grid
!EOC
END MODULE Grid_Mod
