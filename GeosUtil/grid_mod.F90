!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: grid_mod
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

! !REVISION HISTORY:
!  23 Feb 2012 - R. Yantosca - Initial version, based on grid_mod.F
!  01 Mar 2012 - R. Yantosca - Validated for nested grids
!  03 Apr 2012 - M. Payer    - Added ySin for map_a2a regrid (M. Cooper)
!  04 Dec 2012 - R. Yantosca - Modified for GIGC running in ESMF environment
!  26 Feb 2013 - R. Yantosca - Fixed bug in computation of lons & lats when
!                              connecting GEOS-Chem to the GEOS-5 GCM
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
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
!                  GEOS-Chem Global Chemical Transport Model                  !
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
  SUBROUTINE Compute_Grid( am_I_Root,                          &
                           I1, I2, J1,   J2,   JSP,  JNP,      &
                           L1, L2, DLON, DLAT, I_LO, J_LO, RC )
!
! !USES:
!
    USE GIGC_ErrCode_Mod
!
! !INPUT PARAMETERS:
!
    LOGICAL, INTENT(IN)  :: am_I_Root                     ! Root CPU?

    ! Variables with local CPU indices
    INTEGER, INTENT(IN)  :: I1,  I2                       ! Min lon index
    INTEGER, INTENT(IN)  :: J1,  J2                       ! Local lat indices
    INTEGER, INTENT(IN)  :: JSP, JNP                      ! Polar lat indices
    INTEGER, INTENT(IN)  :: L1,  L2                       ! Local lev indices
    REAL*8,  INTENT(IN)  :: DLON(I2-I1+1,J2-J1+1,L2-L1+1) ! Delta lon [deg]
    REAL*8,  INTENT(IN)  :: DLAT(I2-I1+1,J2-J1+1,L2-L1+1) ! Delta lat [deg]

    ! Variables with global CPU indices
    INTEGER, INTENT(IN)  :: I_LO                          ! Min global lon
    INTEGER, INTENT(IN)  :: J_LO                          ! Min global lat
!
! !OUTPUT PARAMETERS:
!  
    INTEGER, INTENT(OUT) :: RC                            ! Success or failure?
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
!  02 Jul 2013 - R. Yantosca - Now compute lon centers properly for GCAP,
!                              which does not have any half-sized polar boxes
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
! 
    ! Scalars
    INTEGER :: I,     J,     L,        IG, JG
    REAL*8  :: SIN_N, SIN_S, SIN_DIFF, TMP

    ! Arrays
    REAL*8  :: IND_X( I1:I2+1 )
    REAL*8  :: IND_Y( J1:J2+1 )

    !======================================================================
    ! Initialization
    !======================================================================

    ! Index array for longitudes
    DO I = I1, I2+1
       IND_X(I) = ( ( I + I0 - 1 ) * 1d0 ) + ( I_LO - 1 )
    ENDDO

    ! Index array for latitudes
    DO J = J1, J2+1
       IND_Y(J) = ( ( J + J0 - 1 ) * 1d0 ) + ( J_LO - 1 )
    ENDDO

    !======================================================================
    ! Compute longitude and latitude arrays
    !======================================================================
    
    ! Loop over levels
    DO L = L1, L2
       
       !-------------------------------------------------------------------
       ! Longitude center and edge arrays
       !-------------------------------------------------------------------

       ! Loop over local latitudes
       DO J = J1, J2

          ! Loop over local longitudes
          DO I = I1, I2

             ! Longitude centers
             XMID(I,J,L)  = ( DLON(I,J,L) * IND_X(I) ) - 180d0
          
             ! Longitude edges
             XEDGE(I,J,L) = XMID(I,J,L) - ( DLON(I,J,L) * 0.5d0 )

             ! Compute the last longitude edge
             IF ( I == I2 ) THEN
                XEDGE(I+1,J,L) = XEDGE(I,J,L) + DLON(I,J,L)
             ENDIF
             
          ENDDO
       ENDDO

       !-------------------------------------------------------------------
       ! Latitude center and edge arrays
       !-------------------------------------------------------------------

       ! Loop over local latitudes
       DO J = J1, J2

          ! Global latitude index
          JG = J + ( J_LO - 1 )

          ! Loop over local longitudes
          DO I = I1, I2

             !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             !%%%%%%   LATITUDE CENTERS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!             ! Lat centers (degrees)
!             YMID(I,J,L)     = ( DLAT(I,J,L) * IND_Y(J) ) - 90d0
          
#if defined( ESMF_ ) || defined( EXTERNAL_GRID ) || defined( EXTERNAL_FORCING )
             !-------------------------------------------------------------
             !  %%%%% CONNECTING TO GEOS-5 GCM via ESMF INTERFACE %%%%%
             !
             ! Do not define half-sized polar boxes (bmy, 3/21/13)
             !-------------------------------------------------------------
#else

# if defined( GCAP )

             !-------------------------------------------------------------
             !             %%%%% TRADITIONAL GEOS-Chem %%%%%
             !
             ! For the GCAP model, there are no half-size polar boxes.
             ! Compute the latitude centers accordingly.  (bmy, 7/2/13)
             !-------------------------------------------------------------

             ! Lat centers (degrees)
             YMID(I,J,L)     = ( DLAT(I,J,L) * IND_Y(J) ) - 88d0

# else

             !-------------------------------------------------------------
             !             %%%%% TRADITIONAL GEOS-Chem %%%%%
             !
             ! Current practice in the standard GEOS-Chem is to make
             ! the polar grid boxes (for non-GCAP global grids only) be
             ! half the size of other grid boxes.  This lets us make +90
             ! degrees and -90 degrees be the edges of the grid.
             ! (bmy, 7/2/13)
             !-------------------------------------------------------------

             ! Lat centers (degrees)
             YMID(I,J,L)     = ( DLAT(I,J,L) * IND_Y(J) ) - 90d0

             IF ( JG == JSP ) THEN
                YMID(I,J,L)  = -90d0 + ( 0.5d0 * DLAT(I,J,L) )   ! S pole
             ELSE IF ( JG == JNP ) THEN
                YMID(I,J,L)  = +90d0 - ( 0.5d0 * DLAT(I,J,L) )   ! N pole
             ENDIF

# endif
#endif
             ! Lat centers (radians)
             YMID_R(I,J,L)   = ( PI_180 * YMID(I,J,L)  )

#if defined( NESTED_CH ) || defined( NESTED_EU ) || defined( NESTED_NA )

             !-------------------------------------------------------------
             !              %%%%% FOR NESTED GRIDS ONLY %%%%%
             !-------------------------------------------------------------

             ! Lat centers (radians), for nested grid window array
             YMID_R_W(I,J,L) = YMID_R(I,J,L)

             ! Compute YMID_R_W at edges of nested region
             IF ( J == J1 ) THEN
                !YMID_R_W(I,J1-1,1) = YMID_R(I,J1,L) - ( DLAT(I,J1,L) * PI_180 )
                YMID_R_W(I,J-1,1) = YMID_R(I,J,L) - ( DLAT(I,J,L) * PI_180 )
             ELSE IF ( J == J2 ) THEN
                !YMID_R_W(I,J2+1,1) = YMID_R(I,J2,L) + ( DLAT(I,J2,L) * PI_180 )
                YMID_R_W(I,J+1,1) = YMID_R(I,J,L) + ( DLAT(I,J,L) * PI_180 )
             ENDIF
                    
#endif

             !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             !%%%%%%   LATITUDE EDGES   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

             ! Lat edges (degrees and radians)
             YEDGE(I,J,L)    = YMID(I,J,L) - ( DLAT(I,J,L) * 0.5d0 )
            
#if defined( ESMF_ ) || defined( EXTERNAL_GRID ) || defined( EXTERNAL_FORCING )
             !-----------------------------------------------------------
             !  %%%%% CONNECTING TO GEOS-5 GCM via ESMF INTERFACE %%%%%
             !
             ! Do not define half-sized polar boxes (bmy, 3/21/13)
             !-----------------------------------------------------------
#else
             !-----------------------------------------------------------
             !            %%%%% TRADITIONAL GEOS-Chem %%%%%
             !
             ! Current practice in the standard GEOS-Chem is to force
             ! the northern edge of grid boxes along the SOUTH POLE to
             ! be -90 degrees latitude. (bmy, 3/21/13)
             !-----------------------------------------------------------
             IF ( JG == JSP ) THEN
                YEDGE(I,J,L) = -90d0                             ! S pole
             ENDIF
#endif          

             ! Lat edges (radians)
             YEDGE_R(I,J,L)  = ( PI_180  * YEDGE(I,J,L) )

             ! mjc - Compute sine of latitude edges (needed for map_a2a regrid)
             YSIN(I,J,L) = SIN ( YEDGE_R(I,J,L) )
             
          ENDDO
       ENDDO

       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       !%%%%%%   LATITUDE EDGES (the last edge)   %%%%%%%%%%%%%%%%%%%%%%%%%
       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
       ! Test for North Pole
       IF ( J2 == JNP ) THEN
          
          ! North pole case (global grids only)
          DO I = I1, I2
#if defined( ESMF_ ) || defined( EXTERNAL_GRID ) || defined( EXTERNAL_FORCING )
             !-----------------------------------------------------------
             !  %%%%% CONNECTING TO GEOS-5 GCM via ESMF INTERFACE %%%%%
             !
             ! Do not define half-sized polar boxes (bmy, 3/21/13)
             !-----------------------------------------------------------
             YEDGE  (I,J2+1,L)   = YEDGE(I,J2,L)   + DLAT(I,J2,L)
#else
             !-----------------------------------------------------------
             !            %%%%% TRADITIONAL GEOS-Chem %%%%%
             !
             ! Current practice in the standard GEOS-Chem is to force
             ! the northern edge of grid boxes along the NORTH POLE to
             ! be +90 degrees latitude. (bmy, 3/21/13)
             !-----------------------------------------------------------
             YEDGE  (I,J2+1,L)   = +90d0
#endif
             YEDGE_R(I,J2+1,L)   = YEDGE(I,J2+1,L) * PI_180

             ! Also compute sine of last latitude edge! (ckeller, 02/13/12)
             YSIN(I,J2+1,L) = SIN ( YEDGE_R(I,J2+1,L) )
          ENDDO
          
       ELSE
          
         !---------------------------------------------------------------
         !                %%%%% FOR NESTED GRIDS ONLY %%%%%
         !---------------------------------------------------------------

          ! No north pole (nested grids only)
          DO I = I1, I2
             YEDGE  (I,J2+1,L)  = YEDGE(I,J2,L  ) + DLAT(I,J2,L)
             YEDGE_R(I,J2+1,L)  = YEDGE(I,J2+1,L) * PI_180

             ! Also compute sine of last latitude edge! (ckeller, 02/13/12)
             YSIN(I,J2+1,L) = SIN ( YEDGE_R(I,J2+1,L) )
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

       ! Loop over local latitudes
       DO J = J1, J2

          ! Global latitude index
          JG  = J + ( J_LO - 1 )
          
          ! Loop over local longitudes
          DO I = I1, I2

             ! Sine of latitudes at N and S edges of grid box (I,J,L)
             SIN_N       = SIN( YEDGE_R(I,J+1,L) )
             SIN_S       = SIN( YEDGE_R(I,J,  L) )

             ! Difference of sin(latitude) at N and S edges of grid box
             SIN_DIFF    = SIN_N - SIN_S

#if defined( ESMF_ ) || defined( EXTERNAL_GRID ) || defined( EXTERNAL_FORCING )
             !-----------------------------------------------------------
             !  %%%%% CONNECTING TO GEOS-5 GCM via ESMF INTERFACE %%%%%
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
             !-----------------------------------------------------------

             ! South pole kludge
             IF ( JG == JSP ) THEN
                SIN_DIFF = 2d0 * ( SIN_N - ( -1d0 ) )
             ENDIF

             ! North pole kludge
             IF ( JG == JNP ) THEN
                SIN_DIFF = 2d0 * ( 1d0 - SIN_S )
             ENDIF
#endif

             ! Grid box surface areas [m2]
             AREA_M2(I,J,L) = ( DLON(I,J,L) * PI_180 ) * ( Re**2 ) * SIN_DIFF

          ENDDO
       ENDDO

    ENDDO

    ! Return successfully
    RC = GIGC_SUCCESS

    !=================================================================
    ! Echo info to stdout
    !=================================================================
    IF ( am_I_Root ) THEN
       WRITE( 6, '(''Nested-Grid X-offset [boxes]:'', i4 )' ) I0
       WRITE( 6, '(''Nested-Grid Y-offset [boxes]:'', i4 )' ) J0
       WRITE( 6, '(a)' )
       WRITE( 6, '(''Grid box longitude centers [degrees]: '')' )
       WRITE( 6, '(8(f8.3,1x))' ) ( XMID(I,1,1),  I=1,I2-I1+1 )
       WRITE( 6, '(a)' )
       WRITE( 6, '(''Grid box longitude edges [degrees]: '')' )
       WRITE( 6, '(8(f8.3,1x))' ) ( XEDGE(I,1,1), I=1,I2+1-I1+1 )
       WRITE( 6, '(a)' )
       WRITE( 6, '(''Grid box latitude centers [degrees]: '')' )
       WRITE( 6, '(8(f8.3,1x))' ) ( YMID(1,J,1),  J=1,J2-J1+1 )
       WRITE( 6, '(a)' )
       WRITE( 6, '(''Grid box latitude edges [degrees]: '')' )
       WRITE( 6, '(8(f8.3,1x))' ) ( YEDGE(1,J,1), J=1,J2+1-J1+1 )
    ENDIF

  END SUBROUTINE Compute_Grid
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
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
!                  GEOS-Chem Global Chemical Transport Model                  !
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
!                  GEOS-Chem Global Chemical Transport Model                  !
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
!                  GEOS-Chem Global Chemical Transport Model                  !
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
!                  GEOS-Chem Global Chemical Transport Model                  !
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
!                  GEOS-Chem Global Chemical Transport Model                  !
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
!                  GEOS-Chem Global Chemical Transport Model                  !
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
!                  GEOS-Chem Global Chemical Transport Model                  !
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
!                  GEOS-Chem Global Chemical Transport Model                  !
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
!  26 Sep 2013 - R. Yantosca - Removed SEAC4RS C-preprocessor switch
!EOP
!------------------------------------------------------------------------------
!BOC

#if defined( NESTED_CH ) || defined( NESTED_EU ) || defined( NESTED_NA )

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
!                  GEOS-Chem Global Chemical Transport Model                  !
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
!                  GEOS-Chem Global Chemical Transport Model                  !
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
!                  GEOS-Chem Global Chemical Transport Model                  !
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
!                  GEOS-Chem Global Chemical Transport Model                  !
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
!                  GEOS-Chem Global Chemical Transport Model                  !
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
!                  GEOS-Chem Global Chemical Transport Model                  !
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
!                  GEOS-Chem Global Chemical Transport Model                  !
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
  SUBROUTINE Init_Grid( am_I_Root, IM, JM, LM, RC )
!
! !USES:
!
    USE GIGC_ErrCode_Mod
!
! !INPUT PARAMETERS:
!
    LOGICAL, INTENT(IN)  :: am_I_Root   ! Are we on the root CPU 
    INTEGER, INTENT(IN)  :: IM          ! # of longitudes on this CPU
    INTEGER, INTENT(IN)  :: JM          ! # of latitudes  on this CPU
    INTEGER, INTENT(IN)  :: LM          ! # of levels     on this CPU
!
! !OUTPUT PARAMETERS:
!
    INTEGER, INTENT(OUT) :: RC          ! Success or failure?
!
! !REVISION HISTORY:
!  24 Feb 2012 - R. Yantosca - Initial version, based on grid_mod.F
!  01 Mar 2012 - R. Yantosca - Now define IS_NESTED based on Cpp flags
!  03 Dec 2012 - R. Yantosca - Add am_I_Root, RC to argument list
!  04 Dec 2012 - R. Yantosca - Now dimension arrays with IM, JM, LM instead 
!                              of I1, J1, L1, I2, J2, L2.  
!  18 Apr 2013 - R. Yantosca - Bug fix: in nested block, use STAT=RC
!  26 Sep 2013 - R. Yantosca - Removed SEAC4RS C-preprocessor switch
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

    ! Assume success
    RC        = GIGC_SUCCESS

    ! First assume that we are doing a global simulation
    IS_NESTED = .FALSE.

    ALLOCATE( XMID   ( IM,   JM,   LM ), STAT=RC )
    IF ( RC /= 0 ) CALL ALLOC_ERR( 'XMID' )
    XMID    = 0d0
    
    ALLOCATE( XEDGE  ( IM+1, JM,   LM ), STAT=RC )
    IF ( RC /= 0 ) CALL ALLOC_ERR( 'XEDGE' )
    XEDGE   = 0d0
    
    ALLOCATE( YMID   ( IM,   JM,   LM ), STAT=RC )
    IF ( RC /= 0 ) CALL ALLOC_ERR( 'YMID' )
    YMID    = 0d0
    
    ALLOCATE( YEDGE  ( IM,   JM+1, LM ), STAT=RC )
    IF ( RC /= 0 ) CALL ALLOC_ERR( 'YEDGE' )
    YEDGE   = 0d0

    ALLOCATE( YSIN   ( IM,   JM+1, LM ), STAT=RC )
    IF ( RC /= 0 ) CALL ALLOC_ERR( 'YSIN' )
    YSIN    = 0d0

    ALLOCATE( YMID_R ( IM,   JM,   LM ), STAT=RC )               
    IF ( RC /= 0 ) CALL ALLOC_ERR( 'YMID_R' )
    YMID_R  = 0d0
   
    ALLOCATE( YEDGE_R( IM,   JM+1, LM ), STAT=RC )
    IF ( RC /= 0 ) CALL ALLOC_ERR( 'YEDGE_R' )
    YEDGE_R = 0d0

    ALLOCATE( AREA_M2( IM,   JM,   LM ), STAT=RC )
    IF ( RC /= 0 ) CALL ALLOC_ERR( 'AREA_M2' )
    AREA_M2 = 0d0

#if defined( NESTED_CH ) || defined( NESTED_EU ) || defined( NESTED_NA )

    !======================================================================
    ! Special settings for nested-grid simulations only
    !======================================================================

    ! Denote that this is a nested-grid simulation
    IS_NESTED = .TRUE.
    
    ! Allocate nested-grid window array of lat centers (radians)
    ALLOCATE( YMID_R_W( IM, 0:JM+1, LM ), STAT=RC ) 
    IF ( RC /= 0 ) CALL ALLOC_ERR( 'YMID_R_W' )
    YMID_R_W = 0d0

#endif

  END SUBROUTINE Init_Grid
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
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
